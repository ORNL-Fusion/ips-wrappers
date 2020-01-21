from netCDF4 import *
import numpy
from scipy.interpolate import griddata
from copy import *
import time
import argparse

def derive_vmec(x, y):
   #calcolo la derivata in stile VMEC
   # dx=x[1]-x[0]
   # yp[i]= ( y[i]-y[i-1] ) / dx

   ny=numpy.size(y)
   dx=x[1]-x[0]

   yp=numpy.zeros(ny, numpy.double)

   for i in range(1, ny):
      yp[i] = (y[i] - y[i - 1])/dx

   yp[0] = -1.0*(yp[2] - yp[1]) + yp[1]        #per uniformare il tutto: fit lineare

   return yp

if __name__ == "__main__":
    #===============================================================================================
    # PRE-PROCESSING before starting a VMEC run
    #===============================================================================================
    # routine to read cariddi matrices, re-arrange data (K1 and C1 matrices) so that they fit to
    # VMEC ordering of points onto the coupling surface. Also compute quantitites required to
    # compute B on VMEC vacuum grid.
    # the following data are stored in eddy.nc file to be used by VMEC
    #  1) save provided I,J,surface data
    #  2) save K1,K2,C1,C2,c,s
    # create lines for VMEC input file and store in input_VMEC.profiles_#

    #all data are read and written in the directory >percorso<
    
    #input data:
    #   time_index

    #required files:
    # matrices: cs_C1.nc,cs_C2.nc,cs_grid.nc,cs_K1.nc,cs_K2.nc,cs_xJ.nc
    # profiles from cariddi: SURFQUANT_#.dat
    # current: J1.nc
    # global file: eddy.nc (created at time step 0)
    
    #output data:
    #   eddy.nc                 (update or create)
    #   input_VMEC.profiles_#   (create)
    
    #===============================================================================================

    parser = argparse.ArgumentParser(description='PRE-processing routine')
    parser.add_argument('-t',
                        '--time_index',
                        required=True,
                        action='store',
                        dest='time_index',
                        help='time index (integer)',
                        type=int,
                        default=0,
                        metavar='TIME_INDEX')
    parser.add_argument('-p',
                        '--matrix_path',
                        required=True,
                        action='store',
                        dest='matrix_path',
                        help='Path to the folder containing the matrices.',
                        metavar='MATRIX_PATH')
    parser.add_argument('-v',
                        '--vmec_current',
                        required=True,
                        action='store',
                        dest='vmec_current',
                        help='Path to the file containing the vmec current',
                        metavar='VMEC_CURRENT')
    parser.add_argument('-i',
                        '--vmec_input',
                        required=True,
                        action='store',
                        dest='vmec_input',
                        help='Path to the file containing the vmec namelist input.',
                        metavar='VMEC_INPUT')

    args = vars(parser.parse_args())

    #  Remove empty arguments
    for key in [key for key in args if args[key] == None]:
        del args[key]

    verbose = 0   #set to 1 to print more output on the screen

    #========================================================

    time_index = args['time_index'] #time index to pre-process. If time_index=ZERO: build netcdf eddy.nc file
    
    #========================================================

    #percorso='/rfx/home/terranova/work/15_F4E-OPE0951/matrices/PrePostProc/'
    percorso = args['matrix_path']

    #matrix files from Cariddi (read)
    filesIn = ['cs_C1.nc','cs_C2.nc','cs_grid.nc','cs_K1.nc','cs_K2.nc','cs_xJ.nc']
    
    #matrix file for VMEC (read and write)
    matrix_nc_file = '{}eddy.nc'.format(percorso)
    
    #current file produced by VMEC (read)
    fileJ = args['vmec_current'] #these contain Jx,Jy,Jz from VMEC
    
    #file with profiles from Cariddi (read)
    file_prof = '{}SURFACEQUANT_{:04d}.dat'.format(percorso, time_index)
    
    #input file for VMEC (write)
    vmec_input_file = args['vmec_input']
    
    #========================================================
    
    print('PRE-PROCESSING time index {:04d}\n'.format(time_index))
    print('reading data from files in {}\n'.format(percorso))
    
    if (time_index == 0):
        #read matrices
        data={}
        for fileIn in filesIn:
            if (verbose == 1):
                print('reading {}{} ...'.format(percorso, fileIn))
            src = Dataset('{}{}'.format(percorso, fileIn))

            # read all variables
            for name, variable in src.variables.items():
                data[name] = numpy.asarray(src.variables[name][:])
                data[name] = data[name].squeeze()
                if (verbose == 1):
                    print('  -> variable = {} - shape = {}'.format(name, data[name].shape))

            src.close()

        #READ DATA FROM TIME=0 ITERATION I0(484) AND J0(1500)
        #I and J at time=0 are to be read externally, but for other times, they are updated in
        #the post-processing routine amd stored in the eddy.nc directly.
        print('reading {} ...'.format(fileJ))
        src = Dataset(fileJ)
        # read all variables
        currents0 = {}
        for name, variable in src.variables.items():
            currents0[name] = numpy.asarray( src.variables[name][:] )
            currents0[name] = currents0[name].squeeze()
            if (verbose==1):
                print('  -> variable = {} - shape = {}'.format(name, currents0[name].shape))

        src.close()

        #copiamo i dati della griglia
        data['xs'] = copy(currents0['x'])
        data['ys'] = copy(currents0['y'])
        data['zs'] = copy(currents0['z'])

        print('matrices in dictionary: data')
        print(data.keys())
        print('currents in dictionary: currents')
        print(currents0.keys())

        print('\nre-ordering matrix and computing vectors for VMEC ...\n')

        #need to use the proper ordering of indexes for the matrix multiplication
        inds = copy(data['is'])
        inds = numpy.asarray( [int(x - 1) for x in inds] )
        indk = numpy.zeros_like(inds)
        J0A = numpy.zeros(1500) #(1500)
        for i in range(500):
            J0A[i] = currents0['k_x'][i]
            J0A[i+500] = currents0['k_y'][i]
            J0A[i+1000] = currents0['k_z'][i]
            #cerchiamo gli indici giusti per K
            indk[i] = numpy.where(inds == (i))[0][0]

        I0A = numpy.zeros(484) #(484)  #eddy currents are ZERO at time=0

        #========================================================
        #read matrices and transpose for python computation
        K1x = copy(data['K1x']).transpose()   #(10201,1500)
        K1y = copy(data['K1y']).transpose()   #(10201,1500)
        K1z = copy(data['K1z']).transpose()   #(10201,1500)
        K2xA = copy(data['K2x']).transpose()   #(10201,484)
        K2yA = copy(data['K2y']).transpose()   #(10201,484)
        K2zA = copy(data['K2z']).transpose()   #(10201,484)
        C1 = copy(data['C1']).transpose()   #(484,1500)
        C2A = copy(data['C2']).transpose()   #(484,484)
        #========================================================

        #re-arrange columns (1500) data:
        K1xA = numpy.zeros_like(K1x)
        K1yA = numpy.zeros_like(K1y)
        K1zA = numpy.zeros_like(K1z)
        C1A = numpy.zeros_like(C1)
        
        inds3 = numpy.append(inds,500+inds)
        inds3 = numpy.append(inds3,1000+inds)
        inds3 = numpy.asarray( [int(x) for x in inds3] )

        indk3 = numpy.append(indk,500+indk)
        indk3 = numpy.append(indk3,1000+indk)
        indk3 = numpy.asarray( [int(x) for x in indk3] )

        nl = numpy.shape(K1xA)[0]
        for n in range(nl):
            K1xA[n,:] = copy(K1x[n,indk3])
            K1yA[n,:] = copy(K1y[n,indk3])
            K1zA[n,:] = copy(K1z[n,indk3])

        for n in range(numpy.shape(C1)[0]):
            C1A[n,:] = copy(C1[n,indk3])

    else:   #times after time=0 when data are all stored in >matrix_nc_file<
        #read all data
        data = {}
        print('reading {} ...'.format(matrix_nc_file))
        src = Dataset(matrix_nc_file)
        # read all variables
        for name, variable in src.variables.items():
            data[name] = numpy.asarray( src.variables[name][:] )
            data[name] = data[name].squeeze()
            if (verbose == 1):
                print('  -> variable = {} - shape = '.format(name, data[name].shape))

        src.close()

        #read matrices and transpose for python computation
        K1xA = copy(data['K1x'])   #(10201,1500)
        K1yA = copy(data['K1y'])   #(10201,1500)
        K1zA = copy(data['K1z'])   #(10201,1500)
        K2xA = copy(data['K2x'])   #(10201,484)
        K2yA = copy(data['K2y'])   #(10201,484)
        K2zA = copy(data['K2z'])   #(10201,484)
        C1A = copy(data['C1'])   #(484,1500)
        C2A = copy(data['C2'])   #(484,484)

        J0A = data['J']
        I0A = data['I']

    #========================================================
    #compute required quantities: c1,d1 from J0 and I0

    c1A = -C1A.dot(J0A) + C2A.dot(I0A)     #(484)
    d1xA = -K1xA.dot(J0A) + K2xA.dot(I0A)  #(10201)
    d1yA = -K1yA.dot(J0A) + K2yA.dot(I0A)  #(10201)
    d1zA = -K1zA.dot(J0A) + K2zA.dot(I0A)  #(10201)

    #========================================================
    #write the single netcdf file containing the following elements:
    #1) K1,C1 re-ordered
    #2) K2,C2
    #3) x,y,z control surface
    #4) D vector (=-K1_dot_J0) and c vector
    #5) I0,J0
    #necessary to run VMEC in the starting time step

    if (time_index == 0):      #things that are to be written only at time=0
        print('creating and writing EDDY netcdf file {} ...'.format(matrix_nc_file))

        rootgrp = Dataset(matrix_nc_file, 'w', format='NETCDF3_CLASSIC')

        rootgrp.description = 'Cariddi matrix re-alligned and vectors for VMEC run time index {:04d}'.format(time_index)
        rootgrp.history = 'Created {}'.format(time.ctime(time.time()))
        rootgrp.source = 'netCDF3 python module'

        Ng = rootgrp.createDimension('Ng', numpy.shape(K1xA)[0])  #number of points on vacuum grid
        Ns = rootgrp.createDimension('Ns', numpy.shape(K1xA)[1])  #3 time sthe number of points on control surface
        Ns3 = rootgrp.createDimension('Ns3', numpy.size(data['xs']))    #Number of points on control surface
        Nw = rootgrp.createDimension('Nw', numpy.shape(C2A)[1])        #number of points on eddy current

        #--------------------------------------------------------------------------------
        K1x = rootgrp.createVariable('K1x', 'f8', ('Ng','Ns'))
        K1y = rootgrp.createVariable('K1y', 'f8', ('Ng','Ns'))
        K1z = rootgrp.createVariable('K1z', 'f8', ('Ng','Ns'))

        K1x[:,:] = K1xA.reshape((numpy.shape(K1xA)[0], numpy.shape(K1xA)[1]), order='F')    #use Fortran ordering in reshaping
        K1y[:,:] = K1yA.reshape((numpy.shape(K1xA)[0], numpy.shape(K1xA)[1]), order='F')    #use Fortran ordering in reshaping
        K1z[:,:] = K1zA.reshape((numpy.shape(K1xA)[0], numpy.shape(K1xA)[1]), order='F')    #use Fortran ordering in reshaping
        #--------------------------------------------------------------------------------
        K2x = rootgrp.createVariable('K2x', 'f8', ('Ng', 'Nw'))
        K2y = rootgrp.createVariable('K2y', 'f8', ('Ng', 'Nw'))
        K2z = rootgrp.createVariable('K2z', 'f8', ('Ng', 'Nw'))

        K2x[:,:] = K2xA.reshape((numpy.shape(K2xA)[0],numpy.shape(K2xA)[1]),order='F')    #use Fortran ordering in reshaping
        K2y[:,:] = K2yA.reshape((numpy.shape(K2xA)[0],numpy.shape(K2xA)[1]),order='F')    #use Fortran ordering in reshaping
        K2z[:,:] = K2zA.reshape((numpy.shape(K2xA)[0],numpy.shape(K2xA)[1]),order='F')    #use Fortran ordering in reshaping
        #--------------------------------------------------------------------------------
        C1 = rootgrp.createVariable('C1', 'f8', ('Nw', 'Ns'))

        C1[:,:] = C1A.reshape((numpy.shape(C1A)[0],numpy.shape(C1A)[1]),order='F')    #use Fortran ordering in reshaping
        #--------------------------------------------------------------------------------
        C2 = rootgrp.createVariable('C2', 'f8', ('Nw', 'Nw'))

        C2[:,:] = C2A.reshape((numpy.shape(C2A)[0],numpy.shape(C2A)[1]),order='F')    #use Fortran ordering in reshaping
        #--------------------------------------------------------------------------------
        xS = rootgrp.createVariable('xS', 'f8', ('Ns3'))
        yS = rootgrp.createVariable('yS', 'f8', ('Ns3'))
        zS = rootgrp.createVariable('zS', 'f8', ('Ns3'))

        xS[:] = data['xs']
        yS[:] = data['ys']
        zS[:] = data['zs']
        #--------------------------------------------------------------------------------
        #I and J at time=0 are to be written, but for other times, they are updated in
        #the post-processing routine.
        #--------------------------------------------------------------------------------
        I = rootgrp.createVariable('I', 'f8', ('Nw'))

        I[:] = I0A
        #--------------------------------------------------------------------------------
        J = rootgrp.createVariable('J', 'f8', ('Ns'))

        J[:] = J0A
        #--------------------------------------------------------------------------------
        dx = rootgrp.createVariable('dx', 'f8', ('Ng'))
        dy = rootgrp.createVariable('dy', 'f8', ('Ng'))
        dz = rootgrp.createVariable('dz', 'f8', ('Ng'))

        dx[:] = d1xA
        dy[:] = d1yA
        dz[:] = d1zA
        #--------------------------------------------------------------------------------
        c = rootgrp.createVariable('c', 'f8', ('Nw'))

        c[:] = c1A
        #--------------------------------------------------------------------------------

        rootgrp.close()
    
    else:    #things that are always to be updated in each pre-processing: d and c
        print('updating EDDY netcdf file {} ...'.format(matrix_nc_file))

        rootgrp = Dataset(matrix_nc_file, 'a', format='NETCDF3_CLASSIC')

        print(rootgrp.description)

        c = rootgrp.variables['c']
        c = c1A

        dx = rootgrp.variables['dx']
        dy = rootgrp.variables['dy']
        dz = rootgrp.variables['dz']

        dx[:] = d1xA
        dy[:] = d1yA
        dz[:] = d1zA

        rootgrp.close()

    #========================================================
    #read data to prepare VMEC input file
    #========================================================

    #read profiles data
    data = numpy.loadtxt(file_prof, skiprows = 1)

    s_pol = data[:,0]
    p = data[:,1]
    F = data[:,2]
    phi_pl = data[:,3]
    phi_tot = data[:,4]
    I_tor = data[:,5]
    q = data[:,6]

    s = phi_tot/numpy.amax(phi_tot)
    ss = numpy.arange(11)/10.

    qq = numpy.interp(ss, s, q)

    II_tor = numpy.interp(ss, s, I_tor)

    pp = numpy.interp(ss, s, p)

    x = numpy.arange(50)/49.
    y = numpy.interp(x, s, I_tor)
    dy = derive_vmec(x, y)
    dII_tor = numpy.interp(ss, x, dy)

    print ('\nWriting profiles data for VMVEC inpout file in {}'.format(vmec_input_file))

    prof_file = open(vmec_input_file,'w')
    prof_file.write('! profiles file written by python pre-processing procedure\n')
    prof_file.write('! time index {:04d}\n'.format(time_index))

    prof_file.write("\npiota_type = 'akima_spline'")
    prof_file.write('\nai_aux_s = '+str(ss)[1:-1])
    prof_file.write('\nai_aux_f = '+str(-1./qq)[1:-1])

    prof_file.write("\n\npcurr_type = 'akima_spline_I'")
    prof_file.write('\nac_aux_s = '+str(ss)[1:-1])
    prof_file.write('\nac_aux_f = '+str(II_tor)[1:-1])

    prof_file.write("\n\npcurr_type = 'akima_spline_Ip'")
    prof_file.write('\nac_aux_s = '+str(ss)[1:-1])
    prof_file.write('\nac_aux_f = '+str(dII_tor)[1:-1])

    prof_file.write("\n\npmass_type = 'akima_spline'")
    prof_file.write('\nam_aux_s = '+str(ss)[1:-1])
    prof_file.write('\nam_aux_f = '+str(pp)[1:-1])

    prof_file.write('\n\nPHIEDGE = '+str(phi_tot[-1]))

    prof_file.write('\nCURTOR = '+str(I_tor[-1]))

    prof_file.close()

    print('Data for time index {} stored in:'.format(time_index))
    print('1) eddy currents -> {}'.format(matrix_nc_file))
    print('2) VMEC input profiles -> {}'.format(vmec_input_file))

    exit()
