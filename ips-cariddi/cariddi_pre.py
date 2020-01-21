from netCDF4 import *
import numpy
from copy import *
import argparse
import json

#-------------------------------------------------------------------------------
#  Extract profile parameters for the time index.
#
#  Profile data is written into a json dictionary file that can be passed to
#  ips-vmec to update namelist input files.
#
#  param[in] percorso   Precurser path for the stored profile quantities.
#  param[in] time_index Current time index.
#  param[in] vmec_input JSON file to write vmec profile parameters to.
#  param[in] verbose    Flag to specifiy verbose output.
#-------------------------------------------------------------------------------
def set_namelist_parameters(percorso, time_index, vmec_input, verbose):
#  File with profiles from Cariddi.
    file_prof = '{}/SURFACEQUANT_{:04d}.dat'.format(percorso, time_index)

#  Read profile data.
    data = numpy.loadtxt(file_prof, skiprows = 1)

    s_pol = data[:,0]
    p = data[:,1]
    phi_tot = data[:,4]
    I_tor = data[:,5]

    s = phi_tot/numpy.amax(phi_tot)
    ss = numpy.arange(11)/10.

    II_tor = numpy.interp(ss, s, I_tor)
    pp = numpy.interp(ss, s, p)

    profile_changes = {}

    profile_changes['pmass_type'] = 'akima_spline'
    profile_changes['pcurr_type'] = 'akima_spline_I'
    for i, s in enumerate(ss):
        profile_changes['ac_aux_s({})'.format(i)] = s
        profile_changes['am_aux_s({})'.format(i)] = s
    for i, c in enumerate(II_tor):
        profile_changes['ac_aux_f({})'.format(i)] = c
    for i, p in enumerate(pp):
        profile_changes['am_aux_f({})'.format(i)] = p

    profile_changes['phiedge'] = phi_tot[-1]
    profile_changes['curtor'] = I_tor[-1]

    with open(vmec_input, 'w') as json_ref:
        json.dump(profile_changes, json_ref, indent=4)

    if verbose:
        print('VMEC profile profile parameters written to {}.'.format(vmec_input))

#-------------------------------------------------------------------------------
#  Update carriddi matricies.
#
#  Routine to read cariddi matrices, re-arrange data (K1 and C1 matrices) so
#  that they fit to VMEC ordering of points onto the coupling surface. Also
#  compute quantitites required to compute B on VMEC vacuum grid. The following
#  data are stored in eddy.nc file to be used by VMEC
#       1) Read J,surface data
#       2) save K1,K2,C1,C2,c,s
#
#  param[in] percorso     Precurser path for the stored profile quantities.
#  param[in] time_index   Current time index.
#  param[in] vmec_current VMEC J0 current file.
#  param[in] verbose      Flag to specifiy verbose output.
#-------------------------------------------------------------------------------
    def set_matrix_file(percorso, time_index, vmec_current, verbose):
#  Matrix file for VMEC (read and write)
        matrix_nc_file = '{}/eddy.nc'.format(percorso)

        if time_index == 0:
#  Read matrices.
            data = {}
            for file in ['cs_C1.nc',
                         'cs_C2.nc',
                         'cs_grid.nc',
                         'cs_K1.nc',
                         'cs_K2.nc',
                         'cs_xJ.nc']:
                if verbose:
                    print('reading {}/{} ...'.format(percorso, file))

#  Read all variables.
                src = Dataset('{}/{}'.format(percorso, file))
                for name, variable in src.variables.items():
                    data[name] = numpy.asarray(src.variables[name][:])
                    data[name] = data[name].squeeze()
                src.close()

#  Need to use the proper ordering of indexes for the matrix multiplication.
            inds = copy(data['is'])
            inds = numpy.asarray( [int(x - 1) for x in inds] )
            indk = numpy.zeros_like(inds)
# FIXME: Size of this matrix is hard coded.
#  Cerchiamo gli indici giusti per K.
            for i in range(500):
                indk[i] = numpy.where(inds == (i))[0][0]

# FIXME: Size of this matrix is hard coded.
            I0A = numpy.zeros(484) #(484)  #eddy currents are ZERO at time=0

#  Read matrices and transpose for python computation.
            K1x = copy(data['K1x']).transpose()   #(10201,1500)
            K1y = copy(data['K1y']).transpose()   #(10201,1500)
            K1z = copy(data['K1z']).transpose()   #(10201,1500)
            K2xA = copy(data['K2x']).transpose()  #(10201,484)
            K2yA = copy(data['K2y']).transpose()  #(10201,484)
            K2zA = copy(data['K2z']).transpose()  #(10201,484)
            C1 = copy(data['C1']).transpose()     #(484,1500)
            C2A = copy(data['C2']).transpose()    #(484,484)

#  Re-arrange columns (1500) data:
            K1xA = numpy.zeros_like(K1x)
            K1yA = numpy.zeros_like(K1y)
            K1zA = numpy.zeros_like(K1z)
            C1A = numpy.zeros_like(C1)

# FIXME: Size of this matrix is hard coded.
            inds3 = numpy.append(inds, 500 + inds)
            inds3 = numpy.append(inds3, 1000 + inds)
            inds3 = numpy.asarray( [int(x) for x in inds3] )

# FIXME: Size of this matrix is hard coded.
            indk3 = numpy.append(indk, 500 + indk)
            indk3 = numpy.append(indk3, 1000 + indk)
            indk3 = numpy.asarray( [int(x) for x in indk3] )

            nl = numpy.shape(K1xA)[0]
            for n in range(nl):
                K1xA[n,:] = copy(K1x[n, indk3])
                K1yA[n,:] = copy(K1y[n, indk3])
                K1zA[n,:] = copy(K1z[n, indk3])

            for n in range(numpy.shape(C1)[0]):
                C1A[n,:] = copy(C1[n, indk3])

#  Times after t=0, all matrix data is stored in the matrix_nc_file.
        else:
            if verbose:
                print('reading {} ...'.format(matrix_nc_file))

#  Read all data form the matrix file.
            data = {}
            src = Dataset(matrix_nc_file)
            for name, variable in src.variables.items():
                data[name] = numpy.asarray( src.variables[name][:] )
                data[name] = data[name].squeeze()
            src.close()

#  Read matrices and transpose for python computation.
            K1xA = copy(data['K1x']) #(10201,1500)
            K1yA = copy(data['K1y']) #(10201,1500)
            K1zA = copy(data['K1z']) #(10201,1500)
            K2xA = copy(data['K2x']) #(10201,484)
            K2yA = copy(data['K2y']) #(10201,484)
            K2zA = copy(data['K2z']) #(10201,484)
            C1A = copy(data['C1'])   #(484,1500)
            C2A = copy(data['C2'])   #(484,484)

#  Load current eddy surface currents J0.
        if verbose:
            print('reading {} ...'.format(vmec_current))

        currents0 = {}
        src = Dataset(vmec_current)
        for name, variable in src.variables.items():
            currents0[name] = numpy.asarray(src.variables[name][:])
            currents0[name] = currents0[name].squeeze()
        src.close()

        if verbose:
            print('\nre-ordering matrix and computing vectors for VMEC ...\n')
        J0A = numpy.zeros(1500) #(1500)  # FIXME: Size of this matrix is hard coded.
        for i in range(500):
            J0A[i]        = currents0['k_x'][i]
            J0A[i + 500]  = currents0['k_y'][i]
            J0A[i + 1000] = currents0['k_z'][i]

        if time_index != 0:
            I0A = C1A.dot(J0A) + C1A

#  Compute required quantities: c1,d1 from J0 and I0.
        c1A = -C1A.dot(J0A) + C2A.dot(I0A)     #(484)
        d1xA = -K1xA.dot(J0A) + K2xA.dot(I0A)  #(10201)
        d1yA = -K1yA.dot(J0A) + K2yA.dot(I0A)  #(10201)
        d1zA = -K1zA.dot(J0A) + K2zA.dot(I0A)  #(10201)

#  Update matrix file. If this is the first time step, create the file.
        if time_index == 0:      #things that are to be written only at time=0
            if verbose:
                print('creating and writing EDDY netcdf file {} ...'.format(matrix_nc_file))

            rootgrp = Dataset(matrix_nc_file, 'w', format='NETCDF3_CLASSIC')

            rootgrp.description = 'Cariddi matrix re-alligned and vectors for VMEC run time index {:04d}'.format(time_index)
            rootgrp.history = 'Created {}'.format(time.ctime(time.time()))
            rootgrp.source = 'netCDF3 python module'

            Ng = rootgrp.createDimension('Ng', numpy.shape(K1xA)[0]) # Number of points on vacuum grid.
            Ns = rootgrp.createDimension('Ns', numpy.shape(K1xA)[1]) # 3 times the number of points on control surface
            Nw = rootgrp.createDimension('Nw', numpy.shape(C2A)[1])  # Number of points on eddy current.

#-------------------------------------------------------------------------------
            K1x = rootgrp.createVariable('K1x', 'f8', ('Ng','Ns'))
            K1y = rootgrp.createVariable('K1y', 'f8', ('Ng','Ns'))
            K1z = rootgrp.createVariable('K1z', 'f8', ('Ng','Ns'))

            shape = (numpy.shape(K1xA)[0], numpy.shape(K1xA)[1])

            K1x[:,:] = K1xA.reshape(shape, order='F') # Use Fortran ordering in reshaping
            K1y[:,:] = K1yA.reshape(shape, order='F') # Use Fortran ordering in reshaping
            K1z[:,:] = K1zA.reshape(shape, order='F') # Use Fortran ordering in reshaping
#-------------------------------------------------------------------------------
            K2x = rootgrp.createVariable('K2x', 'f8', ('Ng', 'Nw'))
            K2y = rootgrp.createVariable('K2y', 'f8', ('Ng', 'Nw'))
            K2z = rootgrp.createVariable('K2z', 'f8', ('Ng', 'Nw'))

            shape = (numpy.shape(K2xA)[0], numpy.shape(K2xA)[1])

            K2x[:,:] = K2xA.reshape(shape,order='F') # Use Fortran ordering in reshaping
            K2y[:,:] = K2yA.reshape(shape,order='F') # Use Fortran ordering in reshaping
            K2z[:,:] = K2zA.reshape(shape,order='F') # Use Fortran ordering in reshaping
#-------------------------------------------------------------------------------
            C1 = rootgrp.createVariable('C1', 'f8', ('Nw', 'Ns'))

            C1[:,:] = C1A.reshape((numpy.shape(C1A)[0],
                                   numpy.shape(C1A)[1]),order='F') # Use Fortran ordering in reshaping
#-------------------------------------------------------------------------------
            C2 = rootgrp.createVariable('C2', 'f8', ('Nw', 'Nw'))

            C2[:,:] = C2A.reshape((numpy.shape(C2A)[0],
                                   numpy.shape(C2A)[1]), order='F') #Use Fortran ordering in reshaping

#-------------------------------------------------------------------------------
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

#  Things that are always to be updated in each pre-processing: d and c
        else:
            if verbose:
                print('updating EDDY netcdf file {} ...'.format(matrix_nc_file))

            rootgrp = Dataset(matrix_nc_file, 'a', format='NETCDF3_CLASSIC')

            if verbose:
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

        if verbose:
            print('Data for time index {} stored in:'.format(time_index))
            print('Eddy current matricies written to {}.'.format(matrix_nc_file))

#-------------------------------------------------------------------------------
#  Main routine.
#
#  This routine runs a number of different tasks depending on the command line
#  arguments set.
#
#  param[in] args Dictionary of commandline arguments.
#-------------------------------------------------------------------------------
def main(**args):
    if args['verbose']:
        print('PRE-PROCESSING time index {:04d}\n'.format(args['time_index']))
        print('reading data from files in {}\n'.format(args['matrix_path']))

#  If the vmec_input flag was used write out the new VMEC profiles.
    if 'vmec_input' in args:
        set_namelist_parameters(args['matrix_path'],
                                args['time_index'],
                                args['vmec_input'],
                                args['verbose'])

#  If the vmec_currents update the matrix file.
    if 'vmec_current' in args:
        set_matrix_file(args['matrix_path'],
                        args['time_index'],
                        args['vmec_current'],
                        args['verbose'])

#-------------------------------------------------------------------------------
#  Parse commandline arguments.
#-------------------------------------------------------------------------------
if __name__ == '__main__':
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
                        required=False,
                        action='store',
                        dest='vmec_current',
                        help='Path to the file containing the vmec current',
                        metavar='VMEC_CURRENT')
    parser.add_argument('-i',
                        '--vmec_input',
                        required=False,
                        action='store',
                        dest='vmec_input',
                        help='Path to the file containing the vmec namelist input.',
                        metavar='VMEC_INPUT')
    parser.add_argument('-o',
                        '--verbose_output',
                        required=False,
                        action='store_true',
                        default=False,
                        dest='verbose',
                        help='Flag to control verbose output.')

    args = vars(parser.parse_args())

    #  Remove empty arguments
    for key in [key for key in args if args[key] == None]:
        del args[key]

    main(**args)
