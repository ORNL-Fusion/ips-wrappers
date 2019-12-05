import netCDF4
import numpy
from scipy.interpolate import griddata
from copy import *
import argparse


#routine to read cariddi matrices and write the field into mgrid file.
command_line_parser = argparse.ArgumentParser()
command_line_parser.add_argument('-m',
                                 '--mgrid_file',
                                 required=True,
                                 action='store',
                                 dest='mgrid_file',
                                 help='Path to the mgrid file.',
                                 metavar='MGIRD_FILE')
command_line_parser.add_argument('-p',
                                 '--matrix_path',
                                 required=True,
                                 action='store',
                                 dest='matrix_path',
                                 help='Path to the folder containing the matrices.',
                                 metavar='MATRIX_PATH')
command_line_parser.add_argument('-v',
                                 '--vmec_current',
                                 action='store',
                                 dest='vmec_current',
                                 help='Path to the file containing the vmec current',
                                 metavar='VMEC_CURRENT')
args = vars(command_line_parser.parse_args())

#  Remove empty arguments
for key in [key for key in args if args[key] == None]:
    del args[key]

write_mgrid = 1   #to write data into mgrid file
mgridfile = args['mgrid_file']

percorso = args['matrix_path']
filesIn=['cs_K1.nc','cs_xJ.nc','cs_grid.nc']
filesJ=['J0.nc']

print('reading data from fields in '+percorso)

if 'vmec_current' in args:
    data={}
    for fileIn in filesIn:
        print('reading '+percorso+'/'+fileIn+' ...')
        src=netCDF4.Dataset(percorso+'/'+fileIn)
        # read all variables
        for name, variable in src.variables.items():
            #data[name] = src.variables[name][:]
            data[name] = numpy.asarray( src.variables[name][:] )
            data[name] = data[name].squeeze()
            print('  -> variable = '+name+' - shape = '+str(data[name].shape))

        src.close()

    cnt=0
    currents={}
    for fileIn in filesJ:
        print('reading '+percorso+'/'+fileIn+' ...')
        src=netCDF4.Dataset(percorso+'/'+fileIn)
        # read all variables
        for name, variable in src.variables.items():
            namecnt=name+str(cnt)
            #currents[name] = src.variables[name][:]
            currents[namecnt] = numpy.asarray( src.variables[name][:] )
            currents[namecnt] = currents[namecnt].squeeze()
            print('  -> variable = '+namecnt+' - shape = '+str(currents[namecnt].shape))

        src.close()
        cnt+=1

    print('reading '+args['vmec_current']+' ...')
    src=netCDF4.Dataset(args['vmec_current'])
    # read all variables
    for name, variable in src.variables.items():
        namecnt=name+str(cnt)
        currents[namecnt] = numpy.asarray( src.variables[name][:] )
        currents[namecnt] = currents[namecnt].squeeze()
        print('  -> variable = '+namecnt+' - shape = '+str(currents[namecnt].shape))
    src.close()

    #copiamo i dati della griglia
    data['xs']=copy(currents['x0'])
    data['ys']=copy(currents['y0'])
    data['zs']=copy(currents['z0'])

    print('matrices in dictionary: data')
    print(data.keys())
    print('currents in dictionary: currents')
    print(currents.keys())

    #compute magnetic field on vaccum grid for VMEC

    print('\ncomputing magnetic field on vaccum grid for VMEC ...\n')

    #need to use the proper ordering of indexes for the matrix multiplication
    inds=copy (data['is'])
    inds=numpy.asarray( [int(x-1) for x in inds] )
    k_vmec = numpy.zeros(1500) #(1500)
    for i in range(500):
        k_vmec[i]      = currents['k_x1'][inds[i]]-currents['k_x0'][inds[i]]
        k_vmec[i+500]  = currents['k_y1'][inds[i]]-currents['k_y0'][inds[i]]
        k_vmec[i+1000] = currents['k_z1'][inds[i]]-currents['k_z0'][inds[i]]

    K1x = copy(data['K1x']).transpose()   #(10201,1500)
    K1y = copy(data['K1y']).transpose()   #(10201,1500)
    K1z = copy(data['K1z']).transpose()   #(10201,1500)

    Bx = K1x.dot(k_vmec)
    By = K1y.dot(k_vmec)
    Bz = K1z.dot(k_vmec)

    x=copy(data['x_grid'])
    y=copy(data['y_grid'])
    z=copy(data['z_grid'])
    points = numpy.array([x,z]).transpose()
    gx,gz  = numpy.mgrid[numpy.amin(x):numpy.amax(x):500j,numpy.amin(z):numpy.amax(z):500j]
    gBx=griddata(points,Bx,(gx,gz),method='cubic')
    gBy=griddata(points,By,(gx,gz),method='cubic')
    gBz=griddata(points,Bz,(gx,gz),method='cubic')

if (write_mgrid==1):
    print('writing data in '+mgridfile)
    #write mgrid data
    print('mgrid file: {}'.format(mgridfile))
    print('reshaping data for mgrid file ...')

    mgriddata = netCDF4.Dataset(mgridfile, "a")
    nr = numpy.array( mgriddata.variables['ir'] )
    nz = numpy.array( mgriddata.variables['jz'] )
    nphi = numpy.array( mgriddata.variables['kp'] )
    nextcur = numpy.array( mgriddata.variables['nextcur'] )

    #names of variables to write
    br_str='br_'+(str(nextcur).zfill(3))
    bp_str='bp_'+(str(nextcur).zfill(3))
    bz_str='bz_'+(str(nextcur).zfill(3))

    br_mgrid=mgriddata.variables[br_str]
    bp_mgrid=mgriddata.variables[bp_str]
    bz_mgrid=mgriddata.variables[bz_str]

    if 'vmec_current' in args:
        br=Bx.reshape((nr,nz),order='F')    #use Fortran ordering in reshaping
        bp=By.reshape((nr,nz),order='F')    #use Fortran ordering in reshaping
        bz=Bz.reshape((nr,nz),order='F')    #use Fortran ordering in reshaping

        print('writing data into mgrid file ...')

        for k in range(nphi):
            br_mgrid[k,:,:]=br
            bp_mgrid[k,:,:]=bp
            bz_mgrid[k,:,:]=bz
    else:
        br_mgrid[:,:,:] = 0
        bp_mgrid[:,:,:] = 0
        bz_mgrid[:,:,:] = 0

    mgriddata.close()

    print('DONE!')
