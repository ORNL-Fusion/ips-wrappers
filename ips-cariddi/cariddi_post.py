from netCDF4 import *
import numpy
#from scipy.interpolate import griddata
from copy import *
#from matplotlib.pyplot import *
#import time
import argparse

if __name__ == "__main__":
    #===============================================================================================
    # POST-PROCESSING after a VMEC run
    #===============================================================================================
    # routine to read cariddi matrices (already stored in eddy.nc and re-alligned) and
    # current density on control surface coming out of the VMEC run
    #  1) save J and computed I

    #all data are read and written in the directory >percorso<

    #input data:
    #   none

    #required files:
    # current: J1.nc
    # global file: eddy.nc

    #output data:
    #   eddy.nc     (update)

    #===============================================================================================

    parser = argparse.ArgumentParser(description='POST-processing routine')
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
    args = vars(parser.parse_args())

    #  Remove empty arguments
    for key in [key for key in args if args[key] == None]:
        del args[key]

    verbose = 0   #set to 1 to print more output on the screen

    #========================================================
    percorso = args['matrix_path']

    #matrix file for VMEC (read and write)
    matrix_nc_file = '{}eddy.nc'.format(precorso)    #for updating
    
    #current file produced by VMEC (read)
    fileJ = args['vmec_current']    #these contain Jx,Jy,Jz from VMEC

    #========================================================
    
    print('POST-PROCESSING DATA\n')
    
    #read all data from eddy.nc (matrices and vectors)
    data={}
    print('reading {} ...'.format(matrix_nc_file))
    src = Dataset(matrix_nc_file)
    # read all variables
    for name, variable in src.variables.items():
        data[name] = numpy.asarray( src.variables[name][:] )
        data[name] = data[name].squeeze()
        if (verbose == 1):
            print('  -> variable = {} - shape = {}'.format(name, data[name].shape))

    src.close()

    #========================================================
    #read matrix and transpose for python computation
    C1A = copy(data['C1'])   #(484,1500)
    #read vector
    c1A = copy(data['c'])
    #========================================================

    #read J from VMEC output file
    print('reading {} ...'.format(fileJ))
    src = Dataset(fileJ)
    # read all variables
    currents1 = {}
    for name, variable in src.variables.items():
        currents1[name] = numpy.asarray( src.variables[name][:] )
        currents1[name] = currents1[name].squeeze()
        if (verbose == 1):
            print('  -> variable = {} - shape = '.format(name, currents1[name].shape))

    src.close()

    J1A = np.zeros(1500) #(1500)
    for i in range(500):
         J1A[i] = currents1['k_x'][i]
         J1A[i + 500] = currents1['k_y'][i]
         J1A[i + 1000] = currents1['k_z'][i]

    #compute required quantities: I1,from c and J
    I1A = C1A.dot(J1A) + c1A            #(484) NOT USED HERE but required for next time step: needs J1 form thie time step VMEC converged equilibrium
    #========================================================

    #update matrix_nc_file with:
    #1) I
    #2) J

    print('updating EDDY netcdf file {} ...'.format(matrix_nc_file))

    rootgrp = Dataset(matrix_nc_file, 'a', format='NETCDF3_CLASSIC')

    #update time_index in eddy.nc
    time_index = int(rootgrp.description.split()[-1])
    rootgrp.description = 'Cariddi matrix re-alligned and vectors for VMEC run time index {:04d}'.format(time_index + 1)

    I = rootgrp.variables['I']
    I[:] = I1A

    J = rootgrp.variables['J']
    J[:] = J1A

    rootgrp.close()
