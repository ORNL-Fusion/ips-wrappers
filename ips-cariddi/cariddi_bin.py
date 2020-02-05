import netCDF4
import numpy
import argparse

#-------------------------------------------------------------------------------
#  Class to represent an carriddi file.
#-------------------------------------------------------------------------------
class carriddi_file:
#*******************************************************************************
#  CONSTRUCTORS
#*******************************************************************************
#-------------------------------------------------------------------------------
#  Initialize a carriddi file.
#
#  param[inout] self      A carriddi file instance.
#  param[in]    file_name Path to the file.
#-------------------------------------------------------------------------------
    def __init__(self, file_name):
        self.data = netCDF4.Dataset(file_name, 'a')

#*******************************************************************************
#  CONSTRUCTORS
#*******************************************************************************
#-------------------------------------------------------------------------------
#  Destruct a cariddi file instance.
#
#  param[inout] self A mgrid file instance.
#-------------------------------------------------------------------------------
    def __del__(self):
        self.data.close()

#-------------------------------------------------------------------------------
#  Class to represent an eddy current file.
#-------------------------------------------------------------------------------
class eddy(carriddi_file):
#*******************************************************************************
#  CONSTRUCTORS
#*******************************************************************************
    def __init__(self, file_name):
        super(eddy, self).__init__(file_name)

        self.dx = self.data.variables['dx'][:].data
        self.dy = self.data.variables['dy'][:].data
        self.dz = self.data.variables['dz'][:].data

        self.K1x = self.data.variables['K1x'][:].data
        self.K1y = self.data.variables['K1y'][:].data
        self.K1z = self.data.variables['K1z'][:].data

#-------------------------------------------------------------------------------
#  Class to represent an vmec current file.
#-------------------------------------------------------------------------------
class vmec_current(carriddi_file):
#*******************************************************************************
#  CONSTRUCTORS
#*******************************************************************************
    def __init__(self, file_name):
        super(vmec_current, self).__init__(file_name)

        k_x = self.data.variables['k_x'][:].data
        k_y = self.data.variables['k_y'][:].data
        k_z = self.data.variables['k_z'][:].data

        k_vmec = numpy.empty(3*len(k_x))
        k_vmec[0:len(k_x)]            = k_x[:]
        k_vmec[len(k_x):2*len(k_x)]   = k_y[:]
        k_vmec[2*len(k_x):3*len(k_x)] = k_z[:]

#-------------------------------------------------------------------------------
#  Class to represent an mgrid file.
#-------------------------------------------------------------------------------
class mgrid(carriddi_file):
#*******************************************************************************
#  CONSTRUCTORS
#*******************************************************************************
#-------------------------------------------------------------------------------
#  Initialize a mgrid file.
#
#  param[inout] self      A mgrid file instance.
#  param[in]    file_name Path to the mgrid file.
#-------------------------------------------------------------------------------
    def __init__(self, file_name):
        super(mgrid, self).__init__(file_name)

        self.nr = self.data.variables['ir'][:].data
        self.nz = self.data.variables['jz'][:].data
        self.nphi = self.data.variables['kp'][:].data
        nextcur = self.data.variables['nextcur'][:].data

        self.br = self.data.variables['br_{:03d}'.format(nextcur)]
        self.bp = self.data.variables['bp_{:03d}'.format(nextcur)]
        self.bz = self.data.variables['bz_{:03d}'.format(nextcur)]

        self.br[:,:,:] = 0
        self.bp[:,:,:] = 0
        self.bz[:,:,:] = 0

        print('writing data in {}'.format(file_name))
        print('reshaping data for mgrid file ...')

#*******************************************************************************
#  SETTERS
#*******************************************************************************
#-------------------------------------------------------------------------------
#  Set eddy current contribution to the VMEC fields.
#
#  param[inout] self         A mgrid file instance.
#  param[in]    matrix_path  Path to the matrix file.
#  param[in]    vmec_current Path to the vmec current file.
#-------------------------------------------------------------------------------
    def set_fields(self, matrix_path, vmec_current):
        print('\ncomputing magnetic field on vaccum grid for VMEC ...\n')
        eddyfile = eddy('{}/eddy.nc'.format(matrix_path))

        j1file = vmec_current(vmec_current)

#  Need to use the proper ordering of indexes for the matrix multiplication.
        k_vmec = numpy.empty(3*len(j1file.k_x))
        k_vmec[0:len(j1file.k_x)]                   = j1file.k_x[:]
        k_vmec[len(j1file.k_x):2*len(j1file.k_x)]   = j1file.k_y[:]
        k_vmec[2*len(j1file.k_x):3*len(j1file.k_x)] = j1file.k_z[:]

#  Use Fortran ordering in reshaping.
        br = (eddyfile.K1x.dot(k_vmec) + eddyfile.dx).reshape((self.nr,self.nz), order='F')
        bp = (eddyfile.K1y.dot(k_vmec) + eddyfile.dy).reshape((self.nr,self.nz), order='F')
        bz = (eddyfile.K1z.dot(k_vmec) + eddyfile.dz).reshape((self.nr,self.nz), order='F')

        print('writing data into mgrid file ...')

#  FIXME: Assumes axisymmetry?
        for k in range(self.nphi):
            self.br[k,:,:] = br
            self.bp[k,:,:] = bp
            self.bz[k,:,:] = bz

#-------------------------------------------------------------------------------
#  Main routine
#
#  This routine runs a number of different tasks based on the command line
#  arguments set.
#
#  param[in] args Dictionary of commandline arguments.
#-------------------------------------------------------------------------------
def main(**args):
    mgridfile = mgrid(args['mgrid_file'])

    if 'current_path' in args:
        mgridfile.set_fields(args['matrix_path'],
                             args['current_path'])

#-------------------------------------------------------------------------------
#  Parse commandline arguments.
#-------------------------------------------------------------------------------
if __name__ == '__main__':
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
                                     action='store',
                                     dest='matrix_path',
                                     help='Path to the folder containing the matrices.',
                                     metavar='MATRIX_PATH')
    command_line_parser.add_argument('-c',
                                     '--vmec_current',
                                     action='store',
                                     dest='vmec_current',
                                     help='Path to the file containing the vmec current.',
                                     metavar='VMEC_CURRENT')
    args = vars(command_line_parser.parse_args())

    #  Remove empty arguments
    for key in [key for key in args if args[key] == None]:
        del args[key]

    if 'current_path' in args and 'matrix_path' not in args:
        print('Faital error: "matrix_path" must be set when using the "current_path" flag.')
        exit(1)

    main(**args)
