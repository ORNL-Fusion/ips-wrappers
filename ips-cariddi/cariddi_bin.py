import netCDF4
import numpy
import argparse
import json

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

        self.k1xx = self.data.variables['k1xx'][:].data
        self.k1xy = self.data.variables['k1xy'][:].data
        self.k1xz = self.data.variables['k1xz'][:].data

        self.k1yx = self.data.variables['k1yx'][:].data
        self.k1yy = self.data.variables['k1yy'][:].data
        self.k1yz = self.data.variables['k1yz'][:].data

        self.k1zx = self.data.variables['k1zx'][:].data
        self.k1zy = self.data.variables['k1zy'][:].data
        self.k1zz = self.data.variables['k1zz'][:].data

        self.bx = self.data.variables['bx'][:].data
        self.by = self.data.variables['by'][:].data
        self.bz = self.data.variables['bz'][:].data

        self.ax = self.data.variables['ax']
        self.ay = self.data.variables['ay']
        self.az = self.data.variables['az']

#-------------------------------------------------------------------------------
#  Class to represent an vmec current file.
#-------------------------------------------------------------------------------
class vmec_current(carriddi_file):
#*******************************************************************************
#  CONSTRUCTORS
#*******************************************************************************
    def __init__(self, file_name):
        super(vmec_current, self).__init__(file_name)

        self.ax = self.data.variables['ax'][:].data
        self.ay = self.data.variables['ay'][:].data
        self.az = self.data.variables['az'][:].data

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

        print('writing data in {}'.format(file_name))
        print('reshaping data for mgrid file ...')

#*******************************************************************************
#  SETTERS
#*******************************************************************************
#-------------------------------------------------------------------------------
#  Zero eddy current contribution to the VMEC fields.
#
#  param[inout] self              A mgrid file instance.
#-------------------------------------------------------------------------------
    def set_zero(self):
        self.br[:,:,:] = 0
        self.bp[:,:,:] = 0
        self.bz[:,:,:] = 0

#-------------------------------------------------------------------------------
#  Set eddy current contribution to the VMEC fields.
#
#  param[inout] self              A mgrid file instance.
#  param[in]    matrix_path       Path to the matrix file.
#  param[in]    vmec_current_path Path to the vmec current file.
#-------------------------------------------------------------------------------
    def set_fields(self, matrix_path, vmec_current_path):
        print('\ncomputing magnetic field on vaccum grid for VMEC ...\n')
        eddyfile = eddy('{}/eddy.nc'.format(matrix_path))

        a1file = vmec_current(vmec_current_path)

        shape = (self.nphi, self.nz, self.nr)
        if self.nz*self.nr == numpy.shape(eddyfile.k1xx)[0]:
            shape = (1, self.nr, self.nz)

#  Use Fortran ordering in reshaping. May need to account for the change in
#  coordinate from
        br = numpy.reshape(numpy.matmul(eddyfile.k1xx, a1file.ax) +
                           numpy.matmul(eddyfile.k1xy, a1file.ay) +
                           numpy.matmul(eddyfile.k1xz, a1file.az) +
                           eddyfile.dx, shape, order='F')
        bp = numpy.reshape(numpy.matmul(eddyfile.k1yx, a1file.ax) +
                           numpy.matmul(eddyfile.k1yy, a1file.ay) +
                           numpy.matmul(eddyfile.k1yz, a1file.az) +
                           eddyfile.dy, shape, order='F')
        bz = numpy.reshape(numpy.matmul(eddyfile.k1zx, a1file.ax) +
                           numpy.matmul(eddyfile.k1zy, a1file.ay) +
                           numpy.matmul(eddyfile.k1zz, a1file.az) +
                           eddyfile.dz, shape, order='F')

#  Save vector potential to compute alpha
        eddyfile.ax[:] = a1file.ax
        eddyfile.ay[:] = a1file.ay
        eddyfile.az[:] = a1file.az

        print('writing data into mgrid file ...')

        if shape[0] == 1:
            for k in range(self.nphi):
                self.br[k,:,:] = br
                self.bp[k,:,:] = bp
                self.bz[k,:,:] = bz

#-------------------------------------------------------------------------------
#  Set eddy current contribution to the VMEC fields.
#
#  param[inout] self              A mgrid file instance.
#  param[in]    matrix_path       Path to the matrix file.
#  param[in]    vmec_current_path Path to the vmec current file.
#-------------------------------------------------------------------------------
    def scale_fields(self, matrix_path, vmec_current_path):
        print('\nscaling magnetic field on vaccum grid for VMEC ...\n')
        eddyfile = eddy('{}/eddy.nc'.format(matrix_path))

        a1file = vmec_current(vmec_current_path)

        shape = (self.nphi, self.nz, self.nr)
        length = (self.nphi*self.nz*self.nr,)
        nphi = self.nphi
        if self.nz*self.nr == numpy.shape(eddyfile.k1xx)[0]:
            shape = (1, self.nr, self.nz)
            length = (self.nz*self.nr,)
            nphi = 1

        if nphi > 1:
            br = self.br[:,:,:]
            bp = self.bp[:,:,:]
            bz = self.bz[:,:,:]
        else:
            br = self.br[1,:,:]
            bp = self.bp[1,:,:]
            bz = self.bz[1,:,:]

        delta_bx = numpy.reshape(br, length, order='F') - eddyfile.bx
        delta_by = numpy.reshape(bp, length, order='F') - eddyfile.by
        delta_bz = numpy.reshape(bz, length, order='F') - eddyfile.bz

        print(eddyfile.bx)

        delta_ax = a1file.ax - eddyfile.ax
        delta_ay = a1file.ay - eddyfile.ay
        delta_az = a1file.az - eddyfile.az

        temp_x = delta_bx - numpy.matmul(eddyfile.k1xx, delta_ax) - numpy.matmul(eddyfile.k1xy, delta_ay) - numpy.matmul(eddyfile.k1xz, delta_az)
        temp_y = delta_by - numpy.matmul(eddyfile.k1yx, delta_ax) - numpy.matmul(eddyfile.k1yy, delta_ay) - numpy.matmul(eddyfile.k1yz, delta_az)
        temp_z = delta_bz - numpy.matmul(eddyfile.k1zx, delta_ax) - numpy.matmul(eddyfile.k1zy, delta_ay) - numpy.matmul(eddyfile.k1zz, delta_az)

        alpha_x = numpy.matmul(numpy.transpose(temp_x), delta_bx)/numpy.matmul(numpy.transpose(temp_x), temp_x)
        alpha_y = numpy.matmul(numpy.transpose(temp_y), delta_by)/numpy.matmul(numpy.transpose(temp_y), temp_y)
        alpha_z = numpy.matmul(numpy.transpose(temp_z), delta_bz)/numpy.matmul(numpy.transpose(temp_z), temp_z)

        br = eddyfile.bx + alpha_x*delta_bx
        bp = eddyfile.by + alpha_y*delta_by
        bz = eddyfile.bz + alpha_z*delta_bz

        if shape[0] == 1:
            for k in range(self.nphi):
                self.br[k,:,:] = numpy.reshape(br, shape, order='F')
                self.bp[k,:,:] = numpy.reshape(bp, shape, order='F')
                self.bz[k,:,:] = numpy.reshape(bz, shape, order='F')

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

    if args['zero_mgrid']:
        mgridfile.set_zero()

    if args['set_fields']:
        if 'vmec_current' not in args:
            print('Faital error: "vmec_current" must be set when using the "set_fields" flag.')
            exit(1)
        if 'matrix_path' not in args:
            print('Faital error: "matrix_path" must be set when using the "set_fields" flag.')
            exit(1)
        mgridfile.set_fields(args['matrix_path'],
                             args['vmec_current'])

    if args['scale_fields']:
        if 'vmec_current' not in args:
            print('Faital error: "vmec_current" must be set when using the "set_fields" flag.')
            exit(1)
        if 'matrix_path' not in args:
            print('Faital error: "matrix_path" must be set when using the "set_fields" flag.')
            exit(1)
        mgridfile.scale_fields(args['matrix_path'],
                               args['vmec_current'])

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
                                     required=False,
                                     action='store',
                                     dest='matrix_path',
                                     help='Path to the folder containing the matrices.',
                                     metavar='MATRIX_PATH')
    command_line_parser.add_argument('-c',
                                     '--vmec_current',
                                     required=False,
                                     action='store',
                                     dest='vmec_current',
                                     help='Path to the file containing the vmec current.',
                                     metavar='VMEC_CURRENT')
    command_line_parser.add_argument('-z',
                                     '--zero_mgrid',
                                     required=False,
                                     action='store_true',
                                     dest='zero_mgrid',
                                     default=False,
                                     help='Set mgrid fields.')
    command_line_parser.add_argument('-st',
                                     '--set_fields',
                                     required=False,
                                     action='store_true',
                                     dest='set_fields',
                                     default=False,
                                     help='Set mgrid fields.')
    command_line_parser.add_argument('-sc',
                                     '--scale_fields',
                                     required=False,
                                     action='store_true',
                                     dest='scale_fields',
                                     default=False,
                                     help='Scale mgrid file.')
    args = vars(command_line_parser.parse_args())

    #  Remove empty arguments
    for key in [key for key in args if args[key] == None]:
        del args[key]

    main(**args)
