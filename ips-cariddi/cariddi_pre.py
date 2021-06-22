#! /usr/bin/env python

import argparse
import scipy.io
import netCDF4
import numpy
import json

#-------------------------------------------------------------------------------
#  Save fields.
#
#  param[in] path      Path to the file.
#  param[in] file_name Name of the mgrid file.
#-------------------------------------------------------------------------------
def save_fields(path, file_name):
    with netCDF4.Dataset('{}/eddy.nc'.format(path), 'a') as eddy:

        with netCDF4.Dataset(file_name, 'a') as mgrid:
#  FIXME: Doesn't handle 3D fields.
            br = mgrid.variables['br_006'][:].data[1,:,:]
            bp = mgrid.variables['bp_006'][:].data[1,:,:]
            bz = mgrid.variables['bz_006'][:].data[1,:,:]

        bx = eddy.variables['bx']
        by = eddy.variables['by']
        bz = eddy.variables['bz']

        bx[:] = numpy.reshape(br, numpy.shape(bx), order='F')
        by[:] = numpy.reshape(bp, numpy.shape(bx), order='F')
        bz[:] = numpy.reshape(bz, numpy.shape(bx), order='F')

    print('Fields saved')
    
#-------------------------------------------------------------------------------
#  Initialize a eddy file.
#
#  param[in] path       Path to the file.
#  param[in] time_index Current time index.
#  param[in] file_name  Name of the file.
#-------------------------------------------------------------------------------
def set_current(path, time_index, file_name):
    with netCDF4.Dataset('{}/eddy.nc'.format(path), 'a') as eddy:

        with netCDF4.Dataset(file_name, 'a') as current:
            ax = current.variables['ax'][:].data
            ay = current.variables['ay'][:].data
            az = current.variables['az'][:].data

        c1x = eddy.variables['c1x'][:].data
        c1y = eddy.variables['c1y'][:].data
        c1z = eddy.variables['c1z'][:].data

        c2 = eddy.variables['c2'][:].data

        c = eddy.variables['c']

        k1xx = eddy.variables['k1xx'][:].data
        k1xy = eddy.variables['k1xy'][:].data
        k1xz = eddy.variables['k1xz'][:].data

        k1yx = eddy.variables['k1yx'][:].data
        k1yy = eddy.variables['k1yy'][:].data
        k1yz = eddy.variables['k1yz'][:].data

        k1zx = eddy.variables['k1zx'][:].data
        k1zy = eddy.variables['k1zy'][:].data
        k1zz = eddy.variables['k1zz'][:].data

        k2x = eddy.variables['k2x'][:].data
        k2y = eddy.variables['k2y'][:].data
        k2z = eddy.variables['k2z'][:].data

        dx = eddy.variables['dx']
        dy = eddy.variables['dy']
        dz = eddy.variables['dz']

        i = numpy.zeros(numpy.shape(c[:].data))
        if time_index > 0:
            i = numpy.matmul(c1x, ax) + numpy.matmul(c1y, ay) + numpy.matmul(c1z, az) + c[:].data

        c[:] = -numpy.matmul(c1x, ax) - numpy.matmul(c1y, ay) - numpy.matmul(c1z, az) + numpy.matmul(c2, i)

        dx[:] = -numpy.matmul(k1xx, ax) - numpy.matmul(k1xy, ay) - numpy.matmul(k1xz, az) + numpy.matmul(k2x, i)
        dy[:] = -numpy.matmul(k1yx, ax) - numpy.matmul(k1yy, ay) - numpy.matmul(k1yz, az) + numpy.matmul(k2y, i)
        dz[:] = -numpy.matmul(k1zx, ax) - numpy.matmul(k1zy, ay) - numpy.matmul(k1zz, az) + numpy.matmul(k2z, i)

    print('eddy.nc file updated')

#-------------------------------------------------------------------------------
#  Class to represent an matrix file.
#-------------------------------------------------------------------------------
class matrix:
#*******************************************************************************
#  CONSTRUCTORS
#*******************************************************************************
#-------------------------------------------------------------------------------
#  Initialize a matrix file.
#
#  param[inout] self      A matrix file instance.
#  param[in]    path      Path to the file.
#  param[in]    file_name Name of the file.
#-------------------------------------------------------------------------------
    def __init__(self, path, file_name):
        mat = scipy.io.loadmat('{}/{}'.format(path, file_name))

        self.path = path

        self.c1x = mat['C1x']
        self.c1y = mat['C1y']
        self.c1z = mat['C1z']
        self.c2 = mat['C2']

        self.k1xx = mat['K1xx']
        self.k1xy = mat['K1xy']
        self.k1xz = mat['K1xz']
        self.k2x = mat['K2x']

        self.k1yx = mat['K1yx']
        self.k1yy = mat['K1yy']
        self.k1yz = mat['K1yz']
        self.k2y = mat['K2y']

        self.k1zx = mat['K1zx']
        self.k1zy = mat['K1zy']
        self.k1zz = mat['K1zz']
        self.k2z = mat['K2z']

#-------------------------------------------------------------------------------
#  Write out the netcdf matrix files.
#
#  param[in] self A matrix file instance.
#-------------------------------------------------------------------------------
    def write(self):
        with netCDF4.Dataset('{}/eddy.nc'.format(self.path), 'w') as eddy:

            nf = eddy.createDimension('nf', numpy.shape(self.k2x)[1])
            ng = eddy.createDimension('ng', numpy.shape(self.k1xx)[1])
            nv = eddy.createDimension('nv', numpy.shape(self.k1xx)[0])
            scalar = eddy.createDimension('s', 1)

            c = eddy.createVariable('c', 'f8', ('nf'))

            c[:] = 0

            dx = eddy.createVariable('dx', 'f8', ('nv'))
            dy = eddy.createVariable('dy', 'f8', ('nv'))
            dz = eddy.createVariable('dz', 'f8', ('nv'))

            dx[:] = 0
            dy[:] = 0
            dz[:] = 0

            k1xx = eddy.createVariable('k1xx', 'f8', ('nv','ng'))
            k1xy = eddy.createVariable('k1xy', 'f8', ('nv','ng'))
            k1xz = eddy.createVariable('k1xz', 'f8', ('nv','ng'))

            k1yx = eddy.createVariable('k1yx', 'f8', ('nv','ng'))
            k1yy = eddy.createVariable('k1yy', 'f8', ('nv','ng'))
            k1yz = eddy.createVariable('k1yz', 'f8', ('nv','ng'))

            k1zx = eddy.createVariable('k1zx', 'f8', ('nv','ng'))
            k1zy = eddy.createVariable('k1zy', 'f8', ('nv','ng'))
            k1zz = eddy.createVariable('k1zz', 'f8', ('nv','ng'))

            k1xx[:,:] = self.k1xx
            k1xy[:,:] = self.k1xy
            k1xz[:,:] = self.k1xz

            k1yx[:,:] = self.k1yx
            k1yy[:,:] = self.k1yy
            k1yz[:,:] = self.k1yz

            k1zx[:,:] = self.k1zx
            k1zy[:,:] = self.k1zy
            k1zz[:,:] = self.k1zz

            k2x = eddy.createVariable('k2x', 'f8', ('nv','nf'))
            k2y = eddy.createVariable('k2y', 'f8', ('nv','nf'))
            k2z = eddy.createVariable('k2z', 'f8', ('nv','nf'))

            k2x[:,:] = self.k2x
            k2y[:,:] = self.k2y
            k2z[:,:] = self.k2z

            c1x = eddy.createVariable('c1x', 'f8', ('nf','ng'))
            c1y = eddy.createVariable('c1y', 'f8', ('nf','ng'))
            c1z = eddy.createVariable('c1z', 'f8', ('nf','ng'))

            c1x[:,:] = self.c1x
            c1y[:,:] = self.c1y
            c1z[:,:] = self.c1z

            c2 = eddy.createVariable('c2', 'f8', ('nf','nf'))

            c2[:,:] = self.c2

            bx = eddy.createVariable('bx', 'f8', 'nv')
            by = eddy.createVariable('by', 'f8', 'nv')
            bz = eddy.createVariable('bz', 'f8', 'nv')

            bx[:] = 0
            by[:] = 0
            bz[:] = 0

            ax = eddy.createVariable('ax', 'f8', ('ng'))
            ay = eddy.createVariable('ay', 'f8', ('ng'))
            az = eddy.createVariable('az', 'f8', ('ng'))

            ax[:] = 0
            ay[:] = 0
            az[:] = 0

#-------------------------------------------------------------------------------
#  Class to represent an surface file.
#-------------------------------------------------------------------------------
class surface:
#*******************************************************************************
#  CONSTRUCTORS
#*******************************************************************************
#-------------------------------------------------------------------------------
#  Initialize a surface file.
#
#  param[inout] self      A surface file instance.
#  param[in]    path      Path to the file.
#  param[in]    file_name Name of the file.
#-------------------------------------------------------------------------------
    def __init__(self, path, file_name):
        mat = scipy.io.loadmat('{}/{}'.format(path, file_name))

        self.x = mat['xel']
        self.y = mat['yel']
        self.z = mat['zel']

#-------------------------------------------------------------------------------
#  Write out the netcdf matrix files.
#
#  param[in] self A matrix file instance.
#-------------------------------------------------------------------------------
    def write(self):
        with netCDF4.Dataset('A1.nc', 'w') as a1:

            ng = a1.createDimension('ng', numpy.shape(self.x)[0])

            ax = a1.createVariable('ax', 'f8', ('ng'))
            ay = a1.createVariable('ay', 'f8', ('ng'))
            az = a1.createVariable('az', 'f8', ('ng'))

            ax[:] = 0
            ay[:] = 0
            az[:] = 0

            x = a1.createVariable('x', 'f8', ('ng'))
            y = a1.createVariable('y', 'f8', ('ng'))
            z = a1.createVariable('z', 'f8', ('ng'))

            x[:] = self.x[:,0]
            y[:] = self.y[:,0]
            z[:] = self.z[:,0]

#-------------------------------------------------------------------------------
#  Convert matrix files to netcdf.
#
#  param[in] path      Path to the file.
#  param[in] file_name Name of the file.
#-------------------------------------------------------------------------------
def set_matrix(path, file_name):
    mat = matrix(path, file_name)
    mat.write()
    print('eddy.nc file created')

#-------------------------------------------------------------------------------
#  Convert surface files to netcdf.
#
#  param[in] path      Path to the file.
#  param[in] file_name Name of the file.
#-------------------------------------------------------------------------------
def set_surface(path, file_name):
    surf = surface(path, file_name)
    surf.write()
    print('A1.nc file created')

#-------------------------------------------------------------------------------
#  Extract profile parameters for the time index.
#
#  Profile data is written into a json dictionary file that can be passed to
#  ips-vmec to update namelist input files.
#
#  param[in] path       Path for the stored profile quantities.
#  param[in] time_index Current time index.
#  param[in] vmec_input JSON file to write vmec profile parameters to.
#-------------------------------------------------------------------------------
def set_namelist_parameters(path, time_index, vmec_input):
#  File with profiles from Cariddi.
    file_prof = '{}/SURFACEQUANT_{:04d}.dat'.format(path, time_index)

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

    profile_changes['vmec__pmass_type'] = 'akima_spline'
    profile_changes['vmec__pcurr_type'] = 'akima_spline_I'
    for i, s in enumerate(ss):
        profile_changes['vmec__ac_aux_s({})'.format(i + 1)] = s
        profile_changes['vmec__am_aux_s({})'.format(i + 1)] = s
    for i, c in enumerate(II_tor):
        profile_changes['vmec__ac_aux_f({})'.format(i + 1)] = c
    for i, p in enumerate(pp):
        profile_changes['vmec__am_aux_f({})'.format(i + 1)] = p

    profile_changes['vmec__phiedge'] = phi_tot[-1]
#  FIXME: Test convergence by decreasing the total toroidal current.
    profile_changes['vmec__curtor'] = I_tor[-1]

    with open(vmec_input, 'w') as json_ref:
        json.dump(profile_changes, json_ref, indent=4)

    print('VMEC profile profile parameters written to {}.'.format(vmec_input))

#-------------------------------------------------------------------------------
#  Main routine.
#
#  This routine runs a number of different tasks depending on the command line
#  arguments set.
#
#  param[in] args Dictionary of commandline arguments.
#-------------------------------------------------------------------------------
def main(**args):
    print('PRE-PROCESSING time index {:04d}\n'.format(args['time_index']))
    print('reading data from files in {}\n'.format(args['matrix_path']))

#  If the matrix_file flag was used create the first eddy file.
    if 'matrix_file' in args:
        set_matrix(args['matrix_path'],
                   args['matrix_file'])

#  If the coordinate_file flag was used create the first A1 file.
    if 'coordinate_file' in args:
        set_surface(args['matrix_path'],
                    args['coordinate_file'])

#  If the vmec_input flag was used write out the new VMEC profiles.
    if 'vmec_input' in args:
        set_namelist_parameters(args['matrix_path'],
                                args['time_index'],
                                args['vmec_input'])

#  If the vmec_current flag was used update the eddy file.
    if 'vmec_current' in args:
        set_current(args['matrix_path'],
                    args['time_index'],
                    args['vmec_current'])

#  If mgrid_file flag was used was used update save the mgird field.
    if 'mgrid_file' in args:
        save_fields(args['matrix_path'],
                    args['mgrid_file'])

#-------------------------------------------------------------------------------
#  Parse commandline arguments.
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Mat to netcdf conversion utility')
    parser.add_argument('-c',
                        '--coordinate_file',
                        action='store',
                        dest='coordinate_file',
                        help='Name of the coordinate data.',
                        metavar='COORDINATE_FILE')
    parser.add_argument('-m',
                        '--matrix_file',
                        action='store',
                        dest='matrix_file',
                        help='Name of the matrix data.',
                        metavar='MATRIX_FILE')
    parser.add_argument('-mg',
                        '--mgrid_file',
                        action='store',
                        dest='mgrid_file',
                        help='Name of the mgrid file.',
                        metavar='MGRID_FILE')
    parser.add_argument('-p',
                        '--matrix_path',
                        action='store',
                        dest='matrix_path',
                        help='Path to the matrix directory.',
                        metavar='MATRIX_PATH')
    parser.add_argument('-i',
                        '--vmec_input',
                        required=False,
                        action='store',
                        dest='vmec_input',
                        help='Path to the file containing the vmec namelist input.',
                        metavar='VMEC_INPUT')
    parser.add_argument('-t',
                        '--time_index',
                        required=False,
                        action='store',
                        dest='time_index',
                        help='time index (integer)',
                        type=int,
                        default=0,
                        metavar='TIME_INDEX')
    parser.add_argument('-v',
                        '--vmec_current',
                        required=False,
                        action='store',
                        dest='vmec_current',
                        help='Path to the file containing the vmec current',
                        metavar='VMEC_CURRENT')

    args = vars(parser.parse_args())

#  Remove empty arguments
    for key in [key for key in args if args[key] == None]:
        del args[key]

    if 'vmec_input' in args and 'time_index' not in args:
        print('Faital error: "time_index" must be set when using the "vmec_input" flag.')
        exit(1)
    if 'vmec_current' in args and 'time_index' not in args:
        print('Faital error: "time_index" must be set when using the "vmec_current" flag.')
        exit(1)

    main(**args)
