import netCDF4
import numpy
import argparse
import json
import time

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
#  DESTRUCTORS
#*******************************************************************************
#-------------------------------------------------------------------------------
#  Destruct a cariddi file instance.
#
#  param[inout] self A mgrid file instance.
#-------------------------------------------------------------------------------
    def __del__(self):
        self.data.close()

#-------------------------------------------------------------------------------
#  Class to represent an vmec current file.
#-------------------------------------------------------------------------------
class vmec_current(carriddi_file):
#*******************************************************************************
#  CONSTRUCTORS
#*******************************************************************************
#-------------------------------------------------------------------------------
#  Initialize a vmec_current file.
#
#  param[inout] self      A carriddi file instance.
#  param[in]    file_name Path to the file.
#-------------------------------------------------------------------------------
    def __init__(self, file_name):
        super(vmec_current, self).__init__(file_name)

        k_x = self.data.variables['k_x'][:].data
        k_y = self.data.variables['k_y'][:].data
        k_z = self.data.variables['k_z'][:].data

        self.k_vmec = numpy.empty(3*len(k_x))
        self.k_vmec[0         :  len(k_x)] = k_x[:]
        self.k_vmec[  len(k_x):2*len(k_x)] = k_y[:]
        self.k_vmec[2*len(k_x):3*len(k_x)] = k_z[:]

#-------------------------------------------------------------------------------
#  Class to represent an xJ Matrix file.
#-------------------------------------------------------------------------------
class xJ_matrix(carriddi_file):
#*******************************************************************************
#  CONSTRUCTORS
#*******************************************************************************
#-------------------------------------------------------------------------------
#  Initialize a K1 file.
#
#  param[inout] self      A k1 file instance.
#  param[in]    file_name Path to the file.
#-------------------------------------------------------------------------------
    def __init__(self, file_name):
        super(xJ_matrix, self).__init__(file_name)

        idx = self.data.variables['is'][:].data[0,:].astype(int) - 1
        self.ind3 = numpy.concatenate((       idx,
                                       500  + idx,
                                       1000 + idx))

#-------------------------------------------------------------------------------
#  Class to represent an K1 Matrix file.
#-------------------------------------------------------------------------------
class k1_matrix(carriddi_file):
#*******************************************************************************
#  CONSTRUCTORS
#*******************************************************************************
#-------------------------------------------------------------------------------
#  Initialize a K1 file.
#
#  param[inout] self      A k1 file instance.
#  param[in]    file_name Path to the file.
#  param[in]    indicies  Ordering indicies for matrix.
#-------------------------------------------------------------------------------
    def __init__(self, file_name, indicies):
        super(k1_matrix, self).__init__(file_name)

        self.x = self.data.variables['K1x'][:].data.transpose()
        self.y = self.data.variables['K1y'][:].data.transpose()
        self.z = self.data.variables['K1z'][:].data.transpose()

        for n in range(numpy.shape(self.x)[0]):
            self.x[n,:] = self.x[n, indicies]
            self.y[n,:] = self.y[n, indicies]
            self.z[n,:] = self.z[n, indicies]

#-------------------------------------------------------------------------------
#  Class to represent an K2 Matrix file.
#-------------------------------------------------------------------------------
class k2_matrix(carriddi_file):
#*******************************************************************************
#  CONSTRUCTORS
#*******************************************************************************
#-------------------------------------------------------------------------------
#  Initialize a K2 file.
#
#  param[inout] self      A k2 file instance.
#  param[in]    file_name Path to the file.
#-------------------------------------------------------------------------------
    def __init__(self, file_name):
        super(k2_matrix, self).__init__(file_name)

        self.x = self.data.variables['K2x'][:].data.transpose()
        self.y = self.data.variables['K2y'][:].data.transpose()
        self.z = self.data.variables['K2z'][:].data.transpose()

#-------------------------------------------------------------------------------
#  Class to represent an C1 Matrix file.
#-------------------------------------------------------------------------------
class c1_matrix(carriddi_file):
#*******************************************************************************
#  CONSTRUCTORS
#*******************************************************************************
#-------------------------------------------------------------------------------
#  Initialize a C1 file.
#
#  param[inout] self      A c1 file instance.
#  param[in]    file_name Path to the file.
#  param[in]    indicies  Ordering indicies for matrix.
#-------------------------------------------------------------------------------
    def __init__(self, file_name, indicies):
        super(c1_matrix, self).__init__(file_name)

        self.c = self.data.variables['C1'][:].data.transpose()

        for n in range(numpy.shape(self.c)[0]):
            self.c[n,:] = self.c[n, indicies].transpose()

#-------------------------------------------------------------------------------
#  Class to represent an C2 Matrix file.
#-------------------------------------------------------------------------------
class c2_matrix(carriddi_file):
#*******************************************************************************
#  CONSTRUCTORS
#*******************************************************************************
#-------------------------------------------------------------------------------
#  Initialize a C2 file.
#
#  param[inout] self      A c2 file instance.
#  param[in]    file_name Path to the file.
#-------------------------------------------------------------------------------
    def __init__(self, file_name):
        super(c2_matrix, self).__init__(file_name)

        self.c = self.data.variables['C2'][:].data.transpose()

#-------------------------------------------------------------------------------
#  Class to represent an vmec current file.
#-------------------------------------------------------------------------------
class eddy:
#*******************************************************************************
#  CONSTRUCTORS
#*******************************************************************************
#-------------------------------------------------------------------------------
#  Initialize a eddy file.
#
#  param[inout] self A eddy file instance.
#  param[in]    data Open netcdf file reference.
#  param[in]    k1x
#  param[in]    k1y
#  param[in]    k1z
#  param[in]    k2x
#  param[in]    k2y
#  param[in]    k2z
#  param[in]    c1
#  param[in]    c2
#  param[in]    dx
#  param[in]    dy
#  param[in]    dz
#  param[in]    c
#-------------------------------------------------------------------------------
    def __init__(self, data,
                 k1x, k1y, k1z, k2x, k2y, k2z,
                 c1, c2, dx, dy, dz, c):

        self.data = data

        self.k1x = k1x
        self.k1y = k1y
        self.k1z = k1z

        self.k2x = k2x
        self.k2y = k2y
        self.k2z = k2z

        self.c1 = c1
        self.c2 = c2

        self.dx = dx
        self.dy = dy
        self.dz = dz

        self.c = c

#-------------------------------------------------------------------------------
#  Create a new eddy file.
#
#  param[in] matrix_path Path to the matrix files.
#  returns An new initialized eddy file object.
#-------------------------------------------------------------------------------
    @classmethod
    def new_eddy(cls, matrix_path):
        data = netCDF4.Dataset('{}/eddy.nc'.format(matrix_path), 'w')

        xJfile = xJ_matrix('{}/cs_xJ.nc'.format(matrix_path))
        c1file = c1_matrix('{}/cs_C1.nc'.format(matrix_path), xJfile.ind3)
        c2file = c2_matrix('{}/cs_C2.nc'.format(matrix_path))
        k1file = k1_matrix('{}/cs_K1.nc'.format(matrix_path), xJfile.ind3)
        k2file = k2_matrix('{}/cs_K2.nc'.format(matrix_path))

        data.description = 'Cariddi matrix re-alligned and vectors for VMEC run time index {:04d}'.format(0)
        data.history = 'Created {}'.format(time.ctime(time.time()))
        data.source = 'netCDF3 python module'

        Ng = data.createDimension('Ng', numpy.shape(k1file.x)[0]) # Number of points on vacuum grid.
        Ns = data.createDimension('Ns', numpy.shape(k1file.x)[1]) # 3 times the number of points on control surface
        Nw = data.createDimension('Nw', numpy.shape(c2file.c)[1]) # Number of points on eddy current.

#-------------------------------------------------------------------------------
        K1x = data.createVariable('K1x', 'f8', ('Ng','Ns'))
        K1y = data.createVariable('K1y', 'f8', ('Ng','Ns'))
        K1z = data.createVariable('K1z', 'f8', ('Ng','Ns'))

        shape = numpy.shape(k1file.x)

        K1x[:,:] = k1file.x.reshape(shape, order='F') # Use Fortran ordering in reshaping
        K1y[:,:] = k1file.y.reshape(shape, order='F') # Use Fortran ordering in reshaping
        K1z[:,:] = k1file.z.reshape(shape, order='F') # Use Fortran ordering in reshaping

#-------------------------------------------------------------------------------
        K2x = data.createVariable('K2x', 'f8', ('Ng', 'Nw'))
        K2y = data.createVariable('K2y', 'f8', ('Ng', 'Nw'))
        K2z = data.createVariable('K2z', 'f8', ('Ng', 'Nw'))

        shape = numpy.shape(k2file.x)

        K2x[:,:] = k2file.x.reshape(shape, order='F') # Use Fortran ordering in reshaping
        K2y[:,:] = k2file.y.reshape(shape, order='F') # Use Fortran ordering in reshaping
        K2z[:,:] = k2file.z.reshape(shape, order='F') # Use Fortran ordering in reshaping

#-------------------------------------------------------------------------------
        C1 = data.createVariable('C1', 'f8', ('Nw', 'Ns'))

        shape = numpy.shape(c1file.c)

        C1[:,:] = c1file.c.reshape(shape, order='F') # Use Fortran ordering in reshaping

#-------------------------------------------------------------------------------
        C2 = data.createVariable('C2', 'f8', ('Nw', 'Nw'))

        shape = numpy.shape(c2file.c)

        C2[:,:] = c2file.c.reshape(shape, order='F') #Use Fortran ordering in reshaping

#-------------------------------------------------------------------------------
        dx = data.createVariable('dx', 'f8', ('Ng'))
        dy = data.createVariable('dy', 'f8', ('Ng'))
        dz = data.createVariable('dz', 'f8', ('Ng'))

#-------------------------------------------------------------------------------
        c = data.createVariable('c', 'f8', ('Nw'))

        return cls(data, K1x, K1y, K1z, K2x, K2y, K2z,
                   C1, C2, dx, dy, dz, c)

#-------------------------------------------------------------------------------
#  Load an existing eddy file.
#
#  param[in] matrix_path Path to the matrix files.
#  returns An open eddy file object.
#-------------------------------------------------------------------------------
    @classmethod
    def load_eddy(cls, matrix_path):
        data = netCDF4.Dataset('{}/eddy.nc'.format(matrix_path), 'a')

        dx = data.variables['dx']
        dy = data.variables['dy']
        dz = data.variables['dz']

        c = data.variables['c']

        return cls(data,
                   data.variables['K1x'][:].data,
                   data.variables['K1y'][:].data,
                   data.variables['K1z'][:].data,
                   data.variables['K2x'][:].data,
                   data.variables['K2y'][:].data,
                   data.variables['K2z'][:].data,
                   data.variables['C1'][:].data,
                   data.variables['C2'][:].data,
                   dx, dy, dz, c)

#*******************************************************************************
#  DESTRUCTORS
#*******************************************************************************
#-------------------------------------------------------------------------------
#  Destruct a cariddi file instance.
#
#  param[inout] self A mgrid file instance.
#-------------------------------------------------------------------------------
    def __del__(self):
        self.data.close()

#*******************************************************************************
#  SETTERS
#*******************************************************************************
#-------------------------------------------------------------------------------
#  Set vmec current.
#
#  This updates the d* and c matricies with values of the new current. At t=0,
#  these values are zero since there is no surface current initially.
#
#  param[inout] self         A eddy file object instance.
#  param[in]    time_index   Current time index.
#  param[in]    current_path Path to the J0 file.
#-------------------------------------------------------------------------------
    def set_current(self, time_index, current_path):
        J0 = vmec_current(current_path)

#  FIXME: Size of this matrix is hard coded.
        I0A = numpy.zeros(484)
        if time_index != 0:
            I0A = self.c1.dot(J0.k_vmec) + self.c

#  Compute updated quantities.
        self.c[:]  = -self.c1.dot(J0.k_vmec) + self.c2.dot(I0A)
        self.dx[:] = -self.k1x.dot(J0.k_vmec) + self.k2x.dot(I0A)
        self.dy[:] = -self.k1y.dot(J0.k_vmec) + self.k2y.dot(I0A)
        self.dz[:] = -self.k1z.dot(J0.k_vmec) + self.k2z.dot(I0A)

#-------------------------------------------------------------------------------
#  Extract profile parameters for the time index.
#
#  Profile data is written into a json dictionary file that can be passed to
#  ips-vmec to update namelist input files.
#
#  param[in] percorso   Precurser path for the stored profile quantities.
#  param[in] time_index Current time index.
#  param[in] vmec_input JSON file to write vmec profile parameters to.
#-------------------------------------------------------------------------------
def set_namelist_parameters(percorso, time_index, vmec_input):
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
#  param[in] matrix_path  Path to the matrix files.
#  param[in] time_index   Current time index.
#  param[in] vmec_current VMEC J0 current file.
#-------------------------------------------------------------------------------
def set_matrix_file(matrix_path, time_index, vmec_current):
    eddyfile = None
    if time_index == 0:
        eddyfile = eddy.new_eddy(matrix_path)
    else:
        eddyfile = eddy.load_eddy(matrix_path)

    eddyfile.set_current(time_index, vmec_current)

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

#  If the vmec_input flag was used write out the new VMEC profiles.
    if 'vmec_input' in args:
        set_namelist_parameters(args['matrix_path'],
                                args['time_index'],
                                args['vmec_input'])

#  If the vmec_currents update the matrix file.
    if 'vmec_current' in args:
        set_matrix_file(args['matrix_path'],
                        args['time_index'],
                        args['vmec_current'])

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

    args = vars(parser.parse_args())

    #  Remove empty arguments
    for key in [key for key in args if args[key] == None]:
        del args[key]

    main(**args)
