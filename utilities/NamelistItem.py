#-------------------------------------------------------------------------------
#
#  This utility parses a key to separate out the indicies. The it sets the file
#  using the OMFIT fortran environment.
#
#-------------------------------------------------------------------------------

from omfit.classes.namelist import fortran_environment

def set(namelist, key, value):
    if '(' in key:
        key, indicies = key.split('(')
        indicies = indicies.split(')')[0]
        with fortran_environment(namelist):
            exec("namelist['{}'][{}]={}".format(key, indicies, value))
    else:
        namelist[key] = value
