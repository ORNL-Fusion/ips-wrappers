#!/usr/bin/env python

## parameters for the grid model depend on the Xolotl executable
## v1 (e.g., master) uses a single line: grid = [value]
## v2 (e.g., tempGrid) uses two lines:  gridType = nonuniform
##                                      gridParam = [value]
## this ends up being quite long, so better in its own file:

## Xolotl v1 needs one line, v2 needs 2 lines --> write_gridModel contains 2 functions              

import sys

def v1(xp_parameters={}, print_test=False):

    rm_gridParam=False
    rm_gridType=False

    if print_test:
        print('\t called handle_gridModel v1: ')
    
    if 'grid' in xp_parameters:
        print('\t for xolotl_v = 1, grid exists: ', xp_parameters['grid'])
        gridValue=xp_parameters['grid']
    else:
        print('\t WARNING: for xolol_v = 1 , could not find parameter grid')
        if 'gridParam' in xp_parameters:
            print('\t \t found gridParam in xolotl parameters instead; use it as grid value: ', xp_parameters['gridParam'])
            gridValue=xp_parameters['gridParam']
        else:
            print('\t \t could not find gridParam either. Use default value 200')
            gridValue=[200, 0.5]
    if 'gridParam' in xp_parameters:
        print('\t found gridParam in xolotl parameters ; delete from dictionary to avoid it in param file')
        rm_gridParam=True
    if 'gridType' in xp_parameters:
        print('\t found gridType in xolotl parameters ; delete from dictionary to avoid it in param file')
        rm_gridType=True
    sys.stdout.flush()

    if print_test:
        print('\t TEST: for v1, handle_gripModel returns gridValue, rm_gridType and rm_gridValue')
        
    return [gridValue, rm_gridType, rm_gridParam]

def v2(xp_parameters={}, print_test=False):
    
    rm_grid=False
    rm_regGrid=False
    
    if print_test:
        print('\t called handle_gridModel v2: ')
    
    if ('gridType' in xp_parameters) and ('gridParam' in xp_parameters):
        print('\t for xolotl_v = 2, gridType exists: ', xp_parameters['gridType'],', and gridParam exists: ', xp_parameters['gridParam'])
        gridType=xp_parameters['gridType']
        gridValue=xp_parameters['gridParam']
    if ('gridType' not in xp_parameters):
        print('\t WARNING: for xolol_v = 2 , could not find parameter gridType. Assume nonuniform')
        gridType='nonuniform'
    if ('gridParam' not in xp_parameters):
        print('\t WARNING: for xolol_v = 2 , could not find parameter gridParam.')
        if 'grid' in xp_parameters:
            print('\t \t found grid in xolotl parameters instead ; use it as gridParam value: ', xp_parameters['grid'])
            gridValue=xp_parameters['grid']
        else:
            print('\t \t could not find grid either. Use default value 200')
            gridValue=[200, 0.5]
    sys.stdout.flush()
    
    if 'grid' in xp_parameters:
        print('\t found "grid" in xolotl parameters ; delete from dictionary to avoid it in param file')
        rm_grid=True

    if 'regularGrid' in xp_parameters:
        print('\t found "regularGrid" in xolotl parameters ; delete from dictionary to avoid it in param file')
        rm_regularGrid=True
        
    if print_test:
        print('\t TEST: for v2, handle_gripModel returns gridType, gridValue and rm_grid')
    sys.stdout.flush()
        
    return [gridType, gridValue, rm_grid, rm_regularGrid]

