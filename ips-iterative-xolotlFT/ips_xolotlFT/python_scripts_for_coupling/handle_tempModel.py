#!/usr/bin/env python

## parameters for the temperature model depend on the Xolotl executable
## v1 (e.g., master) uses a single line: heat = [value], startTemp = [value], etc.
## v2 (e.g., tempGrid) uses two lines:  tempHandler = heat / ...
##                                      heat / tempParam = [value]
## this ends up being quite long, so better in its own file:     
## Xolotl v1 needs one line, v2 needs 2 lines --> write_tempModel contains 2 functions

import shutil
import sys

def v1(xp_parameters={},plasma={},print_test=False):
    print('\t running write_tempModel for Xolotl v1')

    if print_test:
        print('\t with inputs:')
        print('\t \t dictionary xp_parameters =', xp_parameters)
        print('\t \t dictionary plasma =' , plasma)

    print(' ')
    sys.stdout.flush()
    
    #default value
    rm_startTemp=False
    
    if ('heat' in plasma) or ('heat' in xp_parameters) or ('startTemp' in plasma) or ('startTemp' in xp_parameters):
        if 'heat' in plasma:
            mod='heat'
            if (isinstance(plasma['heat'],list) and (len(plasma['heat'])>1)):
                val=[plasma['heat'][0]*1.0e-18,plasma['heat'][1]] 
            else: #no bulkT given in xolotl                                                                                                                                                                                  
                val=[plasma['heat'][0]*1.0e-18,300]
            print('\t \t use heat flux given by PLASMA ', val ,' (here in W/nm2s = 1e-18 W/m2s)')
            if 'startTemp' in xp_parameters:
                rm_startTemp=True
                print('\t \t \t and removed startTemp from xolotl parameters')
        elif 'heat' in xp_parameters:
            mod='heat'
            if (isinstance(xp_parameters['heat'],list) and (len(xp_parameters['heat'])>1)):
                val=[xp_parameters['heat'][0],xp_parameters['heat'][1]]
            else: #no bulkT given in xolotl
                val=[xp_parameters['heat'][0],300]
            print('\t \t no heat flux specified in PLASMA; use values in Xolotl ', val)
            if 'startTemp' in xp_parameters:
                rm_startTemp=True
                print('\t \t \t and removed startTemp from xolotls parameters')
        elif 'startTemp' in plasma:
            mod='startTemp'
            val=plasma['startTemp']
            print('\t \t no heat flux provided by PLASMA or Xolotl')
            print('\t \t use fixed temperature (starTemp) specified by PLASMA ', val)
        elif 'startTemp' in xp_parameters:
            mod='startTemp'
            val=xp_parameters['startTemp']
            print('\t \t no heat flux provided by PLASMA or Xolotl')
            print('\t \t use fixed temperature (startTemp) specified by Xolotl (in config or default): ', val)            
        else:
            print('\t \t WARNING: no heat or temperature defined in PLASMAs output or Xolotls input')
            print('\t \t \t use default value T = 300K')
            mod='startTemp'
            val=300.0

    print('\t ... done running write_tempModel for Xolotl v1')
    sys.stdout.flush()
    return [mod,val,rm_startTemp]

def v2(xp_parameters={},plasma={},print_test=False):

    print('\t running write_tempModel for Xolotl v2')  

    if print_test:
        print('\t with inputs:')
        print('\t \t xp_parameters = ', xp_parameters)
        print('\t \t plasma = ', plasma)
    print(' ')
    sys.stdout.flush()
    
    #default value
    rm_startTemp=False

    if 'tempHandler' in xp_parameters or 'tempHandler' in plasma:
        if 'heat' in plasma:
            print('\t \t use heat defined by PLASMA')
            mod = 'heat'
            if (isinstance(plasma['heat'],list) and (len(plasma['heat'])>1)):
                val=[plasma['heat'][0]*1.0e-18, plasma['heat'][1]] 
            else:
                val=[plasma['heat']*1.0e-18, 300.0] 
            print('\t \t use heat given by PLASMA ', val , '(here used in W/nm2s = 1e-18 W/m2s)')

        elif 'heat' in xp_parameters:
            print('\t \t use heat defined by Xolotl')
            mod = 'heat'
            if (isinstance(xp_parameters['heat'],list) and (len(xp_parameters['heat'])>1)):
                val=[xp_parameters['heat'][0], xp_parameters['heat'][1]]
            else:
                val=[xp_parameters['heat'], 300.0] 
            print('\t \t use heat given by Xolotl ', val) #if given by Xolotl, assume it's in Xolotl's units                                                                                           

        elif 'tempParam' in plasma:
            print('\t \t no heat flux provided')
            print('\t \t use temperature model (tempParam) defined in PLASMAs input:')
            if 'tempHandler' in plasma:
                mod=plasma['tempHandler']
            elif 'tempHandler' in xp_parameters:
                mod=xp_parameters['tempHandler']
            val=plasma['tempParam']
            print('\t tempHandler=', mod, ' and tempParam =', val)

        elif 'tempParam' in xp_parameters:
            print('\t \t no heat flux provided')
            print('\t \t use temperature model (tempParam) definedin FTX config file')
            if 'tempHandler' in xp_parameters:
                mod=xp_parameters['tempHandler']
            elif 'tempHandler' in plasma:
                mod=plasma['tempHandler']
            val=xp_parameters['tempParam']
            print('\t \t tempHandler=', mod, ' and tempParam =', val)

        else:
            print('\t \t WARNING: no heat or other temperature model given in PLASMA or Xolotl. ')
            print('\t \t \t for xolotl version 2, use standard values: ')
            print('\t \t \t tempHandler = constant ; tempParam = 300')
            mod='constant'
            val=300.0
    else:
        print('\t \t WARNING:  no tempHandler provided by xolotl or plasma parameters')
        print('\t \t \t for xolotl version 2, use standard values: ')
        print('\t \t \t tempHandler = constant ; tempParam = 300')
        mod = 'constant'
        val = 300.0

    if 'startTemp' in xp_parameters:
        rm_startTemp=True
        print('\t \t \t and removed startTemp from xolotls parameters')

    print('\t ... done running write_tempModel for Xolotl v2')
    sys.stdout.flush()
    return [mod, val, rm_startTemp]
