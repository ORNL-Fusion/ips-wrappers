#! /usr/bin/env python

from  component import Component
import fileinput
import numpy as np
import solps
import gitr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4
import math
class solps_iter_data_worker(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeStamp=0.0):
        return

    def step(self, timeStamp=0.0):
        self.services.stage_plasma_state()
        
        print('Hello from solps-data-iter worker')
        
        dict = {'task':self.TASK,
                'rmin':float(self.RMIN),
                'rmax':float(self.RMAX),
                'numr':int(self.NUMR),
                'zmin':float(self.ZMIN),
                'zmax':float(self.ZMAX),
                'numz':int(self.NUMZ),
                'solps_geometry':self.SOLPS_GEOMETRY,
                'solps_state':self.SOLPS_STATE,
                'left_target':self.LEFT_TARGET,
                'right_target':self.RIGHT_TARGET,
                'result':self.RESULT,
                'dakota_result':self.DAKOTA_RESULT
                }
        
        dict2 = {
                'rmin_sep':float(self.RMIN_SEP),
                'rmax_sep':float(self.RMAX_SEP),
                'zmin_sep':float(self.ZMIN_SEP),
                'zmax_sep':float(self.ZMAX_SEP)
                }
        
        arguments='-'+' -'.join("%s=%s" % (key,val) for (key,val) in dict.items())
        print('solps-data arguments ', arguments)
        
        dict['solps_equilibrium'] = self.SOLPS_EQUILIBRIUM 
        dict['solps_mesh'] = self.SOLPS_MESH 
        dict['gitr_geometry'] = self.GITR_GEOMETRY 
        print('Regridding SOLPS Data from %s file using solps geometry %s' %(dict['solps_state'], dict['solps_geometry']))
        print('Making use of data from %s equilibrium file' %(dict['solps_equilibrium']))
        
        task_id = self.services.launch_task(self.NPROC,
                        self.services.get_working_dir(),
                        self.EXE,arguments,task_ppn=1,logfile='task.log')
        
        #monitor task until complete
        if (self.services.wait_task(task_id)):
            self.services.error('solps_iter_data_worker: step failed.')
            return
        
        solps.process_solps_output_for_gitr(dakota_filename = dict['dakota_result'], \
                                  nR = dict['numr'], nZ = dict['numz'], plot_variables=1, \
                                  b2fstate_filename = dict['solps_state'])
        
        gitr.make_gitr_geometry_from_solps(gitr_geometry_filename=dict['gitr_geometry'], \
                                  solps_mesh_extra=dict['solps_mesh'], \
                                  solps_geom=dict['solps_geometry'])
        
        solps.make_solps_targ_file(gitr_geom_filename=dict['gitr_geometry'], \
                                   solps_geom  = dict['solps_geometry'], \
                                   right_target_filename = dict['right_target']) 
        
        self.services.update_state()

        return
    
    def finalize(self, timeStamp=0.0):
        return
