#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for VMEC component. This wapper only takes a VMEC input file and
#  runs VMEC.
#
#-------------------------------------------------------------------------------

from component import Component
import os
from omfit.classes.omfit_namelist import OMFITnamelist
import zipfile

#-------------------------------------------------------------------------------
#
#  VMEC Component Constructor
#
#-------------------------------------------------------------------------------
class vmec(Component):
    def __init__(self, services, config):
        print('vmec: Construct')
        Component.__init__(self, services, config)
        self.task_queue = {}

#-------------------------------------------------------------------------------
#
#  VMEC Component init method. This method prepairs the namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0, **keywords):
        print('vmec: init')
        
        self.current_vmec_namelist = self.services.get_config_param('VMEC_NAMELIST_INPUT')
        self.current_state = self.services.get_config_param('CURRENT_VMEC_STATE')
        
#  Stage plasma state.
        self.services.stage_plasma_state()

#  Unzip files from the plasma state. Use mode a so files can be read and
#  written to.
        self.zip_ref = zipfile.ZipFile(self.current_state, 'a')
        self.zip_ref.extract(self.current_vmec_namelist)
        
#  Update parameters in the namelist.
        namelist = OMFITnamelist(self.current_vmec_namelist)

        for key, value in keywords.iteritems():
            if '(' in key :
                key, indices = key.split('(')
                indices, extra = indices.split(')')
                indices = [[int(i) - 1] for i in indices.split(',')]
                namelist['indata'][key][indices] = value
            else:
                namelist['indata'][key] = value
    
        namelist.save()

#-------------------------------------------------------------------------------
#
#  VMEC Component step method. This runs vmec.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('vmec: step')

        self.task_queue['vmec'] = self.services.launch_task(self.NPROC,
                                                            self.services.get_working_dir(),
                                                            self.VMEC_EXE,
                                                            self.current_vmec_namelist,
                                                            logfile = 'vmec.log')
    
#  While vmec is running check the current plasma state for an existing wout
#  file. If such a file exists replace it from the file.
        wout_file = 'wout_{}.nc'.format(self.current_vmec_namelist.replace('input.','',1))
    
        replace = False
        for item in self.zip_ref.infolist():
            if item.filename == wout_file:
                replace = True
                break

        if replace:
            print('No')
            with zipfile.ZipFile('temp.zip', 'w') as zip_new_ref:
                for item in self.zip_ref.infolist():
                    if item.filename != wout_file:
                        zip_new_ref.writestr(item, zip_old_ref.read(item.filename))
                self.zip_ref.close()
        
            os.remove(self.current_vmec_state)
            os.rename('temp.zip', self.current_vmec_state)
            self.zip_ref = zipfile.ZipFile(self.current_state, 'a')

#  Wait for VMEC to finish.
        if (self.services.wait_task(self.task_queue['vmec']) or not os.path.exists(wout_file)):
            self.services.error('vmec: step failed.')
        del self.task_queue['vmec']

#  Add the wout file to the plasma state.
        self.zip_ref.write(wout_file)
        self.zip_ref.close()

        self.services.update_plasma_state()

#-------------------------------------------------------------------------------
#
#  VMEC Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('vmec: finalize')
        self.services.wait_tasklist(self.task_queue.values(), True)
