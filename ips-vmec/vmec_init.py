#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for VMEC init component. This wapper only takes a VMEC input file
#  and runs VMEC.
#
#-------------------------------------------------------------------------------

from component import Component
import zipfile
import os

#-------------------------------------------------------------------------------
#
#  VMEC init Component Constructor
#
#-------------------------------------------------------------------------------
class vmec_init(Component):
    def __init__(self, services, config):
        print('vmec_init: Construct')
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  VMEC init Component init method. This method prepairs the plasma state. Input
#  files can either be a new namelist input, a new plasma state, or both.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('vmec_init: init')

#  Stage input files.
        self.services.stage_input_files(self.INPUT_FILES)
        
#  Get config filenames.
        self.current_vmec_namelist = self.services.get_config_param('VMEC_NAMELIST_INPUT')
        self.current_vmec_state = self.services.get_config_param('CURRENT_VMEC_STATE')
        
#  Create plasma state from files. Input files can either be a new plasma state,
#  namelist input file or both. If both file were staged, replace the namelist
#  input file.
        if not os.path.exists(self.current_vmec_state) and os.path.exists(self.current_vmec_namelist):
            with zipfile.ZipFile(self.current_vmec_state, 'w') as zip_ref:
                zip_ref.write(self.current_vmec_namelist)
        elif os.path.exists(self.current_vmec_namelist) and os.path.exists(self.current_vmec_namelist):
            with zipfile.ZipFile(self.current_vmec_state, 'r') as zip_old_ref:
                with zipfile.ZipFile('temp.zip', 'w') as zip_new_ref:
                    for item in zip_old_ref.infolist():
                        if item.filename != self.current_vmec_namelist:
                            zip_new_ref.writestr(item, zip_old_ref.read(item.filename))
                        else:
                            zip_new_ref.write(self.current_vmec_namelist)

            os.remove(self.current_vmec_state)
            os.rename('temp.zip', self.current_vmec_state)
        else:
            self.services.error('INPUT_FILES {} are not valid VMEC input files.'.format(self.INPUT_FILES))

        self.services.update_plasma_state()

#-------------------------------------------------------------------------------
#
#  VMEC init Component step method. This runs vmec.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('vmec_init: step')

#-------------------------------------------------------------------------------
#
#  VMEC init Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('vmec_init: finalize')

        if os.path.exists(self.current_vmec_namelist):
            os.remove(self.current_vmec_namelist)
        if os.path.exists(self.current_vmec_state):
            os.remove(self.current_vmec_state)
