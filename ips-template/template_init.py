#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for TEMPLATE Init component. This template shows an example of
#  how to set up an ips component init.
#
#-------------------------------------------------------------------------------

from component import Component
import os
import zipfile

#-------------------------------------------------------------------------------
#
#  TEMPLATE Init Class
#
#-------------------------------------------------------------------------------
class template_init(Component):
    
#-------------------------------------------------------------------------------
#
#  TEMPLATE Init Component Constructor
#
#-------------------------------------------------------------------------------
    def __init__(self, services, config):
        print('template_init: Construct')
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  template_init Component init method. This method prepairs the input files.
#  This allows staging the plasma state files and creates the inital state.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('template_init: init')
    
#  Stage input files.
        self.services.stage_input_files(self.INPUT_FILES)
        input_files = self.INPUT_FILES.split(' ')

#  Some times codes expect input files with a specific names. In these cases
#  rename those files.
        os.rename(self.services.get_config_param('INPUT_3'), 'rename.nc')
        input_files = ['rename.nc' if x == self.services.get_config_param('INPUT_3') else x for x in input_files]

#  Create plasma state from files.
        zip_ref = zipfile.ZipFile(self.services.get_config_param('CURRENT_TEMPLATE_STATE'), 'w')

#  Add files to the zip file.
        for file in input_files:
            zip_ref.write(file)

#  Close the file.
        zip_ref.close()

#  Update the plasma state.
        self.services.update_plasma_state()

#-------------------------------------------------------------------------------
#
#  template_init Component step method. This component does nothing and is
#  never called.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('template_init: step')

#-------------------------------------------------------------------------------
#
#  template_init Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('template_init: finalize')
