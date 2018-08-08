#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for SOLPS-ITER_init component. Take the work flow inputs and
#  generates a plasma state.
#
#-------------------------------------------------------------------------------

from component import Component
import os
import zipfile

#-------------------------------------------------------------------------------
#
#  SOLPS-ITER_init Component Constructor
#
#-------------------------------------------------------------------------------
class solps_iter_init(Component):
    def __init__(self, services, config):
        print('solps_iter_init: Construct')
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  SOLPS-ITER_init Component init method. This method prepairs the input files.
#  This allows staging the plasma state files.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('solps_iter_init: init')

        self.services.stage_input_files(self.INPUT_FILES)

#  Rename the eirene input files.
        os.rename(self.services.get_config_param('EIRENE_INPUT_DAT'), 'fort.1')
        os.rename(self.services.get_config_param('EIRENE_NODES'), 'fort.33')
        os.rename(self.services.get_config_param('EIRENE_CELLS'), 'fort.34')
        os.rename(self.services.get_config_param('EIRENE_LINKS'), 'fort.35')

#  Create plasma state zip file.
        zip_ref = zipfile.ZipFile(self.services.get_config_param('CURRENT_SOLPS_STATE'), 'w')

#  b2 files
        zip_ref.write('b2fgmtry')
        zip_ref.write('b2fpardf')
        zip_ref.write('b2frates')
        zip_ref.write('b2fstati')
        zip_ref.write('b2mn.dat')

        zip_ref.write('b2.transport.parameters')
        zip_ref.write('b2.numerics.parameters')
        zip_ref.write('b2.neutrals.parameters')
        zip_ref.write('b2.boundary.parameters')

#  eirene files
        zip_ref.write('fort.1')
        zip_ref.write('fort.33')
        zip_ref.write('fort.34')
        zip_ref.write('fort.35')

        zip_ref.close()

        self.services.update_plasma_state()

#-------------------------------------------------------------------------------
#
#  SOLPS-ITER_init init Component step method. This runs vmec.
#
#-------------------------------------------------------------------------------
        def step(self, timeStamp=0.0):
            print('solps_iter_init: step')
    
#-------------------------------------------------------------------------------
#
#  SOLPS-ITER_init init Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
        def finalize(self, timeStamp=0.0):
            print('solps_iter_init: finalize')
