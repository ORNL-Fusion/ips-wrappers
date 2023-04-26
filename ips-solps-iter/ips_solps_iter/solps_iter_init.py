#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for SOLPS-ITER_init component. Take the work flow inputs and
#  generates a state.
#
#-------------------------------------------------------------------------------

from ipsframework import Component
import os
from utilities import ZipState
from utilities import ScreenWriter

#-------------------------------------------------------------------------------
#
#  SOLPS-ITER_init Component Constructor
#
#-------------------------------------------------------------------------------
class solps_iter_init(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  SOLPS-ITER_init Component init method. This method prepairs the input files.
#  This allows staging the state files.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'solps_iter_init: init')

#  Get config filenames.
        current_solps_state = self.services.get_config_param('CURRENT_SOLPS_STATE')
        eirene_input_dat = self.services.get_config_param('EIRENE_INPUT_DAT')
        eirene_nodes = self.services.get_config_param('EIRENE_NODES')
        eirene_cells = self.services.get_config_param('EIRENE_CELLS')
        eirene_links = self.services.get_config_param('EIRENE_LINKS')
        eirene_grid = self.services.get_config_param('EIRENE_GRID')

#  Remove old inputs. Stage input files.
        for file in os.listdir('.'):
            os.remove(file)
    
        self.services.stage_input_files(self.INPUT_FILES)

#  Create state zip file.
        with ZipState.ZipState(current_solps_state, 'a') as zip_ref:
#  b2 files
            if os.path.exists('b2fgmtry'):
                zip_ref.write('b2fgmtry')
                zip_ref.set_state(state='needs_update')
            if os.path.exists('b2fpardf'):
                zip_ref.write('b2fpardf')
                zip_ref.set_state(state='needs_update')
            if os.path.exists('b2frates'):
                zip_ref.write('b2frates')
                zip_ref.set_state(state='needs_update')
            if os.path.exists('b2fstati'):
                zip_ref.write('b2fstati')
                zip_ref.set_state(state='needs_update')
            if os.path.exists('b2mn.dat'):
                zip_ref.write('b2mn.dat')
                zip_ref.set_state(state='needs_update')

#  Namelist files.
            if os.path.exists('b2.transport.parameters'):
                zip_ref.write('b2.transport.parameters')
                zip_ref.set_state(state='needs_update')
            if os.path.exists('b2.numerics.parameters'):
                zip_ref.write('b2.numerics.parameters')
                zip_ref.set_state(state='needs_update')
            if os.path.exists('b2.neutrals.parameters'):
                zip_ref.write('b2.neutrals.parameters')
                zip_ref.set_state(state='needs_update')
            if os.path.exists('b2.boundary.parameters'):
                zip_ref.write('b2.boundary.parameters')
                zip_ref.set_state(state='needs_update')
            if os.path.exists('b2.transport.inputfile'):
                zip_ref.write('b2.transport.inputfile')
                zip_ref.set_state(state='needs_update')
            
#  eirene files. We need to rename the eirene input files.
            if os.path.exists(eirene_input_dat):
                os.rename(eirene_input_dat, 'fort.1')
                zip_ref.write('fort.1')
                zip_ref.set_state(state='needs_update')
            if os.path.exists(eirene_grid):
                os.rename(eirene_grid, 'fort30')
                zip_ref.write('fort30')
                zip_ref.set_state(state='needs_update')
            if os.path.exists(eirene_nodes):
                os.rename(eirene_nodes, 'fort.33')
                zip_ref.write('fort.33')
                zip_ref.set_state(state='needs_update')
            if os.path.exists(eirene_cells):
                os.rename(eirene_cells, 'fort.34')
                zip_ref.write('fort.34')
                zip_ref.set_state(state='needs_update')
            if os.path.exists(eirene_links):
                os.rename(eirene_links, 'fort.35')
                zip_ref.write('fort.35')
                zip_ref.set_state(state='needs_update')

        self.services.update_state()

#-------------------------------------------------------------------------------
#
#  SOLPS-ITER_init init Component step method. Not used.
#
#-------------------------------------------------------------------------------
        def step(self, timeStamp=0.0):
            ScreenWriter.screen_output(self, 'verbose', 'solps_iter_init: step')
    
#-------------------------------------------------------------------------------
#
#  SOLPS-ITER_init init Component finalize method. This cleans up afterwards.
#  Not used.
#
#-------------------------------------------------------------------------------
        def finalize(self, timeStamp=0.0):
            ScreenWriter.screen_output(self, 'verbose', 'solps_iter_init: finalize')
