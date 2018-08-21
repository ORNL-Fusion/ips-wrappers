#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for SOLPS-ITER component. This wapper runs b2-eierne then
#  computes synthetic signals from the output.
#
#-------------------------------------------------------------------------------

from component import Component
import os
from omfit.classes.omfit_namelist import OMFITnamelist
from utilities import ZipState

#-------------------------------------------------------------------------------
#
#  SOLPS-ITER Component Constructor
#
#-------------------------------------------------------------------------------
class solps_iter(Component):
    def __init__(self, services, config):
        print('solps_iter: Construct')
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  SOLPS-ITER Component init method. This method prepairs the input files.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('solps_iter: init')
        
#  Set the top of the SOLPS source tree as an environment variable.
        if timeStamp = 0.0:
            self.current_solps_state = services.get_config_param('CURRENT_SOLPS_STATE')
            os.environ['SOLPSTOP'] = self.services.get_config_param('SOLPSTOP')
        
        self.services.stage_plasma_state()
        
        self.zip_ref = zipfile.ZipFile(self.current_solps_state, 'r')
        zip_ref.extractall()

#  Rename state file to inital state if it exists. Then delete the state file.

#  Create eirene symbolic links.
        if timeStamp = 0.0:
            eirene_database_path = self.services.get_config_param('EIRENE_DATABASE_PATH')
            os.symlink(os.path.join(eirene_database_path, 'graphite_ext.dat'), 'graphite_ext.dat')
            os.symlink(os.path.join(eirene_database_path, 'mo_ext.dat'), 'mo_ext.dat')
            os.symlink(os.path.join(eirene_database_path, 'AMJUEL'), 'AMJUEL')
            os.symlink(os.path.join(eirene_database_path, 'H2VIBR'), 'H2VIBR')
            os.symlink(os.path.join(eirene_database_path, 'HYDHEL'), 'HYDHEL')
            os.symlink(os.path.join(eirene_database_path, 'METHANE'), 'METHANE')
            os.symlink(os.path.join(eirene_database_path, 'PHOTON'), 'PHOTON')
            os.symlink(os.path.join(eirene_database_path, 'Surfacedata', 'SPUTTER'), 'SPUTTER')
            os.symlink(os.path.join(eirene_database_path, 'Surfacedata', 'TRIM', 'trim.dat'), 'fort.21')
            os.symlink(os.path.join(eirene_database_path, 'Surfacedata', 'TRIM', 'marlow.dat'), 'fort.22')

#  Update parameters in the namelist.
        self.set_namelists(**keywords)

#-------------------------------------------------------------------------------
#
#  SOLPS-ITER step Component step method. This runs b2.eirene and computes the
#  synthetic signals.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0, **keywords):
        print('solps_iter: step')

        flags = self.zip_ref.get_state()

        if 'state' in flags and flags['state'] == 'needs_update':
        
            b2_eirene_task_id = self.services.launch_task(self.NPROC,
                                                          self.services.get_working_dir(),
                                                          self.B2_EIRENE_EXE,
                                                          logfile = 'b2.eirene.log')

#  Update flags.
            self.zip_ref.set_state(state='updated')

#  Wait for SOLPS to finish.
            if self.services.wait_task(b2_eirene_task_id) :
                self.services.error('solps_iter: step failed')



        else:
            self.zip_ref.set_state(state='unchanged')

        if 'task' in keywords and 'reconstruct' in keywords['task'] :
            solps_signals_task_id = self.services.launch_task(self.NPROC,
                                                              self.services.get_working_dir(),
                                                              self.SOLPS_SIGNALS_EXE,
                                                              '-task=get_result',
                                                              '-geometry=b2fgmtry',
                                                              '-state=b2fstate',
                                                              '-diag_geometry=' + self.services.get_config_param('DIAGNOSTIC_GEOMETRY'),
                                                              '-diag_state=' + self.services.get_config_param('DIAGNOSTIC_STATE'),
                                                              '-dakota_result=' + os.path.join(self.services.get_config_param('SIM_ROOT'), 'RESULT'),
                                                              logfile = 'solps_signals.log')
                
            if self.services.wait_task(solps_signals_task_id) :
                self.services.error('solps_iter: step failed')

#-------------------------------------------------------------------------------
#
#  SOLPS-ITER Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('solps_iter: finalize')
                     
#-------------------------------------------------------------------------------
#
#  SOLPS-ITER Component set_namelists method. This sets the namelist input files
#  from the keywords.
#
#-------------------------------------------------------------------------------
    def set_namelists(self, **keywords):
#  Replace parameters in the plasma state.
        transport_nl = OMFITnamelist('b2.transport.parameters')
        numerics_nl = OMFITnamelist('b2.numerics.parameters')
        neutrals_nl = OMFITnamelist('b2.neutrals.parameters')
        boundary_nl = OMFITnamelist('b2.boundary.parameters')
        
        if len(keywords) > 0:
            self.zip_ref.set_state(state='needs_update')

#  Keywords will contain the following syntax.
#
#      file__param
#
#  The file is a two letter indicating which namelist input to update.
#
#    tr : b2.transport.parameters
#    nu : b2.numerics.parameters
#    ne : b2.neutrals.parameters
#    bo : b2.boundary.parameters
#
#  The namelist item is specified verbatium. For arrays, the array index is
#  specifed in the same manner as an entry in the namelist input file would be.
#
#  EXAMPLE:
#
#    name_of_component__dakota_tr_parm_hci(1)
#
            for key, value in keywords.iteritems():
                if 'tr' in key:
                    transport_nl['transport'][key.replace('tr__','',1)] = value
                if 'nu' in key:
                    transport_nl['numerics'][key.replace('nu__','',1)] = value
                if 'ne' in key:
                    transport_nl['NEUTRALS'][key.replace('ne__','',1)] = value
                if 'bo' in file :
                    transport_nl['boundary'][key.replace('bo__','',1)] = value

        transport_nl.save()
        numerics_nl.save()
        neutrals_nl.save()
        boundary_nl.save()
