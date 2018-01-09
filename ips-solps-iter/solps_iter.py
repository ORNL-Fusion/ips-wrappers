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
import zipfile

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
        os.environ['SOLPSTOP'] = self.services.get_config_param('SOLPSTOP')
        
        self.services.stage_plasma_state()
        
        zip_ref = zipfile.ZipFile(self.services.get_config_param('CURRENT_SOLPS_STATE'), 'r')
        zip_ref.extractall()
        zip_ref.close()

#  Replace parameters in the plasma state.
        transport_nl = OMFITnamelist('b2.transport.parameters')
        numerics_nl = OMFITnamelist('b2.numerics.parameters')
        neutrals_nl = OMFITnamelist('b2.neutrals.parameters')
        boundary_nl = OMFITnamelist('b2.boundary.parameters')

#  Dakota parameters use the following syntax. To pass a parameter to a
#  component, Dakota parameters are prefixed with the name of the component
#  followed by two underscores. To quickly identifiy the dakota parameters a
#  dakota prefix followed by a single underscore is appended to the parameter.
#
#      component__dakota_param
#
#  At the component level, everything up to the two underscores is removed. For
#  b2/eirene, there are multiple files. Files will be indicated by a two letter
#  code followed by a single underscore:
#
#    tr      :  b2.transport.parameters
#    nu      :  b2.numerics.parameters
#    ne      :  b2.neutrals.parameters
#    bo      :  b2.boundary.parameters
#
#  The namelist item is specified verbatium. For arrays, the array index is
#  specifed in the same manner as an entry in the namelist input file would be.
#
#  EXAMPLE:
#
#    name_of_component__dakota_tr_parm_hci(1)
#
        for key, value in self.__dict__.iteritems() :
            if 'dakota_' in key :
                extra, key = key.split('_',1)
                file, key = key.split('_',1)

                if '(' in key :
                    key, indices = key.split('(')
                    indices, extra = indices.split(')')
                    indices = [[int(i) - 1] for i in indices.split(',')]

                    if 'tr' in file :
                        transport_nl['transport'][key][indices] = float(value)

                    if 'nu' in file :
                        numerics_nl['numerics'][key][indices] = float(value)

                    if 'ne' in file :
                        neutrals_nl['NEUTRALS'][key][indices] = float(value)

                    if 'bo' in file :
                        boundary_nl['boundary'][key][indices] = float(value)
                else :
                    if 'tr' in file :
                        transport_nl['transport'][key] = float(value)
                    
                    if 'nu' in file :
                        numerics_nl['numerics'][key] = float(value)
                    
                    if 'ne' in file :
                        neutrals_nl['NEUTRALS'][key] = float(value)
                    
                    if 'bo' in file :
                        boundary_nl['boundary'][key] = float(value)

        transport_nl.save()
        numerics_nl.save()
        neutrals_nl.save()
        boundary_nl.save()

#  Create eirene symbolic links.
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
    
#-------------------------------------------------------------------------------
#
#  SOLPS-ITER step Component step method. This runs b2.eirene and computes the
#  synthetic signals.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0, **keywords):
        print('solps_iter: step')

        b2_eirene_task_id = self.services.launch_task(self.NPROC,
                                                      self.services.get_working_dir(),
                                                      self.B2_EIRENE_EXE,
                                                      logfile = 'b2.eirene.log')

        if self.services.wait_task(b2_eirene_task_id) :
            self.services.error('solps_iter: step failed')

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
