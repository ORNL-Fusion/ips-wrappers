#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for the Massive VMEC Runner Init component.
#
#-------------------------------------------------------------------------------

from ipsframework import Component
from ips_components_utilities import ZipState
from ips_components_utilities import ScreenWriter
import os
import adaptive

#-------------------------------------------------------------------------------
#
#  Check if the file already exists. If it doesn't try to extract it from the
#  zip state file.
#
#-------------------------------------------------------------------------------
def extract_if_needed(zip_ref, file):
    if not os.path.exists(file) and file in zip_ref:
        zip_ref.extract(file)

#-------------------------------------------------------------------------------
#
#  Massive VMEC Runner init Component Constructor
#
#-------------------------------------------------------------------------------
class massive_serial_runner_init(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  Massive VMEC Runner init Compoenet init method. This method prepairs the
#  state.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0, **keywords):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner_init: init')

#  Get config filenames.
        if timeStamp == 0.0:
            self.current_state = self.services.get_config_param('CURRENT_MVR_STATE')
            self.current_batch = self.services.get_config_param('CURRENT_BATCH')
            self.batch_size = self.services.get_config_param('BATCH_SIZE')
            self.constraint_path = self.services.get_config_param('MODULE_PATH')
            self.constraint_name = self.services.get_config_param('MODULE_NAME')
            self.model_config = self.services.get_config_param('MODEL_CONFIG')

#  Remove old inputs.
        if os.path.exists(self.current_batch):
            os.remove(self.current_batch)
        
#  Stage input files and setup intial state.
        self.services.stage_input_files(self.INPUT_FILES)

#  Load or create a massive vmec runner zip state.
        with ZipState.ZipState(self.current_state, 'a') as zip_ref:
#  Check if a new batch of data exists. If it doesn't create a new batch file
#  from scan parameters.
            if self.current_batch not in zip_ref:
                zip_ref.extract(self.model_config)

                model = adaptive.no_model(adaptive.load_json(self.model_config),
                                          self.constraint_path,
                                          self.constraint_name)
                model.create_prediction_data(self.batch_size)
                model.save_random_sample(self.current_batch)

                zip_ref.write(self.current_batch)

        self.services.update_state()

#-------------------------------------------------------------------------------
#
#  Massive VMEC Runner init Component step method. Not used.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner_init: step')

#-------------------------------------------------------------------------------
#
#  Massive VMEC Runner init Component finalize method. This cleans up
#  afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner_init: finalize')

