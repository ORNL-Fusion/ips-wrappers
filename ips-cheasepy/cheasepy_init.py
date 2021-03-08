#! /usr/bin/env python

#-------------------------------------------------------------------------------
#  IPS wrapper for CHEASEPY_init component. Take the work flow inputs and generates a state.
#-------------------------------------------------------------------------------

import os

import cheasepy

from component import Component
from utilities import ZipState
from utilities import ScreenWriter

#-------------------------------------------------------------------------------
#  CHEASEPY_init Component Constructor
#-------------------------------------------------------------------------------
class cheasepy_init(Component):
      def __init__(self, services, config):
          Component.__init__(self, services, config)
          print('created %s' % (self.__class__))

#-------------------------------------------------------------------------------
#  CHEASEPY_init Component init method.
#-------------------------------------------------------------------------------
      def init(self, timeStamp=0):
          ScreenWriter.screen_output(self, 'verbose', 'cheasepy_init(): started')

          try:
              init_run = int(self.INIT_RUN)
          except:
              init_run = 0
          print('init_run', init_run)

          if init_run:
             self.step(-1)

          ScreenWriter.screen_output(self, 'verbose', 'cheasepy_init(): done')

#-------------------------------------------------------------------------------
#  CHEASEPY_init init Component step method.
#-------------------------------------------------------------------------------
      def step(self, timeStamp):
          ScreenWriter.screen_output(self, 'verbose', 'cheasepy_step(): started')

          # CODE ENTRY
          services = self.services

          input_dir_id = getattr(self, "INPUT_DIR_ID", "")
          if input_dir_id:
             services.component_ref.config['INPUT_DIR'] = services.component_ref.config['INPUT_DIR']+"_%d"%int(float(input_dir_id))

          services.stage_input_files(self.INPUT_FILES)

          # PLASMA STATE FILE NAMES
          cur_state_file = services.get_config_param('CURRENT_STATE')
          cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')


          # UPDATE PLASMA STATE
          services.update_plasma_state()

          # ARCHIVING OUTPUT FILES
          services.stage_output_files(timeStamp, self.OUTPUT_FILES)

#-------------------------------------------------------------------------------
#  CHEASEPY_init init Component finalize method. This cleans up afterwards.
#-------------------------------------------------------------------------------
      def finalize(self, timeStamp=0.0):
          ScreenWriter.screen_output(self, 'verbose', 'cheasepy_init: finalize')
