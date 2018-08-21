#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for SOLPS-ITER Driver component. This wapper runs b2-eierne then
#  computes synthetic signals from the output.
#
#-------------------------------------------------------------------------------

from component import Component

#-------------------------------------------------------------------------------
#
#  SOLPS-ITER Driver Component Constructor
#
#-------------------------------------------------------------------------------
class solps_iter_driver(Component):
    def __init__(self, services, config):
        print('solps_iter_driver: Construct')
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  SOLPS-ITER Driver Component init method.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0, **keywords):
        print('solps_iter_driver: init')

#  Separate out the vmec keywords.
        solps_keywords = {}
        for key, value in keywords.iteritems():
            if 'solps__' in key:
                solps_keywords[key.replace('solps__','',1)] = value

        self.solps_port = self.services.get_port('SOLPS')
        self.wait = self.services.call_nonblocking(self.solps, 'init',
                                                   timeStamp, **solps_keywords)

#-------------------------------------------------------------------------------
#
#  SOLPS-ITER Driver step Component step method.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0, **keywords):
        print('solps_iter_driver: step')

        current_solps_state = self.services.get_config_param('CURRENT_SOLPS_STATE')

        if not 'task' in keywords:
            keywords['task'] = self.SOLPS_TASK

#  Run SOLPS.
        self.services.wait_call(self.wait, True)
        self.services.call(self.solps_port, 'step', timeStamp, **keywords)

#  Prepare the output files for a super work flow. Need to remove any old output
#  files first before staging the plasma state.
        if os.path.exists(self.OUTPUT_FILES):
            os.remove(self.OUTPUT_FILES)
        self.services.stage_plasma_state()

#  The super flow may need to rename the output file. Check is the current state
#  matches if output file. If it does not rename the plasma state so it can be
#  staged.
        if not os.path.exists(self.OUTPUT_FILES):
            os.rename(current_solps_state, self.OUTPUT_FILES)
    
#-------------------------------------------------------------------------------
#
#  SOLPS-ITER Driver Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('solps_iter_driver: finalize')

        self.services.call(self.solps_port, 'finalize', timeStamp)

