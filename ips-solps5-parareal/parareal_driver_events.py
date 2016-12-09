# imports need to be cleaned up
from component import Component
#from Scientific.IO.NetCDF import *
#import Numeric


class Driver(Component):

    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

# ------------------------------------------------------------------------------
#
# init function
#
# ------------------------------------------------------------------------------

    def init(self, timestamp=0):
        # Driver initialization ? nothing to be done
        return

# ------------------------------------------------------------------------------
#
# step function
#
# ------------------------------------------------------------------------------

    def step(self, timestamp=0):

        services = self.services

       # Get references to the components to run in simulation
        try:
            fine = services.get_port('FINE')
        except Exception:
            self.services.exception('Error accessing fine component')
            raise
        try:                                    
            coarse = services.get_port('COARSE')
        except  Exception:
            self.services.exception('Error accessing coarse component')
            raise 
        try:
            converge = services.get_port('CONVERGE')
        except Exception:
            self.services.exception('Error accessing converge component')
            raise


        # Get timeloop for simulation
        timeloop = services.get_time_loop()
        tlist_str = ['%.3f'%t for t in timeloop]
        t = tlist_str[0]

        # Call init for each component
        services.call(fine,'init', t)
        print 'fine init called'
        print (' ')

        services.call(coarse,'init', t)
        print 'coarse init called'
        print (' ')

        services.call(converge,'init', t)
        print 'converge init called'
        print (' ')
        # get the number of slices for the simulation
        n_slice = int(self.services.get_config_param('NT_SLICE'))
        print ' init sequence complete--ready for time loop'
        # Iterate through the timeloop
        # make an iteration counter so that we don't have to rely
        # on time loop conversion
        print 'Driver:  start iteration = ', t
        services.updateTimeStamp(t)
        print (' ')

        # Call step for each component

        print (' ')
        coars_callid = services.call_nonblocking(coarse,'step', t)

        print (' ')
        fine_callid = services.call_nonblocking(fine,'step', t)

        print (' ')
        print 'calling converge'
        converge_callid = services.call_nonblocking(converge,'step', t)
        
        call_retval = services.wait_call_list([coars_callid, fine_callid, converge_callid])
        
        services.call(coarse , 'finalize')
        services.call(fine, 'finalize')
        services.call(converge, 'finalize')
        return

# ------------------------------------------------------------------------------
#
# finalize function
#
# ------------------------------------------------------------------------------

    def finalize(self, timestamp=0.0):
        # Driver finalize - nothing to be done
        pass
