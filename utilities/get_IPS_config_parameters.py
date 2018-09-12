#! /usr/bin/env python

"""
A simple wrapper to get global and component specific configuration parameters from an
IPS configuration file with one call.  Wraps exception handling, handling of optional
config parameters, and optional printing of parameter name/value.  For optional parameters
on call set optional=True, to turn off printing set verbose=False.  'self' must be passed
as an explicit argument.  e.g.
X = config.get_component_param(self, services, 'X', optional=True, verbose=False)

"""

    # Try to get config parameter - wraps the exception handling for get_config_parameter()
    def get_global_param(self, services, param_name, optional=False):

        try:
            value = services.get_config_param(param_name)
            print param_name, ' = ', value
        except Exception:
            if optional: 
                print 'optional config parameter ', param_name, ' not found'
                value = None
            else:
                message = 'required config parameter ', param_name, ' not found'
                print message
                services.exception(message)
                raise
        
        return value

    # Try to get component specific config parameter - wraps the exception handling
    def get_component_param(self, services, param_name, optional=False):

        if hasattr(self, param_name):
            value = getattr(self, param_name)
            print param_name, ' = ', value
        elif optional:
            print 'optional config parameter ', param_name, ' not found'
            value = None
        else:
            message = 'required component config parameter ', param_name, ' not found'
            print message
            services.exception(message)
            raise
        
        return value
