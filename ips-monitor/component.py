class Component(object):
    """
    A dummy Component class so monitor_comp.py can be used by monitor_comp_stand_alone.py
    without having the whole IPS framework present.  This file must be in the same directory
    with monitor_comp_stand_alone.py, which also must be in the same directory as the IPS 
    monitor component which is being driven by the stand alone code.
    """
    def __init__(self, services, config):
        """
        Set up config values and reference to services.
        """
        self.services = services
        self.config = config

