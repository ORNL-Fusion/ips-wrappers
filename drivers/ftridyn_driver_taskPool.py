#! /usr/bin/env python

from  component import Component

class ftridynDriver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):
        print('ftridyn_driver: init')

        current_ftridyn_namelist = self.services.get_config_param('FTRIDYN_INPUT_FILE')
        #split filenames into a list of strings
        file_list = current_ftridyn_namelist.split()
        #loop over file names and create dummy files in ftridynInit work area
        for index in range(len(file_list)):
            open(file_list[index], 'a').close()

        driver_out = self.services.get_config_param('EA_OUTPUT')
        fid = open(driver_out,'w')
        fid.write('Energy   Angle   Roughness   Sputtering Yield\n')
        fid.close()

        #update plasma state from relevant files in ftridynInit work area
        self.services.update_plasma_state()        

        #Get info from config file about port 'WORKER'
        self.ftridyn_comp = self.services.get_port('WORKER')

    def step(self, timeStamp=0.0):
        print('ftridyn_driver: step')
        energy = [float(i) for i in self.ENERGY.split()]
        angle = [float(i) for i in self.ANGLE.split()]
	roughness = [float(i) for i in self.ROUGHNESS.split()]
        

        #call init method of WORKER ('ftridyn_worker.py')
        self.services.call(self.ftridyn_comp, 'init', timeStamp, eArg =energy,aArg = angle,dArg = roughness)
                    
        #call step method of WORKER ('ftridyn_worker.py')
        self.services.call(self.ftridyn_comp, 'step', timeStamp,eArg =energy,aArg = angle,dArg = roughness )
    
    def finalize(self, timeStamp=0.0):
        print('ftridyn_driver: finalize')
        self.services.call(self.ftridyn_comp, 'finalize', timeStamp)
