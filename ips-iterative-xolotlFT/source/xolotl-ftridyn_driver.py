#! /usr/bin/env python

from  component import Component
import sys
import os
import subprocess
import numpy
import shutil

class xolotlFtridynDriver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):
        print('xolotl-ftridyn_driver: init')
        
        plasma_state_file = self.services.get_config_param('PLASMA_STATE_FILES')
        plasma_state_list = plasma_state_file.split()
        for index in range(len(plasma_state_list)):
            open(plasma_state_list[index], 'a').close()

        self.services.update_plasma_state()

    def step(self, timeStamp=0.0):
        print('xolotl-ftridyn_driver: step')

        ftridyn = self.services.get_port('WORKER')
        xolotl = self.services.get_port('XWORKER')

        self.services.stage_plasma_state() 

        driverInitTime=0.0
        driverEndTime=0.0002
        driverTimeStep=0.0001

        driverTimeString='driverTime = %f' %driverInitTime
        driverTimeStepString='driverTimeStep = %f' %driverTimeStep

        print 'running IPS from t = %f to t=%f, in steps of dt=%f' % (driverInitTime, driverEndTime, driverTimeStep)

        #get config file and write initial state
        modeLine="mode='INIT'"
        xolotl_ftridyn_config_file = self.services.get_config_param('PARAMETER_CONFIG_FILE')
        fid = open(xolotl_ftridyn_config_file, 'w')
        fid.write("%s \n%s \n%s \n" % (modeLine, driverTimeString, driverTimeStepString))
        fid.close()       
        
        self.services.update_plasma_state()
        sys.path.append(os.getcwd())
        import parameterConfig

        for time in numpy.arange(driverInitTime,driverEndTime,driverTimeStep):

            self.services.stage_plasma_state()

            print 'driver time (in loop)  %f' %(time)
            reload(parameterConfig)
            xolotl_ftridyn_tmp_file='configFileRestart.tmp'
            shutil.copyfile(xolotl_ftridyn_config_file,xolotl_ftridyn_tmp_file)
            timeSedString="sed    -e 's/driverTime = [^ ]*/driverTime = %f/' <%s >%s "   % (time, xolotl_ftridyn_tmp_file, xolotl_ftridyn_config_file)
            subprocess.call([timeSedString], shell=True)
            os.remove(xolotl_ftridyn_tmp_file)

            self.services.update_plasma_state()

            self.services.call(ftridyn, 'init', timeStamp)
            self.services.call(ftridyn, 'step', timeStamp)
            
            self.services.call(xolotl, 'init', timeStamp)
            self.services.call(xolotl, 'step', timeStamp)    
            
            self.services.stage_plasma_state() 
            
#            import parameterConfig
            reload(parameterConfig)
            print 'reading parameter config mode %s' %(parameterConfig.mode)

            if (parameterConfig.mode == 'INIT'):
                #test print
                #parameterConfig.mode = 'RESTART'
                #print 'changing parameter config mode to %s ' %(parameterConfig.mode)
                
                xolotl_ftridyn_tmp_file='configFileRestart.tmp'
                shutil.copyfile(xolotl_ftridyn_config_file,xolotl_ftridyn_tmp_file)
                modeSedString="sed    -e 's/mode=[^ ]*/mode="'"RESTART"'"/' <%s >%s "   % (xolotl_ftridyn_tmp_file, xolotl_ftridyn_config_file)
                subprocess.call([modeSedString], shell=True)
                os.remove(xolotl_ftridyn_tmp_file)

                #test 
                #self.services.update_plasma_state()
                #self.services.stage_plasma_state()
                #reload(parameterConfig)
                #print 'parameter config mode changed to: %s' %(parameterConfig.mode)
 
            self.services.update_plasma_state()

    def finalize(self, timeStamp=0.0):
        print('xolotl-ftridyn_driver: finalize')
        self.services.call(ftridyn, 'finalize', timeStamp)
