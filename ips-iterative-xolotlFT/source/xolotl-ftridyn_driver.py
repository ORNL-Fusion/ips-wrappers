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

        #A MORE ELEGANT WAY --  FOR THE FUTURE
        #for file in plasma_state_list:
        #    open(file, 'a').close()

        self.services.update_plasma_state()
        self.services.stage_plasma_state()

        #start first loop from scratch (INIT) or from a previous run (RESTART)
        #RESTART mode requires providing a list of input files, and placing them in a 'restart_files' folder
        #for FTridyn: last_TRIDYN.dat; for Xolotl: params.txt (of the last run), networkfile (networkRestart.h5)
        #the mode is changed to NEUTRAL after the 1st loop
        self.startMode = 'INIT'

        self.initTime=0.0
        self.endTime=0.2
        self.timeStep=0.1

#        driverInitTimeString='driverInitTime = %f' %driverInitTime
#        driverEndTimeString='driverEndTime = %f' %driverEndTime
#        driverTimeString='driverTime = %f' %driverInitTime
#        driverTimeStepString='driverTimeStep = %f' %driverTimeStep

#        if driverStartMode=='INIT':
#            driverStartModeLine="driverStartMode='INIT'"
#        elif driverStartMode=='RESTART':
#            driverStartModeLine="driverStartMode='RESTART'"
#        elif driverStartMode=='NEUTRAL':
#            driverStartModeLine="driverStartMode='NEUTRAL'"

        print 'running IPS from t = %f to t=%f, in steps of dt=%f' % (self.initTime, self.endTime, self.timeStep)
        
        #get config file and write drivers initial state 
#        if driverStartMode=='INIT':
#            modeLine="mode='INIT'"
#        elif driverStartMode=='RESTART':
#            modeLine="mode='RESTART'"

        self.driverMode=self.startMode

#        driver_config_file = self.services.get_config_param('PARAMETER_CONFIG_FILE')
#        xftid = open(driver_config_file, 'w')
#        xftid.write("%s \n%s \n%s \n%s \n%s \n%s \n" % (modeLine, driverStartModeLine, driverInitTimeString, driverEndTimeString, driverTimeString, driverTimeStepString))
#        xftid.close()

        #get config file and write Xolotl's parameter file options
        #write every parameter that will be passed in write_xolotl_paramfile functions 
        #i.e., different from default values (set to reproduce email-coupling of FTridyn-Xolotl)
#        xolotlStartStopLine="start_stop='True'"
#        flux=0.25e4
#        xolotlFluxLine="flux=%e" %flux
#        if driverStartMode=='INIT':
#            xolotlNetworkFileLine="networkFile='networkInit.h5'"
#        elif driverStartMode=='RESTART':
#            xolotlNetworkFileLine="networkFile='networkRestart.h5'"

#        xolotl_config_file = self.services.get_config_param('XOLOTL_PARAMETER_CONFIG_FILE')
#        xid = open(xolotl_config_file, 'w')
#        xid.write("%s \n%s \n%s \n" % (xolotlStartStopLine, xolotlNetworkFileLine,xolotlFluxLine))
#        xid.close()

        self.xolotlStartStop='True'
        self.xolotlFlux=4e4
        if self.startMode=='INIT':
            self.xolotlNetworkFile='networkInit.h5'
        elif self.startMode=='RESTART':
            self.xolotlNetworkFile='networkRestart.h5'

        #ftridyn config file
        #TotalDepth: total substrate depth in [A]; set to 0.0 to use what Xolotl passes to ftridyn (as deep as He exists)
        #InitialTotalDepth: if TotalDepth=0.0, choose an appropriate depth for the irradiation energy in the 1st loop
        #NImpacts: number of impacts (NH in generateInput) ;  InEnergy: impact energy (energy in generateInput, [eV]); initialize SpYield
#        ftridynTotalDepth=0.0
#        ftridynInitialTotalDepth=200.0
#        ftridynNImpacts=1e5
#        ftridynInEnergy=250.0
#        ftridynSpYieldW=0.0

#        driverFtridynTotalDepthString='totalDepth = %f' %ftridynTotalDepth
#        driverFtridynInitialTotalDepthString='initialTotalDepth = %f' %ftridynInitialTotalDepth
#        driverFtridynNImpactsString='nImpacts = %f' %ftridynNImpacts
#        driverFtridynInEnergy='energyIn=%f' %ftridynInEnergy
#        driverFtridynSpYieldString='spYieldW=%f' %ftridynSpYieldW

#        ftridyn_config_file = self.services.get_config_param('FTRIDYN_PARAMETER_CONFIG_FILE')
#        ftid = open(ftridyn_config_file, 'w')
#        ftid.write("%s \n%s \n%s \n%s \n%s \n" % (driverFtridynTotalDepthString,driverFtridynInitialTotalDepthString,driverFtridynNImpactsString,driverFtridynInEnergy,driverFtridynSpYieldString))
#        ftid.close()

        self.ftridynTotalDepth=0.0
        self.ftridynInitialTotalDepth=200.0
        self.ftridynNImpacts=1e5
        self.ftridynInEnergy=200.0
        self.ftridynSpYieldW=0.0

        #if driver start mode is in Restart, copy files to plasma_state
        #MAYBE THIS CAN ALSO BE WRITTEN MORE ELEGANTLY
        if (self.startMode=='RESTART'):
            restart_files = self.services.get_config_param('RESTART_FILES')
            restart_list = restart_files.split()
            for index in range(len(restart_list)):
                filepath='../../restart_files/'+restart_list[index]
                shutil.copyfile(filepath,restart_list[index])

        self.services.update_plasma_state()

    def step(self, timeStamp=0.0):
        print('xolotl-ftridyn_driver: step')

        ftridyn = self.services.get_port('WORKER')
        xolotl = self.services.get_port('XWORKER')

        self.services.stage_plasma_state() 

#        sys.path.append(os.getcwd())

#        driver_config_file = self.services.get_config_param('PARAMETER_CONFIG_FILE')
#        import driverParameterConfig
#        reload(driverParameterConfig)

#        xolotl_config_file = self.services.get_config_param('XOLOTL_PARAMETER_CONFIG_FILE')
#        import xolotlParameterConfig
#        reload(xolotlParameterConfig)

#        ftridyn_config_file = self.services.get_config_param('FTRIDYN_PARAMETER_CONFIG_FILE')
#        import ftridynParameterConfig
#        reload(ftridynParameterConfig)

        for time in numpy.arange(self.initTime,self.endTime,self.timeStep):

            self.services.stage_plasma_state()

            print 'driver time (in loop)  %f' %(time)
#            reload(driverParameterConfig)
#            driver_tmp_file='configFileRestart.tmp'
#            shutil.copyfile(driver_config_file,driver_tmp_file)
#            timeSedString="sed    -e 's/driverTime = [^ ]*/driverTime = %f/' <%s >%s "   % (time, driver_tmp_file, driver_config_file)
#            subprocess.call([timeSedString], shell=True)
#            os.remove(driver_tmp_file)

#            self.services.update_plasma_state()
            
            #component/method calls now include arguments (variables)
            self.services.call(ftridyn, 'init', timeStamp, dMode=self.driverMode, dTime=time, fInitialTotalDepth=self.ftridynInitialTotalDepth, fTotalDepth=self.ftridynTotalDepth, fNImpacts=self.ftridynNImpacts, fEnergyIn=self.ftridynInEnergy)
            self.services.call(ftridyn, 'step', timeStamp, dTime=time, fNImpacts=self.ftridynNImpacts, fSpYieldW=self.ftridynSpYieldW)

            self.services.call(xolotl, 'init', timeStamp, dStartMode=self.startMode, dMode=self.driverMode, dTime=time, dTimeStep=self.timeStep, xNetworkFile=self.xolotlNetworkFile, xStartStop=self.xolotlStartStop, xFlux=self.xolotlFlux, fNImpacts=self.ftridynNImpacts) #fSpYieldW=self.ftridynSpYieldW
            self.services.call(xolotl, 'step', timeStamp, dTime=time)

#            self.services.stage_plasma_state() 
            
#            import parameterConfig
#            reload(driverParameterConfig)
#            print 'reading driver parameter config mode %s' %(driverParameterConfig.mode)

#            if driverParameterConfig.mode == 'INIT':
#                driver_tmp_file='driverConfigFileRestart.tmp'
#                shutil.copyfile(driver_config_file,driver_tmp_file)
#                modeSedString="sed    -e 's/mode=[^ ]*/mode="'"RESTART"'"/' <%s >%s "   % (driver_tmp_file, driver_config_file)
#                subprocess.call([modeSedString], shell=True)
#                os.remove(driver_tmp_file)

            #update driver mode after the 1st loop, from INIT to RESTART
            if self.driverMode == 'INIT':
                self.driverMode = 'RESTART'
                print 'switched driverMode to ', self.driverMode

#            if driverParameterConfig.driverStartMode != 'NEUTRAL':
#                driver_tmp_file='driverConfigFileRestart.tmp'
#                shutil.copyfile(driver_config_file,driver_tmp_file)
#                modeSedString="sed    -e 's/driverStartMode=[^ ]*/driverStartMode="'"NEUTRAL"'"/' <%s >%s "   % (driver_tmp_file, driver_config_file)
#                subprocess.call([modeSedString], shell=True)
#                os.remove(driver_tmp_file)

            if self.startMode != 'NEUTRAL':
                self.startMode = 'NEUTRAL'
                print 'switched startMode to ',  self.startMode 

#            reload(xolotlParameterConfig)
#            print 'updating xolotl parameter config network File from %s' %(xolotlParameterConfig.networkFile)

#            networkFileLine='xolotlStop_%f.h5' %(driverParameterConfig.driverTime)
#            xolotl_tmp_file='xolotlConfigFile.tmp'
#            shutil.copyfile(xolotl_config_file,xolotl_tmp_file)
#            networkFileSedString="sed    -e 's/networkFile=[^ ]*/networkFile="'"%s"'"/' <%s >%s "   %(networkFileLine,xolotl_tmp_file,xolotl_config_file)
#            subprocess.call([networkFileSedString], shell=True)
#            os.remove(xolotl_tmp_file)

            self.xolotlNetworkFile='xolotlStop_%f.h5' %(time)


            #test updating network file
#            self.services.update_plasma_state()
            self.services.stage_plasma_state()
#            reload(xolotlParameterConfig)
#            print '\t \t \t \t to network File %s' %(xolotlParameterConfig.networkFile)

            self.services.update_plasma_state()

    def finalize(self, timeStamp=0.0):
        print('xolotl-ftridyn_driver: finalize')
        ftridyn = self.services.get_port('WORKER')
        xolotl = self.services.get_port('XWORKER')
        self.services.call(ftridyn, 'finalize', timeStamp)
        self.services.call(xolotl, 'finalize', timeStamp)
