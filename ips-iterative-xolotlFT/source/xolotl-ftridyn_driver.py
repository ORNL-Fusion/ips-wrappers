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

        #start first loop from the beginning (INIT) or from a previous run (RESTART)
        #RESTART mode requires providing a list of input files:
        #for FTridyn: last_TRIDYN.dat; for Xolotl: params.txt (of the last run), networkfile (networkRestart.h5)
        #and placing them in the 'restart_files' folder. The mode is changed to NEUTRAL after the 1st loop
        self.startMode = 'INIT'
        self.driverMode=self.startMode

        self.initTime=0.0
        self.endTime=0.2
        self.timeStep=0.02

        print 'running IPS from t = %f to t=%f, in steps of dt=%f' % (self.initTime, self.endTime, self.timeStep)

        #Xolotl parameter:
        #write every parameter that will be used as argumennts in write_xolotl_paramfile function(s) 
        #i.e., anything different from default values (those set to reproduce email-coupling of FTridyn-Xolotl)

        self.xolotlStartStop='True'
        self.xolotlFlux=4e4
        if self.startMode=='INIT':
            self.xolotlNetworkFile='networkInit.h5'
        elif self.startMode=='RESTART':
            self.xolotlNetworkFile='networkRestart.h5'

        #ftridyn parameters:
        #TotalDepth: total substrate depth in [A]; set to 0.0 to use what Xolotl passes to ftridyn (as deep as He exists)
        #InitialTotalDepth: if TotalDepth=0.0, choose an appropriate depth for the irradiation energy in the 1st loop
        #     use TotalDepth=0.0 if startMode is RESTART (not understood why, but a fixed totalDepth doesn't work on the 1st loop)
        #NImpacts: number of impacts (NH in generateInput) ;  InEnergy: impact energy (energy in generateInput, [eV]); initialize SpYield
        #if spYield < 0 -> use calculated value; else, use fixed value, usually [0,1) 
            
        self.ftridynTotalDepth=0.0
        self.ftridynInitialTotalDepth=300.0
        self.ftridynNImpacts=1e5
        self.ftridynInEnergy=250.0
        self.ftridynSpYieldW=-1.0

        if self.ftridynSpYieldW<0:
            self.ftridynSpYieldMode='calculate'
        else:
            self.ftridynSpYieldMode='fixed'

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

        for time in numpy.arange(self.initTime,self.endTime,self.timeStep):

            self.services.stage_plasma_state()
            print 'driver time (in loop)  %f' %(time)
            self.services.update_plasma_state()
            
            #component/method calls now include arguments (variables)
            self.services.call(ftridyn, 'init', timeStamp, dMode=self.driverMode, dTime=time, fInitialTotalDepth=self.ftridynInitialTotalDepth, fTotalDepth=self.ftridynTotalDepth, fNImpacts=self.ftridynNImpacts, fEnergyIn=self.ftridynInEnergy)
            self.services.call(ftridyn, 'step', timeStamp, dTime=time, fNImpacts=self.ftridynNImpacts)#, fSpYieldW=self.ftridynSpYieldW)

            self.services.call(xolotl, 'init', timeStamp, dStartMode=self.startMode, dMode=self.driverMode, dTime=time, dTimeStep=self.timeStep, xNetworkFile=self.xolotlNetworkFile, xStartStop=self.xolotlStartStop, xFlux=self.xolotlFlux, fNImpacts=self.ftridynNImpacts, fSpYieldMode=self.ftridynSpYieldMode, fSpYieldW=self.ftridynSpYieldW)
            self.services.call(xolotl, 'step', timeStamp, dTime=time)

            self.services.stage_plasma_state()

            #update driver mode after the 1st loop, from INIT to RESTART
            if self.driverMode == 'INIT':
                self.driverMode = 'RESTART'
                print 'switched driverMode to ', self.driverMode

            if self.startMode != 'NEUTRAL':
                self.startMode = 'NEUTRAL'
                print 'switched startMode to ',  self.startMode 

            self.xolotlNetworkFile='xolotlStop_%f.h5' %(time)

            self.services.update_plasma_state()

    def finalize(self, timeStamp=0.0):
        print('xolotl-ftridyn_driver: finalize')
        ftridyn = self.services.get_port('WORKER')
        xolotl = self.services.get_port('XWORKER')
        self.services.call(ftridyn, 'finalize', timeStamp)
        self.services.call(xolotl, 'finalize', timeStamp)
