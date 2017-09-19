#! /usr/bin/env python

from  component import Component
import os
import shutil
import subprocess
import glob
#import translate_xolotl_to_ftridyn
import binTRIDYN
import write_xolotl_paramfile
import sys
import numpy as np

class xolotlWorker(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)
        

    def init(self, timeStamp=0.0, **keywords):
        print('xolotl_worker: init')
        self.services.stage_plasma_state()

        print 'check that all arguments are read well by xolotl-init' 
        for (k, v) in keywords.iteritems():
            print '\t', k, " = ", v

        #asign a local variable to arguments used multiple times 
        driverTime=keywords['dTime']
        driverMode=keywords['dMode']
        startStop=keywords['xStartStop']
        networkFile=keywords['xNetworkFile']
#        spYieldW=keywords['fSpYieldW']
        flux=keywords['xFlux']
        runEndTime=driverTime+keywords['dTimeStep']

        fluxFractionW=keywords['gFractionW']

        cwd = self.services.get_working_dir()

        print 'xolotl-init:'
        print '\t driver mode is', driverMode
        print '\t running starts at time',driverTime
        print '\t \t  ends at time', runEndTime
        print '\t driver step is', keywords['dTimeStep']
        print '\n'

        #if fSpYield in driver < 0, fSpYieldMode = calculate -> get sputtering yield FROM spYield.out
        if keywords['fSpYieldMode']=='calculate':
            spYieldsTemp=keywords['spYieldsFile_temp']
            spYieldHe, spYieldW=angleValue, weightAngle = np.loadtxt(cwd+'/'+spYieldsTemp, usecols = (1,2) , unpack=True)
            print '\t W sputtering yields calculated by FTridyn, read from file: spY (by He) = ', spYieldHe, ' spY (by W) = ', spYieldW
        #if fSpYield in driver >= 0, fSpYieldMode = fixed 
        elif keywords['fSpYieldMode']=='fixed':
            spYieldW=keywords['fSpYieldW']
            spYieldHe=keywords['fSpYieldHe']
            print '\t Fixed value of W sputtering yields are: spY (by He) = ', spYieldHe, ' spY (by W) = ', spYieldW
        else:
            print '\t Invalid value of fSpYieldMode, ', keywords['fSpYieldMode'] 
            print '\t set yield to zero!'
            spYieldW=0.0;
            spYieldHe=0.0;
        
        #WEIGHTED SUM OF SPUTTERING YIELDS!
        totalSpYield=spYieldHe+fluxFractionW*spYieldW
        print "\t the effective (weighted by relative flux) sputtering yield in Xolotl is: ", totalSpYield

        if keywords['dStartMode']=='RESTART':
            restartNetworkFile = networkFile
            filepath='../../restart_files/'+restartNetworkFile
            shutil.copyfile(filepath,restartNetworkFile)

        if driverMode == 'INIT':
            #print ('run xolotl preprocessor')
            #print 'with java', os.system('echo $JAVA_EXE')
            #run prepocessor and copy params.txt input file to plasma state
            #since java is handled a little differently accross machines, 
            #a JAVA-XOLOTL environment variables are defined in the machine environment file
            #os.system('$JAVA_XOLOTL_EXE -Djava.library.path=$JAVA_XOLOTL_LIBRARY -cp .:$JAVA_XOLOTL_LIBRARY/*::$XOLOTL_PREPROCESSOR_DIR gov.ornl.xolotl.preprocessor.Main --perfHandler dummy --nxGrid 160 --maxVSize 250 --phaseCut')
            print 'init mode: run parameter file without preprocessor'
            write_xolotl_paramfile.writeXolotlParameterFile_fromTemplate(start_stop=startStop,ts_final_time=runEndTime,sputtering=totalSpYield,flux=flux,initialV=keywords['xInitialV'],nxGrid=keywords['xNGrid'])
                
        else:
            print 'restart mode: run parameter file without preprocessor'
            write_xolotl_paramfile.writeXolotlParameterFile_fromTemplate(start_stop=startStop,ts_final_time=runEndTime,useNetFile=True,networkFile=networkFile,sputtering=totalSpYield,flux=flux, initialV=keywords['xInitialV'],nxGrid=keywords['xNGrid'])
        
        #store xolotls parameter and network files for each loop 
        currentXolotlParamFile='params_%f.txt' %driverTime
        shutil.copyfile('params.txt',currentXolotlParamFile) 
        
        self.services.update_plasma_state()

    def step(self, timeStamp=0.0,**keywords):
        print('xolotl_worker: step')

        self.services.stage_plasma_state()

        #asign a local variable to arguments used multiple times
        driverTime=keywords['dTime']

        print 'check that all arguments are read well by xolotl-step'
        for (k, v) in keywords.iteritems():
            print '\t', k, " = ", v

        #call shell script that runs Xolotl and pipes input file
        task_id = self.services.launch_task(self.NPROC,
                                            self.services.get_working_dir(),
                                            self.XOLOTL_EXE, 'params.txt')
        #monitor task until complete
        if (self.services.wait_task(task_id)):
            self.services.error('xolotl_worker: step failed.')

        newest = max(glob.iglob('TRIDYN_*.dat'), key=os.path.getctime)
        print('newest file ' , newest)
        shutil.copyfile(newest, 'last_TRIDYN_toBin.dat')
        
        #re-bin last_TRIDYN file
        binTRIDYN.binTridyn()

        #store xolotls profile output for each loop (not plasma state)
        currentXolotlOutputFileToBin='last_TRIDYN_toBin_%f.dat' %driverTime
        shutil.copyfile('last_TRIDYN_toBin.dat', currentXolotlOutputFileToBin)
        currentXolotlOutputFile='last_TRIDYN_%f.dat' %driverTime
        shutil.copyfile('last_TRIDYN.dat', currentXolotlOutputFile)


        #append output
        tempfile = open(self.OUTPUT_XOLOTL_TEMP,"r")
        f = open(self.OUTPUT_XOLOTL_FINAL, "a")
        f.write(tempfile.read())
        f.close()
        tempfile.close()

        #save network file with a different name to use in the next time step
        currentXolotlNetworkFile='xolotlStop_%f.h5' %driverTime
        shutil.copyfile('xolotlStop.h5',currentXolotlNetworkFile)

        #updates plasma state Xolotl output files
        self.services.update_plasma_state()
  
    def finalize(self, timeStamp=0.0):
        return
    
