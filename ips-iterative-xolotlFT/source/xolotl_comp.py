#! /usr/bin/env python

from  component import Component
import os
import shutil
import subprocess
import glob
import translate_xolotl_to_ftridyn
import binTRIDYN
import write_xolotl_paramfile
import sys

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

        print 'xolotl-init:'
        print '\t \t driver mode is', driverMode
        print '\t \t running starts at time',driverTime
        print '\t \t \t  ends at time', runEndTime
        print '\t \t driver step is', keywords['dTimeStep']


        #if fSpYield in driver < 0, fSpYieldMode = calculate -> get sputtering yield FROM He_WOUT.DAT
        if keywords['fSpYieldMode']=='calculate':
            ftridynOutFile=open('He_WOUT.DAT',"r")
            ftridynOutData=ftridynOutFile.read().split('\n')
            searchString='PARTICLES(2)'
            for line in ftridynOutData:
                if searchString in line:
                    break
            stringWithEmptyFields=line.strip().split(" ")
            sputteringNparticlesString=[x for x in stringWithEmptyFields if x]
            sputteringNparticles=sputteringNparticlesString[2]
            spYieldW=float(sputteringNparticles)/float(keywords['fNImpacts'])
            print 'calculated in Xolotl-component: W sputtering yield is =', spYieldW
        #if fSpYield in driver >= 0, fSpYieldMode = fixed 
        elif keywords['fSpYieldMode']=='fixed':
            spYieldW=keywords['fSpYieldW']
            print 'Fixed value of W sputtering yield is =', spYieldW
        else:
            print 'Invalid value of fSpYieldMode, ', keywords['fSpYieldMode'] 
            print 'set yield to zero!'
            spYieldW=0.0;

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
            print 'run parameter file without preprocessor'
            write_xolotl_paramfile.writeXolotlParameterFile_fromTemplate(start_stop=startStop,ts_final_time=runEndTime,sputtering=spYieldW,flux=flux)
                
        else:
            write_xolotl_paramfile.writeXolotlParameterFile_fromTemplate(start_stop=startStop,ts_final_time=runEndTime,useNetFile=True,networkFile=networkFile,sputtering=spYieldW,flux=flux)
        
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
    
