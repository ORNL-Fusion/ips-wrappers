#! /usr/bin/env python

from  component import Component
import os
import shutil
import subprocess
import glob
import translate_xolotl_to_ftridyn
import write_xolotl_paramfile
import sys

class xolotlWorker(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)
        

    def init(self, timeStamp=0.0):
        print('xolotl_worker: init')
        self.services.stage_plasma_state()

        sys.path.append(os.getcwd())
        import driverParameterConfig
        reload(driverParameterConfig)

        print 'from drivers parameterConfig file: mode is ', driverParameterConfig.mode
        print '\t \t \t \t this run starts at time ', driverParameterConfig.driverTime
        print '\t \t \t \t driverTimeStep is ', driverParameterConfig.driverTimeStep
        runEndTime=driverParameterConfig.driverTime+driverParameterConfig.driverTimeStep
        print '\t \t \t \t this runs ends at time: ', runEndTime

        import xolotlParameterConfig
        reload(xolotlParameterConfig)

        #print here what's been loaded from xolotlParameterConfig File
        xolotl_config_file = self.services.get_config_param('XOLOTL_PARAMETER_CONFIG_FILE')
        xid = open(xolotl_config_file, 'r')
        #currently only start_stop might be passed to writing the parameter file.
        #if more parameters, it'd require writing them all as a string. something similar to:
        xolotlParamString=""
        for line in xid: 
            print line  
            xolotlParamString=xolotlParamString+line.replace("\n", " , ")

        print 'from xolotls parameterConfig file:', xolotlParamString

        if (driverParameterConfig.mode == 'INIT'):
            print('run xolotl preprocessor')
            #run prepocessor and copy params.txt input file to plasma state
            os.system('java -Djava.library.path=/project/projectdirs/atom/atom-install-edison/xolotl/xolotl-trunk-source/gov.ornl.xolotl.preprocessor/deps -cp .:/project/projectdirs/atom/atom-install-edison/xolotl/xolotl-trunk-source/gov.ornl.xolotl.preprocessor/deps/*:/project/projectdirs/atom/atom-install-edison/xolotl/xolotl-trunk-build/gov.ornl.xolotl.preprocessor/preprocessor/CMakeFiles/xolotlPreprocessor.dir/ gov.ornl.xolotl.preprocessor.Main --nxGrid 160 --maxVSize 250 --phaseCut')        
            write_xolotl_paramfile.writeXolotlParameterFile_fromPreprocessor(start_stop=xolotlParameterConfig.start_stop,ts_final_time=runEndTime,networkFile=xolotlParameterConfig.networkFile)
                
        else:
            write_xolotl_paramfile.writeXolotlParameterFile_fromTemplate(start_stop=xolotlParameterConfig.start_stop,ts_final_time=runEndTime,networkFile=xolotlParameterConfig.networkFile)
        
        #store xolotls parameter and network files for each loop 
        currentXolotlParamFile='params_%f.txt' %driverParameterConfig.driverTime
        shutil.copyfile('params.txt',currentXolotlParamFile) 
        
        #currentXolotlNetworkFile='xolotlStop_%f.h5' %driverParameterConfig.driverTime
        #shutil.copyfile('xolotlStop.h5',currentXolotlNetworkFile)

        self.services.update_plasma_state()

    def step(self, timeStamp=0.0):
        print('xolotl_worker: step')

        sys.path.append(os.getcwd())
        import driverParameterConfig
        reload(driverParameterConfig)

        self.services.stage_plasma_state()
        #call shell script that runs Xolotl and pipes input file
        task_id = self.services.launch_task(self.NPROC,
                                            self.services.get_working_dir(),
                                            self.XOLOTL_EXE, 'params.txt')
        #monitor task until complete
        if (self.services.wait_task(task_id)):
            self.services.error('xolotl_worker: step failed.')

        newest = max(glob.iglob('TRIDYN_*.dat'), key=os.path.getctime)
        print('newest file ' , newest)
        shutil.copyfile(newest, 'TRIDYN_last.dat')

        #append output
        tempfile = open(self.OUTPUT_XOLOTL_TEMP,"r")
        f = open(self.OUTPUT_XOLOTL_FINAL, "a")
        f.write(tempfile.read())
        f.close()
        tempfile.close()

        #save network file with a different name to use in the next time step
        currentXolotlNetworkFile='xolotlStop_%f.h5' %driverParameterConfig.driverTime
        shutil.copyfile('xolotlStop.h5',currentXolotlNetworkFile)

        #updates plasma state Xolotl output files
        self.services.update_plasma_state()
  
    def finalize(self, timeStamp=0.0):
        return
    
