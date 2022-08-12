#! /usr/bin/env python

from ipsframework import Component
import os
import shutil
import subprocess
import glob
import translate_xolotl_to_lay
import write_xolotl_paramfile
import sys

class xolotlWorker(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)
        #get name of FTridyn executable from config file
        #This is actually a shell script that calls FTridyn and pipes
        #the input to the executable
        #self.ftridyn_exe = self.FTRIDYN_EXE

    def init(self, timeStamp=0.0):
        self.services.stage_state()

        sys.path.append(os.getcwd())
        import parameterConfig
        reload(parameterConfig)
        #parameter_config_file = self.services.get_config_param('PARAMETER_CONFIG_FILE')
        print('from parameterConfig file ', parameterConfig.mode)

        if (parameterConfig.mode == 'INIT'):
            print('run xolotl preprocessor')
            #run prepocessor and copy params.txt input file to plasma state
            os.system('java -Djava.library.path=/project/projectdirs/atom/atom-install-edison/xolotl/xolotl-trunk-source/gov.ornl.xolotl.preprocessor/deps -cp .:/project/projectdirs/atom/atom-install-edison/xolotl/xolotl-trunk-source/gov.ornl.xolotl.preprocessor/deps/*:/project/projectdirs/atom/atom-install-edison/xolotl/xolotl-trunk-build/gov.ornl.xolotl.preprocessor/preprocessor/CMakeFiles/xolotlPreprocessor.dir/ gov.ornl.xolotl.preprocessor.Main --nxGrid 160 --maxVSize 250 --phaseCut')
        
            #self.services.stage_input_files(self.INPUT_FILES)
            #shutil.copyfile('params0.txt','params.txt')
            write_xolotl_paramfile.writeXolotlParameterFile_fromPreprocessor()
            shutil.copyfile('params.txt','paramsInit.txt')  #store file to look into network issue; line to be removed
        else:
            write_xolotl_paramfile.writeXolotlParameterFile_fromTemplate(start_stop=0.2,ts_final_time=0.2,networkFile="xolotlStop.h5",sputtering=0.1)
            shutil.copyfile('params.txt','paramsRestart.txt') #store file to look into network issue; line to be removed
        #self.services.stage_state()
        self.services.update_state()

    def step(self, timeStamp=0.0):
        print('xolotl_worker: step')
        self.services.stage_state()
        #call shell script that runs FTridyn and pipes input file
        task_id = self.services.launch_task(self.NPROC,
                                            self.services.get_working_dir(),
                                            self.XOLOTL_EXE, 'params.txt')
        #monitor task until complete
        if (self.services.wait_task(task_id)):
            self.services.error('ftridyn_worker: step failed.')
        #get number of output files from xolotl
        #nTS = subprocess.Popen(["ls TRIDYN_*.dat | awk 'END(print NR-1)'"], stdout=subprocess.PIPE) 
        #print('nTS printed after subprocess call', nTS.communicate())

        newest = max(glob.iglob('TRIDYN_*.dat'), key=os.path.getctime)
        print('newest file ' , newest)
        shutil.copyfile(newest, 'TRIDYN_last.dat')

        #append output
        tempfile = open(self.OUTPUT_XOLOTL_TEMP,"r")
        f = open(self.OUTPUT_XOLOTL_FINAL, "a")
        f.write(tempfile.read())
        f.close()
        tempfile.close()


        #translate_xolotl_to_lay.xolotlToLay()

        #updates plasma state FTridyn output files
        self.services.update_state()
  
    def finalize(self, timeStamp=0.0):
        return
    
