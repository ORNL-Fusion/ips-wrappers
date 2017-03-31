#! /usr/bin/env python

from  component import Component
import os
import shutil
import sys
import translate_xolotl_to_ftridyn
import generateInputIPS

class ftridynWorker(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)
        #get name of FTridyn executable from config file
        #This is actually a shell script that calls FTridyn and pipe the input to the executable
        self.ftridyn_exe = self.FTRIDYN_EXE

    def init(self, timeStamp=0.0):
        print('fridyn_worker: init')
        #stage plasma state files for use on execution of FTridyn
        self.services.stage_plasma_state()
        sys.path.append(os.getcwd())
        import driverParameterConfig
        reload(driverParameterConfig)

        if (driverParameterConfig.mode == 'INIT'):
            print('init mode yes')
            generateInputIPS.main()
        else:
            print('init mode no')
            nDataPts = translate_xolotl_to_ftridyn.xolotlToLay()
            generateInputIPS.main(IQ0=-1,NQX=nDataPts)

        #get name of FTridyn input file from config file to copy newly generated files to           
        current_ftridyn_namelist = self.services.get_config_param('FTRIDYN_INPUT_FILE')
        #this may be more than one file, not sure yet - need to learn more about FTridyn I/O
        from_file_list = self.COPY_FILES.split()
        file_list = current_ftridyn_namelist.split()

        #copy newly generated files to names specified in config file
        for index in range(0,1): #range(len(file_list)):
            print('copying ', from_file_list[index], ' to ', file_list[index])
            shutil.copyfile(from_file_list[index], file_list[index])

        self.services.update_plasma_state()

    def step(self, timeStamp=0.0):
        print('ftridyn_worker: step')
        self.services.stage_plasma_state()
        #call shell script that runs FTridyn and pipes input file
        task_id = self.services.launch_task(self.NPROC,
                                            self.services.get_working_dir(),
                                            self.FTRIDYN_EXE)
        #monitor task until complete
        if (self.services.wait_task(task_id)):
            self.services.error('ftridyn_worker: step failed.')

        os.system(' '.join(['python', self.POSTPROCESSING_SCRIPT]))

        #append output
        tempfile = open(self.OUTPUT_FTRIDYN_TEMP,"r")
        f = open(self.OUTPUT_FTRIDYN_FINAL, "a")
        f.write(tempfile.read())
        f.close()
        tempfile.close()

        #updates plasma state FTridyn output files
        self.services.update_plasma_state()
  
    def finalize(self, timeStamp=0.0):
        return
    
