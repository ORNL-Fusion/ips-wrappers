#! /usr/bin/env python

from  component import Component
import os
import shutil
import re
import sys

class ftridynWorker(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)
        #get name of FTridyn executable from config file
        #This is actually a shell script that calls FTridyn and pipes
        #the input to the executable
        self.ftridyn_exe = self.FTRIDYN_EXE

    def init(self, timeStamp=0.0,**keywords):
        sys.path.append(os.path.expandvars('$FTRIDYN_PATH/utils'))
        import generateInputTY
        #stage plasma state files for use on execution of FTridyn
        self.services.stage_plasma_state()

        energy = keywords["eArg"]
        angle = keywords["aArg"]
        roughness = keywords["dArg"]

        ##get name of FTridyn input file from config file to copy newly generated files to
        #this may be more than one file, not sure yet - need to learn more about FTridyn I/O
        from_file_list = self.COPY_FILES.split()
        #file_list = current_ftridyn_namelist.split()
 
        for i in range(len(energy)):
            for j in range(len(angle)):
                for k in range(len(roughness)):
                    #os.system(' '.join(['python', self.INPUT_SCRIPT,'-R 1 -s 0 -E',str(energy[i]),'-a',str(angle[j]),'-d',str(roughness[k])]))
                    generateInputTY.main(1, 0,energy[i],angle[j], roughness[k])
                    pathString = "FTRIDYN_"+str(energy[i]) + "_"+str(angle[j])+"_" + str(roughness[k])
                    fileString = pathString+".IN"
                    if not os.path.exists(pathString):
                        os.makedirs(pathString)
                    #copy newly generated files to names specified in config file
                    for index in range(0,1): #range(len(file_list)): this may need to be changed
                    #print('copying ', from_file_list[index], ' to ', file_list[index])
                        shutil.copyfile(from_file_list[index], pathString+"/"+fileString)
                    
                    shutil.copyfile('surface.surf', pathString+"/"+'surface.surf')
                    if energy[i] < 0:
                        #copy energy distribution file
                        shutil.copyfile(self.INPUT_DIR+"/dist"+str(int(j))+".dat", pathString+"/"+"He_W0001.ED1")#pathString+".ED1") 
    def step(self, timeStamp=0.0,**keywords):
        print('ftridyn_worker: step (task pool version)')
        
        energy = keywords["eArg"]
        angle = keywords["aArg"]
        roughness = keywords["dArg"]
        ##call shell script that runs FTridyn and pipes input file
        #task_id = self.services.launch_task(self.NPROC,
        #                                    self.services.get_working_dir(),
        #                                    self.FTRIDYN_EXE)
        ##monitor task until complete
        #if (self.services.wait_task(task_id)):
        #    self.services.error('ftridyn_worker: step failed.')
       
        cwd = self.services.get_working_dir()
        pool = self.services.create_task_pool('pool')
        for i in range(len(energy)):
            for j in range(len(angle)):
                for k in range(len(roughness)):
                    pathString = "FTRIDYN_"+str(energy[i]) + "_"+str(angle[j])+"_" + str(roughness[k])
                    fileString = pathString+".IN"
                    self.services.add_task('pool', 'task'+str(energy[i]) + "_"+str(angle[j])+"_" + str(roughness[k]), 1, cwd+"/"+pathString, self.FTRIDYN_EXE, fileString,logfile='task_'+pathString+'.log' )
        
        ret_val = self.services.submit_tasks('pool')
        print 'ret_val = ', ret_val
        exit_status = self.services.get_finished_tasks('pool')
        print exit_status
        self.services.remove_task_pool('pool')


        spyl_file = 'He_WSPYL.DAT'
        driver_out = self.services.get_config_param('EA_OUTPUT')
        fid = open(driver_out,'a')
        for i in range(len(energy)):
            for j in range(len(angle)):
                for k in range(len(roughness)):
                    pathString = "FTRIDYN_"+str(energy[i]) + "_"+str(angle[j])+"_" + str(roughness[k])
                    fileString = pathString+".IN"
                    print('doing analysis on ' , pathString+"/"+spyl_file)
                    fid0 = open(pathString+"/"+spyl_file,'r')
                    lines = fid0.readlines()
                    fid0.close()
                    if len(lines) > 1:
                        print('lines', str(lines[1]))
                        list1 = re.sub(' +',' ',lines[1]).split()
                        print('list1', list1)
                        print(list1[0],list1[1])

                        fid.write(" ".join([str(energy[i]),str(angle[j]),str(roughness[k]),'  ',list1[2],'\n']))
        fid.close()
        #updates plasma state FTridyn output files
        self.services.update_plasma_state()
  
    def finalize(self, timeStamp=0.0):
        return
    
