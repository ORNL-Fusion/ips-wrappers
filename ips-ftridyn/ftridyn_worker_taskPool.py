#! /usr/bin/env python

from  component import Component
import os
import shutil
import re
import sys
import time
import generate_ftridyn_input
import analyze_ftridyn_simulations

class ftridynWorker(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)
        #get name of FTridyn executable from config file
        #This is actually a shell script that calls FTridyn and pipes
        #the input to the executable
        self.ftridyn_exe = self.FTRIDYN_EXE

    def init(self, timeStamp=0.0,**keywords):
        #stage plasma state files for use on execution of FTridyn
        self.services.stage_plasma_state()
        #unpack arguments from init call in driver
        ffilename = keywords["ffilename"]
        beam = keywords["beam"]
        target = keywords["target"]
        nH = keywords["nH"]
        energy = keywords["eArg"]
        angle = keywords["aArg"]
        roughness = keywords["dArg"]

        for i in range(len(energy)):
            for j in range(len(angle)):
                for k in range(len(roughness)):
                    generate_ftridyn_input.beam_and_target(ffilename,beam,target,number_histories=nH,incident_energy=energy[i],incident_angle=angle[j],fractal_dimension=roughness[k])
                    pathString = "FTRIDYN_"+str(energy[i]) + "_"+str(angle[j])+"_" + str(roughness[k])
                    if not os.path.exists(pathString):
                        os.makedirs(pathString)
                    shutil.copyfile(ffilename+"0001.IN", pathString+"/"+ffilename+"0001.IN")
                    
                    surfname = str(1)+'p'+str(int((roughness[k]-0.99999)*1000.0))+'.surf'
                    shutil.copyfile(surfname, pathString+"/"+surfname)
                    if energy[i] < 0:
                        #copy energy distribution file
                        shutil.copyfile(self.INPUT_DIR+"/dist"+str(int(j))+".dat", pathString+"/"+ffilename+"0001.ED1") 
 
    def step(self, timeStamp=0.0,**keywords):
        print('ftridyn_worker: step (task pool version)')
        
        ffilename = keywords["ffilename"]
        beam = keywords["beam"]
        target = keywords["target"]
        nH = keywords["nH"] #number_histories=nH
        energy = keywords["eArg"]
        angle = keywords["aArg"]
        roughness = keywords["dArg"]
       
        cwd = self.services.get_working_dir()
        pool = self.services.create_task_pool('pool')
        for i in range(len(energy)):
            for j in range(len(angle)):
                for k in range(len(roughness)):
                    pathString = "FTRIDYN_"+str(energy[i]) + "_"+str(angle[j])+"_" + str(roughness[k])
                    self.services.add_task('pool', 'task'+str(energy[i]) + "_"+str(angle[j])+"_" + str(roughness[k]), 1, cwd+"/"+pathString, self.FTRIDYN_EXE, ffilename+"0001.IN",logfile='task_'+pathString+'.log' )
        ret_val = self.services.submit_tasks('pool')
        print 'ret_val = ', ret_val
        exit_status = self.services.get_finished_tasks('pool')
        print exit_status
        self.services.remove_task_pool('pool')


        spyl_file = ffilename+'SPYL.DAT'
        driver_out = self.services.get_config_param('EA_OUTPUT')
        fid = open(driver_out,'a')
        for i in range(len(energy)):
            for j in range(len(angle)):
                for k in range(len(roughness)):
                    pathString = "FTRIDYN_"+str(energy[i]) + "_"+str(angle[j])+"_" + str(roughness[k])
                    print('doing analysis on ' , pathString+"/"+spyl_file)
                    WW=analyze_ftridyn_simulations.ftridyn_output_data(pathString+"/"+ffilename,'0001')
                    thisSpyl = WW.calculate_total_sputtering_yield()
                    print('sputtering yield', thisSpyl)
                    nP = analyze_ftridyn_simulations.plot_sputtered_angular_distributions(pathString+"/"+ffilename)
                    analyze_ftridyn_simulations.plot_implantation_profile(pathString+"/"+ffilename)

                    fid.write(" ".join([str(energy[i]),str(angle[j]),str(roughness[k]),'  ',str(thisSpyl),'\n']))
        fid.close()
        #updates plasma state FTridyn output files
        self.services.update_plasma_state()
  
    def finalize(self, timeStamp=0.0):
        return
    
