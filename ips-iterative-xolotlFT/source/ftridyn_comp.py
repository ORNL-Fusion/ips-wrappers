#! /usr/bin/env python

from  component import Component
import os
import shutil
import glob
import sys
import translate_xolotl_to_ftridyn
import generateInputIPS
import numpy as np
import subprocess
import re

class ftridynWorker(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)
        #get name of FTridyn executable from config file
        #This is actually a shell script that calls FTridyn and pipe the input to the executable
        self.ftridyn_exe = self.FTRIDYN_EXE

    def init(self, timeStamp=0.0,**keywords):
        print('fridyn_worker: init')
        #stage plasma state files for use on execution of FTridyn
        self.services.stage_plasma_state()

#        sys.path.append(os.getcwd())
#        import driverParameterConfig
#        reload(driverParameterConfig)
#        import ftridynParameterConfig
#        reload(ftridynParameterConfig)

        print 'check that all arguments are read well by ftridyn-init'
        for (k, v) in keywords.iteritems():
            print '\t', k, " = ", v

        #asign a local variable to arguments used multiple times
        initialTotalDepth=keywords['fInitialTotalDepth']
        totalDepth=keywords['fTotalDepth']
        driverTime=keywords['dTime']
        energyIn=keywords['fEnergyIn']
        nImpacts=keywords['fNImpacts']

        #prepare and save (for each loop) ftridyn input
        if (keywords['dMode'] == 'INIT'):
            print('init mode yes')
            if (totalDepth==0.0):
                nTT=initialTotalDepth
            else:
                nTT=totalDepth
            print 'calling generateInput with TT=%d, NH=%d and Ein=%f' % (nTT,nImpacts,energyIn)
            generateInputIPS.main(TT=nTT,NH=nImpacts,energy=energyIn)

            newestIn = max(glob.iglob('*.IN'), key=os.path.getctime)
            currentFtridynInFile='%s_%f' %(newestIn,driverTime)
            shutil.copyfile(newestIn,  currentFtridynInFile)

        else:
            print('init mode no')
            nDataPts = translate_xolotl_to_ftridyn.xolotlToLay(totalDepth=totalDepth)
            newestLay = max(glob.iglob('*.LAY'), key=os.path.getctime)
            currentFtridynLayFile='%s_%f' %(newestLay,driverTime)
            shutil.copyfile(newestLay, currentFtridynLayFile)
            
            if (totalDepth==0.0):
                nTT=10*np.max(np.loadtxt('last_TRIDYN.dat')[:,0]) 
            else:
                nTT=totalDepth

            print 'calling generateInput with NQX=%d , TT=%d , NH=%d and Ein=%f' % (nDataPts, nTT,nImpacts,energyIn)
            generateInputIPS.main(IQ0=-1,NQX=nDataPts,TT=nTT,NH=nImpacts,energy=energyIn)
            newestIn = max(glob.iglob('*.IN'), key=os.path.getctime)
            currentFtridynInFile='%s_%f' %(newestIn,driverTime)
            shutil.copyfile(newestIn,  currentFtridynInFile)


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

    def step(self, timeStamp=0.0,**keywords):
        print('ftridyn_worker: step')
        self.services.stage_plasma_state()

        print 'check that all arguments are read well by ftridyn-step'
        for (k, v) in keywords.iteritems():
            print '\t', k, " = ", v


        #call shell script that runs FTridyn and pipes input file
        task_id = self.services.launch_task(self.NPROC,
                                            self.services.get_working_dir(),
                                            self.FTRIDYN_EXE)
        #monitor task until complete
        if (self.services.wait_task(task_id)):
            self.services.error('ftridyn_worker: step failed.')

        #post-processing: re-format output (fit function to profile) for Xolotl
        os.system(' '.join(['python', self.POSTPROCESSING_SCRIPT]))

        #get W sputtering yield from FT output: He_WSPYL.DAT or He_WOUT.DAT
#        import ftridynParameterConfig
#        reload(ftridynParameterConfig)

        #FROM He_WSPYL.DAT
        #ftridynYieldOutFile = open('He_WSPYL.DAT', 'rb')
        #allSputteringData=[row.strip().split(" ") for row in ftridynYieldOutFile]
        #sputteringData=[[x for x in allSputteringData[0] if x],[x for x in allSputteringData[1] if x]]
        #sputteringYieldW=sputteringData[1][2]
        #print 'W sputtering yield from file ',ftridynYieldOutFile,' is =', sputteringYieldW

        #FROM He_WOUT.DAT and use NH
        ftridynOutFile=open('He_WOUT.DAT',"r")
        ftridynOutData=ftridynOutFile.read().split('\n')
        searchString='PARTICLES(2)'
        for line in ftridynOutData:
            if searchString in line:
                break
        stringWithEmptyFields=line.strip().split(" ")
        sputteringNparticlesString=[x for x in stringWithEmptyFields if x]
        sputteringNparticles=sputteringNparticlesString[2]
        keywords['fSpYieldW']=float(sputteringNparticles)/float(keywords['fNImpacts'])
        print 'calculated in ftridyn-component: W sputtering yield is =', keywords['fSpYieldW']

        #and replace the value of spYieldW in FTridyn Parameter Config File
#        ftridyn_config_file = self.services.get_config_param('FTRIDYN_PARAMETER_CONFIG_FILE')
#        ftridyn_tmp_file='ftridynConfigFile.tmp'
#        shutil.copyfile(ftridyn_config_file,ftridyn_tmp_file)
#        sputteringYieldSedString="sed    -e 's/spYieldW=[^ ]*/spYieldW=%s/' <%s >%s "   %(sputteringYieldW,ftridyn_tmp_file,ftridyn_config_file)
#        subprocess.call([sputteringYieldSedString], shell=True)
#        os.remove(ftridyn_tmp_file)

        #append output
        tempfile = open(self.OUTPUT_FTRIDYN_TEMP,"r")
        f = open(self.OUTPUT_FTRIDYN_FINAL, "a")
        f.write(tempfile.read())
        f.close()
        tempfile.close()

        #store ftridyns output for each loop (not plasma state)
#        import driverParameterConfig
        currentFtridynOutputFile='He_WDUMPPRJ_%f.dat' %keywords['dTime']
        shutil.copyfile('He_WDUMPPRJ.dat',currentFtridynOutputFile)

        #updates plasma state FTridyn output files
        self.services.update_plasma_state()
  
    def finalize(self, timeStamp=0.0):
        return
    
