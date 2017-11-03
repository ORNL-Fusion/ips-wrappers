#! /usr/bin/env python

from  component import Component
import os
import shutil
import subprocess
import glob
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
        self.driverTime=keywords['dTime']
        driverMode=keywords['dMode']

        flux=keywords['xFlux']
        totalSpYield=keywords['weightedSpYield']
        runEndTime=self.driverTime+keywords['dTimeStep']

        startStop=keywords['xStartStop']
        networkFile=keywords['xNetworkFile']
        paramTemplateFile=keywords['xParamTemplate']

        self.coupling=keywords['xFtCoupling']
        self.zipOutput=keywords['dZipOutput']
        self.petscHeConc=keywords['xHe_conc']
        self.processes=keywords['xProcess']

        self.xNxGrid=keywords['xNxGrid']
        self.xNyGrid=str(keywords['xNyGrid']) 
        self.xDxGrid=keywords['xDxGrid']
        self.xDyGrid=str(keywords['xDyGrid'])

        cwd = self.services.get_working_dir()

        print 'xolotl-init:'
        print '\t driver mode is', driverMode
        print '\t running starts at time',self.driverTime
        print '\t \t  ends at time', runEndTime
        print '\t driver step is', keywords['dTimeStep']
        print '\n'

        
        if keywords['dStartMode']=='RESTART':
            restartNetworkFile = networkFile
            filepath='../../restart_files/'+restartNetworkFile
            shutil.copyfile(filepath,restartNetworkFile)

        if driverMode == 'INIT':
            print 'init mode: run parameter file without preprocessor'
            write_xolotl_paramfile.writeXolotlParameterFile_fromTemplate(dimensions=keywords['xDimensions'], infile=paramTemplateFile, fieldsplit_1_pc_type=keywords['xFieldsplit_1_pc_type'],start_stop=startStop,ts_final_time=runEndTime,sputtering=totalSpYield,flux=flux,initialV=keywords['xInitialV'],nxGrid=self.xNxGrid,nyGrid=self.xNyGrid,dxGrid=self.xDxGrid,dyGrid=self.xDyGrid, he_conc=self.petscHeConc, process=self.processes, voidPortion=keywords['xVoidPortion'])
                
        else:
            print 'restart mode: run parameter file without preprocessor'
            write_xolotl_paramfile.writeXolotlParameterFile_fromTemplate(dimensions=keywords['xDimensions'], infile=paramTemplateFile, fieldsplit_1_pc_type=keywords['xFieldsplit_1_pc_type'],start_stop=startStop,ts_final_time=runEndTime,useNetFile=True,networkFile=networkFile,sputtering=totalSpYield,flux=flux, initialV=keywords['xInitialV'],nxGrid=self.xNxGrid,nyGrid=self.xNyGrid,dxGrid=self.xDxGrid,dyGrid=self.xDyGrid,he_conc=self.petscHeConc, process=self.processes, voidPortion=keywords['xVoidPortion'])
        
        #store xolotls parameter and network files for each loop 
        currentXolotlParamFile='params_%f.txt' %self.driverTime
        shutil.copyfile('params.txt',currentXolotlParamFile) 
        
        self.services.update_plasma_state()

    def step(self, timeStamp=0.0,**keywords):
        print('xolotl_worker: step')

        self.services.stage_plasma_state()

        #asign a local variable to arguments used multiple times

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
        shutil.copyfile(newest, 'last_TRIDYN.dat')


        #save TRIDYN_*.dat files, zipped
        TRIDYNFiles='TRIDYN_*.dat'

        if self.coupling and self.zipOutput:#=='True'):              
            TRIDYNZipped='allTRIDYN_t%f.zip' %self.driverTime
            zip_ouput='zipTridynDatOuput.txt'

            print 'save and zip output: ', TRIDYNFiles
            zipString='zip %s %s >> %s ' %(TRIDYNZipped, TRIDYNFiles, zip_ouput)
            subprocess.call([zipString], shell=True)

            rmString='rm '+ TRIDYNFiles
            subprocess.call([rmString], shell=True)

        elif self.coupling:
            print 'leaving ',  TRIDYNFiles , 'uncompressed'

        else:
             print 'no ', TRIDYNFiles , 'used in this simulation'


        #save helium concentration files, zipped
        heConcFiles='heliumConc_*.dat'

        if self.petscHeConc and self.zipOutput:#=='True'):
            heConcZipped='allHeliumConc_t%f.zip' %self.driverTime
            zip_ouput='zipHeConcOuput.txt'
            
            print 'save and zip output: ', heConcFiles
            zipString='zip %s %s >> %s ' %(heConcZipped, heConcFiles, zip_ouput)
            subprocess.call([zipString], shell=True)

            rmString='rm '+heConcFiles
            subprocess.call([rmString], shell=True)
            
        elif self.petscHeConc :
            print 'leaving ', heConcFiles ,'uncompressed'

        else:
            print 'no ', heConcFiles , ' in this loops output'

        #save network file with a different name to use in the next time step
        currentXolotlNetworkFile='xolotlStop_%f.h5' %self.driverTime
        shutil.copyfile('xolotlStop.h5',currentXolotlNetworkFile)
            
        #updates plasma state Xolotl output files
        self.services.update_plasma_state()
  
    def finalize(self, timeStamp=0.0):
        return
    
