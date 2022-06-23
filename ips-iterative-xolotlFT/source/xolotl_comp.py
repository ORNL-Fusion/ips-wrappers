#! /usr/bin/env python

from  component import Component
import os
import shutil
import subprocess
import glob
import param_handler #write_xolotl_paramfile
import keepLastTS
import sys
import numpy as np

class xolotlWorker(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))
        

    def init(self, timeStamp=0.0, **keywords):

        self.services.stage_plasma_state()
        cwd = self.services.get_working_dir()

        print('\n')

        #asign a local variable to arguments used multiple times 
        self.driverTime=keywords['dTime']
        #self.coupling=keywords['xFtCoupling']

        cwd = self.services.get_working_dir()

        xp = param_handler.xolotl_params()
        xp.parameters=keywords['xParameters'] 
        print_test=keywords['print_test']
        
        if 'output_file' in keywords:
            outFile=keywords['output_file']
            if outFile  is not None:
                print('\t redirect Xolotl:init output')
                print('\t \t of ', cwd )
                print('\t \t to:', outFile)
                outF=open(outFile , 'a')
                sys.stdout = outF
            else:
                print ('\t no log file defined in keywords or config file')
                print ('\t print output of Xolotl:init to default sys.stdout')

        print (' ')
        print('xolotl_worker: init')
        print(' ')
        # if TEST:
        if print_test:
            print('\t TEST: check that all arguments are read well by xolotl-init and write Xolotl input file (from dictionary)')
            for (k, v) in keywords.items():
                print(('\t \t {0} = {1}'.format(k, v)))
            print('\t ... done checking that all arguments are read well by xolotl-init ')
            print(' ')
        
        # if migration parameters are specified, copy them to a file 'migration.txt'
        if 'migration' in xp.parameters.keys():
            with open('migration.txt', 'w') as f:
                for param in xp.parameters['migration']:
                    f.write(f"{param}\n")
            xp.parameters.pop('migration')

        #write and store xolotls parameter for each loop 
        xp.write('params.txt')

        currentXolotlParamFile='params_%f.txt' %self.driverTime
        shutil.copyfile('params.txt',currentXolotlParamFile)


        try:
            #shutil.copyfile(xp.parameters['networkFile'],xp.parameters['networkFile']+'_t'+str(self.driverTime))
            shutil.copyfile(self.NETWORK_FILE,self.NETWORK_FILE+'_t'+str(self.driverTime))
            if print_test:
                print('\t TEST: copied ',self.NETWORK_FILE,' as ' , self.NETWORK_FILE+'_t'+str(self.driverTime))
                print('\t file size of ', self.NETWORK_FILE ,' is: ', os.path.getsize(self.NETWORK_FILE))

        except Exception as e:
            print('\t', e)
            print('\t could not save network file for t = ', str(self.driverTime))            
            
        sys.stdout.flush()
        self.services.update_plasma_state()

    def step(self, timeStamp=0.0,**keywords):

        self.services.stage_plasma_state()
        cwd = self.services.get_working_dir()
        zipOutput=keywords['dZipOutput']
        petscHeConc=keywords['xHe_conc']
        xp_parameters=keywords['xParameters']
        outFile=keywords['output_file']
        n_overgrid_loops=keywords['n_overgrid_loops']
        print_test=keywords['print_test']
        num_tries = keywords['xolotl_num_tries']

        print(' ')
        if 'output_file' in keywords:
            outFile=keywords['output_file']            
            if outFile  is not None:
                print('\t redirect Xolotl:step output')
                print('\t \t of ', cwd )
                print('\t \t to:', outFile)                
                outF=open(outFile , 'a')
                sys.stdout = outF
            else:
                print ('\t no log file defined in keywords or config file')
                print ('\t print output of Xolotl:step to default sys.stdout')

        print(' ')
        print('xolotl_worker: step ')
        print(' ')


        #asign a local variable to arguments used multiple times      
        # if TEST mode:
        if print_test:
            print('\t TEST: checking that all arguments are read well by xolotl-step')
            for (k, v) in keywords.items():
                print('\t \t', k, " = ", v) 
            print('\t ... done checking that all arguments are read well by xolotl-step ')
            print(' ')

        #xolotlLogFile='xolotl_t%f.log' %self.driverTime
        #print '\t Xolotl log file ', xolotlLogFile

        #call shell script that runs Xolotl and pipes input file
        #task_id = self.services.launch_task(self.NPROC,
        #                                    self.services.get_working_dir(),
        #                                    self.XOLOTL_EXE, 'params.txt', 
        #                                    logfile=xolotlLogFile)

        #monitor task until complete
        #if (self.services.wait_task(task_id)):
        #    self.services.error('xolotl_worker: step failed.')


        import time
        os.environ['OMP_NUM_THREADS']=self.THREADS_PER_TASK
        for i in range(num_tries):
            xolotlLogFile='xolotl_t%f_%d.log' %(self.driverTime,i)
            task_id = self.services.launch_task(self.NPROC,self.services.get_working_dir(),
                                                self.XOLOTL_EXE, 'params.txt',task_ppn=self.task_ppn,logfile=xolotlLogFile)
            ret_val = self.services.wait_task(task_id)
            if print_test:
                print('TEST: Xolotl run done, exited with ret_val = ', ret_val)

            statusFile=open(self.EXIT_STATUS, "r")
            exitStatus=statusFile.read().rstrip('\n')
            statusFile.close()
            print("statusFile: ", exitStatus)
            sys.stdout.flush()
            self.services.update_plasma_state()

            #if out of grid space, do not keep on trying
            if exitStatus=='overgrid':
                #print("aborting run, out of void space in grid: fix network file in next restart")
                print("aborting run, out of void space in grid: add grid points and try again")
                self.services.update_plasma_state()
                break

            if (ret_val == 0):
                break
            else:
                self.services.error('xolotl_worker: step failed in trial %d.' %i)
                sys.stdout.flush()
                time.sleep(5)
                #if it failed, save last networkFile before a new try and set newest network file to use in the
                # next try, so that it starts from the last saved time step, not from the beginning of the loop
                shutil.copyfile(self.NETWORK_FILE,'networkFile_%f_%d.h5' %(self.driverTime,i))
                ## use keepLastTS to produce netfile with only info from the last TS
                print('\t produce new network file using xolotlStop:')
                try:
                    iF=cwd+'/xolotlStop.h5'
                    oF= cwd+'/'+self.NETWORK_FILE
                    os.remove(oF) #can not exist & it's copied as w/ time-stamp above
                    print('\t \t run keepLastTS with:')
                    print('\t \t \t inFile = ', iF)
                    print('\t \t \t  outFile = ', oF)
                    keepLastTS.keepLastTS(inFile=iF, outFile=oF, print_test=print_test)
                    print('\t \t ... keepLastTS done')
                #if fails, use old method of copying entire xolotlStop as networkFile
                except Exception as e:
                    print(e)
                    print('\t \t running keepLastTS failed')
                    print('\t \t just copy xolotlStop as networkFile')
                    shutil.copyfile('xolotlStop.h5',self.NETWORK_FILE)
                print('\t done writing a new network file')
                sys.stdout.flush()

        else:
            self.services.error('xolotl_worker: Aborting after %d num_trials trials' %num_trials)
            raise Exception("Aborting simulation after %d failed xolotl runs" % num_trials)

        #UPDATE: ALREADY IN DRIVER save network file with a different name to use in the next time step
        
        ##UPDATE: TRIDYN file in HDF5 format now, 
        #generated by Xolotl when saving concentrations (start-stop)
        #-> pass HDF5 file to binTRIDYN, will output text file
        #-> no need to delete or compress files
        newest = max(glob.iglob('TRIDYN_*.h5'), key=os.path.getctime)  
        print(('\t newest file {} \n'.format(newest)))
        shutil.copyfile(newest, 'last_TRIDYN_toBin.h5')


        #TRIDYNFiles='TRIDYN_*.dat'
        heConcFiles='heliumConc_*.dat'
        
        statusFile=open(self.EXIT_STATUS, "r")
        exitStatus=statusFile.read().rstrip('\n')
        statusFile.close()
        print('')


        if exitStatus=='collapsed':
            print('\t simulation exited loop with status collapse')
            
            ##UPDATE: no TRIDYN files generated by Xolotl (info in HDF5)
            ##no need to compress TRIDYN.dat files
            if petscHeConc:
                heConcZipped='allHeliumConc_t%f.zip' %self.driverTime
                zip_ouput='zipHeConcOuput.txt'
                print(('\t \t also save and zip output: {}, regardless of zipOutput \n'.format(heConcFiles)))
                zipString='zip %s %s >> %s ' %(heConcZipped, heConcFiles, zip_ouput)
                subprocess.call([zipString], shell=True)
                rmString='rm '+heConcFiles
                subprocess.call([rmString], shell=True)
            else:
                print(('\t \t no {} generated in this loop \n'.format(heConcFiles)))
                
            print('\t \t rename output files as collapsed before trying again')

            currentXolotlNetworkFile='xolotlStop_%f.h5' %self.driverTime        
            networkFile_unfinished='xolotlStop_%f_collapsed.h5' %self.driverTime
            os.rename(currentXolotlNetworkFile,networkFile_unfinished)

            retentionFile = self.RET_FILE
            rententionUnfinished = 'retention_t%f_collapsed.out' %self.driverTime
            shutil.copyfile(retentionFile,rententionUnfinished)

            surfaceFile=self.SURFACE_FILE
            surfaceUnfinished='surface_t%f_collapsed.txt' %self.driverTime
            shutil.copyfile(surfaceFile,surfaceUnfinished)
            

        elif exitStatus=='overgrid':
            print('\t simulation exited loop with status overgrid')

            ##UPDATE: no TRIDYN files generated by Xolotl (info in HDF5)
            ##no need to compress TRIDYN.dat files
            if petscHeConc:
                heConcZipped='allHeliumConc_t%f_overgrid_%d.zip' %(self.driverTime,n_overgrid_loops)
                zip_ouput='zipHeConcOuput.txt'
                print(('\t \t also save and zip output: {}, regardless of zipOutput \n'.format(heConcFiles)))
                zipString='zip %s %s >> %s ' %(heConcZipped, heConcFiles, zip_ouput)
                subprocess.call([zipString], shell=True)
                rmString='rm '+heConcFiles
                subprocess.call([rmString], shell=True)
            else:
                print(('\t \t no {} generated in this loop \n'.format(heConcFiles)))

            print('\t \t rename output files as overgrid before exiting')

            
            networkFile_unfinished='xolotlStop_%f_overgrid_%d.h5' %(self.driverTime,n_overgrid_loops)
            shutil.copyfile('xolotlStop.h5',networkFile_unfinished)


            xolotlLogFile_overgrid='xolotl_t%f_overgrid_%d.log' %(self.driverTime,n_overgrid_loops)
            shutil.copyfile(xolotlLogFile,xolotlLogFile_overgrid)
            
            retentionFile = self.RET_FILE
            rententionUnfinished = 'retention_t%f_overgrid_%d.out' %(self.driverTime,n_overgrid_loops)
            shutil.copyfile(retentionFile,rententionUnfinished)

            surfaceFile=self.SURFACE_FILE
            surfaceUnfinished='surface_t%f_overgrid_%d.txt' %(self.driverTime,n_overgrid_loops)
            shutil.copyfile(surfaceFile,surfaceUnfinished)
            
            #self.services.error('xolotl_worker: out of grid space, modify network file to continue')
            #raise Exception("Aborting simulation: run out of grid space in Xolotls grid")
            self.services.update_plasma_state()
        else:
            print('\t simulation exited loop with status good (not collapsed or overgrid)')

            ##save TRIDYN_*.dat files zipped OR delete them 
            ##UPDATE: not TRIDYN files generated by Xolotl (info in HDF5)

            #save helium concentration files, zipped
            heConcFiles='heliumConc_*.dat'

            if petscHeConc:
                if zipOutput=='True':
                    heConcZipped='allHeliumConc_t%f.zip' %self.driverTime
                    zip_ouput='zipHeConcOuput.txt'                
                    print(('\t \t save and zip output: {} \n'.format(heConcFiles)))
                    zipString='zip %s %s >> %s ' %(heConcZipped, heConcFiles, zip_ouput)
                    subprocess.call([zipString], shell=True)

                else:
                    print(('\t \t deleting {} (without saving compressed) \n'.format(heConcFiles)))

                rmString='rm '+heConcFiles
                subprocess.call([rmString], shell=True)

            else:
                print(('\t \t no {} generated in this loop \n'.format(heConcFiles)))

        #updates plasma state Xolotl output files
        sys.stdout.flush()
        self.services.update_plasma_state()

    def finalize(self, timeStamp=0.0):
        return
    
