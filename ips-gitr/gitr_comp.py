#! /usr/bin/env python

from  component import Component
import gitr
import os
import pickle
import shutil
import gitrParticleSource
import hpic

class gitr_comp(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):
        #Set up input deck
        print('input dir and cwd this now', self.BASE_DIR, ' ', os.getcwd())
        gitr.copy_folder(self.GITR_INPUT_DIR,os.getcwd())
        gitr.modifyInputParam(nT=int(self.NT))
        return

    def step(self, timeStamp=0.0):
        print 'Hello from gitr_comp'
        self.services.stage_plasma_state()
	if os.path.exists('bField.nc'):
            shutil.copyfile('bField.nc','input/bField.nc')
	if os.path.exists('ftridynSelf.nc'):
            shutil.copyfile('ftridynSelf.nc','input/ftridynSelf.nc')
	if os.path.exists('profiles.nc'):
            shutil.copyfile('profiles.nc','input/profiles.nc')
        hpic.plot_hpic_ieadDavide(solps_path=str(self.HPIC_SOLPS_FILE),HpicDataFolder=str(self.HPIC_DATA_DIR))
        hpic.computeSputtYld()
        bgFlux,flux_fracs = hpic.printBackgroundDist()
        fluxPerParticle = gitrParticleSource.particleSource(geomFile=str(self.GEOM_FILE),nParticles=int(self.NP)) 
	if os.path.exists('particleSource.nc'):
            shutil.copyfile('particleSource.nc','input/particleSource.nc')
        task_id = self.services.launch_task(self.NPROC,
                                            self.services.get_working_dir(),
                                            self.GITR_EXE,task_ppn=self.TASK_PPN,
                                            logfile='gitr.log') #,ppn=1)
        #monitor task until complete
        if (self.services.wait_task(task_id)):
            self.services.error('gitr_comp: step failed.')
        
        #gitr.piscesProcessing(path=self.BASE_DIR)
	#gitr.iter2dProcessing()
	nLocations = gitr.iter3dProcessingQ4(nParticles=int(self.NP),totalParticleRate=fluxPerParticle,path = '',locRmRs=[-0.1,0.02,0.09,0.2],flux_fracs =flux_fracs,locWidth = 0.02)
        self.services.update_plasma_state()         
        for i in range(nLocations):
            shutil.copyfile('gitrOut'+str(i)+'.txt',self.BASE_DIR+'/gitrOut'+str(i)+'.txt')
        return
    
    def finalize(self, timeStamp=0.0):
        return
    
