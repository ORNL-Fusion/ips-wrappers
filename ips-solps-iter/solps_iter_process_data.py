#! /usr/bin/env python

from  component import Component
import fileinput
import numpy as np
import solps
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4
import periodic
import math
class solps_iter_data_worker(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):
        return

    def step(self, timeStamp=0.0):
        self.services.stage_plasma_state()
	print 'Hello from solps-data-iter worker'
	dict = {'task':self.TASK,'rmin':float(self.RMIN),
	        'rmax':float(self.RMAX), 'numr':int(self.NUMR), 
	        'zmin':float(self.ZMIN), 'zmax':float(self.ZMAX), 
		'numz':int(self.NUMZ), 
		'solps_geometry':self.SOLPS_GEOMETRY, 
		'solps_equilibrium':self.SOLPS_EQUILIBRIUM, 
		'solps_state':self.SOLPS_STATE, 
		'gitr_geometry':self.GITR_GEOMETRY, 
		'left_target':self.LEFT_TARGET, 'right_target':self.RIGHT_TARGET, 
		'result':self.RESULT, 'dakota_result':self.DAKOTA_RESULT
		}
	dict2 = {
		'rmin_sep':float(self.RMIN_SEP),'rmax_sep':float(self.RMAX_SEP),
		'zmin_sep':float(self.ZMIN_SEP),'zmax_sep':float(self.ZMAX_SEP)
		}
	arguments='-'+' -'.join("%s=%s" % (key,val) for (key,val) in dict.iteritems())
	print 'Regridding SOLPS Data from %s file using solps geometry %s' %(dict['solps_state'], dict['solps_geometry'])
	print 'Making use of data from %s equilibrium file and GITR geometry %s' %(dict['solps_equilibrium'], dict['gitr_geometry'])
	task_id = self.services.launch_task(self.NPROC,
		        self.services.get_working_dir(),
			self.EXE,arguments,task_ppn=1,logfile='task.log')
        #monitor task until complete
        if (self.services.wait_task(task_id)):
	    self.services.error('solps_iter_data_worker: step failed.')
	    return

	f = open(dict['solps_state'],'r')
	txt = f.readlines()[:200]
	f.close()
	
	znLine=0
	zn = ''
	for count,line in enumerate(txt):
	    if znLine:
	        if '*' in line:
		    break
		else:    
	            zn=''.join((zn,line)) 
	    if 'zn' in line:
	        znLine = count
	
	zaminLine=0
	zamin = ''
	for count,line in enumerate(txt):
	    if zaminLine:
	        if '*' in line:
		    break
		else:    
	            zamin=''.join((zamin,line)) 
	    if 'zamin' in line:
	        zaminLine = count
	
	amLine=0
	am = ''
	for count,line in enumerate(txt):
	    if amLine:
	        if '*' in line:
		    break
		else:    
	            am=''.join((am,line)) 
	    if 'am ' in line:
	        amLine = count

	zn = zn.split()
	zn = [float(i) for i in zn]

        zamin = zamin.split()
	zamin = [float(i) for i in zamin]
	
	am = am.split()
	am = [float(i) for i in am]
	#Section to manually add in Tritium
        zn.insert(2,1.0)
        zamin.insert(2,0.0)
        am.insert(2,3.0)
        zn.insert(3,1.0)
        zamin.insert(3,1.0)
        am.insert(3,3.0)
	nIonSpecies = len(zn)
	
	print 'nIonSpecies',nIonSpecies
	species_file = open("speciesList.txt", "w")
	species_file.write('SpeciesIndex   Z   Mass   Charge\n')
	print 'Existing species for current SOLPS run:\n'
	print 'SpeciesIndex   Z   Mass   Charge\n'
        for i in range(nIonSpecies):
	    print '%f       %f        %f      %f \n' %(i,zn[i],am[i],zamin[i])
            species_file.write('%f       %f        %f      %f \n' %(i,zn[i],am[i],zamin[i]))

	species_file.close()
	headers = 'RminusRsep R Z Ti'
	for i in range(len(zn)):
	    if(zn[i] == 1.0 and am[i]==2.0):
	        headers = headers + ' n_{D'+str(int(zamin[i]))+'+}'
	    elif(zn[i] == 1.0 and am[i]==3.0):
		headers = headers + ' n_{T'+str(int(zamin[i]))+'+}'
	    else:
	        headers = headers + ' n_{'+periodic.element(zn[i]).symbol+str(int(zamin[i]))+'+}'
	for i in range(len(zn)):
	    if(zn[i] == 1.0 and am[i]==2.0):
	        headers = headers + ' f_{D'+str(int(zamin[i]))+'+}'
	    elif(zn[i] == 1.0 and am[i]==3.0):
		headers = headers + ' f_{T'+str(int(zamin[i]))+'+}'
	    else:
	        headers = headers + ' f_{'+periodic.element(zn[i]).symbol+str(int(zamin[i]))+'+}'
	headers = headers + ' Te'+ ' ne'
	headers = headers + ' bAngle bMag'
        headers = headers.split()
        
	rt = np.loadtxt(dict['right_target']) 
	r=rt[:,0]
	z=rt[:,1]
        rSep,bAngle,bMag=solps.getBfield(r,z,filename=dict['solps_equilibrium'],geometryFile=dict['gitr_geometry'],rmin=dict2['rmin_sep'],rmax=dict2['rmax_sep'],zmin=dict2['zmin_sep'],zmax=dict2['zmax_sep'])
	print "Returned getBfield Successfully"
	#print 'rSep',rSep
	#Manually Add Tritium to target profiles
	print rt.shape
	print rt[1][0:2]

	ud = np.c_[rSep,rt[:,0:3]]
	ud = np.c_[ud,0.5*rt[:,3]]
	ud = np.c_[ud,0.5*rt[:,4]]
	ud = np.c_[ud,0.5*rt[:,3]]
	ud = np.c_[ud,0.5*rt[:,4]]
	ud = np.c_[ud,rt[:,5:25]]
	ud = np.c_[ud,0.5*rt[:,25]]
	ud = np.c_[ud,0.5*rt[:,26]]
	ud = np.c_[ud,0.5*rt[:,25]]
	ud = np.c_[ud,0.5*rt[:,26]]
	ud = np.c_[ud,rt[:,27:50]]

	ud = np.c_[ud,bAngle]
	ud = np.c_[ud,bMag]
	#print ud
	np.savetxt('solpsTarg.txt',ud)
	for line in fileinput.input('solpsTarg.txt', inplace=True):
	    if fileinput.isfirstline():
	        print ' '.join(headers)
            print line
	plt.close()
	plt.plot(ud[:,0],ud[:,3])
	plt.plot(ud[:,0],ud[:,-4])
	plt.title('ITER Divertor Target Plasma Temperature')
	plt.ylabel('T[eV]')
	plt.xlabel('R-Rsep[m]')
	plt.legend(['Ti','Te'])
	plt.savefig('temp.png')
	plt.close()
	fig = plt.figure()
	ax = plt.subplot(111)
	for y_arr, label in zip(ud[:,4:-3-nIonSpecies].T,headers[4:-3-nIonSpecies]):
	    ax.plot(ud[:,0],y_arr.T,label=label)
	ax.plot(ud[:,0],ud[:,4+2*nIonSpecies],label='ne')
	plt.title('plot title')
	ax.legend()
	plt.title('ITER Divertor Target Plasma Density')
	plt.ylabel('n[m-3]')
	plt.xlabel('R-Rsep[m]')
	plt.yscale('log')
	plt.savefig('dens.png')
	plt.close()
	fig = plt.figure()
	ax = plt.subplot(111)
	#for y_arr, label in zip(ud[:,4+nIonSpecies:4+2*nIonSpecies-3].T,headers[4+nIonSpecies:4+2*nIonSpecies-3]):
	#    ax.plot(ud[:,0],y_arr.T,label=label)
	#    print y_arr
	for i in range(4+nIonSpecies,4+2*nIonSpecies):
	    ax.plot(ud[:,0],np.abs(ud[:,i]),label=headers[i])
        plt.axis([ud[:,0].min(), ud[:,0].max(), 1.0, 1.0e23])
	plt.title('ITER Divertor Target Plasma Flux')
	plt.ylabel('Flux[m-2s-1]')
	plt.xlabel('R-Rsep[m]')
	ax.legend()
	plt.yscale('log')
	plt.savefig('flux.png')
        dak = np.loadtxt('dakota')
	dak= np.reshape(dak,(dict['numr']*dict['numz'],-1))
	print 'dak shape', dak.shape
	rdak = np.unique(dak[:,0])
	print 'rdak',rdak
	rdak = np.unique(dak[:,0])
	zdak = np.unique(dak[:,1])
	te = np.reshape(dak[:,2],(dict['numz'],dict['numr']))
	te[te == -1]=0
	print('te shape', te.shape)
        print(te)
	print('rdak shape', rdak.shape)
        print(rdak[1]-rdak[0])
	print('zdak shape', zdak.shape)
        print(zdak[1]-zdak[0])
	[gradientsR, gradientsZ] = np.gradient(np.array(te),zdak[1]-zdak[0],rdak[1]-rdak[0])
	gradTeZ = gradientsZ
	gradTeR = gradientsR
        plt.close()
        plt.pcolor(rdak,zdak,np.reshape(te,(dict['numz'],dict['numr'])))
        plt.colorbar(orientation='vertical')
        plt.savefig('te.png')
        plt.close()
        plt.pcolor(rdak,zdak,gradTeZ)
        plt.colorbar(orientation='vertical')
        plt.savefig('tez.png')
	#plt.pcolor(rdak,zdak,field1)
	#plt.savefig('field1.png')
	ne = np.reshape(dak[:,3],(dict['numz'],dict['numr']))
	ne[ne == -1]=0
	ti = np.reshape(dak[:,4],(dict['numz'],dict['numr']))
	ti[ti == -1]=0
	[gradientsR,gradientsZ] = np.gradient(np.array(ti),zdak[1]-zdak[0],rdak[1]-rdak[0])
	gradTiZ = gradientsZ
	gradTiR = gradientsR
	ni = np.zeros((nIonSpecies,dict['numz'],dict['numr']))
	vr = np.zeros((nIonSpecies,dict['numz'],dict['numr']))
	vp = np.zeros((nIonSpecies,dict['numz'],dict['numr']))
	vz = np.zeros((nIonSpecies,dict['numz'],dict['numr']))
	niTotal = np.zeros((1,dict['numz'],dict['numr']))
	vrTotal = np.zeros((1,dict['numz'],dict['numr']))
	vpTotal = np.zeros((1,dict['numz'],dict['numr']))
	vzTotal = np.zeros((1,dict['numz'],dict['numr']))
	aveMass = np.zeros((1,dict['numz'],dict['numr']))
	aveCharge = np.zeros((1,dict['numz'],dict['numr']))
	offset=5
	nIonSpecies = 21
	print 'starting ionspec loop',nIonSpecies
	for i in range(nIonSpecies):
	    if i==0:
	        ni[0,:,:] = np.reshape(0.5*dak[:,offset+i],(dict['numz'],dict['numr']))
	        ni[2,:,:] = np.reshape(0.5*dak[:,offset+i],(dict['numz'],dict['numr']))
	        vr[0,:,:] = np.reshape(dak[:,offset+nIonSpecies+i],(dict['numz'],dict['numr']))
	        vr[2,:,:] = np.reshape(dak[:,offset+nIonSpecies+i],(dict['numz'],dict['numr']))
	        vp[0,:,:] = np.reshape(dak[:,offset+2*nIonSpecies+i],(dict['numz'],dict['numr']))
	        vp[2,:,:] = np.reshape(dak[:,offset+2*nIonSpecies+i],(dict['numz'],dict['numr']))
	        vz[0,:,:] = np.reshape(dak[:,offset+3*nIonSpecies+i],(dict['numz'],dict['numr']))
	        vz[2,:,:] = np.reshape(dak[:,offset+3*nIonSpecies+i],(dict['numz'],dict['numr']))
	    elif i==1:
	        ni[1,:,:] = np.reshape(0.5*dak[:,offset+i],(dict['numz'],dict['numr']))
	        ni[3,:,:] = np.reshape(0.5*dak[:,offset+i],(dict['numz'],dict['numr']))
	        vr[1,:,:] = np.reshape(dak[:,offset+nIonSpecies+i],(dict['numz'],dict['numr']))
	        vr[3,:,:] = np.reshape(dak[:,offset+nIonSpecies+i],(dict['numz'],dict['numr']))
	        vp[1,:,:] = np.reshape(dak[:,offset+2*nIonSpecies+i],(dict['numz'],dict['numr']))
	        vp[3,:,:] = np.reshape(dak[:,offset+2*nIonSpecies+i],(dict['numz'],dict['numr']))
	        vz[1,:,:] = np.reshape(dak[:,offset+3*nIonSpecies+i],(dict['numz'],dict['numr']))
	        vz[3,:,:] = np.reshape(dak[:,offset+3*nIonSpecies+i],(dict['numz'],dict['numr']))
	    else:
	        ni[i+2,:,:] = np.reshape(dak[:,offset+i],(dict['numz'],dict['numr']))
	        vr[i+2,:,:] = np.reshape(dak[:,offset+nIonSpecies+i],(dict['numz'],dict['numr']))
	        vp[i+2,:,:] = np.reshape(dak[:,offset+2*nIonSpecies+i],(dict['numz'],dict['numr']))
	        vz[i+2,:,:] = np.reshape(dak[:,offset+3*nIonSpecies+i],(dict['numz'],dict['numr']))
	print "finished ionspec loop"    
        #plt.close()
        #plt.figure(1,figsize=(10, 6), dpi=2000)
        #plotsize = math.ceil(nIonSpecies**(0.5))
        #for i in range(nIonSpecies):
        #    plt.subplot(plotsize,plotsize,i+1)
        ##plot2dGeom('../2d/input/iter2dRefinedOuterTarget.cfg')
        ##plt.title("ITER W Impurity Density")
        ##plt.xlabel("r [m]")
        ##plt.ylabel("z [m]")
        #    plt.pcolor(rdak,zdak,np.reshape(ni[i,:,:],(dict['numz'],dict['numr'])))
        #    plt.colorbar(orientation='vertical')
        #plt.savefig('image1.png')
        print('Starting total loop')
	for i in range(nIonSpecies):
	    if zamin[i] > 0.0:
	        #print niTotal.shape
		#print ni[i,:,:].shape
	        niTotal = niTotal + np.reshape(ni[i,:,:],(1,dict['numz'],dict['numr']))
	        aveMass = aveMass+np.reshape(am[i]*ni[i,:,:],(1,dict['numz'],dict['numr']))
	        aveCharge = aveCharge+np.reshape(zamin[i]*ni[i,:,:],(1,dict['numz'],dict['numr']))
	        vrTotal = vrTotal+np.reshape(np.multiply(vr[i,:,:],ni[i,:,:]),(1,dict['numz'],dict['numr']))
	        vpTotal = vpTotal+np.reshape(np.multiply(vp[i,:,:],ni[i,:,:]),(1,dict['numz'],dict['numr']))
	        vzTotal = vzTotal+np.reshape(np.multiply(vz[i,:,:],ni[i,:,:]),(1,dict['numz'],dict['numr']))
        
        print('Finished total loop')
	aveMass = np.divide(aveMass,niTotal)
	aveMass[np.reshape(ni[i,:,:],(1,dict['numz'],dict['numr'])) == -1]=0
	aveCharge= np.divide(aveCharge,niTotal)
	aveCharge[np.reshape(ni[i,:,:],(1,dict['numz'],dict['numr'])) == -1]=0
	niTotal[np.reshape(ni[i,:,:],(1,dict['numz'],dict['numr'])) == -1]=0
	vrTotal = np.divide(vrTotal,niTotal)
	vrTotal[np.reshape(ni[i,:,:],(1,dict['numz'],dict['numr'])) == -1]=0
	vpTotal = np.divide(vpTotal,niTotal)
	vpTotal[np.reshape(ni[i,:,:],(1,dict['numz'],dict['numr'])) == -1]=0
	vzTotal = np.divide(vzTotal,niTotal)
	vzTotal[np.reshape(ni[i,:,:],(1,dict['numz'],dict['numr'])) == -1]=0
        print('Finished divides')
	#print 'finished loop'	
	#print rdak
	#print zdak
	#print niTotal.shape
	#print niTotal
	plt.close()
	plt.pcolor(rdak,zdak,np.reshape(niTotal,(dict['numz'],dict['numr'])))
	plt.colorbar()
	plt.savefig('niTotal.png')
	plt.close()
	plt.pcolor(rdak,zdak,np.reshape(aveMass,(dict['numz'],dict['numr'])))
	plt.colorbar()
	plt.savefig('aveMass.png')
	plt.close()
	plt.pcolor(rdak,zdak,np.reshape(aveCharge,(dict['numz'],dict['numr'])))
	plt.colorbar()
	plt.savefig('aveCharge.png')
        
	    
        print('Finished plots, getting target bfields')
	br = np.reshape(dak[:,offset+4*nIonSpecies],(dict['numz'],dict['numr']))
	br[br == -1]=0
	bphi = np.reshape(dak[:,offset+4*nIonSpecies+1],(dict['numz'],dict['numr']))
	bphi[bphi == -1]=0
	bz = np.reshape(dak[:,offset+4*nIonSpecies+2],(dict['numz'],dict['numr']))
	bz[bz == -1]=0
	pot = np.reshape(dak[:,offset+4*nIonSpecies+3],(dict['numz'],dict['numr']))
	pot[pot == -1]=0
        [gradz,gradr] = np.gradient(np.array(pot),zdak[1]-zdak[0],rdak[1]-rdak[0])
	Ez = -gradz
        Er = -gradr
	rootgrp = netCDF4.Dataset("profiles.nc", "w", format="NETCDF4")
	nrr = rootgrp.createDimension("nR", len(rdak))
	nzz = rootgrp.createDimension("nZ", len(zdak))
	brr = rootgrp.createVariable("br","f8",("nZ","nR"))
	btt = rootgrp.createVariable("bt","f8",("nZ","nR"))
	bzz = rootgrp.createVariable("bz","f8",("nZ","nR"))
	rr = rootgrp.createVariable("r","f8",("nR"))
	zz = rootgrp.createVariable("z","f8",("nZ"))
	brr[:] = br
	btt[:] = bphi
	bzz[:] = bz
	rr[:] = rdak
	zz[:] = zdak
	tee = rootgrp.createVariable("te","f8",("nZ","nR"))
	nee = rootgrp.createVariable("ne","f8",("nZ","nR"))
	tii = rootgrp.createVariable("ti","f8",("nZ","nR"))
	nii = rootgrp.createVariable("ni","f8",("nZ","nR"))
	mass = rootgrp.createVariable("mass","f8",("nZ","nR"))
	charge = rootgrp.createVariable("charge","f8",("nZ","nR"))
	vrr = rootgrp.createVariable("vr","f8",("nZ","nR"))
	vzz = rootgrp.createVariable("vz","f8",("nZ","nR"))
	vpp = rootgrp.createVariable("vp","f8",("nZ","nR"))
	pott = rootgrp.createVariable("pot","f8",("nZ","nR"))
	Err = rootgrp.createVariable("Er","f8",("nZ","nR"))
	Ett = rootgrp.createVariable("Et","f8",("nZ","nR"))
	Ezz = rootgrp.createVariable("Ez","f8",("nZ","nR"))
	teer = rootgrp.createVariable("gradTeR","f8",("nZ","nR"))
	teez = rootgrp.createVariable("gradTeZ","f8",("nZ","nR"))
	teey = rootgrp.createVariable("gradTeY","f8",("nZ","nR"))
	tiir = rootgrp.createVariable("gradTiR","f8",("nZ","nR"))
	tiiz = rootgrp.createVariable("gradTiZ","f8",("nZ","nR"))
	tiiy = rootgrp.createVariable("gradTiY","f8",("nZ","nR"))
	tee[:] = te
	nee[:] = ne
	tii[:] = ti
	nii[:] = niTotal
	mass[:] = aveMass
	charge[:] = aveCharge
	vrr[:] = vrTotal
	vpp[:] = vpTotal
	vzz[:] = vzTotal
	pott[:] = pot
	Err[:] = Er
	Ett[:] = 0*Er
	Ezz[:] = Ez
	teer[:] = gradTeR
	teey[:] = 0*gradTeR
	teez[:] = gradTeZ
	tiir[:] = gradTiR
	tiiy[:] = 0*gradTiR
	tiiz[:] = gradTiZ
	rootgrp.close()
        self.services.update_plasma_state()        
        return
    
    def finalize(self, timeStamp=0.0):
        return
