#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import os
import sys
import h5py
import numpy      as npy
import traceback  as traceback
import mathtools  as mathtools

from scipy.optimize    import curve_fit
from scipy.integrate   import trapz
from scipy.interpolate import interp1d,interp2d
from scipy.interpolate import CubicSpline,RectBivariateSpline


if   sys.version_info.major == 3:
     PYTHON3 = True; PYTHON2 = False
elif sys.version_info.major == 2:
     PYTHON2 = True; PYTHON3 = False

mu0 = 4.0e-7*npy.pi


def read_csv(csvfn):
    with open(csvfn, mode='r') as csvfid:
         csvdata = csv.DictReader(csvfid)
         recordid = 0
         csvdict  = {}
         for record in csvdata:
             headers = record.keys()
             for iheader in sorted(headers):
                 if recordid < len(headers):
                    csvdict[iheader] = []
                    recordid += 1
                 try:
                    csvdict[iheader].append(int(record[iheader]))
                 except ValueError:
                    csvdict[iheader].append(npy.float64(record[iheader]))
         for iheader in csvdict:
             csvdict[iheader] = npy.array(csvdict[iheader])
    return csvdict

def write_csv(csvfn,csvdict):
    headers   = sorted(csvdict.keys())
    nheaders  = npy.size(headers)
    nrows     = npy.size(csvdict[headers[0]])

    with open(csvfn, mode='w') as csvfid:
         writer = csv.DictWriter(csvfid, fieldnames=headers)
         writer.writeheader()
         for irow in range(nrows):
             record = {}
             for iheader in headers:
                 if   type(csvdict[iheader][irow]) == int:
                      record.update({iheader:"%d" % csvdict[iheader][irow]})
                 elif type(csvdict[iheader][irow]) == float:
                      record.update({iheader:"%8.4e" % csvdict[iheader][irow]})
                 elif type(csvdict[iheader][irow]) == npy.float64:
                      record.update({iheader:"%8.4e" % csvdict[iheader][irow]})
                 elif type(csvdict[iheader][irow]) == npy.float128:
                      record.update({iheader:"%8.4e" % csvdict[iheader][irow]})
             writer.writerow(record)
    return 1


def read_chease(fpath,setParam={},**kwargs):
    if os.path.isfile(fpath) == False:
       errorFunc = traceback.extract_stack(limit=2)[-2][3]
       errorLine = traceback.extract_stack(limit=2)[-2][1]
       errorFile = traceback.extract_stack(limit=2)[-2][2]
       errMSG    = 'Call %s line %5d in file %s Failed.\n'
       errMSG   += 'Fatal: file %s not found.'
       raise IOError(errMSG %(errorFunc,errorLine,errorFile,fpath))

    if   h5py.is_hdf5(fpath):
         CHEASEdata = read_chease_hdf(fpath,setParam,**kwargs)
    else:
         CHEASEdata = read_chease_dat(fpath,setParam,**kwargs)

    return CHEASEdata

def read_chease_dat(fpath,setParam={},**kwargs):
    if os.path.isfile(fpath) == False:
       errorFunc = traceback.extract_stack(limit=2)[-2][3]
       errorLine = traceback.extract_stack(limit=2)[-2][1]
       errorFile = traceback.extract_stack(limit=2)[-2][2]
       errMSG    = 'Call %s line %5d in file %s Failed.\n'
       errMSG   += 'Fatal: file %s not found.'
       raise IOError(errMSG %(errorFunc,errorLine,errorFile,fpath))

    if 'Normalized' in setParam: Normalized = setParam['Normalized']
    else:                        Normalized = False

    fh = open(fpath,"r")
    fcontents = fh.readlines()
    fh.close()
    
    CHEASEdata = {}
    
    CHEASEdata['NPSI']   = int(fcontents[1])
    CHEASEdata['NCHI']   = int(fcontents[3])
    CHEASEdata['NRBOX']  = int(fcontents[9])
    CHEASEdata['NZBOX']  = int(fcontents[11])
    CHEASEdata['NBOUND'] = int(fcontents[13])
    
    CHEASEdata['R0EXP']  = float(fcontents[5])
    CHEASEdata['B0EXP']  = float(fcontents[7])
    
    
    NPSI = CHEASEdata['NPSI']
    NCHI = CHEASEdata['NCHI']
    
    NRBOX = CHEASEdata['NRBOX']
    NZBOX = CHEASEdata['NZBOX']
    
    NBOUND = CHEASEdata['NBOUND']
    
    
    fieldIN = []
    fieldIN.extend(['NPSI','NCHI','NRBOX','NZBOX','NBOUND'])
    
    fieldFN = []
    fieldFN.extend(['R0EXP','B0EXP'])
    
    field1D = []
    field1D.extend(['PSI','CHI','Rgeom','ageom','q','dqdpsi','d2qdpsi2','p','dpdpsi'])
    field1D.extend(['f','fdfdpsi','V','rho_t','shear','dsheardpsi','kappa'])
    field1D.extend(['delta_lower','delta_upper','dVdpsi','dpsidrhotor','GDPSI_av'])
    field1D.extend(['radius_av','R_av','TE','DTEDPSI','NE','DNEDPSI'])
    field1D.extend(['TI','DTIDPSI','NI','DNIDPSI','ZEFF','SIGNEO','JBSBAV'])
    
    fieldEX = []
    fieldEX.extend(['Rmesh','Zmesh','RBOUNDplasma','ZBOUNDplasma'])
    fieldEX.extend(['wallposr','wallposz'])
    
    field2D = []
    field2D.extend(['g11','g12','g22','g33'])
    field2D.extend(['R','Z','B','Jacobian','dBdchi','dBdpsi'])
    field2D.extend(['dChidR','dChidZ','dPsidR','dPsidZ','psiRZ','chiRZ'])

    fields = {}
    
    for icount in range(len(fcontents)):
        if fcontents[icount].strip() in fieldIN:
           fields.update({fcontents[icount].strip():icount})
        if fcontents[icount].strip() in fieldFN:
           fields.update({fcontents[icount].strip():icount})
        if fcontents[icount].strip() in field1D:
           fields.update({fcontents[icount].strip():icount})
        if fcontents[icount].strip() in field2D:
           fields.update({fcontents[icount].strip():icount})
        if fcontents[icount].strip() in fieldEX:
           fields.update({fcontents[icount].strip():icount})
    
    for ifield in fields:
        if ifield in field1D:
           ncount = fields[ifield]+1; iloop=0
           if ifield == 'CHI':
              CHEASEdata[ifield] = npy.zeros(NCHI)
              bgnind = ncount
              endind = ncount+int(npy.ceil(NCHI/5))
           else:
              CHEASEdata[ifield] = npy.zeros(NPSI)
              bgnind = ncount
              endind = ncount+int(npy.ceil(NPSI/5))
           for ind in range(bgnind,endind):
               for irecord in fcontents[ind].split():
                   try:
                       CHEASEdata[ifield][iloop] = float(irecord); iloop+=1
                   except ValueError:
                       continue
    
        if ifield in fieldEX:
           ncount = fields[ifield]+1; iloop=0
           if   ifield == 'Rmesh':
                CHEASEdata[ifield] = npy.zeros(NRBOX)
                bgnind = ncount
                endind = ncount+int(npy.ceil(NRBOX/5))
           elif ifield == 'Zmesh':
                CHEASEdata[ifield] = npy.zeros(NZBOX)
                bgnind = ncount
                endind = ncount+int(npy.ceil(NZBOX/5))
           elif ifield in ['RBOUNDplasma','ZBOUNDplasma']:
                CHEASEdata[ifield] = npy.zeros(NBOUND)
                bgnind = ncount
                endind = ncount+int(npy.ceil(NBOUND/5))
           elif ifield in ['wallposr','wallposz']:
                CHEASEdata[ifield] = npy.zeros(5)
                bgnind = ncount
                endind = ncount+int(npy.ceil(5/5))
           for ind in range(bgnind,endind):
               for irecord in fcontents[ind].split():
                   try:
                       CHEASEdata[ifield][iloop] = float(irecord); iloop+=1
                   except ValueError:
                       continue
     
        if ifield in field2D:
           ncount = fields[ifield]+1; iloop=0
           if ifield in ['psiRZ','chiRZ']:
              CHEASEdata[ifield] = npy.zeros(NRBOX*NZBOX)
              bgnind = ncount
              endind = ncount+int(npy.ceil(NRBOX*NZBOX/5))
           else:
              CHEASEdata[ifield] = npy.zeros(NPSI*NCHI)
              bgnind = ncount
              endind = ncount+int(npy.ceil(NPSI*NCHI/5))
           for ind in range(bgnind,endind):
               for irecord in fcontents[ind].split():
                   try:
                       CHEASEdata[ifield][iloop] = float(irecord); iloop+=1
                   except ValueError:
                       continue
           if ifield in ['psiRZ','chiRZ']:
              CHEASEdata[ifield] = npy.reshape(CHEASEdata[ifield],[NZBOX,NRBOX])
           else:
              CHEASEdata[ifield] = npy.reshape(CHEASEdata[ifield],[NCHI,NPSI])


    extendCHI  = False
    reverseCHI = False
    if CHEASEdata['CHI'][0]>CHEASEdata['CHI'][1]:
       CHEASEdata['CHI']      = CHEASEdata["CHI"][::-1]
       reverseCHI             = True
    if min(abs(CHEASEdata['CHI']-2.0*npy.pi))>=min(npy.diff(CHEASEdata['CHI'])):
       CHEASEdata['CHI']      = npy.append(CHEASEdata['CHI'],2.0*npy.pi)
       extendCHI              = True

    if extendCHI:
       CHEASEdata['R']        = npy.vstack((CHEASEdata['R'],CHEASEdata['R'][0,:]))
       CHEASEdata['Z']        = npy.vstack((CHEASEdata['Z'],CHEASEdata['Z'][0,:]))
       CHEASEdata['B']        = npy.vstack((CHEASEdata['B'],CHEASEdata['B'][0,:]))
       CHEASEdata['g11']      = npy.vstack((CHEASEdata['g11'],CHEASEdata['g11'][0,:]))
       CHEASEdata['g22']      = npy.vstack((CHEASEdata['g22'],CHEASEdata['g22'][0,:]))
       CHEASEdata['g33']      = npy.vstack((CHEASEdata['g33'],CHEASEdata['g33'][0,:]))
       CHEASEdata['dBdpsi']   = npy.vstack((CHEASEdata['dBdpsi'],CHEASEdata['dBdpsi'][0,:]))
       CHEASEdata['dBdchi']   = npy.vstack((CHEASEdata['dBdchi'],CHEASEdata['dBdchi'][0,:]))
       CHEASEdata['dChidZ']   = npy.vstack((CHEASEdata['dChidZ'],CHEASEdata['dChidZ'][0,:]))
       CHEASEdata['dPsidZ']   = npy.vstack((CHEASEdata['dPsidZ'],CHEASEdata['dPsidZ'][0,:]))
       CHEASEdata['dChidR']   = npy.vstack((CHEASEdata['dChidR'],CHEASEdata['dChidR'][0,:]))
       CHEASEdata['dPsidR']   = npy.vstack((CHEASEdata['dPsidR'],CHEASEdata['dPsidR'][0,:]))
       CHEASEdata['Jacobian'] = npy.vstack((CHEASEdata['Jacobian'],CHEASEdata['Jacobian'][0,:]))

    CHEASEdata['rhopsi_SI']   = npy.sqrt(CHEASEdata['PSI']/CHEASEdata['PSI'][-1])
    CHEASEdata['rhotor_SI']   = CHEASEdata["rho_t"]

    CHEASEdata['Te']      = CHEASEdata['TE']
    del CHEASEdata['TE']
    CHEASEdata['ne']      = CHEASEdata['NE']
    del CHEASEdata['NE']
    CHEASEdata['Ti']      = CHEASEdata['TI']
    del CHEASEdata['TI']
    CHEASEdata['ni']      = CHEASEdata['NI']
    del CHEASEdata['NI']
    CHEASEdata['Zeff']    = CHEASEdata['ZEFF']
    del CHEASEdata['ZEFF']
    CHEASEdata['TePrime'] = CHEASEdata['DTEDPSI']
    del CHEASEdata['DTEDPSI']
    CHEASEdata['nePrime'] = CHEASEdata['DNEDPSI']
    del CHEASEdata['DNEDPSI']
    CHEASEdata['TiPrime'] = CHEASEdata['DTIDPSI']
    del CHEASEdata['DTIDPSI']
    CHEASEdata['niPrime'] = CHEASEdata['DNIDPSI']
    del CHEASEdata['DNIDPSI']

    CHEASEdata['Zi']          = 1.0
    CHEASEdata['Zz']          = 6.0
    CHEASEdata['nz']          = CHEASEdata['Zeff']*CHEASEdata['ne']
    CHEASEdata['nz']         -= CHEASEdata['ni']*CHEASEdata['Zi']**2
    CHEASEdata['nz']         /= CHEASEdata['Zz']**2

    CHEASEdata['J']        = CHEASEdata['Jacobian']
    del CHEASEdata['Jacobian']
    CHEASEdata['Jbs']      = CHEASEdata['JBSBAV']
    del CHEASEdata['JBSBAV']
    CHEASEdata['signeo']   = CHEASEdata['SIGNEO']
    del CHEASEdata['SIGNEO']
    CHEASEdata['Volume']   = CHEASEdata['V']
    del CHEASEdata['V']
    CHEASEdata['pressure'] = CHEASEdata['p']
    del CHEASEdata['p']
    CHEASEdata['pprime']   = 2.0*npy.pi*CHEASEdata['dpdpsi']
    del CHEASEdata['dpdpsi']
    CHEASEdata['ffprime']  = 2.0*npy.pi*CHEASEdata['fdfdpsi']
    del CHEASEdata['fdfdpsi']
    CHEASEdata['fprime']   = CHEASEdata['ffprime']/CHEASEdata['f']

    CHEASEdata['rmesh']    = CHEASEdata['Rmesh']
    del CHEASEdata['Rmesh']
    CHEASEdata['zmesh']    = CHEASEdata['Zmesh']
    del CHEASEdata['Zmesh']
    CHEASEdata['rbound']    = CHEASEdata['RBOUNDplasma']
    del CHEASEdata['RBOUNDplasma']
    CHEASEdata['zbound']    = CHEASEdata['ZBOUNDplasma']
    del CHEASEdata['ZBOUNDplasma']

    CHEASEdata['C0']       = npy.trapz(y=CHEASEdata['J']/CHEASEdata['R'],                    x=CHEASEdata['CHI'],axis=0)
    CHEASEdata['C1']       = npy.trapz(y=CHEASEdata['J'],                                    x=CHEASEdata['CHI'],axis=0)
    CHEASEdata['C2']       = npy.trapz(y=CHEASEdata['J']/CHEASEdata['R']**2,                 x=CHEASEdata['CHI'],axis=0)
    CHEASEdata['C3']       = npy.trapz(y=CHEASEdata['J']*CHEASEdata['g11']*CHEASEdata['g33'],x=CHEASEdata['CHI'],axis=0)

    CHEASEdata['y1']       = 1.0+CHEASEdata['C3']/CHEASEdata['C2']/CHEASEdata['f']**2/4.0/npy.pi**2

    CHEASEdata['<B2>']     = npy.trapz(y=CHEASEdata['J']*CHEASEdata['B']**2,x=CHEASEdata['CHI'],axis=0)/CHEASEdata['C1']
    CHEASEdata['<JdotB>']  =-CHEASEdata['f']*CHEASEdata['pprime']-CHEASEdata['fprime']*CHEASEdata['<B2>']/mu0
    CHEASEdata['Jprl']     = CHEASEdata['<JdotB>']/CHEASEdata['B0EXP']

    CHEASEdata['<T/R2>']   = npy.trapz(y=CHEASEdata['J']*CHEASEdata['f']*CHEASEdata['g33'],x=CHEASEdata['CHI'],axis=0)/CHEASEdata['C1']
    CHEASEdata['Iprl']     = CHEASEdata['R0EXP']*CHEASEdata['<JdotB>']/CHEASEdata['<T/R2>']

    CHEASEdata['Istr']     =-((CHEASEdata['C2']/CHEASEdata['C0'])*(CHEASEdata['ffprime']/mu0))
    CHEASEdata['Istr']    +=-((CHEASEdata['C1']/CHEASEdata['C0'])*CHEASEdata['pprime'])
    CHEASEdata['Istr']    *= CHEASEdata['R0EXP']**2

    CHEASEdata['Jphi']     =-(CHEASEdata['R']*CHEASEdata['pprime'])-(CHEASEdata['ffprime']/(mu0*CHEASEdata['R']))

    CHEASEdata['Jtor']     = npy.trapz(y=CHEASEdata['Jphi']*CHEASEdata['J']/CHEASEdata['R'],x=CHEASEdata['CHI'],axis=0)
    CHEASEdata['Itor']     = npy.trapz(y=CHEASEdata['Jtor'],x=CHEASEdata['PSI'],axis=0)


    CHEASEdata['Ibs']      = CHEASEdata['R0EXP']*CHEASEdata['Jbs']/CHEASEdata['<T/R2>']
    CHEASEdata['Iohmic']   = CHEASEdata['Iprl']-CHEASEdata['Ibs']
    CHEASEdata['Johmic']   = CHEASEdata['Jprl']-CHEASEdata['Jbs']

    extendPSI    = npy.linspace(CHEASEdata['PSI'][0],CHEASEdata['PSI'][-1],10*npy.size(CHEASEdata['PSI']))
    extendPHI    = npy.empty_like(extendPSI)
    extendPHI[0] = 0.0
    qfunc        = CubicSpline(CHEASEdata['PSI'],CHEASEdata['q'])
    for i in range(1,npy.size(extendPSI)):
        x           = extendPSI[:i+1]
        y           = qfunc(x)
        extendPHI[i]= npy.trapz(y,x)

    CHEASEdata['PHI'] = npy.empty_like(CHEASEdata['PSI'])
    phifunc           = CubicSpline(extendPSI,extendPHI)
    for i in range(npy.size(CHEASEdata['PSI'])):
        CHEASEdata['PHI'][i] = phifunc(CHEASEdata['PSI'][i])

    CHEASEdata['PHIN']     = (CHEASEdata['PHI']-CHEASEdata['PHI'][0])/(CHEASEdata['PHI'][-1]-CHEASEdata['PHI'][0])
    CHEASEdata['PSIN']     = (CHEASEdata['PSI']-CHEASEdata['PSI'][0])/(CHEASEdata['PSI'][-1]-CHEASEdata['PSI'][0])
    CHEASEdata['CHIN']     = (CHEASEdata['CHI']-CHEASEdata['CHI'][0])/(CHEASEdata['CHI'][-1]-CHEASEdata['CHI'][0])

   #CHEASEdata['rhopsi']   = (CHEASEdata['rhopsi_SI']-CHEASEdata['rhopsi_SI'][0])/(CHEASEdata['rhopsi_SI'][-1]-CHEASEdata['rhopsi_SI'][0])
   #CHEASEdata['rhotor']   = (CHEASEdata['rhotor_SI']-CHEASEdata['rhotor_SI'][0])/(CHEASEdata['rhotor_SI'][-1]-CHEASEdata['rhotor_SI'][0])
    CHEASEdata['rhopsi']   = npy.sqrt(CHEASEdata['PSIN'])
    CHEASEdata['rhotor']   = npy.sqrt(CHEASEdata['PHIN'])

    #Implementing Interpolation to EQDSK Grid
    rhopsiflag = False; rhotorflag = False
    if 'nrhomesh' in setParam:
       if   setParam['nrhomesh'] in [0,'rhopsi']: rhopsiflag = True
       elif setParam['nrhomesh'] in [1,'rhotor']: rhotorflag = True
    else:                                         rhopsiflag = True

    eqdskflag    = False
    expeqflag    = False
    interpflag   = False
    importedflag = False
    for key,value in kwargs.items():
        if   key in ['eqdsk','eqdskdata','eqdskfpath']:
             if    type(value)==str and os.path.isfile(value.strip()):
                   eqdskdata = read_eqdsk(fpath=value.strip())
             else: raise IOError('%s file not found!' % value.strip())
             if 'rhopsi' in eqdskdata: rhopsi = eqdskdata['rhopsi'][:]
             if 'rhotor' in eqdskdata: rhotor = eqdskdata['rhotor'][:]
             if 'PSIN'   in eqdskdata: psi    = eqdskdata['PSIN'][:];  interpflag = True
             if 'PHIN'   in eqdskdata: phi    = eqdskdata['PHIN'][:];  interpflag = True
             eqdskflag = True
        elif key in ['imported','external','other']:
             imported = value.copy()
             if 'rhopsi' in imported: rhopsi = imported['rhopsi'][:]
             if 'rhotor' in imported: rhotor = imported['rhotor'][:]
             if 'PSIN'   in imported: psi = imported['PSIN'][:];       interpflag = True
             else:                    psi = importeddata['rhopsi']**2; interpflag = True
             if 'PHIN'   in imported: phi = imported['PHIN'][:];       interpflag = True
             else:                    phi = importeddata['rhotor']**2; interpflag = True
             importedflag = True

    if interpflag:
       if   rhopsiflag:
            CHEASEdata['q']        = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['q'],psi)
            CHEASEdata['f']        = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['f'],psi)
            CHEASEdata['Te']       = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Te'],psi)
            CHEASEdata['Ti']       = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Ti'],psi)
            CHEASEdata['ne']       = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['ne'],psi)
            CHEASEdata['ni']       = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['ni'],psi)
            CHEASEdata['nz']       = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['nz'],psi)
            CHEASEdata['Ibs']      = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Ibs'],psi)
            CHEASEdata['Jbs']      = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Jbs'],psi)
            CHEASEdata['Zeff']     = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Zeff'],psi)
            CHEASEdata['Istr']     = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Istr'],psi)
            CHEASEdata['Iprl']     = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Iprl'],psi)
            CHEASEdata['Jprl']     = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Jprl'],psi)
            CHEASEdata['kappa']    = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['kappa'],psi)
            CHEASEdata['shear']    = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['shear'],psi)
            CHEASEdata['signeo']   = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['signeo'],psi)
            CHEASEdata['fprime']   = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['fprime'],psi)
            CHEASEdata['Iohmic']   = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Iohmic'],psi)
            CHEASEdata['Johmic']   = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Johmic'],psi)
            CHEASEdata['pprime']   = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['pprime'],psi)
            CHEASEdata['ffprime']  = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['ffprime'],psi)
            CHEASEdata['pressure'] = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['pressure'],psi)

            CHEASEdata['PSIN']     = psi[:]
            CHEASEdata['PHIN']     = phi[:]
            CHEASEdata['rhopsi']   = rhopsi[:]
            CHEASEdata['rhotor']   = rhotor[:]

       elif rhotorflag:
            CHEASEdata['q']        = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['q'],psi,phi,phi)
            CHEASEdata['f']        = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['f'],psi,phi,phi)
            CHEASEdata['Te']       = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Te'],psi,phi,phi)
            CHEASEdata['Ti']       = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Ti'],psi,phi,phi)
            CHEASEdata['ne']       = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['ne'],psi,phi,phi)
            CHEASEdata['ni']       = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['ni'],psi,phi,phi)
            CHEASEdata['nz']       = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['nz'],psi,phi,phi)
            CHEASEdata['Ibs']      = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Ibs'],psi,phi,phi)
            CHEASEdata['Jbs']      = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Jbs'],psi,phi,phi)
            CHEASEdata['Zeff']     = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Zeff'],psi,phi,phi)
            CHEASEdata['Istr']     = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Istr'],psi,phi,phi)
            CHEASEdata['Iprl']     = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Iprl'],psi,phi,phi)
            CHEASEdata['Jprl']     = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Jprl'],psi,phi,phi)
            CHEASEdata['kappa']    = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['kappa'],psi,phi,phi)
            CHEASEdata['shear']    = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['shear'],psi,phi,phi)
            CHEASEdata['signeo']   = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['signeo'],psi,phi,phi)
            CHEASEdata['fprime']   = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['fprime'],psi,phi,phi)
            CHEASEdata['Iohmic']   = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Iohmic'],psi,phi,phi)
            CHEASEdata['Johmic']   = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Johmic'],psi,phi,phi)
            CHEASEdata['pprime']   = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['pprime'],psi,phi,phi)
            CHEASEdata['ffprime']  = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['ffprime'],psi,phi,phi)
            CHEASEdata['pressure'] = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['pressure'],psi,phi,phi)

            CHEASEdata['PSIN']     = psi[:]
            CHEASEdata['PHIN']     = phi[:]
            CHEASEdata['rhotor']   = rhopsi[:]
            CHEASEdata['rhopsi']   = rhotor[:]

    if Normalized:
       CHEASEdata['R']        = CHEASEdata['R']/CHEASEdata['R0EXP']
       CHEASEdata['Z']        = CHEASEdata['Z']/CHEASEdata['R0EXP']
       CHEASEdata['B']        = CHEASEdata['B']/CHEASEdata['B0EXP']
       CHEASEdata['Ibs']      = CHEASEdata['Ibs']*mu0/CHEASEdata['R0EXP']/CHEASEdata['B0EXP']
       CHEASEdata['Jbs']      = CHEASEdata['Jbs'] *mu0*CHEASEdata['R0EXP']/CHEASEdata['B0EXP']
       CHEASEdata['Iprl']     = CHEASEdata['Iprl']*mu0/CHEASEdata['R0EXP']/CHEASEdata['B0EXP']
       CHEASEdata['Jprl']     = CHEASEdata['Jprl']*mu0*CHEASEdata['R0EXP']/CHEASEdata['B0EXP']
       CHEASEdata['Istr']     = CHEASEdata['Istr']*mu0/CHEASEdata['R0EXP']/CHEASEdata['B0EXP']
       CHEASEdata['Jphi']     = CHEASEdata['Jphi']*mu0*CHEASEdata['R0EXP']/CHEASEdata['B0EXP']
       CHEASEdata['Itor']     = CHEASEdata['Itor']*mu0/CHEASEdata['R0EXP']/CHEASEdata['B0EXP']
       CHEASEdata['Iohmic']   = CHEASEdata['Iprl']-CHEASEdata['Ibs']
       CHEASEdata['Johmic']   = CHEASEdata['Jprl']-CHEASEdata['Jbs']
       CHEASEdata['pprime']   = CHEASEdata['pprime']*mu0*CHEASEdata['R0EXP']**2/CHEASEdata['B0EXP']
       CHEASEdata['ffprime']  = CHEASEdata['ffprime']/CHEASEdata['B0EXP']
       CHEASEdata['pressure'] = CHEASEdata['pressure']*mu0/CHEASEdata['B0EXP']**2

    return CHEASEdata


def read_chease_hdf(fpath,setParam={},**kwargs):
    if os.path.isfile(fpath) == False:
       errorFunc = traceback.extract_stack(limit=2)[-2][3]
       errorLine = traceback.extract_stack(limit=2)[-2][1]
       errorFile = traceback.extract_stack(limit=2)[-2][2]
       errMSG    = 'Call %s line %5d in file %s Failed.\n'
       errMSG   += 'Fatal: file %s not found.'
       raise IOError(errMSG %(errorFunc,errorLine,errorFile,fpath))

    if 'Normalized' in setParam: Normalized = setParam['Normalized']
    else:                        Normalized = False

    extendCHI  = False
    reverseCHI = False

    hdffh = h5py.File(fpath,'r')
    hdffhKeys = list(hdffh.keys())[0]
    datag = hdffh[hdffhKeys]

    CHEASEdata                = {}

    CHEASEdata['NPSI']        = datag.attrs["NPSI"]
    CHEASEdata['NCHI']        = datag.attrs["NCHI"]
    CHEASEdata['NRBOX']       = datag.attrs["NRBOX"]
    CHEASEdata['NZBOX']       = datag.attrs["NZBOX"]
    CHEASEdata['NBOUND']      = datag.attrs["NBOUND"]
    CHEASEdata['B0EXP']       = datag.attrs["B0EXP"]
    CHEASEdata['R0EXP']       = datag.attrs["R0EXP"]

    CHEASEdata['PSI']         = npy.array(datag["grid"]["PSI"])
    CHEASEdata['rhopsi_SI']   = npy.sqrt(CHEASEdata['PSI']/CHEASEdata['PSI'][-1])
    CHEASEdata['rhotor_SI']   = npy.array(datag["var1d"]["rho_tor"])

    CHEASEdata['CHI']         = npy.array(datag["grid"]["CHI"])
    if CHEASEdata['CHI'][0]>CHEASEdata['CHI'][1]:
       CHEASEdata['CHI']      = CHEASEdata["CHI"][::-1]
       reverseCHI             = True
    if min(abs(CHEASEdata['CHI']-2.0*npy.pi))>=min(npy.diff(CHEASEdata['CHI'])):
       CHEASEdata['CHI']      = npy.append(CHEASEdata['CHI'],2.0*npy.pi)
       extendCHI              = True

    CHEASEdata['Jbs']         = npy.array(datag["var1d"]["jbsBav"])
    CHEASEdata['Zeff']        = npy.array(datag["var1d"]["zeff"])
    CHEASEdata['kappa']       = npy.array(datag["var1d"]["kappa"])
    CHEASEdata['shear']       = npy.array(datag["var1d"]["shear"])
    CHEASEdata['signeo']      = npy.array(datag["var1d"]["signeo"])
    CHEASEdata['pprime']      = 2.0*npy.pi*npy.array(datag["var1d"]["dpdpsi"])

    CHEASEdata['q']           = npy.array(datag["var1d"]["q"])
    CHEASEdata['R_av']        = npy.array(datag["var1d"]["R_av"])
    CHEASEdata['ageom']       = npy.array(datag["var1d"]["ageom"])
    CHEASEdata['Rgeom']       = npy.array(datag["var1d"]["Rgeom"])
    CHEASEdata['Volume']      = npy.array(datag["var1d"]["Volume"])
    CHEASEdata['pressure']    = npy.array(datag["var1d"]["p"])
    CHEASEdata['GDPSI_av']    = npy.array(datag["var1d"]["GDPSI_av"])
    CHEASEdata['radius_av']   = npy.array(datag["var1d"]["radius_av"])

    CHEASEdata['Te']          = npy.array(datag["var1d"]["Te"])
    CHEASEdata['Ti']          = npy.array(datag["var1d"]["Ti"])
    CHEASEdata['ne']          = npy.array(datag["var1d"]["ne"])
    CHEASEdata['ni']          = npy.array(datag["var1d"]["ni"])

    CHEASEdata['Zi']          = 1.0
    CHEASEdata['Zz']          = 6.0
    CHEASEdata['nz']          = CHEASEdata['Zeff']*CHEASEdata['ne']
    CHEASEdata['nz']         -= CHEASEdata['ni']*CHEASEdata['Zi']**2
    CHEASEdata['nz']         /= CHEASEdata['Zz']**2

    CHEASEdata['TePrime']     = npy.array(datag["var1d"]["dTedpsi"])
    CHEASEdata['TiPrime']     = npy.array(datag["var1d"]["dTidpsi"])
    CHEASEdata['nePrime']     = npy.array(datag["var1d"]["dnedpsi"])
    CHEASEdata['niPrime']     = npy.array(datag["var1d"]["dnidpsi"])

    CHEASEdata['f']           = npy.array(datag["var1d"]["f"])
    CHEASEdata['ffprime']     = 2.0*npy.pi*npy.array(datag["var1d"]["fdfdpsi"])
    CHEASEdata['fprime']      = CHEASEdata['ffprime']/CHEASEdata['f']

    CHEASEdata['rmesh']       = npy.array(datag["var1d"]["rmesh"])
    CHEASEdata['zmesh']       = npy.array(datag["var1d"]["zmesh"])
    CHEASEdata['rbound']      = npy.array(datag["var1d"]["rboundplasma"])
    CHEASEdata['zbound']      = npy.array(datag["var1d"]["zboundplasma"])
    CHEASEdata['delta_upper'] = npy.array(datag["var1d"]["delta_upper"])
    CHEASEdata['delta_lower'] = npy.array(datag["var1d"]["delta_lower"])

    CHEASEdata['dpdpsi']      = npy.array(datag["var1d"]["dpdpsi"])
    CHEASEdata['dqdpsi']      = npy.array(datag["var1d"]["dqdpsi"])
    CHEASEdata['dVdpsi']      = npy.array(datag["var1d"]["dVdpsi"])
    CHEASEdata['d2qdpsi2']    = npy.array(datag["var1d"]["d2qdpsi2"])
    CHEASEdata['dsheardpsi']  = npy.array(datag["var1d"]["dsheardpsi"])
    CHEASEdata['dpsidrhotor'] = npy.array(datag["var1d"]["dpsidrhotor"])


   #THE DIMENSION OF ALL THE FOLLOWING QUNATITIES ARE (NCHI,NPSI)
    CHEASEdata['R']           = npy.array(datag["var2d"]["R"])
    CHEASEdata['Z']           = npy.array(datag["var2d"]["Z"])
    CHEASEdata['B']           = npy.array(datag["var2d"]["B"])
    CHEASEdata['J']           = npy.array(datag["var2d"]["Jacobian"])
    CHEASEdata['g11']         = npy.array(datag["var2d"]["g11"])
    CHEASEdata['g22']         = npy.array(datag["var2d"]["g22"])
    CHEASEdata['g33']         = npy.array(datag["var2d"]["g33"])
    CHEASEdata['dBdpsi']      = npy.array(datag["var2d"]["dBdpsi"])
    CHEASEdata['dBdchi']      = npy.array(datag["var2d"]["dBdchi"])
    CHEASEdata['dChidZ']      = npy.array(datag["var2d"]["dChidZ"])
    CHEASEdata['dPsidZ']      = npy.array(datag["var2d"]["dPsidZ"])
    CHEASEdata['dChidR']      = npy.array(datag["var2d"]["dChidR"])
    CHEASEdata['dPsidR']      = npy.array(datag["var2d"]["dPsidR"])

    if extendCHI:
       CHEASEdata['R']      = npy.vstack((CHEASEdata['R'],CHEASEdata['R'][0,:]))
       CHEASEdata['Z']      = npy.vstack((CHEASEdata['Z'],CHEASEdata['Z'][0,:]))
       CHEASEdata['B']      = npy.vstack((CHEASEdata['B'],CHEASEdata['B'][0,:]))
       CHEASEdata['J']      = npy.vstack((CHEASEdata['J'],CHEASEdata['J'][0,:]))
       CHEASEdata['g11']    = npy.vstack((CHEASEdata['g11'],CHEASEdata['g11'][0,:]))
       CHEASEdata['g22']    = npy.vstack((CHEASEdata['g22'],CHEASEdata['g22'][0,:]))
       CHEASEdata['g33']    = npy.vstack((CHEASEdata['g33'],CHEASEdata['g33'][0,:]))
       CHEASEdata['dBdpsi'] = npy.vstack((CHEASEdata['dBdpsi'],CHEASEdata['dBdpsi'][0,:]))
       CHEASEdata['dBdchi'] = npy.vstack((CHEASEdata['dBdchi'],CHEASEdata['dBdchi'][0,:]))
       CHEASEdata['dChidZ'] = npy.vstack((CHEASEdata['dChidZ'],CHEASEdata['dChidZ'][0,:]))
       CHEASEdata['dPsidZ'] = npy.vstack((CHEASEdata['dPsidZ'],CHEASEdata['dPsidZ'][0,:]))
       CHEASEdata['dChidR'] = npy.vstack((CHEASEdata['dChidR'],CHEASEdata['dChidR'][0,:]))
       CHEASEdata['dPsidR'] = npy.vstack((CHEASEdata['dPsidR'],CHEASEdata['dPsidR'][0,:]))

   #THE DIMENSION OF ALL THE FOLLOWING QUNATITIES ARE (NRBOX,NZBOX)
    CHEASEdata['psiRZ']    = npy.array(datag["var2d"]["psiRZ"])
    CHEASEdata['chiRZ']    = npy.array(datag["var2d"]["chiRZ"])

    CHEASEdata['C0']       = npy.trapz(y=CHEASEdata['J']/CHEASEdata['R'],                    x=CHEASEdata['CHI'],axis=0)
    CHEASEdata['C1']       = npy.trapz(y=CHEASEdata['J'],                                    x=CHEASEdata['CHI'],axis=0)
    CHEASEdata['C2']       = npy.trapz(y=CHEASEdata['J']/CHEASEdata['R']**2,                 x=CHEASEdata['CHI'],axis=0)
    CHEASEdata['C3']       = npy.trapz(y=CHEASEdata['J']*CHEASEdata['g11']*CHEASEdata['g33'],x=CHEASEdata['CHI'],axis=0)

    CHEASEdata['y1']       = 1.0+CHEASEdata['C3']/CHEASEdata['C2']/CHEASEdata['f']**2/4.0/npy.pi**2

    CHEASEdata['<B2>']     = npy.trapz(y=CHEASEdata['J']*CHEASEdata['B']**2,x=CHEASEdata['CHI'],axis=0)/CHEASEdata['C1']
    CHEASEdata['<JdotB>']  =-CHEASEdata['f']*CHEASEdata['pprime']-CHEASEdata['fprime']*CHEASEdata['<B2>']/mu0
    CHEASEdata['Jprl']     = CHEASEdata['<JdotB>']/CHEASEdata['B0EXP']

    CHEASEdata['<T/R2>']   = npy.trapz(y=CHEASEdata['J']*CHEASEdata['f']*CHEASEdata['g33'],x=CHEASEdata['CHI'],axis=0)/CHEASEdata['C1']
    CHEASEdata['Iprl']     = CHEASEdata['R0EXP']*CHEASEdata['<JdotB>']/CHEASEdata['<T/R2>']

    CHEASEdata['Istr']     =-((CHEASEdata['C2']/CHEASEdata['C0'])*(CHEASEdata['ffprime']/mu0))
    CHEASEdata['Istr']    +=-((CHEASEdata['C1']/CHEASEdata['C0'])*CHEASEdata['pprime'])
    CHEASEdata['Istr']    *= CHEASEdata['R0EXP']**2

    CHEASEdata['Jphi']     =-(CHEASEdata['R']*CHEASEdata['pprime'])-(CHEASEdata['ffprime']/(mu0*CHEASEdata['R']))

    CHEASEdata['Jtor']     = npy.trapz(y=CHEASEdata['Jphi']*CHEASEdata['J']/CHEASEdata['R'],x=CHEASEdata['CHI'],axis=0)
    CHEASEdata['Itor']     = npy.trapz(y=CHEASEdata['Jtor'],x=CHEASEdata['PSI'],axis=0)


    CHEASEdata['Ibs']      = CHEASEdata['R0EXP']*CHEASEdata['Jbs']/CHEASEdata['<T/R2>']
    CHEASEdata['Iohmic']   = CHEASEdata['Iprl']-CHEASEdata['Ibs']
    CHEASEdata['Johmic']   = CHEASEdata['Jprl']-CHEASEdata['Jbs']

    extendPSI    = npy.linspace(CHEASEdata['PSI'][0],CHEASEdata['PSI'][-1],10*npy.size(CHEASEdata['PSI']))
    extendPHI    = npy.empty_like(extendPSI)
    extendPHI[0] = 0.0
    qfunc        = CubicSpline(CHEASEdata['PSI'],CHEASEdata['q'])
    for i in range(1,npy.size(extendPSI)):
        x           = extendPSI[:i+1]
        y           = qfunc(x)
        extendPHI[i]= npy.trapz(y,x)

    CHEASEdata['PHI'] = npy.empty_like(CHEASEdata['PSI'])
    phifunc           = CubicSpline(extendPSI,extendPHI)
    for i in range(npy.size(CHEASEdata['PSI'])):
        CHEASEdata['PHI'][i] = phifunc(CHEASEdata['PSI'][i])

    CHEASEdata['PHIN']     = (CHEASEdata['PHI']-CHEASEdata['PHI'][0])/(CHEASEdata['PHI'][-1]-CHEASEdata['PHI'][0])
    CHEASEdata['PSIN']     = (CHEASEdata['PSI']-CHEASEdata['PSI'][0])/(CHEASEdata['PSI'][-1]-CHEASEdata['PSI'][0])
    CHEASEdata['CHIN']     = (CHEASEdata['CHI']-CHEASEdata['CHI'][0])/(CHEASEdata['CHI'][-1]-CHEASEdata['CHI'][0])

   #CHEASEdata['rhopsi']   = (CHEASEdata['rhopsi_SI']-CHEASEdata['rhopsi_SI'][0])/(CHEASEdata['rhopsi_SI'][-1]-CHEASEdata['rhopsi_SI'][0])
   #CHEASEdata['rhotor']   = (CHEASEdata['rhotor_SI']-CHEASEdata['rhotor_SI'][0])/(CHEASEdata['rhotor_SI'][-1]-CHEASEdata['rhotor_SI'][0])
    CHEASEdata['rhopsi']   = npy.sqrt(CHEASEdata['PSIN'])
    CHEASEdata['rhotor']   = npy.sqrt(CHEASEdata['PHIN'])

    #Implementing Interpolation to EQDSK Grid
    rhopsiflag = False; rhotorflag = False
    if 'nrhomesh' in setParam:
       if   setParam['nrhomesh'] in [0,'rhopsi']: rhopsiflag = True
       elif setParam['nrhomesh'] in [1,'rhotor']: rhotorflag = True
    else:                                         rhopsiflag = True

    eqdskflag    = False
    expeqflag    = False
    interpflag   = False
    importedflag = False
    for key,value in kwargs.items():
        if   key in ['eqdsk','eqdskdata','eqdskfpath']:
             if    type(value)==str and os.path.isfile(value.strip()):
                   eqdskdata = read_eqdsk(fpath=value.strip())
             else: raise IOError('%s file not found!' % value.strip())
             if 'rhopsi' in eqdskdata: rhopsi = eqdskdata['rhopsi'][:]
             if 'rhotor' in eqdskdata: rhotor = eqdskdata['rhotor'][:]
             if 'PSIN'   in eqdskdata: psi    = eqdskdata['PSIN'][:];  interpflag = True
             if 'PHIN'   in eqdskdata: phi    = eqdskdata['PHIN'][:];  interpflag = True
             eqdskflag = True
        elif key in ['imported','external','other']:
             imported = value.copy()
             if 'rhopsi' in imported: rhopsi = imported['rhopsi'][:]
             if 'rhotor' in imported: rhotor = imported['rhotor'][:]
             if 'PSIN'   in imported: psi = imported['PSIN'][:];       interpflag = True
             else:                    psi = importeddata['rhopsi']**2; interpflag = True
             if 'PHIN'   in imported: phi = imported['PHIN'][:];       interpflag = True
             else:                    phi = importeddata['rhotor']**2; interpflag = True
             importedflag = True

    if interpflag:
       if   rhopsiflag:
            CHEASEdata['q']        = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['q'],psi)
            CHEASEdata['f']        = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['f'],psi)
            CHEASEdata['Te']       = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Te'],psi)
            CHEASEdata['Ti']       = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Ti'],psi)
            CHEASEdata['ne']       = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['ne'],psi)
            CHEASEdata['ni']       = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['ni'],psi)
            CHEASEdata['nz']       = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['nz'],psi)
            CHEASEdata['Ibs']      = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Ibs'],psi)
            CHEASEdata['Jbs']      = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Jbs'],psi)
            CHEASEdata['Zeff']     = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Zeff'],psi)
            CHEASEdata['Istr']     = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Istr'],psi)
            CHEASEdata['Iprl']     = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Iprl'],psi)
            CHEASEdata['Jprl']     = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Jprl'],psi)
            CHEASEdata['kappa']    = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['kappa'],psi)
            CHEASEdata['shear']    = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['shear'],psi)
            CHEASEdata['signeo']   = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['signeo'],psi)
            CHEASEdata['fprime']   = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['fprime'],psi)
            CHEASEdata['Iohmic']   = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Iohmic'],psi)
            CHEASEdata['Johmic']   = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Johmic'],psi)
            CHEASEdata['pprime']   = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['pprime'],psi)
            CHEASEdata['ffprime']  = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['ffprime'],psi)
            CHEASEdata['pressure'] = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['pressure'],psi)

            CHEASEdata['PSIN']     = psi[:]
            CHEASEdata['PHIN']     = phi[:]
            CHEASEdata['rhopsi']   = rhopsi[:]
            CHEASEdata['rhotor']   = rhotor[:]

       elif rhotorflag:
            CHEASEdata['q']        = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['q'],psi,phi,phi)
            CHEASEdata['f']        = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['f'],psi,phi,phi)
            CHEASEdata['Te']       = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Te'],psi,phi,phi)
            CHEASEdata['Ti']       = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Ti'],psi,phi,phi)
            CHEASEdata['ne']       = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['ne'],psi,phi,phi)
            CHEASEdata['ni']       = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['ni'],psi,phi,phi)
            CHEASEdata['nz']       = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['nz'],psi,phi,phi)
            CHEASEdata['Ibs']      = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Ibs'],psi,phi,phi)
            CHEASEdata['Jbs']      = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Jbs'],psi,phi,phi)
            CHEASEdata['Zeff']     = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Zeff'],psi,phi,phi)
            CHEASEdata['Istr']     = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Istr'],psi,phi,phi)
            CHEASEdata['Iprl']     = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Iprl'],psi,phi,phi)
            CHEASEdata['Jprl']     = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Jprl'],psi,phi,phi)
            CHEASEdata['kappa']    = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['kappa'],psi,phi,phi)
            CHEASEdata['shear']    = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['shear'],psi,phi,phi)
            CHEASEdata['signeo']   = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['signeo'],psi,phi,phi)
            CHEASEdata['fprime']   = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['fprime'],psi,phi,phi)
            CHEASEdata['Iohmic']   = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Iohmic'],psi,phi,phi)
            CHEASEdata['Johmic']   = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['Johmic'],psi,phi,phi)
            CHEASEdata['pprime']   = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['pprime'],psi,phi,phi)
            CHEASEdata['ffprime']  = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['ffprime'],psi,phi,phi)
            CHEASEdata['pressure'] = mathtools.interp(CHEASEdata['PSIN'],CHEASEdata['pressure'],psi,phi,phi)

            CHEASEdata['PSIN']     = psi[:]
            CHEASEdata['PHIN']     = phi[:]
            CHEASEdata['rhotor']   = rhopsi[:]
            CHEASEdata['rhopsi']   = rhotor[:]

    if Normalized:
       CHEASEdata['R']        = CHEASEdata['R']/CHEASEdata['R0EXP']
       CHEASEdata['Z']        = CHEASEdata['Z']/CHEASEdata['R0EXP']
       CHEASEdata['B']        = CHEASEdata['B']/CHEASEdata['B0EXP']
       CHEASEdata['Ibs']      = CHEASEdata['Ibs']*mu0/CHEASEdata['R0EXP']/CHEASEdata['B0EXP']
       CHEASEdata['Jbs']      = CHEASEdata['Jbs'] *mu0*CHEASEdata['R0EXP']/CHEASEdata['B0EXP']
       CHEASEdata['Iprl']     = CHEASEdata['Iprl']*mu0/CHEASEdata['R0EXP']/CHEASEdata['B0EXP']
       CHEASEdata['Jprl']     = CHEASEdata['Jprl']*mu0*CHEASEdata['R0EXP']/CHEASEdata['B0EXP']
       CHEASEdata['Istr']     = CHEASEdata['Istr']*mu0/CHEASEdata['R0EXP']/CHEASEdata['B0EXP']
       CHEASEdata['Jphi']     = CHEASEdata['Jphi']*mu0*CHEASEdata['R0EXP']/CHEASEdata['B0EXP']
       CHEASEdata['Itor']     = CHEASEdata['Itor']*mu0/CHEASEdata['R0EXP']/CHEASEdata['B0EXP']
       CHEASEdata['Iohmic']   = CHEASEdata['Iprl']-CHEASEdata['Ibs']
       CHEASEdata['Johmic']   = CHEASEdata['Jprl']-CHEASEdata['Jbs']
       CHEASEdata['pprime']   = CHEASEdata['pprime']*mu0*CHEASEdata['R0EXP']**2/CHEASEdata['B0EXP']
       CHEASEdata['ffprime']  = CHEASEdata['ffprime']/CHEASEdata['B0EXP']
       CHEASEdata['pressure'] = CHEASEdata['pressure']*mu0/CHEASEdata['B0EXP']**2
 
    return CHEASEdata


def read_efit_file(fpath,setParam={}):
   #Developed by Ehab Hassan on 2019-02-27
    if os.path.isfile(fpath) == False:
       errorFunc = traceback.extract_stack(limit=2)[-2][3]
       errorLine = traceback.extract_stack(limit=2)[-2][1]
       errorFile = traceback.extract_stack(limit=2)[-2][2]
       errMSG    = 'Call %s line %5d in file %s Failed.\n'
       errMSG   += 'Fatal: file %s not found.'
       raise IOError(errMSG %(errorFunc,errorLine,errorFile,fpath))

    ofh = open(fpath,'r')
    eqdskdata = {}
    cline = ofh.readline()
    eqdskdata['idum']   = int(cline[48:52])
    eqdskdata['RDIM']   = int(cline[52:56])
    eqdskdata['ZDIM']   = int(cline[56:61])
    cline = ofh.readline()
    eqdskdata['RLEN']   = float(cline[0:16])
    eqdskdata['ZLEN']   = float(cline[16:32])
    eqdskdata['RCTR']   = float(cline[32:48])
    eqdskdata['RLFT']   = float(cline[48:64])
    eqdskdata['ZMID']   = float(cline[64:80])
    cline = ofh.readline()
    eqdskdata['RMAX']   = float(cline[0:16])
    eqdskdata['ZMAX']   = float(cline[16:32])
    eqdskdata['PSIMAX'] = float(cline[32:48])
    eqdskdata['PSIBND'] = float(cline[48:64])
    eqdskdata['BCTR']   = float(cline[64:80])
    cline = ofh.readline()
    eqdskdata['CURNT']  = float(cline[0:16])
    eqdskdata['PSIMAX'] = float(cline[16:32])
    eqdskdata['XDUM']   = float(cline[32:48])
    eqdskdata['RMAX']   = float(cline[48:64])
    eqdskdata['XDUM']   = float(cline[64:80])
    cline = ofh.readline()
    eqdskdata['ZMAX']   = float(cline[0:16])
    eqdskdata['XDUM']   = float(cline[16:32])
    eqdskdata['PSIBND'] = float(cline[32:48])
    eqdskdata['XDUM']   = float(cline[48:64])
    eqdskdata['XDUM']   = float(cline[64:80])

    nlines1D = int(npy.ceil(eqdskdata['RDIM']/5.0))

    eqdskdata['fpol'] = npy.zeros(eqdskdata['RDIM'])
    for iline in range(nlines1D):
        cline = ofh.readline()
        try:
            eqdskdata['fpol'][iline*5+0] = float(cline[0:16])
            eqdskdata['fpol'][iline*5+1] = float(cline[16:32])
            eqdskdata['fpol'][iline*5+2] = float(cline[32:48])
            eqdskdata['fpol'][iline*5+3] = float(cline[48:64])
            eqdskdata['fpol'][iline*5+4] = float(cline[64:80])
        except:
            error = 'empty records'
    eqdskdata['pressure'] = npy.zeros(eqdskdata['RDIM'])
    for iline in range(nlines1D):
        cline = ofh.readline()
        try:
            eqdskdata['pressure'][iline*5+0] = float(cline[0:16])
            eqdskdata['pressure'][iline*5+1] = float(cline[16:32])
            eqdskdata['pressure'][iline*5+2] = float(cline[32:48])
            eqdskdata['pressure'][iline*5+3] = float(cline[48:64])
            eqdskdata['pressure'][iline*5+4] = float(cline[64:80])
        except:
            error = 'empty records'

    eqdskdata['ffprime'] = npy.zeros(eqdskdata['RDIM'])
    for iline in range(nlines1D):
        cline = ofh.readline()
        try:
            eqdskdata['ffprime'][iline*5+0] = float(cline[0:16])
            eqdskdata['ffprime'][iline*5+1] = float(cline[16:32])
            eqdskdata['ffprime'][iline*5+2] = float(cline[32:48])
            eqdskdata['ffprime'][iline*5+3] = float(cline[48:64])
            eqdskdata['ffprime'][iline*5+4] = float(cline[64:80])
        except:
            error = 'empty records'

    eqdskdata['pprime'] = npy.zeros(eqdskdata['RDIM'])
    for iline in range(nlines1D):
        cline = ofh.readline()
        try:
            eqdskdata['pprime'][iline*5+0] = float(cline[0:16])
            eqdskdata['pprime'][iline*5+1] = float(cline[16:32])
            eqdskdata['pprime'][iline*5+2] = float(cline[32:48])
            eqdskdata['pprime'][iline*5+3] = float(cline[48:64])
            eqdskdata['pprime'][iline*5+4] = float(cline[64:80])
        except:
            error = 'empty records'

    nlines2D = int(npy.ceil(eqdskdata['RDIM']*eqdskdata['ZDIM']/5.0))

    eqdskdata['psiRZ'] = npy.zeros(eqdskdata['RDIM']*eqdskdata['ZDIM'])
    for iline in range(nlines2D):
        cline = ofh.readline()
        try:
            eqdskdata['psiRZ'][iline*5+0] = float(cline[0:16])
            eqdskdata['psiRZ'][iline*5+1] = float(cline[16:32])
            eqdskdata['psiRZ'][iline*5+2] = float(cline[32:48])
            eqdskdata['psiRZ'][iline*5+3] = float(cline[48:64])
            eqdskdata['psiRZ'][iline*5+4] = float(cline[64:80])
        except:
            error = 'empty records'
    eqdskdata['psiRZ'] = npy.reshape(eqdskdata['psiRZ'],(eqdskdata['ZDIM'],eqdskdata['RDIM']))
    eqdskdata['qpsi'] = npy.zeros(eqdskdata['RDIM'])
    for iline in range(nlines1D):
        cline = ofh.readline()
        try:
            eqdskdata['qpsi'][iline*5+0] = float(cline[0:16])
            eqdskdata['qpsi'][iline*5+1] = float(cline[16:32])
            eqdskdata['qpsi'][iline*5+2] = float(cline[32:48])
            eqdskdata['qpsi'][iline*5+3] = float(cline[48:64])
            eqdskdata['qpsi'][iline*5+4] = float(cline[64:80])
        except:
            error = 'empty records'

    cline = ofh.readline()
    eqdskdata['nbound'] = int(cline[0:5])
    eqdskdata['nlimit'] = int(cline[5:10])

    if eqdskdata['nbound'] > 0:
       nlines1D = int(npy.ceil(2*eqdskdata['nbound']/5.0))

       Ary1D = npy.zeros(2*eqdskdata['nbound'])
       for iline in range(nlines1D):
           cline = ofh.readline()
           try:
               Ary1D[iline*5+0] = float(cline[0:16])
               Ary1D[iline*5+1] = float(cline[16:32])
               Ary1D[iline*5+2] = float(cline[32:48])
               Ary1D[iline*5+3] = float(cline[48:64])
               Ary1D[iline*5+4] = float(cline[64:80])
           except:
               error = 'empty records'

       eqdskdata['rbound'] = Ary1D[0::2]
       eqdskdata['zbound'] = Ary1D[1::2]


    if eqdskdata['nlimit'] > 0:
       nlines1D = int(npy.ceil(2*eqdskdata['nlimit']/5.0))

       Ary1D = npy.zeros(2*eqdskdata['nlimit'])
       for iline in range(nlines1D):
           cline = ofh.readline()
           try:
               Ary1D[iline*5+0] = float(cline[0:16])
               Ary1D[iline*5+1] = float(cline[16:32])
               Ary1D[iline*5+2] = float(cline[32:48])
               Ary1D[iline*5+3] = float(cline[48:64])
               Ary1D[iline*5+4] = float(cline[64:80])
           except:
               error = 'empty records'

       eqdskdata['rlimit'] = Ary1D[0::2]
       eqdskdata['zlimit'] = Ary1D[1::2]

    eqdskdata['ZR1D']  = npy.arange(eqdskdata['ZDIM'],dtype=float)*eqdskdata['ZLEN']/(eqdskdata['ZDIM']-1.0)
    eqdskdata['ZR1D'] += eqdskdata['ZMID']-eqdskdata['ZMID']/2.0

    eqdskdata['RR1D']  = npy.arange(eqdskdata['RDIM'],dtype=float)*eqdskdata['RLEN']/(eqdskdata['RDIM']-1.0)
    eqdskdata['RR1D'] += eqdskdata['RLFT']

    eqdskdata['psiRZ'] = (eqdskdata['psiRZ']-eqdskdata['PSIMAX'])/(eqdskdata['PSIBND']-eqdskdata['PSIMAX'])

    eqdskdata['PSI']    = (eqdskdata['PSIBND']-eqdskdata['PSIMAX'])*npy.arange(eqdskdata['RDIM'])/(eqdskdata['RDIM']-1.0)
    eqdskdata['PSIN']   = (eqdskdata['PSI']-eqdskdata['PSI'][0])/(eqdskdata['PSI'][-1]-eqdskdata['PSI'][0])
    eqdskdata['rhopsi'] = npy.sqrt(eqdskdata['PSIN'])

    extendPSI    = npy.linspace(eqdskdata['PSI'][0],eqdskdata['PSI'][-1],10*npy.size(eqdskdata['PSI']))
    extendPHI    = npy.empty_like(extendPSI)
    extendPHI[0] = 0.0
    qfunc        = CubicSpline(eqdskdata['PSI'],eqdskdata['qpsi'])
    for i in range(1,npy.size(extendPSI)):
        x           = extendPSI[:i+1]
        y           = qfunc(x)
        extendPHI[i]= npy.trapz(y,x)

    eqdskdata['PHI'] = npy.empty_like(eqdskdata['PSI'])
    phifunc          = CubicSpline(extendPSI,extendPHI)
    for i in range(npy.size(eqdskdata['PSI'])):
        eqdskdata['PHI'][i] = phifunc(eqdskdata['PSI'][i])

    eqdskdata['PHIN']   = (eqdskdata['PHI']-eqdskdata['PHI'][0])/(eqdskdata['PHI'][-1]-eqdskdata['PHI'][0])
    eqdskdata['rhotor'] = npy.sqrt(eqdskdata['PHIN'])

    return eqdskdata


def read_efit(fpath,setParam={},**kwargs):
    EFITdata = read_eqdsk(fpath,setParam,**kwargs)
    return EFITdata


def read_eqdsk(fpath,setParam={},**kwargs):
    if os.path.isfile(fpath) == False:
       errorFunc = traceback.extract_stack(limit=2)[-2][3]
       errorLine = traceback.extract_stack(limit=2)[-2][1]
       errorFile = traceback.extract_stack(limit=2)[-2][2]
       errMSG    = 'Call %s line %5d in file %s Failed.\n'
       errMSG   += 'Fatal: file %s not found.'
       raise IOError(errMSG %(errorFunc,errorLine,errorFile,fpath))
    
    if 'Normalized' in setParam: Normalized = setParam['Normalized']
    else:                        Normalized = False

    rhopsiflag = False; rhotorflag = False
    if 'nrhomesh' in setParam:
        if   setParam['nrhomesh'] in [0,'rhopsi']: rhopsiflag = True
        elif setParam['nrhomesh'] in [1,'rhotor']: rhotorflag = True
    else:                                          rhopsiflag = True

    interpflag   = False
    cheaseflag   = False
    importedflag = False
    for key,value in kwargs.items():
        if   key in ['chease','cheasedata','cheasefpath']:
             if    type(value)==str and os.path.isfile(value.strip()):
                   cheasedata = read_chease(fpath=value.strip())
             else: raise IOError('%s file not found!' % value.strip())
             if 'rhopsi' in cheasedata: rhopsi = cheasedata['rhopsi'][:]
             if 'rhotor' in cheasedata: rhotor = cheasedata['rhotor'][:]
             if 'PSIN'   in cheasedata: psi    = cheasedata['PSIN'][:]; interpflag = True
             if 'PHIN'   in cheasedata: phi    = cheasedata['PHIN'][:]; interpflag = True
             cheaseflag = True
        elif key in ['imported','external','other']:
             imported = value.copy()
             if 'rhopsi' in imported: rhopsi = imported['rhopsi'][:]
             if 'rhotor' in imported: rhotor = imported['rhotor'][:]
             if 'PSIN'   in imported: psi    = imported['PSIN'];      interpflag = True
             else:                    psi    = imported['rhopsi']**2; interpflag = True
             if 'PHIN'   in imported: phi    = imported['PHIN'];      interpflag = True
             else:                    phi    = imported['rhotor']**2; interpflag = True
             importedflag = True

    EQDSKdata = read_efit_file(fpath)

    if interpflag:
       if   rhopsiflag:
            EQDSKdata['q']        = mathtools.interp(EQDSKdata['PSIN'],EQDSKdata['qpsi'],psi) 
            EQDSKdata['f']        = mathtools.interp(EQDSKdata['PSIN'],EQDSKdata['fpol'],psi)
            EQDSKdata['pprime']   = mathtools.interp(EQDSKdata['PSIN'],EQDSKdata['pprime'],psi)
            EQDSKdata['ffprime']  = mathtools.interp(EQDSKdata['PSIN'],EQDSKdata['ffprime'],psi)
            EQDSKdata['pressure'] = mathtools.interp(EQDSKdata['PSIN'],EQDSKdata['pressure'],psi) 
            EQDSKdata['rhopsi']   = rhopsi[:]
            EQDSKdata['rhotor']   = rhotor[:]
       elif rhotorflag:
            EQDSKdata['q']        = mathtools.interp(EQDSKdata['PSIN'],EQDSKdata['qpsi'],psi,phi,phi) 
            EQDSKdata['f']        = mathtools.interp(EQDSKdata['PSIN'],EQDSKdata['fpol'],psi,phi,phi)
            EQDSKdata['pprime']   = mathtools.interp(EQDSKdata['PSIN'],EQDSKdata['pprime'],psi,phi,phi)
            EQDSKdata['ffprime']  = mathtools.interp(EQDSKdata['PSIN'],EQDSKdata['ffprime'],psi,phi,phi)
            EQDSKdata['pressure'] = mathtools.interp(EQDSKdata['PSIN'],EQDSKdata['pressure'],psi,phi,phi)
            EQDSKdata['rhopsi']   = rhopsi[:]
            EQDSKdata['rhotor']   = rhotor[:]
    else:
            EQDSKdata['q']        = EQDSKdata['qpsi'][:]
            EQDSKdata['f']        = EQDSKdata['fpol'][:]

    if Normalized:
       EQDSKdata['RCTR']          = abs(EQDSKdata['RCTR'])
       EQDSKdata['BCTR']          = abs(EQDSKdata['BCTR'])
       EQDSKdata['rbound']        = EQDSKdata['rbound']/EQDSKdata['RCTR']
       EQDSKdata['zbound']        = EQDSKdata['zbound']/EQDSKdata['RCTR']
       if 'rlimit' in EQDSKdata:
          EQDSKdata['rlimit']        = EQDSKdata['rlimit']/EQDSKdata['RCTR']
       if 'zlimit' in EQDSKdata:
          EQDSKdata['zlimit']        = EQDSKdata['zlimit']/EQDSKdata['RCTR']
       EQDSKdata['pressure']      = EQDSKdata['pressure']*mu0/EQDSKdata['BCTR']**2
       EQDSKdata['pprime']        = EQDSKdata['pprime']*mu0*EQDSKdata['RCTR']**2/EQDSKdata['BCTR']
       EQDSKdata['ffprime']       = EQDSKdata['ffprime']/EQDSKdata['BCTR']

    return EQDSKdata

def read_expeq(fpath,setParam={},**kwargs):
    if os.path.isfile(fpath) == False:
       errorFunc = traceback.extract_stack(limit=2)[-2][3]
       errorLine = traceback.extract_stack(limit=2)[-2][1]
       errorFile = traceback.extract_stack(limit=2)[-2][2]
       errMSG    = 'Call %s line %5d in file %s Failed.\n'
       errMSG   += 'Fatal: file %s not found.'
       raise IOError(errMSG %(errorFunc,errorLine,errorFile,fpath))


    rhopsiflag = False; rhotorflag = False
    if 'nrhomesh' in setParam:
        if   setParam['nrhomesh'] in [0,'rhopsi']: rhopsiflag = True
        elif setParam['nrhomesh'] in [1,'rhotor']: rhotorflag = True
    else:                                          rhopsiflag = True

    eqdskflag    = False
    cheaseflag   = False
    interpflag   = False
    importedflag = False
    for key,value in kwargs.items():
        if   key in ['chease','cheasedata','cheasefpath']:
             if    type(value)==str and os.path.isfile(value.strip()):
                   cheasedata = read_chease(fpath=value.strip())
             else: raise IOError('%s file not found!' % value.strip())
             if 'rhopsi' in cheasedata: rhopsi = cheasedata['rhopsi'][:]; interpflag = True
             if 'rhotor' in cheasedata: rhotor = cheasedata['rhotor'][:]; interpflag = True
             cheaseflag = True
        elif key in ['eqdsk','eqdskdata','eqdskfpath']:
             if    type(value)==str and os.path.isfile(value.strip()):
                   eqdskdata = read_eqdsk(fpath=value.strip())
             else: raise IOError('%s file not found!' % value.strip())
             if 'rhopsi' in eqdskdata: rhopsi = eqdskdata['rhopsi'][:]; interpflag = True
             if 'rhotor' in eqdskdata: rhotor = eqdskdata['rhotor'][:]; interpflag = True
             eqdskflag = True
        elif key in ['imported','external','other']:
             imported = value.copy()
             if 'rhopsi' in imported: rhopsi = imported['rhopsi'][:]; interpflag = True
             if 'rhotor' in imported: rhotor = imported['rhotor'][:]; interpflag = True
             importedflag = True

    ofh = open(fpath,'r')
    EXPEQOUT = ofh.readlines()
    ofh.close()

    EXPEQdata                  = {} 
    EXPEQdata['aspect']        = float(EXPEQOUT[0])
    EXPEQdata['zgeom']         = float(EXPEQOUT[1])
    EXPEQdata['pedge']         = float(EXPEQOUT[2])
    nRZmesh                    =   int(EXPEQOUT[3])
    EXPEQdata['nRZmesh']       = nRZmesh
    EXPEQdata['rbound']        = npy.array([irec.split()[0] for irec in EXPEQOUT[4:nRZmesh+4]],dtype=float)
    EXPEQdata['zbound']        = npy.array([irec.split()[1] for irec in EXPEQOUT[4:nRZmesh+4]],dtype=float)
    
    nrhomesh                   = int(EXPEQOUT[nRZmesh+4].split()[0])
    EXPEQdata['nrhomesh']      = nrhomesh
    EXPEQdata['nppfun']        = int(EXPEQOUT[nRZmesh+4].split()[1])
    EXPEQdata['nsttp']         = int(EXPEQOUT[nRZmesh+5].split()[0])
    EXPEQdata['nrhotype']      = int(EXPEQOUT[nRZmesh+5].split()[1])

    if   EXPEQdata['nrhotype']  == 0:
         EXPEQdata['rhopsi']    = npy.array(EXPEQOUT[nRZmesh+6+0*nrhomesh:nRZmesh+6+1*nrhomesh],dtype=float)
    elif EXPEQdata['nrhotype']  == 1:
         EXPEQdata['rhotor']    = npy.array(EXPEQOUT[nRZmesh+6+0*nrhomesh:nRZmesh+6+1*nrhomesh],dtype=float)
    
    if   EXPEQdata['nppfun']    == 4:
         EXPEQdata['pprime']    = npy.array(EXPEQOUT[nRZmesh+6+1*nrhomesh:nRZmesh+6+2*nrhomesh],dtype=float)
    elif EXPEQdata['nppfun']    == 8:
         EXPEQdata['pressure']  = npy.array(EXPEQOUT[nRZmesh+6+1*nrhomesh:nRZmesh+6+2*nrhomesh],dtype=float)
  
    if   EXPEQdata['nsttp']     == 1:
         EXPEQdata['ffprime']   = npy.array(EXPEQOUT[nRZmesh+6+2*nrhomesh:nRZmesh+6+3*nrhomesh],dtype=float)
    elif EXPEQdata['nsttp']     == 2:
         EXPEQdata['Istr']      = npy.array(EXPEQOUT[nRZmesh+6+2*nrhomesh:nRZmesh+6+3*nrhomesh],dtype=float)
    elif EXPEQdata['nsttp']     == 3:
         EXPEQdata['Iprl']      = npy.array(EXPEQOUT[nRZmesh+6+2*nrhomesh:nRZmesh+6+3*nrhomesh],dtype=float)
    elif EXPEQdata['nsttp']     == 4:
         EXPEQdata['Jprl']      = npy.array(EXPEQOUT[nRZmesh+6+2*nrhomesh:nRZmesh+6+3*nrhomesh],dtype=float)
    elif EXPEQdata['nsttp']     == 5:
         EXPEQdata['q']         = npy.array(EXPEQOUT[nRZmesh+6+2*nrhomesh:nRZmesh+6+3*nrhomesh],dtype=float)

    if interpflag:
       if   rhopsiflag and EXPEQdata['nrhotype']==0:
            if   EXPEQdata['nppfun']    == 4:
                 EXPEQdata['pprime']    = mathtools.interp(EXPEQdata['rhopsi'],EXPEQdata['pprime'],rhopsi)
            elif EXPEQdata['nppfun']    == 8:
                 EXPEQdata['pressure']  = mathtools.interp(EXPEQdata['rhopsi'],EXPEQdata['pressure'],rhopsi)
            if   EXPEQdata['nsttp']     == 1:
                 EXPEQdata['ffprime']   = mathtools.interp(EXPEQdata['rhopsi'],EXPEQdata['ffprime'],rhopsi)
            elif EXPEQdata['nsttp']     == 2:
                 EXPEQdata['Istr']      = mathtools.interp(EXPEQdata['rhopsi'],EXPEQdata['Istr'],rhopsi)
            elif EXPEQdata['nsttp']     == 3:
                 EXPEQdata['Iprl']      = mathtools.interp(EXPEQdata['rhopsi'],EXPEQdata['Iprl'],rhopsi)
            elif EXPEQdata['nsttp']     == 4:
                 EXPEQdata['Jprl']      = mathtools.interp(EXPEQdata['rhopsi'],EXPEQdata['Jprl'],rhopsi)
            elif EXPEQdata['nsttp']     == 5:
                 EXPEQdata['q']         = mathtools.interp(EXPEQdata['rhopsi'],EXPEQdata['q'],rhopsi)
            EXPEQdata['rhopsi']         = rhopsi[:]
            EXPEQdata['rhotor']         = rhotor[:]
       elif rhotorflag and EXPEQdata['nrhotype']==1:
            if   EXPEQdata['nppfun']    == 4:
                 EXPEQdata['pprime']    = mathtools.interp(EXPEQdata['rhotor'],EXPEQdata['pprime'],rhotor)
            elif EXPEQdata['nppfun']    == 8:
                 EXPEQdata['pressure']  = mathtools.interp(EXPEQdata['rhotor'],EXPEQdata['pressure'],rhotor)
            if   EXPEQdata['nsttp']     == 1:
                 EXPEQdata['ffprime']   = mathtools.interp(EXPEQdata['rhotor'],EXPEQdata['ffprime'],rhotor)
            elif EXPEQdata['nsttp']     == 2:
                 EXPEQdata['Istr']      = mathtools.interp(EXPEQdata['rhotor'],EXPEQdata['Istr'],rhotor)
            elif EXPEQdata['nsttp']     == 3:
                 EXPEQdata['Iprl']      = mathtools.interp(EXPEQdata['rhotor'],EXPEQdata['Iprl'],rhotor)
            elif EXPEQdata['nsttp']     == 4:
                 EXPEQdata['Jprl']      = mathtools.interp(EXPEQdata['rhotor'],EXPEQdata['Jprl'],rhotor)
            elif EXPEQdata['nsttp']     == 5:
                 EXPEQdata['q']         = mathtools.interp(EXPEQdata['rhotor'],EXPEQdata['q'],rhotor)
            EXPEQdata['rhopsi']         = rhopsi[:]
            EXPEQdata['rhotor']         = rhotor[:]
       elif rhopsiflag and EXPEQdata['nrhotype']==1:
            if   EXPEQdata['nppfun']    == 4:
                 EXPEQdata['pprime']    = mathtools.interp(EXPEQdata['rhotor'],EXPEQdata['pprime'],rhotor,rhopsi,rhopsi)
            elif EXPEQdata['nppfun']    == 8:
                 EXPEQdata['pressure']  = mathtools.interp(EXPEQdata['rhotor'],EXPEQdata['pressure'],rhotor,rhopsi,rhopsi)
            if   EXPEQdata['nsttp']     == 1:
                 EXPEQdata['ffprime']   = mathtools.interp(EXPEQdata['rhotor'],EXPEQdata['ffprime'],rhotor,rhopsi,rhopsi)
            elif EXPEQdata['nsttp']     == 2:
                 EXPEQdata['Istr']      = mathtools.interp(EXPEQdata['rhotor'],EXPEQdata['Istr'],rhotor,rhopsi,rhopsi)
            elif EXPEQdata['nsttp']     == 3:
                 EXPEQdata['Iprl']      = mathtools.interp(EXPEQdata['rhotor'],EXPEQdata['Iprl'],rhotor,rhopsi,rhopsi)
            elif EXPEQdata['nsttp']     == 4:
                 EXPEQdata['Jprl']      = mathtools.interp(EXPEQdata['rhotor'],EXPEQdata['Jprl'],rhotor,rhopsi,rhopsi)
            elif EXPEQdata['nsttp']     == 5:
                 EXPEQdata['q']         = mathtools.interp(EXPEQdata['rhotor'],EXPEQdata['q'],rhotor,rhopsi,rhopsi)
            EXPEQdata['rhopsi']         = rhopsi[:]
            EXPEQdata['rhotor']         = rhotor[:]
       elif rhotorflag and EXPEQdata['nrhotype']==0:
            if   EXPEQdata['nppfun']    == 4:
                 EXPEQdata['pprime']    = mathtools.interp(EXPEQdata['rhopsi'],EXPEQdata['pprime'],rhopsi,rhotor,rhotor)
            elif EXPEQdata['nppfun']    == 8:
                 EXPEQdata['pressure']  = mathtools.interp(EXPEQdata['rhopsi'],EXPEQdata['pressure'],rhopsi,rhotor,rhotor)
            if   EXPEQdata['nsttp']     == 1:
                 EXPEQdata['ffprime']   = mathtools.interp(EXPEQdata['rhopsi'],EXPEQdata['ffprime'],rhopsi,rhotor,rhotor)
            elif EXPEQdata['nsttp']     == 2:
                 EXPEQdata['Istr']      = mathtools.interp(EXPEQdata['rhopsi'],EXPEQdata['Istr'],rhopsi,rhotor,rhotor)
            elif EXPEQdata['nsttp']     == 3:
                 EXPEQdata['Iprl']      = mathtools.interp(EXPEQdata['rhopsi'],EXPEQdata['Iprl'],rhopsi,rhotor,rhotor)
            elif EXPEQdata['nsttp']     == 4:
                 EXPEQdata['Jprl']      = mathtools.interp(EXPEQdata['rhopsi'],EXPEQdata['Jprl'],rhopsi,rhotor,rhotor)
            elif EXPEQdata['nsttp']     == 5:
                 EXPEQdata['q']         = mathtools.interp(EXPEQdata['rhopsi'],EXPEQdata['q'],rhopsi,rhotor,rhotor)
            EXPEQdata['rhopsi']         = rhopsi[:]
            EXPEQdata['rhotor']         = rhotor[:]

    return EXPEQdata


def write_expeq(setParam={},**kwargs):
    '''
    nrhomesh=[rho_type(0:rhopsi,1:rhotor),rho_src(0:chease,1:eqdsk)]
    nppfun  =[pressure_type(8:pressure,4:pprime),pressure_src(0:chease,1:eqdsk,2:expeq,3:exptnz,4:iterdb,5:profiles)]
    nsttp   =[current_type(1:ffprime,2:istar,3:iprl,4:jprl,5:q),current_src(0:chease,1:eqdsk,2:expeq)]
    boundary= boundary_type(0:asis,1,interpolate)
    '''

    eqdskflag    = False
    expeqflag    = False
    interpflag   = False
    cheaseflag   = False
    iterdbflag   = False
    exptnzflag   = False
    importedflag = False
    profilesflag = False
    for key,value in kwargs.items():
        if key in ['chease','cheasedata','cheasefpath']:
           if os.path.isfile(value.strip()):
              cheasepath = value.strip()
              cheasedata = read_chease(fpath=cheasepath)
              cheaseflag = True

        if key in ['eqdsk','eqdskdata','eqdskfpath']:
           if os.path.isfile(value.strip()):
              eqdskpath = value.strip()
              eqdskdata = read_eqdsk(fpath=eqdskpath)
              eqdskflag = True

        if key in ['expeq','expeqdata','expeqfpath']:
           if os.path.isfile(value.strip()):
              expeqpath = value.strip()
              expeqdata = read_expeq(fpath=expeqpath)
              expeqflag = True

        if key in ['exptnz','exptnzdata','exptnzfpath']:
           if os.path.isfile(value.strip()):
              exptnzpath = value.strip()
              exptnzflag = True

        if key in ['profiles','profilesdata','profilesfpath']:
           if os.path.isfile(value.strip()):
              profilespath = value.strip()
              profilesflag = True

        if key in ['iterdb','iterdbdata','iterdbfpath']:
           if os.path.isfile(value.strip()):
              iterdbpath = value.strip()
              iterdbflag = True

        if key in ['imported','external','others']:
              imported = value.copy()
              importedflag = True

    if not (cheaseflag or expeqflag  or eqdskflag or exptnzflag or profilesflag or iterdbflag or importedflag):
       raise IOError('FATAL: NO VALID INPUT PROFILES AVAILABLE. EXIT!')

    if   PYTHON3:
         setParamKeys = list(setParam.keys())
    elif PYTHON2:
         setParamKeys = setParam.keys()

    if 'outfile' in setParam:
        outfile = setParam['outfile']
    else:
        outfile = True

    if 'nrhomesh' in setParam:
       if   type(setParam['nrhomesh'])==list:
            if   type(setParam['nrhomesh'][0])==float: setParam['nrhomesh'][0] = int(setParam['nrhomesh'][0])
            elif type(setParam['nrhomesh'][0])==str:   setParam['nrhomesh'][0] = setParam['nrhomesh'][0].lower()
            if   type(setParam['nrhomesh'][1])==float: setParam['nrhomesh'][1] = int(setParam['nrhomesh'][1])
            elif type(setParam['nrhomesh'][1])==str:   setParam['nrhomesh'][1] = setParam['nrhomesh'][1].lower()
            elif      setParam['nrhomesh'][1] ==None:  setParam['nrhomesh'][1] = None
            nrhotype = setParam['nrhomesh'][:]
       elif type(setParam['nrhomesh']) in [int,str,float]:
            if   type(setParam['nrhomesh'])==float: setParam['nrhomesh'] = int(setParam['nrhomesh'])
            elif type(setParam['nrhomesh'])==str:   setParam['nrhomesh'] = setParam['nrhomesh'].lower()
            if   cheaseflag:   nrhotype=[setParam['nrhomesh'],0]
            elif eqdskflag:    nrhotype=[setParam['nrhomesh'],1]
            elif importedflag: nrhotype=[setParam['nrhomesh'],7]
    else:
            if   cheaseflag:   nrhotype=[0,0]
            elif eqdskflag:    nrhotype=[0,1]
            elif importedflag: nrhotype=[0,7]

    if 'nppfun' in setParam:
       if   type(setParam['nppfun'])==list:
            if   type(setParam['nppfun'][0])==float: setParam['nppfun'][0] = int(setParam['nppfun'][0])
            elif type(setParam['nppfun'][0])==str:   setParam['nppfun'][0] = setParam['nppfun'][0].lower()
            if   type(setParam['nppfun'][1])==float: setParam['nppfun'][1] = int(setParam['nppfun'][1])
            elif type(setParam['nppfun'][1])==str:   setParam['nppfun'][1] = setParam['nppfun'][1].lower()
            nppfun=setParam['nppfun'][:]
       elif type(setParam['nppfun']) in [int,float,str]:
            if   type(setParam['nppfun'])==float: setParam['nppfun'] = int(setParam['nppfun'])
            elif type(setParam['nppfun'])==str:   setParam['nppfun'] = setParam['nppfun'].lower()
            if   cheaseflag:   nppfun=[setParam['nppfun'],0]
            elif eqdskflag:    nppfun=[setParam['nppfun'],1]
            elif expeqflag:    nppfun=[setParam['nppfun'],2]
            elif exptnzflag:   nppfun=[setParam['nppfun'],3]
            elif profilesflag: nppfun=[setParam['nppfun'],4]
            elif iterdbflag:   nppfun=[setParam['nppfun'],5]
            elif importedflag: nppfun=[setParam['nppfun'],7]
    else:
            if   cheaseflag:   nppfun=[4,0]
            elif eqdskflag:    nppfun=[4,1]
            elif expeqflag:    nppfun=[4,2]
            elif exptnzflag:   nppfun=[4,3]
            elif profilesflag: nppfun=[4,4]
            elif iterdbflag:   nppfun=[4,5]
            elif importedflag: nppfun=[4,7]

    if 'nsttp' in setParam:
       if   type(setParam['nsttp'])==list:
            if   type(setParam['nsttp'][0])==float: setParam['nsttp'][0] = int(setParam['nsttp'][0])
            elif type(setParam['nsttp'][0])==str:   setParam['nsttp'][0] = setParam['nsttp'][0].lower()
            if   type(setParam['nsttp'][1])==float: setParam['nsttp'][1] = int(setParam['nsttp'][1])
            elif type(setParam['nsttp'][1])==str:   setParam['nsttp'][1] = setParam['nsttp'][1].lower()
            nsttp=setParam['nsttp'][:]
       elif type(setParam['nsttp']) in [int,float,str]:
            if   type(setParam['nsttp'])==float: setParam['nsttp'] = int(setParam['nsttp'])
            elif type(setParam['nsttp'])==str:   setParam['nsttp'] = setParam['nsttp'].lower()
            if   cheaseflag:   nsttp=[setParam['nsttp'],0]
            elif eqdskflag:    nsttp=[setParam['nsttp'],1]
            elif expeqflag:    nsttp=[setParam['nsttp'],2]
            elif importedflag: nsttp=[setParam['nsttp'],7]
    else:
            if   cheaseflag:   nsttp=[1,0]
            elif eqdskflag:    nsttp=[1,1]
            elif expeqflag:    nsttp=[1,2]
            elif importedflag: nsttp=[1,7]

    if 'geometry' in setParam:
       if type(setParam['geometry']) in [int,float,str]:
            if   type(setParam['geometry'])==float: setParam['geometry'] = int(setParam['geometry'])
            elif type(setParam['geometry'])==str:   setParam['geometry'] = setParam['geometry'].lower()
            if   cheaseflag:   geometry=setParam['geometry']
            elif eqdskflag:    geometry=setParam['geometry']
            elif expeqflag:    geometry=setParam['geometry']
            elif importedflag: geometry=setParam['geometry']
    else:
            if   cheaseflag:   geometry=0
            elif eqdskflag:    geometry=1
            elif expeqflag:    geometry=2
            elif importedflag: geometry=7


    if 'boundary' in setParam:
         if   type(setParam['boundary'])==float:
              boundary = int(setParam['boundary'])
         elif type(setParam['boundary'])==str:
              boundary = setParam['boundary'].lower()
    else:
              boundary = 0

    if   nrhotype[1] in [0,'chease'] and not cheaseflag:
         raise IOError('chease.h5 FILE IS NOT PROVIDED. EXIT!')
    elif nrhotype[1] in [1,'eqdsk'] and not eqdskflag:
         raise IOError('EQDSK FILE IS NOT PROVIDED. EXIT!')
    elif nrhotype[1] in [7,'imported'] and not eqdskflag:
         raise IOError('IMPORTED DATA IS NOT PROVIDED. EXIT!')

    if nrhotype[0] in [1,'rhotor'] and nsttp[0] in [1,'ffprime','ffprimen']:
       raise NameError('FATAL: nrhotype (must) = 0 or rhopsi. Exit!')

    if   nsttp[1] in [0,'chease'] and not cheaseflag:
         raise IOError('chease.h5 FILE IS NOT PROVIDED. EXIT!')
    elif nsttp[1] in [1,'eqdsk'] and not eqdskflag:
         raise IOError('EQSSK FILE IS NOT PROVIDED. EXIT!')
    elif nsttp[1] in [2,'expeq'] and not expeqflag:
         raise IOError('EXPEQ FILE IS NOT PROVIDED. EXIT!')
    elif nsttp[1] in [7,'imported'] and not importedflag:
         raise IOError('IMPORTED DATA IS NOT PROVIDED. EXIT!')

    if   nppfun[1] in [0,'chease'] and not cheaseflag:
         raise IOError('chease.h5 FILE IS NOT PROVIDED. EXIT!')
    elif nppfun[1] in [1,'eqdsk'] and not eqdskflag:
         raise IOError('EQDSK FILE IS NOT PROVIDED. EXIT!')
    elif nppfun[1] in [2,'expeq'] and not expeqflag:
         raise IOError('EXPEQ FILE IS NOT PROVIDED. EXIT!')
    elif nppfun[1] in [3,'exptnz'] and not exptnzflag:
         raise IOError('EXPTNZ FILE IS NOT PROVIDED. EXIT!')
    elif nppfun[1] in [4,'profiles'] and not profilesflag:
         raise IOError('PROFILES FILE IS NOT PROVIDED. EXIT!')
    elif nppfun[1] in [5,'iterdb'] and not iterdbflag:
         raise IOError('ITERDB FILE IS NOT PROVIDED. EXIT!')
    elif nppfun[1] in [7,'imported'] and not importedflag:
         raise IOError('IMPORTED DATA IS NOT PROVIDED. EXIT!')

    currents = []
    currents.append([2,'istr','istrn','istar','istarn'])
    currents.append([3,'iprl','iprln','iparallel','iparalleln'])
    currents.append([4,'jprl','jprln','jparallel','jparalleln'])
    if (nsttp[1] in [1,'eqdsk']) and (nsttp[0] in currents):
       raise IOError("FATAL: eqdsk option is not accepted with nsttp = 2, 3 or 4")

    expeq = {}

    if   geometry in [7,'imported'] and 'imported' in locals():
         if   'R0EXP'  in imported: expeq['R0EXP']  = imported['R0EXP']
         elif 'R0EXP'  in setParam: expeq['R0EXP']  = setParam['R0EXP']
         else:                      expeq['R0EXP']  = 1.0
         if   'B0EXP'  in imported: expeq['B0EXP']  = imported['B0EXP']
         elif 'B0EXP'  in setParam: expeq['B0EXP']  = setParam['B0EXP']
         else:                      expeq['B0EXP']  = 1.0
         if   'ZMAX'   in imported: expeq['zgeom']  = imported['ZMAX']
         elif 'ZMAX'   in setParam: expeq['zgeom']  = setParam['ZMAX']
         else:                      expeq['zgeom']  = 0.0
         if   'aspect' in imported: expeq['aspect'] = imported['aspect']
         elif 'aspect' in setParam: expeq['aspect'] = setParam['aspect']
         else:
              expeq['aspect']  = (max(imported['rbound'])-min(imported['rbound']))
              expeq['aspect'] /= (max(imported['rbound'])+min(imported['rbound']))

         if   'rbound' in imported and 'zbound' in imported:
              expeq['nRZmesh'] = npy.size(imported['rbound'])
              expeq['rbound']  = imported['rbound'][:]
              expeq['zbound']  = imported['zbound'][:]
         elif 'eqdskdata' in locals():
              eqdskParam = {'boundary_type':boundary}
              rbound,zbound    = find_boundary(eqdskdata,setParam=eqdskParam)
              expeq['nRZmesh'] = npy.size(rbound)
              expeq['rbound']  = rbound/expeq['R0EXP']
              expeq['zbound']  = zbound/expeq['R0EXP']

    elif geometry in [0,'chease'] and 'cheasedata' in locals():
         if   'R0EXP'  in imported: expeq['R0EXP']  = imported['R0EXP']
         elif 'R0EXP'  in setParam: expeq['R0EXP']  = setParam['R0EXP']
         else:                      expeq['R0EXP']  = cheasedata['R0EXP']
         if   'B0EXP'  in imported: expeq['B0EXP']  = imported['B0EXP']
         elif 'B0EXP'  in setParam: expeq['B0EXP']  = setParam['B0EXP']
         else:                      expeq['B0EXP']  = cheasedata['B0EXP']
         if   'ZMAX'   in imported: expeq['zgeom']  = imported['ZMAX']
         elif 'ZMAX'   in setParam: expeq['zgeom']  = setParam['ZMAX']
         else:                      expeq['zgeom']  = npy.mean(cheasedata['zmesh'])/expeq['R0EXP']
         if   'aspect' in imported: expeq['aspect'] = imported['aspect']
         elif 'aspect' in setParam: expeq['aspect'] = setParam['aspect']
         else:
              expeq['aspect']  = (max(cheasedata['rbound'])-min(cheasedata['rbound']))
              expeq['aspect'] /= (max(cheasedata['rbound'])+min(cheasedata['rbound']))

         if 'rbound' in imported and 'zbound' in imported:
            expeq['nRZmesh'] = npy.size(imported['rbound'])
            expeq['rbound']  = imported['rbound'][:]
            expeq['zbound']  = imported['zbound'][:]
         else:
            expeq['nRZmesh'] = npy.size(cheasedata['rbound'])
            expeq['rbound']  = cheasedata['rbound']/expeq['R0EXP']
            expeq['zbound']  = cheasedata['zbound']/expeq['R0EXP']

    elif geometry in [2,'expeq'] and 'expeqdata' in locals():
         if   'R0EXP'  in imported: expeq['R0EXP']  = imported['R0EXP']
         elif 'R0EXP'  in setParam: expeq['R0EXP']  = setParam['R0EXP']
         else:                      expeq['R0EXP']  = 1.0
         if   'B0EXP'  in imported: expeq['B0EXP']  = imported['B0EXP']
         elif 'B0EXP'  in setParam: expeq['B0EXP']  = setParam['B0EXP']
         else:                      expeq['B0EXP']  = 1.0
         if   'ZMAX'   in imported: expeq['zgeom']  = imported['ZMAX']
         elif 'ZMAX'   in setParam: expeq['zgeom']  = setParam['ZMAX']
         else:                      expeq['zgeom']  = expeqdata['aspect']
         if   'aspect' in imported: expeq['aspect'] = imported['aspect']
         elif 'aspect' in setParam: expeq['aspect'] = setParam['aspect']
         else:                                        expeqdata['zgeom']

         if 'rbound' in imported and 'zbound' in imported:
            expeq['nRZmesh'] = npy.size(imported['rbound'])
            expeq['rbound']  = imported['rbound'][:]
            expeq['zbound']  = imported['zbound'][:]
         else:
            expeq['nRZmesh'] = npy.size(expeqdata['rbound'])
            expeq['rbound']  = expeqdata['rbound'][:]
            expeq['zbound']  = expeqdata['zbound'][:]

    elif geometry in [1,'eqdsk'] or 'eqdskdata' in locals():
         if   'R0EXP'  in imported: expeq['R0EXP']  = imported['R0EXP']
         elif 'R0EXP'  in setParam: expeq['R0EXP']  = setParam['R0EXP']
         else:                      expeq['R0EXP']  = abs(eqdskdata['RCTR'])
         if   'B0EXP'  in imported: expeq['B0EXP']  = imported['B0EXP']
         elif 'B0EXP'  in setParam: expeq['B0EXP']  = setParam['B0EXP']
         else:                      expeq['B0EXP']  = abs(eqdskdata['BCTR'])
         if   'ZMAX'   in imported: expeq['zgeom']  = imported['ZMAX']
         elif 'ZMAX'   in setParam: expeq['zgeom']  = setParam['ZMAX']
         else:                      expeq['zgeom']  = eqdskdata['ZMAX']/expeq['R0EXP']
         if   'aspect' in imported: expeq['aspect'] = imported['aspect']
         elif 'aspect' in setParam: expeq['aspect'] = setParam['aspect']
         else:
              expeq['aspect']  = (max(rbound)-min(rbound))
              expeq['aspect'] /= (max(rbound)+min(rbound))

         if 'rbound' in imported and 'zbound' in imported:
            expeq['nRZmesh'] = npy.size(imported['rbound'])
            expeq['rbound']  = imported['rbound'][:]
            expeq['zbound']  = imported['zbound'][:]
         else:
            eqdskParam = {'boundary_type':boundary}
            rbound,zbound    = find_boundary(eqdskdata,setParam=eqdskParam)
            expeq['nRZmesh'] = npy.size(rbound)
            expeq['rbound']  = rbound/expeq['R0EXP']
            expeq['zbound']  = zbound/expeq['R0EXP']

    if   nrhotype[0] in [0,'rhopsi','rhopsin']:
         rhopsiflag = True
         rhotorflag = False
    elif nrhotype[0] in [1,'rhotor','rhotorn']:
         rhopsiflag = False
         rhotorflag = True

    SetParam = {'nrhomesh':nrhotype[0]}

    if   nrhotype[1] in [0,'chease']:
         if cheaseflag:
            cheasedata = read_chease(fpath=cheasepath,setParam=SetParam)
         if eqdskflag:
            eqdskdata = read_eqdsk(fpath=eqdskpath,setParam=SetParam,chease=cheasepath)
         if expeqflag:
            expeqdata = read_expeq(fpath=expeqpath,setParam=SetParam,chease=cheasepath)
         if exptnzflag:
            exptnzdata = read_exptnz(fpath=exptnzpath,setParam=SetParam,chease=cheasepath)
         if iterdbflag:
            iterdbdata = read_iterdb(fpath=iterdbpath,setParam=SetParam,chease=cheasepath)
         if profilesflag:
            profilesdata = read_profiles(fpath=profilespath,setParam=SetParam,chease=cheasepath)
         if importedflag:
            importeddata = read_imported(importeddata=imported,setParam=SetParam,chease=cheasepath)
         if   rhopsiflag:
              expeq['rhopsi'] = cheasedata['rhopsi'][:]
              expeq['nrhomesh']= npy.size(expeq['rhopsi'])
         elif rhotorflag:
              expeq['rhotor'] = cheasedata['rhotor'][:]
              expeq['nrhomesh']= npy.size(expeq['rhotor'])
    elif nrhotype[1] in [1,'eqdsk']:
         if cheaseflag:
            cheasedata = read_chease(fpath=cheasepath,setParam=SetParam,eqdsk=eqdskpath)
         if eqdskflag:
            eqdskdata = read_eqdsk(fpath=eqdskpath,setParam=SetParam)
         if expeqflag:
            expeqdata = read_expeq(fpath=expeqpath,setParam=SetParam,eqdsk=eqdskpath)
         if exptnzflag:
            exptnzdata = read_exptnz(fpath=exptnzpath,setParam=SetParam,eqdsk=eqdskpath)
         if iterdbflag:
            iterdbdata = read_iterdb(fpath=iterdbpath,setParam=SetParam,eqdsk=eqdskpath)
         if profilesflag:
            profilesdata = read_profiles(fpath=profilespath,setParam=SetParam,eqdsk=eqdskpath)
         if importedflag:
            importeddata = read_imported(importeddata=imported,setParam=SetParam,eqdsk=eqdskpath)
         if   rhopsiflag:
              expeq['rhopsi'] = eqdskdata['rhopsi'][:]
              expeq['nrhomesh']= npy.size(expeq['rhopsi'])
         elif rhotorflag:
              expeq['rhotor'] = eqdskdata['rhotor'][:]
              expeq['nrhomesh']= npy.size(expeq['rhotor'])
    elif nrhotype[1] in [7,'imported']:
         if cheaseflag:
            cheasedata = read_chease(fpath=cheasepath,setParam=SetParam,imported=imported)
         if eqdskflag:
            eqdskdata = read_eqdsk(fpath=eqdskpath,setParam=SetParam,imported=imported)
         if expeqflag:
            expeqdata = read_expeq(fpath=expeqpath,setParam=SetParam,imported=imported)
         if exptnzflag:
            exptnzdata = read_exptnz(fpath=exptnzpath,setParam=SetParam,imported=imported)
         if iterdbflag:
            iterdbdata = read_iterdb(fpath=iterdbpath,setParam=SetParam,imported=imported)
         if profilesflag:
            profilesdata = read_profiles(fpath=profilespath,setParam=SetParam,imported=imported)
         if importedflag:
            importeddata = read_imported(importeddata=imported,setParam=SetParam)
         if   rhopsiflag:
              expeq['rhopsi'] = imported['rhopsi'][:]
              expeq['nrhomesh']= npy.size(expeq['rhopsi'])
         elif rhotorflag:
              expeq['rhotor'] = imported['rhotor'][:]
              expeq['nrhomesh']= npy.size(expeq['rhotor'])
    elif nrhotype[1] in [None]:
         if cheaseflag:
            cheasedata = read_chease(fpath=cheasepath,setParam=SetParam)
            if   rhopsiflag:
                 expeq['rhopsi'] = cheasedata['rhopsi'][:]
                 expeq['nrhomesh']= npy.size(expeq['rhopsi'])
            elif rhotorflag:
                 expeq['rhotor'] = cheasedata['rhotor'][:]
                 expeq['nrhomesh']= npy.size(expeq['rhotor'])
         if eqdskflag:
            eqdskdata = read_eqdsk(fpath=eqdskpath,setParam=SetParam)
            if   rhopsiflag:
                 expeq['rhopsi'] = eqdskdata['rhopsi'][:]
                 expeq['nrhomesh']= npy.size(expeq['rhopsi'])
            elif rhotorflag:
                 expeq['rhotor'] = eqdskdata['rhotor'][:]
                 expeq['nrhomesh']= npy.size(expeq['rhotor'])
         if expeqflag:
            expeqdata = read_expeq(fpath=expeqpath,setParam=SetParam)
            if   rhopsiflag:
                 expeq['rhopsi'] = expeqdata['rhopsi'][:]
                 expeq['nrhomesh']= npy.size(expeq['rhopsi'])
            elif rhotorflag:
                 expeq['rhotor'] = expeqdata['rhotor'][:]
                 expeq['nrhomesh']= npy.size(expeq['rhotor'])
         if exptnzflag:
            exptnzdata = read_exptnz(fpath=exptnzpath,setParam=SetParam)
            if   rhopsiflag:
                 expeq['rhopsi'] = exptnzdata['rhopsi'][:]
                 expeq['nrhomesh']= npy.size(expeq['rhopsi'])
            elif rhotorflag:
                 expeq['rhotor'] = exptnzdata['rhotor'][:]
                 expeq['nrhomesh']= npy.size(expeq['rhotor'])
         if iterdbflag:
            iterdbdata = read_iterdb(fpath=iterdbpath,setParam=SetParam)
            if   rhopsiflag:
                 expeq['rhopsi'] = iterdbdata['rhopsi'][:]
                 expeq['nrhomesh']= npy.size(expeq['rhopsi'])
            elif rhotorflag:
                 expeq['rhotor'] = iterdbdata['rhotor'][:]
                 expeq['nrhomesh']= npy.size(expeq['rhotor'])
         if profilesflag:
            profilesdata = read_profiles(fpath=profilespath,setParam=SetParam)
            if   rhopsiflag:
                 expeq['rhopsi'] = profilesdata['rhopsi'][:]
                 expeq['nrhomesh']= npy.size(expeq['rhopsi'])
            elif rhotorflag:
                 expeq['rhotor'] = profilesdata['rhotor'][:]
                 expeq['nrhomesh']= npy.size(expeq['rhotor'])
         if importedflag:
            importeddata = read_imported(importeddata=imported,setParam=SetParam)
            if   rhopsiflag:
                 expeq['rhopsi'] = importeddata['rhopsi'][:]
                 expeq['nrhomesh']= npy.size(expeq['rhopsi'])
            elif rhotorflag:
                 expeq['rhotor'] = importeddata['rhotor'][:]
                 expeq['nrhomesh']= npy.size(expeq['rhotor'])


    if   nppfun[1] in [0,'chease'] and cheaseflag:
         if importedflag and 'pressure' in importeddata:
            expeq['pedge'] = mu0*importeddata['pressure'][-1]/expeq['B0EXP']**2
         else:
            expeq['pedge'] = mu0*cheasedata['pressure'][-1]/expeq['B0EXP']**2
         if   nppfun[0] in [4,'pprime','pprimen']:
              if importedflag and 'pprime' in importeddata:
                 expeq['pprime'] = mu0*importeddata['pprime']*expeq['R0EXP']**2/expeq['B0EXP']
              else:
                 expeq['pprime'] = mu0*cheasedata['pprime']*expeq['R0EXP']**2/expeq['B0EXP']
         elif nppfun[0] in [8,'pressure','pressuren','p','pn']:
              if importedflag and 'pressure' in importeddata:
                 expeq['pressure'] = mu0*importeddata['pressure']/expeq['B0EXP']**2
              else:
                 expeq['pressure'] = mu0*cheasedata['pressure']/expeq['B0EXP']**2

    elif nppfun[1] in [1,'eqdsk'] and eqdskflag:
         if importedflag and 'pressure' in importeddata:
            expeq['pedge'] = mu0*importeddata['pressure'][-1]/expeq['B0EXP']**2
         else:
            expeq['pedge'] = mu0*eqdskdata['pressure'][-1]/expeq['B0EXP']**2
         if   nppfun[0] in [4,'pprime','pprimen']:
              if importedflag and 'pprime' in importeddata:
                 expeq['pprime'] = mu0*importeddata['pprime']*expeq['R0EXP']**2/expeq['B0EXP']
              else:
                 expeq['pprime'] = mu0*eqdskdata['pprime']*expeq['R0EXP']**2/expeq['B0EXP']
         elif nppfun[0] in [8,'pressure','pressuren','p','pn']:
              if importedflag and 'pressure' in importeddata:
                 expeq['pressure'] = mu0*importeddata['pressure']/expeq['B0EXP']**2
              else:
                 expeq['pressure'] = mu0*eqdskdata['pressure']/expeq['B0EXP']**2

    elif nppfun[1] in [2,'expeq'] and expeqflag:
         if importedflag and 'pressure' in importeddata:
            expeq['pedge'] = mu0*importeddata['pressure'][-1]/expeq['B0EXP']**2
         else:
            expeq['pedge'] = expeqdata['pedge']
         if   nppfun[0] in [4,'pprime','pprimen']:
              if importedflag and 'pprime' in importeddata:
                 expeq['pprime'] = mu0*importeddata['pprime']*expeq['R0EXP']**2/expeq['B0EXP']
              else:
                 expeq['pprime'] = expeqdata['pprime'][:]
         elif nppfun[0] in [8,'pressure','pressuren','p','pn']:
              if importedflag and 'pressure' in importeddata:
                 expeq['pressure'] = mu0*importeddata['pressure']/expeq['B0EXP']**2
              else:
                 expeq['pressure'] = expeqdata['pressure'][:]

    elif nppfun[1] in [3,'exptnz'] and exptnzflag:
         if importedflag and 'pressure' in importeddata:
            expeq['pedge'] = mu0*importeddata['pressure'][-1]/expeq['B0EXP']**2
         else:
            expeq['pedge'] = mu0*exptnzdata['pressure'][-1]/expeq['B0EXP']**2
         if   nppfun[0] in [4,'pprime','pprimen']:
              if importedflag and 'pprime' in importeddata:
                 expeq['pprime'] = mu0*importeddata['pprime']*expeq['R0EXP']**2/expeq['B0EXP']
              else:
                 expeq['pprime'] = mu0*exptnzdata['pprime']*expeq['R0EXP']**2/expeq['B0EXP']
         elif nppfun[0] in [8,'pressure','pressuren','p','pn']:
              if importedflag and 'pressure' in importeddata:
                 expeq['pressure'] = mu0*importeddata['pressure']/expeq['B0EXP']**2
              else:
                 expeq['pressure'] = mu0*exptnzdata['pressure']/expeq['B0EXP']**2

    elif nppfun[1] in [4,'profiles'] and profilesflag:
         if importedflag and 'pressure' in importeddata:
            expeq['pedge'] = mu0*importeddata['pressure'][-1]/expeq['B0EXP']**2
         else:
            expeq['pedge'] = mu0*profilesdata['pressure'][-1]/expeq['B0EXP']**2
         if   nppfun[0] in [4,'pprime','pprimen']:
              if importedflag and 'pprime' in importeddata:
                 expeq['pprime'] = mu0*importeddata['pprime']*expeq['R0EXP']**2/expeq['B0EXP']
              else:
                 expeq['pprime'] = mu0*profilesdata['pprime']*expeq['R0EXP']**2/expeq['B0EXP']
         elif nppfun[0] in [8,'pressure','pressuren','p','pn']:
              if importedflag and 'pressure' in importeddata:
                 expeq['pressure'] = mu0*importeddata['pressure']/expeq['B0EXP']**2
              else:
                 expeq['pressure'] = mu0*profilesdata['pressure']/expeq['B0EXP']**2

    elif nppfun[1] in [5,'iterdb'] and iterdbflag:
         if importedflag and 'pressure' in importeddata:
            expeq['pedge'] = mu0*importeddata['pressure'][-1]/expeq['B0EXP']**2
         else:
            expeq['pedge'] = mu0*iterdbdata['pressure'][-1]/expeq['B0EXP']**2
         if   nppfun[0] in [4,'pprime','pprimen']:
              if importedflag and 'pprime' in importeddata:
                 expeq['pprime'] = mu0*importeddata['pprime']*expeq['R0EXP']**2/expeq['B0EXP']
              else:
                 expeq['pprime'] = mu0*iterdbdata['pprime']*expeq['R0EXP']**2/expeq['B0EXP']
         elif nppfun[0] in [8,'pressure','pressuren','p','pn']:
              if importedflag and 'pressure' in importeddata:
                 expeq['pressure'] = mu0*importeddata['pressure']/expeq['B0EXP']**2
              else:
                 expeq['pressure'] = mu0*iterdbdata['pressure']/expeq['B0EXP']**2

    elif nppfun[1] in [7,'imported'] and importedflag:
         if 'pressure' in importeddata:
            expeq['pedge'] = mu0*importeddata['pressure'][-1]/expeq['B0EXP']**2
         else:
            expeq['pedge'] = mu0*eqdskdata['pressure'][-1]/expeq['B0EXP']**2
         if   nppfun[0] in [4,'pprime','pprimen']:
              expeq['pprime'] = mu0*importeddata['pprime']*expeq['R0EXP']**2/expeq['B0EXP']
         elif nppfun[0] in [8,'pressure','pressuren','p','pn']:
              expeq['pressure'] = mu0*importeddata['pressure']/expeq['B0EXP']**2


    if   nsttp[1] in [0,'chease'] and cheaseflag:
         if   nsttp[0] in [1,'ffprime','ffprimen']:
              if importedflag and 'ffprime' in importeddata:
                 expeq['ffprime'] = importeddata['ffprime']/expeq['B0EXP']
              else:
                 expeq['ffprime'] = cheasedata['ffprime']/expeq['B0EXP']
         elif nsttp[0] in [2,'istr','istrn','istar','istarn']:
              if importedflag and 'Istr' in importeddata:
                 expeq['Istr'] = importeddata['Istr']*mu0/expeq['R0EXP']/expeq['B0EXP']
              else:
                 expeq['Istr'] = cheasedata['Istr']*mu0/expeq['R0EXP']/expeq['B0EXP']
         elif nsttp[0] in [3,'iprl','iprln','iparallel','iparalleln']:
              if importedflag and 'Iprl' in importeddata:
                 expeq['Iprl'] = importeddata['Iprl']*mu0/expeq['R0EXP']/expeq['B0EXP']
              else:
                 expeq['Iprl'] = cheasedata['Iprl']*mu0/expeq['R0EXP']/expeq['B0EXP']
         elif nsttp[0] in [4,'jprl','jprln','jparallel','jparalleln']:
              if importedflag and 'Jprl' in importeddata:
                 expeq['Jprl'] = importeddata['Jprl']*mu0*expeq['R0EXP']/expeq['B0EXP']
              else:
                 expeq['Jprl'] = cheasedata['Jprl']*mu0*expeq['R0EXP']/expeq['B0EXP']
         elif nsttp[0] in [5,'q','qpsi','qtor']:
              if importedflag and 'q' in importeddata:
                 expeq['q'] = importeddata['q'][:]
              else:
                 expeq['q'] = cheasedata['q'][:]

    elif nsttp[1] in [2,'expeq'] and expeqflag:
         if   nsttp[0] in [1,'ffprime','ffprimen']:
              if importedflag and 'ffprime' in importeddata:
                 expeq['ffprime'] = importeddata['ffprime']/expeq['B0EXP']
              else:
                 expeq['ffprime'] = expeqdata['ffprime'][:]
         elif nsttp[0] in [2,'istr','istrn','istar','istarn']:
              if importedflag and 'Istr' in importeddata:
                 expeq['Istr'] = importeddata['Istr']*mu0/expeq['R0EXP']/expeq['B0EXP']
              else:
                 expeq['Istr'] = expeqdata['Istr'][:]
         elif nsttp[0] in [3,'iprl','iprln','iparallel','iparalleln']:
              if importedflag and 'Iprl' in importeddata:
                 expeq['Iprl'] = importeddata['Iprl']*mu0/expeq['R0EXP']/expeq['B0EXP']
              else:
                 expeq['Iprl'] = expeqdata['Iprl'][:]
         elif nsttp[0] in [4,'jprl','jprln','jparallel','jparalleln']:
              if importedflag and 'Jprl' in importeddata:
                 expeq['Jprl'] = importeddata['Jprl']*mu0*expeq['R0EXP']/expeq['B0EXP']
              else:
                 expeq['Jprl'] = expeqdata['Jprl'][:]
         elif nsttp[0] in [5,'q','qpsi','qtor']:
              if importedflag and 'q' in importeddata:
                 expeq['q'] = importeddata['q']
              else:
                 expeq['q'] = expeqdata['q'][:]

    elif nsttp[1] in [1,'eqdsk'] and eqdskflag:
         if   nsttp[0] in [1,'ffprime','ffprimen']:
              if importedflag and 'ffprime' in importeddata:
                 expeq['ffprime'] = importeddata['ffprime']/expeq['B0EXP']
              else:
                 expeq['ffprime'] = eqdskdata['ffprime']/expeq['B0EXP']
         elif nsttp[0] in [5,'q','qpsi','qtor']:
              if importedflag and 'q' in importeddata:
                 expeq['q'] = importeddata['q']
              else:
                 expeq['q'] = eqdskdata['q']
         else:
              raise ValueError("FATAL: nsttp[0] (must) = 1 or 5.")

    elif nsttp[1] in [7,'imported'] and importedflag:
         if   nsttp[0] in [1,'ffprime','ffprimen']:
              expeq['ffprime'] = importeddata['ffprime']/expeq['B0EXP']
         elif nsttp[0] in [2,'istr','istrn','istar','istarn']:
              expeq['Istr'] = importeddata['Istr']*mu0/expeq['R0EXP']/expeq['B0EXP']
         elif nsttp[0] in [3,'iprl','iprln','iparallel','iparalleln']:
              expeq['Iprl'] = importeddata['Iprl']*mu0/expeq['R0EXP']/expeq['B0EXP']
         elif nsttp[0] in [4,'jprl','jprln','jparallel','jparalleln']:
              expeq['Jprl'] = importeddata['Jprl']*mu0*expeq['R0EXP']/expeq['B0EXP']
         elif nsttp[0] in [5,'q','qpsi','qtor']:
              expeq['q'] = importeddata['q'][:]

    if outfile:
       ofh = open("EXPEQ",'w')
       ofh.write('%18.8E\n'           % expeq['aspect'])
       ofh.write('%18.8E\n'           % expeq['zgeom'])
       ofh.write('%18.8E\n'           % expeq['pedge'])
       ofh.write('%5d\n'              % expeq['nRZmesh'])
       for i in range(expeq['nRZmesh']):
           ofh.write('%18.8E%18.8E\n' % (expeq['rbound'][i],expeq['zbound'][i]))
       ofh.write('%5d%5d\n'           % (expeq['nrhomesh'],nppfun[0]))
       ofh.write('%5d%5d\n'           % (nsttp[0],nrhotype[0]))


       if   nrhotype[0] in [0,'rhopsi']:
            for i in range(expeq['nrhomesh']):
                ofh.write('%18.8E\n'       % expeq['rhopsi'][i])
       elif nrhotype[0] in [1,'rhotor']:
            for i in range(expeq['nrhomesh']):
                ofh.write('%18.8E\n'       % expeq['rhotor'][i])

       if   nppfun[0] in [4,'pprime','pprimen']:
            for i in range(expeq['nrhomesh']):
                ofh.write('%18.8E\n'       % expeq['pprime'][i])
       elif nppfun[0] in [8,'pressure','pressuren','p','pn']:
            for i in range(expeq['nrhomesh']):
                ofh.write('%18.8E\n'       % expeq['pressure'][i])

       if   nsttp[0] in [1,'ffprime','ffprimen']:
            for i in range(expeq['nrhomesh']):
                ofh.write('%18.8E\n'       % expeq['ffprime'][i])
       elif nsttp[0] in [2,'istar','istarn','istr','istrn']:
            for i in range(expeq['nrhomesh']):
                ofh.write('%18.8E\n'       % expeq['Istr'][i])
       elif nsttp[0] in [3,'iprl','iprln','iparallel','iparalleln']:
            for i in range(expeq['nrhomesh']):
                ofh.write('%18.8E\n'       % expeq['Iprl'][i])
       elif nsttp[0] in [4,'jprl','jprln','jparallel','jparalleln']:
            for i in range(expeq['nrhomesh']):
                ofh.write('%18.8E\n'       % expeq['Jprl'][i])
       elif nsttp[0] in [5,'q']:
            for i in range(expeq['nrhomesh']):
                ofh.write('%18.8E\n'       % expeq['q'][i])
       ofh.close()

    return expeq


def read_imported(importeddata,setParam={},**kwargs):
    rhopsiflag = False; rhotorflag = False
    if 'nrhomesh' in setParam:
        if   setParam['nrhomesh'] in [0,'rhopsi']: rhopsiflag = True
        elif setParam['nrhomesh'] in [1,'rhotor']: rhotorflag = True
    else:                                          rhopsiflag = True

    eqdskflag    = False
    cheaseflag   = False
    interpflag   = False
    for key,value in kwargs.items():
        if   key in ['chease','cheasedata','cheasefpath']:
             if    type(value)==str and os.path.isfile(value.strip()):
                   cheasedata = read_chease(fpath=value.strip())
             else: raise IOError('%s file not found!' % value.strip())
             if 'rhopsi' in cheasedata: rhopsi = cheasedata['rhopsi'][:]; interpflag = True
             if 'rhotor' in cheasedata: rhotor = cheasedata['rhotor'][:]; interpflag = True
             cheaseflag = True
        elif key in ['eqdsk','eqdskdata','eqdskfpath']:
             if    type(value)==str and os.path.isfile(value.strip()):
                   eqdskdata = read_eqdsk(fpath=value.strip())
             else: raise IOError('%s file not found!' % value.strip())
             if 'rhopsi' in eqdskdata: rhopsi = eqdskdata['rhopsi'][:]; interpflag = True
             if 'rhotor' in eqdskdata: rhotor = eqdskdata['rhotor'][:]; interpflag = True
             eqdskflag = True

    IMPORTEDdata = {}
    if not importeddata: return IMPORTEDdata

    if 'rhopsi' in importeddata:
       IMPORTEDdata['rhopsi'] = importeddata['rhopsi'][:]
       if   rhopsiflag: rhotype = 'rhopsi'
       elif rhotorflag: rhotype = 'rhotor'

    if 'rhotor' in importeddata:
       IMPORTEDdata['rhotor'] = importeddata['rhotor'][:]
       if   rhopsiflag: rhotype = 'rhopsi'
       elif rhotorflag: rhotype = 'rhotor'


    if 'rhotor' not in importeddata or 'rhopsi' not in importeddata:
       raise IOError('rhopsi and/or rhotor NOT provided. EXIT!')

    if interpflag and 'Te' in importeddata:
       if   rhopsiflag and rhotype=='rhopsi':
            IMPORTEDdata['Te'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['Te'],rhopsi)
       elif rhotorflag and rhotype=='rhotor':
            IMPORTEDdata['Te'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['Te'],rhotor)
       elif rhotorflag and rhotype=='rhopsi':
            IMPORTEDdata['Te'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['Te'],rhopsi,rhotor,rhotor)
       elif rhopsiflag and rhotype=='rhotor':
            IMPORTEDdata['Te'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['Te'],rhotor,rhopsi,rhopsi)
    elif 'Te' in importeddata:
         IMPORTEDdata['Te'] = importeddata['Te']

    if interpflag and 'ne' in importeddata:
       if   rhopsiflag and rhotype=='rhopsi':
            IMPORTEDdata['ne'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['ne'],rhopsi)
       elif rhotorflag and rhotype=='rhotor':
            IMPORTEDdata['ne'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['ne'],rhotor)
       elif rhotorflag and rhotype=='rhopsi':
            IMPORTEDdata['ne'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['ne'],rhopsi,rhotor,rhotor)
       elif rhopsiflag and rhotype=='rhotor':
            IMPORTEDdata['ne'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['ne'],rhotor,rhopsi,rhopsi)
    elif 'ne' in importeddata:
         IMPORTEDdata['ne'] = importeddata['ne']

    if interpflag and 'Ti' in importeddata:
       if   rhopsiflag and rhotype=='rhopsi':
            IMPORTEDdata['Ti'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['Ti'],rhopsi)
       elif rhotorflag and rhotype=='rhotor':
            IMPORTEDdata['Ti'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['Ti'],rhotor)
       elif rhotorflag and rhotype=='rhopsi':
            IMPORTEDdata['Ti'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['Ti'],rhopsi,rhotor,rhotor)
       elif rhopsiflag and rhotype=='rhotor':
            IMPORTEDdata['Ti'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['Ti'],rhotor,rhopsi,rhopsi)
    elif 'Ti' in importeddata:
         IMPORTEDdata['Ti'] = importeddata['Ti']

    if interpflag and 'ni' in importeddata:
       if   rhopsiflag and rhotype=='rhopsi':
            IMPORTEDdata['ni'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['ni'],rhopsi)
       elif rhotorflag and rhotype=='rhotor':
            IMPORTEDdata['ni'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['ni'],rhotor)
       elif rhotorflag and rhotype=='rhopsi':
            IMPORTEDdata['ni'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['ni'],rhopsi,rhotor,rhotor)
       elif rhopsiflag and rhotype=='rhotor':
            IMPORTEDdata['ni'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['ni'],rhotor,rhopsi,rhopsi)
    elif 'ni' in importeddata:
         IMPORTEDdata['ni'] = importeddata['ni']

    if interpflag and 'q' in importeddata:
       if   rhopsiflag and rhotype=='rhopsi':
            IMPORTEDdata['q'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['q'],rhopsi)
       elif rhotorflag and rhotype=='rhotor':
            IMPORTEDdata['q'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['q'],rhotor)
       elif rhotorflag and rhotype=='rhopsi':
            IMPORTEDdata['q'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['q'],rhopsi,rhotor,rhotor)
       elif rhopsiflag and rhotype=='rhotor':
            IMPORTEDdata['q'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['q'],rhotor,rhopsi,rhopsi)
    elif 'q' in importeddata:
         IMPORTEDdata['q'] = importeddata['q']

    if interpflag and 'Zeff' in importeddata:
       if   rhopsiflag and rhotype=='rhopsi':
            IMPORTEDdata['Zeff'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['Zeff'],rhopsi)
       elif rhotorflag and rhotype=='rhotor':
            IMPORTEDdata['Zeff'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['Zeff'],rhotor)
       elif rhotorflag and rhotype=='rhopsi':
            IMPORTEDdata['Zeff'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['Zeff'],rhopsi,rhotor,rhotor)
       elif rhopsiflag and rhotype=='rhotor':
            IMPORTEDdata['Zeff'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['Zeff'],rhotor,rhopsi,rhopsi)
    elif 'Zeff' in importeddata:
         IMPORTEDdata['Zeff'] = importeddata['Zeff']

    if interpflag and 'pressure' in importeddata:
       if   rhopsiflag and rhotype=='rhopsi':
            IMPORTEDdata['pressure'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['pressure'],rhopsi)
       elif rhotorflag and rhotype=='rhotor':
            IMPORTEDdata['pressure'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['pressure'],rhotor)
       elif rhotorflag and rhotype=='rhopsi':
            IMPORTEDdata['pressure'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['pressure'],rhopsi,rhotor,rhotor)
       elif rhopsiflag and rhotype=='rhotor':
            IMPORTEDdata['pressure'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['pressure'],rhotor,rhopsi,rhopsi)
    elif 'pressure' in importeddata:
         IMPORTEDdata['pressure'] = importeddata['pressure']

    if interpflag and 'pprime' in importeddata:
       if   rhopsiflag and rhotype=='rhopsi':
            IMPORTEDdata['pprime'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['pprime'],rhopsi)
       elif rhotorflag and rhotype=='rhotor':
            IMPORTEDdata['pprime'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['pprime'],rhotor)
       elif rhotorflag and rhotype=='rhopsi':
            IMPORTEDdata['pprime'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['pprime'],rhopsi,rhotor,rhotor)
       elif rhopsiflag and rhotype=='rhotor':
            IMPORTEDdata['pprime'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['pprime'],rhotor,rhopsi,rhopsi)
    elif 'pprime' in importeddata:
         IMPORTEDdata['pprime'] = importeddata['pprime']

    if interpflag and 'Istr' in importeddata:
       if   rhopsiflag and rhotype=='rhopsi':
            IMPORTEDdata['Istr'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['Istr'],rhopsi)
       elif rhotorflag and rhotype=='rhotor':
            IMPORTEDdata['Istr'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['Istr'],rhotor)
       elif rhotorflag and rhotype=='rhopsi':
            IMPORTEDdata['Istr'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['Istr'],rhopsi,rhotor,rhotor)
       elif rhopsiflag and rhotype=='rhotor':
            IMPORTEDdata['Istr'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['Istr'],rhotor,rhopsi,rhopsi)
    elif 'Istr' in importeddata:
         IMPORTEDdata['Istr'] = importeddata['Istr']

    if interpflag and 'Iprl' in importeddata:
       if   rhopsiflag and rhotype=='rhopsi':
            IMPORTEDdata['Iprl'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['Iprl'],rhopsi)
       elif rhotorflag and rhotype=='rhotor':
            IMPORTEDdata['Iprl'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['Iprl'],rhotor)
       elif rhotorflag and rhotype=='rhopsi':
            IMPORTEDdata['Iprl'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['Iprl'],rhopsi,rhotor,rhotor)
       elif rhopsiflag and rhotype=='rhotor':
            IMPORTEDdata['Iprl'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['Iprl'],rhotor,rhopsi,rhopsi)
    elif 'Iprl' in importeddata:
         IMPORTEDdata['Iprl'] = importeddata['Iprl']

    if interpflag and 'Jprl' in importeddata:
       if   rhopsiflag and rhotype=='rhopsi':
            IMPORTEDdata['Jprl'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['Jprl'],rhopsi)
       elif rhotorflag and rhotype=='rhotor':
            IMPORTEDdata['Jprl'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['Jprl'],rhotor)
       elif rhotorflag and rhotype=='rhopsi':
            IMPORTEDdata['Jprl'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['Jprl'],rhopsi,rhotor,rhotor)
       elif rhopsiflag and rhotype=='rhotor':
            IMPORTEDdata['Jprl'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['Jprl'],rhotor,rhopsi,rhopsi)
    elif 'Jprl' in importeddata:
         IMPORTEDdata['Jprl'] = importeddata['Jprl']

    if interpflag and 'ffprime' in importeddata:
       if   rhopsiflag and rhotype=='rhopsi':
            IMPORTEDdata['ffprime'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['ffprime'],rhopsi)
       elif rhotorflag and rhotype=='rhotor':
            IMPORTEDdata['ffprime'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['ffprime'],rhotor)
       elif rhotorflag and rhotype=='rhopsi':
            IMPORTEDdata['ffprime'] = mathtools.interp(IMPORTEDdata['rhopsi'],importeddata['ffprime'],rhopsi,rhotor,rhotor)
       elif rhopsiflag and rhotype=='rhotor':
            IMPORTEDdata['ffprime'] = mathtools.interp(IMPORTEDdata['rhotor'],importeddata['ffprime'],rhotor,rhopsi,rhopsi)
    elif 'ffprime' in importeddata:
         IMPORTEDdata['ffprime'] = importeddata['ffprime']

    if interpflag:
       IMPORTEDdata['rhopsi'] = rhopsi[:]
       IMPORTEDdata['rhotor'] = rhotor[:]

    return IMPORTEDdata


def read_exptnz(fpath,setParam={},**kwargs):
    if os.path.isfile(fpath) == False:
       errorFunc = traceback.extract_stack(limit=2)[-2][3]
       errorLine = traceback.extract_stack(limit=2)[-2][1]
       errorFile = traceback.extract_stack(limit=2)[-2][2]
       errMSG    = 'Call %s line %5d in file %s Failed.\n'
       errMSG   += 'Fatal: file %s not found.'
       raise IOError(errMSG %(errorFunc,errorLine,errorFile,fpath))

    rhopsiflag = False; rhotorflag = False
    if 'nrhomesh' in setParam:
        if   setParam['nrhomesh'] in [0,'rhopsi']: rhopsiflag = True
        elif setParam['nrhomesh'] in [1,'rhotor']: rhotorflag = True
    else:                                          rhopsiflag = True

    eqdskflag    = False
    cheaseflag   = False
    interpflag   = False
    importedflag = False
    for key,value in kwargs.items():
        if   key in ['chease','cheasedata','cheasefpath']:
             if    type(value)==str and os.path.isfile(value.strip()):
                   cheasedata = read_chease(fpath=value.strip())
             else: raise IOError('%s file not found!' % value.strip())
             if 'rhopsi' in cheasedata: rhopsi = cheasedata['rhopsi'][:]; interpflag = True
             if 'rhotor' in cheasedata: rhotor = cheasedata['rhotor'][:]; interpflag = True
             cheaseflag = True
        elif key in ['eqdsk','eqdskdata','eqdskfpath']:
             if    type(value)==str and os.path.isfile(value.strip()):
                   eqdskdata = read_eqdsk(fpath=value.strip())
             else: raise IOError('%s file not found!' % value.strip())
             if 'rhopsi' in eqdskdata: rhopsi = eqdskdata['rhopsi'][:]; interpflag = True
             if 'rhotor' in eqdskdata: rhotor = eqdskdata['rhotor'][:]; interpflag = True
             eqdskflag = True
        elif key in ['imported','external','other']:
             imported = value.copy()
             if 'rhopsi' in imported: rhopsi = imported['rhopsi'][:]; interpflag = True
             if 'rhotor' in imported: rhotor = imported['rhotor'][:]; interpflag = True
             importedflag = True

    ofh = open(fpath,'r')
    EXPTNZOUT = ofh.readlines()
    ofh.close()

    n_rho   = int(EXPTNZOUT[0].split()[0])
    rhotype =     EXPTNZOUT[0].split()[1].strip()[0:6]

    EXPTNZdata               = {}

    EXPTNZdata['Zi'] = 1.0
    EXPTNZdata['Zz'] = 6.0

    if   rhotype=='rhopsi':
         EXPTNZdata['rhopsi'] = npy.array(EXPTNZOUT[0*n_rho+1:1*n_rho+1],dtype=float)
    elif rhotype=='rhotor':
         EXPTNZdata['rhotor'] = npy.array(EXPTNZOUT[0*n_rho+1:1*n_rho+1],dtype=float)
    EXPTNZdata['Te']          = npy.array(EXPTNZOUT[1*n_rho+1:2*n_rho+1],dtype=float)
    EXPTNZdata['ne']          = npy.array(EXPTNZOUT[2*n_rho+1:3*n_rho+1],dtype=float)
    EXPTNZdata['Zeff']        = npy.array(EXPTNZOUT[3*n_rho+1:4*n_rho+1],dtype=float)
    EXPTNZdata['Ti']          = npy.array(EXPTNZOUT[4*n_rho+1:5*n_rho+1],dtype=float)
    EXPTNZdata['ni']          = npy.array(EXPTNZOUT[5*n_rho+1:6*n_rho+1],dtype=float)

    if interpflag:
       if   rhopsiflag and rhotype=='rhopsi':
            EXPTNZdata['Te']     = mathtools.interp(EXPTNZdata['rhopsi'],EXPTNZdata['Te'],rhopsi)
            EXPTNZdata['Ti']     = mathtools.interp(EXPTNZdata['rhopsi'],EXPTNZdata['Ti'],rhopsi)
            EXPTNZdata['ne']     = mathtools.interp(EXPTNZdata['rhopsi'],EXPTNZdata['ne'],rhopsi) 
            EXPTNZdata['ni']     = mathtools.interp(EXPTNZdata['rhopsi'],EXPTNZdata['ni'],rhopsi)
            EXPTNZdata['Zeff']   = mathtools.interp(EXPTNZdata['rhopsi'],EXPTNZdata['Zeff'],rhopsi)
            EXPTNZdata['rhopsi'] = rhopsi[:]
            EXPTNZdata['rhotor'] = rhotor[:]
       elif rhotorflag and rhotype=='rhotor':
            EXPTNZdata['Te']     = mathtools.interp(EXPTNZdata['rhotor'],EXPTNZdata['Te'],rhotor)
            EXPTNZdata['Ti']     = mathtools.interp(EXPTNZdata['rhotor'],EXPTNZdata['Ti'],rhotor)
            EXPTNZdata['ne']     = mathtools.interp(EXPTNZdata['rhotor'],EXPTNZdata['ne'],rhotor) 
            EXPTNZdata['ni']     = mathtools.interp(EXPTNZdata['rhotor'],EXPTNZdata['ni'],rhotor)
            EXPTNZdata['Zeff']   = mathtools.interp(EXPTNZdata['rhotor'],EXPTNZdata['Zeff'],rhotor)
            EXPTNZdata['rhopsi'] = rhopsi[:]
            EXPTNZdata['rhotor'] = rhotor[:]
       elif rhotorflag and rhotype=='rhopsi':
            EXPTNZdata['Te']     = mathtools.interp(EXPTNZdata['rhopsi'],EXPTNZdata['Te'],rhopsi,rhotor,rhotor)
            EXPTNZdata['Ti']     = mathtools.interp(EXPTNZdata['rhopsi'],EXPTNZdata['Ti'],rhopsi,rhotor,rhotor)
            EXPTNZdata['ne']     = mathtools.interp(EXPTNZdata['rhopsi'],EXPTNZdata['ne'],rhopsi,rhotor,rhotor) 
            EXPTNZdata['ni']     = mathtools.interp(EXPTNZdata['rhopsi'],EXPTNZdata['ni'],rhopsi,rhotor,rhotor)
            EXPTNZdata['Zeff']   = mathtools.interp(EXPTNZdata['rhopsi'],EXPTNZdata['Zeff'],rhopsi,rhotor,rhotor)
            EXPTNZdata['rhopsi'] = rhopsi[:]
            EXPTNZdata['rhotor'] = rhotor[:]
       elif rhopsiflag and rhotype=='rhotor':
            EXPTNZdata['Te']     = mathtools.interp(EXPTNZdata['rhotor'],EXPTNZdata['Te'],rhotor,rhopsi,rhopsi)
            EXPTNZdata['Ti']     = mathtools.interp(EXPTNZdata['rhotor'],EXPTNZdata['Ti'],rhotor,rhopsi,rhopsi)
            EXPTNZdata['ne']     = mathtools.interp(EXPTNZdata['rhotor'],EXPTNZdata['ne'],rhotor,rhopsi,rhopsi) 
            EXPTNZdata['ni']     = mathtools.interp(EXPTNZdata['rhotor'],EXPTNZdata['ni'],rhotor,rhopsi,rhopsi)
            EXPTNZdata['Zeff']   = mathtools.interp(EXPTNZdata['rhotor'],EXPTNZdata['Zeff'],rhotor,rhopsi,rhopsi)
            EXPTNZdata['rhopsi'] = rhopsi[:]
            EXPTNZdata['rhotor'] = rhotor[:]
    elif rhopsiflag and rhotype=='rhotor':
         print("WARNING: setParam['nrhomesh'] = 0 or rhopsi, but the path to a target rhopsi is not provided.")
         print("         Converting the profiles to poloidal (psi) coordinates could not be done, and")
         print("         all profiles are provided in the toroidal (phi) coordinates.")
    elif rhotorflag and rhotype=='rhopsi':
         print("WARNING: setParam['nrhomesh'] = 1 or rhotor, but the path to a target rhotor is not provided.")
         print("         Converting the profiles to toroidal (phi) coordinates could not be done, and")
         print("         all profiles are provided in the poloidal (psi) coordinates.")

    EXPTNZdata['nz']         = EXPTNZdata['Zeff']*EXPTNZdata['ne']
    EXPTNZdata['nz']        -= EXPTNZdata['ni']*EXPTNZdata['Zi']**2
    EXPTNZdata['nz']        /= EXPTNZdata['Zz']**2

    EXPTNZdata['pressure']   = EXPTNZdata['Te']*EXPTNZdata['ne']
    EXPTNZdata['pressure']  += EXPTNZdata['Ti']*EXPTNZdata['ni']
    EXPTNZdata['pressure']  += EXPTNZdata['Ti']*EXPTNZdata['nz']
    EXPTNZdata['pressure']  *= 1.602e-19

    if   rhopsiflag:
         EXPTNZdata['pprime']= mathtools.derivative(x=EXPTNZdata['rhopsi'],fx=EXPTNZdata['pressure'],method='CubicSpline')
    elif rhotorflag:
         EXPTNZdata['pprime']= mathtools.derivative(x=EXPTNZdata['rhopsi'],fx=EXPTNZdata['pressure'],method='CubicSpline')
    elif rhotype=='rhopsi':
         EXPTNZdata['pprime']= mathtools.derivative(x=EXPTNZdata['rhopsi'],fx=EXPTNZdata['pressure'],method='CubicSpline')
    elif rhotype=='rhotor':
         EXPTNZdata['pprime']= mathtools.derivative(x=EXPTNZdata['rhotor'],fx=EXPTNZdata['pressure'],method='CubicSpline')

    return EXPTNZdata

def write_exptnz(setParam={},**kwargs):
    '''
    nrhomesh=[rho_type(0:rhopsi,1:rhotor),rho_src(0:chease,1:eqdsk)]
    eprofile=[eprofile_src(0:chease,3:exptnz,4:profiles,5:iterdb)]
    iprofile=[iprofile_src(0:chease,3:exptnz,4:profiles,5:iterdb)]
    '''

    eqdskflag    = False
    interpflag   = False
    cheaseflag   = False
    iterdbflag   = False
    exptnzflag   = False
    importedflag = False
    profilesflag = False
    for key,value in kwargs.items():
        if key in ['chease','cheasedata','cheasefpath']:
           if os.path.isfile(value.strip()):
              cheasepath = value.strip()
              cheaseflag = True

        if key in ['eqdsk','eqdskdata','eqdskfpath']:
           if os.path.isfile(value.strip()):
              eqdskpath = value.strip()
              eqdskflag = True

        if key in ['expeq','expeqdata','expeqfpath']:
           if os.path.isfile(value.strip()):
              expeqpath = value.strip()
              expeqflag = True

        if key in ['exptnz','exptnzdata','exptnzfpath']:
           if os.path.isfile(value.strip()):
              exptnzpath = value.strip()
              exptnzflag = True

        if key in ['profiles','profilesdata','profilesfpath']:
           if os.path.isfile(value.strip()):
              profilespath = value.strip()
              profilesflag = True

        if key in ['iterdb','iterdbdata','iterdbfpath']:
           if os.path.isfile(value.strip()):
              iterdbpath = value.strip()
              iterdbflag = True

        if key in ['imported','external','others']:
              imported = value.copy()
              importedflag = True

    if not (cheaseflag or exptnzflag or profilesflag or iterdbflag or importedflag):
       raise IOError('FATAL: NO VALID INPUT PROFILES AVAILABLE. EXIT!')

    if   PYTHON3:
         setParamKeys = list(setParam.keys())
    elif PYTHON2:
         setParamKeys = setParam.keys()

    if 'outfile' in setParam:
        outfile = setParam['outfile']
    else:
        outfile = True

    if 'nrhomesh' in setParam:
       if   type(setParam['nrhomesh'])==list: 
            if   type(setParam['nrhomesh'][0])==float: setParam['nrhomesh'][0] = int(setParam['nrhomesh'][0])
            elif type(setParam['nrhomesh'][0])==str:   setParam['nrhomesh'][0] = setParam['nrhomesh'][0].lower()
            if   type(setParam['nrhomesh'][1])==float: setParam['nrhomesh'][1] = int(setParam['nrhomesh'][1])
            elif type(setParam['nrhomesh'][1])==str:   setParam['nrhomesh'][1] = setParam['nrhomesh'][1].lower()
            elif      setParam['nrhomesh'][1] ==None:  setParam['nrhomesh'][1] = None
            nrhotype = setParam['nrhomesh'][:]
       elif type(setParam['nrhomesh'])==int:
            if   cheaseflag:   nrhotype=[setParam['nrhomesh'],0]
            elif eqdskflag:    nrhotype=[setParam['nrhomesh'],1]
            elif importedflag: nrhotype=[setParam['nrhomesh'],7]
    else:
            if   cheaseflag:   nrhotype=[0,0]
            elif eqdskflag:    nrhotype=[0,1]
            elif importedflag: nrhotype=[0,7]

    if 'eprofile' in setParam: 
            if   type(setParam['eprofile'])==float: setParam['eprofile'] = int(setParam['eprofile'])
            elif type(setParam['eprofile'])==str:   setParam['eprofile'] = setParam['eprofile'].lower()
            eprofile= setParam['eprofile']
    else:
            if   cheaseflag:   eprofile=0
            elif exptnzflag:   eprofile=3
            elif profilesflag: eprofile=4
            elif iterdbflag:   eprofile=5
            elif importedflag: eprofile=7

    if 'iprofile' in setParam:
            if   type(setParam['iprofile'])==float: setParam['iprofile'] = int(setParam['iprofile'])
            elif type(setParam['iprofile'])==str:   setParam['iprofile'] = setParam['iprofile'].lower()
            iprofile= setParam['iprofile']
    else:
            if   cheaseflag:   iprofile=0
            elif exptnzflag:   iprofile=3
            elif profilesflag: iprofile=4
            elif iterdbflag:   iprofile=5
            elif importedflag: iprofile=7

    if   nrhotype[1] in [0,'chease'] and not cheaseflag:
         raise IOError("FATAL: nrhotype=chease and chease.h5 file is not available. EXIT!")
    elif nrhotype[1] in [1,'eqdsk'] and not eqdskflag:
         raise IOError("FATAL: nrhotype=eqdsk and eqdsk file is not available. EXIT!")
    elif nrhotype[1] in [7,'imported'] and not importedflag:
         raise IOError("FATAL: nrhotype=imported and imported data is not available. EXIT!")

    if   eprofile in [0,'chease'] and not cheaseflag:
         raise IOError("FATAL: eprofile=chease and chease.h5 is not available. EXIT!")
    elif eprofile in [3,'exptnz'] and not exptnzflag:
         raise IOError("FATAL: eprofile=exptnz and exptnz is not available. EXIT!")
    elif eprofile in [4,'profiles'] and not profilesflag:
         raise IOError("FATAL: eprofile=profiples and profiles is not available. EXIT!")
    elif eprofile in [5,'iterdb'] and not iterdbflag:
         raise IOError("FATAL: eprofile=iterdb and iterdb is not available. EXIT!")
    elif eprofile in [7,'imported'] and not importedflag:
         raise IOError("FATAL: eprofile=imported and no profiles are imported. EXIT!")

    if   iprofile in [0,'chease'] and not cheaseflag:
         raise IOError("FATAL: iprofile=chease and chease.h5 is not available. EXIT!")
    elif iprofile in [3,'exptnz'] and not exptnzflag:
         raise IOError("FATAL: iprofile=exptnz and exptnz is not available. EXIT!")
    elif iprofile in [4,'profiles'] and not profilesflag:
         raise IOError("FATAL: iprofile=profiples and profiles is not available. EXIT!")
    elif iprofile in [5,'iterdb'] and not iterdbflag:
         raise IOError("FATAL: iprofile=iterdb and iterdb is not available. EXIT!")
    elif iprofile in [5,'iterdb'] and not iterdbflag:
         raise IOError("FATAL: iprofile=iterdb and iterdb is not available. EXIT!")
    elif eprofile in [7,'imported'] and not importedflag:
         raise IOError("FATAL: iprofile=imported and no profiles are imported. EXIT!")

    if   nrhotype[0] in [0,'rhopsi']:
         rhopsiflag = True
         rhotorflag = False
    elif nrhotype[0] in [1,'rhotor']:
         rhopsiflag = False
         rhotorflag = True

    SetParam={'nrhomesh':nrhotype[0]}

    if   nrhotype[1] in [0,'chease']:
         if cheaseflag:
            cheasedata = read_chease(fpath=cheasepath,setParam=SetParam)
         if exptnzflag:
            exptnzdata = read_exptnz(fpath=exptnzpath,setParam=SetParam,chease=cheasepath)
         if iterdbflag:
            iterdbdata = read_iterdb(fpath=iterdbpath,setParam=SetParam,chease=cheasepath)
         if profilesflag:
            profilesdata = read_profiles(fpath=profilespath,setParam=SetParam,chease=cheasepath)
         if importedflag:
            importeddata = read_imported(importeddata=imported,setParam=SetParam,chease=cheasepath)
    elif nrhotype[1] in [1,'eqdsk']:
         if cheaseflag:
            cheasedata = read_chease(fpath=cheasepath,setParam=SetParam,eqdsk=eqdskpath)
         if exptnzflag:
            exptnzdata = read_exptnz(fpath=exptnzpath,setParam=SetParam,eqdsk=eqdskpath)
         if iterdbflag:
            iterdbdata = read_iterdb(fpath=iterdbpath,setParam=SetParam,eqdsk=eqdskpath)
         if profilesflag:
            profilesdata = read_profiles(fpath=profilespath,setParam=SetParam,eqdsk=eqdskpath)
         if importedflag:
            importeddata = read_imported(importeddata=imported,setParam=SetParam,eqdsk=eqdskpath)
    elif nrhotype[1] in [7,'imported']:
         if cheaseflag:
            cheasedata = read_chease(fpath=cheasepath,setParam=SetParam,imported=imported)
         if exptnzflag:
            exptnzdata = read_exptnz(fpath=exptnzpath,setParam=SetParam,imported=imported)
         if iterdbflag:
            iterdbdata = read_iterdb(fpath=iterdbpath,setParam=SetParam,imported=imported)
         if profilesflag:
            profilesdata = read_profiles(fpath=profilespath,setParam=SetParam,imported=imported)
         if importedflag:
            importeddata = read_imported(importeddata=imported,setParam=SetParam)
    elif nrhotype[1] in [None]:
         if cheaseflag:
            cheasedata = read_chease(fpath=cheasepath,setParam=SetParam)
         if exptnzflag:
            exptnzdata = read_exptnz(fpath=exptnzpath,setParam=SetParam)
         if iterdbflag:
            iterdbdata = read_iterdb(fpath=iterdbpath,setParam=SetParam)
         if profilesflag:
            profilesdata = read_profiles(fpath=profilespath,setParam=SetParam)
         if importedflag:
            importeddata = read_imported(importeddata=imported,setParam=SetParam)

    exptnz = {}

    if   eprofile in [0,'chease'] and cheaseflag:
         if   rhopsiflag: exptnz['rhopsi'] = cheasedata['rhopsi']
         elif rhotorflag: exptnz['rhotor'] = cheasedata['rhotor']
         if importedflag and 'Te' in importeddata:
            exptnz['Te'] = importeddata['Te']
         else:
            exptnz['Te'] = cheasedata['Te']
         if importedflag and 'ne' in importeddata:
            exptnz['ne'] = importeddata['ne']
         else:
            exptnz['ne'] = cheasedata['ne']
    elif eprofile in [3,'exptnz'] and exptnzflag:
         if   rhopsiflag: exptnz['rhopsi'] = exptnzdata['rhopsi']
         elif rhotorflag: exptnz['rhotor'] = exptnzdata['rhotor']
         if importedflag and 'Te' in importeddata:
            exptnz['Te'] = importeddata['Te']
         else:
            exptnz['Te'] = exptnzdata['Te']
         if importedflag and 'ne' in importeddata:
            exptnz['ne'] = importeddata['ne']
         else:
            exptnz['ne'] = exptnzdata['ne']
    elif eprofile in [4,'profiles'] and profilesflag:
         if   rhopsiflag: exptnz['rhopsi'] = profilesdata['rhopsi']
         elif rhotorflag: exptnz['rhotor'] = profilesdata['rhotor']
         if importedflag and 'Te' in importeddata:
            exptnz['Te'] = importeddata['Te']
         else:
            exptnz['Te'] = profilesdata['Te']
         if importedflag and 'ne' in importeddata:
            exptnz['ne'] = importeddata['ne']
         else:
            exptnz['ne'] = profilesdata['ne']
    elif eprofile in [5,'iterdb'] and iterdbflag:
         if   rhopsiflag: exptnz['rhopsi'] = iterdbdata['rhopsi']
         elif rhotorflag: exptnz['rhotor'] = iterdbdata['rhotor']
         if importedflag and 'Te' in importeddata:
            exptnz['Te'] = importeddata['Te']
         else:
            exptnz['Te'] = iterdbdata['Te']
         if importedflag and 'ne' in importeddata:
            exptnz['ne'] = importeddata['ne']
         else:
            exptnz['ne'] = iterdbdata['ne']
    elif eprofile in [7,'imported'] and importedflag:
         if   rhopsiflag: exptnz['rhopsi'] = importeddata['rhopsi']
         elif rhotorflag: exptnz['rhotor'] = importeddata['rhotor']
         exptnz['Te'] = importeddata['Te']
         exptnz['ne'] = importeddata['ne']

    if   iprofile in [0,'chease'] and cheaseflag:
         if importedflag and 'Ti' in importeddata:
            exptnz['Ti']   = importeddata['Ti']
         else:
            exptnz['Ti']   = cheasedata['Ti']
         if importedflag and 'ni' in importeddata:
            exptnz['ni']   = importeddata['ni']
         else:
            exptnz['ni']   = cheasedata['ni']
         if importedflag and 'nz' in importeddata:
            exptnz['nz']   = importeddata['nz']
         else:
            exptnz['nz']   = cheasedata['nz']
         if importedflag and 'Zeff' in importeddata:
            exptnz['Zeff'] = importeddata['Zeff']
         else:
            exptnz['Zeff'] = cheasedata['Zeff']
    elif iprofile in [3,'exptnz'] and exptnzflag:
         if importedflag and 'Ti' in importeddata:
            exptnz['Ti']   = importeddata['Ti']
         else:
            exptnz['Ti']   = exptnzdata['Ti']
         if importedflag and 'ni' in importeddata:
            exptnz['ni']   = importeddata['ni']
         else:
            exptnz['ni']   = exptnzdata['ni']
         if importedflag and 'nz' in importeddata:
            exptnz['nz']   = importeddata['nz']
         else:
            exptnz['nz']   = exptnzdata['nz']
         if importedflag and 'Zeff' in importeddata:
            exptnz['Zeff'] = importeddata['Zeff']
         else:
            exptnz['Zeff'] = exptnzdata['Zeff']
    elif iprofile in [4,'profiles'] and profilesflag:
         if importedflag and 'Ti' in importeddata:
            exptnz['Ti']   = importeddata['Ti']
         else:
            exptnz['Ti']   = profilesdata['Ti']
         if importedflag and 'ni' in importeddata:
            exptnz['ni']   = importeddata['ni']
         else:
            exptnz['ni']   = profilesdata['ni']
         if importedflag and 'nz' in importeddata:
            exptnz['nz']   = importeddata['nz']
         else:
            exptnz['nz']   = profilesdata['nz']
         if importedflag and 'Zeff' in importeddata:
            exptnz['Zeff'] = importeddata['Zeff']
         else:
            exptnz['Zeff'] = profilesdata['Zeff']
    elif iprofile in [5,'iterdb'] and iterdbflag:
         if importedflag and 'Ti' in importeddata:
            exptnz['Ti']   = importeddata['Ti']
         else:
            exptnz['Ti']   = iterdbdata['Ti']
         if importedflag and 'ni' in importeddata:
            exptnz['ni']   = importeddata['ni']
         else:
            exptnz['ni']   = iterdbdata['ni']
         if importedflag and 'nz' in importeddata:
            exptnz['nz']   = importeddata['nz']
         else:
            exptnz['nz']   = iterdbdata['nz']
         if importedflag and 'Zeff' in importeddata:
            exptnz['Zeff'] = importeddata['Zeff']
         else:
            exptnz['Zeff'] = iterdbdata['Zeff']
    elif iprofile in [7,'imported'] and importedflag:
         exptnz['Ti']   = importeddata['Ti']
         exptnz['ni']   = importeddata['ni']
         if 'Zeff' in importeddata:
            if type(importeddata['Zeff']) == list:
               exptnz['Zeff'] = importeddata['Zeff']
            else:
               exptnz['Zeff'] = npy.ones(npy.size(importeddata['Ti']))
               exptnz['Zeff']*= importeddata['Zeff']
         else:
               exptnz['Zeff'] = npy.ones(npy.size(importeddata['Ti']))
         if 'nz' in importeddata:
            exptnz['nz']= importeddata['nz']
         else:
            exptnz['nz']= npy.zeros(npy.size(importeddata['Ti']))

    if 'Zi' not in exptnz: exptnz['Zi'] = 1.0
    if 'Zz' not in exptnz: exptnz['Zz'] = 6.0

    if 'nz' not in exptnz:
       exptnz['nz'] = exptnz['Zeff']*exptnz['ne']
       exptnz['nz']-= exptnz['ni']*exptnz['Zi']**2
       exptnz['nz']/= exptnz['Zz']**2

    if   rhopsiflag: rhosize = npy.size(exptnz['rhopsi'])
    elif rhotorflag: rhosize = npy.size(exptnz['rhotor'])

    if outfile:
       ofh = open("EXPTNZ",'w')
       if   nrhotype[0] in [0,'rhopsi']:
            rhosize = npy.size(exptnz['rhopsi'])
            ofh.write('%5d rhopsi,  Te,   ne,   Zeff,   Ti,   ni  profiles\n' % rhosize)
            for i in range(rhosize): ofh.write('%16.6E\n' % exptnz['rhopsi'][i])
       elif nrhotype[0] in [1,'rhotor']:
            rhosize = npy.size(exptnz['rhotor'])
            ofh.write('%5d rhotor,  Te,   ne,   Zeff,   Ti,   ni  profiles\n' % rhosize)
            for i in range(rhosize): ofh.write('%16.6E\n' % exptnz['rhotor'][i])
       for i in range(rhosize): ofh.write('%16.6E\n' % exptnz['Te'][i])
       for i in range(rhosize): ofh.write('%16.6E\n' % exptnz['ne'][i])
       for i in range(rhosize): ofh.write('%16.6E\n' % exptnz['Zeff'][i])
       for i in range(rhosize): ofh.write('%16.6E\n' % exptnz['Ti'][i])
       for i in range(rhosize): ofh.write('%16.6E\n' % exptnz['ni'][i])
       ofh.close()

    return exptnz


def read_profiles_file(pfpath,setParam={}):
   #Developed by Ehab Hassan on 2019-02-17
    from scipy.interpolate import CubicSpline
    if not os.path.isfile(pfpath):
       print('Fatal: file %s not found.' % pfpath)
       sys.exit()

    ofh = open(pfpath,'r')

    profiles = {}
    units    = {}
    while True:
          recs = ofh.readline().split()
          if   len(recs)>4:
               nrec = int(recs[0])
               ary0=npy.zeros(nrec)
               ary1=npy.zeros(nrec)
               ary2=npy.zeros(nrec)
               var0 = str(recs[1]).lower()
               var1 = str(recs[2]).lower()
               var2 = str(recs[3]).lower()
               for i in range(nrec):
                   recs = ofh.readline().split()
                   ary0[i] = float(recs[0])
                   ary1[i] = float(recs[1])
                   ary2[i] = float(recs[2])
               profiles[var0]=ary0
               profiles[var1]=ary1
               profiles[var2]=ary2
          elif len(recs)==4:
               nrec = int(recs[0])
               ary0=npy.zeros(nrec)
               ary1=npy.zeros(nrec)
               ary2=npy.zeros(nrec)
               var0 = str(recs[1]).lower()
               temp = str(recs[2]).lower()
               var1 = temp[:temp.index("(")]
               if   var1.strip() in ['ne','ni','nb','nz1']:
                    powr = int(temp[temp.index("^")+1:temp.index("/")])
                    unit = 'm^{-3}'
               elif var1.strip() in ['te','ti']:
                    unit = 'eV'
               elif var1.strip() in ['pb','ptot']:
                    unit = 'Pa'
               elif var1.strip() in ['vtor1','vpol1']:
                    unit = 'm/s'
               else:
                    unit = temp[temp.index("(")+1:temp.index(")")]
               var2 = str(recs[3]).lower()
               for i in range(nrec):
                   recs    = ofh.readline().split()
                   ary0[i] = float(recs[0])
                   ary1[i] = float(recs[1])
                   ary2[i] = float(recs[2])
               if var0 in profiles.keys():
                  CS = CubicSpline(ary0,ary1)
                  units[var1]    = unit
                  profiles[var1] = CS(profiles[var0])
                  CS = CubicSpline(ary0,ary2)
                  profiles[var2] = CS(profiles[var0])
               else:
                  profiles[var0] = ary0
                  units[var1]    = unit
                  profiles[var1] = ary1
                  profiles[var2] = ary2
               if   var1.strip() in ['ne','ni','nb','nz1']:
                    profiles[var1] *= 10.0**powr
               elif var1.strip() in ['te','ti','pb','ptot','vtor1','vpol1']:
                    profiles[var1] *= 1.0e3
          else:
               break
    ofh.close()

    if 'rhotor' in setParam:
       if 'eqdskfpath' not in setParam:
          eqdskfpath = raw_input('Path to EQDSK file: ')
       elif 'eqdskfpath' in setParam:
          eqdskfpath = setParam['eqdskfpath']

       eqdskdata = read_efit_file(fpath=eqdskfpath)

       qpsifn  = interp1d(eqdskdata['PSIN'],eqdskdata['qpsi'])
       qpsi    = qpsifn(profiles['psinorm'])
       qpsifn  = interp1d(profiles['psinorm'],qpsi)
       psinorm = npy.linspace(profiles['psinorm'][0],profiles['psinorm'][-1],setParam['rhotor'])
       qpsi    = qpsifn(psinorm)
       phinorm = npy.zeros_like(psinorm)
       for i in range(1,npy.size(qpsi)):
           x = psinorm[:i+1]
           y = qpsi[:i+1]
           phinorm[i] = npy.trapz(y,x)
       profiles['rhotor'] = npy.sqrt(phinorm)

    return profiles,units



def read_profiles(fpath,setParam={},**kwargs):
    if os.path.isfile(fpath) == False:
       errorFunc = traceback.extract_stack(limit=2)[-2][3]
       errorLine = traceback.extract_stack(limit=2)[-2][1]
       errorFile = traceback.extract_stack(limit=2)[-2][2]
       errMSG    = 'Call %s line %5d in file %s Failed.\n'
       errMSG   += 'Fatal: file %s not found.'
       raise IOError(errMSG %(errorFunc,errorLine,errorFile,fpath))

    if 'Zeffprofile' in setParam: Zeffprofile = setParam['Zeffprofile']
    else:                         Zeffprofile = True

    rhopsiflag = False; rhotorflag = False
    if 'nrhomesh' in setParam:
        if   setParam['nrhomesh'] in [0,'rhopsi']: rhopsiflag = True
        elif setParam['nrhomesh'] in [1,'rhotor']: rhotorflag = True
    else:                                          rhopsiflag = True

    eqdskflag    = False
    interpflag   = False
    cheaseflag   = False
    importedflag = False
    for key,value in kwargs.items():
        if   key in ['chease','cheasedata','cheasefpath']:
             if    type(value)==str and os.path.isfile(value.strip()):
                   cheasedata = read_chease(fpath=value.strip())
             else: raise IOError('%s file not found!' % value.strip())
             if 'rhopsi' in cheasedata: rhopsi = cheasedata['rhopsi'][:]
             if 'rhotor' in cheasedata: rhotor = cheasedata['rhotor'][:]
             if 'PSIN'   in cheasedata: psi    = cheasedata['PSIN'][:];   interpflag = True
             if 'PHIN'   in cheasedata: phi    = cheasedata['PHIN'][:];   interpflag = True
             cheaseflag = True
        elif key in ['eqdsk','eqdskdata','eqdskfpath']:
             if    type(value)==str and os.path.isfile(value.strip()):
                   eqdskdata = read_eqdsk(fpath=value.strip())
             else: raise IOError('%s file not found!' % value.strip())
             if 'rhopsi' in eqdskdata: rhopsi = eqdskdata['rhopsi'][:]
             if 'rhotor' in eqdskdata: rhotor = eqdskdata['rhotor'][:]
             if 'PSIN'   in eqdskdata: psi    = eqdskdata['PSIN'][:];  interpflag = True
             if 'PHIN'   in eqdskdata: phi    = eqdskdata['PHIN'][:];  interpflag = True
             eqdskflag = True
        elif key in ['imported','external','other']:
             imported = value.copy()
             if 'rhopsi' in imported: rhopsi = imported['rhopsi'][:]
             if 'rhotor' in imported: rhotor = imported['rhotor'][:]
             psi    = importeddata['rhopsi']**2;  interpflag = True
             phi    = importeddata['rhotor']**2;  interpflag = True
             importedflag = True

    profiles,units = read_profiles_file(fpath.strip())

    PROFILESdata            = {}
    PROFILESdata['rhopsi']  = npy.sqrt(profiles['psinorm'])
    PROFILESdata['PSIN']    = profiles['psinorm'][:]
    PROFILESdata['Te']      = profiles['te'][:]
    PROFILESdata['Ti']      = profiles['ti'][:]
    if 'tb' in profiles:
       PROFILESdata['Tb']      = profiles['tb'][:]
    PROFILESdata['ne']      = profiles['ne'][:]
    PROFILESdata['ni']      = profiles['ni'][:]
    PROFILESdata['nb']      = profiles['nb'][:]
    PROFILESdata['nz']      = profiles['nz1'][:]

    PROFILESdata['Pb']      = profiles['pb'][:]
    PROFILESdata['Vpol']    = profiles['vpol1'][:]
    PROFILESdata['Vtor']    = profiles['vtor1'][:]

    PROFILESdata['pressure']= profiles['ptot'][:]

    nrhopsi                 = npy.size(PROFILESdata['rhopsi'])

    PROFILESdata['Zi']      = profiles['z'][1]
    PROFILESdata['Zz']      = profiles['z'][0]

    PROFILESdata['Zeff']    = PROFILESdata['ni']*PROFILESdata['Zi']**2
    PROFILESdata['Zeff']   += PROFILESdata['nz']*PROFILESdata['Zz']**2
    PROFILESdata['Zeff']   /= PROFILESdata['ne']

    if not Zeffprofile:
       Zeff_array = npy.array([1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0])
       Zeff_mean  = npy.mean(PROFILESdata['Zeff'])
       Zeff_diff  = abs(Zeff_array-Zeff_mean)
       Zeff_value = Zeff_array[Zeff_diff==min(Zeff_diff)][0]
       PROFILESdata['Zeff']    = npy.ones(nrhopsi)*Zeff_value

    PROFILESdata['pprime']  = mathtools.derivative(x=PROFILESdata['PSIN'],fx=PROFILESdata['pressure'],method='CubicSpline')

    if   interpflag:
         if   rhopsiflag:
              PROFILESdata['Te']       = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['Te'],psi)
              PROFILESdata['Ti']       = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['Ti'],psi)
              if 'tb' in profiles:
                 PROFILESdata['Tb']    = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['Tb'],psi)
              PROFILESdata['ne']       = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['ne'],psi)
              PROFILESdata['ni']       = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['ni'],psi)
              PROFILESdata['nb']       = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['nb'],psi)
              PROFILESdata['nz']       = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['nz'],psi)
              PROFILESdata['Pb']       = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['Pb'],psi)
              PROFILESdata['Vtor']     = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['Vtor'],psi)
              PROFILESdata['Vpol']     = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['Vpol'],psi) 
              PROFILESdata['Zeff']     = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['Zeff'],psi)
              PROFILESdata['pprime']   = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['pprime'],psi)
              PROFILESdata['pressure'] = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['pressure'],psi)
              PROFILESdata['PSIN']     = psi[:]
              PROFILESdata['PHIN']     = phi[:]
              PROFILESdata['rhopsi']   = rhopsi[:]
              PROFILESdata['rhotor']   = rhotor[:]
         elif rhotorflag:
              PROFILESdata['Te']       = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['Te'],psi,phi,phi)
              PROFILESdata['Ti']       = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['Ti'],psi,phi,phi)
              if 'tb' in profiles:
                 PROFILESdata['Tb']    = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['Tb'],psi,phi,phi)
              PROFILESdata['ne']       = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['ne'],psi,phi,phi)
              PROFILESdata['ni']       = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['ni'],psi,phi,phi)
              PROFILESdata['nb']       = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['nb'],psi,phi,phi)
              PROFILESdata['nz']       = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['nz'],psi,phi,phi)
              PROFILESdata['Pb']       = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['Pb'],psi,phi,phi)
              PROFILESdata['Vtor']     = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['Vtor'],psi,phi,phi)
              PROFILESdata['Vpol']     = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['Vpol'],psi,phi,phi)
              PROFILESdata['Zeff']     = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['Zeff'],psi,phi,phi)
              PROFILESdata['pprime']   = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['pprime'],psi,phi,phi)
              PROFILESdata['pressure'] = mathtools.interp(PROFILESdata['PSIN'],PROFILESdata['pressure'],psi,phi,phi)
              PROFILESdata['PSIN']     = psi[:]
              PROFILESdata['PHIN']     = phi[:]
              PROFILESdata['rhopsi']   = rhopsi[:]
              PROFILESdata['rhotor']   = rhotor[:]
    elif rhotorflag and not interpflag:
         print("WARNING: setParam['nrhomesh'] = 1 or rhotor, but the path to a target rhotor is not provided.")
         print("         Converting the profiles to toroidal (phi) coordinates could not be done, and")
         print("         all profiles are provided in the poloidal (psi) coordinates.")

    return PROFILESdata


def get_next(data_linesplit,lnum,num):
    sec_num_lines = num/6
    if num % 6 != 0:
        sec_num_lines += 1
    keep_going=1
    while keep_going:
        test=re.search('-DEPENDENT VARIABLE LABEL',data_linesplit[lnum])
        if test :
            quantity=data_linesplit[lnum].split()[0]
            units=data_linesplit[lnum].split()[1]
        test=re.search('DATA FOLLOW',data_linesplit[lnum])
        if test:
            keep_going=(1==2)
        lnum=lnum+1

    rhot=npy.empty(0)
    lnum0 = lnum
    for j in range(lnum0,lnum0+int(sec_num_lines)):
        for k in range(6):
            str_temp=data_linesplit[j][1+k*13:1+(k+1)*13]
            if(str_temp):
                temp=npy.array(data_linesplit[j][1+k*13:1+(k+1)*13],dtype='float')
                rhot=npy.append(rhot,temp)
        lnum=lnum+1
    lnum=lnum+1

    arr=npy.empty(0)
    lnum0 = lnum
    for j in range(lnum0,lnum0+int(sec_num_lines)):
        for k in range(6):
            str_temp=data_linesplit[j][1+k*13:1+(k+1)*13]
            if(str_temp):
                temp=npy.array(data_linesplit[j][1+k*13:1+(k+1)*13],dtype='float')
                arr=npy.append(arr,temp)
        lnum=lnum+1

    lnum_out=lnum
    try_again=1
    if len(data_linesplit)-lnum < 10:
        try_again=False
    return lnum_out, try_again,quantity,units,rhot,arr


def read_iterdb_file(filename):
    '''
    This Code is Written by: David R. Hatch
    It reads the iterdb file and returns three dictionaries,
    each diectionary has five quantities:
    electron density (NE) and temperature (TE),
    ion density (NM1) and temperatures (TI),
    impurity density (NM2), if any,
    rotational velocity (VROT).
    The three dictionaries provide the toroidal coordinate (rhotor),
    profiles, and units for each quantity.
    '''
    f=open(filename,'r')
    data_in=f.read()
    data_linesplit=data_in.split('\n')

    keep_going=1
    i=0
    while keep_going:
        test=re.search(';-# OF X PTS',data_linesplit[i])
        if test:
            num=data_linesplit[i].split()[0]
            num=float(num)
            num=int(num)
            keep_going=(1==2)
        if i == len(data_linesplit):
            keep_going=(1==2)
        i=i+1

    lnum=0
    try_again=1
    prof_out = {}
    rhot_out = {}
    units_out = {}
    while try_again:
        lnum,try_again,quantity,units,rhot,arr=get_next(data_linesplit,lnum,num)
        prof_out[quantity]=arr
        units_out[quantity]=units
        rhot_out[quantity]=rhot
    return rhot_out,prof_out,units_out


def read_iterdb(fpath,setParam={},**kwargs):
    if os.path.isfile(fpath) == False:
       errorFunc = traceback.extract_stack(limit=2)[-2][3]
       errorLine = traceback.extract_stack(limit=2)[-2][1]
       errorFile = traceback.extract_stack(limit=2)[-2][2]
       errMSG    = 'Call %s line %5d in file %s Failed.\n'
       errMSG   += 'Fatal: file %s not found.'
       raise IOError(errMSG %(errorFunc,errorLine,errorFile,fpath))

    rhopsiflag = False; rhotorflag = False
    if   'nrhomesh' in setParam:
         if   setParam['nrhomesh'] in [0,'rhopsi']: rhopsiflag = True
         elif setParam['nrhomesh'] in [1,'rhotor']: rhotorflag = True
    elif kwargs.items():
         rhopsiflag = True
    else:
         rhotorflag = True

    rhotors,profiles,units = read_iterdb_file(fpath)

    '''
    Normalizing the rhotor vectors before using them to interpolate the physical quantities
    '''
    if int(rhotors['NE'][-1]) == 1:
       rhotorNE    = rhotors['NE'][:]
       rhotorTE    = rhotors['TE'][:]
       rhotorNM1   = rhotors['NM1'][:]
       rhotorTI    = rhotors['TI'][:]
       if 'NM2' in profiles:
          rhotorNM2  = rhotors['NM2'][:]
       if 'VROT' in profiles:
          rhotorVROT = rhotors['VROT'][:]
    else:
       rhotorNE    = (rhotors['NE']-rhotors['NE'][0])/(rhotors['NE'][-1]-rhotors['NE'][0])
       rhotorTE    = (rhotors['TE']-rhotors['TE'][0])/(rhotors['TE'][-1]-rhotors['TE'][0])
       rhotorNM1   = (rhotors['NM1']-rhotors['NM1'][0])/(rhotors['NM1'][-1]-rhotors['NM1'][0])
       rhotorTI    = (rhotors['TI']-rhotors['TI'][0])/(rhotors['TI'][-1]-rhotors['TI'][0])
       if 'NM2' in profiles:
        rhotorNM2  = (rhotors['NM2']-rhotors['NM2'][0])/(rhotors['NM2'][-1]-rhotors['NM2'][0])
       if 'VROT' in profiles:
        rhotorVROT = (rhotors['VROT']-rhotors['VROT'][0])/(rhotors['VROT'][-1]-rhotors['VROT'][0])

    eqdskflag    = False
    cheaseflag   = False
    interpflag   = False
    importedflag = False
    for key,value in kwargs.items():
        if   key in ['chease','cheasedata','cheasefpath']:
             if    type(value)==str and os.path.isfile(value.strip()):
                   cheasedata = read_chease(fpath=value.strip())
             else: raise IOError('%s file not found!' % value.strip())
             if 'rhopsi' in cheasedata: rhopsi = cheasedata['rhopsi'][:]; interpflag = True
             if 'rhotor' in cheasedata: rhotor = cheasedata['rhotor'][:]; interpflag = True
             cheaseflag = True
        elif key in ['eqdsk','eqdskdata','eqdskfpath']:
             if    type(value)==str and os.path.isfile(value.strip()):
                   eqdskdata = read_eqdsk(fpath=value.strip())
             else: raise IOError('%s file not found!' % value.strip())
             if 'rhopsi' in eqdskdata: rhopsi = eqdskdata['rhopsi'][:]; interpflag = True
             if 'rhotor' in eqdskdata: rhotor = eqdskdata['rhotor'][:]; interpflag = True
             eqdskflag = True
        elif key in ['imported','external','other']:
             imported = value.copy()
             if 'rhopsi' in imported: rhopsi = imported['rhopsi'][:]; interpflag = True
             if 'rhotor' in imported: rhotor = imported['rhotor'][:]; interpflag = True
             importedflag = True
        else:
             rhotor  = rhotorNE[:]

    ITERDBdata       = {}
    ITERDBdata['Zi'] = 1.0
    ITERDBdata['Zz'] = 6.0

    if   rhopsiflag and interpflag:
         ITERDBdata['rhopsi'] = rhopsi
         ITERDBdata['rhotor'] = rhotor
    elif rhotorflag and interpflag:
         ITERDBdata['rhopsi'] = rhopsi
         ITERDBdata['rhotor'] = rhotor
    elif rhopsiflag and not interpflag:
         print("WARNING: setParam['nrhomesh'] = 0 or rhopsi, but the path to a target rhopsi is not provided.")
         print("         Converting the profiles to poloidal (psi) coordinates could not be done, and")
         print("         all profiles are provided in the toroidal (phi) coordinates.")
    else:
         ITERDBdata['rhotor'] = rhotorNE

    nrhosize = npy.size(ITERDBdata['rhotor'])

    if   rhopsiflag and interpflag:
         ITERDBdata['Te']      = mathtools.interp(rhotors['TE'],profiles['TE'],rhotor,rhopsi,rhopsi)
         ITERDBdata['Ti']      = mathtools.interp(rhotors['TI'],profiles['TI'],rhotor,rhopsi,rhopsi)
         ITERDBdata['ne']      = mathtools.interp(rhotors['NE'],profiles['NE'],rhotor,rhopsi,rhopsi)
         ITERDBdata['ni']      = mathtools.interp(rhotors['NM1'],profiles['NM1'],rhotor,rhopsi,rhopsi)
         if 'NM2' in profiles: 
            ITERDBdata['nz']   = mathtools.interp(rhotors['NM2'],profiles['NM2'],rhotor,rhopsi,rhopsi)
         else:
            ITERDBdata['nz']   = npy.zeros(nrhosize)
         if 'VROT' in profiles :
            ITERDBdata['Vrot']   = mathtools.interp(rhotors['VROT'],profiles['VROT'],rhotor,rhopsi,rhopsi)
         else:
            ITERDBdata['Vrot'] = npy.zeros(nrhosize)
    elif rhotorflag and interpflag:
         ITERDBdata['ne']      = mathtools.interp(rhotors['NE'],profiles['NE'],rhotor)
         ITERDBdata['Te']      = mathtools.interp(rhotors['TE'],profiles['TE'],rhotor)
         ITERDBdata['ni']      = mathtools.interp(rhotors['NM1'],profiles['NM1'],rhotor)
         ITERDBdata['Ti']      = mathtools.interp(rhotors['TI'],profiles['TI'],rhotor)
         if 'NM2' in profiles: 
            ITERDBdata['nz']   = mathtools.interp(rhotors['NM2'],profiles['NM2'],rhotor)
         else:
            ITERDBdata['nz']   = npy.zeros(nrhosize)
         if 'VROT' in profiles :
            ITERDBdata['Vrot'] = mathtools.interp(rhotors['VROT'],profiles['VROT'],rhotor)
         else:
            ITERDBdata['Vrot'] = npy.zeros(nrhosize)
    else:
         ITERDBdata['ne']      = mathtools.interp(rhotors['NE'],profiles['NE'],rhotorNE)
         ITERDBdata['Te']      = mathtools.interp(rhotors['TE'],profiles['TE'],rhotorNE)
         ITERDBdata['ni']      = mathtools.interp(rhotors['NM1'],profiles['NM1'],rhotorNE)
         ITERDBdata['Ti']      = mathtools.interp(rhotors['TI'],profiles['TI'],rhotorNE)
         if 'NM2' in profiles: 
            ITERDBdata['nz']   = mathtools.interp(rhotors['NM2'],profiles['NM2'],rhotorNE)
         else:
            ITERDBdata['nz']   = npy.zeros(nrhosize)
         if 'VROT' in profiles :
            ITERDBdata['Vrot'] = mathtools.interp(rhotors['VROT'],profiles['VROT'],rhotorNE)
         else:
            ITERDBdata['Vrot'] = npy.zeros(nrhosize)

    ITERDBdata['Zeff']    = ITERDBdata['ni']*ITERDBdata['Zi']**2
    ITERDBdata['Zeff']   += ITERDBdata['nz']*ITERDBdata['Zz']**2
    ITERDBdata['Zeff']   /= ITERDBdata['ne']

    ITERDBdata['pressure']   = ITERDBdata['Te']*ITERDBdata['ne']
    ITERDBdata['pressure']  += ITERDBdata['Ti']*ITERDBdata['ni']
    ITERDBdata['pressure']  += ITERDBdata['Ti']*ITERDBdata['nz']
    ITERDBdata['pressure']  *= 1.602e-19

    if   rhopsiflag and interpflag:
         ITERDBdata['pprime']= mathtools.derivative(x=ITERDBdata['rhopsi'],fx=ITERDBdata['pressure'],method='CubicSpline')
    else:
         ITERDBdata['pprime']= mathtools.derivative(x=ITERDBdata['rhotor'],fx=ITERDBdata['pressure'],method='CubicSpline')

    return ITERDBdata




