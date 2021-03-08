#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import mathtools
import cheasefiles
import numpy as npy

from glob       import glob
from subprocess import Popen, PIPE
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


if   sys.version_info.major == 3:
     PYTHON3 = True; PYTHON2 = False
elif sys.version_info.major == 2:
     PYTHON2 = True; PYTHON3 = False

mu0 = 4.0e-7*npy.pi

def findall(inlist,item):
    inds = []
    for i,character in enumerate(inlist):
        if character==item: inds.append(i)
    return inds


def create_namelist(setParam={}):
    wfh = open('chease_namelist','w')
    wfh.write('*** for EQDSK file copied onto EXPEQ file \n')
    wfh.write('*** cp this file to "chease_namelist" and run chease \n')
    wfh.write('***  \n')
    wfh.write('***  \n')
    wfh.write(' &EQDATA \n')
    if   'RELAX'     not in setParam: setParam['RELAX']     = 0.0
    if   'NBSEXPQ'   not in setParam: setParam['NBSEXPQ']   = 1111
    if   'NEQDSK'    not in setParam: setParam['NEQDSK']    = 0
    if   'NS'        not in setParam: setParam['NS']        = 256
    if   'NT'        not in setParam: setParam['NT']        = 256
    if   'NISO'      not in setParam: setParam['NISO']      = 256
    if   'NPSI'      not in setParam: setParam['NPSI']      = 1024
    if   'NCHI'      not in setParam: setParam['NCHI']      = 1024
    if   'NRBOX'     not in setParam: setParam['NRBOX']     = 60
    if   'NZBOX'     not in setParam: setParam['NZBOX']     = 60
    if   'NCSCAL'    not in setParam: setParam['NCSCAL']    = 4
    if   'NOPT'      not in setParam: setParam['NOPT']      = 0
    if   'NSURF'     not in setParam: setParam['NSURF']     = 6
    if   'NFUNC'     not in setParam: setParam['NFUNC']     = 4
    if   'NRHOMESH'  not in setParam: setParam['NRHOMESH']  = 0
    if   'NFUNRHO'   not in setParam: setParam['NFUNRHO']   = setParam['NRHOMESH']
    if   'NSTTP'     not in setParam: setParam['NSTTP']     = 1
    if   'NPROPT'    not in setParam: setParam['NPROPT']    = setParam['NSTTP']  if setParam['NPPFUN']==4 else -setParam['NSTTP']
    else:                             setParam['NPROPT']    = setParam['NPROPT'] if setParam['NPPFUN']==4 else -abs(setParam['NPROPT'])
    if   'NPPFUN'    not in setParam: setParam['NPPFUN']    = 8
    if   'NVERBOSE'  not in setParam: setParam['NVERBOSE']  = 4
    if   'NDIAGOP'   not in setParam: setParam['NDIAGOP']   = 1
    if   'NIDEAL'    not in setParam: setParam['NIDEAL']    = 9
    if   'NDIFPS'    not in setParam: setParam['NDIFPS']    = 0
    if   'NDIFT'     not in setParam: setParam['NDIFT']     = 1
    if   'NMESHC'    not in setParam: setParam['NMESHC']    = 1
    if   'NPOIDC'    not in setParam: setParam['NPOIDC']    = 2
    if   'SOLPDC'    not in setParam: setParam['SOLPDC']    = 0.7
    if   'CPLACE'    not in setParam: setParam['CPLACE']    = '0.9,0.95,0.99,1.0'
    if   'CWIDTH'    not in setParam: setParam['CWIDTH']    = '0.1,0.05,0.01,0.025'
    if   'NMESHPOL'  not in setParam: setParam['NMESHPOL']  = 1
    if   'SOLPDPOL'  not in setParam: setParam['SOLPDPOL']  = 0.1
    if   'NTURN'     not in setParam: setParam['NTURN']     = 20
    if   'NBLC0'     not in setParam: setParam['NBLC0']     = 16
    if   'NPPR'      not in setParam: setParam['NPPR']      = 24
    if   'NINMAP'    not in setParam: setParam['NINMAP']    = 40
    if   'NINSCA'    not in setParam: setParam['NINSCA']    = 40
    if   'NSYM'      not in setParam: setParam['NSYM']      = 0
    if   'NEGP'      not in setParam: setParam['NEGP']      = 0
    if   'NER'       not in setParam: setParam['NER']       = 2
    if   'EPSLON'    not in setParam: setParam['EPSLON']    = 1.0E-10
    if   'ETAEI'     not in setParam: setParam['ETAEI']     = 3.0
    if   'RPEOP'     not in setParam: setParam['RPEOP']     = 0.5
    if   'RZION'     not in setParam: setParam['RZION']     = 1.0
    if   'GAMMA'     not in setParam: setParam['GAMMA']     = 1.6666666667
    if   'AT3(1)'    not in setParam: setParam['AT3(1)']    = -0.69
    if   'TENSPROF'  not in setParam: setParam['TENSPROF']  = 0.0
    if   'TENSBND'   not in setParam: setParam['TENSBND']   = 0.0
    if   'cocos_in'  not in setParam: setParam['cocos_in']  = 2
    if   'cocos_out' not in setParam: setParam['cocos_out'] = 12

    if   'RELAX'     in setParam: wfh.write(' RELAX=%2.2f, '      % (float(setParam['RELAX'])))
    if   'NBSEXPQ'   in setParam: wfh.write(' NBSEXPQ=%04d, '     % (int(setParam['NBSEXPQ'])))
    if   'NEQDSK'    in setParam: wfh.write(' NEQDSK=%1d,     \n' % (int(setParam['NEQDSK'])))
    if   'NS'        in setParam: wfh.write(' NS=%4d, '           % (int(setParam['NS'])))
    if   'NT'        in setParam: wfh.write(' NT=%4d,         \n' % (int(setParam['NT'])))
    if   'NPSI'      in setParam: wfh.write(' NPSI=%4d, '         % (int(setParam['NPSI'])))
    if   'NCHI'      in setParam: wfh.write(' NCHI=%4d, '         % (int(setParam['NCHI'])))
    if   'NISO'      in setParam: wfh.write(' NISO=%4d,       \n' % (int(setParam['NISO'])))
    if   'NRBOX'     in setParam: wfh.write(' NRBOX=%4d, '        % (int(setParam['NRBOX'])))
    if   'NZBOX'     in setParam: wfh.write(' NZBOX=%4d,      \n' % (int(setParam['NZBOX'])))
    if   'NCSCAL'    in setParam: wfh.write(' NCSCAL=%1d, '       % (float(setParam['NCSCAL'])))
    if   'NOPT'      in setParam: wfh.write(' NOPT=%1d,       \n' % (int(setParam['NOPT'])))
    if   'NSURF'     in setParam: wfh.write(' NSURF=%1d, '        % (int(setParam['NSURF'])))
    if   'NFUNC'     in setParam: wfh.write(' NFUNC=%1d,      \n' % (int(setParam['NFUNC'])))
    if   'NPPFUN'    in setParam: wfh.write(' NPPFUN=%1d,     \n' % (int(setParam['NPPFUN'])))
    if   'NFUNRHO'   in setParam: wfh.write(' NFUNRHO=%1d, '      % (int(setParam['NFUNRHO'])))
    if   'NRHOMESH'  in setParam: wfh.write(' NRHOMESH=%1d,   \n' % (int(setParam['NRHOMESH'])))
    if   'NSTTP'     in setParam: wfh.write(' NSTTP=%1d, '        % (int(setParam['NSTTP'])))
    if   'NPROPT'    in setParam: wfh.write(' NPROPT=%1d,     \n' % ( int(setParam['NPROPT'])))
    if   'NVERBOSE'  in setParam: wfh.write(' NVERBOSE=%1d,   \n' % (int(setParam['NVERBOSE'])))
    if   'QSPEC'     in setParam: wfh.write(' QSPEC=%3.3f, '      % (float(setParam['QSPEC'])))
    if   'CSSPEC'    in setParam: wfh.write(' CSSPEC=%3.3f,   \n' % (float(setParam['CSSPEC'])))
    if   'TRIANG'    in setParam: wfh.write(' TRIANG=%3.3f,   \n' % (float(setParam['TRIANG'])))
    if   'R0'        in setParam: wfh.write(' R0=%10.8f, '        % (float(setParam['R0'])))
    if   'RZ0'       in setParam: wfh.write(' RZ0=%10.8f,     \n' % (float(setParam['RZ0'])))
    if   'RBOXLEN'   in setParam: wfh.write(' RBOXLEN=%3.3f, '    % (float(setParam['RBOXLEN'])))
    if   'ZBOXLEN'   in setParam: wfh.write(' ZBOXLEN=%3.3f, '    % (float(setParam['ZBOXLEN'])))
    if   'RBOXLFT'   in setParam: wfh.write(' RBOXLFT=%3.3f,  \n' % (float(setParam['RBOXLFT'])))
    if   'R0EXP'     in setParam: wfh.write(' R0EXP=%3.3f, '      % (float(setParam['R0EXP'])))
    if   'B0EXP'     in setParam: wfh.write(' B0EXP=%3.3f,    \n' % (float(setParam['B0EXP'])))
    if   'NDIAGOP'   in setParam: wfh.write(' NDIAGOP=%1d, '      % (int(setParam['NDIAGOP'])))
    if   'NIDEAL'    in setParam: wfh.write(' NIDEAL=%1d,     \n' % (int(setParam['NIDEAL'])))
    if   'NDIFPS'    in setParam: wfh.write(' NDIFPS=%1d, '       % (int(setParam['NDIFPS'])))
    if   'NDIFT'     in setParam: wfh.write(' NDIFT=%1d,      \n' % (int(setParam['NDIFT'])))
    if   'NMESHC'    in setParam: wfh.write(' NMESHC=%1d, '       % (int(setParam['NMESHC'])))
    if   'NPOIDC'    in setParam: wfh.write(' NPOIDC=%1d, '       % (int(setParam['NPOIDC'])))
    if   'SOLPDC'    in setParam: wfh.write(' SOLPDC=%2.2f,   \n' % (float(setParam['SOLPDC'])))
    if   'CPLACE'    in setParam: wfh.write(' CPLACE=%s,      \n' % (str(setParam['CPLACE'])))
    if   'CWIDTH'    in setParam: wfh.write(' CWIDTH=%s,      \n' % (str(setParam['CWIDTH'])))
    if   'NMESHPOL'  in setParam: wfh.write(' NMESHPOL=%4d, '     % (int(setParam['NMESHPOL'])))
    if   'SOLPDPOL'  in setParam: wfh.write(' SOLPDPOL=%2.2f, \n' % (float(setParam['SOLPDPOL'])))
    if   'NTURN'     in setParam: wfh.write(' NTURN=%2d, '        % (int(setParam['NTURN'])))
    if   'NBLC0'     in setParam: wfh.write(' NBLC0=%2d, '        % (int(setParam['NBLC0'])))
    if   'NPPR'      in setParam: wfh.write(' NPPR=%2d,       \n' % (int(setParam['NPPR'])))
    if   'NINMAP'    in setParam: wfh.write(' NINMAP=%2d, '       % (int(setParam['NINMAP'])))
    if   'NINSCA'    in setParam: wfh.write(' NINSCA=%2d,     \n' % (int(setParam['NINSCA'])))
    if   'NSYM'      in setParam: wfh.write(' NSYM=%1d, '         % (int(setParam['NSYM'])))
    if   'NEGP'      in setParam: wfh.write(' NEGP=%1d, '         % (int(setParam['NEGP'])))
    if   'NER'       in setParam: wfh.write(' NER=%1d,        \n' % (int(setParam['NER'])))
    if   'EPSLON'    in setParam: wfh.write(' EPSLON=%6.2E,   \n' % (float(setParam['EPSLON'])))
    if   'ETAEI'     in setParam: wfh.write(' ETAEI=%2.1f, '      % (float(setParam['ETAEI'])))
    if   'RPEOP'     in setParam: wfh.write(' RPEOP=%2.1f, '      % (float(setParam['RPEOP'])))
    if   'RZION'     in setParam: wfh.write(' RZION=%2.1f, '      % (float(setParam['RZION'])))
    if   'GAMMA'     in setParam: wfh.write(' GAMMA=%12.11f,  \n' % (float(setParam['GAMMA'])))
    if   'AT3(1)'    in setParam: wfh.write(' AT3(1)=%2.2f,   \n' % (float(setParam['AT3(1)'])))
    if   'TENSPROF'  in setParam: wfh.write(' TENSPROF=%2.2f, \n' % (float(setParam['TENSPROF'])))
    if   'TENSBND'   in setParam: wfh.write(' TENSBND=%2.2f,  \n' % (float(setParam['TENSBND'])))
    if   'cocos_in'  in setParam: wfh.write(' cocos_in=%2d,   \n' % (int(setParam['cocos_in'])))
    if   'cocos_out' in setParam: wfh.write(' cocos_out=%2d   \n' % (int(setParam['cocos_out'])))

    wfh.write(' &END \n')
    wfh.write('\n')
    wfh.close()

    return setParam


def find_boundary(eqdsk='',setParam={}):
    if eqdsk:
       eqdskflag  = True
       if   type(eqdsk)==str and os.path.isfile(eqdsk.strip()):
                               eqdskdata = cheasefiles.read_eqdsk(fpath=eqdsk.strip())
       elif type(eqdsk)==dict: eqdskdata = eqdsk.copy()
       else:
            eqdskflag = False
    else:
       eqdskflag = False

    if not eqdskflag: raise IOError('FATAL: EQDSK FILE IS NOT PROVIDED. EXIT!')

    asisflag = False; interpflag = False
    if 'boundary_src' in setParam:
       if   setParam['boundary_src'] in [0,'asis']:   asisflag   = True
       elif setParam['boundary_src'] in [1,'interp']: interpflag = True
       else:                                           asisflag   = True

    if   asisflag:
         rbound = eqdskdata['rbound']
         zbound = eqdskdata['zbound']
    elif interpflag:
       rbndtst = int(eqdskdata['RLEN']/(max(eqdskdata['rbound'])-abs(min(eqdskdata['rbound']))))
       zbndtst = int(eqdskdata['ZLEN']/(max(eqdskdata['zbound'])+abs(min(eqdskdata['zbound']))))
       if  rbndtst==1 and zbndtst==1:
           rbound,zbound = magsurf_solvflines(eqdskdata=eqdskdata,psi=0.999,eps=1.0e-16)
       else:
           rbound=npy.zeros(2*len(eqdskdata['rbound'])-1)
           zbound=npy.zeros(2*len(eqdskdata['zbound'])-1)
           rbound[0] = eqdskdata['rbound'][0]
           zbound[0] = eqdskdata['zbound'][0]
           for i in range(1,len(eqdskdata['rbound'])):
               rbound[i]  = eqdskdata['rbound'][i]
               rbound[-i] = eqdskdata['rbound'][i]
               zbound[i]  = eqdskdata['zbound'][i]
               zbound[-i] =-eqdskdata['zbound'][i]

    return rbound,zbound


def plot_chease(cheasefpath,skipfigs=1,**kwargs):
    if   os.path.isdir(cheasefpath):
         srhpath = os.path.join(cheasefpath,'chease*.h5')
         if not os.path.isfile(srhpath):
            srhpath = os.path.join(cheasefpath,'chease*.dat')
            if not os.path.isfile(srhpath):
               srhpath = os.path.join(cheasefpath,'*_CHEASE.dat')
    elif os.path.isfile(cheasefpath):
         srhpath = cheasefpath
         slashloc = findall(cheasefpath,"/")
         if slashloc:
            cheasefpath = cheasefpath[:slashloc[-1]]
         else:
            cheasefpath = "./"
    cheaselist = sorted(glob(srhpath))

    if   'eqdsk' in kwargs:
         srhpath = os.path.join(cheasefpath,kwargs['eqdsk'])
         if not os.path.isfile(srhpath):
            srhpath = os.path.join(kwargs['eqdsk'])
            if not os.path.isfile(srhpath):
               raise IOError('EQDSK FILE NOT FOUND')
    elif 'eqdskfname' in kwargs:
         srhpath = os.path.join(cheasefpath,kwargs['eqdskfname'])
         if not os.path.isfile(srhpath):
            srhpath = os.path.join(kwargs['eqdskfname'])
            if not os.path.isfile(srhpath):
               raise IOError('EQDSK FILE NOT FOUND')
    elif 'eqdskfpath' in kwargs:
         srhpath = os.path.join(cheasefpath,kwargs['eqdskfpath'])
         if not os.path.isfile(srhpath):
            srhpath = os.path.join(kwargs['eqdskfpath'])
            if not os.path.isfile(srhpath):
               raise IOError('EQDSK FILE NOT FOUND')
    elif 'eqdskfpath' not in kwargs:
         srhpath  = os.path.join(cheasefpath,'*_EQDSK')
         if not os.path.isfile(srhpath):
            srhpath = os.path.join(cheasefpath,'EQDSK_iter*')
            if not os.path.isfile(srhpath):
               srhpath = ""
    eqdsklist  = sorted(glob(srhpath))

    if   'profiles' in kwargs:
         srhpath = os.path.join(cheasefpath,kwargs['profiles'])
         if not os.path.isfile(srhpath):
            srhpath = os.path.join(kwargs['profiles'])
            if not os.path.isfile(srhpath):
               raise IOError('PROFILES FILE NOT FOUND')
    elif 'profilesfname' in kwargs:
         srhpath = os.path.join(cheasefpath,kwargs['profilesfname'])
         if not os.path.isfile(srhpath):
            srhpath = os.path.join(kwargs['profilesfname'])
            if not os.path.isfile(srhpath):
               raise IOError('PROFILES FILE NOT FOUND')
    elif 'profilesfpath' in kwargs:
         srhpath = os.path.join(cheasefpath,kwargs['profilesfpath'])
         if not os.path.isfile(srhpath):
            srhpath = os.path.join(kwargs['profilesfpath'])
            if not os.path.isfile(srhpath):
               raise IOError('PROFILES FILE NOT FOUND')
    elif 'profilesfpath' not in kwargs:
         srhpath  = os.path.join(cheasefpath,'*_PROFILES')
         if not os.path.isfile(srhpath):
            srhpath = os.path.join(cheasefpath,'PROFILES_iter*')
            if not os.path.isfile(srhpath):
               srhpath = ""
    profileslist   = sorted(glob(srhpath))

    if   'expeq' in kwargs:
         srhpath = os.path.join(cheasefpath,kwargs['expeq'])
         if not os.path.isfile(srhpath):
            srhpath = os.path.join(kwargs['expeq'])
            if not os.path.isfile(srhpath):
               raise IOError('EXPEQ FILE NOT FOUND')
    elif 'expeqfname' in kwargs:
         srhpath = os.path.join(cheasefpath,kwargs['expeqfname'])
         if not os.path.isfile(srhpath):
            srhpath = os.path.join(kwargs['expeqfname'])
            if not os.path.isfile(srhpath):
               raise IOError('EXPEQ FILE NOT FOUND')
    elif 'expeqfpath' in kwargs:
         srhpath = os.path.join(cheasefpath,kwargs['expeqfpath'])
         if not os.path.isfile(srhpath):
            srhpath = os.path.join(kwargs['expeqfpath'])
            if not os.path.isfile(srhpath):
               raise IOError('EXPEQ FILE NOT FOUND')
    elif 'expeqfpath' not in kwargs:
         srhpath  = os.path.join(cheasefpath,'*_EXPEQ')
         if not os.path.isfile(srhpath):
            srhpath = os.path.join(cheasepath,'EXPEQ_iter*')
            if not os.path.isfile(srhpath):
               srhpath = ""
    expeqlist  = sorted(glob(srhpath))

    if   'exptnz' in kwargs:
         srhpath = os.path.join(cheasefpath,kwargs['exptnz'])
         if not os.path.isfile(srhpath):
            srhpath = os.path.join(kwargs['exptnz'])
            if not os.path.isfile(srhpath):
               raise IOError('EXPTNZ FILE NOT FOUND')
    elif 'exptnzfname' in kwargs:
         srhpath = os.path.join(cheasefpath,kwargs['exptnzfname'])
         if not os.path.isfile(srhpath):
            srhpath = os.path.join(kwargs['exptnzfname'])
            if not os.path.isfile(srhpath):
               raise IOError('EXPTNZ FILE NOT FOUND')
    elif 'exptnzfpath' in kwargs:
         srhpath = os.path.join(cheasefpath,kwargs['exptnzfpath'])
         if not os.path.isfile(srhpath):
            srhpath = os.path.join(kwargs['exptnzfpath'])
            if not os.path.isfile(srhpath):
               raise IOError('EXPTNZ FILE NOT FOUND')
    elif 'exptnzfpath' not in kwargs:
         srhpath  = os.path.join(cheasefpath,'*_EXPTNZ')
         if not os.path.isfile(srhpath):
            srhpath = os.path.join(cheasefpath,'EXPTNZ_iter*')
            if not os.path.isfile(srhpath):
               srhpath = ""
    exptnzlist = sorted(glob(srhpath))

    icounter = 1
    for cheasefid in cheaselist[0::skipfigs+1]:
        print('Plotting CHEASE data in: %s ...' % cheasefid)
        if   cheasefid[13:17] == 'iter':
             caselabel  = cheasefid[13:21]
        elif cheasefid[8:14] in ['KEFITD','MEFITD']:
             caselabel  = cheasefid[8:21]
        else:
             caselabel  = cheasefid

        CHEASEdata = cheasefiles.read_chease(fpath=cheasefid)
        CHEASEdataKeys = CHEASEdata.keys()

        if exptnzlist:
           EXPTNZdata   = cheasefiles.read_exptnz(exptnzlist[cheaselist.index(cheasefid)],eqdsk=eqdsklist[0])
        if expeqlist:
           EXPEQdata    = cheasefiles.read_expeq(expeqlist[cheaselist.index(cheasefid)],eqdsk=eqdsklist[0])
        if eqdsklist:
           EQDSKdata    = cheasefiles.read_eqdsk(eqdsklist[0])
        if profileslist:
           PROFILESdata = cheasefiles.read_profiles(profileslist[0])

        EDENfig = plt.figure("Electron Density")
        plt.plot(CHEASEdata['rhopsi'],CHEASEdata['ne'],linestyle='-',label='CHESAE-'+caselabel[-6:-3])
        if exptnzlist:
           plt.plot(EXPTNZdata['rhopsi'],EXPTNZdata['ne'],linestyle=':',label='EXPTNZ-'+caselabel[-6:-3])
        plt.title('Electron Density Profiles')
        plt.xlabel('$\\rho_{\psi_N}$')
        plt.ylabel('$n_e$')
        plt.legend()
    
        GDNEfig = plt.figure("Electron Density Gradient")
        plt.plot(CHEASEdata['rhopsi'],CHEASEdata['nePrime'],label=caselabel)
        plt.title('Electron Density Gradient Profiles')
        plt.xlabel('$\\rho_{\psi_N}$')
        plt.ylabel('$\\nabla{n_e}$')
        plt.legend()
    
        ETMPfig = plt.figure("Electron Temperature")
        plt.plot(CHEASEdata['rhopsi'],CHEASEdata['Te'],linestyle='-',label='CHESAE-'+caselabel[-6:-3])
        if exptnzlist:
           plt.plot(EXPTNZdata['rhopsi'],EXPTNZdata['Te'],linestyle=':',label='EXPTNZ-'+caselabel[-6:-3])
        plt.title('Electron Temperature Profiles')
        plt.xlabel('$\\rho_{\psi_N}$')
        plt.ylabel('$T_e$')
        plt.legend()

        GDTEfig = plt.figure("Electron Temperature Gradient")
        plt.plot(CHEASEdata['rhopsi'],CHEASEdata['TePrime'],label=caselabel)
        plt.title('Electron Temperature Gradient Profiles')
        plt.xlabel('$\\rho_{\psi_N}$')
        plt.ylabel('$\\nabla{T_e}$')
        plt.legend()
    
        SFTYfig = plt.figure("Safety Factor (q)")
        plt.plot(CHEASEdata['rhopsi'],CHEASEdata['q'],linestyle='-',label='CHESAE-'+caselabel[-6:-3])
        if expeqlist and 'q' in EXPEQdata:
           plt.plot(EXPEQdata['rhopsi'],EXPEQdata['q'],linestyle='--',label='EXPEQ-'+caselabel[-6:-3])
        if eqdsklist:
           plt.plot(EQDSKdata['rhopsi'], EQDSKdata['q'], linestyle=':',label='EQDSK')
        plt.title("Safety Factor Profiles")
        plt.xlabel('$\\rho_{\psi_N}$')
        plt.ylabel("q")
        plt.legend()

        TPPRfig = plt.figure("Plasma Pressure")
        plt.plot(CHEASEdata['rhopsi'],CHEASEdata['pressure'],  linestyle='solid', label='CHESAE-'+caselabel[-6:-3])
        plt.plot(EQDSKdata['rhopsi'], EQDSKdata['pressure'],   linestyle='dotted', label='EQDSK')
        if profileslist and 'pressure' in PROFILESdata:
           plt.plot(PROFILESdata['rhopsi'],PROFILESdata['pressure'],linestyle=(0,(5,1)),label='PROFILES')
        if exptnzlist and 'pressure' in EXPTNZdata:
           plt.plot(EXPTNZdata['rhopsi'],EXPTNZdata['pressure'],  linestyle='dashdot',label='EXPTNZ-'+caselabel[-6:-3])
        if expeqlist and 'pressure' in EXPEQdata:
           EXPEQdata['pressure'] = EXPEQdata['pressure']*CHEASEdata['B0EXP']**2/mu0
           plt.plot(EXPEQdata['rhopsi'],EXPEQdata['pressure'],linestyle='dashed',label='EXPEQ-'+caselabel[-6:-3])
        plt.title('Plasma Pressure Profiles')
        plt.xlabel('$\\rho_{psi_N}$')
        plt.ylabel('$P$')
        plt.legend()
    
        PPRMfig = plt.figure("P'")
        plt.plot(CHEASEdata['rhopsi'],CHEASEdata['pprime'], linestyle='solid', label='CHESAE-'+caselabel[-6:-3])
        if eqdsklist:
           plt.plot(EQDSKdata['rhopsi'], EQDSKdata['pprime'],  linestyle='dotted', label='EQDSK')
        if profileslist and 'pprime' in PROFILESdata:
           plt.plot(PROFILESdata['rhopsi'],PROFILESdata['pprime'],linestyle=(0,(5,1)),label='PROFILES')
        if exptnzlist and 'pprime' in EXPTNZdata:
           plt.plot(EXPTNZdata['rhopsi'],EXPTNZdata['pprime'],linestyle=(0,(5,1)),label='EXPTNZ')
        if expeqlist and 'pprime' in EXPEQdata:
           EXPEQdata['pprime'] = EXPEQdata['pprime']*CHEASEdata['B0EXP']/mu0/CHEASEdata['R0EXP']**2
           plt.plot(EXPEQdata['rhopsi'],EXPEQdata['pprime'],linestyle='dashed',label='EXPEQ-'+caselabel[-6:-3])
        plt.title("P' Profiles")
        plt.xlabel('$\\rho_{\psi_N}$')
        plt.ylabel("P'")
        plt.legend()
    
        TTPMfig = plt.figure("FF'")
        plt.plot(CHEASEdata['rhopsi'],CHEASEdata['ffprime'],linestyle='-',label='CHESAE-'+caselabel[-6:-3])
        if eqdsklist:
           plt.plot( EQDSKdata['rhopsi'], EQDSKdata['ffprime'],linestyle=':',label='EQDSK')
        if expeqlist and 'ffprime' in EXPEQdata:
           EXPEQdata['ffprime'] = EXPEQdata['ffprime']*CHEASEdata['B0EXP']
           plt.plot(EXPEQdata['rhopsi'],EXPEQdata['ffprime'],linestyle='--',label='EXPEQ-'+caselabel[-6:-3])
        plt.title("TT' Profiles")
        plt.xlabel('$\\rho_{\psi_N}$')
        plt.ylabel("TT'")
        plt.legend()
    
        ISTRfig = plt.figure("I*")
        plt.plot(CHEASEdata['rhopsi'],CHEASEdata['Istr'],linestyle='-',label='CHESAE-'+caselabel[-6:-3])
        if 'EXPEQdata' in locals() and 'Istr' in EXPEQdata:
           EXPEQdata['Istr'] = EXPEQdata['Istr']*CHEASEdata['R0EXP']*CHEASEdata['B0EXP']/mu0
           plt.plot(EXPEQdata['rhopsi'],EXPEQdata['Istr'],linestyle='--',label='EXPEQ-'+caselabel[-6:-3])
        plt.title("I* Profiles")
        plt.xlabel('$\\rho_{\psi_N}$')
        plt.ylabel("I*")
        plt.legend()
    
        ICRTfig = plt.figure("Parallel Current")
        plt.plot(CHEASEdata['rhopsi'],CHEASEdata['Iprl'],linestyle='-',label='CHESAE-'+caselabel[-6:-3])
        if expeqlist and 'Iprl' in EXPEQdata:
           EXPEQdata['Iprl'] = EXPEQdata['Iprl']*CHEASEdata['R0EXP']*CHEASEdata['B0EXP']/mu0
           plt.plot(EXPEQdata['rhopsi'],EXPEQdata['Iprl'],linestyle='--',label='EXPEQ-'+caselabel[-6:-3])
        plt.title('Parallel Current Profiles')
        plt.xlabel('$\\rho_{\psi_N}$')
        plt.ylabel('$I_{||}$')
        plt.legend()
    
        JCRTfig = plt.figure("Parallel Current Density")
        plt.plot(CHEASEdata['rhopsi'],CHEASEdata['Jprl'],linestyle='-',label='CHESAE-'+caselabel[-6:-3])
        if expeqlist and 'Jprl' in EXPEQdata:
           EXPEQdata['Jprl'] = EXPEQdata['Jprl']*CHEASEdata['B0EXP']/CHEASEdata['B0EXP']/mu0
           plt.plot(EXPEQdata['rhopsi'],EXPEQdata['Jprl'],linestyle='--',label='EXPEQ-'+caselabel[-6:-3])
        plt.title('Parallel Current Density Profiles')
        plt.xlabel('$\\rho_{\psi_N}$')
        plt.ylabel('$J_{||}$')
        plt.legend()
    
        SCRTfig = plt.figure("Bootstrap Currents")
        plt.plot(CHEASEdata['rhopsi'],CHEASEdata['Jbs'],label='$J_{BS}$-'+caselabel)
        plt.title('Bootstrap Current Density Profiles')
        plt.xlabel('$\\rho_{\psi_N}$')
        plt.ylabel('$J_{BS}$')
        plt.legend()
    
        (CHEASEdata['PSIN2D'],CHEASEdata['CHIN2D']) = npy.meshgrid(CHEASEdata['PSIN'],CHEASEdata['CHIN'])
        BF2Dfig = plt.figure("Magnetic Field, B($\psi$,$\chi$)")
        plt.contour(CHEASEdata['CHIN2D'],CHEASEdata['PSIN2D'],CHEASEdata['B'])
        plt.title('Magnetic Field Profiles')
        plt.xlabel('$\chi$')
        plt.ylabel('$\psi$')
    
        JPHIfig = plt.figure("Toroidal Current")
        plt.contour(CHEASEdata['CHIN2D'],CHEASEdata['PSIN2D'],CHEASEdata['Jphi'],cmap=plt.cm.hot)
        plt.title('Toroidal Current Profiles')
        plt.xlabel('$\\rho_{\psi_N}$')
        plt.ylabel('$J_{\phi}$')
    
        BFRZfig = plt.figure("Magnetic Field, B(R,Z}")
        plt.contour(CHEASEdata['R'],CHEASEdata['Z'],CHEASEdata['B'])
        plt.title('Magnetic Field Profiles, B(R,Z}')
        plt.xlabel('$R$')
        plt.ylabel('$Z$')

        PSRZfig = plt.figure("Magnetic Poloidal Boundary Surface")
        plt.plot(CHEASEdata['rbound'],CHEASEdata['zbound'],label=caselabel)
        plt.title('Magnetic Poloidal Boundary Surface')
        plt.xlabel('$R$')
        plt.ylabel('$Z$')
        plt.legend()

        del(CHEASEdata)

        chsfigs = PdfPages('cheaseresults.pdf')
        chsfigs.savefig(EDENfig)
        chsfigs.savefig(GDNEfig)
        chsfigs.savefig(ETMPfig)
        chsfigs.savefig(GDTEfig)
        chsfigs.savefig(SFTYfig)
        chsfigs.savefig(TPPRfig)
        chsfigs.savefig(PPRMfig)
        chsfigs.savefig(TTPMfig)
        chsfigs.savefig(ISTRfig)
        chsfigs.savefig(ICRTfig)
        chsfigs.savefig(JPHIfig)
        chsfigs.savefig(JCRTfig)
        chsfigs.savefig(SCRTfig)
        chsfigs.savefig(BF2Dfig)
        chsfigs.savefig(BFRZfig)
        chsfigs.savefig(PSRZfig)
        chsfigs.close()

    return 1

def remove_input_files():
    if glob('./*_PROFILES'): os.system('rm ./*_PROFILES')
    if glob('./*_EXPTNZ'):   os.system('rm ./*_EXPTNZ')
    if glob('./*_CHEASE'):   os.system('rm ./*_CHEASE')
    if glob('./*_ITERDB'):   os.system('rm ./*_ITERDB')
    if glob('./*_EQDSK'):    os.system('rm ./*_EQDSK')
    if glob('./*_EXPEQ'):    os.system('rm ./*_EXPEQ')
    return 1

def remove_output_files():
    if glob('./NGA'):              os.system('rm NGA')
    if glob('./NDES'):             os.system('rm NDES')
    if glob('./*OUT*'):            os.system('rm *OUT*')
    if glob('./*.pdf'):            os.system('rm *.pdf')
    if glob('./EXPEQ'):            os.system('rm EXPEQ')
    if glob('./EXPTNZ'):           os.system('rm EXPTNZ')
    if glob('./ogyropsi*'):        os.system('rm ogyropsi*')
    if glob('./chease.log'):       os.system('rm chease.log')
    if glob('./EQDSK_iter*'):      os.system('rm EQDSK_iter*')
    if glob('./chease_iter*'):     os.system('rm chease_iter*')
    if glob('./EXPEQ_EQDSK*'):     os.system('rm EXPEQ_EQDSK*')
    if glob('./EXPEQ_EXPEQ*'):     os.system('rm EXPEQ_EXPEQ*')
    if glob('./EXPEQ_iter*.IN'):   os.system('rm EXPEQ_iter*.IN')
    if glob('./EXPTNZ_iter*.IN'):  os.system('rm EXPTNZ_iter*.IN')
    if glob('./chease_namelist*'): os.system('rm chease_namelist*')
    return 1

def init_chease_inputs(srcVals={},namelistVals={},importedVals={}):
    '''
    "boundary_src  = n" where n in geometry options
    Boundary Surface options:
    n = 0 or 'asis' --- Take the boundary surface 'as is'
    n = 1 or 'interp' - Interpolate the given boundary surface

    "current_src   = n" where n in geometry options
    Geometry Source options:
    n = 0 or 'chease'
    n = 1 or 'eqdsk'
    n = 2 or 'expeq'
    n = 7 or 'imported'

    "rhomesh_src   = n" where n in meshgrid options
    Meshgrid Source options:
    n = 0 or 'chease'
    n = 1 or 'eqdsk'

    "pressure_src  = n" where n in profiles options
    "eprofiles_src = n" where n in profiles options
    "iprofiles_src = n" where n in profiles options
    Profiles Source options:
    n = 0 or 'chease'
    n = 1 or 'eqdsk'
    n = 2 or 'expeq'
    n = 3 or 'exptnz'
    n = 4 or 'profiles'
    n = 5 or 'iterdb'
    n = 7 or 'imported'
    '''

    if 'boundary_src'  in srcVals: boundary_src  = srcVals['boundary_src']
    else:                          boundary_src  = 0

    if 'rhomesh_src'   in srcVals: rhomesh_src   = srcVals['rhomesh_src']
    else:                          rhomesh_src   = None

    if 'current_src'   in srcVals: current_src   = srcVals['current_src']
    else:                          current_src   = None

    if 'pressure_src'  in srcVals: pressure_src  = srcVals['pressure_src']
    else:                          pressure_src  = None

    if 'eprofiles_src' in srcVals: eprofiles_src = srcVals['eprofiles_src']
    else:                          eprofiles_src = None

    if 'iprofiles_src' in srcVals: iprofiles_src = srcVals['iprofiles_src']
    else:                          iprofiles_src = None

    eqdskrequired    = False
    expeqrequired    = False
    cheaserequired   = False
    exptnzrequired   = False
    iterdbrequired   = False
    profilesrequired = False

    if   rhomesh_src  in [0,'chease']:    cheaserequired   = True
    elif rhomesh_src  in [1,'eqdsk']:     eqdskrequired    = True
    elif rhomesh_src  in [7,'imported']:  importedrequired = True

    if   current_src  in [0,'chease']:    cheaserequired   = True
    elif current_src  in [1,'eqdsk']:     eqdskrequired    = True
    elif current_src  in [2,'expeq']:     expeqrequired    = True
    elif current_src  in [7,'imported']:  importedrequired = True

    if   pressure_src in [0,'chease']:    cheaserequired   = True
    elif pressure_src in [1,'eqdsk']:     eqdskrequired    = True
    elif pressure_src in [2,'expeq']:     expeqrequired    = True
    elif pressure_src in [3,'exptnz']:    exptnzrequired   = True
    elif pressure_src in [4,'profiles']:  profilesrequired = True
    elif pressure_src in [5,'iterdb']:    iterdbrequired   = True
    elif pressure_src in [7,'imported']:  importedrequired = True

    if   eprofiles_src in [0,'chease']:   cheaserequired   = True
    elif eprofiles_src in [3,'exptnz']:   exptnzrequired   = True
    elif eprofiles_src in [4,'profiles']: profilesrequired = True
    elif eprofiles_src in [5,'iterdb']:   iterdbrequired   = True
    elif eprofiles_src in [7,'imported']: importedrequired = True

    if   iprofiles_src in [0,'chease']:   cheaserequired   = True
    elif iprofiles_src in [3,'exptnz']:   exptnzrequired   = True
    elif iprofiles_src in [4,'profiles']: profilesrequired = True
    elif iprofiles_src in [5,'iterdb']:   iterdbrequired   = True
    elif iprofiles_src in [7,'imported']: importedrequired = True
    
    namelistParam = {}
    if namelistVals:
       namelistValsKeys = list(namelistVals.keys())
       for ikey in namelistValsKeys:
           namelistParam[ikey] = namelistVals[ikey]

    if 'inputpath'   in srcVals: inputpath = srcVals['inputpath']
    else:                      inputpath = "./"

    if 'gfname'      in srcVals:
       eqdskfname     = srcVals['gfname']
       if os.path.isfile(os.path.join(inputpath,eqdskfname)):    EQDSKexist = True
    if 'pfname'      in srcVals:
       profilesfname  = srcVals['pfname']
       if os.path.isfile(os.path.join(inputpath,profilesfname)): PROFILESexist = True
    if 'iterdbfname' in srcVals:
       iterdbfname    = srcVals['iterdbfname']
       if os.path.isfile(os.path.join(inputpath,iterdbfname)):   ITERDBexist = True
    if 'cheasefname' in srcVals:
       cheasefname    = srcVals['cheasefname']
       if os.path.isfile(os.path.join(inputpath,cheasefname)):   CHEASEexist = True
    if 'exptnzfname' in srcVals:
       exptnzfname    = srcVals['exptnzfname']
       if os.path.isfile(os.path.join(inputpath,exptnzfname)):   EXPTNZexist = True
    if 'expeqfname'  in srcVals:
       expeqfname     = srcVals['expeqfname']
       if os.path.isfile(os.path.join(inputpath,expeqfname)):    EXPEQexist = True

    if eprofiles_src == None:
       if   EXPTNZexist:   eprofiles_src = 3; exptnzrequired    = True
       elif PROFILESexist: eprofiles_src = 4; proflilesrequired = True
       elif ITERDBexist:   eprofiles_src = 5; iterdbrequired    = True
       elif CHEASEexist:   eprofiles_src = 0; cheaserequired    = True

    if iprofiles_src == None:
       if   EXPTNZexist:   iprofiles_src = 3; exptnzrequired    = True
       elif PROFILESexist: iprofiles_src = 4; profilesrequired  = True
       elif ITERDBexist:   iprofiles_src = 5; iterdbrequired    = True
       elif CHEASEexist:   iprofiles_src = 0; cheaserequired    = True

    if pressure_src == None:
       if   EQDSKexist:    pressure_src = 1;  eqdskrequired     = True
       elif EXPEQexist:    pressure_src = 2;  expeqrequired     = True
       elif EXPTNZexist:   pressure_src = 3;  exptnzrequired    = True
       elif PROFILESexist: pressure_src = 4;  profilesrequired  = True
       elif ITERDBexist:   pressure_src = 5;  iterdbrequired    = True
       elif CHEASEexist:   pressure_src = 0;  cheaserequired    = True

    if current_src == None:
       if   EQDSKexist:    current_src = 1;   eqdskrequired     = True
       elif EXPEQexist:    current_src = 2;   expeqrequired     = True
       elif CHEASEexist:   current_src = 0;   cheaserequired    = True

    namelist = create_namelist(setParam=namelistParam)

    if   'NFUNRHO'  in namelist: rhomesh_type = int(namelist['NFUNRHO'])
    elif 'NRHOMESH' in namelist: rhomesh_type = int(namelist['NRHOMESH'])
    else:                        rhomesh_type = 0

    if   'NPPFUN'   in namelist: pressure_type = int(namelist['NPPFUN'])
    else:                        pressure_type = 8

    if   'NSTTP'    in namelist: current_type  = int(namelist['NSTTP'])
    else:                        current_type  = 1

    if rhomesh_src in [7,'imported']:
       if   'rhopsi' not in importedVals and 'rhotor' not in importedVals:
            raise ValueError('importedVals MUST contain rhopsi and rhotor with rhomesh_src = 7 or "imported"')
    if eprofiles_src in [7,'imported']:
       if   'Te' not in importedVals and 'ne' not in importedVals:
            raise ValueError('importedVals MUST contain Te and ne with eprofiles_src = 7 or "imported"')
    if iprofiles_src in [7,'imported']:
       if   'Ti' not in importedVals and 'ni' not in importedVals and 'Zeff' not in importedVals:
            raise ValueError('importedVals MUST contain Ti, ni, and Zeff with iprofiles_src = 7 or "imported"')
    if pressure_src in [7,'imported']:
       if   pressure_type in [8,'pressure'] and 'pressure' not in importedVals:
            raise ValueError('importedVals MUST contain pressure with pressure_src = 7 or "imported"')
       elif pressure_type in [4,'pprime']   and 'pprime'   not in importedVals:
            raise ValueError('importedVals MUST contain pprime with pressure_src = 7 or "imported"')
    if current_src in [7,'imported']:
       if   current_type in [1,'ffprime'] and 'ffprime' not in importedVals:
            raise ValueError('importedVals MUST contain ffprime current_src = 7 or "imported"')
       elif current_type in [2,'istr'] and 'Istr' not in importedVals:
            raise ValueError('importedVals MUST contain Istr current_src = 7 or "imported"')
       elif current_type in [3,'iprl'] and 'Iprl' not in importedVals:
            raise ValueError('importedVals MUST contain Iprl current_src = 7 or "imported"')
       elif current_type in [4,'jprl'] and 'Jprl' not in importedVals:
            raise ValueError('importedVals MUST contain Jprl current_src = 7 or "imported"')
       elif current_type in [5,'q'] and 'q' not in importedVals:
            raise ValueError('importedVals MUST contain q current_src = 7 or "imported"')

    it = max(0,len(glob('./chease_*.dat')) - 1)

    if   'gfname' in srcVals and os.path.isfile(os.path.join(inputpath,eqdskfname)):
         if inputpath != os.getcwd(): os.system('cp %s/%s .' % (inputpath,eqdskfname))
         eqdskfpath  = eqdskfname
    elif os.path.isfile("EQDSK_iter%03d" % it):
         eqdskfpath =   "EQDSK_iter%03d" % it
    elif eqdskrequired:
         raise IOError('EQDSK FILE IS MISSING!')

    if   'cheasefname' in srcVals and os.path.isfile(os.path.join(inputpath,cheasefname)):
         if inputpath != os.getcwd(): os.system('cp %s/%s .' % (inputpath,cheasefname))
         cheasefpath = cheasefname
    elif os.path.isfile("chease_iter%03d" % it):
         cheasefpath =  "chease_iter%03d" % it
    elif cheaserequired:
         raise IOError('CHEASE FILE IS MISSING!')

    if   'expeqfname' in srcVals and os.path.isfile(os.path.join(inputpath,expeqfname)):
         if inputpath != os.getcwd(): os.system('cp %s/%s .' % (inputpath,expeqfname))
         expeqfpath = expeqfname
    elif os.path.isfile("expeq_iter%03d" % it):
         expeqfpath =   "expeq_iter%03d" % it
    elif expeqrequired:
         raise IOError('EXPEQ FILE IS MISSING!')

    if   'exptnzfname' in srcVals and os.path.isfile(os.path.join(inputpath,exptnzfname)):
         if inputpath != os.getcwd(): os.system('cp %s/%s .' % (inputpath,exptnzfname))
         exptnzfpath = exptnzfname
    elif os.path.isfile("exptnz_iter%03d" % it):
         exptnzfpath =  "exptnz_iter%03d" % it
    elif exptnzrequired:
         raise IOError('EXPTNZ FILE IS MISSING!')

    if   'pfname' in srcVals and os.path.isfile(os.path.join(inputpath,profilesfname)):
         if inputpath != os.getcwd(): os.system('cp %s/%s .' % (inputpath,profilesfname))
         profilesfpath = profilesfname
    elif os.path.isfile("profiles_iter%03d" % it):
         profilesfpath ="profiles_iter%03d" % it
    elif profilesrequired:
         raise IOError('PROFILES FILE IS MISSING!')

    if   'iterdbfname' in srcVals and os.path.isfile(os.path.join(inputpath,iterdbfname)):
         if inputpath != os.getcwd(): os.system('cp %s/%s .' % (inputpath,iterdbfname))
         iterdbfpath = iterdbfname
    elif os.path.isfile("iterdb_iter%03d" % it):
         iterdbfpath =  "iterdb_iter%03d" % it
    elif iterdbrequired:
         raise IOError('ITERDB FILE IS MISSING!')


    exptnzParam = {}
    if int(namelist['NBSEXPQ']) != 0:
       exptnzParam['nrhomesh']  = [rhomesh_type,rhomesh_src]
       exptnzParam['eprofiles'] = eprofiles_src
       exptnzParam['iprofiles'] = iprofiles_src

       if   rhomesh_src in [0,'chease']:
            if   eprofiles_src in [0,'chease']   and iprofiles_src in [0,'chease']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,imported=importedVals)
            elif eprofiles_src in [0,'chease']   and iprofiles_src in [3,'exptnz']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,exptnz=exptnzfpath,imported=importedVals)
            elif eprofiles_src in [0,'chease']   and iprofiles_src in [4,'profiles']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,profiles=profilesfpath,imported=importedVals)
            elif eprofiles_src in [0,'chease']   and iprofiles_src in [5,'iterdb']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,iterdb=iterdbfpath,imported=importedVals)
            elif eprofiles_src in [0,'chease']   and iprofiles_src in [7,'imported']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,imported=importedVals)
            elif eprofiles_src in [3,'exptnz']   and iprofiles_src in [0,'chease']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,exptnz=exptnzfpath,imported=importedVals)
            elif eprofiles_src in [3,'exptnz']   and iprofiles_src in [3,'exptnz']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,exptnz=exptnzfpath,imported=importedVals)
            elif eprofiles_src in [3,'exptnz']   and iprofiles_src in [4,'profiles']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,exptnz=exptnzfpath,profiles=profilesfpath,imported=importedVals)
            elif eprofiles_src in [3,'exptnz']   and iprofiles_src in [5,'iterdb']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,exptnz=exptnzfpath,iterdb=iterdbfpath,imported=importedVals)
            elif eprofiles_src in [3,'exptnz']   and iprofiles_src in [7,'imported']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,exptnz=exptnzfpath,imported=importedVals)
            elif eprofiles_src in [4,'profiles'] and iprofiles_src in [0,'chease']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,profiles=profilesfpath,imported=importedVals)
            elif eprofiles_src in [4,'profiles'] and iprofiles_src in [3,'exptnz']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,profiles=profilesfpath,exptnz=exptnzfpath,imported=importedVals)
            elif eprofiles_src in [4,'profiles'] and iprofiles_src in [4,'profiles']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,profiles=profilesfpath,imported=importedVals)
            elif eprofiles_src in [4,'profiles'] and iprofiles_src in [5,'iterdb']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,profiles=profilesfpath,iterdb=iterdbfpath,imported=importedVals)
            elif eprofiles_src in [4,'profiles'] and iprofiles_src in [7,'imported']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,profiles=profilesfpath,imported=importedVals)
            elif eprofiles_src in [5,'iterdb']   and iprofiles_src in [0,'chease']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,iterdb=iterdbfpath,imported=importedVals)
            elif eprofiles_src in [5,'iterdb']   and iprofiles_src in [3,'exptnz']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,iterdb=iterdbfpath,exptnz=exptnzfpath,imported=importedVals)
            elif eprofiles_src in [5,'iterdb']   and iprofiles_src in [4,'profiles']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,iterdb=iterdbfpath,profiles=profilesfpath,imported=importedVals)
            elif eprofiles_src in [5,'iterdb']   and iprofiles_src in [5,'iterdb']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,iterdb=iterdbfpath,imported=importedVals)
            elif eprofiles_src in [5,'iterdb']   and iprofiles_src in [7,'imported']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,iterdb=iterdbfpath,imported=importedVals)
            elif eprofiles_src in [7,'imported'] and iprofiles_src in [0,'chease']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,imported=importedVals)
            elif eprofiles_src in [7,'imported'] and iprofiles_src in [3,'exptnz']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,exptnz=exptnzfpath,imported=importedVals)
            elif eprofiles_src in [7,'imported'] and iprofiles_src in [4,'profiles']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,profiles=profilesfpath,imported=importedVals)
            elif eprofiles_src in [7,'imported'] and iprofiles_src in [5,'iterdb']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,iterdb=iterdbfpath,imported=importedVals)
            elif eprofiles_src in [7,'imported'] and iprofiles_src in [7,'imported']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,imported=importedVals)
       elif rhomesh_src in [1,'eqdsk']:
            if   eprofiles_src in [0,'chease']   and iprofiles_src in [0,'chease']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,chease=cheasefpath,imported=importedVals)
            elif eprofiles_src in [0,'chease']   and iprofiles_src in [3,'exptnz']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,chease=cheasefpath,exptnz=exptnzfpath,imported=importedVals)
            elif eprofiles_src in [0,'chease']   and iprofiles_src in [4,'profiles']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,chease=cheasefpath,profiles=profilesfpath,imported=importedVals)
            elif eprofiles_src in [0,'chease']   and iprofiles_src in [5,'iterdb']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,chease=cheasefpath,iterdb=iterdbfpath,imported=importedVals)
            elif eprofiles_src in [0,'chease']   and iprofiles_src in [7,'imported']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,chease=cheasefpath,imported=importedVals)
            elif eprofiles_src in [3,'exptnz']   and iprofiles_src in [0,'chease']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,chease=cheasefpath,exptnz=exptnzfpath,imported=importedVals)
            elif eprofiles_src in [3,'exptnz']   and iprofiles_src in [3,'exptnz']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,exptnz=exptnzfpath,imported=importedVals)
            elif eprofiles_src in [3,'exptnz']   and iprofiles_src in [4,'profiles']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,exptnz=exptnzfpath,profiles=profilesfpath,imported=importedVals)
            elif eprofiles_src in [3,'exptnz']   and iprofiles_src in [5,'iterdb']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,exptnz=exptnzfpath,iterdb=iterdbfpath,imported=importedVals)
            elif eprofiles_src in [3,'exptnz']   and iprofiles_src in [7,'imported']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,exptnz=exptnzfpath,imported=importedVals)
            elif eprofiles_src in [4,'profiles'] and iprofiles_src in [0,'chease']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,chease=cheasefpath,profiles=profilesfpath,imported=importedVals)
            elif eprofiles_src in [4,'profiles'] and iprofiles_src in [3,'exptnz']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,profiles=profilesfpath,exptnz=exptnzfpath,imported=importedVals)
            elif eprofiles_src in [4,'profiles'] and iprofiles_src in [4,'profiles']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,profiles=profilesfpath,imported=importedVals)
            elif eprofiles_src in [4,'profiles'] and iprofiles_src in [5,'iterdb']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,profiles=profilesfpath,iterdb=iterdbfpath,imported=importedVals)
            elif eprofiles_src in [4,'profiles'] and iprofiles_src in [7,'imported']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,profiles=profilesfpath,imported=importedVals)
            elif eprofiles_src in [5,'iterdb']   and iprofiles_src in [0,'chease']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,chease=cheasefpath,iterdb=iterdbfpath,imported=importedVals)
            elif eprofiles_src in [5,'iterdb']   and iprofiles_src in [3,'exptnz']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,iterdb=iterdbfpath,exptnz=exptnzfpath,imported=importedVals)
            elif eprofiles_src in [5,'iterdb']   and iprofiles_src in [4,'profiles']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,iterdb=iterdbfpath,profiles=profilesfpath,imported=importedVals)
            elif eprofiles_src in [5,'iterdb']   and iprofiles_src in [5,'iterdb']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,iterdb=iterdbfpath,imported=importedVals)
            elif eprofiles_src in [5,'iterdb']   and iprofiles_src in [7,'imported']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,iterdb=iterdbfpath,imported=importedVals)
            elif eprofiles_src in [7,'imported'] and iprofiles_src in [0,'chease']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,chease=cheasefpath,imported=importedVals)
            elif eprofiles_src in [7,'imported'] and iprofiles_src in [3,'exptnz']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,exptnz=exptnzfpath,imported=importedVals)
            elif eprofiles_src in [7,'imported'] and iprofiles_src in [4,'profiles']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,profiles=profilesfpath,imported=importedVals)
            elif eprofiles_src in [7,'imported'] and iprofiles_src in [5,'iterdb']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,iterdb=iterdbfpath,imported=importedVals)
            elif eprofiles_src in [7,'imported'] and iprofiles_src in [7,'imported']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,imported=importedVals)
       elif rhomesh_src in [7,'imported',None]:
            if   eprofiles_src in [0,'chease']   and iprofiles_src in [0,'chease']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,imported=importedVals)
            elif eprofiles_src in [0,'chease']   and iprofiles_src in [3,'exptnz']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,exptnz=exptnzfpath,imported=importedVals)
            elif eprofiles_src in [0,'chease']   and iprofiles_src in [4,'profiles']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,profiles=profilesfpath,imported=importedVals)
            elif eprofiles_src in [0,'chease']   and iprofiles_src in [5,'iterdb']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,iterdb=iterdbfpath,imported=importedVals)
            elif eprofiles_src in [0,'chease']   and iprofiles_src in [7,'imported']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,imported=importedVals)
            elif eprofiles_src in [3,'exptnz']   and iprofiles_src in [0,'chease']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,exptnz=exptnzfpath,imported=importedVals)
            elif eprofiles_src in [3,'exptnz']   and iprofiles_src in [3,'exptnz']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,exptnz=exptnzfpath,imported=importedVals)
            elif eprofiles_src in [3,'exptnz']   and iprofiles_src in [4,'profiles']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,exptnz=exptnzfpath,profiles=profilesfpath,imported=importedVals)
            elif eprofiles_src in [3,'exptnz']   and iprofiles_src in [5,'iterdb']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,exptnz=exptnzfpath,iterdb=iterdbfpath,imported=importedVals)
            elif eprofiles_src in [3,'exptnz']   and iprofiles_src in [7,'imported']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,exptnz=exptnzfpath,imported=importedVals)
            elif eprofiles_src in [4,'profiles'] and iprofiles_src in [0,'chease']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,profiles=profilesfpath,imported=importedVals)
            elif eprofiles_src in [4,'profiles'] and iprofiles_src in [3,'exptnz']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,profiles=profilesfpath,exptnz=exptnzfpath,imported=importedVals)
            elif eprofiles_src in [4,'profiles'] and iprofiles_src in [4,'profiles']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,profiles=profilesfpath,imported=importedVals)
            elif eprofiles_src in [4,'profiles'] and iprofiles_src in [5,'iterdb']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,profiles=profilesfpath,iterdb=iterdbfpath,imported=importedVals)
            elif eprofiles_src in [4,'profiles'] and iprofiles_src in [7,'imported']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,profiles=profilesfpath,imported=importedVals)
            elif eprofiles_src in [5,'iterdb']   and iprofiles_src in [0,'chease']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,iterdb=iterdbfpath,imported=importedVals)
            elif eprofiles_src in [5,'iterdb']   and iprofiles_src in [3,'exptnz']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,iterdb=iterdbfpath,exptnz=exptnzfpath,imported=importedVals)
            elif eprofiles_src in [5,'iterdb']   and iprofiles_src in [4,'profiles']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,iterdb=iterdbfpath,profiles=profilesfpath,imported=importedVals)
            elif eprofiles_src in [5,'iterdb']   and iprofiles_src in [5,'iterdb']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,iterdb=iterdbfpath,imported=importedVals)
            elif eprofiles_src in [5,'iterdb']   and iprofiles_src in [7,'imported']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,iterdb=iterdbfpath,imported=importedVals)
            elif eprofiles_src in [7,'imported'] and iprofiles_src in [0,'chease']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,chease=cheasefpath,imported=importedVals)
            elif eprofiles_src in [7,'imported'] and iprofiles_src in [3,'exptnz']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,exptnz=exptnzfpath,imported=importedVals)
            elif eprofiles_src in [7,'imported'] and iprofiles_src in [4,'profiles']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,profiles=profilesfpath,imported=importedVals)
            elif eprofiles_src in [7,'imported'] and iprofiles_src in [5,'iterdb']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,iterdb=iterdbfpath,imported=importedVals)
            elif eprofiles_src in [7,'imported'] and iprofiles_src in [7,'imported']:
                 cheasefiles.write_exptnz(setParam=exptnzParam,imported=importedVals)

    eqdskdata = cheasefiles.read_eqdsk(fpath=eqdskfpath)
    if 'R0EXP' in namelistVals:
        R0EXP = namelistVals['R0EXP']
    else:
        R0EXP = abs(eqdskdata['RCTR'])
    if 'B0EXP' in namelistVals:
        B0EXP = namelistVals['B0EXP']
    else:
        B0EXP = abs(eqdskdata['BCTR'])
    if 'ITEXP' in namelistVals:
        ITEXP = namelistVals['ITEXP']
    else:
        ITEXP = abs(eqdskdata['CURNT'])

    expeqParam = {}
    if   int(namelist['NEQDSK']) == 1:
         print('Reading from EQDSK file.')
         os.system('cp *_EQDSK  EXPEQ')
    elif int(namelist['NEQDSK']) == 0:
         print('Reading from EXPEQ file.')
         expeqParam['nrhomesh']   = [rhomesh_type,rhomesh_src]
         expeqParam['nppfun']     = [pressure_type,pressure_src]
         expeqParam['nsttp']      = [current_type,current_src]
         expeqParam['boundary']   =  boundary_src
         expeqParam['cheasemode'] =  1
         expeqParam['ITEXP']      =  ITEXP
         expeqParam['R0EXP']      =  R0EXP
         expeqParam['B0EXP']      =  B0EXP

         if not os.path.isfile('EXPEQ'):
            if   rhomesh_src in [0,'chease']:
                 if   pressure_src in [0,'chease']   and current_src in [0,'chease']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,imported=importedVals)
                 elif pressure_src in [0,'chease']   and current_src in [1,'eqdsk']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,eqdsk=eqdskfpath,imported=importedVals)
                 elif pressure_src in [0,'chease']   and current_src in [2,'expeq']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [0,'chease']   and current_src in [7,'imported']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,imported=importedVals)
                 elif pressure_src in [1,'eqdsk']    and current_src in [0,'chease']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,eqdsk=eqdskfpath,imported=importedVals)
                 elif pressure_src in [1,'eqdsk']    and current_src in [1,'eqdsk']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,eqdsk=eqdskfpath,imported=importedVals)
                 elif pressure_src in [1,'eqdsk']    and current_src in [2,'expeq']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,eqdsk=eqdskfpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [1,'eqdsk']    and current_src in [7,'imported']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,eqdsk=eqdskfpath,imported=importedVals)
                 elif pressure_src in [2,'expeq']    and current_src in [0,'chease']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [2,'expeq']    and current_src in [1,'eqdsk']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,eqdsk=eqdskfpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [2,'expeq']    and current_src in [2,'expeq']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [2,'expeq']    and current_src in [7,'imported']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [3,'exptnz']   and current_src in [0,'chease']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,exptnz=exptnzfpath,imported=importedVals)
                 elif pressure_src in [3,'exptnz']   and current_src in [1,'eqdsk']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,exptnz=exptnzfpath,eqdsk=eqdskfpath,imported=importedVals)
                 elif pressure_src in [3,'exptnz']   and current_src in [2,'expeq']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,exptnz=exptnzfpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [3,'exptnz']   and current_src in [7,'imported']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,exptnz=exptnzfpath,imported=importedVals)
                 elif pressure_src in [4,'profiles'] and current_src in [0,'chease']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,profiles=profilesfpath,imported=importedVals)
                 elif pressure_src in [4,'profiles'] and current_src in [1,'eqdsk']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,profiles=profilesfpath,eqdsk=eqdskfpath,imported=importedVals)
                 elif pressure_src in [4,'profiles'] and current_src in [2,'expeq']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,profiles=profilesfpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [4,'profiles'] and current_src in [7,'imported']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,profiles=profilesfpath,imported=importedVals)
                 elif pressure_src in [5,'iterdb']   and current_src in [0,'chease']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,iterdb=iterdbfpath,imported=importedVals)
                 elif pressure_src in [5,'iterdb']   and current_src in [1,'eqdsk']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,iterdb=iterdbfpath,eqdsk=eqdskfpath,imported=importedVals)
                 elif pressure_src in [5,'iterdb']   and current_src in [2,'expeq']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,iterdb=iterdbfpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [5,'iterdb']   and current_src in [7,'imported']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,iterdb=iterdbfpath,imported=importedVals)
                 elif pressure_src in [7,'imported'] and current_src in [0,'chease']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,imported=importedVals)
                 elif pressure_src in [7,'imported'] and current_src in [1,'eqdsk']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,eqdsk=eqdskfpath,imported=importedVals)
                 elif pressure_src in [7,'imported'] and current_src in [2,'expeq']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [7,'imported'] and current_src in [7,'imported']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,imported=importedVals)

            elif rhomesh_src in [1,'eqdsk']:
                 if   pressure_src in [0,'chease']   and current_src in [0,'chease']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,chease=cheasefpath,imported=importedVals)
                 elif pressure_src in [0,'chease']   and current_src in [1,'eqdsk']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,chease=cheasefpath,imported=importedVals)
                 elif pressure_src in [0,'chease']   and current_src in [2,'expeq']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,chease=cheasefpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [0,'chease']   and current_src in [7,'imported']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,chease=cheasefpath,imported=importedVals)
                 elif pressure_src in [1,'eqdsk']    and current_src in [0,'chease']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,chease=cheasefpath,imported=importedVals)
                 elif pressure_src in [1,'eqdsk']    and current_src in [1,'eqdsk']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,imported=importedVals)
                 elif pressure_src in [1,'eqdsk']    and current_src in [2,'expeq']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [1,'eqdsk']    and current_src in [7,'imported']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,imported=importedVals)
                 elif pressure_src in [2,'expeq']    and current_src in [0,'chease']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,chease=cheasefpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [2,'expeq']    and current_src in [1,'eqdsk']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [2,'expeq']    and current_src in [2,'expeq']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [2,'expeq']    and current_src in [7,'imported']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [3,'exptnz']   and current_src in [0,'chease']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,chease=cheasefpath,exptnz=exptnzfpath,imported=importedVals)
                 elif pressure_src in [3,'exptnz']   and current_src in [1,'eqdsk']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,exptnz=exptnzfpath,imported=importedVals)
                 elif pressure_src in [3,'exptnz']   and current_src in [2,'expeq']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,exptnz=exptnzfpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [3,'exptnz']   and current_src in [7,'imported']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,exptnz=exptnzfpath,imported=importedVals)
                 elif pressure_src in [4,'profiles'] and current_src in [0,'chease']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,chease=cheasefpath,profiles=profilesfpath,imported=importedVals)
                 elif pressure_src in [4,'profiles'] and current_src in [1,'eqdsk']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,profiles=profilesfpath,imported=importedVals)
                 elif pressure_src in [4,'profiles'] and current_src in [2,'expeq']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,profiles=profilesfpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [4,'profiles'] and current_src in [7,'imported']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,profiles=profilesfpath,imported=importedVals)
                 elif pressure_src in [5,'iterdb']   and current_src in [0,'chease']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,chease=cheasefpath,iterdb=iterdbfpath,imported=importedVals)
                 elif pressure_src in [5,'iterdb']   and current_src in [1,'eqdsk']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,iterdb=iterdbfpath,imported=importedVals)
                 elif pressure_src in [5,'iterdb']   and current_src in [2,'expeq']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,iterdb=iterdbfpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [5,'iterdb']   and current_src in [7,'imported']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,iterdb=iterdbfpath,imported=importedVals)
                 elif pressure_src in [7,'imported'] and current_src in [0,'chease']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,chease=cheasefpath,imported=importedVals)
                 elif pressure_src in [7,'imported'] and current_src in [1,'eqdsk']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,imported=importedVals)
                 elif pressure_src in [7,'imported'] and current_src in [2,'expeq']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [7,'imported'] and current_src in [7,'imported']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,imported=importedVals) 

            elif rhomesh_src in [7,'imported',None]:
                 if   pressure_src in [0,'chease']   and current_src in [0,'chease']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,imported=importedVals)
                 elif pressure_src in [0,'chease']   and current_src in [1,'eqdsk']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,chease=cheasefpath,imported=importedVals)
                 elif pressure_src in [0,'chease']   and current_src in [2,'expeq']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [0,'chease']   and current_src in [7,'imported']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,imported=importedVals)
                 elif pressure_src in [1,'eqdsk']    and current_src in [0,'chease']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,chease=cheasefpath,imported=importedVals)
                 elif pressure_src in [1,'eqdsk']    and current_src in [1,'eqdsk']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,imported=importedVals)
                 elif pressure_src in [1,'eqdsk']    and current_src in [2,'expeq']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [1,'eqdsk']    and current_src in [7,'imported']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,imported=importedVals)
                 elif pressure_src in [2,'expeq']    and current_src in [0,'chease']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [2,'expeq']    and current_src in [1,'eqdsk']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [2,'expeq']    and current_src in [2,'expeq']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [2,'expeq']    and current_src in [7,'imported']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [3,'exptnz']   and current_src in [0,'chease']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,exptnz=exptnzfpath,imported=importedVals)
                 elif pressure_src in [3,'exptnz']   and current_src in [1,'eqdsk']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,exptnz=exptnzfpath,imported=importedVals)
                 elif pressure_src in [3,'exptnz']   and current_src in [2,'expeq']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,exptnz=exptnzfpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [3,'exptnz']   and current_src in [7,'imported']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,exptnz=exptnzfpath,imported=importedVals)
                 elif pressure_src in [4,'profiles'] and current_src in [0,'chease']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,profiles=profilesfpath,imported=importedVals)
                 elif pressure_src in [4,'profiles'] and current_src in [1,'eqdsk']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,profiles=profilesfpath,imported=importedVals)
                 elif pressure_src in [4,'profiles'] and current_src in [2,'expeq']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,profiles=profilesfpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [4,'profiles'] and current_src in [7,'imported']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,profiles=profilesfpath,imported=importedVals)
                 elif pressure_src in [5,'iterdb']   and current_src in [0,'chease']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,iterdb=iterdbfpath,imported=importedVals)
                 elif pressure_src in [5,'iterdb']   and current_src in [1,'eqdsk']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,iterdb=iterdbfpath,imported=importedVals)
                 elif pressure_src in [5,'iterdb']   and current_src in [2,'expeq']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,iterdb=iterdbfpath,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [5,'iterdb']   and current_src in [7,'imported']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,iterdb=iterdbfpath,imported=importedVals)
                 elif pressure_src in [7,'imported'] and current_src in [0,'chease']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,chease=cheasefpath,imported=importedVals)
                 elif pressure_src in [7,'imported'] and current_src in [1,'eqdsk']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,imported=importedVals)
                 elif pressure_src in [7,'imported'] and current_src in [2,'expeq']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,expeq=expeqfpath,imported=importedVals)
                 elif pressure_src in [7,'imported'] and current_src in [7,'imported']:
                      expeqdata = cheasefiles.write_expeq(setParam=expeqParam,imported=importedVals)

    namelistParam['R0EXP'] = R0EXP
    namelistParam['B0EXP'] = B0EXP

    if   current_src in [0,'chease']:
         if   rhomesh_src in [0,'chease',None]:
              cheasedata = cheasefiles.read_chease(fpath=cheasefpath)
         elif rhomesh_src in [1,'eqdsk']:
              cheasedata = cheasefiles.read_chease(fpath=cheasefpath,eqdsk=eqdskfpath)
         if 'QSPEC' in namelistVals and 'CSSPEC' in namelistVals:
            namelistParam['QSPEC']  = namelistVals['QSPEC']
            namelistParam['CSSPEC'] = namelistVals['CSSPEC']
         else:
            if 'q' in importedVals:
               namelistParam['QSPEC'] = importedVals['q'][0]
            else:
               namelistParam['QSPEC'] = cheasedata['q'][0]
            namelistParam['CSSPEC'] = 0.0
    elif current_src in [1,'eqdsk']:
         if   rhomesh_src in [0,'chease']:
              eqdskdata = cheasefiles.read_eqdsk(fpath=eqdskfpath,chease=cheasefpath)
         elif rhomesh_src in [1,'eqdsk',None]:
              eqdskdata = cheasefiles.read_eqdsk(fpath=eqdskfpath)
         if 'QSPEC' in namelistVals and 'CSSPEC' in namelistVals:
            namelistParam['QSPEC']  = namelistVals['QSPEC']
            namelistParam['CSSPEC'] = namelistVals['CSSPEC']
         else:
            if 'q' in importedVals:
               namelistParam['QSPEC'] = importedVals['q'][0]
            else:
               namelistParam['QSPEC'] = eqdskdata['q'][0]
            namelistParam['CSSPEC'] = 0.0
    elif current_src in [2,'expeq']:
         if   rhomesh_src in [0,'chease']:
              expeqdata = cheasefiles.read_expeq(fpath=expeqfpath,chease=cheasefpath)
         elif rhomesh_src in [1,'eqdsk']:
              expeqdata = cheasefiles.read_expeq(fpath=expeqfpath,eqdsk=eqdskfpath)
         elif rhomesh_src in [None]:
              expeqdata = cheasefiles.read_expeq(fpath=expeqfpath)
         if 'QSPEC' in namelistVals and 'CSSPEC' in namelistVals:
            namelistParam['QSPEC']  = namelistVals['QSPEC']
            namelistParam['CSSPEC'] = namelistVals['CSSPEC']
         else:
            if   'q' in importedVals:
                 namelistParam['QSPEC'] = importedVals['q'][0]
            elif 'q' in expeqdata:
                 namelistParam['QSPEC'] = expeqdata['q'][0]
            namelistParam['CSSPEC'] = 0.0
    namelistParam['NCSCAL'] = 4

    namelist = create_namelist(setParam=namelistParam)
       
    return 1

def update_chease_outputs(**kwargs):
    if 'suffix' in kwargs:
       it = int(kwargs['suffix'])
    else:
       it = max(0,len(glob('./chease_*.dat'))-1)
    if os.path.isfile('./chease_namelist'):        os.system('cp ./chease_namelist ./chease_namelist_iter%03d' % it)
    if os.path.isfile('./ogyropsi.h5'):            os.system('mv ./ogyropsi.h5 chease_iter%03d.h5'             % it)
    if os.path.isfile('./EQDSK'):                  os.system('mv ./EQDSK EQDSK_iter%03d.IN'                    % it)
    if os.path.isfile('./EXPEQ'):                  os.system('mv ./EXPEQ EXPEQ_iter%03d.IN'                    % it)
    if os.path.isfile('./EXPTNZ'):                 os.system('mv ./EXPTNZ EXPTNZ_iter%03d.IN'                  % it)
    if os.path.isfile('./EXPEQ.OUT'):              os.system('mv ./EXPEQ.OUT EXPEQ_iter%03d.OUT'               % it)
    if os.path.isfile('./EXPTNZ.OUT'):             os.system('mv ./EXPTNZ.OUT EXPTNZ_iter%03d.OUT'             % it)
    if os.path.isfile('./EXPEQ_EXPEQ.IN'):         os.system('mv ./EXPEQ_EXPEQ.IN EXPEQ_EXPEQ_iter%03d.IN'     % it)
    if os.path.isfile('./EQDSK_EXPEQ.IN'):         os.system('mv ./EQDSK_EXPEQ.IN EQDSK_EXPEQ_iter%03d.IN'     % it)
    if os.path.isfile('./ogyropsi.dat'):           os.system('mv ./ogyropsi.dat chease_iter%03d.dat'           % it)
    if os.path.isfile('./EQDSK_COCOS_02_POS.OUT'): os.system('mv ./EQDSK_COCOS_02_POS.OUT EQDSK_iter%03d'      % it)
   #   os.system('mv ./EQDSK_COCOS_02_POS.OUT EQDSK_COCOS_02_iter%03d.OUT' % it)
   #if os.path.isfile('./EQDSK_COCOS_12.OUT'):
   #   os.system('mv ./EQDSK_COCOS_12.OUT EQDSK_COCOS_12_iter%03d.OUT'     % it)

    return 1

def runchease(pltVals={}):
    process      = Popen("./bin_chease chease_namelist > chease.log", shell=True, stdout=PIPE)
    (output,err) = process.communicate()
    exit_status  = process.wait()
    if abs(exit_status) > 0:
       print("Exit Status: ",exit_status)
       print("Error Message:\n",err)
       sys.exit()

    return 1

def main():
    if len(sys.argv) > 1:
       if sys.argv[1] == 'remove-inputs':  remove_input_files()
       if sys.argv[1] == 'remove-outputs': remove_output_files()
       if sys.argv[1] == 'remove-all':     remove_input_files(); remove_output_files()

if __name__ == "__main__":
    main()

