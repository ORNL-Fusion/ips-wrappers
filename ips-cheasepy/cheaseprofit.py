#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
import warnings
import cheasefiles

import numpy as npy
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit


if   sys.version_info.major == 3:
     PYTHON3 = True; PYTHON2 = False
elif sys.version_info.major == 2:
     PYTHON2 = True; PYTHON3 = False

def fit_profile(rhotor,in_profile,method='stefanikova',setParam={},fitParam={},fitBounds={}):
    global profile

    if 'alpha'      in fitParam and fitParam['alpha'] != None:
       global set_alpha;      set_alpha      = fitParam['alpha']
    if 'ped_sol'    in fitParam and fitParam['ped_sol'] != None:
       global set_ped_sol;    set_ped_sol    = fitParam['ped_sol']
    if 'ped_mid'    in fitParam and fitParam['ped_mid'] != None:
       global set_ped_mid;    set_ped_mid    = fitParam['ped_mid']
    if 'ped_width'  in fitParam and fitParam['ped_width'] != None:
       global set_ped_width;  set_ped_width  = fitParam['ped_width']
    if 'ped_height' in fitParam and fitParam['ped_height'] != None:
       global set_ped_height; set_ped_height = fitParam['ped_height']
    if 'cor_exp'    in fitParam and fitParam['cor_exp'] != None:
       global set_cor_exp;    set_cor_exp    = fitParam['cor_exp']
    if 'cor_width'  in fitParam and fitParam['cor_width'] != None:
       global set_cor_width;  set_cor_width  = fitParam['cor_width']
    if 'cor_height' in fitParam and fitParam['cor_height'] != None:
       global set_cor_height; set_cor_height = fitParam['cor_height']


    if   method.lower()=='stefanikova':
         lw_bound = [-npy.inf for i in range(8)]
         up_bound = [ npy.inf for i in range(8)]
    elif method.lower()=='groebner':
         lw_bound = [-npy.inf for i in range(5)]
         up_bound = [ npy.inf for i in range(5)]

    if 'alpha'      in fitBounds and fitBounds['alpha'] != None:
       if type(fitBounds['alpha']) in [list,tuple]:
          lw_bound[0] = fitBounds['alpha'][0]
          up_bound[0] = fitBounds['alpha'][1]
       else:
          lw_bound[0] = fitBounds['alpha']
    if 'ped_sol'    in fitBounds and fitBounds['ped_sol'] != None:
       if type(fitBounds['ped_sol']) in [list,tuple]:
          lw_bound[1] = fitBounds['ped_sol'][0]
          up_bound[1] = fitBounds['ped_sol'][1]
       else:
          lw_bound[1] = fitBounds['ped_sol']
    if 'ped_mid'    in fitBounds and fitBounds['ped_mid'] != None:
       if type(fitBounds['ped_mid']) in [list,tuple]:
          lw_bound[2] = fitBounds['ped_mid'][0]
          up_bound[2] = fitBounds['ped_mid'][1]
       else:
          lw_bound[2] = fitBounds['ped_mid']
    if 'ped_width'  in fitBounds and fitBounds['ped_width'] != None:
       if type(fitBounds['ped_width']) in [list,tuple]:
          lw_bound[3] = fitBounds['ped_width'][0]
          up_bound[3] = fitBounds['ped_width'][1]
       else:
          lw_bound[3] = fitBounds['ped_width']
    if 'ped_height' in fitBounds and fitBounds['ped_height'] != None:
       if type(fitBounds['ped_height']) in [list,tuple]:
          lw_bound[4] = fitBounds['ped_height'][0]
          up_bound[4] = fitBounds['ped_height'][1]
       else:
          lw_bound[4] = fitBounds['ped_height']
    if 'cor_exp'    in fitBounds and fitBounds['cor_exp'] != None:
       if type(fitBounds['cor_exp']) in [list,tuple]:
          lw_bound[5] = fitBounds['cor_exp'][0]
          up_bound[5] = fitBounds['cor_exp'][1]
       else:
          lw_bound[5] = fitBounds['cor_exp']
    if 'cor_width'  in fitBounds and fitBounds['cor_width'] != None:
       if type(fitBounds['cor_width']) in [list,tuple]:
          lw_bound[6] = fitBounds['cor_width'][0]
          up_bound[6] = fitBounds['cor_width'][1]
       else:
          lw_bound[6] = fitBounds['cor_width']
    if 'cor_height' in fitBounds and fitBounds['cor_height'] != None:
       if type(fitBounds['cor_height']) in [list,tuple]:
          lw_bound[7] = fitBounds['cor_height'][0]
          up_bound[7] = fitBounds['cor_height'][1]
       else:
          lw_bound[7] = fitBounds['cor_height']
    fit_bounds=(lw_bound,up_bound)


    if 'norm' in setParam:
       norm = setParam['norm']
    else:
       norm = False

    if 'plot' in setParam:
       fit_plot = setParam['plot']
    else:
       fit_plot = True
    print(fit_plot)

    if norm:
       fmin    = npy.min(in_profile)
       fmax    = npy.max(in_profile)
    else:
       fmin    = 0.0
       fmax    = 1.0

    profile = (in_profile-fmin)/(fmax-fmin)
    if 'set_ped_sol'    in globals(): set_ped_sol    = (set_ped_sol-fmin)/(fmax-fmin)
    if 'set_ped_height' in globals(): set_ped_height = (set_ped_height-fmin)/(fmax-fmin)
    if 'set_cor_height' in globals(): set_cor_height = (set_cor_height-fmin)/(fmax-fmin)

    if   method.lower()=='stefanikova':
         try:
            popt,pcov = curve_fit(stefanikova_fit,rhotor,profile,bounds=fit_bounds)
         except RuntimeError:
            result = [npy.nan for i in range(8)]
            popt,pcov = (result,npy.nan)

         fit_parameters               = {}
         if 'set_alpha'      in globals(): fit_parameters['alpha']      = set_alpah
         else:                             fit_parameters['alpha']      = popt[0]
         if 'set_ped_height' in globals(): fit_parameters['ped_height'] = (fmax-fmin)*set_ped_height+fmin
         else:                             fit_parameters['ped_height'] = (fmax-fmin)*popt[1]+fmin
         if 'set_ped_sol'    in globals(): fit_parameters['ped_sol']    = (fmax-fmin)*set_ped_sol+fmin
         else:                             fit_parameters['ped_sol']    = (fmax-fmin)*popt[2]+fmin
         if 'set_ped_width'  in globals(): fit_parameters['ped_width']  = set_ped_width
         else:                             fit_parameters['ped_width']  = popt[3]
         if 'set_ped_mid'    in globals(): fit_parameters['ped_mid']    = set_ped_mid
         else:                             fit_parameters['ped_mid']    = popt[4]
         if 'set_cor_exp'    in globals(): fit_parameters['cor_exp']    = set_cor_exp
         else:                             fit_parameters['cor_exp']    = popt[5]
         if 'set_cor_width'  in globals(): fit_parameters['cor_width']  = set_cor_width
         else:                             fit_parameters['cor_width']  = popt[6]
         if 'set_cor_height' in globals(): fit_parameters['cor_height'] = (fmax-fmin)*set_cor_height+fmin
         else:                            fit_parameters['cor_height'] = (fmax-fmin)*popt[7]+fmin

         if fit_plot:
            print('Plotting Stefanikova Profile.')
            plt.plot(rhotor,stefanikova_fit(rhotor,*popt),'r',label='CRVFIT')
            plt.plot(rhotor,profile,'b--',label='EXP')
            plt.legend()
            plt.savefig(fname="sfitting.png",dpi=100)
            plt.close()

    elif method.lower()=='groebner':
         try:
            popt,pcov = curve_fit(groebner_fit,rhotor,profile,bounds=fit_bounds)
         except RuntimeError:
            result = [npy.nan for i in range(5)]
            popt,pcov = (result,npy.nan)

         fit_parameters               = {}
         if 'set_alpha'      in globals(): fit_parameters['alpha']      = set_alpah
         else:                             fit_parameters['alpha']      = popt[0]
         if 'set_ped_height' in globals(): fit_parameters['ped_height'] = (fmax-fmin)*set_ped_height+fmin
         else:                             fit_parameters['ped_height'] = (fmax-fmin)*popt[1]+fmin
         if 'set_ped_sol'    in globals(): fit_parameters['ped_sol']    = (fmax-fmin)*set_ped_sol+fmin
         else:                             fit_parameters['ped_sol']    = (fmax-fmin)*popt[2]+fmin
         if 'set_ped_width'  in globals(): fit_parameters['ped_width']  = set_ped_width
         else:                             fit_parameters['ped_width']  = popt[3]
         if 'set_ped_mid'    in globals(): fit_parameters['ped_mid']    = set_ped_mid
         else:                             fit_parameters['ped_mid']    = popt[4]

         if fit_plot:
            print('Plotting Groebner Profile.')
            plt.plot(rhotor,groebner_fit(rhotor,*popt),'r',label='CRVFIT')
            plt.plot(rhotor,profile,'b--',label='EXP')
            plt.legend()
            plt.savefig(fname="gfitting.png",dpi=100)
            plt.close()

    if 'set_alpha'      in globals(): del set_alpha
    if 'set_ped_sol'    in globals(): del set_ped_sol
    if 'set_ped_mid'    in globals(): del set_ped_mid
    if 'set_ped_width'  in globals(): del set_ped_width
    if 'set_ped_height' in globals(): del set_ped_height
    if 'set_cor_exp'    in globals(): del set_cor_exp
    if 'set_cor_width'  in globals(): del set_cor_width
    if 'set_cor_height' in globals(): del set_cor_height

    return fit_parameters


def stefanikova_fit(rho,alpha,ped_height,ped_sol,ped_width,ped_mid,cor_exp,cor_width,cor_height):
    profileFlag=False
    if  'profile' in locals() or 'profile' in globals():
        f = profile[:]
        profileFlag=True
    else:
        print("No profile in the locals() or globals()?")
        sys.exit()

    if 'set_alpha'      in globals(): alpha      = set_alpha
    if 'set_ped_sol'    in globals(): ped_sol    = set_ped_sol
    if 'set_ped_mid'    in globals(): ped_mid    = set_ped_mid
    if 'set_ped_width'  in globals(): ped_width  = set_ped_width
    if 'set_ped_height' in globals(): ped_height = set_ped_height
    if 'set_cor_exp'    in globals(): cor_exp    = set_cor_exp
    if 'set_cor_width'  in globals(): cor_width  = set_cor_width
    if 'set_cor_height' in globals(): cor_height = set_cor_height

    arg       = 2.0*(ped_mid-rho)/ped_width
    mtanh     = ((1.0+alpha*arg)*npy.exp(arg)-npy.exp(-arg))/(npy.exp(arg)+npy.exp(-arg))
    fprof_ped = ((ped_height-ped_sol)*(mtanh+1.0)/2.0)+ped_sol
    fprof_cor = (cor_height-fprof_ped)*npy.exp(-npy.power(rho/cor_width,cor_exp))
    fprof_ful = fprof_cor+fprof_ped

    return fprof_ful

def groebner_fit(rho,alpha,ped_height,ped_sol,ped_width,ped_mid):
    profileFlag=False; refprofFlag=False
    if  'profile' in locals() or 'profile' in globals():
        f = profile[:]
        profileFlag=True
    else:
        print("No profile in the locals() or globals()?")
        sys.exit()

    if 'set_alpha'      in globals(): alpha      = set_alpha
    if 'set_ped_sol'    in globals(): ped_sol    = set_ped_sol
    if 'set_ped_mid'    in globals(): ped_mid    = set_ped_mid
    if 'set_ped_width'  in globals(): ped_width  = set_ped_width
    if 'set_ped_height' in globals(): ped_height = set_ped_height

    A = (ped_height-ped_sol)/2.0
    B = (ped_height+ped_sol)/2.0

    arg       = 2.0*(ped_mid-rho)/ped_width
    mtanh     = ((1.0+alpha*arg)*npy.exp(arg)-npy.exp(-arg))/(npy.exp(arg)+npy.exp(-arg))
    fprof_ped = A*mtanh+B

    return fprof_ped


def find_pedestal(profilefpath,profileftype,setParam={},**kwargs):
    pedParam = {}

    eqdskflag    = False
    cheaseflag   = False
    for key,value in kwargs.items():
        if   key in ['chease','cheasedata','cheasefpath']:
             cheaseflag = True
             cheasepath = value.strip()
        elif key in ['eqdsk','eqdskdata','eqdskfpath']:
             eqdskflag = True
             eqdskpath = value.strip()

    if 'nrhomesh' not in setParam:
       setParam = {'nrhomesh':1}

    if os.path.isfile(profilefpath):
       if   eqdskflag:
            if   profileftype.lower()=='iterdb':
                 profiledata = cheasefiles.read_iterdb(fpath=profilefpath,setParam=setParam,eqdsk=eqdskpath)
            elif profileftype.lower()=='chease':
                 profiledata = cheasefiles.read_chease(fpath=profilefpath,setParam=setParam,eqdsk=eqdskpath)
            elif profileftype.lower()=='profiles':
                 profiledata = cheasefiles.read_profiles(fpath=profilefpath,setParam=setParam,eqdsk=eqdskpath)
       elif cheaseflag:
            if   profileftype.lower()=='iterdb':
                 profiledata = cheasefiles.read_iterdb(fpath=profilefpath,setParam=setParam,chease=cheasepath)
            elif profileftype.lower()=='chease':
                 profiledata = cheasefiles.read_chease(fpath=profilefpath)
            elif profileftype.lower()=='profiles':
                 profiledata = cheasefiles.read_profiles(fpath=profilefpath,setParam=setParam,chease=cheasepath)
       else:
            if   profileftype.lower()=='iterdb':
                 profiledata = cheasefiles.read_iterdb(fpath=profilefpath)
            elif profileftype.lower()=='chease':
                 profiledata = cheasefiles.read_chease(fpath=profilefpath)
            elif profileftype.lower()=='profiles':
                 profiledata = cheasefiles.read_profiles(fpath=profilefpath)

    if setParam['nrhomesh'] == 1:
       rho = profiledata['rhotor']
    else:
       rho = profiledata['rhopsi']

    for prof in ('ne','Te','pressure'):
        f        = profiledata[prof]
        df       = npy.gradient(f,edge_order=2)
        d2f      = npy.gradient(df,edge_order=2)

        lowess = sm.nonparametric.lowess(d2f,rho,frac=0.005)
        d2f = lowess[:,1]

        ped_sol_ind = npy.argmin(abs(f-npy.min(f)))
        ped_mid_ind = npy.argmin(abs(df-npy.min(df)))
        ped_top_ind = npy.argmin(abs(d2f-npy.min(d2f)))
        ped_fot_ind = npy.argmin(abs(d2f-npy.max(d2f)))

        ped_top     = rho[ped_top_ind]
        ped_mid     = rho[ped_mid_ind]
        ped_fot     = rho[ped_fot_ind]

        if (ped_mid-ped_top)>(ped_fot-ped_mid):
           ped_hwd     = ped_mid-ped_top
           ped_fot     = ped_mid+ped_hwd
           ped_fot_ind = npy.argmin(abs(rho-ped_fot))
        else:
           ped_hwd     = ped_fot-ped_mid
           ped_top     = ped_mid-ped_hwd
           ped_top_ind = npy.argmin(abs(rho-ped_top))

        pedParam[prof]               = {}
        pedParam[prof]['ped_sol']    = f[ped_sol_ind]
        pedParam[prof]['ped_mid']    = ped_mid
        pedParam[prof]['ped_width']  = ped_fot-ped_top
        pedParam[prof]['ped_slope']  = -1.0/f[ped_mid_ind]*df[ped_mid_ind]
        pedParam[prof]['ped_height'] = f[ped_top_ind]

    return pedParam,profiledata


def get_iter_geometry():
    a =  2.0
    R =  6.2
    B =  5.3
    kappa = 1.8
    delta = 0.4

    epsilon = a/R
    alpha   = math.asin(delta)
    tau     = npy.linspace(0.0,2.0*math.pi,100)

    iter_geometry = {}
    iter_geometry['R0EXP']  = R
    iter_geometry['B0EXP']  = B
    iter_geometry['rbound'] = R*(1.0+epsilon*npy.cos(tau+alpha*npy.sin(tau)))
    iter_geometry['zbound'] = R*(epsilon*kappa*npy.sin(tau))

    return iter_geometry


def get_iter_profiles(iterdbfpath,eqdskfpath,setparam={},plotParam={}):
    global profile
    warnings.filterwarnings("ignore")

    if 'plot' in plotParam:
       plotprof = plotParam['plot']
    else:
       plotprof = False

    iterdbdata = cheasefiles.read_iterdb(fpath=iterdbfpath,setParam={'nrhomesh':1},eqdsk=eqdskfpath)

    if plotprof:
       fg1 = plt.figure('Electron Density')
       ax1 = fg1.add_subplot(111)

       fg2 = plt.figure('Electron Temperature')
       ax2 = fg2.add_subplot(111)

       fg3 = plt.figure('Ion Density')
       ax3 = fg3.add_subplot(111)

       fg4 = plt.figure('Ion Temperature')
       ax4 = fg4.add_subplot(111)

       fg5 = plt.figure('Pressure')
       ax5 = fg5.add_subplot(111)

    iter_profiles = {}

    setparam = {'norm':False,'profiledata':iterdbdata,'fit_plot':False}
    pres_fit_param = fit_profile(iterdbdata['rhotor'],iterdbdata['pressure'],method='groebner',setParam=setparam)

    fitparam = {}
    setparam = {'norm':True,'profiledata':iterdbdata,'fit_plot':False}
    dens_fit_param = fit_profile(iterdbdata['rhotor'],iterdbdata['ne'],method='groebner',setParam=setparam,fitParam=fitparam,fitBounds={})
    fitparam=[]
    fitparam.append(dens_fit_param['alpha'])
    fitparam.append(dens_fit_param['ped_height'])
    fitparam.append(dens_fit_param['ped_sol'])
    fitparam.append(dens_fit_param['ped_width'])
    fitparam.append(dens_fit_param['ped_mid'])
    neProfile = groebner_fit(iterdbdata['rhotor'],*fitparam)
    if plotprof:
       ax1.plot(iterdbdata['rhotor'],neProfile,label='Original DIII-D')
       ax1.plot(iterdbdata['rhotor'],iterdbdata['ne'],label='Fitted DIII-D')
    fitparam[3]=dens_fit_param['ped_width']+0.05
    fitparam[4]=0.9
    neProfile = groebner_fit(iterdbdata['rhotor'],*fitparam)
    if plotprof:
       ax1.plot(iterdbdata['rhotor'],neProfile,label='Wider Pedestal')
    iter_ne_profile = scale_to_iter(neProfile,1.0e20,0.33e20)
    if plotprof:
       ax1.plot(iterdbdata['rhotor'],iter_ne_profile,label='Scale-to_ITER')
       ax1.set_title('Electron Density')
       ax1.set_ylabel('$n_e$ $(m^{-3})$')
       ax1.set_xlabel('$\\rho_{\\phi}$')

    fitparam = {}
    setparam = {'norm':True,'profiledata':iterdbdata,'fit_plot':False}
    temp_fit_param = fit_profile(iterdbdata['rhotor'],iterdbdata['Te'],method='stefanikova',setParam=setparam,fitParam=fitparam,fitBounds={})
    fitparam=[]
    fitparam.append(temp_fit_param['alpha'])
    fitparam.append(temp_fit_param['ped_height'])
    fitparam.append(temp_fit_param['ped_sol'])
    fitparam.append(temp_fit_param['ped_width'])
    fitparam.append(temp_fit_param['ped_mid'])
    fitparam.append(temp_fit_param['cor_exp'])
    fitparam.append(temp_fit_param['cor_width'])
    fitparam.append(temp_fit_param['cor_height'])
    TeProfile = stefanikova_fit(iterdbdata['rhotor'],*fitparam)
    if plotprof:
       ax2.plot(iterdbdata['rhotor'],TeProfile,label='Original DIII-D')
       ax2.plot(iterdbdata['rhotor'],iterdbdata['Te'],label='Fitted DIII-D')
    fitparam[3]=temp_fit_param['ped_width']+0.05
    fitparam[4]=0.9
    TeProfile = stefanikova_fit(iterdbdata['rhotor'],*fitparam)
    if plotprof:
       ax2.plot(iterdbdata['rhotor'],TeProfile,label='Wider Pedestal')
    iter_Te_profile = scale_to_iter(TeProfile,8000.0,100.0)
    if plotprof:
       ax2.plot(iterdbdata['rhotor'],iter_Te_profile,label='Scale-to_ITER')
       ax2.set_title('Electron Temperature')
       ax2.set_ylabel('$T_e$ $(eV)$')
       ax2.set_xlabel('$\\rho_{\\phi}$')
       ax2.legend()

    fitparam = {'ped_mid':pres_fit_param['ped_mid']}
    setparam = {'norm':True,'profiledata':iterdbdata,'fit_plot':False}
    dens_fit_param = fit_profile(iterdbdata['rhotor'],iterdbdata['ni'],method='groebner',setParam=setparam,fitParam=fitparam,fitBounds={})
    fitparam=[]
    fitparam.append(dens_fit_param['alpha'])
    fitparam.append(dens_fit_param['ped_height'])
    fitparam.append(dens_fit_param['ped_sol'])
    fitparam.append(dens_fit_param['ped_width'])
    fitparam.append(dens_fit_param['ped_mid'])
    niProfile = groebner_fit(iterdbdata['rhotor'],*fitparam)
    if plotprof:
       ax3.plot(iterdbdata['rhotor'],niProfile,label='Original DIII-D')
       ax3.plot(iterdbdata['rhotor'],iterdbdata['ni'],label='Fitted DIII-D')
    fitparam[3]=dens_fit_param['ped_width']+0.05
    fitparam[4]=0.9
    niProfile = groebner_fit(iterdbdata['rhotor'],*fitparam)
    if plotprof:
       ax3.plot(iterdbdata['rhotor'],niProfile,label='Wider Pedestal')
    iter_ni_profile = scale_to_iter(niProfile,1.0e20,0.33e20)
    if plotprof:
       ax3.plot(iterdbdata['rhotor'],iter_ni_profile,label='Scale-to-ITER')
       ax3.set_title('Ion Density')
       ax3.set_ylabel('$n_i$ $(m^{-3})$')
       ax3.set_xlabel('$\\rho_{\\phi}$')
       ax3.legend()

    fitparam = {'cor_height':iterdbdata['Ti'][-1],'ped_mid':pres_fit_param['ped_mid']}
    setparam = {'norm':True,'profiledata':iterdbdata,'fit_plot':False}
    temp_fit_param = fit_profile(iterdbdata['rhotor'],iterdbdata['Ti'],method='groebner',setParam=setparam,fitParam=fitparam,fitBounds={})
    fitparam=[]
    fitparam.append(temp_fit_param['alpha'])
    fitparam.append(temp_fit_param['ped_height'])
    fitparam.append(temp_fit_param['ped_sol'])
    fitparam.append(temp_fit_param['ped_width'])
    fitparam.append(temp_fit_param['ped_mid'])
    TiProfile = groebner_fit(iterdbdata['rhotor'],*fitparam)
    if plotprof:
       ax4.plot(iterdbdata['rhotor'],TiProfile,label='Original DIII-D')
       ax4.plot(iterdbdata['rhotor'],iterdbdata['Ti'],label='Fitted DIII-D')
    fitparam[3]=temp_fit_param['ped_width']+0.05
    fitparam[4]=0.9
    TiProfile = groebner_fit(iterdbdata['rhotor'],*fitparam)
    if plotprof:
       ax4.plot(iterdbdata['rhotor'],TiProfile,label='Wider Pedestal')
    iter_Ti_profile = scale_to_iter(TiProfile,8000.0,100.0)
    if plotprof:
       ax4.plot(iterdbdata['rhotor'],iter_Ti_profile,label='Scale-to-ITER')
       ax4.set_title('Ion Temperature')
       ax4.set_ylabel('$T_i$ $(eV)$')
       ax4.set_xlabel('$\\rho_{\\phi}$')
       ax4.legend()

    fitparam=[]
    fitparam.append(pres_fit_param['alpha'])
    fitparam.append(pres_fit_param['ped_height'])
    fitparam.append(pres_fit_param['ped_sol'])
    fitparam.append(pres_fit_param['ped_width'])
    fitparam.append(pres_fit_param['ped_mid'])
    ptProfile = groebner_fit(iterdbdata['rhotor'],*fitparam)
    iter_pt_profile  = iter_ne_profile*iter_Te_profile
    iter_pt_profile += iter_ni_profile*iter_Ti_profile
    iter_pt_profile *= 1.602e-19
    if plotprof:
       ax5.plot(iterdbdata['rhotor'],ptProfile,label='Original DIII-D')
       ax5.plot(iterdbdata['rhotor'],iter_pt_profile,label="Scale-to-ITER")
       ax5.set_title('Pressure')
       ax5.set_ylabel('$P_T$ $(Pa)$')
       ax5.set_xlabel('$\\rho_{\\phi}$')
       ax5.legend()

    iter_profiles['rhotor']   = iterdbdata['rhotor']
    iter_profiles['rhopsi']   = iterdbdata['rhopsi']
    iter_profiles['ne']       = iter_ne_profile
    iter_profiles['ni']       = iter_ni_profile
    iter_profiles['Te']       = iter_Te_profile
    iter_profiles['Ti']       = iter_Ti_profile
    iter_profiles['pressure'] = iter_pt_profile

    if plotprof:
       plt.show()

    return iter_profiles


def scale_to_iter(in_profile,iter_profile_avg,iter_profile_sol):
    in_profile_avg = npy.mean(in_profile)
    in_profile_sol = in_profile[-1]

    profile_norm = (in_profile-in_profile_avg)/(in_profile_avg-in_profile_sol)

    iter_profile = profile_norm*(iter_profile_avg-iter_profile_sol)+iter_profile_avg

    return iter_profile

