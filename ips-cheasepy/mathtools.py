#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import numpy as npy

from scipy.integrate   import quad
from scipy.integrate   import trapz,simps
from scipy.interpolate import splrep,splev
from scipy.interpolate import interp1d,interp2d
from scipy.interpolate import CubicSpline,RectBivariateSpline

if   sys.version_info.major == 3:
     PYTHON3 = True; PYTHON2 = False
elif sys.version_info.major == 2:
     PYTHON2 = True; PYTHON3 = False

def bisection(fx,xmin,xmax,root=0.0,Nmax=100,eps=1.0e-16):
    for it in range(Nmax):
        x = (xmax+xmin)/2.0
        if abs(fx(x)-root) < eps: break
        if   (fx(xmax)-root)*(fx(x)-root) < 0.0: xmin = x
        else:                                    xmax = x
    return x


def rk4(f1,f2,czt,slns,stp):
    k01 = f1(slns[0],slns[1])
    k11 = f2(slns[0],slns[1])
    k02 = f1(slns[0]+stp*k01/2.0,slns[1]+stp*k11/2.0)
    k12 = f2(slns[0]+stp*k01/2.0,slns[1]+stp*k11/2.0)
    k03 = f1(slns[0]+stp*k02/2.0,slns[1]+stp*k12/2.0)
    k13 = f2(slns[0]+stp*k02/2.0,slns[1]+stp*k12/2.0)
    k04 = f1(slns[0]+stp*k03,slns[1]+stp*k13)
    k14 = f2(slns[0]+stp*k03,slns[1]+stp*k13)
    slns[0] = slns[0]+stp*(k01+2*k02+2*k03+k04)/6.0
    slns[1] = slns[1]+stp*(k11+2*k12+2*k13+k14)/6.0
    return slns


def rk5(f1,f2,czt,slns,stp):
    k01 = f1(slns[0],slns[1])
    k11 = f2(slns[0],slns[1])

    k02 = f1(slns[0]+stp*k01/4.0,slns[1]+stp*k11/4.0)
    k12 = f2(slns[0]+stp*k01/4.0,slns[1]+stp*k11/4.0)

    k03 = f1(slns[0]+stp*(3.0*k01+9.0*k02)/32.0,slns[1]+stp*(3.0*k11+9.0*k12)/32.0)
    k13 = f2(slns[0]+stp*(3.0*k01+9.0*k02)/32.0,slns[1]+stp*(3.0*k11+9.0*k12)/32.0)

    k04 = f1(slns[0]+stp*(1932*k01-7200.0*k02+7296.0*k03)/2197.0,slns[1]+stp*(1932*k11-7200.0*k12+7296.0*k13)/2197.0)
    k14 = f2(slns[0]+stp*(1932*k01-7200.0*k02+7296.0*k03)/2197.0,slns[1]+stp*(1932*k11-7200.0*k12+7296.0*k13)/2197.0)

    k05 = f1(slns[0]+stp*(439*k01/216.0-8.0*k02+3680.0*k03/513.0-845.0*k04/4104.0),slns[1]+stp*(439*k11/216.0-8.0*k12+3680.0*k13/513.0-845.0*k14/4104.0))
    k15 = f2(slns[0]+stp*(439*k01/216.0-8.0*k02+3680.0*k03/513.0-845.0*k04/4104.0),slns[1]+stp*(439*k11/216.0-8.0*k12+3680.0*k13/513.0-845.0*k14/4104.0))

    k06 = f1(slns[0]+stp*(-8.0*k01/27.0+2.0*k02-3544.0*k03/2565.0+1859.0*k04/4104.0-11.0*k05/40.0),slns[1]+stp*(-8.0*k11/27.0+2.0*k12-3544.0*k13/2565.0+1859.0*k14/4104.0-11.0*k15/40.0))
    k16 = f2(slns[0]+stp*(-8.0*k01/27.0+2.0*k02-3544.0*k03/2565.0+1859.0*k04/4104.0-11.0*k05/40.0),slns[1]+stp*(-8.0*k11/27.0+2.0*k12-3544.0*k13/2565.0+1859.0*k14/4104.0-11.0*k15/40.0))

    slns[0] = slns[0]+stp*(16.0*k01/135.0+6656.0*k03/12825.0+28561.0*k04/56430.0-9.0*k05/50.0+2*k06/55.0)
    slns[1] = slns[1]+stp*(16.0*k11/135.0+6656.0*k13/12825.0+28561.0*k14/56430.0-9.0*k15/50.0+2*k16/55.0)
    return slns


def findmonotonic(A,kind="increasing"):
    if kind.lower()=="increasing":
       bgnloc=0
       endloc=len(A)-1
       for i in range(len(A)):
           if A[i+1]>A[i]:
              bgnloc=i
              break
       for i in range(bgnloc,len(A)-1):
           if A[i+1]<A[i]:
              endloc=i
              break
    return bgnloc,endloc


def interp(xin,fxin,xnew,ynew=[],yout=[]):
    #splrep returns knots and coefficients for cubic spline
    #Use these knots and coefficients with splev to get a new y

    fknots = splrep(xin,fxin)
    fxnew  = splev(xnew,fknots,der=0)

    if len(ynew)>0 and len(yout)>0:
       yknots = splrep(ynew,fxnew)
       fxout  = splev(yout,yknots,der=0)

       x_knots = splrep(ynew,xnew)
       xout    = splev(yout,x_knots,der=0)
    else:
       fxout = fxnew[:]

    return fxout

def integrate(x,fx,axis=0,method='trapz'):
    fxShape = npy.shape(fx)
    nDim    = len(fxShape)
    if   method in ['trapz','trapzoid']:
         if   nDim == 1:
              intf = trapz(y=fx,x=x,dx=npy.argmin(npy.diff(x)))
         elif nDim == 2:
              if   axis == 0:
                   intf = trapz(y=fx,x=x,dx=npy.argmin(npy.diff(x)),axis=0)
              elif axis == 1:
                   intf = trapz(y=fx,x=x,dx=npy.argmin(npy.diff(x)),axis=1)
    elif method in ['simps','simpson']:
         if   nDim == 1:
              intf = simps(y=fx,x=x,dx=npy.argmin(npy.diff(x)))
         elif nDim == 2:
              if   axis == 0:
                   intf = simps(y=fx,x=x,dx=npy.argmin(npy.diff(x)),axis=0)
              elif axis == 1:
                   intf = simps(y=fx,x=x,dx=npy.argmin(npy.diff(x)),axis=1)
    elif method=='CubicSpline':
         m = fxShape[0]
         n = fxShape[1]
         if   nDim == 1:
                   intf = CubicSpline(x,fx).integrate(x[0],x[-1])
         elif nDim == 2:
              if   axis == 0:
                   try:
                      CS   = CubicSpline(x,fx,axis=0,bc_type='periodic',extrapolate='periodic')
                      intf = CS.integrate(x[0],x[-1],extrapolate='periodic')
                   except ValueError:
                      intf = npy.zeros(m)
                      for j in range(m):
                          func  = interp1d(x,fx[:,j],kind='linear')
                          fyfnc = lambda z: func(z)
                          intf[j] = quad(fyfnc,x[0],x[-1])
              elif axis == 1:
                   try:
                      CS   = CubicSpline(x,fx,axis=1)
                      intf = CS.integrate(x[0],x[-1])
                   except ValueError:
                      intf = npy.zeros(n)
                      for j in range(n):
                          func  = interp1d(x,fx[:,j],kind='linear')
                          fyfnc = lambda z: func(z)
                          intf[j] = quad(fyfnc,x[0],x[-1])
    return intf

def derivative(x,fx,axis=0,dorder=1,method='gradient'):
    fxShape = npy.shape(fx)
    nDim    = len(fxShape)
    if   method=='gradient':
         if   nDim == 1:
              dfdx = npy.gradient(fx,x)
    elif method=='CubicSpline':
         if   nDim == 1:
              CS = CubicSpline(x,fx)
              dfdx = CS(x,dorder)
         elif nDim == 2:
              m = fxShape[0]
              n = fxShape[1]
              dfdx=npy.zeros((m,n))
              if   axis == 0:
                   for j in range(n):
                      CS = CubicSpline(x,fx[:,j])
                      dfdx[:,j] = CS(x,dorder)
              elif axis == 1:
                   for i in range(m):
                      CS = CubicSpline(x,fx[i,:])
                      dfdx[i,:] = CS(x,dorder)
    return dfdx

def get_mat_fd_d1_o4(size,dx,plot_matrix=False):
    """Creates matrix for centered finite difference, first derivative, 4th order.
    size: size of (number of elements in) quantity to be differentiated
    dx: grid spacing (for constant grid)."""

    prefactor=1.0/(12.0*dx)
    mat=npy.zeros((size,size),dtype='float')
    for i in range(size):
        if i-1 >= 0:
            mat[i,i-1]=-8
        if i-2 >= 0:
            mat[i,i-2]=1
        if i+1 <= size-1:
            mat[i,i+1]=8
        if i+2 <= size-1:
            mat[i,i+2]=-1

    mat=prefactor*mat

    if plot_matrix:
        plt.contourf(mat,50)
        plt.colorbar()
        plt.show()

    return mat


def fd_d1_o4(var,grid,mat=False):
    """Centered finite difference, first derivative, 4th order.
    var: quantity to be differentiated.
    grid: grid for var
    mat: matrix for the finite-differencing operator. if mat=False then it is created"""

    if not mat:
        mat=get_mat_fd_d1_o4(len(var),grid[1]-grid[0])

    dvar=-npy.dot(mat,var)
    dvar[0]=0.0
    dvar[1]=0.0
    #dvar[2]=0.0
    dvar[-1]=0.0
    dvar[-2]=0.0
    #dvar[-3]=0.0
    return -dvar

