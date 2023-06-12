#!/opt/bin/python
#upgrading to scipy version .. , changing interface
import numpy as np
import numpy.fft as ft
import scipy.io.netcdf as nc
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib import cm
import os

def cmap_xmap(function,cmap):
    """ Applies function, on the indices of colormap cmap. Beware, function
    should map the [0, 1] segment to itself, or you are in for surprises.
    
    See also cmap_xmap.
    """
    cdict = cmap._segmentdata
    function_to_map = lambda x : (function(x[0]), x[1], x[2])
    for key in ('red','green','blue'):
        cdict[key] = map(function_to_map, cdict[key])
        cdict[key].sort()
        assert (cdict[key][0]<0 or cdict[key][-1]>1),\
            "Resulting indices extend out of the [0, 1] segment."

    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)


class toric_analysis:
    """Class to encapsulate tools used for toric analysis.
    Works with scipy 0.10.0 and numpy 1.6
    """


    def __init__ (self, toric_name='TORICLH.cdf', mode='LH',
        idebug=False, comment='', layout='poster', path="./"):
        import socket
        from time import gmtime

        self.toric_name=toric_name
        self.mode = mode
        self.__version__ = 1.0
        self.idebug = idebug

        self.mylw=1.0
        self.mypt=18.0
        self.fsc=2.0
        self.fw='bold'
        self.set_layout(layout)

        self.path=path

        self.label = True
        self.equigs = {}
        self.toricdict={}


##Open the toric netcdf file read only
        print('reading toric ouput file')
        try:
            self.cdf_hdl = nc.netcdf_file(self.toric_name,'r')
        except IOError:
            print('CRITICAL: ',self.toric_name,' not found.')
            self.cdf_hdl = -1

        try:
            self.qlde_hdl = nc.netcdf_file(path+"toric_qlde.cdf",'r')
        except IOError:
            #print('Non-CRITICAL: ',path+"toric_qlde.cdf",' not found.')
            self.qlde_hdl = -1

#marker sequence
        self.markers = ['o','1','2','s','*','+','H','x']
        self.ls = ['-','--','-.',':','--','-.',':']
        self.bw = False

        print('finished initialization')
        return

    def close (self):
        try:
            self.cdf_hdl.close()
        except IOError:
            print('CRITICAL: ',self.toric_name,' not found.')

        return

    def info( self ):
        "Prints a list of the contents of present TORIC3D output files"

        if (self.cdf_hdl != -1):
            for hdl in [self.cdf_hdl]:
                print ('The toric file, ',self.toric_name,', contains:')
                print ('----------------------------------------------')
                print ("The global attributes: ",hdl.dimensions.keys())        
                print ("File contains the variables: ", hdl.variables.keys())

        if (self.qlde_hdl != -1):
            for hdl in [self.qlde_hdl]:
                print ('The toric file, ',self.toric_name,', contains:')
                print ('----------------------------------------------')
                print ("The global attributes: ",hdl.dimensions.keys())        
                print ("File contains the variables: ", hdl.variables.keys())


        print ('----------------------------------------------')
        
        return

    def plotb0( self, ir=30, db=0 ):
        """Contour plot of bounce averaged Dql coefficient, dB0 """
        if (self.qlde_hdl == -1):
            print ("qlde File not found")
            return


        dqlpsi=self.qlde_hdl.variables['Psi']
        dql_LD=self.qlde_hdl.variables['Qldce_LD']
        nuperp=self.qlde_hdl.dimensions['VelPrpDim']
        nupar=self.qlde_hdl.dimensions['VelParDim']

        umax=1
        umin=-1
        upar=np.arange(nupar)/float(nupar-1)*(umax-umin)+umin
        uperp=np.arange((nupar-1)/2)/float((nupar-1)/2)
        vx,vz=np.meshgrid(uperp,upar)

        fig=plt.figure(figsize=(2.*8.3,2.*3.7))
        plt.axes().set_aspect(1, 'box')
        #plot passing trapped boundary
        if (db>0):
            vpar=np.sqrt(db)*umax
            plt.plot([0,vpar],[0,umax],'k',[0,-vpar],[0,umax],'k',linewidth=2)

        #cd=plt.contourf(vz,vx,np.transpose(np.log(dql_LD[ir,:,:])/np.log(10)),10)
        cd=plt.contourf(vz,vx,np.transpose(np.log10(dql_LD[ir,:,:])+1),10)
        plt.gca().set_ylim(0,umax)
        cbar=plt.colorbar(cd)
        sroa = str(dqlpsi[ir])[0:4]
        plt.title(r'log10 $\lambda$<B> at r/a='+sroa,size=30)
        plt.ylabel(r'$v_{\bot0}/v_{te}$',size=20)
        plt.xlabel(r'$v_{||0}/v_{te}$',size=20)
        plt.draw() #make sure ylimits are updated in displayed plot
        plt.savefig('B0_'+str(ir)+'.eps',format='eps')    

    def plotpower( self, xaxis='tpsi', species=0 ):
        "Plot power profiles versus specified radius for all, or listed species."

        if (self.mode[:2]=='LH'):
            self.__plot1D(xaxis,'S_eld','Power absorbed on electrons')
            plt.xlabel(r'$\sqrt{\psi_{pol}}$')

        cf=plt.gcf()
        cf.subplots_adjust(bottom=0.14)

        return

    def psiplot( self, y ):
        "Plot versus rhopsi. Returns handle on line to modify line style if desired using setp."
        if (self.mode[:2]=='LH'):
            line=self.__plot1D('tpsi',y)
        else:
            line=self.__plot1D('rhopol',y)

        plt.xlabel(r'$\sqrt{\psi_{pol}}$')
        plt.ylabel(y)
        return line

    def plot_1Dfield( self, yvar ):
        "Field versus midplane specified."

        if (self.mode[:2]=='LH'):
            line=self.__plot1D('xeqpl',yvar,'Wave Field component on the midplane',yvar)
            plt.xlabel(r'$X[cm]$')
        return line

    def __plot1D( self, xvar, yvar, ptitle='', plabel='' ):
        "Internal 1D plot"

        x=self.cdf_hdl.variables[xvar].data
        y=self.cdf_hdl.variables[yvar].data
        if (self.mode[:2]=='LH'):
            xname=''
            yname=plabel
        #else:
            #xname=x.attributes.get('long_name','')
            #yname=y.attributes.get('long_name',plabel)
        if (np.size(y) > np.size(x)):
            print ("ToricTools.__plot1D resizing",yvar)
            y=np.array(y)[0:np.size(x)]

        line=plt.plot(x,y)
        plt.title(ptitle)
        if (self.mode[:2]=='LH'):
            plt.xlabel(xname)
            plt.ylabel(yname)

        else:
            plt.xlabel(xvar)
            plt.ylabel(yvar)
        return line

    def __getvar__( self, name ):
        """Internal function to retrieve variable from data file with checking.
        """
        
        try:
            value=self.cdf_hdl.variables[name].data
        except NameError:
            print ('CRITICAL: variable not found')
            raise Exception('Bad variable name in getvar: %s' % name)

        return value


    def fft( self, component='undef',maxr=1. ):

        if (self.mode[:2]=='LH'):
            radius = 'psime'
        if (component=='undef'):
            if (self.mode[:2]=='LH'):
                component='E2d_z_re'
            else:
                component='Re2Ezeta'

        field = self.__getvar__(component)
        rad   = self.__getvar__(radius)
        #field taken to be 2D with shape (ntheta,npsi)
        ntt=field.shape[0]
        nelm=int(field.shape[1]*maxr)
        nlevels=100
        levels=np.arange(nelm/nlevels,nelm-1,nelm/nlevels)
        fftfield = np.zeros((ntt,levels.shape[0]),'complex128')
        i=0
        for ir in levels:
            ffield = (ft.fft(field[:,ir]))
            fftfield[:,i] = ffield
            i=i+1


        return fftfield



    def spectrum( self, component='undef',maxr=1.,cx=0,levels=-1 ):
        """Calculate poloidal spectrum of two dimensional field component.
        """

        if (self.mode[:2]=='LH'):
            radius = 'psime'
        else:
            radius ='PoyFlx_abscissa'
            
        if (component=='undef'):
            if (self.mode[:2]=='LH'):
                component='E2d_z_re'
                componenti='E2d_z_im'
            else:
                component='Re2Eplus'#'Re2Ezeta'
                componenti='Im2Eplus'#'Im2Ezeta'

        f=plt.figure(figsize=(8,6))

        if (component=="power"):
            field = self.get_power2D()
        else:
            field = (self.__getvar__(component))#[:,:]

        if (cx==1):
            fieldi = (self.__getvar__(componenti))#[:,:]
            field=np.array(field)+np.complex(0.,1.)*np.array(fieldi)

        #print(component)
        rad   = self.__getvar__(radius)
        #field taken to be 2D with shape (ntheta,npsi)
        field=field+1.e-20
        ntt=field.shape[0]
        if(self.mode[:]=='LH'):
            nelm=int(field.shape[1]*maxr)
        else:
            nelm=self.cdf_hdl.dimensions['NpsiCurDim']+1
            
        if (levels==-1):
            nlevels=10
            levels=np.arange(nlevels)*int(nelm/nlevels)
        else:
            levels=(np.array(levels)*nelm).astype(int)

        levels=levels[1:nlevels]
        rlevels=rad[levels]
        th = np.arange(ntt)-ntt/2

        ymax = 0.
        ymin = 0.

        colors = cm.viridis_r(np.linspace(0,1,len(levels)))
        i=0
        for ir in levels:  #fft in python isn't normalized to N
            ffield = ft.fftshift(np.log10(abs(ft.fft(field[:,ir]))/float(ntt)+1.e-20))
            ymax = np.max( [ymax, np.max(ffield)] )
            ymin = np.min( [ymin, np.min(ffield)] )
            plabel='%5.2f' % rad[ir]
            if self.bw:
                plt.plot( th, ffield, label=plabel, linestyle=self.ls[i],color='k')
                i=i+1
            else:
                plt.plot( th, ffield, label=plabel,color=colors[i] )
                i=i+1


        ffield = ft.fftshift(np.log10(abs(ft.fft(field[:,nelm-1]))/float(ntt)+1.e-20))
        ymax = np.max( [ymax, np.max(ffield)] )
        ymin = np.min( [ymin, np.min(ffield)] )
#plot antenna spectrum
        plabel='1.0 (ant)'
        #print ("range, levels", rlevels)
        #print ("ymax", ymax,ymin)
        plt.plot( th, ffield,  label=plabel, color='grey',linewidth=2)
        cf=plt.gcf()
        cf.subplots_adjust(left=0.17,right=0.7,top=0.9,bottom=0.15,hspace=0.0)
        plt.axis ('tight')
        plt.axis( xmin=-ntt/4, xmax=ntt/4)
        ticks=np.arange(-ntt/4,ntt/4+1,ntt/16)
        plt.xticks(ticks,fontsize=11)
        plt.axis( ymin=-10)
        lgnd = plt.legend(loc=(1.05,0),title="$\sqrt{\psi_n}$",fontsize=14)
        plt.setp(lgnd.get_title(),fontsize=14)
        plt.xlabel('m')
        plt.ylabel('log10($E_m$)')
        plt.draw()
        return

    def spectrum2( self, component='undef',maxr=1.,cx=0,levels=-1 ):
        """Calculate poloidal spectrum of two dimensional field component.
        """

        if (self.mode[:2]=='LH'):
            radius = 'psime'
        if (component=='undef'):
            if (self.mode[:2]=='LH'):
                component='E2d_z_re'
                componenti='E2d_z_im'
            else:
                component='Re2Ezeta'

        f=plt.figure()

        if (component=="power"):
            field = self.get_power2D()
        else:
            field = (self.__getvar__(component))#[:,:]

        if (cx==1):
            fieldi = (self.__getvar__(componenti))#[:,:]
            field=np.array(field)+np.complex(0.,1.)*np.array(fieldi)

        rad   = self.__getvar__(radius)
        #field taken to be 2D with shape (ntheta,npsi)
        field=field+1.e-20
        ntt=field.shape[0]
        nelm=int(field.shape[1]*maxr)
        if (levels==-1):
            nlevels=30
#            levels=np.arange(nelm/nlevels,nelm-1,nelm/nlevels)
            levels=np.arange(nlevels)*nelm/nlevels
        else:
            levels=(np.array(levels)*nelm).astype(int)

        levels=levels[1:nlevels]
#        levels=nelm-np.arange(1,20,2)
        rlevels=rad[levels]
        th = np.arange(ntt)-ntt/2

        ymax = 0.
        ymin = 0.

        i=0
        m_avg=[]
        for ir in levels:  #fft in python isn't normalized to N
            ffield_factor = ft.fftshift(np.square(abs(ft.fft(field[:,ir])))/sum(np.square(abs(ft.fft(field[:,ir])))))
            #print ("ffield_factor", ffield_factor)
            #print ("th", th)
            m_avg_ir=0
            for i_th in th:
                m_avg_ir=m_avg_ir+th[i_th]*ffield_factor[i_th]
            m_avg.append(m_avg_ir)

#plot antenna spectrum
        plabel='ant'
        #print ("range, levels", rlevels)
        #print ("m_avg", m_avg)
        plt.plot(rlevels,m_avg)
        cf=plt.gcf()
        plt.xlabel('r/a')
        plt.ylabel('average m')
        plt.title('Average of poloidal mode number in terms of radius')
        plt.draw()
        return

    def set_layout( self, layout='poster' ):

        if (layout == 'paper'):
            self.mylw=2.0
            self.mypt=10.0
            self.fsc=1.0
            self.fw='normal'

        params = {
            'axes.linewidth': self.mylw,
            'lines.linewidth': self.mylw,
            'axes.labelsize': self.mypt,
            'font.size': self.mypt,
            'legend.fontsize': self.mypt,
            #'title.fontsize': self.mypt+2.0,
            'xtick.labelsize':self.mypt,
            'ytick.labelsize':self.mypt,
            'font.weight'  : self.fw
            }
        plt.rcParams.update(params)

        return


#note that if plot commands are in the toplevel, they will not return
#to the prompt, but wait to be killed.
    def plot_2Dfield(self, component='E2d_z',logl=0,xunits=1.0,axis=(0.0,0.0),
                     im=False, scaletop=1.0, fig='undef'):
        """
    
        example of using netcdf python modules to plot toric solutions
        requires numpy and matplotlib and netcdf modules for python.

        Note, under windows you need netcdf.dll installed in SYSTEM32 and the file
        system cannot follow symbolic links.  The DLL needs to have executable
        permissions.

        To overplot with limiter, made from efit plotter:
        R.plot_2Dfield(component='Im2Eplus',logl=20,xunits=0.01,axis=maxis,fig=fig1)

        Easier is to plot solution first, then overplot limiter, scaled appropriately:
        p.plot ( rlim*100.-maxis[0], zlim*100.-maxis[1], 'k', linewidth = 2 )
        
        """

        R0=axis[0]
        Z0=axis[1]
        #what should colorbar with be? format=4.1e means 8 characters
        #the bar and title of the bar add about 4 characters.
        #there are 72.27 pt/in
        #12 characters * self.mypt /72.27 pt/in = #in
        legend_frac=12*self.mypt/72.27
        title=component
        if (self.mode[:2]=='LH'):
            xx  = self.cdf_hdl.variables['x_plasma'].data
            yy  = self.cdf_hdl.variables['z_plasma'].data
            if (im):
                im_e2dname=component+'_im'
                title='|'+component+'|'
                component=component+'_re'
            else:
                component=component+'_re'
        else:
            xx  = self.cdf_hdl.variables['Xplasma'].data
            yy  = self.cdf_hdl.variables['Zplasma'].data
            
            if(im):
                pltcomp = [component,component]
                pltcomp[0]=component+'_re'
                pltcomp[1]=component+'_im'
                title='|'+component+'|'
                for i in range(2):
                    if (pltcomp[i]=='Ezeta_re'):
                        pltcomp[i]='Re2Ezeta'
                    if (pltcomp[i]=='Ezeta_im'):
                        pltcomp[i]='Im2Ezeta'
                    if (pltcomp[i]=='Eplus_re'):
                        pltcomp[i]='Re2Eplus'
                    if (pltcomp[i]=='Eplus_im'):
                        pltcomp[i]='Im2Eplus'
            else:
                if (component=='Ezeta_re'):
                    component='Re2Ezeta'
                if (component=='Ezeta_im'):
                    component='Im2Ezeta'
                if (component=='Eplus_re'):
                    component='Re2Eplus'
                if (component=='Eplus_im'):
                    component='Im2Eplus'
    
        if (im and (self.mode[:2]=='LH')):
            e2d = (self.cdf_hdl.variables[component]).data
            im_e2d=(self.cdf_hdl.variables[im_e2dname]).data 
            e2d = abs(e2d+np.complex(0.,1.)*im_e2d)
        elif(im):
            re_e2d = (self.cdf_hdl.variables[pltcomp[0]]).data
            im_e2d = (self.cdf_hdl.variables[pltcomp[0]]).data
            e2d = abs(re_e2d+np.complex(0.,1.)*im_e2d)    
        elif(component=='TDPwE'):
            component='TDPwE'
            e2d = (self.cdf_hdl.variables[component]).data
        elif(component=='TDPwIF'):
            component='TDPwIF'
            e2d = (self.cdf_hdl.variables[component]).data[:,:,-1] #usually minority species
        else:
            e2d = (self.cdf_hdl.variables[component]).data
        #print ("2D Matrix shape:", np.shape(xx))


    #contour with 3 args is confused unless arrays are indexed slices
    #need wrapper to close periodicity in theta direction for this
    #tricky, array indexing different from ncdf slicing
    #this step is needed because periodic dimension is not closed.
    #i.e. its [0,pi) not [0,pi]
        dd=np.shape(xx)
        sx=dd[0]
        sy=dd[1]

        xxx=np.zeros((sx+1,sy),'d')
        xxx[0:sx,:]=xx[:,:]
        xxx[sx,:]=xx[0,:]
        yyy=np.zeros((sx+1,sy),'d')
        yyy[0:sx,:]=yy[:,:]
        yyy[sx,:]=yy[0,:]

        xxx=(xxx+R0)*xunits
        yyy=(yyy+Z0)*xunits

        ee2d=np.zeros((sx+1,sy),'d')
        ee2d[0:sx,:]=e2d[:,:]
        ee2d[sx,:]=e2d[0,:]
        
        emax=ee2d.ravel()[ee2d.argmax()]
        emin=ee2d.ravel()[ee2d.argmin()]

        #contouring levels
        rmax=max([abs(emax),abs(emin)])
        rmin=min([0.,emax,emin])
        val=np.arange(-rmax*1.1,rmax*1.1,(rmax+rmax)/25.,'d')
        #print ("values",val,'xx',rmax,rmin)
        if (im):
            val=np.arange(0.,rmax*1.1,(rmax)/24.,'d')
            #print ("values",val)

        #finally, make the plot
        cwidth=xxx.max()-xxx.min()
        cheight=yyy.max()-yyy.min()
        asp=cheight/cwidth
        #print ("asp:", asp)

        #leave space for bar
        if (fig=='undef'):
            fig=plt.figure(figsize=(self.fsc*3.0+legend_frac,3.0*self.fsc*asp))
            fig.subplots_adjust(left=0.02,bottom=0.15,top=0.90)

        sax=plt.axes().set_aspect(1, 'box') # the right way to control aspect ratio

        maxpsi=xxx.shape[1]-1
        plt.plot(xxx[:,maxpsi],yyy[:,maxpsi],'k-')

        plt.ioff()
        if (logl <= 0):
            if(im or component=='TDPwIF' or component=='TDPwE'):
                CS=plt.contourf(xxx,yyy,ee2d,25,cmap=cm.jet)
            else:
                CS=plt.contourf(xxx,yyy,ee2d,25,cmap=cm.RdBu_r)
        if (logl > 0):
            lee2d=np.log((ee2d)+0.1)/np.log(10)
            if(im):
                CS=plt.contourf(xxx,yyy,lee2d,logl,cmap=cm.jet)
            else:
                CS=plt.contourf(xxx,yyy,lee2d,logl,cmap=cm.RdBu_r)

        cbar=plt.colorbar(CS,format='%4.1e',ax=sax)#,ticks=cbar_tics)

        if(self.mode[:]=='LH'):
            cbar.ax.set_ylabel('$\log_{10}$[V/m]',rotation=270,labelpad=25)

        plt.xlabel('X [cm]')
        plt.ylabel('Z [cm]')
        plt.title(title,fontsize=self.mypt+2.0, loc='left')

        return


#Handling equigs file
    def __get_varname(self, f):
        "Reads next line from file f and returns it, optionally printing it."
        varname=f.readline()
        if self.idebug:
            print (f.name,varname)
        return varname

    
    def read_equigs(self, equigsfile='equigs.data'):
        "Read the equilibrium file created by toric in toricmode='equil',isol=0."
        if self.idebug:
            print ("Using ", equigsfile)
        equigs_hdl=file(equigsfile,'r')

        varname = self.__get_varname(equigs_hdl)
        self.equigs["rtorm"] = np.fromfile(equigs_hdl,sep=" ",count=1,dtype=float)

        varname = self.__get_varname(equigs_hdl)
        self.equigs["raxis"]= np.fromfile(equigs_hdl,sep=" ",count=1,dtype=float)

        varname = self.__get_varname(equigs_hdl)
        self.equigs["bzero"] = np.fromfile(equigs_hdl,sep=" ",count=1,dtype=float)

        varname = self.__get_varname(equigs_hdl)
        self.equigs["torcur"]= np.fromfile(equigs_hdl,sep=" ",count=1,dtype=float)

        varname = self.__get_varname(equigs_hdl)
        self.equigs["imom"] = np.fromfile(equigs_hdl,sep=" ",count=1,dtype=int)
        imom = self.equigs["imom"]

        varname = self.__get_varname(equigs_hdl)
        self.equigs["nmhd"] = np.fromfile(equigs_hdl,sep=" ",count=1,dtype=int)
        nmhd = self.equigs["nmhd"]

        varname = self.__get_varname(equigs_hdl)
        self.equigs["srad"] = np.fromfile(equigs_hdl,sep=" ",count=nmhd,dtype=float)

#this needs to be reshaped or remapped into the R,Z sin cos arrays toric uses
        varname = self.__get_varname(equigs_hdl)
        self.equigs["rzmcs2d"] = np.fromfile(equigs_hdl,sep=" ",
                                                count=2*nmhd+4*nmhd*imom,dtype=float)

        varname = self.__get_varname(equigs_hdl)
        self.equigs["qqf"] = np.fromfile(equigs_hdl,sep=" ",count=nmhd,dtype=float)

#logic checking for "END"
        varname = self.__get_varname(equigs_hdl)
        self.equigs["jcurr"] = np.fromfile(equigs_hdl,sep=" ",count=nmhd,dtype=float)

        varname = self.__get_varname(equigs_hdl)
        self.equigs["gcov"]= np.fromfile(equigs_hdl,sep=" ",count=nmhd,dtype=float)

        varname = self.__get_varname(equigs_hdl)
        self.equigs["rhotor"] = np.fromfile(equigs_hdl,sep=" ",count=nmhd,dtype=float)

        varname = self.__get_varname(equigs_hdl)
        self.equigs["lastpsi"] = np.fromfile(equigs_hdl,sep=" ",count=1,dtype=float)

        equigs_hdl.close()


### user routines using the above, could be in a different module
    def powpoynt( self ):
        fig = plt.figure(figsize=(8,6))
        ax1 = fig.add_subplot(111)
        line1,=self.psiplot('S_eld')
#can use setp(lines, ) to change plot properties.
        plt.setp(line1,color='b',marker='+',label='seld')
        ax1.set_ylabel('Power_e',color='b')

        ax2 = ax1.twinx()
        line2,=self.psiplot('vpoynt')
#set axis floor at 0
        plt.gca().set_ylim(0)
#change color and symbol
        plt.setp(line2,color='r',marker='.',label='vpoynt')
        ax2.set_ylabel('Poynting',color='r')
#make  legend too
        plt.legend( (line1,line2), (r'$P_{eld}$','<ExB>'),loc=2 )
        fig.subplots_adjust(left=0.1,bottom=0.12,top=0.95,right=0.8,hspace=0.32)
        plt.draw()
        return fig

    def powpoynt_ICRF(self,min_indx=-1):
        fig = plt.figure(figsize=(8,6))
        ax1 = fig.add_subplot(111)
        line_e = self.__plot1D('Pw_abscissa','PwE')
        line_eibw = self.__plot1D('Pw_abscissa','PwEIBW')
        
        plt.setp(line_e,color='b',label='e- FW')
        plt.setp(line_eibw,color='b',linestyle='--',label='e- IBW')
        
        r = self.cdf_hdl.variables['Pw_abscissa'].data
        nspec = self.cdf_hdl.dimensions['SpecDim']
        fund = self.cdf_hdl.variables['PwIF'].data
        harm = self.cdf_hdl.variables['PwIH'].data
        colors = ['g','y','m']
        label1 = 'fund min'
        label2 = 'harm bulk'
        line_min = ax1.plot(r,fund[:,min_indx],color=colors[-1],label=label1)
        line_bulk = ax1.plot(r,harm[:,0]+harm[:,1],color=colors[0],linestyle='--',label=label2)

        ax1.set_ylabel('$P$ [MW/m$^3$]')
        ax1.set_xlabel('r/a')
        ax1.set_xlim(0,np.max(1.05*np.max(r)))

        ax2 = ax1.twinx()
        line_pvec=self.__plot1D('PoyFlx_abscissa','PoyFlx')

        #set axis floor at 0
        plt.gca().set_ylim(0)
        ax1.set_ylim(0,1.3*np.max(fund[3:len(r),min_indx]))
        
        #change color and symbol
        plt.setp(line_pvec,color='r',label='vpoynt')
        ax2.set_ylabel('Poynting',color='r')

        #make legend too
        r_patch = mpatches.Patch(color='r', label='<ExB>')
        m_patch = mpatches.Patch(color='m', label='$P_{min}$')
        b_patch = mpatches.Patch(color='b', label='$P_{e}$ (FW/IBW)')
        g_patch = mpatches.Patch(color='g', label='$P_{bulk}$') 
        plt.legend(handles=[r_patch,m_patch,b_patch,g_patch],fontsize=11,loc='upper left')
        plt.tight_layout()
        plt.draw()
        return fig

    def xpsi_map( self ):
        """Return map of x(theta=0)/x(psi=1,theta=0) versus psipol."""
        xmap=1.0

        return xmap

    def get_power2D( self ):
#figure out a sed way of cutting these lines into the file.
#also need to replace '-0.' with ' -0.'
#sed -n -e '/elec/,/,/p' filename | sed -e '/-0\./ -0./g' > torica_2dpower.sol
        try:
            toricsol = file('torica_2dpower.sol','r')
        except IOError:
            print ('CRITICAL: torica_2dpower.sol not found.')
            return

#skip title and max value
        toricsol.readline()
        toricsol.readline()
        nt=self.cdf_hdl.dimensions['ntt']
        nr=self.cdf_hdl.dimensions['mptpsi']

        power=np.fromfile(toricsol,sep=" ",count=nt*nr,dtype=float)
        toricsol.close()

        power=np.transpose(np.reshape(power,(nr,nt)))
        return power

#    def plotsE2D_ICRF( self, prefix='' ):
#  
#        
#        #self.spectrum2(cx=1)
#	#plt.draw()
#        #plt.savefig(prefix+'spectrum_avg.eps',format='eps')
#
#        fig=plt.figure(figsize=(2.*8.3,2.*3.7))
#        self.plot_2Dfield(component='Re2Eplus',logl=0, scaletop=0.3)
#	    plt.draw()
#        plt.savefig(prefix+'E2Dplus.eps',format='eps')
#        
#        fig=plt.figure(figsize=(2.*8.3,2.*3.7))
#        self.plot_2Dfield(component='Re2Eminus',logl=0, scaletop=0.3)
#	    plt.draw()
#        plt.savefig(prefix+'E2Dminus.eps',format='eps')
#
#        fig=plt.figure(figsize=(2.*8.3,2.*3.7))
#        self.plot_2Dfield(component='Re2Ezeta',logl=0, scaletop=0.3)
#	    plt.draw()
#        plt.savefig(prefix+'E2Dzeta.eps',format='eps')
        

        
    def threeplots( self, prefix='',min_indx=-1):
        """Makes and saves the three most commonly used plots. Plots are saved in 
        the current directory. An optional prefix can be used to label them or change
        the save path.
        * Power and poynting flux on one plot as eps.
        * The polodial power spectrum on six flux surfaces for convergence as eps.
        * And the 2D parallel electric field contour plot as a png."""

        print('plotting spectrum')
        self.spectrum(cx=1)
        plt.tight_layout()
        plt.draw()
        plt.savefig(prefix+'spectrum.eps',format='eps',bbox_inches='tight')
        
        #self.spectrum2(cx=1)
	#plt.draw()
        #plt.savefig(prefix+'spectrum_avg.eps',format='eps')

        print('plotting fields')
        if (self.mode[:]=='ICRF'):
            self.plot_2Dfield(component='Ezeta',logl=25,im=True)
            plt.tight_layout()
            plt.draw()
            plt.savefig(prefix+'log10Ezeta2d.eps',format='eps',bbox_inches='tight')
            self.plot_2Dfield(component='Eplus',logl=25,im=True)
            plt.tight_layout()
            plt.draw()
            plt.savefig(prefix+'log10Eplus2d.eps',format='eps',bbox_inches='tight')
            self.plot_2Dfield(component='TDPwIF')
            plt.tight_layout()
            plt.draw()
            plt.savefig(prefix+'log10Pmin2d.eps',format='eps',bbox_inches='tight')
            self.plot_2Dfield(component='TDPwE')
            plt.tight_layout()
            plt.draw()
            plt.savefig(prefix+'log10Pe2d.eps',format='eps',bbox_inches='tight')
        else:
            self.plot_2Dfield(im=True,logl=25)
        
        print('plotting power profiles')
        if (self.mode[:]=='ICRF'):
            self.powpoynt_ICRF(min_indx = min_indx)
        else:
            self.powpoynt()
        plt.savefig(prefix+'powerpoynt.eps',format='eps')

        return


####main block
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import toric_tools
    import sys
    import getopt


# get file name if provided
    iprefix=""
    ifile="TORICLH.cdf"
    try:
        opts, args = getopt.getopt(argv, "hp:f:",["help","prefix=","file="])
    except getopt.GetoptError:
        print ("Accepted flags are help and prefix=")
        sys.exit(2)
	
    for opt, arg in opts:
        if opt in ("-p","--prefix"):
            iprefix=arg
        elif opt in ("-f","--file"):
            ifile=arg
	    
#Load a run
    LHRun=toric_tools.toric_analysis(toric_name=ifile)

    LHRun.threeplots(prefix=iprefix)

#make sequence of plots, ala the old toric idl driver as an option.



