!
!
      program prepare_input

      USE plasma_state_mod

!.......................................................................
!     Prepares new cqlinput_new namelist file for cql3d, using
!     several files from the cql3d source, and the old cqlinput
!     file.
!
!     Compile and load line:
!     Include files in current dir:
!                             param.h name.h frname.h name_decl.h frname_decl.h
!     lf95 --fix -o prepare_input prepare_cql3d_input.f90 aindfpa.f aindflt.f eqindflt.f urfindfl.f frinitl.f
!     lf95 --fix -o prepare_input prepare_cql3d_input.f aindfpa.f aindflt.f eqindflt.f urfindfl.f frinitl.f
!     ifort --fix -o prepare_input prepare_cql3d_input.f aindfpa.f aindflt.f eqindflt.f urfindfl.f frinitl.f
!
!     090924: Running on franklin with Pathscale compiler.
!             Updated from PS1 to PS2 (no ps%nspec,ps%nspec_nonMax).
!             iarg() function replaced by command_argument_count().
!.......................................................................

       use swim_global_data_mod, only :
     1 rspec, ispec,     ! int: kind specification for real (64-bit)
                         ! and integer (32-bit).
!            & swim_string_length,  ! length of strings for names, files, etc.
     1 swim_error         ! error routine
    
!      implicit none
      implicit integer (i-n), real*8 (a-h,o-z)

!.......................................................................
!     NAMELISTS AND THEIR STORAGE, from CQL3D FILES
!.......................................................................

      include 'param.h'
      include 'name_decl.h'
      include 'frname_decl.h'
      include 'name.h'
      include 'frname.h'

!  Temporary dimensioning (BH)
      REAL(KIND=rspec) ns(1:njenea,1:ntotala),Ts(1:njenea,1:ntotala)

      REAL(KIND=rspec) Zeff_local(1:njenea)
      
      REAL(KIND=rspec) delta_t
      character(len = 128) my_mnemonic, delta_t_string, nsteps_string

!.......................................................................
!     set defaults using cql3d subroutines in
!     aindfpa.f, aindflt.f, eqindflt.f, urfindfl, frinitl.f
!.......................................................................

!     ---------------------------------
!     set default values of input data for orginal cql3d namelist:

c     F2003-syntax: command_argument_count()/get_command_argument(,)
      iargs=command_argument_count()
      write(*,*)'iargs=',iargs
      if(iargs .ne. 3) then
        print*, 'prepare_cql3d_input usage: '
        print*, 'nsteps  delta_t  mnemonic_name'
        stop 'incorrect command line arguments'
      end if
      call getarg(1,nsteps_string)
      call getarg(2,delta_t_string)
      call getarg(3,my_mnemonic)
      
      print*, 'prepare_cql3d_input Command line arguments: ',
     . 	nsteps_string,  delta_t_string, my_mnemonic

      call aindfpa
      call aindflt
      call eqindflt
      call urfindfl
      call frinitl

!.......................................................................
!     Read in new namelist data from cql3d namelist input file cqlinput
!.......................................................................

      open(unit=2, file = 'cqlinput', status = 'old',
     1    form = 'formatted', iostat = ierr)
      if(ierr .ne. 0) stop 'cannot open cqlinput cql3d namelist file'

      read(2,setup0,iostat=istat)
      if (istat.ne.0) then
         write(*,*)
         write(*,*)
     .   'No setup0?: Check have new setup0/setup namelist structure.'
         write(*,*)
         STOP
      endif
      read(2,setup)     !Note:  there are two namelist setup sections
      read(2,trsetup)
      read(2,sousetup)
      read(2,eqsetup)
      read(2,rfsetup)
      read(2,frsetup)
      close(2)
      
      read(nsteps_string, '(i10)') nsteps
      nstop=nsteps     !Reset these variables in namelist
      nplot=nsteps
      nplt3d=nsteps
      read(delta_t_string, '(e15.7)') delta_t
      dtr=delta_t
      mnemonic = trim(my_mnemonic)
      rffile= '/home/bobh/cql3d/bonoli/iter_lh/genray/iter5_060828/iter5_060828.nc'

      write(*,*)'prepare_cql3d_input:  finished reading cqlinput nl'


!.......................................................................
!     Make available all variables in the plasma state module,
!     specified in swim_state_spec.dat (or corresponding .txt) in
!     components/state/src/plasma_state.  See also McCune's Plasma_State_Vnn
!     design doc at http://www.cswim.org/componentdescrip/.
!     These variable are available with the prepended "component selector"
!     for the current state, ps%, or for the previous state,  %psp,
!     and a few others.
!.......................................................................

      call ps_get_plasma_state(ierr)
      if(ierr .ne. 0) stop 'cannot get plasma state needed for profiles '


!.......................................................................
!     In general, as many namelist parameters as possible will be
!     specified only through the cqlinput file.  Depending on these
!     settings, additional data will be obtained from the plasma state.
!.......................................................................


!  Input radial profiles of plasma data on normalized radial grid [0. to 1.].
!  Check that cql3d radial coord is compatible with PS:
      if (radcoord.ne.'sqtorflx') stop 'prepare_cql3d_input: check radcoord'

!  Profiles will be input to cql3d using radial arrays (for splines).
!  If cqlinput does not indicate 'spline' input, then don't modify the profile.

      ispline=0   !spline indicator
      if (iprone.eq.'spline' .or. iprote.eq.'spline'
     1 .or. iproti.eq.'spline'.or. iprozeff.eq.'spline'
     1 .or. iprovphi.eq.'spline' .or. ipronn.eq.'spline'
     1 .or. iproelec.eq.'spline') ispline=1
      write(*,*)'prepare_cql3d_input: ispline=',ispline
	
      if (ispline.eq.1) then     !at least 1 profile to be reset

         if (ps%nrho.gt.njenea) stop 'ps%nrho.gt.njenea'
         njene=ps%nrho     !number of bin boundaries
!        Setup cqlinput radial array
         ryain(1:njene)= ps%rho(1:ps%nrho)
         write(*,*)'prepare_cql3d_input: ryain = ', njene,ryain(1:njene)

         if (iprone.eq.'spline') then
!           Setup electron species
!           ntotal= total number of species in cqlinput, 
!                General FP'd(ngen)+Maxwl(nmax)
	    ntotal=ngen+nmax
!           Step through the species to determine which are electrons
            write(*,*)'prepare_cql3d_input: ntotal,fmass(1:ntotal)=',
     1                ntotal,fmass(1:ntotal)
            do k=1,ntotal   !Possibly, electron density given for general and
                         !Maxwellian species.
	    if (fmass(k).le.1.e-27) then    !electron mass = 9.1e-28 grams
	       call zone_check(ps%nrho,enein(1,k),ps%ns(1,0))
               enein(1:njene,k)=1.e-6*enein(1:njene,k)  !convert to cgs
               write(*,*)'prepare_cql3d_input: k,njene,enein(1:njene,k)=',
     1                   k,njene,enein(1:njene,k)
	    endif
!           If iprozeff.eq."parabola" .or. iprozeff.eq."spline", then no
!           ion species profiles will be required (charges given in cqlinput).
!  ps%nspec and ps%nspec_nonMax are only in Plasma State 1:
!            write(*,*)'prepare_cql3d_input: ps%nspec_th,ps%nspec,ps%nspec_nonMax = ',
!     1                ps%nspec_th,ps%nspec,ps%nspec_nonMax
            write(*,*)'prepare_cql3d_input: ps%nspec_th = ',
     1                ps%nspec_th
            if (iprozeff.eq.'disabled') then    !Set up ion profiles
!               write(*,*)'iprozeff.eq.disabled not yet set up'
!               STOP
	    endif
         enddo
         endif
	       

      endif

!     Use zone_check to interpolate ps% zone centered data
!     onto zone bndries for cql3d
!      call zone_check(ps%nrho,ns(1,0),ps%ns(1,0))
      call zone_check(ps%nrho,Ts(1,0),ps%Ts(1,0))
      call zone_check(ps%nrho,Ts(1,1),ps%Ts(1,1))
      call zone_check(ps%nrho,Zeff_local,ps%Zeff(1))

      write(*,*)'ps%Zeff = ',ps%Zeff


!      write(*,*)'prepare_cql3d_input: njene,ns(,0)  = ',
!     1                           njene,ns(1:njene,0) 
      write(*,*)'prepare_cql3d_input: njene,Ts(,0)  = ',
     1                           njene,Ts(1:njene,0)
      write(*,*)'prepare_cql3d_input: njene,Ts(,1)  = ',
     1                           njene,Ts(1:njene,1)

!!*****************************************************************
!!     TEMPORARY, Since the plasma_state I have is giving zeff=0.
!      write(*,*)
!      write(*,*)'NOTE: REMOVE TEMPORARY set of zeff in prep_cql3d_in.'
!      Zeff_local(1:ps%nrho)=2.000
!      write(*,*)
!!*****************************************************************

      write(*,*)'prepare_cql3d_input: njene,Zeff_local  = ',
     1                           njene,Zeff_local(1:njene)


!  This is setup for  DC electric field driven current:
!  Specify ene,zeff,elecin,Te,Ti profiles.   Use the first ion temperature
!  from the PS.   CQL3D will prepare two ion species, one hydrogenic
!  and a higher Z ion to achieve quasineutrality.

!     Use plasma-state magnetic axis, ps%R_axis, to convert V_loop to 
!     cql3d toroidal electric field (volts/cm):
      elecin(1:njene)=ps%V_loop(1:ps%nrho)/(8*atan2(1.,1.)*ps%R_axis*100.)    

      write(*,*)'prepare_cql3d_input: njene,elecin  = ',
     1                           njene,elecin(1:njene)


      ngen=1   !  For now, one FPd species

!      enein(1:njene,1)=ns(1:njene,0)
      write(*,*) 'HERE1'
      tein(1:njene)=Ts(1:njene,0)
      write(*,*) 'HERE2'
      tiin(1:njene)=Ts(1:njene,1)
      write(*,*) 'HERE3'
      zeffin(1:njene)=Zeff_local(1:njene)
      write(*,*) 'HERE4'

!  Set prefix for output files to default,  mnemonic:
!      mnemonic="mnemonic"






!$$$      ! now the species
!$$$
!$$$
!$$$      !first identfy the species, mass, charge
!$$$      i = 1
!$$$      if( i .le. ps%nspec_th) then
!$$$      	amu1 = anint(ps%m_s(i) / ps_mp)
!$$$	z1 = anint(ps%q_s(i) / ps_xe)
!$$$	i = i + 1
!$$$      end if
!$$$      if( i .le. ps%nspec_th) then
!$$$	amu2 = anint(ps%m_s(i) / ps_mp)
!$$$	z2 = anint(ps%q_s(i) / ps_xe)
!$$$	i = i + 1
!$$$      end if
!$$$      if( i .le. ps%nspec_th) then
!$$$      	amu3 = anint(ps%m_s(i)/ ps_mp)
!$$$	z3 = anint(ps%q_s(i) / ps_xe)
!$$$	i = i + 1
!$$$      end if
!$$$      if( i .le. ps%nspec_th) then
!$$$      	amu4 = anint(ps%m_s(i) / ps_mp)
!$$$	z5 = anint(ps%q_s(i) / ps_xe)
!$$$	i = i + 1
!$$$      end if
!$$$      if( i .le. ps%nspec_th) then
!$$$      	amu1 = anint(ps%m_s(i) / ps_mp)
!$$$	z5 = anint(ps%q_s(i) / ps_xe)
!$$$	i = i + 1
!$$$      end if
!$$$      if( i .le. ps%nspec_th) then
!$$$      	amu1 = anint(ps%m_s(i) / ps_mp)
!$$$	z6 = anint(ps%q_s(i) / ps_xe)
!$$$	i = i + 1
!$$$      end if
!$$$
!$$$
!$$$      !
!$$$      s_s_name(0:ps%nspec_th) = (ps%s_name(0:ps%nspec_th))
!$$$      s_m_s(0:ps%nspec_th) = ps%m_s(0:ps%nspec_th)
!$$$      s_q_s(0:ps%nspec_th) = ps%q_s(0:ps%nspec_th)
!$$$
!$$$
!$$$
!$$$
!$$$      !now deal with values of T and n on grid
!$$$      if(.true.)  then !putting numbers on the whole grid
!$$$      	do i = 0, ps%nspec_th
!$$$      	   call  zone_check(ps%nrho, s_n_s(1:ps%nrho,i), ps%ns(:,i))
!$$$	   call  zone_check(ps%nrho, s_t_s(1:ps%nrho,i), ps%ts(:,i))
!$$$
!$$$      	end do
!$$$        s_rho_n_grid(1:s_nrho_n) = ps%rho  !grids are the same for both n and T  FIX
!$$$        s_rho_t_grid(1:s_nrho_t) = ps%rho
!$$$	s_nrho_n = s_nrho_n - 1  ! missing rho = 1.0--don't let aorsa read it
!$$$	s_nrho_t = s_nrho_t - 1
!$$$      else !putting numbers on the half grid
!$$$	s_rho_n_grid(1:s_nrho_n-1) = (ps%rho(2:s_nrho_n) + ps%rho(1:s_nrho_n-1)) / 2.0  !grids are the same for both n and T  FIX
!$$$        s_rho_t_grid(1:s_nrho_t-1) = (ps%rho(2:s_nrho_n) + ps%rho(1:s_nrho_n-1)) / 2.0
!$$$	s_nrho_n = s_nrho_n - 1  ! one less zone-center
!$$$	s_nrho_t = s_nrho_t - 1
!$$$	s_n_s(1:s_nrho_n,0:ps%nspec_th) = ps%ns
!$$$        s_t_s(1:s_nrho_n,0:ps%nspec_th) = ps%ts
!$$$      end if
!$$$
!$$$      iprofile = 5  !use numerical profiles
!$$$      eqdsk = trim(ps%eqdsk_file)
!$$$


      open(unit=64, file='cqlinput_new', delim = 'apostrophe',
     1         status = 'unknown', form = 'formatted')

      write(64, nml = setup0)
      write(64, nml = setup)
      write(64, nml = trsetup)
      write(64, nml = sousetup)
      write(64, nml = eqsetup)
      write(64, nml = rfsetup)
      write(64, nml = frsetup)

      close(64)

      end

          
      subroutine zone_check(nrho, x_out, x_in)
    
!     unravel zone points to boundary points

      use swim_global_data_mod, only :
     1 rspec, ispec,     ! int: kind specification for real and integer
!            & swim_string_length,  ! length of strings for names, files, etc.
     1 swim_error         ! error routine


      integer, intent(in) :: nrho    ! # of pts covering [0:1]
      REAL(KIND=rspec), intent(in) :: x_in(nrho-1)  ! zone values
      REAL(KIND=rspec), intent(out) :: x_out(nrho)  ! boundary values
      
      integer :: irho
      
      x_out(nrho) = 0.5 * (3.0 * x_in(nrho -1) - x_in(nrho-2))
      do irho = nrho - 1, 1, -1
         x_out(irho) = 2.0 * x_in(irho) - x_out(irho + 1)
      end do
      
      end subroutine zone_check
!      end program prepare_input
      

