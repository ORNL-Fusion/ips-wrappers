! SF cleaned up 2023. I have removed a lot of functionality here
! this was mostly because the routine was in really rough shape
! and many features were long deprecated. I am now working on
! overhauling this into something that is a little more robust.
! Current focus is on maintaining E- functionality while also working
! with minority heating and minority+bulk ICRF. May add beams w/
! FREYA too if BH wants to help out with that.
! Older verison is still availablle in prepare_clq3d_input.f90_old

      program prepare_cql3d_input


      ! The following USE gives access to the extensive plasma state data.
      ! (For genray, it is necessary to rename "output" in the plasma state to
      ! avoid conflict with genray namelist /output/. For cql3d, it is benign.)
      USE plasma_state_mod, ps_output_local => output
      
      USE swim_global_data_mod, only : 
     1 rspec, ispec,          ! kind specification for real (64-bit)
     1 swim_error,            ! error routine
     1 swim_filename_length   ! length of filename string

      use strings
      
      implicit integer (i-n), real*8 (a-h,o-z) ! SF ...kill me now

      ! program parameters
      !----------------------------------------------------------------
      integer, parameter :: swim_string_length = 256  !for compatibility LAB
      REAL(KIND=rspec), parameter :: pi=3.1415926535897931_rspec
      REAL(KIND=rspec), parameter :: twopi=6.2831853071795862_rspec
      REAL(KIND=rspec), parameter :: me=9.1094e-31_rspec 
      REAL(KIND=rspec), parameter :: mp=1.6726e-27_rspec
      REAL(KIND=rspec), parameter :: xe=1.6022e-19_rspec 
      REAL(KIND=rspec), parameter :: zero=0.0_rspec 
      REAL(KIND=rspec), parameter :: one=1.0_rspec

      ! CQL3D NAMELISTS AND THEIR STORAGE
      !----------------------------------------------------------------
      include 'param.h'
      include 'name_decl.h'
      include 'frname_decl.h'
      include 'name.h'
      include 'frname.h'

      ! program variables
      !----------------------------------------------------------------
      INTEGER :: inp_unit,ierr

      integer nj,njp1,nrho,l
      integer isp, ipairspec, indx_loop
      REAL(KIND=rspec) :: dryain,rfmin_chargemass,spec_chargemass

      integer, dimension(:), allocatable :: incl_species !indicator for
                                          !species to send to genray.
      integer iargs, isp_min, isp_fus, ncustom, indxmin

      integer, dimension(ngena) :: alla2custom
      
      !characters
      CHARACTER(len=16) :: base_rdc_str
      character(len=2) :: istr
      CHARACTER(len=16) :: rfmistr, fusnstr,state_var
      
      !Unified plasma parameter arrays (in ps%nspec_alla list)
      !abridged thermal species + all fast species combined.
      REAL(KIND=rspec), dimension(:), allocatable :: rho_bdy_rezon
      REAL(KIND=rspec), dimension(:), allocatable :: dense,tempe,zeff
      REAL(KIND=rspec), dimension(:,:), allocatable :: densi,tempi
      REAL(KIND=rspec), dimension(:), allocatable :: wk_eperp,wk_epll,
     +                                               wk_dens
      REAL(KIND=rspec), dimension(:), allocatable :: rho_user
      real(kind=rspec), dimension(:), allocatable :: rya_grid
      integer :: cclist(ps_ccount)  ! State partial copy controls.
                                    ! Ps_ccount is number of known
                                    ! components in the PS.
      
      !Setup some allocatable variables from PS, which will be
      !checked if allocated:
      REAL(KIND=rspec),dimension(:),allocatable :: nbeami,nmini,nfusi

      ! program namelist variables
      !----------------------------------------------------------------
      character(len =swim_string_length) ::
     + cur_state_file = "ips-state.nc"
      CHARACTER(len=16) :: cql3d_specs='e'
      CHARACTER(len=16) :: cql3d_mode='LH'
      CHARACTER(len=16) :: rf_code='toric'
      CHARACTER(len=16) :: restart='disabled'
      INTEGER :: nsteps = 1, norf = 0
      REAL(KIND=rspec) :: arg_deltat = 0.001_rspec
      REAL(KIND=rspec) :: arg_enorm = 1000.0_rspec
      INTEGER :: arg_nsurfFP
      REAL(KIND=rspec) :: arg_rhoFPlo = 0.01_rspec
      REAL(KIND=rspec) :: arg_rhoFPhi = 0.95_rspec
      REAL(KIND=rspec), dimension(ngena) :: pscale = 1.0_rspec
      CHARACTER(len=16), dimension(ngena) :: spec_list = 'NONE'
      
      namelist /cql3d_prepare_nml/
     +    cur_state_file, cql3d_specs, 
     +    cql3d_mode, rf_code, restart,
     +    nsteps, arg_deltat, arg_enorm, arg_nsurfFP,
     +    arg_rhoFPlo, arg_rhoFPhi, pscale, norf, spec_list

      ! ---Begin Program Logic----------
      
      !read program input namelist
      !----------------------------------------------------------------
      OPEN (unit=21, file="cql3d_prepare.nml", 
     +      status='old',form = 'formatted', iostat=ierr)
      IF (ierr .ne. 0) THEN
         CALL SWIM_error ('open', 'prepare_cql3d_input.f90',
     +                    'cql3d_prepare.nml')
         WRITE (*,*) 'prepare_cql3d_input.f90: Cannot open ',
     +    TRIM('cql3d_prepare.nml')
         call exit(1)
      END IF
      READ(21,cql3d_prepare_nml)
      CLOSE (21)
      WRITE (*, cql3d_prepare_nml)

      !trim text inputs
      cur_state_file=trim(cur_state_file)
      cql3d_specs=trim(cql3d_specs)
      cql3d_mode=trim(cql3d_mode)
      rf_code=trim(rf_code)
      restart=trim(restart)     
      
      ! Create default input data
      !----------------------------------------------------------------
      call aindfpa
      call aindflt
      call eqindflt
      call urfindfl
      call frinitl

      ! some new options that are annoying to deal with as they can cause random
      ! run failures if set arbitrarily. these do not have defaults right now
      ! in cql3d
      !read_data = 'disabled'
      !ipxy = 1
      !jpxy = 1
      
      ! Read the cql3d namelists
      !-----------------------------------------------------------------------
      open(unit=22, file = 'cqlinput', delim='apostrophe',
     1     status = 'old', form = 'formatted', iostat = ierr)
      if(ierr .ne. 0) stop 'cannot open cqlinput cql3d namelist file'

      read(22,setup0,iostat=istat)
      if (istat.ne.0) then
         write(*,*)
         write(*,*)
     1   'No setup0?: Check have new setup0/setup namelist structure.'
         write(*,*)
         STOP 1
      endif
      read(22,setup)     !Note:  Now only one namelist setup section
      read(22,trsetup)
      read(22,sousetup)
      read(22,eqsetup)
      read(22,rfsetup)
      read(22,frsetup)
      close(22)

      ! Get the plasma state
      !-----------------------------------------------------------------------
      call ps_get_plasma_state(ierr,cur_state_file)
      if(ierr .ne. 0) stop
     + 'cannot get plasma state needed for profiles'

      ! Reset CQL3D values based on inputs and PS data
      ! ---------------------------------------------------------------
      if (restart .eq. 'enabled') then
         nlrestrt='ncdfdist'
         nlwritf='ncdfdist'
      else
         nlrestrt='disabled'
         nlwritf='ncdfdist'
      endif

      ! Make sure cql3d namelist density and temp scaling turned off,
      ! except introduce any required unit changes here.
      enescal=1.d-6
      tescal=1.d0  !cqlinput is in keV.
      tiscal=1.d0
      zeffscal=1.d0
      vphiscal=ps%r_axis*1.d2  !Change to cm/sec from omega (rad/sec)
      ennscal=1.d-6   !Change to /cm**3
      eqdskin='eqdsk'
      radcoord='sqtorflx'
      CFP_INTEGRALS='disabled'
      pwrscale(:) = 1.d0
      

      !SF saved for later if I ever want to add in beams
!      if (ps%nbeam.ne.1) then
!           write(*,*)'STOP:  Need to adjust prepare_cql3d_input'
!           write(*,*)'       and PS for nbeam.gt.1'
!            stop
!            endif
!            bptor(1)=ps%power_nbi(1)


      !set timestepping
      nstop = nsteps
      dtr = arg_deltat

      !set enorm
      enorm = arg_enorm
      
      !set up rya grid
      lrz = arg_nsurfFP
      lz = lrz
      allocate(rya_grid(arg_nsurfFP+1))
      rya_grid=0.0_rspec
      dryain = (arg_rhoFPhi-arg_rhoFPlo)/(arg_nsurfFP-1)
      do l=1,arg_nsurfFP
         rya_grid(l+1)=arg_rhoFPlo + dryain*(l-1)
      enddo
      WRITE(*,*) dryain, rya_grid
      rya(0:arg_nsurfFP) = rya_grid 
      deallocate(rya_grid)
      
      !interpolate cql3d T/n profiles from plasma state data
      !----------------------------------------------------------------
      
      nj=101    !local variable
      njp1=nj+1   !local variable
      if (njene.ne.101) then
         write(*,*)'Resetting njene for uniform, 101-pt plasma profs'
         njene=101
         if (njene.gt.njenea) then
            write(*,*)'Stop:  Need njenea bigger in param.h'
            stop
         endif
      endif

      dryain=1.0_rspec/(nj-1)
      do l=1,nj
         ryain(l)=(l-1)*dryain
      enddo
      allocate(rho_bdy_rezon(njp1))
      rho_bdy_rezon(1)=ryain(1)
      rho_bdy_rezon(njp1)=ryain(nj)
      rho_bdy_rezon(2:nj) =
     +     .5_rspec*(ryain(1:nj-1)+ryain(2:nj))

      if (cql3d_specs.eq.'E') then

         !set species
         nmax = ps%nspec_tha+1
         ngen = 1

         !set profiles to splines of dataset
         iprone = 'spline'
         iprote = 'spline'
         iproti = 'spline'
         iprozeff = 'spline'
         izeff = 'backgrnd'
         iproelec = 'spline'

         !set up RF diffusion stuff
         if((rf_code.eq.'toric').or.(rf_code.eq.'aorsa'))then
            if(norf.ne.1)then
               rdcmod = 'format1'
            else
               rdcmod = 'disabled'
            endif
            rdcfile(1) = 'du0u0_input'
            nrdc = 1
            nrdcspecies = 1
            rdc_netcdf = 'disabled'
         elseif(rf_code.eq.'genray')then
            rdcmod = 'disabled'
            !SF finish later
         endif
         
         !electron data
         isp_min = ps%rfmin_to_alla(1)
         kspeci(1,1) = 'e'
         kspeci(2,1) = 'general'
         fmass(1)    = me*1.d+3
         bnumb(1)    = -1.d0
         call ps_rho_rezone(rho_bdy_rezon, ps%id_ns(0), enein(1:nj,1),
     +                         ierr, zonesmoo=.TRUE.)

         !loop over ion species if impurity species in plasma state
         !file skip (we are lumping impurities in this case based on Zeff)
         indx_loop = 2
         do isp = 1,ps%nspec_tha
            if(ps%sa_type(isp).ne.2)then !check if impurity species 
               kspeci(1,indx_loop) = trim(ps%alla_name(isp))
               kspeci(2,indx_loop) = 'maxwell'
               fmass(indx_loop) = ps%m_alla(isp)*1.d+3
               bnumb(indx_loop) = NINT(ps%q_alla(isp)/ps_xe)
               call ps_rho_rezone(rho_bdy_rezon, ps%id_ns(isp), enein(1:nj,indx_loop),
     +                         ierr, zonesmoo=.TRUE.)
               indx_loop=indx_loop+1
            else
               nmax= nmax-1
            endif
         enddo

         !set last species as electron bulk
         kspeci(1,nmax+1) = 'e'
         kspeci(2,nmax+1) = 'maxwell'
         fmass(nmax+1)    = me*1.d+3
         bnumb(nmax+1)    = -1.d0
         call ps_rho_rezone(rho_bdy_rezon, ps%id_ns(0), enein(1:nj,nmax),
     +                         ierr, zonesmoo=.TRUE.)
         
         !skip fusion ash and minority species (very small correction in most cases)
         !and CQL3D does not support damping on alphas

         DO i = 1,nmax+ngen
            WRITE(*,*) trim(kspeci(1,i)),' ', trim(kspeci(2,i)), fmass(i), bnumb(i)
         ENDDO

         !ion and electron temperature profiles
         call ps_rho_rezone(rho_bdy_rezon, ps%id_Ts(0), tein(1:nj), ierr,
     +     zonesmoo=.TRUE.)
         call ckerr('ps_rho_rezone (U2)',ierr)
         
         call ps_rho_rezone(rho_bdy_rezon, ps%id_Ts(1), tiin(1:nj), ierr,
     +     zonesmoo=.TRUE.)
         call ckerr('ps_rho_rezone (U2)',ierr)

         !zeff profiles
         call ps_rho_rezone(rho_bdy_rezon, ps%id_Zeff, zeffin(1:nj), ierr,
     +     zonesmoo=.TRUE.)
         call ckerr('ps_rho_rezone (U2)',ierr)

         do l=1, ps%nrho-1
            if (zeffin(l).lt.1.0) then
               write(*,*)'Problem: PS gives zeffin.lt.1: l, zeffin(l)',
     +                l, zeffin(l)
               stop
            endif
         enddo

      elseif (cql3d_specs .eq.'MIN') then

         if((rf_code.ne.'toric').and.(rf_code.ne.'aorsa'))then
         WRITE(*,*)
     +   'Wrong RF code for FPed species needed toric got:', rf_code 
         stop
         endif 
         
         !set species
         nmax=ps%nspec_alla+1 !all ion species+e-
         ngen=1

         !set profiles to splines of dataset
         iprone = 'spline'
         iprote = 'spline'
         iproti = 'spline'
         iprozeff = 'disabled'
         
         !set number of diff coeffs
         if(norf.ne.1)then
            rdcmod = 'format1'
         else
            rdcmod = 'disabled'
         endif
         rdcfile(1) = 'du0u0_input'
         nrdc = 1
         nrdcspecies = 1
         rdc_netcdf = 'disabled'
         rdcscale(:) = 1.d0
         rdcscale(1) = pscale(1)
         
         !minority species data
         isp_min = ps%rfmin_to_alla(1)
         kspeci(1,1) = trim(ps%alla_name(isp_min))
         kspeci(2,1) = 'general'
         fmass(1)    = ps%m_alla(isp_min)*1.d+3
         bnumb(1)    = NINT(ps%q_alla(isp_min)/ps_xe)
         call ps_rho_rezone(rho_bdy_rezon, ps%id_nmini(1),enein(1:nj,1),
     +                      ierr, zonesmoo=.TRUE.)

         isp_min = ps%rfmin_to_alla(1)
         kspeci(1,2) = trim(ps%alla_name(isp_min))
         kspeci(2,2) = 'maxwell'
         fmass(2)    = ps%m_alla(isp_min)*1.d+3
         bnumb(2)    = NINT(ps%q_alla(isp_min)/ps_xe)
         call ps_rho_rezone(rho_bdy_rezon, ps%id_nmini(1),enein(1:nj,2),
     +                      ierr, zonesmoo=.TRUE.)
     
         !ion species data
         indx_loop=3
         do isp = 1,ps%nspec_tha
            kspeci(1,indx_loop) = trim(ps%alla_name(isp))
            kspeci(2,indx_loop) = 'maxwell'
            fmass(indx_loop) = ps%m_alla(isp)*1.d+3
            bnumb(indx_loop) = NINT(ps%q_alla(isp)/ps_xe)
            call ps_rho_rezone(rho_bdy_rezon, ps%id_ns(isp), enein(1:nj,indx_loop),
     +                         ierr, zonesmoo=.TRUE.)
            indx_loop = indx_loop+1
         enddo

         !fusion species data. not going to be the correct energy, but
         !including for quasineutrality/dilution purposes
         if (ps%nspec_fusion>0)then
            isp_fus = ps%sfus_to_alla(1)
            kspeci(1,indx_loop) = trim(ps%alla_name(isp_fus))
            kspeci(2,indx_loop) = 'maxwell'
            fmass(indx_loop)    = ps%m_alla(isp_fus)*1.d+3
            bnumb(indx_loop)    = NINT(ps%q_alla(isp_fus)/ps_xe)
            call ps_rho_rezone(rho_bdy_rezon, ps%id_nfusi(1), enein(1:nj,indx_loop),
     +                         ierr, zonesmoo=.TRUE.)
            indx_loop = indx_loop+1      
         endif
  
         !electron species data
         kspeci(1,indx_loop) = 'e'
         kspeci(2,indx_loop) = 'maxwell'
         fmass(indx_loop) = me*1.d+3
         bnumb(indx_loop) = -1

         DO i = 1,nmax+ngen
            WRITE(*,*) trim(kspeci(1,i)),' ', trim(kspeci(2,i)), fmass(i), bnumb(i)
         ENDDO

         call ps_rho_rezone(rho_bdy_rezon, ps%id_ns(0), enein(1:nj,indx_loop),
     +                         ierr, zonesmoo=.TRUE.)

         call ps_rho_rezone(rho_bdy_rezon, ps%id_Ts(0), tein(1:nj), ierr,
     +     zonesmoo=.TRUE.)
         call ckerr('ps_rho_rezone (U2)',ierr)
         
         call ps_rho_rezone(rho_bdy_rezon, ps%id_Ts(1), tiin(1:nj), ierr,
     +     zonesmoo=.TRUE.)
         call ckerr('ps_rho_rezone (U2)',ierr)

         !zeff profile
         call ps_rho_rezone(rho_bdy_rezon, ps%id_Zeff, zeffin(1:nj), ierr,
     +                      zonesmoo=.TRUE.)
         call ckerr('ps_rho_rezone (U2)',ierr)
     
      else if (cql3d_specs .eq.'MIN+') then
         write(*,*) 'MODE IS MIN+'
         if((rf_code.ne.'toric').and.(rf_code.ne.'aorsa'))then
         WRITE(*,*)
     +   'Wrong RF code for FPed species needed toric got:', rf_code 
         stop
         endif 

         WRITE(*,*) 'nspec_alla', ps%nspec_alla
         nmax=ps%nspec_alla+1
         ngen=2

         !set profiles to splines of dataset
         iprone = 'spline'
         iprote = 'spline'
         iproti = 'spline'
         iprozeff = 'disabled' !must be disabled right now
         
         !set number of diff coeffs
         if(norf.ne.1)then
            rdcmod = 'format1'
         else
            rdcmod = 'disabled'
         endif
         rdcfile(1) = 'du0u0_input_1'
         rdcfile(2) = 'du0u0_input_2'
         nrdc = 2
         nrdcspecies(1) = 1
         nrdcspecies(2) = 2
         rdc_netcdf = 'disabled'
         rdcscale(:) = 1.d0
         rdcscale(1) = pscale(1)
         rdcscale(2) = pscale(2)
         
         !minority species data
         isp_min = ps%rfmin_to_alla(1)
         kspeci(1,1) = trim(ps%alla_name(isp_min))
         kspeci(2,1) = 'general'
         fmass(1)    = ps%m_alla(isp_min)*1.d+3
         bnumb(1)    = NINT(ps%q_alla(isp_min)/ps_xe)
         call ps_rho_rezone(rho_bdy_rezon, ps%id_nmini(1),enein(1:nj,1),
     +                      ierr, zonesmoo=.TRUE.)

     
         !find minorities "pair"
         ipairspec = 0
         rfmin_chargemass = (1.0_rspec*nint(ps%q_ALLA(isp_min)/ps_xe))
     +                    /(1.0_rspec*nint(ps%m_alla(isp_min)/ps_mp))    
         do isp = 1, ps%nspec_tha
            spec_chargemass = (2.0_rspec*nint(ps%q_ALLA(isp)/ps_xe))
     +                    /(1.0_rspec*nint(ps%m_alla(isp)/ps_mp))
            if((abs(rfmin_chargemass - spec_chargemass) < 1.0E-6).and.(nint(ps%q_ALLA(isp)/ps_xe)<6))then
               ipairspec = isp
            endif
         enddo
         
         !pair ion species data
         kspeci(1,2) = trim(ps%alla_name(ipairspec))
         kspeci(2,2) = 'general'
         fmass(2)    = ps%m_alla(ipairspec)*1.d+3
         bnumb(2)    = NINT(ps%q_alla(ipairspec)/ps_xe)
         call ps_rho_rezone(rho_bdy_rezon, ps%id_ns(ipairspec),enein(1:nj,2),
     +                      ierr, zonesmoo=.TRUE.)

         isp_min = ps%rfmin_to_alla(1)
         kspeci(1,3) = trim(ps%alla_name(isp_min))
         kspeci(2,3) = 'maxwell'
         fmass(3)    = ps%m_alla(isp_min)*1.d+3
         bnumb(3)    = NINT(ps%q_alla(isp_min)/ps_xe)
         call ps_rho_rezone(rho_bdy_rezon, ps%id_nmini(1),enein(1:nj,3),
     +                      ierr, zonesmoo=.TRUE.)
     
         indx_loop=4
         do isp = 1,ps%nspec_tha
            kspeci(1,indx_loop) = trim(ps%alla_name(isp))
            kspeci(2,indx_loop) = 'maxwell'
            fmass(indx_loop) = ps%m_alla(isp)*1.d+3
            bnumb(indx_loop) = NINT(ps%q_alla(isp)/ps_xe)
            call ps_rho_rezone(rho_bdy_rezon, ps%id_ns(isp), 
     +              enein(1:nj,indx_loop),ierr, zonesmoo=.TRUE.)
            indx_loop=indx_loop+1
         enddo

         !fusion species data. not going to be the correct energy, but
         !including for quasineutrality/dilution purposes
         if (ps%nspec_fusion>0)then
            isp_fus = ps%sfus_to_alla(1)
            kspeci(1,indx_loop) = trim(ps%alla_name(isp_fus))
            kspeci(2,indx_loop) = 'maxwell'
            fmass(indx_loop)    = ps%m_alla(isp_fus)*1.d+3
            bnumb(indx_loop)    = NINT(ps%q_alla(isp_fus)/ps_xe)
            call ps_rho_rezone(rho_bdy_rezon, ps%id_nfusi(1), enein(1:nj,indx_loop),
     +                         ierr, zonesmoo=.TRUE.)
            indx_loop = indx_loop+1      
         endif
  
         !electron species data
         kspeci(1,indx_loop) = 'e'
         kspeci(2,indx_loop) = 'maxwell'
         fmass(indx_loop) = me*1.d+3
         bnumb(indx_loop) = -1

         DO i = 1,nmax+ngen
            WRITE(*,*) trim(kspeci(1,i)),' ', trim(kspeci(2,i)), fmass(i), bnumb(i)
         ENDDO

         call ps_rho_rezone(rho_bdy_rezon, ps%id_ns(0), enein(1:nj,indx_loop),
     +                         ierr, zonesmoo=.TRUE.)

         call ps_rho_rezone(rho_bdy_rezon, ps%id_Ts(0), tein(1:nj), ierr,
     +     zonesmoo=.TRUE.)
         call ckerr('ps_rho_rezone (U2)',ierr)
         
         call ps_rho_rezone(rho_bdy_rezon, ps%id_Ts(1), tiin(1:nj), ierr,
     +     zonesmoo=.TRUE.)
         call ckerr('ps_rho_rezone (U2)',ierr)

         !zeff profile
         call ps_rho_rezone(rho_bdy_rezon, ps%id_Zeff, zeffin(1:nj), ierr,
     +                      zonesmoo=.TRUE.)
         call ckerr('ps_rho_rezone (U2)',ierr)

      !the custom case that can assign an arbitrary species to the general species
      !first it will loop over all the species in the plasma state and find their
      !corresponding alla index and link them to their general species index specified
      !in cql3d_specs by their ordering   
      elseif(cql3d_specs .eq.'CUSTOM')then
         write(*,*) 'MODE IS CUSTOM'
         write(*,*) 'SPECIES TO BE EVOLVED: ',spec_list
         
         !determine the number of species in the custom list
         ncustom = 0
         indxmin = 0
         indxe = 0
         DO i=1,ngena
            if (spec_list(i) .ne. 'NONE') then
               !convert list to upper case and trim for string compare purposes later
               spec_list(i) = trim(to_upper(spec_list(i)))
               ncustom = ncustom + 1
               if ((spec_list(i) .eq. 'H') .or. (spec_list(i) .eq. 'HE3')) then
                  indxmin = ncustom
               elseif (spec_list(i) .eq. 'E') then
                  indxe = ncustom
               endif
            end if
         END DO

         if (ncustom.eq.0) then
            WRITE(*,*) 'No species found in speclist'
            stop
         endif

         !set up plasma profiles and diffusion
         iprone = 'spline'
         iprote = 'spline'
         iproti = 'spline'
         iprozeff = 'disabled'
         
         !set up general species
         WRITE(*,*) 'nspec_alla', ps%nspec_alla
         nspec_allap1 = ps%nspec_alla+1
         nmax=nspec_allap1
         ngen=ncustom
         if (ngen.gt.nmax)then
            WRITE(*,*) 'more species set to be evolved than exist in state'
            WRITE(*,*) 'nspec in state : ', nspec_allap1
            WRITE(*,*) 'ncustom : ', ncustom
            stop
         endif
         
         !preset some diffusion quantities
         base_rdc_str = 'du0u0_input_'
         nrdc = ncustom
         rdcscale(:) = 1.d0
         if((rf_code.eq.'toric').or.(rf_code.eq.'aorsa'))then
            rdc_netcdf = 'disabled'
            if (norf.ne.1) then
               rdcmod = 'format1'
            else
               rdcmod = 'disabled'
            endif
         elseif(rf_code.eq.'genray')then
            rdcmod = 'disabled'
            !SF finish later
         elseif(norf.ne.1)then
            WRITE(*,*)
     + 'unsupported RF code expected (genray, aorsa, toric) got: ',rf_code
            stop
         endif

         ! loop over general species
         isp_min=-1
         do i=1,NCUSTOM
            ! find the corresponding species in the plasma state
            ! note that the 
            rfmistr = "_rfmi"
            fusnstr = "_fusn"
            isp=-1
            do j=0,ps%nspec_alla
               state_var = to_upper(trim(ps%alla_name(j)))
               if (state_var.eq.spec_list(i))then
                  isp=j
               elseif (state_var.eq.(spec_list(i)//rfmistr))then
                  isp=j
                  isp_min = isp
               elseif (state_var.eq.(spec_list(i)//fusnstr))then
                  isp=j
                  ips_fus = isp
               end if
            enddo
            if (isp.eq.-1) then
               WRITE(*,*) 'Species not found in plasma state: ', spec_list(i)
               stop
            endif
            
            ! set up species 
            kspeci(1,i) = trim(ps%alla_name(isp))
            kspeci(2,i) = 'general'
            fmass(i) = ps%m_alla(isp_min)*1.d+3
            bnumb(i) = NINT(ps%q_alla(isp_min)/ps_xe)

            ! profiles (assumes only one fusion species and minority species)
            if (isp_min.eq.isp)then
               call ps_rho_rezone(rho_bdy_rezon, ps%id_nmini(1),enein(1:nj,i),
     +                      ierr, zonesmoo=.TRUE.)
            elseif (isp_fus.eq.isp)then
               call ps_rho_rezone(rho_bdy_rezon, ps%id_nfusi(1), enein(1:nj,i),
     +                         ierr, zonesmoo=.TRUE.)
            else
               call ps_rho_rezone(rho_bdy_rezon, ps%id_ns(isp),enein(1:nj,i),
     +                      ierr, zonesmoo=.TRUE.)
            end if
            call ckerr('ps_rho_rezone (U2)',ierr)
            
            ! set up rf diffusion for species
            if((rf_code.eq.'toric').or.(rf_code.eq.'aorsa'))then
               write(istr,*) i
               istr=adjustl(istr)
               rdcfile(i) = trim(base_rdc_str)//trim(istr)
               nrdcspecies(i) = i
               rdcscale(i) = pscale(i)
            elseif(rf_code.eq.'genray')then
               !SF finish later   
            endif
         enddo

         !set the maxwellian species
         indx_loop = ncustom
         do i=0,ps%nspec_tha
            kspeci(1,indx_loop) = trim(ps%alla_name(isp))
            kspeci(2,indx_loop) = 'maxwell'
            fmass(indx_loop) = ps%m_alla(isp)*1.d+3
            bnumb(indx_loop) = NINT(ps%q_alla(isp)/ps_xe)
            call ps_rho_rezone(rho_bdy_rezon, ps%id_ns(isp), enein(1:nj,indx_loop),
     +                         ierr, zonesmoo=.TRUE.)
            call ckerr('ps_rho_rezone (U2)',ierr)
            indx_loop = indx_loop+1
         enddo

         !special treatment for rfmin and fusion if they exist
         if (ps%nspec_rfmin>0)then
            isp_min = ps%rfmin_to_alla(1)
            kspeci(1,indx_loop) = trim(ps%alla_name(isp_min))
            kspeci(2,indx_loop) = 'maxwell'
            fmass(indx_loop)    = ps%m_alla(isp_min)*1.d+3
            bnumb(indx_loop)    = NINT(ps%q_alla(isp_min)/ps_xe)
            call ps_rho_rezone(rho_bdy_rezon, ps%id_nmini(1), enein(1:nj,indx_loop),
     +                         ierr, zonesmoo=.TRUE.)
            call ckerr('ps_rho_rezone (U2)',ierr)
            indx_loop = indx_loop+1      
         endif
  
         if (ps%nspec_fusion>0)then
            isp_fus = ps%sfus_to_alla(1)
            kspeci(1,indx_loop) = trim(ps%alla_name(isp_fus))
            kspeci(2,indx_loop) = 'maxwell'
            fmass(indx_loop)    = ps%m_alla(isp_fus)*1.d+3
            bnumb(indx_loop)    = NINT(ps%q_alla(isp_fus)/ps_xe)
            call ps_rho_rezone(rho_bdy_rezon, ps%id_nfusi(1), enein(1:nj,indx_loop),
     +                         ierr, zonesmoo=.TRUE.)
            call ckerr('ps_rho_rezone (U2)',ierr)      
         endif

         !output check for written species
         WRITE(*,*) 'Added following species to CQL3D input:'
         do i=1,nmax + ngen
             WRITE(*,*) trim(kspeci(1,i)),' ', trim(kspeci(2,i)), fmass(i), bnumb(i)
         ENDDO
         
         !temperature profiles (note cql3d does not presently support
         !different ion temperatures across species, this is not hard to fix
         !if willing to modify cql3d)
         call ps_rho_rezone(rho_bdy_rezon, ps%id_Ts(0), tein(1:nj), ierr,
     +     zonesmoo=.TRUE.)
         call ckerr('ps_rho_rezone (U2)',ierr)
         
         call ps_rho_rezone(rho_bdy_rezon, ps%id_Ts(1), tiin(1:nj), ierr,
     +     zonesmoo=.TRUE.)
         call ckerr('ps_rho_rezone (U2)',ierr)

         !zeff profiles
         call ps_rho_rezone(rho_bdy_rezon, ps%id_Zeff, zeffin(1:nj), ierr,
     +                      zonesmoo=.TRUE.)
         call ckerr('ps_rho_rezone (U2)',ierr)

      endif
      
      !loop voltage profiles SF giving error on compile?
      !iproelec = 'spline'
      !call ps_rho_rezone(rho_bdy_rezon, ps%id_V_loop, elecin(1:nj), ierr,
      !+     zonesmoo=.TRUE.)
      !call ckerr('ps_rho_rezone (U2)',ierr)
      !elecin = elecin / (2.*pi*ps%r_axis*1.0d+02)
      
      !write all namelist to  cqlinput file
      !----------------------------------------------------------------

      write(*,*)'before write_all_namelists'

      open(unit=64, file='cqlinput_new', delim = 'apostrophe',
     +         status = 'unknown', form = 'formatted')

      write(64, nml = setup0)
      write(64, nml = setup)
      write(64, nml = trsetup)
      write(64, nml = sousetup)
      write(64, nml = eqsetup)
      write(64, nml = rfsetup)
      write(64, nml = frsetup)

      close(64)
      write(*,*)'after write all namelists'
      write(*,*)'prepare_cql3d_input: Normal end'
      write(*,*)

      end program prepare_cql3d_input

c----------------------------------------------------------------
c helpers
c----------------------------------------------------------------

      SUBROUTINE ckerr(sbrtn,ierr)
      character*(*), intent(in) :: sbrtn
      write(*,*) 'ckerr: ierr=',ierr
      IF(ierr.NE.0) then
         write(6,*) ' ?plasma_state_test: error in call: '//trim(sbrtn)
         stop
      ENDIF
      END SUBROUTINE ckerr
