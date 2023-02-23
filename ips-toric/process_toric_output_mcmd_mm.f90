!IN PROGRESS DO NOT USE 15 NOV 2012 JCW
program process_toric_output_mcmd_mm

  !---------------------------------------------------------------------
  ! processing toric output
  ! 
  !last update on 11/12/12 by creting _mm version to handle multiple 
  !                        toroidal modes
  !last update on 10/12/11 by bonoli to add current profiles
  !     update on 04/02/08 by JCW to move PS setup to do_toric_init and 
  !                        handle rezoning of variable sized power arrays
  !     update on 03/18/08 by JCW to handle remapping of psi_pol->psi_tor
  !John C. Wright 02/28/07 based on L. Berry's process_aorsa_output program
  !Lee A. Berry   12/06/06 based on D. McMcune's test program
  !
  !--------------------------------------------------------------------
  !
 
  ! State elements are traditional f77 integer, floating point, and character
  ! string scalars and arrays all contained within a large container data 
  ! type "plasma_state".  Two instances of this container data type are
  ! declared in the module "plasma_state_mod":
  !
  !   ps -- the current state (timestep now being computed).
  !   psp -- the prior or "committed" state (completed, prior timestep)
  !
  ! Elements of the state can be referenced directly using the standard
  ! f95 syntax-- e.g. ps%nrho is the size of one of the radial flux 
  ! coordinate "rho" grids in the state.
  !
  ! State elements can be directly modified by codes that use the plasma
  ! state-- although this should be done carefully.  Generally, items in 
  ! the state are meant to be shared between multiple components of the 
  ! IPS; conventions for use of the state data will need to evolve.
  !
  ! States will be mapped to NetCDF files.  The module "plasma_state_mod"
  ! defines two variables for these filenames:
  !
  !   CHARACTER*256 :: state_file = 'plasma_state.cdf'
  !   CHARACTER*256 :: prior_file = 'prior_state.cdf'
  !
  ! The assigned default values can of course be modified.
  !-----------------------------------
  ! In addition to direct access to state data elements, the following
  ! subroutines are available (defined in the f95 module plasma_state_mod)
  ! (this is the module's public interface):
  !
  !    integer :: ierr -- status code returned by many routines (0=OK).
  !

  !
  !    SUBROUTINE ps_clear_profs(ierr)
  !       Set all profile data to zero-- this might be desirable when starting
  !       to build a state at a new time.  It will be easier to tell what has
  !       been added, and what not, if quantities not added are zero, rather
  !       than from the prior timestep.  All the prior timestep data is still
  !       accessible in the prior state object psp-- i.e. psp%rho_eq(...), etc.
  !       Scalar data and grids are not affected by this call-- just profiles.
  !
  !    SUBROUTINE ps_update_equilibrium(<g-filename>,ierr)
  !       Update state MHD equilibrium from G-eqdsk file <g-filename>
  !          (could also be an MDSplus G-eqdsk timeslice from an experiment).
  !       Compute state quantities derived from the equilibrium
  !       Arrays such as enclosed volumes (m^3) ps%vol(1:ps%nrho_eq) are 
  !       filled in.
  !
  !    SUBROUTINE ps_store_plasma_state(ierr)
  !       Update interpolation information and store the state (ps) in
  !       the file (state_file).
  !
  !    SUBROUTINE ps_update_plasma_state(ierr)
  !       Update interpolation information but do not write a file.
  !
  !    SUBROUTINE ps_commit_plasma_state(ierr)
  !       Copy current state (ps) to prior state (psp).  Save prior state
  !       in file (prior_state).  The file (state_file) is not modified--
  !       use ps_store_plasma_state for this.
  !
  !    SUBROUTINE ps_get_plasma_state(ierr)
  !       Read the current state (ps) with all interpolation information
  !       from (state_file) --AND-- read the prior state (psp) with all
  !       its interpolation information from the file (prior_state).
  !
  ! Profile IDs:  each state array that is defined over one or more coordinate
  !   grids is assigned an ID variable or array with name ID_<name>.  These IDs
  !   are needed to refer to specific profiles for interpolation or rezoning.
  !
  !   Examples mapping from swim_state_spec.dat to plasma state object "ps":
  !
  !      R|pclin  chi_e(nrho)          ! electron thermal conductivity
  !
  !      -> ps%chi_e(...)  (1d allocated REAL(KIND=rspec) array (1:ps%nrho)) &
  !      -> ps%id_chi_e    INTEGER scalar
  !
  !      R|units=m^-3|step  ns(~nrho,0:nspec_th)       ! thermal specie density
  !
  !      -> ps%ns(...)     (2d allocated REAL(KIND=rspec) array
  !      -> ps%id_ns(...)  (1d allocated INTEGER array)
  !
  !              ps%ns(1:ps%nrho-1,0:ps%nspec_th)
  !              ps%id_ns(0:ps%nspec_th)
  !
  ! Direct interpolation:
  !
  !    SUBROUTINE PS_INTRP_1D(...)  for 1D profiles  LAB based on steps
  !
  !       CALL PS_INTRP_1D( &
  !            x,  &      ! target of interpolation: scalar or 1d vector
  !            id, &      ! profile(s) interpolated: scalar 1d or 2d INT array
  !            ans,&      ! result, dimensioned to match x(...) and id(...)
  !            ierr,   &  ! completion code 0=OK
  !            icur,   &  ! OPTIONAL: "ps_previous" or "ps_current" (default)
  !            ideriv, &  ! OPTIONAL: derivative control for all IDs
  !            ideriv1s,& ! OPTIONAL: derivative control for each ID separately
  !            iccw_th)   ! OPTIONAL: .FALSE. for clockwise poloidal angle
  !                         (theta) coordinate
  !
  !       If x is a vector, the 1st dimension size of ans matches size(x);
  !       subsequent dimensions match sizes of dimensions of id(...).
  !
  !       icur can be used to select the current or prior state; current state
  !       is the default.
  !
  !       ideriv1s is only available if id(...) is an array.  If so, ideriv1s
  !       takes precedence over ideriv, if both are present.  The dimensioning
  !       of INT array ideriv1s must match dimensioning of ID exactly.
  !
  !       If neither ideriv nor ideriv1s are specified, the interpolating 
  !       function value is evaluated with no derivatives.
  !
  !       iccw_th would only be needed if a profile f(theta) is ever defined.
  !
  !    SUBROUTINE PS_INTRP_2D(...)  for 2D profiles
  !
  !       (interface is like PS_INTRP_1D, except that:
  !
  !          x -> x1,x2 -- two interpolation target scalars or vectors must
  !                        be supplied; the coordinate to which they belong
  !                        will match the declaration
  !
  !          similarly, optional derivative control is available separately
  !          for each coordinate:
  !            ideriv -> ideriv1,ideriv2
  !            ideriv1s -> ideriv1s,ideriv2s
  !
  ! Profile rezoning integration:
  !
  !    SUBROUTINE PS_RHO_REZONE(...) for "conservative" rezoning
  !      of profiles f(rho) 
  !      (rezoning to 1d of profiles f(rho,theta) will be done but is
  !      not yet implemented -- DMC 16 Oct 2006).
  !
  !    SUBROUTINE PS_RHOTH_REZONE(...) for "conservative" rezoning
  !      of profiles f(rho,theta) -- not yet implemented DMC 16 Oct 2006
  ! 
  !-----------------------------------
  ! Implementation of interpolation and rezoning-- the "visible" state
  ! elements ps%* and psp%* are not the entire state.  When interpolation
  ! and rezoning operations are carried out, additional hidden information
  ! is accessed which defines the profile interpolation methods.
  !
  ! Interpolation methods for profiles are part of the state specification
  ! "swim_state_spec.dat" and are handled by generated code.
  !-----------------------------------

  ! define state object and interface...

  USE plasma_state_mod
!--------------------------------------------------------------------------

  USE swim_global_data_mod, only : &
            & rspec, ispec, &               ! int: kind specification for real and integer
!            & swim_string_length, &         ! length of strings for names, files, etc.
            & swim_error                    ! error routine
    
  USE netcdf !toric solution file is in nc format
  USE toric_utils_mod    
!--------------------------------------------------------------------------
!
!   Data declarations
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------------
!
!   TORIC data that will be given to the state via the swim_out file
!   the toric output data is in a file called toric.nc
!
!--------------------------------------------------------------------------
   


 
  IMPLICIT NONE
  real(rspec), allocatable :: rhon(:), dvol(:), redotje_int(:),  &
    & redotji_int(:,:), wdote_int(:,:), wdoti_int(:,:,:),tpwi(:,:), &
    & mtrhon(:), mprhon(:), prhon(:), trhon(:), rho(:), &
    & wdot_tmp(:), darea(:), &
    & CdensFW(:), CdensIBW(:), cdicrf_tmp(:), cdicrf_int(:),&
    & nphi_freq(:),ant_nphi_load(:),alpha(:)


  real(rspec):: tpwe(2),tpw, CtotFW, CtotIBW,ant_load
  logical, parameter:: debug=.true.
  logical:: lex
    
  integer, parameter :: swim_string_length = 256  !for compatibility LAB
     


  !------------------------------------
  !  local
  INTEGER :: ierr, n, id_ncdf, imode, nspec, id_temp, i, j, ind_nphi
  integer :: istat, nmhd, inp_unit, tot_num_nphi, nphi_idx, isp
  integer :: nnoderho, iarg
  integer :: isrc,iphi,fidx
  integer, parameter :: stderr=6
  integer, allocatable :: nphi_map(:),nphi_isrc(:),toric_to_alla(:)

  character(len=swim_string_length) :: cur_state_file, ncdfvarname, &
    &                      toric_output_file, nphi_map_file, toric_data_file
  character(len=swim_string_length), allocatable:: nphi_dirname(:)
  character(len=100)::filename

  namelist /ips/ toric_to_alla
  write(*,*) ' -- get (restore) plasma state from file -- '
  ! this call retrieves the plasma state at the present ps%, and
  ! previous time steps psp%
  
  CALL get_arg_count(iarg)
!Args: state_file, nphi_map_directories
  SELECT CASE (iarg)

      case(0)
          cur_state_file="cur_state.cdf"
          nphi_map_file="run_map.txt"

      case(1)
          CALL get_arg(1,cur_state_file)
          nphi_map_file="run_map.txt"

      case(2)
         CALL get_arg(1,cur_state_file)
         CALL get_arg(2,nphi_map_file)

      case DEFAULT
         write(0,*) 'Error. Illegal number of arguments.'
         write(0,*) 'process_toric_output: '
         write(0,*) 'process_toric_output cur_state_file nphi_map_file'
         stop 'incorrect number of command line arguments'

  end select

  write(stderr,*) 'process_toric_output: cur_state_file = ', trim(cur_state_file)
  write(stderr,*) 'process_toric_output: nphi_map_file  = ', trim(nphi_map_file)
  CALL ps_get_plasma_state(ierr, trim(cur_state_file))

  CALL assert( ierr==0,' process toric: ps_get_plasma_state: ierr=',ierr )
  write(*,*) 'stringlen',swim_string_length
  if (allocated(ps%picrf_srcs) .eqv. .TRUE.) then
     if (debug) then
        write(stderr,*) "RF alloc?",allocated(ps%picrf_srcs),size(ps%picrf_srcs)
        write(stderr,*) '   number of icrf sources = ', ps%nicrf_src
     endif
  else
      stop "ERROR, process_toric_output_mcmd_mm: RF sources are not allocated in PS."
  endif

 if (debug) then 
    write(stderr, *) '   number of species(all,thermal) = ', ps%nspec_alla, ps%nspec_th
    if (allocated(ps%rho_icrf) .eqv. .TRUE.) then
       write(stderr,*)'   radial grid points for icrf = ', ps%nrho_icrf
    endif
 endif
 
 imode=NF90_NOWRITE !read only mode(default), dataset id returned in id_ncdf

!Read in nphi mapping file and make directory names
  call getlun(inp_unit,ierr)
  open(unit=inp_unit,file=nphi_map_file,status='old')
  !parse. format is first line in number of nphi values
  !then repeated records of [srcindx,nphi,filename] with each value on a line alone
  read(inp_unit,*) tot_num_nphi
  allocate(character(len=swim_string_length) :: nphi_dirname(tot_num_nphi))
  allocate(nphi_map(tot_num_nphi),nphi_isrc(tot_num_nphi),nphi_freq(tot_num_nphi))
  allocate(toric_to_alla(30))

  do ind_nphi=1,tot_num_nphi
     read(inp_unit,*) nphi_isrc(ind_nphi)
     read(inp_unit,*) nphi_freq(ind_nphi)
     read(inp_unit,*) nphi_idx
     read(inp_unit,*) nphi_map(ind_nphi)
     read(inp_unit,'(A)') nphi_dirname(ind_nphi)
  end do
  close(inp_unit)


!Start loop over different srcs, nphi, spec, [psi is not nested]
!nicrf_src: Number of sources

!! PICRF_NPHI_SRCS (nrho_icrf,num_nphi,nicrf_src,0:nspec:alla) : icrf power
!! PICRF_SRCS      (nrho_icrf,         nicrf_src,0:nspec:alla) : direct icrf power deposition summed over nphi
!! PICRF_TOTALS    (nrho_icrf,                   0:nspec:alla) : total radial power in each species
!! PICTH           (nrho_icrf,                               ) : direct thermal heating
!! PICRF_ABS       (                   nicrf_src             ) : power absorbed in the plasma per source
!! POWER_IC        (                   nicrf_src             ) : power on each ICRF source, :specified
!! wt_nphi_abs(num_nphi, nicrf_src) : nphi antenna weight. power_ic*wt_nphi_abs
!! NPHI(num_nphi(),nicrf_src) : nphi values for antenna for each source 
!!   may need to search for index given nphi from map file

  if (debug) then
     write(stderr,*) 'tot_num_nphi',tot_num_nphi," ,ps%num_nphi",ps%num_nphi
     write(stderr,*) 'nphi array shape',shape(ps%nphi),'values',ps%nphi
     write(stderr,*) 'map file nphi values',nphi_map
     write(stderr,*) 'map file nphi directories',nphi_dirname
     write(stderr,*) 'shape of picrf_nphi_srcs',shape(ps%picrf_nphi_srcs)
  end if
  allocate(alpha(ps%nicrf_src))
  alpha=0.0_rspec

  do isrc=1,ps%nicrf_src
     do iphi=1,ps%num_nphi(isrc)
        !construct index for toric run used in nphi_map_file
        if (isrc.eq.1) then
           fidx=iphi
        else
           fidx=(isrc-1)*ps%num_nphi(isrc-1)+iphi
        end if
        toric_output_file=trim(nphi_dirname(fidx))//'/toric.nc' !'/fort.21' !change to toric.nc when component is updated
        toric_data_file=trim(nphi_dirname(fidx))//'/toric_cfg.nc' !'/fort.9' !change to toric.nc when component is updated

  open(unit=inp_unit, file=trim(nphi_dirname(fidx))//'/torica.inp',status='old',form='formatted')
  INQUIRE(inp_unit, exist=lex)
  IF (lex) THEN
     read(inp_unit, nml = ips)
  ELSE
     write(*,*) &
          'torica.inp does not exist or there was a read error'
  END IF
  close(inp_unit)
!  toric_to_alla=(/1,2,3,4,5,6,7,9/) !uncomment for testing

  ierr = nf90_open(toric_output_file,imode,id_ncdf)
  CALL assert(ierr == NF90_NOERR, &
  &          'fatal, error in nf90_open in process_toric_output', ierr)

!!Get dimension sizes and allocate local variables
!internal toric var is nelpla=nelm-nelvac
  ierr = nf90_inq_dimid(id_ncdf,"PsiPwdDim",id_temp)
  ierr = nf90_inquire_dimension(id_ncdf, id_temp, len = nnoderho)

  if (.not. allocated(ant_nphi_load) ) then
     allocate(ant_nphi_load(nspec))
  end if

!toric number of ion species
  ierr = nf90_inq_dimid(id_ncdf,"SpecDim",id_temp)
  ierr = nf90_inquire_dimension(id_ncdf, id_temp, len = nspec)

!size of mhd rho meshes
  ierr = nf90_inq_dimid(id_ncdf,"NmhdDim",id_temp)
  ierr = nf90_inquire_dimension(id_ncdf, id_temp, len = nmhd)

  write(*,*) "process toric output: "//trim(toric_output_file)//": profile and spec sizes:"
  write(*,*) "toric mesh, PS mesh", nnoderho,ps%nrho_icrf
  write(*,*) "toric species , PS species count",nspec,ps%nspec_th
  write(*,*) "icrf sources",ps%nicrf_src
  write(*,*) 'alla index ',ps%alla_index
  write(*,*) 'alla name ',ps%alla_name
  write(*,*) 'alla type ',ps%alla_type

!wdot should be cell centered quantities of size nnoderho-1, so what are they in toric?
  if (.not. allocated(rhon))then
  allocate(rhon(nnoderho), dvol(nnoderho), prhon(nnoderho+1), &
   &        trhon(nnoderho+1), redotje_int(nnoderho), redotji_int(nspec,nnoderho), &
   &        wdote_int(2,nnoderho), wdoti_int(2,nspec,nnoderho), &
   &        rho(ps%nrho_icrf), wdot_tmp(nnoderho), darea(nnoderho), &
   &        CdensFW(nnoderho), CdensIBW(nnoderho), cdicrf_tmp(nnoderho), &
   &        cdicrf_int(nnoderho),    stat=ierr)
  CALL assert((ierr == 0),'allocation error in process_toric_output ', ierr)
  allocate(mtrhon(nmhd),mprhon(nmhd))
  allocate(tpwi(2,1:nspec))
  endif
  prhon=0. ; trhon=0.

!radial coordinates
 !this is r_{i+1/2} /a. i.e. cell centered normalized minor radius
  ierr = nf90_inq_varid(id_ncdf,"Pw_abscissa",id_temp)
  ierr = nf90_get_var(id_ncdf,id_temp,rhon)

 !this is sqrt(Psi_pol/Psi_pol_edge) on nnoderho mesh, should be uniform
 ![0, (nnoderho-1)/nnoderho]
  ncdfvarname="CurrDens_abscissa"
  ierr = nf90_inq_varid(id_ncdf,ncdfvarname,id_temp)
  call assert(ierr==0,'ncdf variable not found:'//trim(ncdfvarname),ierr)
  ierr = nf90_get_var(id_ncdf,id_temp,prhon(1:nnoderho))


 !this is sqrt(Psi_tor/Psi_tor_edge) mapping on nmhd size mesh
 !this is uniform [0,1]
  ierr = nf90_inq_varid(id_ncdf,"rho_tor",id_temp)
  ierr = nf90_get_var(id_ncdf,id_temp,mtrhon)

 !this is sqrt(Psi_pol/Psi_pol_edge) mapping on nmhd size mesh
 !nonuniform mapping relative to mtrhon [0,1]
  ierr = nf90_inq_varid(id_ncdf,"rho_pol",id_temp)
  ierr = nf90_get_var(id_ncdf,id_temp,mprhon)

!d volume in m^3/dpsi
  ierr = nf90_inq_varid(id_ncdf,"Spec_volume",id_temp)
  ierr = nf90_get_var(id_ncdf,id_temp,dvol)

!FW electron damping
  ierr = nf90_inq_varid(id_ncdf,"PwE",id_temp)
  ierr = nf90_get_var(id_ncdf,id_temp,wdote_int(1,:))

!IBW electron damping
  ierr = nf90_inq_varid(id_ncdf,"PwEIBW",id_temp)
  ierr = nf90_get_var(id_ncdf,id_temp,wdote_int(2,:))

!Ion damping Fundamental
  ierr = nf90_inq_varid(id_ncdf,"PwIF",id_temp)
  ierr = nf90_get_var(id_ncdf,id_temp,wdoti_int(1,:,:))

!Ion damping Harmonic
  ierr = nf90_inq_varid(id_ncdf,"PwIH",id_temp)
  ierr = nf90_get_var(id_ncdf,id_temp,wdoti_int(2,:,:))

!get total power
  ierr = nf90_inq_varid(id_ncdf,"TPwIF",id_temp)
  ierr = nf90_get_var(id_ncdf,id_temp,tpwi(1,:))
  ierr = nf90_inq_varid(id_ncdf,"TPwIH",id_temp)
  ierr = nf90_get_var(id_ncdf,id_temp,tpwi(2,:))
  ierr = nf90_inq_varid(id_ncdf,"TPwEFW",id_temp)
  ierr = nf90_get_var(id_ncdf,id_temp,tpwe(1))
  ierr = nf90_inq_varid(id_ncdf,"TPwEIBW",id_temp)
  ierr = nf90_get_var(id_ncdf,id_temp,tpwe(2))
  tpw = sum(tpwe)+sum(tpwi)

! PTB - begins 

! area in m^2 / dpsi
  ierr = nf90_inq_varid(id_ncdf,"Spec_area",id_temp)
  ierr = nf90_get_var(id_ncdf,id_temp,darea)

! Current density driven by FW's (A/m^2/MW)-> Note the toric.nc file is wrong ! This should be in (A/m^2/W)
  ierr = nf90_inq_varid(id_ncdf,"CurDensFW",id_temp)
  ierr = nf90_get_var(id_ncdf,id_temp,CdensFW) 

! Current density driven by IBW (A/m^2/MW) -> Note the toric.nc file is wrong ! This should be in (A/m^2/W)
  ierr = nf90_inq_varid(id_ncdf,"CurDensIBW",id_temp)
  ierr = nf90_get_var(id_ncdf,id_temp,CdensIBW) 

! Total current driven by FW's (A/MW) -> Note the toric.nc file is wrong ! This should be in (A/W)
  ierr = nf90_inq_varid(id_ncdf,"TotCurrFW",id_temp)
  ierr = nf90_get_var(id_ncdf,id_temp,CtotFW) 

! Total current driven by IBW (A/MW) -> Note the toric.nc file is wrong ! This should be in (A/W)
  ierr = nf90_inq_varid(id_ncdf,"TotCurrIBW",id_temp)
  ierr = nf90_get_var(id_ncdf,id_temp,CtotIBW) 

! PTB - ends

  ierr = nf90_close(id_ncdf)

!get antenna resistance from toric data file
  ierr = nf90_open(toric_data_file,imode,id_ncdf)
  CALL assert(ierr == NF90_NOERR, &
  &          'fatal, error in nf90_open toric data file in process_toric_output', ierr)
  ierr = nf90_inq_varid(id_ncdf,"total_power",id_temp)
  ierr = nf90_get_var(id_ncdf,id_temp,ant_load)
  ierr = nf90_close(id_ncdf)
!  ant_nphi_load(iphi)=ant_load

  print *,'source index ', isrc,nphi_isrc(fidx)
  print *,'nphi index, value: ', iphi, nphi_map(fidx), ps%nphi(iphi,isrc) !also checks both nphi are the same
  print *,'antenna load',ant_load
  print*,'powers in percent'
  print*,'total:',tpw
  print*,'electron total',tpwe
  print*,'ion totals fundamental',tpwi(1,:)
  print*,'ion totals harmonic',tpwi(2,:)
  print*,   '   number of thermal species (ps%nspec_tha) = ', ps%nspec_tha
  print*,   '   number of icrf sources (ps%nicrf_src) = ', ps%nicrf_src
! ps%nspec_alla = ps%nspec_th

  if(allocated(ps%freq_ic)) then
     print*, 'ps%freq_ic = ', ps%freq_ic
  else
     print*, 'ps%freq_ic not allocated'
  end if

  if(allocated(ps%icrf_src_name)) then
     print*, 'ps%icrf_src_name = ', ps%icrf_src_name
  else
     print*, 'ps%icrf_src_name not allocated'
  end if

  if(allocated(ps%power_ic)) then
     print*, 'ps%power_ic = ', ps%power_ic
  else
     print*, 'ps%power_ic not allocated'
  end if

! use mtrhon, but interpolated to nnoderho mesh on which poloidal flux is uniform
! find new toroidal flux mesh, trhon, by mapping mtrhon to corresponding poloidal mesh
! prhon which is the one originally corresponding to the power profiles.
     prhon(1)=0.0 ! ~epsilon from toric
     prhon(nnoderho+1) = 1.0 ! toric doesn't have last mesh point in prhon
     CALL interp_linear( mprhon, mtrhon, prhon, trhon )

     if (1==2) then
        write(0,*) 'trhon', trhon(:)
        write(0,*) 'mtrhon', mtrhon(:)
        write(0,*) 'prhon', prhon(:)
        write(0,*) 'mprhon', mprhon(:)
     endif

!  endif
  ! check for dimension consistancy between state and toric
  ! see plasma_state_definition_mod.f90 for details on state variables
!j  CALL assert( nnoderho == ps%nrho_icrf, 'toric-state power dimension inconsistancy', nnoderho)
!should implement interpolation here in case of failure.
  
  ! check species consistency--not all species will have RF
  CALL assert( ps%nspec_tha <= 30, 'state ion species greater than 30: ', ps%nspec_tha )

  ! load the grid

!Power densities should be cell centered
!also we want integrated total power vs radial point in plasma state
!Toric's are such that P(i) = avg power in between Psi(i) and Psi(i+1) / volume of that cell
!Units of W/m^3/W_inc

!for toric dvol=volume between surfaces/dpsi, so multiply by dpsi
!  do i=2,nnoderho
!     dvol(i) = dvol(i) * ( prhon(i)-prhon(i-1) )
!  end do
!  dvol(1)=dvol(2) !not used, but just to be careful
! 
!dvol = delta V/ delta prhon
!prhon is uniform

  dvol=dvol/real(nnoderho-1,rspec)

  if (debug) print*,'volume',sum(dvol),'m^3'

  ps%picrf_totals = 0.
  ps%picth = 0.

  j=nphi_isrc(fidx)
  if (debug) then
     write(*,*) 'antenna weights ',j,":",ps%nicrf_src
     write(*,*) ps%num_nphi(j)
     write(*,*) ps%wt_nphi_abs(:,j)
     write(*,*) ps%wt_nphi(:,j)
  end if

!!this section should now just apply antenna weighting
  tpw=0.0
  tpw=tpw+sum(wdote_int(1,:)*dvol)+sum(wdote_int(2,:)*dvol)
  do i=1,nspec   ! PTB nspec_th
     tpw=tpw+sum(wdoti_int(1,i,:)*dvol)+sum(wdoti_int(2,i,:)*dvol)
  end do
  if (debug) print*,"summed power",tpw,"MW/MW_inc" !should be 1 if everything is ok
!normalize, but want to use total power absorbed not just for given nphi

  wdote_int=wdote_int*ps%power_ic(j)*ps%wt_nphi(iphi,j)*ant_load!/tpw
  wdoti_int=wdoti_int*ps%power_ic(j)*ps%wt_nphi(iphi,j)*ant_load!/tpw
  if (debug) print *,'fidx',fidx,iphi,ps%wt_nphi(iphi,j),tpw
 print*,'ICRF Power', ps%power_ic(j)

!now we have to rezone to nrho_icrf sizes for trhon and wdot and dvol
!first need to resample rho to larger mesh.
!  CALL ps_user_rezone1(trhon,rho,wdot_tmp,ps%picrf_srcs(:,j, 0), &
!      ierr, nonorm = .True.)

!create new rho mesh which is uniform root toroidal of size ps%nrho_icrf
![0,1]
  do i=1,ps%nrho_icrf
     rho(i) = real(i-1,rspec)
  end do
  rho=rho/real(ps%nrho_icrf-1,rspec)
  ps%rho_icrf = rho

!electrons are in index 0
!picrf has size nicrf-1, it is the cell center value.

!interpolate power in mesh cell to larger mesh
   wdot_tmp= (wdote_int(1,:)+wdote_int(2,:))*dvol(:)

  if (debug) then
      print*,'wdot_tmp',sum(wdot_tmp), size(wdot_tmp)
  end if

!  CALL ps_user_rezone1(trhon,rho,wdot_tmp,ps%picrf_srcs(:,j, 0), &
!      ierr, nonorm = .True.)
!electrons
  CALL ps_user_rezone1(trhon,rho,wdot_tmp,ps%picrf_nphi_srcs(:,iphi,j, 0), &
      ierr, nonorm = .True.)
 
  if (debug) then 
    print*,'wdot_tmp e interp',sum(ps%picrf_nphi_srcs(:,iphi,j, 0)),size(ps%picrf_nphi_srcs(:,iphi,j, 0))
  endif
!  ps%picrf_totals(:,0) = ps%picrf_totals(:,0) + ps%picrf_srcs(:,j, 0)
 print*,'-----------------------------'
  if (debug) print *,"toric spec",0,"alla spec",0,'power',sum(wdot_tmp)

!store icrf power in ions for each source
  do isp=1,nspec
        i=toric_to_alla(isp)
        wdot_tmp= (wdoti_int(1,isp,:)+wdoti_int(2,isp,:))*dvol(:)
        if (debug) print *,"toric spec",isp,"alla spec",i,'power',sum(wdot_tmp)
!        if (debug) print *,"toric spec int",isp,"alla spec",i,'power',sum(ps%picrf_nphi_srcs(:,iphi, j, i)) !these are the same

!rezone and store in plasma state array for correct nphi, species and source
        CALL ps_user_rezone1(trhon,rho,wdot_tmp,ps%picrf_nphi_srcs(:,iphi, j, i), &
            ierr, nonorm = .True.)
     end do
  if (debug) then 
    print*,'wdot each nphi/source interp',sum(ps%picrf_nphi_srcs(:,iphi, j, :))
  endif
  ! Sum up antenna resistance for each source
  alpha(j) = alpha(j) + ps%wt_nphi(iphi,j)*ant_load

! PTB begins
! Create the array of driven current (in A) per zone - cdicrf_tmp
! First get the specific areas correct. 
! Need to multiply by dPsi = 1/(nnoderoho-1) since the toric.nc file stores d(Area) / dPsi.
! Also need to divide by 2*pi since darea comes from dvol which has a 2*pi in it. 

  darea = darea / (2.0_rspec * 3.14159_rspec * real(nnoderho-1,rspec))
  if (debug) print*,'area',sum(darea),'m^2'

! Form the total ICRF current density in A/m^2 from the given injected ICRF power
  cdicrf_tmp(:) = (CdensFW(:) + CdensIBW(:)) * ps%power_ic(j)

  cdicrf_int(:) = cdicrf_tmp(:) * darea(:)
    print*, 'For source ',j,' nphi ', nphi_map(fidx)
    print*, 'Total FW Driven Current (A)', SUM(cdicrf_int(:))
    print*, 'Total ICRF Power (W)', ps%power_ic(j)
    print*, 'FWCD Efficiency (A/W)',SUM(cdicrf_int(:))/ps%power_ic(j)

! Remap from the TORIC radial mesh (trhon) to the radial mesh of the Plasma State (rho)
    CALL ps_user_rezone1(trhon,rho,cdicrf_int,ps%cdicrf_nphi(:,iphi,j), &
         ierr, nonorm = .True.)
    ps%cdicrf(:,j)=ps%cdicrf(:,j)+ps%cdicrf_nphi(:,iphi,j)
!        CALL ps_user_rezone1(trhon,rho,cdicrf_int,ps%curich(:), &
!            ierr, nonorm = .True.)

! PTB end

        
end do
    print*, 'Total FW Driven Current from each source-> Interpolated (A)',j,SUM(ps%cdicrf(:,j))
    ps%curich(:)=ps%curich(:)+ps%cdicrf(:,j)
end do
    print*, 'Total FW Driven Current from all sources -> Interpolated (A)',SUM(ps%curich(:))

!Normalize with respect to loading
do j=1,ps%nicrf_src
ps%picrf_nphi_srcs(:,:, j, :)=ps%picrf_nphi_srcs(:,:, j, :)/alpha(j)
print *,'alpha',j,alpha(j)
enddo
print *,'total power check from ps',sum(ps%picrf_nphi_srcs)

!Cleanup


  print*,   '  --storing toric data in current plasma state'
  ! stores the modified plasma state--doesn't commit it to the previous one is
  ! still around


  do isrc= 1, ps%nicrf_src
     do iphi = 1,ps%num_nphi(isrc)
        ps%wt_nphi_abs(iphi,isrc)=sum(ps%picrf_nphi_srcs(:,iphi,isrc,:))/ps%power_ic(isrc)
        ps%picrf_srcs(:,isrc,:)=ps%picrf_srcs(:,isrc,:)+ps%picrf_nphi_srcs(:,iphi,isrc,:)
     end do
     ps%picrf_totals=ps%picrf_totals+ps%picrf_srcs(:,isrc,:)
     ps%picrf_abs(isrc)=sum(ps%picrf_srcs(:,isrc,:))
  end do
  do isp = 1,ps%nspec_alla !does not include electrons
     ps%picth = ps%picth + ps%picrf_totals(:,isp)
  end do

  if (debug) then
     print *,'power check with power_ic and picrf_abs: ',ps%power_ic, ps%picrf_abs,ps%wt_nphi_abs
  end if
  


    !--------------------------------------------------------------------------    !
    ! Store the data in partial plasma_state file
    !--------------------------------------------------------------------------

	CALL PS_WRITE_UPDATE_FILE('RF_IC_'//cur_state_file, ierr)
	WRITE (*,*) "Stored Partial RF Plasma State"   

!write the state file to optional filename, can also take optional state
  CALL ps_store_plasma_state(ierr , trim(cur_state_file))
  CALL assert( ierr == 0, 'cannot open state in prepare toric output', ierr )

!JCW diagnostic output
  if (debug) then
    print*, ' ps%rho_icrf = ', ps%rho_icrf
    print*, ' '
  end if

    print*, 'electron power, integral'
    print*, SUM(ps%picrf_totals(:,0)),"Watts"
    print*, 'FW ion power, integral'
    print*, SUM(ps%picth(:)),"Watts"

!Cleanup
  deallocate(rhon, dvol, redotje_int, redotji_int, &
       &    wdote_int, wdoti_int, darea, cdicrf_int, &
       &    cdicrf_tmp, CdensFW, CdensIBW )

  contains
    

  subroutine r8vec_bracket ( x, xval, left, right )

    !*****************************************************************************80
    !
    !! R8VEC_BRACKET searches a sorted R8VEC for successive brackets of a value.
    !
    !  Discussion:
    !
    !    An R8VEC is an array of double precision real values.
    !
    !    If the values in the vector are thought of as defining intervals
    !    on the real line, then this routine searches for the interval
    !    nearest to or containing the given value.
    !
    !  Modified:
    !
    !    06 April 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, length of input array.
    !
    !    Input, real ( kind = 8 ) X(N), an array sorted into ascending order.
    !
    !    Input, real ( kind = 8 ) XVAL, a value to be bracketed.
    !
    !    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
    !    Either:
    !      XVAL < X(1), when LEFT = 1, RIGHT = 2;
    !      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
    !    or
    !      X(LEFT) <= XVAL <= X(RIGHT).
    !

    integer,intent(out):: left, right
    real(rspec),intent(in)::     x(:), xval

    integer::  n, i

    n=size(x)
    do i = 2, n - 1

       if ( xval < x(i) ) then
          left = i - 1
          right = i
          return
       end if

    end do

    left = n - 1
    right = n

    return
  end subroutine r8vec_bracket

  subroutine interp_linear ( t_data, p_data, t_interp, p_interp )

    !*****************************************************************************80
    !
    !! INTERP_LINEAR applies piecewise linear interpolation to data.
    !
    !  Discussion:
    !
    !    From a space of DIM_NUM dimensions, we are given a sequence of
    !    DATA_NUM points, which are presumed to be successive samples
    !    from a curve of points P.
    !
    !    We are also given a parameterization of this data, that is,
    !    an associated sequence of DATA_NUM values of a variable T.
    !    The values of T are assumed to be strictly increasing.
    !
    !    Thus, we have a sequence of values P(T), where T is a scalar,
    !    and each value of P is of dimension DIM_NUM.
    !
    !    We are then given INTERP_NUM values of T, for which values P
    !    are to be produced, by linear interpolation of the data we are given.
    !
    !    Note that the user may request extrapolation.  This occurs whenever
    !    a T_INTERP value is less than the minimum T_DATA or greater than the
    !    maximum T_DATA.  In that case, linear extrapolation is used.  
    !
    !  Modified:
    !
    !    03 December 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) T_DATA(DATA_NUM), the value of the
    !    independent variable at the sample points.  The values of T_DATA
    !    must be strictly increasing.
    !
    !    Input, real ( kind = 8 ) P_DATA(DATA_NUM), the value of the
    !    dependent variables at the sample points.
    !
    !    Input, real ( kind = 8 ) T_INTERP(INTERP_NUM), the value of the
    !    independent variable at the interpolation points.
    !
    !    Output, real ( kind = 8 ) P_INTERP(DATA_NUM), the interpolated
    !    values of the dependent variables at the interpolation points.
    !
    implicit none

    integer :: data_num
    integer :: interp_num

    integer :: interp
    integer :: left
    real(rspec)    :: p_data(:) !(dim_num,data_num)
    real(rspec)    :: p_interp(:) !(dim_num,interp_num)
!    logical :: r8vec_ascends_strictly
    integer :: right
    real(rspec)    :: t
    real(rspec)    :: t_data(:) !(data_num)
    real(rspec)    :: t_interp(:) !(interp_num)


    data_num = size(p_data,1)
    interp_num = size(t_interp,1)
!    if ( .not. r8vec_ascends_strictly ( data_num, t_data ) ) then
!       write ( *, '(a)' ) ' '
!       write ( *, '(a)' ) 'INTERP_LINEAR - Fatal error!'
!       write ( *, '(a)' ) '  Independent variable array T_DATA is not strictly increasing.'
!       stop
!    end if

    do interp = 1, interp_num

       t = t_interp(interp)
       !
       !  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
       !  nearest to, TVAL.
       !
       call r8vec_bracket ( t_data, t, left, right )

       p_interp(interp) = &
            ( ( t_data(right) - t                ) * p_data(left)   &
            + (                 t - t_data(left) ) * p_data(right) ) &
            / ( t_data(right)     - t_data(left) )

    end do

    return
  end subroutine interp_linear


end program process_toric_output_mcmd_mm
