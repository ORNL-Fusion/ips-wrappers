program process_toric_output

  !---------------------------------------------------------------------
  ! processing toric output
  !last update on 04/02/08 by JCW to move PS setup to do_toric_init and 
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
    & CdensFW(:), CdensIBW(:), cdicrf_tmp(:), cdicrf_int(:)
  real(rspec):: frequency, tpwe(2),tpw, CtotFW, CtotIBW
  logical, parameter:: debug=.false.
    
  integer, parameter :: swim_string_length = 256  !for compatibility LAB
     


  !------------------------------------
  !  local
  INTEGER :: ierr, n, id_ncdf, imode, nspec, id_temp, i, j
  integer :: istat, nmhd
  integer :: nnoderho, iarg

  character(len =swim_string_length) :: cur_state_file, ncdfvarname

  write(*,*) ' -- get (restore) plasma state from file -- '
  ! this call retrieves the plasma state at the present ps%, and
  ! previous time steps psp%
  
  CALL get_arg_count(iarg)
  SELECT CASE (iarg)

      case(0)
         cur_state_file="cur_state.cdf"

      case(1)
         CALL get_arg(1,cur_state_file)

      case(2:)
         write(0,*) 'Error. Illegal number of arguments.'
         write(0,*) 'process_toric_output: '
	 write(0,*) 'process_toric_output cur_state_file'
	 stop 'incorrect command line arguments'

  end select

  print*, 'process_toric_output: cur_state_file = ', trim(cur_state_file)
  CALL ps_get_plasma_state(ierr, trim(cur_state_file))

  CALL assert( ierr==0,' process toric: ps_get_plasma_state: ierr=',ierr )

  print *,"freq_ic, picrf  alloc?",allocated(ps%freq_ic),allocated(ps%picrf_srcs)
  if (allocated(ps%picrf_srcs) .eqv. .TRUE.) then
      print *,"RF alloc?",allocated(ps%picrf_srcs),size(ps%picrf_srcs)
      print*,   '   number of icrf sources = ', ps%nicrf_src
  endif

  print *,  '   number of species = ', ps%nspec_alla, ps%nspec_th

 if (debug) then 
    if (allocated(ps%rho_icrf) .eqv. .TRUE.) then
       print*,   '   radial grid points for icrf = ', ps%nrho_icrf
       print*,   '   ps%rho_icrf = ', allocated(ps%rho_icrf), size(ps%rho_icrf), ps%rho_icrf
    endif
 endif
 

!##  print *,"RF alloc?",allocated(ps%picrf_srcs),size(ps%picrf_srcs)
!##  print*,   '   number of icrf sources = ', ps%nicrf_src
!##  print *,  '   number of species = ', ps%nspec_alla, ps%nspec_th
!##  print*,   '   radial grid points for icrf = ', ps%nrho_icrf
!##  print*,   '   ps%rho_icrf = ', allocated(ps%rho_icrf), size(ps%rho_icrf), ps%rho_icrf
 

!get some run parameters
  imode=NF90_NOWRITE !read only mode(default), dataset id returned in id_ncdf
  !netcdf does its own file unit testing
  ierr = nf90_open('toric_cfg.nc',imode,id_ncdf)
  CALL assert(ierr == NF90_NOERR, &
  &          'panic, error in nf90_open in process_toric_output', ierr)
  ierr = nf90_inq_varid(id_ncdf,"frequency",id_temp)
  ierr = nf90_get_var(id_ncdf,id_temp,frequency)
  ierr = nf90_close(id_ncdf)

  imode=NF90_NOWRITE !read only mode(default), dataset id returned in id_ncdf
  !netcdf does its own file unit testing
  ierr = nf90_open('toric.nc',imode,id_ncdf)
  CALL assert(ierr == NF90_NOERR, &
  &          'panic, error in nf90_open in process_toric_output', ierr)

!!Get dimension sizes and allocate local variables
!internal toric var is nelpla=nelm-nelvac
  ierr = nf90_inq_dimid(id_ncdf,"PsiPwdDim",id_temp)
  ierr = nf90_inquire_dimension(id_ncdf, id_temp, len = nnoderho)


!toric number of ion species
  ierr = nf90_inq_dimid(id_ncdf,"SpecDim",id_temp)
  ierr = nf90_inquire_dimension(id_ncdf, id_temp, len = nspec)

!size of mhd rho meshes
  ierr = nf90_inq_dimid(id_ncdf,"NmhdDim",id_temp)
  ierr = nf90_inquire_dimension(id_ncdf, id_temp, len = nmhd)

  write(*,*) "process toric output: profile and spec sizes:"
  write(*,*) "toric mesh, PS mesh", nnoderho,ps%nrho_icrf
  write(*,*) "toric species, PS species counts",nspec,ps%nspec_th
!wdot should be cell centered quantities of size nnoderho-1, so what are they in toric?
  allocate(rhon(nnoderho), dvol(nnoderho), prhon(nnoderho+1), &
   &        trhon(nnoderho+1), redotje_int(nnoderho), redotji_int(nspec,nnoderho), &
   &        wdote_int(2,nnoderho), wdoti_int(2,nspec,nnoderho), &
   &        rho(ps%nrho_icrf), wdot_tmp(nnoderho), darea(nnoderho), &
   &        CdensFW(nnoderho), CdensIBW(nnoderho), cdicrf_tmp(nnoderho), &
   &        cdicrf_int(nnoderho),    stat=ierr)
  allocate(mtrhon(nmhd),mprhon(nmhd))
  allocate(tpwi(2,1:nspec))
  CALL assert((ierr == 0),'allocation error in process_toric_output ', ierr)

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
  if (debug) write(0,*) "prhon:",ierr,prhon

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

  print*,'powers in percent'
  print*,'total:',tpw
  print*,'electron total',tpwe
  print*,'ion totals fundamental',tpwi(1,:)
  print*,'ion totals harmonic',tpwi(2,:)
  print*,   '   number of thermal species (ps%nspec_th) = ', ps%nspec_th
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

!j if(allocated(ps%rho_icrf)) then
!j    print*, 'ps%nrho_icrf = ', ps%nrho_icrf
!j else
!j    ps%nrho_icrf = nnoderho !+1 !JCW to allocated
!j    print *,'Allocating RF in prepare_input'
!j    ps%freq_ic(1)=frequency/1.e6

!j    CALL ps_alloc_plasma_state(ierr) 
!j    CALL assert( ierr == 0, trim(cur_state_file)//' ps_alloc_plasma_state: ierr=',ierr )

! use mtrhon, but interpolated to nnoderho mesh on which poloidal flux is uniform
! find new toroidal flux mesh, trhon, by mapping mtrhon to corresponding poloidal mesh
! prhon which is the one originally corresponding to the power profiles.
     prhon(1)=0.0 ! ~epsilon from toric
     prhon(nnoderho+1) = 1.0 ! toric doesn't have last mesh point in prhon
     CALL interp_linear( mprhon, mtrhon, prhon, trhon )

     if (debug) then
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
  
  ! check to species consistency--not all species will have RF
  CALL assert( ps%nspec_th <= 30, 'state ion species greater than 30: ', ps%nspec_th )

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
!there is just one source right now.
  j=1


  tpw=0.0
  tpw=tpw+sum(wdote_int(1,:)*dvol)+sum(wdote_int(2,:)*dvol)
  do i=1,nspec   ! PTB nspec_th
     tpw=tpw+sum(wdoti_int(1,i,:)*dvol)+sum(wdoti_int(2,i,:)*dvol)
  end do
  if (debug) print*,"summed power",tpw,"MW/MW_inc" !should be 1 if everything is ok
  wdote_int=wdote_int*ps%power_ic(j)/tpw
  wdoti_int=wdoti_int*ps%power_ic(j)/tpw

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

!interpolate to larger mesh
   wdot_tmp= (wdote_int(1,:)+wdote_int(2,:))*dvol(:)

  if (debug) then
      print*,'wdot_tmp e',sum(wdot_tmp), size(wdot_tmp)
      print*,wdot_tmp
  end if

  CALL ps_user_rezone1(trhon,rho,wdot_tmp,ps%picrf_srcs(:,j, 0), &
      ierr, nonorm = .True.)
 
  if (debug) then 
    print*,'wdot_tmp e interp',sum(ps%picrf_srcs(:,j, 0)),size(ps%picrf_srcs(:,j, 0))
    print*,ps%picrf_srcs(:,j, 0)
  endif
   ps%picrf_totals(:,0) = ps%picrf_totals(:,0) + ps%picrf_srcs(:,j, 0)


  if(ps%nspec_th >= 1) then
     do i=1,nspec   ! PTB nspec_th 
!store icrf power in ions for each source
        wdot_tmp= (wdoti_int(1,i,:)+wdoti_int(2,i,:))*dvol(:)
!        CALL interp_linear(trhon(:),wdot_tmp,rho(:),ps%picrf_srcs(:,j, i))

        CALL ps_user_rezone1(trhon,rho,wdot_tmp,ps%picrf_srcs(:,j, i), &
            ierr, nonorm = .True.)
 
!sum over sources
        ps%picrf_totals(:,i) = ps%picrf_totals(:,i) +  ps%picrf_srcs(:, j, i)
!sum over species
        ps%picth(:) = ps%picth(:) + ps%picrf_totals(:,i)
     end do
  end if

  print*,   '  --storing toric data in current plasma state'
  ! stores the modified plasma state--doesn't commit it to the previous one is
  ! still around

! PTB begins
! Create the array of driven current (in A) per zone - cdicrf_tmp
! First get the specific areas correct. 
! Need to multiply by dPsi = 1/(nnoderoho-1) since the toric.nc file stores d(Area) / dPsi.
! Also need to divide by 2*pi since darea comes from dvol which has a 2*pi in it. 

  darea = darea / (2.0_rspec * 3.14159_rspec * real(nnoderho-1,rspec))
  if (debug) print*,'area',sum(darea),'m^2'

! Form the total ICRF current density in A/m^2 from the given injected ICRF power
  cdicrf_tmp(:) = (CdensFW(:) + CdensIBW(:)) * ps%power_ic(j)

! Remember we are assuming just one ICRF source so j=1
  cdicrf_int(:) = cdicrf_tmp(:) * darea(:)
    print*, 'Total FW Driven Current (A)', SUM(cdicrf_int(:))
    print*, 'Total ICRF Power (W)', ps%power_ic(j)
    print*, 'FWCD Effciency (A/W)',SUM(cdicrf_int(:))/ps%power_ic(j)
! Remap from the TORIC radial mesh (trhon) to the radial mesh of the Plasma State (rho)
        CALL ps_user_rezone1(trhon,rho,cdicrf_int,ps%cdicrf(:,j), &
            ierr, nonorm = .True.)
        CALL ps_user_rezone1(trhon,rho,cdicrf_int,ps%curich(:), &
            ierr, nonorm = .True.)
    print*, 'Total FW Driven Current from each source-> Interpolated (A)',SUM(ps%cdicrf(:,j))
    print*, 'Total FW Driven Current from all sources -> Interpolated (A)',SUM(ps%curich(:))

! PTB end


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
    
  subroutine assert( lcond, mesg, ivalue )
    logical, intent(in) :: lcond
    character*(*),intent(in) :: mesg
    integer, intent(in) :: ivalue

    if (.not.lcond) then
       write(0,*) mesg,ivalue
       stop '** assertion error **'
    endif
    return
  end subroutine assert


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


end program process_toric_output
