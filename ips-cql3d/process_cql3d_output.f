      program process_cql3d_output

  !- -------------------------------------------------------------------
  ! " processing cql3d output
  !Lee A. Berry 12/6/06  based on D. McMcune's test program
  !BobH, 03/28/07        based on Berry  process_aorsa_output
  !Long-Poe, 02/13/08    Updated from PS1 to PS2
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
  !   CHARACTER*256 :: state_file = 'cur_state.cdf'
  !   CHARACTER*256 :: prior_file = 'prev_state.cdf'
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
  !    SUBROUTINE ps_update_state(ierr)
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
  !      of profiels f(rho,theta) -- not yet implemented DMC 16 Oct 2006
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

      use swim_global_data_mod, only :
     1 rspec, ispec,                 ! int: kind specification for real and integer
!!            & swim_string_length, & ! length of strings for names, files, etc.
     1 swim_error		    ! error routine
    

    
!--------------------------------------------------------------------------
!
!   Data declarations
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------------
!
!   AORSA data that will be given to the state via the swim_out file
!   the aorsa output data is in a file called swim_out
!
!--------------------------------------------------------------------------
   

!BH070328:  Begin by reading in toroidal current density from cql3d
!           output netcdf file (mnemonic.nc) and putting it into the state

 
!      IMPLICIT NONE
      implicit integer (i-n), real*8 (a-h,o-z)

      include 'netcdf.inc'
    
      integer, parameter :: swim_string_length = 256 !for compatibility LAB
      integer :: nnoderho

!     Storage for netcdf file elements and retrieval
      real*8, allocatable, dimension(:)   :: rya
      real*8, allocatable, dimension(:,:) :: curtor  !just printout for now
      real*8, allocatable, dimension(:,:) :: ccurtor !just printout for now
      real*8, allocatable, dimension(:,:) :: denra
      real*8, allocatable, dimension(:,:) :: curra

      integer :: ncid,vid,istatus
      integer :: start(2),count(2)
      integer :: nt_id,lrz_id   !Time step dim id, radial fp bins id
      integer :: ngen_id        !number of general species id
      integer :: nt,lrz,ngen
      character*128 :: name

 
!$$$!     Cubic spline related storage
!$$$      real, allocatable, dimension(:) :: w,work
!$$$      integer, dimension(:) :: iop(2)
!$$$      real, dimension(:) :: tab(3)
!$$$      integer, dimension(:) :: itab(3)

!     Variables which could be put into the plasma state
      real*8, allocatable, dimension(:) :: PScurtor,PSrho_fp
      real*8 :: PStot_curtor

     


  !------------------------------------
  !  local
      INTEGER :: ierr
      INTEGER :: iout = 6

   
      write(iout,*) 'process_cql3d_output: -- restoring plasma state from file -- '

  ! this call retrives the plasma state at the present ps%, and
  ! prvious time steps psp%, plus sets storage for the variables
  
      CALL ps_get_plasma_state(ierr)

      if(ierr.ne.0) then
         write(iout,*) ' process cql3d: ps_get_plasma_state: ierr=',ierr
         stop
      endif

!.......................................................................
!     Open cql3d netcdf file
!.......................................................................
      ncid = ncopn('mnemonic.nc',NCNOWRIT,istatus)
      write(*,*)'after ncopn ncid=',ncid,'istatus',istatus

!     read in dimension IDs
      lrz_id = ncdid(ncid,'rdim',istatus)
!$$$      write(*,*)'proc_cql3d_op: after ncdid lrz_id',lrz_id,'istatus',istatus
      nt_id = ncdid(ncid,'tdim',istatus)
!$$$      write(*,*)'proc_cql3d_op: after ncdid nt_id',nt_id,'istatus',istatus
      ngen_id = ncdid(ncid,'gen_species_dim',istatus)
!$$$      write(*,*)'proc_cql3d_op:after ncdid ngen_id',ngen_id,'istatus',istatus


!.......................................................................
!     inquire about dimension sizes:time steps,grid size,#species---
!.......................................................................
      call ncdinq(ncid,lrz_id,name,lrz,istatus)
      write(*,*)'proc_cql3d_op: after ncdinq, # of rad bins=',lrz,
     .     '  istatus=',istatus
      call ncdinq(ncid,nt_id,name,nt,istatus)
      write(*,*)'proc_cql3d_op: after ncdinq, # of t steps =',nt,
     .     '  istatus=',istatus
      call ncdinq(ncid,ngen_id,name,ngen,istatus)
      write(*,*)'proc_cql3d_op: after ncdinq, #  gen specs =',ngen,
     .     '  istatus=',istatus

!     allocate space for arrays to be read
      allocate (rya(lrz))
      allocate (curtor(lrz,nt))
!$$$      allocate (w(lrz))
!$$$      allocate (work(3*lrz+1))
      allocate (ccurtor(lrz,nt))
      allocate (denra(lrz,nt))
      allocate (curra(lrz,nt))

!     Variables for possible use in Plasma State
!$$$      allocate (PScurtor(ps%nrho-1))

!     Specify reading ranges
      start(1:2)=1
      count(1)=lrz
      count(2)=nt

!.......................................................................
!     read netcdf variables
!.......................................................................
      vid=ncvid(ncid,'rya',istatus)
!$$$      write(*,*)'proc_cql3d_op: after ncvid, rya =',vid,
!$$$    .     '  istatus=',istatus
      call ncvgt(ncid,vid,1,lrz,rya,istatus) !normalized rho, bin centers
!$$$      write(*,*)'proc_cql3d_op: after ncvgt rya, istatus =',istatus
      write(*,*)'proc_cql3d_op: after ncvgt, rya =',rya

      vid=ncvid(ncid,'curtor',istatus)
!$$$      write(*,*)'proc_cql3d_op: after ncvid, curtor =',vid,'istatus=',
!$$$     .     istatus
      call ncvgt(ncid,vid,start,count,curtor,istatus) !tor cur/bin centers

      vid=ncvid(ncid,'ccurtor',istatus)
      call ncvgt(ncid,vid,start,count,ccurtor,istatus) !cum tor cur/bin centers
      write(*,*)'proc_cql3d_op: after ncvgt, tot curtor =',ccurtor(lrz,nt)

      vid=ncvid(ncid,'denra',istatus)
      call ncvgt(ncid,vid,start,count,denra,istatus) !runaway elec density
      write(*,*)'proc_cql3d_op: after ncvgt, denra =',denra(1:lrz,nt)

      vid=ncvid(ncid,'curra',istatus)
      call ncvgt(ncid,vid,start,count,curra,istatus) !runaway current density
      write(*,*)'proc_cql3d_op: after ncvgt, curra =',curra(1:lrz,nt)

      call ncclos(ncid,istatus)

!     CQL3D rho(1:lrz) will generally be different from plasma_state.
!     For current and density, will simply rebin by cubic spline
!     to plasma_state variables, then renormalize to preserve
!     totals from cql3d.  See zcunix.f (from cql3d distn) for spline routines.
!     ps%rho(1:ps%nrho) is on bin boundaries.  Need bin centers.

      PStot_curtor=ccurtor(lrz,nt)

!.......................................................................
!     Set some Plasma State variables
!.......................................................................

!  Necessary PS dimensions (fp user has to set his fp component dims):
!     ps%nspec_nonMax=ngen                        !Usually 1
      ps%nrho_rw=lrz+1                            !Length of bin bndry array
      call ps_alloc_plasma_state(ierr)            !set these PS dims in PS
      if (ierr.ne.0) then
         write(*,*)' proc_cql3d_op:  ps_alloc_plasma_state: ierr=',ierr
         stop
      endif

!     no longer exists--contained in component description
!     ps%fp_codes(1)='cql3d-by-itself'

      allocate (PSrho_fp(lrz+1))
  
      PSrho_fp(1)=0.
      PSrho_fp(lrz+1)=1.
      do ll=2,lrz
         PSrho_fp(ll)=0.5*(rya(ll-1)+rya(ll))
      enddo
      write(*,*)'proc_cql3d_op: PSrho_fp =',PSrho_fp

      write(*,*)'proc_cql3d_op: ps%nrho_fp =',ps%nrho_rw
      write(*,*)'size of ps%rho_fp = ',size(ps%rho_rw)
      ps%rho_rw=PSrho_fp(1:ps%nrho_rw)           !Bin bndries
	  if (allocated(ps%rho_rw)) then 
	     print*, 'ps%rho_rw IS ALLOCTED'
	  else
	     print*, 'ps%rho_rw IS NOTNOTNOT ALLOCTED'
		 stop
	  end if
!      ps%nonMax(1)=7   !  fast electrons = type 7   DOESN'T WORK

! fast electrons...
! I don't think any of this is necessary now LAB
!     ps%nonMax_type(1)=ps_fast_electron ! code for fast electrons
!     CALL ps_label_spectype(-1, -1, -1, ps%nonMax_type(1),
!    .  ps%nonMax_name(1), ps%qatom_nonMax(1), ps%q_nonMax(1),
!    .  ps%m_nonMax(1))

      ps%nrw(1:lrz)=denra(1:lrz,nt)*1.e6 !Runaway density above ucrit
      ps%cur_rw(1:lrz)=curra(1:lrz,nt)*1.e4    !Runaway curr dens above ucrit
      ps%dcur_rw_dvloop(1:lrz)=0.0             !For now, BH070312


!.......................................................................
!      Fill in "all" species list, and Store data in PS
!.......................................................................

!     CALL ps_merge_species_lists(ierr)
!     if(ierr .ne. 0) print *, "ps merge error"
!     print *, "total species : ", ps%nspec_all

      do ii=0,ps%nspec_th
       write(iout,1001) ii,ps%s_type(ii),
     .          trim(ps%s_name(ii)),ps%q_s(ii),ps%m_s(ii)
      enddo	
      if(ps%nspec_fusion .ne. ps_uninit .or. allocated(ps%m_sfus)) then
      do ii=1,ps%nspec_fusion
       write(iout,1001) ii,ps%sfus_type(ii),
     .          trim(ps%sfus_name(ii)),ps%q_sfus(ii),ps%m_sfus(ii)
      enddo	
      endif
      if(ps%nspec_rfmin .ne. ps_uninit .or. allocated(ps%m_rfmin)) then
      do ii=1,ps%nspec_rfmin
       write(iout,1001) ii,ps%rfmin_type(ii),
     .          trim(ps%rfmin_name(ii)),ps%q_rfmin(ii),ps%m_rfmin(ii)
      enddo	
      endif
      if(ps%nspec_beam .ne. ps_uninit .or. allocated(ps%m_snbi)) then
      do ii=1,ps%nspec_beam
       write(iout,1001) ii,ps%snbi_type(ii),
     .          trim(ps%snbi_name(ii)),ps%q_snbi(ii),ps%m_snbi(ii)
      enddo
      endif

 1001 format(' Specie index,type: ',i2,1x,i2,1x,
     .   	   '"',a,'" charge,mass: ',2(1pe12.5,1x))



      print*,   'process_cql3d_output: --storing cql3d data in current PS'
!     stores the modified plasma state--doesn't commit it. The previous one is
!     still around

      call ps_store_plasma_state(ierr)
      if(ierr .ne. 0) then
      write(iout,*) '  cannot open state in process_cql3d_output'
      end if

      
      end program process_cql3d_output
