PROGRAM do_toric_init
  !JCW created APR 2008
  !This is the initialization program called by rf_ic_toric.py init step
  !It sets up the rf component parts of the plasma state.  These have been
  !moved from process_toric_output

!--------------------------------------------------------------------------
! Multi toroidal mode version 0.0 11/15/2012 (Batchelor) 
!
! Modified from do_toric_init_abr.f90  It reads a plasma state file
! ps which has EPA data already initialized.  Then it reads namelist files containing
! ICRF machine description and shot configuration namelists using the plasma
! state functions ps_mdescr_read() and ps_sconfig_read() and stores this data in 
! temporary plasma state file , aux.  It copies the new ICRF data to the current plasma
! state and saves it.  The original section of do_toric_init does further initialization 
! to allocate and fill the profile arrays as is normally done in the component init.  That
! is the data dimensioned by nrho_icrf.
!
!       Don Batchelor
!       ORNL
!       Oak Ridge, TN 37831
!       batchelordb@ornl.gov
!
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! Details:
!
!	1) Machine description namelist file: mdescr_namelist.dat.  mdescr files can be
!	written by plasma state function ps_mdescr_write().  Sample mdescr files are produced
!	by the plasma_state_test routine.  The ICRF relevant items in namelist mdescr are:
!        icrf_src_name  = number & name of ICRF sources
!        ant_model = antenna model filenames (1 per antenna source)
!        nrz_antgeo = number of (R,Z) points, antenna geometries
!        max_nrz_antgeo = Maximum size of variable length enumeration: NRZ_ANTGEO
!        R_antgeo = antenna geo: R pts 
!        Z_antgeo = antenna geo: Z pts 
!        dx_fshield = distance, antenna to Faraday shield
!
!	2) Shot configuration namelist file: sconfig_namelist.dat.  sconfig files can be
!	written by plasma state function ps_sconfig_write(). Sample sconfig files are produced
!	by the plasma_state_test routine.  The ICRF relevant items in namelist sconfig are:
!        RFMIN = ICRF minority species
!        freq_ic = frequency on each ICRF source
!        n_straps = number of straps in the antenna
!        max_n_straps = Maximum size of variable length enumeration: N_STRAPS
!        R_ant = major radius of antenna
!        Z_ant = height of antenna
!        Z_mid_ant = center of antenna relative to TF coil midplane
!        num_nphi_vac = number of non-zero n_phi values
!        max_num_nphi_vac = Maximum size of variable length enumeration: NUM_NPHI_VAC
!        Real_ant_coef = real part of Fourier Coef
!        Imag_ant_coef = imag part of Fourier Coef
!        num_nphi = number of non-zero n_phi values
!        max_num_nphi = Maximum size of variable length enumeration: NUM_NPHI
!        nphi = n_phi wave spectrum from antenna
!        wt_nphi = vacuum spectrum n_phi weight
!
!	3)	The code takes one command-line argument, the name of the current plasma state
!		to be initialized.  The names of the namelist files are assumed to be
!		mdescr_namelist.dat and sconfig_namelist.dat.
!
! One other point: The EPA component is responsible for setting the number of ICRF
! sources, nicrf_src and allocating the arrays it dimensions.  This is because EPA
! controls the power on the different sources.  Note that nicrf_src also 
! appears in the mdescr namelist, although it appears there implicitly as the number
! of elements in the namelist that are dimensioned nicrf_src.  For example icrf_src_name,
! ant_model and dx_fshield. The user must make sure this value is consistent
! with what is set by EPA.  I know, data should only be set in one place, but it
! would be complicated to fix, since the EPA may not be initialized from the same,
! or indeed any, mdescr namelist.  Just be careful at run time, ok.
!
! For that matter many of the arrays in the sconfig namelist are also dimensioned
! nicrf_src, for example freq_ic.  So the dimensons in the two namelist files must be
! consistent.
!
!--------------------------------------------------------------------------

  USE plasma_state_mod

  USE swim_global_data_mod, only : &
       & rspec, ispec, &                ! int: kind specification for real and integer
       & swim_error                     ! error routine

  USE toric_utils_mod


  IMPLICIT NONE

  integer, parameter :: swim_string_length = 256  !for compatibility LAB
  character(len =swim_string_length) :: cur_state_file, program_name
    
  INTEGER :: cclist(ps_ccount)   ! component activation list 
  INTEGER:: inp_unit, ierr, iarg, i, iwarn = 0
  LOGICAL:: lex

  include "toric_namelists.h"

!this is all just to get NELM from toric input namelists or the machine file
  CALL getlun(inp_unit,ierr)
  CALL assert( ierr == 0, 'cannot find free i/o unit', ierr )

  WRITE(*,*) 'do_toric init reading machine.inp'
  OPEN(unit=inp_unit, file='machine.inp', status='old', &
       form='formatted')
  INQUIRE(inp_unit, exist=lex)
  IF (lex) THEN
     READ(inp_unit, nml = toricainp)
     READ(inp_unit, nml = nonthermals)
  ELSE
     WRITE(*,*) &
          'machine.inp does not exist or there was a read error'
  END IF
  CLOSE(inp_unit)


	!------------------------------------------------------------------------------------
	!   Get command line argument
	!------------------------------------------------------------------------------------

      call get_arg_count(iarg)
      if (iarg .ne. 1) then
         print*, 'model_RF_IC usage: '
         print*, ' command line arg = cur_state_file'
         call exit(1)
      end if
      
      call getarg(1,cur_state_file)
      
     WRITE (*,*)
     WRITE (*,*) 'load_mdescr_sconfig'      
     print*, 'cur_state_file = ', trim(cur_state_file)
      
	!------------------------------------------------------------------------------------
	!  Get current plasma state 
	!------------------------------------------------------------------------------------
			
    call ps_get_plasma_state(ierr, trim(cur_state_file))
    if (ierr .ne. 0) then
       print*, 'model_RF_IC:failed to get_plasma_state'
       call exit(1)
    end if


	!---------------------------------------------------------------------------------
	!  Check if nicrf_src dimensioned arrays are already allocated in current plasma
	!  state.  EPA component should have done this in its init.
	!---------------------------------------------------------------------------------
				 
	IF ( allocated(ps%icrf_src_name) ) THEN
		 print*, 'Checking ICRF data in current Plasma State'  
		 
		 !ps%freq_ic = freq_ic(1:ps%nicrf_src) ! EPA doesn't set source frequency
			
		 print*,   '   number of icrf sources = ', ps%nicrf_src
		 print*,   '   names of sources = ', ps%icrf_src_name
		 print*,   '   icrf frequencies = ', ps%freq_ic
		 print*,   '   icrf source powers = ', ps%power_ic
		 print*,   '   number of minority species = ', ps%nspec_rfmin
		 print*,   '   names of minority species = ', ps%RFMIN_name
		 print*,   '   eperp_icrf allocated = ', ALLOCATED(ps%eperp_mini)
		 print*,   '   epll_icrf allocated = ', ALLOCATED(ps%epll_mini)

	!-------------------------------------------------------------------------- 
	! Otherwise quit
	!--------------------------------------------------------------------------
		
	ELSE		  
		 print*, 'RF_IC sources NOT allocated in initial Plasma State'
		 print*, 'Can''t be done in this component without stepping on species lists'
		 print*, 'model_RF_IC_2.f90: I give up'
		 call exit(1)
					
	END IF


	!-------------------------------------------------------------------------- 
	! Get ICRF machine configuration data from machine description
	! namelist file.  Put into aux plasma state
	!--------------------------------------------------------------------------

	! for debugging, remove later
	CALL ps_store_plasma_state(ierr, 'initial_aux_ps.cdf', aux)

	call ps_mdescr_namelist_read(.False., 'mdescr_namelist.dat', ' ',  ' ', aux, ierr)
	IF (ierr .ne. 0) THEN
		print*, 'Could not get namelist mdescr'
		call exit(1)
	END IF
	
	! Check that nicrf_src read from namelist is consistent with that from current
	! plasma state.
	
	IF (ps%nicrf_src .ne. aux%nicrf_src) THEN
		print*, 'nicrf_src read from namelist is inconsistent with that from current&
			& plasma state'
		call exit(1)
	END IF		
		 print*,   '   number of icrf sources = ', ps%nicrf_src

	!-------------------------------------------------------------------------- 
	! Get ICRF shot configuration data from shot configuration description
	! namelist file.  Put into aux plasma state
	!--------------------------------------------------------------------------

	call ps_sconfig_namelist_read(.False., 'sconfig_namelist.dat', ' ', aux, ierr)
	IF (ierr .ne. 0) THEN
		print*, 'Could not get namelist sconfig'
		call exit(1)
	END IF

		 print*,   'aux%icrf_src_name = ', aux%icrf_src_name	  
		 print*,   'aux%ant_model = ', aux%ant_model		  
		 print*,   'aux%nrz_antgeo = ', aux%nrz_antgeo		  
		 print*,   'aux%freq_ic = ', aux%freq_ic		  
		 print*,   'aux%max_num_nphi = ', aux%max_num_nphi		  
		 print*,   'aux%nphi = ', aux%nphi		  
		 print*,   'aux%Real_ant_coef = ', aux%Real_ant_coef		  

	! for debugging, remove later
	CALL ps_store_plasma_state(ierr, 'initial_ps.cdf', ps)
	CALL ps_store_plasma_state(ierr, 'aux_ps.cdf', aux)
	
	!--------------------------------------------------------------------------
	! Copy Zimp1 and Zimp2 from ps to aux to avoid these being clobbered in the
	! plasma state copy from aux to ps. Same with ps%power_ic
	!--------------------------------------------------------------------------
	
	aux%Zimp1 = ps%Zimp1
	aux%Zimp2 = ps%Zimp2
	aux%power_ic = ps%power_ic
	
	!--------------------------------------------------------------------------
	! Copy data in RF_IC sections of 'aux' state to current (i.e. 'ps') state. 
	!--------------------------------------------------------------------------
	 
	 CALL ps_cclist_remove("*", cclist, iwarn)
	 CALL ps_cclist_add("IC" ,cclist, iwarn)
	 
	 CALL PS_COPY_DIMS(aux, ps, 0 ,ierr, cclist, 0)
	 if (ierr .ne. 0) then
		 print*, 'call failed to PS_COPY_DIMS for aux state to ps state'
		 call exit(1)
	 end if
	 CALL PS_COPY_PLASMA_STATE(aux, ps, ierr, cclist, 1, 0)
	 if (ierr .ne. 0) then
		 print*, 'call failed to PS_COPY_PLASMA_STATE for aux state to ps state'
		 call exit(1)
	 end if

		 print*,   'ps%icrf_src_name = ', ps%icrf_src_name		  
		 print*,   'ps%ant_model = ', ps%ant_model		  
		 print*,   'ps%nrz_antgeo = ', ps%nrz_antgeo		  
		 print*,   'ps%freq_ic = ', ps%freq_ic		  
		 print*,   'ps%max_num_nphi = ', ps%max_num_nphi		  
		 print*,   'ps%nphi = ', ps%nphi		  
		 print*,   'ps%Real_ant_coef = ', ps%Real_ant_coef		  
		 print*,   'ps%icrf_source powers = ', ps%power_ic
		

!added by DBB Feb 2008
  !----------------------------------
  !  form combined species list -- some arrays e.g. RF power coupling will
  !    want to be defined over the combined list.
      
  !  print out the species...
  ! PTB LAB removed this because it is no longer needed
  !       
  ! CALL ps_merge_species_lists(ierr)
  ! CALL assert((ierr == 0),' ?swim_state_test: ps_merge_species_lists: ierr= ', ierr)

  !  Print out combined species list...
! Added by PTB on 6-28-2012
! Set the  ICRF minority ion species parameters based on input data from the machine.inp
! file and then modify the "ALL", the "ALLA" and the "RFMIN" species lists.

! General RFMIN parameters
  ps%fracmin = fracmin 
  ps%isThermal = isThermal
  ps%kdens_rfmin = kdens_rfmin

! RFMIN Species list
  ps%m_RFMIN = m_rfmin * ps_mp
  ps%qatom_RFMIN = qatom_rfmin * ps_xe
  ps%q_RFMIN = q_rfmin * ps_xe
  ps%RFMIN_name = trim(rfmin_name)

! write(*,*) "rfmin_name =", rfmin_name
! write(*,*) "ps%RFMIN_name =", ps%RFMIN_name 

! "ALLA" Abridged Species list
  ps%m_ALLA(ps%rfmin_to_alla) = m_rfmin * ps_mp
  ps%qatom_ALLA(ps%rfmin_to_alla) = qatom_rfmin * ps_xe
  ps%q_ALLA(ps%rfmin_to_alla) = q_rfmin * ps_xe
  ps%ALLA_name(ps%rfmin_to_alla) = trim(rfmin_name)

! write(*,*) "rfmin_name =", rfmin_name
! write(*,*) "ps%rfmin_to_alla =",ps%rfmin_to_alla 
! write(*,*) "ps%ALLA_name =", ps%ALLA_name(ps%rfmin_to_alla) 

! "ALL" species list
  ps%m_ALL(ps%rfmin_to_all) = m_rfmin * ps_mp
  ps%qatom_ALL(ps%rfmin_to_all) = qatom_rfmin * ps_xe
  ps%q_ALL(ps%rfmin_to_all) = q_rfmin * ps_xe
  ps%ALL_name(ps%rfmin_to_all)= trim(rfmin_name)

! write(*,*) "rfmin_name =", rfmin_name
! write(*,*) "ps%ALL_name =", ps%ALL_name(ps%rfmin_to_all) 

   do i=0,ps%nspec_alla
!     write(*,1001) i,ps%all_type(i), &
!          trim(ps%all_name(i)),ps%q_all(i),ps%m_all(i)
     write(*,1001) i,ps%alla_type(i), &
          trim(ps%alla_name(i)),ps%q_alla(i),ps%m_alla(i)
  enddo
1001 format(' Species index & type: ',i2,1x,i2,1x, &
             '"',a,'" charge & mass: ',2(1pe12.5,1x))

 !----------------------------------



!create plasma state rf fields
  ps%freq_ic(1)=freqcy  !who sets #icrf_srcs, ps%icrf_src_name?
  print *,"freq_ic, picrf  alloc?",allocated(ps%freq_ic),allocated(ps%picrf_srcs)
  if (allocated(ps%picrf_srcs) .eqv. .TRUE.) then
      print *,"RF alloc?",allocated(ps%picrf_srcs),size(ps%picrf_srcs)
      print*,   '   number of icrf sources = ', ps%nicrf_src
  endif

  print *,  '   number of species = ', ps%nspec_alla, ps%nspec_th

  if (allocated(ps%rho_icrf) .eqv. .TRUE.) then
      print*,   '   radial grid points for icrf = ', ps%nrho_icrf
      print*,   '   ps%rho_icrf = ', allocated(ps%rho_icrf), size(ps%rho_icrf), ps%rho_icrf
  endif

!-------------------------------------------------------------------------- 
! Allocate rho_icrf dimensioned profile arrays in plasma state.
! nrho_icrf is in /static_state_data/ namelist.  For a real RF component
! it would be in standard TORIC or AORSA namelist.
!--------------------------------------------------------------------------
		 
IF ( allocated(ps%rho_icrf) ) THEN

	 print*, "model_RF_IC: rho_icrf(:) already allocated in plasma state &
	 &        Error - I'm supposed to do that"
	 call exit(1)
ELSE
	ps%nrho_icrf = nelm  
	CALL ps_alloc_plasma_state(ierr)
	WRITE (*,*) 'model_RF_IC: Allocated icrf profiles in plasma state'      
END IF

  do i=1,ps%nrho_icrf
     ps%rho_icrf(i) = real(i-1,rspec)
  end do
  ps%rho_icrf=ps%rho_icrf/real(ps%nrho_icrf-1,rspec)
!  ps%picrf_srcs = 0.0_rspec

!write the state file to optional filename, can also take optional state
  CALL ps_store_plasma_state(ierr , trim(cur_state_file))
  CALL assert( ierr == 0, 'cannot open state in do toric init', ierr )

END PROGRAM do_toric_init
