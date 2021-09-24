PROGRAM do_torlh_init
  !JPL adapted from do_toric_init in June 2016
  !This is the initialization program called by rf_ic_torlh.py init step
  !It sets up the rf component parts of the plasma state.  These have been
  !moved from process_torlh_output
  
  ! Notes DBB 8/12/2016
  ! Typically one would use this in an IPS simulation in which the EPA component would have
  ! allocated all the arrays for thermal species and also any RF minorities in the plasma
  ! state.  The TORIC version of the component included the ability to allocate these 
  ! profiles here. However there were some problems due to conflicting use of the variable
  ! 'nspec' pulled out of the torica.inp file, and possible confusions due to this code 
  ! overwriting profile data in the plasma state with minority data pulled out of the 
  ! machine.inp file.  The modification eliminates the reading of torica.inp altogether
  ! but retains the ablility to allocate and initialize the minority data from the 
  ! machine.inp file.  Specifically:
  
  ! 1) A new variable 'nspec_rf_min' is added to the '/nonthermals/' namelist.
  ! 2) If 
  
! Plasma state data for LH component
! 	MACHINE_DESCRIPTION
! 		LHRF_SRC_NAME:  number & name of LHRF sources 
! 		NLHRF_SRC: item list dimension of lhrf_src_name (lhrf_source)
!  
! 	SHOT_CONFIGURATION 
! 		FREQ_LH: frequency on each LHRF source
!  
! 	SIMULATION_INIT 
! 		LH_CODE_INFO: Information: code implementing LH component 
! 		LH_DATA_INFO: Information on source of LHRF power data 
! 		NRHO_LHRF: grid dimension of rho_lhrf (RHO coordinate) 
! 		RHO_LHRF: rho grid -- LHRF
!  
! 	STATE_DATA
! 		 POWER_LH: power on each LHRF source
!  
! 	STATE_PROFILES 
! 		CURLH: LH current drive 
! 		CURLH_SRC: LH current drive (by antenna) 
! 		PELH: electron heating by LH 
! 		PELH_SRC: LH electron heating by antenna 
! 		PILH: ion heating by LH
! 		PILH_SRC: LH ion heating by antenna

! This code is required to set NRHO_LHRF and initialize RHO_LHRF.  The machine description
! and shot configuration data should have already been initialized either by the 
! generic_ps_init component or by the EPA component.

  USE plasma_state_mod

  USE swim_global_data_mod, only : &
       & rspec, ispec, &                ! int: kind specification for real and integer
       & swim_error                     ! error routine

  USE toric_utils_mod


  IMPLICIT NONE

  integer, parameter :: swim_string_length = 256  !for compatibility LAB
  character(len =swim_string_length) :: cur_state_file, program_name

  INTEGER:: inp_unit, ierr, iarg, i, lun
  LOGICAL:: lex

  include "torlh_namelists.h"

!this is all just to get NELM from torlh input namelists or the machine file
  CALL getlun(inp_unit,ierr)
  CALL assert( ierr == 0, 'cannot find free i/o unit', ierr )

  WRITE(*,*) 'do_torlh init reading machine.inp'
  OPEN(unit=inp_unit, file='machine.inp', status='old', &
       form='formatted')
  INQUIRE(inp_unit, exist=lex)
  IF (lex) THEN
     READ(inp_unit, nml = toricainp)
  ELSE
     WRITE(*,*) &
          'machine.inp does not exist or there was a read error'
  END IF
  CLOSE(inp_unit)

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
         write(0,*) 'do_torlh_init: '
    write(0,*) 'do_torlh-init cur_state_file'
    stop 'incorrect command line arguments'

  end select

  print*, 'do_torlh_init: cur_state_file = ', trim(cur_state_file)
  call get_arg(0,program_name)
  print*,'program ',program_name
  CALL ps_get_plasma_state(ierr, trim(cur_state_file))
  CALL assert( ierr==0,' initialize torlh: ps_get_plasma_state: ierr=',ierr )

    if(.not. allocated(ps%freq_lh)) then
       write(*,*)'LRF_SRC not allocated in initial plasma state'
       ps%nlhrf_src = 1
       ps%lhrf_src_name = 'LH_1'
    
       write(*,*)'Allocating LH source in do_torlh_init - default one source'
       CALL ps_alloc_plasma_state(ierr)
    
       if(ierr.ne.0) then
           write(*,*) trim(cur_state_file),' ps_alloc_plasma_state: ierr = ', &
           ierr,  'allocating rf component in do_torlh_init '
           stop
       endif
    
    endif

   do i=0,ps%nspec_alla
     write(*,*) 'Spec Index : ', i
     if(allocated(ps%alla_type)) write(*,*) 'Spec Type : ', ps%alla_type(i)
     if(allocated(ps%alla_name)) write(*,*) 'Spec Name : ', trim(ps%alla_name(i))
     if(allocated(ps%q_alla)) write(*,*) 'Spec q : ', ps%q_alla(i)
     if(allocated(ps%m_alla)) write(*,*) 'Spec m : ', ps%m_alla(i)
  enddo

!create plasma state rf fields
!   ps%freq_lh(1)=freqcy  !who sets #lh_srcs, ps%lh_src_name?
!   print *,"freq_lh, picrf  alloc?",allocated(ps%freq_lh),allocated(ps%pelh_src)
  if (allocated(ps%pelh_src) .eqv. .TRUE.) then
      print *,"RF alloc?",allocated(ps%pelh_src),size(ps%pelh_src)
      print*,   '   number of lower hybrid wave sources = ', ps%nlhrf_src
  endif

  print *,  '   number of species = ', ps%nspec_alla, ps%nspec_th

  if (allocated(ps%rho_lhrf) .eqv. .TRUE.) then
      print*,   '   radial grid points for LH waves = ', ps%nrho_lhrf
      print*,   '   ps%rho_lh = ', allocated(ps%rho_lhrf), size(ps%rho_lhrf), ps%rho_lhrf
  endif

! Check to see if the RF arrays are already allocated
!
  if (allocated(ps%epll_mini) .eqv. .FALSE.) then
     print *,'allcoating plasma state rf arrays'
     ps%nrho_lhrf = nelm
     CALL ps_alloc_plasma_state(ierr) 
     CALL assert( ierr == 0, trim(cur_state_file)//' ps_alloc_plasma_state: ierr=',ierr )
  endif

  do i=1,ps%nrho_lhrf
     ps%rho_lhrf(i) = real(i-1,rspec)
  end do
  ps%rho_lhrf=ps%rho_lhrf/real(ps%nrho_lhrf-1,rspec)
  ps%pelh_src = 0.0_rspec

!write the state file to optional filename, can also take optional state
  CALL ps_store_plasma_state(ierr , trim(cur_state_file))
  CALL assert( ierr == 0, 'cannot open state in do torlh init', ierr )

END PROGRAM do_torlh_init
