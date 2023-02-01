PROGRAM do_toric_init
  !JCW created APR 2008
  !This is the initialization program called by rf_ic_toric.py init step
  !It sets up the rf component parts of the plasma state.  These have been
  !moved from process_toric_output

  !SF 12/2022
  !Removed allocation of minority species here. This should be done in the plasma
  !state setup. Not in the TORIC setup.
  
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
         write(0,*) 'do_toric_init: '
    write(0,*) 'do_toric-init cur_state_file'
    stop 'incorrect command line arguments'

  end select

  print*, 'do_toric_init: cur_state_file = ', trim(cur_state_file)
  call get_arg(0,program_name)
  print*,'program ',program_name
  CALL ps_get_plasma_state(ierr, trim(cur_state_file))
  CALL assert( ierr==0,' initialize toric: ps_get_plasma_state: ierr=',ierr )

    if(.not. allocated(ps%freq_ic)) then
       write(*,*)'RF SRC NOT ALLOCATED'
       ps%nicrf_src = 1
       ps%freq_ic=freqcy  !who sets #icrf_srcs, ps%icrf_src_name?
       write(*,*)'Allocating RF in prepare_input'
       CALL ps_alloc_plasma_state(ierr)
    
       if(ierr.ne.0) then
           write(*,*) trim(cur_state_file),' ps_alloc_plasma_state: ierr = ', &
           ierr,  'allocating rf component in do_toric_init '
           stop
       endif
    
    endif

!create plasma state rf fields
  if (allocated(ps%epll_mini) .eqv. .FALSE.) then
     print *,'allcoating plasma state rf arrays'
     ps%nrho_icrf = nelm
     CALL ps_alloc_plasma_state(ierr) 
     CALL assert( ierr == 0, trim(cur_state_file)//' ps_alloc_plasma_state: ierr=',ierr )
  endif

  do i=1,ps%nrho_icrf
     ps%rho_icrf(i) = real(i-1,rspec)
  end do
  ps%rho_icrf=ps%rho_icrf/real(ps%nrho_icrf-1,rspec)
  ps%picrf_srcs = 0.0_rspec

!write the state file to optional filename, can also take optional state
  CALL ps_store_plasma_state(ierr , trim(cur_state_file))
  CALL assert( ierr == 0, 'cannot open state in do toric init', ierr )

END PROGRAM do_toric_init
