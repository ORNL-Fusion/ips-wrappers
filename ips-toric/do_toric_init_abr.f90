PROGRAM do_toric_init
  !JCW created APR 2008
  !This is the initialization program called by rf_ic_toric.py init step
  !It sets up the rf component parts of the plasma state.  These have been
  !moved from process_toric_output

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
     READ(inp_unit, nml = nonthermals)
  ELSE
     WRITE(*,*) &
          'machine.inp does not exist or there was a read error'
  END IF
  CLOSE(inp_unit)

    ! This is just to read in 'nspec'
    call getlun(lun,ierr)
    call assert(ierr == 0, 'cannot find free i/o unit', ierr )
    open(unit=lun, file='torica.inp', status='old', form='formatted')
    inquire(lun, exist=lex)
    if(lex)then
        write(*,*) 'TORIC INIT : reading torica.inp'
        read(lun, nml=equidata)
    else
        write(*,*)'torica.inp does not exist or could not be read'
    endif
    close(lun)

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
    
       write(*,*)'Allocating RF in prepare_input'
       CALL ps_alloc_plasma_state(ierr)
    
       if(ierr.ne.0) then
           write(*,*) trim(cur_state_file),' ps_alloc_plasma_state: ierr = ', &
           ierr,  'allocating rf component in do_toric_init '
           stop
       endif
    
    endif

    ! vars of size ps%nspec_rfmin
    if(.not.allocated(ps%fracmin))then
        if(nspec.le.0)then
            write(*,*)'TORIC INIT : ERROR - nspec read from torica.inp not valid'
            stop
        endif
        ps%nspec_rfmin = nspec
        write(*,*)'Setting nspec_rfmin = ',nspec
        CALL ps_alloc_plasma_state(ierr)
        CALL assert( ierr==0,' initialize toric A: ps_get_plasma_state: ierr=',ierr )
    endif

    write(*,*)'TORIC INIT : Setting nspec_rfmin sized vars'

    ps%fracmin = fracmin 
    ps%isThermal = isThermal
    ps%m_RFMIN = m_rfmin * ps_mp
    ps%qatom_RFMIN = qatom_rfmin * ps_xe
    ps%q_RFMIN = q_rfmin * ps_xe
    ps%RFMIN_name = trim(rfmin_name)

    write(*,*)'TORIC INIT : Setting kdens_rfmin scalar'

    ps%kdens_rfmin = kdens_rfmin

    ! vars of size ps%nspec_alla
    write(*,*)'TORIC INIT : Setting nspec_alla sized vars'
    if(.not.allocated(ps%m_ALLA))then
        if(nspec.le.0)then
            write(*,*)'TORIC INIT : ERROR - nspec read from torica.inp not valid'
            stop
        endif
 
        ps%nspec_alla = nspec
        write(*,*)'Setting nspec_alla = ',nspec
        CALL ps_alloc_plasma_state(ierr)
        CALL assert( ierr==0,' initialize toric B: ps_get_plasma_state: ierr=',ierr )
    endif
    ps%m_ALLA(ps%rfmin_to_alla) = m_rfmin * ps_mp
    ps%qatom_ALLA(ps%rfmin_to_alla) = qatom_rfmin * ps_xe
    ps%q_ALLA(ps%rfmin_to_alla) = q_rfmin * ps_xe
    ps%ALLA_name(ps%rfmin_to_alla) = trim(rfmin_name)

    ! vars of size ps%nspec_all
    write(*,*)'TORIC INIT : Setting nspec_all sized vars'
    if(.not.allocated(ps%m_ALL))then
        if(nspec.le.0)then
            write(*,*)'TORIC INIT : ERROR - nspec read from torica.inp not valid'
            stop
        endif
 
        ps%nspec_all = nspec
        write(*,*)'Setting nspec_all = ',nspec
        CALL ps_alloc_plasma_state(ierr)
        CALL assert( ierr==0,' initialize toric C: ps_get_plasma_state: ierr=',ierr )
    endif
    ps%m_ALL(ps%rfmin_to_all) = m_rfmin * ps_mp
    ps%qatom_ALL(ps%rfmin_to_all) = qatom_rfmin * ps_xe
    ps%q_ALL(ps%rfmin_to_all) = q_rfmin * ps_xe
    ps%ALL_name(ps%rfmin_to_all)= trim(rfmin_name)

   do i=0,ps%nspec_alla
     write(*,*) 'Spec Index : ', i
     if(allocated(ps%alla_type)) write(*,*) 'Spec Type : ', ps%alla_type(i)
     if(allocated(ps%alla_name)) write(*,*) 'Spec Name : ', trim(ps%alla_name(i))
     if(allocated(ps%q_alla)) write(*,*) 'Spec q : ', ps%q_alla(i)
     if(allocated(ps%m_alla)) write(*,*) 'Spec m : ', ps%m_alla(i)
  enddo

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

!  if (ps%nspec==0) then
!     print*,'Error, nspec not set.  Setting it to value of nspec_th'
!     ps%nspec=ps%nspec_th
!  end if

! Check to see if the RF arrays are already allocated
!
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
