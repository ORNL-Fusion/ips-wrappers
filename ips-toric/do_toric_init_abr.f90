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

  INTEGER:: inp_unit, ierr, iarg, i
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
   write(iout, *)'RF SRC NOT ALLOCATED'
   ps%nicrf_src = 1

   print *,'Allocating RF in prepare_input'
   CALL ps_alloc_plasma_state(ierr)

   if(ierr.ne.0) then
       write(iout, *) trim(cur_state_file),' ps_alloc_plasma_state: ierr = ', &
       ierr,  'allocating rf component in do_toric_init '
       stop
   endif
   !initing the pieces of rf that will be used every time
!   ps%icrf_src_name(1) = 'TORIC'
!   ps%freq_ic(1) = 80.0e6     
!   ps%dist_fun = 'rf_min'!wired for CMOD minority
   !how to do multiple sources needs to be FIXED
endif

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
