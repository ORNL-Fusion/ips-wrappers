      program prepare_abj_input
      USE plasma_state_mod

      use swim_global_data_mod, only : &
           & rspec, ispec, &
           & swim_error

      implicit none
      
      !Parameters
      INTEGER, PARAMETER :: PAR = 1, PERP = 2, LOWER = 1, UPPER = 2
      integer, parameter :: swim_string_length = 256
      INTEGER, PARAMETER:: FNLEN = 50, SPEC_NMLEN = 25
      INTEGER, PARAMETER:: NSPEC_MAX = 4

      !program vars
      character(len =swim_string_length) :: cur_state_file, cur_geq_file
      character *(25):: name_minspec,arg_enorm,arg_nsurf,arg_specs
      double precision:: z_minspec, m_minspec,enorm
      INTEGER :: ierr,isp_min,iarg,inp_unit, out_unit=108
      logical :: lex, enorm_set, surf_set
      
      !abj driver namelist vars
      CHARACTER *(20):: &
      distribution, profile, profile_datafile
  
      
      !abj namelist vars
      CHARACTER *(FNLEN):: cdf_iofn
      CHARACTER *(SPEC_NMLEN),DIMENSION(NSPEC_MAX):: species_name
      INTEGER, DIMENSION(PAR:PERP):: npts
      INTEGER:: ntheta, debug_level
      INTEGER:: ngen, npsi
      DOUBLE PRECISION, DIMENSION(LOWER:UPPER):: v_asymptotic
      DOUBLE PRECISION:: Bmax_by_Bmin, velocity_normalization
      DOUBLE PRECISION, DIMENSION(PAR:PERP, LOWER:UPPER):: vlims
      DOUBLE PRECISION, DIMENSION(NSPEC_MAX):: Z_charge, A_mass 

      !namelists
      NAMELIST / ABJ_driver_nml / profile_datafile, distribution, profile 

      NAMELIST / ABJ_nml / cdf_iofn, ntheta, npts, vlims, v_asymptotic, &
              Bmax_by_Bmin, NPSI, NGEN, Z_Charge, A_mass, species_name, &
              velocity_normalization, debug_level

      !read input args
      call get_arg_count(iarg)
      WRITE(*,*) 'arg count', iarg
      select case (iarg)
      CASE(0)
         enorm_set = .False.
         cur_state_file="cur_state.cdf"
      CASE(1)
         call get_arg(1,cur_state_file)
      CASE(2)
         enorm_set = .True.
         call get_arg(1,cur_state_file)
         call get_arg(2,arg_enorm)
         read(arg_enorm,*) enorm
      CASE(3)
         enorm_set = .True.
         call get_arg(1,cur_state_file)
         call get_arg(2,arg_enorm)
         read(arg_enorm,*) enorm
         call get_arg(3,arg_specs)
      end select
      
      WRITE(*,*) 'cur_state enorm', cur_state_file, enorm
  
      !get plasma state data
      call ps_get_plasma_state(ierr,trim(cur_state_file))
      if(ierr .ne. 0) stop 'cannot get plasma state'

      !read ABJ namelists
      WRITE(*,*) 'Prepare ABJ input reading ABJ namelists'
      open(unit=inp_unit, file='ABJ.inp', status='old', &
           form='formatted')
      inquire(inp_unit,exist=lex)
      if (lex) then
         read(inp_unit, nml=ABJ_nml)
      else
         write(*,*) 'ABJ.inp does not exist or there was a read error'
         stop
      endif
      close(inp_unit)

      open(unit=inp_unit, file='ABJ_driver.inp', status='old', &
           form='formatted')
      inquire(inp_unit,exist=lex)
      if (lex) then
         read(inp_unit, nml=ABJ_driver_nml)
      else
         write(*,*) 'ABJ_driver.inp does not exist or there was a read error'
         stop
      endif
      close(inp_unit)
      
      !set namelist values

      !use plasma state values to set certain namelist params
      isp_min = ps%rfmin_to_alla(1)
      z_minspec = 1.d0*NINT(ps%q_alla(isp_min)/ps_xe)
      m_minspec = 1.d0*NINT(ps%m_alla(isp_min)/ps_mp)
      name_minspec = trim(ps%alla_name(isp_min))
      
      !abj_driver_nml
      profile_datafile = 'cql3d.nc'
      distribution = 'cql3d'
      profile = 'datafile'

      !abj_nml
      if(enorm_set)then
         velocity_normalization = 3.106208817D7*sqrt(0.5*enorm/m_minspec)/25.0
         vlims(PAR,LOWER) = -25.0
         vlims(PAR,UPPER) = 25.0
         vlims(PERP,LOWER) = 0.0
         vlims(PERP,UPPER) = 25.0
         v_asymptotic(LOWER) = -25.0
         v_asymptotic(UPPER) = 25.0
      endif

      debug_level = 0
      Bmax_by_Bmin = 2.0 !SF this is as high as this should ever need to be for an aspect ratio 3 or larger

      if (trim(arg_specs).eq.'MIN')then
         ngen = 1
         A_mass(1) = m_minspec
         Z_charge(1) = z_minspec
         species_name(1) = name_minspec
      endif

      if (trim(arg_specs).eq.'MIN+')then
         ngen = 2
         A_mass(1) = m_minspec
         Z_charge(1) = z_minspec
         species_name(1) = name_minspec
         A_mass(2) = 1.d0*NINT(ps%m_alla(1)/ps_mp)
         Z_charge(2) = 1.d0*NINT(ps%q_alla(1)/ps_xe)
         species_name(2) = trim(ps%alla_name(1))
      endif
         
      !write namelists to ABJ input files
      open(unit=out_unit, file = 'ABJ_driver.inp_new', delim = 'apostrophe', &
          status = 'unknown', form = 'formatted') 

      write(out_unit, nml=ABJ_driver_nml)
      close(out_unit)

      open(unit=out_unit, file = 'ABJ.inp_new', delim = 'apostrophe', &
          status = 'unknown', form = 'formatted') 
      write(out_unit, nml=ABJ_nml)
      close(out_unit)
      
      end program prepare_abj_input
