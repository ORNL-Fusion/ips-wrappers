!SF Cleanup 2023
!
program process_toric_output
  USE PLASMA_STATE_MOD
  USE SWIM_GLOBAL_DATA_MOD, ONLY :  &
    rspec, ispec

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  !internal params
  integer, parameter :: swim_string_length = 256

  !internal vars
  INTEGER :: &
       ierr, i, nedg, nprodt, nproeq
  REAL(rspec):: &
       dpsi
  
  !cdf vars
  INTEGER :: &
       istatus, ncid, vid, &
       npsi_id, nspec_id, nmhd_id, &
       npsi, nspec, nmhd

  REAL(rspec), allocatable, dimension(:) :: &
       psi_mesh, pow_prof_tot_elec, specif_volume, &
       pwr_IC_e, ps_psi
       
  REAL(rspec), allocatable, dimension(:,:) :: &
       pow_prof_ICfund_ions, pow_prof_ICharm_ions, &
       pwr_IC_ions
  
  !namelist vars
  character (len =swim_string_length) :: cur_state_file = "ips-state.nc"
  CHARACTER(len=16) :: toric_output_file = "fort.9"

  !namelist
  namelist /toric_process_nml/ &
       cur_state_file, toric_output_file
  
  ! --- Begin Program Logic ---------

  ! read program input namelist
  !--------------------------------------------------------------------
  open (unit=21, file="toric_process.nml", status='old', &
       form='formatted', iostat=ierr)
  

  ! Open the plasma state
  !--------------------------------------------------------------------

  CALL ps_get_plasma_state(ierr, cur_state_file)
  if(ierr.ne.0) then
     write(*,*)'process_toric_output: ps_get_plasma_state: ierr=' &
                      ,ierr
     stop 2
  end if
  IF (ierr .ne. 0) THEN
     CALL SWIM_error ('open', 'process_toric_output.f90', &
          'toric_process.nml')
     WRITE (*,*) 'process_toric_output.f90: Cannot open ', &
          TRIM('toric_process.nml')
     call exit(1)
  END IF
  READ(21,toric_process_nml)
  CLOSE (21)
  WRITE (*, toric_process_nml)

  ! Open the toric output file and read some values
  !--------------------------------------------------------------------

  istatus = nf_open(trim(toric_output_file),nf_nowrite,ncid)
  write(*,*)'after nf_open ncid=',ncid,'istatus',istatus

  !read dimension ids
  istatus = nf_inq_dimid(ncid,'n_of_rad_pts',npsi_id)
  write(*,*)'after ncdid npsi_id',npsi_id,'istatus',istatus

  istatus = nf_inq_dimid(ncid,'n_of_ionspec',nspec_id)
  write(*,*)'after ncdid nspec_id',nspec_id,'istatus',istatus

  !read dimensions
  istatus = nf_inq_dimlen(ncid, npsi_id, npsi)
  write(*,*)'proc_toric_op: after ncdinq, # of rad bins= ',npsi, &
            '  istatus=',istatus

  istatus = nf_inq_dimlen(ncid, nspec_id, nspec)
  write(*,*)'proc_toric_op: after ncdinq, # of ion spec.= ',nspec, &
            '  istatus=',istatus

  !allocate space for arrays to be read in
  allocate(psi_mesh(npsi))
  allocate(pow_prof_ICfund_ions(nspec,npsi))
  allocate(pow_prof_ICharm_ions(nspec,npsi))
  allocate(pow_prof_tot_elec(npsi))
  allocate(specif_volume(npsi))

  !read in arrays
  istatus = nf_inq_varid(ncid, 'psi_mesh', vid)
  istatus = nf_get_var_double(ncid, vid, psi_mesh)
  WRITE(*,*) 'proc_toric_op: after psi_mesh read. istatus = ', istatus
  
  istatus = nf_inq_varid(ncid, 'pow_prof_ICfund_ions', vid)
  istatus = nf_get_var_double(ncid, vid, pow_prof_ICfund_ions)
  WRITE(*,*) 'proc_toric_op: after pow_ic_fund read. istatus = ', istatus
  
  istatus = nf_inq_varid(ncid, 'pow_prof_ICharm_ions', vid)
  istatus = nf_get_var_double(ncid, vid, pow_prof_ICharm_ions)
  WRITE(*,*) 'proc_toric_op: after pow_ic_harm read. istatus = ', istatus
  
  istatus = nf_inq_varid(ncid, 'pow_prof_tot_elec', vid)
  istatus = nf_get_var_double(ncid, vid, pow_prof_tot_elec)
  WRITE(*,*) 'proc_toric_op: after pow_elec read. istatus = ', istatus
  
  istatus = nf_inq_varid(ncid, 'specif_volume', vid)
  istatus = nf_get_var_double(ncid, vid, specif_volume)
  WRITE(*,*) 'proc_toric_op: after specif_vol read. istatus = ', istatus
  
  !find edge of the the grid need 
  DO i = 1,npsi
     IF(psi_mesh(i) <= 1.0)then
        nedg = i
     ENDIF
  ENDDO
  dpsi = 1.0_rspec/nedg !specif_volume = volume/dpsi (for some reason) so need this

  !derived quantity grids
  allocate(pwr_IC_ions(nspec,nedg))
  allocate(pwr_IC_e(nedg))

  !the toric grid is so fine the difference between zone v grid center very small
  !not going to worry about it here. somebody else can fix this if they really want
  !to in the future
  DO i = 1,nspec
     pwr_IC_ions(i,1:nedg) = dpsi*specif_volume(1:nedg)*(pow_prof_ICfund_ions(i,1:nedg)+ &
          pow_prof_ICharm_ions(i,1:nedg))
  ENDDO
  pwr_IC_e = specif_volume(1:nedg)*pow_prof_tot_elec(1:nedg)

  !get plasma state poloidal flux grid so we can convert toric to rho
  nproeq = size(ps%rho_eq)
  nprodt = size(ps%rho)

  allocate(ps_psi(nprodt))
  ps_psi = 0.0_rspec
  call ps_user_1dintrp_vec(ps%rho, ps%rho_eq, ps%psipol, ps_psi, ierr )
  ps_psi = sqrt(ps_psi/ps_psi(nprodt))

  !rezone power onto grids in plasma state
  psi_mesh(1) = 0.0
  !electrons
  call ps_user_rezone1(psi_mesh(1:nedg), ps_psi, pwr_IC_e, &
       ps%picrf_srcs(:,1,i), ierr, nonorm = .TRUE., zonesmoo = .TRUE.)

  !ions
  do i = 1,ps%nspec_alla
     call ps_user_rezone1(psi_mesh(1:nedg), ps_psi, pwr_IC_ions(i,:), &
          ps%picrf_srcs(:,1,i), ierr, nonorm = .TRUE., zonesmoo = .TRUE.)
  enddo

  
  call ncclos(ncid,istatus)

  CALL ps_store_plasma_state(ierr)

  deallocate(psi_mesh)
  deallocate(pow_prof_ICfund_ions)
  deallocate(pow_prof_ICharm_ions)
  deallocate(pow_prof_tot_elec)
  deallocate(specif_volume)
  deallocate(pwr_IC_ions)
  deallocate(pwr_IC_e)
  deallocate(ps_psi)
   
end program process_toric_output
