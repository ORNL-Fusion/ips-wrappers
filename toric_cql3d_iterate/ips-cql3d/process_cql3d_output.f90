!SF Cleanup 2023
!
      program process_cql3d_output

      USE PLASMA_STATE_MOD
      USE SWIM_GLOBAL_DATA_MOD, ONLY :  
     & rspec, ispec

      IMPLICIT NONE

      INCLUDE 'netcdf.inc'

      integer, parameter :: swim_string_length=256
      
      !error check integers
      integer ::
     & ierr, istatus

      !netcdf vars
      integer ::
     & ncid, vid,
     & lrz_id, nt_id, ngen_id, ntotal_id, r0dim_id, 
     & lrz, nt, ngen, ntotal, r0dim
      
      real(rspec), allocatable, dimension(:) ::
     & rya, darea, dvol

      real(rspec), allocatable, dimension(:,:,:) ::
     & curr, wperp, wpar

      real(rspec), allocatable, dimension(:,:,:,:) ::
     & powers

      character*8 radcoord
     
      !buffer vars for output to plasma state
      integer ::
     & i 
     
      real(rspec), allocatable, dimension(:)::
     & rho_cql,pelh_tmp, curlh_tmp,
     & picrf_tmp,pmine_tmp,pmini_tmp 

     
      ! program namelist variables
      !----------------------------------------------------------------

      character (len =swim_string_length) :: cur_state_file = "ips-state.nc"
      CHARACTER(len=16) :: cql3d_specs='MIN+'
      CHARACTER(len=16) :: cql3d_mode='IC'
      CHARACTER(len=16) :: cql3d_output_file='cql3d.nc'
      
      ! ---Begin Program Logic----------
      
      ! read program input namelist
      !----------------------------------------------------------------
      namelist /cql3d_process_nml/
     &    cur_state_file, cql3d_specs, cql3d_mode, cql3d_output_file 

      
      OPEN (unit=21, file="cql3d_process.nml", 
     &      status='old',form = 'formatted', iostat=ierr)
      IF (ierr .ne. 0) THEN
         CALL SWIM_error ('open', 'process_cql3d_output.f90',
     &                    'cql3d_process.nml')
         WRITE (*,*) 'process_cql3d_output.f90: Cannot open ',
     &    TRIM('cql3d_process.nml')
         call exit(1)
      END IF
      READ(21,cql3d_process_nml)
      CLOSE (21)
      WRITE (*, cql3d_process_nml)

      !trim inputs
      cur_state_file=trim(cur_state_file)
      cql3d_specs=trim(cql3d_specs)
      cql3d_mode=trim(cql3d_mode)


      ! Open the plasma state
      !----------------------------------------------------------------
      
      CALL ps_get_plasma_state(ierr, cur_state_file)
      if(ierr.ne.0) then
         write(*,*)'process_cql3d_output: ps_get_plasma_state: ierr='
     &                 ,ierr
         stop 2
      end if

      ! Open the cql3d output file and read dims/vars
      !----------------------------------------------------------------
      
      istatus = nf_open(trim(cql3d_output_file),nf_nowrite,ncid)
      write(*,*)'after nf_open ncid=',ncid,'istatus',istatus
      
      !read dimension ids
      istatus = nf_inq_dimid(ncid,'rdim',lrz_id)
      write(*,*)'after ncdid lrz_id',lrz_id,'istatus',istatus

      istatus = nf_inq_dimid(ncid,'tdim',nt_id)
      write(*,*)'after ncdid nt_id',nt_id,'istatus',istatus
      
      istatus = nf_inq_dimid(ncid,'gen_species_dim',ngen_id)
       write(*,*)'after ncdid ngen_id',ngen_id,'istatus',istatus
      
      istatus = nf_inq_dimid(ncid,'species_dim',ntotal_id)
      write(*,*)'after ncdid ntotal_id',ntotal_id,'istatus',istatus
      
      istatus =  nf_inq_dimid(ncid,'r0dim',r0dim_id)
      write(*,*)'proc_cql3d_op:after ncdid r0dim_id',r0dim_id,'istatus',istatus

      !read dimensions
      istatus = nf_inq_dimlen(ncid, lrz_id, lrz)
      write(*,*)'proc_cql3d_op: after ncdinq, # of rad bins= ',lrz,
     &     '  istatus=',istatus

      istatus = nf_inq_dimlen(ncid, nt_id, nt)
      write(*,*)'proc_cql3d_op: after ncdinq, # of t steps = ',nt,
     &     '  istatus=',istatus
      
      istatus = nf_inq_dimlen(ncid, ngen_id, ngen)
      write(*,*)'proc_cql3d_op: after ncdinq, #  gen specs = ', ngen,
     &     '  istatus=',istatus

      istatus = nf_inq_dimlen(ncid, ntotal_id, ntotal)
      write(*,*)'proc_cql3d_op: after ncdinq, #  species total = ',
     &    ntotal, '  istatus=',istatus

      istatus = nf_inq_dimlen(ncid, r0dim_id, r0dim)
      write(*,*)'proc_cql3d_op: after ncdinq, # of rad bins= ',r0dim,
     &     '  istatus=',istatus

      !allocate space for arrays to be read in
      allocate (rya(lrz))  !the radial grid
      allocate (darea(lrz),dvol(lrz))
      allocate (powers(lrz, 13, ngen, nt))
      allocate (wperp(lrz,ngen,nt),wpar(lrz,ngen,nt))
      allocate (curr(lrz,ngen,nt))
      
      istatus = nf_inq_varid(ncid, 'radcoord', vid)
      istatus = nf_get_var_text(ncid, vid, radcoord)
      if (radcoord.ne.'sqtorflx') then
         write(*,*)'STOP: problem with cql3d radial coord type, not sqtorflx'
         stop
      endif

      istatus = nf_inq_varid(ncid, 'rya', vid)
      istatus = nf_get_var_double(ncid, vid, rya)

      istatus = nf_inq_varid(ncid, 'darea', vid)
      istatus = nf_get_var_double(ncid, vid, darea)

      istatus = nf_inq_varid(ncid, 'dvol', vid)
      istatus = nf_get_var_double(ncid, vid, dvol)

      istatus = nf_inq_varid(ncid, 'powers', vid)            
      istatus = nf_get_var_double(ncid, vid, powers)

      istatus = nf_inq_varid(ncid, 'wperp', vid)
      istatus = nf_get_var_double(ncid, vid, wperp)

      istatus = nf_inq_varid(ncid, 'wpar', vid)
      istatus = nf_get_var_double(ncid, vid, wpar)

      istatus = nf_inq_varid(ncid, 'curr', vid)
      istatus = nf_get_var_double(ncid, vid, curr)
      
      if (cql3d_mode .eq. 'LH') then
         allocate(pelh_tmp(lrz-1))
         allocate(curlh_tmp(lrz-1))
         allocate(rho_cql(lrz-1))
         pelh_tmp = 0.0_rspec
         curlh_tmp = 0.0_rspec
         rho_cql = 0.0_rspec

         do i=1,lrz-1
            rho_cql(i) = 0.5*(rya(i)+rya(i+1))
            pelh_tmp(i) = 0.5*dvol(i)*(powers(i,5,1,nt)+powers(i+1,5,1,nt))
            curlh_tmp(i) = 0.5*darea(i)*(curr(i,1,nt)+curr(i+1,1,nt))
         enddo
         
         call ps_user_rezone1(rho_cql, ps%rho_lhrf, pelh_tmp,
     &    ps%pelh, ierr, nonorm = .TRUE., zonesmoo = .TRUE.)

         call ps_user_rezone1(rho_cql, ps%rho_lhrf, curlh_tmp,
     &    ps%curlh, ierr, nonorm = .TRUE., zonesmoo = .TRUE.)
         
         deallocate(pelh_tmp)
         deallocate(curlh_tmp)
         deallocate(rho_cql)

      elseif (cql3d_mode .eq. 'IC') then
         
         allocate(picrf_tmp(lrz+1))
         allocate(pmine_tmp(lrz+1))
         allocate(pmini_tmp(lrz+1))
         allocate(rho_cql(lrz+1))
         picrf_tmp = 0.0_rspec
         pmine_tmp = 0.0_rspec
         pmini_tmp = 0.0_rspec
         rho_cql = 0.0_rspec

         rho_cql(1)=0.0
         do i=1,lrz-1
            rho_cql(i+1) = 0.5*(rya(i)+rya(i+1))
         enddo
         rho_cql(lrz+1) = 1.0
         
         !total icrf power and power to minority
         if ((cql3d_specs.eq.'MIN').or.(cql3d_specs.eq.'MIN+')) then
            do i=1,lrz-1
               picrf_tmp(i+1) = 0.5*dvol(i)*(powers(i,5,1,nt)+powers(i+1,5,1,nt)) !icrf power abs
               pmine_tmp(i+1) = -0.5*dvol(i)*(powers(i,1,1,nt)+powers(i+1,1,1,nt)) !power min->e-
               pmini_tmp(i+1) = - 0.5*dvol(i)*(powers(i,2,1,nt)+powers(i+1,2,1,nt)) !power min->i+
            end do
         endif

         if (cql3d_specs.eq.'MIN+')then
            do i=1,lrz-1
               picrf_tmp(i+1) = picrf_tmp(i+1) + 0.5*dvol(i)*(powers(i,5,2,nt)+powers(i+1,5,2,nt))
               pmini_tmp(i+1) = pmini_tmp(i+1) -0.5*dvol(i)*(powers(i,3,1,nt)+powers(i+1,4,1,nt)) !power min->i+gen
            enddo
         endif
            
         do i=1,lrz+1
            if (pmine_tmp(i)<0) pmine_tmp(i)=0.0_rspec
            if (pmini_tmp(i)<0) pmini_tmp(i)=0.0_rspec
         enddo

         !may need to remap to zone center here?
         
         call ps_user_rezone1(rho_cql, ps%rho_icrf, pmine_tmp, 
     &        ps%pmine, ierr, nonorm = .TRUE., zonesmoo = .TRUE.)
         
         call ps_user_rezone1(rho_cql, ps%rho_icrf, pmini_tmp, 
     &        ps%pmini, ierr, nonorm = .TRUE., zonesmoo = .TRUE.)

         
         ps%picrf_abs(1) = sum(picrf_tmp)
         
         !WRITE(*,*) wperp(:,1,nt)
         !WRITE(*,*) wpar(:,1,nt)

         !WRITE(*,*) shape(ps%rho_icrf)
         !WRITE(*,*) shape(rya)
         !WRITE(*,*) shape(wpar(:,1,nt))
         !WRITE(*,*) shape(ps%epll_mini(:,1))
         
         
         !effective temperatures (breaking the heap sometimes for some reason)
!         call ps_user_1dintrp_vec(ps%rho_icrf, rya, wperp(:,1,nt), 
!     &        ps%eperp_mini(:,1), ierr)
!         call ps_user_1dintrp_vec(ps%rho_icrf, rya, wpar(:,1,nt), 
!     &        ps%epll_mini(:,1), ierr)
         
         
         !debug check
         WRITE(*,*) 'pwrrf_ic: ', SUM(ps%picrf_abs)
         WRITE(*,*) 'pmini:    ', SUM(ps%pmini)
         WRITE(*,*) 'pmine:    ', SUM(ps%pmine)
     
      endif

      call ncclos(ncid,istatus)

      CALL ps_store_plasma_state(ierr)

      deallocate(picrf_tmp)
      deallocate(pmine_tmp)
      deallocate(pmini_tmp)
      deallocate (rya)
      deallocate (darea)
      deallocate (dvol)
      deallocate (powers)
      deallocate (wperp)
      deallocate (wpar)
      deallocate (curr)
      
      end program process_cql3d_output
