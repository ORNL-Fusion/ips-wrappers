! Code to check parareal convergence for SOLPS - SD
!Code checks convergence for both pwmxip & pwmxap - SD (9Dec2013)
 program PR_conv_SOLPS
  use netcdf
  implicit none
 ! include <netcdf.inc>
  !Variables for command line inputs
  integer :: n_args
  ! This is the name of the data file we will create and those read in from the command line
  character (len = 1680) :: FILE_NAME_Conv
  character (len = 1680) :: FILE_NAME_Fold
  character (len = 1680):: FILE_NAME_Fnew
  character(len = 20) :: err_tol !error tolerance read as character

 ! When we create netCDF files, variables and dimensions, we get back
  ! an ID for each one.
  integer :: ncid_Fold,ncid_Fnew !ncid is the ID of the file opened
  integer :: es_h_id_Fold, es_h_id_Fnew !ID of the variables
  integer :: es_h_id_Fold2, es_h_id_Fnew2 !ID of the variables
  character (len = 19) :: species_var,time_var,var_name !names of dimensions(probably optional)
  character (len = 19) :: species_var2,time_var2,var_name2 !names of dimensions(probably optional)
  integer:: species_len, time_len !length of the dimensions
  integer:: species_len2, time_len2 !length of the dimensions
  integer::xtype, ndims, xtype2, ndims2 !type and number of dimensions in energy variable in each file
  integer:: dimids(5), dimids2(5)

  double precision, allocatable, dimension(:,:):: flux_in_old, flux_in_new !Variables used from file
  double precision, allocatable, dimension(:,:):: flux_in_old2, flux_in_new2 !Variables used from file
  double precision, allocatable, dimension(:) :: err,err_av_time,err2,err_av_time2
  double precision :: err_tolerance,err_av,err_av2,err_av_final
  integer:: i,time

!Get command line arguments
  n_args = command_argument_count() !number of arguments from command line
      call get_command_argument(1,FILE_NAME_Fold)
      call get_command_argument(2,FILE_NAME_Fnew)
      call get_command_argument(3,FILE_NAME_Conv)
      call get_command_argument(4, err_tol)

 !Convert characters to integers- SD
  read(err_tol,'(f24.20)') err_tolerance

!Open output file:
  open(143,file= trim(FILE_NAME_Conv),status='unknown')  

!Open existing file only.NF90_WRITE gives read-write access.-SD
  call check( nf90_open(FILE_NAME_Fold, NF90_WRITE, ncid_Fold) )
  call check( nf90_open(FILE_NAME_Fnew, NF90_WRITE, ncid_Fnew) )

 !Get the varid of the data variable from F run, based on its name. -SD
  call check( nf90_inq_varid(ncid_Fold, "pwmxip", es_h_id_Fold) )
  call check( nf90_inq_varid(ncid_Fnew, "pwmxip", es_h_id_Fnew) )
  call check( nf90_inq_varid(ncid_Fold, "pwmxap", es_h_id_Fold2) )
  call check( nf90_inq_varid(ncid_Fnew, "pwmxap", es_h_id_Fnew2) )

  write(6,*) "pwmxip id found"
!Get details of pwmxip & pwmxap variables from just one input file. The dimensions will be same for other input file.
  call check(nf90_inquire_variable(ncid_Fold, es_h_id_Fold, var_name,xtype, ndims, dimids))
  call check(nf90_inquire_variable(ncid_Fold, es_h_id_Fold2, var_name2,xtype2, ndims2, dimids2))
  write(6,*) "var_name",var_name,"xtype",xtype,"ndims",ndims
!Note: It is known from ncdump of data that it is a ** dimensional array. So, dimids(5) is redundant - could just have left it as dimids(2)
  write(6,*) "IDs of dimensions of flux variable",dimids(1),"&",dimids(2)

!Get length of dimension of data
  call check(nf90_inquire_dimension(ncid_Fold,dimids(1),species_var,species_len ))
  call check(nf90_inquire_dimension(ncid_Fold,dimids(2),time_var,time_len ))
  call check(nf90_inquire_dimension(ncid_Fold,dimids2(1),species_var2,species_len2 ))
  call check(nf90_inquire_dimension(ncid_Fold,dimids2(2),time_var2,time_len2 ))

   write(6,*) "species & time_index lengths:",species_var,":",species_len,time_var,":",time_len

   !allocate array
   allocate(flux_in_old(species_len,time_len),flux_in_new(species_len,time_len))
   allocate(err(species_len),err_av_time(species_len))
   allocate(flux_in_old2(species_len2,time_len2),flux_in_new2(species_len2,time_len2))
   allocate(err2(species_len2),err_av_time2(species_len2))

  ! Read the data from F run.
   call check( nf90_get_var(ncid_Fold,es_h_id_Fold,flux_in_old))
   call check( nf90_get_var(ncid_Fnew,es_h_id_Fnew,flux_in_new))
   call check( nf90_get_var(ncid_Fold,es_h_id_Fold2,flux_in_old2))
   call check( nf90_get_var(ncid_Fnew,es_h_id_Fnew2,flux_in_new2))
   err_av=0.
   do i=1,species_len
        err(i)=0.
        do time=1,time_len
!		write(6,*) "heat_flux at time_index",time,"=",energy_tmte_old(i,time)
!		write(6,*) "heat_flux at time_index",time,"=",energy_tmte_new(i,time)
                err(i) = err(i) + ((flux_in_new(i,time) - flux_in_old(i,time))/(flux_in_old(i,time)+1.0E-32))             
	end do
         err_av_time(i)=err(i)/time_len
         err_av = err_av + err_av_time(i)
   end do
   err_av = abs(err_av/species_len)

   err_av2=0.
   do i=1,species_len2
        err2(i)=0.
        do time=1,time_len2
                err2(i) = err2(i) + ((flux_in_new2(i,time) - flux_in_old2(i,time))/(flux_in_old2(i,time)+1.0E-32))             
	end do
         err_av_time2(i)=err2(i)/time_len2
         err_av2 = err_av2 + err_av_time2(i)
   end do
   err_av2 = abs(err_av2/species_len2)
   err_av_final=(err_av+err_av2)/2.0
!We know species_len=1 in this case. If there were more, all would be taken into account
   write(6,*) "Error = ",err_av_final,"tolerance",err_tolerance

!Write to output file
  if (err_av_final <= err_tolerance) then
     write(143,*) 1,err_av_final !=1 if convergence
     else
     write(143,*) 0,err_av_final !=0 if convergence
  end if


contains
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check
end program PR_conv_SOLPS

  
