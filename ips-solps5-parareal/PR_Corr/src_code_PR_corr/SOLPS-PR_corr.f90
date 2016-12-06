program main
      use b2mod_types
      use b2mod_version


!*****************************************************************************80
!
! code to apply parareal correction to SOLPS - 23 Oct 2013 - SD
  implicit none

  real ( kind = 8 ) time_init
  real ( kind = 8 ) time_old
  integer time_step
  integer time_step_num
  integer inx,iny,is
  !***Variables read in from or related to SOLPS file - SD
  !From old G run
  character*10 version_plasma_Gold
  character :: lblcp_Gold*120
  real (kind=R8), allocatable :: zamin_Gold(:), zamax_Gold(:), zn_Gold(:), am_Gold(:)
  real (kind=R8), allocatable :: na_Gold(:,:,:), ne_Gold(:,:),ua_Gold(:,:,:),uadia_Gold(:,:,:,:), te_Gold(:,:),ti_Gold(:,:)
  real (kind=R8), allocatable :: po_Gold(:,:),fna_Gold(:,:,:,:), fhe_Gold(:,:,:), fhi_Gold(:,:,:), fch_Gold(:,:,:)
  real (kind=R8), allocatable :: fch_32_Gold(:,:,:), fch_52_Gold(:,:,:), kinrgy_Gold(:,:,:),fch_p_Gold(:,:)
  real (kind=R8) :: time_p_Gold(1)  

  integer:: nx,ny,ns,num_i_Gold(0:2)
  !From new G run
  character*10 version_plasma_Gnew
  character :: lblcp_Gnew*120
  real (kind=R8), allocatable :: zamin_Gnew(:), zamax_Gnew(:), zn_Gnew(:), am_Gnew(:)
  real (kind=R8), allocatable :: na_Gnew(:,:,:), ne_Gnew(:,:),ua_Gnew(:,:,:),uadia_Gnew(:,:,:,:), te_Gnew(:,:),ti_Gnew(:,:)
  real (kind=R8), allocatable :: po_Gnew(:,:),fna_Gnew(:,:,:,:), fhe_Gnew(:,:,:), fhi_Gnew(:,:,:), fch_Gnew(:,:,:)
  real (kind=R8), allocatable :: fch_32_Gnew(:,:,:), fch_52_Gnew(:,:,:), kinrgy_Gnew(:,:,:),fch_p_Gnew(:,:)
  real (kind=R8) :: time_p_Gnew(1)  

  integer:: num_i_Gnew(0:2)

  !From old F run
  character*10 version_plasma_F
  character :: lblcp_F*120
  real (kind=R8), allocatable :: zamin_F(:), zamax_F(:), zn_F(:), am_F(:)
  real (kind=R8), allocatable :: na_F(:,:,:), ne_F(:,:),ua_F(:,:,:),uadia_F(:,:,:,:), te_F(:,:),ti_F(:,:)
  real (kind=R8), allocatable :: po_F(:,:),fna_F(:,:,:,:), fhe_F(:,:,:), fhi_F(:,:,:), fch_F(:,:,:)
  real (kind=R8), allocatable :: fch_32_F(:,:,:), fch_52_F(:,:,:), kinrgy_F(:,:,:),fch_p_F(:,:)
  real (kind=R8) :: time_p_F(1)  

  integer:: num_i_F(0:2)

  !Final corrected value to be written out
  character*10 version_plasma
  character :: lblcp*120
  real (kind=R8), allocatable :: zamin(:), zamax(:), zn(:), am(:)
  real (kind=R8), allocatable :: na(:,:,:), ne(:,:),ua(:,:,:),uadia(:,:,:,:), te(:,:),ti(:,:)
  real (kind=R8), allocatable :: po(:,:),fna(:,:,:,:), fhe(:,:,:), fhi(:,:,:), fch(:,:,:)
  real (kind=R8), allocatable :: fch_32(:,:,:), fch_52(:,:,:), kinrgy(:,:,:),fch_p(:,:)
  real (kind=R8) :: time_p(1)  

  integer:: num_i(0:2)
  

  !************************************
  
  !Variables for PR-SD
  integer*8 :: n_args, ival !number of command line variables for parareal PYTHON framework -SD
  real ( kind = 8 ) :: time_init_tmp, time_final_tmp
  real (kind = 8) :: time1,time2
  character(len=128) c_time_final,c_time_init, c_pe_slice, c_time_step_num !Variable for time setup for slice in PR-SD
      !corr_result: file for result of PR correction
      character(len = 1680) :: file_values_oldG, file_values_newG,file_values_F
      character(len = 1680) :: file_corr_result !all variables for reading in  file names  -SD

!  call timestamp ( )

  call CPU_TIME(time1) !calculate time-SD
!Command line entries for PR-SD
 n_args = command_argument_count() !number of arguments from command line for parareal PYTHON framework -SD
      call get_command_argument(1, file_values_F)
      call get_command_argument(2, file_values_newG)
      call get_command_argument(3, file_values_oldG)
      call get_command_argument(4, file_corr_result)
     !****************open files**********

!      call cfopen (52,'b2fstati','old','formatted')
!      call cfopen (52,'b2fstati','old','un*formatted')
      call cfopen (52,trim(file_values_oldG),'old','formatted')
      call cfopen (62,trim(file_values_newG),'old','formatted')
      call cfopen (72,trim(file_values_F),'old','formatted')
      call cfopen (82,trim(file_corr_result),'new','formatted')
!      open(141,file= trim(file_values_old),status='unknown')
      !*******************************************

     !Read everything in from old F file
     !1st line of b2fstati gives version
     call cfverr (72,version_plasma_F)
      !Get dimensions:
     call cfruin (72, 3, num_i_F, 'nx,ny,ns')


     nx=num_i_F(0)!150
     ny=num_i_F(1)!36
     ns=num_i_F(2)!9
      allocate(zamin_F(0:ns-1), zamax_F(0:ns-1), zn_F(0:ns-1), am_F(0:ns-1))
      allocate(na_F(-1:nx,-1:ny,0:ns-1), ne_F(-1:nx,-1:ny), ua_F(-1:nx,-1:ny,0:ns-1))
      allocate(uadia_F(-1:nx,-1:ny,0:1,0:ns-1),te_F(-1:nx,-1:ny),ti_F(-1:nx,-1:ny))
      allocate (po_F(-1:nx,-1:ny),fna_F(-1:nx,-1:ny,0:1,0:ns-1),fhe_F(-1:nx,-1:ny,0:1),fhi_F(-1:nx,-1:ny,0:1))
      allocate (fch_F(-1:nx,-1:ny,0:1),fch_32_F(-1:nx,-1:ny,0:1),fch_52_F(-1:nx,-1:ny,0:1))
      allocate (kinrgy_F(-1:nx,-1:ny,0:ns-1),fch_p_F(-1:nx,-1:ny))
     !Read label
     call cfruch (72, 120, lblcp_F, 'label')

     !Read data
     call cfrure (72, ns, zamin_F, 'zamin')
     call cfrure (72, ns, zamax_F, 'zamax')
     call cfrure (72, ns, zn_F, 'zn')
     call cfrure (72, ns, am_F, 'am')
     call cfrure (72, size(na_F), na_F, 'na')
     call cfrure (72, size(ne_F), ne_F, 'ne')
     call cfrure (72, size(ua_F), ua_F, 'ua')
     call cfrure (72, size(uadia_F), uadia_F, 'uadia')
     call cfrure (72, size(te_F), te_F, 'te')
     call cfrure (72, size(ti_F), ti_F, 'ti')
     call cfrure (72, size(po_F), po_F, 'po')
     call cfrure (72, size(fna_F), fna_F, 'fna')
     call cfrure (72, size(fhe_F), fhe_F, 'fhe')
     call cfrure (72, size(fhi_F), fhi_F, 'fhi')
     call cfrure (72, size(fch_F), fch_F, 'fch')
     call cfrure (72, size(fch_32_F), fch_32_F, 'fch_32')
     call cfrure (72, size(fch_52_F), fch_52_F, 'fch_52')
     call cfrure (72, size(kinrgy_F), kinrgy_F, 'kinrgy')
     call cfrure (72, 1, time_p_F, 'time')
     call cfrure (72, size(fch_p_F), fch_p_F, 'fch_p')
     
     

     !Read everything in from new G file
     !1st line of b2fstati gives version
     call cfverr (62,version_plasma_Gnew)
      !Get dimensions:
     call cfruin (62, 3, num_i_Gnew, 'nx,ny,ns')


     nx=num_i_Gnew(0)!150
     ny=num_i_Gnew(1)!36
     ns=num_i_Gnew(2)!9
      allocate(zamin_Gnew(0:ns-1), zamax_Gnew(0:ns-1), zn_Gnew(0:ns-1), am_Gnew(0:ns-1))
      allocate(na_Gnew(-1:nx,-1:ny,0:ns-1), ne_Gnew(-1:nx,-1:ny), ua_Gnew(-1:nx,-1:ny,0:ns-1))
      allocate(uadia_Gnew(-1:nx,-1:ny,0:1,0:ns-1),te_Gnew(-1:nx,-1:ny),ti_Gnew(-1:nx,-1:ny))
      allocate (po_Gnew(-1:nx,-1:ny),fna_Gnew(-1:nx,-1:ny,0:1,0:ns-1),fhe_Gnew(-1:nx,-1:ny,0:1),fhi_Gnew(-1:nx,-1:ny,0:1))
      allocate (fch_Gnew(-1:nx,-1:ny,0:1),fch_32_Gnew(-1:nx,-1:ny,0:1),fch_52_Gnew(-1:nx,-1:ny,0:1))
      allocate (kinrgy_Gnew(-1:nx,-1:ny,0:ns-1),fch_p_Gnew(-1:nx,-1:ny))
     !Read label
     call cfruch (62, 120, lblcp_Gnew, 'label')

     !Read data
     call cfrure (62, ns, zamin_Gnew, 'zamin')
     call cfrure (62, ns, zamax_Gnew, 'zamax')
     call cfrure (62, ns, zn_Gnew, 'zn')
     call cfrure (62, ns, am_Gnew, 'am')
     call cfrure (62, size(na_Gnew), na_Gnew, 'na')
     call cfrure (62, size(ne_Gnew), ne_Gnew, 'ne')
     call cfrure (62, size(ua_Gnew), ua_Gnew, 'ua')
     call cfrure (62, size(uadia_Gnew), uadia_Gnew, 'uadia')
     call cfrure (62, size(te_Gnew), te_Gnew, 'te')
     call cfrure (62, size(ti_Gnew), ti_Gnew, 'ti')
     call cfrure (62, size(po_Gnew), po_Gnew, 'po')
     call cfrure (62, size(fna_Gnew), fna_Gnew, 'fna')
     call cfrure (62, size(fhe_Gnew), fhe_Gnew, 'fhe')
     call cfrure (62, size(fhi_Gnew), fhi_Gnew, 'fhi')
     call cfrure (62, size(fch_Gnew), fch_Gnew, 'fch')
     call cfrure (62, size(fch_32_Gnew), fch_32_Gnew, 'fch_32')
     call cfrure (62, size(fch_52_Gnew), fch_52_Gnew, 'fch_52')
     call cfrure (62, size(kinrgy_Gnew), kinrgy_Gnew, 'kinrgy')
     call cfrure (62, 1, time_p_Gnew, 'time')
     call cfrure (62, size(fch_p_Gnew), fch_p_Gnew, 'fch_p')




     !Read everything in from old G file
     !1st line of b2fstati gives version
     call cfverr (52,version_plasma_Gold)
      !Get dimensions:
     call cfruin (52, 3, num_i_Gold, 'nx,ny,ns')


     nx=num_i_Gold(0)!150
     ny=num_i_Gold(1)!36
     ns=num_i_Gold(2)!9
      allocate(zamin_Gold(0:ns-1), zamax_Gold(0:ns-1), zn_Gold(0:ns-1), am_Gold(0:ns-1))
      allocate(na_Gold(-1:nx,-1:ny,0:ns-1), ne_Gold(-1:nx,-1:ny), ua_Gold(-1:nx,-1:ny,0:ns-1))
      allocate(uadia_Gold(-1:nx,-1:ny,0:1,0:ns-1),te_Gold(-1:nx,-1:ny),ti_Gold(-1:nx,-1:ny))
      allocate (po_Gold(-1:nx,-1:ny),fna_Gold(-1:nx,-1:ny,0:1,0:ns-1),fhe_Gold(-1:nx,-1:ny,0:1),fhi_Gold(-1:nx,-1:ny,0:1))
      allocate (fch_Gold(-1:nx,-1:ny,0:1),fch_32_Gold(-1:nx,-1:ny,0:1),fch_52_Gold(-1:nx,-1:ny,0:1))
      allocate (kinrgy_Gold(-1:nx,-1:ny,0:ns-1),fch_p_Gold(-1:nx,-1:ny))
     !Read label
     call cfruch (52, 120, lblcp_Gold, 'label')

     !Read data
     call cfrure (52, ns, zamin_Gold, 'zamin')
     call cfrure (52, ns, zamax_Gold, 'zamax')
     call cfrure (52, ns, zn_Gold, 'zn')
     call cfrure (52, ns, am_Gold, 'am')
     call cfrure (52, size(na_Gold), na_Gold, 'na')
     call cfrure (52, size(ne_Gold), ne_Gold, 'ne')
     call cfrure (52, size(ua_Gold), ua_Gold, 'ua')
     call cfrure (52, size(uadia_Gold), uadia_Gold, 'uadia')
     call cfrure (52, size(te_Gold), te_Gold, 'te')
     call cfrure (52, size(ti_Gold), ti_Gold, 'ti')
     call cfrure (52, size(po_Gold), po_Gold, 'po')
     call cfrure (52, size(fna_Gold), fna_Gold, 'fna')
     call cfrure (52, size(fhe_Gold), fhe_Gold, 'fhe')
     call cfrure (52, size(fhi_Gold), fhi_Gold, 'fhi')
     call cfrure (52, size(fch_Gold), fch_Gold, 'fch')
     call cfrure (52, size(fch_32_Gold), fch_32_Gold, 'fch_32')
     call cfrure (52, size(fch_52_Gold), fch_52_Gold, 'fch_52')
     call cfrure (52, size(kinrgy_Gold), kinrgy_Gold, 'kinrgy')
     call cfrure (52, 1, time_p_Gold, 'time')
     call cfrure (52, size(fch_p_Gold), fch_p_Gold, 'fch_p')


  write(6,*) "label",lblcp_Gnew
  write (6,*) "ver:",version_plasma_F,"na is",na_F(1,1,0),"ne is",ne_F(1,1)
  write(6,*) "nx,ny,ns are",nx,ny,ns

     !Perform parareal correction. 
     !First prepare output  variables.
  
      allocate(zamin(0:ns-1), zamax(0:ns-1), zn(0:ns-1), am(0:ns-1))
      allocate(na(-1:nx,-1:ny,0:ns-1), ne(-1:nx,-1:ny), ua(-1:nx,-1:ny,0:ns-1))
      allocate(uadia(-1:nx,-1:ny,0:1,0:ns-1),te(-1:nx,-1:ny),ti(-1:nx,-1:ny))
      allocate (po(-1:nx,-1:ny),fna(-1:nx,-1:ny,0:1,0:ns-1),fhe(-1:nx,-1:ny,0:1),fhi(-1:nx,-1:ny,0:1))
      allocate (fch(-1:nx,-1:ny,0:1),fch_32(-1:nx,-1:ny,0:1),fch_52(-1:nx,-1:ny,0:1))
      allocate (kinrgy(-1:nx,-1:ny,0:ns-1),fch_p(-1:nx,-1:ny))

      !First equate all output with F run before applying PR correction:
      zamin=zamin_F
      zamax=zamax_F
      zn=zn_F
      am=am_F
      na=na_F      
      ne=ne_F
      ua=ua_F
      uadia=uadia_F
      te=te_F
      ti=ti_F
      po=po_F
      fna=fna_F
      fhe=fhe_F
      fhi=fhi_F
      fch=fch_F
      fch_32=fch_32_F
      fch_52=fch_52_F
      kinrgy=kinrgy_F
      time_p=time_p_F
      fch_p=fch_p_F

      !Now, apply parareal correction:
      !na(is), ua(is), te,  ti, po     
      na = na_F + na_Gnew - na_Gold
      ua = ua_F + ua_Gnew - ua_Gold
      te = te_F + te_Gnew - te_Gold
      ti = ti_F + ti_Gnew - ti_Gold
      po = po_F + po_Gnew - po_Gold
      ne = ne_F + ne_Gnew - ne_Gold

     !Check condition to ensure na, te, ti are .ge. 0, as suggested in preprocessing/b2aidr.F 

     do is=0,ns-1
        do inx=-1,nx
               do iny=-1,ny
                   !No correction to neutrals if following if loop is NOT commented out
                   !if (zamin(is).le.0) then
                           !write(6,*) "as zamin is 0 at is=",is,"no corr"
                           !na(inx,iny,is) = na_F(inx,iny,is)
                           !ua(inx,iny,is) = ua_F(inx,iny,is)
                   !end if
                   if (na(inx,iny,is) .lt. 0) then
                          write(6,*) "na left at na_F at ix,iy,is:",inx,iny,is,na(inx,iny,is), na_F(inx,iny,is),na_Gnew(inx,iny,is),na_Gold(inx,iny,is)
                          na(inx,iny,is) = na_F(inx,iny,is)
                   end if
                end do
        end do
     end do
     do inx=-1,nx
         do iny=-1,ny
             if (te(inx,iny) .lt. 0) then
   write(6,*) "te left at te_F at ix,iy:",inx,iny,te(inx,iny),te_F(inx,iny),te_Gnew(inx,iny), te_Gold(inx,iny)
                  te(inx,iny) = te_F(inx,iny)
             end if
             if (ti(inx,iny) .lt. 0) then
    write(6,*) "ti left at ti_F at ix,iy:",inx,iny,ti(inx,iny),ti_F(inx,iny),ti_Gnew(inx,iny), ti_Gold(inx,iny)
                  ti(inx,iny) = ti_F(inx,iny)
             end if
             if (ne(inx,iny) .lt. 0) then
    write(6,*) "ne left at ne_F at ix,iy:",inx,iny,ne(inx,iny),ne_F(inx,iny),ne_Gnew(inx,iny), ne_Gold(inx,iny)
                  ne(inx,iny) = ne_F(inx,iny)
             end if
         end do
     end do

     
      !Write to file
      call cfverw (82, version_plasma_F)
      !Put dimensions:
      call cfwuin (82, 3, num_i_F, 'nx,ny,ns')
      !Put label same as F run
      call cfwuch (82, 120, lblcp_F, 'label')

      !Write remaining data
     call cfwure (82, ns, zamin, 'zamin')
     call cfwure (82, ns, zamax, 'zamax')
     call cfwure (82, ns, zn, 'zn')
     call cfwure (82, ns, am, 'am')
     call cfwure (82, size(na), na, 'na')
     call cfwure (82, size(ne), ne, 'ne')
     call cfwure (82, size(ua), ua, 'ua')
     call cfwure (82, size(uadia), uadia, 'uadia')
     call cfwure (82, size(te), te, 'te')
     call cfwure (82, size(ti), ti, 'ti')
     call cfwure (82, size(po), po, 'po')
     call cfwure (82, size(fna), fna, 'fna')
     call cfwure (82, size(fhe), fhe, 'fhe')
     call cfwure (82, size(fhi), fhi, 'fhi')
     call cfwure (82, size(fch), fch, 'fch')
     call cfwure (82, size(fch_32), fch_32, 'fch_32')
     call cfwure (82, size(fch_52), fch_52, 'fch_52')
     call cfwure (82, size(kinrgy), kinrgy, 'kinrgy')
     call cfwure (82, 1, time_p, 'time')
     call cfwure (82, size(fch_p), fch_p, 'fch_p')
            


  write ( *, '(a)' ) '  Normal end of execution.'

!  call timestamp ( )

  call CPU_TIME(time2)
  write(6,*) "Time for run is",time2-time1
  !stop
end
