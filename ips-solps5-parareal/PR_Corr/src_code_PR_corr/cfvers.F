      subroutine cfvers
      implicit none
      character*10 version_in, newversion
      integer nget, nput
      character chf*12,label*7
      logical streql
      external streql

      entry cfverr(nget, newversion)
      inquire (nget,form=chf,err=10)
      if (streql(chf,'formatted')) then
        read(nget,'(a,a)',err=10) label,version_in
      else
        read(nget,err=10) label,version_in
      endif
      if(label.eq.'VERSION') then
	newversion=version_in
        write(*,*) 'VERSION set to: ',newversion
        return
      endif
   10 continue
      rewind(nget)
      write(*,*) 'VERSION defaulted to : ',newversion
      return

      entry cfverw(nput, newversion)
      inquire (nput,form=chf,err=20)
      if (streql(chf,'formatted')) then
        write(nput,'(a,a)') 'VERSION',newversion
      else
        write(nput) 'VERSION',newversion
      endif
      write(*,*) 'VERSION: ',newversion
   20 continue
      return

      end
