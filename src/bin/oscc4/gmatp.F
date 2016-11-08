      subroutine gmatp(mat,sf11,nlamda,nsym)
      real*8 mat(*), sf11(*)
      integer nlamda(*)
      ii=0
      do 2034 is=1,nsym
         nis=nlamda(is)
         nis2=nis*(nis+1)/2
         do 2033 i=1,nis2
            ii=ii+1
            sf11(i)=mat(ii)
 2033    continue
         write(6,*)
         write(6,*)'irrep   ', is
         call print(sf11,nis2,nis,6)
 2034 continue
      return
      end
