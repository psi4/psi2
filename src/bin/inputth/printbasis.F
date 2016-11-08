c Utilities to print out information the basis set.
      subroutine intpb(c,maxcor)
#ifdef sun
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
#include <error.h>
      integer maxcor
      integer c(maxcor)

#include <basis.h>

      call intpb1(c(pexpon),c(pconco),c(ppshc),c(pnprim),
     &            c(ppprim),c(ptypes),c(pnshc),
     &            nshell,nprimt,nucent)


      return
      end
c---------------------------------------------------------------------
      subroutine intpb1(expon,conco,pshc,nprim,pprim,types,nshc,
     &                        nshell,nprimt,nucent)
#ifdef sun
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
      integer nshell,nprimt,nucent
      real*8 expon(nprimt),conco(nprimt)
      integer pshc(nucent),nshc(nucent),nprim(nshell),pprim(nshell)
      integer types(nshell)

      character*3 c1
      character*1 c2,c3,c4
      character*10 ctype
      integer i,j,k,ish,ipr

      integer chrlen

c     write(6,*) 'The number of unique centers = ',nucent
c     write(6,*) 'The total number of unique shells = ',nshell
c     write(6,*) 'The total number of unique primitives = ',nprimt

      ish = 0
      ipr = 0

      do 101 i = 1, nucent
        write(6,*)
        write(6,'(A,I6,A)') ' Basis set on unique center ',i,':'
c       write(6,*) 'Number of shells on this center = ',nshc(i)
        do 201 j = 1, nshc(i)
          ish = ish + 1
c         write(6,*) 'NPRIM(',ish,') =',nprim(ish)
          call inttyp('I_TO_C',ctype,types(ish))
c         write(6,1) '  Shell number = ', j,
c    &               ', type = ', ctype(1:chrlen(ctype)),
c    &               ', number of primitives = ', nprim(ish)
c  1      format(a,i5,a,a,a,i5)
          do 301 k = 1, nprim(ish)
            ipr = ipr + 1

            if ((j.eq.1).and.(k.eq.1)) then
              c1 = '  ('
              c2 = '('
            else if (k.eq.1) then
              c1 = '   '
              c2 = '('
            else
              c1 = '   '
              c2 = ' '
              endif
            if ((j.eq.nshc(i)).and.(k.eq.nprim(ish))) then
              c3 = ')'
              c4 = ')'
            else if (k.eq.nprim(ish)) then
              c3 = ')'
              c4 = ' '
            else
              c3 = ' '
              c4 = ' '
              endif
            
            write(6,3) c1,c2,ctype(1:max0(chrlen(ctype),4)),
     &                 '    ( ',expon(ipr),conco(ipr),' )',c3,c4
   3        format(a,a,a,a,f17.8,1x,f17.8,a,a,a)
c           write(6,2) '    Primitive number = ', k,
c    &                 ', exponent = ',expon(ipr),
c    &                 ', con. coef. = ',conco(ipr)
c  2        format(a,i5,a,f20.10,a,f20.10)
            ctype = ' '
 301        continue
 201      continue
 101    continue

      return
      end
