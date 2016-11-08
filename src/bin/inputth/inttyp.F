c This converts a character type to an integer type or the reverse
c conv is 'I_TO_C' or 'C_TO_I' to specify which conversion is to be
c done.
      subroutine inttyp(conv,ctype,itype)
#ifdef UNDEFINED
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
#include <error.h>
      character*(*) conv,ctype
      integer itype

      integer i,mxtype
      parameter (mxtype=11)
      character*3 ct(mxtype)
      data ct/'S','P','D','F','G','D6','D5','F10','F7','G15','G9'/

c     write(6,*) 'Entered inttyp:',conv,'  ',ctype,itype

      if (conv.eq.'I_TO_C') then
        if ((itype.le.0).or.(itype.gt.mxtype)) then
          write(6,*) 'ERROR: bad itype given to inttyp'
          write(6,*) ' itype = ',itype
          call qabort
          endif
        ctype = ct(itype)
        if (ctype .ne. ct(itype)) then
          write(6,*) 'ERROR: truncation of ctype in inttyp'
          call qabort
          endif
c       write(6,*) '      Leaving inttyp:','  ',ctype,itype
        return
      else if (conv.eq.'C_TO_I') then
        do 101 i = 1, mxtype
          if (ctype.eq.ct(i)) then
            itype = i
c           write(6,*) '      Leaving inttyp:','  ',ctype,itype
            return
            endif
  101     continue
        write(6,*) 'ERROR: bad ctype given to inttyp'
        write(6,*) ' ctype = ',ctype
        write(6,*) ' possible values of ctype are:'
        write(6,'(20(1X,A))') ct
        call qabort
      else
        write(6,*) 'ERROR: bad conv given to inttyp'
        call qabort
        endif 

      return
      end
