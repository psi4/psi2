      subroutine consym(symm,nord)
      character*(*) symm
      integer nord

      if (symm.eq.'C1') then
        nord = 0
      else if (symm.eq.'CS') then
        nord = 0
      else if (symm.eq.'CI') then
        nord = 0
      else if (symm.eq.'C2') then
        symm = 'CN'
        nord = 2
      else if (symm.eq.'C2H') then
        symm = 'CNH'
        nord = 2
      else if (symm.eq.'C2V') then
        symm = 'CNV'
        nord = 2
      else if (symm.eq.'D2') then
        symm = 'DN'
        nord = 2
      else if (symm.eq.'D2H') then
        symm = 'DNH'
        nord = 2
      else
        write(6,*) 'ERROR: bad value for ''SYMMETRY'''
        call fioflu(6)
        call qabort
        endif
      return
      end
