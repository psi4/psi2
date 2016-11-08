

      SUBROUTINE MABORT

      write(6,*) 'mabort: ERROR, aborting'
      write(*,*) 'mabort: ERROR, aborting'

      write(0,*) 'mabort: ERROR, aborting'

      call exit(3)

      stop 876
C
      END
