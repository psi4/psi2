      SUBROUTINE BERROR 
      COMMON /INOUT/  IOUT,ICHECK,ITAP30
C                                      
      entry msread
      write(6,*) 'mips f77 2.0 requires statement here'
      entry n21g
      write(6,*) 'mips f77 2.0 requires statement here'
      entry n31g
      write(6,*) 'mips f77 2.0 requires statement here'
      entry dijgen
      write(6,*) 'mips f77 2.0 requires statement here'
C                      
      WRITE(IOUT,1)                   
    1 FORMAT(//,' ILLEGAL SUBROUTINE CALLED:BASIS ERROR?',//)
      CALL mabort
C                                                         
      END                                                  
