      SUBROUTINE SKIPBL
C
      IMPLICIT INTEGER (A-Z)
C
      COMMON /CONST/  BLANK, DIGIT(10)
      CHARACTER*1 BLANK,DIGIT,CHAR,GETC
      character*1 tmp
C
clj080889  Modified to work around MIPS 2.0 f77 compiler bugs
    1 continue
        tmp = getc(char)
        IF (tmp.EQ.BLANK) GO TO 1
      CALL BKSPAC
      RETURN
      END
