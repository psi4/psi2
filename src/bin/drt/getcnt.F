clj070289 0ALP1 is now really 0 alpha's of sym 1, etc.
      INTEGER FUNCTION GETCNT()
C
      IMPLICIT INTEGER (A-Z)
      CHARACTER*1 BLANK,DIGIT,JUNK,CHAR,NEXTC
      logical found
C
      COMMON /CONST/  BLANK, DIGIT(10)
C
      GETCNT=0
      found = .false.
    1 JUNK=NEXTC(CHAR)
      DO 2 N=1,10
         IF (CHAR.EQ.DIGIT(N)) GO TO 3
    2 CONTINUE
      CALL BKSPAC
      IF (.not.found) GETCNT=1
      RETURN
C
    3 GETCNT=GETCNT*10 + (N-1)
      found = .true.
      GO TO 1
      END
