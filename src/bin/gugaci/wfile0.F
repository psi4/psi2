      SUBROUTINE WFILE0(FILE,NSECT,IA,LENGTH)
C
C10-13-90
C
      IMPLICIT INTEGER (A-Z)
C
      IF (SEC2I(NSECT).GT.LENGTH) THEN
      N = LENGTH/SEC2I(1)
      NPASS = NSECT/N
C10-13-90   ADD 1 LINE BELOW --- NEXSEC is now set if "IF" below False
      NEXSEC = N
      IF(NPASS*N.NE.NSECT) THEN
C10-13-90  NPASS = NPASS + 1
C10-13-90  NEXSEC = NPASS*N - NSECT   this gives wrong remainder
C10-13-90  switch order of lines and switch order in subtraction
      NEXSEC =  NSECT - NPASS*N
      NPASS = NPASS + 1
      END IF
      NLEN = SEC2I(N)
C
      DO 1 I = 1,NPASS
      IF(I.EQ.NPASS) NLEN = SEC2I(NEXSEC)
      CALL SWRIT(FILE,IA,NLEN)
   1  CONTINUE
C
      ELSE
      CALL SWRIT(FILE,IA,SEC2I(NSECT))
      END IF
C
C
      RETURN
      END
