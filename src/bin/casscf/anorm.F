      SUBROUTINE ANORM(CNORM,C,NNC,NC)                                  SMM21560
      IMPLICIT REAL*8 (A-H,O-Z)                                         SMM21570
      DIMENSION C(NNC)                                                  SMM21580
      DATA A00 / 0.0D+00 /                                              SMM21590
C                                                                       SMM21600
      CNORM=A00                                                         SMM21610
      DO 101 I=1,NC                                                     SMM21620
      CNORM=CNORM+C(I)*C(I)                                             SMM21630
  101 CONTINUE                                                          SMM21640
      CNORM=DSQRT(CNORM)                                                SMM21650
C                                                                       SMM21660
      RETURN                                                            SMM21670
      END                                                               SMM21680
