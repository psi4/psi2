      SUBROUTINE NORMLZ(V,NDIM,N,ISGN)                                  SMM23930
      IMPLICIT REAL*8 (A-H,O-Z)                                         SMM23940
      DIMENSION V(NDIM)                                                 SMM23950
    1 FORMAT(//,2X,'***ERROR IN NORMALIZATION***'/)                     SMM23960
      DATA RLIM / 1.0D-38 /                                             SMM23970
      DATA A00 / 0.0D+00 /                                              SMM23980
C                                                                       SMM23990
      SQL=A00                                                           SMM24000
      DO 101 I=1,N                                                      SMM24010
  101 SQL=SQL+V(I)*V(I)                                                 SMM24020
      R=DSQRT(SQL)                                                      SMM24030
      IF(R.LE.RLIM) GO TO 202                                           SMM24040
C                                                                       SMM24050
      IF(ISGN.LE.0) GO TO 201                                           SMM24060
      IF(V(ISGN).LT.A00) R=-R                                           SMM24070
C                                                                       SMM24080
  201 DO 102 I=1,N                                                      SMM24090
  102 V(I)=V(I)/R                                                       SMM24100
      GO TO 203                                                         SMM24110
  202 WRITE(6,1)                                                        SMM24120
      STOP                                                              SMM24130
C                                                                       SMM24140
  203 RETURN                                                            SMM24150
      END                                                               SMM24160
