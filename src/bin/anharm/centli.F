      SUBROUTINE CENTLI(BXY,CXY,F3Q,EE)
C   THE CALCULATION OF CENTRIFUGAL DISTORTION CONSTANTS
C   FOR A LINEAR POLYATOMIC MOLECULE
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION BXY(3,3,N3N),CXY(3,3,N3N)
      DIMENSION F3Q(N3N,N3N,N3N),EE(N3N,N3N)
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)
      COMMON/VIB203/IOFF(150),IPRNT
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)
      COMMON/VIB205/IFREQ,NFRQ(150)
      COMMON/VIB207/TABCD(3,3,3,3),PABC(3,3,3)
      COMMON/VIB208/TAABB(3,3),TABAB(3,3)
      DATA WLIMIT / 1.0D+00 /
      DATA CL / 2.99792458D+00 /
      DATA A00,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /
    1 FORMAT(//,2X,' ******************************************'/
     1          2X,' ***THE CENTRIFUGAL DISTORTION CONSTANTS***'/
     2          2X,' ******************************************')
    2 FORMAT(//,2X,' F3(QI*QJ*QK) MATRIX (IN CM-1), IVIB = ',I5/)
    3 FORMAT(//,2X,' TAABB MATRIX FROM BXY MATRICES (IN CM-1)'/)
    4 FORMAT(//,2X,' TABAB MATRIX FROM BXY MATRICES (IN CM-1)'/)
    5 FORMAT(//,2X,' FOR A LINEAR POLYATOMIC MOLECULE'/
     1          2X,' DJ = ',F20.10,' D-10 CM-1',
     2         10X,' DJ = ',F20.10,' MHZ'/)
    6 FORMAT(2X,' HJ = ',F20.10,' D-10 CM-1',
     1         10X,' HJ = ',F20.10,' HZ'/)
C
      WRITE(6,1)
C   ELEMENTS OF F3Q ARE IN CM-1
      IF(IPRNT.LE.2) GO TO 201
      DO 102 III=1,NVIB
      IVIB=NFRQ(III)
      WRITE(6,2) IVIB
      DO 101 I=1,NVIB
      II=NFRQ(I)
      DO 101 J=1,NVIB
      JJ=NFRQ(J)
      EE(I,J)=F3Q(II,JJ,IVIB)
  101 CONTINUE
      CALL ANHOUT(EE,N3N,N3N,NVIB,NVIB,6)
  102 CONTINUE
C
  201 CONTINUE
C   FOR A LINEAR POLYATOMIC MOLECULES
C   QUARTIC CENTRIFUGAL CONSTANTS
C
C***********************************************************************
C***FOLLOWING CODES ARE FOR ALTERNATE EXPRESSIONS                      *
C***THEY ARE TESTED TO GIVE THE SAME RESULTS                           *
C***********************************************************************
C*    DO 103 I=1,3                                                     *
C*    DO 103 J=1,3                                                     *
C*    TAABB(I,J)=TABCD(I,I,J,J)                                        *
C*    TABAB(I,J)=TABCD(I,J,I,J)                                        *
C*103 CONTINUE                                                         *
C*    IF(IPRNT.LE.2) GO TO 202                                         *
C*    WRITE(6,3)                                                       *
C*    CALL MATOUT(TAABB,3,3,3,3,6)                                     *
C*    WRITE(6,4)                                                       *
C*    CALL MATOUT(TABAB,3,3,3,3,6)                                     *
C                                                                      *
C*202 CONTINUE                                                         *
C*    IXX=1                                                            *
C*    DJ=-TAABB(IXX,IXX)*0.25D+00                                      *
C*    DJCM=DJ*1.0D+10                                                  *
C*    DJMH=DJ*CL*1.0D+4                                                *
C*    WRITE(6,5) DJCM,DJMH                                             *
C***********************************************************************
C*    DJ=A00                                                           *
C*    DO 104 IVIB=1,NVIB                                               *
C*    IF(NDEG(IVIB).NE.1) GO TO 104                                    *
C*    FRQI=FREQ(IVIB)                                                  *
C*    IF(DABS(FRQI).LE.WLIMIT) GO TO 104                               *
C*    FACT=BXY(IXX,IXX,IVIB)*BXY(IXX,IXX,IVIB)/FRQI                    *
C*    DJ=DJ+FACT/TWO                                                   *
C*104 CONTINUE                                                         *
C*    DJCM=DJ*1.0D+10                                                  *
C*    DJMH=DJ*CL*1.0D+4                                                *
C*    WRITE(6,5) DJCM,DJMH                                             *
C***********************************************************************
C
      IXX=1
      DJ=A00
      DO 105 IVIB=1,NVIB
      IF(NDEG(IVIB).NE.1) GO TO 105
      FRQI=FREQ(IVIB)
      IF(DABS(FRQI).LE.WLIMIT) GO TO 105
      FACT=CXY(IXX,IXX,IVIB)*CXY(IXX,IXX,IVIB)*FRQI
      DJ=DJ+FACT/TWO
  105 CONTINUE
      DJCM=DJ*1.0D+10
      DJMH=DJ*CL*1.0D+4
      WRITE(6,5) DJCM,DJMH
C
C   SEXTIC CENTRIFUGAL DISTORTION CONSTANTS
C***********************************************************************
C***FOLLWOING CODES ARE FOR AN ALTERNATE EXPRESSION                    *
C***THEY ARE TESTED TO GIVE THE SAME RESUTLS                           *
C***********************************************************************
C*    IXX=1                                                            *
C*    VALU1=(DJ*DJ*4.0D+00)/ROTCM(IXX)                                 *
C*    VALU2=A00                                                        *
C*    VALU3=A00                                                        *
C*    ROTSQ2=ROTCM(IXX)*ROTCM(IXX)*TWO                                 *
C*    DO 110 IVIB=1,NVIB                                               *
C*    IF(NDEG(IVIB).NE.1) GO TO 110                                    *
C*    FRQI=FREQ(IVIB)                                                  *
C*    FRQI3=FRQI**3                                                    *
C*    BI=BXY(IXX,IXX,IVIB)                                             *
C*    VALU2=VALU2-ROTSQ2*BXY(IXX,IXX,IVIB)*BXY(IXX,IXX,IVIB)/FRQI3     *
C*    DO 109 JVIB=1,NVIB                                               *
C*    IF(NDEG(JVIB).NE.1) GO TO 109                                    *
C*    FRQJ=FREQ(JVIB)                                                  *
C*    BJ=BXY(IXX,IXX,JVIB)                                             *
C*    DO 108 KVIB=1,NVIB                                               *
C*    IF(NDEG(KVIB).NE.1) GO TO 108                                    *
C*    FRQK=FREQ(KVIB)                                                  *
C*    BK=BXY(IXX,IXX,KVIB)                                             *
C*    FACT=BI*BJ*BK                                                    *
C*    VALU3=VALU3-FACT*F3Q(IVIB,JVIB,KVIB)/(FRQI*FRQJ*FRQK)            *
C*108 CONTINUE                                                         *
C*109 CONTINUE                                                         *
C*110 CONTINUE                                                         *
C*    HJ=VALU1+VALU2+VALU3/6.0D+00                                     *
C*    HJCM=HJ*1.0D+10                                                  *
C*    HJHZ=HJ*CL*1.0D+10                                               *
C*    WRITE(6,6) HJCM,HJHZ                                             *
C***********************************************************************
C
      IXX=1
      VALU1=(DJ*DJ*4.0D+00)/ROTCM(IXX)
      VALU2=A00
      VALU3=A00
      ROTSQ2=ROTCM(IXX)*ROTCM(IXX)*TWO
      DO 115 IVIB=1,NVIB
      IF(NDEG(IVIB).NE.1) GO TO 115
      FRQI=FREQ(IVIB)
      CI=CXY(IXX,IXX,IVIB)
      VALU2=VALU2-ROTSQ2*CXY(IXX,IXX,IVIB)*CXY(IXX,IXX,IVIB)/FRQI
      DO 114 JVIB=1,NVIB
      IF(NDEG(JVIB).NE.1) GO TO 114
      CJ=CXY(IXX,IXX,JVIB)
      DO 113 KVIB=1,NVIB
      IF(NDEG(KVIB).NE.1) GO TO 113
      CK=CXY(IXX,IXX,KVIB)
      FACT=CI*CJ*CK
      VALU3=VALU3+FACT*F3Q(IVIB,JVIB,KVIB)
  113 CONTINUE
  114 CONTINUE
  115 CONTINUE
      HJ=VALU1+VALU2+VALU3/6.0D+00
C
      HJCM=HJ*1.0D+10
      HJHZ=HJ*CL*1.0D+10
      WRITE(6,6) HJCM,HJHZ
      RETURN
      END
