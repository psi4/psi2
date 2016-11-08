C     PROGRAM ANHLIN
C**********************************************************
c Moved to PSI area on 2/4/89 by clj.
C**********************************************************
C***LAST UPDATED ON FEBRUARY 11, 1988 BY YUKIO YAMAGUCHI***
C**********************************************************
      SUBROUTINE LTYPLI(ZETA,CXY,F3Q,CLTYP)
C   THE CALCULATION OF L-TYPE DOUBLING CONSTANTS
C   FOR A LINEAR POLYATOMIC MOLECULE
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ZETA(N3N,N3N,3),CXY(3,3,N3N),F3Q(N3N,N3N,N3N)
      DIMENSION CLTYP(N3N,3,3)
      DIMENSION CCC(3,3)
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)
      COMMON/VIB109/CLIMIT,FLIM1,FLIM2
      COMMON/VIB203/IOFF(150),IPRNT
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)
      COMMON/VIB205/IFREQ,NFRQ(150)
      DATA A00,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /
      DATA WLIMIT / 1.0D-02 /
      DATA DLIMIT / 1.0D-02 /
      DATA PI,CL / 3.1415926536D+00 , 2.99792458D+00 /
    1 FORMAT(//,2X,' ***********************************'/
     1          2X,' ***THE L-TYPE DOUBLING CONSTANTS***'/
     2          2X,' ***********************************')
    2 FORMAT(//,2X,' ZETA MATRIX, IABC = ',I5/)
    3 FORMAT(2X,I5,4F10.5)
    4 FORMAT(//,2X,' QE(T) VALU'/)
    5 FORMAT(2X,' QVIB',5X,' IVIB',7X,' FREQENCY',9X,' HARMONIC',
     1      14X,' ANHARMONIC',12X,' TOTAL',14X,' TOTAL'/
     2      24X,' IN CM-1',10X,' IN CM-1',15X,' IN CM-1',
     3      15X,' IN CM-1',12X,' IN MHZ'/)
    6 FORMAT(2X,I5,5X,I5,8X,F10.5,F20.10,3X,F20.10,3X,F20.10,3X,F20.10)
    7 FORMAT(//,2X,' DE = ',F20.10/
     1          2X,' BE = ',F20.10/)
    8 FORMAT(//,2X,' QET = ',F20.10/
     1          2X,' AET = ',F20.10/)
    9 FORMAT(//,2X,' VALU4H = ',F20.10/
     1          2X,' VALU4A = ',F20.10/
     2          2X,' VALU4  = ',F20.10/
     3          2X,' VALU4X = ',F20.10/)
   10 FORMAT(//,2X,' VALU8T = ',F20.10/
     1          2X,' VALU8S = ',F20.10/
     2          2X,' VALU8  = ',F20.10/)
   11 FORMAT(//,2X,' VALU9T = ',F20.10/
     1          2X,' VALU9S = ',F20.10/
     2          2X,' VALU9SP= ',F20.10/
     3          2X,' VALU9  = ',F20.10/)
   12 FORMAT(//,2X,' IVIB   =  ',I5/
     1          2X,' VALU1  = ',F20.10/
     2          2X,' VALU2  = ',F20.10/
     3          2X,' VALU3  = ',F20.10/
     4          2X,' VALU4  = ',F20.10/
     5          2X,' VALU5  = ',F20.10/
     6          2X,' VALU6  = ',F20.10/
     7          2X,' VALU7  = ',F20.10/
     8          2X,' VALU8  = ',F20.10/
     9          2X,' VALU9  = ',F20.10/
     A          2X,' VALU10 = ',F20.10/
     B          2X,' QJT    = ',F20.10/
     C          2X,' VALUH  = ',F20.10/
     D          2X,' VALUA  = ',F20.10/
     E          2X,' VALUT  = ',F20.10/)
   13 FORMAT(//,2X,' QJ(T) VALU'/)
   14 FORMAT(2X,' QVIB',5X,' IVIB',7X,' FREQENCY',9X,' HARMONIC',
     1      14X,' ANHARMONIC',12X,' TOTAL',14X,' TOTAL'/
     2      24X,' IN CM-1',10X,' IN D-07 CM-1',10X,' IN D-07 CM-1',
     3      10X,' IN D-07 CM-1',7X,' IN KHZ'/)
   15 FORMAT(//,2X,' QK(T) VALU'/)
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%  HIGHER-ORDER L-TYPE DOUBLING CONSTANTS FOR LINEAR MOLEUCLES        %
C%  REFERENCE : J.K.G. WATSON, J. MOL. SPECTROSC. 101,83 (1983).       %
C%                                                                     %
C%  NOTE :                                                             %
C%         ZETA(S,T)=-ZETA(T,S)                                        %
C%         IN THIS SUBROUTINE ZETA(S,T) HAS BEEN USED.                 %
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      WRITE(6,1)
C
C   NOTE : THE MOLECULAR AXIS IS CHOSEN AS THE Z AXIS
      IXX=1
      IYY=2
      IZZ=3
      IF(IPRNT.LE.2) GO TO 201
      DO 101 IABC=1,3
      WRITE(6,2) IABC
      CALL MATOUT(ZETA(1,1,IABC),N3N,N3N,NVIB,NVIB,6)
  101 CONTINUE
      IF(IPRNT.GT.2) THEN
        DO 150 IABC=1,NVIB
        WRITE(6,*) ' IABC = ',IABC
        DO 149 I=1,3
        DO 149 J=1,3
        CCC(I,J)=CXY(I,J,IABC)*1.0D+07
  149   CONTINUE
        CALL MATOUT(CCC,3,3,3,3,6)
  150   CONTINUE
      END IF
C
C:::::::::::::::::::::::::::::::
C:::ONLY FOR DEGENERATE MODES:::
C:::::::::::::::::::::::::::::::
  201 CONTINUE
C   CALCULATION OF Q(E)
      CALL ZERO(CLTYP,N3N*9)
C***********************************************************************
C***FOLLOWINGS ARE FOR AN ALTERNATE EXPRESSION                         *
C***THEY ARE TESTED TO GIVE THE SAME RESULTS                           *
C***********************************************************************
C*    DO 103 IVIB=1,NVIB                                               *
C*    IF(NDEG(IVIB).EQ.1) GO TO 103                                    *
C*    FRQI=FREQ(IVIB)                                                  *
C*    FRQI2=FRQI*FRQI                                                  *
C                                                                      *
C   SUM OVER NON-DEGENERATE MODES                                      *
C*    VALU2=A00                                                        *
C*    DO 102 JVIB=1,NVIB                                               *
C*    IF(NDEG(JVIB).NE.1) GO TO 102                                    *
C*    FRQJ=FREQ(JVIB)                                                  *
C*    FRQJ2=FRQJ*FRQJ                                                  *
C*    FACT2=(FRQI2*3.0D+00+FRQJ2)/(FRQI2-FRQJ2)                        *
C*    VALU2=VALU2+(ZETA(IVIB,JVIB,IXX)**2+ZETA(IVIB,JVIB,IYY)**2)*FACT2*
C*102 CONTINUE                                                         *
C*    VALUH=-(VALU2*ROTCM(IXX)*ROTCM(IXX)*TWO)/FRQI                    *
C*    VALUA=A00                                                        *
C*    VALUT=VALUH+VALUA                                                *
C*    CLTYP(IVIB,1,1)=VALUH                                            *
C*    CLTYP(IVIB,2,1)=VALUA                                            *
C*    CLTYP(IVIB,3,1)=VALUT                                            *
C*    IF(IPRNT.LE.2) GO TO 103                                         *
C*    WRITE(6,3) IVIB,FRQI,VALUH,VALUA,VALUT                           *
C*103 CONTINUE                                                         *
C***********************************************************************
C
      BE=ROTCM(IXX)
      BE2=BE*BE
      DO 105 IVIB=1,NVIB
      IF(NDEG(IVIB).EQ.1) GO TO 105
      FRQI=FREQ(IVIB)
      FRQI2=FRQI*FRQI
C
C   SUM OVER NON-DEGENERATE MODES
      VALU2=ONE
      DO 104 JVIB=1,NVIB
      IF(NDEG(JVIB).NE.1) GO TO 104
      FRQJ=FREQ(JVIB)
      FRQJ2=FRQJ*FRQJ
      FACT2=(FRQI2*4.0D+00)/(FRQJ2-FRQI2)
      VALU2=VALU2+(ZETA(IVIB,JVIB,IXX)**2+ZETA(IVIB,JVIB,IYY)**2)*FACT2
  104 CONTINUE
      VALUH=(VALU2*BE2*TWO)/FRQI
      VALUA=A00
      VALUT=VALUH+VALUA
      CLTYP(IVIB,1,1)=VALUH
      CLTYP(IVIB,2,1)=VALUA
      CLTYP(IVIB,3,1)=VALUT
      IF(IPRNT.LE.2) GO TO 105
      WRITE(6,3) IVIB,FRQI,VALUH,VALUA,VALUT
  105 CONTINUE
C
      WRITE(6,4)
      WRITE(6,5)
      DO 106 II=1,NVIB
      IVIB=NFRQ(II)
      IF(NDEG(IVIB).LE.1) GO TO 106
      FRQI=FREQ(IVIB)
      VALUH=CLTYP(IVIB,1,1)
      VALUA=CLTYP(IVIB,2,1)
      VALUT=CLTYP(IVIB,3,1)
      VALUM=VALUT*CL*1.0D+04
      WRITE(6,6) II,IVIB,FRQI,VALUH,VALUA,VALUT,VALUM
  106 CONTINUE
C
C   HIGHER-ORDER L-TYPE DOUBLING CONSTANTS
      VALU0=A00
      DO 107 IVIB=1,NVIB
      IF(NDEG(IVIB).NE.1) GO TO 107
      FRQI=FREQ(IVIB)
      FACT=CXY(IXX,IXX,IVIB)*CXY(IXX,IXX,IVIB)
      VALU0=VALU0+FACT*FRQI
  107 CONTINUE
      DE=VALU0*HALF
      BE=ROTCM(IXX)
      BE2=BE*BE
      BE4=BE2*BE2
      IF(IPRNT.GT.2) THEN
        WRITE(6,7) DE,BE
      END IF
C
      DO 130 IVIB=1,NVIB
      IF(NDEG(IVIB).LE.1) GO TO 130
      FRQI=FREQ(IVIB)
      FRQI2=FRQI*FRQI
      FRQI3=FRQI2*FRQI
      FRQI4=FRQI2*FRQI2
C
C   CALCULATE QE(T) AND AE(T)
C   1ST , 2ND AND 3RD TERMS
      VAL1=A00
      VAL2=A00
      DO 111 JVIB=1,NVIB
      IF(NDEG(JVIB).NE.1) GO TO 111
      FRQJ=FREQ(JVIB)
      FRQJ2=FRQJ*FRQJ
      DNUM=(FRQI2*3.0D+00+FRQJ2)
      FACT=(ZETA(JVIB,IVIB,IXX)**2+ZETA(JVIB,IVIB,IYY)**2)
      DNOM=(FRQJ2-FRQI2)
      VAL1=VAL1+(FACT*DNUM)/DNOM
      VAL2=VAL2+CXY(IXX,IXX,JVIB)*F3Q(JVIB,IVIB,IVIB)
  111 CONTINUE
      QET=VAL1*BE2*TWO/FRQI
      AET=QET*HALF-VAL2*HALF
      IF(IPRNT.GT.2) THEN
        WRITE(6,8) QET,AET
      END IF
C
      VALU1=-(4.0D+00*DE*QET)/BE
      VALU2=(TWO*BE4)/FRQI3
      VALU3=BE2*(AET*TWO-QET)/FRQI2
C
C   4TH , 5TH AND 6TH TERMS
      VAL4=A00
      VAL5=A00
      VAL6=A00
      DO 112 JVIB=1,NVIB
      IF(NDEG(JVIB).NE.1) GO TO 112
      FRQJ=FREQ(JVIB)
      FRQJ2=FRQJ*FRQJ
      FRQJH=DSQRT(FRQJ)
      FACT4=(ZETA(JVIB,IVIB,IXX)**2+ZETA(JVIB,IVIB,IYY)**2)
      DNUM4=FACT4*(FRQJ2+FRQI2)
      DNOM4=(FRQJ2-FRQI2)**2
      VAL4=VAL4+DNUM4/DNOM4
      DNUM5=CXY(IXX,IXX,JVIB)*(ZETA(JVIB,IVIB,IXX)+ZETA(JVIB,IVIB,IYY))
     1     *FRQJ**(1.50D+00)
      DNOM5=(FRQJ2-FRQI2)
      VAL5=VAL5+DNUM5/DNOM5
      VAL6=VAL6+CXY(IXX,IXX,JVIB)*(ZETA(JVIB,IVIB,IXX)
     1                            +ZETA(JVIB,IVIB,IYY))/FRQJH
  112 CONTINUE
      VALU4H=-(8.0D+00*BE2*QET)*VAL4
      VALU4A=(4.0D+00*BE2*VAL2)*VAL4
      VALU4=VALU4H+VALU4A
      VALU5=-(12.0D+00*BE*FRQI)*(VAL5*VAL5)
      VALU6=-(8.0D+00*BE*FRQI)*(VAL5*VAL6)
      IF(IPRNT.GT.2) THEN
        VALU4X=-(4.0D+00*BE2*(TWO*AET+QET))*VAL4
        WRITE(6,9) VALU4H,VALU4A,VALU4,VALU4X
      END IF
C
C   7TH TERM
      VAL7=A00
      DO 115 JVIB=1,NVIB
      IF(NDEG(JVIB).NE.1) GO TO 115
      FRQJ=FREQ(JVIB)
      FRQJ2=FRQJ*FRQJ
      DO 114 KVIB=1,NVIB
      IF(NDEG(KVIB).NE.1) GO TO 114
      FRQK=FREQ(KVIB)
      FRQK2=FRQK*FRQK
      FRQKH=DSQRT(FRQK)
      DO 113 LVIB=1,NVIB
      IF(NDEG(LVIB).NE.1) GO TO 113
      FRQL=FREQ(LVIB)
      FRQL2=FRQL*FRQL
      FRQLH=DSQRT(FRQL)
      DNUM=CXY(IXX,IXX,JVIB)*F3Q(JVIB,KVIB,LVIB)
     1    *(ZETA(KVIB,IVIB,IXX)+ZETA(KVIB,IVIB,IYY))
     2    *(ZETA(LVIB,IVIB,IXX)+ZETA(LVIB,IVIB,IYY))* FRQKH*FRQLH
      DNOM=(FRQK2-FRQI2)*(FRQL2-FRQI2)
      VAL7=VAL7+DNUM/DNOM
  113 CONTINUE
  114 CONTINUE
  115 CONTINUE
      VALU7=-(VAL7*BE2*FRQI)*8.0D+00
C
C   8TH TERM
      VAL8T=A00
      VAL8S=A00
      DO 117 JVIB=1,NVIB
      IF(NDEG(JVIB).NE.2) GO TO 117
      FRQJ=FREQ(JVIB)
      FRQJ2=FRQJ*FRQJ
      IF(DABS(FRQI-FRQJ).LE.DLIMIT) GO TO 117
      VAL8T=VAL8T+(FRQI2*3.0D+00+FRQJ2)/(FRQI2-FRQJ2)
      DO 116 KVIB=1,NVIB
      IF(NDEG(KVIB).NE.1) GO TO 116
      FRQK=FREQ(KVIB)
      FRQK2=FRQK*FRQK
      DNOMS=(FRQK2-FRQI2)
      VAL8S=VAL8S+(ZETA(KVIB,IVIB,IXX)+ZETA(KVIB,IVIB,IYY))
     1           *(ZETA(KVIB,JVIB,IXX)+ZETA(KVIB,JVIB,IYY))/DNOMS
  116 CONTINUE
  117 CONTINUE
      VAL8T=VAL8T*HALF
      VAL8S=VAL8S*HALF
      VALU8=-(BE4*FRQI*VAL8T*VAL8S*VAL8S)*32.0D+00
      IF(IPRNT.GT.2) THEN
        WRITE(6,10) VAL8T,VAL8S,VALU8
      END IF
C
C   9TH TERM
      VAL9T=A00
      VAL9S=A00
      VAL9SP=A00
      DO 120 JVIB=1,NVIB
      IF(NDEG(JVIB).NE.2) GO TO 120
      FRQJ=FREQ(JVIB)
      FRQJ2=FRQJ*FRQJ
      IF(DABS(FRQI-FRQJ).LE.DLIMIT) GO TO 120
      VAL9T=VAL9T+DSQRT(FRQI3*FRQJ)/(FRQI2-FRQJ2)
      DO 118 KVIB=1,NVIB
      IF(NDEG(KVIB).NE.1) GO TO 118
      FRQK=FREQ(KVIB)
      FRQK2=FRQK*FRQK
      DNOMS=FRQK2-FRQI2
      VAL9S=VAL9S+(ZETA(KVIB,IVIB,IXX)+ZETA(KVIB,IVIB,IYY))
     1           *(ZETA(KVIB,JVIB,IXX)+ZETA(KVIB,JVIB,IYY))/DNOMS
  118 CONTINUE
      DO 119 LVIB=1,NVIB
      IF(NDEG(LVIB).NE.1) GO TO 119
      VAL9SP=VAL9SP+CXY(IXX,IXX,LVIB)*F3Q(LVIB,IVIB,JVIB)
  119 CONTINUE
  120 CONTINUE
      VAL9T=VAL9T*HALF
      VAL9S=VAL9S*HALF
      VAL9SP=VAL9SP*HALF
      VALU9=(BE2*VAL9T*VAL9S*VAL9SP)*16.0D+00
      IF(IPRNT.GT.2) THEN
        WRITE(6,11) VAL9T,VAL9S,VAL9SP,VALU9
      END IF
C
      VAL10=A00
      DO 121 JVIB=1,NVIB
      IF(NDEG(JVIB).NE.1) GO TO 121
      FRQJ=FREQ(JVIB)
      FRQJ2=FRQJ*FRQJ
      FRQJ4=FRQJ2*FRQJ2
      DNOM10=FRQI*(FRQJ2-FRQI2)**3
      DNUM10=(ZETA(JVIB,IVIB,IXX)**2+ZETA(JVIB,IVIB,IYY)**2)
     1      *(FRQJ4+FRQJ2*FRQI2*10.0D+00+FRQI4*5.0D+00)
      VAL10=VAL10+DNUM10/DNOM10
  121 CONTINUE
      VALU10=-(VAL10*BE4)*8.0D+00
C
      VALUH=VALU1+VALU2+VALU4H+VALU5+VALU6+VALU8
      VALUA=VALU3+VALU4A+VALU7+VALU9
      QJT=VALU1+VALU2+VALU3+VALU4+VALU5+VALU6+VALU7+VALU8+VALU9
      IF(IPRNT.GT.2) THEN
        VALU1=VALU1*1.0D+08
        VALU2=VALU2*1.0D+08
        VALU3=VALU3*1.0D+08
        VALU4=VALU4*1.0D+08
        VALU5=VALU5*1.0D+08
        VALU6=VALU6*1.0D+08
        VALU7=VALU7*1.0D+08
        VALU8=VALU8*1.0D+08
        VALU9=VALU9*1.0D+08
        VALU10=VALU10*1.0D+08
        VALUH=VALUH*1.0D+08
        VALUA=VALUA*1.0D+08
C
        VALUTX=VALUH+VALUA
        WRITE(6,12) IVIB,VALU1,VALU2,VALU3,VALU4,VALU5,VALU6,VALU7,
     1              VALU8,VALU9,VALU10,QJT,VALUH,VALUA,VALUTX
      END IF
      CLTYP(IVIB,1,2)=VALUH
      CLTYP(IVIB,2,2)=VALUA
      CLTYP(IVIB,3,2)=QJT
      QKT=-QJT+VALU10
      CLTYP(IVIB,1,3)=-VALUH+VALU10
      CLTYP(IVIB,2,3)=-VALUA
      CLTYP(IVIB,3,3)=QKT
C
  130 CONTINUE
C
      WRITE(6,13)
      WRITE(6,14)
      DO 132 II=1,NVIB
      IVIB=NFRQ(II)
      IF(NDEG(IVIB).LE.1) GO TO 132
      FRQI=FREQ(IVIB)
      VALUH=CLTYP(IVIB,1,2)*1.0D+07
      VALUA=CLTYP(IVIB,2,2)*1.0D+07
      VALUT=CLTYP(IVIB,3,2)*1.0D+07
      VALUM=VALUT*CL
      WRITE(6,6) II,IVIB,FRQI,VALUH,VALUA,VALUT,VALUM
  132 CONTINUE
      WRITE(6,15)
      WRITE(6,14)
      DO 133 II=1,NVIB
      IVIB=NFRQ(II)
      IF(NDEG(IVIB).LE.1) GO TO 133
      FRQI=FREQ(IVIB)
      VALUH=CLTYP(IVIB,1,3)*1.0D+07
      VALUA=CLTYP(IVIB,2,3)*1.0D+07
      VALUT=CLTYP(IVIB,3,3)*1.0D+07
      VALUM=VALUT*CL
      WRITE(6,6) II,IVIB,FRQI,VALUH,VALUA,VALUT,VALUM
  133 CONTINUE
C
      RETURN
      END
