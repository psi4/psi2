      SUBROUTINE ASYTOP(ZETA,F3Q,F4Q,XIJ)
C   THE ANHARMONICITY CONTSTANTS FOR AN ASYMMETRIC TOP
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ZETA(N3N,N3N,3)
      DIMENSION F3Q(N3N,N3N,N3N),F4Q(N3N,N3N,N3N,N3N)
      DIMENSION XIJ(N3N,N3N)
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)
      COMMON/VIB109/CLIMIT,FLIM1,FLIM2
      COMMON/VIB203/IOFF(150),IPRNT
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)
      COMMON/VIB205/IFREQ,NFRQ(150)
      DATA A00,ONE / 0.0D+00 , 1.0D+00 /
      DATA WLIMIT / 1.0D+00 /
      DATA PLIMIT / 1.0D-05 /
    1 FORMAT(//,2X,' ********************************************'
     1        /,2X,' ***THE VIBRATIONAL ANHARMONICTY CONSTANTS***'
     2        /,2X,' ********************************************')
    2 FORMAT(//,2X,' ANHARMONIC CONSTANTS IN CM-1'/
     1          2X,' IVIB',3X,' JVIB',8X,' FREQ(I)',5X,' FREQ(J)',
     2          7X,' VALU1',7X,' VALU2',7X,' VALU3',7X,' VALU4',
     3          7X,' XIJ'/)
    3 FORMAT(2X,I5,3X,I5,3X,7(F13.3))
    4 FORMAT(//,2X,' ANHARMONIC CONSTANTS IN CM-1'/
     1          4X,' QI',5X,' QJ',8X,' FREQ(I)',5X,' FREQ(J)',
     2          7X,' XIJ'/)
    5 FORMAT(2X,I5,3X,I5,3X,3(F13.3))
    6 FORMAT(//,2X,' THEORETICAL FREQUENCIES IN CM-1'/
     1          2X,' NO.',9X,' HARMONIC',14X,' DELTA',
     2         13X,' FUNDAMENTAL'/)
    7 FORMAT(2X,I5,4X,F13.3,7X,F13.3,12X,F13.3)
    8 FORMAT(//,2X,' :::::::::::::::::::::::::::::::::::::::::::::::'/
     1          2X,' :::ANHARMONIC RESONANCE (TYPE I) IS OBSERVED:::'/
     2          2X,' :::::::::::::::::::::::::::::::::::::::::::::::'/)
    9 FORMAT(//,2X,' ::::::::::::::::::::::::::::::::::::::::::::::::'/
     1          2X,' :::ANHARMONIC RESONANCE (TYPE II) IS OBSERVED:::'/
     2          2X,' ::::::::::::::::::::::::::::::::::::::::::::::::'/)
   10 FORMAT(2X,' NUMBER OF OCCURENCES IS ',I5/
     1       2X,' RESONANCE LIMIT IS SET TO ',F10.3,' CM-1'/)
C
C:::::::::::::::::::::::::::::::::::::::
C:::CALCULATE ANHARMONICITY CONSTANTS:::
C:::::::::::::::::::::::::::::::::::::::
      IANH1=0
      IANH2=0
C
      WRITE(6,1)
      IF(IPRNT.GT.2)
     *WRITE(6,2)
      DO 110 IVIB=1,NVIB
      FRQI=FREQ(IVIB)
      FRQI2=FRQI*FRQI
      DO 110 JVIB=1,NVIB
      FRQJ=FREQ(JVIB)
      FRQJ2=FRQJ*FRQJ
C
      VALU1=A00
      VALU2=A00
      VALU3=A00
      VALU4=A00
      XRS=A00
      IF(DABS(FRQI).LE.WLIMIT) GO TO 210
      IF(DABS(FRQJ).LE.WLIMIT) GO TO 210
      IF(IVIB.NE.JVIB) GO TO 205
C
C++++++++++++++++++++++++
C+++DIAGONAL CONSTANTS+++
C++++++++++++++++++++++++
      VALU1=F4Q(IVIB,IVIB,IVIB,IVIB)
      DO 101 KVIB=1,NVIB
      FRQK=FREQ(KVIB)
      IF(DABS(FRQK).LE.WLIMIT) GO TO 101
      FRQK2=FRQK*FRQK
      FAC=F3Q(IVIB,IVIB,KVIB)
      FACT=FAC*FAC
      IF(DABS(FRQI+FRQI-FRQK).LE.FLIM1) GO TO 201
      DNUM=FRQI2*8.0D+00-FRQK2*3.0D+00
      DENM=FRQI2*4.0D+00-FRQK2
      VALU2=VALU2+(DNUM*FACT)/(FRQK*DENM)
C***********************************************************************
C***FOLLWOING CODES ARE FOR AN ALTERNATE EXPRESSION                    *
C***THEY ARE TESTED TO GIVE THE SAME RESUTLS                           *
C***********************************************************************
C*    DNUM=ONE/((FRQI+FRQI+FRQK)*2.0D+00)-ONE/((FRQI+FRQI-FRQK)*2.0D+00)
C*   1    +2.0D+00/FRQK                                                *
C*    VALU2=VALU2+(DNUM*FACT)                                          *
C***********************************************************************
      GO TO 101
C/////////////////////////////////
C///FERMI RESONANCE CORRECTIONS///
C/////////////////////////////////
  201 CONTINUE
      DNUM=ONE/((FRQI+FRQI+FRQK)*2.0D+00)+2.0D+00/FRQK
      IF(DABS(FAC).GT.PLIMIT) THEN
        IANH1=IANH1+1
      END IF
      VALU2=VALU2+FACT*DNUM
C////////////////////////////
C///END OF FERMI RESONANCE///
C////////////////////////////
  101 CONTINUE
      XRS=(VALU1-VALU2)/16.0D+00
      GO TO 210
C
C++++++++++++++++++++++++++++
C+++OFF-DIAGONAL CONSTANTS+++
C++++++++++++++++++++++++++++
  205 CONTINUE
      VALU1=F4Q(IVIB,IVIB,JVIB,JVIB)
      DO 102 KVIB=1,NVIB
      FRQK=FREQ(KVIB)
      IF(DABS(FRQK).LE.WLIMIT) GO TO 102
      FRQK2=FRQK*FRQK
      WPPP=FRQI+FRQJ+FRQK
      WPPM=FRQI+FRQJ-FRQK
      WPMP=FRQI-FRQJ+FRQK
      WPMM=FRQI-FRQJ-FRQK
      DIJK=WPPP*WPPM*WPMP*WPMM
      DNUM2=F3Q(IVIB,IVIB,KVIB)*F3Q(KVIB,JVIB,JVIB)
      VALU2=VALU2+DNUM2/FRQK
      FAC=F3Q(IVIB,JVIB,KVIB)
      FACT=FAC*FAC
      IF(DABS(WPPM).LE.FLIM2) GO TO 206
      IF(DABS(WPMP).LE.FLIM2) GO TO 207
      IF(DABS(WPMM).LE.FLIM2) GO TO 208
      DNUM3=(FRQK2-FRQI2-FRQJ2)*FRQK
      VALU3=VALU3+(DNUM3*FACT)/DIJK
C***********************************************************************
C***FOLLOWING CODES ARE FOR AN ALTERNATE EXPRESSION                    *
C***THEY ARE TESTED TO GIVE THE SAME RESULTS                           *
C***********************************************************************
C*    DNUM3=ONE/WPPP-ONE/WPPM+ONE/WPMP-ONE/WPMM                        *
C*    VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                 *
C***********************************************************************
      GO TO 102
C/////////////////////////////////
C///FERMI RESONANCE CORRECTIONS///
C/////////////////////////////////
C   WI+WJ=WK
  206 CONTINUE
      DNUM3=ONE/WPPP+ONE/WPMP-ONE/WPMM
      IF(DABS(FAC).GT.PLIMIT) THEN
        IANH2=IANH2+1
      END IF
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00
      GO TO 102
C   WI+WK=WJ
  207 CONTINUE
      DNUM3=ONE/WPPP-ONE/WPPM-ONE/WPMM
      IF(DABS(FAC).GT.PLIMIT) THEN
        IANH2=IANH2+1
      END IF
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00
      GO TO 102
C   WJ+WK=WI
  208 CONTINUE
      DNUM3=ONE/WPPP-ONE/WPPM+ONE/WPMP
      IF(DABS(FAC).GT.PLIMIT) THEN
        IANH2=IANH2+1
      END IF
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00
C////////////////////////////
C///END OF FERMI RESONANCE///
C////////////////////////////
  102 CONTINUE
      ZFACT=A00
      DO 103 IABC=1,3
      II=IAXIS(IABC)
      ZET=ZETA(IVIB,JVIB,II)
      ZFACT=ZFACT+ZET*ZET*ROTCM(IABC)
  103 CONTINUE
      VALU4=(FRQI/FRQJ+FRQJ/FRQI)*ZFACT
      XRS=(VALU1-VALU2)/4.0D+00-VALU3/2.0D+00+VALU4
C
  210 CONTINUE
      XIJ(IVIB,JVIB)=XRS
      IF(IPRNT.LE.2) GO TO 110
      WRITE(6,3) IVIB,JVIB,FRQI,FRQJ,VALU1,VALU2,VALU3,VALU4,XRS
  110 CONTINUE
      IF(IANH1.NE.0) THEN
        WRITE(6,8)
        WRITE(6,10) IANH1,FLIM1
      END IF
      IF(IANH2.NE.0) THEN
        WRITE(6,8)
        WRITE(6,10) IANH2,FLIM2
      END IF
C
      IF(IPRNT.GT.2) GO TO 211
      WRITE(6,4)
      DO 112 II=1,NVIB
      IVIB=NFRQ(II)
      FRQI=FREQ(IVIB)
      DO 112 JJ=1,II
      JVIB=NFRQ(JJ)
      FRQJ=FREQ(JVIB)
      XRS=XIJ(IVIB,JVIB)
      WRITE(6,5) II,JJ,FRQI,FRQJ,XRS
  112 CONTINUE
C
C:::::::::::::::::::::::::::::::::::::::
C:::CALCULATE FUNDAMENTAL FREQUENCIES:::
C:::::::::::::::::::::::::::::::::::::::
  211 CONTINUE
      WRITE(6,6)
      DO 115 II=1,NVIB
      IVIB=NFRQ(II)
      VALU1=FREQ(IVIB)
      VALU2=XIJ(IVIB,IVIB)*2.0D+00
      VALU3=A00
      DO 114 JJ=1,NVIB
      JVIB=NFRQ(JJ)
      IF(IVIB.EQ.JVIB) GO TO 114
      VALU3=VALU3+XIJ(IVIB,JVIB)
  114 CONTINUE
      DELTA=VALU2+VALU3*0.5D+00
      VALUT=VALU1+DELTA
      WRITE(6,7) II,VALU1,DELTA,VALUT
  115 CONTINUE
C
      RETURN
      end
