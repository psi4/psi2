      SUBROUTINE LINEAR(ZETA,F3Q,F4Q,XIJ,GIJ)
C   THE ANHARMONICITY CONTSTANTS FOR A LINEAR POLYATOMIC MOLECULE
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F3Q(N3N,N3N,N3N),F4Q(N3N,N3N,N3N,N3N),ZETA(N3N,N3N,3)
      DIMENSION XIJ(N3N,N3N),GIJ(N3N,N3N)
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)
      COMMON/VIB109/CLIMIT,FLIM1,FLIM2
      COMMON/VIB203/IOFF(150),IPRNT
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)
      COMMON/VIB205/IFREQ,NFRQ(150)
      DATA A00,ONE / 0.0D+00 , 1.0D+00 /
      DATA WLIMIT / 1.0D+00 /
      DATA DLIMIT / 1.0D-01 /
      DATA PLIMIT / 1.0D-05 /
    1 FORMAT(//,2X,' *********************************************'/
     1          2X,' ***THE VIBRATIONAL ANHARMONICITY CONSTANTS***'/
     2          2X,' *********************************************'/)
    2 FORMAT(//,2X,' ANHARMONIC CONSTANTS IN CM-1'/
     1          2X,' IVIB',3X,' JVIB',8X,' FREQ(I)',5X,' FREQ(J)',
     2          7X,' VALU1',7X,' VALU2',7X,' VALU3',7X,' VALU4',
     3          7X,' XIJ',9X,' GIJ'/)
    3 FORMAT(2X,I5,3X,I5,3X,8(F13.3))
    4 FORMAT(//,2X,' ANHARMONIC CONSTANTS IN CM-1'/
     1          4X,' QI',5X,' QJ',8X,' FREQ(I)',5X,' FREQ(J)',
     2          7X,' XIJ',9X,' GIJ'/)
    5 FORMAT(2X,I5,3X,I5,3X,4(F13.3))
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
   10 FORMAT(2X,' RESONANCE LIMIT IS SET TO ',F10.3,' CM-1'/)
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
      IAA=IAXIS(1)
      IBB=IAXIS(2)
      ICC=IAXIS(3)
C
      DO 120 IVIB=1,N3N
      FRQI=FREQ(IVIB)
      FRQI2=FRQI*FRQI
      IVA=NDAB(IVIB,1)
      IVB=NDAB(IVIB,2)
      DO 120 JVIB=1,N3N
      FRQJ=FREQ(JVIB)
      FRQJ2=FRQJ*FRQJ
      JVA=NDAB(JVIB,1)
      JVB=NDAB(JVIB,2)
C
      VALU1=A00
      VALU2=A00
      VALU3=A00
      VALU4=A00
      VALU5=A00
      GVAL1=A00
      GVAL2=A00
      XRS=A00
      GRS=A00
      IF(DABS(FRQI).LE.WLIMIT) GO TO 230
      IF(DABS(FRQJ).LE.WLIMIT) GO TO 230
      IF(IVIB.NE.JVIB) GO TO 210
C
C+++++++++++++++++++++++++++++++++++++++++
C+++DIAGONAL CONSTANTS (NON-DEGENERATE)+++
C+++++++++++++++++++++++++++++++++++++++++
      IF(NDEG(IVIB).NE.1) GO TO 205
      VALU1=F4Q(IVIB,IVIB,IVIB,IVIB)
      DO 101 KVIB=1,N3N
C   SUM OVER ONLY NON-DEGENERATE MODES
      FRQK=FREQ(KVIB)
      IF(DABS(FRQK).LE.WLIMIT) GO TO 101
      IF(NDEG(KVIB).NE.1) GO TO 101
      FRQK2=FRQK*FRQK
      FAC=F3Q(IVIB,IVIB,KVIB)
      FACT=FAC*FAC
      IF(DABS(FRQI+FRQI-FRQK).LE.FLIM1) GO TO 201
      DNUM=FRQI2*8.0D+00-FRQK2*3.0D+00
      DENM=FRQI2*4.0D+00-FRQK2
      VALU2=VALU2+(DNUM*FACT)/(FRQK*DENM)
C***********************************************************************
C***FOLLOWING CODES ARE FOR AN ALTERNATE EXPRESSION                    *
C***THEY ARE TESTED TO GIVE THE SAME RESULTS                           *
C***********************************************************************
C*    DNUM=ONE/((FRQI+FRQI+FRQK)*2.0D+00)-ONE/((FRQI+FRQI-FRQK)*2.0D+00)
C*   *    +2.0D+00/FRQK                                                *
C     VALU2=VALU2+(DNUM*FACT)                                          *
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
      GO TO 230
C
C+++++++++++++++++++++++++++++++++++++
C+++DIAGONAL CONSTANTS (DEGENERATE)+++
C+++++++++++++++++++++++++++++++++++++
  205 CONTINUE
      VALU1=F4Q(IVIB,IVIB,IVIB,IVIB)
      GVAL1=-VALU1
      DO 102 KVIB=1,N3N
C   SUM OVER NON-DEGENERATE MODES
      FRQK=FREQ(KVIB)
      IF(DABS(FRQK).LE.WLIMIT) GO TO 102
      IF(NDEG(KVIB).NE.1) GO TO 102
      FRQK2=FRQK*FRQK
      FAC=F3Q(IVIB,IVIB,KVIB)
      FACT=FAC*FAC
      IF(DABS(FRQI+FRQI-FRQK).LE.FLIM1) GO TO 206
      DNUM=FRQI2*8.0D+00-FRQK2*3.0D+00
      DENM=FRQI2*4.0D+00-FRQK2
      VALU2=VALU2+(DNUM*FACT)/(FRQK*DENM)
      GVAL2=GVAL2+(FACT*FRQK)/DENM
      GO TO 102
C/////////////////////////////////
C///FERMI RESONANCE CORRECTIONS///
C/////////////////////////////////
  206 CONTINUE
      DNUM=ONE/((FRQI+FRQI+FRQK)*2.0D+00)+2.0D+00/FRQK
      IF(DABS(FAC).GT.PLIMIT) THEN
        IANH1=IANH1+1
      END IF
      VALU2=VALU2+FACT*DNUM
      GNUM=ONE/((FRQI+FRQI+FRQK)*2.0D+00)
      GVAL2=GVAL2-FACT*GNUM
C////////////////////////////
C///END OF FERMI RESONANCE///
C////////////////////////////
  102 CONTINUE
      XRS=(VALU1-VALU2)/16.0D+00
      GRS=(GVAL1/3.0D+00-GVAL2)/16.0D+00
      GO TO 230
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++OFF-DIAGONAL CONSTANTS (NON-DEGENERATE & NON-DEGENERATE)+++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  210 CONTINUE
      IF(NDEG(IVIB).NE.1) GO TO 215
      IF(NDEG(JVIB).NE.1) GO TO 215
      VALU1=F4Q(IVIB,IVIB,JVIB,JVIB)
      DO 103 KVIB=1,NVIB
C   SUM OVER ONLY NON-DEGENERATE MODES
      FRQK=FREQ(KVIB)
      IF(DABS(FRQK).LE.WLIMIT) GO TO 103
      IF(NDEG(KVIB).NE.1) GO TO 103
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
      IF(DABS(WPPM).LE.FLIM2) GO TO 211
      IF(DABS(WPMP).LE.FLIM2) GO TO 212
      IF(DABS(WPMM).LE.FLIM2) GO TO 213
      DNUM3=(FRQK2-FRQI2-FRQJ2)*FRQK
      VALU3=VALU3+(DNUM3*FACT)/DIJK
C***********************************************************************
C***FOLLOWING CODES ARE FOR AN ALTERNATIVE EXPRESSION                  *
C***THEY ARE TESTED TO GIVE THE SAME RESULTS                           *
C***********************************************************************
C*    DNUM3=ONE/WPPP-ONE/WPPM+ONE/WPMP-ONE/WPMM                        *
C*    VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                 *
C***********************************************************************
      GO TO 103
C/////////////////////////////////
C///FERMI RESONANCE CORRECTIONS///
C/////////////////////////////////
C   WI+WJ=WK
  211 CONTINUE
      DNUM3=ONE/WPPP+ONE/WPMP-ONE/WPMM
      IF(DABS(FAC).GT.PLIMIT) THEN
        IANH2=IANH2+1
      END IF
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00
      GO TO 103
C   WI+WK=WJ
  212 CONTINUE
      DNUM3=ONE/WPPP-ONE/WPPM-ONE/WPMM
      IF(DABS(FAC).GT.PLIMIT) THEN
        IANH2=IANH2+1
      END IF
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00
      GO TO 103
C   WJ+WK=WI
  213 CONTINUE
      DNUM3=ONE/WPPP-ONE/WPPM+ONE/WPMP
      IF(DABS(FAC).GT.PLIMIT) THEN
        IANH2=IANH2+1
      END IF
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00
C////////////////////////////
C///END OF FERMI RESONANCE///
C////////////////////////////
  103 CONTINUE
      XRS=(VALU1-VALU2-VALU3*2.0D+00)/4.0D+00
      GO TO 230
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++OFF-DIAGONAL CONSTANTS (DEGENERATE & NON-DEGENERATE)+++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   NOTE : RELATIONSHIPS AMONG CUBIC FORCE CONSTANTS
C     (S,TA,TA)=(S,TB,TB)
C     (S,TA,TB)=0
C     (S,TA,TA')=(S,TB,TB')
C     (S,TA,TB')=0
  215 CONTINUE
      IF(NDEG(IVIB).NE.1.AND.NDEG(JVIB).NE.1) GO TO 225
      VALU1=F4Q(IVIB,IVIB,JVIB,JVIB)
C   SUM OVER NON-DEGENERATE MODES
      DO 105 KVIB=1,NVIB
      FRQK=FREQ(KVIB)
      IF(DABS(FRQK).LE.WLIMIT) GO TO 105
      IF(NDEG(KVIB).NE.1) GO TO 105
      DNUM2=F3Q(IVIB,IVIB,KVIB)*F3Q(KVIB,JVIB,JVIB)
      VALU2=VALU2+DNUM2/FRQK
  105 CONTINUE
C   SUM OVER DEGENERATE MODES
      DO 106 KVIB=1,NVIB
      FRQK=FREQ(KVIB)
      IF(DABS(FRQK).LE.WLIMIT) GO TO 106
      FRQK2=FRQK*FRQK
      IF(NDEG(KVIB).EQ.1) GO TO 106
      WPPP=FRQI+FRQJ+FRQK
      WPPM=FRQI+FRQJ-FRQK
      WPMP=FRQI-FRQJ+FRQK
      WPMM=FRQI-FRQJ-FRQK
      DIJK=WPPP*WPPM*WPMP*WPMM
      DNUM3=(FRQK2-FRQI2-FRQJ2)*FRQK
      FAC=F3Q(IVIB,JVIB,KVIB)
      FACT=FAC*FAC
      IF(DABS(WPPM).LE.FLIM2) GO TO 217
      IF(DABS(WPMP).LE.FLIM2) GO TO 218
      IF(DABS(WPMM).LE.FLIM2) GO TO 219
      VALU3=VALU3+(DNUM3*FACT)/DIJK
      GO TO 106
C/////////////////////////////////
C///FERMI RESONANCE CORRECTIONS///
C/////////////////////////////////
C   WI+WJ=WK
  217 CONTINUE
      DNUM3=ONE/WPPP+ONE/WPMP-ONE/WPMM
      IF(DABS(FAC).GT.PLIMIT) THEN
        IANH2=IANH2+1
      END IF
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00
      GO TO 106
C   WI+WK=WJ
  218 CONTINUE
      DNUM3=ONE/WPPP-ONE/WPPM-ONE/WPMM
      IF(DABS(FAC).GT.PLIMIT) THEN
        IANH2=IANH2+1
      END IF
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00
      GO TO 106
C   WJ+WK=WI
  219 CONTINUE
      DNUM3=ONE/WPPP-ONE/WPPM+ONE/WPMP
      IF(DABS(FAC).GT.PLIMIT) THEN
        IANH2=IANH2+1
      END IF
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00
C////////////////////////////
C///END OF FERMI RESONANCE///
C////////////////////////////
  106 CONTINUE
      ZFACT=(ZETA(IVIB,JVIB,IAA)**2+ZETA(IVIB,JVIB,IBB)**2)*ROTCM(IBB)
      VALU4=(FRQI/FRQJ+FRQJ/FRQI)*ZFACT
      XRS=(VALU1-VALU2-VALU3*2.0D+00)/4.0D+00+VALU4
      GO TO 230
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++OFF-DIAGONAL CONSTANTS (DEGENERATE & DEGENERATE)+++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++
  225 CONTINUE
      IF(DABS(FRQI-FRQJ).LE.DLIMIT) GO TO 230
      VALU1=(F4Q(IVA,IVA,JVA,JVA)+F4Q(IVA,IVA,JVB,JVB))*0.5D+00
      DO 109 KVIB=1,NVIB
C   SUM OVER NON-DEGENERATE MODES
      FRQK=FREQ(KVIB)
      IF(DABS(FRQK).LE.WLIMIT) GO TO 109
      IF(NDEG(KVIB).NE.1) GO TO 109
      FRQK2=FRQK*FRQK
      WPPP=FRQI+FRQJ+FRQK
      WPPM=FRQI+FRQJ-FRQK
      WPMP=FRQI-FRQJ+FRQK
      WPMM=FRQI-FRQJ-FRQK
      DIJK=WPPP*WPPM*WPMP*WPMM
      DNUM2=F3Q(IVIB,IVIB,KVIB)*F3Q(KVIB,JVIB,JVIB)
      DNUM3=(FRQK2-FRQI2-FRQJ2)*FRQK
      FAC=F3Q(IVA,JVA,KVIB)
      FACT=FAC*FAC
      VALU2=VALU2+DNUM2/FRQK
      IF(DABS(WPPM).LE.FLIM2) GO TO 226
      IF(DABS(WPMP).LE.FLIM2) GO TO 227
      IF(DABS(WPMM).LE.FLIM2) GO TO 228
      VALU3=VALU3+(DNUM3*FACT)/DIJK
      GVAL1=GVAL1+(FRQI*FRQJ*FRQK*FACT)/DIJK
      GO TO 109
C/////////////////////////////////
C///FERMI RESONANCE CORRECTIONS///
C/////////////////////////////////
C   WI+WJ=WK
  226 CONTINUE
      DNUM3=ONE/WPPP+ONE/WPMP-ONE/WPMM
      IF(DABS(FAC).GT.PLIMIT) THEN
        IANH2=IANH2+1
      END IF
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00
      GFACT=ONE/WPPP-ONE/WPMP+ONE/WPMM
      GVAL1=GVAL1+(GFACT*FACT)/8.0D+00
      GO TO 109
C   WI+WK=WJ
  227 CONTINUE
      DNUM3=ONE/WPPP-ONE/WPPM-ONE/WPMM
      IF(DABS(FAC).GT.PLIMIT) THEN
        IANH2=IANH2+1
      END IF
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00
      GFACT=ONE/WPPP-ONE/WPPM+ONE/WPMM
      GVAL1=GVAL1+(GFACT*FACT)/8.0D+00
      GO TO 109
C   WJ+WK=WI
  228 CONTINUE
      DNUM3=ONE/WPPP-ONE/WPPM+ONE/WPMP
      IF(DABS(FAC).GT.PLIMIT) THEN
        IANH2=IANH2+1
      END IF
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00
      GFACT=ONE/WPPP-ONE/WPPM-ONE/WPMP
      GVAL1=GVAL1+(GFACT*FACT)/8.0D+00
      GO TO 109
C////////////////////////////
C///END OF FERMI RESONANCE///
C////////////////////////////
  109 CONTINUE
      XRS=(VALU1-VALU2-VALU3)/4.0D+00
      GRS=GVAL1*0.5D+00
C
  230 CONTINUE
      XIJ(IVIB,JVIB)=XRS
      GIJ(IVIB,JVIB)=GRS
      IF(IPRNT.LE.2) GO TO 120
      WRITE(6,3) IVIB,JVIB,FRQI,FRQJ,VALU1,VALU2,VALU3,VALU4,XRS,GRS
  120 CONTINUE
C
      IF(IANH1.NE.0) THEN
        WRITE(6,8)
        WRITE(6,10) FLIM1
      END IF
      IF(IANH2.NE.0) THEN
        WRITE(6,9)
        WRITE(6,10) FLIM2
      END IF
C
C*    IF(IPRNT.GT.2) GO TO 231
      WRITE(6,4)
      DO 122 II=1,NVIB
      IVIB=NFRQ(II)
      FRQI=FREQ(IVIB)
C*    DO 122 JJ=1,II
      DO 122 JJ=1,NVIB
      JVIB=NFRQ(JJ)
      FRQJ=FREQ(JVIB)
      XRS=XIJ(IVIB,JVIB)
      GRS=GIJ(IVIB,JVIB)
      WRITE(6,5) II,JJ,FRQI,FRQJ,XRS,GRS
  122 CONTINUE
C
C:::::::::::::::::::::::::::::::::::::::
C:::CALCULATE FUNDAMENTAL FREQUENCIES:::
C:::::::::::::::::::::::::::::::::::::::
  231 CONTINUE
      WRITE(6,6)
      DO 125 II=1,NVIB
      IVIB=NFRQ(II)
      DEGI=NDEG(IVIB)
      VALU1=FREQ(IVIB)
      VALU2=XIJ(IVIB,IVIB)*(ONE+DEGI)
      VALU3=A00
      VALU4=GIJ(IVIB,IVIB)
      DO 124 JJ=1,NVIB
      JVIB=NFRQ(JJ)
      IF(IVIB.EQ.JVIB) GO TO 124
      VALU3=VALU3+XIJ(IVIB,JVIB)
  124 CONTINUE
      DELTA=VALU2+VALU3*0.5D+00+VALU4
      VALUT=VALU1+DELTA
      WRITE(6,7) II,VALU1,DELTA,VALUT
  125 CONTINUE
C
      RETURN
      end
