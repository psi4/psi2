C***********************************************************************
C* For array size changes from/to  (50625) or (10000) search C11-14-88 *
C* (50625) is for running in 15Mb while (10000) is for 7Mb.    RBR     *
C***********************************************************************
C* THIS SMALL CORE VERSION FOR THE 9370 WILL ONLY DO UP TO F FUNCTIONS.*
C* CHANGED ABDENS(50625) TO (10000), derint(50625,9) to (10000,9)      *
C***********************************************************************
C* LAST UPDATED NOVEMBER 11, 1988 BY TRACY HAMILTON FOR D5,D6,F7,F10   *
C* G9,G15 FUNCTIONS.  SEARCH FOR CTPH.                                 *
C***********************************************************************
C* LAST UPDATED NOVEMBER 11, 1987 BY RICHARD REMINGTON                 *
C* VARIABLE DLIM CHANGED FROM 1.0E-10 TO 1.0E-16 FOR BETTER PRECISION; *
C* ALSO PUT DLIM IN COMMON TOLER                                       *
C***********************************************************************
C* LAST UPDATED ON APRIL 25, 1985 BY YUKIO YAMAGUCHI                   *
C***********************************************************************
      SUBROUTINE tder2
ctph********************************************************************
c     may need to have h and w expanded for 2nd derivatives of g's     *
c     done november 11, 1988.                                          *
ctph********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/DSTVG/T,X0,Y0,Z0,NI,NJ,CX,CY,CZ,TAI,TAJ,TAA,UU
      COMMON/DSTVR/VX,VY,VZ,DVX,DVY,DVZ,DDVX,DDVY,DDVZ,WIVX,WIVY,WIVZ
     #,            WIWIVX,WIWIVY,WIWIVZ,WIWJVX,WIWJVY,WIWJVZ
     #,            WJWJVX,WJWJVY,WJWJVZ,WJVX,WJVY,WJVZ
     #,            WIXDVX,WIYDVY,WIZDVZ,WJXDVX,WJYDVY,WJZDVZ
     #,            DEXDXI,DEXDYI,DEXDZI,DEXDXJ,DEXDYJ,DEXDZJ
     #,            DEXDII,DEXDIJ,DEXDJJ,DEXDIC,DEXDJC
      COMMON/IJATM/IATOM,JATOM,ICATOM
      COMMON/IJXYZ/XI,YI,ZI,XJ,YJ,ZJ
      COMMON/NDERS/NDERIV
      COMMON/TOLER/DLIM
ctph  DIMENSION H(21),W(21),LMIN(6),LMAX(6)
      DIMENSION H(28),W(28),LMIN(7),LMAX(7)
      DIMENSION COORD0(3),COORDI(3)
      DIMENSION COORDJ(3),COORDC(3),V(3),DV(3),DDV(3),WIV(3),WIWIV(3)
      DIMENSION WIWJV(3),WJWJV(3),WIDV(3),DEXDI(3),DEXDJ(3),WJV(3)
      DIMENSION WJDV(3)
      EQUIVALENCE (V,VX),(DV,DVX),(DDV,DDVX),(WIV,WIVX),(WIWIV,WIWIVX)
      EQUIVALENCE (WIWJV,WIWJVX),            (WJV,WJVX)
      EQUIVALENCE (WJWJV,WJWJVX),(WIDV,WIXDVX),(WJDV,WJXDVX)
      EQUIVALENCE (COORD0,X0),(COORDI,XI),(COORDJ,XJ),(COORDC,CX)
      EQUIVALENCE (DEXDI,DEXDXI),(DEXDJ,DEXDXJ)
      DATA ZERO,ONE,TWO / 0.0D+00 , 1.0D+00 , 2.0D+00 /
ctph  DATA LMIN /1,  2,  4,  7, 11, 16/
ctph  DATA LMAX /1,  3,  6, 10, 15, 21/
      DATA LMIN /1,  2,  4,  7, 11, 16, 22/
      DATA LMAX /1,  3,  6, 10, 15, 21, 28/
CTPH  TABLE 25.10 OF HANDBOOK OF MATHEMATICAL FUNCTIONS, ED.
CTPH                 MILTON ABRAMOWITZ AND IRENE A. STEGUN
      DATA H   /0.0D+00
     1,         -.707106781186548D+00,  0.707106781186548D+00
     2,         -1.22474487139159D+00,0.0D+00,1.22474487139159D+00
     3,         -1.65068012388578D+00, -0.524647623275290D+00
     4,          0.524647623275290D+00, 1.65068012388578D+00
     5,   -2.02018287045609D+00,-0.958572464613819D+00,0.0D+00
     6,          0.958572464613819D+00, 2.02018287045609D+00
     7,         -2.350604973674D+00  , -1.335849074014D+00
     8,         -0.436077411928D+00  ,  0.436077411928D+00
     9,          1.335849074014D+00  ,  2.350604973674D+00
     A,         -2.651961356835233D+00 , -1.673551628767471D+00
     B,         -0.816287882858965D+00 , 0.0D00, 0.816287882858965D+00
     C,          1.673551628767471D+00 , 2.651961356835233D+00/
      DATA W   /1.77245385090552D+00
     1,         0.8862269254528D+00  ,  0.8862269254528D+00
     2,         0.2954089751509D+00  ,  1.181635900604D+00
     3,         0.2954089751509D+00
     4,         8.131283544725D-02   ,  8.049140900055D-01
     5,         8.049140900055D-01   ,  8.131283544725D-02
     6,         1.995324205905D-02   ,  3.936193231522D-01
     7,         9.453087204829D-01   ,  3.936193231522D-01
     8,         1.995324205905D-02
     9,         4.530009905509D-03   ,  1.570673203229D-01
     A,         7.246295952244D-01   ,  7.246295952244D-01
     B,         1.570673203229D-01   ,  4.530009905509D-03
     C,         9.717812450995D-04   ,  5.451558281913D-02
     D,         4.256072526101D-01   ,  8.102646175568D-01
     E,         4.256072526101D-01   ,  5.451558281913D-02
     F,         9.717812450995D-04/
C
cges  ENTRY TDER2
C
      NPTS=(NI+NJ+NDERIV)/2+1
      TWOAI=TAI
      TWOAJ=TAJ
      NA=NI-1
      NB=NJ-1
      DO 106 IXYZ=1,3
      V(IXYZ)=ZERO
      DV(IXYZ)=ZERO
      WIV(IXYZ)=ZERO
      WIDV(IXYZ)=ZERO
      IF(NDERIV.LE.1) GO TO 203
      WIWIV(IXYZ)=ZERO
      DDV(IXYZ)=ZERO
C
  203 CONTINUE
      DO 105 IPT=LMIN(NPTS),LMAX(NPTS)
      PT=H(IPT)/T+COORD0(IXYZ)
      A=PT-COORDI(IXYZ)
      B=PT-COORDJ(IXYZ)
      WT=W(IPT)
      IF(NA.GT.2) THEN
      ANAM2=A**(NA-2)
      ANAM1=ANAM2*A
      ANA  =ANAM1*A
      ANAP1=ANA*A
      ELSE IF(NA.EQ.2) THEN
      ANAM2=ONE
      ANAM1=A
      ANA  =A*A
      ANAP1=ANA*A
      ELSE IF(NA.EQ.1) THEN
      ANAM2=ZERO
      ANAM1=ONE
      ANA  =A
      ANAP1=A*A
      ELSE IF(NA.EQ.0) THEN
      ANAM2=ZERO
      ANAM1=ZERO
      ANA  =ONE
      ANAP1=A
      ELSE
      STOP 'NEGATIVE VALUE OF NA IN TDER'
      END IF
      IF(NB.GT.2) THEN
      BNBM2=B**(NB-2)
      BNBM1=BNBM2*B
      BNB  =BNBM1*B
      BNBP1=BNB*B
      ELSE IF(NB.EQ.2) THEN
      BNBM2=ONE
      BNBM1=B
      BNB  =B*B
      BNBP1=BNB*B
      ELSE IF(NB.EQ.1) THEN
      BNBM2=ZERO
      BNBM1=ONE
      BNB  =B
      BNBP1=B*B
      ELSE IF(NB.EQ.0) THEN
      BNBM2=ZERO
      BNBM1=ZERO
      BNB  =ONE
      BNBP1=B
      ELSE
      STOP 'NEGATIVE VALUE OF NB IN TDER'
      END IF
      ANAP2=ANAP1*A
      BNBP2=BNBP1*B
      TWOAJ2=TWOAJ*TWOAJ
      TWOAI2=TWOAI*TWOAI
      DELSQB=NB*(NB-1)*BNBM2-TWOAJ*(2*NB+1)*BNB+TWOAJ2*BNBP2
C
      V(IXYZ)=V(IXYZ) + ANA*BNB*WT
      DV(IXYZ)=DV(IXYZ) + ANA*DELSQB*WT
      WIV(IXYZ)=WIV(IXYZ) + (TWOAI*ANAP1-NA*ANAM1)*BNB*WT
      WIDV(IXYZ)=WIDV(IXYZ) + (TWOAI*ANAP1-NA*ANAM1)*DELSQB*WT
      IF(NDERIV.LE.1) GO TO 105
      WIWIV(IXYZ)=WIWIV(IXYZ) + (NA*(NA-1)*ANAM2-TWOAI*(2*NA+1)*
     *      ANA+                 TWOAI2*ANAP2)*BNB*WT
      DDV(IXYZ)=DDV(IXYZ) + (NA*(NA-1)*ANAM2-TWOAI*(2*NA+1)*ANA
     #      +                    TWOAI2*ANAP2)*DELSQB*WT
  105 CONTINUE
  106 CONTINUE
      RETURN
      END
