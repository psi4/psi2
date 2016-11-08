C     ////////////////////////////////////////////////////////////
      SUBROUTINE EXTRL(NEXT,NMAT,SF,SF1,SF2,IPRCT,DAMPQ,IFEXT)
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C
      COMMON/IFDMP/DAMPSV,INFLG,MICMX,IDFIX
      DIMENSION SF(1),SF1(1),SF2(1)
      DAMP = DAMPQ
      IF(NEXT.EQ.1)GO TO 20
      IF(IPRCT.EQ.0) GO TO 1
      IF(NEXT.LT.3)GO TO 1
      IF((DAMP.EQ.0.0D0).AND.(DAMPSV.NE.0.0D0)) GO TO 1
      IF (MOD(NEXT+1,IPRCT) .LE. 1)  DAMP = DAMPSV
      IF(MOD((NEXT),IPRCT).EQ.0) GO TO 5
    1 DO 2 I=1,NMAT
         SF2(I)=SF1(I)
         SF(I)=(SF(I)+DAMP*SF2(I))/(1.0D0+DAMP)
         SF1(I)=SF(I)
    2 CONTINUE
      DAMPSV = DAMP
      DAMPQ = DAMP
      RETURN
    5 IFEXT=1
      CUTOFF=0.5D0
      IF(CUTOFF.GE.1.0D0) CUTOFF=0.99D0
      IF(CUTOFF.LE.0.0D0) CUTOFF=0.0D0
      DO 30 I=1,NMAT
         SF(I)=(SF(I)+DAMP*SF1(I))/(1.0D0+DAMP)
         SFI=SF(I)
         IF (ABS( SF1(I)-SF2(I)).EQ.0.0D0) GO TO 32
         EPS=(SF(I)-SF1(I))/(SF1(I)-SF2(I))
         IF(EPS.GT.CUTOFF)GO TO 31
         IF(EPS.LT.-1.0D0) GO TO 33
         IF(ABS((SF(I)-2.0D0*SF1(I)+SF2(I))*(SF(I)-SF1(I))).LT.
     X   1.0D-8*SF(I)*SF(I)) GO TO 32
         SF(I)=(SF(I)*SF2(I)-SF1(I)*SF1(I))/(SF(I)-2.0D0*SF1(I)+SF2(I)
     *   )
         GO TO 32
   31    SF(I)=SFI+(SF1(I)-SF2(I))*(1.0D0/(1.0D0-CUTOFF)-1.0D0-CUTOFF)
         GO TO 32
   33    SF(I)=(SFI+DAMP*SF1(I))/(1.0D0+DAMP)
   32    CONTINUE
         SF2(I)=SF1(I)
         SF1(I)=SF(I)
   30 CONTINUE
      DAMPSV = DAMPQ
      RETURN
   20 DO 21 I=1,NMAT
   21 SF1(I)=SF(I)
      DAMPSV = DAMP
      RETURN
      END
