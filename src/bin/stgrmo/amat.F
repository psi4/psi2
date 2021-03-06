      SUBROUTINE AMAT(NIJ,C,ZETA,IND,LBLI,BUFI)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NIJ(NTRI,2),C(NDTRI)
      DIMENSION ZETA(NTRI,NTYPEP)
      DIMENSION IND(NTRI)
      DIMENSION LBLI(MAXBF4),BUFI(MAXBF2)
      COMMON/BASIS/NBASIS,NTRI
      COMMON/FUNCS/NTYPES,NTYPEP
      COMMON/GVBAB/A(10,10),B(10,10)
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4
      COMMON/SIGNS/IOFF(32897),IPRNT
      COMMON/START/FOCC(10),NSORB(10),NSTR(10),NEND(10),MOTYP(256)
      COMMON/CI102/NIND,NDEP,NDTRI
      DATA A00 / 0.0D+00 /
    1 FORMAT(//,2X,' ZETA MATRIX, ITYP = ',I5/)
    2 FORMAT(//,2X,' THE A MATRIX (TWO-ELECTRON PART)'/)
    3 FORMAT(//,2X,' THE TOTAL A MATRIX'/)
C
C   READ IN THE GENERALIZED LAGRANGIAN MATRICES
      CALL ZERO(ZETA,NTRI*NTYPEP)
      DO 101 ITYP=1,NTYPES
      CALL MREAD(ZETA(1,ITYP),26+ITYP)
      IF(IPRNT.LE.3) GO TO 101
      WRITE(6,1) ITYP
      CALL PRINT(ZETA(1,ITYP),NTRI,NBASIS,6)
  101 CONTINUE
C
C   FORM THE A MATRIX
      CALL IZERO(IND,NTRI)
      DO 102 II=1,NIND
      I=NIJ(II,1)
      J=NIJ(II,2)
      IJ=IOFF(I)+J
      IND(IJ)=II
  102 CONTINUE
C
C   CALCULATE TWO-ELECTRON PART FOR INDEPENDENT PAIRS
      CALL PQMAT(C,IND,LBLI,BUFI,NDTRI)
      IF(IPRNT.LE.2) GO TO 301
      WRITE(6,2)
      CALL PRINT(C,NDTRI,NIND,6)
C
C   ADD THE TERMS DUE TO DIFFERENCE IN THE GENERALIZED LAGRANGIANS
  301 CONTINUE
      DO 105 II=1,NIND
      I=NIJ(II,1)
      J=NIJ(II,2)
      ITYP=MOTYP(I)
      JTYP=MOTYP(J)
      DO 104 JJ=1,II
      K=NIJ(JJ,1)
      L=NIJ(JJ,2)
      KTYP=MOTYP(K)
      LTYP=MOTYP(L)
      ZVAL1=A00
      ZVAL2=A00
      ZVAL3=A00
      ZVAL4=A00
      IF(J.NE.K) GO TO 201
      IL=I21D(I,L)
      ZVAL1=ZETA(IL,ITYP)-ZETA(IL,JTYP)
  201 IF(I.NE.K) GO TO 202
      JL=I21D(J,L)
      ZVAL2=-ZETA(JL,JTYP)+ZETA(JL,ITYP)
  202 IF(J.NE.L) GO TO 203
      IK=I21D(I,K)
      ZVAL3=-ZETA(IK,ITYP)+ZETA(IK,JTYP)
  203 IF(I.NE.L) GO TO 204
      JK=I21D(J,K)
      ZVAL4=ZETA(JK,JTYP)-ZETA(JK,ITYP)
  204 ZVALT=ZVAL1+ZVAL2+ZVAL3+ZVAL4
C
      IIJJ=IOFF(II)+JJ
      C(IIJJ)=C(IIJJ)+ZVALT
  104 CONTINUE
  105 CONTINUE
C
      IF(IPRNT.LE.4) GO TO 302
      WRITE(6,3)
      CALL PRINT(C,NDTRI,NIND,6)
C
  302 CONTINUE
      RETURN
      END
