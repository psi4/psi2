      SUBROUTINE POLAD2(ESO,UA,UF,SSAF,SSFG,U,T,DP,DQ,DM,BUFI,LBLI)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ESO(NBFAO,NBASIS)
      DIMENSION UA(NBASIS,NBASIS,N3N),UF(NBASIS,NBASIS,3)
      DIMENSION SSAF(NTRI,N3N,3),SSFG(NTRI,3,3)
      DIMENSION U(NBASIS,NBASIS),T(NBASIS,NBASIS)
      DIMENSION DP(NTRI,N3NP3,NTYPEP),DQ(NTRI,N3NP3,NTYPEP)
      DIMENSION DM(NTRI,N3NP3,NTYPEP)
      DIMENSION BUFI(MAXBF2),LBLI(MAXBF4)
      COMMON/BASIS/NBASIS,NBFAO,NTRI,NTRI2,NBATRI,NBASQ
      COMMON/EAFGS/AF1(3,3,150),AF2(3,3,150),AF3(3,3,150),AFT(3,3,150)
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4
      COMMON/SIGNS/IOFF(466),IPRNT
      COMMON/START/FOCC(10),NSORB(10),NSTR(10),NEND(10),MOTYP(256)
      DATA ZERO,TWO,FOUR / 0.0D+00 , 2.0D+00 , 4.0D+00 /
    1 FORMAT(//,2X,' DP MATRIX IN POLAD2'/)
    2 FORMAT(2X,3I5,6F10.6)
    3 FORMAT(//,2X,' AFT MATRIX, IABC = ',I5/)
C
C   FORM DENSITY LIKE MATRIX IN SO BASIS
      CALL PAMAT(UA,UF,DP,DQ,ESO,N3NP3)
C
C   FORM HALF-TRANSFORMED DENSITY LIKE MATRIX
      CALL PQMAT(DP,DQ,DM,N3NP3*NTYPES,BUFI,LBLI)
C
C   COMPLETE FINAL TRANSFORMATION
      DO 101 IABC=1,N3NP3
      DO 101 ITYP=1,NTYPES
  101 CALL MOCONS(DM(1,IABC,ITYP),NTRI,DP(1,IABC,ITYP),NTRI,ESO,U,T)
      IF(IPRNT.LE.3) GO TO 201
      WRITE(6,1)
      CALL MATOUT(DP,NTRI,N3NP3*NTYPEP,NTRI,N3NP3*NTYPEP,6)
C
  201 CONTINUE
C:::A COORDINATE:::
      DO 130 IABC=1,N3N
C
C:::F COORDINATE:::
      DO 120 JXYZ=1,3
      JPQR=JXYZ+N3N
C
C:::G COORDINATE:::
      DO 110 KXYZ=1,3
      KPQR=KXYZ+N3N
C
      VALU1=ZERO
      VALU2=ZERO
C
      DO 103 I=1,NBASIS
      ITYP=MOTYP(I)
      DO 103 J=1,NBASIS
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)
      VALU1=VALU1+SSAF(IJ,IABC,JXYZ)*DP(IJ,KPQR,ITYP)
     1           +SSFG(IJ,JXYZ,KXYZ)*DP(IJ,IABC,ITYP)
     2           +SSAF(IJ,IABC,KXYZ)*DP(IJ,JPQR,ITYP)
      DO 102 K=1,NBASIS
      KTYP=MOTYP(K)
      VALU2=VALU2+UF(I,K,JXYZ)*UF(J,K,KXYZ)*DP(IJ,IABC,KTYP)
     1           +UF(I,K,KXYZ)*UA(J,K,IABC)*DP(IJ,JPQR,KTYP)
     2           +UA(I,K,IABC)*UF(J,K,JXYZ)*DP(IJ,KPQR,KTYP)
  102 CONTINUE
  103 CONTINUE
C
      VALU1=VALU1*TWO
      VALU2=-VALU2*FOUR
      VALUZ=VALU1+VALU2
      VALUX=AF1(JXYZ,KXYZ,IABC)
      VALUY=AF2(JXYZ,KXYZ,IABC)
      VALUT=VALUX+VALUY+VALUZ
      AFT(JXYZ,KXYZ,IABC)=VALUT
      IF(IPRNT.LE.2) GO TO 110
      WRITE(6,2) IABC,JXYZ,KXYZ,VALU1,VALU2,VALUZ,VALUX,VALUY,VALUT
C
  110 CONTINUE
  120 CONTINUE
  130 CONTINUE
C
      IF(IPRNT.LE.2) GO TO 202
      DO 140 IABC=1,N3N
      WRITE(6,3) IABC
      CALL MATOUT(AFT(1,1,IABC),3,3,3,3,6)
  140 CONTINUE
C
  202 CONTINUE
      RETURN
      END
