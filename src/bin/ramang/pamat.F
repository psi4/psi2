      SUBROUTINE PAMAT(UA,UF,DP,DQ,ESO,NABC)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DP(NTRI,NABC,NTYPEP),DQ(NTRI,NABC,NTYPEP)
      DIMENSION UA(NBASIS,NBASIS,N3N),UF(NBASIS,NBASIS,3)
      DIMENSION ESO(NBFAO,NBASIS)
      COMMON/BASIS/NBASIS,NBFAO,NTRI,NTRI2,NBATRI,NBASQ
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3
      COMMON/GVBAB/A(10,10),B(10,10)
      COMMON/SIGNS/IOFF(466),IPRNT
      COMMON/START/FOCC(10),NSORB(10),NSTR(10),NEND(10),MOTYP(256)
      DATA ZERO / 0.0D+00 /
    1 FORMAT(//,2X,' DP MATRIX'/)
    2 FORMAT(//,2X,' DQ MATRIX'/)
C
      DO 101 IABC=1,NABC
      DO 101 ITYP=1,NTYPEP
      DO 101 I=1,NTRI
      DP(I,IABC,ITYP)=ZERO
      DQ(I,IABC,ITYP)=ZERO
  101 CONTINUE
C
      IJ=0
      DO 106 I=1,NBASIS
      DO 106 J=1,I
      IJ=IJ+1
      DO 105 K=1,NBASIS
      DO 105 L=1,NBASIS
      KTYP=MOTYP(K)
      LTYP=MOTYP(L)
      FAC=ESO(I,K)*ESO(J,L)+ESO(I,L)*ESO(J,K)
      DO 104 ITYP=1,NTYPES
      AVAL=A(ITYP,LTYP)
      BVAL=B(ITYP,LTYP)
      FACA=AVAL*FAC
      FACB=BVAL*FAC
      DO 102 IABC=1,N3N
      DP(IJ,IABC,ITYP)=DP(IJ,IABC,ITYP)+UA(K,L,IABC)*FACA
      DQ(IJ,IABC,ITYP)=DQ(IJ,IABC,ITYP)+UA(K,L,IABC)*FACB
  102 CONTINUE
      DO 103 IXYZ=1,3
      IPQR=IXYZ+N3N
      DP(IJ,IPQR,ITYP)=DP(IJ,IPQR,ITYP)+UF(K,L,IXYZ)*FACA
      DQ(IJ,IPQR,ITYP)=DQ(IJ,IPQR,ITYP)+UF(K,L,IXYZ)*FACB
  103 CONTINUE
  104 CONTINUE
  105 CONTINUE
  106 CONTINUE
      IF(IPRNT.LE.4) GO TO 201
      WRITE(6,1)
      CALL MATOUT(DP,NTRI,NABC*NTYPEP,NTRI,NABC*NTYPES,6)
      WRITE(6,2)
      CALL MATOUT(DQ,NTRI,NABC*NTYPEP,NTRI,NABC*NTYPES,6)
C
  201 CONTINUE
      RETURN
      END
