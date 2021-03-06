      SUBROUTINE mbmat2(b,d2,m2,nn2,iset2)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION D(6),D2(6),DC(6),M(6),M2(6),MC(6)
      COMMON/BASIS/NBASIS,NBFAO,NBATRI,NBATR2
      COMMON/CENTR/NCENTR
CTPH  COMMON/DERIN/DERINT(1296,3,3)
C11-14-88  (For  7Mb) COMMON/DERIN/DERINT(10000,9)
C11-14-88  (For 15Mb) COMMON/DERIN/DERINT(50625,9)
      COMMON/DERIN/DERINT(10000,3,3)
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOM,N3N,NATRI,NSYM
      COMMON/PARA3/NATMS(4),NPASS
      COMMON/TOLER/DLIM
      DIMENSION B(NBATRI,3,NATOM,NTREAD)
      DATA ZERO / 0.0D+00 /
C
cges  ENTRY MBMAT2(B,D2,M2,NN2,ISET2)
C
      DO 108 IXYZ=1,3
      TEMP=ZERO
      DO 106 I=1,NCENTR-1
      T=DERINT(NN2,IXYZ,I)
      TEMP=TEMP+T
      IATOM=NATMS(I)
      DO 105 J=3,6
      B(M2(J),IXYZ,IATOM,ISET2)=B(M2(J),IXYZ,IATOM,ISET2)+T*D2(J)
  105 CONTINUE
  106 CONTINUE
C
C   USE TRANSLATIONAL INVARIANCE FOR LAST ATOM
      IATOM=NATMS(NCENTR)
      DO 107 J=3,6
      B(M2(J),IXYZ,IATOM,ISET2)=B(M2(J),IXYZ,IATOM,ISET2)-TEMP*D2(J)
  107 CONTINUE
C
  108 CONTINUE
C
      RETURN
      end
