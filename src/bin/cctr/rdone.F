C
C                         *********************
C-------------------------******  RDONE  ******-----------------------
C                         *********************
C
      SUBROUTINE RDONE(S,H,T,V,LBLI,STVI,INBUFS,ITAP34,JOUT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION LBLI(INBUFS*2),STVI(INBUFS)
      COMMON/BASIS/NBASIS,MTRI,NSYM
      COMMON/SIGNS/IOFF(256),IPRNT
      DIMENSION S(MTRI),H(MTRI),T(MTRI),V(MTRI)
    1 FORMAT(//,2X,' S MATRIX'/)
    2 FORMAT(//,2X,' T MATRIX'/)
    3 FORMAT(//,2X,' V MATRIX'/)
    4 FORMAT(//,2X,' H MATRIX'/)
C
      IPRNT = 0
C   READ IN ONE ELECTRON INTEGRALS
C
C   S INTEGRALS
      CALL ONEI(S,LBLI,STVI,INBUFS,ITAP34)
      IF(IPRNT.LE.3) GO TO 201
      WRITE(JOUT,1)
      CALL PRINT(S,MTRI,NBASIS,JOUT)
C
C   T INTEGRALS
  201 CONTINUE
      CALL ONEI(T,LBLI,STVI,INBUFS,ITAP34)
      IF(IPRNT.LE.4) GO TO 202
      WRITE(JOUT,2)
      CALL PRINT(T,MTRI,NBASIS,JOUT)
C
C   V INTEGRALS
  202 CONTINUE
      CALL ONEI(V,LBLI,STVI,INBUFS,ITAP34)
      IF(IPRNT.LE.4) GO TO 203
      WRITE(JOUT,3)
      CALL PRINT(V,MTRI,NBASIS,JOUT)
C
C   TOTAL ONE ELECTRON INTEGRALS
  203 CONTINUE
      DO 101 I=1,MTRI
  101 H(I)=T(I)+V(I)
      IF(IPRNT.LE.3) GO TO 204
      WRITE(JOUT,4)
      CALL PRINT(H,MTRI,NBASIS,JOUT)
C
  204 CONTINUE
CGES
      ITAP67=67
      CALL RFILE(ITAP67)
      CALL SWRIT(ITAP67,S,INTOWP(MTRI))
      CALL SWRIT(ITAP67,H,INTOWP(MTRI))
      CALL RCLOSE (ITAP67,3)
CGES
C     DO 1000 I = 1,NBASIS
C     DO 1000 J=1,I
C     IJ=I*(I-1)/2+J
C1000 WRITE(14,123) I,J,S(IJ)
C 123 FORMAT(2I5,F20.10)
C     DO 2000 I = 1,NBASIS
C     DO 2000 J=1,I
C     IJ=I*(I-1)/2+J
C2000 WRITE(14,123) I,J,H(IJ)
C
      RETURN
      END
