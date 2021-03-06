      SUBROUTINE RDONE(S,H,T,V,LBLI,STVI,INBUFS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION S(NTRI),H(NTRI),T(NTRI),V(NTRI)
      DIMENSION LBLI(INBUFS*2),STVI(INBUFS)
      COMMON/BASIS/NBASIS,NTRI,NST,NSYM,NBFAO,NBFSO,NBATRI
      COMMON/SIGNS/IOFF(256),IPRNT
    1 FORMAT(//,2X,' S MATRIX'/)
    2 FORMAT(//,2X,' T MATRIX'/)
    3 FORMAT(//,2X,' V MATRIX'/)
    4 FORMAT(//,2X,' H MATRIX'/)
C
C   READ IN ONE ELECTRON INTEGRALS
C
C***  ITAPE8=8
C***  REWIND ITAPE8
C   S INTEGRALS
      CALL ONEI(S,LBLI,STVI,INBUFS)
      IF(IPRNT.LE.3) GO TO 201
      WRITE(6,1)
      CALL PRINT(S,NTRI,NBASIS,6)
C
C   T INTEGRALS
  201 CONTINUE
      CALL ONEI(T,LBLI,STVI,INBUFS)
      IF(IPRNT.LE.4) GO TO 202
      WRITE(6,2)
      CALL PRINT(T,NTRI,NBASIS,6)
C
C   V INTEGRALS
  202 CONTINUE
      CALL ONEI(V,LBLI,STVI,INBUFS)
      IF(IPRNT.LE.4) GO TO 203
      WRITE(6,3)
      CALL PRINT(V,NTRI,NBASIS,6)
C
C   TOTAL ONE ELECTRON INTEGRALS
  203 CONTINUE
      DO 101 I=1,NTRI
  101 H(I)=T(I)+V(I)
      IF(IPRNT.LE.3) GO TO 204
      WRITE(6,4)
      CALL PRINT(H,NTRI,NBASIS,6)
C
  204 CONTINUE
CCC   WRITE(6,*) ' ONE ELECTRON INTEGRALS---S , T, V, H'
CCC   DO 102 I=1,NTRI
CCC   WRITE(6,5) I,S(I),T(I),V(I),H(I)
CC102 CONTINUE
CCC 5 FORMAT(2X,I5,5X,4F20.10)
C
      RETURN
      END
