      SUBROUTINE MAGCHK(EIGVAL,NVAR)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 EIGVAL(1)
      LOGICAL PRNT,NRFLAG,DUMP
      COMMON/IO/IN,IOUT
      COMMON/OPTEF/
     $             ES,DMAX,DD,CONVF,CNVFMX,CONVX,CNVXMX,
     $             EIGMAX,EIGMIN,vmode(50),
     $             X(50),XNAME(50),OLDF(50),D(50),HESS(50,50),
     $             IC(50),MODE,NSTEP,MXSTEP,IHESS,IS,
     $             NEGREQ,IUPDAT,
     $             MXHESS,PRNT,NRFLAG,DUMP,IDUM
C
      DATA ZERO/0.0D0/
C
C  FIRST CHECK IF THE MAGNITUDES OF THE LOWEST EIGENVALUES
C  EXCEED EIGMIN
C
      DO 10 I=1,NVAR
      IF(EIGVAL(I).GT.EIGMIN) GO TO 11
      IF(DABS(EIGVAL(I)).LT.EIGMIN) THEN
      IF(EIGVAL(I).LT.ZERO) EIGVAL(I)=-EIGMIN
      IF(EIGVAL(I).GE.ZERO) EIGVAL(I)= EIGMIN
      IF(PRNT) WRITE(IOUT,1000) I,EIGVAL(I)
      ENDIF
   10 CONTINUE
   11 CONTINUE
C
C  NOW CHECK IF THE MAGNITUDES OF THE HIGHEST EIGENVALUES
C  EXCEED EIGMAX
C
      DO 20 I=1,NVAR
      II=NVAR-I+1
      IF(EIGVAL(II).LT.EIGMAX) GO TO 21
      IF(DABS(EIGVAL(II)).GT.EIGMAX) THEN
      IF(EIGVAL(II).LT.ZERO) EIGVAL(II)=-EIGMAX
      IF(EIGVAL(II).GE.ZERO) EIGVAL(II)= EIGMAX
      IF(PRNT) WRITE(IOUT,2000) II,EIGVAL(II)
      ENDIF
   20 CONTINUE
   21 CONTINUE
C
      RETURN
C
1000  FORMAT(' WARNING! EIGENVALUE ',I2,' TOO SMALL. REPLACED',
     $       ' BY ',F12.6)
C
2000  FORMAT(' WARNING! EIGENVALUE ',I2,' TOO LARGE. REPLACED',
     $       ' BY ',F12.6)
C
      END
