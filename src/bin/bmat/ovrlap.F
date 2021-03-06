      SUBROUTINE OVRLAP(U,NMODE,NVAR)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 U(50,50)
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
C  ON THE FIRST STEP DETERMINE WHICH MODE TO FOLLOW
C
      IF(NSTEP.EQ.1) THEN
C
C  (A) FOLLOWING A GIVEN MODE
C
      IF(MODE.GT.NVAR) THEN
      IT=MODE-NVAR
      GO TO 10
      ENDIF
C
C  FIND THE HESSIAN EIGENVECTOR WITH THE LARGEST MAGNITUDE
C  COMPONENT IN THE POSITION GIVEN BY MODE (READ IN VIA IC ARRAY)
C
      IT=1
      TOVLP=DABS(U(MODE,1))
      DO 5 I=2,NVAR
      IF(DABS(U(MODE,I)).GT.TOVLP) THEN
      TOVLP=DABS(U(MODE,I))
      IT=I
      ENDIF
    5 CONTINUE
C
   10 MODE=IT
      IF(PRNT) WRITE(IOUT,1000) MODE
C
      ELSE
C
C  ON SUBSEQUENT STEPS DETERMINE WHICH HESSIAN EIGENVECTOR HAS
C  THE GREATEST OVERLAP WITH THE MODE WE ARE FOLLOWING
C
      IT=1
      TOVLP=DVTV(U(1,1),VMODE,NVAR)
      TOVLP=DABS(TOVLP)
      DO 6 I=2,NVAR
      OVLP=DVTV(U(1,I),VMODE,NVAR)
      OVLP=DABS(OVLP)
      IF(OVLP.GT.TOVLP) THEN
      TOVLP=OVLP
      IT=I
      ENDIF
    6 CONTINUE
C
      IF(PRNT) WRITE(IOUT,1010) TOVLP
      ENDIF
C
C  SAVE THE EIGENVECTOR IN VMODE
C
      DO 7 I=1,NVAR
      VMODE(I)=U(I,IT)
    7 CONTINUE
C
      NMODE=IT
      RETURN
C
1000  FORMAT(' HESSIAN MODE FOLLOWING SWITCHED ON'/
     $       ' FOLLOWING MODE ',I3)
C
1010  FORMAT(' OVERLAP OF CURRENT MODE WITH PREVIOUS MODE IS ',F12.6)
C
      END
