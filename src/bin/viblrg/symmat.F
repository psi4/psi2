      SUBROUTINE SYMMAT(SMAT,K,TYPE,INUM,NSYM1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SMAT(3,3)
      INTEGER TYPE,SYM
      REAL*8 MDYNE
       COMMON /PARAM/VECT(60),ANGS,FRCS,MDYNE
     1,SYM(60),ICEN,N,ISYM,ITYPE(20),IATNO(20),ITAPE5
      DIMENSION TMAT(3,3),TMATT(3,3),PROD(3,3)
      ONE=1.0D+00
      ZERO=0.0D+00
      TWO=2.0D+00
      IF(TYPE.GT.1) GO TO 20
      DO 10 I=1,3
      DO 10 J=1,3
   10   SMAT(I,J)=ZERO
      SMAT(1,1)=ONE
      SMAT(2,2)=ONE
      SMAT(3,3)=-(ONE)
      GO TO 30
   20 PI=3.141592653589793D+00
      THETA=TWO*PI/INUM
      SMAT(1,1)=COS(THETA)
      SMAT(2,1)=SIN(THETA)
      SMAT(3,1)=ZERO
      SMAT(1,2)=-SIN(THETA)
      SMAT(2,2)=COS(THETA)
      SMAT(3,2)=ZERO
      SMAT(1,3)=ZERO
      SMAT(2,3)=ZERO
      SMAT(3,3)=ONE
      IF(TYPE.GT.2) SMAT(3,3)=-(ONE)
   30 NSYMS=3*(NSYM1-1)+1
      SUM=ZERO
      DO 33 I=0,2
        SUM=SUM+VECT(NSYMS+I)**2
   33 CONTINUE
      IF(SUM.NE.ZERO) GO TO 35
      WRITE(6,34) SYM(K+1)
   34 FORMAT( 'NORM OF ',I2,'TH SYMMETRY VECTOR IS ZERO')
      ISTOP=0
      IF(ISTOP.EQ.0)STOP
   35 DO 37 I=0,2
        VECT(NSYMS+I)=VECT(NSYMS+I)/DSQRT(SUM)
   37 CONTINUE
      COSPHI=VECT(NSYMS+2)
      IF(COSPHI.NE.ONE) GO TO 45
      DO 40 I=1,3
      DO 40 J=1,3
        TMAT(I,J)=ZERO
        IF(I.EQ.J) TMAT(I,J)=ONE
   40   TMATT(J,I)=TMAT(I,J)
      WRITE(6,1001)
1001  FORMAT(' REACHED 40')
      GO TO 60
   45 SINPHI=DSQRT(ONE-COSPHI*COSPHI)
      COSTHT=VECT(NSYMS)/SINPHI
      SINTHT=VECT(NSYMS+1)/SINPHI
      TMAT(1,3)=-VECT(NSYMS)
      TMAT(2,3)=-VECT(NSYMS+1)
      TMAT(3,1)=VECT(NSYMS)
      TMAT(3,2)=VECT(NSYMS+1)
      TMAT(3,3)=COSPHI
      TMAT(1,1)=COSTHT*COSTHT*COSPHI+SINTHT*SINTHT
      TMAT(2,2)=SINTHT*SINTHT*COSPHI+COSTHT*COSTHT
      TMAT(2,1)=COSTHT*SINTHT*(COSPHI-ONE)
      TMAT(1,2)=TMAT(2,1)
      DO 50 I=1,3
      DO 50 J=1,3
   50   TMATT(I,J)=TMAT(J,I)
   60 CALL MMULT(SMAT,TMAT,PROD,3,3,3)
      CALL MMULT(TMATT,PROD,SMAT,3,3,3)
      K=K+TYPE+1
      RETURN
      END
