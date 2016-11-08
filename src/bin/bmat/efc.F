C**********************************************************************
C     NOTIFICATION OF PROGRAM CHANGE:                                 *
C     BY:     TRACY HAMILTON                                          *
C     DATE:    JULY 18,1988                                           *
C     REASON:  SCALE FORCES CORRECTLY IN DO 20 LOOP                   *
C**********************************************************************
      SUBROUTINE EFC(NUMIT,NQ,F1,QQ,FI,C,MODEBM,QQ1,IFSTRE,GEC,vmode1)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C  THIS ROUTINE SETS UP ALL THE COMMON BLOCKS NEEDED PRIOR TO CALLING
C  LINK 110.
C
      real*8 vmode1(nq)
      LOGICAL PRNT,NRFLAG,DUMP
      COMMON/IOP/IOP(50)
      COMMON/IO/IN,IOUT,INP2,IPUN,IGMUP
      COMMON/GRDNT/ENERGY,F(50),FRCNST(1275),NVAR,IGETFC
      COMMON/OPTEF/
     $             ES,DMAX,DD,CONVF,CNVFMX,CONVX,CNVXMX,
     $             EIGMAX,EIGMIN,vmode(50),
     $             X(50),XNAME(50),OLDF(50),D(50),HESS(50,50),
     $             IC(50),MODE,NSTEP,MXSTEP,IHESS,IS,
     $             NEGREQ,IUPDAT,
     $             MXHESS,PRNT,NRFLAG,DUMP,IDUM
      COMMON/ASTRE/STRE(50)
      DIMENSION F1(1),QQ(1),FI(1),C(1),QQ1(1),CORE(1),IFSTRE(1),GEC(1)
      LOGICAL STRE,IFSTRE
      DATA TOANG/0.52917706D0/, HARTRE/4.359813653D0/
      CONV1 = TOANG/HARTRE
      CONV2 = TOANG*CONV1
      DO 5 I=1,NQ
    5 STRE(I) = IFSTRE(I)
C
      NSTEP = NUMIT
      MODE = MODEBM
C
C     NUMBER OF SYMMETRY INTERNAL COORDINATES:
      NVAR = NQ
C
C     VALUES OF SYMMETRY INTERNAL COORDINATES:
C     CONVERT  ANG --> BOHR
c
      DO 10 I=1,NVAR
        X(I) = QQ(I)
        IF(STRE(I)) X(I) = X(I) / TOANG
   10 CONTINUE
 
C     FORCES (CONVERT  MDYNE --> HARTREE/BOHR):

      DO 20 I=1,NVAR
         IF(STRE(I)) THEN
            F(I) = FI(I) * CONV1
         ELSE
            F(I) = FI(I) / HARTRE
         ENDIF
   20 CONTINUE
 
C     HESSIAN:
C              MDYNE/ANG --> HARTREE/BOHR**2
C              MDYNE/RAD --> HARTREE/BOHR/RAD
C              MDYNE*ANG/RAD**2 --> HARTREE/RAD**2

      IJ = 0
      DO 31 I=1,NVAR
        DO 30 J=1,NVAR
          IJ = IJ+1
          HESS(I,J) = C(IJ)
          IF(STRE(I).AND.STRE(J)) THEN
            HESS(I,J) = HESS(I,J) * CONV2
          ELSE IF(STRE(I).OR.STRE(J)) THEN
            HESS(I,J) = HESS(I,J) * CONV1
          ELSE
            HESS(I,J) = HESS(I,J) / HARTRE
          ENDIF
   30   CONTINUE
   31 CONTINUE
C
      IF(NSTEP.GT.0) THEN
C
C     OLD GRADIENT (NEGATIVE OF OLD FORCES) (CONVERT TO HARTREE/BOHR):
        DO 40 I=1,NVAR
   40   OLDF(I) = -F1(I) * CONV1

        DD = 0.0D0
        DO 45 I=1,NVAR
           D(I) = GEC(I)
           IF(STRE(I)) D(I) = D(I) / TOANG
           DD = DD + D(I)**2
   45   CONTINUE
        DD = DSQRT(DD)
 
C     POSSIBLY RECOVER OLD HESSIAN MODE:
        IF(MODE.GT.0) THEN
          DO 46 I=1,NVAR
             READ(INP2,1020) VMODE(I)
   46     CONTINUE
          IF(MODE.EQ.1) MODE=0
        ENDIF
      endif
 
      CALL OPTEFC(CORE)
 
C     PUT HESSIAN BACK INTO BMAT FORM:
      MODEBM = MODE
      DO 51 I=1,NVAR
         VMODE1(I) = VMODE(I)
  51  CONTINUE

C         HARTREE/BOHR**2 --> MDYNE/ANG
C         HARTREE/BOHR/RAD --> MDYNE/RAD
C         HARTREE/RAD**2 --> MDYNE*ANG/RAD**2
      IJ = 0
      DO 70 I=1,NVAR
        DO 70 J=1,NVAR
          IJ = IJ+1
          IF(STRE(I).AND.STRE(J)) THEN
            HESS(I,J) = HESS(I,J) / CONV2
          ELSE IF(STRE(I).OR.STRE(J)) THEN
            HESS(I,J) = HESS(I,J) / CONV1
          ELSE
            HESS(I,J) = HESS(I,J) * HARTRE
        ENDIF
   70 C(IJ) = HESS(I,J)
C
C     COPY DISPLACEMENTS INTO BMAT ARRAY (AND CONVERT BOHR --> ANG):
      DO 80 I=1,NVAR
      QQ1(I) = D(I)
      IF(STRE(I)) QQ1(I) = QQ1(I) * TOANG
   80 CONTINUE
C
      RETURN
C
 1000 FORMAT(F10.6)
 1010 FORMAT(2I5)
 1020 FORMAT(8F10.6)
C
      END
