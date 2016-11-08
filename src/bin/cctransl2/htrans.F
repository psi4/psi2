C                          ********************
C -------------------------****** HTRANS ******------------------------
C                          ********************
C
      SUBROUTINE HTRANS(X,XTP,X1,RINTS,TINTS,ITINTS,NDIM,LDIM,XDIM,
     .                  EIV,NBFAO,NBF,JOUT,NTRI1,NTRI2,IFLAG,ISORT,
     .                  ITAP78,NTRIMX,IBKT,ITRI,JBUF,NBKT,ORBSYM,FLOV,
     .                  NIRRED,NLAMDA,NSYMHF,KLPNT,IJPNT,NTRI,PTOCC,
     .                  IOFF1,IOFF2,NO,nfzo,nfzv)
C
C THIS SUBROUTINE PERFORMS A HALF TRANSFORMATION OF ATOMIC (SYMMETRY)
C ORBITALS TO MOLECULAR ORBITALS.  RINTS IS THE REAL*8 ARRAY CONTAINING
C THE INTEGRALS TO BE TRANSFORMED, TINTS IS THE REAL*8 ARRAY CONTAINING
C THE TRANSFORMED INTS, X & XTP ARE BUFFER ARRAYS, AND NDIM & LDIM ARE
C THE DIMENSIONS OF THE INTEGRAL ARRAYS.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER XDIM,FK,FL,FMU,FNU,FI,FJ,FFJ,nfzo(8),nfzv(8)
      REAL*8 RINTS(NDIM),TINTS(LDIM),XTP(XDIM,XDIM),X1(XDIM,XDIM),
     .          X(XDIM,XDIM),EIV(NBFAO,NBF)
      INTEGER ITINTS(2*LDIM),ORBSYM(NBF),FLOV(NIRRED,2),NLAMDA(NSYMHF),
     .        KLPNT(NTRI),IJPNT(NTRI),PTOCC(NBF),IOFF1(NTRI),
     .        IOFF2(NTRI)
      COMMON/PAKQ/I1,I2,I3
      COMMON/SIGNS/IOFF(256),IPRNT
      DATA ITEMP /255/
    1 FORMAT ('                                                ')
    2 FORMAT ('************************************************')
   29 FORMAT (2I3,F20.12)
   39 FORMAT (4I3,F20.12)
      KOUT=8
      R0=1.0D-15
C     WRITE(JOUT,*)
C     WRITE(JOUT,*) ' THE SCF EIGENVECTORS IN HTRANS '
C     WRITE(JOUT,*)
C     CALL MATOUT (EIV,NBFAO,NBF,NBFAO,NBF,JOUT)
C     WRITE (JOUT,*) ' LDIM=',LDIM
C     WRITE (JOUT,7) ISORT
C     WRITE (JOUT,17) IFLAG
      CALL ZERO(TINTS,LDIM)
C     WRITE (JOUT,7) ISORT
C     WRITE (JOUT,17) IFLAG
    7 FORMAT (' ISORT=',I2)
   17 FORMAT (' IFLAG=',I2)
   18 FORMAT (' ITINTS(',I4,')=',I8,'(',I3,',',I3,')  TINTS(',I4,')=',
     #           E20.10)
   19 FORMAT (2I5,F20.12)
   72 FORMAT (F20.12)
C     WRITE (JOUT,*) ' LDIM=',LDIM
      INTLEN=(INTOWP(LDIM)-2)/INTOWP(1)
      MAXINT=INTOWP(INTLEN)/(1+INTOWP(1))
C     WRITE (JOUT,*) ' MAXINT=',MAXINT
CTJL  IVOFF=(MAXINT+2)/INTOWP(1)
      IVOFF=(MAXINT+3)/INTOWP(1)
      ICNTR=1
C     WRITE (JOUT,*) ' ICNTR=',ICNTR
      ITRUF=ITRI*(IBKT-1)+1
      FI=IOFF1(ITRUF)
      FJ=IOFF2(ITRUF)
      ITRUL=ITRUF+NTRIMX-1
      LI=IOFF1(ITRUL)
      LJ=IOFF2(ITRUL)
C     INBF=NBFAO
C     IF (IFLAG.EQ.1) INBF=NBF
      JTRI=0
C     WRITE (JOUT,*) ' FI,LI',FI,LI
      DO 10 I=FI,LI
      ISYM=ORBSYM(I)
      FFJ=1
      LLJ=I
      IF (I.EQ.FI) FFJ=FJ
      IF (I.EQ.LI) LLJ=LJ
C     WRITE (JOUT,*) ' FFJ,LLJ',FFJ,LLJ
      DO 10 J=FFJ,LLJ
         JTRI=JTRI+1
C        WRITE (JOUT,*) ' I,J,JTRI',I,J,JTRI
C        IF (JTRI.GT.NTRIMX) GOTO 500
         JSYM=ORBSYM(J)
         IJSYM=IEOR(ISYM,JSYM)
         INDR=NTRI1*(JTRI-1)+1
C        WRITE (JOUT,*) ' JTRI=',JTRI
         CALL SQUR(RINTS(INDR),X,XDIM)
C        WRITE (JOUT,*) ' SQUARED RINTS INTO X'
         CALL ZERO(XTP,XDIM*XDIM)
C        WRITE (JOUT,*) ' ZEROED XTP'
         CALL ZERO(X1,XDIM*XDIM)
C        WRITE (JOUT,*) ' ZEROED X1'
C        IF (JTRI.EQ.1) THEN
C        WRITE (JOUT,*) ' ******************  X  ****************'
C        CALL MATOUT(X,XDIM,XDIM,XDIM,XDIM,JOUT)
C        WRITE (JOUT,*) ' ******************  XTP  ****************'
C        CALL MATOUT(XTP,XDIM,XDIM,XDIM,XDIM,JOUT)
C        WRITE (JOUT,*) ' ******************  EIV  ****************'
C        CALL MATOUT(EIV,NBFAO,NBF,NBFAO,NBF,JOUT)
C        END IF
C        WRITE (JOUT,*) ' NIRRED=',NIRRED
         DO 13 KIRRED=1,NIRRED
c           WRITE (JOUT,*) ' KIRRED=',KIRRED
            KSYM=KIRRED-1
            FK=FLOV(KSYM+1,1)
c           WRITE (JOUT,*) ' FK=',FK
            LK=FLOV(KSYM+1,2)
c           WRITE (JOUT,*) ' LK=',LK
            NUMK=LK-FK+1
C           WRITE (JOUT,*) ' NUMK=',NUMK
c           NUMMU=NUMK
            FMU=FK-nfzo(kirred)
c           WRITE (JOUT,*) ' FMU=',FMU
            LMU=LK+nfzv(kirred)
c           WRITE (JOUT,*) ' LMU=',LMU
            nummu=lmu-fmu+1
            LSYM=IEOR(KSYM,IJSYM)
            LIRRED=LSYM+1
C           WRITE (JOUT,*) ' LIRRED=',LIRRED
            FL=FLOV(LSYM+1,1)
c           WRITE (JOUT,*) ' FL=',FL
            LL=FLOV(LSYM+1,2)
c           WRITE (JOUT,*) ' LL=',LL
            NUML=LL-FL+1
C           WRITE (JOUT,*) ' NUML=',NUML
c           NUMNU=NUML
            FNU=FL-nfzo(lirred)
c           WRITE (JOUT,*) ' FNU=',FNU
            LNU=LL+nfzv(lirred)
c           WRITE (JOUT,*) ' LNU=',LNU
            numnu=lnu-fnu+1
            IF (FK.EQ.0.OR.FL.EQ.0) GOTO 13
            CALL ZERO(XTP,XDIM*XDIM)
            CALL MXMB(EIV(FMU,FK),NBFAO,1,X(FMU,FNU),1,
     .                XDIM,XTP(FK,FNU),1,XDIM,NUMK,NUMMU,NUMNU)
C           WRITE (JOUT,*) ' MATRIX MULTIPLIED EIV * X'
C        IF (JTRI.EQ.1) THEN
C        WRITE (JOUT,*) ' ******************  XTP  ****************'
C        CALL MATOUT(XTP,XDIM,XDIM,XDIM,XDIM,JOUT)
C        END IF
            CALL MXMB(XTP(FK,FNU),1,XDIM,EIV(FNU,FL),1,NBFAO,
     .                X1(FK,FL),1,XDIM,NUMK,NUMNU,NUML)
C           WRITE (JOUT,*) ' MATRIX MULTIPLIED XTP * EIV'
C        IF (JTRI.EQ.1) THEN
C        WRITE (JOUT,*) ' ******************  X1  ****************'
C        CALL MATOUT(X1,XDIM,XDIM,XDIM,XDIM,JOUT)
C        END IF
   13    CONTINUE
C        IF (JTRI.EQ.1) CALL MATOUT(X1,XDIM,XDIM,XDIM,XDIM,JOUT)
         KL=0
C        WRITE (JOUT,*) ' BEFORE WRITING INTS'
         IF (ISORT.EQ.1) THEN
         IF (IFLAG.EQ.2) THEN
         DO 20 K=1,NBF
            DO 30 L=1,K
               KL=KL+1
               INDT2=NTRI2*(KL-1)+JTRI
               TINTS(INDT2)=X1(K,L)
   30       CONTINUE
   20    CONTINUE
         ELSE
C        WRITE (JOUT,*) 'IJ=',JTRI
         DO 40 K=1,NBF
            DO 50 L=1,K
               KL=KL+1
               IF (DABS(X1(K,L)).GT.RO) THEN
               IF (ICNTR.GT.MAXINT) THEN
                  CALL WRINTS(ITAP78,TINTS,LDIM,ISCTR)
                  JBUF=JBUF+1
                  CALL ZERO (TINTS,LDIM)
                  ICNTR=1
C                 WRITE (JOUT,*) ' ICNTR=',ICNTR
               END IF
               KCC=PTOCC(K)
               LCC=PTOCC(L)
               ITMPX=MAX0(KCC,LCC)
               ITMPN=MIN0(KCC,LCC)
               KCC=ITMPX
               LCC=ITMPN
               LBL=IOR(LCC,ISHFT(KCC,8))
C              INDT1=NTRI2*(JTRI-1)+KL
               ITINTS(1)=0
               ITINTS(2)=ICNTR
               ITINTS(ICNTR+2)=LBL
               TINTS(IVOFF+ICNTR)=X1(K,L)
               ICNTR=ICNTR+1
C              WRITE (JOUT,*) ' ICNTR=',ICNTR
            END IF
   50       CONTINUE
   40    CONTINUE
         ITINTS(1)=-1
         CALL WRINTS(ITAP78,TINTS,LDIM,ISCTR)
C        CALL PRINT(TINTS(INDR),NTRI2,NBF,JOUT)
         END IF
         ELSE
         IF (IFLAG.EQ.2) THEN
C        ILOOP=1
         DO 22 K=1,NBF
            DO 32 L=1,K
               KL=KL+1
               IF (DABS(X1(K,L)).GT.R0) THEN
               IF (ICNTR.GT.MAXINT) THEN
                  CALL RGETSA(ITAP78,IAD)
C                 WRITE (JOUT,1)
C                 WRITE (JOUT,2)
C                 WRITE (JOUT,*)' ******   TAPE78 DUMP   ******'
C                 WRITE (JOUT,2)
C                 WRITE (JOUT,1)
C                 IADHLD=IAD
C                 WRITE(JOUT,*) ' IAD=',IAD
                  CALL WRINTS (ITAP78,TINTS,LDIM,ISECTR)
                  JBUF=JBUF+1
                  CALL ZERO (TINTS,LDIM)
C                 WRITE (JOUT,1)
C                 WRITE (JOUT,2)
C                 WRITE (JOUT,*)' ******   TINTS DUMP   ******'
C                 WRITE (JOUT,2)
C                 WRITE (JOUT,1)
C                 CALL RSETSA(ITAP78,IADHLD)
C                 WRITE (JOUT,*) ' IADHLD=',IADHLD
C                 CALL SREAD(ITAP78,TINTS,INTOWP(LDIM))
C                 DO 53 JQ=1,MAXINT
C                    I3=ITINTS(JQ)
C                    CALL UNPAK
C                    QINT=TINTS(JQ+IVOFF)
C                    WRITE (JOUT,18) JQ,I3,I1,I2,JQ+IVOFF,QINT
C  53             CONTINUE
                  ICNTR=1
C                 WRITE (JOUT,*) ' ICNTR=',ICNTR
C                 ILOOP=ILOOP+1
               END IF
               I1=KL
               I2=JTRI+ITRI*(IBKT-1)
               CALL PAK
               ITINTS(1)=0
               ITINTS(2)=ICNTR
               ITINTS(ICNTR+2)=I3
               TINTS(IVOFF+ICNTR)=X1(K,L)
               ICNTR=ICNTR+1
C              WRITE (JOUT,*) ' ICNTR=',ICNTR
C              WRITE (JOUT,18) ICNTR,ITINTS(ICNTR+2),I1,I2,IVOFF+ICNTR,
C    .                TINTS(IVOFF+ICNTR)
               END IF
   32       CONTINUE
   22    CONTINUE
         ELSE
C        WRITE (JOUT,*) ' BEFORE 42'
         DO 42 K=1,NBF
            DO 52 L=1,K
C              WRITE (JOUT,*) ' K,L=',K,L
               KL=KL+1
               IF (DABS(X1(K,L)).GT.R0) THEN
               IF (ICNTR.GT.MAXINT) THEN
C                 CALL RGETSA(ITAP78,IAD)
C                 WRITE (JOUT,*) ' IAD=',IAD
C                 WRITE (JOUT,1)
C                 WRITE (JOUT,2)
C                 WRITE (JOUT,*)' ******   TINTS DUMP   ******'
C                 WRITE (JOUT,2)
C                 WRITE (JOUT,1)
C                 DO 63 JQ=1,MAXINT
C                    I3=ITINTS(JQ+2)
C                    CALL UNPAK
C                    PINT=TINTS(JQ+IVOFF)
C                    WRITE (JOUT,18) JQ,I3,I1,I2,JQ+IVOFF,PINT
C  63             CONTINUE
                  CALL WRINTS (ITAP78,TINTS,LDIM,ISECTR)
                  JBUF=JBUF+1
                  CALL ZERO (TINTS,LDIM)
                  ICNTR=1
C                 WRITE (JOUT,*) ' ICNTR=',ICNTR
               END IF
               I1=KL
               I2=JTRI+ITRI*(IBKT-1)
               IF (I1.GE.I2) THEN
                  ICC=PTOCC(I)
                  JCC=PTOCC(J)
                  KCC=PTOCC(K)
                  LCC=PTOCC(L)
                  ITMPX=MAX0(ICC,JCC)
                  ITMPN=MIN0(ICC,JCC)
                  ICC=ITMPX
                  JCC=ITMPN
                  ITMPX=MAX0(KCC,LCC)
                  ITMPN=MIN0(KCC,LCC)
                  KCC=ITMPX
                  LCC=ITMPN
                  IJCC=IOFF(ICC)+JCC
                  KLCC=IOFF(KCC)+LCC
                  ITMPX=MAX0(IJCC,KLCC)
                  ITMPN=MIN0(IJCC,KLCC)
                  IJCC=ITMPX
                  KLCC=ITMPN
                  ICC=IOFF1(IJCC)
                  JCC=IOFF2(IJCC)
                  KCC=IOFF1(KLCC)
                  LCC=IOFF2(KLCC)
                  IJKL=IOR(JCC,ISHFT(ICC,8))
                  IJKL=IOR(KCC,ISHFT(IJKL,8))
                  IJKL=IOR(LCC,ISHFT(IJKL,8))
C                 CALL PAK
                  ITINTS(1)=0
                  ITINTS(2)=ICNTR
                  ITINTS(ICNTR+2)=IJKL
                  TINTS(IVOFF+ICNTR)=X1(K,L)
C                 IF (KCC.GT.NO.AND.JCC.LE.NO.AND.LCC.LE.NO)
C    .            WRITE (JOUT,39) ICC,JCC,KCC,LCC,TINTS(IVOFF+ICNTR)
C                 WRITE (JOUT,29) I1,I2,TINTS(IVOFF+ICNTR)
                  ICNTR=ICNTR+1
C                 WRITE (JOUT,*) ' ICNTR=',ICNTR
               END IF
               END IF
   52       CONTINUE
   42    CONTINUE
         END IF
         END IF
   10 CONTINUE
  500 IF (ISORT.EQ.2) THEN
         IF ((IBKT.EQ.NBKT).AND.(IFLAG.EQ.1)) ITINTS(1)=-1
         CALL WRINTS (ITAP78,TINTS,LDIM,ISCTR)
C        WRITE (JOUT,*) ' ICNTR=',ITINTS(2)
         JBUF=JBUF+1
C        IF ((IBKT.EQ.NBKT).AND.(IFLAG.EQ.1)) THEN
C        CALL SREW(ITAP78)
C        DO 33 IQ=1,JBUF
C           WRITE (JOUT,*) ' JBUF=',IQ
C           CALL ZERO(TINTS,LDIM)
C           CALL SREAD(ITAP78,TINTS,INTOWP(LDIM))
C           WRITE (JOUT,*) ' IFLG=',ITINTS(1)
C           CALL SETMBF(MBUF,ITINTS(2))
C           DO 43 JQ=1,MBUF
C              IJKL=ITINTS(JQ+2)
C              I=ISHFT(IJKL,-24)
C              J=IAND(ITEMP,ISHFT(IJKL,-16))
C              K=IAND(ITEMP,ISHFT(IJKL,-8))
C              L=IAND(ITEMP,IJKL)
C              QINT=TINTS(IVOFF+JQ)
C              WRITE (JOUT,*) I,J,K,L,QINT
C  43      CONTINUE
C  33    CONTINUE
C     END IF
      END IF
      RETURN
      END
