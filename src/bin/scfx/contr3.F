C     ///////////////////////////////////////////////////////
      SUBROUTINE CONTR3 (VALUE,BUFFER,DT,DQ,DR,DAB,P,DQP,DRP,DQQ,DRQ,
     1DQRP,NMAT)
C
C
C
      IMPLICIT INTEGER (A-Z)
C
C
      REAL*8 VALUE(1),DT(1),DQ(1),DR(1),DAB(1),P(1),DQP(1),
     1DRP(1),DQQ(1),DRQ(1),DQRP(1)
      REAL*8 ALPHA1,ALPHA2,ALPHA3,QTEMP
      INTEGER PRINT
      INTEGER START,END
      REAL*8 RP1
      LOGICAL TWOCON
      REAL*8 AVG
      DIMENSION BUFFER(1)
      COMMON /timex/ITIME,ICONT0,JCOUNT,ITIME2,ICONT2,ITIME3,ICONT3,
     1ITIME4,ICONT4,ITIME1,ICONT1
      COMMON /RROB/ BTAG,POINTR,TAGMAX,RDCNT
      COMMON /T102/ NOUT,HALF,THIRD,MAXNO,OUTBUF,KINADR
      COMMON /ABOPEN/ NOST
      COMMON /PRINT/ PRINT
      COMMON /GVB/ ALPHA1,ALPHA2,ALPHA3,IJU,DUMMY
      COMMON /TAPES/ NIN,IR,IOUT,ICHECK,IFLAG3,ITAPE,ITERML
C
      ICONT4=ICONT4+1
      LIM=MIN(MAXNO,TAGMAX)
      NUMBUF=0
      NOREAD=MAXNO-TAGMAX
      DO 2 I=1,NMAT
        P(I)=0.0D+00
        DRP(I)=0.0D+00
        DQP(I)=0.0D+00
        DQQ(I)=0.0D+00
        DRQ(I)=0.0D+00
        DQRP(I)=0.0D+00
    2 CONTINUE
C
      DO 8 NUMBUF=1,MAXNO
        CASE=BUFFER(POINTR)
        NOINBF=BUFFER(POINTR+1)
        IF (CASE.EQ.4) THEN
C     CALL CASE4(BUFFER(POINTR+2),DT,P,NOINBF)
          LOCAT=(POINTR+1)/2+HALF+1
C
          CALL CASE4(VALUE(LOCAT),BUFFER(POINTR+2),DT,P,NOINBF)
          CALL CASE4(VALUE(LOCAT),BUFFER(POINTR+2),DQ,DQP,NOINBF)
          CALL CASE4(VALUE(LOCAT),BUFFER(POINTR+2),DR,DRP,NOINBF)
         CALL CASE4(VALUE(LOCAT),BUFFER(POINTR+2),DAB,DQRP,NOINBF)
        ELSE IF (CASE.EQ.3) THEN
          N3I=POINTR+1
          N3RB=N3I/2+THIRD
          N3RC=N3RB+THIRD
          DO 4 JUNK=1,NOINBF
            N3I=N3I+2
            N3RB=N3RB+1
            N3RC=N3RC+1
            KL=BUFFER(N3I-1)
            IJ=BUFFER(N3I)
            DQRP(IJ)=DQRP(IJ)+VALUE(N3RB)*DAB(KL)
            DQRP(KL)=DQRP(KL)+VALUE(N3RB)*DAB(IJ)
            P(IJ)=P(IJ)+VALUE(N3RB)*DT(KL)
            P(KL)=P(KL)+VALUE(N3RB)*DT(IJ)
            DQP(IJ)=DQP(IJ)+VALUE(N3RB)*DQ(KL)
            DQP(KL)=DQP(KL)+VALUE(N3RB)*DQ(IJ)
            DRP(IJ)=DRP(IJ)+VALUE(N3RB)*DR(KL)
            DRP(KL)=DRP(KL)+VALUE(N3RB)*DR(IJ)
            DQQ(IJ)=DQQ(IJ)+VALUE(N3RC)*DQ(KL)
            DQQ(KL)=DQQ(KL)+VALUE(N3RC)*DQ(IJ)
            DRQ(IJ)=DRQ(IJ)+VALUE(N3RC)*DR(KL)
            DRQ(KL)=DRQ(KL)+VALUE(N3RC)*DR(IJ)
C     IF ((PRINT.AND.128).NE.0) PRINT '(2(' P(',I4,')=',F10.7,' DT(',I4,
C    1                    ')=',F10.7),/,2(' Q(',I4,'2)=',F10.7,' DO(',I4
C    *                    ,')=',F10.7))', IJ, P(IJ),IJ,DT(IJ),KL,P(KL),K
C    *L,DT(KL),IJ,Q(IJ),IJ,DO(IJ),KL,Q(KL),   KL,DO(KL)
    4     CONTINUE
        ELSE IF (CASE.EQ.2) THEN
          N2I=POINTR+1
          N2R=N2I/2+HALF
          DO 5 JUNK=1,NOINBF
            N2I=N2I+1
            N2R=N2R+1
            IJ=BUFFER(N2I)
            P(IJ)=P(IJ)+VALUE(N2R)*DT(IJ)
            DQRP(IJ)=DQRP(IJ)+VALUE(N2R)*DAB(IJ)
            DQP(IJ)=DQP(IJ)+VALUE(N2R)*DQ(IJ)
            DRP(IJ)=DRP(IJ)+VALUE(N2R)*DR(IJ)
C     IF ((PRINT.AND.128).NE.0) PRINT '(' P(',I4,')=',F10.7,' DT(',I4,
C    1               ')=',F10.7)', IJ,P(IJ),IJ,DT(IJ)
    5     CONTINUE
        ELSE IF (CASE.EQ.1) THEN
          N1I=POINTR+1
          N1RA=N1I/2+THIRD
          N1RB=N1RA+THIRD
          DO 7 JUNK=1,NOINBF
            N1I=N1I+1
            N1RA=N1RA+1
            N1RB=N1RB+1
            IJ=BUFFER(N1I)
            P(IJ)=P(IJ)+VALUE(N1RA)*DT(IJ)
            DQRP(IJ)=DQRP(IJ)+VALUE(N1RA)*DAB(IJ)
            DQP(IJ)=DQP(IJ)+VALUE(N1RA)*DQ(IJ)
            DRP(IJ)=DRP(IJ)+VALUE(N1RA)*DR(IJ)
            DQQ(IJ)=DQQ(IJ)+VALUE(N1RB)*DQ(IJ)
            DRQ(IJ)=DRQ(IJ)+VALUE(N1RB)*DR(IJ)
C     IF ((PRINT.AND.128).NE.0) PRINT '(' P(',I4,')=',F10.7,' DT(',I4,
C    1                    ')=',F10.7,' Q(',I4,')=',F10.7,' DO(',I4,')=',
C    *                    F10.7)', IJ,P(IJ),  IJ,DT(IJ),IJ,Q(IJ),IJ,DO(I
    7     CONTINUE
        ELSE
          WRITE (ICHECK,910) CASE
  910     FORMAT (' ERROR IN SCF(CONTR2): CASE IS',I8)
          STOP
        END IF
        IF (NOREAD .GT. 0) THEN
C         CALL AWAIT(NOUT)
          IF (RDCNT .GE. MAXNO) THEN
            CALL SREW(NOUT)
            RDCNT=0
          END IF
          CALL SREAD(NOUT,BUFFER(POINTR),OUTBUF)
          NOREAD=NOREAD-1
          RDCNT=RDCNT+1
        END IF
        BTAG=BTAG+1
        POINTR=POINTR+OUTBUF
        IF (BTAG .GT. LIM) THEN
          BTAG=1
          POINTR=1
        END IF
    8 CONTINUE
C     AVG=FLOAT(SUM1)/SUMCNT
C     WRITE(6,'(/'  AVERAGE NUMBER OF CONSECITIVE IJ IS ',F10.3/)')AVG
C     WRITE(6,'('    VALUE RANGE             NUMBER ')')
C     WRITE(6,'('        X < 0.10          ',I7)')CNT3
C     WRITE(6,'('  0.10< X < 1.00          ',I7)')CNT2
C     WRITE(6,'('  1.00< X                 ',I7)')CNT1
C     STOP
      RETURN
      END
