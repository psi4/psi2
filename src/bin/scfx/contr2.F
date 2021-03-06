C     //////////////////////////////////////////////////////////////
      SUBROUTINE CONTR2 (VALUE,BUFFER,DT,DO,P,Q,NMAT,TWOCON)
C
C
C
      IMPLICIT INTEGER (A-Z)
C
C
      REAL*8 VALUE(1),P(1),Q(1),DT(1),DO(1)
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
C
      ICONT4=ICONT4+1
      LIM=MIN(MAXNO,TAGMAX)
      NUMBUF=0
      NOREAD=MAXNO-TAGMAX
      IF (NOST.EQ.0) THEN
         DO 1 I=1,NMAT
            P(I)=0.0D+00
    1    CONTINUE
      ELSE
         DO 2 I=1,NMAT
            P(I)=0.0D+00
            Q(I)=0.0D+00
    2    CONTINUE
      END IF
C
      DO 8 NUMBUF=1,MAXNO
         CASE=BUFFER(POINTR)
         NOINBF=BUFFER(POINTR+1)
         IF (CASE.EQ.4) THEN
C     CALL CASE4(BUFFER(POINTR+2),DT,P,NOINBF)
            LOCAT=(POINTR+1)/2+HALF+1
C
            CALL CASE4(VALUE(LOCAT),BUFFER(POINTR+2),DT,P,NOINBF)
         ELSE IF (CASE.EQ.3) THEN
            N3I=POINTR+1
            N3RB=N3I/2+THIRD
            N3RC=N3RB+THIRD
            IF (TWOCON) THEN
               DO 3 JUNK=1,NOINBF
                  N3I=N3I+2
                  N3RB=N3RB+1
                  N3RC=N3RC+1
                  KL=BUFFER(N3I-1)
                  IJ=BUFFER(N3I)
CD    IF ((PRINT.AND.128).NE.0) PRINT '(' N3I=',I5,' N3RB=',I5,' N3RC=',
CD   1                    I5,' IJ=',I5,' KL=',I5)', N3I,N3RB,N3RC,IJ,KL
                  IF (IJ.LE.IJU) THEN
                     Q(IJ)=Q(IJ)+ALPHA1*VALUE(N3RB)*DO(KL)
                     Q(KL)=Q(KL)+ALPHA1*VALUE(N3RB)*DO(IJ)
                  ELSE
                     IF (KL.LE.IJU) THEN
                        QTEMP=VALUE(N3RB)+ALPHA2*VALUE(N3RC)
                        Q(IJ)=Q(IJ)+QTEMP*DO(KL)
                        Q(KL)=Q(KL)+QTEMP*DO(IJ)
                     ELSE
                        Q(IJ)=Q(IJ)+ALPHA3*VALUE(N3RB)*DO(KL
     *                  )
                        Q(KL)=Q(KL)+ALPHA3*VALUE(N3RB)*DO(IJ
     *                  )
                     END IF
                  END IF
                  P(IJ)=P(IJ)+VALUE(N3RB)*DT(KL)
                  P(KL)=P(KL)+VALUE(N3RB)*DT(IJ)
CD    IF ((PRINT.AND.128).NE.0) PRINT '(2(' P(',I4,')=',F10.7,' DT(',I4,
CD   1                    ')=',F10.7),/,2(' Q(',I4,'2)=',F10.7,' DO(',I4
CD   *                    ,')=',F10.7))', IJ, P(IJ),IJ,DT(IJ),KL,P(KL),K
CD   *L,DT(KL),IJ,Q(IJ),IJ,DO(IJ),KL,Q(KL),   KL,DO(KL)
    3          CONTINUE
            ELSE
               DO 4 JUNK=1,NOINBF
                  N3I=N3I+2
                  N3RB=N3RB+1
                  N3RC=N3RC+1
                  KL=BUFFER(N3I-1)
                  IJ=BUFFER(N3I)
                  P(IJ)=P(IJ)+VALUE(N3RB)*DT(KL)
                  P(KL)=P(KL)+VALUE(N3RB)*DT(IJ)
                  Q(IJ)=Q(IJ)+VALUE(N3RC)*DO(KL)
                  Q(KL)=Q(KL)+VALUE(N3RC)*DO(IJ)
CD    IF ((PRINT.AND.128).NE.0) PRINT '(2(' P(',I4,')=',F10.7,' DT(',I4,
CD   1                    ')=',F10.7),/,2(' Q(',I4,'2)=',F10.7,' DO(',I4
CD   *                    ,')=',F10.7))', IJ, P(IJ),IJ,DT(IJ),KL,P(KL),K
CD   *L,DT(KL),IJ,Q(IJ),IJ,DO(IJ),KL,Q(KL),   KL,DO(KL)
    4          CONTINUE
            END IF
         ELSE IF (CASE.EQ.2) THEN
            N2I=POINTR+1
            N2R=N2I/2+HALF
            DO 5 JUNK=1,NOINBF
               N2I=N2I+1
               N2R=N2R+1
               IJ=BUFFER(N2I)
               P(IJ)=P(IJ)+VALUE(N2R)*DT(IJ)
CD    IF ((PRINT.AND.128).NE.0) PRINT '(' P(',I4,')=',F10.7,' DT(',I4,
CD   1               ')=',F10.7)', IJ,P(IJ),IJ,DT(IJ)
    5       CONTINUE
         ELSE IF (CASE.EQ.1) THEN
            N1I=POINTR+1
            N1RA=N1I/2+THIRD
            IF (TWOCON) THEN
               DO 6 JUNK=1,NOINBF
                  N1I=N1I+1
                  N1RA=N1RA+1
                  IJ=BUFFER(N1I)
                  IF (IJ.LE.IJU) THEN
                     Q(IJ)=Q(IJ)+ALPHA1*VALUE(N1RA)*DO(IJ)
                  ELSE
                     Q(IJ)=Q(IJ)+ALPHA3*VALUE(N1RA)*DO(IJ)
                  END IF
                  P(IJ)=P(IJ)+VALUE(N1RA)*DT(IJ)
CD    IF ((PRINT.AND.128).NE.0) PRINT '(' P(',I4,')=',F10.7,' DT(',I4,
CD   1                    ')=',F10.7,' Q(',I4,')=',F10.7,' DO(',I4,')=',
CD   *                    F10.7)', IJ,P(IJ),  IJ,DT(IJ),IJ,Q(IJ),IJ,DO(I
    6          CONTINUE
            ELSE
               N1RB=N1RA+THIRD
               DO 7 JUNK=1,NOINBF
                  N1I=N1I+1
                  N1RA=N1RA+1
                  N1RB=N1RB+1
                  IJ=BUFFER(N1I)
                  P(IJ)=P(IJ)+VALUE(N1RA)*DT(IJ)
                  Q(IJ)=Q(IJ)+VALUE(N1RB)*DO(IJ)
CD    IF ((PRINT.AND.128).NE.0) PRINT '(' P(',I4,')=',F10.7,' DT(',I4,
CD   1                    ')=',F10.7,' Q(',I4,')=',F10.7,' DO(',I4,')=',
CD   *                    F10.7)', IJ,P(IJ),  IJ,DT(IJ),IJ,Q(IJ),IJ,DO(I
    7          CONTINUE
            END IF
         ELSE
            WRITE (ICHECK,910) CASE
  910       FORMAT (' ERROR IN SCF(CONTR2): CASE IS',I8)
            STOP
         END IF
         IF (NOREAD .GT. 0) THEN
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
      ITIME4=ITIME4+END-START
C     AVG=FLOAT(SUM1)/SUMCNT
C     WRITE(ICHECK,'(/'  AVERAGE NUMBER OF CONSECITIVE IJ IS ',
C    +F10.3/)')AVG
C     WRITE(ICHECK,'('    VALUE RANGE             NUMBER ')')
C     WRITE(ICHECK,'('        X < 0.10          ',I7)')CNT3
C     WRITE(ICHECK,'('  0.10< X < 1.00          ',I7)')CNT2
C     WRITE(ICHECK,'('  1.00< X                 ',I7)')CNT1
C     STOP
      RETURN
      END
