C
C----------------------------------------------------------------------
C
      SUBROUTINE SORT1(NT,ITAPIN,ITAPOU,NIRRED,FLOV,ORBSYM,NCNTIN,
     .                 IPQ,IOFF,IADD,BUCK,CORE,LEFT,FZT,NFZO,NFZV,NIJ)
      IMPLICIT INTEGER (A-Z)
      REAL*8 CORE(LEFT),BUCK(NT*NT)
      INTEGER ORBSYM(NT),FLOV(NIRRED,2),IOFF(*),IADD(*),FZT(*),
     .        IPQ(*),NCNTIN(*),NFZO(*),NFZV(*)
C
C >>> THIS SORT LIMITED TO 500 LOOPS
C
      INTEGER FPAIR(500),LPAIR(500),FINT(500),LINT(500),NPQ(8),NIJ(8)
      DATA ITEMP /255/
C
      CALL RFILE(ITAPIN)
      CALL RFILE(ITAPOU)
C
      WRITE(6,*)
      WRITE(6,*) '>>>>>>   IN SORT1 <<<<<<<<<'
      WRITE(6,*)
C
C >>> CONSTRUCT SYMMETRY-PACKING ARRAYS FOR EXIT FILE
C
      TOFF=0
      SOFF=0
      PQPAIR=0
      LOOP=1
      FPAIR(1)=1
      LPAIR(1)=NT*(NT+1)/2
      FINT(1)=1
      LINT(1)=0
      DO 190 TSYM=1,NIRRED
         NOFF=0
         NPQ(TSYM)=0
         DO 180 P=1,NT
            PSYM=ORBSYM(P)
            DO 170 Q=1,P
               QSYM=ORBSYM(Q)
               PQSYM=IEOR(PSYM,QSYM)
               IF(PQSYM.NE.TSYM-1)GO TO 170
               PQ=IPQ(P)+Q
               IOFF(PQ)=TOFF
               NPQ(TSYM)=NPQ(TSYM)+1
               PQPAIR=PQPAIR+1
               NIJ(TSYM)=0
               NADD=0
               DO 150 I=1,NT
                  IF(FZT(I).EQ.1) GO TO 150
                  ISYM=ORBSYM(I)
                  JSYM=IEOR(ISYM,PQSYM)+1
                  FJ=FLOV(JSYM,1)
                  LJ=FLOV(JSYM,2)
                  IF(LJ.GT.I)LJ=I
                  DO 140 J=FJ,LJ
                     IJ=IPQ(I)+J
                     NOFF=NOFF+1
                     NADD=NADD+1
                     IADD(IJ)=NADD
                     TOFF=TOFF+1
                     SOFF=SOFF+1
                     TOT=IOFF(PQ)+IADD(IJ)
                     NIJ(TSYM)=NIJ(TSYM)+1
C               WRITE(6,634)P,Q,PQ,IOFF(PQ),I,J,IJ,IADD(IJ),TOT,TOFF
  634 FORMAT(10(1X,I6))
  140             CONTINUE
  150          CONTINUE
               IF(SOFF.LT.LEFT)THEN
                  LINT(LOOP)=TOFF
                  LPAIR(LOOP)=PQPAIR
               ELSE
                  LOOP=LOOP+1
                  FINT(LOOP)=LINT(LOOP-1)+1
                  FPAIR(LOOP)=LPAIR(LOOP-1)+1
                  LINT(LOOP)=TOFF
                  LPAIR(LOOP)=PQPAIR
                  SOFF=NADD
               ENDIF
  170       CONTINUE
  180    CONTINUE
  190 CONTINUE
C
      WRITE(6,*)
      WRITE(6,*)LOOP,' LOOPS ARE NEEDED TO RESORT FILE',ITAPIN
      WRITE(6,*)
C
      IF(LOOP.GT.500)THEN
         WRITE(6,*)' TOO MANY LOOPS !!   GET LOST! '
         STOP 'TOO MANY LOOPS IN SORT1'
      ENDIF
C
      WRITE(6,*  )'   I,   FPAIR(I),    LPAIR(I),    FINT(I),   LINT(I)'
      DO 231 I=1,LOOP
      WRITE(6,654)I,FPAIR(I),LPAIR(I),FINT(I),LINT(I)
  654 FORMAT(5(2X,I8))
  231 CONTINUE
      WRITE(6,*)
      WRITE(6,*  )'   ISYM        NRS(ISYM)      NIJ(ISYM)'
      TPAR=0
      SIZEOU=0
      DO 232 I=1,NIRRED
      WRITE(6,652)I,NPQ(I),NIJ(I)
      TPAR=TPAR+NPQ(I)
      SIZEOU=SIZEOU+NPQ(I)*NIJ(I)
  652 FORMAT(3(2X,I8))
  232 CONTINUE
      WRITE(6,*)
      WRITE(6,*)'NUMBER OF ELEMENTS IN DISK IS =',SIZEOU
      WRITE(6,*)
C
      PNTOU=1
      DO 1000 ILOOP=1,LOOP
      CALL ZERO(BUCK,NT*NT)
      CALL ZERO(CORE,LEFT)
      FFINT =FINT (ILOOP)
      LLINT =LINT (ILOOP)
      FFPAIR=FPAIR(ILOOP)
      LLPAIR=LPAIR(ILOOP)
      SHIFT=FFINT-1
      WRITE(6,*)'LOOP=',ILOOP,' SHIFT=',SHIFT
C
      PNTIN=1
      DO 290 TSYM=1,NIRRED
      DO 280 I=1,NT
      IF(FZT(I).EQ.1)GO TO 280
      ISYM=ORBSYM(I)
      DO 270 J=1,I
      IF(FZT(J).EQ.1)GO TO 270
      JSYM=ORBSYM(J)
      IJSYM=IEOR(ISYM,JSYM)
      IF(IJSYM.NE.TSYM-1)GO TO 270
      IJ=IPQ(I)+J
      LENGTH=NCNTIN(TSYM)
C     WRITE(6,650)I,J,LENGTH
  650 FORMAT(' READING PAIR',' I=',I3,' J=',I3,' OF LENGTH ',I6)
      CALL WREADW(ITAPIN,BUCK,INTOWP(LENGTH),PNTIN,PNTIN)
      NRS=0
      DO 250 R=1,NT
      RSYM=ORBSYM(R)
      SSYM=IEOR(RSYM,IJSYM)+1
      FS=FLOV(SSYM,1)-NFZO(SSYM)
      LS=FLOV(SSYM,2)+NFZV(SSYM)
      IF(LS.GT.R)LS=R
      DO 240 S=FS,LS
      RS=IPQ(R)+S
      NRS=NRS+1
      RSIJ=IOFF(RS)+IADD(IJ)
      IF(RSIJ.GE.FFINT.AND.RSIJ.LE.LLINT) THEN
         RSIJ=RSIJ-SHIFT
         CORE(RSIJ)=BUCK(NRS)
      ENDIF
  240 CONTINUE
  250 CONTINUE
  270 CONTINUE
  280 CONTINUE
  290 CONTINUE
C
      TOFF=0
      XRS=0
      DO 390 TSYM=1,NIRRED
      DO 380 R=1,NT
      RSYM=ORBSYM(R)
      DO 370 S=1,R
      SSYM=ORBSYM(S)
      RSSYM=IEOR(RSYM,SSYM)
      IF(RSSYM.NE.TSYM-1)GO TO 370
      XRS=XRS+1
      IF(XRS.LT.FFPAIR.OR.XRS.GT.LLPAIR)GO TO 365
      NADD=0
      DO 350 I=1,NT
      IF(FZT(I).EQ.1)GO TO 350
      ISYM=ORBSYM(I)
      JSYM=IEOR(ISYM,RSSYM)+1
      FJ=FLOV(JSYM,1)
      LJ=FLOV(JSYM,2)
      IF(LJ.GT.I)LJ=I
      DO 340 J=FJ,LJ
      NADD=NADD+1
      TOFF=TOFF+1
      BUCK(NADD)=CORE(TOFF)
  340 CONTINUE
  350 CONTINUE
      LENGTH=NIJ(TSYM)
C     WRITE(6,651)XRS,R,S,LENGTH
  651 FORMAT(' WRITING PAIR',I4,' R=',I3,' S=',I3,' OF LENGTH ',I6)
      CALL WWRITW(ITAPOU,BUCK,INTOWP(LENGTH),PNTOU,PNTOU)
  365 CONTINUE
  370 CONTINUE
  380 CONTINUE
  390 CONTINUE
C
 1000 CONTINUE
C
      CALL RCLOSE(ITAPIN,4)
      CALL RCLOSE(ITAPOU,3)
C
      RETURN
      END
