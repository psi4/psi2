C
C----------------------------------------------------------------------
C
      SUBROUTINE SORT2(NT,NTAO,ITAPIN,ITAPOU,IPQ,IOFF,IADD,NT2,NTAO2,
     .                 BUCK,CORE,LEFT,KLOC,AOPS,NSHELL,ORBSYM,NIRRED)
      IMPLICIT INTEGER (A-Z)
      REAL*8 CORE(LEFT),BUCK(NT*NT)
      INTEGER IOFF(*),IADD(*),IPQ(*),KLOC(*),AOPS(*),ORBSYM(*)
C
C >>> THIS SORT LIMITED TO 500 LOOPS
C
      INTEGER FPAIR(500),LPAIR(500),FINT(500),LINT(500)
      DATA ITEMP /255/
C
      CALL RFILE(ITAPIN)
      CALL RFILE(ITAPOU)
C
      WRITE(6,*)
      WRITE(6,*) '>>>>>>   IN SORT2 <<<<<<<<<'
      WRITE(6,*)
C
C >>> CONSTRUCT SYMMETRY-PACKING ARRAYS FOR EXIT FILE
C
      TOFF=0
      SOFF=0
      PQPAIR=0
      LOOP=1
      FPAIR(1)=1
      LPAIR(1)=NTAO*(NTAO+1)/2
      FINT(1)=1
      LINT(1)=0
      DO 195 MUSHEL=1,NSHELL
         FMU=KLOC(MUSHEL)
         LMU=FMU+AOPS(MUSHEL)-1
         DO 185 NUSHEL=1,MUSHEL
            FNU=KLOC(NUSHEL)
            LNU=FNU+AOPS(NUSHEL)-1
      DO 190 MU=FMU,LMU
         IF(NUSHEL.EQ.MUSHEL)LNU=MU
         DO 180 NU=FNU,LNU
               MUNU=IPQ(MU)+NU
               IOFF(MUNU)=TOFF
               PQPAIR=PQPAIR+1
               NADD=0
               DO 150 R=1,NT
                  DO 140 S=1,R
                     RS=IPQ(R)+S
                     NADD=NADD+1
                     IADD(RS)=NADD
                     TOFF=TOFF+1
                     SOFF=SOFF+1
                     TOT=IOFF(MUNU)+IADD(RS)
C         WRITE(6,634)MU,NU,MUNU,IOFF(MUNU),R,S,RS,IADD(RS),TOT,TOFF
  634     FORMAT(10(1X,I6))
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
  180    CONTINUE
  190    CONTINUE
  185 CONTINUE
  195 CONTINUE
C
      WRITE(6,*)
      WRITE(6,*)LOOP,' LOOPS ARE NEEDED TO RESORT FILE',ITAPIN
      WRITE(6,*)
C
      IF(LOOP.GT.500)THEN
         WRITE(6,*)' TOO MANY LOOPS !!   GET LOST! '
         STOP 'TOO MANY LOOPS IN SORT2'
      ENDIF
C
      WRITE(6,*  )'   I,   FPAIR(I),    LPAIR(I),    FINT(I),   LINT(I)'
      DO 231 I=1,LOOP
      WRITE(6,654)I,FPAIR(I),LPAIR(I),FINT(I),LINT(I)
  654 FORMAT(5(2X,I8))
  231 CONTINUE
      WRITE(6,*)
      WRITE(6,*)'NUMBER OF ELEMENTS IN DISK IS =',TOT
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
      ISYM=ORBSYM(I)
      DO 270 J=1,I
      JSYM=ORBSYM(J)
      IJSYM=IEOR(ISYM,JSYM)
      IF(IJSYM.NE.TSYM-1)GO TO 270
      IJ=IPQ(I)+J
C     WRITE(6,650)I,J,NTAO2
  650 FORMAT(' READING PAIR',' I=',I3,' J=',I3,' OF LENGTH ',I6)
      CALL WREADW(ITAPIN,BUCK,INTOWP(NTAO2),PNTIN,PNTIN)
      NRS=0
      DO 250 R=1,NTAO
      DO 240 S=1,R
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
      DO 395 MUSHEL=1,NSHELL
         FMU=KLOC(MUSHEL)
         LMU=FMU+AOPS(MUSHEL)-1
         DO 385 NUSHEL=1,MUSHEL
            FNU=KLOC(NUSHEL)
            LNU=FNU+AOPS(NUSHEL)-1
      DO 390 MU=FMU,LMU
         IF(NUSHEL.EQ.MUSHEL)LNU=MU
         DO 380 NU=FNU,LNU
               MUNU=IPQ(MU)+NU
      XRS=XRS+1
      IF(XRS.LT.FFPAIR.OR.XRS.GT.LLPAIR)GO TO 365
      NADD=0
      DO 350 I=1,NT
      DO 340 J=1,I
      NADD=NADD+1
      TOFF=TOFF+1
      BUCK(NADD)=CORE(TOFF)
  340 CONTINUE
  350 CONTINUE
C     WRITE(6,651)XRS,MU,NU,NT2
  651 FORMAT(' WRITING PAIR',I4,' MU=',I3,' NU=',I3,' OF LENGTH ',I6)
      CALL WWRITW(ITAPOU,BUCK,INTOWP(NT2),PNTOU,PNTOU)
  365 CONTINUE
  380 CONTINUE
  390 CONTINUE
  385 CONTINUE
  395 CONTINUE
C
 1000 CONTINUE
C
      CALL RCLOSE(ITAPIN,4)
      CALL RCLOSE(ITAPOU,3)
C
      RETURN
      END
