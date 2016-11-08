C
C----------------------------------------------------------------------
C
      SUBROUTINE SORTB(BUF,IBUF,LENGTH,NBF,
     .                 NO,ITAP61,ITAP71,JOUT,NIRRED,FLOV,ORBSYM,
     .                 IPQ,ITV,NV,BOFF,BADD,LEFT,CORE,GBUCK)
      IMPLICIT INTEGER (A-Z)
      REAL*8 BUF(LENGTH),RINT,CORE(LEFT),GBUCK(NV*NV)
      INTEGER IBUF(LENGTH*2)
      INTEGER ORBSYM(NBF),FLOV(NIRRED,4),BOFF(NV*NV),BADD(NV*NV),
     .        ITV(NBF),IPQ(NBF),
     .        NAC(8),NBD(8),
     .        FPAIR(100),LPAIR(100),FINT(100),LINT(100)
      integer prcntr
C
      DATA ITEMP /255/
      DATA MAXLP /100/
C
C >>> CONSTRUCT SYMMETRY-PACKING ARRAYS
C
      DO 50 I=1,NBF
         IPQ(I)=(I*(I-1))/2
         ITV(I)=(I-1)*NV
   50 CONTINUE
C
      TOFF=0
      SOFF=0
      ACPR=0
      LOOP=1
      FPAIR(1)=1
      LPAIR(1)=NV*(NV+1)/2
      FINT(1)=1
      LINT(1)=0
      DO 190 TSYM=1,NIRRED
         NOFF=0
         NAC(TSYM)=0
         DO 180 A=1,NV
            ASYM=ORBSYM(A+NO)
            DO 170 C=1,A
               CSYM=ORBSYM(C+NO)
               ACSYM=IEOR(ASYM,CSYM)
               IF(ACSYM.NE.TSYM-1)GO TO 170
               AC=IPQ(A)+C
               BOFF(AC)=TOFF
               NAC(TSYM)=NAC(TSYM)+1
               ACPR=ACPR+1
               NBD(TSYM)=0
               NADD=0
               DO 150 B=1,NV
                  BSYM=ORBSYM(B+NO)
                  DSYM=IEOR(BSYM,ACSYM)+1
                  FD=FLOV(DSYM,3)-NO
                  LD=FLOV(DSYM,4)-NO
                  DO 140 D=FD,LD
                     BD=ITV(B)+D
                     NOFF=NOFF+1
                     NADD=NADD+1
                     BADD(BD)=NADD
                     TOFF=TOFF+1
                     SOFF=SOFF+1
                     TOT=BOFF(AC)+BADD(BD)
                     NBD(TSYM)=NBD(TSYM)+1
C     WRITE(6,634)A,C,AC,BOFF(AC),B,D,BD,BADD(BD),TOT,TOFF
  634 FORMAT(10(1X,I6))
  140             CONTINUE
  150          CONTINUE
               IF(SOFF.LT.LEFT)THEN
                  LINT(LOOP)=TOFF
                  LPAIR(LOOP)=ACPR
               ELSE
                  LOOP=LOOP+1
                  FINT(LOOP)=LINT(LOOP-1)+1
                  FPAIR(LOOP)=LPAIR(LOOP-1)+1
                  LINT(LOOP)=TOFF
                  LPAIR(LOOP)=ACPR
                  SOFF=NADD
               ENDIF
  170       CONTINUE
  180    CONTINUE
  190 CONTINUE
C
      if (prcntr('IS_ON BRIEF').eq.0) then
        WRITE(6,*)
        WRITE(6,*)LOOP,' LOOPS ARE NEEDED TO RESORT BINTS'
        WRITE(6,*)
        WRITE(6,*  )'      I     FPAIR(I)   LPAIR(I)   FINT(I) LINT(I)'
        endif
      if (prcntr('IS_ON BRIEF').eq.0) then
        DO 231 I=1,LOOP
          WRITE(6,654)I,FPAIR(I),LPAIR(I),FINT(I),LINT(I)
  654     FORMAT(5(I8,2X))
  231     CONTINUE
        WRITE(6,*)
        WRITE(6,*  )'    ISYM     NAC(ISYM)  NBD(ISYM)'
        endif
      TPAR=0
      SIZE71=0
      DO 232 I=1,NIRRED
         if (prcntr('IS_ON BRIEF').eq.0) WRITE(6,652)I,NAC(I),NBD(I)
         TPAR=TPAR+NAC(I)
         SIZE71=SIZE71+NAC(I)*NBD(I)
  652 FORMAT(3(I8,2X))
  232 CONTINUE
      if (prcntr('IS_ON BRIEF').eq.0) then
        WRITE(6,*)
        WRITE(6,*)'NV2 =', NV*(NV+1)/2
        WRITE(6,*)'TPAR=', TPAR
        WRITE(6,*)
        WRITE(6,*)'SIZE71=',SIZE71
        WRITE(6,*)
        endif
C
      IF(LOOP.GT.MAXLP)THEN
         WRITE(6,*)LOOP,' loops are needed to resort bints'
         write(6,*)maxlp,' is the maximum number of loops'
         WRITE(6,*)'ERROR: number of required loops exceeds maximum'
         WRITE(6,*)' increase maxloop in nccsort'
         call qabort
      ENDIF
C
      INTLEN=(INTOWP(LENGTH)-2)/INTOWP(1)
      MAXVAL=INTOWP(INTLEN)/(1+INTOWP(1))
      IVOFF=(MAXVAL+3)/INTOWP(1)
C
      CALL RFILE(ITAP71)
      POINT=1
      DO 1000 ILOOP=1,LOOP
         CALL ZERO(GBUCK,NV*NV)
         CALL ZERO(CORE ,LEFT )
         FFINT =FINT (ILOOP)
         LLINT =LINT (ILOOP)
         FFPAIR=FPAIR(ILOOP)
         LLPAIR=LPAIR(ILOOP)
         SHIFT=FFINT-1
         if (prcntr('IS_ON BRIEF').eq.0)
     &     WRITE(6,*)'LOOP=',ILOOP,' SHIFT=',SHIFT
         CALL SREW(ITAP61)
  111    CALL ZERO(BUF,LENGTH)
         CALL SREAD(ITAP61,IBUF,INTOWP(LENGTH))
         IFLG=IBUF(1)
         MBUF=IBUF(2)
C
      DO 101 II=1,MBUF
         ABCD=IBUF(2+II)
         A=ISHFT(ABCD,-24)
         B=IAND(ITEMP,ISHFT(ABCD,-16))
         C=IAND(ITEMP,ISHFT(ABCD,-8))
         D=IAND(ITEMP,ABCD)
         A=A-NO
         B=B-NO
         C=C-NO
         D=D-NO
         RINT=BUF(IVOFF+II)
C        WRITE (6,19) A,B,C,D,RINT
   19 FORMAT (4I3,F20.12,3X,'**')
C
      AC=IPQ(A)+C
      BD=ITV(B)+D
      ACBD=BOFF(AC)+BADD(BD)
      IF(ACBD.GE.FFINT.AND.ACBD.LE.LLINT) THEN
         ACBD=ACBD-SHIFT
         CORE(ACBD)=RINT
      ENDIF
C
      IF(B.GE.C) THEN
         BC=IPQ(B)+C
         AD=ITV(A)+D
         BCAD=BOFF(BC)+BADD(AD)
         IF(BCAD.GE.FFINT.AND.BCAD.LE.LLINT) THEN
            BCAD=BCAD-SHIFT
            CORE(BCAD)=RINT
         ENDIF
      ENDIF
C
      AD=IPQ(A)+D
      BC=ITV(B)+C
      ADBC=BOFF(AD)+BADD(BC)
      IF(ADBC.GE.FFINT.AND.ADBC.LE.LLINT) THEN
         ADBC=ADBC-SHIFT
         CORE(ADBC)=RINT
      ENDIF
C
      IF(B.GE.D) THEN
         BD=IPQ(B)+D
         AC=ITV(A)+C
         BDAC=BOFF(BD)+BADD(AC)
         IF(BDAC.GE.FFINT.AND.BDAC.LE.LLINT) THEN
            BDAC=BDAC-SHIFT
            CORE(BDAC)=RINT
         ENDIF
      ENDIF
C
      IF(C.GE.A) THEN
         CA=IPQ(C)+A
         DB=ITV(D)+B
         CADB=BOFF(CA)+BADD(DB)
         IF(CADB.GE.FFINT.AND.CADB.LE.LLINT) THEN
            CADB=CADB-SHIFT
            CORE(CADB)=RINT
         ENDIF
      ENDIF
C
      IF(C.GE.B) THEN
         CB=IPQ(C)+B
         DA=ITV(D)+A
         CBDA=BOFF(CB)+BADD(DA)
         IF(CBDA.GE.FFINT.AND.CBDA.LE.LLINT) THEN
            CBDA=CBDA-SHIFT
            CORE(CBDA)=RINT
         ENDIF
      ENDIF
C
      IF(D.GE.A) THEN
         DA=IPQ(D)+A
         CB=ITV(C)+B
         DACB=BOFF(DA)+BADD(CB)
         IF(DACB.GE.FFINT.AND.DACB.LE.LLINT) THEN
            DACB=DACB-SHIFT
            CORE(DACB)=RINT
         ENDIF
      ENDIF
C
      IF(D.GE.B) THEN
         DB=IPQ(D)+B
         CA=ITV(C)+A
         DBCA=BOFF(DB)+BADD(CA)
         IF(DBCA.GE.FFINT.AND.DBCA.LE.LLINT) THEN
            DBCA=DBCA-SHIFT
            CORE(DBCA)=RINT
         ENDIF
      ENDIF
C
  101 CONTINUE
      IF (IFLG.EQ.0) GOTO 111
C
      TOFF=0
      IAC=0
      DO 390 TSYM=1,NIRRED
         DO 380 A=1,NV
            ASYM=ORBSYM(A+NO)
            DO 370 C=1,A
               CSYM=ORBSYM(C+NO)
               ACSYM=IEOR(ASYM,CSYM)
               IF(ACSYM.NE.TSYM-1)GO TO 370
               IAC=IAC+1
               IF(IAC.LT.FFPAIR.OR.IAC.GT.LLPAIR)GO TO 365
               NADD=0
               DO 350 B=1,NV
                  BSYM=ORBSYM(B+NO)
                  DSYM=IEOR(BSYM,ACSYM)+1
                  FD=FLOV(DSYM,3)-NO
                  LD=FLOV(DSYM,4)-NO
                  DO 340 D=FD,LD
                     NADD=NADD+1
                     TOFF=TOFF+1
                     GBUCK(NADD)=CORE(TOFF)
  340             CONTINUE
  350          CONTINUE
               LBD=NBD(TSYM)
C              WRITE(6,651)IAC,A,C,NBD(TSYM)
  651 FORMAT(' WRITING PAIR',I4,' A=',I3,' C=',I3,' OF LENGTH ',I6)
               CALL WWRITW(ITAP71,GBUCK,INTOWP(LBD),POINT,POINT)
  365       CONTINUE
  370    CONTINUE
  380 CONTINUE
  390 CONTINUE
C
 1000 CONTINUE
      CALL RCLOSE(ITAP71,3)
C
      RETURN
      END
