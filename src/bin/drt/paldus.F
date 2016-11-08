      SUBROUTINE PALDUS(BFSYM,BFCODE,ORBTBF,A,B,S
     #,                 LEVPT,LEVNR,ARC,NLWKS,X,NELECS,SPC)
C
C***********************************************************************
C     THIS SUBROUTINE GENERATES THE DISTINCT ROW TABLE BY SEARCHING    *
C     FROM THE TOP OF THE GRAPH DOWN POSSIBLE PATHS. NOTE THAT THERE   *
C     IS A PROGRAM LOOP OVER ORBITALS RUNNING FROM STATEMENT LABELED   *
C     '4'. THE IDEA IS FOR EACH LEVEL, SEARCH ACROSS ALL POSSIBLE      *
C     A AND B COMBINATIONS (ATEST AND BTEST) THEN ACROSS ALL POINTS    *
C     ON THE ROW ABOVE (ROW) TESTING WHETHER A PARTICULAR ARC (CASE)   *
C     GETS DOWN TO THE POINTS GIVEN BY ATEST,BTEST. THE COMPLICATION   *
C     LIES IN COMPUTING THE EXCITATION LEVEL WHEN ELECTRONS ARE        *
C     EXCITED INTO OCCUPIED OR MULTI-REFERENCE ORBITALS. THIS IS       *
C     ACCOMPLISHED USING THE X ARRAY WHICH CONTAINS THE NUMBER OF      *
C     ELECTRONS EXCITED INTO ORBITALS ABOVE THE PRESENT POINT FOR      *
C     ANY WALK LEADING TO THE PARTICULAR POINT IN QUESTION. THEREFORE, *
C     POINTS MAY BE THE SAME IN ALL RESPECTS EXCEPT FOR X VALUE AND    *
C     FINALLY, AT THE FERMI LEVEL THE X VALUE DETERMINES IF THE WALK   *
C     IS POSSIBLE OR NOT. FOR SIMPLE CASES SUCH AS HIGH SPIN OPEN      *
C     SHELLS THE INTERACTING SPACE IS IMPLIMENTED BY COUNTING SPIN-    *
C     FLIPS AS EXCITATIONS INTO ORBITALS.                              *
C                                                                      *
C     THE VARIABLES AMAX, AMIN, BMAX, AND BMIN LIMIT THE PORTION       *
C     OF THE SEARCH TO ONLY THOSE POINTS THAT CAN BE REACHED FROM      *
C     THE PRESENT ROW. SPECIAL ORBITALS ARE HANDLED BY INPUT           *
C     EXPLICITLY THE X VALUES FOR ALL POSSIBLE PARTIAL WALKS IN THE    *
C     SPECIAL ORBITALS. ORIGINALLY THE EXCITATION LEVEL IS SET         *
C     TO THE SUM OF THE GENERAL VALUE AND THE VALUE IN THE VALENCE     *
C     (%) SPACE. UPON LEAVING THE VALENCE SPACE IT IS RESET TO THE     *
C     GENERAL VALUE, THUS GIVING SAY ALL SINGLES AND DOUBLES           *
C     FROM REFERENCES OF ALL SINGLES AND DOUBLES IN THE VALENCE        *
C     SPACE.                                                           *
C                                                                      *
C     MAXB KEEPS TRACK OF THE LARGEST B VALUE ENCOUNTERED. THIS        *
C     CAN THEN BE USED IN THE CI PROGRAM TO COMPUTE SUFFICIENT         *
C     COEFFICIENTS FOR SEGMENT VALUES. IF THIS IS DONE, THERE IS       *
C     NO LIMIT TO THE VALUE OF THE SPIN POSSIBLE. FINALLY THE LAST     *
C     PORTION OF THE ROUTINE ELIMINATES ALL POINTS AND ARCS FROM       *
C     WALKS THAT DONT MAKE IT FROM HEAD TO TAIL.                       *
C***********************************************************************
C
      IMPLICIT INTEGER (A-Z)
      INTEGER NUMINT
      LOGICAL SPINI,SYMI
      CHARACTER*1 MULTRF,VALENC,KEY,BFKEY
      CHARACTER*3 CODES,WORDS
      CHARACTER*4 LABEL
C
      COMMON /CBFKEY/ BFKEY(1000)
      COMMON /TAPES/  OUT,ERROUT,INPUT,DRTTAP
      COMMON /DRTINF/ NA,NB,NS,NESPEC,MAXB,LEVFRM,LEVVAL,LEVOPN,LEVMUL
     #,               LEVOCC,SPEC,SSPESH,VAL
      COMMON /DRTCOD/ NCODES,DELA(9),DELB(9),DELELE(9)
     #,               NTYPES,VIRTUL,OCCUPD,VALOCC,RESCOR,RESVIR,FROZEN
     #,               VALVIR,OPENSH,MULTI,SPESHL,MULTRF,VALENC
     #,               WORDS(6,9),CODES(9)
      COMMON /CODE/  FZC, FZV, COR, VIR, DOC, UOC, ALP, BET, SPE
      COMMON /DIMENS/ NBF,NSYM,NORBS,NROWSP,NROW4P,NROWS,NROWS4
     #,               NLEVS,NREFS,NROWOC,NROW4O,NWKS,NWKSOC,NLEVOC
     #,               ORBFRM,SYMORB,NUMIJ,NGROUP,NUMINT,NMAX,NSPC,NVREF
     #,               NIJVIR
      COMMON /OPTION/ OPTION(16),IDENT(26),LABEL(26)
C
      DIMENSION BFSYM(NBF),BFCODE(NREFS,NBF),ORBTBF(NORBS)
      DIMENSION A(NROWSP),B(NROWSP),S(NROWSP),NLWKS(NROWSP)
      DIMENSION LEVPT(NLEVS),LEVNR(NLEVS),ARC(NROW4P)
      DIMENSION X(NREFS,NROWSP),NELECS(NREFS),SPC(NSPC)
C
      DATA SPCNT /0/
C
C
      SPINI=IAND(OPTION(4),1).EQ.0
      SYMI =IAND(OPTION(4),2).NE.0
C
      MAXB=-999999
C     IF (OPTION(2).EQ.0) OPTION(2)=2
      OPTION(3)=OPTION(2)+OPTION(3)
      EXCITA=OPTION(3)
      DO 1 REF=1,NREFS
         NELECS(REF)=2*NA+NB
         X(REF,1)=0
    1 CONTINUE
      LEVM1=NLEVS
      LEVPT(LEVM1)=0
      LEVNR(LEVM1)=1
      A(1)=NA
      B(1)=NB
      S(1)=NS
      X(1,1)=0
      RS1=NS
C
C     ----- HERE BEGINS THE LOOP OVER ORBITALS (LEVELS) -----
C
    4 IF (LEVM1.LE.1) GO TO 25
      LEV=LEVM1
      LEVM1=LEV-1
      IF (LEVM1.LE.LEVOCC) EXCITA=OPTION(3)
      IF (LEVM1.LT.LEVVAL) EXCITA=OPTION(2)
      NROWM1=0
      BF=ORBTBF(LEVM1)
      KEY=BFKEY(BF)
      CODE=BFCODE(1,BF)
      NROWLV=LEVNR(LEV)
      PONTLV=LEVPT(LEV)
      PONTM1=PONTLV+NROWLV
      LEVPT(LEVM1)=PONTM1
      IF (CODE.NE.SPE) GO TO 5
      IF (SPCNT.EQ.0) NELECS(1)=NELECS(1)-NESPEC
      SPCNT=SPCNT+1
      GO TO 11
    5 CONTINUE
      IF (DELELE(BFCODE(1,BF)).EQ.1) RS1=IEOR(RS1,BFSYM(BF))
      DO 2 REF=1,NREFS
         NELECS(REF)=NELECS(REF)-DELELE(BFCODE(REF,BF))
    2 CONTINUE
   11 CONTINUE
C
      AMAX=0
      AMIN=999999
      BMAX=0
      BMIN=999999
      DO 3 ROW=PONTLV+1,PONTLV+NROWLV
         IF (A(ROW).GT.AMAX) AMAX=A(ROW)
         IF (A(ROW).LT.AMIN) AMIN=A(ROW)
         IF (B(ROW).GT.BMAX) BMAX=B(ROW)
         IF (B(ROW).LT.BMIN) BMIN=B(ROW)
    3 CONTINUE
      IF (AMIN.GT.0) AMIN=AMIN-1
      IF (BMIN.GT.0) BMIN=BMIN-1
      IF (BMAX.GT.MAXB) MAXB=BMAX
      BMAX=BMAX+1
C
      DO 23 AJUNK=1,AMAX-AMIN+1
         ATEST=AMAX-AJUNK+1
         DO 22 BJUNK=1,BMAX-BMIN+1
            BTEST=BMAX-BJUNK+1
            DO 21 ROW=PONTLV+1,PONTLV+NROWLV
               DO 20 CASE=1,4
                  IF (CODE.EQ.COR.AND.CASE.NE.4) GO TO 19
                  IF (CODE.EQ.VIR.AND.CASE.NE.1) GO TO 19
                  IA=A(ROW)
                  IB=B(ROW)
                  IS=S(ROW)
                  IC=LEVM1-IA-IB
                  GO TO (6,7,8,9),CASE
C
                  WRITE (OUT,18) CASE
                  WRITE (ERROUT,18) CASE
   18             FORMAT (//,' IMPOSSIBLE CASE VALUE IN PALDUS:',I5)
                  CALL mabort
C
    6             IC=IC-1
                  DELE=0
                  GO TO 10
C
    7             IB=IB-1
                  IS=IEOR(IS,BFSYM(BF))
                  DELE=1
                  GO TO 10
C
    8             IA=IA-1
                  IB=IB+1
                  IC=IC-1
                  IS=IEOR(IS,BFSYM(BF))
                  DELE=1
                  IF (CODE.EQ.ALP.AND.SPINI) DELE=2
                  GO TO 10
C
    9             CONTINUE
                  IA=IA-1
                  DELE=2
C
   10             CONTINUE
                  IF (IA.LT.0.OR.IB.LT.0.OR.IC.LT.0) GO TO 19
                  IF (IA.NE.ATEST.OR.IB.NE.BTEST) GO TO 19
                  IF (2*IA+IB.GT.EXCITA.AND.LEVM1.LE.LEVFRM) GO TO 19
                  POSSBL=0
                  DO 60 REF=1,NREFS
                     IF ((2*IA+IB.LE.NELECS(REF)+EXCITA).AND
     #               .   (2*IA+IB.GE.NELECS(REF)-EXCITA)) POSSBL=1
   60             CONTINUE
                  IF (POSSBL.EQ.0.) GO TO 19
C
                  DO 13 ROWM1=PONTM1+1,PONTM1+NROWM1
                     IF (IA.NE.A(ROWM1).OR.IB.NE.B(ROWM1).OR.IS.NE.S(
     *               ROWM1)) GO TO 13
                     IF (LEVM1.LT.LEVFRM) GO TO 17
                     DIFF=0
                     POSSBL=0
                     IF (LEVVAL-LEVFRM.EQ.1) POSSBL=1
                     IF (CODE.NE.SPE) GO TO 62
                     IF (SPCNT.EQ.SPEC) GO TO 61
                     POSSBL=1
                     IX=4*(X(1,ROW)+CASE-1)
                     IF (IX.NE.X(1,ROWM1)) DIFF=1
                     GO TO 63
   61                CONTINUE
                     REF=1
                     IX=SPC(X(1,ROW)+CASE)
                     FRMX=0
                     IF (LEVM1.EQ.LEVFRM) FRMX=IX
                     IF ((2*IA+IB.LE.NELECS(REF)+EXCITA-IX).AND
     #               .   (2*IA+IB.GE.NELECS(REF)-EXCITA+FRMX)) POSSBL=1
                     IF (LEVM1.EQ.LEVFRM) IX=0
                     IF (IX.NE.X(1,ROWM1)) DIFF=1
                     GO TO 63
   62                CONTINUE
                     DO 40 REF=1,NREFS
                        IX=X(REF,ROW)
                        IF (DELE.GT.DELELE(BFCODE(REF,BF)).AND.KEY.NE.
     *                  VALENC)IX=IX+DELE-DELELE(BFCODE(REF,BF))
                        IF (LEV.NE.LEVVAL) GO TO 70
                        IXSV=NELECS(REF)-(2*IA+IB)
                        IX=IX+IXSV
                        IF (IX.LT.0) THEN
                           IF (SPINI) IX=(IX+IB)/2
                        ELSE
                           IF (SPINI.AND.IB.GT.IX) IX=IX+(IB-IX)/2
                        END IF
                        IF (IX.LT.0) IX=0
                        IF (SYMI.AND.IX.EQ.0.AND.IS.NE.RS1.AND.IXSV.EQ.0
     *                  ) IX=1
                        IF (2*IA+IB.LE.EXCITA-IX) IX=0
   70                   CONTINUE
                        FRMX=0
                        IF (LEVM1.EQ.LEVFRM) FRMX=IX
                        IF ((2*IA+IB.LE.NELECS(REF)+EXCITA-IX).AND
     #                  .   (2*IA+IB.GE.NELECS(REF)-EXCITA+FRMX))
     *                   POSSBL=1
                        IF (LEVM1.EQ.LEVFRM) IX=0
   40                IF (IX.NE.X(REF,ROWM1)) DIFF=1
   63                CONTINUE
                     IF (POSSBL.EQ.0) GO TO 19
                     IF (DIFF.EQ.0) GO TO 17
   13             CONTINUE
C
C     ----- CHECK THAT THIS IS INDEED A POSSIBLE NEW POINT,     -----
C           ESPECIALLY, IF GOT HERE BECAUSE A, B OR S DIFFERENT
C
                  POSSBL=0
                  IF (LEVVAL-LEVFRM.EQ.1) POSSBL=1
                  IF (CODE.NE.SPE) GO TO 162
                  IF (SPCNT.EQ.SPEC) GO TO 161
                  POSSBL=1
                  GO TO 163
  161             CONTINUE
                  REF=1
                  IX=SPC(X(1,ROW)+CASE)
                  FRMX=0
                  IF (LEVM1.EQ.LEVFRM) FRMX=IX
                  IF ((2*IA+IB.LE.NELECS(REF)+EXCITA-IX).AND
     #            .   (2*IA+IB.GE.NELECS(REF)-EXCITA+FRMX)) POSSBL=1
                  GO TO 163
  162             CONTINUE
                  DO 140 REF=1,NREFS
                  IX=X(REF,ROW)
                  IF (DELE.GT.DELELE(BFCODE(REF,BF)).AND.KEY.NE.VALENC)
     #            IX=IX+DELE-DELELE(BFCODE(REF,BF))
                  IF (LEV.NE.LEVVAL) GO TO 170
                  IXSV=NELECS(REF)-(2*IA+IB)
                  IX=IX+IXSV
                  IF (IX.LT.0) THEN
                     IF (SPINI) IX=(IX+IB)/2
                  ELSE
                     IF (SPINI.AND.IB.GT.IX) IX=IX+(IB-IX)/2
                  END IF
                  IF (IX.LT.0) IX=0
                  IF (SYMI.AND.IX.EQ.0.AND.IS.NE.RS1.AND.IXSV.EQ.0) IX=1
                  IF (2*IA+IB.LE.EXCITA-IX) IX=0
  170             CONTINUE
                  FRMX=0
                  IF (LEVM1.EQ.LEVFRM) FRMX=IX
                  IF ((2*IA+IB.LE.NELECS(REF)+EXCITA-IX).AND
     #            .   (2*IA+IB.GE.NELECS(REF)-EXCITA+FRMX)) POSSBL=1
  140             CONTINUE
  163             CONTINUE
                  IF (POSSBL.EQ.0) GO TO 19
C
                  NROWM1=NROWM1+1
                  ROWM1=PONTM1+NROWM1
                  IF (ROWM1.LT.NROWSP) GO TO 16
                  WRITE (OUT,15) NROWSP,LEVM1
                  WRITE (ERROUT,15) NROWSP,LEVM1
   15             FORMAT (//,' NOT ENOUGH SPACE TO MAKE DRT, HAVE ONLY R
     *OOM FOR',           I6,' ROWS AND ARE EXCEEDING THAT AT LEVEL',I4)
                  CALL mabort
   16             A(ROWM1)=IA
                  B(ROWM1)=IB
                  S(ROWM1)=IS
C
                  IF (LEVM1.LE.LEVFRM) GO TO 50
                  IF (LEV.NE.LEVVAL) GO TO 42
                  DO 77 REF=1,NREFS
                  IXSV=NELECS(REF)-(2*IA+IB)
                  IX=IXSV
                  IF (IX.LT.0) THEN
                     IF (SPINI) IX=(IX+IB)/2
                  ELSE
                     IF (SPINI.AND.IB.GT.IX) IX=IX+(IB-IX)/2
                  END IF
                  IF (IX.LT.0) IX=0
                  IF (SYMI.AND.IX.EQ.0.AND.IS.NE.RS1.AND.IXSV.EQ.0) IX=1
                  IF (2*IA+IB.LE.EXCITA-IX) IX=0
                  X(REF,ROWM1)=X(REF,ROW)+IX
   77             CONTINUE
                  GO TO 17
C
   42             CONTINUE
                  IF (CODE.NE.SPE) GO TO 44
                  IF (SPCNT.EQ.SPEC) GO TO 43
                  X(1,ROWM1)=4*(X(1,ROW)+CASE-1)
                  GO TO 17
   43             CONTINUE
                  X(1,ROWM1)=SPC(X(1,ROW)+CASE)
                  GO TO 17
   44             CONTINUE
                  DO 41 REF=1,NREFS
                     IX=X(REF,ROW)
                     IF (DELE.GT.DELELE(BFCODE(REF,BF)).AND.KEY.NE.
     *               VALENC)IX=IX+DELE-DELELE(BFCODE(REF,BF))
                     IF (2*IA+IB.LE.EXCITA-IX) IX=0
                     X(REF,ROWM1)=IX
   41             CONTINUE
                  GO TO 17
   50             CONTINUE
                  DO 51 REF=1,NREFS
                     X(REF,ROWM1)=0
   51             CONTINUE
   17             CONTINUE
                  ARC((ROW-1)*4+CASE)=ROWM1-PONTM1
   19             CONTINUE
   20          CONTINUE
   21       CONTINUE
   22    CONTINUE
   23 CONTINUE
      LEVNR(LEVM1)=NROWM1
      GO TO 4
   25 CONTINUE
C
C     ----- SET WEIGHT OF BOTTOM OF SHAVITT GRAPH TO ONE, -----
C                  ELIMINATE ALL OTHER BOTTOMS
C
      NROOT=0
      DO 27 ROOT=LEVPT(1)+1,LEVPT(1)+LEVNR(1)
      IF (A(ROOT).NE.0.OR.B(ROOT).NE.0.OR.S(ROOT).NE.0) GO TO 26
      NROOT=NROOT+1
      NLWKS(ROOT)=1
   26 CONTINUE
   27 CONTINUE
C
      IF (NROOT.EQ.1) GO TO 29
      WRITE (OUT,28) NROOT
      WRITE (ERROUT,28) NROOT
   28 FORMAT (//,' INVALID NUMBER OF BOTTOMS TO GRAPH:',I4,//)
      CALL mabort
   29 CONTINUE
C
C     ----- GENERATE THE WEIGHTS (NLWKS) OF ALL ROWS -----
C
      IF (LEV.LT.2) GO TO 91
      DO 32 LEV=2,NLEVS
         LEVM1=LEV-1
         PONTM1=LEVPT(LEVM1)
         DO 31 ROW=LEVPT(LEV)+1,LEVPT(LEV)+LEVNR(LEV)
            NLWK=0
            DO 30 CASE=1,4
               ARCPT=(ROW-1)*4+CASE
               IF (ARC(ARCPT).GT.0) NLWK=NLWK+NLWKS(ARC(ARCPT)+PONTM1)
   30       CONTINUE
            NLWKS(ROW)=NLWK
   31    CONTINUE
   32 CONTINUE
C
      NWKS=NLWKS(1)
      IF (.NOT.SPINI) GO TO 93
      WRITE (OUT,92) NWKS
   92 FORMAT (//,T25,'***** NUMBER OF INTERACTING WALKS *****',/,T25,'*'
     #,      T63,'*',/,T25,'*',T35,I12,T63,'*',/,T25,'*',T63,'*',/
     #,          T25,'***************************************')
      GO TO 95
   93 CONTINUE
      WRITE (OUT,94) NWKS
   94 FORMAT (//,T25,'******** TOTAL NUMBER OF WALKS ********',/T25,'*'
     #,       T63,'*',/,T25,'*',T35,I12,T63,'*',/,T25,'*',T63,'*',/
     #,          T25,'***************************************')
   95 CONTINUE
C
C     ----- REMOVE ALL ROWS WITH ZERO WEIGHTS -----
C
      PONT=1
      DO 37 JUNK=1,NLEVS-1
         LEV=NLEVS-JUNK
         NROWLV=LEVNR(LEV)
         PONTLV=LEVPT(LEV)
         LEVPT(LEV)=PONT
         DO 36 ROW=PONTLV+1,PONTLV+NROWLV
         IF (NLWKS(ROW).EQ.0) GO TO 35
         PONT=PONT+1
         A(PONT)=A(ROW)
         B(PONT)=B(ROW)
         S(PONT)=S(ROW)
         NLWKS(PONT)=NLWKS(ROW)
         DO 34 CASE=1,4
            ARC((PONT-1)*4+CASE)=ARC((ROW-1)*4+CASE)
   34    CONTINUE
         DO 55 ROWP1=LEVPT(LEV+1)+1,LEVPT(LEV+1)+LEVNR(LEV+1)
         DO 54 CASE=(ROWP1-1)*4+1,(ROWP1-1)*4+4
            IF (ARC(CASE).EQ.ROW-PONTLV) ARC(CASE)=PONT-LEVPT(LEV)
   54    CONTINUE
   55    CONTINUE
         GO TO 36
   35    CONTINUE
         DO 53 ROWP1=LEVPT(LEV+1)+1,LEVPT(LEV+1)+LEVNR(LEV+1)
            DO 52 CASE=(ROWP1-1)*4+1,(ROWP1-1)*4+4
               IF (ARC(CASE).EQ.ROW-PONTLV) ARC(CASE)=0
   52       CONTINUE
   53    CONTINUE
   36    CONTINUE
         LEVNR(LEV)=PONT-LEVPT(LEV)
   37 CONTINUE
C
   91 CONTINUE
      NWKS=NLWKS(1)
      NROWS=PONT
      NROWS4=NROWS*4
      NROWOC=0
      DO 38 LEV=LEVFRM,NLEVS
         NROWOC=NROWOC+LEVNR(LEV)
   38 CONTINUE
      NROW4O=NROWOC*4
      RETURN
      END
