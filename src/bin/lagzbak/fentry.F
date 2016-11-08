cges  PROGRAM LGZBK
      subroutine fentry (cc,ic,maxcor)
C
C >>>   LAST UPDATE: november 22, 1989   <<<
C
C >> WRITTEN BY GUSTAVO E. SCUSERIA, SEPTEMBER 27, 1988.
C         UPDATE: NOVEMBER 21, 1988.
C
C >> THIS CODE REQUIRES:
C    1-E MO INTS IN FILE79 AS WRITTEN BY ccsr
C    FOCK MATRIX IN FILE79 AS WRITTEN BY cceg
C
C >> FILE79 CONTENT:                                       INTOWP(...)
C    1. 1E-MO INTS IN CC ORDERING (WRITTEN BY ccsr)          NT2
C    2. FOCK MATRIX   CC ORDERING (WRITTEN BY cceg)          NT2
C    3. TOTAL CC ENERGY           (WRITTEN BY cceg)          1
c    :  written by lagzbak:
C    4. TOTAL ONEPDM:  Q + Q (Z) IN SCF ORDER  MO BASIS      NT2
C    5. TOTAL LAGRANGIAN:  X + X(Z) (THIS IS SYMMETRIC!)     NT*NT
C    6. TOTAL ONEPDM IN AO BASIS                             NTAO2
C    7. TOTAL (SYMMETRIC) LAGRANGIAN IN AO BASIS             NTAO2
C    8. Z VECTOR                                             NO*NV
C
C
C >> LAGRANGIAN:
C    THIS CODE READS FI78 (MO INTS) AND FI68 (1,2PDM) WITH CC ORDERING
C    LABELS. USING CCTSCF THESE ARE CONVERTED TO SCF ORDERING.
C    IF CORE IS ENOUGH LAGRANGIAN IS CONSTRUCTED IN CORE.
C    IF NOT, SORTS TO FILE82 AND FILE83.
C    NOTE THAT LAGRANGIAN IS OBTAINED IN SCF ORDERING.
C
C >> Z-CPHF:
C    REQUIRES FILE79 (ONE-E MO INTS & FOCK MATRIX) IN CC ORDERING.
C    CALCULATION IS PERFORMED IN CC ORDERING. LAGRANGIAN NEEDS TO BE
C    CONVERTED FROM SCF TO CC ORDERING.
C
C >> BACK TRANSFORMATION:
C    IS DONE FROM THE MO TO THE SO BASIS FIRST, USING SYMMETRY.
C    THE SO TO AO TRANSFORMATION IS ONLY N**4.
C    IF MEMORY IS ENOUGH BACK TRANSFORMATION IS DONE IN CORE.
C    IF NOT, 2 SORTS ARE REQUIRED. AUXILIARY FILES ARE 95, 96, 97 & 98.
C
c-------------------------------------------------------------------------
c     Notice of Program Modification
c
c     by Roger Grev (crsg)
c     fix gradient bug associated with frozen core orbitals (core and vir)
c     changes are in rdcore.f, sort.f, and lagr.f
c     July 25, 1994
c
c     by Edward Seidl (cets030491)
c     add new style input
c     march 4, 1991
c
c     by Ges, fix bug. so --> ao back transformation needs to be
c     done even in c1 if using 5d, 7f or 9g !
c     july 20, 1990    
c
c     by Ges, fix bug in dummy cycles (AOSTRT) and add >>> <<<
c     november 22, 1989
c
c     by Ges, rclose 4  itaps 82  83  68  78
c     november 21, 1989
c
c     by Ges, skip so --> ao back transformation if symm is c1  
c     november 6, 1989
c
c     by Ges, change escf to eref in lagr (using info in tape69)
c     october 27, 1989
c-------------------------------------------------------------------------
c
      IMPLICIT INTEGER (A-Z)
      REAL*8 CC(maxcor)
      INTEGER IC(*)
      REAL*8 TIMLIM,ENUC,ESCF,ECORR,ETOT69
      REAL*4 EECPU,EETIME,EETCPU
      CHARACTER*32 DATTIM
      COMMON/TIEMPO/EECPU,EETIME,EETCPU
      CHARACTER*4 LABEL(20),OPTION,PRTS,RSTR,CHAR,SRTOUT
      INTEGER NLAMDA(8),ITEMP(100),NC(8),ITYP(8),FLOV(8*2),
     .        NFZO(8),NFZV(8),NDOC(8),NUOC(8),NPQ(8),NRS(8)
C
C
      call psinit('APPEND')
      CALL TSTART(6)
C      CALL DATETM(DATTIM,EECPU,EETIME,EETCPU)
C
      WRITE(6,6000)
 6000 FORMAT(//,49('-'),/,49('-'),/,
     .'                    LAGZBAK                          ',/,
     .'LAGRANGIAN , Z-CPHF AND 2PDM BACK TRANSFORMATION',/,
     .'         WRITTEN BY GUSTAVO E. SCUSERIA  ',/,
     .'          LAGZBAK VS 0.1  JUN 14, 1989',/,49('-'),/,49('-'),//)
C
      INPUT=5
      JOUT=6
      ITAP30=30
      ITAP55=55
      ITAP68=68
      ITAP78=78
      ITAP79=79
      ITAP82=82
      ITAP83=83
      ITAP95=95
      ITAP96=96
      ITAP97=97
      ITAP98=98
C
      CALL RFILE(ITAP55)
      CALL RFILE(ITAP68)
      CALL RFILE(ITAP78)
      CALL RFILE(ITAP79)
C
C >>> BUFFER SIZES
C
      DMBUF =SEC2I(100)/INTOWP(1)
      LENBUF=SEC2I(100)/INTOWP(1)
      WRITE (JOUT,*) ' DMBUF IN REAL WORDS =',DMBUF
      CALL SREAD(ITAP78,INTBUF,1)
      WRITE (JOUT,*) 'INTBUF IN REAL WORDS =',INTBUF
      WRITE (JOUT,*) '  '
      MAXBUF=MAX0(DMBUF,INTBUF)
C
C >>> GET TAPE30 INFO, CONSTRUCT ORBSYM AND FLOV ARRAYS
C
      CALL RFILE(ITAP30)
      CALL WREADW(ITAP30,IC,200,101,JEND)
      MPOINT = IC(2)
      MCONST = IC(3)
      MCALCS = IC(4)
      NCALCS = IC(5)
      NT     = IC(18)
      NAT    = IC(19)
      NTAO   = IC(22)
      NSHELL = IC(27)
      NIRRED = IC(28)
      NSYMHF = IC(41)
      MXCOEF = IC(42)
C
      WRITE(6,*)
      WRITE(6,*)' FILE 30  PARAMETERS'
      WRITE(6,*)
      WRITE(6,*)' NT     = ',NT
      WRITE(6,*)' NTAO   = ',NTAO
      WRITE(6,*)' NSYMHF = ',NSYMHF
      WRITE(6,*)' NIRRED = ',NIRRED
C
C >>> ZERO OUT THE CORE
C
      CALL ZERO(CC,MAXCOR)
C
C >>> CORE ALLOCATION
C
      WTEMP = 1
      COORD = WTEMP + NT
      ZAN   = COORD + 3*NAT
      EIG   = ZAN + NAT
      XVEC  = EIG + NT
      DVEC  = XVEC + NTAO*NT
      WVEC  = DVEC + NT*NT
      TOP   = WVEC + NTAO*NT
      SCFTCC = WPADTI(TOP)
      CCTSCF= SCFTCC + NT
      ORBSYM= CCTSCF + NT
      KATOM = ORBSYM + NT
      KLOC  = KATOM + NSHELL
      KTYPE = KLOC  + NSHELL
      ILAST = KTYPE + NSHELL
C
C >>> READ FILE30
C
      CALL NRED30(ITAP30,MPOINT,MCONST,MCALCS,NCALCS,NT,NTAO,NSYMHF,NAT,
     .            MXCOEF,CC(EIG),NLAMDA,ITEMP,6,NC,NO,
     .            CC(WTEMP),IC(SCFTCC),ITYP,IC(ORBSYM),FLOV,NIRRED,
     .            ENUC,ESCF,CHAR,NORD,CC(ZAN),CC(COORD),IC(CCTSCF),
     .            IC(KATOM),IC(KTYPE),IC(KLOC),NSHELL,
     .            CC(XVEC),CC(DVEC),CC(WVEC),ic1)
C
      NV=NT-NO
      NM=NO*NV
      NT2=NT*(NT+1)/2
      NTAO2=NTAO*(NTAO+1)/2
      NOV=MAX0(NO,NV)
      NV2=NV*(NV+1)/2
      NO2=NO*(NO+1)/2
      NST=NSHELL*(NSHELL+1)/2
C
C >>> READ CCSD INPUT AND IPRINT
C
      CALL INDATA(OPTION,NFZO,NFZV,LABEL,NIRRED,IPRINT)
C
C >>> WRITE OUT INPUT INFO
C
      WRITE(JOUT,6003)OPTION,IPRINT
 6003 FORMAT(//,2X,'>>>  CCSD PARAMETERS   <<<',/,
     .        /,2X,'OPTION =',1X,A4,
     .        /,2X,'IPRINT =',I5)
C
C >>> MORE CORE ALLOCATION
C
      FZT   = ILAST
      IPQ   = FZT + NT
      ITR   = IPQ + NTAO2
      ITV   = ITR + NOV
      IFZO  = ITV + NOV
      IFZV  = IFZO + NO*NO
      AOFF  = IFZV + NV*NV
      AADD  = AOFF + NO*NO
      COFF  = AADD + NO*NO
      CADD  = COFF + NV*NV
      DOFF  = CADD + NO*NO
      DADD  = DOFF + NO*NV
      EOFF  = DADD + NO*NV
      EADD  = EOFF + NO*NV
      IOFF  = EADD + NO*NO
      IADD  = IOFF  + NTAO2
      I1OFF = IADD  + NTAO2
      I1ADD = I1OFF + NT2
      C2OFF = I1ADD + NT2
      C2ADD = C2OFF + NT2
      D2OFF = C2ADD + NT2
      D2ADD = D2OFF + NT2
      E3OFF = D2ADD + NT2
      E3ADD = E3OFF + NT2
      JOFF  = E3ADD + NT2
      JADD  = JOFF  + NT2
      AOPS  = JADD  + NT2
      LNSR  = AOPS  + NSHELL
      KLOFF = LNSR  + NST
      IJOFF = KLOFF + NST
      MAPIJ = IJOFF + NST
      ILAST = MAPIJ + NTAO2
C
C >>> SET FROZEN ARRAYS AND ACTIVE ORBITALS PER SYMMETRY
C
      CALL FROZEN (FLOV,NLAMDA,NC,NFZO,NFZV,NDOC,NUOC,IC(FZT),
     .             ITYP,NIRRED,NO,NV,JOUT,NOACT,NVACT,NSYMHF)
      LFZO=NO-NOACT
      LFZV=NV-NVACT
C
C >>> SET ADDRESSING ARRAYS
C
      CALL SYMCNT(IC(ORBSYM),FLOV,NIRRED,NO,NV,NT,IC(FZT),
     .            NSGOO,NSGOV,NSHOV,NSLOV,NSINT,N1INT,NC2INT,ND2INT,
     .            IC(IPQ),IC(ITR),IC(ITV),IC(CCTSCF),IC(AOFF),IC(AADD),
     .            IC(COFF),IC(CADD),IC(DOFF),IC(DADD),IC(EOFF),IC(EADD),
     .            NT2,NTAO2,NOV,NV2,NO2,NM,IC(IOFF),IC(IADD),IC(IFZO),
     .            IC(IFZV),IC(I1OFF),IC(I1ADD),IC(C2OFF),IC(C2ADD),
     .            IC(D2OFF),IC(D2ADD),IC(E3OFF),IC(E3ADD),NE3INT)
C
      WRITE(6,*)'       '
      WRITE(6,*)' NSGOO=',NSGOO
      WRITE(6,*)' NSGOV=',NSGOV
      WRITE(6,*)' NSHOV=',NSHOV
      WRITE(6,*)' NSLOV=',NSLOV
      WRITE(6,*)' NSINT=',NSINT
      WRITE(6,*)' N1INT=',N1INT
      WRITE(6,*)' N2INT=',NC2INT+ND2INT
      WRITE(6,*)' N3INT=',NE3INT
      WRITE(6,*)'       '
C
c  note that proper equivalence buf and ibuf are as follows
c   ibuf= ...  & buf=iadtwp(ibuf) will not work because iadtwp
c   gives the NEXT available real word. (I+3)/2.
c   when ibuf is odd this will not work.
c
      BUF   = IADTWP(ILAST)
      IBUF  = WPADTI(BUF)
      XLAG  = BUF + MAXBUF
      HI    = XLAG + NT*NT
      QDM   = HI + NT2
      QDMT  = QDM + NT2
      TOP   = QDMT + NT2
C
C >>> READ TAPE79 JUST TO GET THE TOTAL ENERGY
C
      IADR=1+INTOWP(NT2)+INTOWP(NT2)
      CALL WREADW(ITAP79,ETOT69,INTOWP(1),IADR,JUNK)
C
      DELX  = TOP
      ZETA  = DELX + NO*NV
      ZETO  = ZETA + NO*NV
      ZETV  = ZETO + NO*LFZO
      EX    = ZETV + NV*LFZV
      FOCK  = EX   + NO*NV
      AUXOO = FOCK + NT2
      AUX1  = AUXOO + NO*NO
      AUX2  = AUX1 + NTAO*NTAO
      AUX3  = AUX2 + NTAO*NTAO
      BZET  = AUX3 + NTAO*NTAO
      QMAT  = BZET + NTAO*NTAO
      QCOR  = QMAT + NTAO*NTAO
      QACT  = QCOR + NTAO*NTAO
      TOPINT= QACT + NTAO*NTAO
      AIN   = TOPINT
      CIN   = AIN + NSGOO
      DIN   = CIN + NSGOV
      EIN   = DIN + NSHOV
      X1INT = EIN + NSLOV
      C2INT = X1INT + N1INT
      D2INT = C2INT + NC2INT
      E3INT = D2INT + ND2INT
      TOP   = E3INT + NE3INT
C
      LEFT  =MAXCOR-TOP-1
      LEFTIN=MAXCOR-TOPINT-1
C
      IF(TOP.GE.MAXCOR)THEN
         WRITE(6,*)'NOT ENOUGH CORE, REQUIRES',TOP, '  REAL*8 WORDS'
         call qabort
         call mabort
      ENDIF
C
C >>> DO THE CALCULATION IN CORE IF MEMORY IS ENOUGH
C
      II    = TOP
      GDM   = II + NSINT
      CORE1 = GDM + NSINT
      CORE2 = CORE1+ NT*NT*NT*NT
      TOPCOR= CORE2+ NTAO*NTAO*NT*NT
C
      IF(TOPCOR.GT.MAXCOR)THEN
C
C >>> SORT FILE78 TO 82 AND FILE68 TO 83 (ONLY 4 ACTIVE INDICES)
C
      WRITE(6,*)
      WRITE(6,*)' AVAILABLE MEMORY IS ',MAXCOR,'  REAL*8 WORDS'
      WRITE(6,*)' REQUIRED  MEMORY IS ',TOP   ,'  REAL*8 WORDS'
      WRITE(6,*)' CORE FOR SORTING IS ',LEFT  ,'  REAL*8 WORDS'
      WRITE(6,*)' CORE FOR 2ND SRT IS ',LEFTIN,'  REAL*8 WORDS'
      WRITE(6,*)
C
      CALL SORT(CC(BUF),IC(IBUF),INTBUF,IC(CCTSCF),IC(SCFTCC),NPQ,NRS,
     .          NO,NV,NT,ITAP78,ITAP82,JOUT,
     .          NIRRED,FLOV,IC(ORBSYM),IC(IPQ),IC(ITV),
     .          IC(JOFF),IC(JADD),CC(AUX1),CC(TOP),LEFT,IC(FZT),CC(HI),
     .          IC(I1OFF),IC(I1ADD),CC(X1INT),
     .          IC(C2OFF),IC(C2ADD),CC(C2INT),
     .          IC(D2OFF),IC(D2ADD),CC(D2INT),
     .          IC(E3OFF),IC(E3ADD),CC(E3INT))
      WRITE(6,*)
      WRITE(6,*)'  >>>   SORTING OF INTEGRALS  COMPLETED   <<<'
      WRITE(6,*)
c     CALL RCLOSE(ITAP78,3)
      CALL RCLOSE(ITAP78,4)
C
      CALL SORT(CC(BUF),IC(IBUF), DMBUF,IC(CCTSCF),IC(SCFTCC),NPQ,NRS,
     .          NO,NV,NT,ITAP68,ITAP83,JOUT,
     .          NIRRED,FLOV,IC(ORBSYM),IC(IPQ),IC(ITV),
     .          IC(JOFF),IC(JADD),CC(AUX1),CC(TOP),LEFT,IC(FZT),CC(QDM),
     .          IC(I1OFF),IC(I1ADD),CC(X1INT),
     .          IC(C2OFF),IC(C2ADD),CC(C2INT),
     .          IC(D2OFF),IC(D2ADD),CC(D2INT),
     .          IC(E3OFF),IC(E3ADD),CC(E3INT))
      SRTOUT=' YES'
      WRITE(6,*)
      WRITE(6,*)'  >>>   SORTING OF  TWO-PDM   COMPLETED   <<<'
      WRITE(6,*)
c     CALL RCLOSE(ITAP68,3)
      CALL RCLOSE(ITAP68,4)
C
      ELSE
C
      TOPINT=TOPCOR
      LEFTIN=MAXCOR-TOPINT-1
      WRITE(6,*)' '
      WRITE(6,*)' ENOUGH MEMORY TO CONSTRUCT THE LAGRANGIAN IN CORE'
      WRITE(6,*)' REQUIRED  MEMORY IS ',TOPCOR,'  REAL*8 WORDS'
      WRITE(6,*)' AVAILABLE MEMORY IS ',MAXCOR,'  REAL*8 WORDS'
      WRITE(6,*)' '
C
      CALL RDCORE(ITAP78,CC(II),NSINT,6,CC(BUF),IC(IBUF),INTBUF,
     .            IC(IPQ),NT2,IC(CCTSCF),NO,NT,CC(HI),IC(IOFF),
     .            IC(IADD),IC(FZT),
     .            IC(I1OFF),IC(I1ADD),CC(X1INT),
     .            IC(C2OFF),IC(C2ADD),CC(C2INT),
     .            IC(D2OFF),IC(D2ADD),CC(D2INT),
     .            IC(E3OFF),IC(E3ADD),CC(E3INT))
      CALL RCLOSE(ITAP78,3)
c     CALL RCLOSE(ITAP78,4)
C
      CALL RDCORE(ITAP68,CC(GDM),NSINT,6,CC(BUF),IC(IBUF),DMBUF,
     .            IC(IPQ),NT2,IC(CCTSCF),NO,NT,CC(QDM),IC(IOFF),
     .            IC(IADD),IC(FZT),
     .            IC(I1OFF),IC(I1ADD),CC(X1INT),
     .            IC(C2OFF),IC(C2ADD),CC(C2INT),
     .            IC(D2OFF),IC(D2ADD),CC(D2INT),
     .            IC(E3OFF),IC(E3ADD),CC(E3INT))
      SRTOUT='  NO'
      CALL RCLOSE(ITAP68,3)
c     CALL RCLOSE(ITAP68,4)
C
      ENDIF
C
C
C >>> CALCULATE ACTIVE LAGRANGIAN AND ADD HF CONTRIBUTION
C
      CALL LAGR  (CC(II),CC(HI),CC(GDM),CC(QDM),IC(IPQ),CC(XLAG),
     .            NO,NV,NT,IC(ORBSYM),FLOV,NIRRED,IC(IOFF),IC(IADD),
     .            IC(SCFTCC),IC(CCTSCF),CC(FOCK),ENUC,ESCF,ETOT69,
     .            ITAP79,IC(FZT),CC(ZETO),CC(ZETV),CC(EX),IC(ITR),
     .            LFZO,LFZV,IC(IFZO),IC(IFZV),IC(ITV),NRS,SRTOUT,
     .            ITAP82,ITAP83,CC(QMAT),CC(QCOR),CC(AUX1),CC(AUX2),
     .            IC(I1OFF),IC(I1ADD),CC(X1INT),
     .            IC(C2OFF),IC(C2ADD),CC(C2INT),
     .            IC(D2OFF),IC(D2ADD),CC(D2INT),
     .            CC(AIN),CC(CIN),CC(DIN),CC(EIN),
     .            IC(AOFF),IC(AADD),IC(COFF),IC(CADD),IC(DOFF),IC(DADD),
     .            IC(EOFF),IC(EADD),NSGOO,NSGOV,NSHOV,NSLOV,
     .            CC(BUF),IC(IBUF),LENBUF,IC(E3OFF),IC(E3ADD),CC(E3INT))
C
C >>> SOLVE CPHF
C
      CALL CPHF  (CC(XLAG),CC(CIN),CC(DIN),CC(ZETA),CC(DELX),NO,NV,NT,
     .            IC(COFF),IC(CADD),IC(DOFF),IC(DADD),
     .            IC(ORBSYM),NIRRED,FLOV,
     .            IC(ITR),IC(ITV),IC(IPQ),
     .            CC(BUF),IC(IBUF),LENBUF,NSGOV,NSHOV,IC(CCTSCF),
     .            CC(AUXOO),CC(FOCK),CC(EX))
C
C >>> ADD Z-CPHF AND HARTREE-FOCK CONTRIBUTIONS TO ONE-PDM
C
      CALL ADDZ  (CC(QDMT),CC(QDM),
     .            IC(ORBSYM),NIRRED,FLOV,
     .            IC(ITR),IC(ITV),IC(IPQ),NO,NV,NT,
     .            IC(CCTSCF),IC(SCFTCC),CC(ZETA),
     .            CC(ZETO),CC(ZETV),IC(FZT),IC(IFZO),IC(IFZV))
C
C >>> ADD Z-CPHF CONTRIBUTION TO LAGRANGIAN
C
      CALL ZLAG  (CC(XLAG),CC(FOCK),CC(ZETA),
     .            CC(AIN),CC(CIN),CC(DIN),CC(EIN),
     .            IC(AOFF),IC(AADD),IC(COFF),IC(CADD),IC(DOFF),IC(DADD),
     .            IC(EOFF),IC(EADD),NSGOO,NSGOV,NSHOV,NSLOV,
     .            CC(BUF),IC(IBUF),LENBUF,
     .            IC(ORBSYM),NIRRED,FLOV,
     .            IC(ITR),IC(ITV),IC(IPQ),NO,NV,NT,
     .            IC(CCTSCF),IC(SCFTCC),
     .            CC(ZETO),CC(ZETV),IC(IFZO),IC(IFZV),IC(FZT),
     .            CC(AUX1),NSYMHF,ITYP,NFZO,NFZV)
C
C >>> BACK TRANSFORM ONEPDM AND LAGRANGIAN TO AO BASIS AND WRITE TO 79
C
      CALL BACK1 (CC(QDMT),CC(ZETA),CC(XLAG),CC(WVEC),CC(AUX1),CC(AUX2),
     .            ITAP79,IC(ORBSYM),NIRRED,FLOV,IC(IPQ),NO,NV,NT,NT2,
     .            CC(AUX1),IC(ITR),IC(CCTSCF),NTAO,NTAO2)
C
C >>> BACK TRANSFORM THE 2PDM
C
      CALL SET55 (IC(AOPS),IC(IJOFF),IC(KLOFF),IC(LNSR),IC(MAPIJ),
     .            IC(KTYPE),IC(KLOC),IC(IPQ),
     .            NSHELL,NST,NTAO,NTAO2,ITAP55,AOSTRT,CC(TOPINT),LEFTIN)
C
      CALL BACK2 (CC(DVEC),CC(XVEC),CC(AUX1),CC(AUX2),CC(GDM),
     .            IC(IOFF),IC(IADD),IC(ORBSYM),NIRRED,FLOV,
     .            IC(ITR),IC(ITV),IC(IPQ),NO,NV,NT,NT2,NTAO,NTAO2,
     .            IC(AOPS),IC(IJOFF),IC(KLOFF),IC(LNSR),IC(MAPIJ),
     .            IC(KTYPE),IC(KLOC),NSHELL,ITAP55,AOSTRT,
     .            IC(CCTSCF),IC(SCFTCC),CC(ZETA),CC(CORE1),CC(CORE2),
     .            CC(WVEC),CC(BZET),CC(QMAT),IC(FZT),
     .            CC(QCOR),CC(QACT),CC(QDM),NFZO,NFZV,
     .            CC(ZETO),CC(ZETV),IC(IFZO),IC(IFZV),ITAP83,SRTOUT,
     .            NRS,CC(AUX3),CC(TOPINT),LEFTIN,ic1)
C
      CALL RCLOSE(ITAP30,3)
      CALL RCLOSE(ITAP55,3)
      CALL RCLOSE(ITAP79,3)
c
      write(6,*)
      write(6,*)'    >>>   Lagrangian  Calculation  Completed   <<< '
      write(6,*)
      CALL TSTOP(6)
      call psdone()
      stop
C
      END
