cges  PROGRAM CCEG
c     program cceg
c     real*8 cc(6500000)
c     integer ic(1)
c     equivalence (cc(1),ic(1))
c     call fentry(cc,ic,6500000)
c     end
      subroutine fentry(cc,ic,mxcor)
C
C   >>> LAST UPDATE jan 26, 1990  <<<
C
C   ******************************************************
C   **  THE NEW CLOSED-SHELL COUPLED  CLUSTER  PROGRAM  **
C   **             JULY  17  , 1988  BY GES             **
C   **      CISD IS 1/4 O2 V4 + 2 O3 V3 + 1/4 O4 V2     **
C   **      CCSD IS 1/4 O2 V4 + 4 O3 V3 + 1/2 O4 V2     **
C   **      VERSION 9, FINAL FORM SEP 2, 1988           **
C   ******************************************************
C
C     NOTES:
C    > PROGRAM COMPATIBLE ONLY WITH CCSR DUE TO NEW ORBSYM ORDER.
C    > FILE79 IS WRITTEN (USED BY LAGZBAK).
C
C
C***********************************************************************
C  NOTICE OF PROGRAM MODIFICATIONS
C***********************************************************************
C      BY:  David Sherrill                        search: cds2-11-93
C    date:  February 11, 1993
C  reason:  update to use new-style input
C
C***********************************************************************
C      BY:  GUSTAVO E. SCUSERIA
c
c    date:  february 25, 1991
c  reason:  add reference to output and beautify eigenvalue output
c
c    date:  october 15, 1990
c  reason:  if option=FSDT and level=0 stop calculation.
c           fix problem in may 23 modification.
c           add energies to energy.dat
c
c    date:  june 27, 1990
c  reason:  speed-up ncount to get #t3 walks; use option FSDTCHEK
c
c    date:  may 23, 1990
c  reason:  add all non-scf terms to full ccsdt    
c
c    date:  april 27, 1990
c  reason:  when printing # walks including triples
c           change ndimt2 to nshov
c
c    date:  january 31, 1990
c  reason:  implement full ccsdt                  
c
c    date:  december 6, 1989
c  reason:  implement t3 in core option           
c
c    date:  november 8, 1989
c  reason:  change t1 diagnostic to Tim's version.
c
c    date:  october 25, 1989
c  reason:  include ccd option.
c
c    date:  september 11, 1989
c  reason:  create term.dat with real time info.
c
c    date:  september 8, 1989
c  reason:  file81 is required by ccqg. include option SDT1GRAD to create it
c
C    date:  august 31,1989.
C  reason:  get rid of file81 in ccsdt1 subroutine.                     
C
C    DATE:  JUNE 16, 1989.
C  REASON:  FIX FOCK2 AND RDINTS ROUTINES FOR NON-SCF AND FROZEN CALCS
C           If frozen orbitals are present, frozen part of the fock
c              matrix is constructed while reading the 2-e ints.
c           Only the active 2e ints are kept in core although all of
c              them are read and needed for the Fock matrix.
C
C    DATE:  MAY 10, 1989.
C  REASON:  CHANGE NAME TO CCEG
C
C    DATE:  MAY 1, 1989.
C  REASON:  IMPLEMENT TWRIT FOR CISD AND CCSD.
C
C    DATE:  APRIL 20, 1989.
C  REASON:  BUG FOUND! ITR SHOULD BE DIMENSION MAX(NO,NV)
C           (THIS WAS CAUSING PROBLEMS FOR NO.GT.NV )
C
C    DATE:  MARCH 31, 1989.
C  REASON:  CHANGE TDIFF TO MINIMAL SET AND INCLUDE TWRIT
C
C    DATE:  MARCH 27, 1989.
C  REASON:  FINISH CLEANUP OF READ30, NCOUNT, SYMARR, FROZEN & MAIN
C           GET RID OF RFI30 ROUTINE
C           CHANGE ORDERING OF ORBITAL SYMMETRIES (DONE IN NCCSORT TOO)
C
C    DATE:  MARCH 26, 1989.
C  REASON:  SPLIT CLUS9 INTO SUBROUTINES
C           CREATE SUBROUTINE CISD9
C
C    DATE:  MARCH 22, 1989.
C  REASON:  GET RID OF F*T1 & E*T1 ---> T2 TERMS
C           FIX II POINTER IN CALL CLUS & CALL PLUS
C           FIX DIMENSION OF ITEMP TO 100
C           GET RID OF SUBROUTINES RINTS AND PAK-UNPAK (FILE78)
C           FIX DIMENSIONS OF AUX ARRAYS IN CLUS
C
C    DATE:  NOVEMBER 16, 1988
C  REASON:  CALCULATE FOCK OPERATOR IN MO BASIS READING ONE-E MO INTS
C           FROM FILE79. WRITE OUT FOCK MATRIX TO FILE79.
C
C    DATE:  OCTOBER 3, 1988
C  REASON:  INCLUDE REFORMULATED CCSDT-1A OPTION.
C
C    DATE:  SEPTEMBER 16, 1988
C  REASON:  FORCE NBFAO=NBF FOR NEW GCINTS CODE.
C
C***********************************************************************
C
      IMPLICIT INTEGER (A-Z)
      REAL*8 TIMLIM,CC,ENUC,ESCF,ECORR,ETOT
      REAL*4 EECPU,EETIME,EETCPU
      CHARACTER*32 DATTIM
      CHARACTER*4 LABEL(20),CASE,OPTION,RSTR,CHAR,ITYP(8),fl81,t3core
      COMMON/TIEMPO/EECPU,EETIME,EETCPU
      COMMON/RESTAR/TIMLIM
      COMMON/PARAM/ ESCF,ENUC,ECORR,
     .              EPSI,NGO,NMIN ,MINDIM,MAXDIM,MAXIT,CONVI,
     .              NO,NV,NT,NOV,NM,NV2,NT2,IMX,NONO,NVNV,NDIMT2,NTIN,
     .              NSYMHF,NIRRED,MAXOO,MAXVV,MAXOV,
     .              LEVEL,INTBUF,DIISFL,ITAP69,ITAP71,ITAP98,ITAP99,RSTR
      COMMON/INTS/NSGOO,NSGVV,NSGOV,NSHOV,NSLOV,NSLVO
      INTEGER NLAMDA(8),NCL(8),FLOV(8,4),NOP(8),
     .        NOSYM(8),NVSYM(8),NFZO(8),NFZV(8),NDOC(8),NUOC(8),NT3(8)
      DIMENSION CC(mxcor),IC(*)
cges  EQUIVALENCE (CC,IC)
C
c     CALL NOUNFL
c     mxcor=1500000
C
      call psinit('APPEND')
      call ffile(4,'term',0)
      call ffile(10,'energy',1)
      write(10,'(a1)')'*'
C      CALL DATETM(DATTIM,EECPU,EETIME,EETCPU)
      CALL TSTART(6)
C
 6000 FORMAT(//,42('-'),/,42('-'),/,
     .' the full ccsdt coupled cluster program',/,
     .'       written by Gustavo Scuseria       ',/,
     .'        vs fcceg 1.0 feb  02, 1990     ',/,42('-'),/,42('-'),///,
     .'  Reference:',/,
     .'> full ccsdt energy closed-shell',/,
     .'  G.E. Scuseria and H.F. Schaefer,',/,
     .'  Chem. Phys. Lett. 152, 382 (1988). ',//)
C
      INPUT=5
      IW = 6
      ITAP30=30
      ITAP60=60
      ITAP61=61
      ITAP62=62
      ITAP63=63
      ITAP64=64
      ITAP65=65
      ITAP66=66
      ITAP69=69
      ITAP71=71
      ITAP78=78
      ITAP79=79
      ITAP90=90
      ITAP98=98
      ITAP99=99
      WRITE(IW,6000)
      INTBUF=SEC2I(100)/INTOWP(1)
      CALL RFILE(ITAP78)
      CALL RFILE(ITAP60)
      CALL RFILE(ITAP61)
      CALL RFILE(ITAP62)
      CALL RFILE(ITAP63)
      CALL RFILE(ITAP64)
      CALL RFILE(ITAP65)
      CALL RFILE(ITAP66)
      CALL RFILE(ITAP79)
      CALL RFILE(ITAP98)
      CALL RFILE(ITAP99)
      CALL RFILE(ITAP90)
C
C >>> READ FILE30
C
      CALL RFILE(ITAP30)
      CALL WREADW(ITAP30,IC,200,101,LWORD)
      MPOINT = IC(2)
      MCONST = IC(3)
      MCALCS = IC(4)
      NCALCS = IC(5)
      NT     = IC(18)
      NATOM  = IC(19)
      NTAO   = IC(22)
      NSHELL = IC(27)
      NIRRED = IC(28)
      NPRIM  = IC(32)
      NSYMHF = IC(41)
      MXCOEF = IC(42)
      IOPEN  = IC(43)
C
      WRITE(6,*)
      WRITE(6,*)' >>>  PARAMETERS FROM FILE30  <<<'
      WRITE(6,*)
      WRITE (6,*) ' MPOINT=',MPOINT
      WRITE (6,*) ' MCONST=',MCONST
      WRITE (6,*) ' MCALCS=',MCALCS
      WRITE (6,*) ' NCALCS=',NCALCS
      WRITE (6,*) ' NT    =',NT
      WRITE (6,*) ' NATOM =',NATOM
      WRITE (6,*) ' NTAO  =',NTAO
      WRITE (6,*) ' NSHELL=',NSHELL
      WRITE (6,*) ' NIRRED=',NIRRED
      WRITE (6,*) ' NPRIM =',NPRIM
      WRITE (6,*) ' NSYMHF=',NSYMHF
      WRITE (6,*) ' MXCOEF=',MXCOEF
      WRITE (6,*) ' IOPEN =',IOPEN
C
C >>> ZERO OUT THE CORE
C
      CALL ZERO(CC,MXCOR)
C
C >>> DYNAMICAL ALLOCATION
C
      WTEMP=1
      COORD=WTEMP+NT
      ZAN=COORD+3*NATOM
      EEE=ZAN+NATOM
      SCFTCC=WPADTI(EEE+NT)
      ORBSYM=SCFTCC+NT
      FZO=ORBSYM+NT
C
C >>> READ FILE30
C
      CALL READ30(ITAP30,MPOINT,MCONST,MCALCS,NCALCS,NT,
     .            NSYMHF,MXCOEF,CC(EEE),NLAMDA,6,NCL,NO,
     .            CC(WTEMP),IC(SCFTCC),ITYP,IC(ORBSYM),FLOV,NIRRED,
     .            ENUC,ESCF,NATOM,NSHELL,NPRIM,
     .            CHAR,NORD,CC(ZAN),CC(COORD),
     .            IOPEN,NOP)
C
      NV=NT-NO
C
C >>> MORE DYNAMICAL ALLOCATION
C
      FZV=FZO+NO
      TOP=FZV+NV
C
C
C >>> READ INPUT FILE
cds2-11-93 new function reads from old OR new input styles
      call rdinpt(epsi,ngo,nmin,mindim,maxdim,maxit,convi,level,
     &            ncsf,fldiis,option,fl81,label,rstr,nfzo,nfzv,
     &            nirred)

cds2-11-93      CALL LOCATE (INPUT,'# CCSD ###',IERR)
cds2-11-93      IF(IERR.NE.0) THEN
cds2-11-93         WRITE (IW,6007)
cds2-11-93         WRITE ( 4,6007)
cds2-11-93         EPSI=5
cds2-11-93         NGO=2
cds2-11-93         NMIN =2
cds2-11-93         MINDIM=2
cds2-11-93         MAXDIM=8
cds2-11-93         MAXIT=20
cds2-11-93         CONVI=7
cds2-11-93         LEVEL=0
cds2-11-93         RSTR='NO  '
cds2-11-93         CASE='CS  '
cds2-11-93         NCSF=20
cds2-11-93         DO 622 I=1,NIRRED
cds2-11-93            NFZO(I)=0
cds2-11-93            NFZV(I)=0
cds2-11-93  622    CONTINUE
cds2-11-93      ELSE
cds2-11-93         READ(5,6004)LABEL
cds2-11-93         WRITE(6,*)
cds2-11-93         WRITE(6,*)
cds2-11-93         WRITE(6,*)' >>>  CCSD INPUT LABEL  <<<'
cds2-11-93         WRITE(6,*)
cds2-11-93         WRITE(6,6024)LABEL
cds2-11-93         READ(5,6002)EPSI,NGO,NMIN ,MINDIM,MAXDIM,FLDIIS,IRSTR,NCSF
cds2-11-93         IF(EPSI  .EQ.0)EPSI  =5
cds2-11-93         IF(NGO   .EQ.0)NGO   =2
cds2-11-93         IF(NMIN  .EQ.0)NMIN  =2
cds2-11-93         IF(MINDIM.EQ.0)MINDIM=2
cds2-11-93         IF(MAXDIM.EQ.0)MAXDIM=8
cds2-11-93         IF(NCSF  .EQ.0)NCSF  =20
cds2-11-93         IF(IRSTR.EQ.1)THEN
cds2-11-93            RSTR='YES '
cds2-11-93         ELSE
cds2-11-93            RSTR='NO '
cds2-11-93         ENDIF
cds2-11-93         READ(5,6002)CONVI,MAXIT,LEVEL
cds2-11-93         IF(CONVI .EQ.0)CONVI =7
cds2-11-93         IF(MAXIT .EQ.0)MAXIT =20
cds2-11-93         DO 623 I=1,NIRRED
cds2-11-93            READ(INPUT,6051)NFZO(I),NFZV(I)
cds2-11-93  623    CONTINUE
cds2-11-93         READ(5,6004)OPTION,fl81
cds2-11-93      ENDIF
cds2-11-93      if (option.ne.'FSDT') then
cds2-11-93        write(6,*) 'ERROR: Use cceg instead of fcceg.'
cds2-11-93        call qabort
cds2-11-93        endif
cds2-11-93      if(level.eq.0.and.option.eq.'FSDT')then
cds2-11-93         write(6,*)'you are using option=FSDT and level=0'
cds2-11-93         write(6,*)'either set level=1-5 or option=CCSD'  
cds2-11-93         call qabort
cds2-11-93      endif
cds2-11-93C
cds2-11-93 6002 FORMAT(16I5)
cds2-11-93 6004 FORMAT(20A4)
cds2-11-93 6006 FORMAT(A4,1X,2I5)
cds2-11-93  6007 FORMAT(//,2X,' WARNING! NO INPUT TO CCSD HAS BEEN FOUND.',/,
cds2-11-93      .          2X,'          DEFAULT PARAMETERS WILL BE USED.')
cds2-11-93  6024 FORMAT(2X,20A4)
cds2-11-93  6051 FORMAT(I2,1X,I2)

C
C >>> WRITE OUT INPUT INFO
C
      WRITE(IW,6003)EPSI,NGO,NMIN ,MINDIM,MAXDIM,FLDIIS,
     .              CONVI,MAXIT,LEVEL,OPTION,fl81,RSTR,NCSF
 6003 FORMAT(//,2X,'>>>  DIIS PARAMETERS  <<< ',/,
     .        /,2X,'EPSI   =',I5,
     .        /,2X,'NGO    =',I5,
     .        /,2X,'NMIN   =',I5,
     .        /,2X,'MINDIM =',I5,
     .        /,2X,'MAXDIM =',I5,
     .        /,2X,'FLDIIS =',I5,
     .      ///,2X,'>>>  CCSD PARAMETERS   <<<',/,
     .        /,2X,'CONVI  =',I5,
     .        /,2X,'MAXIT  =',I5,
     .        /,2X,'LEVEL  =',I5,
     .        /,2X,'OPTION =',1X,A4,
     .        /,2X,'fl81   =',1X,A4,
     .        /,2X,'RSTR   =',1X,A4,
     .        /,2X,'NCSF   =',I5)
  749 CONTINUE
C
C >>> SET FROZEN ARRAYS AND ACTIVE ORBITALS PER SYMMETRY
C
      CALL FROZEN (FLOV,NOSYM,NVSYM,NFZO,NFZV,NDOC,NUOC,IC(FZO),IC(FZV),
     .             ITYP,NIRRED,NO,NV,IW,NOACT,NVACT,NSYMHF)
C
C >>> COUNT INTEGRALS AND VECTORS ONLY FOR ALLOCATION PURPOSES
C
      CALL NCOUNT(IC(ORBSYM),FLOV,NIRRED,NO,NV,NT,IC(FZO),IC(FZV),
     .            NDIMT1,NDIMT2,NSGOO,NSGVV,NSGOV,NSHOV,NSLOV,NSLVO,
     .            NT3,NTAU,OPTION,MAXOO,MAXVV,MAXOV,DIMAR,NDIMT3,
     .            MAXO,MAXV)
C
C >>> CALCULATE BUNCH OF USEFUL PARAMETERS
C
      NO2=(NO*(NO+1))/2
      NV2=(NV*(NV+1))/2
      NT2=(NT*(NT+1))/2
      NONO=NO*NO
      NVNV=NV*NV
      NM =NO*NV
      NMNM=NM*NM
      NTNT=NT*NT
      NM2=(NM*(NM+1))/2
      NGOO=(NO2*(NO2+1))/2
      NGVV=(NV2*(NV2+1))/2
      NGOV=(NM*(NO+1)*(NV+1))/4
      NHOV=(NM*(NM+1))/2
      NLOV=NO2*NM
      NLVO=NV2*NM
      NINT=NGOO+NGVV+NGOV+NHOV+NLOV+NLVO
      NSINT=NSGOO+NSGVV+NSGOV+NSHOV+NSLOV+NSLVO
C     MAXAR=NSGOO
C     IF(NSGVV.GT.MAXAR)MAXAR=NSGVV
C     IF(NSGOV.GT.MAXAR)MAXAR=NSGOV
C     IF(NSHOV.GT.MAXAR)MAXAR=NSHOV
C     IF(NSLOV.GT.MAXAR)MAXAR=NSLOV
C     IF(NSLVO.GT.MAXAR)MAXAR=NSLVO
C     NSTIN= MAXAR
      NSTIN=NSGOO+      NSGOV+NSHOV+NSLOV+NSLVO
      NONO2=NO*NO2
      NONV2=NO*NV2
      NONM=NO*NM
      NVNO2=NV*NO2
      NVNV2=NV*NV2
      NVNM=NV*NM
      NOV2=MAX0(NO2,NV2)
      NOV =MAX0(NO ,NV )
      IMX=MAX0(NM,NV2)
      NOVNOV=NOV*NOV
C
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)' >>>  NUMBER OF INTEGRALS   <<<'
      WRITE(6,*)
      WRITE(6,*)' NGOO =',NGOO,'   NSGOO =',NSGOO
      WRITE(6,*)' NGVV =',NGVV,'   NSGVV =',NSGVV
      WRITE(6,*)' NGOV =',NGOV,'   NSGOV =',NSGOV
      WRITE(6,*)' NHOV =',NHOV,'   NSHOV =',NSHOV
      WRITE(6,*)' NLOV =',NLOV,'   NSLOV =',NSLOV
      WRITE(6,*)' NLVO =',NLVO,'   NSLVO =',NSLVO
      WRITE(6,*)' NINT =',NINT,'   NSINT =',NSINT
      WRITE(6,*)
      WRITE(6,*)' REQUIRED CORE FOR INTEGRALS',NSTIN
C
      WRITE(IW,6042)
      WRITE(IW,6013)NDIMT1
      WRITE(IW,6014)NSHOV
      NDIMTT=1+NDIMT1+NSHOV
      IF(OPTION.EQ.'SDT1'.or.option.eq.'FSDT')THEN
         WRITE(IW,6016)NDIMT3
         NDIMTT=1+NDIMT1+nshov+NDIMT3
      ENDIF
      WRITE(IW,6015)NDIMTT
 6042 FORMAT(//,1X,' >>>  NUMBER OF WALKS   <<<<<<<<<<<',/)
 6013 FORMAT(2X,'THE NUMBER OF SINGLE EXCITATIONS IS',I12)
 6014 FORMAT(2X,'THE NUMBER OF DOUBLE EXCITATIONS IS',I12)
 6016 FORMAT(2X,'THE NUMBER OF TRIPLE EXCITATIONS IS',I12)
 6015 FORMAT(2X,'THE TOTAL NUMBER OF S & D WALKS  IS',I12)
C
C >>> FINISH CORE ALLOCATION
C
      IPV=TOP
      IPO=IPV+NM
      ITRV=IPO+NM
      ITRT=ITRV+NM
      ITR=ITRT+NM
      IPL=ITR+NOV
      IPQ=IPL+NT
      ITRO=IPQ+NT2
      ITV=ITRO+IMX
      UOFF=ITV+NV
      VADD=UOFF+NO*NO*2
      ZLX= VADD+NV*NV
      AOFF=ZLX+NV*NV
      AADD=AOFF+NONO
      BOFF=AADD+NONO
      BADD=BOFF+NVNV
      COFF=BADD+NVNV
      CADD=COFF+NVNV
      DOFF=CADD+NONO
      DADD=DOFF+NM
      EOFF=DADD+NM
      EADD=EOFF+NM
      FOFF=EADD+NONO
      FADD=FOFF+NVNV
      IVAL=FADD+NM
      ININ=IVAL+NCSF
c
      IF(OPTION.EQ.'SDT1'.or.option.eq.'FSDT')THEN
         T3OFF=ININ
         T3ADD=T3OFF+NO*NO*2*NIRRED
         i1off=T3ADD+NV*NV*NIRRED
         i2off=i1off+no*nv
         i3off=i2off+no*nv*nirred
         trsym=i3off+no*nv
         inin =trsym+no*nv
      ELSE
         T3OFF=1
         T3ADD=1
         i1off=1
         i2off=1
         i3off=1
      ENDIF
c
      AIN=IADTWP(ININ)
      BIN=AIN+NSGOO
      CIN=BIN+NV*NV
      DIN=CIN+NSGOV
      EIN=DIN+NSHOV
      FIN=EIN+NSLOV
      BUF=FIN+NSLVO
      IBUF=WPADTI(BUF)
      TOPBUF=BUF+INTBUF
      T2O=TOPBUF
      T2=T2O+NDIMT2
      FOCK=T2+NDIMT2
      HH  =FOCK+NT2
      ICC =HH  +NT2
      IBB=ICC+MAXDIM
      IB2=IBB+(MAXDIM+1)*(MAXDIM+2)
      TVAL=IB2+(MAXDIM+1)*(MAXDIM+2)
      ILAST=TVAL+NCSF
C
C     ALLOCATE CORE FOR NEW CCSD INTERMEDIATES
C
      AR1=ILAST
      AR3=AR1
      AR5=AR1
      AR2=AR1+DIMAR
      AR4=AR2
      AR6=AR2
      GO=AR2+DIMAR
      GV=GO+NO*NO
      HO =GV+NV*NV
      HIA= HO+NO*NO
      HV =HIA+NO*NV
      T1O= HV+NV*NV
      T1=T1O+NO*NV
      FFO=T1+NO*NV
      SFOO=FFO
      FFV=FFO+NO*NO
      SFVV=FFV
      AUX1=FFV+NV*NV
      AUX2=AUX1+NO*NV
      AUX3=AUX2+NO*NO
      AUX4=AUX3+NV*NV
      JAR =AUX4+NV*NO
      KAR =JAR +NO*NV
      SUX1=KAR +NO*NV
      SUX2=SUX1+MAX0(NO*NO,NO*NV)
      SUX3=SUX2+MAX0(NO*NO,NO*NV)
      XVAL=SUX3+NV*NV
      VVAL=XVAL+NT
      VPAL=VVAL+NT
      BFT=VPAL+NT
      ILAST=BFT+NV*NV
c
c >>> allocate extra core for ccsdt-1 and check if in-core is possible
c
      IF(OPTION.EQ.'SDT1')THEN
         EXTRA=NTAU+NO*NO*NIRRED+NV*NV*NIRRED/2
         WRITE(6,*)
         WRITE(6,*)' ALLOCATING ',EXTRA,' EXTRA WORDS FOR CCSDT-1A'
         WRITE(6,*)
         TAU=ILAST
         ILAST=TAU+NTAU
         if(ilast+ndimt3.lt.mxcor)then
            write(6,*)
            write(6,*)' enough memory to hold t3 in core !!'
            write(6,*)
            t3arr=ilast
            ilast=t3arr+ndimt3
            t3core='YES '
          endif
      ELSE
         TAU=1
         t3arr=1
         t3core='NO  '
      ENDIF
c
C     ALLOCATE CORE FOR FULL CCSDT
C
      if(option.eq.'FSDT')then
      write(6,*)
      write(6,*)' pointers to full ccsdt arrays'
      write(6,*)
      CHI=ILAST
      F1=CHI+NT*NT
      write(6,*)' f1 ',f1
      F2=F1+NO*NV*NV*NV
      write(6,*)' f2 ',f2
      F3=F2+NO*NO*NO*NV
      write(6,*)' f3 ',f3
      F4=F3+NO*NV*NO*NV
      write(6,*)' f4 ',f4
      F5=F4+NO*NV*NO*NV
      write(6,*)' f5 ',f5
      F6=F5+NO*NV*NO*NV
      write(6,*)' f6 ',f6
      CHI1=F6+NO*NV*NO*NV
      write(6,*)'chi1',chi1
      CHI2=CHI1+NO*NV*NV*NV
      write(6,*)'chi2',chi2
      CHIO=CHI2+NO*NO*NO*NV
      write(6,*)'chio',chio
      CHIV=CHIO+NO*NO
      write(6,*)'chiv',chiv
      FFF =CHIV+NV*NV
      write(6,*)' fff', fff
      t3=FFF+NO*NV
      write(6,*)'  t3',  t3
      t3o=t3+NDIMT3
      write(6,*)' t3o', t3o
      tau  =t3o+NDIMT3
      write(6,*)' tau',  tau
      bint=tau+ntau
      write(6,*)'bint',bint
      ilast=bint+nsgvv
      write(6,*)'ilast',ilast
      write(6,*)
      else
         CHI=1   
         F1=1
         F2=1
         F3=1
         F4=1
         F5=1
         F6=1
         CHI1=1
         CHI2=1
         CHIO=1
         CHIV=1
         FFF =1
         t3=1
         t3o=1
         tau=1
         bint=1
      endif
C
      IMCOR=MXCOR*8/1000
      REQCOR=MAX0(TOPBUF,ILAST)
      IF(REQCOR.GT.MXCOR) THEN
         ICOR =REQCOR *8/1000
         WRITE(IW,6001)IMCOR,MXCOR,ICOR,REQCOR
         STOP 'NOT ENOUGH CORE '
      END IF
C
      IF(FLDIIS.EQ.2)THEN
         DIISFL=2
      ELSE
         DIISFL=0
      ENDIF
C
      ICOR =REQCOR *8/1000
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)' >>>  CORE MEMORY REQUIREMENTS  <<<'
      WRITE(6,*)
      WRITE(IW,6001)IMCOR,MXCOR,ICOR,REQCOR
c    
      if((option.eq.'SDT1'.or.option.eq.'FSDT')
     .    .and.fl81.eq.'CHEK')then
         write(6,*)'memory check finished'
         call tstop(6)
         call qabort
      endif
C
C >>> FORM ADDRESING ARRAYS
C
      DO 140  I=1,NM
      I1=I-1
      IC(IPV+I1)=NO+1+I1/NO
      IC(IPO+I1)=I-(IC(IPV+I1)-NO-1)*NO
      IC(ITRV+I1)=I1*NV2
  140 IC(ITRT+I1)=I1*NM
      DO 141 I=1,NOV
      I1=I-1
  141 IC(ITR+I1)=I1*NO
      DO 142 I=1,NT
      I1=I-1
  142 IC(IPL+I1)=I1*NT
      IC(IPQ)=0
      DO 143 I=1,NT2-1
  143 IC(IPQ+I)=IC(IPQ+I-1)+I
      DO 144 I=1,IMX
      I1=I-1
  144 IC(ITRO+I1)=I1*NO2
      DO 145 I=1,NV
      I1=I-1
  145 IC(ITV+I1)=I1*NV
C
C >>> FORM SYMMETRY-PACKING ARRAYS
C
      CALL SYMARR(IC(ORBSYM),FLOV,NIRRED,NO,NV,NT,NOV,NM,
     .            IC(FZO),IC(FZV),NFZO,NFZV,
     .            IC(UOFF),IC(VADD),IC(ZLX),IC(ITR),IC(ITV),NONO,NVNV,
     .            IC(AOFF),IC(AADD),IC(BOFF),IC(BADD),IC(COFF),IC(CADD),
     .            IC(DOFF),IC(DADD),IC(EOFF),IC(EADD),IC(FOFF),IC(FADD),
     .            IC(T3OFF),IC(T3ADD),OPTION,
     .            IC(i1off),IC(i2off),ic(i3off))
C
C >>> CISD CALCULATION
C
      IF(OPTION.EQ.'CISD'.OR.OPTION.EQ.'BOTH')THEN
      CALL CISD9(IC(IPQ),IC(IPL),IC(ITR),IC(ITRO),IC(ITRV),
     .          IC(ITRT),IC(ITV),CC(EEE),IC(ORBSYM),FLOV,
     .          IC(UOFF),IC(VADD),IC(ZLX),CC(T2O),CC(T2),
     .          CC(ICC),CC(IBB),CC(IB2),
     .          IC(FZO),IC(FZV),
     .          CC(HH),CC(FOCK),CC(FFO),CC(FFV),CC(SFOO),CC(SFVV),
     .                  CC(AR1),CC(AR2),CC(AR3),CC(AR4),CC(AR5),CC(AR6),
     .          CC(JAR),CC(KAR),CC(GO),CC(GV),CC(T1O),CC(T1),
     .          CC(AUX1),CC(AUX2),CC(AUX3),CC(XVAL),CC(VVAL),CC(VPAL),
     .          CC(SUX1),CC(SUX2),CC(SUX3),
     .          IC(AOFF),IC(AADD),IC(BOFF),IC(BADD),IC(COFF),IC(CADD),
     .          IC(DOFF),IC(DADD),IC(EOFF),IC(EADD),IC(FOFF),IC(FADD),
     .          CC(BUF),IC(IBUF),
     .          CC(AIN),CC(BIN),CC(CIN),CC(DIN),CC(EIN),CC(FIN),CC(BFT),
     .          NIT,I79,
     .          ITYP,IC(IVAL),CC(TVAL),NCSF,NDIMT1)
C
      WRITE(IW,6074)
      WRITE(IW,6071) ESCF
      WRITE( 4,6071) ESCF
      WRITE(IW,6072) ECORR
      WRITE( 4,6072) ECORR
      ETOT=ESCF+ECORR
      WRITE(IW,6073) ETOT
      WRITE( 4,6073) ETOT
C
      ENDIF
C
C >>> CCSD CALCULATION
C
      IF(OPTION.EQ.'CCSD'.OR.OPTION.EQ.'SDT1'.OR.OPTION.EQ.'BOTH'
     .   .or.option.eq.'CCD '.or.option.eq.'FSDT')THEN
      CALL CLUS9(IC(IPQ),IC(IPL),IC(ITR),IC(ITRO),IC(ITRV),
     .          IC(ITRT),IC(ITV),CC(EEE),IC(ORBSYM),FLOV,
     .          IC(UOFF),IC(VADD),IC(ZLX),CC(T2O),CC(T2),
     .          CC(ICC),CC(IBB),CC(IB2),
     .          IC(FZO),IC(FZV),
     .          CC(HH),CC(FOCK),CC(FFO),CC(FFV),CC(SFOO),CC(SFVV),
     .          CC(AR1),CC(AR2),CC(AR3),CC(AR4),CC(AR5),CC(AR6),
     .          CC(JAR),CC(KAR),CC(GO),CC(GV),CC(T1O),CC(T1),
     .          CC(AUX1),CC(AUX2),CC(AUX3),CC(AUX4),
     .          CC(XVAL),CC(VVAL),CC(VPAL),
     .          CC(SUX1),CC(SUX2),CC(SUX3),
     .          IC(AOFF),IC(AADD),IC(BOFF),IC(BADD),IC(COFF),IC(CADD),
     .          IC(DOFF),IC(DADD),IC(EOFF),IC(EADD),IC(FOFF),IC(FADD),
     .          CC(BUF),IC(IBUF),CC(HIA),
     .          CC(AIN),CC(BIN),CC(CIN),CC(DIN),CC(EIN),CC(FIN),
     .          CC(HO),CC(HV),CC(BFT),NIT,IC(T3OFF),IC(T3ADD),CC(TAU),
     .          OPTION,NT3,CC(AR1),CC(AR2),NOACT,NVACT,MAXO,MAXV,I79,
     .          ITYP,IC(IVAL),CC(TVAL),NCSF,NDIMT1,NCL,fl81,
     .          t3core,cc(t3arr),ic(i1off),ic(i2off),ic(i3off),ndimt3,
     .          ic(trsym),
     .          cc(F1),cc(F2),cc(F3),cc(F4),cc(F5),cc(F6),
     .          cc(CHI),cc(CHI1),cc(CHI2),cc(CHIO),cc(CHIV),cc(FFF),
     .          cc(t3),cc(t3o),cc(bint))
C
      WRITE(IW,6074)
      WRITE(IW,6071) ESCF
      WRITE( 4,6071) ESCF
      WRITE(IW,6072) ECORR
      WRITE( 4,6072) ECORR
      ETOT=ESCF+ECORR
      WRITE(IW,6073) ETOT
      WRITE( 4,6073) ETOT
 6071 FORMAT('  SCF  ENERGY = ',F22.12)
 6072 FORMAT(' CORR. ENERGY = ',F22.12)
 6073 FORMAT(' TOTAL ENERGY = ',F22.12)
 6074 FORMAT(//,'CALCULATION RESULTS',/)
      write(10,6076) escf 
      write(10,6077) level,etot
 6076 FORMAT('SCF ',F22.12)
 6077 FORMAT('FSDT',i1,F22.12)
C
      ENDIF
C
      CALL RCLOSE(ITAP60,3)
      CALL RCLOSE(ITAP61,3)
      CALL RCLOSE(ITAP62,3)
      CALL RCLOSE(ITAP63,3)
      CALL RCLOSE(ITAP64,3)
      CALL RCLOSE(ITAP65,3)
      CALL RCLOSE(ITAP66,3)
C
      CALL WWRITW(ITAP79,ETOT,INTOWP(1),I79,JUNK)
      CALL RCLOSE(ITAP79,3)
c
      CALL RCLOSE(ITAP90,4)
C
      IF(NIT.GT.1.AND.DIISFL.NE.2)THEN
         CALL RCLOSE(ITAP98,4)
         CALL RCLOSE(ITAP99,4)
      ELSE
         CALL RCLOSE(ITAP98,3)
         CALL RCLOSE(ITAP99,3)
      ENDIF
C
      CALL TSTOP (6)
      call psdone()
      stop
C
 6001 FORMAT(   '  MAXIMUM  CORE = ',I7,2X,'K BYTES  OR  ',I9 ,' REAL',
     .         ' WORDS',/,
     .          '  REQUIRED CORE = ',I7,2X,'K BYTES  OR  ',I9 ,' REAL',
     .         ' WORDS',///)
      END
