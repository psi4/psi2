cges  PROGRAM CCEG
      subroutine fentry(cc,ic,mxcor)
C
C***********************************************************************
C  NOTICE OF PROGRAM MODIFICATIONS
C      BY:  GUSTAVO E. SCUSERIA (unless otherwise noted)
c
c    date:  Mar. 6, 1995
c  reason:  Made wfn restart the default
c     who:  Daniel Crawford
c   
c    date:  Thu Feb 28 08:18:44 EST 1991
c  reason:  Conversion to new style input format
c     who:  Curtis L Janssen
c   where:  U of Georgia
c
c    date:  june 12, 1990   
c  reason:  implement non-scf mos for cisd                    
c           debugged for h2dzp(1doc/9uoc) and hfsto(5occ/1uoc)
c
c    date:  june 7, 1990   
c  reason:  replace w-matrix  formation with faster procedure
c            both for energy and gradient CCT
c
c    date:  june 5, 1990   
c  reason:  implement restart option for CCSD(T) energy
c
c    date:  may 29, 1990   
c  reason:  do not rfile itap78
c
c    date:  may 25, 1990   
c  reason:  re-check cisd option & change t1 norm printout to diag 
c
c    date:  may 18, 1990   
c  reason:  new version of ccsd(t) gradients for i>j>k; a>b>c       
c
c    date:  april 27, 1990
c  reason:  include in this version modifications done to energy.dat
c
c    date:  april 27, 1990
c  reason:  when printing # walks including triples
c           change ndimt2 to nshov
c
c    date:  april 17, 1990
c  reason:  include ccsd(t) gradients!
c
c    date:  march 28, 1990
c  reason:  add subroutine sdt for ccsd(t), option=CCT 
c           and write out unrecognizable option
c
c    date:  march 16, 1990
c  reason:  add file10 containing geometry and energies
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
      IMPLICIT INTEGER (A-Z)
      REAL*8 TIMLIM,CC,ENUC,ESCF,ECORR,ETOT,xet,xet4
      REAL*4 EECPU,EETIME,EETCPU
      CHARACTER*32 DATTIM
      CHARACTER*4 CASE,OPTION,RSTR,CHAR,ITYP(8),fl81,trstr
      character*80 label
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
c
#include <error.h>
      integer frdi, frdc, frdboo, errcod, bRSTR, frdia1, cout, fcnt
      integer prdup, prcntr
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
      errcod = prcntr('ON BRIEF')
      errcod = prdup()
      if (prcntr('IS_ON DEBUG').ne.0) then
        errcod = prcntr('OFF BRIEF')
        endif
clj Initialize the print routines.
C
clj Save paper in the output.  The references have been moved to the
clj manual page.
 6000 FORMAT(/,
     .' THE CLOSED-SHELL COUPLED CLUSTER PROGRAM',/,
     .'      REFORMULATED VERSION JUNE 1988     ',/,
     .'       WRITTEN BY GUSTAVO SCUSERIA       ',/,
     .'        VS CCEG 2.0  jun  05, 1990     ')
c6000 FORMAT(//,42('-'),/,42('-'),/,
c    .' THE CLOSED-SHELL COUPLED CLUSTER PROGRAM',/,
c    .'      REFORMULATED VERSION JUNE 1988     ',/,
c    .'       WRITTEN BY GUSTAVO SCUSERIA       ',/,
c    .'        VS CCEG 2.0  jun  05, 1990     ',/,42('-'),/,42('-'),///,
c    .' References :',/,
c    .'',/,
c    .'> ccsd energy',/,
c    .'  G.E.Scuseria, A.C.Scheiner, T.J.Lee, J.E.Rice and H.F.Schaefer,
c    .',/,
c    .'  J. Chem. Phys. 86, 2881 (1986).',/,
c    .'  G.E.Scuseria, C.L.Janssen and H.F.Schaefer,',/,
c    .'  J. Chem. Phys. 89, 7382 (1988).',/,
c    .'',/,
c    .'> ccsd energy gradient',/,
c    .'  A.C.Scheiner, G.E.Scuseria, J.E.Rice, T.J.Lee and H.F.Schaefer,
c    .',/,
c    .'  J. Chem. Phys. 87, 5361 (1987).',/,
c    .'',/,
c    .'> ccsdt-1 energy and gradient',/,
c    .'  G.E.Scuseria and H.F.Schaefer,',/,
c    .'  Chem. Phys. Lett. 146, 23 (1988).',/,
c    .'',/,
c    .'> ccsd(t) closed-shell energy ',/,
c    .'  G.E. Scuseria and T.J. Lee',/,
c    .'  J. Chem. Phys. 93, 5851 (1990).',/,
c    .'',/,
c    .'> ccsd(t) closed-shell gradient ',/,
c    .'  G.E. Scuseria',/,
c    .'  J. Chem. Phys. 94, 442 (1991).',///)
C
      kopt=0
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
      ITAP79=79
      ITAP98=98
      ITAP99=99
      WRITE(IW,6000)
      INTBUF=SEC2I(100)/INTOWP(1)
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
      iphase = ic(69)
C
      if (prcntr('IS_ON BRIEF').eq.0) then
        WRITE(6,*)
        WRITE(6,*)' >>>  PARAMETERS FROM FILE30  <<<'
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
        endif
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
      dist=zan+natom
      EEE=dist+natom*(natom-1)/2
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
     .            IOPEN,NOP,cc(dist))
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
C
      call rdinpt(epsi,ngo,nmin,mindim,maxdim,maxit,convi,
     &            ncsf,fldiis,option,fl81,label,rstr,trstr,nfzo,nfzv,
     &            nirred,iphase)

c Make sure the the option returned by rdinpt is acceptable
      if (     (option.ne.'CCSD')
     &    .and.(option.ne.'CCT')
     &    .and.(option.ne.'CCD')
     &    .and.(option.ne.'SDT1')
     &    .and.(option.ne.'CISD')
     &    .and.(option.ne.'BOTH')) then
        write(6,*) 'ERROR: Illegal ''WFN'''
        call qabort
        endif
c
c >>> DONE READING INPUT
c
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
     .        NT3,NTAU,OPTION,MAXOO,MAXVV,MAXOV,DIMAR,NDIMT3,MAXO,MAXV)
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
      if (prcntr('IS_ON BRIEF').eq.0) then
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
        endif
      WRITE(6,*)
      WRITE(6,*)' REQUIRED MEMORY FOR INTEGRALS',NSTIN
C
      WRITE(IW,6042)
      WRITE(IW,6013)NDIMT1
      WRITE(IW,6014)NSHOV
      NDIMTT=1+NDIMT1+NSHOV
      IF(OPTION.EQ.'SDT1')THEN
         WRITE(IW,6016)NDIMT3
         NDIMTT=1+NDIMT1+nshov+NDIMT3
      ENDIF
      WRITE(IW,6015)NDIMTT
 6042 FORMAT(/,1X,' >>>  NUMBER OF WALKS   <<<<<<<<<<<')
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
      tra=IVAL+NCSF
      tva=tra+no*no
      inin=tva+nv*nv
      IF(OPTION.EQ.'SDT1')THEN
         T3OFF=ININ
         T3ADD=T3OFF+NO*NO*2*NIRRED
         ININ=T3ADD+NV*NV*NIRRED
      ELSE
         T3OFF=1
         T3ADD=1
      ENDIF
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
c For the (T) correction each thread of execution needs a copy of
c w and v.  Nthrd is used to allocate storage for these arrays.
#ifdef SGI_PARALLEL
      nthrd = mp_numthreads()
#else
      nthrd = 1
#endif
      IF(OPTION.EQ.'SDT1')THEN
         EXTRA=NTAU+NO*NO*NIRRED+NV*NV*NIRRED/2
         WRITE(6,*)
         WRITE(6,*)' ALLOCATING ',EXTRA,' EXTRA WORDS FOR CCSDT-1A'
         WRITE(6,*)
         TAU=ILAST
         ILAST=TAU+NTAU
      ELSE if(option.eq.'CCT '.and.fl81.eq.'GRAD')then
         ww=ilast
         vv=ww+no*no*no*nthrd
         w1=vv+no*no*no*nthrd
         v1=w1+nv*nv*nv*nthrd
         j1=v1+nv*nv*nv*nthrd
         j2=j1+no*nv*nthrd
         moo=j2+ndimt2*nthrd
         nvv=moo+no*no*nthrd
         dden=nvv+nv*nv*nthrd
         eden=dden+nshov*nthrd
         fden=eden+nslov*nthrd
         ffden=fden+no*nv*nv
         ilast=ffden+nslvo*nthrd
         extra=ilast-ww  
         write(6,*)
         write(6,*)'allocating',extra,' extra words for ccsd(t) grad'
         write(6,*)
      ELSE if(option.eq.'CCT '.and.fl81.ne.'GRAD')then
         ww=ilast
         vv=ww+no*no*no*nthrd
         ilast=vv+no*no*no*nthrd
         extra=ilast-ww  
         write(6,*)
         write(6,*)'allocating',extra,' extra words for ccsd(t) energy'
         write(6,*)
      else
         TAU=1
         vv=1
         ww=1
         j1=1          
         j2=1
         moo=1
         nvv=1
         dden=1
         eden=1
         fden=1
      ENDIF
C
      IMCOR=MXCOR*8/1000
      REQCOR=MAX0(TOPBUF,ILAST)
      IF(REQCOR.GT.MXCOR) THEN
         ICOR =REQCOR *8/1000
         WRITE(IW,6001)IMCOR,MXCOR,ICOR,REQCOR
         write(6,*) 'ERROR: NOT ENOUGH MEMORY HAS BEEN ALLOCATED'
         call qabort
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
      WRITE(6,*)' >>>  MEMORY REQUIREMENTS  <<<'
      WRITE(IW,6001)IMCOR,MXCOR,ICOR,REQCOR
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
     .            IC(T3OFF),IC(T3ADD),OPTION,ic(tra),ic(tva))
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
      kopt=1
      WRITE(IW,6074)
      WRITE(IW,6071) ESCF
      WRITE( 4,6071) ESCF
      WRITE(IW,6072) ECORR
      WRITE( 4,6072) ECORR
      ETOT=ESCF+ECORR
      WRITE(IW,6073) ETOT
      WRITE( 4,6073) ETOT
      write(10,6076) escf 
      write(10,6079) option,etot
      go to 1000
      ENDIF
C
C >>> CC CALCULATION
C
      IF(OPTION.EQ.'CCSD'.OR.OPTION.EQ.'SDT1'.OR.OPTION.EQ.'BOTH'
     .   .or.option.eq.'CCD '.or.option.eq.'CCT '
     .   .or.option(1:3).eq.'QCI')THEN
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
     .          ITYP,IC(IVAL),CC(TVAL),NCSF,NDIMT1,NCL,fl81,xet,xet4,
     .          cc(ww),cc(vv),ic(tra),ic(tva),
     .          cc(j1),cc(j2),cc(moo),cc(nvv),cc(dden),cc(eden),
     .          cc(fden),cc(ffden),cc(w1),cc(v1),trstr)
C
      kopt=1
      if(option.eq.'CCT ')then
      WRITE(IW,6074)
      WRITE(IW,6171) ESCF
      WRITE( 4,6171) ESCF
      write(10,6076) escf 
      WRITE(IW,6172) ECORR
      WRITE( 4,6172) ECORR
      ETOT=ESCF+ECORR
      WRITE(IW,6173) ETOT
      WRITE( 4,6173) ETOT
      write(10,6077) etot
      write( 6,6175) xet
      write( 4,6175) xet
      etot=etot+xet
      WRITE(IW,6177) ETOT
      WRITE( 4,6177) ETOT
      write(10,6078) etot
 6076 FORMAT('SCF ',F22.12)
 6077 FORMAT('CCSD',F22.12)
 6078 FORMAT('CCT ',F22.12)
 6171 FORMAT('  scf  energy = ',F22.12)
 6172 FORMAT(' ccsd  corr.  = ',F22.12)
 6173 FORMAT(' ccsd  energy = ',F22.12)
 6175 FORMAT(' (t)   corr.  = ',F22.12)
 6177 FORMAT('ccsd(t)energy = ',F22.12)
      else 
      WRITE(IW,6074)
      WRITE(IW,6071) ESCF
      WRITE( 4,6071) ESCF
      WRITE(IW,6072) ECORR
      WRITE( 4,6072) ECORR
      ETOT=ESCF+ECORR
      WRITE(IW,6073) ETOT
      WRITE( 4,6073) ETOT
      write(10,6076) escf 
      write(10,6079) option,etot
 6079 FORMAT(A4,    F22.12)
 6071 FORMAT('  SCF  ENERGY = ',F22.12)
 6072 FORMAT(' CORR. ENERGY = ',F22.12)
 6073 FORMAT(' TOTAL ENERGY = ',F22.12)
 6074 FORMAT(/,'CALCULATION RESULTS')
      endif
      go to 1000
      ENDIF
C
      if(kopt.eq.0) then
         write(6,*)
         write(6,*)' cceg: unrecognizable option',option
         write(6,*)
         call qabort
      endif
c
 1000 continue
c
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
#if defined(AIX)
      stop 'end'
#else
      stop
#endif
C
 6001 FORMAT(   '  MAXIMUM  CORE = ',I7,2X,'K BYTES  OR  ',I9 ,' REAL',
     .         ' WORDS',/,
     .          '  REQUIRED CORE = ',I7,2X,'K BYTES  OR  ',I9 ,' REAL',
     .         ' WORDS')
      END
