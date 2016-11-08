c     program oscc
      subroutine fentry(cc,ic,mxcor)
C
C***********************************************************************
C  NOTICE OF PROGRAM MODIFICATIONS
c
c      by: Daniel Crawford
c    date: March 6, 1995
c  reason: Made wfn restart the default.
c
c      by: Daniel Crawford
c    date: january 30, 1995
c  reason: Added Tim Lee's new T1 Diagnostic as described in
c          Jayatilaka and Lee, JCP 98, 9734(1993).  See tdiff() for
c          details.
c
C      BY:  GUSTAVO E. SCUSERIA
c
c    date:  february 25, 1991 
c  reason:  add reference to output
c
c    date:  september 17, 1990
c  reason:  clean out output and add note concerning vs 3.0
c
c    date:  august 23, 1990
c  reason:  work out cct option
C
c    date:  august 21, 1990
c  reason:  add energy.dat (file10) containing geometry and energies
C
C***********************************************************************
C
      IMPLICIT INTEGER (A-Z)
      REAL*8 TIMLIM,CC,ENUC,ESCF,ECORR,ETOT
      REAL*4 EECPU,EETIME,EETCPU
      CHARACTER*32 DATTIM
      CHARACTER*4 LABEL(20),CASE,OPTION,RSTR,CHAR,ITYP(8),fl81
      COMMON/TIEMPO/EECPU,EETIME,EETCPU
      COMMON/RESTAR/TIMLIM
      COMMON/PARAM/ ESCF,ENUC,ECORR,
     .              EPSI,NGO,NMIN ,MINDIM,MAXDIM,MAXIT,CONVI,
     .              NO,NV,NT,NOV,NM,NV2,NT2,IMX,NONO,NVNV,NDIMT2,NTIN,
     .              NSYMHF,NIRRED,MAXOO,MAXVV,MAXOVa,maxovb,
     .              LEVEL,INTBUF,DIISFL,ITAP69,ITAP71,ITAP98,ITAP99,
     .              RSTR,ntmo,moint
      COMMON/INTS/NSGOO,NSGVV,NSGOV,NSHOV,NSLOV,NSLVO
      INTEGER NLAMDA(8),NCL(8),FLOV(8,2,4),NOP(8),ssov(2,4),nasym(8),
     .        NOSYM(8),NVSYM(8),NFZO(8),NFZV(8),NDOC(8),NUOC(8),NT3(8),
     .        xspin(3),yspin(3)
      DIMENSION CC(mxcor),IC(*)
      real*8 xet,xet4
      character*4 trstr
C
      call psinit('APPEND')
      call ffile(4,'term',0)
      call ffile(10,'energy',1)
C      CALL DATETM(DATTIM,EECPU,EETIME,EETCPU)
      CALL TSTART(6)
C
 6000 FORMAT(//,42('-'),/,42('-'),/,
     .' The open - shell coupled cluster program',/,
     .'       written by Gustavo Scuseria       ',/,
     .'        vs oscc 2.0   23  aug  1990     ',/,42('-'),/,42('-'),///,
     .'  Reference:',/,
     .'> ccsd(t) energy open-shell  ',/,
     .'  G. E. Scuseria',/,
     .'  Chem. Phys. Lett. 176, 27 (1991).',//)
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
      ITAP98=98
      ITAP99=99
C
      call izero(nop,8)
      xspin(1)=1
      yspin(1)=1
      xspin(2)=2
      yspin(2)=1
      xspin(3)=2
      yspin(3)=2
c        
      WRITE(IW,6000)
      INTBUF=SEC2I(100)/INTOWP(1)
      CALL RFILE(ITAP60)
      CALL RFILE(ITAP61)
      CALL RFILE(ITAP62)
      CALL RFILE(ITAP63)
      CALL RFILE(ITAP64)
      CALL RFILE(ITAP65)
      CALL RFILE(ITAP66)
      CALL RFILE(ITAP78)
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
      NTmo   = IC(18)
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
      WRITE(6,*)
      WRITE(6,*)' >>>  PARAMETERS FROM FILE30  <<<'
      WRITE(6,*)
      WRITE (6,*) ' MPOINT=',MPOINT
      WRITE (6,*) ' MCONST=',MCONST
      WRITE (6,*) ' MCALCS=',MCALCS
      WRITE (6,*) ' NCALCS=',NCALCS
      WRITE (6,*) ' NTmo  =',NTmo
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
      nt=ntmo*2
c
      WTEMP=1
      COORD=WTEMP+NT
      ZAN=COORD+3*NATOM
      dist=ZAN+NATOM
      eee=dist+natom*(natom-1)/2
      SCFTCC=WPADTI(EEE+NT)
      ORBSYM=SCFTCC+NT
      cctscf=ORBSYM+NT
      map=cctscf+nt 
      spin=map+nt
      mofz=spin+nt
      mosym=mofz+nt
      fzo=mosym+nt
C
C >>> READ FILE30
C
      CALL READ30(ITAP30,MPOINT,MCONST,MCALCS,NCALCS,NT,
     .            NSYMHF,MXCOEF,CC(EEE),NLAMDA,6,NCL,NO,
     .            CC(WTEMP),IC(SCFTCC),ITYP,IC(ORBSYM),FLOV,NIRRED,
     .            ENUC,ESCF,NATOM,NSHELL,NPRIM,
     .            CHAR,NORD,CC(ZAN),CC(COORD),ic(cctscf),
     .            IOPEN,NOP,ic(map),ic(spin),ntmo,ic(mosym),cc(dist))
C
c   note that no, nv and nt are the number of occupied, virtual and 
c   total spin orbitals respectively
c
      NV=NT-NO
C
C >>> MORE DYNAMICAL ALLOCATION
C
      FZV=FZO+NO
      TOP=FZV+NV
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

C
 6002 FORMAT(16I5)
 6004 FORMAT(20A4)
 6006 FORMAT(A4,1X,2I5)
 6007 FORMAT(//,2X,' WARNING! NO INPUT TO CCSD HAS BEEN FOUND.',/,
     .          2X,'          DEFAULT PARAMETERS WILL BE USED.')
 6024 FORMAT(2X,20A4)
 6051 FORMAT(I2,1X,I2)
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
     .             ITYP,NIRRED,NO,NV,IW,NOACT,NVACT,NSYMHF,nasym,
     .             ic(mofz))
C
      ssov(1,1)=flov(1,1,1)
      ssov(1,2)=flov(nsymhf,1,2)
      ssov(1,3)=flov(1,1,3)
      ssov(1,4)=flov(nsymhf,1,4)
      ssov(2,1)=flov(1,2,1)
      ssov(2,2)=flov(nsymhf,2,2)
      ssov(2,3)=flov(1,2,3)
      ssov(2,4)=flov(nsymhf,2,4)
c     
c     write(6,*)'ssov,alpha',ssov(1,1),ssov(1,2)
c     write(6,*)'ssov, beta',ssov(2,1),ssov(2,2)
c     write(6,*)'ssov,alpha',ssov(1,3),ssov(1,4)
c     write(6,*)'ssov, beta',ssov(2,3),ssov(2,4)
C
C >>> COUNT INTEGRALS AND VECTORS ONLY FOR ALLOCATION PURPOSES
C
      CALL NCOUNT(IC(ORBSYM),FLOV,NIRRED,NO,NV,NT,IC(FZO),IC(FZV),
     .            NDIMT1,NDIMT2,NSGOO,NSGVV,NSGOV,NSHOV,NSLOV,NSLVO,
     .            NT3,NTAU,OPTION,MAXOO,MAXVV,DIMAR,NDIMT3,
     .            MAXO,MAXV,ic(spin),ssov,ic(mofz),ic(mosym),
     .            ntmo,moint,maxova,maxovb)
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
c     WRITE(6,*)
c     WRITE(6,*)
c     WRITE(6,*)' >>>  NUMBER OF INTEGRALS   <<<'
c     WRITE(6,*)
c     WRITE(6,*)' NGOO =',NGOO,'   NSGOO =',NSGOO
c     WRITE(6,*)' NGVV =',NGVV,'   NSGVV =',NSGVV
c     WRITE(6,*)' NGOV =',NGOV,'   NSGOV =',NSGOV
c     WRITE(6,*)' NHOV =',NHOV,'   NSHOV =',NSHOV
c     WRITE(6,*)' NLOV =',NLOV,'   NSLOV =',NSLOV
c     WRITE(6,*)' NLVO =',NLVO,'   NSLVO =',NSLVO
c     WRITE(6,*)' NINT =',NINT,'   NSINT =',NSINT
c     WRITE(6,*)
c     WRITE(6,*)' REQUIRED CORE FOR INTEGRALS',NSTIN
      WRITE(6,*)' REQUIRED CORE FOR INTEGRALS',moint
C
      WRITE(IW,6042)
      WRITE(IW,6013)NDIMT1
      WRITE(IW,6014)Ndimt2
      NDIMTT=1+NDIMT1+Ndimt2
      IF(OPTION.EQ.'SDT1')THEN
         WRITE(IW,6016)NDIMT3
         NDIMTT=1+NDIMT1+NDIMT2+NDIMT3
      ENDIF
      WRITE(IW,6015)NDIMTT
 6042 FORMAT(//,1X,' >>>  NUMBER OF WALKS   <<<<<<<<<<<',/)
 6013 FORMAT(2X,'The number of single amplitudes is ',I12)
 6014 FORMAT(2X,'The number of double amplitudes is ',I12)
 6016 FORMAT(2X,'The number of triple amplitudes is ',I12)
 6015 FORMAT(2X,'The TOTAL NUMBER  of amplitudes is ',I12)
C
C >>> FINISH CORE ALLOCATION
C
      IPV=TOP
      IPO=IPV+no*nv
      ITRV=IPO+no*nv
      ITRT=ITRV+no*nv
      ITR=ITRT+no*nv
      IPL=ITR+no*nv
      IPQ=IPL+NT
      ITRO=IPQ+NT2
      ITV=ITRO+IMX
      UOFF=ITV+NV
      VADD=UOFF+no*no
      ZLX= VADD+nv*nv
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
      ioff=FADD+NM
      iadd=ioff+nt2
      IVAL=iadd+nt2
      ININ=IVAL+NCSF
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
      xin=FIN+NSLVO
      BUF=xin+moint
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
c     AR2=AR1+DIMAR
      AR2=AR1
      AR4=AR2
      AR6=AR2
      GO=AR2+DIMAR
      GV=GO+NO*NO
      HO =GV+NV*NV
      HIA= HO+NO*NO
      HV =HIA+NO*NV
      T1O= HV+NV*NV
      T1=T1O+no*nv
      FFO=T1+no*nv
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
      SUX2=SUX1+MAX0(no*no,no*nv)
      SUX3=SUX2+MAX0(no*no,no*nv)
      XVAL=SUX3+nv*nv
      VVAL=XVAL+NT
      VPAL=VVAL+NT
      BFT=VPAL+NT
      ILAST=BFT+NV*NV
      IF(OPTION.EQ.'SDT1')THEN
         EXTRA=NTAU+NO*NO*NIRRED+NV*NV*NIRRED/2
         WRITE(6,*)
         WRITE(6,*)' ALLOCATING ',EXTRA,' EXTRA WORDS FOR CCSDT-1A'
         WRITE(6,*)
         TAU=ILAST
         ILAST=TAU+NTAU
      ELSE
         TAU=1
      ENDIF
C
      IMCOR=MXCOR*8/1000
      REQCOR=MAX0(TOPBUF,ILAST)
      IF(REQCOR.GT.MXCOR) THEN
         ICOR =REQCOR *8/1000
         WRITE(IW,6001)IMCOR,MXCOR,ICOR,REQCOR
         write(6,*) 'NOT ENOUGH CORE '
         write(6,*) ' try using version 3.0 if required core is huge'
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
      WRITE(6,*)
      WRITE(6,*)' >>>  CORE MEMORY REQUIREMENTS  <<<'
      WRITE(6,*)
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
     .            IC(T3OFF),IC(T3ADD),OPTION,
     .            ntmo,ic(spin),ssov,ic(mofz),ic(mosym),
     .            ic(ioff),ic(iadd),ic(ipq),ic(map))
C
C >>> CISD CALCULATION
C
c     IF(OPTION.EQ.'CISD'.OR.OPTION.EQ.'BOTH')THEN
c     ENDIF
C
C >>> CCSD CALCULATION
C
      IF(OPTION.EQ.'CCSD'.OR.OPTION.EQ.'SDT1'.OR.OPTION.EQ.'BOTH'
     .   .or.option.eq.'CCD '.or.option.eq.'CCT '
     .   .or.option(1:3).eq.'QCI')THEN
c     IF(OPTION.EQ.'CCSD'.OR.OPTION.EQ.'SDT1'.OR.OPTION.EQ.'BOTH')THEN
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
     .          cc(xin),ic(ioff),ic(iadd),ic(mofz),ic(mosym),ic(spin),
     .          ic(map),ssov,ic(cctscf),nfzo,nfzv,xet,xet4,trstr,
     .          nop)
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
         write(10,6075) option,etot
 6075 FORMAT(a4    ,F22.12)
 6071 FORMAT('  SCF  ENERGY = ',F22.12)
 6072 FORMAT(' CORR. ENERGY = ',F22.12)
 6073 FORMAT(' TOTAL ENERGY = ',F22.12)
 6074 FORMAT(//,'CALCULATION RESULTS',/)
      endif
C
      ENDIF
C
      if(kopt.eq.0) then
         write(6,*)
         write(6,*)' unrecognizable option',option
         write(6,*)
         call qabort
      endif
C
      CALL RCLOSE(ITAP60,3)
      CALL RCLOSE(ITAP61,3)
      CALL RCLOSE(ITAP62,3)
      CALL RCLOSE(ITAP63,3)
      CALL RCLOSE(ITAP64,3)
      CALL RCLOSE(ITAP65,3)
      CALL RCLOSE(ITAP66,3)
      CALL RCLOSE(ITAP78,3)
C
c     CALL WWRITW(ITAP79,ETOT,INTOWP(1),I79,JUNK)
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
c530
      stop 'end'
C
 6001 FORMAT(   '  MAXIMUM  CORE = ',I7,2X,'K BYTES  OR  ',I9 ,' REAL',
     .         ' WORDS',/,
     .          '  REQUIRED CORE = ',I7,2X,'K BYTES  OR  ',I9 ,' REAL',
     .         ' WORDS',///)
      END
