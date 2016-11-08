      SUBROUTINE CLUS9(IPQ,IPL,ITR,ITRO,ITRV,ITRT,ITV,
     .                 EIG,ORBSYM,FLOV,TOFF,TADD,ZLX,T2O,T2,
     .                 CCC,BBB,BB2,FZO,FZV,
     .                 HH,FOCK,FFO,FFV,SFOO,SFVV,
     .                 AR1,AR2,AR3,AR4,AR5,AR6,JAR,KAR,GO,GV,
     .                 T1O,T1,AUX1,AUX2,AUX3,AUX4,XVAL,VVAL,VPAL,
     .                 SUX1,SUX2,SUX3,
     .                 AOFF,AADD,BOFF,BADD,COFF,CADD,
     .                 DOFF,DADD,EOFF,EADD,FOFF,FADD,
     .                 BUF,IBUF,HIA,
     .                 AIN,BIN,CIN,DIN,EIN,FIN,HO,HV,BUFINT,NIT,
     .                 T3OFF,T3ADD,TAU,OPTION,NT3,
     .                 XAR1,XAR2,NOACT,NVACT,MAXO,MAXV,I79,
     .                 ITYP,IVAL,TVAL,NCSF,NDIMT1,NCL,fl81,
     .                 xin,ioff,iadd,mofz,mosym,spin,map,ssov,cctscf,
     .                 nfzo,nfzv,xet,xet4,trstr,nocs,ipont,omap,f1,g1,
     .                 nalp,nvcs)
      IMPLICIT INTEGER (A-Z)
      COMMON/INTS/NSGOO,NSGVV,NSGOV,NSHOV,NSLOV,NSLVO
      COMMON/PARAM/ ESCF,ENUC,NXE,
     .              EPSI,NGO,NMIN ,MINDIM,MAXDIM,MAXIT,CONVI,
     .              NO,NV,NT,NOV,NM,NV2,NT2,IMX,NONO,NVNV,NDIMT2,NTIN,
     .              NSYMHF,NIRRED,MAXOO,MAXVV,MAXOVa,maxovb,
     .              LEVEL,INTBUF,DIISFL,ITAP69,ITAP71,ITAP98,ITAP99,
     .              RSTR,ntmo,moint
      CHARACTER*4 OPTION,RSTR,fl81,trstr
      REAL*8  EREF,ESCF,ENUC,NXE,NXEO,
     .        VAL,VAL1,VAL2,DELTA,TMP,TMV,T1NORM,CONV,
     .        DELT,VALJ,VALK,VALP,VALM,xet,xet4
      REAL*8 AR1(MAXOVa,MAXOVa),AR2(MAXOVb,MAXOVb),
     .       AR3(MAXVV,MAXOO),AR4(MAXVV,MAXOO),
     .       AR5(MAXOO,MAXVV),AR6(MAXOO,MAXVV),
     .       JAR(NV*NO),KAR(NV*NO),GO(NO,NO),GV(NV,NV),
     .       T1O(NO,NV),T1(NO,NV),
     .       AUX1(NV,NO),AUX2(NO,NO),AUX3(NV,NV),AUX4(NO,NV),
     .       HO(NO,NO),HV(NV,NV),HIA(NO,NV),
     .       SUX1(NV*NO),SUX2(NV*NO),SUX3(NV*NV),HH(NT2),
     .       FOCK(NT2),FFO(NO,NO),FFV(NV,NV),SFOO(NO*NO),SFVV(NV*NV),
     .       VVAL(NT),VPAL(NT),XVAL(NT),TAU(*),XAR1(*),XAR2(*),
     .       EIG(NT),BUF(INTBUF),
     .       AIN(NSGOO),BIN(NV,NV),CIN(NSGOV),
     .       DIN(NSHOV),EIN(NSLOV),FIN(NSLVO),
     .       T2O(NDIMT2),T2(NDIMT2),BUFINT(NV*NV),
     .       CCC(MAXDIM),BBB(MAXDIM+1,MAXDIM+2),BB2(MAXDIM+1,MAXDIM+2),
     .       TVAL(*),
     .       xin(*)
      INTEGER IBUF(INTBUF*2),NAB(8),NT3(8),
     .        AOFF(*),BOFF(*),COFF(*),DOFF(*),EOFF(*),FOFF(*),
     .        AADD(*),BADD(*),CADD(*),DADD(*),EADD(*),FADD(*),
     .        IPQ(NT2),IPL(NT),ITR(NOV),ITRO(IMX),
     .        ITRV(NM),ITRT(NM),ITV(NV),
     .        TOFF(NO,NO),TADD(NV,NV),ZLX(NV,NV),
     .        ORBSYM(NT),FLOV(8,2,4),FZO(NO),FZV(NV),NCL(*),
     .        T3OFF(NO,NO,2,NIRRED),T3ADD(NV,NV,NIRRED),ITYP(*),IVAL(*),
     .        ioff(*),iadd(*),mofz(*),mosym(*),spin(*),map(*),ssov(2,4),
     .        cctscf(*),nfzo(*),nfzv(*),ipont(nv,nv),omap(nt),f1(nt),
     .        g1(nt)
C
      NIT=0
      CONV=10.D0**(-CONVI)
      NXEO=0.0D0
      IT=0
      IFLAG=0
      NCOUNT=0
      itap91=91
      itap78=78
C
C >>> CONSTRUCT FOCK MATRIX IN MO BASIS USING FILE78 AND 79 MO INTS
C
      CALL ZERO (FOCK,NT2)
C
C >>> READ SYMMETRY-PACKED INTEGRALS
C
c     CALL RDINTS(60,AIN,NSGOO,BUF,IBUF,INTBUF,NO,NV,
c    .            AOFF,AADD,ITR,ITR,FOCK)
c     CALL RDINTS(61,BIN,NSGVV,BUF,IBUF,INTBUF,NO,NV,
c    .            BOFF,BADD,ITV,ITV,FOCK)
c     CALL RDINTS(62,CIN,NSGOV,BUF,IBUF,INTBUF,NO,NV,
c    .            COFF,CADD,ITV,ITR,FOCK)
c     CALL RDINTS(63,DIN,NSHOV,BUF,IBUF,INTBUF,NO,NV,
c    .            DOFF,DADD,ITR,ITR,FOCK)
c     CALL RDINTS(64,EIN,NSLOV,BUF,IBUF,INTBUF,NO,NV,
c    .            EOFF,EADD,ITR,ITR,FOCK)
c     CALL RDINTS(65,FIN,NSLVO,BUF,IBUF,INTBUF,NO,NV,
c    .            FOFF,FADD,ITV,ITR,FOCK)
C
      CALL RDINTS(itap78,xIN,moint,BUF,IBUF,INTBUF,NO,NV,
     .            iOFF,iADD,ipq,ipq,FOCK,mofz,nvcs,g1,nalp)
c
      call rclose(itap78,3)
c     
      CALL FOCK2 (NT,HH,FOCK,EREF,ENUC,EIG,
     .            IPQ,ITR,ITV,NO,NV,xIN,
     .            iOFF,iADD,
     .            ORBSYM,FLOV,NIRRED,I79,FZO,FZV,
     .            ntmo,spin,map,cctscf,nfzo,nfzv,omap)
      ESCF=EREF
C
C >>> SET INITIAL VALUES OF T1 AND T2
C
      IF(RSTR.EQ.'YES ')THEN
         CALL RFILE(ITAP69)
         CALL SREAD(ITAP69,JU,1)
         CALL SREAD(ITAP69,JU,1)
         CALL SREAD(ITAP69,T1,INTOWP(NO*NV))
         CALL SREAD(ITAP69,T2,INTOWP(NDIMT2))
         CALL RCLOSE(ITAP69,3)
      ELSE
         CALL INITT (T1,T2,xIN,FOCK,
     .               ORBSYM,FLOV,iOFF,iADD,TOFF,TADD,IPQ,ITR,FZO,FZV,
     .               NO,NV,NIRRED,map,spin,ssov)
      ENDIF
C
C >>> MP2 (OR RESTART) CORRELATION ENERGY
C
c     call outn(t1,t2,no,nv,map,spin,ssov,flov,orbsym,toff,tadd,mosym) 
      WRITE(6,*)' '
      write(4,'(a16,f20.12)')'reference energy',escf
      CALL ECORR(T2,T1,xIN,NXE,TOFF,TADD,IPQ,ORBSYM,FLOV,
     .           NIT,NO,NV,NT,NT2,NTIN,NDIMT2,NIRRED,FZO,FZV,FOCK,
     .           iOFF,iADD,NM,ITR,map,spin,ssov)
      write(4,'(a16,f20.12)')'   mp2    energy',escf+nxe
C
C >>> BEGIN A NEW ITERATION
C
   50 CONTINUE
C
      NIT=NIT+1
      NXEO=NXE
C
C >>> GET THE TIME
C
C      CALL TIMIT(4,6)
      close(unit=4,err=51)
      go to 52
   51 continue
      write(6,*)'error closing unit 4'
   52 continue
C
C >>> TRANSFER T1 AND T2 TO BEGIN NEW ITERATION
C
      DO 180 U=1,NO
         USYM=ORBSYM(U)
         uspn=spin(u)
         FBE=FLOV(USYM+1,uspn,3)-NO
         LBE=FLOV(USYM+1,uspn,4)-NO
         DO 170 BE=FBE,LBE
            T1O(U,BE)=T1(U,BE)
  170    CONTINUE
  180 CONTINUE
      DO 190 I=1,NDIMT2
         T2O(I)=T2(I)
  190 CONTINUE
C
      CALL ZERO(T1,NO*NV)
      CALL ZERO(T2,NDIMT2)
C
C >>> HO INTERMEDIATE
C
      CALL HOINT (T1O,T2O,xIN,FOCK,HO,XVAL,
     .            ORBSYM,FLOV,iOFF,iADD,TOFF,TADD,IPQ,ITR,FZO,FZV,
     .            ZLX,NO,NV,NIRRED,map,spin,ssov)
C
C >>> HV INTERMEDIATE
C
      CALL HVINT (T1O,T2O,xIN,FOCK,HV,XVAL,
     .            ORBSYM,FLOV,iOFF,iADD,TOFF,TADD,IPQ,ITR,FZO,FZV,
     .            ZLX,NO,NV,NIRRED,map,spin,ssov)
C
C >>> HIA INTERMEDIATE
C
      CALL HIAINT(T1O,T2O,xIN,FOCK,HIA,XVAL,
     .            ORBSYM,FLOV,iOFF,iADD,TOFF,TADD,IPQ,ITR,FZO,FZV,
     .            ZLX,NO,NV,NIRRED,map,spin,ssov)
C
C >>> (HO, HV, HIA) * T1  --->  T1
C
      CALL HT1   (T1,T1O,T2O,xIN,xIN,FOCK,HO,HV,HIA,
     .            ORBSYM,FLOV,iOFF,iADD,iOFF,iADD,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED,map,spin,ssov)
C
C >>> E * T2  --->  T1
C
      CALL ET2   (T1O,T2O,xIN,AUX1,XVAL,
     .            ORBSYM,FLOV,iOFF,iADD,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX, 
     .            NO,NV,NIRRED,map,spin,ssov)
c
      CALL TRNPS(AUX4,AUX1,NV,NO)
C
C >>> F * T2  --->  T1
C
      CALL FT2   (T1,T1O,T2O,xIN,AUX4,XVAL,
     .            ORBSYM,FLOV,iOFF,iADD,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED,map,spin,ssov)
C
C >>> GO INTERMEDIATE
C
      CALL GOINT (T1O,T2O,xIN,FOCK,HO,GO,
     .            ORBSYM,FLOV,iOFF,iADD,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED,map,spin,ssov)
c     write(6,*)'goo after goint'
c     do 293 i=1,no
c     ispn=spin(i)
c     do 292 j=1,no
c     jspn=spin(j)
c     if(ispn.ne.jspn)go to 292
c     write(6,673)i,j,go(i,j)
c 673 format(2i6,f20.12)
c 292 continue
c 293 continue
C
C >>> GV INTERMEDIATE
C
      CALL GVINT (T1O,T2O,xIN,FOCK,HV,GV,
     .            ORBSYM,FLOV,iOFF,iADD,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED,map,spin,ssov)
c
c >>> A * TAU   --->  T2
C
      CALL ATAU  (T1O,T2O,T2,xIN,xIN,xIN,AR3,AR4,AR5,AR6,
     .            SUX1,SFOO,SUX3,SFVV,AUX1,
     .            ORBSYM,FLOV,iOFF,iADD,iOFF,iADD,iOFF,iADD,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED,MAXOO,MAXVV,ITAP91,
     .            map,spin,ssov)
c     write(6,*)'after atau'
c     call outn(t1,t2,no,nv,map,spin,ssov,flov,orbsym,toff,tadd,mosym) 
c     write(6,*)'goo after atau'
c     do 393 i=1,no
c     ispn=spin(i)
c     do 392 j=1,no
c     jspn=spin(j)
c     if(ispn.ne.jspn)go to 392
c     write(6,673)i,j,go(i,j)
c 392 continue
c 393 continue
C
C >>> B * TAU  --->  T2
C
      CALL BTAU  (T1O,T2O,T2,xIN,AR3,BIN,BUFINT,
     .            SUX1,SUX3,
     .            ORBSYM,FLOV,iOFF,iADD,NAB,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,
     .            NO,NV,NIRRED,MAXOO,MAXVV,ITAP71,
     .            map,spin,ssov,ndimt2,nocs,ipont,omap)
c     write(6,*)'after btau'
c     call outn(t1,t2,no,nv,map,spin,ssov,flov,orbsym,toff,tadd,mosym) 
c
C >>> D * T2  --->  T2
C     E * T1  --->  T2
C     F * T1  --->  T2
C
      CALL JKT2a (T1O,T2O,T2,xIN,xIN,xIN,xIN,AR1,AR2,
     .            SUX1,SUX2,JAR,KAR,AUX1,
     .            ORBSYM,FLOV,iOFF,iADD,iOFF,iADD,
     .            iOFF,iADD,iOFF,iADD,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED,MAXOVa,ITAP91,
     .            map,spin,ssov)
c     write(6,*)'after jkt2a'
c     call outn(t1,t2,no,nv,map,spin,ssov,flov,orbsym,toff,tadd,mosym) 
c     write(6,*)'goo after jkt2a'
c     do 493 i=1,no
c     ispn=spin(i)
c     do 492 j=1,no
c     jspn=spin(j)
c     if(ispn.ne.jspn)go to 492
c     write(6,673)i,j,go(i,j)
c 492 continue
c 493 continue
c
c
      CALL JKT2b (T1O,T2O,T2,xIN,xIN,xIN,xIN,AR1,AR2,
     .            SUX1,SUX2,JAR,KAR,AUX1,
     .            ORBSYM,FLOV,iOFF,iADD,iOFF,iADD,
     .            iOFF,iADD,iOFF,iADD,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED,MAXOVb,ITAP91,
     .            map,spin,ssov)
c     write(6,*)'after jkt2b'
c     call outn(t1,t2,no,nv,map,spin,ssov,flov,orbsym,toff,tadd,mosym) 
c     write(6,*)'goo after jkt2b'
c     do 593 i=1,no
c     ispn=spin(i)
c     do 592 j=1,no
c     jspn=spin(j)
c     if(ispn.ne.jspn)go to 592
c     write(6,673)i,j,go(i,j)
c 592 continue
c 593 continue
c
c >>> D  --->  T2
C
      CALL ADDD  (T2,xIN,XVAL,
     .            ORBSYM,FLOV,iOFF,iADD,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,
     .            NO,NV,NIRRED,map,spin,ssov)
c     write(6,*)'after addd'
c     call outn(t1,t2,no,nv,map,spin,ssov,flov,orbsym,toff,tadd,mosym) 
C
C >>> C * T1**2  --->  T2
C     D * T1**2  --->  T2
C
      CALL CDT1  (T1O,T2,xIN,xIN,XVAL,
     .            ORBSYM,FLOV,iOFF,iADD,iOFF,iADD,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED,map,spin,ssov)
c     write(6,*)'after cdt1'
c     call outn(t1,t2,no,nv,map,spin,ssov,flov,orbsym,toff,tadd,mosym) 
c
c
C >>> GO * T2  --->  T2
C
      CALL GOT2  (T2,T2O,GO,XVAL,AUX2,
     .            ORBSYM,FLOV,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED,map,spin,ssov)
c     write(6,*)'after got2'
c     call outn(t1,t2,no,nv,map,spin,ssov,flov,orbsym,toff,tadd,mosym) 
c
C >>> GV * T2  --->  T2
C
      CALL GVT2  (T2,T2O,GV,XVAL,AUX3,
     .            ORBSYM,FLOV,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED,map,spin,ssov)
c     write(6,*)'after gvt2'
c     call outn(t1,t2,no,nv,map,spin,ssov,flov,orbsym,toff,tadd,mosym) 
c
C >>> ADD CCSDT-1A TERMS
C
c     IF(OPTION.EQ.'SDT1')THEN
c        CALL CCSDT1 (T2,T1,T2O,T1O,TAU,NO,NV,NDIMT2,NIT,
c    .                ORBSYM,FLOV,NIRRED,T3OFF,T3ADD,ZLX,FZO,FZV,
c    .                TOFF,TADD,OPTION,RSTR,NT3,DOFF,DADD,EOFF,EADD,
c    .                FOFF,FADD,ITV,ITR,IPQ,DIN,EIN,FIN,sort,FOCK,
c    .                AUX1,AUX2,AUX3,SUX1,SUX2,SUX3,
c    .                XAR1,XAR2,NOACT,NVACT,MAXO,MAXV,MAXOV,
c    .                AR1,AR2,SFVV,AUX3,AR1,AR2,AUX3,AUX1)
c     ENDIF
C
C >>> DIVIDE BY MO ENERGIES
C
      CALL DIVIDE(T1,T2,FOCK,ORBSYM,FLOV,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED,spin)
C
C >>> CORRELATION ENERGY
C
c     open(unit=4,file='term.dat',status='old')
      call ffile(4,'term',1)
      WRITE(6,*)' '
      CALL ECORR(T2,T1,xIN,NXE,TOFF,TADD,IPQ,ORBSYM,FLOV,
     .           NIT,NO,NV,NT,NT2,NTIN,NDIMT2,NIRRED,FZO,FZV,FOCK,
     .           iOFF,iADD,NM,ITR,map,spin,ssov)
c     write(6,*)'t.s used for above energy'
c     call outn(t1,t2,no,nv,map,spin,ssov,flov,orbsym,toff,tadd,mosym) 
C
C >>> WRITE OUT T1,T2 TO FILE69
C
      CALL RFILE(ITAP69)
      CALL SWRIT(ITAP69,NO*NV,1)
      CALL SWRIT(ITAP69,NDIMT2,1)
      CALL SWRIT(ITAP69,T1,INTOWP(NO*NV))
      CALL SWRIT(ITAP69,T2,INTOWP(NDIMT2))
      CALL SWRIT(ITAP69,ESCF,INTOWP(1))
      CALL SWRIT(ITAP69,NXE,INTOWP(1))
      CALL RCLOSE(ITAP69,3)
C
C >>> CHECK CONVERGENCE
C
      CALL TDIFF(T2,T1,T2O,T1O,T1NORM,DELT,NO,NV,NDIMT2,spin,ssov,
     .           ORBSYM,FLOV,NIRRED,TOFF,TADD,FZO,FZV,NXE,NXEO)
C
      WRITE(6,601)T1NORM,DELT
      WRITE(4,601)T1NORM,DELT
C
c >>> stop calculation if convergence has been achieved
c
      IF(DELT.LT.CONV) THEN
         WRITE(6,602) CONV
         NT1=NDIMT1
         NT2=NSHOV
c        CALL TWRIT (T1,T2,ORBSYM,FLOV,NIRRED,TOFF,TADD,IVAL,TVAL,
c    .               ZLX,FZO,FZV,NO,NV,NSYMHF,ITYP,NT1,NT2,NCSF,NCL)
c
c >>> if ccsdt1 gradients need to create file81 for ccqg
c
c        IF(OPTION.EQ.'SDT1'.AND.fl81.EQ.'GRAD')THEN
c           write(6,*)   
c           write(6,*)' fl81 = grad , file81 with t3 is created'
c           write(6,*)   
c           CALL CCSDT1  (T2,T1,T2O,T1O,TAU,NO,NV,NDIMT2,NIT,
c    .                    ORBSYM,FLOV,NIRRED,T3OFF,T3ADD,ZLX,FZO,FZV,
c    .                    TOFF,TADD,OPTION,RSTR,NT3,DOFF,DADD,EOFF,EADD,
c    .                    FOFF,FADD,ITV,ITR,IPQ,DIN,EIN,FIN,fl81,FOCK,
c    .                    AUX1,AUX2,AUX3,SUX1,SUX2,SUX3,
c    .                    XAR1,XAR2,NOACT,NVACT,MAXO,MAXV,MAXOV,
c    .                    AR1,AR2,SFVV,AUX3,AR1,AR2,AUX3,AUX1)
c        ENDIF
c
c >>> ccsd(t) option
c
      if(option.eq.'CCT ')then
      call sdt(t1,t2,xin,xin,xin,fock,
     .         orbsym,flov,fzo,fzv,toff,tadd,zlx,ipq,
     .         itr,itv,ioff,iadd,ioff,iadd,ioff,iadd,
     .         no,nv,nirred,xet,xet4,trstr,
     .         map,spin,ssov,omap)
      endif
         RETURN
      ENDIF
C
      IF(NIT.GE.MAXIT)THEN
         WRITE(6,603)
         RETURN
      ENDIF
C
C >>> DIIS EXTRAPOLATION
C
      IF(DIISFL.NE.2) THEN
         CALL DIISD (T1O,T1,T2O,T2,NO,NV,NDIMT2,NIT,IT,IFLAG,
     .               EPSI,NGO,NDIIS,MINDIM,MAXDIM,
     .               CCC,BBB,ITAP98,ITAP99,ITC,BB2)
         CALL ECORR(T2,T1,xIN,NXE,TOFF,TADD,IPQ,ORBSYM,FLOV,
     .              NIT,NO,NV,NT,NT2,NTIN,NDIMT2,NIRRED,FZO,FZV,FOCK,
     .              iOFF,iADD,NM,ITR,map,spin,ssov)
      ENDIF
C
      GO TO 50
C
  601 FORMAT(13X,'t1 diag =',f15.8,5X,'DELT=',D15.8)
  602 FORMAT(///,32('*'),/,'* CORRELATION ENERGY CONVERGED *',3X,D12.2,
     . /,32('*'))
  603 FORMAT(///,32('*'),/,'* ITERATIONS EXHAUSTED *',/,32('*'))
C
      END
c
      subroutine outn(t1,t2,no,nv,map,spin,ssov,flov,orbsym,toff,tadd,
     .                mosym) 
      implicit integer(a-z)
      real*8 t1(no,nv),t2(*),val1,val2,val3
      integer map(*),spin(*),ssov(2,4),flov(8,2,4),orbsym(*),
     .        toff(no,no),tadd(nv,nv),mosym(*)
      real*8 x2(3,5,5,9,9)
c
      do 400 i=1,no
         xi=map(i)
         ispn=spin(i)
         isym=orbsym(i)
            fa=flov(isym+1,ispn,3)-no
            la=flov(isym+1,ispn,4)-no
            do 420 a=fa,la
               xa=map(a+no)
               aspn=spin(a+no)
                  write(6,677)xi,ispn,xa,aspn,t1(i,a) 
  420       continue
  400 continue
  677 format(2(i4,i3),f20.12)
c
      do 300 i=1,no
         xi=map(i)
         ispn=spin(i)
         isym=orbsym(i)
         do 310 j=1,i-1
            xj=map(j)
            jspn=spin(j)
            jsym=orbsym(j)
            ijsym=ieor(isym,jsym)
            fa=ssov(ispn,3)-no
            la=ssov(ispn,4)-no
            do 320 a=fa,la
               xa=map(a+no)
               aspn=spin(a+no)
               asym=orbsym(a+no)
               bsym=ieor(asym,ijsym)
               fb=flov(bsym+1,jspn,3)-no
               lb=flov(bsym+1,jspn,4)-no
               if(lb.ge.a)lb=a-1
               do 330 b=fb,lb
                  xb=map(b+no)
                  bspn=spin(b+no)
                  ijab=toff(i,j)+tadd(a,b)
                  write(6,678)xi,ispn,xj,jspn,xa,aspn,xb,bspn,t2(ijab)
  330          continue
  320       continue
  310    continue
  300 continue
c
  678 format(4(i4,i3),f20.12)
c
      return
      end
c
      subroutine outt(t1,t2,no,nv,map,spin,ssov,flov,orbsym,toff,tadd,
     .                mosym) 
      implicit integer(a-z)
      real*8 t1(no,nv),t2(*),val1,val2,val3
      integer map(*),spin(*),ssov(2,4),flov(8,2,4),orbsym(*),
     .        toff(no,no),tadd(nv,nv),mosym(*)
      real*8 x2(3,5,5,9,9)
c
      no2=no/2
      nv2=nv/2
      do 300 i=1,no
         xi=map(i)
         ispn=spin(i)
         isym=orbsym(i)
         do 310 j=1,i-1
            xj=map(j)
            jspn=spin(j)
            jsym=orbsym(j)
            ijsym=ieor(isym,jsym)
            fa=ssov(ispn,3)-no
            la=ssov(ispn,4)-no
            do 320 a=fa,la
               xa=map(a+no)
               asym=orbsym(a+no)
               bsym=ieor(asym,ijsym)
               fb=flov(bsym+1,jspn,3)-no
               lb=flov(bsym+1,jspn,4)-no
               if(lb.ge.a)lb=a-1
               do 330 b=fb,lb
                  xb=map(b+no)
                  ijab=toff(i,j)+tadd(a,b)
                  if(ispn.eq.1.and.ispn.eq.jspn)then
                  x2(1,xi,xj,xa-no2,xb-no2)=t2(ijab)
                  x2(1,xj,xi,xb-no2,xa-no2)=t2(ijab)
                  else if(ispn.eq.2.and.ispn.eq.jspn)then
                  x2(3,xi,xj,xa-no2,xb-no2)=t2(ijab)
                  x2(3,xj,xi,xb-no2,xa-no2)=t2(ijab)
                  else
                  x2(2,xi,xj,xa-no2,xb-no2)=t2(ijab)
                  endif
c                 write(6,678)ispn,jspn,xi,xj,xa,xb,t2(ijab)
  330          continue
  320       continue
  310    continue
  300 continue
c
  678 format(2i2,3x,4i3,f20.12)
c
      do 400 i=1,no2
         isym=mosym(i)
         do 410 a=1,nv2 
            asym=mosym(a+no2)
            iasym=ieor(isym,asym)
            do 420 j=1,i  
               jsym=mosym(j)
               lb=nv2
               if(i.eq.j)lb=a
               do 430 b=1,lb 
                  bsym=mosym(b+no2)
                  jbsym=ieor(jsym,bsym)
                  if(iasym.ne.jbsym)go to 430
                  val1=x2(1,i,j,a,b)
c                 val2=x2(2,i,j,a,b)-x2(2,i,j,b,a)
                  val2=x2(2,i,j,a,b)
                  val3=x2(3,i,j,a,b)
            write(6,674)i,j,a+no2,b+no2,val1,val2,val3
  674   format(4i3,1x,3f20.12)
  430          continue
  420       continue
  410    continue
  400 continue
c
      return
      end
