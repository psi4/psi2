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
     .                 t3core,t3arr,i1off,i2off,i3off,ndimt3,trsym,
     .                 F1,F2,F3,F4,F5,F6,
     .                 CHI,CHI1,CHI2,CHIO,CHIV,FFF,
     .                 t3,t3o,bint)
      IMPLICIT INTEGER (A-Z)
      COMMON/INTS/NSGOO,NSGVV,NSGOV,NSHOV,NSLOV,NSLVO
      COMMON/PARAM/ ESCF,ENUC,NXE,
     .              EPSI,NGO,NMIN ,MINDIM,MAXDIM,MAXIT,CONVI,
     .              NO,NV,NT,NOV,NM,NV2,NT2,IMX,NONO,NVNV,NDIMT2,NTIN,
     .              NSYMHF,NIRRED,MAXOO,MAXVV,MAXOV,
     .              LEVEL,INTBUF,DIISFL,ITAP69,ITAP71,ITAP98,ITAP99,RSTR
      CHARACTER*4 OPTION,RSTR,fl81,t3core
      real*8 t3(*),t3o(*),bint(*)
      REAL*8  EREF,ESCF,ENUC,NXE,NXEO,
     .        VAL,VAL1,VAL2,DELTA,TMP,TMV,T1NORM,CONV,
     .        DELT,VALJ,VALK,VALP,VALM
      REAL*8 AR1(MAXOV,MAXOV),AR2(MAXOV,MAXOV),
     .       AR3(MAXVV,MAXOO),AR4(MAXVV,MAXOO),
     .       AR5(MAXOO,MAXVV),AR6(MAXOO,MAXVV),
     .       JAR(NV*NO),KAR(NV*NO),GO(NO,NO),GV(NV,NV),
     .       T1O(NO,NV),T1(NO,NV),
     .       AUX1(NV,NO),AUX2(NO,NO),AUX3(NV,NV),AUX4(NO,NV),
     .       HO(NO,NO),HV(NV,NV),HIA(NO,NV),
     .       SUX1(NV*NO),SUX2(NV*NO),SUX3(NV*NV),HH(NT2),
     .       FOCK(NT2),FFO(NO,NO),FFV(NV,NV),SFOO(NO*NO),SFVV(NV*NV),
     .       VVAL(NT),VPAL(NT),XVAL(NT),TAU(*),XAR1(*),XAR2(*),
     .       EIG(NT),BUF(INTBUF),t3arr(*),
     .       AIN(NSGOO),BIN(NV,NV),CIN(NSGOV),
     .       DIN(NSHOV),EIN(NSLOV),FIN(NSLVO),
     .       T2O(NDIMT2),T2(NDIMT2),BUFINT(NV*NV),
     .       CCC(MAXDIM),BBB(MAXDIM+1,MAXDIM+2),BB2(MAXDIM+1,MAXDIM+2),
     .       TVAL(*)
c these are passing arrays
      real *8 F1(*),F2(*),F3(*),F4(*),F5(*),F6(*),
     .        CHI(*),CHI1(*),CHI2(*),CHIO(*),CHIV(*),FFF(*) 
      INTEGER IBUF(INTBUF*2),NAB(8),NT3(8),
     .        AOFF(*),BOFF(*),COFF(*),DOFF(*),EOFF(*),FOFF(*),
     .        AADD(*),BADD(*),CADD(*),DADD(*),EADD(*),FADD(*),
     .        IPQ(NT2),IPL(NT),ITR(NOV),ITRO(IMX),
     .        ITRV(NM),ITRT(NM),ITV(NV),
     .        TOFF(NO,NO,2),TADD(NV,NV),ZLX(NV,NV),
     .        ORBSYM(NT),FLOV(NIRRED,4),FZO(NO),FZV(NV),NCL(*),
     .        T3OFF(NO,NO,2,NIRRED),T3ADD(NV,NV,NIRRED),ITYP(*),
     .        IVAL(*),i1off(*),i2off(no*nv,nirred),i3off(*),trsym(*)
C
      NIT=0
      CONV=10.D0**(-CONVI)
      NXEO=0.0D0
      IT=0
      IFLAG=0
      NCOUNT=0
C
C >>> CONSTRUCT FOCK MATRIX IN MO BASIS USING FILE78 AND 79 MO INTS
C
      CALL ZERO (FOCK,NT2)
C
C   TO CONSTRUCT FOCK OPERATOR ALL A, C, D, & E INTS ARE NEEDED
C   THIS IS PRESENTLY INCOMPATIBLE WITH ADDRESSING SCHEME FOR INTS
C   WHEN CORS AND VIRS ARE PRESENT.
C   TO DEBUG CORS AND VIRS, I'M FUDGING FOCK MATRIX WITH EIGENVALUES
C   THIS WILL NEED TO BE REVISITED FOR NON-SCF MOS CALCULATIONS.
C
C
C >>> READ SYMMETRY-PACKED INTEGRALS
C
      CALL RDINTS(60,AIN,NSGOO,BUF,IBUF,INTBUF,NO,NV,
     .            AOFF,AADD,ITR,ITR,FOCK)
      if(option.eq.'FSDT')then
         CALL RDINTS(61,bint,NSGVV,BUF,IBUF,INTBUF,NO,NV,
     .               BOFF,BADD,ITV,ITV,FOCK)
      endif
      CALL RDINTS(62,CIN,NSGOV,BUF,IBUF,INTBUF,NO,NV,
     .            COFF,CADD,ITV,ITR,FOCK)
      CALL RDINTS(63,DIN,NSHOV,BUF,IBUF,INTBUF,NO,NV,
     .            DOFF,DADD,ITR,ITR,FOCK)
      CALL RDINTS(64,EIN,NSLOV,BUF,IBUF,INTBUF,NO,NV,
     .            EOFF,EADD,ITR,ITR,FOCK)
      CALL RDINTS(65,FIN,NSLVO,BUF,IBUF,INTBUF,NO,NV,
     .            FOFF,FADD,ITV,ITR,FOCK)
C
      CALL FOCK2 (NT,HH,FOCK,EREF,ENUC,EIG,
     .            IPQ,ITR,ITV,NO,NV,AIN,EIN,CIN,DIN,
     .            AOFF,AADD,EOFF,EADD,COFF,CADD,DOFF,DADD,
     .            ORBSYM,FLOV,NIRRED,I79,FZO,FZV)
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
         if(option.eq.'FSDT')call sread(itap69,t3,intowp(ndimt3))
         CALL RCLOSE(ITAP69,3)
      ELSE
         CALL INITT (T1,T2,DIN,FOCK,
     .               ORBSYM,FLOV,DOFF,DADD,TOFF,TADD,IPQ,ITR,FZO,FZV,
     .               NO,NV,NIRRED)
         if(option.eq.'FSDT')call zero(t3,ndimt3)
      ENDIF
c
      if(option.eq.'CCD ')call zero(t1,no*nv)
C
C >>> MP2 (OR RESTART) CORRELATION ENERGY
C
      WRITE(6,*)' '
      CALL ECORR(T2,T1,DIN,NXE,TOFF,TADD,IPQ,ORBSYM,FLOV,
     .           NIT,NO,NV,NT,NT2,NTIN,NDIMT2,NIRRED,FZO,FZV,FOCK,
     .           DOFF,DADD,NM,ITR)
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
c
      DO 180 U=1,NO
         USYM=ORBSYM(U)
         FBE=FLOV(USYM+1,3)-NO
         LBE=FLOV(USYM+1,4)-NO
         DO 170 BE=FBE,LBE
            T1O(U,BE)=T1(U,BE)
  170    CONTINUE
  180 CONTINUE
      DO 190 I=1,NDIMT2
         T2O(I)=T2(I)
  190 CONTINUE
      if(option.eq.'FSDT')then
      do 195 i=1,ndimt3
         t3o(i)=t3(i)
  195 continue
      endif
C
c >> zero out t's 
c
      CALL ZERO(T1,NO*NV)
      CALL ZERO(T2,NDIMT2)
      if(option.eq.'FSDT')call zero(t3,ndimt3)
C
C >>> HO INTERMEDIATE
C
      CALL HOINT (T1O,T2O,DIN,FOCK,HO,XVAL,
     .            ORBSYM,FLOV,DOFF,DADD,TOFF,TADD,IPQ,ITR,FZO,FZV,
     .            ZLX,NO,NV,NIRRED)
C
C >>> HV INTERMEDIATE
C
      CALL HVINT (T1O,T2O,DIN,FOCK,HV,XVAL,
     .            ORBSYM,FLOV,DOFF,DADD,TOFF,TADD,IPQ,ITR,FZO,FZV,
     .            ZLX,NO,NV,NIRRED)
C
C >>> HIA INTERMEDIATE
C
      CALL HIAINT(T1O,T2O,DIN,FOCK,HIA,XVAL,
     .            ORBSYM,FLOV,DOFF,DADD,TOFF,TADD,IPQ,ITR,FZO,FZV,
     .            ZLX,NO,NV,NIRRED)
C
C >>> (HO, HV, HIA) * T1  --->  T1
C
      CALL HT1   (T1,T1O,T2O,CIN,DIN,FOCK,HO,HV,HIA,
     .            ORBSYM,FLOV,COFF,CADD,DOFF,DADD,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED)
C
C >>> E * T2  --->  T1
C
      CALL ET2   (T1O,T2O,EIN,AUX1,XVAL,
     .            ORBSYM,FLOV,EOFF,EADD,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED)
C
      CALL TRNPS(AUX4,AUX1,NV,NO)
C
C >>> F * T2  --->  T1
C
      CALL FT2   (T1,T1O,T2O,FIN,AUX4,XVAL,
     .            ORBSYM,FLOV,FOFF,FADD,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED)
C
C >>> GO INTERMEDIATE
C
      CALL GOINT (T1O,T2O,EIN,FOCK,HO,GO,
     .            ORBSYM,FLOV,EOFF,EADD,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED)
C
C >>> GV INTERMEDIATE
C
      CALL GVINT (T1O,T2O,FIN,FOCK,HV,GV,
     .            ORBSYM,FLOV,FOFF,FADD,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED)
C
C >>> B * TAU  --->  T2
C
      CALL BTAU  (T1O,T2O,T2,FIN,AR3,AR4,BIN,BUFINT,
     .            SUX1,SUX2,SUX3,SFVV,
     .            ORBSYM,FLOV,FOFF,FADD,NAB,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED,MAXOO,MAXVV,ITAP71)
C
C >>> D * T2  --->  T2
C     E * T1  --->  T2
C     F * T1  --->  T2
C
      CALL JKT2  (T1O,T2O,T2,CIN,DIN,EIN,FIN,AR1,AR2,
     .            SUX1,SUX2,JAR,KAR,AUX1,
     .            ORBSYM,FLOV,COFF,CADD,DOFF,DADD,
     .            EOFF,EADD,FOFF,FADD,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED,MAXOV,ITAP91)
C
C >>> A * TAU   --->  T2
C
      CALL ATAU  (T1O,T2O,T2,AIN,DIN,EIN,AR3,AR4,AR5,AR6,
     .            SUX1,SFOO,SUX3,SFVV,AUX1,
     .            ORBSYM,FLOV,AOFF,AADD,DOFF,DADD,EOFF,EADD,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED,MAXOO,MAXVV,ITAP91)
C
C >>> D  --->  T2
C
      CALL ADDD  (T2,DIN,XVAL,
     .            ORBSYM,FLOV,DOFF,DADD,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED)
C
C >>> C * T1**2  --->  T2
C     D * T1**2  --->  T2
C
      CALL CDT1  (T1O,T2,CIN,DIN,XVAL,
     .            ORBSYM,FLOV,COFF,CADD,DOFF,DADD,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED)
c
C >>> GO * T2  --->  T2
C
      CALL GOT2  (T2,T2O,GO,XVAL,AUX2,
     .            ORBSYM,FLOV,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED)
C
C >>> GV * T2  --->  T2
C
      CALL GVT2  (T2,T2O,GV,XVAL,AUX3,
     .            ORBSYM,FLOV,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED)
C
C >>> ADD CCSDT-1A TERMS
C
c     call flush(6)
      IF(OPTION.EQ.'SDT1')THEN
         if(t3core.eq.'YES ')then
         CALL coret3 (T2,T1,T2O,T1O,TAU,NO,NV,NDIMT2,NIT,
     .                ORBSYM,FLOV,NIRRED,T3OFF,T3ADD,ZLX,FZO,FZV,
     .                TOFF,TADD,OPTION,RSTR,NT3,DOFF,DADD,EOFF,EADD,
     .                FOFF,FADD,ITV,ITR,IPQ,DIN,EIN,FIN,sort,FOCK,
     .                AUX1,AUX2,AUX3,SUX1,SUX2,SUX3,
     .                XAR1,XAR2,NOACT,NVACT,MAXO,MAXV,MAXOV,
     .                AR1,AR2,SFVV,AUX3,AR1,AR2,AUX3,AUX1, 
     .                t3arr,i1off,i2off,i3off,ndimt3,trsym)
         else
         CALL CCSDT1 (T2,T1,T2O,T1O,TAU,NO,NV,NDIMT2,NIT,
     .                ORBSYM,FLOV,NIRRED,T3OFF,T3ADD,ZLX,FZO,FZV,
     .                TOFF,TADD,OPTION,RSTR,NT3,DOFF,DADD,EOFF,EADD,
     .                FOFF,FADD,ITV,ITR,IPQ,DIN,EIN,FIN,sort,FOCK,
     .                AUX1,AUX2,AUX3,SUX1,SUX2,SUX3,
     .                XAR1,XAR2,NOACT,NVACT,MAXO,MAXV,MAXOV,
     .                AR1,AR2,SFVV,AUX3,AR1,AR2,AUX3,AUX1)
         endif
      ENDIF
C
C >>> full ccsdt
C
      IF(OPTION.EQ.'FSDT')THEN
            call fccsdt(IPQ,ITR,itv,ain,bint,cin,din,ein,fin,
     .                 AOFF,AADD,BOFF,BADD,COFF,CADD,
     .                 DOFF,DADD,EOFF,EADD,FOFF,FADD,
     .                  NO,NV,NT,NT2,NDIMT2,ndimt3,nirred,level,
     .                  t1o,t1,t2o,t2,tau,t3o,t3,
     .                  ORBSYM,FLOV,tOFF,tADD,ZLX,FZO,FZV,FOCK,
     .                  F1,F2,F3,F4,F5,F6,
     .                  CHI,CHI1,CHI2,CHIO,CHIV,FFF,
     .                  i1off,i2off,i3off,trsym)
      ENDIF
C
C >>> DIVIDE BY MO ENERGIES (t3 is done in fccsdt)
C
      CALL DIVIDE(T1,T2,FOCK,ORBSYM,FLOV,
     .            TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .            NO,NV,NIRRED)
C
C >>> CORRELATION ENERGY
C
c     open(unit=4,file='term.dat',status='old')
      call ffile(4,'term',1)
c     call flush(4)
c
      if(option.eq.'CCD ')call zero(t1,no*nv)
c
      WRITE(6,*)' '
      CALL ECORR(T2,T1,DIN,NXE,TOFF,TADD,IPQ,ORBSYM,FLOV,
     .           NIT,NO,NV,NT,NT2,NTIN,NDIMT2,NIRRED,FZO,FZV,FOCK,
     .           DOFF,DADD,NM,ITR)
C
C >>> WRITE OUT T1,T2 TO FILE69
C
      CALL RFILE(ITAP69)
      CALL SWRIT(ITAP69,NO*NV,1)
      CALL SWRIT(ITAP69,NDIMT2,1)
      CALL SWRIT(ITAP69,T1,INTOWP(NO*NV))
      CALL SWRIT(ITAP69,T2,INTOWP(NDIMT2))
      if(option.eq.'FSDT')call swrit(itap69,t3,intowp(ndimt3))
      CALL SWRIT(ITAP69,ESCF,INTOWP(1))
      CALL SWRIT(ITAP69,NXE,INTOWP(1))
      CALL RCLOSE(ITAP69,3)
C
C >>> CHECK CONVERGENCE
C
      CALL TDIFF(T2,T1,T2O,T1O,T1NORM,DELT,NO,NV,NDIMT2,
     .           ORBSYM,FLOV,NIRRED,TOFF,TADD,ZLX,FZO,FZV,NXE,NXEO)
C
c     CALL TDIFF(S,t1,T2,T1O,T1NORM,DELT,NO,NV,NDIMT2,
c    .           ORBSYM,FLOV,NIRRED,UOFF,VADD,ZLX,FZO,FZV,
c    .           t3,t3o,NDIMT3,X3OFF,X3ADD)
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
         if(option.ne.'FSDT')then
         CALL TWRIT (T1,T2,ORBSYM,FLOV,NIRRED,TOFF,TADD,IVAL,TVAL,
     .               ZLX,FZO,FZV,NO,NV,NSYMHF,ITYP,NT1,NT2,NCSF,NCL)
C
         else
         CALL TWRIT3(t1,t2,ORBSYM,FLOV,NIRRED,tOFF,tADD,IVAL,TVAL,
     .               ZLX,FZO,FZV,NO,NV,NIRRED,ITYP,NT1,NT2,NCSF,
     .               t3,ndimt3,i1off,i2off,i3off,trsym)
         endif
c
c >>> if ccsdt1 gradients need to create file81 for ccqg
c
         IF(OPTION.EQ.'SDT1'.AND.fl81.EQ.'GRAD')THEN
            write(6,*)   
            write(6,*)' fl81 = grad , file81 with t3 is created'
            write(6,*)   
            CALL CCSDT1  (T2,T1,T2O,T1O,TAU,NO,NV,NDIMT2,NIT,
     .                    ORBSYM,FLOV,NIRRED,T3OFF,T3ADD,ZLX,FZO,FZV,
     .                    TOFF,TADD,OPTION,RSTR,NT3,DOFF,DADD,EOFF,EADD,
     .                    FOFF,FADD,ITV,ITR,IPQ,DIN,EIN,FIN,fl81,FOCK,
     .                    AUX1,AUX2,AUX3,SUX1,SUX2,SUX3,
     .                    XAR1,XAR2,NOACT,NVACT,MAXO,MAXV,MAXOV,
     .                    AR1,AR2,SFVV,AUX3,AR1,AR2,AUX3,AUX1)
         ENDIF
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
      if(option.ne.'FSDT')then
         CALL DIISD (T1O,T1,T2O,T2,NO,NV,NDIMT2,NIT,IT,IFLAG,
     .               EPSI,NGO,NDIIS,MINDIM,MAXDIM,
     .               CCC,BBB,ITAP98,ITAP99,ITC,BB2)
c
      else
         CALL DIISD3(T1O,t1,T2o,t2,NO,NV,NDIMT2,IT,MAXDIM,
     .               CCC,BBB,ITAP98,ITAP99,ITC,BB2,
     .               t3,t3o,NDIMT3)
      endif
c
         CALL ECORR(T2,T1,DIN,NXE,TOFF,TADD,IPQ,ORBSYM,FLOV,
     .              NIT,NO,NV,NT,NT2,NTIN,NDIMT2,NIRRED,FZO,FZV,FOCK,
     .              DOFF,DADD,NM,ITR)
      ENDIF
C
      GO TO 50
C
  601 FORMAT(13X,'t1 diag =',f12.5,5X,'DELT=',D15.8)
  602 FORMAT(///,32('*'),/,'* CORRELATION ENERGY CONVERGED *',3X,D12.2,
     . /,32('*'))
  603 FORMAT(///,32('*'),/,'* ITERATIONS EXHAUSTED *',/,32('*'))
C
      END
