      SUBROUTINE GTINTS(A,ETA,ZETA,F,ALPHA,BETA,IJORD,KADD,LADD,
     1 IJADD,IJGRP,ORBSYM,INTS,S,EPA,JT1C,IOFF,EX,RS,LG,IJDEP,BOX,
     2 EAVF,LSTEPA,aa,doff,dadd,idim,idima,e11,e12,e22)
      IMPLICIT INTEGER (A-Z)
      LOGICAL IJEQ,IJKEQ,JKEQ
      INTEGER JULIA(70),IJDEP(NUMIJ)
      REAL*8 A0,A1,A2,A1S2,A4,FOCC,BOX(NDF),A3S2
      REAL*8 EAVF(NUMIJ),LSTEPA(NUMIJ)
      real*8 aa(idima)
      integer doff(idim),dadd(idim)
      DIMENSION IJORD(NUMIJ),KADD(SYMORB),LADD(SYMORB),JT1C(1)
      DIMENSION ORBSYM(NORBS),IJADD(NUMIJ),IJGRP(NUMIJ)
      DIMENSION IOFF(NUMIJ)
      real*8 ci1,ci2,h11,h12,h22,e11(nbf,nbf),e12(nbf,nbf),e22(nbf,nbf)
      real*8 derc1,derc2,one
      REAL*8 A(IJIND,NDF),F(NORBS),INTS(NMAX),LG(NORBS,NORBS)
     1,S(NUMIJ,NDF),EPA(NORBS,NORBS),EX(IJIND),RS(NUMDEP,NDF)
cets121890      REAL*8 ZETA(NUMIJ,NORBS),ALPHA(NUMIJ),BETA(NUMIJ)
      REAL*8 zeta(ntypes+1,numij),ALPHA(NUMIJ),BETA(NUMIJ)
      REAL*8 C1A,C2A,C1,C2
      REAL*8 VAL1,VAL2,VAL3,ETA(NUMIJ),ESCF,REPNUC
      COMMON /CIA/ C1A(45),C2A(45),C1,C2,ESCF
C!!!  COMMON /TAPES/     ,   ,NOUT,IT2,    ,JTAP1,     ,  OLD TAPES USED
C!!!  COMMON /TAPES/ IT22,IT8,NOUT,IT2,IT42,JTAP1,IUTAP,JTAP2
      COMMON /TAPES/ NOUT,ICHK,IT30,IT42,IT47,IT52,IT58,IT94,IT96,IT98
      COMMON /DIMSG/ SYMORB,NUMIJ,NBF,NMAX,NBFORB,NSYM,NDF,NORBS
      COMMON /POINT/ FOCC(5),IPOINT(5),JPOINT(5),NSORB(5)
      COMMON /FUNCS/ IOCC,JOCC,KOCC,NTYPES,NBSET,NAT
      COMMON /OCCS/ NSPE,NSA,NSB,NALP,NOP,NDOC,LFERM,NUOC,IJIND,IJIND3
      COMMON/OCCS2/NUMCOR,NCOR,FDOC,CORDOC,NUMVIR,FVIR,LUOC,NUMDEP,VUOC
      COMMON /LOC22/ IA22(100),IBF2AT(284)
      COMMON /IEXT/ IPRINT,IGRP,ICIDIP
      COMMON /REN/ REPNUC
      common /tcscf/ derc1,derc2,ci1,ci2,h11,h12,h22,one,itcscf
C
 3001 FORMAT(10I5)
 1100 FORMAT(/'  IN THE CONSTRUCTION OF THE A MATRIX, THE CASE ',A7
     * /'  IS NOT POSSIBLE IN THE ',A11,'  SET OF LOOPS.'/'  ERROR
     * IN LOOP STRUCTURE SOMEWHERE.')
   2  FORMAT(/'  IN GTINT2, IJKLT UNKNOWN TYPE ',I5/'  ORBITALS: ',4I5
     * /'  BASIS FUNCTIONS: ',4I5/,'  SYMMETRIES: ',4I5)
 2005 FORMAT(/'  IN GTINT2, LAD TOO BIG: ',I6)
      IJIND2=IJIND
      A0=0.0D+00
      A1=1.0D+00
      A2=2.0D+00
      A4=4.00D+00
      A1S2=0.50D+00
      A3S2=1.50D+00
      h11=a0
      h12=a0
      h22=a0
      one=a0
      derc1=a0
      derc2=a0
      call zero(aa,idima)
      if(itcscf.eq.1) then
       call zero(e11,nbf*nbf)
       call zero(e12,nbf*nbf)
       call zero(e22,nbf*nbf)
      endif
C
C THIS IS THE BEGINNING OF THE SECOND PART OF GETINT TJL 11/24/84
C
      IF (NSPE.EQ.0) RETURN
      DO 3400 IOR=NSB,NSA,-1
        IO=IOFF(IOR)
        II=IO+IOR
        ISM=ORBSYM(IOR)
        ITYP=1
        DO 3410 JOR=IOR,NSA,-1
          JO=IOFF(JOR)
          IJ=IO+JOR
          JJ=JO+JOR
          IJND=IJORD(IJ)
          JSM=ORBSYM(JOR)
          IJEQ=IOR.EQ.JOR
          IJT=ITYP
          IF (IJEQ) IJT=IJT+2
          IJSYM=IEOR((ISM-1),(JSM-1))
          KOFF=IJSYM*NORBS
          JAD=IJADD(IJ)
          IF (IGRP .NE. IJGRP(IJ)) THEN
            CALL SREAD(IT52,INTS,NMAX*2)
            SRIN=SRIN+NMAX*2/112
            IGRP=IJGRP(IJ)
          END IF
          DO 3420 KOR=1,NUOC
            KO=IOFF(KOR)
            IK=IO+KOR
            JK=JO+KOR
            KK=KO+KOR
            IKND=IJORD(IK)
            JKND=IJORD(JK)
            KSM=ORBSYM(KOR)
            iksym=ieor((ism-1),(ksm-1))
            jksym=ieor((jsm-1),(ksm-1))
            JKEQ=JOR.EQ.KOR
            IJKT=IJT
            IF (JKEQ) IJKT=IJKT+1
            IJKEQ=IJEQ.AND.JKEQ
            IJKSYM=IEOR(IJSYM,(KSM-1))
            LOFF=IJKSYM*NORBS
            KAD=JAD+KADD(KOFF+KOR)
CTJL
C      I1FL = 0
C      IF(I1FL.EQ.1) GOTO 3431
            DO 3430 LOR=1,KOR
         IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO 303
              LSM=ORBSYM(LOR)
              ilsym=ieor((ism-1),(lsm-1))
              jlsym=ieor((jsm-1),(lsm-1))
              klsym=ieor((ksm-1),(lsm-1))
              IF ((LSM-1) .NE. IJKSYM) GOTO 303
              LO=IOFF(LOR)
              IL=IO+LOR
              JL=JO+LOR
              KL=KO+LOR
              LL=LO+LOR
              ILND=IJORD(IL)
              JLND=IJORD(JL)
              KLND=IJORD(KL)
              IJKLT=IJKT
              IF (IJKEQ .AND. KOR.GT.LOR) IJKLT=IJKLT+1
              IF (KOR.EQ.LOR) IJKLT=IJKLT+3
              LAD=KAD+LADD(LOFF+LOR)
              IF (IJKLT.GT.7 .OR. IJKLT.LT.1) GOTO 9903
              IF (LAD .GT. NMAX) GOTO 9902
CTJL
         KLVU = IJDEP(KL)
CTJL
C
C     SPE-SPE-UOC-UOC
C
              GOTO (13,23,33,43,53,63,73)IJKLT
 13           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
cets121890
      if(itcscf.eq.1) then
        IF(IOR.EQ.NSB) CALL DCICOF(DERC1,KK,LL,LOR,KOR,VAL2,VAL1,
     +     rs,ex,eavf,lg,s)
        IF(IOR.EQ.NSA) CALL DCICOF(DERC2,KK,LL,LOR,KOR,VAL2,VAL1,
     +     rs,ex,eavf,lg,s)
      endif
cets
C
      IF(IKSYM.EQ.0.AND.JLSYM.EQ.0)
     .   CALL GA1(IK  ,JL  ,IL,KL,IJ,JK,VAL1,VAL2,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
      IF(ILSYM.EQ.0.AND.JKSYM.EQ.0)
     .   CALL GA1(IL  ,JK  ,IK,KL,IJ,JL,VAL3,VAL2,VAL1,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB1(IKND,IL,KL,IJ,JK,VAL1,VAL2,VAL3,JL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JLND,JK,KL,IJ,IL,VAL1,VAL2,VAL3,IK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(ILND,IJ,JL,IK,KL,VAL3,VAL2,VAL1,JK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JKND,IJ,IK,JL,KL,VAL3,VAL2,VAL1,IL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
      IF(KLVU.NE.0) CALL EX1(IJND,LL,KK,LOR,KOR,VAL2,VAL1,VAL3,
     1                       rs,ex,eavf,lg,s)
      IF(KLVU.NE.0) CALL RS1(KLVU,A1,VAL2,VAL1,VAL3,IJ,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 303
 23           CONTINUE
              WRITE(NOUT,2012)
 2012         FORMAT(/     '  IN CONSTRUCTON OF A, J=K IS NOT POSSIBLE I
     *N THE FIRST SET OF LOOPS'/'  ERROR IN LOOPS SOMEWHERE')
              CALL MABORT
 33           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      IF(IKSYM.EQ.0.AND.ILSYM.EQ.0)
     .   CALL GA3(IK  ,IL  ,IL,KL,II,IK,VAL1,VAL2,IOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB3(IKND,IL,KL,II,IK,VAL1,VAL2,IOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB3(ILND,IK,KL,II,IL,VAL1,VAL2,IOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
      IF(KLVU.NE.0) CALL RS1(KLVU,A1S2,VAL2,VAL1,VAL1,IJ,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 303
 43           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      IF(IKSYM.EQ.0.AND.JKSYM.EQ.0)
     .   CALL GA4(IK  ,JK  ,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB4(IKND,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB4(JKND,JK,KK,IJ,IK,VAL1,VAL2,KOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 303
 53           CONTINUE
              WRITE(NOUT,2013)
 2013         FORMAT(/     '  IN CONSTRUCTON OF THE A MATRIX, I=J=K IS N
     *OT POSSIBLE DO THE FIRST SET OF LOOPS'/'  ERROR IN LOOPS SOMEWHERE
     *')
              CALL MABORT
 63           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      IF(IKSYM.EQ.0)
     .   CALL GA6(IK  ,IK,KK,II,VAL1,VAL2,IOR,KOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB6(IKND,II,KK,VAL1,VAL2,IK,IOR,KOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 303
 73           CONTINUE
              WRITE(NOUT,2014)
 2014         FORMAT(/     '  IN CONSTRUCTON OF THE A MATRIX,I=J=K=L IS
     *NOT POSSIBLE DO THE FIRST SET OF LOOPS'/'  ERROR IN LOOPS SOMEWHER
     *E')
              CALL MABORT
 303          CONTINUE
 3430       CONTINUE
C3431 CONTINUE
 3420     CONTINUE
          DO 3480 KOR=LFERM,NDOC
            KO=IOFF(KOR)
            IK=IO+KOR
            JK=JO+KOR
            KK=KO+KOR
            IKND=IJORD(IK)
            JKND=IJORD(JK)
            KSM=ORBSYM(KOR)
            iksym=ieor((ism-1),(ksm-1))
            jksym=ieor((jsm-1),(ksm-1))
            JKEQ=JOR.EQ.KOR
            IJKT=IJT
            IF (JKEQ) IJKT=IJKT+1
            IJKEQ=IJEQ.AND.JKEQ
            IJKSYM=IEOR(IJSYM,(KSM-1))
            LOFF=IJKSYM*NORBS
            KAD=JAD+KADD(KOFF+KOR)
CTJL
       I2FL = 0
       IF(I2FL.EQ.1) GOTO 3491
            DO 3490 LOR=1,NUOC
         IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO 304
              LSM=ORBSYM(LOR)
              ilsym=ieor((ism-1),(lsm-1))
              jlsym=ieor((jsm-1),(lsm-1))
              klsym=ieor((ksm-1),(lsm-1))
              IF ((LSM-1) .NE. IJKSYM) GOTO 304
              LO=IOFF(LOR)
              IL=IO+LOR
              JL=JO+LOR
              KL=KO+LOR
              LL=LO+LOR
              ILND=IJORD(IL)
              JLND=IJORD(JL)
              KLND=IJORD(KL)
              IJKLT=IJKT
              IF (IJKEQ .AND. KOR.GT.LOR) IJKLT=IJKLT+1
              IF (KOR.EQ.LOR) IJKLT=IJKLT+3
              LAD=KAD+LADD(LOFF+LOR)
              IF (IJKLT.GT.7 .OR. IJKLT.LT.1) GOTO 9903
              IF (LAD .GT. NMAX) GOTO 9902
C
C   SPE-SPE-DOC-UOC
C
              GOTO (14,24,34,44,54,64,74)IJKLT
 14           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
C
      if(ijsym.eq.0.and.klsym.eq.0)
     1    CALL GA1(IJ  ,KL  ,IL,JL,IK,JK,VAL2,VAL1,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
      if(iksym.eq.0.and.jlsym.eq.0)
     1    CALL GA1(IK  ,JL  ,IL,KL,IJ,JK,VAL1,VAL2,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
      if(ilsym.eq.0.and.jksym.eq.0)
     1    CALL GA1(IL  ,JK  ,IK,KL,IJ,JL,VAL3,VAL2,VAL1,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB1(IJND,IL,JL,IK,JK,VAL2,VAL1,VAL3,KL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(IKND,IL,KL,IJ,JK,VAL1,VAL2,VAL3,JL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JLND,JK,KL,IJ,IL,VAL1,VAL2,VAL3,IK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(ILND,IJ,JL,IK,KL,VAL3,VAL2,VAL1,JK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(KLND,IK,IL,JK,JL,VAL2,VAL1,VAL3,IJ,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JKND,IJ,IK,JL,KL,VAL3,VAL2,VAL1,IL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 304
 24           CONTINUE
              WRITE(NOUT,1100) 'J=K    ','SECOND     '
              CALL MABORT
  34          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
cets121890
      if(itcscf.eq.1) then
        IF(IOR.EQ.NSB) E11(LOR,KOR) = E11(LOR,KOR) + A2*VAL2 - VAL1
        IF(IOR.EQ.NSA) E22(LOR,KOR) = E22(LOR,KOR) + A2*VAL2 - VAL1
      endif
cets
      if(iksym.eq.0.and.ilsym.eq.0)
     1    CALL GA3(IK  ,IL  ,IL,KL,II,IK,VAL1,VAL2,IOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GBD(KLND,IK,IL,VAL2,VAL1,II,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB3(IKND,IL,KL,II,IK,VAL1,VAL2,IOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB3(ILND,IK,KL,II,IL,VAL1,VAL2,IOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 304
 44           CONTINUE
              WRITE(NOUT,1100) 'K=L    ','SECOND     '
              CALL MABORT
 54           CONTINUE
              WRITE(NOUT,1100) 'I=J=K  ','SECOND     '
              CALL MABORT
 64           CONTINUE
              WRITE(NOUT,1100) 'I=J;K=L','SECOND     '
              CALL MABORT
 74           CONTINUE
              WRITE(NOUT,1100) 'I=J=K=L','SECOND     '
              CALL MABORT
 304          CONTINUE
 3490       CONTINUE
 3491 CONTINUE
CTJL
       I3FL = 0
       IF(I3FL.EQ.1) GOTO 3521
            DO 3520 LOR=LFERM,KOR
         IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO 305
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO 305
              ilsym=ieor((ism-1),(lsm-1))
              jlsym=ieor((jsm-1),(lsm-1))
              klsym=ieor((ksm-1),(lsm-1))
              LO=IOFF(LOR)
              IL=IO+LOR
              JL=JO+LOR
              KL=KO+LOR
              LL=LO+LOR
              ILND=IJORD(IL)
              JLND=IJORD(JL)
              KLND=IJORD(KL)
              IJKLT=IJKT
              IF (IJKEQ .AND. KOR.GT.LOR) IJKLT=IJKLT+1
              IF (KOR.EQ.LOR) IJKLT=IJKLT+3
              LAD=KAD+LADD(LOFF+LOR)
              IF (IJKLT.GT.7 .OR. IJKLT.LT.1) GOTO 9903
              IF (LAD .GT. NMAX) GOTO 9902
CTJL
         KLDC = IJDEP(KL)
CTJL
C
C    SPE-SPE-DOC-DOC
C
              GOTO (15,25,35,45,55,65,75)IJKLT
 15           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
C
      if(iksym.eq.0.and.jlsym.eq.0)
     1    CALL GA1(IK  ,JL  ,IL,KL,IJ,JK,VAL1,VAL2,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
      if(ilsym.eq.0.and.jksym.eq.0)
     1    CALL GA1(IL  ,JK  ,IK,KL,IJ,JL,VAL3,VAL2,VAL1,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB1(IJND,IL,JL,IK,JK,VAL2,VAL1,VAL3,KL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(IKND,IL,KL,IJ,JK,VAL1,VAL2,VAL3,JL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JLND,JK,KL,IJ,IL,VAL1,VAL2,VAL3,IK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(ILND,IJ,JL,IK,KL,VAL3,VAL2,VAL1,JK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JKND,IJ,IK,JL,KL,VAL3,VAL2,VAL1,IL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
      IF(KLDC.NE.0) CALL EX1(IJND,LL,KK,LOR,KOR,VAL2,VAL1,VAL3,
     1                       rs,ex,eavf,lg,s)
      IF(KLDC.NE.0) CALL RS1(KLDC,A1,VAL2,VAL1,VAL3,IJ,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 305
 25           CONTINUE
              WRITE(NOUT,1100) 'J=K    ','THIRD      '
              CALL MABORT
 35           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
cets121890
      if(itcscf.eq.1) then
        IF(IOR.EQ.NSB) THEN
           E11(KOR,LOR) = E11(KOR,LOR) + A2*VAL2 - VAL1
           E11(LOR,KOR) = E11(LOR,KOR) + A2*VAL2 - VAL1
        CALL DCICOF(DERC1,KK,LL,LOR,KOR,VAL2,VAL1,
     +     rs,ex,eavf,lg,s)
        END IF
        IF(IOR.EQ.NSA) THEN
           E22(KOR,LOR) = E22(KOR,LOR) + A2*VAL2 - VAL1
           E22(LOR,KOR) = E22(LOR,KOR) + A2*VAL2 - VAL1
        CALL DCICOF(DERC2,KK,LL,LOR,KOR,VAL2,VAL1,
     +     rs,ex,eavf,lg,s)
        END IF
      endif
cets
C
      if(iksym.eq.0.and.ilsym.eq.0)
     1   CALL GA3(IK  ,IL  ,IL,KL,II,IK,VAL1,VAL2,IOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB3(IKND,IL,KL,II,IK,VAL1,VAL2,IOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB3(ILND,IK,KL,II,IL,VAL1,VAL2,IOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
      IF(KLDC.NE.0) CALL RS1(KLDC,A1S2,VAL2,VAL1,VAL1,IJ,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 305
 45           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
cets121890
      if(itcscf.eq.1) then
        IF(IOR.EQ.NSB) THEN
           E11(JOR,IOR) = E11(JOR,IOR) + A2*VAL2 - VAL1
        END IF
        IF(JOR.EQ.NSA) THEN
           E22(IOR,JOR) = E22(IOR,JOR) + A2*VAL2 - VAL1
        END IF
      endif
cets
C
      if(iksym.eq.0.and.jksym.eq.0)
     1   CALL GA4(IK  ,JK  ,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GBD(IJND,IK,JK,VAL2,VAL1,KK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB4(IKND,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB4(JKND,JK,KK,IJ,IK,VAL1,VAL2,KOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 305
 55           CONTINUE
              WRITE(NOUT,1100) 'I=J=K  ','THIRD     '
              CALL MABORT
 65           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
cets121890
      if(itcscf.eq.1) then
       IF(IOR.EQ.NSB) E11(KOR,KOR) = E11(KOR,KOR) + A2*VAL2 - VAL1
       IF(IOR.EQ.NSA) E22(KOR,KOR) = E22(KOR,KOR) + A2*VAL2 - VAL1
       IF(IOR.EQ.NSB) E11(IOR,IOR) = E11(IOR,IOR) + A2*VAL2 - VAL1
       IF(IOR.EQ.NSA) E22(IOR,IOR) = E22(IOR,IOR) + A2*VAL2 - VAL1
C
       IF(IOR.EQ.NSB) H11 = H11 + A4*VAL2 - A2*VAL1
       IF(IOR.EQ.NSA) H22 = H22 + A4*VAL2 - A2*VAL1
      endif
cets
C
      if(iksym.eq.0)
     1   CALL GA6(IK  ,IK,KK,II,VAL1,VAL2,IOR,KOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB6(IKND,II,KK,VAL1,VAL2,IK,IOR,KOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 305
 75           CONTINUE
              WRITE(NOUT,1100) 'I=J=K=L','THIRD     '
              CALL MABORT
 305          CONTINUE
 3520       CONTINUE
 3521 CONTINUE
 3480     CONTINUE
          IF(IAND(IPRINT,1024).NE.0) WRITE(NOUT,1505)
 1505     FORMAT(/'  ABOUT TO JUMP TO 520')
C
          IF (NALP.EQ.0) GOTO 520
C
          DO 3570 KOR=NOP,NALP
            KO=IOFF(KOR)
            IK=IO+KOR
            JK=JO+KOR
            KK=KO+KOR
            IKND=IJORD(IK)
            JKND=IJORD(JK)
            KSM=ORBSYM(KOR)
            iksym=ieor((ism-1),(ksm-1))
            jksym=ieor((jsm-1),(ksm-1))
            JKEQ=JOR.EQ.KOR
            IJKT=IJT
            IF (JKEQ) IJKT=IJKT+1
            IJKEQ=IJEQ.AND.JKEQ
            IJKSYM=IEOR(IJSYM,(KSM-1))
            LOFF=IJKSYM*NORBS
            KAD=JAD+KADD(KOFF+KOR)
            DO 3580 LOR=1,NDOC
         IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO 306
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO 306
              ilsym=ieor((ism-1),(lsm-1))
              jlsym=ieor((jsm-1),(lsm-1))
              klsym=ieor((ksm-1),(lsm-1))
              LO=IOFF(LOR)
              IL=IO+LOR
              JL=JO+LOR
              KL=KO+LOR
              LL=LO+LOR
              ILND=IJORD(IL)
              JLND=IJORD(JL)
              KLND=IJORD(KL)
              IJKLT=IJKT
              IF (IJKEQ .AND. KOR.GT.LOR) IJKLT=IJKLT+1
              IF (KOR.EQ.LOR) IJKLT=IJKLT+3
              LAD=KAD+LADD(LOFF+LOR)
              IF (IJKLT.GT.7 .OR. IJKLT.LT.1) GOTO 9903
              IF (LAD .GT. NMAX) GOTO 9902
C
C    SPE-SPE-ALP-UOC/DOC
C
              GOTO (16,26,36,46,56,66,76)IJKLT
 16           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
C
      if(ijsym.eq.0.and.klsym.eq.0)
     1   CALL GA1(IJ  ,KL  ,IL,JL,IK,JK,VAL2,VAL1,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
      if(iksym.eq.0.and.jlsym.eq.0)
     1   CALL GA1(IK  ,JL  ,IL,KL,IJ,JK,VAL1,VAL2,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
      if(ilsym.eq.0.and.jksym.eq.0)
     1   CALL GA1(IL  ,JK  ,IK,KL,IJ,JL,VAL3,VAL2,VAL1,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB1(IJND,IL,JL,IK,JK,VAL2,VAL1,VAL3,KL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(KLND,IK,IL,JK,JL,VAL2,VAL1,VAL3,IJ,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(IKND,IL,KL,IJ,JK,VAL1,VAL2,VAL3,JL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JLND,JK,KL,IJ,IL,VAL1,VAL2,VAL3,IK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(ILND,IJ,JL,IK,KL,VAL3,VAL2,VAL1,JK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JKND,IJ,IK,JL,KL,VAL3,VAL2,VAL1,IL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 306
 26           CONTINUE
              WRITE(NOUT,1100) 'J=K    ','FOURTH     '
              CALL MABORT
  36          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      if(iksym.eq.0.and.ilsym.eq.0)
     1   CALL GA3(IK  ,IL  ,IL,KL,II,IK,VAL1,VAL2,IOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GBD(KLND,IK,IL,VAL2,VAL1,II,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB3(IKND,IL,KL,II,IK,VAL1,VAL2,IOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB3(ILND,IK,KL,II,IL,VAL1,VAL2,IOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 306
 46           CONTINUE
              WRITE(NOUT,1100) 'K=L    ','FOURTH     '
              CALL MABORT
 56           CONTINUE
              WRITE(NOUT,1100) 'I=J=K  ','FOURTH     '
              CALL MABORT
 66           CONTINUE
              WRITE(NOUT,1100) 'I=J;K=L','FOURTH     '
              CALL MABORT
 76           CONTINUE
              WRITE(NOUT,1100) 'I=J=K=L','FOURTH     '
              CALL MABORT
 306          CONTINUE
 3580       CONTINUE
            DO 3610 LOR=NOP,KOR
         IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO 307
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO 307
              ilsym=ieor((ism-1),(lsm-1))
              jlsym=ieor((jsm-1),(lsm-1))
              klsym=ieor((ksm-1),(lsm-1))
              LO=IOFF(LOR)
              IL=IO+LOR
              JL=JO+LOR
              KL=KO+LOR
              LL=LO+LOR
              ILND=IJORD(IL)
              JLND=IJORD(JL)
              KLND=IJORD(KL)
              IJKLT=IJKT
              IF (IJKEQ .AND. KOR.GT.LOR) IJKLT=IJKLT+1
              IF (KOR.EQ.LOR) IJKLT=IJKLT+3
              LAD=KAD+LADD(LOFF+LOR)
              IF (IJKLT.GT.7 .OR. IJKLT.LT.1) GOTO 9903
              IF (LAD .GT. NMAX) GOTO 9902
C
C     SPE-SPE-ALP-ALP
C
              GOTO (17,27,37,47,57,67,77)IJKLT
 17           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
C
      if(iksym.eq.0.and.jlsym.eq.0)
     1   CALL GA1(IK  ,JL  ,IL,KL,IJ,JK,VAL1,VAL2,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
      if(ilsym.eq.0.and.jksym.eq.0)
     1   CALL GA1(IL  ,JK  ,IK,KL,IJ,JL,VAL3,VAL2,VAL1,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB1(IJND,IL,JL,IK,JK,VAL2,VAL1,VAL3,KL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(IKND,IL,KL,IJ,JK,VAL1,VAL2,VAL3,JL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JLND,JK,KL,IJ,IL,VAL1,VAL2,VAL3,IK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(ILND,IJ,JL,IK,KL,VAL3,VAL2,VAL1,JK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JKND,IJ,IK,JL,KL,VAL3,VAL2,VAL1,IL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 307
 27           CONTINUE
              WRITE(NOUT,1100) 'J=K    ','FIFTH      '
              CALL MABORT
 37           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      if(iksym.eq.0.and.ilsym.eq.0)
     1   CALL GA3(IK  ,IL  ,IL,KL,II,IK,VAL1,VAL2,IOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB3(IKND,IL,KL,II,IK,VAL1,VAL2,IOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB3(ILND,IK,KL,II,IL,VAL1,VAL2,IOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 307
 47           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      if(iksym.eq.0.and.jksym.eq.0)
     1   CALL GA4(IK  ,JK  ,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GBD(IJND,IK,JK,VAL2,VAL1,KK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB4(IKND,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB4(JKND,JK,KK,IJ,IK,VAL1,VAL2,KOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 307
 57           CONTINUE
              WRITE(NOUT,1100) 'I=J=K  ','FIFTH     '
              CALL MABORT
 67           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      if(iksym.eq.0)
     1  CALL GA6(IK  ,IK,KK,II,VAL1,VAL2,IOR,KOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB6(IK  ,II,KK,VAL1,VAL2,IK,IOR,KOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 307
 77           CONTINUE
              WRITE(NOUT,1100) 'I=J=K=L','FIFTH     '
              CALL MABORT
 307          CONTINUE
 3610       CONTINUE
 3570     CONTINUE
C
 520      CONTINUE
          DO 3660 KOR=NSA,JOR
            KO=IOFF(KOR)
            IK=IO+KOR
            JK=JO+KOR
            KK=KO+KOR
            IKND=IJORD(IK)
            JKND=IJORD(JK)
            KSM=ORBSYM(KOR)
            iksym=ieor((ism-1),(ksm-1))
            jksym=ieor((jsm-1),(ksm-1))
            JKEQ=JOR.EQ.KOR
            IJKT=IJT
            IF (JKEQ) IJKT=IJKT+1
            IJKEQ=IJEQ.AND.JKEQ
            IJKSYM=IEOR(IJSYM,(KSM-1))
            LOFF=IJKSYM*NORBS
            KAD=JAD+KADD(KOFF+KOR)
CTJL
       I6FL = 0
       IF(I6FL.EQ.1) GOTO 3671
            DO 3670 LOR=1,KOR
         IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO 508
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO 508
              ilsym=ieor((ism-1),(lsm-1))
              jlsym=ieor((jsm-1),(lsm-1))
              klsym=ieor((ksm-1),(lsm-1))
              LO=IOFF(LOR)
              IL=IO+LOR
              JL=JO+LOR
              KL=KO+LOR
              LL=LO+LOR
              ILND=IJORD(IL)
              JLND=IJORD(JL)
              KLND=IJORD(KL)
              IJKLT=IJKT
              IF (IJKEQ .AND. KOR.GT.LOR) IJKLT=IJKLT+1
              IF (KOR.EQ.LOR) IJKLT=IJKLT+3
              LAD=KAD+LADD(LOFF+LOR)
              IF (IJKLT.GT.7 .OR. IJKLT.LT.1) GOTO 9903
              IF (LAD .GT. NMAX) GOTO 9902
C
C       SPE-SPE-SPE-UOC/DOC/ALP/SPE
C
              GOTO (518,528,538,548,558,568,578)IJKLT
518           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
C
      if(ijsym.eq.0.and.klsym.eq.0)
     1  CALL GA1(IJ  ,KL  ,IL,JL,IK,JK,VAL2,VAL1,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
      if(iksym.eq.0.and.jlsym.eq.0)
     1  CALL GA1(IK  ,JL  ,IL,KL,IJ,JK,VAL1,VAL2,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
      if(ilsym.eq.0.and.jksym.eq.0)
     1  CALL GA1(IL  ,JK  ,IK,KL,IJ,JL,VAL3,VAL2,VAL1,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB1(IJND,IL,JL,IK,JK,VAL2,VAL1,VAL3,KL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(KLND,IK,IL,JK,JL,VAL2,VAL1,VAL3,IJ,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(IKND,IL,KL,IJ,JK,VAL1,VAL2,VAL3,JL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JLND,JK,KL,IJ,IL,VAL1,VAL2,VAL3,IK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(ILND,IJ,JL,IK,KL,VAL3,VAL2,VAL1,JK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JKND,IJ,IK,JL,KL,VAL3,VAL2,VAL1,IL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 508
528           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
cets121890
      if(itcscf.eq.1) then
      E12(LOR,IOR) =  A1S2*VAL1
C
        IF(LOR.GT.NUOC.AND.LOR.LT.NSA) THEN
          IF(JOR.EQ.NSB) E11(IOR,LOR) = E11(IOR,LOR) + A2*VAL2 - VAL1
          IF(JOR.EQ.NSA) E22(IOR,LOR) = E22(IOR,LOR) + A2*VAL2 - VAL1
        END IF
      endif
cets
C
      if(ijsym.eq.0.and.jlsym.eq.0)
     1  CALL GA2(IJ  ,JL  ,IL,JL,IJ,JJ,VAL1,VAL2,JOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GBD(ILND,IJ,JL,VAL2,VAL1,JJ,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB22(IJND,IL,JL,IJ,JJ,VAL1,VAL2,JOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB21(JLND,JJ,JL,IJ,IL,VAL1,VAL2,JOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 508
538           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
cets121890
      if(itcscf.eq.1) then
        E12(LOR,KOR) =  A1S2*VAL1
C
        IF(LOR.GT.NUOC.AND.LOR.LT.NSA) THEN
          IF(IOR.EQ.NSB) E11(KOR,LOR) = E11(KOR,LOR) + A2*VAL2 - VAL1
          IF(IOR.EQ.NSA) E22(KOR,LOR) = E22(KOR,LOR) + A2*VAL2 - VAL1
        END IF
      endif
cets
C
      if(iksym.eq.0.and.ilsym.eq.0)
     1  CALL GA3(IK  ,IL  ,IL,KL,II,IK,VAL1,VAL2,IOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GBD(KLND,IK,IL,VAL2,VAL1,II,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB3(IKND,IL,KL,II,IK,VAL1,VAL2,IOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB3(ILND,IK,KL,II,IL,VAL1,VAL2,IOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 508
548           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      if(iksym.eq.0.and.jksym.eq.0)
     1  CALL GA4(IK  ,JK  ,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GBD(IJND,IK,JK,VAL2,VAL1,KK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB4(IKND,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB4(JKND,JK,KK,IJ,IK,VAL1,VAL2,KOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 508
558           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
cets121890
      if(itcscf.eq.1) then
              VAL3=INTS(LAD+3)
C
      IF(LOR.EQ.NSA.OR.LOR.EQ.NSB) THEN
      E12(IOR,LOR) =  A1S2*VAL1
      E12(LOR,IOR) =  A1S2*VAL2
      END IF
C
      IF(LOR.EQ.NSA) THEN
      E11(LOR,IOR) =  VAL1 + VAL3
      E22(IOR,LOR) =  VAL2 + VAL3
      END IF
C
      IF(LOR.GT.NUOC.AND.LOR.LT.NSA) THEN
        IF(IOR.EQ.NSB) E11(KOR,LOR) = E11(KOR,LOR) + VAL1 + VAL3 + VAL2
        IF(IOR.EQ.NSA) E22(KOR,LOR) = E22(KOR,LOR) + VAL1 + VAL3 + VAL2
        IF(IOR.EQ.NSA) E11(IOR,LOR) = E11(IOR,LOR) + VAL3 + VAL2
        IF(IOR.EQ.NSB) E22(IOR,LOR) = E22(IOR,LOR) + VAL3 + VAL2
        IF(IOR.EQ.NSB) E11(LOR,KOR) = E11(LOR,KOR) + VAL1 + VAL3 + VAL2
        IF(IOR.EQ.NSA) E22(LOR,KOR) = E22(LOR,KOR) + VAL1 + VAL3 + VAL2
      END IF
C
      IF(LOR.LT.LFERM) THEN
CTJL  WRITE(*,*) ' IN S-S-S-U ; VAL1,VAL3 ',VAL1,VAL3
        IF(IOR.EQ.NSB) E11(LOR,KOR) = E11(LOR,KOR) + VAL1 + VAL3
        IF(IOR.EQ.NSA) E22(LOR,KOR) = E22(LOR,KOR) + VAL1 + VAL3
      END IF
      endif
cets121890
      CALL GBD(ILND,II,IL,VAL1,VAL1,II,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      IF(LOR.GT.NUOC) CALL GBD(ILND,IL,LL,VAL2,VAL2,LL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 508
568           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
cets121890
      if(itcscf.eq.1) then
        E12(LOR,LOR) =  A1S2*VAL1
        E12(IOR,IOR) =  A1S2*VAL1
C
        IF(IOR.EQ.NSB.AND.LOR.EQ.NSA) H12 = VAL1
      endif
cets
C
      if(iksym.eq.0)
     1   CALL GA6(IK  ,IK,KK,II,VAL1,VAL2,IOR,KOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB6(IKND,II,KK,VAL1,VAL2,IK,IOR,KOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 508
578           CONTINUE
cets121890
              VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
      if(itcscf.eq.1) then
        IF(IOR.EQ.NSB) ONE = ONE + A2*CI1*CI1*VAL2
        IF(IOR.EQ.NSA) ONE = ONE + A2*CI2*CI2*VAL2
C
        IF(IOR.EQ.NSB) E11(IOR,IOR) = E11(IOR,IOR) + VAL1 + VAL2
        IF(IOR.EQ.NSA) E22(IOR,IOR) = E22(IOR,IOR) + VAL1 + VAL2
CTJL  WRITE(*,*) ' IN SPE-SPE-SPE-SPE PART; IOR = ',IOR
        IF(IOR.EQ.NSB) H11 = H11 + VAL1 + A2*VAL2
        IF(IOR.EQ.NSA) H22 = H22 + VAL1 + A2*VAL2
      endif
cets
C  THIS TYPE (I=J=K=L) DOES NOT CONTRIBUTE
 508          CONTINUE
 3670       CONTINUE
 3671 CONTINUE
 3660     CONTINUE
 3410   CONTINUE
C
        IF (NALP.EQ.0) GOTO 530
        DO 3740 JOR=NALP,NOP,-1
          JO=IOFF(JOR)
          IJ=IO+JOR
          JJ=JO+JOR
          IJND=IJORD(IJ)
          JSM=ORBSYM(JOR)
          IJEQ=IOR.EQ.JOR
          IJT=ITYP
          IF (IJEQ) IJT=IJT+2
          IJSYM=IEOR((ISM-1),(JSM-1))
          KOFF=IJSYM*NORBS
          JAD=IJADD(IJ)
          IF (IGRP .NE. IJGRP(IJ)) THEN
            CALL SREAD(IT52,INTS,NMAX*2)
            SRIN=SRIN+NMAX*2/112
            IGRP=IJGRP(IJ)
          END IF
          DO 3750 KOR=1,NUOC
            KO=IOFF(KOR)
            IK=IO+KOR
            JK=JO+KOR
            KK=KO+KOR
            IKND=IJORD(IK)
            JKND=IJORD(JK)
            KSM=ORBSYM(KOR)
            iksym=ieor((ism-1),(ksm-1))
            jksym=ieor((jsm-1),(ksm-1))
            JKEQ=JOR.EQ.KOR
            IJKT=IJT
            IF (JKEQ) IJKT=IJKT+1
            IJKEQ=IJEQ.AND.JKEQ
            IJKSYM=IEOR(IJSYM,(KSM-1))
            LOFF=IJKSYM*NORBS
            KAD=JAD+KADD(KOFF+KOR)
            DO 3760 LOR=1,KOR
         IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO 309
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO 309
              ilsym=ieor((ism-1),(lsm-1))
              jlsym=ieor((jsm-1),(lsm-1))
              klsym=ieor((ksm-1),(lsm-1))
              LO=IOFF(LOR)
              IL=IO+LOR
              JL=JO+LOR
              KL=KO+LOR
              LL=LO+LOR
              ILND=IJORD(IL)
              JLND=IJORD(JL)
              KLND=IJORD(KL)
              IJKLT=IJKT
              IF (IJKEQ .AND. KOR.GT.LOR) IJKLT=IJKLT+1
              IF (KOR.EQ.LOR) IJKLT=IJKLT+3
              LAD=KAD+LADD(LOFF+LOR)
              IF (IJKLT.GT.7 .OR. IJKLT.LT.1) GOTO 9903
              IF (LAD .GT. NMAX) GOTO 9902
CTJL
         KLVU = IJDEP(KL)
CTJL
C
C    SPE-ALP-UOC-UOC
C
              GOTO (19,29,39,49,59,69,79)IJKLT
 19           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
C
      if(iksym.eq.0.and.jlsym.eq.0)
     1  CALL GA1(IK  ,JL  ,IL,KL,IJ,JK,VAL1,VAL2,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
      if(ilsym.eq.0.and.jksym.eq.0)
     1  CALL GA1(IL  ,JK  ,IK,KL,IJ,JL,VAL3,VAL2,VAL1,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB1(IKND,IL,KL,IJ,JK,VAL1,VAL2,VAL3,JL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JLND,JK,KL,IJ,IL,VAL1,VAL2,VAL3,IK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(ILND,IJ,JL,IK,KL,VAL3,VAL2,VAL1,JK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JKND,IJ,IK,JL,KL,VAL3,VAL2,VAL1,IL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
      IF(KLVU.NE.0) CALL EX1(IJND,LL,KK,LOR,KOR,VAL2,VAL1,VAL3,
     1                       rs,ex,eavf,lg,s)
      IF(KLVU.NE.0) CALL RS1(KLVU,A1,VAL2,VAL1,VAL3,IJ,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 309
 29           CONTINUE
              WRITE(NOUT,1100) 'J=K    ','SEVENTH    '
              CALL MABORT
 39           CONTINUE
              WRITE(NOUT,1100) 'I=J    ','SEVENTH    '
              CALL MABORT
 49           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      if(iksym.eq.0.and.jksym.eq.0)
     1  CALL GA4(IK  ,JK  ,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB4(IKND,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB4(JKND,JK,KK,IJ,IK,VAL1,VAL2,KOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 309
 59           CONTINUE
              WRITE(NOUT,1100) 'I=J=K  ','SEVENTH    '
              CALL MABORT
 69           CONTINUE
              WRITE(NOUT,1100) 'I=J;K=L','SEVENTH    '
              CALL MABORT
 79           CONTINUE
              WRITE(NOUT,1100) 'I=J=K=L','SEVENTH    '
              CALL MABORT
 309          CONTINUE
 3760       CONTINUE
 3750     CONTINUE
          DO 3790 KOR=LFERM,NDOC
            KO=IOFF(KOR)
            IK=IO+KOR
            JK=JO+KOR
            KK=KO+KOR
            IKND=IJORD(IK)
            JKND=IJORD(JK)
            KSM=ORBSYM(KOR)
            iksym=ieor((ism-1),(ksm-1))
            jksym=ieor((jsm-1),(ksm-1))
            JKEQ=JOR.EQ.KOR
            IJKT=IJT
            IF (JKEQ) IJKT=IJKT+1
            IJKEQ=IJEQ.AND.JKEQ
            IJKSYM=IEOR(IJSYM,(KSM-1))
            LOFF=IJKSYM*NORBS
            KAD=JAD+KADD(KOFF+KOR)
            DO 3800 LOR=1,NUOC
         IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO 401
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO 401
              ilsym=ieor((ism-1),(lsm-1))
              jlsym=ieor((jsm-1),(lsm-1))
              klsym=ieor((ksm-1),(lsm-1))
              LO=IOFF(LOR)
              IL=IO+LOR
              JL=JO+LOR
              KL=KO+LOR
              LL=LO+LOR
              ILND=IJORD(IL)
              JLND=IJORD(JL)
              KLND=IJORD(KL)
              IJKLT=IJKT
              IF (IJKEQ .AND. KOR.GT.LOR) IJKLT=IJKLT+1
              IF (KOR.EQ.LOR) IJKLT=IJKLT+3
              LAD=KAD+LADD(LOFF+LOR)
              IF (IJKLT.GT.7 .OR. IJKLT.LT.1) GOTO 9903
              IF (LAD .GT. NMAX) GOTO 9902
C
C    SPE-ALP-DOC-UOC
C
              GOTO (111,121,131,141,151,161,171)IJKLT
 111          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
C
      if(ijsym.eq.0.and.klsym.eq.0)
     1   CALL GA1(IJ  ,KL  ,IL,JL,IK,JK,VAL2,VAL1,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
      if(iksym.eq.0.and.jlsym.eq.0)
     1   CALL GA1(IK  ,JL  ,IL,KL,IJ,JK,VAL1,VAL2,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
      if(ilsym.eq.0.and.jksym.eq.0)
     1   CALL GA1(IL  ,JK  ,IK,KL,IJ,JL,VAL3,VAL2,VAL1,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB1(IJND,IL,JL,IK,JK,VAL2,VAL1,VAL3,KL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(KLND,IK,IL,JK,JL,VAL2,VAL1,VAL3,IJ,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(IKND,IL,KL,IJ,JK,VAL1,VAL2,VAL3,JL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JLND,JK,KL,IJ,IL,VAL1,VAL2,VAL3,IK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(ILND,IJ,JL,IK,KL,VAL3,VAL2,VAL1,JK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JKND,IJ,IK,JL,KL,VAL3,VAL2,VAL1,IL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 401
 121          CONTINUE
              WRITE(NOUT,1100) 'J=K    ','EIGHTH     '
              CALL MABORT
 131          CONTINUE
              WRITE(NOUT,1100) 'I=J    ','EIGHTH     '
              CALL MABORT
 141          CONTINUE
              WRITE(NOUT,1100) 'K=L    ','EIGHTH     '
              CALL MABORT
 151          CONTINUE
              WRITE(NOUT,1100) 'I=J=K  ','EIGHTH     '
              CALL MABORT
 161          CONTINUE
              WRITE(NOUT,1100) 'I=J;K=L','EIGHTH     '
              CALL MABORT
 171          CONTINUE
              WRITE(NOUT,1100) 'I=J=K=L','EIGHTH     '
              CALL MABORT
 401          CONTINUE
 3800       CONTINUE
            DO 3820 LOR=LFERM,KOR
         IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO 402
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO 402
              ilsym=ieor((ism-1),(lsm-1))
              jlsym=ieor((jsm-1),(lsm-1))
              klsym=ieor((ksm-1),(lsm-1))
              LO=IOFF(LOR)
              IL=IO+LOR
              JL=JO+LOR
              KL=KO+LOR
              LL=LO+LOR
              ILND=IJORD(IL)
              JLND=IJORD(JL)
              KLND=IJORD(KL)
              IJKLT=IJKT
              IF (IJKEQ .AND. KOR.GT.LOR) IJKLT=IJKLT+1
              IF (KOR.EQ.LOR) IJKLT=IJKLT+3
              LAD=KAD+LADD(LOFF+LOR)
              IF (IJKLT.GT.7 .OR. IJKLT.LT.1) GOTO 9903
              IF (LAD .GT. NMAX) GOTO 9902
CTJL
         KLDC = IJDEP(KL)
CTJL
C
C   SPE-ALP-DOC-DOC
C
              GOTO (112,122,132,142,152,162,172)IJKLT
 112          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
C
      if(iksym.eq.0.and.jlsym.eq.0)
     1  CALL GA1(IK  ,JL  ,IL,KL,IJ,JK,VAL1,VAL2,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
      if(ilsym.eq.0.and.jksym.eq.0)
     1  CALL GA1(IL  ,JK  ,IK,KL,IJ,JL,VAL3,VAL2,VAL1,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB1(IJND,IL,JL,IK,JK,VAL2,VAL1,VAL3,KL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(IKND,IL,KL,IJ,JK,VAL1,VAL2,VAL3,JL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JLND,JK,KL,IJ,IL,VAL1,VAL2,VAL3,IK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(ILND,IJ,JL,IK,KL,VAL3,VAL2,VAL1,JK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JKND,IJ,IK,JL,KL,VAL3,VAL2,VAL1,IL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
      IF(KLDC.NE.0) CALL EX1(IJND,LL,KK,LOR,KOR,VAL2,VAL1,VAL3,
     1                       rs,ex,eavf,lg,s)
      IF(KLDC.NE.0) CALL RS1(KLDC,A1,VAL2,VAL1,VAL3,IJ,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 402
 122          CONTINUE
              WRITE(NOUT,1100) 'J=K    ','NINTH      '
              CALL MABORT
 132          CONTINUE
              WRITE(NOUT,1100) 'I=J    ','NINTH      '
              CALL MABORT
 142          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      if(iksym.eq.0.and.jksym.eq.0)
     1   CALL GA4(IK  ,JK  ,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GBD(IJND,IK,JK,VAL2,VAL1,KK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB4(IKND,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB4(JKND,JK,KK,IJ,IK,VAL1,VAL2,KOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 402
152           CONTINUE
              WRITE(NOUT,1100) 'I=J=K  ','NINTH      '
              CALL MABORT
162           CONTINUE
              WRITE(NOUT,1100) 'I=J;K=L','NINTH      '
              CALL MABORT
172           CONTINUE
              WRITE(NOUT,1100) 'I=J=K=L','NINTH      '
              CALL MABORT
 402          CONTINUE
 3820       CONTINUE
 3790     CONTINUE
C
          DO 3850 KOR=NOP,JOR
            KO=IOFF(KOR)
            IK=IO+KOR
            JK=JO+KOR
            KK=KO+KOR
            IKND=IJORD(IK)
            JKND=IJORD(JK)
            KSM=ORBSYM(KOR)
            iksym=ieor((ism-1),(ksm-1))
            jksym=ieor((jsm-1),(ksm-1))
            JKEQ=JOR.EQ.KOR
            IJKT=IJT
            IF (JKEQ) IJKT=IJKT+1
            IJKEQ=IJEQ.AND.JKEQ
            IJKSYM=IEOR(IJSYM,(KSM-1))
            LOFF=IJKSYM*NORBS
            KAD=JAD+KADD(KOFF+KOR)
            DO 3860 LOR=1,NDOC
         IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO 403
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO 403
              ilsym=ieor((ism-1),(lsm-1))
              jlsym=ieor((jsm-1),(lsm-1))
              klsym=ieor((ksm-1),(lsm-1))
              LO=IOFF(LOR)
              IL=IO+LOR
              JL=JO+LOR
              KL=KO+LOR
              LL=LO+LOR
              ILND=IJORD(IL)
              JLND=IJORD(JL)
              KLND=IJORD(KL)
              IJKLT=IJKT
              IF (IJKEQ .AND. KOR.GT.LOR) IJKLT=IJKLT+1
              IF (KOR.EQ.LOR) IJKLT=IJKLT+3
              LAD=KAD+LADD(LOFF+LOR)
              IF (IJKLT.GT.7 .OR. IJKLT.LT.1) GOTO 9903
              IF (LAD .GT. NMAX) GOTO 9902
C
C     SPE-ALP-ALP-UOC/DOC
C
              GOTO (113,123,133,143,153,163,173)IJKLT
 113          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
C
      if(ijsym.eq.0.and.klsym.eq.0)
     1  CALL GA1(IJ  ,KL  ,IL,JL,IK,JK,VAL2,VAL1,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
      if(iksym.eq.0.and.jlsym.eq.0)
     1  CALL GA1(IK  ,JL  ,IL,KL,IJ,JK,VAL1,VAL2,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB1(IJND,IL,JL,IK,JK,VAL2,VAL1,VAL3,KL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(IKND,IL,KL,IJ,JK,VAL1,VAL2,VAL3,JL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JLND,JK,KL,IJ,IL,VAL1,VAL2,VAL3,IK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(ILND,IJ,JL,IK,KL,VAL3,VAL2,VAL1,JK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(KLND,IK,IL,JK,JL,VAL2,VAL1,VAL3,IJ,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 403
 123          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      if(ijsym.eq.0.and.jlsym.eq.0)
     1   CALL GA2(IJ  ,JL  ,IL,JL,IJ,JJ,VAL1,VAL2,JOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GBD(ILND,IJ,JL,VAL2,VAL1,JJ,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB22(IJND,IL,JL,IJ,JJ,VAL1,VAL2,JOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB21(JLND,JJ,JL,IJ,IL,VAL1,VAL2,JOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 403
 133          CONTINUE
              WRITE(NOUT,1100) 'I=J    ','TENTH      '
              CALL MABORT
 143          CONTINUE
              WRITE(NOUT,1100) 'K=L    ','TENTH      '
              CALL MABORT
 153          CONTINUE
              WRITE(NOUT,1100) 'I=J=K  ','TENTH      '
              CALL MABORT
 163          CONTINUE
              WRITE(NOUT,1100) 'I=J;K=L','TENTH      '
              CALL MABORT
 173          CONTINUE
              WRITE(NOUT,1100) 'I=J=K=L','TENTH      '
              CALL MABORT
 403          CONTINUE
 3860       CONTINUE
C
            DO 3890 LOR=NOP,KOR
         IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO 603
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO 603
              ilsym=ieor((ism-1),(lsm-1))
              jlsym=ieor((jsm-1),(lsm-1))
              klsym=ieor((ksm-1),(lsm-1))
              LO=IOFF(LOR)
              IL=IO+LOR
              JL=JO+LOR
              KL=KO+LOR
              LL=LO+LOR
              ILND=IJORD(IL)
              JLND=IJORD(JL)
              KLND=IJORD(KL)
              IJKLT=IJKT
              IF (IJKEQ .AND. KOR.GT.LOR) IJKLT=IJKLT+1
              IF (KOR.EQ.LOR) IJKLT=IJKLT+3
              LAD=KAD+LADD(LOFF+LOR)
              IF (IJKLT.GT.7 .OR. IJKLT.LT.1) GOTO 9903
              IF (LAD .GT. NMAX) GOTO 9902
C
C    SPE-ALP-ALP-ALP
C
              GOTO (613,623,633,643,653,663,673)IJKLT
 613          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
      CALL GB1(IJND,IL,JL,IK,JK,VAL2,VAL1,VAL3,KL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(IKND,IL,KL,IJ,JK,VAL1,VAL2,VAL3,JL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(ILND,IJ,JL,IK,KL,VAL3,VAL2,VAL1,JK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 603
 623          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
      CALL GBD(ILND,IJ,JL,VAL2,VAL1,JJ,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB22(IJND,IL,JL,IJ,JJ,VAL1,VAL2,JOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 603
 633          CONTINUE
              WRITE(NOUT,1100) 'I=J    ','TENTH+1    '
              CALL MABORT
 643          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
      CALL GBD(IJND,IK,JK,VAL2,VAL1,KK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB4(IKND,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 603
 653          CONTINUE
              WRITE(NOUT,1100) 'I=J=K  ','TENTH+1    '
              CALL MABORT
 663          CONTINUE
              WRITE(NOUT,1100) 'I=J;K=L','TENTH+1    '
              CALL MABORT
 673          CONTINUE
              WRITE(NOUT,1100) 'I=J=K=L','TENTH+1    '
              CALL MABORT
 603          CONTINUE
 3890       CONTINUE
 3850     CONTINUE
 3740   CONTINUE
C
 530    CONTINUE
C
        DO 3930 JOR=NDOC,LFERM,-1
          JO=IOFF(JOR)
          IJ=IO+JOR
          JJ=JO+JOR
          IJND=IJORD(IJ)
          JSM=ORBSYM(JOR)
          IJEQ=IOR.EQ.JOR
          IJT=ITYP
          IF (IJEQ) IJT=IJT+2
          IJSYM=IEOR((ISM-1),(JSM-1))
          KOFF=IJSYM*NORBS
          JAD=IJADD(IJ)
          IF (IGRP .NE. IJGRP(IJ)) THEN
            CALL SREAD(IT52,INTS,NMAX*2)
            SRIN=SRIN+NMAX*2/112
            IGRP=IJGRP(IJ)
          END IF
          DO 3940 KOR=1,NUOC
            KO=IOFF(KOR)
            IK=IO+KOR
            JK=JO+KOR
            KK=KO+KOR
            IKND=IJORD(IK)
            JKND=IJORD(JK)
            KSM=ORBSYM(KOR)
            iksym=ieor((ism-1),(ksm-1))
            jksym=ieor((jsm-1),(ksm-1))
            JKEQ=JOR.EQ.KOR
            IJKT=IJT
            IF (JKEQ) IJKT=IJKT+1
            IJKEQ=IJEQ.AND.JKEQ
            IJKSYM=IEOR(IJSYM,(KSM-1))
            LOFF=IJKSYM*NORBS
            KAD=JAD+KADD(KOFF+KOR)
CTJL
       I12FL = 0
       IF(I12FL.EQ.1) GOTO 3951
            DO 3950 LOR=1,KOR
         IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO 404
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO 404
              ilsym=ieor((ism-1),(lsm-1))
              jlsym=ieor((jsm-1),(lsm-1))
              klsym=ieor((ksm-1),(lsm-1))
              LO=IOFF(LOR)
              IL=IO+LOR
              JL=JO+LOR
              KL=KO+LOR
              LL=LO+LOR
              ILND=IJORD(IL)
              JLND=IJORD(JL)
              KLND=IJORD(KL)
              IJKLT=IJKT
              IF (IJKEQ .AND. KOR.GT.LOR) IJKLT=IJKLT+1
              IF (KOR.EQ.LOR) IJKLT=IJKLT+3
              LAD=KAD+LADD(LOFF+LOR)
              IF (IJKLT.GT.7 .OR. IJKLT.LT.1) GOTO 9903
              IF (LAD .GT. NMAX) GOTO 9902
CTJL
         KLVU = IJDEP(KL)
CTJL
C
C   SPE-DOC-UOC-UOC
C
              GOTO (114,124,134,144,154,164,174)IJKLT
 114          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
C
      if(iksym.eq.0.and.jlsym.eq.0)
     1   CALL GA1(IK  ,JL  ,IL,KL,IJ,JK,VAL1,VAL2,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
      if(ilsym.eq.0.and.jksym.eq.0)
     1  CALL GA1(IL  ,JK  ,IK,KL,IJ,JL,VAL3,VAL2,VAL1,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB1(IKND,IL,KL,IJ,JK,VAL1,VAL2,VAL3,JL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JLND,JK,KL,IJ,IL,VAL1,VAL2,VAL3,IK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(ILND,IJ,JL,IK,KL,VAL3,VAL2,VAL1,JK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JKND,IJ,IK,JL,KL,VAL3,VAL2,VAL1,IL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
      IF(KLVU.NE.0) CALL EX1(IJND,LL,KK,LOR,KOR,VAL2,VAL1,VAL3,
     1                       rs,ex,eavf,lg,s)
      IF(KLVU.NE.0) CALL RS1(KLVU,A3S2,VAL2,VAL1,VAL3,IJ,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 404
 124          CONTINUE
              WRITE(NOUT,1100) 'J=K    ','ELEVENTH   '
              CALL MABORT
 134          CONTINUE
              WRITE(NOUT,1100) 'I=J    ','ELEVENTH   '
              CALL MABORT
 144          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      if(iksym.eq.0.and.jksym.eq.0)
     1  CALL GA4(IK  ,JK  ,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB4(IKND,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB4(JKND,JK,KK,IJ,IK,VAL1,VAL2,KOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 404
 154          CONTINUE
              WRITE(NOUT,1100) 'I=J=K  ','ELEVENTH   '
              CALL MABORT
 164          CONTINUE
              WRITE(NOUT,1100) 'I=J;K=L','ELEVENTH   '
              CALL MABORT
 174          CONTINUE
              WRITE(NOUT,1100) 'I=J=K=L','ELEVENTH   '
              CALL MABORT
 404          CONTINUE
 3950       CONTINUE
 3951 CONTINUE
 3940     CONTINUE
          DO 3980 KOR=LFERM,JOR
            KO=IOFF(KOR)
            IK=IO+KOR
            JK=JO+KOR
            KK=KO+KOR
            IKND=IJORD(IK)
            JKND=IJORD(JK)
            KSM=ORBSYM(KOR)
            iksym=ieor((ism-1),(ksm-1))
            jksym=ieor((jsm-1),(ksm-1))
            JKEQ=JOR.EQ.KOR
            IJKT=IJT
            IF (JKEQ) IJKT=IJKT+1
            IJKEQ=IJEQ.AND.JKEQ
            IJKSYM=IEOR(IJSYM,(KSM-1))
            LOFF=IJKSYM*NORBS
            KAD=JAD+KADD(KOFF+KOR)
CTJL
         JKDC = IJDEP(JK)
CTJL
       I13FL = 0
       IF(I13FL.EQ.1) GOTO 3991
            DO 3990 LOR=1,NUOC
         IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO 405
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO 405
              ilsym=ieor((ism-1),(lsm-1))
              jlsym=ieor((jsm-1),(lsm-1))
              klsym=ieor((ksm-1),(lsm-1))
              LO=IOFF(LOR)
              IL=IO+LOR
              JL=JO+LOR
              KL=KO+LOR
              LL=LO+LOR
              ILND=IJORD(IL)
              JLND=IJORD(JL)
              KLND=IJORD(KL)
              IJKLT=IJKT
              IF (IJKEQ .AND. KOR.GT.LOR) IJKLT=IJKLT+1
              IF (KOR.EQ.LOR) IJKLT=IJKLT+3
              LAD=KAD+LADD(LOFF+LOR)
              IF (IJKLT.GT.7 .OR. IJKLT.LT.1) GOTO 9903
              IF (LAD .GT. NMAX) GOTO 9902
C
C   SPE-DOC-DOC-UOC
C
              GOTO (115,125,135,145,155,165,175)IJKLT
 115          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
C
      if(ijsym.eq.0.and.klsym.eq.0)
     1   CALL GA1(IJ  ,KL  ,IL,JL,IK,JK,VAL2,VAL1,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
      if(iksym.eq.0.and.jlsym.eq.0)
     1   CALL GA1(IK  ,JL  ,IL,KL,IJ,JK,VAL1,VAL2,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB1(IJND,IL,JL,IK,JK,VAL2,VAL1,VAL3,KL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(IKND,IL,KL,IJ,JK,VAL1,VAL2,VAL3,JL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JLND,JK,KL,IJ,IL,VAL1,VAL2,VAL3,IK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(ILND,IJ,JL,IK,KL,VAL3,VAL2,VAL1,JK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(KLND,IK,IL,JK,JL,VAL2,VAL1,VAL3,IJ,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
      IF(JKDC.NE.0) CALL RS1(JKDC,A1S2,VAL3,VAL1,VAL2,IL,
     1                       rs,ex,eavf,lg,s)
      IF(JKDC.NE.0) CALL EX1(ILND,KK,JJ,KOR,JOR,VAL3,VAL1,VAL2,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 405
 125          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
cets121890
      if(itcscf.eq.1) then
        IF(IOR.EQ.NSB) E11(LOR,IOR) = E11(LOR,IOR) + A2*VAL2 - VAL1
        IF(IOR.EQ.NSA) E22(LOR,IOR) = E22(LOR,IOR) + A2*VAL2 - VAL1
      endif
cets
C
      if(ijsym.eq.0.and.jlsym.eq.0)
     1   CALL GA2(IJ  ,JL  ,IL,JL,IJ,JJ,VAL1,VAL2,JOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GBD(ILND,IJ,JL,VAL2,VAL1,JJ,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB21(JLND,JJ,JL,IJ,IL,VAL1,VAL2,JOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB22(IJND,IL,JL,IJ,JJ,VAL1,VAL2,JOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
              GOTO 405
 135          CONTINUE
              WRITE(NOUT,1100) 'I=J    ','TWELVTH    '
              CALL MABORT
 145          CONTINUE
              WRITE(NOUT,1100) 'K=L    ','TWELVTH    '
              CALL MABORT
 155          CONTINUE
              WRITE(NOUT,1100) 'I=J=K  ','TWELVTH    '
              CALL MABORT
 165          CONTINUE
              WRITE(NOUT,1100) 'I=J;K=L','TWELVTH    '
              CALL MABORT
 175          CONTINUE
              WRITE(NOUT,1100) 'I=J=K=L','TWELVTH    '
              CALL MABORT
 405          CONTINUE
 3990       CONTINUE
 3991 CONTINUE
CTJL
       I14FL = 0
       IF(I14FL.EQ.1) GOTO 4021
            DO 4020 LOR=LFERM,KOR
         IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO 605
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO 605
              LO=IOFF(LOR)
              IL=IO+LOR
              JL=JO+LOR
              KL=KO+LOR
              LL=LO+LOR
              ILND=IJORD(IL)
              JLND=IJORD(JL)
              KLND=IJORD(KL)
              IJKLT=IJKT
              IF (IJKEQ .AND. KOR.GT.LOR) IJKLT=IJKLT+1
              IF (KOR.EQ.LOR) IJKLT=IJKLT+3
              LAD=KAD+LADD(LOFF+LOR)
              IF (IJKLT.GT.7 .OR. IJKLT.LT.1) GOTO 9903
              IF (LAD .GT. NMAX) GOTO 9902
CTJL
         JLDC = IJDEP(JL)
         KLDC = IJDEP(KL)
CTJL
C
C    SPE-DOC-DOC-DOC
C
              GOTO (615,625,635,645,655,665,675)IJKLT
 615          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
      CALL GB1(IJND,IL,JL,IK,JK,VAL2,VAL1,VAL3,KL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(IKND,IL,KL,IJ,JK,VAL1,VAL2,VAL3,JL,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(ILND,IJ,JL,IK,KL,VAL3,VAL2,VAL1,JK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
      IF(JKDC.NE.0) CALL RS1(JKDC,A2,VAL3,VAL1,VAL2,IL,
     1                       rs,ex,eavf,lg,s)
      IF(JLDC.NE.0) CALL RS1(JLDC,A2,VAL1,VAL1,VAL2,IK,
     1                       rs,ex,eavf,lg,s)
      IF(KLDC.NE.0) CALL RS1(KLDC,A2,VAL2,VAL1,VAL3,IJ,
     1                       rs,ex,eavf,lg,s)
      IF(JKDC.NE.0) CALL EX1(ILND,KK,JJ,KOR,JOR,VAL3,VAL1,VAL2,
     1                       rs,ex,eavf,lg,s)
      IF(JLDC.NE.0) CALL EX1(IKND,LL,JJ,LOR,JOR,VAL1,VAL3,VAL2,
     1                       rs,ex,eavf,lg,s)
      IF(KLDC.NE.0) CALL EX1(IJND,LL,KK,LOR,KOR,VAL2,VAL1,VAL3,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 605
 625          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
cets121890
      if(itcscf.eq.1) then
        E11(IOR,LOR) = E11(IOR,LOR) + A2*VAL2 - VAL1
        E22(IOR,LOR) = E22(IOR,LOR) + A2*VAL2 - VAL1
        IF(IOR.EQ.NSB) E11(LOR,IOR) = E11(LOR,IOR) + A2*VAL2 - VAL1
        IF(IOR.EQ.NSA) E22(LOR,IOR) = E22(LOR,IOR) + A2*VAL2 - VAL1
      endif
cets
      CALL GBD(ILND,IJ,JL,VAL2,VAL1,JJ,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB22(IJND,IL,JL,IJ,JJ,VAL1,VAL2,JOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
      IF(JLDC.NE.0) CALL RS1(JLDC,A1,VAL1,VAL1,VAL2,IK,
     1                       rs,ex,eavf,lg,s)
      IF(JLDC.NE.0) CALL EX1(IKND,LL,JJ,LOR,JOR,VAL1,VAL1,VAL2,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 605
 635          CONTINUE
              WRITE(NOUT,1100) 'I=J    ','TWELVTH+1  '
              CALL MABORT
 645          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
cets121890
      if(itcscf.eq.1) then
        E11(IOR,JOR) = E11(IOR,JOR) + A2*VAL2 - VAL1
        E22(IOR,JOR) = E22(IOR,JOR) + A2*VAL2 - VAL1
        IF(IOR.EQ.NSB) E11(JOR,IOR) = E11(JOR,IOR) + A2*VAL2 - VAL1
        IF(IOR.EQ.NSA) E22(JOR,IOR) = E22(JOR,IOR) + A2*VAL2 - VAL1
      endif
cets
      CALL GBD(IJND,IK,JK,VAL2,VAL1,KK,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB4(IKND,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     1   a,eta,zeta,alpha,beta,s,ntypes+1)
C
      IF(JKDC.NE.0) CALL RS1(JKDC,A1,VAL1,VAL1,VAL2,IL,
     1                       rs,ex,eavf,lg,s)
ctph  IF(JLDC.NE.0) CALL EX1(IKND,LL,JJ,LOR,JOR,VAL1,VAL1,VAL2,
      IF(JKDC.NE.0) CALL EX1(ILND,KK,JJ,KOR,JOR,VAL1,VAL1,VAL2,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 605
 655       CONTINUE
C THERE IS NO CONTRIBUTION HERE; DONE IN SPE-SPE-SPE-ALL SECTION
              WRITE(NOUT,1100) 'I=J=K  ','TWELVTH+1  '
              CALL MABORT
              GOTO 605
 665          CONTINUE
              WRITE(NOUT,1100) 'I=J;K=L','TWELVTH+1  '
              CALL MABORT
 675          CONTINUE
              WRITE(NOUT,1100) 'I=J=K=L','TWELVTH    '
              CALL MABORT
 605          CONTINUE
 4020       CONTINUE
 4021 CONTINUE
 3980     CONTINUE
 3930   CONTINUE
C     SKIP I,UOC,UOC,UOC PART
        DO 4060 JOR=NUOC,1,-1
          JO=IOFF(JOR)
          IJ=IO+JOR
          JJ=JO+JOR
          IJND=IJORD(IJ)
          JSM=ORBSYM(JOR)
          IJEQ=IOR.EQ.JOR
          IJT=ITYP
          IF (IJEQ) IJT=IJT+2
          IJSYM=IEOR((ISM-1),(JSM-1))
          KOFF=IJSYM*NORBS
          JAD=IJADD(IJ)
CTJL
       I20FL = 0
       IF(I20FL.EQ.1) GOTO 8411
C
          IF(IGRP.NE.IJGRP(IJ)) THEN
            CALL SREAD(IT52,INTS,NMAX*2)
            SRIN=SRIN+NMAX*2/112
            IGRP=IJGRP(IJ)
          END IF
C
          DO 8410 KOR=1,JOR
            KO=IOFF(KOR)
            IK=IO+KOR
            JK=JO+KOR
            KK=KO+KOR
            IKND=IJORD(IK)
            JKND=IJORD(JK)
            KSM=ORBSYM(KOR)
            JKEQ=JOR.EQ.KOR
            IJKT=IJT
            IF (JKEQ) IJKT=IJKT+1
            IJKEQ=IJEQ.AND.JKEQ
            IJKSYM=IEOR(IJSYM,(KSM-1))
            LOFF=IJKSYM*NORBS
            KAD=JAD+KADD(KOFF+KOR)
CTJL
            JKVU = IJDEP(JK)
CTJL
C
C     SPE-UOC-UOC-UOC
C
            DO 8420 LOR=1,KOR
          IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO  8430
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO 8430
              LO=IOFF(LOR)
              IL=IO+LOR
              JL=JO+LOR
              KL=KO+LOR
              LL=LO+LOR
              ILND=IJORD(IL)
              JLND=IJORD(JL)
              KLND=IJORD(KL)
              IJKLT=IJKT
              IF (IJKEQ .AND. KOR.GT.LOR) IJKLT=IJKLT+1
              IF (KOR.EQ.LOR) IJKLT=IJKLT+3
              LAD=KAD+LADD(LOFF+LOR)
              IF (IJKLT.GT.7 .OR. IJKLT.LT.1) GOTO 9903
              IF (LAD .GT. NMAX) GOTO 9902
              JULIA(IJKLT+56)=JULIA(IJKLT+56)+1
CTJL
          JLVU = IJDEP(JL)
          KLVU = IJDEP(KL)
CTJL
              GOTO (8412,8422,8432,8442,8452,8462,8472)IJKLT
 8412         VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
CTJL
      IF(JKVU.NE.0) CALL EX1(ILND,KK,JJ,KOR,JOR,VAL3,VAL1,VAL2,
     1                       rs,ex,eavf,lg,s)
      IF(JLVU.NE.0) CALL EX1(IKND,LL,JJ,LOR,JOR,VAL1,VAL3,VAL2,
     1                       rs,ex,eavf,lg,s)
      IF(KLVU.NE.0) CALL EX1(IJND,LL,KK,LOR,KOR,VAL2,VAL1,VAL3,
     1                       rs,ex,eavf,lg,s)
C
      IF(JKVU.NE.0) CALL RS1(JKVU,A1S2,VAL3,VAL2,VAL1,IL,
     1                       rs,ex,eavf,lg,s)
      IF(JLVU.NE.0) CALL RS1(JLVU,A1S2,VAL1,VAL2,VAL3,IK,
     1                       rs,ex,eavf,lg,s)
      IF(KLVU.NE.0) CALL RS1(KLVU,A1S2,VAL2,VAL1,VAL3,IJ,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 8430
 8422         VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
      IF(JLVU.NE.0) CALL EX1(IKND,LL,JJ,LOR,JOR,VAL1,VAL1,VAL2,
     1                       rs,ex,eavf,lg,s)
C
      IF(JLVU.NE.0) CALL RS1(JLVU,A1S2,VAL1,VAL1,VAL2,IK,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 8430
 8432      WRITE(*,*) ' HIT I=J IN D-U-U-U ERROR IN LOOPS'
              GOTO 8430
 8442         VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
      IF(JKVU.NE.0) CALL EX1(ILND,KK,JJ,KOR,JOR,VAL1,VAL1,VAL2,
     1                       rs,ex,eavf,lg,s)
C
      IF(JKVU.NE.0) CALL RS1(JKVU,A1S2,VAL1,VAL1,VAL2,IL,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 8430
 8452      WRITE(*,*) ' HIT CASE 5 IN D-U-U-U ERROR IN LOOPS'
CTJL
CTJL THIS TERM , EITHER I=J=K OR J=K=L, DOES NOT CONTRIBUTE.
CTJL
              GOTO 8430
 8462      WRITE(*,*) ' HIT CASE 6 IN D-U-U-U ERROR IN LOOPS'
              GOTO 8430
 8472         CONTINUE
      WRITE(*,*) ' I=K=J=L ERROR IN LOOPS',IOR,JOR,KOR,LOR
CTJL          WRITE(NOUT,1100) 'I=J=K=L','NINTEENTH+1'
CTJL   I=J=J=L DOES NOT CONTRIBUTE
CTJL          CALL MABORT
 8430         CONTINUE
 8420       CONTINUE
 8410       CONTINUE
 8411   CONTINUE
 4060   CONTINUE
 3400 CONTINUE
CTJL
      RETURN
CTJL
 9902 CONTINUE
      WRITE(NOUT,2005)LAD
      CALL MABORT
 9903 CONTINUE
      WRITE(NOUT,2)IJKLT,IOR,JOR,KOR,LOR,IBF,JBF,KBF,LBF,ISM,JSM,KSM,LSM
      CALL MABORT
CTJL
      END
