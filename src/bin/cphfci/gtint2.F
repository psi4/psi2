C**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
C**********************************************************************
C***********************************************************************
c       by: Edward Seidl (cets121890)
c     date: December 18, 1990
c   reason: 1) finish changes to gtints started by tph
c           2) removed all entry points
c           3) changed dimension of zeta from numij*norbs to
c                numij*(ntypes+1) and changed getint, ga* and gb*
c                subroutines
c           4) added ability to do tcscf-ci (maybe)
C***********************************************************************
C*  BY:  TRACY HAMILTON                                               *
C*  DATE:  July 2, 1990                                               *
C*  REASON: CHANGE SUBROUTINE CALLS TO GA1, GA2, GA3, GA4 and GA6     *
C*          IN SUBROUTINE GTINTS TO BE THE SAME AS IN GTINT2          *
C*          SO OPEN-SHELL SINGLET CI GRADIENTS WILL WORK              *
C**********************************************************************
C*  BY:  RICHARD REMINGTON                         search:  c5-22-88  *
C*  DATE:  MAY   22,   1988                                           *
C*  REASON: DECREASE CORE TO RUN IN 7MB ON 9370                       *
C**********************************************************************
C     PROGRAM LAST UPDATED OCTOBER 6, 1987                            *
C     BY: GES                                                         *
C     REASON: INCLUDE OPEN-SHELL AND TCSCF.                           *
C**********************************************************************
C     PROGRAM LAST UPDATED JULY  7, 1987                              *
C     BY: GES                                                         *
C     REASON: VERSION FOR CPHFL  (LOOK FOR CGES)                      *
C**********************************************************************
C     PROGRAM LAST UPDATED APRIL 8, 1986                              *
C     BY: TJL                                                         *
C     REASON: FIX BUG IN D-D-D-D CASE 5.                              *
C     PROGRAM LAST UPDATED JULY 18, 1985                              *
C     BY: TJL                                                         *
C**********************************************************************
C
      SUBROUTINE GTINT2(A,ETA,ZETA,F,ALPHA,BETA,IJORD,KADD,LADD,
     1 IJADD,IJGRP,ORBSYM,INTS,S,EPA,JT1C,IOFF,EX,RS,LG,IJDEP,BOX,
     2 EAVF,LSTEPA,EX2,aa,doff,dadd,idim,idima,
     3 e11,e12,e22,ha11,ha22,ha12,nma1)
      IMPLICIT INTEGER (A-Z)
cets072091
#include <error.h>
      integer frdboo, errcod, idl52
      LOGICAL IJEQ,IJKEQ,JKEQ
      INTEGER JULIA(70),IJDEP(NUMIJ)
      REAL*8 A0,A1,A2,A1S2,A4,FOCC,EFAC,BOX(NDF),XFAC,SFAC,A3S2
      REAL*8 EAVF(NUMIJ),LSTEPA(NUMIJ),EX2(IJIND)
      real*8 aa(nma1,nma1+1)
      integer doff(idim),dadd(idim)
      DIMENSION IJORD(NUMIJ),KADD(SYMORB),LADD(SYMORB),JT1C(1)
      DIMENSION ORBSYM(NORBS),IJADD(NUMIJ),IJGRP(NUMIJ)
      DIMENSION IOFF(NUMIJ)
      REAL*8 A(IJIND,NDF),F(NORBS),INTS(NMAX),LG(NORBS,NORBS)
     1,S(NUMIJ,NDF),EPA(NORBS,NORBS),EX(IJIND),RS(NUMDEP,NDF)
cets121890*8 ZETA(NUMIJ,NORBS),ALPHA(NUMIJ),BETA(NUMIJ)
      REAL*8 zeta(ntypes+1,numij),ALPHA(NUMIJ),BETA(NUMIJ)
      real*8 ha11(ndf),ha22(ndf),ha12(ndf),eh11,eh22,eh12,etcsc1
      real*8 ci1,ci2,h11,h12,h22,e11(nbf,nbf),e12(nbf,nbf),e22(nbf,nbf)
      real*8 escf2,etcscf,derc1,derc2,one
      REAL*8 C1A,C2A,C1,C2
      real*8 bi1,bi2,saij1,saij2,etest,ea
      REAL*8 VAL1,VAL2,VAL3,ETA(NUMIJ),ESCF,REPNUC
      COMMON /CIA/ C1A(45),C2A(45),C1,C2,ESCF
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
 1100 FORMAT(/'  IN THE CONSTRUCTION OF THE A MATRIX, THE CASE ',A7
     * /'  IS NOT POSSIBLE IN THE ',A11,'  SET OF LOOPS.'/'  ERROR
     * IN LOOP STRUCTURE SOMEWHERE.')
   2  FORMAT(/'  IN GTINT2, IJKLT UNKNOWN TYPE ',I5/'  ORBITALS: ',4I5
     * /'  BASIS FUNCTIONS: ',4I5/,'  SYMMETRIES: ',4I5)
 2005 FORMAT(/'  IN GTINT2, LAD TOO BIG: ',I6)
C
      IJIND2 = IJIND
      ITEST = 0
      A0 = 0.0D+00
      A1 = 1.0D+00
      A2 = 2.0D+00
      A4 = 4.0D+00
      A1S2 = 0.50D+00
      A3S2 = 1.5D+00
      escf2=a0 
CTJL
C     WRITE(*,*) ' THE RS MATRIX '
C     CALL MATOUT(RS,NUMDEP,NDF,NUMDEP,NDF,NOUT)
C     WRITE(*,*) ' IN GTINT2 ; NOP,NBF ',NOP,NBF
C     WRITE(*,*) ' THE AVERAGED FOCK MATRIX '
C     CALL PRINT(EAVF,NUMIJ,NBF,NOUT)
CTJL
C
C THIS IS THE BEGINNING OF THE THIRD PART OF GETINT TJL 7/09/84
C
COUT  WRITE(*,*)'IN SUBROUTINE GTINT2,  NALP=',NALP
      IF (NALP.EQ.0) GOTO 540
C
      DO 4070 IOR=NALP,NOP,-1
        IO=IOFF(IOR)
        II=IO+IOR
        ISM=ORBSYM(IOR)
        ITYP=1
        DO 4080 JOR=IOR,NOP,-1
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
C
          IF (IGRP.NE.IJGRP(IJ)) THEN
            CALL SREAD(IT52,INTS,NMAX*2)
            SRIN=SRIN+NMAX*2/112
            IGRP=IJGRP(IJ)
          END IF
C
          DO 4090 KOR=1,NUOC
            KO=IOFF(KOR)
            IK=IO+KOR
            JK=JO+KOR
            KK=KO+KOR
            IKND=IJORD(IK)
            JKND=IJORD(JK)
            KSM=ORBSYM(KOR)
            IKSYM=IEOR((ISM-1),(KSM-1))
            JKSYM=IEOR((JSM-1),(KSM-1))
            JKEQ=JOR.EQ.KOR
            IJKT=IJT
            IF (JKEQ) IJKT=IJKT+1
            IJKEQ=IJEQ.AND.JKEQ
            IJKSYM=IEOR(IJSYM,(KSM-1))
            LOFF=IJKSYM*NORBS
            KAD=JAD+KADD(KOFF+KOR)
C
C  ALP-ALP-UOC-UOC
C
            DO 4100 LOR=1,KOR
         IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO  406
              LSM=ORBSYM(LOR)
              ILSYM=IEOR((ISM-1),(LSM-1))
              JLSYM=IEOR((JSM-1),(LSM-1))
              IF ((LSM-1) .NE. IJKSYM) GOTO 406
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
              JULIA(IJKLT)=JULIA(IJKLT)+1
CTJL
         KLVU = IJDEP(KL)
CTJL
              GOTO (116,126,136,146,156,166,176)IJKLT
 116          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
C
      IF(IKSYM.EQ.0.AND.JLSYM.EQ.0)
     .   CALL GA1(IK  ,JL  ,IL,KL,IJ,JK,VAL1,VAL2,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
      IF(ILSYM.EQ.0.AND.JKSYM.EQ.0)
     .   CALL GA1(IL  ,JK  ,IK,KL,IJ,JL,VAL3,VAL2,VAL1,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB1(IKND,IL,KL,IJ,JK,VAL1,VAL2,VAL3,JL,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JLND,JK,KL,IJ,IL,VAL1,VAL2,VAL3,IK,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(ILND,IJ,JL,IK,KL,VAL3,VAL2,VAL1,JK,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JKND,IJ,IK,JL,KL,VAL3,VAL2,VAL1,IL,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
C
      IF(KLVU.NE.0) CALL RS1(KLVU,A1,VAL2,VAL1,VAL3,IJ,
     1                       rs,ex,eavf,lg,s)
C
 4110      CONTINUE
           GO TO 406
 126          CONTINUE
              WRITE(NOUT,1100) 'J=K    ','THIRTEENTH '
              CALL MABORT
 136          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      IF(IKSYM.EQ.0.AND.ILSYM.EQ.0)
     .   CALL GA3(IK  ,IL  ,IL,KL,II,IK,VAL1,VAL2,IOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB3(IKND,IL,KL,II,IK,VAL1,VAL2,IOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB3(ILND,IK,KL,II,IL,VAL1,VAL2,IOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
C
      IF(KLVU.NE.0) CALL RS1(KLVU,A1S2,VAL2,VAL1,VAL1,IJ,
     1                       rs,ex,eavf,lg,s)
C
 4120      CONTINUE
              GOTO 406
 146          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      IF(IKSYM.EQ.0.AND.JKSYM.EQ.0)
     .   CALL GA4(IK  ,JK  ,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB4(IKND,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB4(JKND,JK,KK,IJ,IK,VAL1,VAL2,KOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
              GOTO 406
 156          CONTINUE
              WRITE(NOUT,1100) 'I=J=K  ','THIRTEENTH '
              CALL MABORT
 166          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      IF(IKSYM.EQ.0)
     .   CALL GA6(IK  ,IK,KK,II,VAL1,VAL2,IOR,KOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB6(IKND,II,KK,VAL1,VAL2,IK,IOR,KOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
              GOTO 406
 176          CONTINUE
              WRITE(NOUT,1100) 'I=J=K=L','THIRTEENTH '
              CALL MABORT
 406          CONTINUE
 4100       CONTINUE
 4090     CONTINUE
          DO 4150 KOR=LFERM,NDOC
            KO=IOFF(KOR)
            IK=IO+KOR
            JK=JO+KOR
            KK=KO+KOR
            IKND=IJORD(IK)
            JKND=IJORD(JK)
            KSM=ORBSYM(KOR)
            IKSYM=IEOR((ISM-1),(KSM-1))
            JKSYM=IEOR((JSM-1),(KSM-1))
            JKEQ=JOR.EQ.KOR
            IJKT=IJT
            IF (JKEQ) IJKT=IJKT+1
            IJKEQ=IJEQ.AND.JKEQ
            IJKSYM=IEOR(IJSYM,(KSM-1))
            LOFF=IJKSYM*NORBS
            KAD=JAD+KADD(KOFF+KOR)
C
C  ALP-ALP-DOC-UOC
C
            DO 4160 LOR=1,NUOC
         IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO 407
              LSM=ORBSYM(LOR)
              ILSYM=IEOR((ISM-1),(LSM-1))
              JLSYM=IEOR((JSM-1),(LSM-1))
              IF ((LSM-1) .NE. IJKSYM) GOTO 407
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
              JULIA(IJKLT+7)=JULIA(IJKLT+7)+1
              GOTO (117,127,137,147,157,167,177)IJKLT
 117          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
C
      IF(IKSYM.EQ.0.AND.JLSYM.EQ.0)
     .   CALL GA1(IK  ,JL  ,IL,KL,IJ,JK,VAL1,VAL2,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
      IF(ILSYM.EQ.0.AND.JKSYM.EQ.0)
     .   CALL GA1(IL  ,JK  ,IK,KL,IJ,JL,VAL3,VAL2,VAL1,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB1(KLND,IK,IL,JK,JL,VAL2,VAL1,VAL3,IJ,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(IKND,IL,KL,IJ,JK,VAL1,VAL2,VAL3,JL,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JLND,JK,KL,IJ,IL,VAL1,VAL2,VAL3,IK,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(ILND,IJ,JL,IK,KL,VAL3,VAL2,VAL1,JK,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JKND,IJ,IK,JL,KL,VAL3,VAL2,VAL1,IL,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
              GOTO 407
 127          CONTINUE
              WRITE(NOUT,1100) 'J=K    ','FOURTEENTH '
              CALL MABORT
 137          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      IF(IKSYM.EQ.0.AND.ILSYM.EQ.0)
     .   CALL GA3(IK  ,IL  ,IL,KL,II,IK,VAL1,VAL2,IOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GBD(KLND,IK,IL,VAL2,VAL1,II,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB3(IKND,IL,KL,II,IK,VAL1,VAL2,IOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB3(ILND,IK,KL,II,IL,VAL1,VAL2,IOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
              GOTO 407
 147          CONTINUE
              WRITE(NOUT,1100) 'K=L    ','FOURTEENTH '
              CALL MABORT
 157          CONTINUE
              WRITE(NOUT,1100) 'I=J=K  ','FOURTEENTH '
              CALL MABORT
 167          CONTINUE
              WRITE(NOUT,1100) 'I=J;K=L','FOURTEENTH '
              CALL MABORT
 177          CONTINUE
              WRITE(NOUT,1100) 'I=J=K=L','FOURTEENTH '
              CALL MABORT
 407          CONTINUE
 4160       CONTINUE
C
C ALP-ALP-DOC-DOC
C
            DO 4190 LOR=LFERM,KOR
         IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO 408
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO 408
              ILSYM=IEOR((ISM-1),(LSM-1))
              JLSYM=IEOR((JSM-1),(LSM-1))
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
              JULIA(IJKLT+14)=JULIA(IJKLT+14)+1
CTJL
            KLDC = IJDEP(KL)
CTJL
              GOTO (118,128,138,148,158,168,178)IJKLT
 118          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
C
      IF(IKSYM.EQ.0.AND.JLSYM.EQ.0)
     .   CALL GA1(IK  ,JL  ,IL,KL,IJ,JK,VAL1,VAL2,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
      IF(ILSYM.EQ.0.AND.JKSYM.EQ.0)
     .   CALL GA1(IL  ,JK  ,IK,KL,IJ,JL,VAL3,VAL2,VAL1,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB1(IKND,IL,KL,IJ,JK,VAL1,VAL2,VAL3,JL,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JLND,JK,KL,IJ,IL,VAL1,VAL2,VAL3,IK,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(ILND,IJ,JL,IK,KL,VAL3,VAL2,VAL1,JK,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JKND,IJ,IK,JL,KL,VAL3,VAL2,VAL1,IL,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
C
      IF(KLDC.NE.0) CALL RS1(KLDC,A1,VAL2,VAL1,VAL3,IJ,
     1                       rs,ex,eavf,lg,s)
C
 4200      CONTINUE
              GOTO 408
 128          CONTINUE
              WRITE(NOUT,1100) 'J=K    ','FIFTEENTH  '
              CALL MABORT
 138          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      IF(IKSYM.EQ.0.AND.ILSYM.EQ.0)
     .   CALL GA3(IK  ,IL  ,IL,KL,II,IK,VAL1,VAL2,IOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB3(IKND,IL,KL,II,IK,VAL1,VAL2,IOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB3(ILND,IK,KL,II,IL,VAL1,VAL2,IOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
C
      IF(KLDC.NE.0) CALL RS1(KLDC,A1S2,VAL2,VAL1,VAL1,IJ,
     1                       rs,ex,eavf,lg,s)
C
 4210      CONTINUE
              GOTO 408
 148          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      IF(IKSYM.EQ.0.AND.JKSYM.EQ.0)
     .   CALL GA4(IK  ,JK  ,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB4(IKND,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB4(JKND,JK,KK,IJ,IK,VAL1,VAL2,KOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
              GOTO 408
 158          CONTINUE
              WRITE(NOUT,1100) 'I=J=K  ','FIFTEENTH '
              CALL MABORT
 168          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      IF(IKSYM.EQ.0)
     .   CALL GA6(IK  ,IK,KK,II,VAL1,VAL2,IOR,KOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB6(IKND,II,KK,VAL1,VAL2,IK,IOR,KOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
              GOTO 408
 178          CONTINUE
              WRITE(NOUT,1100) 'I=J=K=L','FIFTEENTH '
              CALL MABORT
 408          CONTINUE
 4190       CONTINUE
 4150     CONTINUE
          DO 4240 KOR=NOP,JOR
            KO=IOFF(KOR)
            IK=IO+KOR
            JK=JO+KOR
            KK=KO+KOR
            KSM=ORBSYM(KOR)
            JKEQ=JOR.EQ.KOR
            IJKT=IJT
            IF (JKEQ) IJKT=IJKT+1
            IJKEQ=IJEQ.AND.JKEQ
            IJKSYM=IEOR(IJSYM,(KSM-1))
            LOFF=IJKSYM*NORBS
            KAD=JAD+KADD(KOFF+KOR)
C
C  ALP-ALP-ALP-UOC/DOC
C
            DO 4250 LOR=1,NDOC
         IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO 607
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO 607
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
              JULIA(IJKLT+21)=JULIA(IJKLT+21)+1
              GOTO (717,727,737,747,757,767,777)IJKLT
 717          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
      CALL GB1(KLND,IK,IL,JK,JL,VAL2,VAL1,VAL3,IJ,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JLND,JK,KL,IJ,IL,VAL1,VAL2,VAL3,IK,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(ILND,IJ,JL,IK,KL,VAL3,VAL2,VAL1,JK,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
              GOTO 607
 727          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
      CALL GBD(ILND,IJ,JL,VAL2,VAL1,JJ,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB21(JLND,JJ,JL,IJ,IL,VAL1,VAL2,JOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
              GOTO 607
 737          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
      CALL GBD(KLND,IK,IL,VAL2,VAL1,II,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB3(ILND,IK,KL,II,IL,VAL1,VAL2,IOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
              GOTO 607
 747          CONTINUE
              WRITE(NOUT,1100) 'K=L    ','FIFTEENTH+1'
              CALL MABORT
 757          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
      CALL GBD(ILND,II,IL,VAL1,VAL1,II,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      IF(LOR.GT.NUOC) CALL GBD(ILND,IL,LL,VAL2,VAL2,LL,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
CTJL  IF(LOR.GT.NUOC) CALL GBD(ILND,II,IL,VAL2,VAL2,LL)
              GOTO 607
 767          CONTINUE
              WRITE(NOUT,1100) 'I=J=K=L','FIFTEENTH+1'
              CALL MABORT
 777          CONTINUE
              WRITE(NOUT,1100) 'I=J=K=L','FIFTEENTH+1'
              CALL MABORT
 607          CONTINUE
 4250       CONTINUE
 4240     CONTINUE
 4080   CONTINUE
C
        DO 4300 JOR=NDOC,LFERM,-1
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
C
          IF (IGRP.NE.IJGRP(IJ)) THEN
            CALL SREAD(IT52,INTS,NMAX*2)
            SRIN=SRIN+NMAX*2/112
            IGRP=IJGRP(IJ)
          END IF
C
          DO 4310 KOR=1,NUOC
            KO=IOFF(KOR)
            IK=IO+KOR
            JK=JO+KOR
            KK=KO+KOR
            IKND=IJORD(IK)
            JKND=IJORD(JK)
            KSM=ORBSYM(KOR)
            IKSYM=IEOR((ISM-1),(KSM-1))
            JKSYM=IEOR((JSM-1),(KSM-1))
            JKEQ=JOR.EQ.KOR
            IJKT=IJT
            IF (JKEQ) IJKT=IJKT+1
            IJKEQ=IJEQ.AND.JKEQ
            IJKSYM=IEOR(IJSYM,(KSM-1))
            LOFF=IJKSYM*NORBS
            KAD=JAD+KADD(KOFF+KOR)
C
C  ALP-DOC-UOC-UOC
C
            DO 4320 LOR=1,KOR
         IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO 409
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO 409
              ILSYM=IEOR((ISM-1),(LSM-1))
              JLSYM=IEOR((JSM-1),(LSM-1))
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
              JULIA(IJKLT+28)=JULIA(IJKLT+28)+1
CTJL
         KLVU = IJDEP(KL)
CTJL
              GOTO (119,129,139,149,159,169,179)IJKLT
 119          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
CTJL
      IF(KLVU.NE.0) CALL EX1(IJND,LL,KK,LOR,KOR,VAL2,VAL1,VAL3,
     1                       rs,ex,eavf,lg,s)
C
      IF(IKSYM.EQ.0.AND.JLSYM.EQ.0)
     .   CALL GA1(IK  ,JL  ,IL,KL,IJ,JK,VAL1,VAL2,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
      IF(ILSYM.EQ.0.AND.JKSYM.EQ.0)
     .   CALL GA1(IL  ,JK  ,IK,KL,IJ,JL,VAL3,VAL2,VAL1,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB1(IKND,IL,KL,IJ,JK,VAL1,VAL2,VAL3,JL,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JLND,JK,KL,IJ,IL,VAL1,VAL2,VAL3,IK,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(ILND,IJ,JL,IK,KL,VAL3,VAL2,VAL1,JK,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JKND,IJ,IK,JL,KL,VAL3,VAL2,VAL1,IL,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
C
      IF(KLVU.NE.0) CALL RS1(KLVU,A3S2,VAL2,VAL1,VAL3,IJ,
     1                       rs,ex,eavf,lg,s)
C
 4330     CONTINUE
              GOTO 409
 129          CONTINUE
              WRITE(NOUT,1100) 'J=K    ','SIXTEENTH  '
              CALL MABORT
 139          CONTINUE
              WRITE(NOUT,1100) 'I=J    ','SIXTEENTH  '
              CALL MABORT
 149          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      IF(IKSYM.EQ.0.AND.JKSYM.EQ.0)
     .   CALL GA4(IK  ,JK  ,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB4(IKND,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB4(JKND,JK,KK,IJ,IK,VAL1,VAL2,KOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
              GOTO 409
 159          CONTINUE
              WRITE(NOUT,1100) 'I=J=K  ','SIXTEENTH  '
              CALL MABORT
 169          CONTINUE
              WRITE(NOUT,1100) 'I=J;K=L','SIXTEENTH  '
              CALL MABORT
 179          CONTINUE
              WRITE(NOUT,1100) 'I=J=K=L','SIXTEENTH  '
              CALL MABORT
 409          CONTINUE
 4320       CONTINUE
 4310     CONTINUE
          DO 4350 KOR=LFERM,JOR
            KO=IOFF(KOR)
            IK=IO+KOR
            JK=JO+KOR
            KK=KO+KOR
            IKND=IJORD(IK)
            JKND=IJORD(JK)
            KSM=ORBSYM(KOR)
            IKSYM=IEOR((ISM-1),(KSM-1))
            JKSYM=IEOR((JSM-1),(KSM-1))
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
C
C  ALP-DOC-DOC-UOC
C
            DO 4360 LOR=1,NUOC
         IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO 410
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO 410
              ILSYM=IEOR((ISM-1),(LSM-1))
              JLSYM=IEOR((JSM-1),(LSM-1))
              KLSYM=IEOR((KSM-1),(LSM-1))
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
              JULIA(IJKLT+35)=JULIA(IJKLT+35)+1
              GOTO (211,221,231,241,251,261,271)IJKLT
 211          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
C
      IF(JKDC.NE.0) CALL EX1(ILND,KK,JJ,KOR,JOR,VAL3,VAL1,VAL2,
     1                       rs,ex,eavf,lg,s)
C
      IF(IJSYM.EQ.0.AND.KLSYM.EQ.0)
     .   CALL GA1(IJ  ,KL  ,IL,JL,IK,JK,VAL2,VAL1,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
      IF(IKSYM.EQ.0.AND.JLSYM.EQ.0)
     .   CALL GA1(IK  ,JL  ,IL,KL,IJ,JK,VAL1,VAL2,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB1(IJND,IL,JL,IK,JK,VAL2,VAL1,VAL3,KL,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(IKND,IL,KL,IJ,JK,VAL1,VAL2,VAL3,JL,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JLND,JK,KL,IJ,IL,VAL1,VAL2,VAL3,IK,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(ILND,IJ,JL,IK,KL,VAL3,VAL2,VAL1,JK,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(KLND,IK,IL,JK,JL,VAL2,VAL1,VAL3,IJ,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
C
      IF(JKDC.NE.0) CALL RS1(JKDC,A1S2,VAL3,VAL1,VAL2,IL,
     1                       rs,ex,eavf,lg,s)
C
 4370      CONTINUE
              GOTO 410
 221          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      IF(IJSYM.EQ.0.AND.JLSYM.EQ.0)
     .   CALL GA2(IJ  ,JL  ,IL,JL,IJ,JJ,VAL1,VAL2,JOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GBD(ILND,IJ,JL,VAL2,VAL1,JJ,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB22(IJND,IL,JL,IJ,JJ,VAL1,VAL2,JOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB21(JLND,JJ,JL,IJ,IL,VAL1,VAL2,JOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
              GOTO 410
 231          CONTINUE
              WRITE(NOUT,1100) 'I=J    ','SEVENTEENTH'
              CALL MABORT
 241          CONTINUE
              WRITE(NOUT,1100) 'K=L    ','SEVENTEENTH'
              CALL MABORT
 251          CONTINUE
              WRITE(NOUT,1100) 'I=J=K  ','SEVENTEENTH'
              CALL MABORT
 261          CONTINUE
              WRITE(NOUT,1100) 'I=J;K=L','SEVENTEENTH'
              CALL MABORT
 271          CONTINUE
              WRITE(NOUT,1100) 'I=J=K=L','SEVENTEENTH'
              CALL MABORT
 410          CONTINUE
 4360       CONTINUE
C
C ALP-DOC-DOC-DOC
C
            DO 4390 LOR=LFERM,KOR
         IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO 710
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO 710
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
              JULIA(IJKLT+42)=JULIA(IJKLT+42)+1
CTJL
            JLDC = IJDEP(JL)
            KLDC = IJDEP(KL)
CTJL
              GOTO (711,721,731,741,751,761,771)IJKLT
 711          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
      CALL GB1(IJND,IL,JL,IK,JK,VAL2,VAL1,VAL3,KL,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(IKND,IL,KL,IJ,JK,VAL1,VAL2,VAL3,JL,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(ILND,IJ,JL,IK,KL,VAL3,VAL2,VAL1,JK,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
CTJL
      IF(JKDC.NE.0) CALL EX1(ILND,KK,JJ,KOR,JOR,VAL3,VAL1,VAL2,
     1                       rs,ex,eavf,lg,s)
      IF(JLDC.NE.0) CALL EX1(IKND,LL,JJ,LOR,JOR,VAL1,VAL3,VAL2,
     1                       rs,ex,eavf,lg,s)
      IF(KLDC.NE.0) CALL EX1(IJND,LL,KK,LOR,KOR,VAL2,VAL1,VAL3,
     1                       rs,ex,eavf,lg,s)
C
      IF(JKDC.NE.0) CALL RS1(JKDC,A3S2,VAL3,VAL1,VAL2,IL,
     1                       rs,ex,eavf,lg,s)
      IF(JLDC.NE.0) CALL RS1(JLDC,A3S2,VAL1,VAL3,VAL2,IK,
     1                       rs,ex,eavf,lg,s)
      IF(KLDC.NE.0) CALL RS1(KLDC,A3S2,VAL2,VAL1,VAL3,IJ,
     1                       rs,ex,eavf,lg,s)
C
 4400      CONTINUE
              GOTO 710
 721          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
      CALL GBD(ILND,IJ,JL,VAL2,VAL1,JJ,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB22(IJND,IL,JL,IJ,JJ,VAL1,VAL2,JOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
CTJL
      IF(JLDC.NE.0) CALL EX1(IKND,LL,JJ,LOR,JOR,VAL1,VAL1,VAL2,
     1                       rs,ex,eavf,lg,s)
C
      IF(JLDC.NE.0) CALL RS1(JLDC,A3S2,VAL1,VAL1,VAL2,IK,
     1                       rs,ex,eavf,lg,s)
 4410      CONTINUE
              GOTO 710
 731          CONTINUE
              WRITE(NOUT,1100) 'I=J    ','SEVENTNTH+1'
              CALL MABORT
 741          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
      CALL GBD(IJND,IK,JK,VAL2,VAL1,KK,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB4(IKND,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
C
      IF(JLDC.NE.0) CALL EX1(IKND,LL,JJ,LOR,JOR,VAL1,VAL1,VAL2,
     1                       rs,ex,eavf,lg,s)
C
      IF(JLDC.NE.0) CALL RS1(JLDC,A3S2,VAL1,VAL1,VAL2,IK,
     1                       rs,ex,eavf,lg,s)
C
 4420      CONTINUE
              GOTO 710
 751          CONTINUE
              WRITE(NOUT,1100) 'I=J=K  ','SEVENTNTH+1'
              CALL MABORT
 761          CONTINUE
              WRITE(NOUT,1100) 'I=J;K=L','SEVENTNTH+1'
              CALL MABORT
 771          CONTINUE
              WRITE(NOUT,1100) 'I=J=K=L','SEVENTNTH+1'
              CALL MABORT
 710          CONTINUE
 4390       CONTINUE
 4350     CONTINUE
 4300   CONTINUE
        DO 4430 JOR=NUOC,1,-1
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
C
          IF (IGRP.NE.IJGRP(IJ)) THEN
            CALL SREAD(IT52,INTS,NMAX*2)
            SRIN=SRIN+NMAX*2/112
            IGRP=IJGRP(IJ)
          END IF
C
       ISKIP3 = 0
       IF(ISKIP3.EQ.1) GOTO 8400
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
C     ALP-UOC-UOC-UOC
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
 8400       CONTINUE
CTJL
 4430   CONTINUE
 4070 CONTINUE
C
 540  CONTINUE
C
      DO 4440 IOR=NDOC,LFERM,-1
        IO=IOFF(IOR)
        II=IO+IOR
        ISM=ORBSYM(IOR)
        ITYP=1
        DO 4450 JOR=IOR,LFERM,-1
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
C
          IF (IGRP.NE.IJGRP(IJ)) THEN
            CALL SREAD(IT52,INTS,NMAX*2)
            SRIN=SRIN+NMAX*2/112
            IGRP=IJGRP(IJ)
          END IF
C
CTJL
        IJDC = IJDEP(IJ)
CTJL
          DO 4460 KOR=1,NUOC
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
            IKSYM=IEOR((ISM-1),(KSM-1))
            JKSYM=IEOR((JSM-1),(KSM-1))
            LOFF=IJKSYM*NORBS
            KAD=JAD+KADD(KOFF+KOR)
C
C  DOC-DOC-UOC-UOC
C
            DO 4470 LOR=1,KOR
         IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO 411
              LSM=ORBSYM(LOR)
              ILSYM=IEOR((ISM-1),(LSM-1))
              JLSYM=IEOR((JSM-1),(LSM-1))
              IF ((LSM-1) .NE. IJKSYM) GOTO 411
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
              JULIA(IJKLT+49)=JULIA(IJKLT+49)+1
CTJL
             KLVU = IJDEP(KL)
CTJL
              GOTO (212,222,232,242,252,262,272)IJKLT
 212          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
C
CGES
C     CONSTRUCT A MATRIX ONLY IF BOTH PAIRS ARE TOTALLY SYMMETRIC
C
CGES  CALL GA1(IKND,JLND,IL,KL,IJ,JK,VAL1,VAL2,VAL3)
CGES  CALL GA1(ILND,JKND,IK,KL,IJ,JL,VAL3,VAL2,VAL1)
C
      IF(IKSYM.EQ.0.AND.JLSYM.EQ.0)
     .   CALL GA1(IK  ,JL  ,IL,KL,IJ,JK,VAL1,VAL2,VAL3,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
      IF(ILSYM.EQ.0.AND.JKSYM.EQ.0)
     .   CALL GA1(IL  ,JK  ,IK,KL,IJ,JL,VAL3,VAL2,VAL1,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB1(IKND,IL,KL,IJ,JK,VAL1,VAL2,VAL3,JL,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JLND,JK,KL,IJ,IL,VAL1,VAL2,VAL3,IK,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(ILND,IJ,JL,IK,KL,VAL3,VAL2,VAL1,JK,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JKND,IJ,IK,JL,KL,VAL3,VAL2,VAL1,IL,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
C
      IF(KLVU.NE.0) CALL RS1(KLVU,A2,VAL2,VAL1,VAL3,IJ,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 411
 222          CONTINUE
              WRITE(NOUT,1100) 'J=K    ','EIGHTEENTH '
              CALL MABORT
 232          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      IF(IKSYM.EQ.0.AND.ILSYM.EQ.0)
     .   CALL GA3(IK  ,IL  ,IL,KL,II,IK,VAL1,VAL2,IOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB3(IKND,IL,KL,II,IK,VAL1,VAL2,IOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB3(ILND,IK,KL,II,IL,VAL1,VAL2,IOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
C
      IF(KLVU.NE.0) CALL RS1(KLVU,A1,VAL2,VAL1,VAL1,IJ,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 411
 242          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      IF(IKSYM.EQ.0.AND.JKSYM.EQ.0)
     .   CALL GA4(IK  ,JK  ,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB4(IKND,IK,KK,IJ,JK,VAL1,VAL2,KOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB4(JKND,JK,KK,IJ,IK,VAL1,VAL2,KOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
              GOTO 411
 252          CONTINUE
              WRITE(NOUT,1100) 'I=J=K  ','EIGHTEENTH '
              CALL MABORT
 262          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
      IF(IKSYM.EQ.0)
     .   CALL GA6(IK  ,IK,KK,II,VAL1,VAL2,IOR,KOR,
     2      aa,eta,zeta,alpha,beta,doff,dadd,idim,idima,ntypes+1)
C
      CALL GB6(IKND,II,KK,VAL1,VAL2,IK,IOR,KOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
              GOTO 411
 272          CONTINUE
              WRITE(NOUT,1100) 'I=J=K=L','EIGHTEENTH '
              CALL MABORT
 411          CONTINUE
 4470       CONTINUE
 4460     CONTINUE
          DO 4520 KOR=LFERM,JOR
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
            IKDC = IJDEP(IK)
            JKDC = IJDEP(JK)
CTJL
C
C  DOC-DOC-DOC-DOC
C
       ISKIP = 0
       IF(ISKIP.EQ.1) GOTO 9220
CTJL
            DO 9530 LOR=LFERM,KOR
            IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO 9811
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO 9811
              LO=IOFF(LOR)
              IL=IO+LOR
              JL=JO+LOR
              KL=KO+LOR
              LL=LO+LOR
              IJKLT=IJKT
              IF (IJKEQ .AND. KOR.GT.LOR) IJKLT=IJKLT+1
              IF (KOR.EQ.LOR) IJKLT=IJKLT+3
              LAD=KAD+LADD(LOFF+LOR)
              IF (IJKLT.GT.7 .OR. IJKLT.LT.1) GOTO 9903
              IF (LAD .GT. NMAX) GOTO 9902
              JULIA(IJKLT+56)=JULIA(IJKLT+56)+1
CTJL
          ILDC = IJDEP(IL)
          JLDC = IJDEP(JL)
          KLDC = IJDEP(KL)
CTJL
              GOTO (9712,9722,9732,9742,9752,9762,9772)IJKLT
 9712         VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
CTJL
C
      IF(IJDC.NE.0) CALL RS1(IJDC,A2,VAL2,VAL1,VAL3,KL,
     1                       rs,ex,eavf,lg,s)
      IF(IKDC.NE.0) CALL RS1(IKDC,A2,VAL1,VAL3,VAL2,JL,
     1                       rs,ex,eavf,lg,s)
      IF(ILDC.NE.0) CALL RS1(ILDC,A2,VAL3,VAL1,VAL2,JK,
     1                       rs,ex,eavf,lg,s)
      IF(JKDC.NE.0) CALL RS1(JKDC,A2,VAL3,VAL1,VAL2,IL,
     1                       rs,ex,eavf,lg,s)
      IF(JLDC.NE.0) CALL RS1(JLDC,A2,VAL1,VAL3,VAL2,IK,
     1                       rs,ex,eavf,lg,s)
      IF(KLDC.NE.0) CALL RS1(KLDC,A2,VAL2,VAL1,VAL3,IJ,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 9811
 9722         VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
CTJL
cets121890
      if(itcscf.eq.1) then
         e11(ior,lor)=e11(ior,lor)+a2*val2-val1
         e22(ior,lor)=e22(ior,lor)+a2*val2-val1
         e11(lor,ior)=e11(lor,ior)+a2*val2-val1
         e22(lor,ior)=e22(lor,ior)+a2*val2-val1
      end if
cets
C
      IF(IJDC.NE.0) CALL RS1(IJDC,A2,VAL1,VAL1,VAL2,KL,
     1                       rs,ex,eavf,lg,s)
      IF(ILDC.NE.0) CALL RS1(ILDC,A1,VAL2,VAL1,VAL1,JK,
     1                       rs,ex,eavf,lg,s)
      IF(JLDC.NE.0) CALL RS1(JLDC,A2,VAL1,VAL1,VAL2,IK,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 9811
 9732         VAL1 = INTS(LAD+1)
              VAL2 = INTS(LAD+2)
CTJL
cets121890
      if(itcscf.eq.1) then
         e11(kor,lor)=e11(kor,lor)+a2*val2-val1
         e22(kor,lor)=e22(kor,lor)+a2*val2-val1
         e11(lor,kor)=e11(lor,kor)+a2*val2-val1
         e22(lor,kor)=e22(lor,kor)+a2*val2-val1
      end if
cets
C
      IF(IKDC.NE.0) CALL RS1(IKDC,A2,VAL1,VAL1,VAL2,JL,
     1                       rs,ex,eavf,lg,s)
      IF(ILDC.NE.0) CALL RS1(ILDC,A2,VAL1,VAL1,VAL2,JK,
     1                       rs,ex,eavf,lg,s)
      IF(KLDC.NE.0) CALL RS1(KLDC,A1,VAL2,VAL1,VAL1,IJ,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 9811
 9742         VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
CTJL
cets121890
      if(itcscf.eq.1) then
         e11(ior,jor)=e11(ior,jor)+a2*val2-val1
         e22(ior,jor)=e22(ior,jor)+a2*val2-val1
         e11(jor,ior)=e11(jor,ior)+a2*val2-val1
         e22(jor,ior)=e22(jor,ior)+a2*val2-val1
      end if
cets
C
      IF(IJDC.NE.0) CALL RS1(IJDC,A1,VAL2,VAL1,VAL1,KL,
     1                       rs,ex,eavf,lg,s)
      IF(IKDC.NE.0) CALL RS1(IKDC,A2,VAL1,VAL1,VAL2,JL,
     1                       rs,ex,eavf,lg,s)
      IF(JKDC.NE.0) CALL RS1(JKDC,A2,VAL1,VAL1,VAL2,IL,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 9811
 9752         VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              val3=ints(lad+3)
CTJL
CTJL THIS TERM , EITHER I=J=K OR J=K=L, DOES NOT CONTRIBUTE.
CTJL
cets121890
      if(itcscf.eq.1) then
         e11(ior,lor)=e11(ior,lor)+val2+val1+val3
         e22(ior,lor)=e22(ior,lor)+val2+val1+val3
         e11(lor,ior)=e11(lor,ior)+val2+val1+val3
         e22(lor,ior)=e22(lor,ior)+val2+val1+val3
      end if
cets
C
cets  IF(ILDC.NE.0) CALL RS1(ILDC,A1,VAL2,VAL2,VAL2,LL,
cets 1                       rs,ex,eavf,lg,s)
      IF(IJDC.NE.0) CALL RS1(IJDC,A1,VAL2,VAL2,VAL2,KL,
     1                       rs,ex,eavf,lg,s)
      IF(ILDC.NE.0) CALL RS1(ILDC,A1,VAL1,VAL1,VAL1,JK,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 9811
 9762         VAL1 = INTS(LAD+1)
              VAL2 = INTS(LAD+2)
CTJL
              h11=h11+ a4*val2 - a2*val1
              h22=h22+ a4*val2 - a2*val1
cets121890
      if(itcscf.eq.1) then
         e11(ior,ior)=e11(ior,ior)+a2*val2-val1
         e22(ior,ior)=e22(ior,ior)+a2*val2-val1
         e11(lor,lor)=e11(lor,lor)+a2*val2-val1
         e22(lor,lor)=e22(lor,lor)+a2*val2-val1
      end if
cets
C
      IF(IKDC.NE.0) CALL RS1(IKDC,A2,VAL1,VAL1,VAL2,JL,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 9811
 9772         CONTINUE
cets121890
              val1=ints(lad+1)
              val2=ints(lad+2)
      if(itcscf.eq.1) then
         e11(ior,ior)=e11(ior,ior)+val2+val1
         e22(ior,ior)=e22(ior,ior)+val2+val1
         h11 = h11 +val1 + a2*val2
         h22 = h22 +val1 + a2*val2
         one = one +a2*val2
      endif
cets
CTJL          WRITE(NOUT,1100) 'I=J=K=L','NINTEENTH+1'
CTJL   I=J=J=L DOES NOT CONTRIBUTE
CTJL          CALL MABORT
 9811         CONTINUE
 9530       CONTINUE
 9220       CONTINUE
C
C  DOC-DOC-DOC-UOC
C
            DO 4530 LOR=1,NUOC
         IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO 811
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO 811
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
              JULIA(IJKLT+63)=JULIA(IJKLT+63)+1
              GOTO (712,722,732,742,752,762,772)IJKLT
 712          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
C
      IF(IJDC.NE.0) CALL EX1(KLND,JJ,II,JOR,IOR,VAL2,VAL1,VAL3,
     1                       rs,ex,eavf,lg,s)
      IF(IKDC.NE.0) CALL EX1(JLND,KK,II,KOR,IOR,VAL1,VAL3,VAL2,
     1                       rs,ex,eavf,lg,s)
      IF(JKDC.NE.0) CALL EX1(ILND,KK,JJ,KOR,JOR,VAL3,VAL1,VAL2,
     1                       rs,ex,eavf,lg,s)
C
      IF(IJDC.NE.0) CALL RS1(IJDC,A1,VAL2,VAL1,VAL3,KL,
     1                       rs,ex,eavf,lg,s)
      IF(IKDC.NE.0) CALL RS1(IKDC,A1,VAL1,VAL3,VAL2,JL,
     1                       rs,ex,eavf,lg,s)
      IF(JKDC.NE.0) CALL RS1(JKDC,A1,VAL3,VAL1,VAL2,IL,
     1                       rs,ex,eavf,lg,s)
C
      CALL GB1(KLND,IK,IL,JK,JL,VAL2,VAL1,VAL3,IJ,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(JLND,JK,KL,IJ,IL,VAL1,VAL2,VAL3,IK,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB1(ILND,IJ,JL,IK,KL,VAL3,VAL2,VAL1,JK,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
              GOTO 811
 722          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
cets121890
      if(itcscf.eq.1) then
         e11(lor,ior)=e11(lor,ior)+a2*val2-val1
         e22(lor,ior)=e22(lor,ior)+a2*val2-val1
      endif
c
      IF(IJDC.NE.0) CALL EX1(KLND,JJ,II,JOR,IOR,VAL1,VAL1,VAL2,
     1                       rs,ex,eavf,lg,s)
C
      IF(IJDC.NE.0) CALL RS1(IJDC,A1,VAL1,VAL2,VAL1,KL,
     1                       rs,ex,eavf,lg,s)
C
      CALL GBD(ILND,IJ,JL,VAL2,VAL1,JJ,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB21(JLND,JJ,JL,IJ,IL,VAL1,VAL2,JOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
              GOTO 811
 732          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
C
cets121890
      if(itcscf.eq.1) then
         e11(lor,kor)=e11(lor,kor)+a2*val2-val1
         e22(lor,kor)=e22(lor,kor)+a2*val2-val1
      endif
c
      IF(IKDC.NE.0) CALL EX1(JLND,KK,II,KOR,IOR,VAL1,VAL1,VAL2,
     1                       rs,ex,eavf,lg,s)
C
      IF(IKDC.NE.0) CALL RS1(IKDC,A1,VAL1,VAL1,VAL2,JL,
     1                       rs,ex,eavf,lg,s)
C
      CALL GBD(KLND,IK,IL,VAL2,VAL1,II,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
      CALL GB3(ILND,IK,KL,II,IL,VAL1,VAL2,IOR,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
              GOTO 811
 742          CONTINUE
              WRITE(NOUT,1100) 'K=L    ','NINTEENTH  '
              CALL MABORT
 752          VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
cets121890
      if(itcscf.eq.1) then
         val3=ints(lad+3)
         e11(lor,ior)=e11(lor,ior)+val1+val3
         e22(lor,ior)=e22(lor,ior)+val1+val3
      endif
c
      CALL GBD(ILND,II,IL,VAL1,VAL1,II,
     1         a,eta,zeta,alpha,beta,s,ntypes+1)
              GOTO 811
 762          CONTINUE
              WRITE(NOUT,1100) 'I=J;K=L','NINETEENTH '
              CALL MABORT
 772          CONTINUE
              WRITE(NOUT,1100) 'I=J=K=L','NINETEENTH '
              CALL MABORT
 811          CONTINUE
 4530       CONTINUE
 4520     CONTINUE
 4450   CONTINUE
        DO 4580 JOR=NUOC,1,-1
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
C
          IF (IGRP.NE.IJGRP(IJ)) THEN
            CALL SREAD(IT52,INTS,NMAX*2)
            SRIN=SRIN+NMAX*2/112
            IGRP=IJGRP(IJ)
          END IF
C
C     DOC-UOC-UOC-UOC
C
       ISKIP2 = 0
       IF(ISKIP2.EQ.1) GOTO 9400
          DO 9410 KOR=1,JOR
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
            DO 9420 LOR=1,KOR
          IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO  9430
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO 9430
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
          JLVU = IJDEP(JL)
          KLVU = IJDEP(KL)
CTJL
              GOTO (9412,9422,9432,9442,9452,9462,9472)IJKLT
 9412         VAL1=INTS(LAD+1)
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
      IF(JKVU.NE.0) CALL RS1(JKVU,A1,VAL3,VAL1,VAL2,IL,
     1                       rs,ex,eavf,lg,s)
      IF(JLVU.NE.0) CALL RS1(JLVU,A1,VAL1,VAL3,VAL2,IK,
     1                       rs,ex,eavf,lg,s)
      IF(KLVU.NE.0) CALL RS1(KLVU,A1,VAL2,VAL1,VAL3,IJ,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 9430
 9422         VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
      IF(JLVU.NE.0) CALL EX1(IKND,LL,JJ,LOR,JOR,VAL1,VAL1,VAL2,
     1                       rs,ex,eavf,lg,s)
C
      IF(JLVU.NE.0) CALL RS1(JLVU,A1,VAL1,VAL2,VAL1,IK,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 9430
 9432      WRITE(*,*) ' HIT I=J IN D-U-U-U ERROR IN LOOPS'
              GOTO 9430
 9442         VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
      IF(JKVU.NE.0) CALL EX1(ILND,KK,JJ,KOR,JOR,VAL1,VAL1,VAL2,
     1                       rs,ex,eavf,lg,s)
C
      IF(JKVU.NE.0) CALL RS1(JKVU,A1,VAL1,VAL2,VAL1,IL,
     1                       rs,ex,eavf,lg,s)
C
              GOTO 9430
 9452      WRITE(*,*) ' HIT CASE 5 IN D-U-U-U ERROR IN LOOPS'
CTJL
CTJL THIS TERM , EITHER I=J=K OR J=K=L, DOES NOT CONTRIBUTE.
CTJL
              GOTO 9430
 9462      WRITE(*,*) ' HIT CASE 6 IN D-U-U-U ERROR IN LOOPS'
CTJL
              GOTO 9430
 9472         CONTINUE
      WRITE(*,*) ' I=K=J=L ERROR IN LOOPS',IOR,JOR,KOR,LOR
CTJL          WRITE(NOUT,1100) 'I=J=K=L','NINTEENTH+1'
CTJL   I=J=J=L DOES NOT CONTRIBUTE
CTJL          CALL MABORT
 9430         CONTINUE
 9420       CONTINUE
 9410       CONTINUE
 9400       CONTINUE
CTJL
 4580   CONTINUE
 4440 CONTINUE
CC    JT4=104
CC    CALL SREW(JT4)
c     WRITE(*,*) ' THE A AND B MATRICES BEFORE ADDING EPA '
c     CALL MATOUT(aa,nma1,nma1+1,nma1,nma1+1,NOUT)
c     CALL MATOUT(A,IJIND,NDF,IJIND,NDF,NOUT)
      IEA=2*NDF+NBSET*NDF+1
      ISEA=JT1C(IEA)
      DO 4590 II=1,NDF
        CALL WREADW(IT94,EPA,INTOWP(NORBS*NORBS),ISEA,ISEA)
        IF (IAND(IPRINT,256) .NE. 0) THEN
          WRITE(NOUT,3014) II,ISEA
 3014     FORMAT(/'  EPA FOR DEGREE ',I5,'  AS READ FROM ADDRESS ',I8)
          CALL MATOUT(EPA,NORBS,NORBS,NORBS,NORBS,NOUT)
        END IF
        IJ=0
        DO 4600 I=1,NORBS
          DO 4610 J=1,I
            IJ=IJ+1
            INDEX=IJORD(IJ)
            IF(INDEX.NE.0)THEN
CGES     A(INDEX,II+IJIND2)=A(INDEX,II+IJIND2)-A2*(EPA(I,J)-EPA(J,I))
         A(INDEX,II       )=A(INDEX,II       )-A2*(EPA(I,J)-EPA(J,I))
            END IF
 4610     CONTINUE
 4600   CONTINUE
 4590 CONTINUE
CTJL
CTJLNEW
C       WRITE(*,*) ' THE AVERAGED FOCK OPERATOR FOCK MATRIX '
C       CALL PRINT(EAVF,NUMIJ,NORBS,NOUT)
        DO 9750 IDF = 1,NDF
        DO 9760 I = 1,NBF
         IO = IOFF(I)
         II = IO + I
         DO 9770 J = 1,I
         JO = IOFF(J)
          JJ = JO + J
          IJ = IO + J
          IJDC = IJDEP(IJ)
          IF(IJDC.EQ.0) GOTO 9770
C         DO 9780 M = 1,NUOC
C         MI = IO + M
C         MJ = JO + M
C         RS(IJDC,IDF) = RS(IJDC,IDF) + (S(MI,IDF)*EAVF(MJ) +
C    1  S(MJ,IDF)*EAVF(MI))
C         WRITE(*,*) ' IN VIRTUAL PART; I,J,M,MI,MJ,F1,F2 '
C         WRITE(*,*) I,J,M,MI,MJ,EAVF(MI),EAVF(MJ)
 9780     CONTINUE
          IF(NOP.EQ.0) GOTO 9771
          DO 9790 M = NOP,NBF
          MO = IOFF(M)
          MI = MO + I
          MJ = MO + J
          RS(IJDC,IDF) = RS(IJDC,IDF) + (S(MI,IDF)*EAVF(MJ) +
     1  S(MJ,IDF)*EAVF(MI))
 9790     CONTINUE
 9771     CONTINUE
          IF(NSPE.EQ.0) GOTO 9770
          DO 9791 M = NSA,NBF
          MO = IOFF(M)
          MI = MO + I
          MJ = MO + J
          RS(IJDC,IDF) = RS(IJDC,IDF) + (S(MI,IDF)*EAVF(MJ) +
     1  S(MJ,IDF)*EAVF(MI))
 9791     CONTINUE
 9770     CONTINUE
 9760     CONTINUE
 9750     CONTINUE
CTJLOLD
C     WRITE(*,*) ' BEFORE END ',IJIND
      ISEA = 1
      DO 5590 IDF=1,NDF
        CALL WREADW(IT96,LSTEPA,INTOWP(NUMIJ),ISEA,ISEA)
CTJL      WRITE(*,*) ' IT96,ISEA ',IT96,ISEA
CTJL    CALL PRINT(LSTEPA,NUMIJ,NORBS,NOUT)
        DO 9600 I = 1,NBF
         IO = IOFF(I)
         II = IO + I
         DO 9610 J = 1,I
         JO = IOFF(J)
          JJ = JO + J
          IJ = IO + J
          IJDC = IJDEP(IJ)
          IF(IJDC.EQ.0) GOTO 9610
          EFAC = A1/(EAVF(JJ)-EAVF(II))
          SFAC = (EAVF(II)+EAVF(JJ))*(-EFAC)
          XFAC = (LG(J,I)-LG(I,J))
          BOX(IDF)=BOX(IDF)+XFAC*((A2*LSTEPA(IJ)-RS(IJDC,IDF))*EFAC +
     1 SFAC*S(IJ,IDF))
CTJL      WRITE(*,*) ' IT96,ISEA ',IT96,ISEA
 9650     CONTINUE
 9610     CONTINUE
 9600     CONTINUE
 5590 CONTINUE
CTJL
CTJL
      CALL ZERO(EX2,IJIND)
      I = 0
      IF(I.EQ.1) RETURN
        DO 9660 I = 1,NBF
         IO = IOFF(I)
         II = IO + I
         DO 9670 J = 1,I
         JO = IOFF(J)
          JJ = JO + J
          IJ = IO + J
          IJDC = IJDEP(IJ)
          IF(IJDC.EQ.0) GOTO 9670
          EFAC = A1/(EAVF(JJ)-EAVF(II))
          XFAC = (LG(J,I)-LG(I,J))*EFAC
CTJL      WRITE(*,*) ' IN ALPHA PART; EFAC,XFAC ',EFAC,XFAC
cets121890
          if(itcscf.eq.1) then
            DO 9680 M = 1,NUOC
            MI = IO + M
            MJ = JO + M
            MJND = IJORD(MJ)
            MIND = IJORD(MI)
            EX2(MJND) = EX2(MJND) - XFAC*EAVF(MI)
            EX2(MIND) = EX2(MIND) - XFAC*EAVF(MJ)
CTJL      WRITE(*,*) ' IN VIRTUAL PART; M,EX2 ',M,EX2(MJND),EX2(MIND)
 9680     CONTINUE
          endif
cets
CTJL      WRITE(*,*) ' IN EX2 PART ; NOP,NBF ',NOP,NBF
          IF(NOP.EQ.0) GOTO 9671
          DO 9690 M = NOP,NBF
          MO = IOFF(M)
          MI = MO + I
          MJ = MO + J
          MJND = IJORD(MJ)
          MIND = IJORD(MI)
          EX2(MJND) = EX2(MJND) + XFAC*EAVF(MI)
          EX2(MIND) = EX2(MIND) + XFAC*EAVF(MJ)
CTJL      WRITE(*,*) ' IN ALPHA PART; M,EX2 ',M,EX2(MJND),EX2(MIND)
 9690     CONTINUE
 9671     CONTINUE
          IF(NSPE.EQ.0) GOTO 9670
          DO 9691 M = NSA,NBF
          MO = IOFF(M)
          MI = MO + I
          MJ = MO + J
          MJND = IJORD(MJ)
          MIND = IJORD(MI)
          EX2(MJND) = EX2(MJND) + XFAC*EAVF(MI)
          EX2(MIND) = EX2(MIND) + XFAC*EAVF(MJ)
CTJL      WRITE(*,*) ' IN ALPHA PART; M,EX2 ',M,EX2(MJND),EX2(MIND)
 9691     CONTINUE
 9670     CONTINUE
 9660     CONTINUE
C
C  EXTRA STUFF FOR TCSCF-CI GRADIENTS
C
      IF(ITCSCF.ne.0) THEN
      ETCSCF = CI1*CI1*H11 + CI2*CI2*H22 + A2*CI1*CI2*H12
      ETCSC1 = ETCSCF
      ETCSCF = ETCSCF + REPNUC
      WRITE(*,*) ' THE TCSCF ENERGY ',ETCSCF
c     write(6,*) 'ci1,ci2 again',ci1,ci2
c     write(6,*) 'h11,h12,h22',h11,h12,h22
c     write(6,*) ' THE e11,e12,e22 matrices'
c     call matout(e11,nbf,nbf,nbf,nbf,6)
c     write(6,*)
c     call matout(e12,nbf,nbf,nbf,nbf,6)
c     write(6,*)
c     call matout(e22,nbf,nbf,nbf,nbf,6)
c     write(6,*)
      EH12 = A0
      EH11 = A0
      EH22 = A0
      DO 9379 I = 1,NBF
C
      IF(I.EQ.NSA.OR.I.EQ.NSB) EH12 = EH12 + A2*CI1*CI2*E12(I,I)
      IF(I.LT.NSA) EH11 = EH11 + CI1*CI1*E11(I,I)
      IF(I.LT.NSA) EH22 = EH22 + CI2*CI2*E22(I,I)
C
      ESCF2=ESCF2+CI1*CI1*E11(I,I)+CI2*CI2*E22(I,I)+A2*CI1*CI2*E12(I,I)
 9379 CONTINUE
c
      IJ = 0
      ij1=0
      DO 19 I = 1,NBF
        ism=orbsym(i)
        DO 29 J = 1,I
          jsm=orbsym(j)
          ijsym=ieor((ism-1),(jsm-1))
          IJ = IJ + 1
          ijnd=ijord(ij)
          IF(IJND.EQ.0.or.ijsym.ne.0) GOTO 29
          ij1=ij1+1
          aa(nma1-1,ij1) =-A4*(CI1*(E11(I,J) - E11(J,I))
     1     + CI2*(E12(I,J) - E12(J,I)))
          aa(nma1,ij1) =-A4*(CI2*(E22(I,J) - E22(J,I))
     1     + CI1*(E12(I,J) - E12(J,I)))
          aa(ij1,nma1-1)=aa(nma1-1,ij1)
          aa(ij1,nma1)=aa(nma1,ij1)
   29   CONTINUE
   19 CONTINUE
c
      aa(nma1,nma1) = -A2*(H22 - ETCSC1 + CI2*CI2)
      aa(nma1-1,nma1-1) = -A2*(H11 - ETCSC1 + CI1*CI1)
      aa(nma1,nma1-1) = -A2*(H12 + CI1*CI2)
      aa(nma1-1,nma1) = -A2*(H12 + CI1*CI2)
C
C NOW FORM THE BAI TERMS.
C
      DO 18 IDF = 1,NDF
C
      SAIJ1 = A0
      SAIJ2 = A0
      ij=0
      ij1=0
      DO 68 I = 1,nbf
        ism=orbsym(i)
        do 67 j=1,i
          jsm=orbsym(j)
          ijsym=ieor((ism-1),(jsm-1))
          ij=ij+1
          IJND = IJORD(IJ)
          IF(IJND.EQ.0.or.ijsym.ne.0) GOTO 67
          ij1=ij1+1
          SAIJ1 = SAIJ1 + S(IJ,IDF)*aa(nma1-1,ij1)
          SAIJ2 = SAIJ2 + S(IJ,IDF)*aa(nma1,ij1)
  67    continue
   68 CONTINUE
      EA = A0
      DO 28 I = LFERM,NBF
      IO = IOFF(I)
      II = IO + I
      DO 38 J = LFERM,I-1
      IJ = IO + J
      EA = EA + ETA(IJ)*S(IJ,IDF)*A2
   38 CONTINUE
      EA = EA + ETA(II)*S(II,IDF)
   28 CONTINUE
C
      BI1 = A1S2*CI1*(CI1*CI1*HA11(IDF) + A2*CI1*CI2*HA12(IDF)
     1      + CI2*CI2*HA22(IDF) - A2*EA)
      BI2 = A1S2*CI2*(CI1*CI1*HA11(IDF) + A2*CI1*CI2*HA12(IDF)
     1      + CI2*CI2*HA22(IDF) - A2*EA)
      ETEST = (A2*BI1)/CI1
      ETEST = (A2*BI2)/CI2
      BI1 = BI1 - A1S2*(CI1*HA11(IDF) + CI2*HA12(IDF))
      BI2 = BI2 - A1S2*(CI1*HA12(IDF) + CI2*HA22(IDF))
      DO 48 I = 1,NBF
      IO = IOFF(I)
      II = IO + I
      DO 58 J = 1,I-1
      IJ = IO + J
      BI1 = BI1 + A2*(CI1*E11(J,I) + CI2*E12(J,I))*S(IJ,IDF)
      BI2 = BI2 + A2*(CI1*E12(J,I) + CI2*E22(J,I))*S(IJ,IDF)
   58 CONTINUE
      BI1 = BI1 + (CI1*E11(I,I) + CI2*E12(I,I))*S(II,IDF)
      BI2 = BI2 + (CI1*E12(I,I) + CI2*E22(I,I))*S(II,IDF)
   48 CONTINUE
cets121890      A(IJIND-1,IJIND+IDF) = -(BI1*A4 - SAIJ1)
cets121890      A(IJIND,IJIND+IDF) = -(BI2*A4 - SAIJ2)
      A(IJIND-1,IDF) = -(BI1*A4 - SAIJ1)
      A(IJIND,IDF) = -(BI2*A4 - SAIJ2)
   18 CONTINUE
      END IF

C     WRITE(*,*) ' THE A AND B MATRICES AS CONSTRUCTED IN GETINT '
C     CALL MATOUT(A,IJIND2,IJIND3,IJIND2,IJIND3,NOUT)
CTJL  CALL MATOUT(B0,IJIND,NDF,IJIND,NDF,NOUT)
      idl52=1
      errcod = frdboo('DELETE52',idl52)
      if(idl52.eq.1) then
        call rclose(IT52,4)
      else
        call rclose(IT52,3)
      endif
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
