C**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
C**********************************************************************
C*  BY:  RICHARD REMINGTON                         search:  c5-22-88  *
C*  DATE:  MAY   22,   1988                                           *
C*  REASON: DECREASE CORE TO RUN IN 7MB ON 9370                       *
C**********************************************************************
C     PROGRAM LAST UPDATED JULY  07, 1987                             *
C     BY: GES                                                         *
C     REASON: COMMENT LINES ZEROING A MATRIX FOR CPHFL (LOOK FOR CGES)*
C     THIS SUBROUTINE CONSTRUCTS ETA,ZETA,BETA, ETC. AND WAS ONLY     *
C     ZEROING OUT A.
C**********************************************************************
C     PROGRAM LAST UPDATED JULY  12, 1985                             *
C     BY: TJL                                                         *
C**********************************************************************
C
      SUBROUTINE GETINT(ETA,ZETA,F,ALPHA,BETA,KADD,LADD,
     1IJADD,IJGRP,ORBSYM,INTS,IOFF,EAVF,CCSD)
      IMPLICIT INTEGER (A-Z)
CIBM  EXTENDED DUMMY A,ETA,ZETA,F,ALPHA,BETA,INTS,BUF,EIG,S,EPA,W,
CIBM 1E11,E22,E12
      LOGICAL IJEQ,IJKEQ,JKEQ,CCSD
      INTEGER N112(112),JULIA(70)
      REAL*8 A0,A1,A2,A1S2,A4,FOCC,EAVF(NUMIJ)
      DIMENSION KADD(SYMORB),LADD(SYMORB)
      DIMENSION ORBSYM(NORBS),IJADD(NUMIJ),IJGRP(NUMIJ)
      DIMENSION IOFF(NUMIJ)
      REAL*8 F(NORBS),INTS(NMAX)
cets121890      REAL*8 ZETA(NUMIJ,NORBS),ALPHA(NUMIJ),BETA(NUMIJ)
      REAL*8 zeta(ntypes+1,numij),ALPHA(NUMIJ),BETA(NUMIJ)
      REAL*8 C1A,C2A,C1,C2
      REAL*8 VAL1,VAL2,VAL3,ETA(NUMIJ),ESCF,REPNUC
      COMMON /CIA/ C1A(45),C2A(45),C1,C2,ESCF
C!!!  COMMON /TAPES/ IT22,IT8,NOUT,IT2,IT42,JTAP1,IUTAP,JTAP2
      COMMON /TAPES/ NOUT,ICHK,IT30,IT42,IT47,IT52,IT58,IT94,IT96,IT98
      COMMON /DIMSG/ SYMORB,NUMIJ,NBF,NMAX,NBFORB,NSYM,NDF,NORBS
      COMMON /POINT/ FOCC(5),IPOINT(5),JPOINT(5),NSORB(5)
      common /ptrs/ motyp(500),mopos(5)
      COMMON /FUNCS/ IOCC,JOCC,KOCC,NTYPES,NBSET,NAT
      COMMON /OCCS/ NSPE,NSA,NSB,NALP,NOP,NDOC,LFERM,NUOC,IJIND,IJIND3
      COMMON /LOC22/ IA22(100),IBF2AT(284)
      COMMON /IEXT/ IPRINT,IGRP,ICIDIP
      COMMON /REN/ REPNUC
C
      WRITE(ICHK,2222)
 2222 FORMAT(/'  IN GETINT')
 3001 FORMAT(10I5)
      DO 3010 IJUL=1,70
        JULIA(IJUL)=0
 3010 CONTINUE
      WRITE(ICHK,*) '   THE SPECIAL PARAMETERS '
      WRITE(ICHK,3001) NSPE,NSA,NSB,NALP,NOP,NDOC,
     +                 LFERM,NUOC,IJIND,IJIND3
CTJL  WRITE(ICHK,*) '   THE IJORD ARRAY:'
CTJL  WRITE(ICHK,3001) (IJORD(III), III=1,NUMIJ)
      IJIND2=IJIND
      A0=0.0D+00
      A1=1.0D+00
      A2=2.0D+00
      A4=4.00D+00
      A1S2=0.50D+00
      DO 3040 I=1,NORBS
        F(I)=F(I)*A1S2
 3040 CONTINUE
      ESCF=A0
C
C  NEW WAY OF OBTAINING INTSRT
      CALL RFILE(IT52)
      WRITE(ICHK,2223)
 2223 FORMAT(/'  OPENED IT52')
CTJL  CALL SREW(IT52)   THIS IS DONE IN RFILE
      IF(.NOT.CCSD)CALL SREAD(IT52,N112,112)
      CALL RGETSA(IT52,INTSRT)
CTJL  WRITE(ICHK,2224) INTSRT
C2224 FORMAT(/'  CALLED RGETSA',I8)
CIBM  CALL RSETSA(IT52,INTSRT)
CTJL  WRITE(ICHK,2225)  INTSRT
CTJL  WRITE(NOUT,2225) INTSRT
C2225 FORMAT(/'  CALLED RSETSA',I8)
      IGRP=0
C
      IIJJ=0
      DO 3080 II=1,NORBS
        DO 3090 JJ=1,II
          IIJJ=IIJJ+1
          ETA(IIJJ)=A0
          EAVF(IIJJ)=A0
          DO 3100 III=1,ntypes+1
            ZETA(iii,iijj)=A0
 3100     CONTINUE
 3090   CONTINUE
 3080 CONTINUE
c
c construct motyp and mopos arrays
      nclosd=nsorb(1)
      nopen1=0
      nopen2=0
      mopos(ntypes+1)=1
      if(ntypes.eq.1) then
         nvir=nsorb(2)
         mopos(1)=nvir+1
      else if(ntypes.eq.2) then
         nopen1=nsorb(2)
         nvir=nsorb(3)
         mopos(1)=nvir+1
         mopos(2)=nvir+nclosd+1
      else if(ntypes.eq.3) then
         nopen1=1
         nopen2=1
         nvir=nsorb(4)
         mopos(1)=nvir+1
         mopos(2)=nvir+nclosd+1
         mopos(3)=nvir+nclosd+2
      else
         write(6,*) 'ntypes = ',ntypes
         write(6,*) 'what are you doing???'
         call mabort
      endif

      do 4000 i=1,nvir
         motyp(i)=ntypes+1
 4000 continue
      if(ntypes.eq.1) then
         do 4010 i=nvir+1,norbs
            motyp(i)=1
 4010    continue
      else if(ntypes.eq.2) then
         do 4020 i=nvir+1,nvir+nclosd
            motyp(i)=1
 4020    continue
         do 4021 i=nvir+nclosd+1,norbs
            motyp(i)=2
 4021    continue
      else
         do 4030 i=nvir+1,nvir+nclosd
            motyp(i)=1
 4030    continue
         motyp(norbs-1)=2
         motyp(norbs)=3
      endif
C
C
C  CONSTRUCT ETA AND ZETA MATRICIES
C
      DO 3110 IOR=NORBS,LFERM,-1
        IO=IOFF(IOR)
        II=IO+IOR
        ISM=ORBSYM(IOR)
        ITYP=1
        itype=motyp(ior)
        ipos=mopos(itype)
        DO 3120 JOR=IOR,LFERM,-1
          JO=IOFF(JOR)
          IJ=IO+JOR
          JJ=JO+JOR
          jtype=motyp(jor)
          jpos=mopos(jtype)
          IF (IGRP .NE. IJGRP(IOFF(IOR)+JOR)) THEN
            WRITE(ICHK,3002) IOR,JOR,IGRP,IJGRP(IOFF(IOR)+JOR)
 3002       FORMAT('  IOR=',I5,' JOR=',I5,' IGRP=',I5,' NEW GROUP=',I5)
CIBM        CALL RGETSA(IT52,IX)
CIBM        WRITE(NOUT,2001) IX
C2001       FORMAT('  READING FROM IT52 AT ADR ',I8)
            CALL SREAD(IT52,INTS,NMAX*2)
CCC            WRITE(NOUT,1901)
CCC            WRITE(NOUT,1902) (INTS(III), III=1,NMAX)
 1901       FORMAT(/'  HERE ARE A FEW INTEGRALS FROM IT52:')
 1902       FORMAT(10F10.6)
            SRIN=SRIN+NMAX*2/112
            IGRP=IJGRP(IOFF(IOR)+JOR)
            WRITE(ICHK,3003) IGRP
 3003       FORMAT('  NOW IGRP=',I5)
          END IF
          JSM=ORBSYM(JOR)
          IJEQ=IOR.EQ.JOR
          IJT=ITYP
          IF (IJEQ) IJT=IJT+2
          IJSYM=IEOR((ISM-1),(JSM-1))
          KOFF=IJSYM*NORBS
          JAD=IJADD(IJ)
          DO 3130 KOR=1,NUOC
            KO=IOFF(KOR)
            IK=IO+KOR
            JK=JO+KOR
            KK=KO+KOR
            KSM=ORBSYM(KOR)
            ktype=motyp(kor)
            kpos=mopos(ktype)
            JKEQ=JOR.EQ.KOR
            IJKT=IJT
            IF (JKEQ) IJKT=IJKT+1
            IJKEQ=IJEQ.AND.JKEQ
            IJKSYM=IEOR(IJSYM,(KSM-1))
            LOFF=IJKSYM*NORBS
            KAD=JAD+KADD(KOFF+KOR)
            DO 3140 LOR=1,KOR
              LO=IOFF(LOR)
              IL=IO+LOR
              JL=JO+LOR
              KL=KO+LOR
              LL=LO+LOR
              ltype=motyp(lor)
              lpos=mopos(ltype)
CTJL          IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO
CTJL *        300
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO 300
              IJKLT=IJKT
              IF (IJKEQ .AND. KOR.GT.LOR) IJKLT=IJKLT+1
              IF (KOR.EQ.LOR) IJKLT=IJKLT+3
              LAD=KAD+LADD(LOFF+LOR)
              IF (IJKLT.GT.7 .OR. IJKLT.LT.1) THEN
                WRITE(NOUT,2004) IJKLT,IOR,JOR,KOR,LOR,IBF,JBF,KBF,
     *          LBF,ISM,JSM,KSM,LSM
 2004           FORMAT(/'  IN GETDM, IJKLT UNKNOWN TYPE ',I5/'  ORBITALS
     * :',4I5/'  BASIS FUNCTIONS :',4I5/'   SYMMETRIES :',4I5)
                CALL ABORTS
              END IF
              IF (LAD .GT. NMAX) THEN
                WRITE(NOUT,2005) LAD
 2005           FORMAT(/'  IN GETDM, LAD TOO BIG :',I6)
                CALL ABORTS
              END IF
              GOTO (10,20,30,40,50,60,70)IJKLT
 10           CONTINUE
C  THIS TYPE (ALL DIFFERENT) DOES NOT CONTRIBUTE
              GOTO 300
 20           CONTINUE
              WRITE(NOUT,2006)
 2006         FORMAT(/     '  IN CONSTRUCTON OF ETA AND ZETA, J=K IS NOT
     * POSSIBLE DO THE FIRST SET OF LOOPS'/'  ERROR IN LOOPS SOMEWHERE')
              CALL ABORTS
 30           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              ETA(KL)=ETA(KL)+ALPHA(IK)*VAL2+BETA(IK)*VAL1
              EAVF(KL)=EAVF(KL)+F(IOR)*(A2*VAL2-VAL1)
cets121890              DO 3150 J2=1,NORBS
cets121890                IJ2=IOFF(MAX(IOR,J2))+MIN(IOR,J2)
cets121890                ZETA(KL,J2)=ZETA(KL,J2)+ALPHA(IJ2)*VAL2+BETA(IJ2)*
cets121890     *          VAL1
              DO 3150 J2=1,ntypes
                IJ2=IOFF(MAX(ipos,mopos(j2)))+MIN(ipos,mopos(j2))
                ZETA(j2,kl)=ZETA(j2,kl)+ALPHA(IJ2)*VAL2+BETA(IJ2)*
     *          VAL1
 3150         CONTINUE
              GOTO 300
 40           CONTINUE
C  THIS TYPE (K=L) DOES NOT CONTRIBUTE
              GOTO 300
 50           CONTINUE
              WRITE(NOUT,2007)
 2007         FORMAT(/     '  IN CONSTRUCTON OF ETA AND ZETA, I=J=K IS N
     *OT POSSIBLE DO THE FIRST SET OF LOOPS'/'  ERROR IN LOOPS SOMEWHERE
     *')
              CALL ABORTS
 60           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              ETA(LL)=ETA(LL)+ALPHA(IL)*VAL2+BETA(IL)*VAL1
              EAVF(LL)=EAVF(LL)+F(IOR)*(A2*VAL2-VAL1)
cets121890              DO 3160 J2=1,NORBS
cets121890                IJ2=IOFF(MAX(IOR,J2))+MIN(IOR,J2)
cets121890                ZETA(LL,J2)=ZETA(LL,J2)+ALPHA(IJ2)*VAL2+BETA(IJ2)*
cets121890     *          VAL1
              DO 3160 J2=1,ntypes
                IJ2=IOFF(MAX(ipos,mopos(j2)))+MIN(ipos,mopos(j2))
                ZETA(j2,ll)=ZETA(j2,ll)+ALPHA(IJ2)*VAL2+BETA(IJ2)*
     *          VAL1
 3160         CONTINUE
              GOTO 300
 70           CONTINUE
              WRITE(NOUT,2008)
 2008         FORMAT(/     '  IN CONSTRUCTON OF ETA AND ZETA,I=J=K=L IS
     *NOT POSSIBLE DO THE FIRST SET OF LOOPS'/'  ERROR IN LOOPS SOMEWHER
     *E')
              CALL ABORTS
 300          CONTINUE
 3140       CONTINUE
 3130     CONTINUE
          DO 3170 KOR=LFERM,JOR
            KSM=ORBSYM(KOR)
            KO=IOFF(KOR)
            IK=IO+KOR
            JK=JO+KOR
            KK=KO+KOR
            ktype=motyp(kor)
            kpos=mopos(ktype)
            JKEQ=JOR.EQ.KOR
            IJKT=IJT
            IF (JKEQ) IJKT=IJKT+1
            IJKEQ=IJEQ.AND.JKEQ
            IJKSYM=IEOR(IJSYM,(KSM-1))
            LOFF=IJKSYM*NORBS
            KAD=JAD+KADD(KOFF+KOR)
            DO 3180 LOR=1,NUOC
CTJL          IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO
CTJL *        301
              LSM=ORBSYM(LOR)
              ltype=motyp(lor)
              lpos=mopos(ltype) 
              IF ((LSM-1) .NE. IJKSYM) GOTO 301
              LO=IOFF(LOR)
              IL=IO+LOR
              JL=JO+LOR
              KL=KO+LOR
              LL=LO+LOR
              IJKLT=IJKT
              IF (IJKEQ .AND. KOR.GT.LOR) IJKLT=IJKLT+1
              IF (KOR.EQ.LOR) IJKLT=IJKLT+3
              LAD=KAD+LADD(LOFF+LOR)
              IF (IJKLT.GT.7 .OR. IJKLT.LT.1) THEN
                WRITE(NOUT,2004)IJKLT,IOR,JOR,
     *          KOR,LOR,IBF,JBF,KBF,LBF,ISM,JSM,KSM,LSM
                CALL ABORTS
              END IF
              IF (LAD .GT. NMAX) THEN
                WRITE(NOUT,2005) LAD
                CALL ABORTS
              END IF
              GOTO (11,21,31,41,51,61,71)IJKLT
 11           CONTINUE
C  THIS TYPE (ALL DIFFERENT) DOES NOT CONTRIBUTE
              GOTO 301
 21           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              ETA(IL)=ETA(IL)+ALPHA(IJ)*VAL2+BETA(IJ)*VAL1
              EAVF(IL)=EAVF(IL)+F(JOR)*(A2*VAL2-VAL1)
cets121890              DO 3190 K2=1,NORBS
cets121890                KL2=IOFF(MAX(JOR,K2))+MIN(JOR,K2)
cets121890                ZETA(IL,K2)=ZETA(IL,K2)+ALPHA(KL2)*VAL2+BETA(KL2)*
cets121890     *          VAL1
              DO 3190 K2=1,ntypes
                KL2=IOFF(MAX(jpos,mopos(k2)))+MIN(jpos,mopos(k2))
                ZETA(k2,il)=ZETA(k2,il)+ALPHA(KL2)*VAL2+BETA(KL2)*
     *          VAL1
 3190         CONTINUE
              GOTO 301
 31           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              ETA(KL)=ETA(KL)+ALPHA(IK)*VAL2+BETA(IK)*VAL1
              EAVF(KL)=EAVF(KL)+F(IOR)*(A2*VAL2-VAL1)
cets121890              DO 3200 J2=1,NORBS
cets121890                IJ2=IOFF(MAX(IOR,J2))+MIN(IOR,J2)
cets121890                ZETA(KL,J2)=ZETA(KL,J2)+ALPHA(IJ2)*VAL2+BETA(IJ2)*
cets121890     *          VAL1
              DO 3200 J2=1,ntypes
                IJ2=IOFF(MAX(ipos,mopos(j2)))+MIN(ipos,mopos(j2))
                ZETA(j2,kl)=ZETA(j2,kl)+ALPHA(IJ2)*VAL2+BETA(IJ2)*
     *          VAL1
 3200         CONTINUE
              GOTO 301
 41           CONTINUE
              WRITE(NOUT,2009)
 2009         FORMAT(/     '  IN CONSTRUCTON OF ETA AND ZETA,  K=L  IS N
     *OT POSSIBLE DO THE 2ND SET OF LOOPS'/'  ERROR IN LOOPS SOMEWHERE')
              CALL ABORTS
              GOTO 301
 51           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
              ETA(IL)=ETA(IL)+ALPHA(II)*VAL1+BETA(II)*VAL1+F(IOR)*
     *        VAL3
              EAVF(IL)=EAVF(IL)+F(IOR)*(A2*VAL1-VAL1) + VAL3
cets121890              DO 3210 J2=1,NORBS
cets121890                IJ2=IOFF(MAX(IOR,J2))+MIN(IOR,J2)
cets121890                ZETA(IL,J2)=ZETA(IL,J2)+ALPHA(IJ2)*VAL1+BETA(IJ2)*
cets121890     *          VAL1+F(J2)*VAL3
              DO 3210 J2=1,ntypes
                IJ2=IOFF(MAX(ipos,mopos(j2)))+MIN(ipos,mopos(j2))
                ZETA(j2,il)=ZETA(j2,il)+ALPHA(IJ2)*VAL1+BETA(IJ2)*
     *          VAL1+F(mopos(j2))*VAL3
 3210         CONTINUE
              GOTO 301
 61           CONTINUE
              WRITE(NOUT,2010)
 2010         FORMAT(/     '  IN CONSTRUCTON OF ETA AND ZETA,I=J=K  IS N
     *OT POSSIBLE DO THE 2ND SET OF LOOPS'/'  ERROR IN LOOPS SOMEWHERE')
              CALL ABORTS
 71           CONTINUE
              WRITE(NOUT,2011)
 2011         FORMAT(/     '  IN CONSTRUCTON OF ETA AND ZETA,I=J=K=L IS
     *NOT POSSIBLE DO THE 2ND SET OF LOOPS'/'  ERROR IN LOOPS SOMEWHERE'
     *        )
              CALL ABORTS
 301          CONTINUE
 3180       CONTINUE
            DO 3220 LOR=LFERM,KOR
              ltype=motyp(lor)
              lpos=mopos(ltype)
              IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO
     *        302
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO 302
              LO=IOFF(LOR)
              IL=IO+LOR
              JL=JO+LOR
              KL=KO+LOR
              LL=LO+LOR
              IJKLT=IJKT
              IF (IJKEQ .AND. KOR.GT.LOR) IJKLT=IJKLT+1
              IF (KOR.EQ.LOR) IJKLT=IJKLT+3
              LAD=KAD+LADD(LOFF+LOR)
              IF (IJKLT.GT.7 .OR. IJKLT.LT.1) THEN
                WRITE(NOUT,2004) IJKLT,IOR,JOR,
     *          KOR,LOR,IBF,JBF,KBF,LBF,ISM,JSM,KSM,LSM
                CALL ABORTS
              END IF
              IF (LAD .GT. NMAX) THEN
                WRITE(NOUT,2005) LAD
                CALL ABORTS
              END IF
              GOTO (12,22,32,42,52,62,72)IJKLT
 12           CONTINUE
C  THIS TYPE (ALL DIFFERENT) DOES NOT CONTRIBUTE
              GOTO 302
 22           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              ETA(IL)=ETA(IL)+ALPHA(IJ)*VAL2+BETA(IJ)*VAL1
              EAVF(IL)=EAVF(IL)+F(JOR)*(A2*VAL2-VAL1)
cets121890              DO 3230 K2=1,NORBS
cets121890                KL2=IOFF(MAX(K2,JOR))+MIN(K2,JOR)
cets121890                ZETA(IL,K2)=ZETA(IL,K2)+ALPHA(KL2)*VAL2+BETA(KL2)*
cets121890     *          VAL1
              DO 3230 K2=1,ntypes
                KL2=IOFF(MAX(mopos(k2),jpos))+MIN(mopos(k2),jpos)
                ZETA(k2,il)=ZETA(k2,il)+ALPHA(KL2)*VAL2+BETA(KL2)*
     *          VAL1
 3230         CONTINUE
              GOTO 302
 32           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              ETA(KL)=ETA(KL)+ALPHA(IK)*VAL2+BETA(IK)*VAL1
              EAVF(KL)=EAVF(KL)+F(IOR)*(A2*VAL2-VAL1)
cets121890              DO 3240 J2=1,NORBS
cets121890                IJ2=IOFF(MAX(IOR,J2))+MIN(IOR,J2)
cets121890                ZETA(KL,J2)=ZETA(KL,J2)+ALPHA(IJ2)*VAL2+BETA(IJ2)*
cets121890     *          VAL1
              DO 3240 J2=1,ntypes
                IJ2=IOFF(MAX(ipos,mopos(j2)))+MIN(ipos,mopos(j2))
                ZETA(j2,kl)=ZETA(j2,kl)+ALPHA(IJ2)*VAL2+BETA(IJ2)*
     *          VAL1
 3240         CONTINUE
              GOTO 302
 42           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              ETA(IJ)=ETA(IJ)+ALPHA(IL)*VAL2+BETA(IL)*VAL1
              EAVF(IJ)=EAVF(IJ)+F(KOR)*(A2*VAL2-VAL1)
cets121890              DO 3250 K2=1,NORBS
cets121890                KL2=IOFF(MAX(K2,LOR))+MIN(K2,LOR)
cets121890                ZETA(IJ,K2)=ZETA(IJ,K2)+ALPHA(KL2)*VAL2+BETA(KL2)*
cets121890     *          VAL1
              DO 3250 K2=1,ntypes
                KL2=IOFF(MAX(mopos(k2),lpos))+MIN(mopos(k2),lpos)
                ZETA(k2,ij)=ZETA(k2,ij)+ALPHA(KL2)*VAL2+BETA(KL2)*
     *          VAL1
 3250         CONTINUE
              GOTO 302
 52           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              VAL3=INTS(LAD+3)
              ETA(IL)=ETA(IL)+(ALPHA(II)+BETA(II))*VAL1
     1        +(ALPHA(IL)+BETA(IL))*VAL2+F(IOR)*VAL3
              EAVF(IL)=EAVF(IL)+F(LOR)*(VAL2) + VAL3
              EAVF(IL)=EAVF(IL)+F(IOR)*(VAL1)
cets121890              DO 3260 J2=1,NORBS
cets121890                IJ2=IOFF(MAX(IOR,J2))+MIN(IOR,J2)
cets121890                JL2=IOFF(MAX(LOR,J2))+MIN(LOR,J2)
cets121890                ZETA(IL,J2)=ZETA(IL,J2)+(ALPHA(IJ2)+BETA(IJ2))*VAL1
cets121890     1          +(ALPHA(JL2)+BETA(JL2))*VAL2+F(J2)*VAL3
              DO 3260 J2=1,ntypes
                IJ2=IOFF(MAX(ipos,mopos(j2)))+MIN(ipos,mopos(j2))
                JL2=IOFF(MAX(lpos,mopos(j2)))+MIN(lpos,mopos(j2))
                ZETA(j2,il)=ZETA(j2,il)+(ALPHA(IJ2)+BETA(IJ2))*VAL1
     1          +(ALPHA(JL2)+BETA(JL2))*VAL2+F(mopos(j2))*VAL3
 3260         CONTINUE
              GOTO 302
 62           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              ETA(II)=ETA(II)+ALPHA(IL)*VAL2+BETA(IL)*VAL1
              EAVF(II)=EAVF(II)+F(KOR)*(A2*VAL2-VAL1)
              ETA(LL)=ETA(LL)+ALPHA(IL)*VAL2+BETA(IL)*VAL1
              EAVF(LL)=EAVF(LL)+F(IOR)*(A2*VAL2-VAL1)
              ESCF=ESCF+ALPHA(IL)*VAL2*2+BETA(IL)*VAL1*2
C     WRITE(*,*) ' I,L,A,B,V2,V1,ESCF '
C     WRITE(*,*) IOR,LOR,ALPHA(IL),BETA(IL),VAL2,VAL1,ESCF
cets121890              DO 3270 J2=1,NORBS
cets121890                JL2=IOFF(MAX(J2,LOR))+MIN(J2,LOR)
cets121890                ZETA(II,J2)=ZETA(II,J2)+ALPHA(JL2)*VAL2+BETA(JL2)*
cets121890     *          VAL1
              DO 3270 J2=1,ntypes
                JL2=IOFF(MAX(mopos(j2),lpos))+MIN(mopos(j2),lpos)
                ZETA(j2,ii)=ZETA(j2,ii)+ALPHA(JL2)*VAL2+BETA(JL2)*
     *          VAL1
 3270         CONTINUE
cets121890              DO 3280 J2=1,NORBS
cets121890                IJ2=IOFF(MAX(IOR,J2))+MIN(IOR,J2)
cets121890                ZETA(LL,J2)=ZETA(LL,J2)+ALPHA(IJ2)*VAL2+BETA(IJ2)*
cets121890     *          VAL1
              DO 3280 J2=1,ntypes
                IJ2=IOFF(MAX(ipos,mopos(j2)))+MIN(ipos,mopos(j2))
                ZETA(j2,ll)=ZETA(j2,ll)+ALPHA(IJ2)*VAL2+BETA(IJ2)*
     *          VAL1
 3280         CONTINUE
              GOTO 302
 72           VAL1=INTS(LAD+1)
              VAL2=INTS(LAD+2)
              ESCF=ESCF+(ALPHA(II)+BETA(II))*VAL1+2*F(IOR)*VAL2
C     WRITE(*,*) ' I,A,B,V1,ESCF,F,V2 '
C     WRITE(*,*) IOR,ALPHA(II),BETA(II),VAL1,ESCF,F(IOR),VAL2
              ETA(II)=ETA(II)+ALPHA(II)*VAL1+BETA(II)*VAL1+F(IOR)*
     *        VAL2
              EAVF(II)=EAVF(II)+F(IOR)*(A2*VAL1-VAL1) + VAL2
cets121890              DO 3290 J2=1,NORBS
cets121890                IJ2=IOFF(MAX(IOR,J2))+MIN(IOR,J2)
cets121890                ZETA(II,J2)=ZETA(II,J2)+ALPHA(IJ2)*VAL1+BETA(IJ2)*
cets121890     *          VAL1+F(J2)*VAL2
              DO 3290 J2=1,ntypes
                IJ2=IOFF(MAX(ipos,mopos(j2)))+MIN(ipos,mopos(j2))
                ZETA(j2,ii)=ZETA(j2,ii)+ALPHA(IJ2)*VAL1+BETA(IJ2)*
     *          VAL1+F(mopos(j2))*VAL2
 3290         CONTINUE
 302          CONTINUE
 3220       CONTINUE
 3170     CONTINUE
 3120   CONTINUE
C  SKIP OVER THE  OCC,UOC / UOC,UOC  TERMS
        DO 3300 JOR=NUOC,1,-1
          IF (IGRP .NE. IJGRP(IOFF(IOR)+JOR)) THEN
CIBM        CALL RGETSA(IT52,IX)
CIBM        WRITE(NOUT,2001) IX
            CALL SREAD(IT52,INTS,NMAX*2)
            SRIN=SRIN+NMAX*2/112
            IGRP=IJGRP(IOFF(IOR)+JOR)
          END IF
 3300   CONTINUE
 3110 CONTINUE
C  GET THE TWO ELECTRON CONTRIBUTIONS TO ETA AND ZETA FROM THE EXTERNAL
      DO 3310 IOR=NUOC,1,-1
        IO=IOFF(IOR)
        II=IO+IOR
        ISM=ORBSYM(IOR)
        ITYP=1
        DO 3320 JOR=IOR,1,-1
          JO=IOFF(JOR)
          IJ=IO+JOR
          JJ=JO+JOR
          IF (IGRP .NE. IJGRP(IOFF(IOR)+JOR)) THEN
CIBM        CALL RGETSA(IT52,IX)
CIBM        WRITE(NOUT,2001) IX
            CALL SREAD(IT52,INTS,NMAX*2)
            SRIN=SRIN+NMAX*2/112
            IGRP=IJGRP(IOFF(IOR)+JOR)
          END IF
          JSM=ORBSYM(JOR)
          IJEQ=IOR.EQ.JOR
          IJT=ITYP
          IF (IJEQ) IJT=IJT+2
          IJSYM=IEOR((ISM-1),(JSM-1))
          KOFF=IJSYM*NORBS
          JAD=IJADD(IJ)
          DO 3330 KOR=1,JOR
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
            DO 3340 LOR=1,KOR
              LO=IOFF(LOR)
              IL=IO+LOR
              JL=JO+LOR
              KL=KO+LOR
              LL=LO+LOR
              IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GOTO
     *        9
              LSM=ORBSYM(LOR)
              IF ((LSM-1) .NE. IJKSYM) GOTO  9
              IJKLT=IJKT
              IF (IJKEQ .AND. KOR.GT.LOR) IJKLT=IJKLT+1
              IF (KOR.EQ.LOR) IJKLT=IJKLT+3
              LAD=KAD+LADD(LOFF+LOR)
              IF (IJKLT.GT.7 .OR. IJKLT.LT.1) THEN
                WRITE(NOUT,2004) IJKLT,IOR,JOR,
     *          KOR,LOR,IBF,JBF,KBF,LBF,ISM,JSM,KSM,LSM
                CALL ABORTS
              END IF
              IF (LAD .GT. NMAX) THEN
                WRITE(NOUT,2005) LAD
                CALL ABORTS
              END IF
              GOTO (9,9,9,9,5,9,7)IJKLT
C  ONLY 2 E TERMS, I.E., H(IJ), CONTRIBUTE
 5            VAL3=INTS(LAD+3)
              ETA(IL)=ETA(IL)+F(IOR)*VAL3
              EAVF(IL)=EAVF(IL) + VAL3
cets121890              DO 3350 J2=1,NORBS
cets121890                ZETA(IL,J2)=ZETA(IL,J2)+F(J2)*VAL3
cets121890 3350         CONTINUE
              DO 3350 J2=1,ntypes
                ZETA(j2,il)=ZETA(j2,il)+F(mopos(j2))*VAL3
 3350         CONTINUE
              GOTO 9
 7            VAL2=INTS(LAD+2)
              ETA(II)=ETA(II)+F(IOR)*VAL2
              EAVF(II)=EAVF(II) + VAL2
cets121890              DO 3360 J2=1,NORBS
cets121890                ZETA(II,J2)=ZETA(II,J2)+F(J2)*VAL2
cets121890 3360         CONTINUE
              DO 3360 J2=1,ntypes
                ZETA(j2,ii)=ZETA(j2,ii)+F(mopos(j2))*VAL2
 3360         CONTINUE
 9            CONTINUE
 3340       CONTINUE
 3330     CONTINUE
 3320   CONTINUE
 3310 CONTINUE
      WRITE(NOUT,2002) ESCF,REPNUC,ESCF+REPNUC
CTJL  WRITE(*,2003)
 2002 FORMAT(/'  ELECTRONIC ENERGY:',F16.10
     1       /'     NUCLEAR ENERGY:',F16.10
     2       /'   TOTAL SCF ENERGY:',F16.10)
 2003 FORMAT(/'  ETA AND ZETA ARE COMPLETE')
      IF (IAND(IPRINT,128) .NE. 0) THEN
        WRITE(NOUT,1904)
 1904   FORMAT(/'  ETA MATRIX:')
        CALL PRINT(ETA,NUMIJ,NORBS,NOUT)
        WRITE(NOUT,*) ' EAVF MATRIX '
        CALL PRINT(EAVF,NUMIJ,NORBS,NOUT)
CTJL  IF (IAND(IPRINT,128) .NE. 0) THEN
cets        DO 3370 IORB=1,norbs
cets          WRITE(NOUT,1905) IORB
cets 1905     FORMAT(/'  THE ZETA MATRIX FOR IORB=',I5)
cets          CALL PRINT(ZETA(1,IORB),NUMIJ,NORBS,NOUT)
cets 3370   CONTINUE
        call matout(zeta,ntypes+1,numij,ntypes,numij,6)
      END IF
C
C  CONSTRUCT A AND B0 MATRICIES; B0 IS IN RIGHT SIDE OF A MATRIX
C
 1100 FORMAT(/'  IN THE CONSRTUCTION OF THE A MATRIX, THE CASE ',A7
     1/'  IS NOT POSSIBLE IN THE ',A11,' SET OF LOOPS.'/'  ERROR IN LOOP
     * STRUCTURE SOMEWHERE.')
      CALL SREW(IT52)
      IF(.NOT.CCSD)CALL SREAD(IT52,N112,112)
      INTSRT = 0
      CALL RGETSA(IT52,INTSRT)
      WRITE(ICHK,*) ' IN GETINT AND INTSRT = ',INTSRT
CIBM
CIBM  CALL RSETSA(IT52,INTSRT)
      IGRP=0
C
C     COMMENT OUT THESE LINES FOR CPHFL
C
CGES  DO 3380 I=1,IJIND2
CGES    DO 3390 J=1,IJIND3
C         A(I,J)=A0
C3390   CONTINUE
C3380 CONTINUE
C
      RETURN
      END
