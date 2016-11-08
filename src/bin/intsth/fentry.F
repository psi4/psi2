cges  PROGRAM INTS
      subroutine fentry(ra,ia,maxcor)
C            (PINTS)
C**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
C**********************************************************************
c      by: Edward Seidl (cets071991)
c    date: July 19, 1991
c  reason: make mblu = maxcor/2 rather than hardwiring it
c
c      by: Edward Seidl (cets052691)
c    date: may 5, 1991
c  reason: add boolean keywords oneint and twoint to allow the skipping
c          of the actual calculation of oei's and tei's a la 
c          wesley allen's option
C**********************************************************************
c      by: ges
c    date: aug 22 1989
c  reason: fix bug for 5d,7f,9g when using machine other than ibm
c          subroutine fixed is togcs. see comment there.
C**********************************************************************
C      BY: RICHARD REMINGTON                    SEARCH: C04-19-89     *
C    DATE: April     19, 1989                                         *
C      BY: TPH                                                        *
C    DATE: MARCH 24, 1989                                             *
C  REASON: Increase MPRU, MBLU, and MCCU to the values in the old     *
C          INTS program.  Also make the dimensioning of ICB           *
C          correct.                                                   *
C**********************************************************************
C      BY: RICHARD REMINGTON                    SEARCH: CMachineSize  *
C    DATE: March     14, 1989                                         *
C  REASON: CHANGE CORE TO RUN IN 7MB ON 9370 OR 15MB ON 3090          *
C**********************************************************************
C    DATE: OCTOBER 25, 1988                                           *
C      BY: TPH                                                        *
C  REASON: ARRAY NPERSH TO TELL IF D5 OR D6, ETC.  ALL SHELLS MUST    *
C          BE D5 OR D6, F7 OR F10,ETC.  THIS IS FOR COMPARISON WORK   *
C          with old calculations or GAUSSIAN.  NPERSH(5)              *
c          Default is now d6, f10, g15.                               *
C**********************************************************************
C    DATE: OCTOBER 19, 1988                                           *
C      BY: TPH                                                        *
C  REASON: THIS PROGRAM WILL NOW DO GENERAL CONTRACTIONS AND ALSO     *
C          D5, F7, AND G9 FUNCTIONS.  THE LIMITATION ON THE NUMBER    *
C          OF GENERAL CONTRACTIONS IS 6. PAST THE INTEGRALS CODE,     *
C          YOU ONLY NEED TO KEEP STRAIGHT THE QUANTITIES NBFAO AND    *
C          NBFSO FOR ENERGY CALCULATION.  THIS VERSION HAS 10 ROOTS   *
C          IN COMMON ROOT IN ANTICIPATION OF CHANGES NEEDED FOR 2ND   *
C          DERIVATIVES OF G FUNCTIONS.  RIGHT NOW A SPECIAL VERSION   *
C          OF IOMR IS NEEDED.  INPUTTH$ AND INTSTH$ ARE THE EXECS     *
C          NEEDED.                                                    *
C**********************************************************************
C THIS PROGRAM IS THE INTEGRALS PACKAGE WHICH SHOULD HANDLE UP TO AND *
C INCLUDING F & G  FUNCTIONS.  UNFORTUNATELY, THE INFORMATION RELATING*
C ANY PROGRAM CHANGES HAS NOT BEEN UPDATED HERE BY ANY OF THOSE WHO   *
C DID THE PROGRAM MODIFICATIONS (WHO WAS THAT MASKED PROGRAMER?)...   *
C   ANYWAY, TO USE THIS PROGRAM REQUIRES THAT A DIFFERENT VERSION OF  *
C   INPUT BE RUN;   " PINPUT$ "                                       *
C   CURRENTLY THE PROGRAMMING FOR F FUNCTIONS STOPS AT THE SCF LEVEL  *
C   AND THE NORMAL SCF PROGRAM MAY BE RUN.                            *
C**********************************************************************
C    DATE: FEBRUARY 12, 1987                                          *
C      BY: RBR                                                        *
C  REASON: CLEAN UP OUTPUT A LITTLE AND ADD THE ABOVE ABOUT "F'S"     *
C**********************************************************************
C    DATE: MARCH 19, 1985                                             *
C      BY: WDA ??????????                                             *
C  REASON: ADD OPTION TO WRITE OUT TAPE30 INFO AND INTEGRAL TAPE      *
C          HEADER WITHOUT COMPUTING ANY INTEGRALS.                    *
C**********************************************************************
C    DATE: FEBRUARY 5, 1985       1:00 PM                             *
C      BY: RBR                                                        *
C  REASON: CHANGE TO NEW UNIT NUMBERS ! ! ! ! ! ! !!!!!!!!!!!!!!      *
C**********************************************************************
C    DATE: NOVEMBER 8, 1984      10:00 PM                             *
C      BY: TJL                                                        *
C  REASON: CONVERT TO NEW IO SYSTEM                                   *
C**********************************************************************
C
C  THIS PROGRAM CALCULATES ONE- AND TWO-ELECTRON INTEGRALS OVER
C  SYMMETRY ORBITALS OF GENERALLY CONTRACTED GAUSSIAN ORBITALS
C  AND FORMS THE SYMMETRY ORBITAL INTEGRALS ORDERED BY PK SHELL-BLOCKS.
C  THE POINT GROUPS ARE LIMITED TO D2H AND ITS SUBGROUPS.
C           BY RUSSELL M. PITZER
C
      IMPLICIT INTEGER (A-Z)
CMachineSize  (For either 7Mb or 15Mb)   PARAMETER (CRTP=600000)
cges  PARAMETER (CRTP=600000)
#include <error.h>
      integer errcod, frdboo, mk1ei, mk2ei
c     REAL*8 RA(CRTP),SEC
      REAL*8 RA(maxcor),SEC
      dimension ia(1)
      INTEGER FCB10(16)
C
      COMMON /INPUT/  I24,NGEN,NAORDS,NCONS,ITOL,ICUT,LABEL(26)
      COMMON /OUT/    IOUT,ICHECK,IFLAG3,NSECT,NINTS
      common /parmi/ mblu, mccu, mconu, mcu, mcxu,
     1   mrcru, mstu, msu, nblu, ng, ns, nst
      common /ntgr/ kaords, mconsu, mgcsu, mru, mcsu, mctu, mcru,
     1   msfu, mgu, msfru, mnru, mnsfu, mpru, ngcs, nu, mcxu2,
     2   lxyzir(3), inam(5), nnam, mdum(30)
      COMMON /DIM21/  IPQ(256)
      COMMON /TAPE/   NTAPE
      COMMON /SALCS/  NIRRED
      COMMON /SYMTRY/ NTRANS
      COMMON /PSIZE/  MSQ,KAORDQ,MCONSQ,MGCSQ,MSTQ,MRQ,MCSQ,MCTQ,MCONQ,
     1                MCQ,MSFQ,MGQ,MSFRQ,MNSFQ,MPRQ,MCXQ,MCCQ,MBLQ,MRCRQ
      COMMON /INFOA/  NAT,NUM,NQ,NX,NE,NA,NB
      COMMON /NSHEL/  NSHELL
      COMMON /MAXDIM/ MXAT,MXTRN,MXPR,MXSH,MXSYM,MXCF,MXT,MXISO,MXPSYM
      COMMON /TAPE30/ I10(200),ITAP30
C
C***********************************************************************
C*                                                                     *
cges  COMMON /CORE/   RA
C
cges  DIMENSION IA(1)
cges  EQUIVALENCE (RA,IA)
C*                                                                     *
C*    ----- NOTE ABOVE EQUIVALENCE, CRUCIAL TO DIVIDING CORE.          *
C*          IA SHOULD BE DIMENSIONED TWICE AS LARGE AS RA, AND         *
C*          CORTOP SHOULD BE DIMENSION OF RA                    -----  *
C***********************************************************************
C
C
      call psinit('APPEND')
C      CALL TSTART(3)
      CALL TSTART(6)
C
      NASK=190
      CALL IZERO(I24,32)
c      CALL IZERO(KAORDS,61)
      call izero(kaords,55)
      NTAPE=0
      NIRRED=0
      NTRANS=0
c      CALL IZERO(MSQ,18)
      call izero(msq,19)
c      CALL IZERO(NAT,9)
      call izero(nat,7)
      NSHELL=0
      CALL IZERO(MXAT,9)
      CALL IZERO(I10,201)
C
      IFLAG3 = 1
      ICHECK = 3
      IOUT   = 6
CWA   LINES ADDED BY WA TO ALLOW EXIT AFTER WRITING INPUT.
      ITAP5=5
      CALL LOCATE (ITAP5,'# INTS ###',IERLOC)
      IF(IERLOC.NE.0) THEN
          NWRTFG=0
          GO TO 3333
      END IF
      READ(ITAP5,3332) NWRTFG
 3332 FORMAT(I5)
 3333 CONTINUE
CENDWA
cets052691
      mk1ei=1
      mk2ei=1

      errcod = frdboo('MK_1EI',mk1ei)
      errcod = frdboo('MK_2EI',mk2ei)
cets
C
      NINTS=0
      NSECT=0
C
C     ----- READ IN HONDO ARRAYS -----
C
C     REWIND (10)
C     READ (10) MSQ,KAORDQ,MCONSQ,MGCSQ,MSTQ,MRQ,MCSQ,MCTQ,MCONQ,
C    #           MCQ,MSFQ,MGQ,MSFRQ,NAT,NUM,NQ,NX,NE,NB,
C    #           NSHELL,NIRRED,NTRANS,MXAT,MXTRN,MXPR,MXSH,MXSYM,
C    #           MXCF,MXT,MXISO,MXPSYM
      ITAP30=30
      CALL RFILE(ITAP30)
      CALL WREADW(ITAP30,I10,200,101,JUNK)
C
      MSQ=I10(6)
      KAORDQ=I10(7)
      MCONSQ=I10(8)
      MGCSQ=I10(9)
      MSTQ=I10(10)
      MRQ=I10(11)
      MCSQ=I10(12)
      MCTQ=I10(13)
      MCONQ=I10(14)
      MCQ=I10(15)
      MSFQ=I10(16)
      MGQ=I10(17)
      MSFRQ=I10(18)
      NAT=I10(19)
      MUL=I10(21)
      NUM=I10(22)
      NQ=I10(23)
      NX=I10(24)
      NE=I10(25)
      NB=I10(26)
      NSHELL=I10(27)
      NIRRED=I10(28)
      NTRANS=I10(29)
      MXAT=NAT
      MXTRN=NTRANS
      MXPR=I10(32)
      MXSH=NSHELL
      MXSYM=NIRRED
      MXCF=I10(35)
      MXT=I10(36)
      MXISO=I10(40)
      MXPSYM=I10(39)
      MRCRQ=I10(44)
C
      WRITE(ICHECK,997)MSQ,KAORDQ,MCONSQ,MGCSQ,MSTQ,MRQ,MCSQ,MCTQ,
     1  MCONQ,MCQ,MSFQ,MGQ,MSFRQ,NAT,NUM,NQ,NX,NE,NB,NSHELL,
     2  NIRRED,NTRANS,MXAT,MXTRN,MXPR,MXSH,MXSYM,MXCF,MXT,MXISO,
     3  MXPSYM,MRCRQ
  997 FORMAT(10I5)
C
      MSU=MSQ
      KAORDS=KAORDQ
      MCONSU=MCONSQ
      MGCSU=MGCSQ
      MSTU=MSTQ
      MRU=MRQ
      MCSU=MCSQ
      MCTU=MCTQ
      MCONU=MCONQ
      MCU=MCQ
      MSFU=MSFQ
      MGU=MGQ
      MSFRU=MSFRQ
      MNSFU=MSFU*(MSFU+1)/2
C04-19-89  MPRU=11000
cets071991      MPRU  = 22000
      mpru = 66000
C04-19-89  MCCU=182
      MCCU  = 364
CMachineSize  (For  7Mb)  MBLU  = 300000
CMachineSize  (For 15Mb)  MBLU  = 400000
cges change mblu from 300000 to 400000 oct 1 89.
cges change mblu from 400000 to 530000 22 dec 89.
c     MBLU  = 400000
cets071991      MBLU  = 530000
      mblu = maxcor/2
      MRCRU = MRCRQ
      MNRU  = 1000
C
      NS=MSU
      NAORDS=KAORDS
      NCONS=MCONSU
      NGCS=MGCSU
      NST=MSTU
      NGEN=MGU
      NG=MGU
C
C
      IF (IFLAG3.GT.0) WRITE (ICHECK,104) NGEN,NS,NAORDS,NCONS,NGCS,ITOL
     #,                              ICUT,NTAPE,MCONU,MCTU,MSFU,MSFRU
CCCC  WRITE (IOUT,104) NGEN,NS,NAORDS,NCONS,NGCS,ITOL,ICUT,NTAPE,MCONU
CCCC #,                MCTU,MSFU,MSFRU
 
      ITOL=20
      NCUT=9
      NTAPE=34
C
C     ----- ARRANGE CORE STORAGE BASED ON ABOVE PARAMETERS -----
C
CSEL  CALL MAXCOR(CORLEN)
CSEL  CORLEN=(CORLEN-1000)
CSEL  CALL GETCOR(CORLEN,CORADD)
CSEL  CORLEN=CORLEN/INTOWP(1)
CSEL  CALL LOCF(RA,RADD)
CCSEL WRITE (IOUT,'(//,' CORADD=',I9,' RADD=',I9)') CORADD,RADDA
CSEL  CORADD=CORADD/INTOWP(1)
CSEL  RADD=RADD/INTOWP(1)
CSEL  CORBOT=CORADD-RADD+1
CSEL  CORTOP=CORBOT+CORLEN-1
cges  CORTOP=CRTP
      cortop=maxcor
C^^^^^^^^^^^^^^^^^^^^^^
C
C     ----- ZERO CORE NOW -----
C
      CALL ZERO(RA,CORTOP)
C
      NF=1
      NC=NF+MSU
      NIR=NC+MSU
      NCON=NIR+KAORDS
      NRCR=NCON+MCONSU
      LMNP1=NRCR+MCONSU
      MAORDS=LMNP1+MCONSU
      ICXSV1=MAORDS+MGCSU
      ICXSV2=ICXSV1+MGCSU
      NSO=ICXSV2+MGCSU**2
      ND=NSO+MSTU
      NSOPR=ND+MSTU
      NBLPR=NSOPR+MSTU
      LA=NBLPR+MSTU
      NTL=LA+MRU*KAORDS
      NTU=NTL+MSFU
      NT=NTU+MSFU
      NCT=NT+MSFU
      MCONS=NCT+MSFU
      MGCS=MCONS+MSFU
      ICA=MGCS+MSFU
      LB=ICA+MCU*MSU*MGU
      MS=LB+MSFRU
      MNL=MS+MSFRU
      IPRST=MNL+MSFRU
      NPAIR=IPRST+MNSFU
      NFCT=NPAIR+2*MNSFU
      NGW=NFCT+MCCU
      ICB=NGW+MCCU
ctph  NPRIR=ICB+4*MGU*MCCU
      NPRIR=ICB+4*24*MCCU
      IDP=NPRIR+2*MSTU*MNSFU
      MTYPE=IDP+MSTU**3
      ITYP=MTYPE+2*MSU
      IL=ITYP+MSTU
      ICXST=IL+MPRU
      MAU=ICXST+2*MNSFU
      MICS=MAU+MNRU
      MJCS=MICS+MNRU
      NIRL=MJCS+MNRU
      INTTOP=NIRL+KAORDS
C
C     ----- NOW WORK FOR REAL ARRAYS, STARTING AT REALST (JUST ABOVE
C                                            INTEGER ARRAYS) -----
      REALST=IADTWP(INTTOP)
      IGC=REALST
      CHG=IGC+MCSU*MCTU*MGCSU
      ZET=CHG+MSU
      ETA=ZET+MCONU*MCONSU
      X=ETA+MRCRU*MCONU*MCONSU
      Y=X+MCU*MSU
      Z=Y+MCU*MSU
      AJ=Z+MCU*MSU
      CX=AJ+MBLU
C
C     ----- CX CAN NOW TAKE REST OF CORE, SO -----
C
      MCXU=CORTOP-CX-2
      IF (IFLAG3.GT.0) WRITE(ICHECK,105) INTTOP,REALST,CX,MCXU
CCCC  WRITE (IOUT,105) INTTOP,REALST,CX,MCXU
      IF (MCXU.LE.0) WRITE (IOUT,106) CORTOP,CX
      IF (MCXU.LE.0) THEN
        WRITE(ICHECK,*) ' NOT ENOUGH CORE SPACE TO CONTINUE, CX'
        call qabort
      END IF
C
C        SET UP IPQ(I)=(I*(I-1))/2 ARRAY
C
      IPQ(1)=0
      DO 80 I=1,255
   80 IPQ(I+1)=IPQ(I)+I
C
C     ----- ALLOCATE SPACE FOR HONDO ARRAYS. NOTE OVERLAP
C                AJ AND CX
C
C     ----- REAL ARRAYS FIRST -----
C
      EX=AJ
      CSPD=EX+MXPR
      CC=CSPD+MRCRU*MXPR*5
      ZAN=CC+MXCF
      C=ZAN+MXAT
      USOTAO=C+3*MXAT
C
C     ---- INTEGER EQUIVALENTS -----
C
      REX=WPADTI(EX)
      RCSPD=WPADTI(CSPD)
      RCC=WPADTI(CC)
      RZAN=WPADTI(ZAN)
      RC=WPADTI(C)
      RSOTAO=WPADTI(USOTAO)
C
C     ----- INTEGER ARRAYS -----
C
      ICT=WPADTI(USOTAO+NUM*MSFRU)
      NSQ=ICT+MXAT*MXTRN
      KS=NSQ+MXAT
      NGNCN=KS+MXAT
      KSTART=NGNCN+MXSH
      KATOM=KSTART+MXSH
      KTYPE=KATOM+MXSH
      KNG=KTYPE+MXSH
      KMIN=KNG+MXSH
      KMAX=KMIN+MXSH
      NPSH=KMAX+MXSH
      ISC=NPSH+5
      IPC=ISC+MXSH*MXPSYM
      SLAB=IPC+MXSH*MXPSYM
      LAB=SLAB+MXSYM
      LOC=LAB+MXCF
      LOC2=LOC+MXAT*5
      KLASS=LOC2+MXAT*5
      LAMBDA=KLASS+MXAT
      IADDR=LAMBDA+MXSYM
      TITLE=IADDR+MXSYM
      LBLAT=TITLE+20
      ITRNAO=LBLAT+2*MXAT
ctph  NSOT=ITRNAO+MSFRU
      NSOT=ITRNAO+NUM
      NEED=NSOT+MSTU
C
      IF (NEED.GT.INTOWP(CORTOP)) THEN
        WRITE(ICHECK,*) ' NOT ENOUGH CORE FOR HONDO'
        call qabort
      END IF
C
      CALL IN(IA(ICT),IA(NSQ),IA(KS),IA(REX),IA(RCSPD),IA(NGNCN),
     1  IA(KSTART),IA(KATOM),IA(KTYPE),IA(KNG),IA(KMIN),IA(KMAX),
     2  IA(NPSH),IA(ISC),IA(IPC),IA(SLAB),IA(RCC),IA(LAB),IA(LOC),
     3  IA(LOC2),IA(KLASS),IA(LAMBDA),IA(IADDR),IA(RZAN),IA(RC),
     4  IA(TITLE),IA(LBLAT))
C
C     ----- FORM PITZER ARRAYS FROM HONDO INFO -----
C
ctph  write(6,*) ' ex before topitz ',(ra(ex+i-1),i=1,mxpr)
      CALL TOPITZ(IA(NF),IA(NC),IA(NIR),IA(NCON),IA(NRCR),IA(LMNP1),
     1  IA(NSO),IA(ND),IA(NSOPR),IA(NBLPR),IA(LA),RA(IGC),IA(NTL),
     2  IA(NTU),IA(NT),IA(NCT),IA(MCONS),IA(MAORDS),IA(MGCS),IA(ICA),
     3  IA(LB),IA(MS),IA(MNL),IA(IDP),IA(MTYPE),IA(ITYP),RA(CHG),
     4  RA(ZET),RA(ETA),RA(X),IA(NIRL),RA(ZAN),RA(C),IA(ICT),IA(NSQ),
     5  IA(KS),RA(EX),RA(CSPD),IA(NGNCN),IA(KSTART),IA(KATOM),
     6  IA(KTYPE),IA(KNG),IA(KMIN),IA(KMAX),IA(ISC),IA(IPC),IA(SLAB),
     7  RA(CC),IA(LAB),IA(LOC),IA(LOC2),IA(KLASS),IA(LAMBDA),
     8  IA(IADDR),IA(LBLAT),IA(NPSH))
c
c  input processing
c
      CALL SYMINP(RA(IGC),RA(CHG),RA(ETA),IA(ICA),IA(IDP),IA(ITYP),
     1  IA(LA),IA(LB),IA(LMNP1),IA(MAORDS),IA(MCONS),IA(MGCS),
     2  IA(MNL),IA(MS),IA(MTYPE),IA(NBLPR),IA(NC),IA(NCON),IA(NCT),
     3  IA(ND),IA(NF),IA(NIR),IA(NRCR),IA(NSO),IA(NSOPR),IA(NT),
     4  IA(NTL),IA(NTU),RA(X),RA(Y),RA(Z),RA(ZET),IA(TITLE))
C
C     ----- FORM SO-TO-AO TRANSFORMATION MATRIX -----
C
      CALL SOTOAO(IA(NC),IA(NF),IA(MCONS),IA(ITRNAO),IA(NSOT),
     1  IA(NSO),RA(USOTAO),IA(MGCS),IA(MAORDS),IA(NIR),IA(NCT),
     2  IA(NRCR),IA(LA),RA(IGC),IA(LMNP1),IA(RSOTAO))
C
C     ----- WRITE OUT INTER-ATOMIC DISTANCES -----
C
      CALL INTR(RA(C),MXAT,IOUT)
c
c  symmetry coefficient products
c
      CALL SOCFPD(RA(IGC),RA(CX),IA(ICXST),IA(ICXSV1),IA(ICXSV2),
     1  IA(IDP),IA(IL),IA(IPRST),IA(LA),IA(LB),IA(MAORDS),IA(MAU),
     2  IA(MCONS),IA(MGCS),IA(MICS),IA(MJCS),
     3  IA(NBLPR),IA(NC),IA(NCT),IA(ND),IA(NF),IA(NIR),IA(NPAIR),
     4  IA(NPRIR),IA(NRCR),IA(NT),IA(NTL),IA(NTU))
C     PRINT '(//' RETURNED FROM SOCFPD'//)'
C
CWA  SKIP COMPUTATION OF INTEGRALS IF DESIRED.
      IF(NWRTFG.EQ.1) GO TO 3334
CENDWA
c
      do 20 i=1,3
   20 inam(i)=i
c
c  one-electron integrals
c
cets052691
      if(mk1ei.eq.1) then
        DO 10 N=1,3
        nnam=inam(n)
   10 CALL ONEINT(RA(AJ),RA(CHG),RA(CX),RA(ETA),IA(ICA),IA(ICB),
     1  IA(ICXST),IA(IL),IA(IPRST),IA(LMNP1),IA(MCONS),IA(NC),
     2  IA(NCON),IA(NF),IA(NFCT),IA(NGW),IA(NPAIR),IA(NPRIR),
     3  IA(NRCR),IA(NT),IA(NTL),IA(NTU),RA(X),RA(Y),RA(Z),RA(ZET))
      else
        write(6,*)' skipping calculation of one-electron integrals'
      endif
C
      if(mk2ei.eq.1) then
        CALL TWOINT(RA(AJ),RA(CX),RA(ETA),IA(ICA),IA(ICB),IA(ICXST),
     1  IA(IL),IA(IPRST),IA(LMNP1),IA(MCONS),IA(NC),IA(NCON),IA(NF),
     2  IA(NFCT),IA(NGW),IA(NPAIR),IA(NPRIR),IA(NRCR),IA(NT),IA(NTL),
     3  IA(NTU),RA(X),RA(Y),RA(Z),RA(ZET))
      else
        write(6,*)' skipping calculation of two-electron integrals'
      endif
C
      WRITE(IOUT,222)NINTS,NSECT,NTAPE
      IF (IFLAG3.GT.0) WRITE(ICHECK,222)NINTS,NSECT,NTAPE
  222 FORMAT(/' TOTAL OF ',I10,'  INTEGRALS IN',I10,
     1  ' SECTORS ARE WRITTEN ON TAPE',I5)
C     CALL ETIME2(ENDD)
 3334 CONTINUE
      ITAPE1 = 93
      JCODE = 3
      CALL RCLOSE(ITAP30,JCODE)
CTJL  CALL RCLOSE(ITAPE1,JCODE)
      CALL RCLOSE(NTAPE,JCODE)
C
C      CALL TSTOP(3)
      CALL TSTOP(6)
      call psdone()
C
      stop
C
C
  101 FORMAT (26A3)
  102 FORMAT (1H1,26A3)
  103 FORMAT (26I3)
  104 FORMAT (1X,26I3)
  105 FORMAT (//,1X,I6,' WORDS USED FOR INTEGER ARRAYS, SO REAL ARRAYS '
     #,       'START AT',I7,/,7X,' CX STARTS AT',I7,' AND HAS SPACE FOR'
     #,       I7,' ELEMENTS')
  106 FORMAT (//,' ***** NOT ENOUGH CORE TO STORE ARRAYS (!). CORTOP:'
     #,       I7,' AND CX:',I7,//)
C
      END
