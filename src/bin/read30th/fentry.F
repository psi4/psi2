c $Id: fentry.F,v 1.8 1994/07/22 18:33:18 crawdad Exp $
c
c $Log: fentry.F,v $
c Revision 1.8  1994/07/22 18:33:18  crawdad
c Added changes by Roger Grev such that read30th can now read general
c contractions from file30.
c
c Revision 1.7  1991/10/08  20:15:07  psi
c no change
c
c Revision 1.6  1991/09/11  00:53:33  seidl
c test
c
c Revision 1.5  1991/08/23  08:16:20  seidl
c *** empty log message ***
c
c Revision 1.4  1991/08/21  21:29:16  psi
c *** empty log message ***
c
cRevision 1.3  1991/08/21  21:25:47  psi
cinitial revision
c
cges  PROGRAM READ30
      subroutine fentry (cc,ic,maxcor)
C     THIS PROGRAM CHECKS FILE30
c
C**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
C**********************************************************************
C   BY:  Roger Grev                                search:  crsg
C   DATE: December 15,1991
C   REASON: read file30 for general contractions from inputth
C**********************************************************************
C*  BY:  TRACY HAMILTON                            SEARCH:  CTPH      *
C*  DATE:  NOVEMBER 11,1988                                           *
C*  REASON: READ FILE30 FROM INPUTTH                                  *
C**********************************************************************
C*  BY:  RICHARD REMINGTON                         SEARCH:  C3-26-88  *
C*  DATE:  MARCH  26,  1988                                           *
C*  REASON: MIGRATION TO RUN IN 7MB ON 9370                           *
C**********************************************************************
C*  LAST UPDATED MARCH 3, 1985 BY ROBERT MACLAGAN                     *
C**********************************************************************
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
c
crsg  REAL*8 ZAN(200),C(3,200),EX(200),CSPD(200,5),FTR(10,10,200)
crsg
crsg  cspd dimension changed to accomodate generalized contractions
crsg
      REAL*8 ZAN(200),C(3,200),EX(200),CSPD(10,200,5),FTR(10,10,200)
      REAL*8 CC(200),T(200),PTR(3,3,200),DTR(6,6,200),TEMP(100)
      REAL*8 SU(5000),EIG(100),WMOT(5000),USOTAO(100,100),DE(300)
      REAL*8 ALPHA(10),BETA(10),GTR(15,15,200),NGNCN(200)
      INTEGER ICT(200,200),NS(200),KS(200),KSTART(200)
      INTEGER KATOM(200),KTYPE(200),KNG(200),KLOC(200)
      INTEGER KMIN(200),KMAX(200),ISC(200,200)
      INTEGER IPC(200,200),LAB(200),LOC(200,5),LOC2(200,5)
      INTEGER KLASS(200),LAMBDA(200),IADDR(200),IPP(200)
      INTEGER INVT(200),ISO(200,200),TITLE(20),ITEMP(200)
      INTEGER SB(200),LBLAT(2,200),NPERSH(5)
      INTEGER ITYP(200),NLAMDA(200),NC(200),NO(200)
      EQUIVALENCE (ITEMP,TEMP)
C
      call psinit('APPEND')
      ITAP30=30
      WRITE(6,*)' CALLING RFILE -- ITAP30= ',ITAP30
      CALL RFILE(ITAP30)
C
      WRITE(6,*)' RETURNED FROM RFILE -- ITAP30 = ',ITAP30
      IEND=1
      MPOINT=200
      MCONST=200
      MCALCS=100
C
      DO 1 I=1,MPOINT
         ITEMP(I)=0
    1 CONTINUE
C
      CALL WREADW (ITAP30,ITEMP,MCONST,101,JUNK)
C
C
      ITOTAL= ITEMP(1)
C     MPOINT = ITEMP(2)
C     MCONST = ITEMP(3)
C     MCALCS = ITEMP(4)
      MSU = ITEMP(6)
      KAORDS = ITEMP(7)
      MCONSU = ITEMP(8)
      MGCSU = ITEMP(9)
      MSTU = ITEMP(10)
      MRU = ITEMP(11)
      MCSU = ITEMP(12)
      MCTU = ITEMP(13)
      MCONU = ITEMP(14)
      MCU = ITEMP(15)
      MSFU = ITEMP(16)
      MGU = ITEMP(17)
      MSFRU = ITEMP(18)
      NAT = ITEMP(19)
      ICH = ITEMP(20)
      MUL = ITEMP(21)
      NUM = ITEMP(22)
      NQ = ITEMP(23)
      NX = ITEMP(24)
      NE = ITEMP(25)
      NB = ITEMP(26)
      NSHELL = ITEMP(27)
      NIRRED = ITEMP(28)
      NT = ITEMP(29)
      MXAT = ITEMP(30)
      MXTRN = ITEMP(31)
      MXPR = ITEMP(32)
      MXSH = ITEMP(33)
      MXSYM = ITEMP(34)
      MXCF = ITEMP(35)
      MXT = ITEMP(36)
      MXISO = ITEMP(37)
      MXPSYM = ITEMP(38)
      NPSYM = ITEMP(39)
      NISO = ITEMP(40)
      NSYMHF = ITEMP(41)
      MXCOEF = ITEMP(42)
      IOPEN = ITEMP(43)
C
crsg need to read in the max number of general contractions, mrcru
      mrcru = itemp(44)
crsg
      WRITE (6,*) '          '
      WRITE (6,*) '****TITLE*******************************************'
      WRITE (6,*) '          '
      DO 2 I=1,5
         CALL WREADW (ITAP30,ITEMP,20,IEND,IEND)
         WRITE (6,21) I,(ITEMP(J),J=1,20)
   21    FORMAT (I5,5X,20A4)
    2 CONTINUE
C
C     ----- ARRAY FOR CONSTANTS -----
C
      CALL WREADW (ITAP30,ITEMP,MCONST,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****CONSTANTS***************************************'
      WRITE (6,*) '          '
crsg  WRITE (6,22) (ITEMP(J),J=1,43)
      WRITE (6,22) (ITEMP(J),J=1,44)
   22 FORMAT (10I5)
C
C     ----- ARRAY FOR POINTERS -----
C
      CALL WREADW (ITAP30,ITEMP,MPOINT,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****POINTERS TO BASIS SET INFO ETC.*****************'
      WRITE (6,*) '          '
      WRITE (6,24) (ITEMP(I),I=1,MPOINT)
   24 FORMAT (10I5)
      ITMP29 = ITEMP(29)
      ITMP30 = ITEMP(30)
C
C     ----- POINTERS TO CALCULATIONS -----
C
      CALL WREADW (ITAP30,ITEMP,MCALCS,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****POINTERS TO SUBDIRECTORIES FOR CALCS************'
      WRITE (6,*) '          '
      WRITE (6,25) (ITEMP(I),I=1,MCALCS)
   25 FORMAT (10I5)
C
C
C     ----- START WRITING HONDO ARRAYS TO FILE, PUTTING ADDRESSES
C               IN ITEMP AS WE GO
C
      CALL WREADW(ITAP30,ZAN,INTOWP(NAT),IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****NUCLEAR CHARGES - ZAN(I)************************'
      WRITE (6,*) '          '
      WRITE (6,26) (ZAN(I),I=1,NAT)
   26 FORMAT (3F20.10)
C
C     ----- TRANSFORMATION TABLE OF ATOMS -----
C
      DO 3 I=1,NT
      CALL WREADW(ITAP30,ICT(1,I),NAT,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****TRANSFORMATION TABLE OF ATOMS - ICT(J,',I,')****'
      WRITE (6,*) '          '
      WRITE (6,27) (ICT(J,I),J=1,NAT)
   27 FORMAT (10I5)
    3 CONTINUE
C
C    ----- NUMBER OF SHELLS PER ATOM -- NS(NAT) -----
C
      CALL WREADW(ITAP30,NS,NAT,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****NUMBER OF SHELLS PER ATOM - NS(I)***************'
      WRITE (6,*) '          '
      WRITE (6,28) (NS(I),I=1,NAT)
   28 FORMAT (10I5)
C
C     ----- POINTER TO SHELLS ON ATOM -- KS(NAT) -----
C
      CALL WREADW(ITAP30,KS,NAT,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****POINTER TO SHELLS ON ATOMS - KS(I)**************'
      WRITE (6,*) '          '
      WRITE (6,29) (KS(I),I=1,NAT)
   29 FORMAT (10I5)
C
C     ----- EXPONENTS OF PRIMITIVE GAUSSIANS -- EX(MXPR) -----
C
      CALL WREADW(ITAP30,EX,INTOWP(MXPR),IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****PRIMITIVE GAUSSIAN EXPONENTS - EX(I)************'
      WRITE (6,*) '          '
      WRITE (6,30) (EX(I),I=1,MXPR)
   30 FORMAT (3F20.10)
crsg
crsg  ----- S,P AND D CONTRACTION COEFFICIENTS -- CSPD(MXPR,3) -----
crsg
crsg  WRITE (6,*) '          '
crsg  WRITE (6,*) '****S,P,D,F,G CONTRACTION COEFFS - CSPD(J,I)********'
crsg  DO 4 I=1,5
crsg  CALL WREADW(ITAP30,CSPD(1,I),INTOWP(MXPR),IEND,IEND)
crsg  WRITE (6,*) '          '
crsg  WRITE (6,*) 'I =',I
crsg  WRITE (6,*) '          '
crsg  WRITE (6,31) (CSPD(J,I),J=1,MXPR)
crsg31 FORMAT (3F20.10)
crsg4 CONTINUE
crsg
crsg  this is changed to work for generalized contractions
C
C     ----- S,P AND D CONTRACTION COEFFICIENTS -- CSPD(mrcru,MXPR,3) -----
C
      WRITE (6,*) '          '
      WRITE (6,*) '****S,P,D,F,G CONTRACTION COEFFS - CSPD(k,J,I)******'
      DO 4 I=1,5
crsg    CALL WREADW(ITAP30,CSPD(1,1,I),INTOWP(mrcru*MXPR),IEND,IEND)
        WRITE (6,*) '          '
        WRITE (6,*) 'I =',I
        WRITE (6,*) '          '
        call cx(cspdfg,mrcru,mxpr,i,iend)
    4 CONTINUE
c2    DO 2344 I=1,5
crsg    CALL WREADW(ITAP30,CSPD(1,1,I),INTOWP(mrcru*MXPR),IEND,IEND)
c2      WRITE (6,*) '          '
c2      WRITE (6,*) 'I =',I
c2      WRITE (6,*) '          '
c       if(mrcru.eq.1) then
c           write(6,31) (cspd(J,1,I), J = 1,mxpr)
c       else
c2        do 2345 k=1,mrcru
c2          WRITE (6,*) '          '
c2          WRITE (6,*) '    k =',k
c2          WRITE (6,*) '          '
crsg        WRITE (6,31) (k,j,i,CSPD(k,j,I),J=1,MXPR)
crsg2       WRITE (6,31) (CSPD(k,J,I),J=1,MXPR)
c2          WRITE (6,31) (CSPD(j,k,I),J=1,MXPR)
c2345     continue
c       endif
c2 31 FORMAT (3F20.10)
crsg31 format(3(3i3,f20.10))
c2344 CONTINUE
crsg *****************************************************************
crsg ***** END CHANGES FOR GENERAL CONTRACTIONS BY RSG HERE **********
crsg *****************************************************************
C
C     ----- POINTER TO PRIMITIVES FOR A SHELL -- KSTART(NSHELL) -----
C
      CALL WREADW(ITAP30,KSTART,NSHELL,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****POINTER TO PRIMITIVES IN SHELL - KSTART(I)******'
      WRITE (6,*) '          '
      WRITE (6,32) (KSTART(I),I=1,NSHELL)
   32 FORMAT (10I5)
C
C     ----- ATOM ON WHICH NTH SHELL IS CENTRED -- KATOM(NSHELL) -----
C
      CALL WREADW(ITAP30,KATOM,NSHELL,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****ATOM ON WHICH SHELL IS CENTRED - KATOM(I)*******'
      WRITE (6,*) '          '
      WRITE (6,33) (KATOM(I),I=1,NSHELL)
   33 FORMAT (10I5)
C
C     ----- TYPE (S,P,D...) OF NTH SHELL -- KTYPE(NSHELL) -----
C
      CALL WREADW(ITAP30,KTYPE,NSHELL,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****TYPE(S,P,D ETC) OF SHELLS - KTYPE(I)************'
      WRITE (6,*) '          '
      WRITE (6,34) (KTYPE(I),I=1,NSHELL)
   34 FORMAT (10I5)
C
C     ----- NUMBER OF CONTRACTED FUNCTIONS IN NTH SHELL -- KNG(NSHELL))
C
      CALL WREADW(ITAP30,KNG,NSHELL,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****NUMBER OF PRIMITIVES PER SHELL - KNG(I)*********'
      WRITE (6,*) '          '
      WRITE (6,35) (KNG(I),I=1,NSHELL)
   35 FORMAT (10I5)
C
C     ----- POINTER TO FIRST BASIS FNCT IN SHELL --- KLOC(NSHELL) ----
C
      CALL WREADW(ITAP30,KLOC,NSHELL,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****PTR TO FIRST BASIS FNCT IN SHELL - KLOC(I)******'
      WRITE (6,*) '          '
      WRITE (6,36) (KLOC(I),I=1,NSHELL)
   36 FORMAT (10I5)
C
C     ----- FIRST FUNCTION (S,X,Y,Z,XX...) IN SHELL -- KMIN(NSHELL) --
C
      CALL WREADW(ITAP30,KMIN,NSHELL,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****1ST FN (S,X,XX,XXX...) IN A SHELL - KMIN(I)*****'
      WRITE (6,*) '          '
      WRITE (6,37) (KMIN(I),I=1,NSHELL)
   37 FORMAT (10I5)
C
C     ----- LAST FUNCTION IN SHELL -- KMAX(NSHELL) -----
C
      CALL WREADW(ITAP30,KMAX,NSHELL,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****LAST FN (S,X,XX,XXX...) IN A SHELL - KMAX(I)****'
      WRITE (6,*) '          '
      WRITE (6,38) (KMAX(I),I=1,NSHELL)
   38 FORMAT (10I5)
C
C     ----- PACKED SO TRANSFORMATION MATRICES -- ISC(NSHELL,NPSYM)
C
      NPSYM=(NIRRED+3)/4
      DO 5 I=1,NPSYM
      CALL WREADW(ITAP30,ISC(1,I),NSHELL,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****PACKED SO TRANSFORMATION MATRICES - ISC(J,',I,')'
      WRITE (6,*) '          '
      WRITE (6,39) (ISC(J,I),J=1,NSHELL)
   39 FORMAT (5I10)
    5 CONTINUE
C
C     ----- PACKED SO TRANSFORMATION INFORMATION -- IPC(NSHELL,NPSYM)
C
      DO 6 I=1,NPSYM
      CALL WREADW(ITAP30,IPC(1,I),NSHELL,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****PACKED SO TRANSFORMATION INFO - IPC(J,',I,')****'
      WRITE (6,*) '          '
      WRITE (6,39) (IPC(J,I),J=1,NSHELL)
    6 CONTINUE
C
C     ----- LABELS OF IRREPS -- SB(NIRRED) -----
C
      CALL WREADW(ITAP30,SB,NIRRED,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****LABELS OF IRREDUCIBLE REPRS - SB(I)*************'
      WRITE (6,*) '          '
      WRITE (6,20) (SB(I),I=1,NIRRED)
C
C     ----- COEFFICIENTS OF SO'S -- CC(MXCF) -----
C
      CALL WREADW(ITAP30,CC,INTOWP(MXCF),IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****COEFFS FOR SYMMETRY ORBITALS - CC(I)************'
      WRITE (6,*) '          '
      WRITE (6,26) (CC(I),I=1,MXCF)
C
C     ----- LABELS OF SO MATRICES -- LAB(MXCF) -----
C
      CALL WREADW(ITAP30,LAB,MXCF,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****PACKED LABELS FOR SYMMETRY ORBITALS - LAB(I)****'
      WRITE (6,*) '          '
      WRITE (6,20) (LAB(I),I=1,MXCF)
   20 FORMAT (20A4)
C
C     ----- POINTER TO BEGINNING OF SO MATRIX -- LOC(NAT,5) -----
C
      WRITE (6,*) '          '
      WRITE (6,*) '****PTR TO BEGINNING OF SO MATRIX - LOC(J,I)********'
      DO 7 I=1,5
      CALL WREADW(ITAP30,LOC(1,I),NAT,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) 'I =',I
      WRITE (6,*) '          '
      WRITE (6,24) (LOC(J,I),J=1,NAT)
    7 CONTINUE
C
C     ----- POINTER TO END OF SO MATRIX -- LOC2(NAT,5) -----
C
      WRITE (6,*) '          '
      WRITE (6,*) '****PTR TO END OF SO MATRIX - LOC2(J,I)*************'
      DO 8 I=1,5
      CALL WREADW(ITAP30,LOC2(1,I),NAT,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) 'I =',I
      WRITE (6,*) '          '
      WRITE (6,24) (LOC2(J,I),J=1,NAT)
    8 CONTINUE
C
C     ----- CLASS AN ATOM BELONGS TO -- KLASS(NAT) -----
C
      CALL WREADW(ITAP30,KLASS,NAT,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****CLASS AN ATOM BELONGS TO - KLASS(I)*************'
      WRITE (6,*) '          '
      WRITE (6,27) (KLASS(I),I=1,NAT)
C
C     ----- DEGENERACY OF NTH IRREP -- LAMBDA(NIRRED) -----
C
      CALL WREADW(ITAP30,LAMBDA,NIRRED,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****DEGENERACY OF THE IRREPS - LAMBDA(I)************'
      WRITE (6,*) '          '
      WRITE (6,27) (LAMBDA(I),I=1,NIRRED)
C
C     ----- POINTER TO START OF NTH SOS ? -- IADDR(NIRRED) -----
C
      CALL WREADW(ITAP30,IADDR,NIRRED,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****PTR TO START OF NTH SOS ? - IADDR(I)************'
      WRITE (6,*) '          '
      WRITE (6,27) (IADDR(I),I=1,NIRRED)
C
C     ----- NUMBER OF SOS OF NTH SYMMETRY - IPP(NIRRED) -----
C
      CALL WREADW(ITAP30,IPP,NIRRED,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****NO. OF AO''S OF EACH SYMMETRY - IPP(I)**********'
      WRITE (6,*) '          '
      WRITE (6,27) (IPP(I),I=1,NIRRED)
C
C     ----- TRANSFORMATION MATRICES -- T(MXT) -----
C
      MXT=9*NT
      CALL WREADW(ITAP30,T,INTOWP(MXT),IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****TRANFORMATION MATRICES OF COORDINATES - T(I)****'
      WRITE (6,*) '          '
      WRITE (6,26) (T(I),I=1,MXT)
C
C     ----- INVERSE OF SYMMETRY OPERATION -- INVT(NT)
C
      CALL WREADW(ITAP30,INVT,NT,IEND,IEND)
C
C     ----- PACKED TRANSFORMATION MATRIX FOR SHELLS -- ISO(NSHELL,MXISO)
C
      NISO=(NT+3)/4
      DO 9 I=1,NISO
      CALL WREADW(ITAP30,ISO(1,I),NSHELL,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****PACKED TRANSFORMATION TABLE FOR SHELLS - ISO(J,',
     *I,')'
      WRITE (6,*) '          '
      WRITE (6,39) (ISO(J,I),J=1,NSHELL)
    9 CONTINUE
C
C     ----- LABEL OF ATOMS -- LBLAT(2,NAT)
C
      CALL WREADW(ITAP30,LBLAT,2*NAT,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****LABELS OF THE ATOMS -LBLAT(I,J)*****************'
      WRITE (6,*) '          '
      DO 227 I=1,NAT
      WRITE (6,20) (LBLAT(J,I),J=1,2)
  227 CONTINUE
C
C     ----- AO TO SO TRANSFORMATION MATRIX -- USOTAO(NUM,MSFRU)
C
      CALL WREADW(ITAP30,USOTAO,INTOWP(NUM*MSFRU),ITMP29,JUNK)
      WRITE (6,*) '          '
      WRITE (6,*) '****AO TO SO TRANSFORMATION MATRIX - USOTAO(NUM,MSFRU
     *)'
      WRITE (6,*) '          '
      CALL MATOUT(USOTAO,NUM,MSFRU,NUM,MSFRU,6)
C
C     ----- TRANSFORMATION MATRICES FOR P FUNCTIONS -- PTR(3,3,NT)
C
      CALL WREADW(ITAP30,PTR,INTOWP(3*3*NT),IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****TRANSFORMATION MATRICES FOR P FNS - PTR(I,J,K)*'
      WRITE (6,*) '          '
      DO 225 K=1,NT
      DO 226 J=1,3
      WRITE (6,26) (PTR(I,J,K),I=1,3)
  226 CONTINUE
  225 CONTINUE
C
C     ----- TRANSFORMATION MATRICES FOR D FUNCTIONS -- DTR(6,6,NT)
C
      CALL WREADW(ITAP30,DTR,INTOWP(6*6*NT),IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '******TRANSFORMATION MATICES FOR D FNS - DTR(I,J,K)*'
      WRITE (6,*) '          '
      DO 223 K=1,NT
      DO 224 J=1,6
      WRITE (6,26) (DTR(I,J,K),I=1,6)
  224 CONTINUE
  223 CONTINUE
C
C     ----- TRANSFORMATION MATRICES FOR F FUNCTIONS -- FTR(10,10,NT)
C
      CALL WREADW(ITAP30,FTR,INTOWP(10*10*NT),IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '******TRANSFORMATION MATICES FOR F FNS - FTR(I,J,K)*'
      WRITE (6,*) '          '
      DO 233 K=1,NT
      DO 234 J=1,10
      WRITE (6,26) (FTR(I,J,K),I=1,10)
  234 CONTINUE
  233 CONTINUE
C
C     ----- TRANSFORMATION MATRICES FOR G FUNCTIONS -- GTR(15,15,NT)
C
      CALL WREADW(ITAP30,GTR,INTOWP(15*15*NT),IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '******TRANSFORMATION MATICES FOR G FNS - GTR(I,J,K)*'
      WRITE (6,*) '          '
      DO 235 K=1,NT
      DO 236 J=1,15
      WRITE (6,26) (GTR(I,J,K),I=1,15)
  236 CONTINUE
  235 CONTINUE
C
C     ----- NUMBER OF GENERAL CONTRACTIONS PER SET OF PRIMITIVES -----
C
      WRITE (6,*) '          '
      WRITE (6,*) '*****NGNCN ARRAY************************************'
      WRITE (6,*) '          '
      CALL WREADW (ITAP30,NGNCN,NSHELL,IEND,IEND)
      WRITE (6,237) (NGNCN(J),J=1,NSHELL)
  237 FORMAT(10I5)
C
C     ----- WHETHER PURE SETS OF S,P,D,F,G FUNCTIONS  -  NPERSH  -----
C
      WRITE (6,*) '          '
      WRITE (6,*) '*****NPERSH ARRAY***********************************'
      WRITE (6,*) '          '
      CALL WREADW (ITAP30,NPERSH,5,IEND,IEND)
      WRITE (6,238) (NPERSH(J),J=1,5)
  238 FORMAT(5I5)
C
C     ----- WRITE OUT DIRECTORY FOR FIRST CALCULATION -----
C
      WRITE (6,*) '          '
      WRITE (6,*) '****TITLE OF INDIVIDUAL CALCULATION*****************'
      WRITE (6,*) '          '
      CALL WREADW (ITAP30,TITLE,20,IEND,IEND)
      WRITE (6,21) I,(TITLE(J),J=1,20)
      DO 12 I=1,2
      CALL WREADW (ITAP30,ITEMP,20,IEND,IEND)
      WRITE (6,21) I,(ITEMP(J),J=1,20)
   12 CONTINUE
C
C     ----- ZERO POINTERS TO VECTORS ETC. -----
C
      CALL WREADW (ITAP30,ITEMP,20,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****POINTERS TO VECTORS, EIGENVALUES ETC.***********'
      WRITE (6,*) '          '
      WRITE (6,27) (ITEMP(I),I=1,20)
C
C     ----- WRITE OUT GEOMETRY -----
C
      CALL WREADW(ITAP30,C,INTOWP(3*NAT),IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****GEOMETRY****************************************'
      WRITE (6,*) '          '
      DO 222 I=1,NAT
      WRITE (6,45) (C(J,I),J=1,3)
  222 CONTINUE
   45 FORMAT(3F20.15)
C
C     ----- WRITE OUT NUCLEAR REPULSION ENERGY, SCF ENERGY
C
      CALL WREADW(ITAP30,T,INTOWP(2),IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****NUCLEAR REPULSION ENERGY************************'
      WRITE (6,*) '          '
      WRITE (6,26) T(1)
      WRITE (6,*) '          '
      WRITE (6,*) '****SCF ENERGY**************************************'
      WRITE (6,*) '          '
      WRITE (6,26) T(2)
C
      IF(ITEMP(1).EQ.0) GO TO 13
C
C     ----- WRITE SYMMETRY BLOCKED SCF VECTOR
C
      CALL WREADW(ITAP30,SU,INTOWP(MXCOEF),ITEMP(1),IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****SYMMETRY BLOCKED SCF VECTOR*********************'
      WRITE (6,*) '          '
      WRITE (6,26) (SU(J),J=1,MXCOEF)
C
C     ----- WRITE SCF ORBITAL ENERGIES
C
      CALL WREADW(ITAP30,EIG,INTOWP(MSFRU),IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****SCF ORBITAL EIGENVALUES*************************'
      WRITE (6,*) '          '
      WRITE (6,26) (EIG(J),J=1,MSFRU)
C
C     ----- WRITE LABELS OF IRREDUCIBLE REPS IN SCF
C
      CALL WREADW(ITAP30,ITYP,NSYMHF,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****LABELS OF IRREDUCIBLE REPS IN SCF***************'
      WRITE (6,*) '          '
      WRITE (6,20) (ITYP(J),J=1,NSYMHF)
C
C     ----- WRITE NUMBER OF ORBITALS PER IRREP
C
      CALL WREADW(ITAP30,NLAMDA,NSYMHF,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****NUMBERS OF ORBITALS PER IRREP*******************'
      WRITE (6,*) '          '
      WRITE (6,22) (NLAMDA(J),J=1,NSYMHF)
C
C     ----- WRITE NUMBER CLOSED SHELL ORBITALS PER IRREP
C
      CALL WREADW(ITAP30,NC,NSYMHF,IEND,IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****NUMBER CLOSED SHELL ORBITALS PER IRREP**********'
      WRITE (6,*) '          '
      WRITE (6,22) (NC(J),J=1,NSYMHF)
C
      IF(IOPEN.NE.0) THEN
C
C      ----- WRITE NUMBER OPEN SHELL ORBITALS PER IRREP
C
       CALL WREADW(ITAP30,NO,NSYMHF,IEND,IEND)
       WRITE (6,*) '          '
       WRITE (6,*) '****NUMBER CLOSED SHELL ORBITALS PER IRREP*********'
       WRITE (6,*) '          '
       WRITE (6,22) (NO(J),J=1,NSYMHF)
C
C      ----- WRITE ALPHA ENERGY COEFFICIENTS
C
       CALL WREADW(ITAP30,ALPHA,INTOWP(IABS(IOPEN)),IEND,IEND)
       WRITE (6,*) '          '
       WRITE (6,*) '****ALPHA ENERGY COEFFICIENTS**********************'
       WRITE (6,*) '          '
       WRITE (6,26) (ALPHA(J),J=1,IABS(IOPEN))
C
C      ----- WRITE BETA ENERGY COEFFICIENTS
C
       CALL WREADW(ITAP30,BETA,INTOWP(IABS(IOPEN)),IEND,IEND)
       WRITE (6,*) '          '
       WRITE (6,*) '******BETA ENERGY COEFFICIENTS*********************'
       WRITE (6,*) '          '
       WRITE (6,26) (BETA(J),J=1,IABS(IOPEN))
C
      END IF
C
C      ----- WRITE LAGRANGIAN COEFFICIENTS
C
       CALL WREADW(ITAP30,WMOT,INTOWP(NX),IEND,IEND)
       WRITE (6,*) '          '
       WRITE (6,*) '******LAGRANGIAN COEFFICIENTS**********************'
       WRITE (6,*) '          '
       CALL PRINT(WMOT,NX,MSFRU,6)
C
   13 IF(ITEMP(2).EQ.0) GO TO 14
C
      CALL WREADW(ITAP30,DE,INTOWP(3*NAT),ITEMP(2),IEND)
      WRITE (6,*) '          '
      WRITE (6,*) '****SCF GRADIENTS***********************************'
      WRITE (6,*) '          '
      WRITE (6,40) (DE(I),I=1,3*NAT)
   40 FORMAT (3F20.12)
C
   14 WRITE (6,90) IEND,I2SEC(IEND)
   90 FORMAT (//,' TAPE30 USED ',I8,' WORDS OR ',I5,' SECTORS')
C
      call psdone()
      stop
      RETURN
      END
