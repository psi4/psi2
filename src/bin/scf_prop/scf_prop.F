C**********************************************************************
      subroutine fentry(cc,ia,maxcor)
C**********************************************************************
c     PROGRAM SCF_PROP                                                *
C   THE MULLIKEN POPULATON AND DIPOLE MOMENT PROGRAM                  *
C**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
C**********************************************************************
c   Last updated on September 01, 2005 by Dr. Yukio Yamaguchi         *
C**********************************************************************
c   Last updated on November 15, 2002 by Dr. Yukio Yamaguchi          *
c   Calculate contribution of each occupied MO to the dipole moment   *
c   for the SCF wave function                                         *
C**********************************************************************
c   Last updated on December 20, 1999 by Dr. Yukio Yamaguchi          *
C**********************************************************************
c   Last updated on April 29, 1997 by Dr. Yukio Yamaguchi             *
c   Modification for the second moments                               *
C**********************************************************************
c   Last updated on January 18, 1997 by Dr. Yukio Yamaguchi           *
c   Modification for a power2-ibm-aix4.1 system                       *
c   Index has been changes to Indxx                                   *
C**********************************************************************
c   Last updated on May 02, 1994  by Dr. Yukio Yamaguchi              *
c   Modification for an SCF wavefunction                              *
C**********************************************************************
c   by: Roger Grev, December 20, 1991             search:  crsg       *
c       changes to basset.f and popul.f so it works for 5d7f9g.       *
c       basically changing NBASIS to nbfao.                           *
c                                                                     *
c       by: Edward Seidl (cets072891)                                 *
c     date: July 28, 1991                                             *
c   reason: change dimension of ex,ey and ez in common exyzs and      *
c           contr in basset to allow for f&g functions                *
c                                                                     *
c       by: Edward Seidl (cets030291)                                 *
c     date: March 2, 1991                                             *
c   reason: replace drum with psinit and use new input format         *
c                                                                     *
c       by: Edward Seidl (cets)                                       *
c     date: August 11, 1988                                           *
c   reason: change to fentry to allow dynamic core allocation         *
c**********************************************************************
C*  BY:  RICHARD REMINGTON                         search:  c3-24-88  *
C*  DATE:  MARCH  24,  1988                                           *
C*  REASON: DECREASE CORE TO RUN IN 7MB ON 9370                       *
C**********************************************************************
C*  LAST UPDATED ON MARCH 20, 1985 BY TJL                             *
C**********************************************************************
C*  LAST UPDATED ON MAY 27, 1985  BY RICHARD REMINGTON                *
C*   CALL TO RFILE(ITAP30) WAS REPLACED BY SREW AND THE               *
C*   RFILE WAS MOVED TO AFTER ITAP30=30, ALSO DEFAULTS                *
C*   WERE ADDED IN CASE NO INPUT TO PROPER IS FOUND.                  *
C**********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include <error.h>
      integer frdc, frdi, frdf, fwkclr, fwkadd, frdboo, frdpos
      integer fcnt, frdia1
      integer errcod
      character*6 refctr
C
      dimension cc(maxcor),ia(maxcor*2)
      DIMENSION I30(200),A30(100)
      COMMON/BASIS/NBASIS,NBFAO,NBFSO,NBATRI,NTRI
      common/calif/lpara(1024),apara(1024)
      COMMON/COORD/CHGS(112),X(3,112),W(112)
      common/dipct/icent
      common/dipnl/dnx,dny,dnz
      COMMON/FLAGS/MPOINT,MCONST
      COMMON/FUNCS/NATOM,N3N,NSYM
      COMMON/GAUSS/ZS(200),CS(200)
      COMMON/NORBS/NAST(112),NAED(112)
      COMMON/NSHEL/KNUC(200),KTYPE(200),KPRIM(200),KLOC(200),
     1 KMIN(200),KMAX(200),KNG(200),NSHELL,NPRIM
      COMMON/PARAS/MXCOEF,LOCVEC
      COMMON/QVALS/QORB(256),QATM(112),QNET(112)
      common/saven/escf
      common/savlb/label(20)
      COMMON/SIGNS/IOFF(500),IPRNT
      common/symmt/mosym(256)
      COMMON/CI101/IOPEN,IOCC,JOCC,KOCC
      EQUIVALENCE (I30,A30)
    1 FORMAT(//,2X,' THE MULLIKEN POPULATION AND DIPOLE MOMENT PROGRAM',
     1        /,2X,'                F & G function version            ',
     2        /,2X,'                (Version of 09/01/2005)           ')
    2 FORMAT(3I5)
    3 FORMAT(/,2X,' PARAMETERS'/
     * 2X,' ISCFCI = ',I8/
     * 2X,' ICENT  = ',I8/
     * 2X,' IPRNT  = ',I8/
     * 2X,' NBASIS = ',I8/
     * 2X,' NBFAO  = ',I8/
     * 2X,' NBFSO  = ',I8/
     * 2X,' NBATRI = ',I8/
     * 2X,' NTRI   = ',I8/
     * 2X,' NATOM  = ',I8/
     * 2X,' N3N    = ',I8/
     * 2X,' NSYM   = ',I8/
     * 2X,' IOPEN  = ',I8/
     * 2X,' NSHELL = ',I8/
     * 2X,' NPRIM  = ',I8/)
    4 FORMAT(/,3X,20A4)
    5 FORMAT(//,2X,' ENUC   = ',F20.10/
     1          2X,' ESCF   = ',F20.10/)
    6 FORMAT(//,2X,' REQUIRED MEMORY EXCEEDS MAXCOR'/
     1          2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)
    7 FORMAT(//,2X,' WARNING! NO INPUT TO PROPER HAS BEEN FOUND. '/
     1          2X,'          DEFAULT PARAMETERS WILL BE USED. ',/)
C
      call psinit('APPEND')
      CALL TSTART(6)
c
      ITAPE3=3
      INPUT=5
      ITAPE6=6
      ITAP30=30
c
      CALL RFILE(ITAP30)
C
      WRITE(6,1)
      WRITE(3,1)
      write(6,*) '  '
      write(6,*) '  maxcor = ',maxcor
C
C   ISCFCI IS A PARAMETER FOR TYPE OF WAVE FUNCTION
C     ISCFCI = 0 OR 1  FOR SCF
C     ISCFCI = 2       FOR CI
C     ISCFCI = 3       FOR MCSCF
C   ICENT IS A PARAMETER FOR REFERENCE COORDINATE OF DIPOLE MOMENT
C     ICENT  = 1       FOR CENTER OF MASS
C     ICENT  = 2       FOR ORIGIN OF SPACE FIXED COORDINATE
C   IPRINT IS A PARAMETER FOR SIZE OF OUTPUT
C
C   SET UP DEFAULT VALUES (SCF) IF INPUT TO PROPER IS NOT FOUND
C
      ISCFCI = 1
C
      errcod = fwkclr()
      errcod = fwkadd(5,'SCF_PROP')
      errcod = fwkadd(5,'DEFAULT')
      errcod = frdc('CENTER',refctr)
      icent = 1 
      if(refctr.eq.'CMS   ') then
        icent = 1
      end if
      if(refctr.eq.'ORIGIN') then
        icent = 2
      end if
      iprnt = 0
      errcod = frdi('PRINT',iprnt)
C
      IOFF(1)=0
      DO 101 I=1,499
      IOFF(I+1)=IOFF(I)+I
  101 continue
C
      CALL SREW(ITAP30)
      CALL WREADW(ITAP30,I30,200,101,JUNK)
      MPOINT=I30(2)
      MCONST=I30(3)
      NCALCS=I30(5)
      NBFSO=I30(18)
      NATOM=I30(19)
      NBFAO=I30(22)
      NSHELL=I30(27)
      NPRIM=I30(32)
      NSYM=I30(41)
      MXCOEF=I30(42)
      IOPEN=I30(43)
      NBASIS=NBFSO
      N3N=NATOM*3
      NBATRI=IOFF(NBFAO+1)
      NTRI=IOFF(NBASIS+1)
C
      WRITE(6,3) ISCFCI,ICENT,IPRNT,NBASIS,NBFAO,NBFSO,NBATRI,NTRI,
     1           NATOM,N3N,NSYM,IOPEN,NSHELL,NPRIM
C
C   READ IN LABELS
      JUNK=101+MCONST+MPOINT+NCALCS-1
      CALL WREADW(ITAP30,LOCCAL,1,JUNK,JUNK)
      CALL WREADW(ITAP30,I30,60,LOCCAL,LOCCAL)
      WRITE(6,4) (I30(I),I=1,20)
      do 102 i=1,20
      label(i)=i30(i)
  102 continue
      CALL WREADW(ITAP30,I30,20,LOCCAL,JUNK)
      LOCVEC=I30(1)
C
C   READ IN GEOMETRY
      LOCCAL=LOCCAL+20
      CALL WREADW(ITAP30,X,NATOM*6,LOCCAL,LOCCAL)
C
C   READ IN ENUC AND ESCF
      CALL WREADW(ITAP30,A30,20,LOCCAL,LOCCAL)
      ENUC=A30(1)
      ESCF=A30(2)
      write(6,5) enuc,escf
c
c*******************************************************
c   Dynamic allocation
c*******************************************************
c     IC1 : Uas(nbfao,nbasis)
c     IC2 : Eig(nbasis)
c     IC3 : Occ(nbasis)
c     IC4 : Eso(nbfso,nbasis)
c     IC5 : Eao(nbfao,nbasis)
c     IC6 : Pao(nbatri)
c*******************************************************
      IC1=1
      IC2=IC1+NBFAO*NBFSO
      IC3=IC2+NBASIS
      IC4=IC3+NBASIS
      IC5=IC4+NBFSO*NBFSO
      IC6=IC5+NBFAO*NBFSO
      ICSAV=IC6+nbatri
c*******************************************************
C
C   READ IN BASIS SET INFORMATION
CC    WRITE(6,21)
      WRITE(3,21)
   21 FORMAT(//,2X,' NOW YOU ARE IN BASSET'/)
C.................UAS.....
      CALL BASSET(CC(IC1))
C
c=============================================================
c   Read in eigenvectors for an SCF wavefunction
c=============================================================
CC    WRITE(6,22)
      WRITE(3,22)
   22 FORMAT(//,2X,'NOW YOU ARE IN EIGVEC'/)
      ic7=icsav
      icmax=ic7+nbfao*nbfso
      IF(ICMAX.GT.MAXCOR) GO TO 399
C................UAS     EIG     OCC     ESO     EAO     EE......
      CALL EIVEC(CC(IC1),CC(IC2),CC(IC3),CC(IC4),CC(IC5),CC(IC7))
C
c=============================================================
c   Form a density matrix for an SCF wavefunction
c=============================================================
CC    WRITE(6,23)
      WRITE(3,23)
   23 FORMAT(//,2X,' NOW YOU ARE IN SCFDEN'/)
      icmax=icsav
      IF(ICMAX.GT.MAXCOR) GO TO 399
C.................OCC     EAO     P.......
      CALL SCFDEN(CC(IC3),CC(IC5),CC(IC6))
c
C   CALCULATE ONE ELECTRON INTEGRALS FOR A TEST
CC    WRITE(6,24)
      WRITE(3,24)
   24 FORMAT(//,2X,' NOW YOU ARE IN ONEHGT'/)
      IC7=icsav
      IC8=IC7+NBATRI
      IC9=IC8+NBATRI
      IC10=IC9+NBATRI
      ICMAX=IC10+NBATRI
cc    write(6,* ) ' ic1 = ',ic1
cc    write(6,* ) ' ic2 = ',ic2
cc    write(6,* ) ' ic3 = ',ic3
cc    write(6,* ) ' ic4 = ',ic4
cc    write(6,* ) ' ic5 = ',ic5
cc    write(6,* ) ' ic6 = ',ic6
cc    write(6,* ) ' ic7 = ',ic7
cc    write(6,* ) ' ic8 = ',ic8
cc    write(6,* ) ' ic9 = ',ic9
cc    write(6,* ) ' ic10= ',ic10
      IF(ICMAX.GT.MAXCOR) GO TO 399
C.................P       SS      TT      VV      HH.......
C
      CALL ONEHGT(CC(IC6),CC(IC7),CC(IC8),CC(IC9),CC(IC10))
C
CC    WRITE(6,25)
      WRITE(3,25)
   25 FORMAT(//,2X,' NOW YOU ARE IN ONEHGT'/)
      IC7=icsav
      IC8=IC7+NBATRI
      IC9=IC8+NBATRI
      IC10=IC9+NBATRI
      IC11=IC10+NBATRI
      IC12=IC11+NTRI
      IC13=IC12+NBFAO*NBFAO
      ICMAX=IC13+NBFAO*NBFAO
C................OCC     EAO     HH       HMO      U        T........
      CALL ENONE(CC(IC3),CC(IC5),CC(IC10),CC(IC11),CC(IC12),CC(IC13))
C
C   CALCULATE MULLIKEN POPULATION
CC    WRITE(6,26)
      WRITE(3,26)
   26 FORMAT(//,2X,' NOW YOU ARE IN POPUL'/)
      IC7=icsav
      IC8=IC7+NBATRI
      IC9=IC8+NBATRI
      ICMAX=IC9+NATOM*NATOM
      IF(ICMAX.GT.MAXCOR) GO TO 399
C.................P       S       AOPOP   BATM....
      CALL  POPUL(CC(IC6),CC(IC7),CC(IC8),CC(IC9))
C
C   FIND A REFERENCE CENTER IF NECESSARY
CC    WRITE(6,27)
      WRITE(3,27)
   27 FORMAT(//,2X,' NOW YOU ARE IN CENTER'/)
      CALL CENTER
C
C   CALCULATE DIPOLE MOMENT IN THE AO BASIS
CC    WRITE(6,28)
      WRITE(3,28)
   28 FORMAT(//,2X,' NOW YOU ARE IN DIPAO'/)
      IC7=icsav
      IC8=IC7+NBATRI
      IC9=IC8+NBATRI
      IC10=IC9+NBATRI
      ICMAX=IC10+NBATRI
      IF(ICMAX.GT.MAXCOR) GO TO 399
C................P       DIPX    DIPY    DIPZ.....
      CALL DIPAO(CC(IC6),CC(IC8),CC(IC9),CC(IC10))
C
C   CALCULATE DIPOLE MOMENT IN THE MO BASIS
CC    WRITE(6,29)
      WRITE(3,29)
   29 FORMAT(//,2X,' NOW YOU ARE IN DIPMO'/)
      IC7=icsav
      IC8=IC7+NBATRI
      IC9=IC8+NBATRI
      IC10=IC9+NBATRI
      IC11=IC10+NBATRI
      IC12=IC11+NTRI
      IC13=IC12+NTRI
      IC14=IC13+NTRI
      IC15=IC14+NBFAO*NBFAO
      ICMAX=IC15+NBFAO*NBFAO
      IF(ICMAX.GT.MAXCOR) GO TO 399
C................OCC     EAO     DIPAX   DIPAY   DIPAZ....
      CALL DIPMO(CC(IC3),CC(IC5),CC(IC8),CC(IC9),CC(IC10),
C................DIMX     DIPMY    DIPMZ    U        T........
     *           CC(IC11),CC(IC12),CC(IC13),CC(IC14),CC(IC15))
      go to 400
C
  399 continue
      WRITE(6,6) ICMAX,MAXCOR
  400 CONTINUE
      CALL TSTOP(6)
C
      STOP
      END
      SUBROUTINE BASSET(UAS)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer nbfao,nbfso
      DIMENSION UAS(NBFAO,NBFSO)
      DIMENSION CONTR(1000),I30(200)
      COMMON/BASIS/NBASIS,NBFAO,NBFSO,NBATRI,NTRI
      COMMON/COORD/CHGS(112),XX(3,112),W(112)
      COMMON/FLAGS/MPOINT,MCONST
      COMMON/FUNCS/NATOM,N3N,NSYM
      COMMON/GAUSS/ZS(200),CS(200)
      COMMON/NORBS/NAST(112),NAED(112)
      COMMON/NSHEL/KNUC(200),KTYPE(200),KPRIM(200),KLOC(200),
     1 KMIN(200),KMAX(200),KNG(200),NSHELL,NPRIM
      COMMON/SIGNS/IPQ(500),IPRNT
      DATA a00 / 0.0D+00 /
    1 FORMAT(//,2X,' NPRIM  IS ',I5)
    2 FORMAT(2X,I5,2F20.10)
    3 FORMAT(//,2X,' NSHELL IS ',I5)
    4 FORMAT(2X,8I5)
    5 FORMAT(//,2X,' NATOM IS ',I5)
    6 FORMAT(2X,I5,4F20.10,5X,2I5)
C
      ITAP30=30
      JUNK=101+MCONST
      CALL WREADW(ITAP30,I30,MPOINT,JUNK,JUNK)
C
C   READ IN BASIS SET
      CALL WREADW(ITAP30,CHGS,NATOM*2,I30(1),JUNK)
      CALL WREADW(ITAP30,ZS,NPRIM*2,I30(5),JUNK)
      call wreadw(itap30,contr,intowp(nprim)*5,i30(6),junk)
      CALL WREADW(ITAP30,KPRIM,NSHELL,I30(7),JUNK)
      CALL WREADW(ITAP30,KNUC,NSHELL,I30(8),JUNK)
      CALL WREADW(ITAP30,KTYPE,NSHELL,I30(9),JUNK)
      CALL WREADW(ITAP30,KNG,NSHELL,I30(10),JUNK)
      CALL WREADW(ITAP30,KLOC,NSHELL,I30(11),JUNK)
      CALL WREADW(ITAP30,KMIN,NSHELL,I30(12),JUNK)
      CALL WREADW(ITAP30,KMAX,NSHELL,I30(13),JUNK)
      CALL WREADW(ITAP30,UAS,NBFAO*NBFSO*2,I30(29),JUNK)
C
      IJ=0
cets072891      DO 102 I=1,3
      do 102 i=1,5
        DO 101 J=1,NPRIM
          IJ=IJ+1
          VALU=CONTR(IJ)
          IF(VALU.EQ.a00) GO TO 101
          CS(J)=VALU
  101   CONTINUE
  102 CONTINUE
C
      IF(IPRNT.LE.3) GO TO 201
      WRITE(6,1) NPRIM
      DO 103 I=1,NPRIM
      WRITE(6,2) I,ZS(I),CS(I)
  103 continue
      WRITE(6,3) NSHELL
      DO 104 I=1,NSHELL
      WRITE(6,4) I,KNUC(I),KTYPE(I),KPRIM(I),KLOC(I),KMIN(I),KMAX(I),
     * KNG(I)
  104 continue
C
C   FORM ARRAIES FOR AO'S
  201 CONTINUE
      NAST(1)=1
      IATOM=1
      DO 105 I=1,NSHELL
      II=KNUC(I)
      IF(II.EQ.IATOM) GO TO 105
      IATOM=IATOM+1
      NAED(IATOM-1)=KLOC(I)-1
      NAST(IATOM)=KLOC(I)
  105 CONTINUE
c
crsg change nbasis to nbfao for 5d7f9g
crsg  NAED(NATOM)=NBASIS
      NAED(NATOM)=nbfao
      IF(IPRNT.LE.3) GO TO 202
      WRITE(6,5) NATOM
      DO 106 I=1,NATOM
      WRITE(6,6) I,XX(1,I),XX(2,I),XX(3,I),CHGS(I),NAST(I),NAED(I)
  106 continue
  202 CONTINUE
C
      RETURN
      END
      SUBROUTINE EIVEC(UAS,EIG,OCC,ESO,EAO,EE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION UAS(NBFAO,NBFSO),EIG(NBASIS),OCC(NBASIS)
      DIMENSION ESO(NBFSO,NBFSO),EAO(NBFAO,NBFSO)
cets080989 the mips compiler can't handle the next line
cets080989      DIMENSION EE(NBFAO*NBFSO)
      dimension ee(1)
      DIMENSION NLAMDA(10),NCLOSD(10),NOPEN(10),OCT(2)
      DIMENSION ALPT(15),BETT(15)
      COMMON/BASIS/NBASIS,NBFAO,NBFSO,NBATRI,NTRI
      COMMON/FUNCS/NATOM,N3N,NSYM
      COMMON/SIGNS/IPQ(500),IPRNT
      common/symmt/mosym(256)
      COMMON/PARAS/MXCOEF,LOCVEC
      COMMON/CI101/IOPEN,IOCX,JOCX,KOCX
      DATA a00,ONE,TWO / 0.0D+00 , 1.0D+00 , 2.0D+00 /
    1 FORMAT(//,2X,' UAS MATRIX'/)
    2 FORMAT(//,2X,' NOPEN = ',10I5/)
    3 FORMAT(//,2X,' ESO MATRIX'/)
    4 FORMAT(//,2X,' EAO MATRIX'/)
C
      ITAP30=30
      CALL WREADW(ITAP30,EE,MXCOEF*2,LOCVEC,LOCVEC)
      CALL WREADW(ITAP30,EIG,NBFSO*2,LOCVEC,LOCVEC)
      IF(IPRNT.LE.3) GO TO 201
      WRITE(6,1)
      CALL MATOUT(UAS,NBFAO,NBFSO,NBFAO,NBFSO,6)
C
C   READ IN SO-MO EIGENVECTORS
  201 CONTINUE
      LOCVEC=LOCVEC+NSYM
      CALL WREADW(ITAP30,NLAMDA,NSYM,LOCVEC,LOCVEC)
      CALL WREADW(ITAP30,NCLOSD,NSYM,LOCVEC,LOCVEC)
      WRITE(3,28) (NLAMDA(I),I=1,NSYM)
   28 FORMAT(2X,' NLAMDA = ',10I5/)
      WRITE(3,29) (NCLOSD(I),I=1,NSYM)
   29 FORMAT(2X,' NCLOSD = ',10I5/)
C
C   READ IN COUPLING CONSTANTS FOR OPEN SHELLS
      WRITE(3,30) IOPEN
      IF(IPRNT.GT.2)
     *WRITE(6,30) IOPEN
   30 FORMAT(2X,' IOPEN = ',I5/)
      IF(IOPEN.NE.0) THEN
      KOPEN=IABS(IOPEN)
      CALL WREADW(ITAP30,NOPEN,NSYM,LOCVEC,LOCVEC)
      WRITE(3,2) (NOPEN(I),I=1,NSYM)
      IF(IPRNT.GT.2)
     *WRITE(6,2) (NOPEN(I),I=1,NSYM)
      CALL WREADW(ITAP30,ALPT,KOPEN*2,LOCVEC,LOCVEC)
      WRITE(3,31) (ALPT(I),I=1,KOPEN)
      IF(IPRNT.GT.2)
     *WRITE(6,31) (ALPT(I),I=1,KOPEN)
   31 FORMAT(2X,' ALPT = ',10F10.5/)
      CALL WREADW(ITAP30,BETT,KOPEN*2,LOCVEC,LOCVEC)
      WRITE(3,32) (BETT(I),I=1,KOPEN)
      IF(IPRNT.GT.2)
     *WRITE(6,32) (BETT(I),I=1,KOPEN)
   32 FORMAT(2X,' BETT = ',10F10.5/)
      END IF
      IF(IOPEN.GE.0) GO TO 202
      OCT(1)=TWO/(ONE-ALPT(1))
      OCT(2)=TWO/(ONE-ALPT(3))
C
C   FORM OCCUPATION VECTOR
  202 CONTINUE
      DO 101 I=1,NBASIS
      OCC(I)=a00
  101 continue
      ITC=0
      I=0
      DO 105 ISM=1,NSYM
      DO 103 J=I+1,I+NCLOSD(ISM)
      OCC(J)=TWO
  103 CONTINUE
      K=I+NCLOSD(ISM)
      IF(IOPEN.LT.0) GO TO 203
      IF(IOPEN.EQ.0) GO TO 204
      DO 104 J=K+1,K+NOPEN(ISM)
      OCC(J)=ONE
  104 CONTINUE
      GO TO 204
  203 CONTINUE
      ITC=ITC+1
      OCC(K+1)=OCT(ITC)
  204 CONTINUE
      I=I+NLAMDA(ISM)
  105 CONTINUE
C
C   FORM SO-MO EIGENVECTOR MATRIX
      DO 106 I=1,NBFSO
      DO 106 J=1,NBFSO
      ESO(I,J)=a00
  106 CONTINUE
      I=0
      J=0
      DO 108 ISYM=1,NSYM
      N=NLAMDA(ISYM)
      DO 107 K=J+1,J+N
      DO 107 L=J+1,J+N
      I=I+1
      ESO(L,K)=EE(I)
  107 CONTINUE
      J=J+N
  108 CONTINUE
C
      IF(IPRNT.LE.3) GO TO 205
      WRITE(6,3)
      CALL EIGOUT(ESO,EIG,OCC,NBFSO,NBFSO,NBASIS,NBASIS,6)
C
  205 CONTINUE
      call izero(mosym,nbasis)
      II=0
      DO 110 ISYM=1,NSYM
      N=NLAMDA(ISYM)
      DO 109 I=1,N
      II=II+1
      MOSYM(II)=ISYM
  109 CONTINUE
  110 CONTINUE
      WRITE(3,33) (MOSYM(I),I=1,NBASIS)
   33 FORMAT(2X,' MOSYM = ',10I5/)
C
C   TRANSFORM EIGENVECTORS FROM SO TO AO BASIS
      DO 112 I=1,NBFAO
      DO 112 J=1,NBASIS
      VALU=a00
      DO 111 K=1,NBFSO
      VALU=VALU+UAS(I,K)*ESO(K,J)
  111 CONTINUE
      EAO(I,J)=VALU
  112 CONTINUE
      IF(IPRNT.LE.3) GO TO 206
      WRITE(6,4)
      CALL EIGOUT(EAO,EIG,OCC,NBFAO,NBFAO,NBASIS,NBASIS,6)
  206 CONTINUE
C
      itape8 = 8
      call ffile(itape8,'file8',0)
c     rewind itape8
      iuas = 0
      do 115 i=1,nbfao
      do 115 j=1,nbasis
      valu = uas(i,j)
      if(dabs(valu).lt.1.0d-6) go to 115
      iuas = iuas + 1
      write(itape8,20) i,j,valu
      if(iprnt.ge.5) then
        write(6,20) i,j,valu
      end if
  115 continue
   20 format(2i5,f20.10)
      if(iprnt.ge.5) then
        write(6,*) ' iuas = ',iuas
      end if
C
      RETURN
      END
      SUBROUTINE SCFDEN(OCC,EAO,P)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer nbasis,nbfao,nbfso,nbatri
      DIMENSION OCC(NBASIS),EAO(NBFAO,NBFSO),P(NBATRI)
      COMMON/BASIS/NBASIS,NBFAO,NBFSO,NBATRI,NTRI
      COMMON/SIGNS/IOFF(500),IPRNT
      DATA a00 / 0.0D+00 /
    1 FORMAT(//,2X,' THE SCF DENSITY MATRIX'/)
C
C   FORM THE DENSITY MATRIX
      DO 101 I=1,NBATRI
      P(I)=a00
  101 continue
      IJ=0
      DO 103 I=1,NBFAO
      DO 103 J=1,I
      IJ=IJ+1
      VALU=a00
      DO 102 K=1,NBASIS
      IF(OCC(K).LE.a00) GO TO 102
      VALU=VALU+EAO(I,K)*EAO(J,K)*OCC(K)
  102 CONTINUE
      P(IJ)=VALU
  103 CONTINUE
      IF(IPRNT.LE.3) GO TO 201
      WRITE(6,1)
      CALL PRINT(P,NBATRI,NBFAO,6)
C
  201 CONTINUE
      RETURN
      END
      SUBROUTINE ONEHGT(P,SS,TT,VV,HH)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL IANDJ,EQUAL
      integer nbatri
      DIMENSION P(NBATRI),SS(NBATRI),TT(NBATRI),VV(NBATRI),HH(NBATRI)
      DIMENSION VX0(36),VY0(36),VZ0(36)
      COMMON/BASIS/NBASIS,NBFAO,NBFSO,NBATRI,NTRI
      COMMON/COORD/CHGS(112),C(3,112),W(112)
cets072891      COMMON/EXYZS/EX(7,7,13),EY(7,7,13),EZ(7,7,13)
      common/exyzs/ex(7,7,14),ey(7,7,14),ez(7,7,14)
      COMMON/FUNCS/NATOM,N3N,NSYM
      COMMON/GAUSS/ZS(200),CS(200)
      COMMON/IJVEC/IX(84),IY(84),IZ(84),JX(84),JY(84),JZ(84)
      COMMON/INDXX/IJX(225),IJY(225),IJZ(225),IJPOS(225)
      COMMON/NSHEL/KNUC(200),KTYPE(200),KSTART(200),KLOC(200),
     1 KMIN(200),KMAX(200),KNG(200),NSHELL,NPRIM
      COMMON/ROOT /XX,U(9),WT(9),NROOTS
      COMMON/SIGNS/IPQ(500),IPRNT
      COMMON/TEMP0/SX0(36),SY0(36),SZ0(36),TX0(36),TY0(36),TZ0(36)
      EQUIVALENCE (SX0,VX0),(SY0,VY0),(SZ0,VZ0)
      DATA ITOL / 10 /
      DATA PI212 / 1.1283791670955D+00 /
      DATA PI32 / 5.56832799683170D+00 /
      DATA a00,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /
      DATA NQLIM / 7 /
    1 FORMAT(//,2X,' S MATRIX'/)
    2 FORMAT(//,2X,' T MATRIX'/)
    3 FORMAT(//,2X,' V MATRIX'/)
    4 FORMAT(//,2X,' H MATRIX'/)
    5 FORMAT(/,2X,' EONE (calculated in the AO basis) = ',F20.10/)
C
      TOL=ITOL*2.30258D+00
C
C     ZERO ALL THE CORE THAT WE NEED FOR S AND ONE ELECTRON INTEGRALS
C
cc    write(6,*) ' in onehgt'
cc    write(6,*) ' nbatri = ',nbatri
      DO 101 I=1,NBATRI
      SS(I)=a00
      TT(I)=a00
      VV(I)=a00
  101 CONTINUE
C
C
C:::::::::::::::::::::::::::::
C:::::----- I SHELL -----:::::
C:::::::::::::::::::::::::::::
C
      DO 6000 ISH=1,NSHELL
C
      IATOM=KNUC(ISH)
      AX=C(1,IATOM)
      AY=C(2,IATOM)
      AZ=C(3,IATOM)
      IGMIN=KSTART(ISH)
      IGMAX=IGMIN+KNG(ISH)-1
      LIT=KTYPE(ISH)
      MINI=KMIN(ISH)
      MAXI=KMAX(ISH)
      LOCI=KLOC(ISH)
C
C:::::::::::::::::::::::::::::
C:::::----- J SHELL -----:::::
C:::::::::::::::::::::::::::::
C
      DO 5000 JSH=1,ISH
C
      JATOM=KNUC(JSH)
      BX=C(1,JATOM)
      BY=C(2,JATOM)
      BZ=C(3,JATOM)
      JGMIN=KSTART(JSH)
      JGMAX=JGMIN+KNG(JSH)-1
      LJT=KTYPE(JSH)
      MINJ=KMIN(JSH)
      MAXJ=KMAX(JSH)
      LOCJ=KLOC(JSH)
      IANDJ=ISH.EQ.JSH
      EQUAL=IATOM.EQ.JATOM
C
C     ----- WORK OUT INDICES FOR COMBINING 2D INTEGRALS -----
C
      MIJ=0
      JMAX=MAXJ
      II=-1
      DO 103 I=MINI,MAXI
      II=II+1
      JJ=-1
      NX=IX(I)
      NY=IY(I)
      NZ=IZ(I)
      IF(IANDJ) JMAX=I
      DO 102 J=MINJ,JMAX
      MIJ=MIJ+1
      JJ=JJ+1
      IJX(MIJ)=NX+JX(J)
      IJY(MIJ)=NY+JY(J)
      IJZ(MIJ)=NZ+JZ(J)
      IDXI=LOCI+II
      IDXJ=LOCJ+JJ
      IDXS=IPQ(MAX0(IDXI,IDXJ))+MIN0(IDXI,IDXJ)
      IJPOS(MIJ)=IDXS
  102 CONTINUE
  103 CONTINUE
C
      XAB=AX-BX
      YAB=AY-BY
      ZAB=AZ-BZ
      RAB=XAB*XAB+YAB*YAB+ZAB*ZAB
CIBM  ZEROTH ORDER
      NROOTS=(LIT+LJT-2)/2+1
CIBM  ZEROTH ORDER
      MAXIJ=MAX0(LIT,LJT+2)
C
      DO 4000 IG=IGMIN,IGMAX
      CSI=CS(IG)*PI32
      AI=ZS(IG)
      AXI=AI*AX
      AYI=AI*AY
      AZI=AI*AZ
C
      IF(IANDJ) JGMAX=IG
      DO 3000 JG=JGMIN,JGMAX
      CSJ=CS(JG)
      BJ=ZS(JG)
      BXJ=BJ*BX
      BYJ=BJ*BY
      BZJ=BJ*BZ
C
      PP=AI+BJ
      TPP=ONE/PP
      TEMP=(AI*BJ*RAB)*TPP
C
C     ----- IGNORE INTEGRAL IF PRE-EXPONENTIAL TOO SMALL -----
C
      IF(TEMP.GT.TOL) GO TO 3000
C
C     ----- DENSITY FACTOR -----
C
      TP=DSQRT(TPP)
      DIJ=CSI*CSJ*DEXP(-TEMP)*TPP
      IF(IANDJ.AND.IG.NE.JG) DIJ=DIJ+DIJ
      PFAC=DIJ*TP
      PFACH=PFAC*HALF
      DIJPI=DIJ*PI212
      PTWO=TPP*HALF
C
C     FORM SOME PRODUCTS NEEDED IN INTEGRAL EVALUATION
C
      PX=(AXI+BXJ)*TPP
      PY=(AYI+BYJ)*TPP
      PZ=(AZI+BZJ)*TPP
      IF(EQUAL) GO TO 201
      PAX=PX-AX
      PAY=PY-AY
      PAZ=PZ-AZ
      PBX=PX-BX
      PBY=PY-BY
      PBZ=PZ-BZ
      GO TO 202
  201 CONTINUE
      PAX=a00
      PAY=a00
      PAZ=a00
      PBX=a00
      PBY=a00
      PBZ=a00
C
  202 CONTINUE
      CALL ECAL(PAX,PAY,PAZ,PBX,PBY,PBZ,PTWO,MAXIJ)
C
      BJ2=BJ+BJ
      B23=BJ2*BJ2
      IN=-NQLIM
      DO 105 NI=1,LIT
      IN=IN+NQLIM
      DO 104 NJ=1,LJT
      JN=IN+NJ
      XNJ=FLOAT(NJ)
      B22=BJ2*(XNJ+XNJ-ONE)
C
C   NORMAL OVERLAP INTEGRALS
      SX0(JN)=EX(NI,NJ,1)
      SY0(JN)=EY(NI,NJ,1)
      SZ0(JN)=EZ(NI,NJ,1)
C
C   NORMAL KINETIC INTEGRALS
      TX00=EX(NI,NJ+2,1)*B23
      TY00=EY(NI,NJ+2,1)*B23
      TZ00=EZ(NI,NJ+2,1)*B23
      IF(NJ.LE.2) GO TO 203
      B21=(XNJ-ONE)*(XNJ-TWO)
      TX00=TX00+EX(NI,NJ-2,1)*B21
      TY00=TY00+EY(NI,NJ-2,1)*B21
      TZ00=TZ00+EZ(NI,NJ-2,1)*B21
  203 CONTINUE
      TX00=TX00-EX(NI,NJ,1)*B22
      TY00=TY00-EY(NI,NJ,1)*B22
      TZ00=TZ00-EZ(NI,NJ,1)*B22
      TX0(JN)=TX00
      TY0(JN)=TY00
      TZ0(JN)=TZ00
  104 CONTINUE
  105 CONTINUE
C
      DO 106 I=1,MIJ
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      X=SX0(NX)
      Y=SY0(NY)
      Z=SZ0(NZ)
      IJ=IJPOS(I)
      SS(IJ)=SS(IJ)+X*Y*Z*PFAC
      TT(IJ)=TT(IJ)-(TX0(NX)*Y*Z+TY0(NY)*X*Z+TZ0(NZ)*X*Y)*PFACH
  106 CONTINUE
C
C     ----- LOOP OVER ATOMS IN POTENTIAL TERM -----
C
      DO 2000 IC=1,NATOM
      ZNUC=-CHGS(IC)*DIJPI
      CX=C(1,IC)
      CY=C(2,IC)
      CZ=C(3,IC)
      PCX=PX-CX
      PCY=PY-CY
      PCZ=PZ-CZ
      XX=PP*(PCX*PCX+PCY*PCY+PCZ*PCZ)
      IF(NROOTS.LE.3) THEN
      CALL RT123
      ELSE IF(NROOTS.EQ.4) THEN
      CALL ROOT4
      ELSE IF(NROOTS.EQ.5) THEN
      CALL ROOT5
      ELSE IF(NROOTS.LE.9) THEN
      CALL DROOT
      ELSE
      STOP ' NOT ENOUGH ROOTS IN ONEINT'
      END IF
C
      DO 1000 IROOT=1,NROOTS
      ROOT2=U(IROOT)/(ONE+U(IROOT))
      WW=WT(IROOT)*ZNUC
      AROOT2=-ROOT2*(PP+PP)
C
      IN=-NQLIM
      DO 110 NI=1,LIT
      IN=IN+NQLIM
      DO 109 NJ=1,LJT
      MAXT=NI+NJ-1
      JN=IN+NJ
      CALL GCAL(PCX,PCY,PCZ,VXX,VYY,VZZ,AROOT2,NI,NJ,MAXT)
      VX0(JN)=VXX
      VY0(JN)=VYY
      VZ0(JN)=VZZ
  109 CONTINUE
  110 CONTINUE
C
      DO 112 I=1,MIJ
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      X=VX0(NX)
      Y=VY0(NY)
      Z=VZ0(NZ)
      IJ=IJPOS(I)
      VV(IJ)=VV(IJ)+X*Y*Z*WW
  112 CONTINUE
C
 1000 CONTINUE
 2000 CONTINUE
C
 3000 CONTINUE
 4000 CONTINUE
 5000 CONTINUE
 6000 CONTINUE
C
      DO 115 I=1,NBATRI
      HH(I)=TT(I)+VV(I)
  115 continue
C
      IF(IPRNT.LE.2) GO TO 210
      WRITE(6,1)
      CALL PRINT(SS,NBATRI,NBFAO,6)
      WRITE(6,2)
      CALL PRINT(TT,NBATRI,NBFAO,6)
      WRITE(6,3)
      CALL PRINT(VV,NBATRI,NBFAO,6)
      WRITE(6,4)
      CALL PRINT(HH,NBATRI,NBFAO,6)
  210 CONTINUE
C
C   CALCULATE ONE ELECTRON ENERGY FOR A TEST
      EONE=a00
      IJ=0
      DO 120 I=1,NBFAO
      DO 120 J=1,I
      IJ=IJ+1
      FAC=TWO
      IF(I.EQ.J) FAC=ONE
      EONE=EONE+HH(IJ)*P(IJ)*FAC
  120 CONTINUE
      WRITE(6,5) EONE
C
      RETURN
      END
      SUBROUTINE ENONE(OCC,EAO,HAO,HMO,U,T)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION OCC(NBASIS)
      DIMENSION HAO(NBATRI),HMO(NTRI)
      DIMENSION EAO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)
      COMMON/BASIS/NBASIS,NBFAO,NBFSO,NBATRI,NTRI
      COMMON/SIGNS/IOFF(500),IPRNT
      DATA a00 / 0.0D+00 /                                              CPC05240
    1 FORMAT(/,2X,' HMO Matrix')
    2 FORMAT(/,2X,' EONE (calculated in the MO basis) = ',F20.10/)
C                                                                       CPC05250
C   TRANSFORM INTEGRALS FROM AO TO MO BASIS                             CPC05260
      call moconv(eao,hao,hmo,u,t)
      if(iprnt.ge.3) then
        write(6,1)
        call print(hmo,ntri,nbasis,6)
      end if
C
C   CALCULATE ONE ELECTRON ENERGY FOR A TEST
      EONE=a00
      do 101 i=1,nbasis
      ii=ioff(i+1)
      eone = eone + hmo(ii)*occ(i)
  101 continue
      write(6,2) eone
C
      RETURN
      END
      SUBROUTINE MOCONV(EAO,SA,SM,U,T)                        
      IMPLICIT REAL*8 (A-H,O-Z)                                      
      DIMENSION SA(NBATRI),SM(NTRI)                                     
      DIMENSION EAO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)    
      COMMON/BASIS/NBASIS,NBFAO,NBFSO,NBATRI,NTRI
      COMMON/SIGNS/IOFF(500),IPRNT      
      DATA a00 / 0.0D+00 /                                              CPC05240
C                                                                       CPC05250
C   TRANSFORM INTEGRALS FROM AO TO MO BASIS                             CPC05260
      DO 101 II=1,NBFAO                                                 CPC05270
      DO 101 JJ=1,NBFAO                                                 CPC05280
      IIJJ=IOFF(MAX0(II,JJ))+MIN0(II,JJ)                                CPC05290
      T(II,JJ)=SA(IIJJ)                                                 CPC05300
  101 CONTINUE                                                          CPC05310
      DO 103 II=1,NBFAO                                                 CPC05320
      DO 103 J=1,NBASIS                                                 CPC05330
      SUM=a00                                                           CPC05340
      DO 102 JJ=1,NBFAO                                                 CPC05350
      SUM=SUM+EAO(JJ,J)*T(II,JJ)                                        CPC05360
  102 CONTINUE                                                          CPC05370
      U(II,J)=SUM                                                       CPC05380
  103 CONTINUE                                                          CPC05390
      DO 105 I=1,NBASIS                                                 CPC05400
      DO 105 J=1,I                                                      CPC05410
      IJ=IOFF(I)+J                                                      CPC05420
      SUM=a00                                                           CPC05430
      DO 104 II=1,NBFAO                                                 CPC05440
      SUM=SUM+EAO(II,I)*U(II,J)                                         CPC05450
  104 CONTINUE                                                          CPC05460
      SM(IJ)=SUM                                                        CPC05470
  105 CONTINUE                                                          CPC05480
C                                                                       CPC05490
      RETURN                                                            CPC05500
      END                                                               CPC05510
      SUBROUTINE POPUL(P,SS,AOPOP,BATM)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL IANDJ,EQUAL
      integer nbatri,natom
      DIMENSION P(NBATRI),SS(NBATRI),AOPOP(NBATRI),BATM(NATOM,NATOM)
      COMMON/BASIS/NBASIS,NBFAO,NBFSO,NBATRI,NTRI
      COMMON/COORD/CHGS(112),X(3,112),W(112)
      COMMON/FUNCS/NATOM,N3N,NSYM
      COMMON/NORBS/NAST(112),NAED(112)
      COMMON/QVALS/QORB(256),QATM(112),QNET(112)
      COMMON/SIGNS/IPQ(500),IPRNT
      DATA a00 / 0.0D+00 /
    1 FORMAT(/,2X,' ATOMIC ORBITAL POPULATION'/)
    2 FORMAT(//,2X,' GROSS ORBITAL POPULATION'/5X,'AO',9X,'Q(AO)'/)
    3 FORMAT(3X,I4,5X,F10.5)
    4 FORMAT(//,2X,' ATOMIC BOND POPULATION'/)
    5 FORMAT(//,2X,' GROSS ATOMIC POPULATION AND NET CHARGE'/
     1 5X,'ATOM',7X,'Q(ATOM)',7X,' QNET'/)
    6 FORMAT(3X,I4,5X,F10.5,5X,F10.5)
C
C   CALCULATE THE ATOMIC ORBITAL POPULATION
      DO 101 I=1,NBATRI
      AOPOP(I)=P(I)*SS(I)
  101 CONTINUE
      IF(IPRNT.LE.0) GO TO 201
      WRITE(6,1)
c
crsg change nbasis to nbfao for 5d7f9g
crsg  CALL PRINT(AOPOP,NBATRI,NBASIS,6)
      CALL PRINT(AOPOP,NBATRI,nbfao,6)
C
C   CALCULATE THE GROSS ORBITAL POPULATION
  201 CONTINUE
c
crsg  change range of do loops from nbasis to nbfao for 5d7f9g
crsg
crsg  DO 103 I=1,NBASIS
crsg  VAL=a00
crsg  DO 102 J=1,NBASIS
      DO 103 I=1,nbfao
      VAL=a00
      DO 102 J=1,nbfao
      IJ=IPQ(MAX0(I,J))+MIN0(I,J)
      VAL=VAL+AOPOP(IJ)
  102 CONTINUE
      QORB(I)=VAL
  103 CONTINUE
C
      if(iprnt.le.1) go to 202
crsg change NBASIS to nbfao for 5d7f9g
crsg  DO 104 I=1,NBASIS
     
      WRITE(6,2)
      DO 104 I=1,nbfao 
      WRITE(6,3) I,QORB(I)
  104 CONTINUE
C
C   CALCULATE THE ATOMIC BOND POPULATION AND GROSS ATOMIC POPULATION
  202 continue
      DO 105 I=1,NATOM
      QATM(I)=a00
      DO 105 J=1,NATOM
      BATM(I,J)=a00
  105 CONTINUE
      DO 110 IATOM=1,NATOM
      IST=NAST(IATOM)
      IED=NAED(IATOM)
      DO 106 I=IST,IED
      QATM(IATOM)=QATM(IATOM)+QORB(I)
  106 continue
      DO 109 JATOM=1,NATOM
      JST=NAST(JATOM)
      JED=NAED(JATOM)
      DO 107 I=IST,IED
      DO 107 J=JST,JED
      IJ=IPQ(MAX0(I,J))+MIN0(I,J)
      BATM(IATOM,JATOM)=BATM(IATOM,JATOM)+AOPOP(IJ)
  107 CONTINUE
  109 CONTINUE
  110 CONTINUE
C
      WRITE(6,4)
      CALL MATOUT(BATM,NATOM,NATOM,NATOM,NATOM,6)
      WRITE(6,5)
      DO 112 I=1,NATOM
      QNET(I)=CHGS(I)-QATM(I)
      WRITE(6,6) I,QATM(I),QNET(I)
  112 CONTINUE
C
      RETURN
      END
      SUBROUTINE CENTER
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ATM(92)
      COMMON/COORD/CHGS(112),X(3,112),W(112)
      common/dipct/icent
      COMMON/FUNCS/NATOM,N3N,NSYM
      COMMON/QVALS/QORB(256),QATM(112),QNET(112)
      common/savgm/schg(112),savx(112),savy(112),savz(112)
      COMMON/SIGNS/IPQ(500),IPRNT
      DATA a00 / 0.0D+00 /
      DATA DLIM / 1.0D-10 /
      DATA ATM / 1.007825D+00 , 4.00260D+00 , 7.01600D+00 ,
     *           9.01218D+00  , 11.00931D+00 , 12.00000D+00 ,
     *           14.00307D+00 , 15.99491D+00 , 18.99840D+00 ,
     *           19.99244D+00 , 22.98977D+00 , 23.98504D+00 ,
     *           26.98154D+00 , 27.97693D+00 , 30.97376D+00 ,
     *           31.97207D+00 , 34.96885D+00 , 39.96238D+00 ,
     *           38.96371D+00 , 39.96259D+00 , 44.95591D+00 ,
     *           47.94795D+00 , 50.94396D+00 , 51.94051D+00 ,
     *           54.93805D+00 , 55.93494D+00 , 58.93320D+00 ,
     *           57.93535D+00 , 62.92960D+00 , 63.92915D+00 ,
     *           68.92558D+00 , 73.92118D+00 , 74.92160D+00 ,
     *           79.91652D+00 , 78.91834D+00 , 83.91151D+00 ,
Crawdad.7.18.93.  Adding Atomic masses through Uranium (92)
     + 84.9117996d+0, 87.9056249d+0, 88.9058560d+0, 89.9047080d+0,
     + 92.9063780d+0, 97.9054050d+0, 0.0d+0, 101.9043475d+0,
     + 102.905503d+0, 105.903475d+0, 106.905095d+0, 113.9033607d+0,
     + 114.903875d+0, 119.9021990d+0, 120.9038237d+0, 129.906229d+0,
     + 126.904477d+0, 131.904148d+0, 132.905433d+0, 137.905236d+0,
     + 138.906355d+0, 139.905442d+0, 140.907657d+0, 141.907731d+0,
     + 0.0d+0, 151.919741d+0, 152.921243d+0, 157.924111d+0,
     + 158.925350d+0, 163.929183d+0, 164.930332d+0, 165.930305d+0,
     + 168.934225d+0, 173.938873d+0, 174.940785d+0, 179.946561d+0,
     + 180.948014d+0, 183.950953d+0, 186.955765d+0, 191.961487d+0,
     + 192.962942d+0, 194.964785d+0, 196.966560d+0, 201.970632d+0,
     + 204.974410d+0, 207.976641d+0, 208.980388d+0, 0.0d+0,
     + 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 232.038053805d+0,
     + 0.0d+0, 238.050785782d+0/
    1 FORMAT(//,2X,' CARTESIAN COORDINATES'/
     1 5X,' NO.',13X,' X',18X,' Y',18X,' Z',18X,' NUC'/)
    2 FORMAT(2X,I7,5X,4F20.10)
    3 FORMAT(//,2X,' CENTER OF MASS'/
     1 22X,' CMX',16X,' CMY',16X,' CMZ'//14X,3F20.10/)
    4 FORMAT(//,2X,' CARTESIAN COORDINATES W.R.T. CENTER OF MASS'/
     1 5X,' NO.',13X,' X',18X,' Y',18X,' Z',18X,' W'/)
    5 FORMAT(//,2X,' CARTESIAN COORDINATES W.R.T. ORIGIN'/
     1 5X,' NO.',13X,' X',18X,' Y',18X,' Z',18X,' NUC'/)
    6 FORMAT(//,2X,' CENTER OF ELECTRONIC CHARGE'/
     1 22X,' CCX',16X,' CCY',16X,' CCZ'//14X,3F20.10/)
C
      WRITE(6,1)
      DO 101 I=1,NATOM
      schg(i)=chgs(i)
      savx(i)=x(1,i)
      savy(i)=x(2,i)
      savz(i)=x(3,i)
      WRITE(6,2) I,X(1,I),X(2,I),X(3,I),CHGS(I)
  101 continue
C
C   CALCULATE THE CENTER OF MASS
      if(icent.ne.1) go to 301
      SUMW=a00
      SUMWX=a00
      SUMWY=a00
      SUMWZ=a00
      DO 102 I=1,NATOM
      II=CHGS(I)
      W(I)=ATM(II)
      SUMW=SUMW+W(I)
      SUMWX=SUMWX+X(1,I)*W(I)
      SUMWY=SUMWY+X(2,I)*W(I)
      SUMWZ=SUMWZ+X(3,I)*W(I)
  102 CONTINUE
      CMX=SUMWX/SUMW
      CMY=SUMWY/SUMW
      CMZ=SUMWZ/SUMW
      WRITE(6,3) CMX,CMY,CMZ
C
C   NOW DETERMINE THE COORDINATES FOR A CENTER OF MASS FRAME
      WRITE(6,4)
      DO 103 I=1,NATOM
      X(1,I)=X(1,I)-CMX
      X(2,I)=X(2,I)-CMY
      X(3,I)=X(3,I)-CMZ
      WRITE(6,2) I,X(1,I),X(2,I),X(3,I),W(I)
  103 CONTINUE
      GO TO 302
C
C   REFERENCE POINT IS ORIGIN OF SPACED FIXED COORDINATE
  301 CONTINUE
      WRITE(6,5)
      DO 104 I=1,NATOM
      WRITE(6,2) I,X(1,I),X(2,I),X(3,I),CHGS(I)
  104 CONTINUE
C
  302 CONTINUE
      RETURN
      END
      SUBROUTINE DIPAO(P,DIPX,DIPY,DIPZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL IANDJ
      CHARACTER*22 CENT(2)
      integer nbatri
      DIMENSION P(NBATRI),DIPX(NBATRI),DIPY(NBATRI),DIPZ(NBATRI)
      COMMON/BASIS/NBASIS,NBFAO,NBFSO,NBATRI,NTRI
      COMMON/COORD/CHGS(112),C(3,112),W(112)
cets072891      COMMON/EXYZS/EX(7,7,13),EY(7,7,13),EZ(7,7,13)
      common/dipct/icent
      common/dipnl/dnx,dny,dnz
      common/exyzs/ex(7,7,14),ey(7,7,14),ez(7,7,14)
      COMMON/FUNCS/NATOM,N3N,NSYM
      COMMON/GAUSS/ZS(200),CS(200)
      COMMON/IJVEC/IX(84),IY(84),IZ(84),JX(84),JY(84),JZ(84)
      COMMON/INDXX/IJX(225),IJY(225),IJZ(225),IJPOS(225)
      COMMON/NSHEL/KNUC(200),KTYPE(200),KSTART(200),KLOC(200),
     1 KMIN(200),KMAX(200),KNG(200),NSHELL,NPRIM
      COMMON/ROOT /XX,U(9),WT(9),NROOTS
      common/savdp/dtx,dty,dtz
      COMMON/SIGNS/IPQ(500),IPRNT
      COMMON/TEMP0/SX0(36),SY0(36),SZ0(36),DX0(36),DY0(36),DZ0(36)
      DATA ITOL / 10 /
      DATA PI212 / 1.1283791670955D+00 /
      DATA PI32 / 5.56832799683170D+00 /
      DATA DEBYE / 2.541765480D+00 /
      DATA a00,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /
      DATA NQLIM / 7 /
      DATA CENT(1) / ' W.R.T. CENTER OF MASS'/
      DATA CENT(2) / ' W.R.T. ORIGIN        '/
    1 FORMAT(//,2X,' CALCULATION OF DIPOLE MOMENT'/)
    2 FORMAT(//,2X,' P MATRIX'/)
    3 FORMAT(//,2X,' DIPX MATRIX (AO BASIS)'/)
    4 FORMAT(//,2X,' DIPY MATRIX (AO BASIS)'/)
    5 FORMAT(//,2X,' DIPZ MATRIX (AO BASIS)'/)
    6 FORMAT(//,2X,' ======================================'/,
     *          2X,' DIPOLE MOMENT ANALYSIS IN THE AO BASIS'/,
     *          2X,' ======================================')
    7 FORMAT(/,2X,' DIPOLE MOMENTS',A22,' IN DEBYE'/
     1       /,2X,' DEX = ',F20.10,5X,' DNX = ',F20.10,5X,
     2            ' DTX = ',F20.10
     3       /,2X,' DEY = ',F20.10,5X,' DNY = ',F20.10,5X,
     4            ' DTY = ',F20.10
     5       /,2X,' DEZ = ',F20.10,5X,' DNZ = ',F20.10,5X,
     6            ' DTZ = ',F20.10
     7       /)
    8 FORMAT(2X,' TOTAL DIPOLE MOMENT',A22,F20.10,' DEBYE')
C
      TOL=ITOL*2.30258D+00
      IF(IPRNT.LE.2) GO TO 301
      WRITE(6,1)
      WRITE(6,2)
      CALL PRINT(P,NBATRI,NBASIS,6)
C
C     ZERO ALL THE CORE THAT WE NEED FOR S AND ONE ELECTRON INTEGRALS
C
  301 CONTINUE
      DO 101 I=1,NBATRI
      DIPX(I)=a00
      DIPY(I)=a00
      DIPZ(I)=a00
  101 CONTINUE
C
C:::::::::::::::::::::::::::::
C:::::----- I SHELL -----:::::
C:::::::::::::::::::::::::::::
C
      DO 6000 ISH=1,NSHELL
C
      IATOM=KNUC(ISH)
      AX=C(1,IATOM)
      AY=C(2,IATOM)
      AZ=C(3,IATOM)
      IGMIN=KSTART(ISH)
      IGMAX=IGMIN+KNG(ISH)-1
      LIT=KTYPE(ISH)
      MINI=KMIN(ISH)
      MAXI=KMAX(ISH)
      LOCI=KLOC(ISH)
C
C:::::::::::::::::::::::::::::
C:::::----- J SHELL -----:::::
C:::::::::::::::::::::::::::::
C
      DO 5000 JSH=1,ISH
C
      JATOM=KNUC(JSH)
      BX=C(1,JATOM)
      BY=C(2,JATOM)
      BZ=C(3,JATOM)
      JGMIN=KSTART(JSH)
      JGMAX=JGMIN+KNG(JSH)-1
      LJT=KTYPE(JSH)
      MINJ=KMIN(JSH)
      MAXJ=KMAX(JSH)
      LOCJ=KLOC(JSH)
      IANDJ=ISH.EQ.JSH
C
C     ----- WORK OUT INDICES FOR COMBINING 2D INTEGRALS -----
C
      MIJ=0
      JMAX=MAXJ
      II=-1
      DO 104 I=MINI,MAXI
      II=II+1
      JJ=-1
      NX=IX(I)
      NY=IY(I)
      NZ=IZ(I)
      IF(IANDJ) JMAX=I
      DO 103 J=MINJ,JMAX
      MIJ=MIJ+1
      JJ=JJ+1
      IJX(MIJ)=NX+JX(J)
      IJY(MIJ)=NY+JY(J)
      IJZ(MIJ)=NZ+JZ(J)
      IDXI=LOCI+II
      IDXJ=LOCJ+JJ
      IDXS=IPQ(MAX0(IDXI,IDXJ))+MIN0(IDXI,IDXJ)
      IJPOS(MIJ)=IDXS
  103 CONTINUE
  104 CONTINUE
C
      XAB=AX-BX
      YAB=AY-BY
      ZAB=AZ-BZ
      RAB=XAB*XAB+YAB*YAB+ZAB*ZAB
CIBM  ZEROTH ORDER
      MAXIJ=MAX0(LIT+2,LJT+2)
C
      DO 4000 IG=IGMIN,IGMAX
      CSI=CS(IG)*PI32
      AI=ZS(IG)
      AXI=AI*AX
      AYI=AI*AY
      AZI=AI*AZ
C
      IF(IANDJ) JGMAX=IG
      DO 3000 JG=JGMIN,JGMAX
      CSJ=CS(JG)
      BJ=ZS(JG)
      BXJ=BJ*BX
      BYJ=BJ*BY
      BZJ=BJ*BZ
C
      PP=AI+BJ
      TPP=ONE/PP
      TEMP=(AI*BJ*RAB)*TPP
C
C     ----- IGNORE INTEGRAL IF PRE-EXPONENTIAL TOO SMALL -----
C
      IF(TEMP.GT.TOL) GO TO 3000
C
C     ----- DENSITY FACTOR -----
C
      TP=DSQRT(TPP)
      DIJ=CSI*CSJ*DEXP(-TEMP)*TPP
      IF(IANDJ.AND.IG.NE.JG) DIJ=DIJ+DIJ
      PFAC=DIJ*TP
      PTWO=TPP*HALF
C
C     FORM SOME PRODUCTS NEEDED IN INTEGRAL EVALUATION
C
      PX=(AXI+BXJ)*TPP
      PY=(AYI+BYJ)*TPP
      PZ=(AZI+BZJ)*TPP
      PAX=PX-AX
      PAY=PY-AY
      PAZ=PZ-AZ
      PBX=PX-BX
      PBY=PY-BY
      PBZ=PZ-BZ
C
      CALL ECAL(PAX,PAY,PAZ,PBX,PBY,PBZ,PTWO,MAXIJ)
C
      CX=a00
      CY=a00
      CZ=a00
      PCX=PX-CX
      PCY=PY-CY
      PCZ=PZ-CZ
C
      A12=AI+AI
      B12=BJ+BJ
C
      IN=-NQLIM
      DO 106 NI=1,LIT
      IN=IN+NQLIM
      A11=FLOAT(NI-1)
      DO 105 NJ=1,LJT
      JN=IN+NJ
      B11=FLOAT(NJ-1)
C
C   NORMAL OVERLAP INTEGRALS
      SX0(JN)=EX(NI,NJ,1)
      SY0(JN)=EY(NI,NJ,1)
      SZ0(JN)=EZ(NI,NJ,1)
C
C   NORMAL DIPOLE INTEGRALS
      DX0(JN)=EX(NI,NJ,2)+EX(NI,NJ,1)*PCX
      DY0(JN)=EY(NI,NJ,2)+EY(NI,NJ,1)*PCY
      DZ0(JN)=EZ(NI,NJ,2)+EZ(NI,NJ,1)*PCZ
C
C
  105 CONTINUE
  106 CONTINUE
C
      DO 107 I=1,MIJ
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      X=SX0(NX)
      Y=SY0(NY)
      Z=SZ0(NZ)
      IJ=IJPOS(I)
      XFAC=X*PFAC
      YFAC=Y*PFAC
      ZFAC=Z*PFAC
      YZFAC=Y*Z*PFAC
      XZFAC=X*Z*PFAC
      XYFAC=X*Y*PFAC
      DIPX(IJ)=DIPX(IJ)+DX0(NX)*YZFAC
      DIPY(IJ)=DIPY(IJ)+DY0(NY)*XZFAC
      DIPZ(IJ)=DIPZ(IJ)+DZ0(NZ)*XYFAC
  107 CONTINUE
C
 3000 CONTINUE
 4000 CONTINUE
 5000 CONTINUE
 6000 CONTINUE
C
      IF(IPRNT.LE.2) GO TO 302
      WRITE(6,3)
      CALL PRINT(DIPX,NBATRI,NBFAO,6)
      WRITE(6,4)
      CALL PRINT(DIPY,NBATRI,NBFAO,6)
      WRITE(6,5)
      CALL PRINT(DIPZ,NBATRI,NBFAO,6)
C
C   CALCULATE DIPOLE MOMENTS
  302 CONTINUE
C
      DEX=a00
      DEY=a00
      DEZ=a00
      IJ=0
      DO 110 I=1,NBFAO
      DO 110 J=1,I
      IJ=IJ+1
      FAC=TWO
      IF(I.EQ.J) FAC=ONE
      DEX=DEX+DIPX(IJ)*P(IJ)*FAC
      DEY=DEY+DIPY(IJ)*P(IJ)*FAC
      DEZ=DEZ+DIPZ(IJ)*P(IJ)*FAC
  110 CONTINUE
      DEX=-DEX*DEBYE
      DEY=-DEY*DEBYE
      DEZ=-DEZ*DEBYE
      DNX=a00
      DNY=a00
      DNZ=a00
      DO 112 I=1,NATOM
      DNX=DNX+CHGS(I)*C(1,I)
      DNY=DNY+CHGS(I)*C(2,I)
      DNZ=DNZ+CHGS(I)*C(3,I)
  112 CONTINUE
C
      dnx00 = dnx
      dny00 = dny
      dnz00 = dnz
      DNX=DNX*DEBYE
      DNY=DNY*DEBYE
      DNZ=DNZ*DEBYE
      DTX=DNX+DEX
      DTY=DNY+DEY
      DTZ=DNZ+DEZ
C
      WRITE(6,6)
      WRITE(6,7) CENT(ICENT),DEX,DNX,DTX,DEY,DNY,DTY,DEZ,DNZ,DTZ
C
C   CALCULATE THE TOTAL DIPOLE MOMENT
      DD=DTX*DTX+DTY*DTY+DTZ*DTZ
      DD=DSQRT(DD)
      WRITE(6,8) CENT(ICENT),DD
C
cc    itape8 = 8
cc    call ffile(itape8,'file8',0)
cc    rewind itape8
      itape8 = 8
      write(itape8,21) dnx00,dny00,dnz00
   20 format(2i5,f20.10)
   21 format(3f20.10)
      IJ=0
      idpx = 0
      DO 115 I=1,NBFAO
      DO 115 J=1,I
      IJ=IJ+1
      valu = dipx(ij)
      if(dabs(valu).lt.1.0d-06) go to 115
      idpx = idpx + 1
      write(itape8,20) i,j,valu
  115 CONTINUE
      if(iprnt.ge.5) then
        write(6,*) ' idpx  for dipx = ',idpx
      end if
C
      IJ=0
      idpy = 0
      DO 116 I=1,NBFAO
      DO 116 J=1,I
      IJ=IJ+1
      valu = dipy(ij)
      if(dabs(valu).lt.1.0d-06) go to 116
      idpy = idpy + 1
      write(itape8,20) i,j,valu
  116 CONTINUE
      if(iprnt.ge.5) then
        write(6,*) ' idpy  for dipy = ',idpy
      end if
C
      IJ=0
      idpz = 0
      DO 117 I=1,NBFAO
      DO 117 J=1,I
      IJ=IJ+1
      valu = dipz(ij)
      if(dabs(valu).lt.1.0d-06) go to 117
      idpz = idpz + 1
      write(itape8,20) i,j,valu
  117 CONTINUE
      if(iprnt.ge.5) then
        write(6,*) ' idpz  for dipx = ',idpz
      end if
C
c     rewind itape8
C
      RETURN
      END
      SUBROUTINE DIPMO(OCC,EAO,DIPAX,DIPAY,DIPAZ,DIPMX,DIPMY,DIPMZ,
     *                 U,T)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*22 CENT(2)
      DIMENSION OCC(NBASIS)
      DIMENSION DIPAX(NBATRI),DIPAY(NBATRI),DIPAZ(NBATRI)
      DIMENSION DIPMX(NTRI),DIPMY(NTRI),DIPMZ(NTRI)
      DIMENSION EAO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)
      COMMON/BASIS/NBASIS,NBFAO,NBFSO,NBATRI,NTRI
      common/dipct/icent
      common/dipnl/dnx,dny,dnz
      COMMON/SIGNS/IOFF(500),IPRNT
      common/symmt/mosym(256)
      DATA a00 / 0.0D+00 /   
      DATA DEBYE / 2.541765480D+00 /
      DATA DLIMIT / 1.0D-06 /
      DATA CENT(1) / ' W.R.T. CENTER OF MASS'/
      DATA CENT(2) / ' W.R.T. ORIGIN        '/
    1 FORMAT(/,2X,' DIPMX Matrix')
    2 FORMAT(/,2X,' DIPMY Matrix')
    3 FORMAT(/,2X,' DIPMZ Matrix')
    4 FORMAT(//,2X,' ======================================'/,
     *          2X,' DIPOLE MOMENT ANALYSIS IN THE MO BASIS'/,
     *          2X,' ======================================')
    5 FORMAT(/,2X,' DIPOLE CONTRIBUTION IN DEBYE'/,
     1         4X,' NO      MO     SYM     OCC            DEX',
     2         2X,'         DEY           DEZ',/)
    6 FORMAT(2X,I5,3X,I5,5X,I3,2X,F10.5,2X,F12.5,2X,F12.5,2X,F12.5)
    7 FORMAT(/,2X,' TOTAL ELECTRONIC',18X,F12.5,2X,F12.5,2X,F12.5,/)
    8 FORMAT(/,2X,' DIPOLE MOMENTS',A22,' IN DEBYE'/
     1       /,2X,' DEX = ',F20.10,5X,' DNX = ',F20.10,5X,
     2            ' DTX = ',F20.10
     3       /,2X,' DEY = ',F20.10,5X,' DNY = ',F20.10,5X,
     4            ' DTY = ',F20.10
     5       /,2X,' DEZ = ',F20.10,5X,' DNZ = ',F20.10,5X,
     6            ' DTZ = ',F20.10
     7       /)
    9 FORMAT(2X,' TOTAL DIPOLE MOMENT',A25,F20.10,' DEBYE'/)
C
C   TRANSFORM INTEGRALS FROM AO TO MO BASIS
      call moconv(eao,dipax,dipmx,u,t)
      if(iprnt.ge.3) then
        write(6,1)
        call print(dipmx,ntri,nbasis,6)
      end if
      call moconv(eao,dipay,dipmy,u,t)
      if(iprnt.ge.3) then
        write(6,2)
        call print(dipmy,ntri,nbasis,6)
      end if
      call moconv(eao,dipaz,dipmz,u,t)
      if(iprnt.ge.3) then
        write(6,3)
        call print(dipmz,ntri,nbasis,6)
      end if
C
C   CALCULATE ONE ELECTRON ENERGY FOR A TEST
      dex=a00
      dey=a00
      dez=a00
      write(6,4)
      write(6,5)
      imo = 0
      do 101 i=1,nbasis
      ii=ioff(i+1)
      fac = occ(i) * debye
      if(dabs(fac).lt.dlimit) go to 101
      imo = imo + 1
      dx = - dipmx(ii) * fac
      dy = - dipmy(ii) * fac
      dz = - dipmz(ii) * fac
      dex = dex + dx
      dey = dey + dy
      dez = dez + dz
      isym=mosym(i)
      write(6,6) imo,i,isym,occ(i),dx,dy,dz
  101 continue
      write(6,7) dex,dey,dez
C
      dtx = dex + dnx
      dty = dey + dny
      dtz = dez + dnz
C
      WRITE(6,8) CENT(ICENT),DEX,DNX,DTX,DEY,DNY,DTY,DEZ,DNZ,DTZ
C
C   CALCULATE THE TOTAL DIPOLE MOMENT
      DD=DTX*DTX+DTY*DTY+DTZ*DTZ
      DD=DSQRT(DD)
      WRITE(6,9) CENT(ICENT),DD
C
      RETURN
      END
      SUBROUTINE ECAL(PAX,PAY,PAZ,PBX,PBY,PBZ,PTWO,MAXIJ)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL EQUAL
cets072891      COMMON/EXYZS/EX(7,7,13),EY(7,7,13),EZ(7,7,13)
      common/exyzs/ex(7,7,14),ey(7,7,14),ez(7,7,14)
      DATA a00,ONE / 0.0D+00 , 1.0D+00 /
C
      EX(1,1,1)=ONE
      EY(1,1,1)=ONE
      EZ(1,1,1)=ONE
C
      DO 110 II=2,MAXIJ
      DO 109 JJ=1,II
      EQUAL=II.EQ.JJ
      IJTOT=(II-2)+(JJ-1)
      KMAX=II+JJ-1
      DO 108 KK=1,KMAX
      PXX=a00
      PYY=a00
      PZZ=a00
      QXX=a00
      QYY=a00
      QZZ=a00
C
      KM1=KK-1
      IF(KM1.LT.1) GO TO 201
      PXX=PXX+EX(II-1,JJ,KM1)*PTWO
      PYY=PYY+EY(II-1,JJ,KM1)*PTWO
      PZZ=PZZ+EZ(II-1,JJ,KM1)*PTWO
      IF(EQUAL) GO TO 201
      QXX=QXX+EX(JJ,II-1,KM1)*PTWO
      QYY=QYY+EY(JJ,II-1,KM1)*PTWO
      QZZ=QZZ+EZ(JJ,II-1,KM1)*PTWO
C
  201 CONTINUE
      IF(KM1.GT.IJTOT) GO TO 202
      PXX=PXX+EX(II-1,JJ,KK)*PAX
      PYY=PYY+EY(II-1,JJ,KK)*PAY
      PZZ=PZZ+EZ(II-1,JJ,KK)*PAZ
      IF(EQUAL) GO TO 202
      QXX=QXX+EX(JJ,II-1,KK)*PBX
      QYY=QYY+EY(JJ,II-1,KK)*PBY
      QZZ=QZZ+EZ(JJ,II-1,KK)*PBZ
C
  202 CONTINUE
      IF(KK.GT.IJTOT) GO TO 203
      KP1=KK+1
      T=FLOAT(KK)
      PXX=PXX+EX(II-1,JJ,KP1)*T
      PYY=PYY+EY(II-1,JJ,KP1)*T
      PZZ=PZZ+EZ(II-1,JJ,KP1)*T
      IF(EQUAL) GO TO 203
      QXX=QXX+EX(JJ,II-1,KP1)*T
      QYY=QYY+EY(JJ,II-1,KP1)*T
      QZZ=QZZ+EZ(JJ,II-1,KP1)*T
C
  203 CONTINUE
      EX(II,JJ,KK)=PXX
      EY(II,JJ,KK)=PYY
      EZ(II,JJ,KK)=PZZ
      IF(EQUAL) GO TO 205
      EX(JJ,II,KK)=QXX
      EY(JJ,II,KK)=QYY
      EZ(JJ,II,KK)=QZZ
C
  205 CONTINUE
  108 CONTINUE
  109 CONTINUE
  110 CONTINUE
C
      RETURN
      END
      SUBROUTINE GCAL(PCX,PCY,PCZ,VX,VY,VZ,A22,NI,NJ,MAXT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION GX(10),GY(10),GZ(10)
cets072891      COMMON/EXYZS/EX(7,7,13),EY(7,7,13),EZ(7,7,13)
      common/exyzs/ex(7,7,14),ey(7,7,14),ez(7,7,14)
      DATA a00,ONE / 0.0D+00 , 1.0D+00 /
C
      GX(1)=ONE
      GY(1)=ONE
      GZ(1)=ONE
      GX(2)=PCX*A22
      GY(2)=PCY*A22
      GZ(2)=PCZ*A22
C
      IF(MAXT.LE.2) GO TO 201
      DO 101 I=2,MAXT-1
      XI1=FLOAT(I-1)
      GX(I+1)=(PCX*GX(I)+XI1*GX(I-1))*A22
      GY(I+1)=(PCY*GY(I)+XI1*GY(I-1))*A22
      GZ(I+1)=(PCZ*GZ(I)+XI1*GZ(I-1))*A22
  101 CONTINUE
C
  201 CONTINUE
      VX=a00
      VY=a00
      VZ=a00
      DO 102 I=1,MAXT
      VX=VX+EX(NI,NJ,I)*GX(I)
      VY=VY+EY(NI,NJ,I)*GY(I)
      VZ=VZ+EZ(NI,NJ,I)*GZ(I)
  102 CONTINUE
C
      RETURN
      END
      BLOCK DATA
      COMMON/IJVEC/IX(84),IY(84),IZ(84),JX(84),JY(84),JZ(84)
C
C FIRST ROW CONTAINS S,P, AND D INDICES
C SECOND ROW CONTAINS F
C THIRD ROW CONTAINS G
C FOURTH AND FIFTH ROWS CONTAIN  H
C SIXTH AND SEVENTH ROWS CONTAIN  I
C
      DATA IX / 1, 8, 1, 1,15, 1, 1, 8, 8, 1,
     *         22, 1, 1,15,15, 8, 1, 8, 1, 8,
     *         29, 1, 1,22,22, 8, 1, 8, 1,15,15, 1,15, 8, 8,
     *         36, 1, 1,29,29, 8, 1, 8, 1,22,22,15, 1,15, 1,22, 8, 8,15,
     *          8,15,
     *         43, 1, 1,36,36, 8, 1, 8, 1,29,29,15, 1,15, 1,22,22, 1,29,
     *          8, 8,22,22,15, 8,15, 8,15/
      DATA JX / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0,
     *          3, 0, 0, 2, 2, 1, 0, 1, 0, 1,
     *          4, 0, 0, 3, 3, 1, 0, 1, 0, 2, 2, 0, 2, 1, 1,
     *          5, 0, 0, 4, 4, 1, 0, 1, 0, 3, 3, 2, 0, 2, 0, 3, 1, 1, 2,
     *          1, 2,
     *          6, 0, 0, 5, 5, 1, 0, 1, 0, 4, 4, 2, 0, 2, 0, 3, 3, 0, 5,
     *          1, 1, 3, 3, 2, 1, 2, 1, 2/
      DATA IY / 1, 1, 8, 1, 1,15, 1, 8, 1, 8,
     *          1,22, 1, 8, 1,15,15, 1, 8, 8,
     *          1,29, 1, 8, 1,22,22, 1, 8,15, 1,15, 8,15, 8,
     *          1,36, 1, 8, 1,29,29, 1, 8,15, 1,22,22, 1,15, 8,22, 8,15,
     *         15, 8,
     *          1,43, 1, 8, 1,36,36, 1, 8,15, 1,29,29, 1,15,22, 1,22, 8,
     *         29, 8,15, 8,22,22, 8,15,15/
      DATA JY / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1,
     *          0, 3, 0, 1, 0, 2, 2, 0, 1, 1,
     *          0, 4, 0, 1, 0, 3, 3, 0, 1, 2, 0, 2, 1, 2, 1,
     *          0, 5, 0, 1, 0, 4, 4, 0, 1, 2, 0, 3, 3, 0, 2, 1, 3, 1, 2,
     *          2, 1,
     *          0, 6, 0, 1, 0, 5, 5, 0, 1, 2, 0, 4, 4, 0, 2, 3, 0, 3, 1,
     *          5, 1, 2, 1, 3, 3, 1, 2, 2/
      DATA IZ / 1, 1, 1, 8, 1, 1,15, 1, 8, 8,
     *          1, 1,22, 1, 8, 1, 8,15,15, 8,
     *          1, 1,29, 1, 8, 1, 8,22,22, 1,15,15, 8, 8,15,
     *          1, 1,36, 1, 8, 1, 8,29,29, 1,15, 1,15,22,22, 8, 8,22, 8,
     *         15,15,
     *          1, 1,43, 1, 8, 1, 8,36,36, 1,15, 1,15,29,29, 1,22,22, 8,
     *          8,29, 8,15, 8,15,22,22,15/
      DATA JZ / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1,
     *          0, 0, 3, 0, 1, 0, 1, 2, 2, 1,
     *          0, 0, 4, 0, 1, 0, 1, 3, 3, 0, 2, 2, 1, 1, 2,
     *          0, 0, 5, 0, 1, 0, 1, 4, 4, 0, 2, 0, 2, 3, 3, 1, 1, 3, 1,
     *          2, 2,
     *          0, 0, 6, 0, 1, 0, 1, 5, 5, 0, 2, 0, 2, 4, 4, 0, 3, 3, 1,
     *          1, 4, 1, 2, 1, 2, 3, 3, 2/
C
      END
