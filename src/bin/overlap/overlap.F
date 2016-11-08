C**********************************************************************
      subroutine fentry(cc,ia,maxcor)
C**********************************************************************
c     PROGRAM OVERLAP                                                 *
C   A program to calculate overlap integrals                          *
c     and to examine linearity of the basis set                       *
C**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
C**********************************************************************
c   Last updated on August 04, 2005 by Dr. Yukio Yamaguchi            *
C**********************************************************************
c   Last updated on January 18, 1997 by Dr. Yukio Yamaguchi           *
c   Modification for a power2-ibm-aix4.1 system                       *
c   Index has been changes to Indxx                                   *
C**********************************************************************
c   by: Roger Grev, December 20, 1991             search:  crsg       *
c       changes to basset.f and popul.f so it works for 5d7f9g.       *
c       basically changing NBASIS to nbfao.                           *
C**********************************************************************
c       by: Edward Seidl (cets072891)                                 *
c     date: July 28, 1991                                             *
c   reason: change dimension of ex,ey and ez in common exyzs and      *
c           contr in basset to allow for f&g functions                *
c**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
cets030291
c#include <error.h>
      integer errcod, frdi
c
      dimension cc(maxcor),ia(maxcor*2)
      DIMENSION I30(200),A30(100)
      COMMON/BASIS/NBASIS,NBFAO,NBFSO,NBATRI
      COMMON/COORD/CHGS(112),X(3,112),W(112)
      COMMON/FLAGS/MPOINT,MCONST
      COMMON/FUNCS/NATOM,N3N,NSYM
      COMMON/GAUSS/ZS(200),CS(200)
      COMMON/NORBS/NAST(112),NAED(112)
      COMMON/NSHEL/KNUC(200),KTYPE(200),KPRIM(200),KLOC(200),
     1 KMIN(200),KMAX(200),KNG(200),NSHELL,NPRIM
      COMMON/SIGNS/IOFF(500),IPRNT
      common/limit/slimit
      data minovl / 7 /
      EQUIVALENCE (I30,A30)
    1 FORMAT(//,2X,' THE PROGRAM FOR ATOMIC OVERLAP INTEGRALS'/
     1             '             (Version of 08/05/2005)')
    2 FORMAT(2I5)
    3 FORMAT(/,2X,' PARAMETERS'/
     * 2X,' LIMITS = ',I8/
     * 2X,' IPRNT  = ',I8/
     * 2X,' NBASIS = ',I8/
     * 2X,' NBFAO  = ',I8/
     * 2X,' NBFSO  = ',I8/
     * 2X,' NBATRI = ',I8/
     * 2X,' NATOM  = ',I8/
     * 2X,' N3N    = ',I8/
     * 2X,' NSYM   = ',I8/
     * 2X,' NSHELL = ',I8/
     * 2X,' NPRIM  = ',I8/)
    4 FORMAT(/,3X,20A4)
    5 FORMAT(//,2X,' ENUC   = ',F20.10/
     1          2X,' ESCF   = ',F20.10/)
    6 FORMAT(//,2X,' REQUIRED MEMORY EXCEEDS MAXCOR'/
     1          2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)
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
      write(6,*) '       F & G function version '
      WRITE(3,1)
      write(6,*) '  '
      write(6,*) '  maxcor = ',maxcor
C
      CALL LOCATE(INPUT,'# OVERLAP ',IERR)
C
      if(ierr.ne.0) then
        errcod = frdi('PRINT',iprnt)
        errcod = frdi('LIMITS',limits)
      else
        READ(5,2) LIMITS,IPRNT
      endif
      if(limits.eq.0) then
        limits = minovl
      end if
      slimit = (10.0d+00)**(-limits)
      write(3,*) ' slimit = ',slimit
c
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
      NBASIS=NBFSO
      N3N=NATOM*3
      NBATRI=IOFF(NBFAO+1)
C
      WRITE(6,3) LIMITS,IPRNT,NBASIS,NBFAO,NBFSO,NBATRI,NATOM,
     1           N3N,NSYM,NSHELL,NPRIM
C
C   READ IN LABELS
      JUNK=101+MCONST+MPOINT+NCALCS-1
      CALL WREADW(ITAP30,LOCCAL,1,JUNK,JUNK)
      CALL WREADW(ITAP30,I30,60,LOCCAL,LOCCAL)
      WRITE(6,4) (I30(I),I=1,20)
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
C   READ IN BASIS SET INFORMATION
      WRITE(3,21)
   21 FORMAT(//,2X,' NOW YOU ARE IN BASSET'/)
c
      CALL BASSET
c
C   CALCULATE ONE ELECTRON INTEGRALS
      WRITE(3,22)
   22 FORMAT(//,2X,' NOW YOU ARE IN ONEHGT'/)
      IC1=1
      IC2=IC1+NBATRI
      IC3=IC2+NBATRI
      IC4=IC3+NBATRI
      ICMAX=IC4+NBATRI
      write(3,* ) ' ic1 = ',ic1
      write(3,* ) ' ic2 = ',ic2
      write(3,* ) ' ic3 = ',ic3
      write(3,* ) ' ic4 = ',ic4
      if(icmax.gt.maxcor) go to 399
C.................SS      TT      VV      HH......
      CALL ONEHGT(CC(IC1),CC(IC2),CC(IC3),CC(IC4))
c
c   Diagonalize the overlap matrix
      write(3,23)
   23 format(//,2x,' Now you are in sdiag'/)
      ic1=1
      ic2=ic1+nbatri
      ic3=ic2+nbfao
      ic4=ic3+nbfao*nbfao
      ic5=ic4+nbfao
      icmax=ic5+nbfao
      if(icmax.gt.maxcor) go to 399
c................SS      Eigs    Eivs    Fv1     Fv2.....
      call sdiag(cc(ic1),cc(ic2),cc(ic3),cc(ic4),cc(ic5))
      GO TO 400
c
  399 WRITE(6,6) ICMAX,MAXCOR
  400 CONTINUE
      CALL TSTOP(6)
C
      STOP
      END
      SUBROUTINE BASSET
      IMPLICIT REAL*8 (A-H,O-Z)
      integer nbfao,nbfso
      DIMENSION CONTR(1000),I30(200)
      COMMON/BASIS/NBASIS,NBFAO,NBFSO,NBATRI
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
  103 WRITE(6,2) I,ZS(I),CS(I)
      WRITE(6,3) NSHELL
      DO 104 I=1,NSHELL
  104 WRITE(6,4) I,KNUC(I),KTYPE(I),KPRIM(I),KLOC(I),KMIN(I),KMAX(I),
     1 KNG(I)
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
      IF(IPRNT.LE.5) GO TO 202
      WRITE(6,5) NATOM
      DO 106 I=1,NATOM
  106 WRITE(6,6) I,XX(1,I),XX(2,I),XX(3,I),CHGS(I),NAST(I),NAED(I)
  202 CONTINUE
C
      RETURN
      END
      SUBROUTINE ONEHGT(SS,TT,VV,HH)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL IANDJ,EQUAL
      integer nbatri
      DIMENSION SS(NBATRI),TT(NBATRI),VV(NBATRI),HH(NBATRI)
      DIMENSION VX0(36),VY0(36),VZ0(36)
      COMMON/BASIS/NBASIS,NBFAO,NBFSO,NBATRI
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
    1 FORMAT(//,2X,' S matrix'/)
    2 FORMAT(//,2X,' T matrix'/)
    3 FORMAT(//,2X,' V matrix'/)
    4 FORMAT(//,2X,' H matrix'/)
C
      TOL=ITOL*2.30258D+00
C
C     ZERO ALL THE CORE THAT WE NEED FOR S AND ONE ELECTRON INTEGRALS
C
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
      if(iprnt.ge.1) then
        WRITE(6,1)
        CALL PRINT(SS,NBATRI,NBFAO,6)
      end if
      if(iprnt.gt.4) then
        WRITE(6,2)
        CALL PRINT(TT,NBATRI,NBFAO,6)
        WRITE(6,3)
        CALL PRINT(VV,NBATRI,NBFAO,6)
        WRITE(6,4)
        CALL PRINT(HH,NBATRI,NBFAO,6)
      end if
C
      RETURN
      END
      subroutine sdiag(s,eigs,eivs,fv1,fv2)
      implicit real*8 (a-h,o-z)
      dimension s(nbatri),eigs(nbfao),eivs(nbfao,nbfao)
      dimension fv1(nbfao),fv2(nbfao)
      common/basis/nbasis,nbfao,nbfso,nbatri
      common/signs/ioff(500),iprnt
      common/limit/slimit
c
c   Diagonalize the S matrix
      call rsp(nbfao,nbfao,nbatri,s,eigs,1,eivs,fv1,fv2)
      if(iprnt.gt.2) then
        write(6,*) '   '
        write(6,*) '   '
        write(6,*) '  Eigenvalues and Eigenvectors of the S matrix'
        write(6,*) '   '
        call eivout(eivs,eigs,nbfao,nbfao,nbfao,nbfao,6)
        write(6,*) '   '
      end if
c
c   Find the smallest eigenvalue
      absmin=10.0d+00
      do 101 i=1,nbfao
      svalu=eigs(i)
      if(dabs(svalu).lt.absmin) then
        absmin=svalu
      end if
  101 continue
      write(6,*) '   '
      write(6,*) '  The smallest eigenvalue of the S matrix = ',absmin
      write(6,*) '   '
      if(dabs(absmin).lt.slimit) then
        write(6,*) '   '
        write(6,*) '  WARNING!!!   WARNING!!!   WARNING!!!   WARNING!!!'
        write(6,*) '  *************************************************'
        write(6,*) '  ** Possible linear dependency in the basis set **'
        write(6,*) '  *************************************************'
        write(6,*) '  WARNING!!!   WARNING!!!   WARNING!!!   WARNING!!!'
        write(6,*) '   '
        write(6,*) '   '
      end if
c
      return
      end
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
