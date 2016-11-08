C**********************************************************************
      subroutine fentry(cc,ia,maxcor)
C**********************************************************************
c     PROGRAM NATURAL                                                 *
c   Constructs natural orbitals for correlated wave functions         *
C**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
C**********************************************************************
c   Last updated on September 07, 2005 by Dr. Yukio Yamaguchi         *
C**********************************************************************
c   Last updated on January 23, 2003 by Dr. Yukio Yamaguchi           *
C**********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include <error.h>
      integer frdc, frdi, frdf, fwkclr, fwkadd, frdboo, frdpos
      integer fcnt, frdia1
      integer errcod
      character*6 wfn
      character*6 mos
      character*3 upd
C
      dimension cc(maxcor),ia(maxcor*2)
      dimension label(20)
      DIMENSION I30(200),A30(100)
      COMMON/BASIS/NBASIS,NBFAO,NBFSO,NBATRI,NTRI
      common/calif/lpara(1024),apara(1024)
      COMMON/COORD/CHGS(112),X(3,112),W(112)
      COMMON/FLAGS/MPOINT,MCONST
      COMMON/FUNCS/NATOM,N3N,NSYM
      COMMON/GAUSS/ZS(200),CS(200)
      common/mobas/iaoso,naoso
      COMMON/NORBS/NAST(112),NAED(112)
      COMMON/NSHEL/KNUC(200),KTYPE(200),KPRIM(200),KLOC(200),
     1 KMIN(200),KMAX(200),KNG(200),NSHELL,NPRIM
      COMMON/PARAS/MXCOEF,LOCVEC
      COMMON/SIGNS/IOFF(500),IPRNT
      common/symmt/mosym(256),nlamda(10)
      COMMON/CI101/IOPEN,IOCC,JOCC,KOCC
      EQUIVALENCE (I30,A30)
    1 FORMAT(//,2X,' THE NATURAL ORBITAL PROGRAM',/
     1          2X,'   (Version of 01/23/2003)  ',/)
    2 FORMAT(3I5)
    3 FORMAT(/,2X,' PARAMETERS'/
     * 2X,' WFN     = ',A6/
     * 2X,' MOS     = ',A6/
     * 2X,' UPD     = ',A3,3x/
     * 2X,' IWFN    = ',I8/
     * 2X,' IAOSO   = ',I8/
     * 2X,' IUPD    = ',I8/
     * 2X,' IPRNT   = ',I8/
     * 2X,' NBASIS  = ',I8/
     * 2X,' NBFAO   = ',I8/
     * 2X,' NBFSO   = ',I8/
     * 2X,' NAOSO   = ',I8/
     * 2X,' NBATRI  = ',I8/
     * 2X,' NTRI    = ',I8/
     * 2X,' NATOM   = ',I8/
     * 2X,' N3N     = ',I8/
     * 2X,' NSYM    = ',I8/
     * 2X,' IOPEN   = ',I8/
     * 2X,' NSHELL  = ',I8/
     * 2X,' NPRIM   = ',I8/)
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
      write(6,*) '       F & G function version '
      WRITE(3,1)
      write(6,*) '  '
      write(6,*) '  maxcor = ',maxcor
C
C   IPRINT IS A PARAMETER FOR SIZE OF OUTPUT
C   ALLOWED WAVE FUNCTIONS
C     SCF     : WFN = 'SCF   '
C     CI      : WFN = 'CI    '
C     MCSCF   : WFN = 'MCSCF '
C     CCSD    : WFN = 'CCSD  '
C     CCSD(T) : WFN = 'CCSD_T'
C
      errcod = fwkclr()
      errcod = fwkadd(5,'NATURAL')
      errcod = fwkadd(5,'DEFAULT')
      wfn = 'SCF   '
      errcod = frdc('WFN',wfn)
      if(wfn.eq.'SCF   ') then
        iwfn = 1
      end if
      if(wfn.eq.'CI    ') then
        iwfn = 2
      end if
      if(wfn.eq.'MCSCF ') then
        iwfn = 3
      end if
      if(wfn.eq.'CCSD  ') then
        iwfn = 4
      end if
      if(wfn.eq.'CCSD_T') then
        iwfn = 5
      end if
      mos = 'SO_MO '
      errcod = frdc('MOS',mos)
      if(mos.eq.'SO_MO ') then
        iaoso = 1
      end if
      if(mos.eq.'AO_MO ') then
        iaoso = 2
      end if
      upd = 'YES'
      errcod = frdc('UPDATE',upd)
      if(upd.eq.'YES') then
        iupd = 1
      end if
      if(upd.eq.'NO ') then
        iupd = 2
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
      if(iaoso.eq.1) then
        naoso = nbfso
      end if
      if(iaoso.eq.2) then
        naoso = nbfao
      end if
C                                                                       CPC05250`
C
      WRITE(6,3) WFN,MOS,UPD,IWFN,IAOSO,IUPD,IPRNT,NBASIS,NBFAO,NBFSO,
     *           NAOSO,NBATRI,NTRI,
     *           NATOM,N3N,NSYM,IOPEN,NSHELL,NPRIM
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
C
c*******************************************************
c   Dynamic allocation
c*******************************************************
c     IC1  : Uas(nbfao,nbasis)
c     IC2  : Eig(nbasis)
c     IC3  : Occ(nbasis)
c     IC4  : Eso(nbfso,nbasis)
c     IC5  : Eao(nbfao,nbasis)
c     IC6  : Eno(naoso,nbasis)
c     IC7  : Pao(nbatri)
c     IC8  : Pmo(ntri)
c     IC9  : Pno(ntri)
c     IC10 : Sao(nbatri)
c     IC11 : Hao(nbatri)
c     IC12 : Hmo(ntri)
c*******************************************************
      IC1=1
      IC2=IC1+NBFAO*NBFSO
      IC3=IC2+NBASIS
      IC4=IC3+NBASIS
      IC5=IC4+NBFSO*NBFSO
      IC6=IC5+NBFAO*NBFSO
      IC7=IC6+NAOSO*NBFSO
      IC8=IC7+nbatri
      IC9=IC8+ntri
      IC10=IC9+ntri
      IC11=IC10+nbatri
      IC12=IC11+nbatri
      ICSAV=IC12+ntri
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
      ic13=icsav
      icmax=ic13+nbfao*nbfso
      IF(ICMAX.GT.MAXCOR) GO TO 399
C................UAS     EIG     OCC     ESO     EAO     EE.......
      CALL EIVEC(CC(IC1),CC(IC2),CC(IC3),CC(IC4),CC(IC5),CC(IC13))
      call rclose(itap30,3)
C
c=============================================================
c   Read in the density matrix for a CI wavefunction
c=============================================================
      IF(WFN.EQ.'CI    '.OR.WFN.EQ.'MCSCF ') THEN
CC      WRITE(6,23)
        WRITE(3,23)
   23   FORMAT(//,2X,' NOW YOU ARE IN CIDENS'/)
        ic13=icsav
        ic14=ic13+ntri
        ic15=ic14+nbfao*nbfao
        icmax=ic15+nbfao*nbfao
        IF(ICMAX.GT.MAXCOR) GO TO 399
C...................EAO     PAO     PMO     TEMP.....
        CALL CIDENS(CC(IC5),CC(IC7),CC(IC8),CC(IC13),
C...................U        T........
     *              CC(IC14),CC(IC15))
       END IF
C
c=================================================================
c   Read in the density matrix for a CCSD or CCSD(T) wavefunction
c=================================================================
      IF(WFN.EQ.'CCSD  '.OR.WFN.EQ.'CCSD_T') THEN
CC      WRITE(6,24)
        WRITE(3,24)
   24   FORMAT(//,2X,' NOW YOU ARE IN CCDENS'/)
        icmax=icsav
        IF(ICMAX.GT.MAXCOR) GO TO 399
C...................PAO     PMO.....
        CALL CCDENS(CC(IC7),CC(IC8))
      END IF
C
CC    WRITE(6,25)
      WRITE(3,25)
   25 FORMAT(//,2X,' NOW YOU ARE IN NATUR_ORB'/)
      ic13=icsav
      ic14=ic13+nbasis
      ic15=ic14+nbasis
      ic16=ic15+nbasis
      ic17=ic16+nbasis
      ic18=ic17+nbasis
      ic19=ic18+nbasis
      ic20=ic19+ntri
      ic21=ic20+ntri
      ic22=ic21+ntri
      ic23=ic22+naoso*nbasis
      ic24=ic23+naoso*nbasis
      ic25=ic24+naoso*nbasis
      ic26=ic25+naoso*nbasis
      icmax=ic26+naoso*nbasis
      if(iprnt.ge.5) then
        write(6,* ) ' ic1   = ',ic1
        write(6,* ) ' ic2   = ',ic2
        write(6,* ) ' ic3   = ',ic3
        write(6,* ) ' ic4   = ',ic4
        write(6,* ) ' ic5   = ',ic5
        write(6,* ) ' ic6   = ',ic6
        write(6,* ) ' ic7   = ',ic7
        write(6,* ) ' ic8   = ',ic8
        write(6,* ) ' ic9   = ',ic9
        write(6,* ) ' ic10  = ',ic10
        write(6,* ) ' ic11  = ',ic11
        write(6,* ) ' ic12  = ',ic12
        write(6,* ) ' ic13  = ',ic13
        write(6,* ) ' ic14  = ',ic14
        write(6,* ) ' ic15  = ',ic15
        write(6,* ) ' ic16  = ',ic16
        write(6,* ) ' ic17  = ',ic17
        write(6,* ) ' ic18  = ',ic18
        write(6,* ) ' ic19  = ',ic19
        write(6,* ) ' ic20  = ',ic20
        write(6,* ) ' ic21  = ',ic21
        write(6,* ) ' ic22  = ',ic22
        write(6,* ) ' ic23  = ',ic23
        write(6,* ) ' ic24  = ',ic24
        write(6,* ) ' ic25  = ',ic25
        write(6,* ) ' ic26  = ',ic26
        write(6,* ) ' icmax = ',icmax
      end if
      IF(ICMAX.GT.MAXCOR) GO TO 399
C                                                                       CPC05250
C.................ESO     EAO     ENO     PMO     PNO     OCCNO....
      CALL NORBIT(CC(IC4),CC(IC5),CC(IC6),CC(IC8),CC(IC9),CC(IC13),
C.................PROOT    QROOT    RROOT    FV1      FV2......
     *            CC(IC14),CC(IC15),CC(IC16),CC(IC17),CC(IC18),
C.................PTEMP    QTEMP    RTEMP....  
     *            CC(IC19),CC(IC20),CC(IC21),
C.................PEIV     QEIV     REIV     SEIV     ENEIV....
     *            CC(IC22),CC(IC23),CC(IC24),CC(IC25),CC(IC26))
c
C   CALCULATE ONE ELECTRON INTEGRALS FOR A TEST
CC    WRITE(6,26)
      WRITE(3,26)
   26 FORMAT(//,2X,' NOW YOU ARE IN ONEHGT'/)
      IC13=icsav
      IC14=IC13+NBATRI
      ICMAX=IC14+NBATRI
      if(iprnt.ge.5) then
        write(6,* ) ' ic1   = ',ic1
        write(6,* ) ' ic2   = ',ic2
        write(6,* ) ' ic3   = ',ic3
        write(6,* ) ' ic4   = ',ic4
        write(6,* ) ' ic5   = ',ic5
        write(6,* ) ' ic6   = ',ic6
        write(6,* ) ' ic7   = ',ic7
        write(6,* ) ' ic8   = ',ic8
        write(6,* ) ' ic9   = ',ic9
        write(6,* ) ' ic10  = ',ic10
        write(6,* ) ' ic11  = ',ic11
        write(6,* ) ' ic12  = ',ic12
        write(6,* ) ' ic13  = ',ic13
        write(6,* ) ' ic14  = ',ic14
        write(6,* ) ' icmax = ',icmax
      end if
      IF(ICMAX.GT.MAXCOR) GO TO 399
C.................PAO     SAO      HAO      TT       VV.......
      CALL ONEHGT(CC(IC7),CC(IC10),CC(IC11),CC(IC13),CC(IC14))
C
C   CALCULATE MULLIKEN POPULATION
CC    WRITE(6,27)
      WRITE(3,27)
   27 FORMAT(//,2X,' NOW YOU ARE IN POPUL'/)
      IC13=icsav
      IC14=IC13+NBATRI
      ICMAX=IC14+NATOM*NATOM
      IF(ICMAX.GT.MAXCOR) GO TO 399
C.................PAO     SAO      AOPOP    BATM.....
      CALL  POPUL(CC(IC7),CC(IC10),CC(IC13),CC(IC14))
C
C=====================================================================
C  Correlated wave function
C  Using a standard one particle density matrix
C=====================================================================
      write(6,*) ' '
      write(6,*) '  **************************************'
      write(6,*) '  ANALYSIS USING SCF ORBITALS           '
      write(6,*) '  **************************************'
      write(6,*) ' '
CC    WRITE(6,28)
      WRITE(3,28)
   28 FORMAT(//,2X,' NOW YOU ARE IN TRANS_MO'/)
      IC13=icsav
      IC14=IC13+NBFAO*NBFAO
      ICMAX=IC14+NBFAO*NBFAO
      if(iprnt.ge.5) then
        write(6,* ) ' ic1   = ',ic1
        write(6,* ) ' ic2   = ',ic2
        write(6,* ) ' ic3   = ',ic3
        write(6,* ) ' ic4   = ',ic4
        write(6,* ) ' ic5   = ',ic5
        write(6,* ) ' ic6   = ',ic6
        write(6,* ) ' ic7   = ',ic7
        write(6,* ) ' ic8   = ',ic8
        write(6,* ) ' ic9   = ',ic9
        write(6,* ) ' ic10  = ',ic10
        write(6,* ) ' ic11  = ',ic11
        write(6,* ) ' ic12  = ',ic12
        write(6,* ) ' ic13  = ',ic13
        write(6,* ) ' ic14  = ',ic14
        write(6,* ) ' icmax = ',icmax
      end if
      IF(ICMAX.GT.MAXCOR) GO TO 399
C...................EAO     HAO      HMO      U        T........
      CALL TRANS_MO(CC(IC5),CC(IC11),CC(IC12),CC(IC13),CC(IC14))
C
CC    WRITE(6,29)
      WRITE(3,29)
   29 FORMAT(//,2X,' NOW YOU ARE IN EONE_CORR'/)
C....................PMO     HMO......
      CALL EONE_CORR(CC(IC8),CC(IC12))
C
C=====================================================================
C  Correlated wave function
C  Using a natural one particle density matrix
C=====================================================================
      write(6,*) ' '
      write(6,*) '  **************************************'
      write(6,*) '  ANALYSIS USING NATURAL ORBITALS       '
      write(6,*) '  **************************************'
      write(6,*) ' '
C
      if(iaoso.eq.1) then
        ic13=icsav
        ICMAX=IC13+NBFAO*NBASIS
        if(iprnt.ge.5) then
          write(6,* ) ' before sotoao '
          write(6,* ) ' iaoso = ',iaoso
          write(6,* ) ' nbfao = ',nbfao
          write(6,* ) ' nbfso = ',nbfso
          write(6,* ) ' nbasis = ',nbasis
          write(6,* ) ' iprnt  = ',iprnt
          write(6,* ) ' ic1   = ',ic1
          write(6,* ) ' ic2   = ',ic2
          write(6,* ) ' ic3   = ',ic3
          write(6,* ) ' ic4   = ',ic4
          write(6,* ) ' ic5   = ',ic5
          write(6,* ) ' ic6   = ',ic6
          write(6,* ) ' ic7   = ',ic7
          write(6,* ) ' ic8   = ',ic8
          write(6,* ) ' ic9   = ',ic9
          write(6,* ) ' ic10  = ',ic10
          write(6,* ) ' ic11  = ',ic11
          write(6,* ) ' ic12  = ',ic12
          write(6,* ) ' ic13  = ',ic13
          write(6,* ) ' icmax = ',icmax
        end if
C...................UAS     ENO     ENAO.....
        call sotoao(CC(IC1),CC(IC6),CC(IC13))
C
CC      WRITE(6,30)
        WRITE(3,30)
   30   FORMAT(//,2X,' NOW YOU ARE IN TRANS_MO'/)
        IC13=icsav
        IC14=IC13+NBFAO*NBASIS
        IC15=IC14+NBFAO*NBFAO
        ICMAX=IC14+NBFAO*NBFAO
        IF(ICMAX.GT.MAXCOR) GO TO 399
C.....................ENAO     HAO      HMO      U        T........
        CALL TRANS_MO(CC(IC13),CC(IC11),CC(IC12),CC(IC14),CC(IC15))
C
        if(iupd.eq.1) then
CC        WRITE(6,31)
          WRITE(3,31)
   31     FORMAT(//,2X,' NOW YOU ARE IN UPDATE'/)
          IC13=icsav
          IC14=IC13+NBASIS
          IC15=IC14+NBASIS*NBASIS
          ICMAX=IC14+NBASIS*NBASIS
          IF(ICMAX.GT.MAXCOR) GO TO 399
C.....................ENO     EIG      EE       ESS......
          CALL UPDATE(CC(IC6),CC(IC13),CC(IC14),CC(IC15))
        end if
      end if
C
      if(iaoso.eq.2) then
CC      WRITE(6,32)
        WRITE(3,32)
   32   FORMAT(//,2X,' NOW YOU ARE IN TRANS_MO'/)
        IC13=icsav
        IC14=IC13+NBFAO*NBFAO
        ICMAX=IC14+NBFAO*NBFAO
        IF(ICMAX.GT.MAXCOR) GO TO 399
C.....................ENO     HAO      HMO      U        T........
        CALL TRANS_MO(CC(IC6),CC(IC11),CC(IC12),CC(IC13),CC(IC14))
      end if
C
CC    WRITE(6,33)
      WRITE(3,33)
   33 FORMAT(//,2X,' NOW YOU ARE IN EONE_CORR'/)
C......................PNO     HMO......
      CALL EONE_CORR(CC(IC9),CC(IC12))
C
C=====================================================================
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
      DIMENSION NCLOSD(10),NOPEN(10),OCT(2)
      DIMENSION ALPT(15),BETT(15)
      COMMON/BASIS/NBASIS,NBFAO,NBFSO,NBATRI,NTRI
      COMMON/FUNCS/NATOM,N3N,NSYM
      COMMON/PARAS/MXCOEF,LOCVEC
      COMMON/SIGNS/IPQ(500),IPRNT
      common/symmt/mosym(256),nlamda(10)
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
      RETURN
      END
      SUBROUTINE CIDENS(EAO,PAO,PMO,TEMP,U,T)
C=======================================================================
C  Modification for LINUX symstem
C  Subroutine Mwrit has been commented out
C  on September 07, 2005 by Dr. Yukio Yamaguchi
C=======================================================================
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension iorscf(256)
      DIMENSION PAO(NBATRI),PMO(NTRI)
      dimension eao(nbfao,nbasis)
      dimension temp(ntri),u(nbfao,nbfao),t(nbfao,nbfao)
      COMMON/BASIS/NBASIS,NBFAO,NBFSO,NBATRI,NTRI
      COMMON/SIGNS/IOFF(500),IPRNT
      DATA a00 / 0.0D+00 /
    1 FORMAT(//,2X,' THE CI DENSITY MATRIX (MO BASIS)'/)
    2 FORMAT(//,2X,' THE CI DENSITY MATRIX (AO BASIS)'/)
C
C   READ IN THE DENSITY MATRIX IN THE MO BASIS
      call initmf(1)
      call mread(iorscf,4)
c*    write(6,*) ' nbasis = ',nbasis
c*    write(6,*) ' ntri   = ',ntri
c*    write(6,*) (iorscf(i),i=1,nbasis)
cyy9-07-05  call mread(temp,35)
C
C================================================
C ITAP59 has been used for pass 'one-pdm'
C================================================
      itap59 = 59
      call rfile(itap59)
      call sread(itap59,temp,ntri*2)
      call rclose(itap59,3)
C================================================
C
      if(iprnt.ge.3) then
        WRITE(6,1)
        CALL PRINT(temp,NTRI,NBASIS,6)
      end if
      ij=0
      do 101 i=1,nbasis
      ii=iorscf(i)
      do 101 j=1,i
      jj=iorscf(j)
      ij=ij+1
      iijj=ioff(max0(ii,jj)) + min0(ii,jj)
      pmo(iijj)=temp(ij)
  101 continue
      if(iprnt.ge.3) then
        WRITE(6,1)
        CALL PRINT(PMO,NTRI,NBASIS,6)
      end if
C
c   Transform one-pdm from mo to ao basis
      do 103 i=1,nbasis
      do 103 jj=1,nbfao
      valu=a00
      do 102 j=1,nbasis
      ij=ioff(max0(i,j))+min0(i,j)
      valu=valu+eao(jj,j)*pmo(ij)
  102 continue
      u(i,jj)=valu
  103 continue
c
      do 105 ii=1,nbfao
      do 105 jj=1,nbfao
      valu=a00
      do 104 i=1,nbasis
      valu=valu+eao(ii,i)*u(i,jj)
  104 continue
      t(ii,jj)=valu
  105 continue
      if(iprnt.gt.4) then
        write(6,*) '  '
        write(6,*) '  One-PDM matrix in AO basis'
        write(6,*) '  '
        call matout(t,nbfao,nbfao,nbfao,nbfao,6)
      end if
c
      ij=0
      do 107 i=1,nbfao
      do 107 j=1,i
      ij=ij+1
      pao(ij)=t(i,j)
  107 continue
      if(iprnt.gt.4) then
        write(6,*) '  '
        write(6,*) '  One-PDM in AO basis (triang.)'
        write(6,*) '  '
        call print(pao,nbatri,nbfao,6)
      end if
c
cyy9-07-05      call mwrit(pao,34)
c*    write(6,*) ' nbfao  = ',nbfao
c*    write(6,*) ' nbatri = ',nbatri
      if(iprnt.ge.3) then
        WRITE(6,2)
        CALL PRINT(PAO,NBATRI,NBFAO,6)
      end if
C
      RETURN
      END
      SUBROUTINE CCDENS(PAO,PMO)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PAO(NBATRI),PMO(NTRI)
      COMMON/BASIS/NBASIS,NBFAO,NBFSO,NBATRI,NTRI
      COMMON/SIGNS/IOFF(500),IPRNT
      DATA a00 / 0.0D+00 /
    1 FORMAT(//,2X,' THE CC DENSITY MATRIX (MO BASIS)'/)
    2 FORMAT(//,2X,' THE CC DENSITY MATRIX (AO BASIS)'/)
C
C   READ IN THE DENSITY MATRIX
      nt = nbasis
      nt2 = ioff(nt+1)
      ntao2 = ioff(nbfao+1)
      itap79 = 79
      call rfile(itap79)
C
C   READ IN THE DENSITY MATRIX IN THE MO BASIS
      ipt = 1 + intowp(nt2+nt2+1)
      call wreadw(itap79,pmo,intowp(nt2),ipt,ipt)
      if(iprnt.ge.3) then
        WRITE(6,1)
        CALL PRINT(PMO,NTRI,NBASIS,6)
      end if
C
C   READ IN THE DENSITY MATRIX IN THE AO BASIS
      ipt = 1 + intowp(nt2+nt2+1+nt2+nt*nt)
      call wreadw(itap79,pao,intowp(ntao2),ipt,ipt)
      if(iprnt.ge.3) then
        WRITE(6,2)
        CALL PRINT(PAO,NBATRI,NBFAO,6)
      end if
      call rclose(itap79,3)
C
      RETURN
      END
      SUBROUTINE EONE_CORR(PMO,HMO)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PMO(NTRI),HMO(NTRI)
      COMMON/BASIS/NBASIS,NBFAO,NBFSO,NBATRI,NTRI
      COMMON/SIGNS/IOFF(500),IPRNT
      DATA a00,one,two / 0.0d+00 , 1.0d+00, 2.0d+00 /
    1 FORMAT(/,2X,' EONE_CORR (calculated in the MO basis) = ',F20.10/)
C                                                                       CPC05250
C   CALCULATE ONE ELECTRON ENERGY FOR A TEST
      eone=a00
      ij=0
      do 101 i=1,nbasis
      do 101 j=1,i
      ij=ij+1
      fac=2.0d+00
      if(i.eq.j) fac=1.0d+00
      eone = eone + hmo(ij)*pmo(ij)*fac
  101 continue
      write(6,1) eone
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
      SUBROUTINE POPUL(PAO,SAO,AOPOP,BATM)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL IANDJ,EQUAL
      integer nbatri,natom
      DIMENSION PAO(NBATRI),SAO(NBATRI),AOPOP(NBATRI),BATM(NATOM,NATOM)
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
      AOPOP(I)=PAO(I)*SAO(I)
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
      SUBROUTINE SOTOAO(UAS,ESO,EAO)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION UAS(NBFAO,NBASIS),ESO(NBFSO,NBASIS),EAO(NBFAO,NBASIS)
      COMMON/BASIS/NBASIS,NBFAO,NBFSO,NBATRI,NTRI
      COMMON/SIGNS/IOFF(500),IPRNT
      data a00 / 0.0d+00 /
    1 FORMAT(/,2X,' HMO Matrix')
C
C   TRANSFORM EIGENVECTORS FROM SO TO AO BASIS
      write(6,*) '  in sotoao '
      write(6,*) '  nbfao  =  ',nbfao
      write(6,*) '  nbfso  =  ',nbfso
      write(6,*) '  nbasis =  ',nbasis
      write(6,*) '  iprnt  =  ',iprnt
      do 102 i=1,nbfao
      do 102 j=1,nbasis
      valu=a00
      do 101 k=1,nbfso
      valu=valu+uas(i,k)*eso(k,j)
  101 CONTINUE
      eao(i,j)=valu
  102 continue
C
      if(iprnt.ge.3) then
        write(6,1)
        call matout(eao,nbfao,nbasis,nbfao,nbasis,6)
      end if
C
      return
      end
      SUBROUTINE TRANS_MO(EAO,HAO,HMO,U,T)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION HAO(NBATRI),HMO(NTRI)
      DIMENSION EAO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)
      COMMON/BASIS/NBASIS,NBFAO,NBFSO,NBATRI,NTRI
      COMMON/SIGNS/IOFF(500),IPRNT
    1 FORMAT(/,2X,' HMO Matrix')
C
C   TRANSFORM ONE ELECTRON INTEGRALS FROM AO TO MO BASIS
      call moconv(eao,hao,hmo,u,t)
      if(iprnt.ge.3) then
        write(6,1)
        call print(hmo,ntri,nbasis,6)
      end if
C
      return
      end
      SUBROUTINE VECOUT(A,B,NAD,NBD,M,N,IOUT)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*6 LINE(21)
      DIMENSION A(NAD,NBD),B(NBD)
      DATA LINE / 21*'------' /
    1 FORMAT(2X,10(7X,I5))
    2 FORMAT(2X,21A6)
    3 FORMAT(2X,I2,2X,10F12.7)
    4 FORMAT(/,6X,10F12.7)
    5 FORMAT(/)
C
      II=0
      JJ=0
  200 II=II+1
      JJ=JJ+1
      KK=10*JJ
      NN=N
      IF(N.GT.KK) NN=KK
      LL=2*(NN-II+1)+1
      WRITE(IOUT,1) (I,I=II,NN)
      WRITE(IOUT,2) (LINE(I),I=1,LL)
      DO 101 I=1,M
  101 WRITE(IOUT,3) I,(A(I,J),J=II,NN)
      WRITE(IOUT,4) (B(J),J=II,NN)
      IF(N.LE.KK) GO TO 201
      WRITE(IOUT,5)
      II=KK
      GO TO 200
  201 RETURN
      END
