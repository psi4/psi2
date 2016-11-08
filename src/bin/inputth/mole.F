      SUBROUTINE MOLE(ZAN,C,EX,C1,CSPD,NGNCN,NPERSH,KSTART,KATOM,KTYPE,
     1  KNG,KLOC,KMIN,KMAX,NEWSH,TITLE,T,INVT,ISO,ICT,NS,KS,LBLAT,CSINP,
     &  stpool,maxcor)
CChoas 12-03-92 added call to enuc to calc. and write nuc. rep. energy.

C
      IMPLICIT REAL*8 (A-H,O-Z)
      integer maxcor
      integer stpool(maxcor)
      COMMON /MAXDIM/ MXAT,MXTRN,MXPR,MXSH,MXSYM,MXCF,MXT,MXISO,MXPSYM
      COMMON /PSIZE/  MSU,KAORDS,MCONSU,MGCSU,MSTU,MRU,MCSU,MCTU,MCONU,
     1                MCU,MSFU,MGU,MSFRU,MNSFU,MPRU,MCXU,MCCU,MBLU,MRCRU
      COMMON /OUTPUT/ NPRINT,NGEOM,NORMF,NORMP
c06-09-89  COMMON /RESTAR/ TIMLIM
      COMMON /INFOA/  NAT,NUM,NQW,NX
      COMMON /NSHEL/  NSHELL,MSHELL
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV
      COMMON /INOUT/  IOUT,ICHECK,ITAP30
clj031091
      integer errcod, bool, frdboo, frdc
      integer inold
      common /instyl/ inold
      INTEGER INVT(MXTRN),ISO(MXSH,MXISO),ICT(MXAT,MXTRN),NPERSH(5)
c530  INTEGER TITLE(20),NS(MXAT),KS(MXAT),LBLAT(2,MXAT)
clj032690 The following change by Gus(?) makes no sense
clj      INTEGER TITLE(17),NS(MXAT),KS(MXAT),LBLAT(2,MXAT)
      INTEGER TITLE(20),NS(MXAT),KS(MXAT),LBLAT(2,MXAT)
      character*80 label
      integer etitle(20)
      REAL*8 ZAN(MXAT),C(3,MXAT),EX(MXPR),C1(MRCRU),
     1  CSPD(MRCRU,MXPR,5),T(MXT), CSINP(MRCRU,MXPR,5), nucener
      INTEGER NGNCN(MXSH),KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),
     1  KNG(MXSH),KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NEWSH(MXSH,MXTRN)
      DIMENSION JFIX(9)
      equivalence(etitle,label)
      DATA ZERO /0.0D+00/
C
      IR=5
      IW=6
clj031091
      if (inold.eq.1) then
C***********************************************************************
C           FIND INPUT RECORDS                                         *
C***********************************************************************
      CALL LOCATE(IR,'# INPUT ##',IERROR)
      IF (IERROR.NE.0) THEN
        WRITE (IP,90)
        WRITE (IW,90)
   90   FORMAT (//,' # INPUT ##  ERROR--NO INPUT FOUND',//)
        CALL qabort
      END IF
C***********************************************************************
C           READ TITLE AND OPTIONS (ONLY NPRINT IS USED)               *
C***********************************************************************
      READ(IR,7774) TITLE
c06-09-89 READ(IR,7773) TIMLIM,NGEOM,NPRINT,NORMF,NORMP
      READ(IR,7773) NGEOM,NPRINT,NORMF,NORMP
      IF(NPRINT.EQ.7) WRITE(IP,7774) TITLE
      IF(NPRINT.EQ.7) WRITE(IP,7776) NPRINT,NGEOM,NORMF,NORMP
      IF(NPRINT.EQ.6) GO TO 1
      WRITE(IW,7777) TITLE
    1 CONTINUE
      IF(NORMF.NE.1) NORMF=0
      IF(NORMP.NE.1) NORMP=0
      WRITE(IW,7775) NPRINT,NGEOM,NORMF,NORMP

clj021091 New style input:
      else
        write(6,*)
        label = ' '
        errcod = frdc('LABEL',label)
        write(6,'(A9,A71)') ' LABEL = ',label
        errcod = frdboo('NORMBASIS',BOOL)
        if (errcod.ne.0) BOOL = 1
        write(6,*) 'NORMBASIS =',bool
        if (bool.eq.0) then
          normf = 1
        else
          normf = 0
          endif
        errcod = frdboo('PRIMNORM',BOOL)
        if (errcod.ne.0) BOOL = 1
        write(6,*) 'PRIMNORM =',bool
        if (bool.eq.0) then
          normp = 1
        else
          normp = 0
          endif
        errcod = frdboo('PRESERVE',BOOL)
        if (errcod.ne.0) BOOL = 0
        write(6,*) 'PRESERVE =',bool
        if (bool.eq.0) then
          ngeom = 0
        else
          ngeom = 1
          endif
        nprint = 0
        do 101 i = 1, 20
          title(i) = etitle(i)
 101      continue
        endif

C***********************************************************************
C        READ SYMMETRY POINT GROUP OF THE MOLECULE,                    *
C        AND GENERATE THE TRANSFORMATION MATRICES.                     *
C***********************************************************************
      CALL PTGRP(T,T,T,INVT)
C***********************************************************************
C        READ THE BASIS SET FOR THE UNIQUE CENTERS,                    *
C        AND GENERATE THE MOLECULAR BASIS SET.                         *
C***********************************************************************
      if (inold.eq.1) then
      CALL ATOMS(EX,C1,CSPD(1,1,1),CSPD(1,1,2),CSPD(1,1,3),
     1 CSPD(1,1,4),CSPD(1,1,5),NGNCN,KSTART,KATOM,KTYPE,KNG,KLOC,
     2 KMIN,KMAX,NPERSH,T,ISO,NEWSH,ZAN,C,INVT,ICT,NS,KS,LBLAT,
     3 CSINP(1,1,1),CSINP(1,1,2),CSINP(1,1,3),CSINP(1,1,4),CSINP(1,1,5))
      else
      CALL PATOMS(EX,C1,CSPD(1,1,1),CSPD(1,1,2),CSPD(1,1,3),
     1 CSPD(1,1,4),CSPD(1,1,5),NGNCN,KSTART,KATOM,KTYPE,KNG,KLOC,
     2 KMIN,KMAX,NPERSH,T,ISO,NEWSH,ZAN,C,INVT,ICT,NS,KS,LBLAT,
     3 CSINP(1,1,1),CSINP(1,1,2),CSINP(1,1,3),CSINP(1,1,4),CSINP(1,1,5),
     & stpool,maxcor)
        endif
CChaos 12-03-92
       nucener = ENUC(ZAN,C)
       write (6, 7800) nucener
 7800  FORMAT ('Nuclear Repulsion Energy = ',F20.10)
C
      RETURN
c06-09-89 7773 FORMAT(F10.0,4I5)
 7773 FORMAT(10X,4I5)
 7774 FORMAT(20A4)
 7775 FORMAT(/,
c06-09-89  1 ' TIME LIMIT                             (NOT USED)',/,
     1 ' OUTPUT OPTION                               = ',I5,/,
     1 ' OPTION TO PRESERVE M.O. COEFFICIENTS        = ',I5,/,
     1 ' OPTION FOR NORMALIZATION OF BASIS FUNCTIONS = ',I5,/,
     1 ' OPTION FOR UNNORMALIZATION OF PRIMITIVES    = ',I5,/)
 7776 FORMAT(' NPRINT  NGEOM  NORMF  NORMP OPTIONS ARE',/,4I5)
 7777 FORMAT(1H ,///,72(1H*),/,1H*,70X,1H*,/,
     1 1H*,11X,'  INPUT PROGRAM FOR HFS PROGRAMS (HONDO-LIKE)  ',
     2 12X,1H*,/,1H*,70X,1H*,/,1H*,10X,14X,
     3 ' KING-DUPUIS-RYS-SAXE ',10X,14X,1H*,/,1H*,70X,1H*,/,
     4 1H*,1X,20A4,1X,1H*,/,1H*,70X,1H*,/,72(1H*),/)
      END
