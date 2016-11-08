      SUBROUTINE OUT(ZAN,C,ICT,NS,KS,EX,CSPD,NGNCN,KSTART,KATOM,KTYPE,
     1  KNG,KLOC,KMIN,KMAX,NPERSH,ISC,IPC,CC,LAB,LOC,LOC2,KLASS,LAMBDA,
     2  IADDR,IPP,T,INVT,ISO,TITLE,LBLAT,SB,PTR,DTR,MXPRD,FTR,GTR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      EQUIVALENCE (ITEMP,TEMP)
      DIMENSION ZZZ(3),ICRUD(200)
      COMMON /OUTPUT/ NPRINT,NGEOM,NORMF,NORMP
      COMMON /INOUT/  IOUT,ICHECK,ITAP30
      COMMON /SALCS/  NIRRED
      COMMON /SYMTRY/ NT
      COMMON /PSIZE/  MSU,KAORDS,MCONSU,MGCSU,MSTU,MRU,MCSU,MCTU,MCONU,
     1                MCU,MSFU,MGU,MSFRU,MNSFU,MPRU,MCXU,MCCU,MBLU,MRCRU
      COMMON /MAXDIM/ MXAT,MXTRN,MXPR,MXSH,MXSYM,MXCF,MXT,MXISO,MXPSYM
      COMMON /INFOA/  NAT,NUM,NQ,NX
      COMMON /NSHEL/  NSHELL,MSHELL
      REAL*8 ZAN(MXAT),C(3,MXAT),EX(MXPR),CSPD(MRCRU,MXPR,5)
      REAL*8 CC(MXCF),T(MXT),PTR(3,3,NT),DTR(6,6,NT),FTR(10,10,NT)
      REAL*8 GTR(15,15,NT)
      INTEGER ICT(MXAT,MXTRN),NS(MXAT),KS(MXAT),NGNCN(MXSH),KSTART(MXSH)
      INTEGER KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),KLOC(MXSH),NPERSH(5)
      INTEGER KMIN(MXSH),KMAX(MXSH),ISC(MXSH,MXPSYM)
      INTEGER IPC(MXSH,MXPSYM),LAB(MXCF),LOC(MXAT,5),LOC2(MXAT,5)
      INTEGER KLASS(MXAT),LAMBDA(MXSYM),IADDR(MXSYM),IPP(MXSYM)
      INTEGER INVT(MXTRN),ISO(MXSH,MXISO),TITLE(20)
      INTEGER SB(MXSYM),LBLAT(2,MXAT)
      INTEGER FCB(16),ITEMP(200)
      REAL*8   TEMP(100)
      integer isizof
C***********************************************************************
C     FILE30 IS INITIALLY MADE AS FOLLOWS:                             *
C       1 - 100 TITLE AND 80 ZEROS                                     *
C     101 - 300 CONSTANTS                                              *
C     301 - 500 POINTERS TO FILE30 FOR VARIOUS ARRAYS THAT ARE STORED  *
C     501 - 600 DIRECTORY TO CALCULATIONS, WITH END OF ARRAYS POINTER  *
C               IN 501 AND THE REST ZEROS                              *
C     601 - END OF ARRAYS STORED                                       *
C     NEXT 60   DIRECTORY FOR FIRST CALCULATION (TITLE AND 40 ZEROS)   *
C     NEXT 20   POINTERS TO VECTORS (ZEROS)                            *
C     NEXT 6*NATOMS  GEOMETRY                                          *
C     NEXT 2*10  ENERGIES FOR FUTURE USE                               *
C***********************************************************************
      ITAP30=30
      CALL RFILE(ITAP30)
      IEND=1
      MPOINT=200
      MCONST=200
      MCALCS=100
      NCALCS=0
C***********************************************************************
C     ZERO FIRST 200 ELEMENTS OF ITEMP.                                *
C***********************************************************************
      DO 1 I=1,MPOINT
        ITEMP(I)=0
    1 CONTINUE
C***********************************************************************
C     WRITE OUT TITLE, THEN 80 ZEROS IN FILE30.                        *
C***********************************************************************
      CALL WWRITW(ITAP30,TITLE,20,IEND,IEND)
      DO 2 I=2,5
        CALL WWRITW (ITAP30,ITEMP,20,IEND,IEND)
    2 CONTINUE
C***********************************************************************
C           ZERO ARRAY FOR CONSTANTS,POINTERS, AND CALCULATIONS        *
C     IF NGEOM = 1, PRESERVE INFO FROM FILE30 BEFORE WRITING BACK OUT  *
C***********************************************************************
      IF(NGEOM.EQ.1) THEN
        CALL WREADW(ITAP30,ITEMP,43,101,JUNK)
        IEND0  = ITEMP(1)
        NSYMHF = ITEMP(41)
        MXCOEF = ITEMP(42)
        IOPEN  = ITEMP(43)
      ENDIF
      CALL WWRITW (ITAP30,ITEMP,MCONST,IEND,IEND)
C
      IF(NGEOM.EQ.1) THEN
        CALL WREADW(ITAP30,ITEMP,MPOINT,IEND,JUNK)
        I29 = ITEMP(29)
        I30 = ITEMP(30)
      ENDIF
      CALL WWRITW (ITAP30,ITEMP,MPOINT,IEND,IEND)
C
      IF(NGEOM.EQ.0) THEN
        CALL WWRITW (ITAP30,ITEMP,MCALCS,IEND,IEND)
      ELSE
        IEND = IEND + NCALCS
      ENDIF
C***********************************************************************
C           START WRITING HONDO ARRAYS TO FILE, PUTTING ADDRESSES      *
C               IN ITEMP AS WE GO. FIRST IS NUCLEAR CHARGES.           *
C***********************************************************************
      ITEMP(1)=IEND
      CALL WWRITW(ITAP30,ZAN,INTOWP(NAT),IEND,IEND)
C***********************************************************************
C           TRANSFORMATION TABLE OF ATOMS                              *
C***********************************************************************
      ITEMP(2)=IEND
      DO 3 I=1,NT
        CALL WWRITW(ITAP30,ICT(1,I),NAT,IEND,IEND)
    3 CONTINUE
C***********************************************************************
C          NUMBER OF SHELLS PER ATOM -- NS(NAT)                        *
C***********************************************************************
      ITEMP(3)=IEND
      CALL WWRITW(ITAP30,NS,NAT,IEND,IEND)
C***********************************************************************
C           POINTER TO SHELLS ON ATOM -- KS(NAT)                       *
C***********************************************************************
      ITEMP(4)=IEND
      CALL WWRITW(ITAP30,KS,NAT,IEND,IEND)
C***********************************************************************
C           EXPONENTS OF PRIMITIVE GAUSSIANS -- EX(MXPR)               *
C***********************************************************************
      ITEMP(5)=IEND
      CALL WWRITW(ITAP30,EX,INTOWP(MXPR),IEND,IEND)
C***********************************************************************
C     S,P,D,F,G CONTRACTION COEFFICIENTS -- CSPD(MRCRU,MXPR,5)         *
C***********************************************************************
      ITEMP(6)=IEND
      DO 4 I=1,5
        CALL WWRITW(ITAP30,CSPD(1,1,I),INTOWP(MRCRU*MXPR),IEND,IEND)
    4 CONTINUE
C***********************************************************************
C           POINTER TO PRIMITIVES FOR A SHELL -- KSTART(NSHELL)        *
C***********************************************************************
      ITEMP(7)=IEND
      CALL WWRITW(ITAP30,KSTART,NSHELL,IEND,IEND)
C***********************************************************************
C           ATOM ON WHICH NTH SHELL IS CENTRED -- KATOM(NSHELL)        *
C***********************************************************************
      ITEMP(8)=IEND
      CALL WWRITW(ITAP30,KATOM,NSHELL,IEND,IEND)
C***********************************************************************
C           TYPE (S,P,D...) OF NTH SHELL -- KTYPE(NSHELL)              *
C***********************************************************************
      ITEMP(9)=IEND
      CALL WWRITW(ITAP30,KTYPE,NSHELL,IEND,IEND)
C***********************************************************************
C     NUMBER OF CONTRACTED FUNCTIONS IN NTH SHELL -- KNG(NSHELL)       *
C***********************************************************************
      ITEMP(10)=IEND
      CALL WWRITW(ITAP30,KNG,NSHELL,IEND,IEND)
C***********************************************************************
C      POINTER TO FIRST BASIS FNCT IN SHELL --- KLOC(NSHELL)           *
C***********************************************************************
      ITEMP(11)=IEND
      CALL WWRITW(ITAP30,KLOC,NSHELL,IEND,IEND)
C***********************************************************************
C       FIRST FUNCTION (S,X,Y,Z,XX...) IN SHELL -- KMIN(NSHELL)        *
C***********************************************************************
      ITEMP(12)=IEND
      CALL WWRITW(ITAP30,KMIN,NSHELL,IEND,IEND)
C***********************************************************************
C           LAST FUNCTION IN SHELL -- KMAX(NSHELL)                     *
C***********************************************************************
      ITEMP(13)=IEND
      CALL WWRITW(ITAP30,KMAX,NSHELL,IEND,IEND)
C***********************************************************************
C         PACKED SO TRANSFORMATION MATRICES -- ISC(NSHELL,NPSYM)       *
C***********************************************************************
      ITEMP(14)=IEND
      NPSYM=(NIRRED+3)/4
      DO 5 I=1,NPSYM
        CALL WWRITW(ITAP30,ISC(1,I),NSHELL,IEND,IEND)
    5 CONTINUE
C***********************************************************************
C       PACKED SO TRANSFORMATION INFORMATION -- IPC(NSHELL,NPSYM)      *
C***********************************************************************
      ITEMP(15)=IEND
      DO 6 I=1,NPSYM
        CALL WWRITW(ITAP30,IPC(1,I),NSHELL,IEND,IEND)
    6 CONTINUE
C***********************************************************************
C           LABELS OF IRREPS -- SB(NIRRED)                             *
C***********************************************************************
      ITEMP(16)=IEND
      CALL WWRITW(ITAP30,SB,NIRRED,IEND,IEND)
C***********************************************************************
C           COEFFICIENTS OF SO'S -- CC(MXCF)                           *
C***********************************************************************
      ITEMP(17)=IEND
      CALL WWRITW(ITAP30,CC,INTOWP(MXCF),IEND,IEND)
c     write(6,*)'cc written to file30'
c     write(6,*)(cc(i),i=1,mxcf)
C***********************************************************************
C           LABELS OF SO MATRICES -- LAB(MXCF)                         *
C***********************************************************************
      ITEMP(18)=IEND
      CALL WWRITW(ITAP30,LAB,MXCF,IEND,IEND)
C***********************************************************************
C           POINTER TO BEGINNING OF SO MATRIX -- LOC(NAT,3)            *
C***********************************************************************
      ITEMP(19)=IEND
      DO 7 I=1,5
        CALL WWRITW(ITAP30,LOC(1,I),NAT,IEND,IEND)
    7 CONTINUE
C***********************************************************************
C           POINTER TO END OF SO MATRIX -- LOC2(NAT,3)                 *
C***********************************************************************
      ITEMP(20)=IEND
      DO 8 I=1,5
        CALL WWRITW(ITAP30,LOC2(1,I),NAT,IEND,IEND)
    8 CONTINUE
C***********************************************************************
C           CLASS AN ATOM BELONGS TO -- KLASS(NAT)                     *
C***********************************************************************
      ITEMP(21)=IEND
      CALL WWRITW(ITAP30,KLASS,NAT,IEND,IEND)
C***********************************************************************
C           DEGENERACY OF NTH IRREP -- LAMBDA(NIRRED)                  *
C***********************************************************************
      ITEMP(22)=IEND
      CALL WWRITW(ITAP30,LAMBDA,NIRRED,IEND,IEND)
C***********************************************************************
C           POINTER TO START OF NTH SOS ? -- IADDR(NIRRED)             *
C***********************************************************************
      ITEMP(23)=IEND
      CALL WWRITW(ITAP30,IADDR,NIRRED,IEND,IEND)
C***********************************************************************
C           NUMBER OF SOS OF NTH SYMMETRY ? -- IPP(NIRRED)             *
C***********************************************************************
      ITEMP(24)=IEND
      CALL WWRITW(ITAP30,IPP,NIRRED,IEND,IEND)
C***********************************************************************
C           TRANSFORMATION MATRICES -- T(MXT)                          *
C***********************************************************************
      MXT=9*NT
      ITEMP(25)=IEND
      CALL WWRITW(ITAP30,T,INTOWP(MXT),IEND,IEND)
C***********************************************************************
C           INVERSE OF SYMMETRY OPERATION -- INVT(NT)                  *
C***********************************************************************
      ITEMP(26)=IEND
      CALL WWRITW(ITAP30,INVT,NT,IEND,IEND)
C***********************************************************************
C      PACKED TRANSFORMATION MATRIX FOR SHELLS -- ISO(NSHELL,MXISO)    *
C***********************************************************************
      ITEMP(27)=IEND
      NISO=(NT+3)/4
      DO 9 I=1,NISO
        CALL WWRITW(ITAP30,ISO(1,I),NSHELL,IEND,IEND)
    9 CONTINUE
C***********************************************************************
C           LABEL OF ATOMS -- LBLAT(2,NAT)                             *
C***********************************************************************
      ITEMP(28)=IEND
      CALL WWRITW(ITAP30,LBLAT,2*NAT,IEND,IEND)
C***********************************************************************
C       TRANSFORMATION MATRICES FOR P FUNCTIONS -- PTR(3,3,NT)         *
C***********************************************************************
      ITEMP(31)=IEND
      CALL WWRITW(ITAP30,PTR,INTOWP(3*3*NT),IEND,IEND)
C***********************************************************************
C       TRANSFORMATION MATRICES FOR D FUNCTIONS -- DTR(6,6,NT)         *
C***********************************************************************
      ITEMP(32)=IEND
      CALL WWRITW(ITAP30,DTR,INTOWP(6*6*NT),IEND,IEND)
C***********************************************************************
C       TRANSFORMATION MATRICES FOR F FUNCTIONS -- FTR(10,10,NT)       *
C***********************************************************************
      ITEMP(33)=IEND
      CALL WWRITW(ITAP30,FTR,INTOWP(10*10*NT),IEND,IEND)
C***********************************************************************
C       TRANSFORMATION MATRICES FOR G FUNCTIONS -- GTR(15,15,NT)       *
C***********************************************************************
      ITEMP(34)=IEND
      CALL WWRITW(ITAP30,GTR,INTOWP(15*15*NT),IEND,IEND)
C***********************************************************************
C          NUMBER OF CONTRACTIONS IN EACH SHELL -- NGNCN               *
C***********************************************************************
      ITEMP(35)=IEND
      CALL WWRITW(ITAP30,NGNCN,NSHELL,IEND,IEND)
C***********************************************************************
C          WHETHER THE SHELL IS PURE D5,F7,G9 OR NOT - NPERSH -        *
C***********************************************************************
      ITEMP(36)=IEND
      CALL WWRITW(ITAP30,NPERSH,5,IEND,IEND)
C***********************************************************************
C           WRITE POINTERS TO ARRAYS ON FILE                           *
C***********************************************************************
      IF(NGEOM.EQ.1) THEN
        ITEMP(29) = I29
        ITEMP(30) = I30
      ENDIF
      CALL WWRITW (ITAP30,ITEMP,MPOINT,301,JUNK)
C***********************************************************************
C           WRITE OUT DIRECTORY TO CALCULATIONS                        *
C***********************************************************************
      IF(NGEOM.EQ.0) THEN
        DO 14 I=1,MCALCS
          ITEMP(I)=0
   14   CONTINUE
        ITEMP(1)=IEND
        CALL WWRITW (ITAP30,ITEMP,MCALCS,501,JUNK)
C***********************************************************************
C        WRITE OUT DIRECTORY FOR FIRST CALCULATION AFTER ARRAYS        *
C***********************************************************************
        DO 11 I=1,20
          ITEMP(I)=0
   11   CONTINUE
      ENDIF
      CALL WWRITW (ITAP30,TITLE,20,IEND,IEND)
      DO 12 I=1,2
        CALL WWRITW (ITAP30,ITEMP,20,IEND,IEND)
   12 CONTINUE
C***********************************************************************
C           ZERO POINTERS TO VECTORS ETC.                              *
C***********************************************************************
      IF(NGEOM.EQ.0) THEN
        CALL WWRITW (ITAP30,ITEMP,20,IEND,IEND)
      ELSE
        IEND = IEND + 20
      ENDIF
C***********************************************************************
C           WRITE OUT GEOMETRY                                         *
C***********************************************************************
c The following line is required to make the geometry nonzero on SGI
c machines:
      write(6,*) 'to write geometry: length ',intowp(3*nat)
      CALL WWRITW(ITAP30,C,intowp(3*nat),IEND,IEND)
C***********************************************************************
C           ZERO ENERGIES FOR LATER USE                                *
C***********************************************************************
      IF(NGEOM.EQ.0) THEN
        DO 13 I=1,10
          TEMP(I)=0.0D+00
   13   CONTINUE
        CALL WWRITW (ITAP30,TEMP,INTOWP(10),IEND,IEND)
      ELSE
        IEND = IEND + INTOWP(10)
      ENDIF
C***********************************************************************
C           PUT ALL CONSTANTS INTO ITEMP TO PUT ON FILE                *
C***********************************************************************
      DO 10 I=1,MCONST
        ITEMP(I)=0
   10 CONTINUE
      ITEMP(1)  = IEND
      ITEMP(2)  = MPOINT
      ITEMP(3)  = MCONST
      ITEMP(4)  = MCALCS
      ITEMP(5)  = 1
      ITEMP(6)  = MSU
      ITEMP(7)  = KAORDS
      ITEMP(8)  = MCONSU
      ITEMP(9)  = MGCSU
      ITEMP(10) = MSTU
      ITEMP(11) = MRU
      ITEMP(12) = MCSU
      ITEMP(13) = MCTU
      ITEMP(14) = MCONU
      ITEMP(15) = MCU
      ITEMP(16) = MSFU
      ITEMP(17) = MGU
      ITEMP(18) = MSFRU
      ITEMP(19) = NAT
      ITEMP(20) = 0
      ITEMP(21) = 0
      ITEMP(22) = NUM
      ITEMP(23) = NQ
      ITEMP(24) = NX
      ITEMP(25) = 0
      ITEMP(26) = 0
      ITEMP(27) = NSHELL
      ITEMP(28) = NIRRED
      ITEMP(29) = NT
      ITEMP(30) = MXAT
      ITEMP(31) = MXTRN
      ITEMP(32) = MXPR
      ITEMP(33) = MXSH
      ITEMP(34) = MXSYM
      ITEMP(35) = MXCF
      ITEMP(36) = MXT
      ITEMP(37) = MXISO
      ITEMP(38) = MXPSYM
      ITEMP(39) = NPSYM
      ITEMP(40) = NISO
      IF(NGEOM.EQ.1) THEN
        ITEMP(1) = IEND0
        ITEMP(41) = NSYMHF
        ITEMP(42) = MXCOEF
        ITEMP(43) = IOPEN
      ELSE
        ITEMP(41) = 0
        ITEMP(42) = 0
        ITEMP(43) = 0
      ENDIF
      ITEMP(44) = MRCRU
      CALL WWRITW (ITAP30,ITEMP,44,101,JUNK)
C
      WRITE(IOUT,90) IEND*isizof('INT_IN_BYTES')
   90 FORMAT (/,' FILE30 USED ',I8,' BYTES' )
C
      RETURN
      END
