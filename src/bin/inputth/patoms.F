      SUBROUTINE PATOMS(EX,C1,CS,CP,CD,CF,CG,NGNCN,KSTART,KATOM,
     1  KTYPE,KNG,KLOC,KMIN,KMAX,NPERSH,T,ISO,NEWSH,ZAN,C,INVT,ICT,NS,
     2  KS,LBLAT,CSINP,CPINP,CDINP,CFINP,CGINP,stpool,maxcor)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      integer maxcor
      integer stpool(maxcor)
      COMMON /MAXDIM/ MXAT,MXTRN,MXPR,MXSH,MXSYM,MXCF,MXT,MXISO,MXPSYM
      COMMON /PSIZE/  MSU,KAORDS,MCONSU,MGCSU,MSTU,MRU,MCSU,MCTU,MCONU,
     1                MCU,MSFU,MGU,MSFRU,MNSFU,MPRU,MCXU,MCCU,MBLU,MRCRU
      REAL*8 EX(MXPR),C1(MRCRU),CS(MRCRU,MXPR),CP(MRCRU,MXPR),
     1  CD(MRCRU,MXPR),CF(MRCRU,MXPR),CG(MRCRU,MXPR),T(MXT),
     2  ZAN(MXAT),C(3,MXAT),CSINP(MRCRU,MXPR),CPINP(MRCRU,MXPR),
     3  CDINP(MRCRU,MXPR),CFINP(MRCRU,MXPR),CGINP(MRCRU,MXPR)
      INTEGER NGNCN(MXSH),KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),
     1  KNG(MXSH),KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),ISO(MXSH,MXISO),
     2  NEWSH(MXSH,MXTRN),INVT(MXTRN),ICT(MXAT,MXTRN),NS(MXAT),
     3  KS(MXAT),LBLAT(2,MXAT),NPERSH(5)
#include <error.h>
#include <basis.h>
      integer errcod,fwkadd,prcntr,frdboo,frdc,gencon
      character*10 cunits
      integer puream

c See if pure angular momentum is desired (default = FALSE)
      errcod = frdboo('PUREAM',puream)
      if (errcod.ne.EOK) puream = 0
      write(6,'(2X,A,I2)') 'PUREAM = ',puream

c See if we want to look for generalized contractions (default = FALSE)
      errcod = frdboo('GENCON',gencon)
      if (errcod.ne.EOK) gencon = 0
      write(6,'(2X,A,I2)') 'GENCON = ',gencon

c Get the units
      errcod = frdc('UNITS',cunits)
      if (errcod.ne.EOK) then
        cunits = 'BOHR'
      else if (cunits.eq.'ANGSTROM') then
        cunits = 'ANGSTROMS'
      else if (cunits.ne.'BOHR'.and.cunits.ne.'ANGSTROMS') then
        write(6,*) 'ERROR: bad ''UNITS'''
        call pabort
        call qabort
        endif
      write(6,'(2X,2A)') 'UNITS = ',cunits

c Initialize memory
      call coinit(stpool,maxcor)

c Read the basis set
      call intrb(stpool,maxcor)

c Read the charges
      call intrc(stpool,maxcor)

c Read the geometry
      call intrg(stpool,maxcor)

c Print the basis set
      write(6,*)
      call intpb(stpool,maxcor)

c Print the charges
      write(6,*)
      write(6,*) 'Charges on unique centers (a.u.):'
      call intpc(stpool,maxcor)

c Print the geometry
      write(6,*)
      write(6,*) 'Geometry of unique centers (a.u.):'
      call intpg(stpool,maxcor,'CARTESIAN','BOHR')

clj031091 Convert the data read in above to the internal data needed
c        for the input program.
      call PATOMC(EX,C1,CS,CP,CD,CF,CG,NGNCN,KSTART,KATOM,
     &  KTYPE,KNG,KLOC,KMIN,KMAX,NPERSH,T,ISO,NEWSH,ZAN,C,INVT,ICT,NS,
     &  KS,LBLAT,CSINP,CPINP,CDINP,CFINP,CGINP,
     &  nshell,nprimt,nprim2,nucent,
     &  stpool(pexpon),stpool(pconco),stpool(ptypes),stpool(ppshc),
     &  stpool(pnshc),stpool(pnprim),stpool(ppprim),stpool(pcharg),
     &  stpool(pgeom),cunits,puream,gencon
     &  )


      return
      end
clj031091 ---------------------------------------------------------
c    Note that nshell -> nshll happened in the function call to
c    avoid a naming conflict.
      SUBROUTINE PATOMc(EX,C1,CS,CP,CD,CF,CG,NGNCN,KSTART,KATOM,
     1  KTYPE,KNG,KLOC,KMIN,KMAX,NPERSH,T,ISO,NEWSH,ZAN,C,INVT,ICT,NS,
     2  KS,LBLAT,CSINP,CPINP,CDINP,CFINP,CGINP,
     &  nshll,nprimt,nprim2,nucent,
     &  expon,conco,types,pshc,
     &  nshc,nprim,pprim,charg,
     &  geom,cunits,puream,gencon
     &  )
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include <error.h>
      integer errcod, frdc1
      integer nshll,nprimt,nprim2,nucent
      real*8 expon(nprimt),conco(nprimt)
      integer types(nshll)
      integer pshc(nucent),nshc(nucent),nprim(nshll)
      integer pprim(nshll)
      real*8 charg(nucent),geom(3,nucent)
      integer puream,gencon
      character*(*) cunits
      COMMON /MAXDIM/ MXAT,MXTRN,MXPR,MXSH,MXSYM,MXCF,MXT,MXISO,MXPSYM
      COMMON /PSIZE/  MSU,KAORDS,MCONSU,MGCSU,MSTU,MRU,MCSU,MCTU,MCONU,
     1                MCU,MSFU,MGU,MSFRU,MNSFU,MPRU,MCXU,MCCU,MBLU,MRCRU
c06-09-89  COMMON /RESTAR/ TIMLIM
      COMMON /OUTPUT/ NPRINT,NGEOM,NORMF,NORMP
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV
      COMMON /INOUT/  IOUT,ICHECK,ITAP30
      COMMON /INFOA/  NAT,NUM,NQW,NX
      COMMON /NSHEL/  NSHELL,MSHELL
      COMMON /SYMTRY/ NT
      COMMON /TRANSF/ XSMAL,YSMAL,ZSMAL,XNEW,YNEW,ZNEW,XP,YP,ZP
      COMMON /FRAME/  U1,U2,U3,V1,V2,V3,W1,W2,W3,X0,Y0,Z0
      COMMON /ISOPAC/ INDIN(120),INDOUT(30)
      DIMENSION NBFS(21),MINF(21),MAXF(21),LABEL(21),NANGM(21),NBFSP(21)
      DIMENSION EXX(6),CSS(6),CPP(6),CDD(6),JLABEL(10)
      integer ishc,iushll
      character*4 cjlab(10),clab(21),ctype
cges dec 12 1989  increase istack to 30 primitives
cges may 24 1989  increase istack to 35 primitives
c     DIMENSION INTYP(200),ISTACK(20),JNTYP(3),MS(200)
c     DIMENSION INTYP(200),ISTACK(30),JNTYP(3),MS(200)
      DIMENSION INTYP(200),ISTACK(35),JNTYP(3),MS(200)
      REAL*8 EX(MXPR),C1(MRCRU),CS(MRCRU,MXPR),CP(MRCRU,MXPR),
     1  CD(MRCRU,MXPR),CF(MRCRU,MXPR),CG(MRCRU,MXPR),T(MXT),
     2  ZAN(MXAT),C(3,MXAT),CSINP(MRCRU,MXPR),CPINP(MRCRU,MXPR),
     3  CDINP(MRCRU,MXPR),CFINP(MRCRU,MXPR),CGINP(MRCRU,MXPR)
      INTEGER NGNCN(MXSH),KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),
     1  KNG(MXSH),KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),ISO(MXSH,MXISO),
     2  NEWSH(MXSH,MXTRN),INVT(MXTRN),ICT(MXAT,MXTRN),NS(MXAT),
     3  KS(MXAT),LBLAT(2,MXAT),NPERSH(5)
      integer prcntr
      integer iush2,ipcnt
      integer ilbl(2)
      character*8 clbl
      character*40 albl
      equivalence (ilbl,clbl)
      CHARACTER*80 GET
      REAL*4 ATNUM
      REAL*8 LINE,LINEV,LEESF
      PARAMETER (ZERO=0.0D0, PT125=0.125D0, PT3S16=0.1875D0,
     1  PT25=0.25D0, PT5=0.5D0, PT0625=0.0625D0)
      DATA LABEL /4H  1S,4H  2S,4H  2P,4H 2SP,4H  3S,4H  3P,4H 3SP,
     1  4H  3D,4H  4S,4H  4P,4H 4SP,4H  4D,4H  5S,4H  5P,4H  5D,
     2  4H   S,4H   P,4H   D,4H   F,4H   G,4H   H/
      DATA JLABEL /4H   S,4H   P,4H  D5,4H  F7,4H  G9,4H   S,4H   P,
     1  4H  D6,4H F10,4H G15/
      DATA CJLAB  /'   S','   P','  D5','  F7','  G9','   S','   P',
     1  '  D6',' F10',' G15'/
      DATA NBFS /1,1,3,4,1,3,4,6,1,3,4,6,1,3,6,1,3,6,10,15,21/
      DATA NBFSP/1,1,3,4,1,3,4,6,1,3,4,6,1,3,6,1,3,5,7,9,11/
      DATA MINF /1,1,2,1,1,2,1,5,1,2,1,5,1,2,5,1,2,5,11,21,36/
      DATA MAXF /1,1,4,4,1,4,4,10,1,4,4,10,1,4,10,1,4,10,20,35,56/
      DATA NANGM /1,1,2,2,1,2,2,3,1,2,2,3,1,2,3,1,2,3,4,5,6/
      DATA PI32 /5.5683279968317D+00/
      DATA TOL /1.0D-10/
      DATA LINE,LINEV /8H * * * *,8H       */
      DATA IBLANK /4H    /
C**********************************************************************
C     READ IN MOLECULE INFORMATION (DUMMY,DUMMY,UNITS)                *
C     READ IN UNIQUE CENTERS AND ATOMIC BASIS SETS GROUPED IN SHELLS  *
C     GENERATE ALL NEW CENTERS                                        *
C     SET TABLE ( CENTERS VERSUS TRANSFORMATIONS )                    *
C     SET TABLE ( SHELLS VERSUS TRANSFORMATIONS )                     *
C**********************************************************************
      if (gencon.ne.0) then
        write(6,*)
        endif
      MSU=0
      MCONSU=0
      MCU=0
      NAT=0
c06-09-89 NE=0
      MSHELL=0
      NSHELL=0
      iushll=0
      iprimt = 1
      LOC=0
      LOC1=0
      NGAUSS=0
      DO 5 I=1,5
        NPERSH(I)=0
    5 CONTINUE
      UNITS=1.0D+00/0.529177249D+00
      DO 20 I=1,MXPR
        EX(I)=ZERO
        DO 10 J=1,MRCRU
          CS(J,I)=ZERO
          CP(J,I)=ZERO
          CD(J,I)=ZERO
          CF(J,I)=ZERO
          CG(J,I)=ZERO
          CSINP(J,I)=ZERO
          CPINP(J,I)=ZERO
          CDINP(J,I)=ZERO
          CFINP(J,I)=ZERO
          CGINP(J,I)=ZERO
   10   CONTINUE
   20 CONTINUE

C --------------------------- Convert the basis set data to the
c                             internal format used by input.
  100 CONTINUE
C**********************************************************************
C     ANOTHER SYMMETRY UNIQUE ATOM, INCREMENT MSU.            *
C**********************************************************************
      if (msu.eq.nucent) goto 3000
      MSU=MSU+1
      znuc = charg(msu)
      x = geom(1,msu)
      y = geom(2,msu)
      z = geom(3,msu)
      errcod = frdc1('ATOMS',albl,msu)
      if (errcod.ne.0) then
        write(6,*) 'ERROR: unexpected problem getting one of'
     &             // ' the ''ATOMS'''
        call pabort
        call qabort
        endif
c This puts the label in the ilbl array due to the equivanlence
      clbl = albl
C**********************************************************************
C     CONVERT ANGSTROMS TO BOHR if needed.                            *
C**********************************************************************
      IF(cunits.eq.'ANGSTROMS') then
        X=X*UNITS
        Y=Y*UNITS
        Z=Z*UNITS
        endif
C**********************************************************************
C     CHECK FOR REDUNDANCY OF CENTERS.                                *
C**********************************************************************
      IF(NAT.EQ.0) GO TO 121
      DO 120 IAT=1,NAT
        TEST=(X-C(1,IAT))**2+(Y-C(2,IAT))**2+(Z-C(3,IAT))**2
        IF(TEST.GE.TOL) GO TO 120
        WRITE(IW,8879)
        call qabort
  120 CONTINUE
  121 NAT=NAT+1
      LBLAT(1,NAT)= ilbl(1)
      LBLAT(2,NAT)= ilbl(2)
      NS(NAT)=0
      MS(NAT)=0
      KS(NAT)=NSHELL+1
      ZAN(NAT)=ZNUC
      C(1,NAT)=X
      C(2,NAT)=Y
      C(3,NAT)=Z
C**********************************************************************
C           READ IN BASIS SET ON EACH ATOM.                           *
C**********************************************************************
      ishc = 0
  200 CONTINUE
      if (ishc.eq.nshc(msu)) goto 2000
      ishc = ishc + 1
C
C**********************************************************************
C     READ IN SHELL NUMBER (NOT USED), SHELL TYPE, NUMBER OF          *
C          PRIMITIVES, AND NUMBER OF GENERAL CONTRACTIONS (DEFAULT 1) *
C**********************************************************************
c NSHELL counts the total number of shells
      NSHELL=NSHELL+1
c iushll counts the number of unique shells
      iushll = iushll + 1
      call inttyp('I_TO_C',ctype,types(iushll))
c the number of primitives
      igauss = nprim(iushll)
c no generalized contractions yet
      igncn = 0
clj ityp was originally established with a comparison with CLABEL
clj ipersh was originally read in from the input
      if (ctype.eq.'S') then
        ityp = 16
        ipersh = 0
      else if (ctype.eq.'P') then
        ityp = 17
        ipersh = 0
      else if (ctype.eq.'D') then
        ityp = 18
        ipersh = puream
      else if (ctype.eq.'F') then
        ityp = 19
        ipersh = puream
      else if (ctype.eq.'G') then
        ityp = 20
        ipersh = puream
      else if (ctype.eq.'D6') then
        ityp = 18
        ipersh = 0
      else if (ctype.eq.'F10') then
        ityp = 19
        ipersh = 0
      else if (ctype.eq.'G15') then
        ityp = 20
        ipersh = 0
      else if (ctype.eq.'D5') then
        ityp = 18
        ipersh = 1
      else if (ctype.eq.'F7') then
        ityp = 19
        ipersh = 1
      else if (ctype.eq.'G9') then
        ityp = 20
        ipersh = 1
      else
        write(6,*) 'ERROR: cannot handle angular momentum of ',ctype
        call qabort
        endif
cges  change max # primitives to 35 above
      IF(IGAUSS.GT.35) THEN
        write(6,*) 'ERROR: too many primitives in a contraction'
        call qabort
        END IF
      MCONSU=MCONSU+1
      MCONU=MAX(MCONU,IGAUSS)
      KMIN(NSHELL)=MINF(ITYP)
      KMAX(NSHELL)=MAXF(ITYP)
      KSTART(NSHELL)=NGAUSS+1
      KATOM(NSHELL)=NAT
      KTYPE(NSHELL)=NANGM(ITYP)
      INTYP(NSHELL)=ITYP
      KNG(NSHELL)=IGAUSS
      NGNCN(NSHELL)=IGNCN
CTPH
      IF(IPERSH.EQ.1) THEN
         NTYP=ITYP-15
         NPERSH(NTYP)=1
      ENDIF
clj Search the following contractions to see if generalized contractions
clj are possible (the maximum number of contractions in a generalized
clj contraction is mrcru) the search is done even if gencon == 0. This
clj is to find obvious linear dependencies in the basis set.
      igncnu = 1
clj   iush2 is a temporary counter for the number of the unique shell
      iush2 = iushll
clj   k is the number of the shell on this center
c     write(6,*) 'mrcru:',mrcru
      do 1101 i = ishc+1,nshc(msu)
        iush2 = iush2 + 1
c       write(6,*) 'iushll,iush2:',iushll,iush2
c       write(6,*) '  nprim:',nprim(iushll),nprim(iush2)
c       write(6,*) '  types:',types(iushll),types(iush2)
        if (nprim(iushll).ne.nprim(iush2)) goto 9101
        if (types(iushll).ne.types(iush2)) goto 9101
clj     ipcnt is a temporary counter for the number of the primitive
        ipcnt = iprimt + (i-ishc)*nprim(iushll)
        do 1201 j = 1, nprim(iushll)
c         write(6,*) '  iprimt+j-1,ipcnt+j-1:',iprimt+j-1,ipcnt+j-1
c         write(6,*) '    expon:',expon(iprimt+j-1),expon(ipcnt+j-1)
          if (expon(iprimt+j-1).ne.expon(ipcnt+j-1)) goto 9101
 1201     continue
c We found a generalized contraction
        igncnu = igncnu + 1
c Break out if we have reached the upper bound for generalized contractions
        if (igncnu.eq.mrcru) goto 9101
 1101   continue
c     Break out of the above do loop to here if no more gen cons can be found
 9101 continue
      if (igncnu.gt.nprim(iushll)) then
        write(6,*) 'ERROR: the basis set appears to have'
     &             // ' linear dependencies'
        write(6,*) 'unique center = ',msu
        write(6,*) 'shell on this center = ',ishc
        write(6,*) 'number of primitives = ',nprim(iushll)
        write(6,*) 'number of times these primitives appear = ',igncnu
        call qabort
        endif
clj   reset igncnu to 1 if gencon is 0
      if (gencon.eq.0) then
        igncnu = 1
      else if (igncnu.gt.1) then
        write(6,1100)
     &     'Using a generalized contraction on center ',
     &     msu,' for the shells ',ishc,'-',ishc+igncnu-1,'.'
 1100   format(1X,A,I5,A,I5,A,I5,A)
        endif
      ngncn(nshell) = igncnu
      MSHELL=MSHELL+igncnu
      NGAUSS=NGAUSS+IGAUSS
      NS(NAT)=NS(NAT)+1
      MS(NAT)=MS(NAT)+igncnu
C**********************************************************************
C     GENERAL BASIS SET: READ IN EXPONENTS AND CONTRACTION            *
C     COEFFICIENTS                                                    *
C**********************************************************************
      IF(ITYP.LE.15) CALL BERROR(3)
      K1=KSTART(NSHELL)
      K2=K1+KNG(NSHELL)-1
      DO 230 K=K1,K2
        ex(k) = expon(iprimt)
clj set up the contraction coefficients for this exponent
        do 1102 i = 1, igncnu
          c1(i) = conco(iprimt+(i-1)*nprim(iushll))
 1102     continue
        iprimt = iprimt + 1
        IF(ITYP.EQ.16) THEN
          DO 220 I=1,IGNCNU
            CSINP(I,K)=C1(I)
            CS(I,K)=C1(I)
  220     CONTINUE
        ELSE IF(ITYP.EQ.17) THEN
          DO 222 I=1,IGNCNU
            CPINP(I,K)=C1(I)
            CP(I,K)=C1(I)
  222     CONTINUE
        ELSE IF(ITYP.EQ.18) THEN
          DO 224 I=1,IGNCNU
            CDINP(I,K)=C1(I)
            CD(I,K)=C1(I)
  224     CONTINUE
        ELSE IF(ITYP.EQ.19) THEN
          DO 226 I=1,IGNCNU
            CFINP(I,K)=C1(I)
            CF(I,K)=C1(I)
  226     CONTINUE
        ELSE IF(ITYP.EQ.20) THEN
          DO 228 I=1,IGNCNU
            CGINP(I,K)=C1(I)
            CG(I,K)=C1(I)
  228     CONTINUE
        END IF
  230 CONTINUE

clj If we had a generalized contraction, then the iushll, ishc, and iprimt
clj counters are now wrong, correct these:
      iprimt = iprimt + (igncnu - 1)*nprim(iushll)
      iushll = iushll + igncnu - 1
      ishc = ishc + igncnu - 1
C**********************************************************************
C        CHANGING NORMALIZATION OF D FUNCTIONS TO SUIT PITZER.        *
C        USE PT25 INSTEAD OF PT75 FOR D FUNCTIONS, FOR EXAMPLE.       *
C        PRIMITIVE FUNCTIONS ARE ALWAYS NORMALIZED NOW.               *
C**********************************************************************
      DO 1010 IG=K1,K2
        EE=EX(IG)+EX(IG)
        FACS=PI32/(EE*DSQRT(EE))
        FACP=PT5*FACS/EE
        FACD=PT25*FACS/(EE*EE)
        FACF=PT125*FACS/(EE*EE*EE)
        FACG=PT0625*FACS/(EE*EE*EE*EE)
C
        DO 1005 I=1,IGNCNU
          CS(I,IG)=CS(I,IG)/DSQRT(FACS)
          CP(I,IG)=CP(I,IG)/DSQRT(FACP)
          CD(I,IG)=CD(I,IG)/DSQRT(FACD)
          CF(I,IG)=CF(I,IG)/DSQRT(FACF)
          CG(I,IG)=CG(I,IG)/DSQRT(FACG)
 1005   CONTINUE
 1010 CONTINUE
C**********************************************************************
C        NORMALIZE THE CONTRACTED BASIS FUNCTIONS IF(NORMF.NE.1)      *
C**********************************************************************
      IF(NORMF.EQ.1) GO TO 200
      DO 1060 I=1,IGNCNU
      FACS=ZERO
      FACP=ZERO
      FACD=ZERO
      FACF=ZERO
      FACG=ZERO
      DO 1040 IG=K1,K2
        DO 1040 JG=K1,IG
          EE=EX(IG)+EX(JG)
          FAC=EE*DSQRT(EE)
          DUMS=CS(I,IG)*CS(I,JG)/FAC
          DUMP=PT5*CP(I,IG)*CP(I,JG)/(EE*FAC)
          DUMD=PT25*CD(I,IG)*CD(I,JG)/(EE**2*FAC)
          DUMF=PT125*CF(I,IG)*CF(I,JG)/(EE**3*FAC)
          DUMG=PT0625*CG(I,IG)*CG(I,JG)/(EE**4*FAC)
          IF(IG.EQ.JG) GO TO 1030
          DUMS=DUMS+DUMS
          DUMP=DUMP+DUMP
          DUMD=DUMD+DUMD
          DUMF=DUMF+DUMF
          DUMG=DUMG+DUMG
 1030     FACS=FACS+DUMS
          FACP=FACP+DUMP
          FACD=FACD+DUMD
          FACF=FACF+DUMF
          FACG=FACG+DUMG
 1040 CONTINUE
      DO 1050 IG=K1,K2
        IF(FACS.GT.TOL) CS(I,IG)=CS(I,IG)/DSQRT(FACS*PI32)
        IF(FACP.GT.TOL) CP(I,IG)=CP(I,IG)/DSQRT(FACP*PI32)
        IF(FACD.GT.TOL) CD(I,IG)=CD(I,IG)/DSQRT(FACD*PI32)
        IF(FACF.GT.TOL) CF(I,IG)=CF(I,IG)/DSQRT(FACF*PI32)
        IF(FACG.GT.TOL) CG(I,IG)=CG(I,IG)/DSQRT(FACG*PI32)
 1050 CONTINUE
 1060 CONTINUE
C**********************************************************************
C        GO TO 200 TO TRY TO READ IN ANOTHER SET OF PRIMITIVES.       *
C**********************************************************************
      GO TO 200
 2000 CONTINUE
      ICU=1
      IF (NT.EQ.1) GO TO 100
C**********************************************************************
C     GENERATE EQUIVALENT CENTERS                                     *
C**********************************************************************
      CALL LOCAL(X,Y,Z,XS,YS,ZS)
      XSMAL=XS
      YSMAL=YS
      ZSMAL=ZS
      NAT0=NAT
      DO 2030 IT=1,NT
        IF(IT.EQ.1) GO TO 2030
        NN=9*(IT-1)
        CALL TRANS(T,NN)
        CALL ROT
        DO 2010 IAT=1,NAT
          TEST=(XP-C(1,IAT))**2+(YP-C(2,IAT))**2+(ZP-C(3,IAT))**2
          IF(TEST.LE.TOL) GO TO 2030
 2010   CONTINUE
        ICU=ICU+1
        NUAT=NAT
        NAT=NAT+1
        C(1,NAT)=C(1,NUAT)
        C(2,NAT)=C(2,NUAT)
        C(3,NAT)=C(3,NUAT)
        C(1,NUAT)=XP
        C(2,NUAT)=YP
        C(3,NUAT)=ZP
        NS(NAT)=NS(NAT0)
        MS(NAT)=MS(NAT0)
        KS(NAT)=KS(NAT-1)+NS(NAT0)
        LBLAT(1,NAT)=LBLAT(1,NAT0)
        LBLAT(2,NAT)=LBLAT(2,NAT0)
        ZAN(NAT)=ZAN(NAT0)
c06-09-89  NE=NE+ZAN(NAT)
        NSHELL=NSHELL+NS(NAT)
        MSHELL=MSHELL+MS(NAT)
        NS1=NS(NAT)
        DO 2020 K=1,NS1
          J=KS(NAT0)+K-1
          JJ=KS(NAT)+K-1
          KMIN(JJ)=KMIN(J)
          KMAX(JJ)=KMAX(J)
          KSTART(JJ)=KSTART(J)
          KTYPE(JJ)=KTYPE(J)
          NGNCN(JJ)=NGNCN(J)
CTPH      MSHELL=MSHELL+NGNCN(JJ)
          KNG(JJ)=KNG(J)
          KATOM(JJ)=NAT
          INTYP(JJ)=INTYP(J)
 2020   CONTINUE
 2030 CONTINUE
      IF (ICU.GT.MCU) MCU=ICU
C**********************************************************************
C        GO TO 100 TO TRY TO READ IN ANOTHER ATOM.                    *
C**********************************************************************
      GO TO 100
C**********************************************************************
C        NO MORE ATOMS.                                               *
C**********************************************************************
 3000 CONTINUE
C**********************************************************************
C           SET MCU FOR CASE WHERE NO ATOMS ARE RELATED               *
C**********************************************************************
      IF (MCU.LE.0) MCU=1
      if ((normp.ne.1).or.(normf.ne.1)) write(6,*)
      IF(NORMP.NE.1) WRITE(IW,8864)
      IF(NORMF.NE.1) WRITE(IW,8850)
C**********************************************************************
C     FORM TRANSFORMATION TABLES FOR ATOMS AND SHELLS.                *
C**********************************************************************
      DO 3050 IAT=1,NAT
         X=C(1,IAT)
         Y=C(2,IAT)
         Z=C(3,IAT)
         NS1=KS(IAT)-1
         NS2=NS(IAT)
         CALL LOCAL(X,Y,Z,XS,YS,ZS)
         XSMAL=XS
         YSMAL=YS
         ZSMAL=ZS
         DO 3040 IT=1,NT
            NN=9*(IT-1)
            CALL TRANS(T,NN)
            CALL ROT
            DO 3010 I=1,NAT
               TEST=(XP-C(1,I))**2+(YP-C(2,I))**2+(ZP-C(3,I))**2
               IF(TEST.GT.TOL) GO TO 3010
               ICTR=I
               GO TO 3020
 3010       CONTINUE
 3020       ICT(IAT,IT)=ICTR
            NS3=KS(ICTR)-1
            DO 3030 ISH=1,NS2
 3030       NEWSH(NS1+ISH,IT)=NS3+ISH
 3040    CONTINUE
 3050 CONTINUE
      NTWD=(NT+3)/4
      IF (NTWD.GT.30) then
        write(6,*) 'ERROR: NTWD TOO BIG...PWS'
        call qabort
        endif
      DO 3090 I=1,NSHELL
         DO 3060 IT=1,NTWD
 3060    INDOUT(IT)=0
         DO 3070 IT=1,NT
            INDIN(IT)=0
 3070    INDIN(IT)=NEWSH(I,IT)
         CALL ISOIN
         DO 3080 IT=1,NTWD
 3080    ISO(I,IT)=INDOUT(IT)
 3090 CONTINUE
C**********************************************************************
C     NPRINT=7   ... PUNCH BASIS SET INPUT DATA + MO'S                *
C     NPRINT=1   ... EXTRA PRINTING FOR BASIS SET + SYMMETRY.         *
C     NPRINT=0   ... NORMAL PRINTING                                  *
C**********************************************************************
      if (prcntr('IS_ON SHELLS').eq.1) then
      WRITE(IW,8884)
      DO 4000 IAT=1,NAT
         WRITE(IW,8883) LBLAT(1,IAT),LBLAT(2,IAT),ZAN(IAT),C(1,IAT),
     1   C(2,IAT),C(3,IAT),MS(IAT)
         NS1=KS(IAT)
         NS2=NS1+NS(IAT)-1
CTPH     WRITE(IW,8871) (LABEL(INTYP(ISH)),ISH=NS1,NS2)
         I1=0
         DO 3998 ISH=NS1,NS2
           DO 3997 IGCN=1,NGNCN(ISH)
             I1=I1+1
CTPH   IF IMPURE  (IPERSH=1) ADD 6 TO INTYP)
CTPH         JNTYP(I1)=INTYP(ISH)+5*NPERSH(INTYP(ISH)-15)-15
             JNTYP(I1)=INTYP(ISH)+5*(1-NPERSH(INTYP(ISH)-15))-15
             IF(I1.EQ.3) THEN
               WRITE(IW,8871) (JLABEL(JNTYP(I)),I=1,3)
               I1=0
             ENDIF
 3997      CONTINUE
 3998    CONTINUE
         IF(I1.GT.0) THEN
           WRITE(IW,8871) (JLABEL(JNTYP(I)),I=1,I1)
         ENDIF
 4000 CONTINUE
      endif
      if (prcntr('IS_ON BRIEF').eq.0) then
      WRITE(IW,8882)
      DO 4090 IAT=1,NAT
        DO 4010 IT=1,NT
          IF(ICT(IAT,IT).GT.IAT) GO TO 4090
 4010   CONTINUE
        WRITE(IW,8848) LBLAT(1,IAT),LBLAT(2,IAT)
        NS1=KS(IAT)
        NS2=NS1+NS(IAT)-1
        ISH1=0
        DO 4080 ISH=NS1,NS2
          WRITE(IW,8855)
          I1=KSTART(ISH)
          I2=I1+KNG(ISH)-1
          IGNCNU=NGNCN(ISH)
          ISH1=ISH1+IGNCNU
          ITYP=INTYP(ISH)
          DO 4070 IG=I1,I2
            GO TO (4020,4020,4030,4050,4020,4030,4050,4040,4020,4030,
     1        4050,4040,4020,4030,4040,4020,4030,4040,4055,4058,4040),
     2        ITYP
 4020       DO 4024 I=1,IGNCNU
 4024       C1(I)=CSINP(I,IG)
            GO TO 4060
 4030       DO 4034 I=1,IGNCNU
 4034       C1(I)=CPINP(I,IG)
            GO TO 4060
 4040       DO 4044 I=1,IGNCNU
 4044       C1(I)=CDINP(I,IG)
            GO TO 4060
 4050       C1(1)=CSINP(1,IG)
            C1(2)=CPINP(1,IG)
            WRITE(IW,8849) ISH1,LABEL(ITYP),IG,EX(IG),C1(1),C1(2)
            GO TO 4070
 4055       DO 4057 I=1,IGNCNU
 4057       C1(I)=CFINP(I,IG)
            GO TO 4060
 4058       DO 4059 I=1,IGNCNU
 4059       C1(I)=CGINP(I,IG)
            GO TO 4060
 4060       WRITE(IW,8881) ISH1,LABEL(ITYP),IG,EX(IG),(C1(I),I=1,IGNCNU)
 4070     CONTINUE
 4080   CONTINUE
 4090 CONTINUE
      endif
      if (prcntr('IS_ON DEBUG').eq.0) goto 5000
      WRITE(IW,8860)
      WRITE(IW,8870)
      IMAX=0
 4100 IMIN=IMAX+1
      IMAX=IMAX+15
      IF(IMAX.GT.NT) IMAX=NT
      IMAX1=IMAX+1
      WRITE(IW,8869)
      WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
      WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
      WRITE(IW,8867) LINEV,(I,I=IMIN,IMAX)
      WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
      WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
      WRITE(IW,8868) (LINE ,I=IMIN,IMAX1)
      DO 4200 IAT=1,NAT
        WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
        WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
        WRITE(IW,8866) IAT,(ICT(IAT,I),I=IMIN,IMAX)
        WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
        WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
        WRITE(IW,8868) (LINE ,I=IMIN,IMAX1)
 4200 CONTINUE
      IF(IMAX.LT.NT) GO TO 4100
      WRITE(IW,8860)
      WRITE(IW,8865)
      IMAX=0
 4300 IMIN=IMAX+1
      IMAX=IMAX+15
      IF(IMAX.GT.NT) IMAX=NT
      IMAX1=IMAX+1
      WRITE(IW,8869)
      WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
      WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
      WRITE(IW,8867) LINEV,(I,I=IMIN,IMAX)
      WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
      WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
      WRITE(IW,8868) (LINE ,I=IMIN,IMAX1)
      DO 4400 ISH=1,NSHELL
        WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
        WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
        WRITE(IW,8866) ISH,(NEWSH(ISH,I),I=IMIN,IMAX)
        WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
        WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
        WRITE(IW,8868) (LINE ,I=IMIN,IMAX1)
 4400 CONTINUE
      IF(IMAX.LT.NT) GO TO 4300
      WRITE(IW,8860)
      WRITE(IW,8851)
      IMAX=0
 4500 IMIN=IMAX+1
      IMAX=IMAX+15
      IF(IMAX.GT.NT) IMAX=NT
      IMAX1=IMAX+1
      WRITE(IW,8869)
      WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
      WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
      WRITE(IW,8867) LINEV,(I,I=IMIN,IMAX)
      WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
      WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
      WRITE(IW,8868) (LINE ,I=IMIN,IMAX1)
      WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
      WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
      WRITE(IW,8867) LINEV,(INVT(I),I=IMIN,IMAX)
      WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
      WRITE(IW,8868) (LINEV,I=IMIN,IMAX1)
      WRITE(IW,8868) (LINE ,I=IMIN,IMAX1)
      IF(IMAX.LT.NT) GO TO 4500
 5000 CONTINUE
      IF(NAT.LE.MXAT) GO TO 5100
      WRITE(IW,8876)
      call qabort
 5100 IF(NSHELL.LE.MXSH) GO TO 5200
      WRITE(IW,8875)
      call qabort
 5200 CONTINUE
C5200 IF(LOC.LE.100) GO TO 5300
C5200 IF(LOC.LE.120) GO TO 5300
C     WRITE(IW,8874)
C     STOP
 5300 IF(NGAUSS.GT.MXPR) THEN
        WRITE(IW,8873)
        call qabort
      END IF
C
C   NOW COUNT NUMBER OF AOS = NBFAO = LOC,  SOS = NBASIS = NBFSO = LOC1
C
      DO 5400 I=1,NSHELL
        ITYP=INTYP(I)
        KLOC(I)=LOC+1
        LOC=LOC+NBFS(ITYP)*NGNCN(I)
        NTYP=ITYP-15
        IF(NPERSH(NTYP).EQ.1) THEN
          LOC1=LOC1+NBFSP(ITYP)*NGNCN(I)
        ELSE
          LOC1=LOC1+NBFS(ITYP)*NGNCN(I)
        ENDIF
 5400 CONTINUE
      MSFU=MCONSU
      NUM=LOC
      MSFRU=LOC1
      NX=NUM*(NUM+1)/2
      RETURN
 8848 FORMAT(/,2A4)
 8849 FORMAT(10X,I3,1X,A4,3X,I3,1X,2F15.6,3H  (,F12.6,3H  ),
     1  F15.6,3H  (,F12.6,3H  ))
 8850 FORMAT(' THE CONTRACTED BASIS FUNCTIONS ARE NOW NORMALIZED TO',
     1  ' UNITY')
 8851 FORMAT(/,1X,5H*****,' INVERSE TRANSFORMATIONS ',5H*****,/)
 8855 FORMAT(/)
 8859 FORMAT(10X,I5)
 8860 FORMAT(1H1)
 8864 FORMAT(' THE CONTRACTED PRIMITIVE FUNCTIONS HAVE BEEN ',
     1  'UNNORMALIZED')
 8865 FORMAT(/' ***** TRANSFORMATION TABLE OF SHELLS *****'/
     1  30X,' ROWS ARE SHELLS'/30X,' COLUMNS ARE SYMMETRY OPERATIONS')
 8866 FORMAT(1X,16(I4,3X,1H*))
 8867 FORMAT(1X,A8,15(I4,3X,1H*))
 8868 FORMAT(1X,16A8)
 8869 FORMAT(//)
 8870 FORMAT(/' ***** TRANSFORMATION TABLE OF ATOMS *****'/
     1  30X,' ROWS ARE ATOMS'/30X,' COLUMNS ARE SYMMETRY OPERATIONS')
 8871 FORMAT(60X,3A4)
 8872 FORMAT(' EXCESSIVE CONTRACTION INTO SHELL')
 8873 FORMAT(' EXCESSIVE NUMBER OF CONTRACTED PRIMITIVE FUNCTIONS')
 8874 FORMAT(' EXCESSIVE NUMBER OF BASIS FUNCTIONS')
 8875 FORMAT(' EXCESSIVE NUMBER OF SHELLS')
 8876 FORMAT(' EXCESSIVE NUMBER OF ATOMS')
 8878 FORMAT(38H UNRECOGNIZABLE BASIS FUNCTION TYPE = ,A4)
 8879 FORMAT(50H HONDO DOESN'T TAKE THE DANGER OUT OF BEING CLOSE./
     1  19H REDUNDANT CENTERS.)
cges change format to allow 10 contractions
c8881 FORMAT(I16,1X,A4,I6,1X,F15.6,6F12.6)
 8881 FORMAT(I16,1X,A4,I6,1X,F15.6,10F12.6)
 8882 FORMAT(/,78(1H*),//,8X,' CONTRACTED PRIMITIVE FUNCTIONS',//,
     1  '  ATOM       ','SHELL TYPE PRIM      EXPONENTS   ',
     2  1X,'CONTRACTION COEFFICIENTS')
 8883 FORMAT(2A4,2X,F5.1,3(F12.7,3X),I3)
 8884 FORMAT(/,38(2H**),/,'  ATOM    ','ATOMIC',16X,
     1  'COORDINATES',16X,'NUMBER OF',/,10X,'CHARGE',7X,'X',13X,
     1  'Y',14X,'Z',6X,'SHELLS',/,38(2H**))
 8886 FORMAT(5X,E15.9,3E20.12)
 8896 FORMAT(E10.6,7E10.6)
 8887 FORMAT(2A4,2X,F5.1,3F20.12)
 8888 FORMAT(6X,A4,3I5)
 8889 FORMAT(20X,3E20.12)
      END
