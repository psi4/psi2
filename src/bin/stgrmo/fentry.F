      subroutine fentry(cc,ia,maxcor)
C**********************************************************************
C* 10/09/91        PROGRAM STGRMO
C**********************************************************************
C   AN SCF STABILITY PROGRAM FOR OPEN-SHELL SCF SYSTEMS
C**********************************************************************
C   Last updated on January 10, 1995 by Evan Bolton                   *
C   Reason:  To put in PSI 2.0 input parsing.  See comments Ceeb      *
C**********************************************************************
C   Last updated on August 25, 1993 by Dr. Yukio Yamaguchi           *
C   Dimension of `istor', `cstor' and `eigsv' were increased to 50   *
C*********************************************************************
C   Last updated on August 21, 1993 by Dr. Yukio Yamaguchi           *
C   Dimension of `ioff' was increased to (256)*(256+1)/2 + 1 = 32897 *
C*********************************************************************
C   Last updated on July 30, 1993 by Dr. Yukio Yamaguchi             *
C   ``call tstart'' and ``call initmf'' were moved                   *
C     after ``call psinit''                                          *
C*********************************************************************
C   Last updated on October 09, 1991 by Dr. Yukio Yamaguchi          *
C   Modification for the Unix Version                                *
C*********************************************************************
C***LAST UPDATED ON MARCH 28, 1990   BY YUKIO YAMAGUCHI              *
C*********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension cc(maxcor),ia(maxcor*2)
      CHARACTER*10 TITLE
      CHARACTER*3 ISMAT
      DIMENSION NDEG(8),LORD(256),IPZSYM(256),ISYM(8)
Ceeb Needed For Psi2.0 Input
      integer errcod, inold
      integer fwkclr, fwkadd, frdi, frdboo, frdc
      common /instyl/ inold
      character*80    method, symm
Ceeb
      COMMON/BASIS/NBASIS,NTRI
      COMMON/CALIF/LPARA(1024),APARA(1024)
      COMMON/CONVS/IQCONV,IDCONV,IECONV
      COMMON/COUPL/ALPC(15),BETC(15)
      COMMON/ENRGY/ESCF,ENUC
      COMMON/FUNCS/NTYPES,NTYPEP
      COMMON/GVBAB/A(10,10),B(10,10)
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4
      COMMON/MFSEC/MASTER,NSECT
      COMMON/POSIT/IABP(500)
      COMMON/ROOTS/NROOT,NPAIR
      COMMON/SDATA/ISMAT(8,8)
      COMMON/SIGNS/IOFF(32897),IPRNT
      COMMON/START/FOCC(10),NSORB(10),NSTR(10),NEND(10),MOTYP(256)
      common/stori/istor(50),maxstr
      common/storr/cstor(50),eigsv(50)
      COMMON/SYMTR/ISTYPE,MOSYM(256),MONUM(256)
      COMMON/TAPES/ITAP35,ITAP91,ITAP92
      COMMON/CI101/IOCC,JOCC,KOCC
      COMMON/CI102/NIND,NDEP,NDTRI
      COMMON/CI103/MXITER,MXVEC,MXTRI
      COMMON/CI104/NVEC
      DATA ONE,QURT,HALF / 1.0D+00 , 0.25D+00 , 0.5D+00 /
    1 FORMAT(//,2X,' STGRMO: An SCF Stability Program ',
     1 'in the MO basis  (Version of 01/14/95)',/,11X,
     2 'Written by Dr. Yukio ',
     3 'Yamaguchi',//,11X, 'Modifications for PSI 2.0 use by ',
     4 'Evan Bolton'///)
    2 FORMAT(14I5)
    3 FORMAT(//,2X,' PARAMETERS'/
     * 2X,' NBASIS = ',I8/
     * 2X,' NTRI   = ',I8/
     * 2X,' IOCC   = ',I8/
     * 2X,' JOCC   = ',I8/
     * 2X,' KOCC   = ',I8/
     * 2X,' NTYPES = ',I8/
     * 2X,' NTYPEP = ',I8/
     * 2X,' NROOT  = ',I8/
     * 2X,' NPAIR  = ',I8/
     * 2X,' MXITER = ',I8/
     * 2X,' MXVEC  = ',I8/
     * 2X,' IQCONV = ',I8/
     * 2X,' IDCONV = ',I8/
     * 2X,' IECONV = ',I8/
     * 2X,' IAUTO  = ',I8/
     * 2X,' IDIAG  = ',I8/
     * 2X,' IDAVL  = ',I8/
     * 2X,' IDAVX  = ',I8/
     * 2X,' IBPRT  = ',I8/
     * 2X,' IPRNT  = ',I8/)
    4 FORMAT(//,2X,' ALPA MATRIX'/)
    5 FORMAT(//,2X,' BETA MATRIX'/)
    6 FORMAT(A10)
    7 FORMAT(//,2X,' NO SYMMETRY INPUT ---> STOP !'/)
    8 FORMAT(//,2X,' SCF EIGENVECTOR (SO BASIS) MATRIX'/)
    9 FORMAT(//,2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)
   10 FORMAT(//,2X,' REQUIRED MEMORY EXCEEDS MAXCOR'/
     1          2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)
   11 FORMAT(//,2X,' NIND EXCEEDS MAXSIZ'/
     1          2X,' NIND  = ',I10,5X,' MAXSIZ = ',I10/)
C
      call psinit('APPEND')
      CALL TSTART(6)
      CALL INITMF(1)
C::::::::::::::::::::::::::::::::::::::::
C:::DEFINE MACHINE DEPENDENT VARIABLES:::
C::::::::::::::::::::::::::::::::::::::::
      MAXBAS=256
      MAXSIZ=MAXBAS*(MAXBAS+1)/2
      MAXVEC=500
      maxstr=50
C
C   ITAP35 = TWO ELECTRON MO INTEGRALS
C   ITAP91 = WORKING TAPE
C   ITAP92 = WORKING TAPE
C
      ITAPE3=3
      INPUT=5
      ITAPE6=6
      ITAP35=35
      ITAP91=91
      ITAP92=92
      CALL LOCATE(INPUT,'# STABIL #',IERR)

Ceeb  Check to see if using old style input
      if (ierr.ne.0) inold = 0
      if (ierr.eq.0) inold = 1

      if (inold.eq.1) then
        write(6,*) 'WARNING: using old style ''# INPUT ##'' and ''
     1# STABIL ##'''
      endif
Ceeb
C
      WRITE(6,1)
      WRITE(3,1)
      write(6,*) '  maxcor = ',maxcor
C
C   Parameters from the master file
      nbasis=lpara(3)
      ntri=lpara(4)
      if(nbasis.gt.maxbas) then
        write(6,*) '  Nbasis exceeds Maxbas -----> Stop ! '
        write(6,*) '  Nbasis = ',nbasis,' Maxbas = ',maxbas
        stop
      end if
      if(ntri.gt.maxsiz) then
        write(6,*) '  Ntri exceeds Maxsiz -----> Stop ! '
        write(6,*) '  Ntri = ',ntri,' Maxsiz = ',maxsiz
        stop
      end if
C
c   Form the IOFF array
      IOFF(1)=0
      DO 101 I=1,NTRI
      IOFF(I+1)=IOFF(I)+I
  101 CONTINUE

C
Ceeb  Read in Old style input
      if (inold.eq.1) then
C
C   READIN PARAMETERS
      READ(5,2) NROOT,NPAIR,MXITER,MXVEC,IQCONV,IDCONV,IECONV,
     1          IDIAG,IDAVL,IDAVX,IDUM1,IDUM2,IBPRT,IPRNT
      endif
Ceeb  Finished reading in Old Style input

Ceeb  Read in New Style input
      if (inold.eq.0) then
        errcod = fwkclr()
        errcod = fwkadd(5,'STABILITY')
        errcod = fwkadd(5,'DEFAULT')

        errcod = frdi('NROOT',nroot)
        errcod = frdi('NPAIR',npair)
        errcod = frdi('MAXITER',mxiter)
        errcod = frdi('MAXVEC',mxvec)

        errcod = frdi('QCONVERGENCE',iqconv)
        errcod = frdi('DCONVERGENCE',idconv)
        errcod = frdi('ECONVERGENCE',ieconv)

        errcod = frdi('NVECTOR',ibprt)
        errcod = frdi('PRINT',iprnt)

        errcod = frdboo('SINGLET',isngl)
        errcod = frdboo('TRIPLET',itrpl)

Ceeb  For different types of Diagonalization Routines
        errcod = frdc('DMETHOD', method)

        if (method.eq.'RSP') then
           idiag = 1
        else if (method.eq.'DAVIDSON') then
           idavl = 1
        else if (method.eq.'DLARGE') then
           idavx = 1
        endif

Ceeb  This determines whether any method was requested
Ceeb    If no method requested we will automatically set it internally

        if (idiag.eq.0.and.idavl.eq.0.and.idavx.eq.0)
     +    iauto = 1

Ceeb  For now both Singlet and Triplet formalisms are not to be used
        isngl = 0
        itrpl = 0
      endif
Ceeb  Finished reading in New Style input

Ceeb  This sets up the Defaults properly after reading the input
      IF(NROOT.LE.0) NROOT=10
C*    NROOT=MIN0(NROOT,20)
      IF(NPAIR.LE.0) NPAIR=20
      IF(MXITER.LE.0) MXITER=NROOT*20
      IF(MXVEC.LE.0) MXVEC=NROOT*20
      if(mxvec.gt.maxvec) then
        write(6,*) '  Mxvec exceeds Maxvec -----> Stop ! '
        write(6,*) '  Increase dimension of iabp(maxvec) array '
        write(6,*) '  Mxvec = ',mxvec,' Maxvec = ',maxvec
        stop
      end if
       if(nroot.gt.maxstr) then
        write(6,*) '  Nroot exceeds Maxstr -----> Stop ! '
        write(6,*) '  Number of initial guess vectors '
        write(6,*) '    is limited to maxstr '
        write(6,*) '  Nroot = ',nroot,' Maxstr = ',maxstr
      end if
      MXTRI=MXVEC*(MXVEC+1)/2
      IF(IQCONV.EQ.0) IQCONV=6
      IF(IDCONV.EQ.0) IDCONV=6
      IF(IECONV.EQ.0) IECONV=6

Ceeb  Commented out in favor of auto setting it later
Ceeb     after NIND is known
Ceeb
Ceeb  IDIAG=1
Ceeb  IF(IDAVL.NE.0) IDIAG=0
Ceeb  IF(IDAVX.NE.0) IDIAG=0

C
C   READIN SCF INFORMATION
      IOCC=LPARA(7)
      KOCC=LPARA(8)
      JOCC=LPARA(9)
      NTYPES=LPARA(18)
      NST=LPARA(30)
      NTYPEP=NTYPES+1
      ENUC=APARA(1)
      ESCF=APARA(2)
      DO 102 I=1,10
      NSORB(I)=LPARA(I+40)
      FOCC(I)=APARA(I+40)
  102 CONTINUE
      DO 103 I=1,15
      ALPC(I)=APARA(I+10)
      BETC(I)=APARA(I+25)
  103 CONTINUE
      CALL ZERO(A,100)
      CALL ZERO(B,100)
      IJ=0
      DO 104 I=1,NTYPES
      DO 104 J=1,I
      IJ=IJ+1
      A(I,J)= (ONE-APARA(10+IJ))*FOCC(I)*FOCC(J)*HALF
      B(I,J)=-(ONE-APARA(25+IJ))*FOCC(I)*FOCC(J)*QURT
      IF(I.EQ.0) GO TO 104
      A(J,I)=A(I,J)
      B(J,I)=B(I,J)
  104 CONTINUE
      IF(IBPRT.LE.0) IBPRT=MIN0(NBASIS,JOCC+10)
      WRITE(6,3) NBASIS,NTRI,IOCC,JOCC,KOCC,NTYPES,NTYPEP,NROOT,NPAIR,
     1 MXITER,MXVEC,IQCONV,IDCONV,IECONV,IAUTO,IDIAG,IDAVL,IDAVX,
     2 IBPRT,IPRNT
      WRITE(6,4)
C*    CALL PRINT(ALPC,15,NTYPEP,6)
      CALL MATOUT(A,10,10,NTYPEP,NTYPEP,6)
      WRITE(6,5)
C*    CALL PRINT(BETC,15,NTYPEP,6)
      CALL MATOUT(B,10,10,NTYPEP,NTYPEP,6)
C
C   FIND THE SYMMETRY TYPE
      if (inold.eq.1) then
        CALL LOCATE(INPUT,'# INPUT ##',IERR)
        READ(INPUT,6) TITLE
        READ(INPUT,6) TITLE
        READ(INPUT,6) TITLE
        ISTYPE=0
        IF(TITLE.EQ.'C1        ')  ISTYPE=1
        IF(TITLE.EQ.'CS        ')  ISTYPE=2
        IF(TITLE.EQ.'CI        ')  ISTYPE=3
        IF(TITLE.EQ.'CN       2')  ISTYPE=4
        IF(TITLE.EQ.'CNH      2')  ISTYPE=5
        IF(TITLE.EQ.'CNV      2')  ISTYPE=6
        IF(TITLE.EQ.'DN       2')  ISTYPE=7
        IF(TITLE.EQ.'DNH      2')  ISTYPE=8

        IF(ISTYPE.EQ.0) THEN
          WRITE(6,7)
          STOP
        END IF
      endif
Ceeb  End of Read in Symmetry for old Style input

Ceeb  Read in symmetry for new style input
      if (inold.eq.0) then
        errcod = frdc('SYMMETRY',symm)
        if (errcod.ne.0) then
          write(6,*) 'ERROR: couldn''t read ''SYMMETRY'' from input'
          write(6,*) ' '
          call mabort
        endif

        if (symm.eq.'C1') then
          istype = 1
        else if (symm.eq.'CS') then
          istype = 2
        else if (symm.eq.'CI') then
          istype = 3
        else if (symm.eq.'C2') then
          istype = 4
        else if (symm.eq.'C2H') then
          istype = 5
        else if (symm.eq.'C2V') then
          istype = 6
        else if (symm.eq.'D2') then
          istype = 7
        else if (symm.eq.'D2H') then
          istype = 8
        else
          write(6,*) 'ERROR: bad value for ''symmetry'''
          call qabort
        endif
      endif
Ceeb  End of Read in Symmetry for new Style input

      WRITE(3,*) ' ISTYPE = ',ISTYPE
C
      DO 105 I=1,NST
      NDEG(I)=LPARA(I+30)
  105 CONTINUE
      write(3,*) '  Ndeg'
      WRITE(3,*) (NDEG(I),I=1,NST)
      CALL MREAD(LORD,4)
      write(3,*) '  Lord'
      WRITE(3,*) (LORD(I),I=1,NBASIS)
      MM=0
      DO 107 I=1,NST
      NN=NDEG(I)
      IF(NN.EQ.0) GO TO 107
      DO 106 J=1,NN
      MM=MM+1
      IPZSYM(MM)=I
  106 CONTINUE
  107 CONTINUE
      write(3,*) ' Ipzsym'
      WRITE(3,*) (IPZSYM(I),I=1,NBASIS)
      CALL IZERO(ISYM,8)
      CALL IZERO(MOSYM,nbasis)
      CALL IZERO(MONUM,nbasis)
      DO 108 I=1,NBASIS
      MM=LORD(I)
      NN=IPZSYM(MM)
      ISYM(NN)=ISYM(NN)+1
      MOSYM(I)=NN
      MONUM(I)=ISYM(NN)
  108 CONTINUE
      write(3,*) '  Mosym'
      WRITE(3,*) (MOSYM(I),I=1,NBASIS)
      write(3,*) '  Monum'
      WRITE(3,*) (MONUM(I),I=1,NBASIS)
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C  DYNAMIC ALLOCATION
C  1 : NIJ       (NTRI)
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IC1=1
      IA1=1
      ICKEEP=IC1+NTRI
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C   READ IN AO-MO EIGENVECTORS
      IC2=ICKEEP
      IC3=IC2+NBASIS
      IC4=IC3+NBASIS
      ICMAX=IC4+NBASIS*NBASIS
      CALL MREAD(CC(IC2),17)
      CALL MREAD(CC(IC3),18)
      CALL MREAD(CC(IC4),19)
C*    IF(IPRNT.LE.2) GO TO 301
      WRITE(6,8)
      CALL EIGOUT(CC(IC4),CC(IC2),CC(IC3),NBASIS,NBASIS,NBASIS,IBPRT,6)
C
  301 CONTINUE
C   FORM REGISTER ARRAYS
      WRITE(3,21)
   21 FORMAT(//,2X,' NOW YOU ARE IN REGIST'/)
C.................NIJ.....
      CALL REGIST(IA(IA1))
      IF(NIND.GT.MAXSIZ) GO TO 499

Ceeb  Time to automatically set the type of Diagonalization routine

581   format(//,2x,'DEFAULT METHOD OF DIAGONALIZATION  =  RSP')
582   format(//,2x,'DEFAULT METHOD OF DIAGONALIZATION  =  DAVIDSON')
583   format(//,2x,'DEFAULT METHOD OF DIAGONALIZATION  =  DLARGE')

      if (iauto.eq.1.and.nind.le.600) then
        idiag = 1
        write(6,581)
        write(3,581)
      else if (iauto.eq.1.and.nind.gt.600.and.nind.le.2500) then
        idavl = 1
        write(6,582)
        write(3,582)
      else if (iauto.eq.1.and.nind.gt.2500) then
        idavx = 1
        write(6,583)
        write(3,583)
      endif

C
C   CALCULATE SCF ENERGY FOR A TEST
      WRITE(3,22)
   22 FORMAT(//,2X,' NOW YOU ARE IN E0CAL'/)
      IC2=ICKEEP
      IC3=IC2+NTRI
      IA3=IC3+IC3-1
      ICMAX=IC3+MAXBF2
      maxuse=icmax
      WRITE(3,9) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 399
C................H       LBLI    BUFI....
      CALL E0CAL(CC(IC2),IA(IA3),CC(IC3))

Ceeb  Flush io buffer before calculation of stability
      call fioflu(6)
      call fioflu(3)
C
C   FORM THE A MATRIX
      WRITE(3,23)
   23 FORMAT(//,2X,' NOW YOU ARE IN AMAT'/)
      IC2=ICKEEP
      IC3=IC2+NDTRI
      IC4=IC3+NTRI*NTYPEP
      IA4=IC4+IC4-1
      IC5=IC4+NTRI
      IA5=IC5+IC5-1
      ICMAX=IC5+MAXBF2
      maxuse=max0(icmax,maxuse)
      WRITE(3,9) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 399
C...............NIJ     C       ZETA    IND     LBLI    BUFI....
      CALL AMAT(IA(IA1),CC(IC2),CC(IC3),IA(IA4),IA(IA5),CC(IC5))
C
C   DIAGONALIZE THE MO HESSIAN DIRECTLY
C   USING THE SUBROUTINE "RSP"
      IF(IDIAG.EQ.0) GO TO 302
      WRITE(3,24)
   24 FORMAT(//,2X,' NOW YOU ARE IN STDIAG'/)
      IC2=ICKEEP
      IC3=IC2+NDTRI
      IC4=IC3+NIND
      IC5=IC4+NIND*NIND
      IC6=IC5+NIND
      IC7=IC6+NIND
      IA7=IC7+IC7-1
      IC8=IC7+MAXBF2
      IC9=IC8+NIND
      IA9=IC9+IC9-1
      IC10=IC9+NIND
      IA10=IC10+IC10-1
      ICMAX=IC10+NIND
      maxuse=max0(icmax,maxuse)
      WRITE(3,9) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 303
C.................NIJ     C       EIGC    EIVC    FV1     FV2.....
      CALL STDIAG(IA(IA1),CC(IC2),CC(IC3),CC(IC4),CC(IC5),CC(IC6),
C.................LBLI    BUFI    EIVT    NEVT    NORD.....
     1            IA(IA7),CC(IC7),CC(IC8),IA(IA9),IA(IA10))
      GO TO 500
C
C   DIAGONALIZE THE MO HESSIAN FOR SEVERAL ROOTS
C   USING THE DAVIDSON-LIU METHOD WITHOUT WORKING TAPES
  302 CONTINUE
      IF(IDAVL.EQ.0) GO TO 304
  303 CONTINUE
      WRITE(3,25)
   25 FORMAT(//,2X,' NOW YOU ARE IN DAVLIU'/)
      IC2=ICKEEP
      IC3=IC2+NDTRI
      IC4=IC3+MXVEC
      IC5=IC4+NIND*MXVEC
      IC6=IC5+MXVEC*MXVEC
      IC7=IC6+MXTRI
      IC8=IC7+MXVEC
      IC9=IC8+MXVEC
      IC10=IC9+MXTRI
      IC11=IC10+MXVEC
      IC12=IC11+NIND*MXVEC
      IC13=IC12+NIND
      IC14=IC13+NIND
      ICMAX=IC14+NIND
      maxuse=max0(icmax,maxuse)
      WRITE(3,9) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 305
C.................C       EIG     B       VEC     G       FV1.....
      CALL DAVLIU(CC(IC2),CC(IC3),CC(IC4),CC(IC5),CC(IC6),CC(IC7),
C.................FV2     CC      ALPA     AB       Q        Z........
     1            CC(IC8),CC(IC9),CC(IC10),CC(IC11),CC(IC12),CC(IC13),
C.................D........
     2            CC(IC14))
C
      WRITE(3,26)
   26 FORMAT(//,2X,' NOW YOU ARE IN STABIL'/)
      IC2=ICKEEP
      IC3=IC2+NDTRI
      IC4=IC3+MXVEC
      IC5=IC4+NIND*MXVEC
      IC6=IC5+MXVEC*MXVEC
      IC7=IC6+NIND*NROOT
      IC8=IC7+NIND
      IA8=IC8+IC8-1
      IC9=IC8+NIND
      IA9=IC9+IC9-1
      ICMAX=IC9+NIND
      maxuse=max0(icmax,maxuse)
      WRITE(3,9) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 305
C.................NIJ     EIG     B       VEC     EVEC    EIVT....
      CALL STABIL(IA(IA1),CC(IC3),CC(IC4),CC(IC5),CC(IC6),CC(IC7),
C.................NEVT    NORD....
     1            IA(IA8),IA(IA9))
      GO TO 500
C
C   DIAGONALIZE THE MO HESSIAN FOR SEVERAL ROOTS
C   USING THE DAVIDSON-LIU METHOD WITH WORKING TAPES
  304 CONTINUE
      IF(IDAVX.EQ.0) GO TO 310
  305 CONTINUE
      WRITE(3,27)
   27 FORMAT(//,2X,' NOW YOU ARE IN DAVLX'/)
      IC2=ICKEEP
      IC3=IC2+NDTRI
      IC4=IC3+MXVEC
      IC5=IC4+MXVEC*MXVEC
      IC6=IC5+MXTRI
      IC7=IC6+MXVEC
      IC8=IC7+MXVEC
      IC9=IC8+MXTRI
      IC10=IC9+MXVEC
      IC11=IC10+NIND
      IC12=IC11+NIND
      IC13=IC12+NIND
      IC14=IC13+NIND
      ICMAX=IC14+NIND
      maxuse=max0(icmax,maxuse)
      WRITE(3,9) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 399
C.................C       EIG     VEC     G       FV1     FV2.....
      CALL DAVLIX(CC(IC2),CC(IC3),CC(IC4),CC(IC5),CC(IC6),CC(IC7),
C.................CC      ALPA    B        AB       Q        Z........
     1            CC(IC8),CC(IC9),CC(IC10),CC(IC11),CC(IC12),CC(IC13),
C.................D........
     2            CC(IC14))
C
      WRITE(3,28)
   28 FORMAT(//,2X,' NOW YOU ARE IN STABIX'/)
      IC2=ICKEEP
      IC3=IC2+NDTRI
      IC4=IC3+MXVEC
      IC5=IC4+MXVEC*MXVEC
      IC6=IC5+NIND
      IC7=IC6+NIND*NROOT
      IC8=IC7+NIND
      IA8=IC8+IC8-1
      IC9=IC8+NIND
      IA9=IC9+IC9-1
      ICMAX=IC9+NIND
      maxuse=max0(icmax,maxuse)
      WRITE(3,9) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 399
C.................NIJ     EIG     VEC     B       EVEC    EIVT....
      CALL STABIX(IA(IA1),CC(IC3),CC(IC4),CC(IC5),CC(IC6),CC(IC7),
C.................NEVT    NORD....
     1            IA(IA8),IA(IA9))
  310 CONTINUE
      GO TO 500
C
  399 CONTINUE
      WRITE(6,10) ICMAX,MAXCOR
      WRITE(3,10) ICMAX,MAXCOR
      GO TO 500
  499 CONTINUE
      WRITE(6,11) NIND,MAXSIZ
      WRITE(3,11) NIND,MAXSIZ
  500 CONTINUE
C
      write(6,*) '  '
      write(6,*) '  largest core used = ',maxuse
      CALL RCLOSE(MASTER,3)
      CALL TSTOP(6)
C
      STOP
      END
