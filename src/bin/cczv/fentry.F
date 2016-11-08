cges  PROGRAM CCZV
      subroutine fentry(cc,ic,maxcor)
C
C  >>>   LAST UPDATE   april 17 , 1990  <<<
C
C***********************************************************************
C  NOTICE OF PROGRAM MODIFICATION
C      BY:  GUSTAVO E. SCUSERIA (unless otherwise noted)
c
c    date:  Thu Feb 28 12:08:27 EST 1991
c      by:  Curtis L Janssen
c  reason:  converted to new style input
c
c    date:  april 17, 1990
c  reason:  finish implementing ccsd(t) gradients
c           only changes are in initz.f where j1 & j2 are read in
c           every iteration.
C
C    DATE:  november 29, 1989
C  REASON:  split ysdt1.f into ydterm, yfterm & yeterm (opt O2)     
C
C    DATE:  november 28, 1989
C  REASON:  install reformulated z-ccsdt-1 subroutine ysdt1 based on
c           cceg/ccsdt1.f (old version is y3.f)
C
C    DATE:  november 26, 1989
C  REASON:  create term.dat with updated timing information
C
C    DATE:  november 6, 1989
C  REASON:  if option = SDT1 and not SDT1GRAD stops giving an error msg
C
C    DATE:  october 26, 1989
C  REASON:  include CCD option
C
C    DATE:  august 22, 1989
C  REASON:  implement cors and virs in subroutine y3
c
C    DATE:  MAY 10, 1989.
C  REASON:  CHANGE NAME TO CCZV
C
C    DATE:  MAY 01, 1989.
C  REASON:  IMPLEMENT TWRIT
C
C    DATE:  APRIL 28, 1989
C  REASON:  FIX DIMENSION OF ABINTS TO NOV-NOV IN SUBROUTINE ABTERM
C           THIS CAUSED PROBLEMS FOR NO.GT.NV.AND.NV.GT.1
C
C    DATE:  MARCH 31, 1989
C  REASON:  INCLUDE TWRIT
C
C    DATE:  MARCH 28, 1989
C  REASON:  UPDATE SUBROUTINE READ30 (NEW ORDER OF SYMMETRIES)
C           UPDATE SUBROUTINE FROZEN
C           UPDATE SUBROUTINE RDINTS
C           GET RID OF ROUTINE FFINTS
C
C    DATE:  FEBRUARY 23,1989
C  REASON:  IMPLEMENT CORS AND VIRS
C
C    DATE:  NOVEMBER 15,1988
C  REASON:  SUBROUTINE NFINTS: TRADE OFF IO FOR D AND F INTS IN CORE.
C
C    DATE:  OCTOBER 20, 1988
C  REASON:  SYMMETRY-PACKED EMAT.
C           REDUCE CORE REQUIREMENTS BY SPLITTING INTERMEDIATE AND
C           ITERATIVE PARTS (CCSDT1 WILL REQUIRE MORE CORE THAN CCSD)
C
C    DATE:  OCTOBER 24, 1988
C  REASON:  CCSDT-1 WORKING BUT NOT OPTIMIZED
C
C***********************************************************************
C
      IMPLICIT INTEGER (A-Z)
      REAL*8 CC(maxcor),DELZ,ZCONV,TIMLIM,DELZ2,ENUC,ESCF
      REAL*4 EECPU,EETIME,EETCPU
      CHARACTER*32 DATTIM
      CHARACTER*4 fl82,ZOPT,CHAR,OPTION,SORT,DIISFL
      character*80 label
      character*40 corrty, dertyp
      integer i, bDIIS, bRSTR, frdc, frdi, fcnt, errcod
      INTEGER NOSYM(8),NVSYM(8),NFZO(8),NFZV(8),NDOC(8),NUOC(8)
      COMMON/TIEMPO/EECPU,EETIME,EETCPU
      COMMON/INTS/NSGOO,NSGVV,NSGOV,NSHOV,NSLOV,NSLVO
      COMMON/RESTAR/TIMLIM
      COMMON/TAPES/TPBMAT,TPFMAT,TPHMAT,TPRMAT,TPQMAT,TPJMAT,ITAP71,
     .             TPAMAT
C     DIMENSION IC(3400000),I30(200)
cges  DIMENSION IC( 600000),I30(200)
      DIMENSION IC(1),I30(200)
#include <error.h>
cges  EQUIVALENCE (CC,IC)
C
cges  CALL NOUNFL
C
c     MAXCOR= 300000
C
      call psinit('APPEND')
      call ffile(4,'term',0)
c     call flush(4)
c
C      CALL DATETM(DATTIM,EECPU,EETIME,EETCPU)
      CALL TSTART(6)
C
      WRITE(6,6000)
 6000 FORMAT(   56('-'),/,56('-'),/,
     .' THE  CLOSED-SHELL  CCSD  LINEAR  RESPONSE  PROGRAM       ',/,
     .'           WRITTEN  BY  GUSTAVO SCUSERIA'     ,/,
     .'           REFORMULATED VERSION JULY 1988',/,
     .'            VS CCZV 2.0 april 17,  1990  ',/,56('-'),/,56('-') )
C
 
      INPUT = 5
      JOUT  = 6
      FILE7 = 7
      FILE9 = 9
      ITAP30=30
      ITAP60=60
      ITAP61=61
      ITAP62=62
      ITAP63=63
      ITAP64=64
      ITAP65=65
      ITAP66=66
      ITAP68=86
      ITAP69=69
      ITAP71=71
      ITAP82=82
      ITAP83=83
      TPAMAT=90
      TPBMAT=91
      TPFMAT=92
      TPHMAT=93
      TPRMAT=94
      TPQMAT=95
      TPJMAT=96
      ITAP98=98
      ITAP97=97
      ITAP99=99
      fl82='NO  '
C
      NIOBF=100
      LENBUF=NIOBF*(SEC2I(1)/INTOWP(1))
C
      CALL RFILE(ITAP60)
c     CALL RFILE(ITAP61)
      CALL RFILE(ITAP62)
      CALL RFILE(ITAP63)
      CALL RFILE(ITAP64)
      CALL RFILE(ITAP65)
      CALL RFILE(ITAP66)
      CALL RFILE(ITAP82)
      CALL RFILE(ITAP83)
      CALL RFILE(ITAP97)
      CALL RFILE(ITAP98)
      CALL RFILE(ITAP99)
C
C >>> START READING FILE30
C
      CALL RFILE(ITAP30)
      CALL WREADW(ITAP30,I30,200,101,LWORD)
      MPOINT = I30( 2)
      MCONST = I30( 3)
      MCALCS = I30( 4)
      NCALCS = I30( 5)
      NBF    = I30(18)
      NAT    = I30(19)
      NBFAO  = I30(22)
      NSHELL = I30(27)
      NIRRED = I30(28)
      NPRIM  = I30(32)
      NSYMHF = I30(41)
      MXCOEF = I30(42)
      IOPEN  = I30(43)
C
      NTR=(NBF*(NBF+1))/2
      NTRAO=(NBFAO*(NBFAO+1))/2
      NT=NBF
C
C >>> ZERO OUT THE CORE
C
      CALL ZERO(CC,MAXCOR)
C
C >>> DYNAMICAL ALLOCATION
C
      EIGVAL = 1
      FOCK   = EIGVAL + NBF
      ZAN    = FOCK + NTR
      COORD  = ZAN + NBF
      NLAMDA = WPADTI(COORD + NAT*3)
      SCFTCC = NLAMDA + NSYMHF
      ITYP   = SCFTCC + NBF
      FLOV   = ITYP + NSYMHF
      ORBSYM = FLOV + 4*NIRRED
      FZO    = ORBSYM + NBF
      FZV    = FZO + NBF
      NT3    = FZV + NBF
      ITEMP  = NT3 + NIRRED
      NCL    = ITEMP + MCALCS
      NOP    = NCL + NSYMHF
      WTEMP  = IADTWP(NOP + NSYMHF)
      TOP    = WTEMP + NBF
C
C >>> READ FILE30
C
      CALL READ30(ITAP30,MPOINT,MCONST,MCALCS,NCALCS,NBF,
     .            NSYMHF,MXCOEF,CC(EIGVAL),IC(NLAMDA),6,IC(NCL),NO,
     .            CC(WTEMP),IC(SCFTCC),IC(ITYP),IC(ORBSYM),
     .            IC(FLOV),NIRRED,ENUC,ESCF,NAT,NSHELL,NPRIM,
     .            CHAR,NORD,CC(ZAN),CC(COORD),
     .            IOPEN,NOP)
C
      NV=NBF-NO
      NTRO=NO*(NO+1)/2
      NTRV=NV*(NV+1)/2
      NO2=NO*NO
      NV2=NV*NV
      NOV=MAX0(NO,NV)
      N2OV=MAX0(NTRV,NO*NV)
      NM=NO*NV
C
C >>> CALCULATE FOCK MATRIX
C
      CALL ZERO(CC(FOCK),NTR)
c     DO 42 I=1,NBF
c     IJ=(I*(I+1))/2
c  42 CC(FOCK+IJ-1)=CC(EIGVAL+I-1)
c     EIGVAL=FOCK
c     
      itap79=79
      call rfile(itap79)
      i79=1+intowp(ntr)
      call wreadw(itap79,cc(fock),intowp(ntr),i79,i79)
      call rclose(itap79,3)
      write(6,*)
      write(6,*)'fock matrix read from file79'
c     call print(cc(fock),ntr,nbf,6)
      eigval=fock
c
C >>>>>>>>>>>>>>>>>>>> Begin input processing.
C

c Make sure that the user isn't trying to use the old style input
      CALL LOCATE (INPUT,'# CCSD ###',IERR)
      if (ierr.eq.0) then
        write(6,*)'ERROR: old style ''# CCSD ###'' input has been found'
        call qabort
        endif
      CALL LOCATE (INPUT,'# ZCCSD ##',IERR)
      if (ierr.eq.0) then
        write(6,*)'ERROR: old style ''# ZCCSD ##'' input has been found'
        call qabort
        endif

      errcod = frdc('WFN',corrty)
      if (errcod.ne.0) then
        write(6,*)
        write(6,*) 'ERROR: required input ''wfn'' not found'
        call qabort
        endif

      errcod = fcnt('FROZEN_DOCC',cout)
      if ((errcod.eq.EOK) .and. (cout.ne.NIRRED)) then
         write(6,*) 'ERROR: an incorrect number of FROZEN_DOCC irreps'
     &             // ' has been given'
         call pabort
         call qabort
         endif
       errcod = frdia1('FROZEN_DOCC',NFZO,NIRRED)
       if (errcod.eq.EEOD) then
         DO 622 I=1,NIRRED
            NFZO(I)=0
  622    CONTINUE
       else if (errcod.ne.EOK) then
         write(6,*) 'ERROR: problem reading ''frozen_docc'''
         call pabort
         call qabort
         endif

       errcod = fcnt('FROZEN_UOCC',cout)
       if ((errcod.eq.EOK) .and. (cout.ne.NIRRED)) then
         write(6,*) 'ERROR: an incorrect number of FROZEN_UOCC irreps'
     &             // ' has been given'
         call pabort
         call qabort
         endif
       errcod = frdia1('FROZEN_UOCC',NFZV,NIRRED)
       if (errcod.eq.EEOD) then
         DO 1622 I=1,NIRRED
            NFZV(I)=0
 1622    CONTINUE
       else if (errcod.ne.EOK) then
         write(6,*) 'ERROR: problem reading ''frozen_uocc'''
         call pabort
         call qabort
         endif

c Set up the defaults for the rest of the input parameters
      label = ' '
      dertyp = ' '
      EPSI  =5
      NGO   =2
      MINDIM=2
      MAXDIM=8
      NCSF  =20
      ICONV =10
      ZMAXIT =30
      bDIIS = 1

c Ignore errors in the rest of the input.
      errcod = frdc('LABEL',LABEL)
      errcod = frdc('DERTYPE',dertyp)
      errcod = frdi('NCSF',NCSF)
      errcod = frdi('MAXITER',ZMAXIT)
      errcod = frdi('CONVERGENCE',ICONV)
      errcod = frdboo('RESTART',bRSTR)

c The more obscure diis parameters are prefixed with ccdiis
      errcod = frdi('CCDIIS:EPSILON',EPSI)
      errcod = frdi('CCDIIS:NGO',NGO)
      errcod = frdi('CCDIIS:MINDIM',MINDIM)
      errcod = frdi('CCDIIS:MAXDIM',MAXDIM)
      errcod = frdboo('CCDIIS:DIIS',bDIIS)


      write(6,*)
      write(6,*) 'Parameters to be used:'
      write(6,*)    '  label        = ',LABEL
      write(6,*)    '  wfn          = ',corrty
      write(6,*)    '  dertype      = ',dertyp
      write(6,*)    '  ncsf         = ',NCSF
      write(6,*)    '  maxiter      = ',ZMAXIT
      write(6,*)    '  convergence  = ',ICONV
      write(6,*)    '  restart      = ',bRSTR
      write(6,6101) '  frozen_docc  = ',(nfzo(i),i=1,nirred)
      write(6,6101) '  frozen_uocc  = ',(nfzv(i),i=1,nirred)
      write(6,*)    '  ccdiis:epsilon = ',EPSI
      write(6,*)    '  ccdiis:ngo     = ',NGO
      write(6,*)    '  ccdiis:mindim  = ',MINDIM
      write(6,*)    '  ccdiis:maxdim  = ',MAXDIM
      write(6,*)    '  ccdiis:diis    = ',bDIIS
      write(6,*)
6101  format(1X,A,8(1X,I3,1X))

c >>>>>>>>> Check the input and convert to the internal flags required by cczv.
      if (corrty.eq.'CCSD') then
        option = corrty
      else if (corrty.eq.'CCD') then
        option = corrty
      else if (corrty.eq.'CCSDT1') then
        option = 'SDT1'
      else if (corrty.eq.'CCSD_T') then
        option = 'CCT'
      else
        write(6,*) 'ERROR: An invalid ''wfn'' has been specified.'
        call qabort
        endif

      if (dertyp.ne.'FIRST') then
        write(6,*) 'ERROR: An invalid ''dertype'' has been specified'
     &          // ' or ''dertype'' is missing.'
        call qabort
        endif

      if (bDIIS.eq.0) then
        diisfl = 'NO'
      else
        diisfl = 'YES'
        endif

      if (bRSTR.eq.0) then
        ZOPT = 'NORM'
      else
        ZOPT = 'RSTR'
        endif
c
c >>>>>>>>>>>>>>>>>>>>>> Done with input processing
C
C >>> FORCE SORT TO 'NO'
C
      SORT='NO  '
C
      CALL FROZEN (IC(FLOV),NOSYM,NVSYM,NFZO,NFZV,NDOC,NUOC,
     .             IC(FZO),IC(FZV),IC(ITYP),NIRRED,NO,NV,JOUT,
     .             NOACT,NVACT,NSYMHF)
C
      CALL NCOUNT(IC(ORBSYM),IC(FLOV),NIRRED,NO,NV,NBF,
     .            IC(FZO),IC(FZV),NT1,NT2,
     .            NSGOO,NSGVV,NSGOV,NSHOV,NSLOV,NSLVO,
     .            IC(NT3),NTAU,OPTION,DIMAR,maxo,maxv)
C
      LENA=NTRO*(NTRO+1)/2
      LENB=NTRV*(NTRV+1)/2
      LENC=NTRO*NTRV
      LEND=NO*NV*(NO*NV+1)/2
      LENE=NTRO*NO*NV
      LENF=NTRV*NO*NV
      NINT=LENA+LENB+LENC+LEND+LENE+LENF
      NSINT=NSGOO+NSGVV+NSGOV+NSHOV+NSLOV+NSLVO
C
C ~~~~~~~~~~~~  ALLOCATE CC MEMORY  ~~~~~~~~~~~~~~~~
C
C   NOTE: TOFF AND TADD ARE THE "UOFF AND VADD" POINTERS OF BCCSD.
C   ZLX IS ALSO CALCULATED IN SUBROUTINE OFFSET AND NO LONGER IN SYMARR.
C   THE BCCSD ITR(NOV) AND ITV(NV) ARE INO AND INV CALCULATED IN OFFSET.
C   OJO NO2, NONO, NV2 AND NVNV!
C
      IOFF  = WPADTI(TOP)
      INO   = IOFF + NTR
      INV   = INO + NOV
      INTRO = INV + NOV
      INTRV = INTRO + N2OV
      ZLX   = INTRV + NO*NV
      UOFF  = ZLX + NV2
      ZLXO  = UOFF + 2*NO2
      TOFF  = ZLXO + NO2
      TADD  = TOFF + NO*NO*2
      AOFF  = TADD + NV*NV
      AADD  = AOFF + NO*NO
      BOFF  = AADD + NO*NO
      BADD  = BOFF + NV*NV
      COFF  = BADD + NV*NV
      CADD  = COFF + NV*NV
      DOFF  = CADD + NO*NO
      DADD  = DOFF + NM
      EOFF  = DADD + NM
      EADD  = EOFF + NM
      EMOFF = EADD + NO*NO
      EMADD = EMOFF + NM
      FOFF  = EMADD + NO*NO
      FADD  = FOFF + NV*NV
      UOFFV = FADD + NM
      IVAL  = UOFFV + 2*NV2
      ITOP  = IVAL + NCSF
      T1    = IADTWP( ITOP )
      TOP   = T1 + NO*NV
C
      IF(OPTION.EQ.'SDT1')THEN
         TAU   = TOP
         T3OFF = WPADTI(TAU + NTAU)
         WRITE(6,*)'NTAU',NTAU
         T3ADD = T3OFF + NO*NO*NIRRED*2
         TOP   = IADTWP(T3ADD + NV*NV*NIRRED)
      ELSE
         TAU   = 1
         T3OFF = 1
         T3ADD = 1
      ENDIF
C
      CALL OFFSET(IC(INO),IC(INV),IC(INTRO),IC(INTRV),
     .            IC(ZLX),IC(UOFF),NBF,NO,NV,NOV,NTRO,NTRV,N2OV,
     .            JOUT,IC(ZLXO),IC(UOFFV),IC(IOFF),NTR)
C
      CALL SYMARR(IC(ORBSYM),IC(FLOV),NIRRED,NO,NV,NBF,NOV,NM,
     .            IC(FZO),IC(FZV),NFZO,NFZV,
     .            IC(TOFF),IC(TADD),IC(INO),IC(INV),
     .            NO*NO,NV*NV,
     .            IC(AOFF),IC(AADD),IC(BOFF),IC(BADD),
     .            IC(COFF),IC(CADD),IC(DOFF),IC(DADD),
     .            IC(EOFF),IC(EADD),IC(FOFF),IC(FADD),
     .            IC(T3OFF),IC(T3ADD),OPTION,
     .            IC(EMOFF),IC(EMADD),LEMAT)
C
C >>> ALLOCATE CORE FOR ARRAYS NEEDED IN REFORMULATION
C
      AMAT  = TOP
      GOMAT = AMAT + NO*NO
      GVMAT = GOMAT + NO2
      I1MAT = GVMAT + NV2
      TZV   = I1MAT + NV*NO
      TZO   = TZV + NV2
      EMAT  = TZO + NO2
      CCC   = EMAT + LEMAT
      BB    = CCC + MAXDIM
      BB2   = BB + (MAXDIM+1)*(MAXDIM+2)
      LAST  = BB2 + (MAXDIM+1)*(MAXDIM+2)-1
      AR1   = LAST
      AR3   = AR1
      AR5   = AR1
      AR2   = AR1 + DIMAR
      AR4   = AR2
      AR6   = AR2
      AUX1  = AR2 + DIMAR
      AUX2  = AUX1 + MAX0(NV*NV,NO*NO)
      AUX3  = AUX2 + MAX0(NV*NV,NO*NO)
      SFV   = AUX3 + max0(NV*NV,no*no)
      ABIN  = SFV  + max0(NV*NV,no*no)
      XMAT  = ABIN + MAX0(NV*NV,NO*NO)
      LXMAT = MAX0(NO*NO,NV*NV)
      BMAT  = XMAT
      FMAT  = XMAT
      HMAT  = XMAT
      JMAT  = XMAT
      RMAT  = XMAT
      QMAT  = XMAT + LXMAT
      OMG1  = QMAT + NO*NV
      OMG2  = OMG1 + NO*NOV
      OMG3  = OMG2 + NO*NV
      TVAL  = OMG3 + NV2
      TOP   = TVAL + NCSF
      IF (TOP.GT.MAXCOR) THEN
         WRITE(*,*)'TOP',TOP
         WRITE (JOUT,*) 'ERROR: BEFORE READ OF INTS : NOT ENOUGH MEMORY'
         call qabort
      ENDIF
C
      WRITE(*,*)
      WRITE(*,*)'LENGTH OF INTEGRALS ARRAYS AND T AMPLITUDES'
      WRITE(*,*)
      WRITE(*,*)'NGOO =',LENA,'NSGOO =',NSGOO
      WRITE(*,*)'NGVV =',LENB,'NSGVV =',NSGVV
      WRITE(*,*)'NGOV =',LENC,'NSGOV =',NSGOV
      WRITE(*,*)'NHOV =',LEND,'NSHOV =',NSHOV
      WRITE(*,*)'NLOV =',LENE,'NSLOV =',NSLOV
      WRITE(*,*)'NLVO =',LENF,'NSLVO =',NSLVO
      WRITE(*,*)'NINT =',NINT,'NSINT =',NSINT
      WRITE(*,*)'NT1  =',NT1 ,'NT2   =',NT2
      WRITE(*,*)'EMAT =',LEMAT
      WRITE(6,*)
C
      MAXINT=NSGOO
C     IF(NSGVV.GT.MAXINT)MAXINT=NSGVV
      IF(NSGOV.GT.MAXINT)MAXINT=NSGOV
      IF(NSHOV.GT.MAXINT)MAXINT=NSHOV
      IF(NSLOV.GT.MAXINT)MAXINT=NSLOV
      IF(NSLVO.GT.MAXINT)MAXINT=NSLVO
      WRITE(6,*)' '
      WRITE(6,*)' INTEGRAL ARRAY IS ',MAXINT
      WRITE(6,*)' '
C
C >>> ALLOCATE CORE FOR INTERMEDIATE PART
C
      DINTS = TOP
      MOINTS= TOP
      BUF   = MOINTS + MAXINT
      IBUF  = WPADTI(BUF)
      T2    = BUF + LENBUF
      TOPINT= T2 + NT2
C
C >>> ALLOCATE CORE FOR ITERATIVE PART
C
      Z2O   = DINTS + NSHOV
      Z2    = Z2O + NT2
      Z1O   = Z2 + NT2
      Z1    = Z1O + NO*NV
      TOPITR= Z1 + NO*NV
C
      IF(OPTION.EQ.'SDT1')THEN
         DINTS = TOPINT
         ein   = DINTS + NSHOV
         fin   = ein   + nslov
         z2o   = fin   + nslvo
         Z2    = Z2O + NT2
         Z1O   = Z2 + NT2
         Z1    = Z1O + NO*NV
         TOPITR= Z1 + NO*NV
         
      ENDIF
C
      ICRINT = TOPINT*8/1024
      ICRITR = TOPITR*8/1024
      IMCOR  = MAXCOR*8/1024
      WRITE(JOUT,6001)IMCOR,MAXCOR,ICRINT,TOPINT,ICRITR,TOPITR
 6001 FORMAT(//,' MAXIMUM  CORE  = ',I7,2X,'K BYTES  OR  ',I9 ,' REAL',
     .         ' WORDS',/,
     .          ' REQUIRED CORE1 = ',I7,2X,'K BYTES  OR  ',I9 ,' REAL',
     .         ' WORDS',/,
     .          ' REQUIRED CORE2 = ',I7,2X,'K BYTES  OR  ',I9 ,' REAL',
     .         ' WORDS',//)
C
      TOP=MAX0(TOPINT,TOPITR)
      IF (TOP.GT.MAXCOR) THEN
         WRITE(6,*)'   '
         WRITE(6,*)'ERROR: NOT ENOUGH CORE. PROGRAM QUITS ...'
         call qabort
      ENDIF
C
C >>> READ TAPE69
C
      CALL RFILE(ITAP69)
      CALL SREAD(ITAP69,NDIMT1,1)
      CALL SREAD(ITAP69,NDIMT2,1)
      IF(NDIMT2.NE.NT2)THEN
         WRITE(6,*)' PROBLEMS ... NT2.NE.NDIMT2'
      ENDIF
      CALL READ69(CC(T1),CC(T2),NO,NV,NDIMT2,ITAP69,FILE7,JOUT,
     .            IT269,IZ169)
C
C>>>  FORM INTEGRAL DRIVEN INTERMEDIATES AND STORE THEM IN FILES 90'S
C
C
C     A INTS
C
      CALL RDINTS(ITAP60,CC(MOINTS),NSGOO,
     .            CC(BUF),IC(IBUF),LENBUF,NO,NV,
     .            IC(INO),IC(INO),IC(AOFF),IC(AADD))
      CALL AINTS(CC(AMAT),CC(MOINTS),NTRO,NO2,NSGOO,NO,NTR,
     .           IC(IOFF),JOUT,IC(INO),NOV,
     .           IC(AOFF),NO*NO,IC(AADD),NO*NO,
     .           IC(FLOV),NIRRED,IC(ORBSYM),IC(FZO),IC(FZV))
C
C     B INTS
C
      CALL BINTS (CC(BMAT),CC(AUX3),NV,NO,IC(FLOV),IC(ORBSYM),NIRRED,
     .            IC(FZO),IC(FZV))
C
C     C INTS
C
      CALL RDINTS(ITAP62,CC(MOINTS),NSGOV,
     .            CC(BUF),IC(IBUF),LENBUF,NO,NV,
     .            IC(INV),IC(INO),IC(COFF),IC(CADD))
      CALL CINTS(CC(EMAT),CC(FMAT),CC(MOINTS),NO,NV,
     .           CC(T1),NSGOV,NTRO,NTR,IC(IOFF),JOUT,IC(INO),
     .           IC(INTRO),IC(INTRV),IC(INV),NOV,N2OV,
     .           IC(FLOV),IC(ORBSYM),NIRRED,
     .           IC(COFF),NV*NV,IC(CADD),NO*NO,
     .           IC(EMOFF),IC(EMADD),LEMAT,IC(FZO),IC(FZV))
      CALL RCINTS(CC(RMAT),CC(MOINTS),NO,NV,
     .            NSGOV,NTRO,NTR,IC(IOFF),JOUT,IC(INO),
     .            IC(INTRO),IC(INTRV),IC(INV),NOV,N2OV,
     .            IC(FLOV),IC(ORBSYM),NIRRED,
     .            IC(COFF),NV*NV,IC(CADD),NO*NO,IC(FZO),IC(FZV))
C
C     E INTS
C
      CALL RDINTS(ITAP64,CC(MOINTS),NSLOV,
     .            CC(BUF),IC(IBUF),LENBUF,NO,NV,
     .            IC(INO),IC(INO),IC(EOFF),IC(EADD))
      CALL EINTS(CC(JMAT),CC(GOMAT),CC(EMAT),
     .           CC(FMAT),CC(AMAT),CC(MOINTS),NO,NV,
     .           CC(T1),CC(T2),NSLOV,NO2,NV2,NTRO,NTRV,NTR,
     .           IC(IOFF),JOUT,CC(EIGVAL),IC(INO),IC(INTRO),
     .           IC(INTRV),NOV,IC(INV),N2OV,IC(ZLX),
     .           IC(UOFF),IC(FLOV),IC(ORBSYM),NIRRED,
     .           IC(EOFF),NM   ,IC(EADD),NO*NO,
     .           IC(TOFF),IC(TADD),NDIMT2,
     .           IC(EMOFF),IC(EMADD),LEMAT,IC(FZO),IC(FZV))
      CALL RQEINT(CC(RMAT),CC(QMAT),CC(MOINTS),CC(T1),
     .            NSLOV,NO,NV,NO2,NV2,NTRO,NTRV,NTR,
     .            IC(IOFF),JOUT,IC(INO),IC(INTRO),
     .            IC(INTRV),NOV,IC(INV),N2OV,IC(FLOV),
     .            IC(ORBSYM),NIRRED,
     .            IC(EOFF),NM,IC(EADD),NO*NO,IC(FZO),IC(FZV))
C
C     F INTS
C
      CALL RDINTS(ITAP65,CC(MOINTS),NSLVO,
     .            CC(BUF),IC(IBUF),LENBUF,NO,NV,
     .            IC(INV),IC(INO),IC(FOFF),IC(FADD))
      CALL FINTS(CC(GVMAT),CC(EMAT),
     .           CC(BMAT),CC(MOINTS),NO,NV,
     .           CC(T1),CC(T2),NSLVO,NO2,NV2,NTRO,NTRV,NTR,
     .           IC(IOFF),JOUT,CC(FOCK),IC(INO),IC(INTRO),
     .           IC(INTRV),NOV,IC(INV),N2OV,IC(ZLX),
     .           IC(UOFF),IC(FLOV),IC(ORBSYM),NIRRED,
     .           IC(FOFF),NV*NV,IC(FADD),NM,
     .           IC(TOFF),IC(TADD),NDIMT2,
     .           CC(AR3),CC(AR4),CC(AUX1),CC(AUX2),
     .           CC(AUX3),CC(SFV),CC(OMG1),
     .           IC(EMOFF),IC(EMADD),LEMAT,IC(FZO),IC(FZV))
C
C  NFINTS REQUIRES BOTH MOINTS AND DINTS. THIS IS THE ONLY PLACE.
C
      CALL NFINTS(CC(FMAT),CC(MOINTS),CC(T1),CC(T2),NO2,NV2,
     .            NTRO,NTRV,NO,NV,NSLVO,NTR,IC(IOFF),JOUT,IC(INO),
     .            IC(INTRO),IC(INTRV),NOV,IC(INV),N2OV,
     .            IC(ZLX),IC(UOFF),IC(FLOV),IC(ORBSYM),
     .            NIRRED,
     .            IC(FOFF),NV*NV,IC(FADD),NM,
     .            IC(TOFF),IC(TADD),NDIMT2,
     .            CC(DINTS),IC(DOFF),IC(DADD),
     .            CC(AR1),CC(AR2),CC(OMG1),CC(OMG2),
     .            CC(AUX1),CC(AUX2),CC(AUX1),
     .            CC(BUF),IC(IBUF),LENBUF,NSHOV,NSLVO,IC(FZO),IC(FZV))
      CALL HFINTS(CC(HMAT),CC(MOINTS),NO2,NV2,
     .            NTRO,NTRV,NO,NV,NSLVO,NTR,IC(IOFF),JOUT,IC(INO),
     .            IC(INTRO),IC(INTRV),NOV,IC(INV),N2OV,
     .            IC(FLOV),IC(ORBSYM),NIRRED,
     .            IC(FOFF),NV*NV,IC(FADD),NM,IC(FZO),IC(FZV))
      CALL RQFINT(CC(QMAT),CC(RMAT),CC(MOINTS),CC(T1),
     .            NSLVO,NO,NV,NO2,NV2,
     .            NTRO,NTRV,NTR,IC(IOFF),JOUT,IC(INO),
     .            IC(INTRO),IC(INTRV),NOV,IC(INV),N2OV,
     .            IC(FLOV),IC(ORBSYM),NIRRED,
     .            IC(FOFF),NV*NV,IC(FADD),NM,IC(FZO),IC(FZV))
C
C     D INTS
C
      CALL RDINTS(ITAP63,CC(DINTS),NSHOV,
     .            CC(BUF),IC(IBUF),LENBUF,NO,NV,
     .            IC(INO),IC(INO),IC(DOFF),IC(DADD))
      CALL DINTSS(CC(I1MAT),CC(GVMAT),CC(GOMAT),
     .            CC(JMAT),CC(AMAT),CC(OMG3),
     .            CC(FMAT),CC(EMAT),CC(T1),CC(T2),CC(OMG1),
     .            CC(OMG2),CC( DINTS),NSHOV,NO,NV,NTRO,NTRV,NO2,
     .            NV2,NTR,IC(IOFF),JOUT,IC(INO),IC(INTRO),
     .            IC(INTRV),NOV,IC(INV),N2OV,IC(ZLX),
     .            IC(UOFF),IC(FLOV),IC(ORBSYM),NIRRED,
     .            IC(DOFF),NM,IC(DADD),NM,
     .            IC(TOFF),IC(TADD),NDIMT2,
     .            CC(AR3),CC(AR4),CC(AUX1),
     .            CC(AUX2),CC(SFV),CC(AUX3),
     .            IC(EMOFF),IC(EMADD),LEMAT,IC(FZO),IC(FZV))
      CALL HDINTS(CC(HMAT),CC( DINTS),NO2,NV2,CC(T1),
     .            NTRO,NTRV,NO,NV,NSHOV,NTR,IC(IOFF),JOUT,IC(INO),
     .            IC(INTRO),IC(INTRV),NOV,IC(INV),N2OV,
     .            IC(FLOV),IC(ORBSYM),NIRRED,
     .            IC(DOFF),NM,IC(DADD),NM,IC(FZO),IC(FZV))
      CALL BDINTS(CC(BMAT),CC( DINTS),NO2,NV2,CC(T1),CC(T2),
     .            NTRO,NTRV,NO,NV,NSHOV,NTR,IC(IOFF),JOUT,
     .            IC(INO),IC(INTRO),IC(INTRV),NOV,IC(INV),
     .            N2OV,IC(FLOV),IC(ORBSYM),NIRRED,IC(ZLX),
     .            IC(UOFF),
     .            IC(DOFF),NM   ,IC(DADD),NM,
     .            IC(TOFF),IC(TADD),NDIMT2,
     .            CC(AR5),CC(AR6),CC(AUX1),CC(AUX2),
     .            CC(AUX3),CC(SFV),CC(OMG3),IC(FZO),IC(FZV))
      CALL RQDINT(CC(QMAT),CC(RMAT),CC( DINTS),
     .            NO2,NV2,CC(T1),CC(T2),
     .            NTRO,NTRV,NO,NV,NSHOV,NTR,IC(IOFF),JOUT,
     .            IC(INO),IC(INTRO),IC(INTRV),NOV,IC(INV),
     .            N2OV,IC(ZLX),IC(UOFF),IC(FLOV),IC(
     .            ORBSYM),NIRRED,
     .            IC(DOFF),NM   ,IC(DADD),NM,
     .            IC(TOFF),IC(TADD),NDIMT2,
     .            CC(AR1),CC(AR2),CC(AUX1),CC(AUX2),
     .            CC(OMG1),CC(OMG2),IC(FZO),IC(FZV))
C
C >>> CHECK THE CORRELATION ENERGY
C
      CALL ECORR(CC(T2),CC(T1),CC( DINTS),NSHOV,NO,NV,NO2,NV2,
     .           NTRO,NTRV,IC(IOFF),NTR,JOUT,IC(INO),IC(INTRO),
     .           IC(INTRV),NOV,IC(INV),N2OV,IC(ZLX),
     .           IC(UOFF),IC(FLOV),IC(ORBSYM),NIRRED,
     .           IC(DOFF),NM,IC(DADD),NM,
     .           IC(TOFF),IC(TADD),NDIMT2,IC(FZO),IC(FZV))
C
      CALL ADDE(CC(EMAT),CC(I1MAT),CC(JMAT),CC(T1),
     .          CC(T2),IC(IOFF),NO,NV,NTRO,NTRV,NO2,NV2,NTR,JOUT,
     .          IC(INO),IC(INTRO),IC(INTRV),NOV,IC(INV),
     .          N2OV,IC(ZLX),IC(ZLXO),IC(UOFF),IC(FLOV),
     .          IC(ORBSYM),NIRRED,
     .          IC(TOFF),IC(TADD),NDIMT2,
     .          IC(EMOFF),IC(EMADD),LEMAT,IC(FZO),IC(FZV))
C
      WRITE(6,*)
      WRITE(6,*)'ALL INTERMEDIATES HAVE BEEN FORMED'
      WRITE(4,*)'ALL INTERMEDIATES HAVE BEEN FORMED'
C      CALL TIMIT(4,JOUT)
      close(unit=4,err=51)
   51 continue
      call ffile(4,'term',1)
c
C
      IF (ZOPT.EQ.'RSTR') THEN
         CALL RSETSA(ITAP69,IZ169)
         CALL SREAD(ITAP69,CC(Z1O),INTOWP(NO*NV))
         CALL SREAD(ITAP69,CC(Z2O),INTOWP(NDIMT2))
      ELSE
         IFLG=0
         CALL INITZ(CC(Z1O),CC(Z2O),CC(I1MAT),CC( DINTS),NSHOV,
     .              CC(EIGVAL),IC(IOFF),NO,NV,NTRO,NTRV,
     .              NO2,NV2,NTR,NBF,JOUT,IFLG,IC(INO),IC(INTRO),
     .              IC(INTRV),NOV,IC(INV),N2OV,IC(FLOV),IC(
     .              ORBSYM),NIRRED,
     .              IC(DOFF),NM,IC(DADD),NM,
     .              IC(ZLX),IC(TOFF),IC(TADD),NDIMT2,IC(FZO),IC(FZV),
     .              option)
      ENDIF
      IFLG=1
c
      if(option.eq.'CCD ')call zero(cc(z1),no*nv)
C
      NIT=0
      DELZ=2.0D+00
      DELZ2=2.0D+00
      ZCONV=10.D0 ** (-ICONV)
C
C >>> START ITERATIVE PROCEDURE
C
 1001 CONTINUE
C
C     THIS IS THE ONLY PLACE IN THE ITERATIVE PART WHERE T2 IS NEEDED.
C     IT IS READ INTO Z2 AND PASS TO T2Z2.
C
      CALL RSETSA(ITAP69,IT269)
      CALL SREAD(ITAP69,CC(Z2),INTOWP(NDIMT2))
      CALL T2Z2(CC(TZV),CC(TZO),CC(Z2O),CC(Z2),NO,NV,NO2,NV2,
     .          NTRO,NTRV,NTR,IC(IOFF),JOUT,IC(INO),
     .          IC(INTRO),IC(INTRV),NOV,IC(INV),N2OV,IC(ZLX),
     .          IC(UOFF),IC(FLOV),IC(ORBSYM),NIRRED,
     .          IC(TOFF),IC(TADD),NDIMT2,IC(FZO),IC(FZV))
C
      CALL INITZ(CC(Z1),CC(Z2),CC(I1MAT),CC( DINTS),NSHOV,
     .           CC(EIGVAL),IC(IOFF),NO,NV,NTRO,NTRV,
     .           NO2,NV2,NTR,NBF,JOUT,IFLG,IC(INO),IC(INTRO),
     .           IC(INTRV),NOV,IC(INV),N2OV,IC(FLOV),
     .           IC(ORBSYM),NIRRED,
     .           IC(DOFF),NM   ,IC(DADD),NM,
     .           IC(ZLX),IC(TOFF),IC(TADD),NDIMT2,IC(FZO),IC(FZV),
     .           option)
      CALL RQTERM(CC(Z1),CC(Z2),CC(Z1O),CC(Z2O),
     .            CC(RMAT),CC(QMAT),
     .            IC(IOFF),NO,NV,NO2,NV2,NTRO,NTRV,NTR,JOUT,
     .            IC(INO),IC(INTRO),IC(INTRV),NOV,IC(INV),
     .            N2OV,IC(ZLX),IC(UOFF),
     .            IC(FLOV),IC(ORBSYM),NIRRED,
     .            CC(AUX1),CC(AUX2),CC(AR1),CC(AR2),
     .            IC(TOFF),IC(TADD),NDIMT2,IC(FZO),IC(FZV))
      CALL HJTERM(CC(Z1),CC(Z2),CC(Z1O),CC(TZO),CC(TZV),
     .            CC(HMAT),CC(JMAT),NO,NV,NO2,NV2,NTRO,NTRV,NTR,
     .            IC(IOFF),JOUT,CC(FOCK),IC(INO),
     .            IC(INTRO),IC(INTRV),NOV,IC(INV),N2OV,
     .            IC(FLOV),IC(ORBSYM),NIRRED,
     .            IC(ZLX),IC(TOFF),IC(TADD),NDIMT2,IC(FZO),IC(FZV))
      CALL ABTERM(CC(Z2),CC(Z2O),CC(AMAT),CC(BMAT),NO,NV,NO2,
     .            NV2,NTRO,NTRV,NTR,IC(IOFF),JOUT,IC(INO),
     .            IC(INTRO),IC(INTRV),NOV,IC(INV),N2OV,
     .            IC(ZLX),IC(UOFF),IC(ZLXO),IC(UOFFV),
     .            IC(FLOV),IC(ORBSYM),NIRRED,
     .            CC(Z1),CC(T1),CC(AR3),CC(AR4),
     .            CC(AR5),CC(AR6),CC(ABIN),CC(SFV),
     .            CC(AUX1),CC(AUX2),CC(AUX3),
     .            IC(TOFF),IC(TADD),NDIMT2,IC(FZO),IC(FZV))
      CALL ETERM(CC(Z1),CC(Z2O),CC(EMAT),NO,NV,NO2,NV2,NTRO,
     .           NTRV,NTR,IC(IOFF),JOUT,IC(INO),IC(INTRO),
     .           IC(INTRV),NOV,IC(INV),N2OV,IC(ZLX),
     .           IC(UOFF),IC(FLOV),IC(ORBSYM),NIRRED,
     .           IC(TOFF),IC(TADD),NDIMT2,
     .           IC(EMOFF),IC(EMADD),LEMAT,IC(FZO),IC(FZV))
      CALL FTERM(CC(Z1),CC(Z2O),CC(FMAT),NO,NV,NO2,NV2,NTRO,
     .           NTRV,NTR,IC(IOFF),JOUT,IC(INO),IC(INTRO),
     .           IC(INTRV),NOV,IC(INV),N2OV,IC(ZLX),
     .           IC(UOFF),IC(FLOV),IC(ORBSYM),NIRRED,
     .           IC(TOFF),IC(TADD),NDIMT2,IC(FZO),IC(FZV))
      CALL MTERM(CC(Z2),CC(TZO),CC(TZV),NO,NV,NO2,NV2,NTRO,
     .           NTRV,NTR,IC(IOFF),CC( DINTS),NSHOV,JOUT,
     .           IC(INO),IC(INTRO),IC(INTRV),NOV,
     .           IC(INV),N2OV,
     .           IC(FLOV),IC(ORBSYM),NIRRED,
     .           IC(DOFF),NM,IC(DADD),NM,
     .           IC(TOFF),IC(TADD),NDIMT2,IC(FZO),IC(FZV))
      CALL GTERM(CC(Z2),CC(Z1),CC(Z1O),CC(Z2O),CC(GOMAT),
     .           CC(GVMAT),NO,NV,NO2,NV2,NTRO,NTRV,NTR,IC(IOFF),
     .           JOUT,CC(FOCK),CC(T1),IC(INO),IC(INTRO),
     .           IC(INTRV),NOV,IC(INV),N2OV,IC(ZLX),
     .           IC(UOFF),IC(FLOV),IC(ORBSYM),NIRRED,
     .           IC(TOFF),IC(TADD),NDIMT2,IC(FZO),IC(FZV))
      CALL ITERM(CC(Z2),CC(Z1O),CC(I1MAT),NO,NV,NO2,NV2,NTRO,
     .           NTRV,NTR,IC(IOFF),JOUT,IC(INO),IC(INTRO),
     .           IC(INTRV),NOV,IC(INV),N2OV,
     .           IC(ZLX),IC(FLOV),IC(ORBSYM),NIRRED,
     .           IC(TOFF),IC(TADD),NDIMT2,IC(FZO),IC(FZV))
C
C      R3 MATRIX FOR CCSDT
C
C  T2 IS NEEDED IN Y3 ONLY IF SORT='YES '. SORT IS FORCED TO 'NO' ABOVE
C
      IF(OPTION.EQ.'SDT1')THEN
c        CALL Y3(CC(MOINTS),NSHOV,NSLVO,NSLOV,MAXINT,NO,NV,
c    .           CC(BUF),IC(IBUF),LENBUF,
c    .           IC(DOFF),NM   ,IC(DADD),NM,
c    .           IC(EOFF),NM,IC(EADD),NO*NO,
c    .           IC(FOFF),NV*NV,IC(FADD),NM,
c    .           IC(FLOV),IC(ORBSYM),IC(ZLX),NIRRED,
c    .           CC(TAU),IC(TOFF),IC(TADD),NDIMT2,NTAU,
c    .           IC(INO),NOV,IC(INV),NV,
c    .           CC(Z1O),CC(Z2O),CC(Z1),CC(Z2),NO2,NTRV,ITAP97,
c    .           IC(T3OFF),IC(T3ADD),IC(NT3),
c    .           IC(IOFF),NTR,IC(UOFF),CC(FOCK),
c    .           ITAP82,ITAP83,CC(T2),NNT2,SORT,NIT,ic(fzo),ic(fzv))
c
cges this is new reformulted/optimized version of y3. 
c    file82 is created in the last iteration.
c
      call YSDT1(cc(z2),cc(z1),cc(z2o),cc(z1o),cc(TAU),
     .           NO,NV,NDIMT2,NIT,
     .           ic(ORBSYM),ic(FLOV),NIRRED,ic(T3OFF),
     .           ic(T3ADD),ic(ZLX),ic(FZO),ic(FZV),
     .           ic(tOFF),ic(tADD),ic(NT3),
     .         IC(DOFF),IC(DADD),IC(EOFF),IC(EADD),IC(FOFF),IC(FADD),
     .           ic(inv),ic(ino),ic(ioff),
     .           cc(DINts),cc(EIN),cc(FIN),fl82,cc(FOCK),
     .           cc(AUX1),cc(AUX3),cc(omg1),
     .           cc(AR1),cc(AR2),MAXO,MAXV,
     .           cc(AR1),cc(omg2),cc(AR1),cc(AR2),
     .           cc(AUX2),cc(sfv),cc(ar1),cc(ar2),cc(ar1),cc(ar2),
     .           cc(buf),ic(ibuf),lenbuf,nshov,nslov,nslvo)
      ENDIF
c
      if(option.eq.'CCD ')call zero(cc(z1),no*nv)
c
      CALL NORMZ(CC(Z2),CC(Z1),CC(EIGVAL),NO,NV,NBF,NTR,
     .           NO2,NTRV,IC(IOFF),JOUT,IC(INO),IC(INTRO),
     .           IC(INTRV),NOV,IC(INV),N2OV,
     .           IC(ZLX),IC(FLOV),IC(ORBSYM),NIRRED,
     .           IC(TOFF),IC(TADD),NDIMT2,IC(FZO),IC(FZV))
      CALL ZDIFF(CC(Z2),CC(Z1),CC(Z2O),CC(Z1O),DELZ,NO,NV,
     .           NTR,NO2,NTRV,IC(IOFF),JOUT,IC(INO),
     .           IC(INTRO),IC(INTRV),NOV,IC(INV),N2OV,DELZ2,
     .           IC(ZLX),IC(FLOV),IC(ORBSYM),NIRRED,
     .           IC(TOFF),IC(TADD),NDIMT2,IC(FZO),IC(FZV))
C
      IF(DIISFL.EQ.'YES ') THEN
      CALL DIISD(CC(Z1O),CC(Z1),CC(Z2O),CC(Z2),NO,NV,
     .           NDIMT2  ,NIT,IT,IECORR,EPSI,NGO,NDIIS,MINDIM,MAXDIM,
     .           CC(CCC),CC(BB),ITAP98,ITAP99,ITC,CC(BB2))
      ENDIF
C
      CALL NEWZ (CC(Z2),CC(Z1),CC(Z2O),CC(Z1O),NO,NV,
     .           NTR,NO2,NTRV,IC(IOFF),JOUT,IC(INO),IC(INTRO),
     .           IC(INTRV),NOV,IC(INV),N2OV,
     .           IC(ORBSYM),IC(FLOV),NIRRED,
     .           IC(ZLX),IC(TOFF),IC(TADD),NDIMT2,IC(FZO),IC(FZV))
C
      NIT=NIT+1
      WRITE (JOUT,6057) NIT,DELZ,DELZ2
      WRITE (4   ,6057) NIT,DELZ,DELZ2
 6057 FORMAT (/,2X,'NIT =',I3,' DELTA Z =',F20.12,10X,F20.12)
C
C      IF(NIT.LE.5)CALL TIMIT(4,JOUT)
C
      CALL RSETSA(ITAP69,IZ169)
      CALL SWRIT(ITAP69,CC(Z1),INTOWP(NO*NV))
      CALL SWRIT(ITAP69,CC(Z2),INTOWP(NDIMT2))
C
C >>> CHECK CONVERGENCE AND MAXIMUM NUMBER OF ITERATIONS
C
      close(unit=4,err=52)
   52 continue
      call ffile(4,'term',1)
c
      IF (NIT.Lt.ZMAXIT.AND.DELZ2.GT.ZCONV) GO TO 1001
c
c >>> create file82 for ccqg 
c
      if (option.eq.'SDT1')then
      fl82='YES '
      call YSDT1(cc(z2),cc(z1),cc(z2o),cc(z1o),cc(TAU),
     .           NO,NV,NDIMT2,NIT,
     .           ic(ORBSYM),ic(FLOV),NIRRED,ic(T3OFF),
     .           ic(T3ADD),ic(ZLX),ic(FZO),ic(FZV),
     .           ic(tOFF),ic(tADD),ic(NT3),
     .         IC(DOFF),IC(DADD),IC(EOFF),IC(EADD),IC(FOFF),IC(FADD),
     .           ic(inv),ic(ino),ic(ioff),
     .           cc(DINts),cc(EIN),cc(FIN),fl82,cc(FOCK),
     .           cc(AUX1),cc(AUX3),cc(omg1),
     .           cc(AR1),cc(AR2),MAXO,MAXV,
     .           cc(AR1),cc(omg2),cc(AR1),cc(AR2),
     .           cc(AUX2),cc(sfv),cc(ar1),cc(ar2),cc(ar1),cc(ar2),
     .           cc(buf),ic(ibuf),lenbuf,nshov,nslov,nslvo)
      ENDIF
C
C >>> EXIT NZCCSD
C
      CALL TWRIT (CC(Z1),CC(Z2),IC(ORBSYM),IC(FLOV),NIRRED,
     .            IC(TOFF),IC(TADD),IC(IVAL),CC(TVAL),
     .            IC(ZLX),IC(FZO),IC(FZV),NO,NV,NSYMHF,
     .            IC(ITYP),NT1,NSHOV,NCSF,IC(NCL))
C
C     CALL RSETSA(ITAP69,IZ169)
C     CALL SWRIT(ITAP69,CC(Z1),INTOWP(NO*NV))
C     CALL SWRIT(ITAP69,CC(Z2O),INTOWP(NDIMT2))
C
      CALL RCLOSE(ITAP60,3)
c     CALL RCLOSE(ITAP61,3)
c     CALL RCLOSE(ITAP61,4)
      CALL RCLOSE(ITAP62,3)
      CALL RCLOSE(ITAP63,3)
      CALL RCLOSE(ITAP64,3)
      CALL RCLOSE(ITAP65,3)
c     CALL RCLOSE(ITAP65,4)
      CALL RCLOSE(ITAP66,3)
c     CALL RCLOSE(ITAP66,4)
      CALL RCLOSE(ITAP69,3)
      CALL RCLOSE(ITAP83,3)
      IF(OPTION.EQ.'SDT1')THEN
         CALL RCLOSE(ITAP97,3)
      ENDIF
      CALL RFILE (TPAMAT)
      CALL RFILE (TPBMAT)
      CALL RFILE (TPFMAT)
      CALL RFILE (TPHMAT)
      CALL RFILE (TPRMAT)
      CALL RFILE (TPQMAT)
      CALL RFILE (TPJMAT)
      CALL RCLOSE(TPAMAT,4)
      CALL RCLOSE(TPBMAT,4)
      CALL RCLOSE(TPFMAT,4)
      CALL RCLOSE(TPHMAT,4)
      CALL RCLOSE(TPRMAT,4)
      CALL RCLOSE(TPQMAT,4)
      CALL RCLOSE(TPJMAT,4)
      IF(DIISFL.EQ.'YES ')THEN
         CALL RCLOSE(ITAP98,4)
         CALL RCLOSE(ITAP99,4)
      ENDIF
C
      CALL TSTOP(6)
      call psdone()
      stop
      END
