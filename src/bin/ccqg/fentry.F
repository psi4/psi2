cges  PROGRAM CCQG
      subroutine fentry(cc,ic,maxcor)
C
C    >>>>  LAST UPDATE  april 17, 1990   <<<
C
C**********************************************************************
C    NOTICE OF PROGRAM MODIFICATION                                   *
C         BY:  GUSTAVO E. SCUSERIA                                    *
c
c         by:  Edward Seidl (cets030491)
c       date:  March 4, 1991
c     reason:  add new style input, program dies if old style input
c              is found
C                                                                     *
c    updated:  april 17, 1990                                         *
c     reason:  include CCT option for ccsd(t) gradients               *
c
C    UPDATED:  november 12,1989                                       *
C     REASON:  rclose 4 itaps 81 82 83                                *
C                                                                     *
C    UPDATED:  MAY   10, 1989                                         *
C     REASON:  CHANGE NAME TO CCQG                                    *
C                                                                     *
C    UPDATED:  MARCH 28, 1989                                         *
C     REASON:  UPDATE SUBROUTINE READ30 (NEW SYMMETRY ORDERING)       *
C              INCLUDE SUBROUTINE FROZEN                              *
C              IMPLEMENT CORS AND VIRS                                *
C                                                                     *
C    UPDATED:  OCTOBER 18, 1988                                       *
C     REASON:  FORMGD: ELIMINATE 5 O3 V3 LOOPS                        *
C              FORMGA: ADD SYMMETRY TO 4 LAST LOOPS                   *
C                                                                     *
C    UPDATED:  JUNE 26 , 1988.                                        *
C     REASON:  INCLUDE D2 POINT GROUP CONSTANTS.                      *
C                                                                     *
C    UPDATED:  DECEMBER 17, 1987                                      *
C     REASON:  ALL LOOPS SYMMETRY ADAPTED                             *
C           :  Z2, T2 AND TAU SYMMETRY-PACKED                         *
C                                                                     *
C    UPDATED:  MARCH 7 , 1988.                                        *
C     REASON:  TAKE OUT FORMGC TO COMPILE IT LEVEL 2                  *
C                                                                     *
C   VERSION :  CCDMAT  -ORIGINAL- WRITTEN BY ANDY SCHEINER            *
C      DATE :  MAY  25   , 1987                                       *
C**********************************************************************
C
      IMPLICIT INTEGER (A-Z)
cets030491
#include <error.h>
      character*40 corrty,dertyp
      CHARACTER*4 OPTION,CLABEL(20),RSTR,CHAR
      character*80 label
      REAL*8 CC(maxcor),ECOR,ETOT,ENUC,ESCF
      INTEGER NOSYM(8),NVSYM(8),NFZO(8),NFZV(8),NDOC(8),NUOC(8)
      COMMON/TAPES/ITAP81,ITAP82,ITAP83,ITAP99,itap76,i76,itap74
C9370 DIMENSION IC(1000000),I30(200)
C3090 DIMENSION IC(3480000),I30(200)
      DIMENSION IC(*),I30(200)
cges  EQUIVALENCE (CC,IC)
C
cges  CALL NOUNFL
C
C9370 MAXCOR= 500000
C3090 MAXCOR=1740000
C
      call psinit('APPEND')
      CALL TSTART(6)
C
      INPUT=5
      JOUT=6
      WRITE(JOUT,6000)
 6000 FORMAT(   46('-'),/,46('-'),/,
     .' EFFECTIVE ONE- & TWO PDM FOR COUPLED CLUSTER',/,
     .'     WRITTEN BY GUSTAVO E. SCUSERIA'     ,/,
     .'      VS CCQG 2.0   april 17 1990 ',/,46('-'),/,46('-')   )
C
C >>> ZERO OUT THE CORE
C
      CALL ZERO(CC,MAXCOR)
C
      ITAP30=30
      ITAP68=68
      ITAP69=69
      itap74=74
      itap76=76
      ITAP81=81
      ITAP82=82
      ITAP83=83
      ITAP99=99
C
      CALL RFILE(ITAP68)
      CALL RFILE(ITAP69)
      CALL RFILE(ITAP81)
      CALL RFILE(ITAP82)
      CALL RFILE(ITAP83)
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
C
      EIGVAL = 1
      ZAN    = EIGVAL + NBF
      COORD  = ZAN + NBF
      NLAMDA = WPADTI(COORD + 3*NAT)
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
     .            IOPEN,IC(NOP))
C
      NV=NBF-NO
      NTRO=NO*(NO+1)/2
      NTRV=NV*(NV+1)/2
      NO2=NO*NO
      NV2=NV*NV
      NOV=MAX0(NO,NV)
      N2OV=MAX0(NO*NV,NTRV)
C
C >>> READ CCSD INPUT
C
c Make sure that the user isn't trying to use the old style input

      CALL LOCATE (INPUT,'# CCSD ###',IERR)
      if (ierr.eq.0) then
        write(6,*)'ERROR: old style ''# CCSD ###'' input has been found'
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

      errcod = frdc('LABEL',label)
      dertyp = 'FIRST '
      errcod = frdc('DERTYPE',dertyp)

      write(6,*)
      write(6,*) 'Parameters to be used:'
      write(6,*)    '  label        = ',label
      write(6,*)    '  wfn          = ',corrty
      write(6,*)    '  dertype      = ',dertyp
      write(6,6101) '  frozen_docc  = ',(nfzo(i),i=1,nirred)
      write(6,6101) '  frozen_uocc  = ',(nfzv(i),i=1,nirred)
      write(6,*)
6101  format(1X,A,8(1X,I3,1X))

c >>>>>>>>> Check the input and convert to the internal flags required by ccqg.
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


C
 6002 FORMAT(16I5)
 6004 FORMAT(20A4)
 6007 FORMAT(//,2X,' WARNING! NO INPUT TO CCSD HAS BEEN FOUND.',/,
     .          2X,'          DEFAULT PARAMETERS WILL BE USED.',//)
 6051 FORMAT(I2,1X,I2)
 6052 FORMAT(A4,1X,5X,A4)
C
C >>> WRITE OUT INPUT OPTIONS
C
      WRITE(JOUT,6022)label,OPTION,CHAR,NORD
 6022 FORMAT(//,2X,'***** CC PARAMETERS *****',
     .        /,2X,'LABEL  = ',A80,
     .        /,2X,'OPTION = ',A4,
     .        /,2X,'SYMM   = ',A4,I1)
C
      CALL FROZEN(IC(FLOV),NOSYM,NVSYM,NFZO,NFZV,NDOC,NUOC,
     .            IC(FZO),IC(FZV),IC(ITYP),NIRRED,NO,NV,JOUT,
     .            NOACT,NVACT,NSYMHF)
C
      CALL NCOUNT(IC(ORBSYM),IC(FLOV),NIRRED,NO,NV,NBF,
     .            IC(FZO),IC(FZV),NT1,NT2,
     .            IC(NT3),NTAU,OPTION)
C
C >>> DYNAMICAL ALLOCATION
C
      TOFF  = WPADTI(TOP)
      TADD  = TOFF + NO*NO*2
      IOFF  = TADD + NV*NV
      INO   = IOFF + NTR
      INV   = INO + NOV
      INTRO = INV + NOV
      INTRV = INTRO + N2OV
      ZLX   = INTRV + NO*NV
      AOFF  = ZLX + NV*NV
      AADD  = AOFF + NO*NO
      COFF  = AADD + NO*NO
      CADD  = COFF + NV*NV
      DOFF  = CADD + NO*NO
      DADD  = DOFF + NO*NV
      EOFF  = DADD + NO*NV
      EADD  = EOFF + NO*NV
      T1    = IADTWP(EADD + NO*NO)
      Z1    = T1 + NO*NV
      TOP   = Z1 + NO*NV
      IF(OPTION.EQ.'SDT1')THEN
         Y3    = TOP
         W3    = Y3 + NTAU
         TOP   = W3 + NTAU
         T3OFF = WPADTI(TOP)
         T3ADD = T3OFF + NO*NO*NIRRED*2
         MOO   = IADTWP(T3ADD + NV*NV*NIRRED)
         NVV   = MOO + NO*NO
         TOP   = NVV + NV*NV
      else if(option.eq.'CCT ')then
         write(6,*)'allocating moo and nvv for ccsd(t)'
         moo=top
         nvv=moo+no*no
         top =nvv+nv*nv
      ELSE
         Y3    = 1
         W3    = 1
         T3OFF = 1
         T3ADD = 1
         MOO   = 1
         NVV   = 1
      ENDIF
C
C >>> FORM SYMMETRY-PACKING ARRAYS AND OFFSETS
C
      CALL OFFSET(IC(INO),IC(INV),IC(INTRO),IC(INTRV),
     .            IC(ZLX),IC(IOFF),
     .            NBF,NO,NV,NOV,NTRO,NTRV,N2OV,NTR,JOUT)
C
      CALL SYMARR(IC(ORBSYM),IC(FLOV),NIRRED,NO,NV,NBF,NOV,NM,
     .            IC(FZO),IC(FZV),NFZO,NFZV,
     .            IC(TOFF),IC(TADD),IC(INO),IC(INV),NO2,NV2,
     .            IC(T3OFF),IC(T3ADD),OPTION,
     .            NSGOO,NSGOV,NSHOV,NSLOV,
     .            IC(AOFF),IC(AADD),IC(COFF),IC(CADD),
     .            IC(DOFF),IC(DADD),IC(EOFF),IC(EADD))
C
      WRITE(JOUT,6428)NSGOO,NSGOV,NSHOV,NSLOV
 6428 FORMAT(/,'LENGTH OF SYMMETRY-PACKED TWO-PDM',
     .       /,' GA TERMS = ',I8,
     .       /,' GC TERMS = ',I8,
     .       /,' GD TERMS = ',I8,
     .       /,' GE TERMS = ',I8)
      MAXPDM=NSGOO
      IF(NSGOV.GT.MAXPDM)MAXPDM=NSGOV
      IF(NSHOV.GT.MAXPDM)MAXPDM=NSHOV
      IF(NSLOV.GT.MAXPDM)MAXPDM=NSLOV
c
c >>> read tape 76. NOTE that order is important! moo,nvv,dden,eden
c
      if(option.eq.'CCT ')then
         call rfile(itap76)
         call wreadw(itap76,cc(moo),intowp(no*no),  1,i76)
c        write(6,*)'reading moo',no*no,'  long'
c        write(6,'(8f10.6)')(cc(moo+i-1),i=1,no*no)
c
         call wreadw(itap76,cc(nvv),intowp(nv*nv),i76,i76)
c        write(6,*)'reading nvv',nv*nv,'  long'
c        write(6,'(8f10.6)')(cc(nvv+i-1),i=1,nv*nv)
      endif
         
C
C >>> READ TAPE 69
C
      CALL SREAD(ITAP69,NM,1)
      CALL SREAD(ITAP69,NDIMT2,1)
C
C >>> MORE DYNAMICAL ALLOCATION
C
      T2  = TOP
      Z2  = T2 + NDIMT2
      TAU = Z2 + NDIMT2
      TOP = TAU + NDIMT2
      IF (TOP.GT.MAXCOR) THEN
         WRITE (JOUT,*) ' BEFORE READ TAPE 69 : NOT ENOUGH CC'
         WRITE (JOUT,*) ' TOP ',TOP,'  MAXCOR ',MAXCOR
         call qabort
      ENDIF
C
      CALL RGETSA(ITAP69,IT169)
      CALL SREAD(ITAP69,CC(T1),INTOWP(NO*NV))
C
      CALL RGETSA(ITAP69,IT269)
      CALL SREAD(ITAP69,CC(T2),INTOWP(NDIMT2))
      CALL FRMTAU(IC(ORBSYM),IC(FLOV),NIRRED,NO,NV,
     .            CC(TAU),CC(T2),NDIMT2,CC(T1),IC(FZO),IC(FZV))
C
      CALL RGETSA(ITAP69,IETOT)
      CALL SREAD(ITAP69,ETOT,INTOWP(1))
C
      CALL RGETSA(ITAP69,IEC)
      CALL SREAD(ITAP69,ECOR,INTOWP(1))
C
      CALL RGETSA(ITAP69,IZ169)
      CALL SREAD(ITAP69,CC(Z1),INTOWP(NO*NV))
C
      CALL RGETSA(ITAP69,IZ269)
      CALL SREAD(ITAP69,CC(Z2),INTOWP(NDIMT2))
C
      WRITE(6,*)
      WRITE(6,*)'ECOR FROM FILE69   ',ECOR
      WRITE(6,*)'ETOT FROM FILE69   ',ETOT
C
C >>> FINISH DYNAMICAL ALLOCATION
C
      QOO    = TOP
      QVV    = QOO + NTRO
      QOV    = QVV + NTRV
      GBUF   = SEC2I(100)/INTOWP(1)
      BUF    = QOV + NO*NV
      IBUF   = WPADTI(BUF)
      GMAT   = BUF + GBUF
      GAMAT  = GMAT
      GCMAT  = GMAT
      GDMAT  = GMAT
      GEMAT  = GMAT
      VALV2  = GMAT + MAXPDM
      VALOV  = VALV2 + NTRV
      T3Z2   = VALOV + NV*NO
      VALOV2 = T3Z2 + NO*NV*NV
      VALOV3 = VALOV2 + NV2*NO
      VALO3V = VALOV3 + NV2*NO
      Q1     = VALO3V + NV*NO2*NO
      Q2     = Q1 + NO*NV
      Q3     = Q2 + NO*NV
      TZ2    = Q3 + NO*NV
      VALO   = TZ2 + NO2
      VALV   = VALO + NO
      Q1A    = VALV + NV
      TZA    = Q1A + NO*NO
      TOP    = TZA + NO*NO
C
      ICOR   = TOP*8/1024
      IMCOR  = MAXCOR*8/1024
      WRITE(JOUT,6001)IMCOR,MAXCOR,ICOR,TOP
 6001 FORMAT(//,' MAXIMUM  CORE = ',I7,2X,'K BYTES  OR  ',I9 ,' REAL',
     .         ' WORDS',/,
     .          ' REQUIRED CORE = ',I7,2X,'K BYTES  OR  ',I9 ,' REAL',
     .         ' WORDS',//)
      IF (TOP.GT.MAXCOR) THEN
         WRITE (JOUT,*) ' NOT ENOUGH CORE'
         WRITE (JOUT,*) ' TOP ',TOP,'  MAXCOR ',MAXCOR
         call qabort
      ENDIF
C
      CALL FORMGD(CC(GDMAT),CC(TAU),CC(T2),CC(T1),CC(Z2),
     .            CC(Z1),NO,NV,NO2,NV2,NTRO,NTRV,IC(IOFF),NTR,
     .            CC(Q1),CC(Q2),CC(Q3),CC(TZ2),CC(VALO),
     .            CC(VALV),JOUT,OPTION,
     .            IC(ORBSYM),IC(FLOV),IC(ZLX),IC(INO),IC(INV),
     .            IC(NT3),NIRRED,NTAU,NOV,IC(T3OFF),IC(T3ADD),CC(W3),
     .            CC(MOO),CC(NVV),CC(Y3),NDIMT2,IC(TOFF),IC(TADD),
     .            IC(DOFF),IC(DADD),NSHOV,CC(Q1A),CC(TZA),
     .            IC(FZO),IC(FZV))
      CALL WRTGD (CC(BUF),IC(IBUF),GBUF,NBF,NO,NV,ITAP68,
     .            CC(GDMAT),NSHOV,IC(IOFF),NTR,JOUT,ICNT,
     .            IC(DOFF),IC(DADD),IC(INO),
     .            IC(FLOV),IC(ORBSYM),NIRRED,IC(FZO),IC(FZV))
C
      CALL FORMGA(CC(GAMAT),CC(TAU),CC(T2),CC(T1),CC(Z2),
     .            CC(Z1),NO,NV,NO2,NV2,NTRO,NTRV,IC(IOFF),NTR,JOUT,
     .            IC(ZLX),IC(FLOV),IC(ORBSYM),NIRRED,
     .            CC(Q1A),CC(TZA),OPTION,CC(MOO),
     .            NDIMT2,IC(TOFF),IC(TADD),
     .            IC(AOFF),IC(AADD),IC(INO),NSGOO,IC(FZO),IC(FZV))
      CALL WRTGA (CC(BUF),IC(IBUF),GBUF,NBF,NO,NV,ITAP68,
     .            CC(GAMAT),NSGOO,IC(IOFF),NTR,JOUT,ICNT,
     .            IC(AOFF),IC(AADD),IC(INO),
     .            IC(FLOV),IC(ORBSYM),NIRRED,IC(FZO),IC(FZV))
C
      CALL FORMGF(CC(BUF),IC(IBUF),GBUF,CC(TAU),CC(T2),
     .            CC(T1),CC(Z2),CC(Z1),NO,CC(VALOV),CC(VALOV2),
     .            NV,NO2,NV2,NTRO,NTRV,ICNT,IC(ZLX),CC(VALV2),
     .            CC(VALOV3),CC(VALO3V),IC(IOFF),NTR,CC(VALV),
     .            JOUT,ITAP68,OPTION,
     .            IC(ORBSYM),IC(FLOV),CC(T3Z2),IC(INO),IC(INV),
     .            IC(NT3),NIRRED,NTAU,NOV,IC(T3OFF),IC(T3ADD),
     .            CC(W3),CC(Y3),NDIMT2,IC(TOFF),IC(TADD),
     .            IC(FZO),IC(FZV))
C
      CALL FORMGB(CC(BUF),IC(IBUF),GBUF,CC(TAU),CC(Z2),
     .            NO,NV,NO2,NV2,NTRO,
     .            NTRV,ICNT,IC(IOFF),NTR,IC(ZLX),
     .            JOUT,ITAP68,IC(FLOV),IC(ORBSYM),NIRRED,
     .            NDIMT2,IC(TOFF),IC(TADD),IC(FZO),IC(FZV))
C
      CALL FORMGC(CC(GCMAT),CC(TAU),CC(T2),CC(T1),CC(Z2),
     .            CC(Z1),NO,NV,NO2,NV2,NTRO,NTRV,IC(IOFF),NTR,JOUT,
     .            IC(ZLX),IC(FLOV),IC(ORBSYM),NIRRED,
     .            OPTION,CC(NVV),NDIMT2,IC(TOFF),IC(TADD),
     .            IC(COFF),IC(CADD),NSGOV,IC(INO),IC(INV),NOV,
     .            IC(FZO),IC(FZV))
      CALL WRTGC (CC(BUF),IC(IBUF),GBUF,NBF,NO,NV,ITAP68,
     .            CC(GCMAT),NSGOV,IC(IOFF),NTR,JOUT,ICNT,
     .            IC(COFF),IC(CADD),IC(INO),IC(INV),NOV,
     .            IC(FLOV),IC(ORBSYM),NIRRED,IC(FZO),IC(FZV))
C
      CALL FORMGE(CC(GEMAT),CC(TAU),CC(T2),CC(T1),CC(Z2),
     .            CC(Z1),NO,NV,NO2,NV2,NTRO,NTRV,IC(IOFF),NTR,
     .            CC(Q1),CC(Q2),CC(Q3),CC(VALO),JOUT,OPTION,
     .            IC(ORBSYM),IC(FLOV),IC(ZLX),IC(INO),IC(INV),
     .            IC(NT3),NIRRED,NTAU,NOV,IC(T3OFF),IC(T3ADD),CC(W3),
     .            CC(Y3),NDIMT2,IC(TOFF),IC(TADD),
     .            IC(EOFF),IC(EADD),NSLOV,IC(FZO),IC(FZV))
      CALL WRTGE (CC(BUF),IC(IBUF),GBUF,NBF,NO,NV,ITAP68,
     .            CC(GEMAT),NSLOV,IC(IOFF),NTR,JOUT,ICNT,
     .            IC(EOFF),IC(EADD),IC(INO),NOV,
     .            IC(FLOV),IC(ORBSYM),NIRRED,IC(FZO),IC(FZV))
C
      CALL FORMQ (CC(QOO),CC(QVV),CC(QOV),CC(TAU),CC(T2),
     .            CC(T1),CC(Z2),CC(Z1),NO,NV,NO2,NV2,NTRO,NTRV,
     .            IC(IOFF),NTR,JOUT,OPTION,
     .            IC(ORBSYM),IC(FLOV),IC(ZLX),IC(INO),IC(INV),
     .            IC(NT3),NIRRED,NTAU,NOV,IC(T3OFF),IC(T3ADD),CC(W3),
     .            CC(MOO),CC(NVV),NDIMT2,IC(TOFF),IC(TADD),
     .            IC(FZO),IC(FZV))
      CALL WRTQQ (CC(BUF),IC(IBUF),GBUF,NBF,NO,NV,ITAP68,
     .            CC(QOO),NTRO,CC(QVV),NTRV,CC(QOV),IC(IOFF),
     .            NTR,JOUT,ICNT,IC(FZO),IC(FZV))
C
C
      WRITE(JOUT,6059)
 6059 FORMAT(//,18X,' >>> CONSTRUCTION OF ONE & TWO-PDM COMPLETED <<<')
c     CALL RCLOSE(ITAP30,3)
      CALL RCLOSE(ITAP68,3)
      CALL RCLOSE(ITAP69,3)
      if(option.eq.'CCT ')call rclose(itap74,3)
      if(option.eq.'CCT ')call rclose(itap76,3)
clean CALL RCLOSE(ITAP81,3)
clean CALL RCLOSE(ITAP82,3)
clean CALL RCLOSE(ITAP83,3)
      CALL RCLOSE(ITAP81,4)
      CALL RCLOSE(ITAP82,4)
      CALL RCLOSE(ITAP83,4)
      CALL RCLOSE(ITAP99,4)
      CALL TSTOP(6)
      call psdone()
      stop
C
      END
