cges  PROGRA CCTR
      subroutine fentry(cc,ic,maxcor)
C
C >>>  LAST UPDATE MAY 10, 1989  <<<
C
C***********************************************************************
C  NOTICE OF PROGRAM MODIFICATIONS                                     *
c     by:  Daniel Crawford		search: Crawdad.1.25.94
c   date:  1/25/94
c reason:  Make sure cctr dies when nirred != nsymhf for open-shells.
c
c     by:  Edward Seidl                  search: cets070991
c   date:  July 9,1991
c reason:  hardwire to use <= 8MB of core.  option LARGE_CORE added
c          will use amount of core specified by MEMORY if LARGE_CORE is
c          true
c
C      BY:  GUSTAVO E. SCUSERIA                                        *
c
c    date:  november 16, 1990
c  reason:  under new 2.11 compiler (installed with os 4.51) this code *
c           does not work if -O2 is used. Calculations between nov 9
c           17:27 and november 16 1:00 pm are screwed up. VIVA mips!   *
c
c    date:  august 23, 1990 
c  reason:  symorb in read30 should be integer (caused problems on 530)*
c
c    date:  january 17, 1990
c  reason:  add warning concerning nirred.ne.nsymhf                    *
c
c    date:  october 27, 1989
c  reason:  add # CCTR ### input to read code34                        *
C                                                                      *
C    DATE:  MAY 10,  1989.                                             *
C  REASON:  CHANGE NAME TO CCTR                                        *
C                                                                      *
C    DATE:  APRIL 13, 1989                                             *
C  REASON:  CHANGE READ30 ROUTINE TO NEW VERSION (NEW SYMMETRY ORDER)  *
C                                                                      *
C    DATE:  SEPTEMBER 16, 1988                                         *
C  REASON:  FORCE NBFAO=NBF FOR NEW GCINTS CODE.                       *
C           BECAUSE NBFAO IS ACTUALLY NBFSO !
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include <error.h>
      integer icore
      REAL*8 CC(maxcor) 
      INTEGER IC(*),I30(200),code34
      INTEGER WPADTI,SEC2I,SPACE,SPACE2
      COMMON/BASIS/NBASIS,MTRI,NSYM
      integer frdboo, errcod, d34, prcntr, prdup
      call psinit('APPEND')
c     call vadvis(1,cc(1),cc(maxcor))
      errcod = prcntr('ON BRIEF')
      errcod = prdup()
      if (prcntr('IS_ON DEBUG').eq.1) then
        errcod = prcntr('OFF BRIEF')
        endif
      CALL TSTART(6)
C
      WRITE(6,6000)
clj Conserve paper
c6000 FORMAT(//,49('-'),/,49('-'),/,
c    .' CCTR(ANSL) : INTEGRAL TRANSFORMATION FOR CC CODES',/,
c    .'        WRITTEN BY ANDY SCHEINER     ',/,
c    .'     LAST UPDATE MAY 10, 1989 BY GES ',/,49('-'),/,49('-'),//)
 6000 FORMAT(/,
     .' CCTR(ANSL) : INTEGRAL TRANSFORMATION FOR CC CODES',/,
     .'        WRITTEN BY ANDY SCHEINER     ',/,
     .'     LAST UPDATE MAY 10, 1989 BY GES ',/)
C
    1 FORMAT ('                                                ')
    2 FORMAT ('************************************************')
  121 FORMAT (' NOT ENOUGH CC TO DO FINAL 2-ELECTRON HALF TRANSFORM')
C
      INPUT=5
      JOUT=6
      KOUT=8
      ITAP30=30
      ITAP71=71
      ITAP72=72
      ITAP77=77
      ITAP78=78
      ITAP34=34
      IBUFHF=4471
      IFLAG=1
      ISORT=1
      JBUF=0
      INTBUF=SEC2I(100)/INTOWP(1)
      MAXBKT=120
C
      CALL RFILE(ITAP71)
      CALL RFILE(ITAP72)
      CALL RFILE(ITAP77)
      CALL RFILE(ITAP78)
      CALL RFILE(ITAP30)
C
C >>> READ POINTERS IN FILE30
C
      CALL WREADW(ITAP30,I30,200,101,LWORD)
      MPOINT = I30( 2)
      MCONST = I30( 3)
      MCALCS = I30( 4)
      NCALCS = I30( 5)
      NBF    = I30(18)
      NAT    = I30(19)
      NTAO   = I30(22)
      NSHELL = I30(27)
      NIRRED = I30(28)
      NPRIM  = I30(32)
      NSYMHF = I30(41)
      MXCOEF = I30(42)
      IOPEN  = I30(43)
      if (prcntr('IS_ON BRIEF').eq.0) then
      write(6,*)
      write(6,*) 'nbf =',nbf
      write(6,*) 'nat =',nat
      write(6,*) 'nshell =',nshell
      write(6,*) 'nirred =',nirred
      write(6,*) 'nsymhf =',nsymhf
      write(6,*) 'nprim =',nprim
      write(6,*)
      endif
c
      if (nirred.ne.nsymhf)then
         write(6,*)
         write(6,*)' WARNING ...      WARNING !!!!'
         write(6,*)' nirred not equal nsymhf'
         write(6,*)' there is an irrep with no orbitals'
Crawdad.1.25.94         write(6,*)' this has produced problems in D2h'
Crawdad.1.25.94         write(6,*)' lower your symmetry or take the risk!'
         write(6,*)' this has produced problems in high symmetry'
         write(6,*)' lower your symmetry if this is open-shell!'
       endif
       if((iopen.ne.0).and.(nirred.ne.nsymhf)) then
         write(6,*)
	 write(6,*) ' Cannot continue for open-shells! '
         write(6,*) ' Error: exiting cctr '
         write(6,*)
         call fioflu(6)
	 call qabort
       endif
Crawdad.1.25.94
C
C >>> FORCE NBFAO = NBF FOR 5D, 7F AND 9G
C
      NBFAO=NBF
      NTRI=(NBF*(NBF+1))/2
      NTRIAO=(NBFAO*(NBFAO+1))/2
      MTRI=NTRIAO
c
cets070991
c    test to see if using more than 8mb of core
c    if MEMORY is > 8mb, then reduce it unless LARGE_CORE is true
c
#ifdef AIX
      icore=0
      errcod = frdboo('LARGE_CORE',icore)
      if(maxcor.gt.1000000.and.icore.eq.0) then
        write(6,*) ' '
        write(6,*) 'You have requested more than 8 MB of memory'
        write(6,*) 'which will cause poor performance on the'
        write(6,*) 'rs6000s, thus I will only use 8 MB of core.'
        write(6,*) 'If this is not enough, may I suggest using'
        write(6,*) 'the LARGE_CORE option (see psi manual)'
        write(6,*) ' '
        maxcor=1000000
      endif
#endif
cets070991
C
C >>> ZERO OUT THE CC
C
      CALL ZERO(CC,MAXCOR)
C
C >>> DYNAMICAL ALLOCATION
C
      K1     = 1
      K2     = K1 + NBF*NBFAO
      K3     = WPADTI(K2 + NBF)
      KTPO   = K3 + NSYMHF
      KTYP   = KTPO + NBF
      KFLOV  = KTYP + NSYMHF
      KORBSM = KFLOV + 2*NIRRED
      KNMSM  = KORBSM + NBF
      KIWRK  = KNMSM + NIRRED
      KKLPT  = KIWRK + NIRRED
      KIJPT  = KKLPT + NTRI
      KOFF1  = KIJPT + NTRI
      KOFF2  = KOFF1 + NTRI
      K4     = KOFF2 + NTRI
      K5     = K4 + MCALCS
      K6     = IADTWP(K5 + NSYMHF)
      K7     = K6 + MXCOEF
      KTOP   = K7 + NBF*NBF - 1
C
c     write(6,*)'calling read30'
c     call flush(6)
      CALL READ30 (ITAP30,MPOINT,MCONST,MCALCS,NCALCS,NBF,NBFAO,NSYMHF,
     .             MXCOEF,CC(K1),CC(K2),IC(K3),IC(K4),CC(K6)
     .             ,JOUT,IC(K5),NO,IC(KTPO),CC(K7),IC(KTYP),
     .             NIRRED,IC(KFLOV),IC(KORBSM),NTRI,IC(KNMSM),
     .             IC(KIWRK),IC(KKLPT),IC(KIJPT),IC(KOFF1),
     .             IC(KOFF2))
c     write(6,*)'back from read30'
c     call flush(6)
      NV=NBF-NO
c
      call locate(5,'# CCTR ####',ierr)
      if(ierr.ne.0) then
        errcod = frdboo('DELETE34',d34)
        if (errcod.ne.EOK) d34=1
        write(6,*)
        if (d34.eq.1) then
          write(6,*) 'delete34 = YES'
        else
          write(6,*) 'delete34 = NO'
          endif
        write(6,*)
        if (d34.eq.1) code34=4
        if (d34.ne.1) code34=3
      else
        write(6,*) 'WARNING: using old style ''# CCTR ####'' input'
        code34=0
        read(5,621)code34
  621   format(16i5)
        if(code34.ne.3.and.code34.ne.4)code34=4
        write(6,*)
        write(6,*)'code34 =',code34
        write(6,*)
        endif
c
c     WRITE(*,*)  '  NV = ',NV
c     WRITE(JOUT,*)
c     WRITE(JOUT,*) ' THE SCF EIGENVECTORS '
c     WRITE(JOUT,*)
c     CALL MATOUT (CC(K1),NBFAO,NBF,NBFAO,NBF,JOUT)
c     WRITE (JOUT,1)
C     WRITE (JOUT,2)
C     WRITE (JOUT,*) ' ************  EIGENVECTOR  *************'
C     WRITE (JOUT,2)
C     WRITE (JOUT,1)
C     DO 4 I=1,196
C        WRITE (JOUT,*) I,CC(I)
C   4 CONTINUE
c     WRITE (JOUT,1)
c     WRITE (JOUT,2)
c     WRITE (JOUT,*) ' ************  EIGENVALUES  *************'
c     WRITE (JOUT,2)
c     WRITE (JOUT,1)
c     call flush(6)
C     DO 5 I=197,210
C        WRITE (JOUT,*) I,CC(I)
C   5 CONTINUE
C     WRITE (JOUT,1)
C     WRITE (JOUT,2)
C     WRITE (JOUT,*) ' ************  NLAMDA  *************'
C     WRITE (JOUT,2)
C     WRITE (JOUT,1)
C     WRITE (JOUT,*) IC (421)
C     WRITE (JOUT,1)
C     WRITE (JOUT,2)
C     WRITE (JOUT,*) ' ************  ITEMP  *************'
C     WRITE (JOUT,2)
C     WRITE (JOUT,1)
C     DO 6 I=422,521
C        WRITE (JOUT,*) I,IC(I)
C   6 CONTINUE
C     WRITE (JOUT,1)
C     WRITE (JOUT,2)
C     WRITE (JOUT,*) ' ************  XVEC  *************'
C     WRITE (JOUT,2)
C     WRITE (JOUT,1)
C     DO 7 I=262,457
C        WRITE (JOUT,*) I,CC(I)
C   7 CONTINUE
C
      J1   = IADTWP(K4)
      J2   = J1 + NTRIAO
      J3   = J2 + NTRIAO
      J4   = J3 + NTRIAO
      J5   = J4 + NTRIAO
      J6   = WPADTI(J5)
      JTOP = J5 + IBUFHF-1
      IF (JTOP.GT.MAXCOR) THEN
         WRITE (JOUT,*) ' ERROR: out of memory'
         WRITE (JOUT,*) ' JTOP=',JTOP
         WRITE (JOUT,*) ' MAXCOR=',MAXCOR
         WRITE (JOUT,*) ' increase size of ''MEMORY'''
         call fioflu(6)
         call qabort
      END IF
C
c     write(6,*)'calling read34'
c     call flush(6)
      CALL READ34 (CC(J1),CC(J2),CC(J3),CC(J4),CC(J5),IC(J6
     .            ),IBUFHF,ITAP34,JOUT)
C     WRITE (JOUT,1)
C     WRITE (JOUT,2)
C     WRITE (JOUT,*) ' ************  H MATRIX  *************'
C     WRITE (JOUT,2)
C     WRITE (JOUT,1)
C     DO 8 I=212,316
C        WRITE (JOUT,*) I,CC(I)
C   8 CONTINUE
C     WRITE (JOUT,1)
C     WRITE (JOUT,2)
C     WRITE (JOUT,*) ' ************  S MATRIX  *************'
C     WRITE (JOUT,2)
C     WRITE (JOUT,1)
C     DO 9 I=317,421
C        WRITE (JOUT,*) I,CC(I)
C   9 CONTINUE
C     WRITE (JOUT,1)
C     WRITE (JOUT,2)
C     WRITE (JOUT,*) ' ************  T MATRIX  *************'
C     WRITE (JOUT,2)
C     WRITE (JOUT,1)
C     DO 11 I=422,526
C        WRITE (JOUT,*) I,CC(I)
C  11 CONTINUE
C     WRITE (JOUT,1)
C     WRITE (JOUT,2)
C     WRITE (JOUT,*) ' ************  V MATRIX  *************'
C     WRITE (JOUT,2)
C     WRITE (JOUT,1)
C     DO 12 I=527,631
C        WRITE (JOUT,*) I,CC(I)
C  12 CONTINUE
CC    WRITE (JOUT,1)
CC    WRITE (JOUT,2)
CC    WRITE (JOUT,*) ' ************  XVEC  *************'
CC    WRITE (JOUT,2)
CC    WRITE (JOUT,1)
CC    DO 7 I=262,457
CC       WRITE (JOUT,*) I,CC(I)
CC  7 CONTINUE
C
      LENINT=MAXCOR-IBUFHF-J2-3*NBFAO*NBFAO-(MAXBKT+1)/INTOWP(1)
C     LENINT=MAXCOR-INTBUF-J2-3*NBFAO*NBFAO-(MAXBKT+1)/INTOWP(1)
      if (prcntr('IS_ON BRIEF').eq.0)
     &  WRITE (JOUT,*) ' LENINT=',LENINT
CGES
      IF(LENINT.LE.0)THEN
         write(6,*)'ERROR: out of memory'
         WRITE(6,*)' LENINT.LE.0'
         WRITE(6,*)' increase size of ''MEMORY'''
         call fioflu(6)
         call qabort
      ENDIF
CGES
      ITRI=LENINT/NTRIAO
C     WRITE (JOUT,*) ' ITRI=',ITRI
C     IF (ITRI.GE.NTRIAO*2) THEN
C        WRITE(*,*)'DOING TRANSFORMATION IN CC'
C        L1=I2
C        L2=L1+NTRIAO*NTRIAO
C        L3=WPADTI(L2)
C        LTOP=IADTWP(L3+IBUFHF)-1
C        CALL ZERO(CC(L1),MAXCOR-L1)
C        WRITE(*,*) '  ZEROING OUT CC BEFORE RDTWO '
C
C     CALL RDTWO(CC(L1),CC(L2),IC(L3),IBUFHF,NTRIAO,ITAP34,JOUT)
C
C        M1=L2
C        IM1=WPADTI(M1)
C        M2=M1+NTRIAO*NTRI
C        M3=M2+NBFAO*NBFAO
C        M4=M3+NBFAO*NBFAO
C        MTOP=M4+NBFAO*NBFAO-1
C        IF (MTOP.GT.MAXCOR) THEN
C           WRITE(JOUT,*) ' MTOP EXCEEDS MAXCOR'
C        END IF
C        IFLAG=2
C
C        WRITE (JOUT,*) ' NTRIAO*NTRI=',NTRIAO*NTRI
C        CALL HTRANS(CC(M2),CC(M3),CC(M4),CC(L1),
C    .               CC(M1),IC(IM1),
C    .               NTRIAO*NTRIAO,NTRIAO*NTRI,NBFAO,CC(K1),NBFAO,NBF,
C    .               JOUT,NTRIAO,NTRI,IFLAG,ISORT,ITAP78,NTRIAO,1,1,
C    .               JBUF,1,IC(KORBSM),IC(KFLOV),NIRRED,IC(K3),
C    .               NSYMHF,IC(KKLPT),IC(KIJPT),NTRI,IC(KTPO),
C    .               IC(KOFF1),IC(KOFF2),NO)
C
C        M5=L1
C        IM5=WPADTI(M5)
C        IF (MTOP.GT.MAXCOR) THEN
C           WRITE(JOUT,*) ' MTOP EXCEEDS MAXCOR'
C        END IF
C        IFLAG=1
C
C        WRITE (JOUT,*) ' NTRI*NTRI=',NTRI*NTRI
C      CALL HTRANS(CC(M2),CC(M3),CC(M4),CC(M1),
C    .             CC(M5),IC(IM5),NTRIAO
C    .             *NTRI,NTRI*NTRI,NBFAO,CC(K1),NBFAO,NBF,JOUT,NTRI,
C    .              NTRI,IFLAG,ISORT,ITAP78,NTRI,1,1,JBUF,1,
C    .              IC(KORBSM),IC(KFLOV),NIRRED,IC(K3),NSYMHF,
C    .              IC(KKLPT),IC(KIJPT),NTRI,IC(KTPO),
C    .              IC(KOFF1),IC(KOFF2),NO)
C
C
C        CALL WRINTS(ITAP71,CC(M4),NTRI*NTRI,IAD)
C        CALL SREW(ITAP71)
C        CALL RCLOSE(ITAP71,3)
C        CALL SREW(ITAP30)
C        CALL RCLOSE(ITAP30,3)
C        WRITE (JOUT,*) ' IAD=',IAD
C     ELSE
         ISORT=2
         IF (ITRI.GT.NTRIAO) ITRI=NTRIAO
         NBUCK=NTRIAO/ITRI
         MODTRI=MOD(NTRIAO,ITRI)
         IF (ITRI*NBUCK.NE.NTRIAO) NBUCK=NBUCK+1
C        WRITE (JOUT,*) ' NBUCK=',NBUCK
         IBKTSP=MAXCOR-IBUFHF-J2
     .            -(MAXBKT+1)/INTOWP(1)
C        WRITE (JOUT,*) ' IBKTSP=',IBKTSP
         LENGTH=IBKTSP/NBUCK
         NIOBF=INTOWP(LENGTH)/SEC2I(1)
         IF (NIOBF.GT.100) NIOBF=100
         LENGTH=NIOBF*(SEC2I(1)/INTOWP(1))
         if (prcntr('IS_ON BRIEF').eq.0)
     &     WRITE (JOUT,*) ' LENGTH=',LENGTH
         SPACE=MAXCOR-LENGTH-3*NBFAO*NBFAO-J2
     .         -(MAXBKT+1)/INTOWP(1)
         IF (NTRIAO*ITRI.GT.SPACE) THEN
            LENGTH=LENGTH-(NTRIAO*ITRI-SPACE)
            NIOBF=INTOWP(LENGTH)/SEC2I(1)
            LENGTH=NIOBF*(SEC2I(1)/INTOWP(1))
C        WRITE (JOUT,*) ' LENGTH=',LENGTH
         SPACE=MAXCOR-LENGTH-3*NBFAO*NBFAO-J2
     .         -(MAXBKT+1)/INTOWP(1)
         END IF
         IF (SPACE.GT.NTRIAO*NTRIAO) SPACE=NTRIAO*NTRIAO
         if (prcntr('IS_ON BRIEF').eq.0) WRITE (JOUT,*) ' SPACE=',SPACE
CGES
      IF(SPACE.LE.0)THEN
         write(6,*) 'ERROR: out of memory'
         WRITE(6,*)' SPACE.LE.0'
         WRITE(6,*)' increase size of ''MEMORY'' or decrease niobf'
         call fioflu(6)
         call qabort
      ENDIF
CGES
C
         N1=J2
         IN1=WPADTI(N1)
         N2=IADTWP(IN1+MAXBKT)
         IN2=WPADTI(N2)
         N3=N2+NBUCK*LENGTH
         IN3=WPADTI(N3)
         NTOP=N3+IBUFHF-1
C        WRITE (JOUT,*) ' ENTERING SORT 1'
         CALL ZERO(CC(N1),MAXCOR-N1)
C
c       write(6,*)'calling sort'
c       call flush(6)
        CALL SORT(CC(N2),IC(IN2),IC(IN3),CC(N3),NBUCK*LENGTH,
     .     IBUFHF,NBUCK,JOUT,ITAP34,ITRI,ITAP77,IC(IN1),MAXBKT,NIOBF
     .            ,LENGTH,MAXVAL,IWRIT)
C
C        WRITE (JOUT,*) ' EXITING SORT 1'
c        CALL RCLOSE(ITAP34,3)
         CALL RCLOSE(ITAP34,code34)
c        CALL RCLOSE(ITAP34,4)
C
         N6=N2
         N7=N6+SPACE
         IN7=WPADTI(N7)
         N8=N7+LENGTH
         N9=N8+NBFAO*NBFAO
         N10=N9+NBFAO*NBFAO
         NTOP=N10+NBFAO*NBFAO-1
         CALL ZERO(CC(N6),MAXCOR-N6)
         JBUF=0
         DO 101 I=1,NBUCK
C           WRITE (JOUT,*) ' ENTERING BCKCHN 1'
C           WRITE (JOUT,*) ' NBUCK=',I
C           WRITE (JOUT,*) ' LENGTH=',LENGTH
C           WRITE (JOUT,*) ' SPACE=',SPACE
C           WRITE (JOUT,*) ' MAXBKT=',MAXBKT
            CALL BCKCHN(IC(IN1),I,ITAP77,CC(N7),IC(IN7),LENGTH,
     .                  MAXBKT,MAXVAL,CC(N6),SPACE,NTRIAO,JOUT,ITRI)
            NTRIMX=ITRI
            IF ((I.EQ.NBUCK).AND.(MODTRI.NE.0)) NTRIMX=MODTRI
            IFLAG=2
            CALL HTRANS(CC(N8),CC(N9),CC(N10),CC(N6),
     .               CC(N7),IC(IN7),
     .               SPACE,LENGTH,NBFAO,CC(K1),NBFAO,NBF,JOUT,NTRIAO
     .               ,NTRI,IFLAG,ISORT,ITAP78,NTRIMX,I,ITRI,JBUF,
     .               NBUCK,IC(KORBSM),IC(KFLOV),NIRRED,IC(K3),
     .               NSYMHF,IC(KKLPT),IC(KIJPT),NTRI,IC(KTPO),
     .               IC(KOFF1),IC(KOFF2),NO)
            CALL ZERO(CC(N6),MAXCOR-N6)
  101    CONTINUE
         CALL RCLOSE(ITAP77,4)
         CALL RFILE(ITAP77)
         CALL ZERO(CC(J2),MAXCOR-J2)
         NBUCK=NTRI/ITRI
         MODTRI=MOD(NTRI,ITRI)
         IF (ITRI*NBUCK.NE.NTRI) NBUCK=NBUCK+1
         IBKTSP=MAXCOR-LENGTH-J2
     .            -(MAXBKT+1)/INTOWP(1)
C        WRITE (JOUT,*) ' IBKTSP=',IBKTSP
         LEN2=IBKTSP/NBUCK
         NIOBF=INTOWP(LEN2)/SEC2I(1)
         IF (NIOBF.GT.100) NIOBF=100
         LEN2=NIOBF*(SEC2I(1)/INTOWP(1))
C        WRITE (JOUT,*) '   *******  BEFORE SORT2 *******'
C        WRITE (JOUT,*) ' LEN2=',LEN2
C        WRITE (JOUT,*) ' LENGTH=',LENGTH
C        WRITE (JOUT,*) ' NBUCK=',NBUCK
         SPACE2=MAXCOR-LEN2-3*NBF*NBFAO-J2
     .         -(MAXBKT+1)/INTOWP(1)
         IF (NTRIAO*ITRI.GT.SPACE2) THEN
            LEN2=LEN2-(NTRIAO*ITRI-SPACE2)
            NIOBF=INTOWP(LEN2)/SEC2I(1)
            LEN2=NIOBF*(SEC2I(1)/INTOWP(1))
            if (prcntr('IS_ON BRIEF').eq.0) WRITE (JOUT,*) ' LEN2=',LEN2
         SPACE2=MAXCOR-LEN2-3*NBF*NBFAO-J2
     .         -(MAXBKT+1)/INTOWP(1)
         END IF
         IF (SPACE2.GT.NTRI*NTRIAO) SPACE2=NTRI*NTRIAO
         if (prcntr('IS_ON BRIEF').eq.0) WRITE (JOUT,*) ' SPACE2=',SPACE2
CGES
      IF(SPACE2.LE.0)THEN
         write(6,*)'ERROR: out of memory'
         WRITE(6,*)' SPACE2.LE.0'
         WRITE(6,*)' increase size of ''MEMORY'' or decrease niobf'
         call fioflu(6)
         call qabort
      ENDIF
CGES
C        IF (INTBUF.GT.LEN2) THEN
C           WRITE (JOUT,121)
C           STOP
C        END IF
 
C
         N1=J2
         IN1=WPADTI(N1)
         N2=IADTWP(IN1+MAXBKT)
         IN2=WPADTI(N2)
         N3=N2+LEN2*NBUCK
         IN3=WPADTI(N3)
         NTOP=N3+LENGTH-1
         CALL ZERO(CC(N1),MAXCOR-N1)
C
         CALL SORT2(ITAP77,ITAP78,CC(N3),IC(IN3),LENGTH,CC(N2),
     .              IC(IN2),NBUCK*LEN2,ITRI,IC(IN1),MAXBKT,
     .              JOUT,NBUCK,NIOBF,LEN2,JBUF,MAXVAL)
C
         CALL RCLOSE(ITAP78,4)
         CALL RFILE(ITAP78)
         CALL SWRIT(ITAP78,LEN2,1)
         CALL SREW(ITAP78)
         CALL SREAD(ITAP78,LEN2,1)
         if (prcntr('IS_ON BRIEF').eq.0) WRITE (JOUT,*) ' LEN2=',LEN2
         N6=N2
         N7=N6+SPACE2
         IN7=WPADTI(N7)
         N8=N7+LEN2
         N9=N8+NBFAO*NBFAO
         N10=N9+NBFAO*NBFAO
         NTOP=N10+NBFAO*NBFAO-1
         CALL ZERO(CC(N6),MAXCOR-N6)
         JBUF=0
         DO 111 I=1,NBUCK
C           WRITE (JOUT,*) ' LEN2=',LENGTH
C           WRITE (JOUT,*) ' SPACE2=',SPACE2
C           WRITE (JOUT,*) ' MAXBKT=',MAXBKT
            CALL BCKCHN(IC(IN1),I,ITAP77,CC(N7),IC(IN7),LEN2,
     .                  MAXBKT,MAXVAL,CC(N6),SPACE2,NTRIAO,JOUT,ITRI)
            NTRIMX=ITRI
            IF ((I.EQ.NBUCK).AND.(MODTRI.NE.0)) NTRIMX=MODTRI
            IFLAG=1
            CALL SETMBF(N8,N7+LEN2)
            CALL SETMBF(N9,N8+NBFAO*NBFAO)
            CALL SETMBF(N10,N9+NBFAO*NBFAO)
            CALL SETMBF(NTOP,N10+NBFAO*NBFAO)
C
            CALL HTRANS(CC(N8),CC(N9),CC(N10),CC(N6),
     .               CC(N7),IC(IN7),
     .               SPACE2,LEN2,NBFAO,CC(K1),NBFAO,NBF,JOUT,NTRI
     .               ,NTRI,IFLAG,ISORT,ITAP78,NTRIMX,I,ITRI,JBUF,
     .               NBUCK,IC(KORBSM),IC(KFLOV),NIRRED,IC(K3),
     .               NSYMHF,IC(KKLPT),IC(KIJPT),NTRI,IC(KTPO),
     .               IC(KOFF1),IC(KOFF2),NO)
            CALL ZERO(CC(N6),MAXCOR-N6)
  111    CONTINUE
         CALL RCLOSE(ITAP77,4)
C     END IF
      I1=J2
      II1=WPADTI(I1)
      I2=I1+LEN2
      I3=I2+NBFAO*NBFAO
      I4=I3+NBFAO*NBFAO
      ITOP=I4+NBFAO*NBFAO-1
C
      CALL ZERO(CC(I1),MAXCOR-I1)
      ISORT=1
      IFLAG=1
      CALL HTRANS(CC(I2),CC(I3),CC(I4),CC(J1),
     .            CC(I1),IC(II1),NTRIAO
     .            *1,LEN2,NBFAO,CC(K1),NBFAO,NBF,JOUT,1,1,IFLAG,
     .            ISORT,ITAP78,1,1,1,JBUF,1,IC(KORBSM),IC(KFLOV),
     .           NIRRED,IC(K3),NSYMHF,IC(KKLPT),IC(KIJPT),NTRI,
     .            IC(KTPO),IC(KOFF1),IC(KOFF2),NO)
C     WRITE (JOUT,1)
C     WRITE (JOUT,2)
C     WRITE (JOUT,*) ' ************  TINTS  *************'
C     WRITE (JOUT,2)
C     WRITE (JOUT,1)
C     CALL PRINT(CC(I1),NTRI,NBF,JOUT)
C     DO 13 I=317,421
C        WRITE (JOUT,*) I,CC(I)
C  13 CONTINUE
C
C     CALL WRINTS(ITAP71,CC(I1),NTRI,IAD)
C     WRITE (JOUT,*) ' IAD=',IAD
C     CALL REINTS(ITAP71,CC(ITOP),NTRI,JOUT,IAD)
C
      CALL RCLOSE(ITAP30,3)
      CALL RCLOSE(ITAP78,3)
C     WRITE (JOUT,*) ' CALLED RCLOSE TAPE 78'
      CALL TSTOP(6)
      call psdone()
c
      stop
C
      END
