C bmatv2     FINAL MODIFIED VERSION OF PULAY'S BMATIN PROGRAM
      subroutine fentry(ra,ia,maxcor)
C***********************************************************************
C   NOTICE OF PROGRAM MODIFICATION                                     *
C***********************************************************************
c      by:  Edward Seidl (cets072991)
c    date:  July 29, 1991
c  reason:  replace SYMM vector with SYMM and ASYMM vectors
c           also, read intco information from input.dat first
c           keep intco.dat for now for compatability
c
c      by:  Edward Seidl (cets050591)
c    date:  May 5, 1991
c  reason:  merge ibm and unix versions
c           1) totally redo input to use new input parsing routines
c              written by Curt Janssen.
c              this requires massive changes to machb and fifth
c              and the addition of routines rdsimp and rdsymm
c
c           2) read geometry from file11.dat, write new geom
c              to geom.dat
c           3) read initial force constants from fconst.dat
c           4) read gradient from file11.dat so don't need bmwrta
c              any more
c           5) get internal coordinates from intco.dat
c           6) dynamically allocate most arrays so that there will
c              no longer be a hard limit on the number of atoms
c           7) write to output.dat rather than resul1.dat
c
C      BY:  Richard Remington                      Search: c010990     *
C    DATE:  JANUARY  9, 1990                                           *
C  REASON:  CHANGE "TPI" 2 pi to 6.0D0 as directed by Tracy Hamilton.  *
C           Done to assure convergence for torsion problems.           *
C***********************************************************************
C      BY:  TPH                                                        *
C    DATE:  JANUARY 12, 1989                                           *
C  REASON:  CHANGE DIMENSIONS AS PER YATES SUGGESTION TO 50 ATOMS      *
C           CHANGE TORSIONS IN MACHB SO THAT DISCONTINUITY HANDLED.    *
C           CHANGE FORMAT TO ACCOMODATE FORCE CONSTANTS OF 999 TO -99  *
C           CHANGE DUMMY ATOM PRINT TO PRETTY FORM.                    *
C***********************************************************************
C      BY:  TPH                                                        *
C    DATE:  DECEMBER 19, 1988                                          *
C  REASON:  READ IN NUMBER OF TIMES TO ATTEMPT TO SCALE DISPLACEMENTS  *
C           FOR GEOMETRY OPTIMIZATIONS.  DEFAULT IS 5.  OPTION IS NTRY *
C           THIS WAS TO ACCOUNT FOR CASES WHERE THE INTERNAL COORD.    *
C           CONVERGED TO A VALUE DIFFERENT FROM THAT EXPECTED, BUT     *
C           ONLY DIFFERING BY NUMERICAL ERROR.                         *
C***********************************************************************
C      BY:  TPH                                                        *
C    DATE:  DECEMBER 12, 1988                                          *
C  REASON:  CORRECT PROBLEM IN DISPLACEMENTS FOR TORSIONS.             *
C***********************************************************************
C      BY:  TPH                                                        *
C    DATE:  SEPT 9, 1988                                               *
C  REASON:  CHANGE SCALING OF DISPLACEMENT TO HALVE 5 TIMES. DISP=TRUE *
C***********************************************************************
C      BY:  TPH                                                        *
C    DATE:  JUNE 13,1988                                               *
C  REASON:  PUT IN INVERSE BOND LENGTH AS A COORDINATE  -  INVR        *
C***********************************************************************
C      BY:  TPH                                                        *
C    DATE:  JUNE 3,1988                                                *
C  REASON:  CHANGE SO THAT DUMMY ATOMS CAN BE USED FOR LINEAR BEND     *
C           DISPLACEMENTS.  THE PROGRAM ZOOMS THE DUMMY ATOMS TO A     *
C           DISTANCE OF 1 BILLION SO THAT LIN1 AND LIN2 ARE EQUIVALENT *
C           FOR LINEAR MOLECULES.                                      *
C           ALSO SCALE DISPLACEMENT ONCE IN GEOMETRY OPTIMIZATIONS     *
C***********************************************************************
C      BY:  TPH                                                        *
C    DATE:  APRIL 25, 1988                                             *
C  REASON:  CHANGE SOME OF THE DEFAULTS. NO LONGER HAVE OLDF, OLDH,    *
C           AND EXFI.  THE PROGRAM WILL READ THEM AND NOT CRASH.       *
C           NOUP AND NOEX MUST BE SPECIFIED TO SUPRESS HESSIAN UPDATE  *
C           AND READING IN FROM INP2, WHICH ARE NOW DEFAULT.           *
C           MURT MUST BE USED TO NOT HAVE POWELL OR BFGS UPDATES       *
C***********************************************************************
C      BY:  TPH                                                        *
C    DATE:  APRIL 23, 1988                                             *
C  REASON:  CHANGE MACHB TO USE INTERNAL COORDINATES SPECIFICATION     *
C           THE SAME AS IN INTDER. (BEND,OUT,LIN1,LIN2)                *
C***********************************************************************
C      BY:  TPH                                                        *
C    DATE:  FEBRUARY 26, 1988                                          *
C  REASON:  PUT PARTS BRIAN YATES CHANGED IN HIS VERSION,              *
C           KEEPING JON BAKER'S OPTEFC NEARLY THE SAME AS IN           *
C           GAUSSIAN82, AND USING AN INTERFACE CALLED EFC.             *
C***********************************************************************
C      BY:  TPH                                                        *
C    DATE:  FEBRUARY 8, 1988                                           *
C  REASON:  PUT COMMENTS IN PROGRAM FIFTH, RENUMBER                    *
C           INDICATING WHERE CHANGES IN THE NEAR FUTURE WILL BE MADE   *
C***********************************************************************
C      BY:  TPH                                                        *
C    DATE:  OCTOBER 26, 1987                                           *
C  REASON:  PUT IN EIGENVECTOR FOLLOWING ROUTINE OF JON BAKER,         *
C           POWELL UPDATE OF THE HESSIAN, AND DIIS FOR GEOMETRY OPT.   *
C           BESSLER PREVIOUSLY PUT IN BFGS AND DAVIDON-FLETCHER-POWELL *
C           UPDATES.  FILE INP4 HAS INTERNAL FORCES PREVIOUS STEPS.    *
C           FILE INP3 HAS THE EIGENVECTOR FROM PREVIOUS STEP IN OPTEFC.*
C***********************************************************************
C***********************************************************************
C   FUTURE NEEDS FOR BMAT: TAKE CARE OF DISCONTINUITY IN ANGULAR       *
C           COORDINATES (FOR NOW CHANGING ANGLE DEFINITIONS MAY HELP), *
C           FIX GDIIS,   PRINT RMS GRADIENT,                           *
C           AUTOMATIC ROTATION SO THAT ONE COORDINATE IS AT ORIGIN     *
C           OR SYMMETRY PLANE PRESERVED ETC.,                          *
C           CALCULATION OF # OF TOTALLY SYMMETRIC MODES FOR WARNINGS   *
C           STEP LENGTH CONTROL FOR REGULAR BMAT                       *
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include <error.h>
      integer errcod, frdi, fcnt, fwkadd, fwkclr
      integer nstre,nbnd,ntor,nlin1,nout,ninv,nlin2
      integer coallc,isizof
      real*8 ra(maxcor)
      integer ia(1)
      character junk
      CHARACTER AI*80,WKEY1*4,WKEY2*4
      INTEGER INP,IOUT,INP2,IPUN,IGMUP
      COMMON /IOP/ IOP(50)
      COMMON /GANZ/ NA,LOPT(27)
      parameter(inp=30,iout=6,inp2=32,ipun=33,igmup=8)
cets050591      PARAMETER(INP=30,IOUT=31,INP2=32,IPUN=33,IGMUP=5)
      PARAMETER(WKEY1='BMAT',WKEY2='STOP')
C....................................................................
C     THIS IS WRITTEN SO THAT SEVERAL BMAT SECTIONS IN THE BMAT     .
C     INPUT FILE CAN BE RUN, SHOULD THE NEED ARISE                  .
C....................................................................
cets050591      CALL NOUNFL
      call psinit('APPEND')
      call ffile(inp,'intco',0)
      call ffile(igmup,'geom',0)
      call ffile(11,' ',1)
      call ffile(inp2,'resul2',0)
      call ffile(ipun,'resul3',0)

      call tstart(iout)

      write(iout,*) ' '
      write(iout,101)
      write(iout,*) ' '
 101  format(22x,'BMAT: The geometry updating program')

      call coinit(ia,maxcor*isizof('REAL'))
      
      errcod = fwkclr()
      errcod = fwkadd(5,'BMAT')
      errcod = fwkadd(5,'DEFAULT')

      ndum=0
      errcod = fcnt('DUMMY_ATOMS',ndum)

      ndisp=0
      errcod = fcnt('DISPLACEMENTS',ndisp)

      nfix=0
      errcod = fcnt('FIX_COORD',nfix)

c read in internal coordinates

      errcod = fwkclr()
      errcod = fwkadd(5,'DUMMY')
      errcod = fwkadd(inp,'DUMMY')

      nstre=0
      nbnd=0
      ntor=0
      nout=0
      ninv=0
      nlin1=0
      nlin2=0

      errcod = fcnt(':INTCO:STRE',nstre)
      errcod = fcnt(':INTCO:BEND',nbnd)
      errcod = fcnt(':INTCO:TORS',ntor)
      errcod = fcnt(':INTCO:OUT',nout)
      errcod = fcnt(':INTCO:LIN1',nlin1)
      errcod = fcnt(':INTCO:LIN2',nlin2)
      errcod = fcnt(':INTCO:INVR',ninv)
      
      nsim=nstre+nbnd+ntor+nout+nlin1+nlin2+ninv

      istype = coallc(ia,nsim*isizof('INT'))
      isatom = coallc(ia,4*nsim*isizof('INT'))

      call rdsimp(ia(istype),ia(isatom),nsim,nstre,nbnd,
     1 ntor,nlin1,nlin2,ninv,nout)

      errcod = fcnt(':INTCO:SYMM',nsym)

c if doing displacements, read in ASYMM vector as well
      nasym=0
      if(ndisp.ne.0) then
        errcod = fcnt(':INTCO:ASYMM',nasym)
        nsym=nsym+nasym
      endif

      isymms = coallc(ia,101*nsym*isizof('INT'))
      isymmc = coallc(ia,100*nsym*isizof('REAL'))
      ilbls = coallc(ia,10*nsym*isizof('REAL'))

      call rdsymm(ia(isymms),ia(isymmc),nsym,nasym,ia(ilbls))
      
      errcod = fwkclr()
      errcod = fwkadd(8,'DUMMY')
      errcod = fcnt(':GEOMETRY',na)
      if(errcod.ne.EOK) then
        rewind(11)
        read(11,100) junk
        read(11,'(i5)') na
        if(na.eq.0) then
          write(6,*) 'cannot find GEOMETRY'
          call pabort
          call qabort
        endif
        call bottom(11)
      endif

      napdm = na+ndum

  100 format(a1)

      errcod = fwkclr()
      errcod = fwkadd(5,'BMAT')
      errcod = fwkadd(5,'DEFAULT')

      DO 20 I=1,50
         IOP(I)=0
   20 CONTINUE
      DO 30 I=1,27
         LOPT(I)=0
   30 CONTINUE
      LOPT(27)=5

      iccin = coallc(ia,4*ndisp*isizof('REAL'))
      iccc = coallc(ia,nsym*isizof('REAL'))
      isymb = coallc(ia,na*isizof('REAL'))
      iqqq = coallc(ia,nsym*isizof('REAL'))
      ifii = coallc(ia,nsym*isizof('REAL'))
      ibmat = coallc(ia,3*napdm*nsym*isizof('REAL'))
      iww = ibmat
      il = coallc(ia,nsym*isizof('INT'))
      im = coallc(ia,nsym*isizof('INT'))
      iia = coallc(ia,na*isizof('INT'))
      iifix = coallc(ia,nfix*isizof('INT'))
      iibohr = coallc(ia,nstre*isizof('INT'))
      ivmode = coallc(ia,nsym*isizof('REAL'))
      ic = coallc(ia,nsym*nsym*isizof('REAL'))
      icc = coallc(ia,nsym*isizof('REAL'))
      if = coallc(ia,3*napdm*isizof('REAL'))
      ifi = coallc(ia,nsym*isizof('REAL'))
      ixa = coallc(ia,napdm*isizof('REAL'))
      iya = coallc(ia,napdm*isizof('REAL'))
      iza = coallc(ia,napdm*isizof('REAL'))
      ixm = coallc(ia,3*napdm*isizof('REAL'))
      ixy = coallc(ia,3*napdm*isizof('REAL'))
      iqq = coallc(ia,nsym*isizof('REAL'))
      iqq1 = coallc(ia,nsym*isizof('REAL'))
      inv = coallc(ia,nsym*isizof('REAL'))
      ihe = coallc(ia,nsym*nsym*isizof('REAL'))
      ihn = coallc(ia,nsym*nsym*isizof('REAL'))
      in2v = coallc(ia,nsym*isizof('REAL'))
      idv = coallc(ia,nsym*isizof('REAL'))
      in3v = coallc(ia,nsym*isizof('REAL'))
      igrc = coallc(ia,nsym*isizof('REAL'))
      igec = coallc(ia,nsym*isizof('REAL'))
      iqd = coallc(ia,nsym*nsym*isizof('REAL'))
      ifd = coallc(ia,nsym*nsym*isizof('REAL'))
      iaa = coallc(ia,nsym*nsym*isizof('REAL'))
      iindi = coallc(ia,4*ndisp*isizof('INT'))
      iifst = coallc(ia,nsym*isizof('INT'))
      igm = coallc(ia,na*3*isizof('REAL'))
      igeobf = coallc(ia,20*2000*isizof('INT'))

      call fifth (inp,iout,inp2,ipun,igmup,nsim,nsym,nstre,nbnd,
     1 ntor,nlin1,nlin2,ninv,nout,ndisp,ndum,napdm,nfix,
     2 ia(istype),ia(isatom),ia(isymms),ia(isymmc),ia(iccin),ia(iccc),
     3 ia(isymb),ia(iqqq),ia(ifii),
     4 ia(ibmat),ia(iww),ia(il),ia(im),ia(iia),ia(iifix),ia(iibohr),
     5 ia(ivmode),ia(ic),ia(icc),ia(if),ia(ifi),ia(ixa),ia(iya),ia(iza),
     6 ia(ixm),ia(ixy),ia(iqq),ia(iqq1),ia(inv),ia(ihe),ia(ihn),
     7 ia(in2v),ia(idv),ia(in3v),ia(igrc),ia(igec),ia(iqd),ia(ifd),
     8 ia(iaa),ia(iindi),ia(iifst),ia(igm),ia(igeobf),ia,ia(ilbls))

      call tstop(iout)
      stop
      END
