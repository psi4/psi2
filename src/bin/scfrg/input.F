      SUBROUTINE INPUT (NC,NO,ALPHA,BETA,SU,ITYP,NST,SU1,SU2,NSO)
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
cets022791
#include <error.h>
C
      integer nst,nsym,mxcoef
      INTEGER NC(NST),NO(NST),ITYP(NSYM),IELAST(1),FCB10(16),NSO(8)
clj080889  The mips compiler cannot handle the following dimension declarator
clj080889      REAL*8 ALPHA(NSYM*(NSYM+1)/2),BETA(NSYM*(NSYM+1)/2),REPNUC
      real*8 alpha(*),beta(*),repnuc
      REAL*8 THRSH,EPS,FMT
      REAL*8 SU(MXCOEF),SU1(MXCOEF),SU2(MXCOEF)
      INTEGER ALABEL,BLABEL,ALBL2,BLBL2
      INTEGER PRINT,HALF,THIRD,OUTBUF
      INTEGER IORDER(300)
      CHARACTER*80 LINE
      LOGICAL TWOCON,IPOPLE
      DIMENSION MAP(21), IXTI(40),ALBL2(26),BLBL2(26)
C
cets022791
      character*8 opentp
      character*80 alabl2
      integer frdc, frdboo, frdi, frdf, frdfa1, frdia1, frdca1
      integer fcnt
      integer errcod
      logical oldinp
      common /oldin/ oldinp
      common /damp/ dam1,damp
cets030291
      integer infl30
      character*6 dertyp,cityp
clj032691
      integer chrend
c
      COMMON /T102/ NOUT,HALF,THIRD,MAXNO,OUTBUF,KINADR
      COMMON /POP/ IPOPLE
      COMMON /GVB/ ALPHA1,ALPHA2,ALPHA3,IJUP,TWOCON
      COMMON /TAPES/ NIN,IR,IOUT,ICHECK,IFLAG3,ITAPE,I30
      COMMON /PRINTQ/ PRINT
      COMMON /CORE/ NCOR
      COMMON /IFDMP/ DAMPSV,INFLG,MICMX,IDFIX
      COMMON /SCFIN/ REPNUC,IOPEN,ISCF,IPRCT,IVECT,NBFLD,ICUT
      COMMON /FLAGS/ IDEBUG,ISTO,IPUNCH,ISAVE,IRORDR
      COMMON /INPT/  NSYM,IPMAX,NMAT,ITMAX,NSFMAX,NTOTAL,MXCOEF,MSFRU
      COMMON /NEWRUN/ ITPREP,ITCON
      COMMON /ABOPEN/ MM
      COMMON /HOLL/ BLABEL(20),ALABEL(26),FMT(10),JL
      COMMON /DIM21/ IPQ(256)
      COMMON /TAPE30/ I10(200),ITAP30
      COMMON /NEXTFG/ NXTFLG
CTJL
      COMMON /TJL/ ASHIFO,DIISCT,AERRO,AFAC
CTJL
C
clj This is a mess.  Characters as integers should be removed everywhere
      character*80 clabel
      equivalence (blabel(1),clabel)
      EQUIVALENCE (IELAST,ELAST)
C
      DATA NCOL/10/,A0 /0.0D0/
C
C        SET UP IPQ(I)=(I*(I-1))/2 ARRAY
      IPQ(1)=0
      DO 5 I=1,255
    5 IPQ(I+1)=IPQ(I)+I
C
C     ----- FIND INPUT RECORDS -----
C
      ipople=.false.
      twocon=.false.
      oldinp=.true.
      do 5001 i=1,nst
        nc(i)=0
        no(i)=0
 5001 continue
c
      CALL LOCATE(IR,'# SCF ####',IERROR)
      IF (IERROR.NE.0) THEN
cets022791         WRITE (IOUT,900)
cets022791  900    FORMAT (//,' # SCF #### ERROR: NO INPUT ###########',//)
cets022791         STOP
cets022791 add fancy-schmancy janssen-style input
        oldinp = .false.

        iopen=0
cets022791 these should always be zero
        ivect=0
        isave=0
        maxno=0
        micmx=0
        ncor=0

        errcod = frdc('LABEL',alabl2)
clj     write(iout,6502) blabel,alabl2
        write(6,*)
        write(iout,1101) clabel(1:chrend(clabel))
 1101   format(10X,A)
        if (clabel.ne.alabl2) then
          write(iout,1101) alabl2(1:chrend(alabl2))
          endif
        write(6,*)

        ijunk=0
        errcod = frdboo('POPLE',ijunk)
        if(ijunk.eq.1) ipople=.true.

        inflg=1
        errcod = frdboo('RESTART',inflg)

        ipunch=0
        errcod = frdboo('REORDER',ipunch)

        iprct=0
        errcod = frdi('EXTRAP',iprct)

        nxtflg=0
        errcod = frdi('DIISSTART',nxtflg)

        itras=40
        errcod = frdi('MAXITER',itras)

        print=0
        errcod = frdi('IPRINT',print)

        isto=20
        errcod = frdi('CUTOFF',isto)

        dam1 = 1.0d0
        errcod = frdf('LEVELSHIFT',dam1)

        damp = 0.0d0
        errcod = frdf('DAMPING',damp)

        dertyp = 'FIRST '
        cityp = 'SCF   '
        errcod = frdc('DERTYPE',dertyp)
        errcod = frdc('WFN',cityp)

cets figure out if openshell, and then get occupations and coupling 
c    coefs if necessary
c
        opentp = 'NONE    '
        errcod = frdc('OPENTYPE',opentp)
        if(opentp.ne.'NONE    ') then
          errcod = fcnt('SOCC',irtest)
          if(errcod.eq.EOK.and.irtest.ne.nst) then
           write(6,*) ' wrong number of irreps in SOCC'
           call pabort
           call qabort
          endif

          errcod = frdia1('SOCC',no,nst)
          if(errcod.ne.EOK) then
            write(6,*) ' no SOCC input found in input'
            write(6,*) ' OPENTYPE = ',opentp
            call pabort
            call qabort
          endif

          iopen=1
          mm=0
          do 5002 i=1,nst
              if(no(i).ne.0) mm=mm+1
 5002         continue
            jl=(mm*(mm+1))/2

          if(opentp.eq.'SINGLET ') then
            if(mm.eq.2) then
              alpha(1)=0.0d0
              alpha(2)=0.0d0
              alpha(3)=0.0d0
              beta(1)=1.0d0
              beta(2)=-3.0d0
              beta(3)=1.0d0
            else
              write(iout,*) ' '
              write(iout,*) 'this program cannot handle same symmetry'
              write(iout,*) 'open-shell singlets. try scfx'
              call qabort
              endif
          else if(opentp.eq.'TWOCON  ') then
            twocon = .true.
            if(mm.eq.2) then
              alpha(1)=0.0d0
              alpha(2)=0.0d0
              alpha(3)=0.0d0
              beta(1)=0.0d0
              beta(2)=1.0d0
              beta(3)=0.0d0
            else
              write(iout,*) ' '
              write(iout,*) 'this program cannot handle same symmetry'
              write(iout,*) ' tcscf. try scfx'
              call qabort
              endif
          else if(opentp.eq.'HIGHSPIN') then
            do 5020 i=1,jl
              alpha(i)=0.0d0
              beta(i)=1.0d0
 5020       continue
          else if(opentp.eq.'SPECIAL ') then
            errcod = frdfa1('ALPHA',alpha,jl)
            if(errcod.ne.EOK) then
              write(iout,*) ' '
              write(iout,*) 'no ALPHA coupling coefficients found'
              call pabort
              call qabort
            endif
            errcod = frdfa1('BETA',beta,jl)
            if(errcod.ne.EOK) then
              write(iout,*) ' '
              write(iout,*) 'no BETA coupling coefficients found'
              call pabort
              call qabort
            endif

            do 5003 i=1,jl
              beta(i) = alpha(i)-beta(i)
 5003       continue
          endif
        endif

        errcod = fcnt('DOCC',irtest)
        if(errcod.eq.EOK.and.irtest.ne.nst) then
         write(6,*) ' wrong number of irreps in DOCC'
         call pabort
         call qabort
        endif

        errcod = frdia1('DOCC',nc,nst)
        if(errcod.ne.EOK.and.iopen.ne.1) then
          write(iout,*) ' '
          write(iout,*) 'hey, add some electrons buddy!'
          write(iout,*) 'use DOCC and/or SOCC'
          call qabort
        endif

c figure out default convergence criterion
        iscf=10
        if(dertyp.eq.'NONE  '.and.cityp.eq.'SCF   ') then
          iscf=7
        else if(opentp.ne.'NONE  '.and.cityp.eq.'SCF   ') then
          iscf=8
        endif
        
        errcod = frdi('CONVERGENCE',iscf)

c now get new mo order if desired
        if(ipunch.ne.0) then
          ipunch=0
          norder=1
          errcod = frdia1('MOORDER',iorder,msfru)
          if(errcod.ne.EOK) then
            write(iout,*) ' '
            write(iout,*) ' no new mo ordering found'
            write(iout,*) ' calculation continuing'
            write(iout,*) ' '
            norder=0
          endif
        endif
c
        if(twocon) ipople=.false.
C
      else
C
   10   READ (IR,1000) ALABEL
        WRITE(IOUT,6500) BLABEL,ALABEL
        IF(IFLAG3.GT.0) WRITE (ICHECK,6500) BLABEL, ALABEL
        READ (IR,6415) IPRCT,ISCF,IOPEN,INFLG,IVECT,IPUNCH,PRINT,nxtflg
     1,    ISAVE,ITRAS,MAXNO,ISTO,MICMX,NCOR
      END IF
C
C   NXTFLG=ITERATION TO BEGIN USING DIIS
cets022791      NXTFLG=IDEBUG
      IDEBUG=0
C
      IF (IPUNCH .EQ.0) GOTO 201
      NORDER=1
      LOFF=0
      DO 202 L=1,NST
         IF (NSO(L) .LE. 0) GOTO 202
         ISTART=LOFF+1
         IEND=LOFF+NSO(L)
         READ (IR,9410) (IORDER(II) ,II=ISTART,IEND)
         WRITE(IOUT,9410) (IORDER(II) ,II=ISTART,IEND)
 202  LOFF=LOFF+NSO(L)
 201  CONTINUE
 9410 FORMAT(14I5)
      if(oldinp) TWOCON=.FALSE.
      ICUT=IPRCT
      IF (ICUT.EQ.0) ICUT=1
      IF (ICUT.GT.20) ICUT=1
      IF (IOPEN.EQ.2) THEN
         TWOCON=.TRUE.
         IOPEN=1
         IPRCT=0
      END IF
      IF (IPRCT .LT. 0) THEN
         IPRCT=0
         IPOPLE=.TRUE.
         IDFIX=1
      END IF
      IF(IPUNCH.EQ.0) IPUNCH=47
      IF(ISTO.LE.0) ISTO = 20
      THRSH=10.0D0**(-ISTO)
      IF (ISCF .LE. 0) ISCF=5
      EPS=10.0D0**(-ISCF)
      ITMAX = 40
      IF(ITRAS.GT.0) ITMAX = ITRAS
      IWFG=IOUT
      DO 4 JUNK=1,2
         WRITE(IWFG,6460)
 6460    FORMAT(5X,4HCNVG,4X,4HINPT,4X,4HPNCH,4X,4HDIIS,4X,4HITER,4X,
     *   4HMCMX,4X,4HDAMP,4X,4HITCN)
         WRITE(IWFG,6470)
 6470    FORMAT(5H EXTP,4X,4HOPEN,4X,4HVECT,4X,4HPRNT,4X,4HSAVE,4X,
     *   4HTHRS,4X,4HNCOR,4X,4HITPR)
         WRITE(IWFG,6410) IPRCT,ISCF,IOPEN,INFLG,IVECT,IPUNCH,PRINT,
     *   NXTFLG,    ISAVE,ITRAS,ISTO,MICMX,NCOR,IDFIX,ITPREP,ITCON
         IF (NCOR.NE.0) WRITE (IWFG,6480)
 6480    FORMAT (/,' CORE POTENTIALS USED ')
         WRITE(IWFG,6600) THRSH
 6600    FORMAT(/,49H THRESHOLD FOR ELIMINATION OF BASIS FUNCTIONS  = ,
     *   4X,1PE9.1)
         WRITE(IWFG,5100)EPS,IPRCT
         IF (IFLAG3.LE.0) GO TO 6
         IWFG=3
    4 CONTINUE
    6 CONTINUE
      IRORDR=0
      IF (IVECT .LT. 0) THEN
         IF (IAND(ABS(IVECT),4) .NE. 0) THEN
            WRITE (IOUT,6610)
            IF (IFLAG3.GT.0) WRITE (IFLAG3,6610)
 6610       FORMAT (/,' VECTOR WILL NOT BE REORDERED')
            IRORDR=1
            IVECT=IVECT+4
         END IF
      END IF
      IF(ISAVE.GT.0) GO TO 271
      IXCOUN=0
      IPMAX=0
      ICMAX=0
      IDMAX=0
      MXFBI=0
      NTOTAL=0
      WRITE(IOUT,6450) REPNUC
      IF(IFLAG3.GT.0) WRITE(ICHECK,6450) REPNUC
cets022791      MM =0
      III=1
C
C     ----- READ IN HOW MANY OF EACH SYMMETRY TO OCCUPY -----
C
      if(oldinp) then
        MM=0
        DO 106 L=1,NST
          READ (IR,'(2I5)') NC(L),NO(L)
          IF (NO(L).GT.0) MM=MM+1
  106   CONTINUE
        JL=MM*(MM+1)/2
        endif
      IF (INFLG.GE.0) THEN
C
C     ----- FETCH CONSTANTS FROM TAPE30 -----
C
         ITAP30=30
C        CALL RFILE(ITAP30)
         CALL WREADW(ITAP30,I10,200,101,JUNK)
         MPOINT=I10(2)
         MCONST=I10(3)
         MCALCS=I10(4)
         NCALCS=I10(5)
cets030291 add infl30 at position 45.  this is set to 1 in
c          outptx when the vector is written out, so
c          you can leave RESTART set to 'yes' all the time
c
         infl30=i10(45)
C
C     ----- IF INFLG IS 0, GET LAST VECTOR, UNLESS FIRST -----
C           CALCULATION, IN WHICH CASE, USE CORE HAMILTONIAN
C
cets030291         IF (INFLG.EQ.0.AND.NCALCS.EQ.1) GO TO 24
cets030291         IF (INFLG.EQ.0) INFLG=NCALCS-1
         if(inflg.eq.0.or.infl30.eq.0) go to 24
C
         NAT=I10(19)
         MXCOEF=I10(42)
         JL=ABS(I10(43))
C
C    ----- FETCH POINTERS TO CALCULATIONS -----
C
         JUNK=101+MCONST+MPOINT
         CALL WREADW(ITAP30,I10,MCALCS,JUNK,JUNK)
         LOCCAL=I10(INFLG)
C
C    ----- FETCH LABEL FROM OLD VECTOR -----
C
         DO 11 I=1,3
            CALL WREADW(ITAP30,BLABEL,20,LOCCAL,LOCCAL)
   11    CONTINUE
CRBR     WRITE (IOUT,901) INFLG,(BLABEL(IQ),IQ=1,20)
CR901    FORMAT (/,' USING VECTOR SET ',I3,' AS INITIAL GUESS '
CRBR #   ,            '-- LABEL IS:',/,(1X,20A4))
         WRITE (IOUT,901) INFLG
  901    FORMAT (/,' USING VECTOR SET ',I3,' AS INITIAL GUESS ',/)
C
C     ----- FETCH POINTERS TO VECTORS, ETC -----
C
         CALL WREADW(ITAP30,I10,20,LOCCAL,LOCCAL)
C
         LOCVEC=I10(1)
C
         LOCCAL=LOCCAL+INTOWP(3*NAT)
C
C     ----- FETCH OLD ENERGY FROM FILE -----
C
         JUNK=LOCCAL+INTOWP(1)
         CALL WREADW(ITAP30,IELAST,INTOWP(1),JUNK,JUNK)
         WRITE (IOUT,902) ELAST
         WRITE (ICHECK,902) ELAST
  902    FORMAT (' ENERGY FROM OLD VECTOR:',F14.8,/)
C
C     ----- FETCH VECTOR -----
C
         CALL WREADW(ITAP30,SU,INTOWP(MXCOEF),LOCVEC,LOCVEC)
         IF (NORDER .EQ. 0) GOTO 205
         LOFF=0
         LOFSQ=0
         DO 208 L=1,NST
            IF (NSO(L) .LE. 0) GOTO 208
            DO 209 J=1,NSO(L)
               ISTART=(J-1)*NSO(L)+LOFSQ
               JNEW=IORDER(J+LOFF)
               IEND=(JNEW-1)*NSO(L)+LOFSQ
               DO 209 I=1,NSO(L)
                  ISTART=ISTART+1
                  IEND=IEND+1
 209        SU2(IEND)=SU1(ISTART)
            LOFF=LOFF+NSO(L)
 208     LOFSQ=LOFSQ+NSO(L)*NSO(L)
         DO 210 L=1,MXCOEF
 210     SU1(L)=SU2(L)
         LOFF=0
         DO 211 L=1,NST
            WRITE(IOUT,333) L
            LSYM=NSO(L)
            CALL MATOUT(SU1(LOFF+1),LSYM,LSYM,LSYM,LSYM,IOUT)
 211     LOFF=LOFF+LSYM*LSYM
 333     FORMAT(/'  RE=ORDERED VECTOR FOR SYMMETRY TYPE ',I5/)
 205     CONTINUE
C
C
C
      END IF
C
C
C
      GO TO 273
C
C     ----- JUMP HERE IF WISH CORE HAMILTONIAN GUESS -----
C
   24 CONTINUE
      WRITE (IOUT,903)
  903 FORMAT (/,' FIRST RUN, SO DEFAULTING TO CORE-HAMILTONIAN GUESS',/)
      INFLG=-1
C
      GO TO 273
C
  271 PRINT 5800
  273 CONTINUE
c
cets022791
      if(.not.oldinp) then
        if(iopen.eq.1) then
          jl=(mm*(mm+1))/2
          write(iout,5600)
          mm1=1
          mm2=1
          do 5009 i=1,jl
            write (iout,5500) mm1,mm2,alpha(i),alpha(i)-beta(i)
            mm2=mm2+1
            if (mm2.gt.mm1) then
              mm1=mm1+1
              mm2=1
              end if
 5009       continue
          write(iout,*) ' '
          endif
        return
        endif
c
      IF (IOPEN.EQ.0) GO TO 1603
      IF(IFLAG3.GT.0) WRITE(ICHECK,5600)
      WRITE(IOUT,5600)
      JL=(MM*(MM+1))/2
      IF (TWOCON) THEN
         ALPHA(1)=A0
         ALPHA(2)=A0
         ALPHA(3)=A0
         BETA(1)=A0
         BETA(2)=-1.0D0
         BETA(3)=A0
         MM1=1
         MM2=1
         DO 9 I=1,JL
            IF (IFLAG3.GT.0) WRITE (ICHECK,5500) MM1,MM2,ALPHA(I)
     #      ,                                       BETA(I)
            WRITE (IOUT,5500) MM1,MM2,ALPHA(I),BETA(I)
            READ (IR,6420) DUM1,DUM2
            MM2=MM2+1
            IF (MM2.GT.MM1) THEN
               MM1=MM1+1
               MM2=1
            END IF
    9    CONTINUE
C
         BETA(2)=-BETA(2)
      ELSE
         MM1=1
         MM2=1
         DO 8 I=1,JL
            READ (IR,6420) ALPHA(I),BETA(I)
            IF (IFLAG3.GT.0) WRITE (ICHECK,5500) MM1,MM2,ALPHA(I)
     #      ,                                       BETA(I)
            WRITE (IOUT,5500) MM1,MM2,ALPHA(I),BETA(I)
            BETA(I)=ALPHA(I)-BETA(I)
            MM2=MM2+1
            IF (MM2.GT.MM1) THEN
               MM1=MM1+1
               MM2=1
            END IF
    8    CONTINUE
C
      END IF
 1603 CONTINUE
      IF(IR.NE.5) CALL SREW(IR)
      IR=5
      RETURN
C                        FORMAT STATEMENTS
 1000 FORMAT(26A3)
 1001 FORMAT(A3,1X)
 2000 FORMAT(/1H ,6X,6HIRREP ,A3,1X,19X,20HORBITAL COEFFICIENTS)
 2500 FORMAT(30H NUMBER OF ELECTRONS OCCUPYING/
     124H CLOSED SHELL ORBITALS =,F6.1,9F11.1)
 3500 FORMAT(22X,10F11.6)
 4200 FORMAT(24H MATRIX MAP - LAMBDA,P,Q)
 4300 FORMAT(7(5X,2(I2,1H,),I2,5X))
 5000 FORMAT(30H NUMBER OF ELECTRONS OCCUPYING/
     122H OPEN SHELL ORBITALS =,F8.1,9F11.1)
 5100 FORMAT(49H DENSITY MATRIX CONVERGENCE CRITERION          = ,4X,
     11PE9.1/49H ITERATIONS BEFORE EXTRAPOLATION               = ,I7)
 5200 FORMAT(21X,10(6X,2HC(,I2,1H)))
 5500 FORMAT(5X,I2,1H,,I2,14X,F12.8,5X,F12.8)
 5600 FORMAT(31H OPEN-SHELL ENERGY COEFFICIENTS/16H OPEN SHELL PAIR,12X,
     15HALPHA,12X,4HBETA)
 5800 FORMAT(31H NEW STATE - SAME CONFIGURATION)
 6410 FORMAT(20I4)
 6415 FORMAT(14I5)
 6420 FORMAT(2F20.10)
 6430 FORMAT(20F4.0)
 6450 FORMAT(1X,'NUCLEAR REPULSION ENERGY                      =',
     +F22.13)
 6500 FORMAT(' ',/,T10,20A4,/,T10,26A3,//)
 6501 FORMAT(/10X,26A3)
c6502 FORMAT(' ',/,T10,20A4,/,T10,A80,//)
      END
