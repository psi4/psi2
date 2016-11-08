C     ////////////////////////////////////////////////////////
      SUBROUTINE INPUT (NC,NO,ALPHA,BETA,SU,ITYP,NST,SU1,SU2,NSO
     1             ,QALP,RALP)
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
Crawdad.12.21.92
#include <error.h>
C
      INTEGER QALP,RALP
      integer nst,nsym,mxcoef
      INTEGER NC(NST),NO(NST),ITYP(NSYM),IELAST(1),FCB10(16),NSO(8)
cets090589 mips compiler can't handle the following
cets090589      REAL*8 ALPHA(NSYM*(NSYM+1)/2),BETA(NSYM*(NSYM+1)/2),REPNUC
      real*8 alpha(3),beta(4),repnuc
      REAL*8 THRSH,EPS,FMT
      REAL*8 SU(MXCOEF),SU1(MXCOEF),SU2(MXCOEF)
      INTEGER ALABEL,BLABEL,ALBL2,BLBL2
      INTEGER PRINT,HALF,THIRD,OUTBUF
      INTEGER IORDER(300)
      CHARACTER*80 LINE
      LOGICAL TWOCON,IPOPLE,CHOSE,GSTATE,TERMIN,QRDAMP,MOFIX
      DIMENSION MAP(21), IXTI(40),ALBL2(26),BLBL2(26)
C
Crawdad.12.21.92
      character*80 alabl2,clabel,dlabel
      integer frdc,frdboo,frdi,frdia1
      integer errcod
      logical oldinp
      common/oldin/ oldinp
      common/eneg/ egmin,egmax,c2init
      integer chrend
      integer nocc(2),qalp2(4),ralp2(4)
C
      COMMON/CHSR/CHOSE,GSTATE,TERMIN,QRDAMP,MOFIX
      COMMON/QDATA/ISYMQ,NQ,NVIRQ,NDOC,NMATQ,IAVRQ
      COMMON /T102/ NOUT,HALF,THIRD,MAXNO,OUTBUF,KINADR
      COMMON /POP/ IPOPLE
      COMMON /GVB/ ALPHA1,ALPHA2,ALPHA3,IJUP,TWOCON
      COMMON /TAPES/ NIN,IR,IOUT,ICHECK,IFLAG3,ITAPE,ITERML
      COMMON /PRINT/ PRINT
      COMMON /CORE/ NCOR
      COMMON /IFDMP/ DAMPSV,INFLG,MICMX,IDFIX
      COMMON /SCFIN/ REPNUC,IOPEN,ISCF,IPRCT,IVECT,NBFLD,ICUT
      COMMON /FLAGS/ IDEBUG,ISTO,IPUNCH,ISAVE,IRORDR
      COMMON /INPT/  NSYM,IPMAX,NMAT,ITMAX,NSFMAX,NTOTAL,MXCOEF,MSFRU
      COMMON /NEWRUN/ ITPREP
      COMMON /ABOPEN/ MM
      COMMON /HOLL/ BLABEL(20),ALABEL(26),FMT(10),JL
      COMMON /DIM21/ IPQ(256)
      COMMON /TAPE30/ I10(200),ITAP30
C
      EQUIVALENCE (IELAST,ELAST)
C
Crawdad.12.21.92
      equivalence (blabel(1),clabel)
      equivalence (alabel(1),dlabel)
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
      oldinp = .true.
      do 5001 i=1,nst
        nc(i) = 0
        no(i) = 0
 5001 continue
c
      CALL LOCATE(IR,'# SCFEX ##',IERROR)
      IF (IERROR.NE.0) THEN
Crawdad.12.21.92          WRITE (IOUT,900)
Crawdad.12.21.92   900    FORMAT (//,' # SCFEX ## ERROR: NO INPUT ###########',//)
Crawdad.12.21.92          STOP

	oldinp = .false.
Crawdad.12.21.92   The following are set to zero always.
	ivect = 0
	iterm = 0
	maxno = 0
	micmx = 0
	ncor = 0
      
	errcod = frdc('LABEL',dlabel)
	write(iout,1101) clabel(1:chrend(clabel))
	write(iout,1101) dlabel(1:chrend(dlabel))
1101	format(10X,A)
	
	iavrq = 10
	errcod = frdi('DAMPING',iavrq)
      
	iscf = 5
	errcod = frdi('CONVERGENCE',iscf)

	ijunk = 1
	errcod = frdboo('EXCITED',ijunk)
	if(ijunk.eq.0) then
		iopen = 0
	else
		iopen = 1
        endif

	inflg = -1
	errcod = frdi('GUESS', inflg)

	ipunch = 0
	errcod = frdi('REORDER',ipunch)

	print = 0
	errcod = frdi('PRINT',print)

	ichsr = 0
	errcod = frdi('CHOOSE',ichsr)

	itras = 100
	errcod = frdi('MAXITER',itras)

	isto = 20
	errcod = frdi('CUTOFF',isto)

	if(ipunch.ne.0) then
		ipunch = 0
		norder = 1
		errcod = frdia1('MOORDER',iorder,msfru)
		if(errcod.ne.EOK) then
		  write(iout,*) ' '
		  write(iout,*) 'No New MO Ordering Found'
		  write(iout,*) 'Continuing Anyway...'
		  write(iout,*) ' '
		  norder = 0
		endif
	endif

        mm = 0

	errcod = frdia1('DOCC',nc,nst)
	errcod = frdia1('SOCC',no,nst)
	errcod = frdia1('SOCC_ORB',nocc,2)

	do j=1,nst
		if(no(j).ne.0) mm=mm+1
	enddo

	qalp = nocc(1)
	ralp = nocc(2)

	egmin = 0.0
	egmax = 0.0
	c2init = 0.0
	errcod = frdf('EGMIN',egmin)
	errcod = frdf('EGMAX',egmax)
	errcod = frdf('C2INIT',c2init)

      else
C
C
   10 READ (IR,1000) ALABEL
      WRITE(IOUT,6500) BLABEL,ALABEL
      IF(IFLAG3.GT.0) WRITE (ICHECK,6500) BLABEL, ALABEL
      READ (IR,6415) IAVRQ,ISCF,IOPEN,INFLG,IVECT,IPUNCH,PRINT,ICHSR
     1,    ITERM,ITRAS,MAXNO,ISTO,MICMX,NCOR

      end if

CWA   DIFFERENT FOR EXCITED STATES.  INITIALIZE IT48 OR PREPARE TO
CWA   REORDER OPEN SHELLS.
CWA   ITCON=IDFIX=ITPREP=IDEBUG=ISAVE=0
      IPRCT=0
      IT48=48
      QRDAMP=.FALSE.
      MOFIX=.FALSE.
      IF(IAVRQ.NE.0) QRDAMP=.TRUE.
CWA   ALLOW DUMMY RUN WITH CONVERGED MO'S TO SET UP TFOCK.
      IF(ITRAS.LT.0) THEN
         ITRAS=-ITRAS
         MOFIX=.TRUE.
      END IF
      GSTATE=.FALSE.
      IF(IOPEN.EQ.0) THEN
          GSTATE=.TRUE.
          IOPEN=1
      END IF
      IF(INFLG.EQ.-1.OR.INFLG.EQ.-2) THEN
          IN48=IABS(INFLG)
          CALL SWRIT(IT48,IN48,1)
          CALL SREW(IT48)
          INFLG=1
      END IF
      IF(IPUNCH.EQ.-1) THEN
          IN48=-5
          CALL SWRIT(IT48,IN48,1)
          CALL SREW(IT48)
          IPUNCH=0
      END IF
      CHOSE=.FALSE.
      IF(ICHSR.LT.0) CHOSE=.TRUE.
      TERMIN=.FALSE.
      IF(ITERM.EQ.1) TERMIN=.TRUE.
C
      IF (IPUNCH .EQ.0) GOTO 201
      IF(oldinp) then
      NORDER=1
      LOFF=0
      DO 202 L=1,NST
         IF (NSO(L) .LE. 0) GOTO 202
         ISTART=LOFF+1
         IEND=LOFF+NSO(L)
         READ (IR,9410) (IORDER(II) ,II=ISTART,IEND)
         WRITE(IOUT,9410) (IORDER(II) ,II=ISTART,IEND)
 202  LOFF=LOFF+NSO(L)
      else
	 goto 201
      endif
 201  CONTINUE
 9410 FORMAT(14I5)
      TWOCON=.FALSE.
      ICUT=IPRCT
      IF (ICUT.EQ.0) ICUT=1
      IF (ICUT.GT.20) ICUT=1
      IF (IOPEN.EQ.2) THEN
         TWOCON=.TRUE.
         IOPEN=1
         IPRCT=0
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
 6460    FORMAT(5X,4HCNVG,4X,4HINPT,4X,4HPNCH,4X,4HCHSR,4X,4HITER)
         WRITE(IWFG,6470)
 6470    FORMAT(5H AVRG,4X,4HOPEN,4X,4HVECT,4X,4HPRNT,4X,4HTERM,4X,
     *   4HTHRS)
         WRITE(IWFG,6410) IAVRQ,ISCF,IOPEN,INFLG,IVECT,IPUNCH,PRINT,
     *   ICHSR,ITERM,ITRAS
         IF (NCOR.NE.0) WRITE (IOUT,6480)
 6480    FORMAT (/,' CORE POTENTIALS USED ')
         WRITE(IWFG,6600) THRSH
 6600    FORMAT(/,49H THRESHOLD FOR ELIMINATION OF BASIS FUNCTIONS  = ,
     *   4X,1PE9.1)
         WRITE(IWFG,5100)EPS,IPRCT
         IF (IFLAG3.LE.0) GO TO 6
         IWFG=ICHECK
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
CWA   IF(ISAVE.GT.0) GO TO 271
      IXCOUN=0
      IPMAX=0
      ICMAX=0
      IDMAX=0
      MXFBI=0
      NTOTAL=0
      WRITE(IOUT,6450) REPNUC
      IF(IFLAG3.GT.0) WRITE(ICHECK,6450) REPNUC
Crawdad.122292  ARGH! I hate it when a variable is just _randomly_ set to zero!!
Crawdad.122292      MM =0
      III=1
C
C     ----- READ IN HOW MANY OF EACH SYMMETRY TO OCCUPY -----
C
CWA   CHANGES FOR EXCITED STATE SCF.
      if(oldinp) then
      MM=0
C
      DO 106 L=1,NST
         READ (IR,'(4I5)') NC(L),NO(L),IQ1,IR1
         IF(IQ1.NE.0) QALP=IQ1
         IF(IR1.NE.0) RALP=IR1
         IF (NO(L).GT.0) MM=MM+1
  106 CONTINUE
      endif
      JL=MM*(MM+1)/2
      IF(QALP.EQ.0.OR.RALP.EQ.0) THEN
         WRITE(ICHECK,801) QALP,RALP
 801  FORMAT(/'  MISSING OPEN SHELL SPECIFIER :  ',2I5)
      STOP
      END IF
      IF (INFLG.GE.0) THEN
C
C     ----- FETCH CONSTANTS FROM TAPE30 -----
C
         ITAP30=30
         CALL WREADW(ITAP30,I10,200,101,JUNK)
         MPOINT=I10(2)
         MCONST=I10(3)
         MCALCS=I10(4)
         NCALCS=I10(5)
C
C     ----- IF INFLG IS 0, GET LAST VECTOR, UNLESS FIRST -----
C           CALCULATION, IN WHICH CASE, USE CORE HAMILTONIAN
C
         IF (INFLG.EQ.0.AND.NCALCS.EQ.1) GO TO 24
         IF (INFLG.EQ.0) INFLG=NCALCS-1
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
C        WRITE (IOUT,901) INFLG,(BLABEL(IQ),IQ=1,20)
C 901    FORMAT (/,' USING VECTOR SET ',I3,' AS INITIAL GUESS '
C    #   ,            '-- LABEL IS:',/,(1X,20A4))
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
C 271 PRINT 5800
  273 CONTINUE
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
cets053193 still do the read so that the file pointer is in the right place
            if(oldinp) READ (IR,6420) DUM1,DUM2
Crawdad.122292  Setting Dummy Var.  to standards with no option to change..
Crawdad.122292  
            dum1 = 0.0
            dum2 = -1.0
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
cets053193 still do the read so that the file pointer is in the right place
            if(oldinp) READ (IR,6420) ALPHA(I),BETA(I)
Crawdad.122292 Setting Dummy Var.  to standards with no option to change..
Crawdad.122292  
            alpha(i) = 0.0
            beta(i) = -1.0
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
      END
