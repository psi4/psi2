      SUBROUTINE PRN11(DER1,NAT,CCSD,CIDT)
C
      IMPLICIT REAL*8 (A-H,O-Z)
cets030691
#include <error.h>
      integer errcod, frdc, frdia1, fcnt
      integer numfzc, numfzv, numir
      character*58 label
      character*8 wfn, dertyp
      integer i30(20)
C
      REAL*8 DER1(3,NAT),A30(100),CIDT(3)
      INTEGER IA30(400),IPTR30(400)
      INTEGER IALBEL(20),IBLBEL(20),LBLDRT(26)
      CHARACTER*52 TITLE,TITLC
CCCC  CHARACTER*8  SCFTYP,CALTYP
      CHARACTER*1 LET(3)
      LOGICAL OPEN,CCSD
C
      DIMENSION ZAN(50),GEOM(3,50)
C
      COMMON /TAPES/ NOUT,ICHK,IT30,IT42,IT47,IT52,IT58,IT94,IT96,IT98
      COMMON /CALCJ/ ICALT,IDERT,ICIT
      COMMON /FILES/ IT54
C
      EQUIVALENCE (IA30,A30)
C
      DATA LET   /'X','Y','Z'/
      DATA TITLE /'---------- CI GRADIENTS PLUS CPHF TERMS ----------'/
      DATA TITLC /'---------- CC GRADIENTS PLUS CPHF TERMS ----------'/
C
      ITAP11 = 11
      ITAP69 = 69
Ctdc101195  ffile() on file11 has already been called by mainin().  Don't 
Ctdc101195 re-call it or things die with the new IBM compilers.
Ctdc101195      call ffile(itap11,' ',1)
C
CC      IF (CALTYP.NE.'SCF  '.AND.CALTYP.NE.'CI   '.AND.
CC   #  CALTYP.NE.'MCSCF') THEN
CC        WRITE (NOUT,902) CALTYP,LEVEL,IPRINT
CC902     FORMAT (//,' UNKNOWN CALCULATION TYPE: ',A5,1X,A5,I5,//)
CC        STOP
CC      END IF
C
C     ----- GET CONSTANTS FROM TAPE30 -----
C
      CALL SREW(IT30)
      CALL WREADW(IT30,IA30,200,101,JUNK)
C
      IEND   =  IA30(1)
      MPOINT =  IA30(2)
      MCONST =  IA30(3)
      MCALCS =  IA30(4)
      NCALCS =  IA30(5)
      NAT    = IA30(19)
      NUM    = IA30(22)
      NNP    = IA30(23)
      NSHELL = IA30(27)
      NIRRED = IA30(28)
      NT     = IA30(29)
      NPRIM  = IA30(32)
      MXT    = IA30(36)
      NISO   = IA30(40)
C
C     ----- READ POINTERS FROM TAPE30 -----
C
      JUNK   = 101+MCONST
      CALL WREADW(IT30,IPTR30,MPOINT,JUNK,JUNK)
C
C     ----- NOW GO FOR THE GUSTO -----
C
      CALL WREADW(IT30,   ZAN,    INTOWP(NAT),IPTR30(1),JUNK)
C
C     ----- AND GET LOCATION OF INDIVIDUAL CALCULATION TO GET GEOM
C
      JUNK   = 101+MCONST+MPOINT+NCALCS-1
      CALL WREADW(IT30, LOCCAL, 1,JUNK,JUNK)
C
      CALL WREADW(IT30, IA30, 60,LOCCAL,LOCCAL)
C
      DO 500  ILBL=1,20
         IALBEL(ILBL) =  IA30(ILBL)
CCCC     CALL ITC(IA30(ILBL),0,ALABEL(ILBL))
  500 CONTINUE
C
      WRITE(3,*) ' ICALT =',ICALT,' IDERT = ',IDERT,' ICIT = ',ICIT
CCCC  IF(CITYP.NE.'CI    ' .AND. DERTYP.EQ.'FIRST ') THEN
CC    WRITE (ITAP11,899) (IA30(I),I=1,20)
C 899 FORMAT (20A4,/)
CC    WRITE (NOUT,900) (IA30(I),I=1,20)
C 900 FORMAT (/,(1X,20A4)/)
CCCC  END IF
C
      LOCCAL=LOCCAL+20
      CALL WREADW(IT30,  GEOM,  INTOWP(3*NAT),LOCCAL,LOCCAL)
C
C     ----- GET ENERGIES AND WRITE OUT -----
C
      CALL WREADW(IT30,   IA30,     INTOWP(10),LOCCAL,LOCCAL)
C
      ESCF = A30(2)
C
C    GET CI OR CC ENERGY
C
      IF(.NOT.CCSD)THEN
      CALL RFILE(IT54)
      CALL DCDLBL(IT54,LBLDRT,LBLDRT,IX,IX,IX,EGUESS,ECI,CZERO)
      ELSE
      CALL RFILE(ITAP69)
      CALL SREAD(ITAP69,NM,1)
      CALL SREAD(ITAP69,NDIMT2,1)
      CALL RGETSA(ITAP69,ISECT1)
      ISECT2=I2SEC(INTOWP(NM))+ISECT1
      ISECE=I2SEC(INTOWP(NDIMT2))+ISECT2
      CALL RSETSA(ITAP69,ISECE)
      CALL SREAD(ITAP69,EREF,INTOWP(1))
      CALL SREAD(ITAP69,XE,INTOWP(1))
      CALL RCLOSE(ITAP69,3)
      ECI=EREF+XE
      ENDIF
C
      WRITE (NOUT,901) A30(1)
  901 FORMAT (/,1X,' NUCLEAR REPULSION ENERGY =',F20.10)
      WRITE (NOUT,902) ESCF
  902 FORMAT (  1X,' SCF ENERGY               =',F20.10)
      WRITE (NOUT,903) ECI
  903 FORMAT( 1X,'  + CORRELATION ENERGY    =',F20.10,/)
C
C     -------  NOW WRITE CI GRADIENTS TO FILE11  --------
C     -------       AND DIPOLE MOMENT            --------
C
  111 FORMAT (20A4)
C 111 FORMAT (20A4,/) GET RID OF BLANK LINE AFTER TITLE. GES. JAN-7-88.
  112 FORMAT (I5,F20.10,F15.10,5X,F15.10,5X,F15.10)
C 112 FORMAT (I5,F20.10)
  113 FORMAT (4F20.10)
  114 FORMAT (20X,3F20.10)
C
cets030691 write a little more info to file11
      numir=0
      numfzc=0
      numfzv=0
      errcod = frdc('WFN',wfn)
      errcod = frdc('DERTYPE',dertyp)
      errcod = frdc('LABEL',label)

      errcod = fcnt('FROZEN_DOCC',numir)
      if(errcod.eq.EOK.and.numir.ne.0) then
        errcod = frdia1('FROZEN_DOCC',i30,numir)
        do 201 i=1,numir
  201     numfzc = numfzc+i30(i)
      endif

      errcod = fcnt('FROZEN_UOCC',numir)
      if(errcod.eq.EOK.and.numir.ne.0) then
        errcod = frdia1('FROZEN_UOCC',i30,numir)
        do 202 i=1,numir
  202     numfzv = numfzv+i30(i)
      endif

cets030691      WRITE (11,111) (IALBEL(I),I=1,20)
      write(11,115) label,wfn,dertyp,numfzc,numfzv
 115  format(a43,3x,2(a8,2x),i3,' fzc',i3,' fzv')
c
      WRITE (11,112) NAT,ECI,CIDT(1),CIDT(2),CIDT(3)
      DO 100  NG = 1,NAT
        WRITE (11,113)  ZAN(NG),GEOM(1,NG),GEOM(2,NG),GEOM(3,NG)
  100 CONTINUE
      DO 200  ND = 1,NAT
        WRITE (11,114)  DER1(1,ND),DER1(2,ND),DER1(3,ND)
  200 CONTINUE
C
C     -----  NOW WRITE CI GRADIENTS TO OUTPUT FILE  -----
C
      IF(.NOT.CCSD)THEN
      WRITE (NOUT,61) TITLE
      ELSE
      WRITE (NOUT,61) TITLC
      ENDIF
   61 FORMAT (//,T11,A50,//,' CENTER',T15,'X',T35,'Y',T55,'Z')
      DO 63  IC=1,NAT
        WRITE (NOUT,62) IC,DER1(1,IC),DER1(2,IC),DER1(3,IC)
   62   FORMAT (/,1X,I4,3F20.12)
   63 CONTINUE
C
C     -------  NOW WRITE CI GRADIENTS TO FILE30  --------
C
      JUNK=101+MCONST+MPOINT+NCALCS-1
      CALL WREADW(IT30,LOCCAL,1,JUNK,JUNK)
C
      LOCCAL=LOCCAL+61
      CALL WREADW(IT30,LOCDER,1,LOCCAL,JUNK)
      IF (LOCDER.LE.0) THEN
        LOCDER=IEND
        IEND=IEND+INTOWP(3*NAT)
        CALL WWRITW(IT30,LOCDER,1,LOCCAL,JUNK)
        CALL WWRITW(IT30,IEND,1,101,JUNK)
      END IF
C
      CALL WWRITW(IT30,DER1,INTOWP(3*NAT),LOCDER,JUNK)
C
      CALL RCLOSE(IT30,3)
C
C
      RETURN
      END
