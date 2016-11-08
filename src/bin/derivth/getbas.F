      SUBROUTINE GETBAS(NPERSH)
      IMPLICIT REAL*8 (A-H,O-Z)
cets030691
#include <error.h>
      integer ip, prcntr
ctph  DIMENSION I30(200),A30(100),CSPD(600),LISO(1600),LICT(400)
      DIMENSION I30(200),A30(100),CSPD(1000),LISO(1600),LICT(4000)
      DIMENSION NPERSH(5)
      COMMON/BASIS/NBASIS,NBFAO,NBATRI,NBATR2,NTRI,NTRI2
      COMMON/COORD/ZAN(50),X(3,50)
      COMMON/ENRGY/ENUC,ESCF
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOM,N3N,NATRI,NSYM
      COMMON/GAUSA/EX(200),CS(200)
      COMMON/GAUSB/KATOM(200),KTYPE(200),KSTART(200),KLOC(200),
     1             KMIN(200),KMAX(200),KNG(200)
CTPH  COMMON/GAUSC/T(36),PTR(3,3,8),DTR(6,6,8)
      COMMON/GAUSC/T(225),PTR(3,3,8),DTR(6,6,8),FTR(10,10,8),
     1             GTR(15,15,8)
      COMMON/GAUSD/INVT(8),ISO(200,8),LBLAT(2,50),ICT(50,8)
      COMMON/PARA1/NSHELL,NT,MXT,NISO
      COMMON/PARA2/NPRIM
      COMMON/SIGNS/IOFF(256),IPRNT
      COMMON/TOLER/DLIM
      integer chrend
      character*80 label
      EQUIVALENCE (I30,A30,label)
C
    1 FORMAT(/,3X,20A4/)
    2 FORMAT(2X,' NUCLEAR REPULSION ENERGY:',F20.10)
    3 FORMAT(2X,'               SCF ENERGY:',F20.10)
    4 FORMAT(/,2X,' NPRIM  = ',I5/)
    5 FORMAT(2X,I5,2F20.10)
    6 FORMAT(/,2X,' NSHELL = ',I5/)
    7 FORMAT(2X,8I5)
    8 FORMAT(/,2X,' NATOM = ',I5/)
    9 FORMAT(2X,I5,4F20.10)
C
czat
c variable ZERO was not initialized ! 22 may 1990
c
      ZERO=0.0d0
c
C   GET CONSTANTS FROM TAPE30
      ITAP30=30
      CALL SREW(ITAP30)
      CALL WREADW(ITAP30,I30,200,101,JUNK)
C
      IEND    = I30(1)
      MPOINT  = I30(2)
      MCONST  = I30(3)
      MCALCS  = I30(4)
      NCALCS  = I30(5)
      NBASIS  = I30(18)
      NATOM   = I30(19)
      NBFAO   = I30(22)
      NSHELL  = I30(27)
      NT      = I30(29)
      NPRIM   = I30(32)
      MXT     = I30(36)
      NISO    = I30(40)
      NSYM    = I30(41)
      if (prcntr('IS_ON DEBUG').ne.0) then
        write(6,*)' '
        write(6,*)'file30 parameters'
        write(6,*)' '
        write(6,*)' iend   =',iend
        write(6,*)' mpoint =',mpoint
        write(6,*)' mconst =',mconst
        write(6,*)' mcalcs =',mcalcs
        write(6,*)' ncalcs =',ncalcs
        write(6,*)' nbasis =',nbasis
        write(6,*)' natom  =',natom 
        write(6,*)' nbfao  =',nbfao 
        write(6,*)' nshell =',nshell
        write(6,*)' nt     =',nt    
        write(6,*)' nprim  =',nprim 
        write(6,*)' mxt    =',mxt   
        write(6,*)' niso   =',niso   
        write(6,*)' nsym   =',nsym
        write(6,*)' '
        endif
      NBATRI=IOFF(NBFAO+1)
      NTRI=IOFF(NBASIS+1)
      NBATR2=NBATRI*2
      NTRI2=NTRI*2
      N3N=NATOM*3
      NATRI=IOFF(N3N+1)
C
C   READ POINTERS FROM TAPE30
      IPOS=101+MCONST
      CALL WREADW(ITAP30,I30,MPOINT,IPOS,JUNK)
C
C   READ BASIS SET INFO ETC.
      CALL WREADW(ITAP30,ZAN,NATOM*2,I30(1),JUNK)
      CALL WREADW(ITAP30,LICT,NATOM*NT,I30(2),JUNK)
      CALL WREADW(ITAP30,EX,NPRIM*2,I30(5),JUNK)
ctph  CALL WREADW(ITAP30,CSPD,NPRIM*6,I30(6),JUNK)
      CALL WREADW(ITAP30,CSPD,NPRIM*10,I30(6),JUNK)
      CALL WREADW(ITAP30,KSTART,NSHELL,I30(7),JUNK)
      CALL WREADW(ITAP30,KATOM,NSHELL,I30(8),JUNK)
      CALL WREADW(ITAP30,KTYPE,NSHELL,I30(9),JUNK)
      CALL WREADW(ITAP30,KNG,NSHELL,I30(10),JUNK)
      CALL WREADW(ITAP30,KLOC,NSHELL,I30(11),JUNK)
      CALL WREADW(ITAP30,KMIN,NSHELL,I30(12),JUNK)
      CALL WREADW(ITAP30,KMAX,NSHELL,I30(13),JUNK)
      CALL WREADW(ITAP30,T,MXT*2,I30(25),JUNK)
      CALL WREADW(ITAP30,INVT,NT,I30(26),JUNK)
      CALL WREADW(ITAP30,LISO,NSHELL*NISO,I30(27),JUNK)
      CALL WREADW(ITAP30,LBLAT,2*NATOM,I30(28),JUNK)
      CALL WREADW(ITAP30,PTR,3*3*NT*2,I30(31),JUNK)
      CALL WREADW(ITAP30,DTR,6*6*NT*2,I30(32),JUNK)
CTPH
      CALL WREADW(ITAP30,FTR,10*10*NT*2,I30(33),JUNK)
      CALL WREADW(ITAP30,GTR,15*15*NT*2,I30(34),JUNK)
      CALL WREADW(ITAP30,NPERSH,5,I30(36),JUNK)
C   AND GET LOCATION OF INDIVIDUAL CALCULATION TO GET GEOMETRY
      IPOS=101+MCONST+MPOINT+NCALCS-1
      CALL WREADW(ITAP30,LOCCAL,1,IPOS,JUNK)
      CALL WREADW(ITAP30,I30,60,LOCCAL,JUNK)
clj   WRITE(6,1) (I30(I),I=1,20)
      write(6,'(/,3X,A,/)') label(1:chrend(label))
C
      IPOS=LOCCAL+80
      CALL WREADW(ITAP30,X,NATOM*6,IPOS,JUNK)
C
C   GET ENERGIES AND WRITE OUT
      IPOS=LOCCAL+80+NATOM*6
      CALL WREADW(ITAP30,I30,20,IPOS,JUNK)
      ENUC=A30(1)
      ESCF=A30(2)
      WRITE(6,2) ENUC
      WRITE(6,3) ESCF
C
C   COLLAPSE CONTRACTION COEFFICIENTS INTO ONE-D ARRAY
      IJ=0
ctph  DO 102 I=1,3
      DO 102 I=1,5
      DO 101 J=1,NPRIM
      IJ=IJ+1
      VALU=CSPD(IJ)
      IF(VALU.EQ.ZERO) GO TO 101
      CS(J)=VALU
CMachineSize  (Remove following for large core 15Mb in order to run "G's
c     IF(I.EQ.5) THEN
c       WRITE(6,*) 'ERROR: SMALL CORE VERSION WILL NOT RUN G FUNCTIONS'
c       WRITE(3,*) 'ERROR: SMALL CORE VERSION WILL NOT RUN G FUNCTIONS'
c       STOP
c     ENDIF
CMachineSize  (Remove above     for large core 15Mb in order to run "G's
  101 CONTINUE
  102 CONTINUE
C
      DO 104 I=1,NISO
      II=(I-1)*NSHELL
      DO 103 ISHELL=1,NSHELL
      II=II+1
      ISO(ISHELL,I)=LISO(II)
  103 CONTINUE
  104 CONTINUE
      DO 106 I=1,NT
      II=(I-1)*NATOM
      DO 105 IATOM=1,NATOM
      II=II+1
      ICT(IATOM,I)=LICT(II)
  105 CONTINUE
  106 CONTINUE
C
cets030691      IF(IPRNT.LE.3) GO TO 205
      ip = prcntr('IS_ON DEBUG')
      if(ip.eq.1) then
        WRITE(6,4) NPRIM
        DO 107 I=1,NPRIM
          WRITE(6,5) I,EX(I),CS(I)
  107   CONTINUE
        WRITE(6,6) NSHELL
        DO 108 I=1,NSHELL
          WRITE(6,7) I,KATOM(I),KTYPE(I),KSTART(I),KLOC(I),
     1           KMIN(I),KMAX(I),KNG(I)
  108   CONTINUE
        WRITE(6,8) NATOM
        DO 109 I=1,NATOM
  109   WRITE(6,9) I,ZAN(I),X(1,I),X(2,I),X(3,I)
      endif
C
  205 CONTINUE
      CALL RCLOSE(ITAP30,3)
      RETURN
      END
