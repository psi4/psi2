      SUBROUTINE SCFDER(D1N,D1O,D1F,D1W,D1T,D1TOT,
     1                  D2N,D2O,D2F,D2W,D2T,D2TOT,DENS,DENT,WAO,DEN1,
     2                  IA,CC,MAXLFT,NPERSH)
      IMPLICIT REAL*8 (A-H,O-Z)
cets030691
#include <error.h>
      integer ip, prcntr
      CHARACTER*52 TITLE
      COMMON/BASIS/NBASIS,NBFAO,NBATRI,NBATR2,NTRI,NTRI2
      COMMON/CALCS/ISTYPE,ICTYPE,IRAMAN,IAFOCK,IDFOCK,IDINT
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOM,N3N,NATRI,NSYM
      COMMON/NDERS/NDERIV
      COMMON/PARA1/NSHELL,NT,MXT,NISO
      COMMON/TOLER/DLIM
      DIMENSION D1N(3,NATOM),D1O(3,NATOM),D1F(3,NATOM),D1W(3,NATOM)
      DIMENSION D1T(3,NATOM),D1TOT(3,NATOM)
      DIMENSION D2N(N3N,N3N),D2O(N3N,N3N),D2F(N3N,N3N),D2W(N3N,N3N)
      DIMENSION D2T(N3N,N3N),D2TOT(N3N,N3N)
      DIMENSION DENS(NBATRI,NTYPEP),DENT(NBATRI),WAO(NBATRI)
      DIMENSION DEN1(NBATRI)
      DIMENSION IA(1),CC(MAXLFT),NPERSH(5)
C
C   CALCULATE DERIVATIVES OF ONE-ELECTRON AND OVERLAP TERMS
C
      IC1=1
      IC2=IC1+NBATRI*N3N
      IC3=IC2+NBATRI*3
      IC4=IC3+NBATRI*NATRI
      ICMAX=IC4+NBATRI*6
      if (prcntr('IS_ON DEBUG').ne.0) write(6,*)'in scfder icmax=',icmax
C
C12-1-88  Branch for SCF or CI
C
      IF (ICTYPE.EQ.1) THEN
        ip=prcntr('IS_ON DEBUG')
        if(ip.eq.1) then
          write(6,*)'in scfder, dent matrix to be used in makone'
          write(6,*)' this is an scf calculation                '
          call print (dent,nbatri,nbasis,6)
        endif
C...................                                 DFH     DFS.....
        CALL MAKONE(DENT,WAO,D1O,D1F,D1W,D2O,D2F,D2W,CC(IC1),CC(IC2),
C...................DSH     DSS.....
     1              CC(IC3),CC(IC4))
      ELSE
        ip=prcntr('IS_ON DEBUG')
        if(ip.eq.1) then
          write(6,*)'in scfder, den1 matrix to be used in makone'
          write(6,*)' this is a correlated calculation          '
          call print (den1,nbatri,nbasis,6)
        endif
C...................                                 DFH     DFS.....
        CALL MAKONE(DEN1,WAO,D1O,D1F,D1W,D2O,D2F,D2W,CC(IC1),CC(IC2),
C...................DSH     DSS.....
     1              CC(IC3),CC(IC4))
      ENDIF
C
      if(ictype.eq.2) then
        ip=prcntr('IS_ON DEBUG')
        TITLE='HELLMAN-FEYNMAN CONTRIBUTION TO FIRST DERIVATIVES '
        if(ip.eq.0) CALL PRNONE(TITLE,D1F)
        TITLE='ONE-ELECTRON CONTRIBUTION TO FIRST DERIVATIVES    '
        if(ip.eq.0) CALL PRNONE(TITLE,D1O)
      endif
cges
c     TITLE='this is d1o                                       '
c     CALL PRNONE(TITLE,D1O)
c     TITLE='this is d1f                                       '
c     CALL PRNONE(TITLE,D1f)
cges
      CALL VADD(D1O,D1O,D1F,N3N)
      TITLE='ONE-ELECTRON CONTRIBUTION TO FIRST DERIVATIVES    '
cets030691
      ip = prcntr('IS_ON BRIEF')
      if (ip.eq.0) CALL PRNONE(TITLE,D1O)
C
      IF(NDERIV.GE.2) THEN
        CALL FOLD(D2F,N3N)
C       TITLE='HELLMAN-FEYNMAN CONTRIBUTION TO SECOND DERIVATIVES'
C       CALL PRNTWO(TITLE,D2F)
        CALL FOLD(D2O,N3N)
        CALL VADD(D2O,D2O,D2F,N3N*N3N)
        TITLE='ONE-ELECTRON CONTRIBUTION TO SECOND DERIVATIVES   '
cets030691
        ip = prcntr('IS_ON BRIEF')
        if(ip.eq.0) CALL PRNTWO(TITLE,D2O)
      END IF
CGES
      IF (ICTYPE.EQ.1.OR.ICTYPE.EQ.9) THEN
        TITLE='OVERLAP CONTRIBUTION TO FIRST DERIVATIVES         '
cets030691
        ip = prcntr('IS_ON BRIEF')
        if(ip.eq.0) CALL PRNONE(TITLE,D1W)
      ENDIF
C
      IF(NDERIV.GE.2) THEN
        CALL FOLD(D2W,N3N)
        TITLE='OVERLAP CONTRIBUTION TO SECOND DERIVATIVES        '
cets030691
        ip = prcntr('IS_ON BRIEF')
        if(ip.eq.0) CALL PRNTWO(TITLE,D2W)
      END IF
C
C   CALCULATE TWO-ELECTRON CONTRIBUTIONS
      IC1=1
      IC2=IC1+NBATRI*N3N
      IC3=IC2+NBATRI*3
      IC4=IC3+NBATRI*NATRI
      IC5=IC4+NBATRI*6
      IC6=IC5+NBATRI*N3N*NTREAD
      NST=NSHELL*(NSHELL+1)/2
CGES  IF (ICTYPE.EQ.2) THEN
      IF (ICTYPE.EQ.2.OR.ICTYPE.EQ.9) THEN
        IC7=IC6+NST
        IC8=IC7+NST
        IC9=IC8+NST
        IC10=IC9+NSHELL
CMachineSize (For  7Mb)  IFANDG= 10000  =  10**4
CMachineSize (For 15Mb)  IFANDG= 50625  =  15**4
cges    IFANDG= 10000
        IFANDG= 15*15*15*15
        IC11=IC10+IFANDG
        IC12=IC11+NBFAO
        ICMAX=IC12+IFANDG
      ELSE
        IC7=IC6
        IC8=IC7
        IC9=IC8
        IC10=IC9
        IC11=IC10
        IC12=IC11
        ICMAX=IC12
      ENDIF
ctph  CALL JK2DER(DENS,DENT,D1T,D2T,CC(IC5))
C.................                  B1MAT...IJOFF...KLOFF...KLSIZ...
      CALL JK2DER(DENS,DENT,D1T,D2T,CC(IC5),CC(IC6),CC(IC7),CC(IC8),
C.................SHIZ....RDMBF....IMOSHL...ABTEMP...
     1            CC(IC9),CC(IC10),CC(IC11),CC(IC12),NST)
C
C   SYMMETRIZE TWO-ELECTRON FIRST DERIVATIVES
      CALL SYMD(D1T)
      TITLE='TWO-ELECTRON CONTRIBUTION TO FIRST DERIVATIVES    '
cets030691
      ip = prcntr('IS_ON BRIEF')
      if(ip.eq.0) CALL PRNONE(TITLE,D1T)
C
C   SYMMETRIZE TWO-ELECTRON SECOND DERIVATIVES
      IF(NDERIV.GE.2) THEN
        CALL FOLD(D2T,N3N)
        CALL SYMH(D2T)
        TITLE='TWO-ELECTRON CONTRIBUTION TO SECOND DERIVATIVES   '
cets030691
        ip = prcntr('IS_ON BRIEF')
        if(ip.eq.0) CALL PRNTWO(TITLE,D2T)
      END IF
C
C   SUM UP ALL CONTRIBUTIONS FOR FIRST DERIVATIVES
C12-1-88   CHANGES MADE IN SUMMATION AND PRINTING
      DO 20 J=1,NATOM
        DO 10 I=1,3
          D1TOT(I,J) = D1O(I,J)
   10   CONTINUE
   20 CONTINUE
      CALL VADD(D1TOT,D1T,D1TOT,N3N)
      CALL VADD(D1TOT,D1N,D1TOT,N3N)
      TITLE='-----CI FIRST DERIVATIVES LESS CPHF TERMS-----    '
      IF (ICTYPE.EQ.1) THEN
        CALL VADD(D1TOT,D1TOT,D1W,N3N)
        TITLE='-----TOTAL SCF FIRST DERIVATIVES-----            '
      ENDIF
CGES
      IF (ICTYPE.EQ.9) THEN
        CALL VADD(D1TOT,D1TOT,D1W,N3N)
        TITLE='-----TOTAL CC  FIRST DERIVATIVES-----            '
      ENDIF
C
      CALL PRNONE(TITLE,D1TOT)
C
C   SUM UP ALL CONTRIBUTIONS FOR SECOND DERIVATIVES
      IF(NDERIV.GE.2) THEN
        CALL VADD(D2TOT,D2O,D2W,N3N*N3N)
        CALL VADD(D2TOT,D2T,D2TOT,N3N*N3N)
        CALL VADD(D2TOT,D2N,D2TOT,N3N*N3N)
        TITLE='-----TOTAL SCF SECOND DERIVATIVES (AO PART)-----  '
cets030691
        ip = prcntr('IS_ON BRIEF')
        if(ip.eq.0) CALL PRNTWO(TITLE,D2TOT)
      END IF
C
      IF(IDFOCK.GE.1) THEN
C   SYMMETRIZE B1 MATRIX ELEMENTS
C................B1MAT...
      CALL SYMF1(CC(IC5))
C
C   FORM TAPE42
c
cges don't need this crap for ccsd gradients
c
      IC1=1
      IC2=IC1+NBATRI*N3N
      IC3=IC2+NBATRI*3
      IC4=IC3+NBATRI*NATRI
      IC5=IC4+NBATRI*6
      IC6=IC5+NBATRI*N3N*NTREAD
      ICMAX=IC6+NBATRI
C.................            DFH     DFS     B1MAT   SS
      CALL DERFK1(D1TOT,D2TOT,CC(IC1),CC(IC2),CC(IC5),CC(IC6),NPERSH)
C
      END IF
C
C     IF(IDFOCK.GE.2) THEN
C   SYMMETRIZE B2 MATRIX ELEMENTS
C................B2MAT...
C     CALL SYMF2(CC(IC5))
C
C   FORM TAPE44
C     IC1=1
C     IC2=IC1+NBATRI*N3N
C     IC3=IC2+NBATRI*3
C     IC4=IC3+NBATRI*NATRI
C     IC5=IC4+NBATRI*6
C     IC6=IC5+NBATRI*N3N*NTREAD
C     ICMAX=IC6+NBATRI
C.................DSH     DSS     B2MAT    SS.....
C     CALL DERFK2(CC(IC3),CC(IC4),CC(IC6),CC(IC7))
C
C     END IF
C
C   SAVE FIRST DERIVATIVES ON TAPE11 AND TAPE30
      CALL WRDER(D1TOT)
C
      RETURN
      END
