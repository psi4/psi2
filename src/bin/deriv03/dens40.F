      SUBROUTINE DENS40(DENS,DENT,WAO,DEN1,EIG,OCC,EAO,U,T,iopen)
      IMPLICIT REAL*8 (A-H,O-Z)
cets030691
#include <error.h>
      integer ip, prcntr
      integer errcod, frdc
      character*8 opentp
      COMMON/BASIS/NBASIS,NBFAO,NBATRI,NBATR2,NTRI,NTRI2
      COMMON/CALCS/ISTYPE,ICTYPE,IRAMAN,IAFOCK,IDFOCK,IDINT
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOM,N3N,NATRI,NSYM
      COMMON/SIGNS/IOFF(256),IPRNT
      COMMON/START/FOCC(10),NSORB(10),NSTR(10),NEND(10),MOTYP(256)
      COMMON/TOLER/DLIM
      DIMENSION DENS(NBATRI,NTYPEP),DENT(NBATRI),WAO(NBATRI)
      DIMENSION U(NBFAO,NBASIS),T(NBASIS,NBASIS)
      DIMENSION EIG(NBASIS),OCC(NBASIS),EAO(NBFAO,NBASIS),DEN1(NBATRI)
      DATA ZERO,TWO / 0.0D+00 , 2.0D+00 /
    1 FORMAT(//,2X,' SCF EIGENVECTOR (AO BASIS) MATRIX'/)
    2 FORMAT(//,2X,' SCF LAGRANGIAN (AO BASIS) MATRIX'/)
    3 FORMAT(//,2X,' TOTAL DENSITY MATRIX'/)
    4 FORMAT(//,2X,' SHELL DENSITY MATRIX, ITYP = ',I5/)
    5 FORMAT(//,2X,' AVERAGE DENSITY MATRIX'/)
    6 FORMAT(//,2X,' CI ONE PDM'/)
C
C   READ IN AO-MO EIGENVECTORS
C
C12-1-88  EIGENVECTOR SHOULD ALREADY BE ROTATED BY TFOCK SUBROUTINE
C
      CALL MREAD(EIG,17)
      CALL MREAD(OCC,18)
      CALL MREAD(EAO,20)
cets030691      IF(IPRNT.LE.2) GO TO 201
      ip = prcntr('IS_ON DEBUG')
      if(ip.eq.1) then
        WRITE(6,1)
        CALL EIGOUT(EAO,EIG,OCC,NBFAO,NBASIS,NBFAO,NBASIS,6)
      endif
C
C   READ IN THE LAGRANGIAN MATRIX
C
C12-1-88  READ CI LAGRANGIAN AND ONEPDM FROM YUKIO'S MASTER FILE 40
C
  201 CONTINUE
      IF(ICTYPE.EQ.1) THEN
        CALL MREAD(WAO,24)
      ELSE
        CALL MREAD(WAO,32)
        CALL MREAD(DEN1,34)
      ENDIF
        DO 101 I=1,NBATRI
  101     WAO(I)=WAO(I)*TWO
cets030691        IF(IPRNT.LE.2) GO TO 202
        ip = prcntr('IS_ON DEBUG')
        if(ip.eq.1) then
          WRITE(6,2)
          CALL PRINT(WAO,NBATRI,NBFAO,6)
          WRITE(6,6)
          CALL PRINT(DEN1,NBATRI,NBFAO,6)
        endif
C
C   FORM THE TOTAL DENSITY MATRIX
  202 CONTINUE
        IJ=0
        DO 103 I=1,NBFAO
          DO 103 J=1,I
            IJ=IJ+1
            VALU=ZERO
            DO 102 K=1,NBASIS
              IF(OCC(K).LE.ZERO) GO TO 102
              VALU=VALU+EAO(I,K)*EAO(J,K)*OCC(K)
  102       CONTINUE
            DENT(IJ)=VALU
  103   CONTINUE
cets030691        IF(IPRNT.LE.3) GO TO 203
        ip = prcntr('IS_ON DEBUG')
        if(ip.eq.1) then
          WRITE(6,3)
          CALL PRINT(DENT,NBATRI,NBFAO,6)
        endif
C
C   FORM THE SHELL DENSITY MATRICES
  203 CONTINUE
      DO 106 ITYP=1,NTYPES
        IJ=0
        COEF=TWO
        IF(ISTYPE.EQ.4.AND.ITYP.GT.1) COEF=COEF*FOCC(ITYP)
        IST=NSTR(ITYP)
        IED=NEND(ITYP)
        DO 105 I=1,NBFAO
        DO 105 J=1,I
          IJ=IJ+1
          VALU=ZERO
          DO 104 K=IST,IED
            VALU=VALU+EAO(I,K)*EAO(J,K)
  104     CONTINUE
          DENS(IJ,ITYP)=VALU*COEF
  105   CONTINUE
cets030691      IF(IPRNT.LE.3) GO TO 106
        ip = prcntr('IS_ON DEBUG')
        if(ip.eq.1) then
          WRITE(6,4) ITYP
          CALL PRINT(DENS(1,ITYP),NBATRI,NBFAO,6)
        endif
  106 CONTINUE
C
cets030691 not sure about this, but it makes triplet methylene work
      errcod = frdc('OPENTYPE',opentp)
      if(ictype.eq.2.and.iopen.ne.0) then
        if(opentp.ne.'TWOCON  ') then
          call rfile(49)
          call sread(49,dens,4)
          do 300 i=1,nbatri
  300       dent(i)=0.0d0
          do 113 i=1,ntypes
            call sread(49,dens(1,i),intowp(nbatri))
            do 112 j=1,nbatri
              dens(j,i)=dens(j,i)*2.0d0
              dent(j)=dent(j)+focc(i)*dens(j,i)*0.25d0
  112       continue
  113     continue
          call rclose(49,3)
        else
c         write(6,*) ' hmm, tcscf-ci huh.  Well I sure cant do that'
c         call qabort
          do 301 i=1,nbatri
            dent(i)=dent(i)*1.0d0
  301     continue
        endif
c
        ip = prcntr('IS_ON DEBUG')
        if(ip.eq.1) then
          WRITE(6,*) 'dent after redoing'
          CALL PRINT(DENT,NBATRI,NBFAO,6)
        endif
      endif
c
C   FORM THE DENSITY MATRIX FOR THE AVERAGED FOCK
      IF(IAFOCK.EQ.0) GO TO 205
      IJ=0
      DO 108 I=1,NBFAO
      DO 108 J=1,I
      IJ=IJ+1
      VALU=ZERO
      DO 107 K=1,NBASIS
      KTYP=MOTYP(K)
      FAC=FOCC(KTYP)
      VALU=VALU+EAO(I,K)*EAO(J,K)*FAC
  107 CONTINUE
      DENS(IJ,NTYPEP)=VALU
  108 CONTINUE
cets030691      IF(IPRNT.LE.3) GO TO 205
      ip = prcntr('IS_ON DEBUG')
      if(ip.eq.1) then
        WRITE(6,5)
        CALL PRINT(DENS(1,NTYPEP),NBATRI,NBFAO,6)
ctph    CALL PRINT(DENS(1,NTYPEP),NBATRI,NBASIS,6)
      endif
C
  205 CONTINUE
      RETURN
      END
