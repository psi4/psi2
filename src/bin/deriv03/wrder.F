      SUBROUTINE WRDER(D1TOT)
      IMPLICIT REAL*8 (A-H,O-Z)
cets030691
#include <error.h>
      integer frdia1, fcnt, frdc, errcod
      integer numfzc, numfzv, numir
      character*58 label
      character*8 wfn, dertyp
      DIMENSION I30(200),A30(100)
      COMMON/CALCS/ISTYPE,ICTYPE,IRAMAN,IAFOCK,IDFOCK,IDINT
      COMMON/COORD/ZAN(50),X(3,50)
      COMMON/ENRGY/ENUC,ESCF
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOM,N3N,NATRI,NSYM
      COMMON/TOLER/DLIM
      DIMENSION D1TOT(3,NATOM)
      EQUIVALENCE (I30,A30)
cets    1 FORMAT(20A4)
    1 format(a43,3x,2(a8,2x),i3,' fzc',i3,' fzv')
    5 format(a58,2x,2(a8,2x))
    2 FORMAT(I5,F20.10)
    3 FORMAT(4F20.10)
    4 FORMAT(20X,3F20.10)
C
C::::::::::::::::::::::::::::::::::::::::::
C:::SAVE FIRST DERIVATIVES ON THE TAPE30:::
C::::::::::::::::::::::::::::::::::::::::::
C
C   GET CONSTANTS FROM TAPE30
      ITAP30=30
      CALL RFILE(ITAP30)
      CALL WREADW(ITAP30,I30,200,101,JUNK)
      IEND    = I30(1)
      MPOINT  = I30(2)
      MCONST  = I30(3)
      MCALCS  = I30(4)
      NCALCS  = I30(5)
C
C   AND GET LOCATION OF INDIVIDUAL CALCULATION TO SAVE GRADIENTS
      IPOS=101+MCONST+MPOINT+NCALCS-1
      CALL WREADW(ITAP30,LOCCAL,1,IPOS,JUNK)
      CALL WREADW(ITAP30,I30,60,LOCCAL,JUNK)
      LOCCAL=LOCCAL+61
      CALL WREADW(ITAP30,LOCDER,1,LOCCAL,JUNK)
      IF(LOCDER.LE.0) THEN
        LOCDER=IEND
        IEND=IEND+N3N*2
        CALL WWRITW(ITAP30,LOCDER,1,LOCCAL,JUNK)
        CALL WWRITW(ITAP30,IEND,1,101,JUNK)
      END IF
      CALL WWRITW(ITAP30,D1TOT,N3N*2,LOCDER,JUNK)
      CALL RCLOSE(ITAP30,3)
C
C::::::::::::::::::::::::::::::::::::::::::::::::::
C:::SAVE FIRST DERIVATIVES ON THE TAPE11 FOR SCF:::
C::::::::::::::::::::::::::::::::::::::::::::::::::
C12-1-88
CGES  IF(ICTYPE.EQ.1) THEN
      IF(ICTYPE.EQ.1.OR.ICTYPE.EQ.9) THEN
cets030691 add more junk to file11
        numir=0
        numfzc=0
        numfzv=0
        errcod = frdc('WFN',wfn)
        errcod = frdc('DERTYPE',dertyp)
        errcod = frdc('LABEL',label)

        if(ictype.eq.9) then
          errcod = fcnt('FROZEN_DOCC',numir)
          if(errcod.eq.EOK.and.numir.ne.0) then
            errcod = frdia1('FROZEN_DOCC',i30,numir)
            do 201 i=1,numir
  201         numfzc = numfzc+i30(i)
          endif

          errcod = fcnt('FROZEN_UOCC',numir)
          if(errcod.eq.EOK.and.numir.ne.0) then
            errcod = frdia1('FROZEN_UOCC',i30,numir)
            do 202 i=1,numir
  202         numfzv = numfzv+i30(i)
          endif
        endif

        ITAP11=11
        call ffile(itap11,' ',1)
        if(ictype.eq.9) then
          write(itap11,1) label,wfn,dertyp,numfzc,numfzv
        else
          write(itap11,5) label,wfn,dertyp
        endif
c
        WRITE(ITAP11,2) NATOM,ESCF
        DO 101 I=1,NATOM
          WRITE(ITAP11,3) ZAN(I),X(1,I),X(2,I),X(3,I)
  101   CONTINUE
        DO 102 I=1,NATOM
          WRITE(ITAP11,4) D1TOT(1,I),D1TOT(2,I),D1TOT(3,I)
  102   CONTINUE
        REWIND ITAP11
      ENDIF
C
      RETURN
      END
