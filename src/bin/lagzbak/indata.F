C
C----------------------------------------------------------------------
C
      SUBROUTINE INDATA(OPTION,NFZO,NFZV,LABEL,NIRRED,IPRINT)
#include <error.h>
      character*40 corrty,dertyp
      character*80 labl
      integer cout,errcod
      integer frdc,fcnt,frdia1
      CHARACTER*4 OPTION,LABEL(20)
      INTEGER NFZO(NIRRED),NFZV(NIRRED)
      INPUT=5
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

      errcod = frdc('LABEL',labl)
      dertyp = 'FIRST '
      errcod = frdc('DERTYPE',dertyp)

      write(6,*)
      write(6,*) 'Parameters to be used:'
      write(6,*)    '  label        = ',labl
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
 6006 FORMAT(A4,1X,2I5)
 6007 FORMAT(//,2X,' WARNING! NO INPUT TO CCSD HAS BEEN FOUND.',/,
     .          2X,'          DEFAULT PARAMETERS WILL BE USED.')
 6024 FORMAT(2X,20A4)
 6051 FORMAT(I2,1X,I2)
      RETURN
      END
