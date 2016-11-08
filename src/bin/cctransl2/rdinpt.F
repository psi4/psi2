      subroutine rdinpt(label,nfzo,nfzv,nirred)
      IMPLICIT INTEGER (A-Z)
      character*(*) label
      INTEGER NFZO(8),NFZV(8)
c
#include <error.h>
      integer frdi, frdc, frdboo, errcod, cout, fcnt

      label = ' '

      CALL LOCATE (5,'# CCSD ###',IERR)
      IF(IERR.NE.0) THEN

         errcod = frdc('LABEL',LABEL)

         errcod = fcnt('FROZEN_DOCC',cout)
         if ((errcod.eq.EOK) .and. (cout.ne.NIRRED)) then
           write(6,*) 'ERROR: an incorrect number of FROZEN_DOCC irreps'
     &               // ' has been given'
           call pabort
           call qabort
           endif
         errcod = frdia1('FROZEN_DOCC',NFZO,NIRRED)
         if (errcod.eq.EEOD) then
           DO 622 I=1,NIRRED
              NFZO(I)=0
  622      CONTINUE
         else if (errcod.ne.EOK) then
           write(6,*) 'ERROR: problem reading ''frozen_docc'''
           call pabort
           call qabort
           endif

         errcod = fcnt('FROZEN_UOCC',cout)
         if ((errcod.eq.EOK) .and. (cout.ne.NIRRED)) then
           write(6,*) 'ERROR: an incorrect number of FROZEN_UOCC irreps'
     &               // ' has been given'
           call pabort
           call qabort
           endif
         errcod = frdia1('FROZEN_UOCC',NFZV,NIRRED)
         if (errcod.eq.EEOD) then
           DO 1622 I=1,NIRRED
              NFZV(I)=0
 1622      CONTINUE
         else if (errcod.ne.EOK) then
           write(6,*) 'ERROR: problem reading ''frozen_uocc'''
           call pabort
           call qabort
           endif

c Write input parameters out.
        write(6,*)
        write(6,6100) 'Parameters to be used:'
        write(6,6100) '  LABEL          = ',LABEL
        write(6,6101) '  FROZEN_DOCC    = ',(nfzo(i),i=1,nirred)
        write(6,6101) '  FROZEN_UOCC    = ',(nfzv(i),i=1,nirred)
6100    format(1X,A,A)
6101    format(1X,A,8(1X,I4))
6102    format(1X,A,I5)

      ELSE
         write(6,*)
         write(6,*) 'WARNING: using old style ''# CCSD ###'' input'
         READ(5,6014)LABEL
 6014    FORMAT(A)
         READ(5,6002)EPSI,NGO,NMIN ,MINDIM,MAXDIM,FLDIIS,IRSTR,NCSF
 6002    FORMAT(16I5)
         READ(5,6002)CONVI,MAXIT,LEVEL
         DO 623 I=1,NIRRED
            READ(5,6051)NFZO(I),NFZV(I)
 6051       FORMAT(I2,1X,I2)

  623    CONTINUE
c        READ(5,6004)OPTION,fl81
 6004    FORMAT(20A4)

        if (label.ne.' ') then
           WRITE(6,*)
           WRITE(6,*)
           WRITE(6,*)' >>>  CCSD INPUT LABEL  <<<'
           WRITE(6,*)
           WRITE(6,6024)LABEL
           endif
C
 6024   FORMAT(2X,A80)
        endif

      return
      end
