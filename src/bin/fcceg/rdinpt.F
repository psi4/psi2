      subroutine rdinpt(epsi,ngo,nmin,mindim,maxdim,maxit,convi,level,
     &            ncsf,fldiis,option,fl81,label,rstr,nfzo,nfzv,
     &            nirred)
      IMPLICIT INTEGER (A-Z)
      character*(*) option,fl81,label(20),rstr
      INTEGER NFZO(8),NFZV(8)
      character*40 corrty
c
#include <error.h>
      integer frdi, frdc, frdboo, errcod, frdia1, cout, fcnt

C
C
C >>> READ INPUT FILE
C
      CALL LOCATE (5,'# CCSD ###',IERR)
      IF(IERR.NE.0) THEN
         errcod = frdc('LABEL',LABEL) 
         errcod = frdi('CCDIIS:EPSILON',EPSI)
         errcod = frdi('CCDIIS:NGO',NGO)
         errcod = frdi('CCDIIS:NMIN',NMIN)
         errcod = frdi('CCDIIS:MINDIM',MINDIM)
         errcod = frdi('CCDIIS:MAXDIM',MAXDIM)
         errcod = frdi('CCDIIS:FLDIIS',FLDIIS)
         errcod = frdboo('RESTART',IRSTR)
         if (errcod.ne.EOK) then
            RSTR='NO  '
         else 
            if (IRSTR .ne. 0) then
               RSTR = 'YES '
            else
               RSTR = 'NO  '
            endif
         endif
         errcod = frdi('NCSF',NCSF)
         IF(EPSI  .EQ.0)EPSI  =5
         IF(NGO   .EQ.0)NGO   =2
         IF(NMIN  .EQ.0)NMIN  =2
         IF(MINDIM.EQ.0)MINDIM=2
         IF(MAXDIM.EQ.0)MAXDIM=8
         IF(NCSF  .EQ.0)NCSF  =20
         errcod = frdi('CONVERGENCE',CONVI)
         errcod = frdi('MAXITER',MAXIT)
         errcod = frdi('LEVEL',LEVEL)
         IF(LEVEL.EQ.0)LEVEL  =5
         IF(CONVI .EQ.0)CONVI =7
         IF(MAXIT .EQ.0)MAXIT =20

         errcod = fcnt('FROZEN_DOCC',cout)
         if ((errcod.eq.EOK) .and. (cout.ne.NIRRED)) then
           write(6,*) 'ERROR: an incorrect number of FROZEN_DOCC irreps'
     &               // ' has been given'
           call pabort
           call qabort
         endif
         errcod = frdia1('FROZEN_DOCC',NFZO,NIRRED)
         if (errcod.eq.EEOD) then
           DO 100 I=1,NIRRED
              NFZO(I)=0
  100      CONTINUE
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
           DO 200 I=1,NIRRED
              NFZV(I)=0
 200       CONTINUE
         else if (errcod.ne.EOK) then
           write(6,*) 'ERROR: problem reading ''frozen_uocc'''
           call pabort
           call qabort
         endif
         errcod = frdc('FL81',fl81) 
         errcod = frdc('WFN',corrty)
c Convert WFN to option
         if (corrty.eq.'CCD') then
           option = corrty
         else if (corrty.eq.'CCSD') then
           option = corrty
         else if (corrty.eq.'CCSD_T') then
           option = 'CCT'
         else if (corrty.eq.'CCSDT1') then
           option = 'SDT1'
         else if (corrty.eq.'CCSDT2') then
           option = 'SDT2'
         else if (corrty.eq.'CCSDT3') then
           option = 'SDT3'
         else if (corrty.eq.'CCSDT4') then
           option = 'SDT4'
         else if (corrty.eq.'CCSDT') then
           option = 'FSDT'
         else if (corrty.eq.'CC_CISD') then
           option = 'CISD'
         else if (corrty.eq.'CC_BOTH') then
           option = 'BOTH'
         else
           write(6,*)
           write(6,*) 'ERROR: Illegal ''WFN'''
           call qabort
           endif
      ELSE
         write(6,*) '/READING OLD STYLE INPUT'
         READ(5,6004)LABEL
         WRITE(6,*)
         WRITE(6,*)
         WRITE(6,*)' >>>  CCSD INPUT LABEL  <<<'
         WRITE(6,*)
         WRITE(6,6024)LABEL
         READ(5,6002)EPSI,NGO,NMIN ,MINDIM,MAXDIM,FLDIIS,IRSTR,NCSF
         IF(EPSI  .EQ.0)EPSI  =5
         IF(NGO   .EQ.0)NGO   =2
         IF(NMIN  .EQ.0)NMIN  =2
         IF(MINDIM.EQ.0)MINDIM=2
         IF(MAXDIM.EQ.0)MAXDIM=8
         IF(NCSF  .EQ.0)NCSF  =20
         IF(IRSTR.EQ.1)THEN
            RSTR='YES '
         ELSE
            RSTR='NO '
         ENDIF
         READ(5,6002)CONVI,MAXIT,LEVEL
         IF(CONVI .EQ.0)CONVI =7
         IF(MAXIT .EQ.0)MAXIT =20
         DO 623 I=1,NIRRED
            READ(5,6051)NFZO(I),NFZV(I)
  623    CONTINUE
         READ(5,6004)OPTION,fl81
      ENDIF
      if (option.ne.'FSDT') then
        write(6,*) 'ERROR: Use cceg instead of fcceg.'
        call qabort
        endif
      if(level.eq.0.and.option.eq.'FSDT')then
         write(6,*)'you are using option=FSDT and level=0'
         write(6,*)'either set level=1-5 or option=CCSD'
         call qabort
      endif
C
 6002 FORMAT(16I5)
 6004 FORMAT(20A4)
 6024 FORMAT(2X,20A4)
 6051 FORMAT(I2,1X,I2)
C

      return
      end

