      subroutine rdinpt(epsi,ngo,nmin,mindim,maxdim,maxit,convi,
     &            ncsf,fldiis,option,fl81,label,rstr,trstr,nfzo,nfzv,
     &            nirred)
      IMPLICIT INTEGER (A-Z)
      character*(*) option,fl81,label,rstr,trstr
      INTEGER NFZO(8),NFZV(8)
      character*40 corrty,dertyp
c
#include <error.h>
      integer frdi, frdc, frdboo, errcod, bRSTR, frdia1, cout, fcnt
      integer bTRSTR
      integer rddiis

c Initialize input parameters
      epsi = 0
      ngo = 0
      nmin = 0
      mindim = 0
      maxdim = 0
      maxit = 0
      convi = 0
      ncsf = 0
      fldiis = 0
      dertyp = 'NONE'
      corrty = ' '
      label = ' '
c     option = ' '
c     fl81 = ' '
c     rstr = ' '
c     trstr = ' '

      CALL LOCATE (5,'# CCSD ###',IERR)
      IF(IERR.NE.0) THEN

         errcod = frdc('LABEL',LABEL)
         errcod = frdi('NCSF',NCSF)
         errcod = frdi('MAXITER',MAXIT)
         errcod = frdi('CONVERGENCE',CONVI)
         errcod = frdc('DERTYPE',dertyp)
         errcod = frdc('WFN',corrty)
c        errcod = frdc('OPTION',OPTION)
c        errcod = frdc('FL81',FL81)
c        errcod = frdi('LEVEL',LEVEL)

c DIIS stuff
         rddiis = 0
         errcod = frdi('CCDIIS:EPSILON',EPSI)
         if (errcod.eq.0) rddiis = 1
         errcod = frdi('CCDIIS:NGO',NGO)
         if (errcod.eq.0) rddiis = 1
         errcod = frdi('CCDIIS:NMIN',NMIN)
         if (errcod.eq.0) rddiis = 1
         errcod = frdi('CCDIIS:MINDIM',MINDIM)
         if (errcod.eq.0) rddiis = 1
         errcod = frdi('CCDIIS:MAXDIM',MAXDIM)
         if (errcod.eq.0) rddiis = 1
         errcod = frdi('CCDIIS:FLDIIS',FLDIIS)
         if (errcod.eq.0) rddiis = 1

         errcod = frdboo('RESTART',bRSTR)
         if (errcod.ne.EOK) then
           RSTR='NO  '
         else
           if (bRSTR .ne. 0) then
             RSTR = 'YES '
           else
             RSTR = 'NO  '
             endif
           endif

         errcod = frdboo('T_RESTART',bTRSTR)
         if (errcod.ne.EOK) then
           TRSTR='NO  '
         else
           if (bTRSTR .ne. 0) then
             RSTR = 'YES '
             TRSTR = 'YES '
           else
             TRSTR = 'NO  '
             endif
           endif

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

c Convert DERTYPE to fl81
        if (dertyp.eq.'FIRST') then
          fl81 = 'GRAD'
        else
          fl81 = ' '
          endif

c Convert DERTYPE to CONVI, if CONVI is set to the default
        if (CONVI.eq.0) then
          if (dertyp.eq.'FIRST') then
            convi = 10
          else
            convi = 7
            endif
          endif

c Convert DERTYPE to MAXIT, if MAXIT is set to the default
        if (MAXIT.eq.0) then
          if (dertyp.eq.'FIRST') then
            MAXIT = 60
          else
            MAXIT = 50
            endif
          endif

c Set up defaults for the rest of the unspecified input parameters
        IF(EPSI  .EQ.0)EPSI  =5
        IF(NGO   .EQ.0)NGO   =2
        IF(NMIN  .EQ.0)NMIN  =2
        IF(MINDIM.EQ.0)MINDIM=2
        IF(MAXDIM.EQ.0)MAXDIM=8
        IF(NCSF  .EQ.0)NCSF  =20

c Write input parameters out.
        write(6,*)
        write(6,6100) 'Parameters to be used:'
        write(6,6100) '  LABEL          = ',LABEL
        write(6,6100) '  WFN            = ',corrty
        write(6,6100) '  DERTYPE        = ',dertyp
        write(6,6102) '  NCSF           = ',NCSF
        write(6,6102) '  MAXITER        = ',MAXIT
        write(6,6102) '  CONVERGENCE    = ',CONVI
        write(6,6102) '  RESTART        = ',bRSTR
        write(6,6102) '  T_RESTART      = ',bTRSTR
        write(6,6101) '  FROZEN_DOCC    = ',(nfzo(i),i=1,nirred)
        write(6,6101) '  FROZEN_UOCC    = ',(nfzv(i),i=1,nirred)
        if (rddiis.eq.1) then
          write(6,6102) '  CCDIIS:EPSILON = ',EPSI
          write(6,6102) '  CCDIIS:NGO     = ',NGO
          write(6,6102) '  CCDIIS:MINDIM  = ',MINDIM
          write(6,6102) '  CCDIIS:MAXDIM  = ',MAXDIM
          write(6,6102) '  CCDIIS:FLDIIS  = ',FLDIIS
          endif
6100    format(1X,A,A)
6101    format(1X,A,8(1X,I4))
6102    format(1X,A,I5)


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
         else if (corrty.eq.'FULL_CCSDT') then
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
         write(6,*)
         write(6,*) 'WARNING: using old style ''# CCSD ###'' input'
         READ(5,6014)LABEL
 6014    FORMAT(A)
         READ(5,6002)EPSI,NGO,NMIN ,MINDIM,MAXDIM,FLDIIS,IRSTR,NCSF
 6002    FORMAT(16I5)
         IF(IRSTR.EQ.1)THEN
            RSTR='YES '
            RSTR='NO  '
         else if(irstr.eq.3)then
            rstr='YES '
            trstr='YES '
         ELSE
            RSTR='NO '
            tRSTR='NO '
         ENDIF
         READ(5,6002)CONVI,MAXIT,LEVEL
         DO 623 I=1,NIRRED
            READ(5,6051)NFZO(I),NFZV(I)
 6051       FORMAT(I2,1X,I2)

  623    CONTINUE
         READ(5,6004)OPTION,fl81
 6004    FORMAT(20A4)

c Set up defaults for unspecified input parameters
        IF(EPSI  .EQ.0)EPSI  =5
        IF(NGO   .EQ.0)NGO   =2
        IF(NMIN  .EQ.0)NMIN  =2
        IF(MINDIM.EQ.0)MINDIM=2
        IF(MAXDIM.EQ.0)MAXDIM=8
        IF(NCSF  .EQ.0)NCSF  =20
        IF(CONVI .EQ.0)CONVI =7
        IF(MAXIT .EQ.0)MAXIT =50

        if (label.ne.' ') then
           WRITE(6,*)
           WRITE(6,*)
           WRITE(6,*)' >>>  CCSD INPUT LABEL  <<<'
           WRITE(6,*)
           WRITE(6,6024)LABEL
           endif
C
 6024   FORMAT(2X,A80)


C
C >>> WRITE OUT INPUT INFO
C
        WRITE(6,6003)EPSI,NGO,NMIN ,MINDIM,MAXDIM,FLDIIS,
     .                CONVI,MAXIT,LEVEL,OPTION,fl81,RSTR,NCSF
 6003   FORMAT(//,2X,'>>>  DIIS PARAMETERS  <<< ',/,
     .          /,2X,'EPSI   =',I5,
     .          /,2X,'NGO    =',I5,
     .          /,2X,'NMIN   =',I5,
     .          /,2X,'MINDIM =',I5,
     .          /,2X,'MAXDIM =',I5,
     .          /,2X,'FLDIIS =',I5,
     .        ///,2X,'>>>  CCSD PARAMETERS   <<<',/,
     .          /,2X,'CONVI  =',I5,
     .          /,2X,'MAXIT  =',I5,
     .          /,2X,'LEVEL  =',I5,
     .          /,2X,'OPTION =',1X,A4,
     .          /,2X,'fl81   =',1X,A4,
     .          /,2X,'RSTR   =',1X,A4,
     .          /,2X,'NCSF   =',I5)
        endif

      return
      end
