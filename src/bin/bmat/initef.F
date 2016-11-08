      SUBROUTINE INITEF(CORE)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 CORE(1)
      LOGICAL PRNT,NRFLAG,DUMP
      COMMON/IOP/IOP(50)
      COMMON/IO/IN,IOUT
      COMMON/GRDNT/ENERGY,F(50),FRCNST(1275),NVAR,IGETFC
      COMMON/OPTEF/
     $             ES,DMAX,DD,CONVF,CNVFMX,CONVX,CNVXMX,
     $             EIGMAX,EIGMIN,vmode(50),
     $             X(50),XNAME(50),OLDF(50),D(50),HESS(50,50),
     $             IC(50),MODE,NSTEP,MXSTEP,IHESS,IS,
     $             NEGREQ,IUPDAT,
     $             MXHESS,PRNT,NRFLAG,DUMP,IDUM
      DATA DELTA/0.005D0/, ZERO/0.0D0/
C
C
C  SET OPTIMIZATION CONTROL PARAMETERS
C
      NEGREQ=1
      MXSTEP=MIN0(NVAR+20,40)
      EIGMAX=25.0D0
      EIGMIN=0.0001D0
      IUPDAT=0
      NRFLAG= IOP(19).NE.0
      PRNT = IOP(33).EQ.0
      DUMP = IOP(34).NE.0
      DMAX=0.3D0
      CONVF=0.0003D0
C
C  CHANGE DEFAULT VALUES IF INSTRUCTED BY OPTIONS
C
      IF(IOP(5).NE.0) NEGREQ=0
      IF(IOP(6).NE.0) MXSTEP=IOP(6)
      IF(IOP(13).NE.0) IUPDAT=IOP(13)
      IF(IOP(16).NE.0) EIGMAX=0.1D0*IOP(16)
      IF(IOP(17).NE.0) EIGMIN=1.0D0/IOP(17)
      IF(IOP(8).NE.0) DMAX=0.01D0*IOP(8)
      IF(IOP(7).NE.0) CONVF=0.001D0/IOP(7)
C
C  CALCULATE REMAINING CONVERGENCE FACTORS
C
      CNVFMX=CONVF*1.5D0
      CONVX=CONVF*4.0D0
      CNVXMX=CONVX*1.5D0
C
C  CHECK HESSIAN MODE FOLLOWING
C
cets071891      IF(MODE.GT.0 .and. nstep.lt.1) THEN
      if(mode.gt.0) then
        IF(MODE.GT.NVAR) THEN
          IHOP=MODE-NVAR
          WRITE(IOUT,1025) IHOP
        ELSE
          WRITE(IOUT,1030) MODE
        ENDIF
        IF(MODE.GT.2*NVAR) WRITE(IOUT,1040)
        IF(MODE.GT.2*NVAR) stop
      ENDIF
C
      RETURN
C
 500  FORMAT(8F10.6)
 510  FORMAT(2I3,F20.0)
1000  FORMAT(//' ************************************************'/
     $         ' ** ERROR IN INITEF. NUMBER OF VARIABLES (',I3,') **'/
     $         ' **   INCORRECT (SHOULD BE BETWEEN 1 AND 50)   **'/
     $         ' ************************************************'//)
1005  FORMAT(//' CALL FROM INITEF:  BAD IC ENTRY FOR VARIABLE ',I3)
1010  FORMAT(//' CALL FROM INITEF:  OPTION 10 OUT OF RANGE ',
     $       ' IOP(10) = ',I4)
1020  FORMAT(//' ***************************************************'/
     $         ' ** O P T I M I Z A T I O N   P A R A M E T E R S **'/
     $         ' **  MAXIMUM NO. OF ITERATIONS      :  ',I10,' **'/
     $         ' **  MAXIMUM ALLOWED STEP LENGTH    :  ',F10.6,' **'/
     $         ' **  STEPS BEFORE HESSIAN REEVALUATED :  ',I8,' **'/
     $         ' ***************************************************'//)
1025  FORMAT(//' HESSIAN MODE FOLLOWING SWITCHED ON',/
     $       ' FOLLOWING MODE ',I3)
1030  FORMAT(//' HESSIAN MODE FOLLOWING SWITCHED ON',/
     $       ' FOLLOWING MODE WITH LARGEST MAGNITUDE COMPONENT',
     $       ' FOR VARIABLE ',I3)
1040  FORMAT(//' CALL FROM INITEF:  MODE FOLLOWING OUT OF RANGE')

      END
