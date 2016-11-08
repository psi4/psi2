      SUBROUTINE MOGRAD(X,NIJ,GMO)                                      SMM15310
      IMPLICIT REAL*8 (A-H,O-Z)                                         SMM15320
      DIMENSION X(NBASIS,NBASIS),NIJ(NTRI,2),GMO(NTRI)                  SMM15330
      COMMON/BASIS/NBASIS,NTRI,NTOT,NSYM
      COMMON/GRADS/GLIM,GMNORM
      COMMON/SIGNS/IOFF(32897),IPRNT 
      COMMON/CI102/NIND,NINTRI
      DATA A00 / 0.0D+00 /                                              SMM15380
    1 FORMAT(//,2X, ' INDEPENDENT PAIRS'/)                              SMM15400
    2 FORMAT(2X,I5,3X,2I5,3X,F20.10)                                    SMM15410
    3 FORMAT(//,2X,' NO. OF I.P. = ',I5,5X,' GNORM = ',F20.10,5X,
     *             ' GMAX = ',F20.10)
C                                                                       SMM15430
C   FORM THE GRADIENT VECTOR FOR MO ROTATIONS                           SMM15440
c=====================================================================
c   Note:
c     Pairs with the same symmetry of the Lagrangian matrix
c     are non-vanishing
c=====================================================================
      IF(IPRNT.LE.2) GO TO 201                                          SMM15450
      WRITE(6,1)                                                        SMM15460
  201 CONTINUE                                                          SMM15470
      II=0                                                              SMM15480
      imax=0
      gtmax=a00
      DO 101 I=2,NBASIS                                                 SMM15490
      DO 101 J=1,I-1                                                    SMM15500
      XGRAD=X(I,J)-X(J,I)                                               SMM15510
      IF(DABS(XGRAD).LT.GLIM) GO TO 101                                 SMM15520
      II=II+1                                                           SMM15530
      GMO(II)=XGRAD                                                     SMM15540
      if(dabs(xgrad).gt.gtmax) then
        imax=ii
        gtmax=dabs(xgrad)
      end if
      NIJ(II,1)=I                                                       SMM15550
      NIJ(II,2)=J                                                       SMM15560
      IF(IPRNT.LE.2) GO TO 101                                          SMM15570
      WRITE(6,2) II,I,J,XGRAD                                           SMM15580
  101 CONTINUE                                                          SMM15590
C                                                                       SMM15600
      NIND=II                                                           SMM15610
      GNORM=A00                                                         SMM15620
      IF(NIND.EQ.0) GO TO 202                                           SMM15630
      call anorm(gnorm,gmo,nind,nind)
      denom=dsqrt(dfloat(ntri))
      gnorm=gnorm/denom
  202 CONTINUE                                                          SMM15680
      GMNORM=GNORM                                                      SMM15690
      if(imax.eq.0) then
        gmax=a00
      else
        gmax=gmo(imax)
      end if
c     IF(IPRNT.LE.2) GO TO 203                              
      WRITE(6,3) NIND,GNORM,gmax
C                                                                       SMM15720
  203 CONTINUE                                                          SMM15730
      RETURN                                                            SMM15740
      END                                                               SMM15750
