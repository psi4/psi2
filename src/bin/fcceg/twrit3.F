      SUBROUTINE TWRIT3(T1,T2,ORBSYM,FLOV,NIRRED,TOFF,TADD,IVAL,TVAL,
     .                  ZLX,FZO,FZV,NO,NV,NSYMHF,ITYP,NT1,NT2,NCSF,
     .                  t3,ndimt3,i1off,i2off,i3off,trsym)
      IMPLICIT INTEGER(A-Z)
      REAL*8 T2(*),T1(NO,NV),TVAL(*),TMAX,T3(*)
      INTEGER ORBSYM(NO+NV),FLOV(NIRRED,4),TOFF(NO,NO,2),TADD(NV,NV),
     .        ZLX(NV,NV),FZO(NO),FZV(NV),IVAL(*),ITYP(NSYMHF),
     .        i1off(*),i2off(no*nv,nirred),i3off(*),trsym(*) 
 
C
      NCSF1=NCSF
      IF(NCSF.GT.NT1)NCSF1=NT1
      WRITE(6,*)
      WRITE(6,*)' T1 VECTOR '
      WRITE(6,*)
      WRITE(6,*)'    I     T1( I )     CSF '
      WRITE(6,*)
C
      DO 880 CSF=1,NCSF1
         TMAX = 0.0D0
         IMAX = 1
C
      DO 790 USYM=1,NSYMHF
         FU=FLOV(USYM,1)
         LU=FLOV(USYM,2)
         SU=ITYP(USYM)
         DO 780 U=FU,LU
            DO 770 BESYM=1,NSYMHF
               UBESYM=IEOR(USYM-1,BESYM-1)+1
               IF(UBESYM.NE.1)GO TO 770
               FBE=FLOV(BESYM,3)-NO
               LBE=FLOV(BESYM,4)-NO
               SBE=ITYP(BESYM)
               DO 760 BE=FBE,LBE
                  XX=(U-1)*NV+BE
            FLAG=0
            DO 740 J=1,CSF-1
               IF(XX.EQ.IVAL(J))FLAG=1
  740       CONTINUE
            IF(FLAG.EQ.1)GO TO 760
C              WRITE(6,'(4I6,F12.8)')CSF,U,BE,XX,T1(U,BE)
C
            IF(DABS(T1(U,BE)).GT.DABS(TMAX))THEN
               TMAX=T1(U,BE)
               XMAX=XX
               XU=U
               XB=BE+NO
               SSU=SU
               SSB=SBE
            ENDIF
C
  760          CONTINUE
  770       CONTINUE
  780    CONTINUE
  790 CONTINUE
C
      TVAL(CSF)=TMAX
      IVAL(CSF)=XMAX
      WRITE(6,612)CSF,TMAX,XMAX,XU,SSU,XB,SSB
  880 CONTINUE
  612 FORMAT(I6,F12.6,'  (',I7,')',5X,(I3,A4,1X),'  --- >  ',(I3,A4,1X))
C
C
      NCSF2=NCSF
      IF(NCSF.GT.NT2)NCSF2=NT2
      WRITE(6,*)
      WRITE(6,*)' T2 VECTOR '
      WRITE(6,*)
      WRITE(6,*)'    I     T2( I )     CSF '
      WRITE(6,*)
C
      DO 180 CSF=1,NCSF2
         TMAX = 0.0D0
         IMAX = 1
C
      DO 290 USYM=1,NSYMHF
         FU=FLOV(USYM,1)
         LU=FLOV(USYM,2)
         SU=ITYP(USYM)
         DO 280 U=FU,LU
            DO 270 BESYM=1,NSYMHF
               UBESYM=IEOR(USYM-1,BESYM-1)
               FBE=FLOV(BESYM,3)-NO
               LBE=FLOV(BESYM,4)-NO
               SBE=ITYP(BESYM)
               DO 260 BE=FBE,LBE
                  DO 250 VSYM=1,USYM
                     TSYM=IEOR(VSYM-1,UBESYM)+1
                     FV=FLOV(VSYM,1)
                     LV=FLOV(VSYM,2)
                     SV=ITYP(VSYM)
                     IF(VSYM.EQ.USYM)LV=U
                     DO 240 V=FV,LV
                        DO 230 GASYM=1,NSYMHF
                           IF(TSYM.NE.GASYM)GO TO 230
                           SGA=ITYP(GASYM)
                           FGA=FLOV(GASYM,3)-NO
                           LGA=FLOV(GASYM,4)-NO
                           IF(U.EQ.V)LGA=BE
                           DO 220 GA=FGA,LGA
                              ZL=ZLX(BE,GA)
                              XX=TOFF(U,V,ZL)+TADD(BE,GA)
            FLAG=0
            DO 140 J=1,CSF-1
               IF(XX.EQ.IVAL(J))FLAG=1
  140       CONTINUE
            IF(FLAG.EQ.1)GO TO 220
C              WRITE(6,'(6I6,F12.8)')CSF,U,V,BE,GA,XX,T2(XX)
C
            IF(DABS(T2(XX)).GT.DABS(TMAX))THEN
               TMAX=T2(XX)
               XMAX=XX
               XU=U
               XV=V
               XB=BE+NO
               XG=GA+NO
               SSU=SU
               SSV=SV
               SSB=SBE
               SSG=SGA
            ENDIF
C
  220                      CONTINUE
  230                   CONTINUE
  240                CONTINUE
  250             CONTINUE
  260          CONTINUE
  270       CONTINUE
  280    CONTINUE
  290 CONTINUE
C
      TVAL(CSF)=TMAX
      IVAL(CSF)=XMAX
      WRITE(6,602)CSF,TMAX,XMAX,XU,SSU,XV,SSV,XB,SSB,XG,SSG
  180 CONTINUE
  602 FORMAT(I6,F12.6,'  (',I7,')',2(I3,A4,1X),'  --- >  ',2(I3,A4,1X))
C
C
      NCSF3=NCSF
      IF(NCSF.GT.NT3)NCSF3=NT3
      WRITE(6,*)
      WRITE(6,*)' T3 VECTOR '
      WRITE(6,*)
      WRITE(6,*)'    I     T3( I )     CSF '
      WRITE(6,*)
C
      DO 495 CSF=1,NCSF3
         TMAX = 0.0D0
         IMAX = 1
C
      DO 490 USYM=1,NSYMHF
         FU=FLOV(USYM,1)
         LU=FLOV(USYM,2)
         SU=ITYP(USYM)
         DO 480 U=FU,LU
            DO 470 BESYM=1,NSYMHF
               UBESYM=IEOR(USYM-1,BESYM-1)
               FBE=FLOV(BESYM,3)-NO
               LBE=FLOV(BESYM,4)-NO
               SBE=ITYP(BESYM)
               DO 460 BE=FBE,LBE
                  beu=(be-1)*no+u
                  DO 450 VSYM=1,NSYMHF
                     UVBSYM=IEOR(VSYM-1,UBESYM)
                     FV=FLOV(VSYM,1)
                     LV=FLOV(VSYM,2)
                     SV=ITYP(VSYM)
                     IF(VSYM.EQ.USYM)LV=U
                     DO 440 V=FV,LV
                        DO 430 GASYM=1,NSYMHF
                           TSYM=IEOR(UVBSYM,GASYM-1)
                           SGA=ITYP(GASYM)
                           FGA=FLOV(GASYM,3)-NO
                           LGA=FLOV(GASYM,4)-NO
C                          IF(U.EQ.V)LGA=BE
                           DO 420 GA=FGA,LGA
                              gav=(ga-1)*no+v
                  DO 415 WSYM=1,NSYMHF
                     FW=FLOV(WSYM,1)
                     LW=FLOV(WSYM,2)
                     SW=ITYP(WSYM)
C                    IF(WSYM.EQ.VSYM)LW=V
                     DO 413 W=FW,LW
                        DO 411 EPSYM=1,NSYMHF
                           WEPSYM=IEOR(WSYM-1,EPSYM-1)
                           IF(WEPSYM.NE.TSYM)GO TO 411
                           SEP=ITYP(EPSYM)
                           FEP=FLOV(EPSYM,3)-NO
                           LEP=FLOV(EPSYM,4)-NO
C                          IF(W.EQ.V)LEP=GA
                           DO 409 EP=FEP,LEP
                              epw=(ep-1)*no+w
c                             XX=T3OFF(U,V,W)+T3ADD(BE,GA,EP)
                           k1=max0(beu,gav,epw)
                           k3=min0(beu,gav,epw)
                           k3s=trsym(k3)
                           k2=beu+gav+epw-k1-k3
                           xx=i1off(k1)+i2off(k2,k3s)+i3off(k3)
            FLAG=0
            DO 404 J=1,CSF-1
               IF(DABS(T3(XX)-TVAL(J)).LT.1.D-10)FLAG=1
  404       CONTINUE
            IF(FLAG.EQ.1)GO TO 409
C              WRITE(6,'(6I6,F12.8)')CSF,U,V,BE,GA,XX,T2(XX)
C
            IF(DABS(T3(XX)).GT.DABS(TMAX))THEN
               TMAX=T3(XX)
               XMAX=XX
               XU=U
               XV=V
               XW=W
               XB=BE+NO
               XG=GA+NO
               XE=EP+NO
               SSU=SU
               SSV=SV
               SSW=SW
               SSB=SBE
               SSG=SGA
               SSE=SEP
            ENDIF
C
  409                      CONTINUE
  411                      CONTINUE
  413                   CONTINUE
  415                CONTINUE
  420                      CONTINUE
  430                   CONTINUE
  440                CONTINUE
  450             CONTINUE
  460          CONTINUE
  470       CONTINUE
  480    CONTINUE
  490 CONTINUE
C
      TVAL(CSF)=TMAX
      IVAL(CSF)=XMAX
      WRITE(6,632)CSF,TMAX,XMAX,XU,SSU,XV,SSV,XW,SSW,
     .                          XB,SSB,XG,SSG,XE,SSE
  495 CONTINUE
  632 FORMAT(I6,F12.6,'  (',I7,')',3(I3,A4,1X),'  --- >  ',3(I3,A4,1X))
C
      RETURN
      END
