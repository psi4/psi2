C
      SUBROUTINE TWRIT (T1,T2,ORBSYM,FLOV,NIRRED,TOFF,TADD,IVAL,TVAL,
     .                  ZLX,FZO,FZV,NO,NV,NSYMHF,ITYP,NT1,NT2,NCSF,NCL)
      IMPLICIT INTEGER(A-Z)
      REAL*8 T2(*),T1(NO,NV),TVAL(*),TMAX
      INTEGER ORBSYM(NO+NV),FLOV(8,2,4),TOFF(NO,NO),TADD(NV,NV),
     .        ZLX(NV,NV),FZO(NO),FZV(NV),IVAL(*),ITYP(NSYMHF),NCL(*)
      character*6 spen(2)
      integer xspin(4),yspin(4)
 
C
      NCSF1=NCSF
      IF(NCSF1.GT.NT1)NCSF1=NT1
c
      WRITE(6,*)
      WRITE(6,*)' T1 VECTOR alpha'
      WRITE(6,*)
      WRITE(6,*)'    I     T1( I )     CSF '
      WRITE(6,*)
C
      DO 880 CSF=1,NCSF1
         TMAX = 0.0D0
         IMAX = 1
C
      DO 790 USYM=1,NSYMHF
         FU=FLOV(USYM,1,1)
         LU=FLOV(USYM,1,2)
         SU=ITYP(USYM)
         DO 780 U=FU,LU
            DO 770 BESYM=1,NSYMHF
               UBESYM=IEOR(USYM-1,BESYM-1)+1
               IF(UBESYM.NE.1)GO TO 770
               FBE=FLOV(BESYM,1,3)-NO
               LBE=FLOV(BESYM,1,4)-NO
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
               XU=U-FU+1
               XB=BE-FBE+1+NCL(BESYM)
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
c
      WRITE(6,*)
      WRITE(6,*)' T1 VECTOR beta'
      WRITE(6,*)
      WRITE(6,*)'    I     T1( I )     CSF '
      WRITE(6,*)
C
      DO 981 CSF=1,NCSF1
         TMAX = 0.0D0
         IMAX = 1
C
      DO 990 USYM=1,NSYMHF
         FU=FLOV(USYM,2,1)
         LU=FLOV(USYM,2,2)
         SU=ITYP(USYM)
         DO 980 U=FU,LU
            DO 970 BESYM=1,NSYMHF
               UBESYM=IEOR(USYM-1,BESYM-1)+1
               IF(UBESYM.NE.1)GO TO 970
               FBE=FLOV(BESYM,2,3)-NO
               LBE=FLOV(BESYM,2,4)-NO
               SBE=ITYP(BESYM)
               DO 960 BE=FBE,LBE
                  XX=(U-1)*NV+BE
            FLAG=0
            DO 940 J=1,CSF-1
               IF(XX.EQ.IVAL(J))FLAG=1
  940       CONTINUE
            IF(FLAG.EQ.1)GO TO 960
C              WRITE(6,'(4I6,F12.8)')CSF,U,BE,XX,T1(U,BE)
C
            IF(DABS(T1(U,BE)).GT.DABS(TMAX))THEN
               TMAX=T1(U,BE)
               XMAX=XX
               XU=U-FU+1
               XB=BE-FBE+1+NCL(BESYM)
               SSU=SU
               SSB=SBE
            ENDIF
C
  960          CONTINUE
  970       CONTINUE
  980    CONTINUE
  990 CONTINUE
C
      TVAL(CSF)=TMAX
      IVAL(CSF)=XMAX
      WRITE(6,612)CSF,TMAX,XMAX,XU,SSU,XB,SSB
  981 CONTINUE
c
  612 FORMAT(I6,F12.6,'  (',I7,')',5X,(I3,A4,1X),'  --- >  ',(I3,A4,1X))
C
C
      spen(1)=' alpha'
      spen(2)=' beta '
      xspin(1)=1
      yspin(1)=1
      xspin(2)=2
      yspin(2)=1
      xspin(3)=2
      yspin(3)=2
c
      IF(NCSF.GT.NT2)NCSF=NT2
      do 1000 case=1,3
      spin1=xspin(case)
      spin2=yspin(case)
      WRITE(6,*)
      WRITE(6,638)spen(spin1),spen(spin1),spen(spin2),spen(spin2)
  638 format(' T2 VECTOR ',4(3x,a6))
      WRITE(6,*)
      WRITE(6,*)'    I     T2( I )     CSF '
      WRITE(6,*)
C
      DO 180 CSF=1,NCSF
         TMAX = 0.0D0
         IMAX = 1
C
      DO 290 USYM=1,NSYMHF
         FU=FLOV(USYM,spin1,1)
         LU=FLOV(USYM,spin1,2)
         SU=ITYP(USYM)
         DO 280 U=FU,LU
            DO 270 BESYM=1,NSYMHF
               UBESYM=IEOR(USYM-1,BESYM-1)
               FBE=FLOV(BESYM,spin1,3)-NO
               LBE=FLOV(BESYM,spin1,4)-NO
               SBE=ITYP(BESYM)
               DO 260 BE=FBE,LBE
                  DO 250 VSYM=1,USYM
                     TSYM=IEOR(VSYM-1,UBESYM)+1
                     FV=FLOV(VSYM,spin2,1)
                     LV=FLOV(VSYM,spin2,2)
                     SV=ITYP(VSYM)
                     IF(VSYM.EQ.USYM)LV=U
                     DO 240 V=FV,LV
                        DO 230 GASYM=1,NSYMHF
                           IF(TSYM.NE.GASYM)GO TO 230
                           SGA=ITYP(GASYM)
                           FGA=FLOV(GASYM,spin2,3)-NO
                           LGA=FLOV(GASYM,spin2,4)-NO
                           IF(U.EQ.V)LGA=BE
                           DO 220 GA=FGA,LGA
                              ZL=ZLX(BE,GA)
                              XX=TOFF(U,V)+TADD(BE,GA)
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
               XU=U-FU+1
               XV=V-FV+1
               XB=BE-FBE+1+NCL(BESYM)
               XG=GA-FGA+1+NCL(GASYM)
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
      WRITE(6,602)CSF,TMAX,XMAX,XV,SSV,XU,SSU,XG,SSG,XB,SSB
C     WRITE(6,602)CSF,TMAX,XMAX,XU,SSU,XV,SSV,XB,SSB,XG,SSG
  180 CONTINUE
 1000 continue
  602 FORMAT(I6,F12.6,'  (',I7,')',2(I3,A4,1X),'  --- >  ',2(I3,A4,1X))
C
      RETURN
      END
