      SUBROUTINE REDUN2(II,JJ,KK,LL,Q4,KATOM)
C     ----- THIS VERSION HALVES Q4 TO ACCOUNT FOR DUPLICATION -----
C                     OF [II,JJ] OR [IJ,IJ]
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION KATOM(200)
      LOGICAL IANDJ,IANDK,IANDL,JANDK,JANDL,KANDL
      COMMON/PARA3/NATMS(4),NPASS
      COMMON/SHELL/ISH,JSH,KSH,LSH
      COMMON/TOLER/DLIM
      DATA HALF / 0.5D+00 /
C
      DO 101 I=1,4
  101 NATMS(I)=0
C
      IAT=KATOM(II)
      JAT=KATOM(JJ)
      KAT=KATOM(KK)
      LAT=KATOM(LL)
      IANDJ=IAT.EQ.JAT
      IANDK=IAT.EQ.KAT
      IANDL=IAT.EQ.LAT
      JANDK=JAT.EQ.KAT
      JANDL=JAT.EQ.LAT
      KANDL=KAT.EQ.LAT
      IF(.NOT.IANDJ) GO TO 300
      IF(.NOT.IANDK) GO TO 200
      IF(.NOT.IANDL) GO TO 100
C
C     ----- IAT = JAT = KAT = LAT -----
C
      NPASS=0
      GO TO 1000
C
C     ----- IAT = JAT = KAT _ LAT ----
C
  100 CONTINUE
      NATMS(1)=LAT
      NATMS(2)=IAT
      ISH=LL
      JSH=KK
      KSH=JJ
      LSH=II
      NPASS=4
      GO TO 1000
C
C     ----- IAT = JAT _ KAT = LAT ----
C
  200 IF(.NOT.KANDL) GO TO 250
      IF(KK.EQ.LL) Q4=Q4*HALF
      NATMS(1)=KAT
      NATMS(2)=LAT
      NATMS(3)=IAT
      ISH=KK
      JSH=LL
      KSH=II
      LSH=JJ
      NPASS=2
      GO TO 1000
C
C     ----- IAT = JAT _ KAT _ LAT -----
C
  250 CONTINUE
      NATMS(1)=KAT
      NATMS(2)=LAT
      NATMS(3)=IAT
      ISH=KK
      JSH=LL
      KSH=II
      LSH=JJ
      NPASS=2
      GO TO 1000
C
C     ----- IAT = KAT = LAT _ JAT -----
C
  300 IF(.NOT.IANDK) GO TO 500
      IF(.NOT.IANDL) GO TO 400
      NATMS(1)=JAT
      NATMS(2)=IAT
      ISH=JJ
      JSH=II
      KSH=KK
      LSH=LL
      NPASS=4
      GO TO 1000
C
C     ----- IAT = KAT _ JAT = LAT -----
C     ----- IAT = KAT _ JAT _ LAT -----
C
  400 CONTINUE
      IF(II.EQ.KK.AND.JJ.EQ.LL) Q4=Q4*HALF
      NATMS(1)=JAT
      NATMS(2)=LAT
      NATMS(3)=IAT
      ISH=JJ
      JSH=II
      KSH=LL
      LSH=KK
      NPASS=3
      GO TO 1000
C
C     ----- JAT = KAT = LAT _ IAT -----
C
  500 IF(.NOT.JANDK) GO TO 700
      IF(.NOT.JANDL) GO TO 600
      NATMS(1)=IAT
      NATMS(2)=JAT
      ISH=II
      JSH=JJ
      KSH=KK
      LSH=LL
      NPASS=4
      GO TO 1000
C
C     ----- JAT = KAT _ IAT _ LAT -----
C
  600 CONTINUE
      NATMS(1)=IAT
      NATMS(2)=LAT
      NATMS(3)=JAT
      ISH=II
      JSH=JJ
      KSH=LL
      LSH=KK
      NPASS=3
      GO TO 1000
C
C     ----- JAT = LAT _ IAT _ KAT -----
C
  700 IF(.NOT.JANDL) GO TO 800
      NATMS(1)=IAT
      NATMS(2)=KAT
      NATMS(3)=JAT
      ISH=II
      JSH=JJ
      KSH=KK
      LSH=LL
      NPASS=3
      GO TO 1000
C
C     ----- KAT = LAT _ IAT _ JAT ----
C
  800 IF(.NOT.KANDL) GO TO 900
      NATMS(1)=IAT
      NATMS(2)=JAT
      NATMS(3)=KAT
      ISH=II
      JSH=JJ
      KSH=KK
      LSH=LL
      NPASS=2
      GO TO 1000
C
C     ----- IAT _ JAT _ KAT _ LAT ----
C
  900 CONTINUE
      NATMS(1)=IAT
      NATMS(2)=JAT
      NATMS(3)=KAT
      NATMS(4)=LAT
      ISH=II
      JSH=JJ
      KSH=KK
      LSH=LL
      NPASS=1
C
C     -----
C
 1000 CONTINUE
      RETURN
      END
