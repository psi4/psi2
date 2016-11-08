      SUBROUTINE SYMBOL(S,IGROUP,NAXIS)
      IMPLICIT INTEGER (A-Z)
C
      DIMENSION S(14)
C
      DATA A,AP,APP,A1,A2,AG,AU,A1P,A2P,A1PP,A2PP,A1G,A2G,A1U,A2U
     1/4HA   ,4HA'  ,4HA"  ,4HA1  ,4HA2  ,4HAG  ,4HAU  ,4HA1' ,4HA2'
     2,4HA1" ,4HA2" ,4HA1G ,4HA2G ,4HA1U ,4HA2U /
C
      DATA B,BG,BU,B1,B2,B3,B1G,B2G,B3G,B1U,B2U,B3U
     1/4HB   ,4HBG  ,4HBU  ,4HB1  ,4HB2  ,4HB3  ,4HB1G ,4HB2G ,4HB3G
     2,4HB1U ,4HB2U ,4HB3U /
C
      DATA E,E1,E2,E3,EG,EU,EP,EPP,E1P,E2P,E1PP,E2PP,E1G,E2G,E1U,E2U
     1,E4,E5,E3G,E3U
     2/4HE   ,4HE1  ,4HE2  ,4HE3  ,4HEG  ,4HEU  ,4HE'  ,4HE"  ,4HE1'
     3,4HE2' ,4HE1" ,4HE2" ,4HE1G ,4HE2G ,4HE1U ,4HE2U ,4HE4  ,4HE5
     4,4HE3G ,4HE3U /
C
      DATA F,FG,FU,F1,F2,F1G,F1U,F2G,F2U,G,H,GG,GU,HG,HU/4HF   ,
     14HFG  ,4HFU  ,4HF1  ,4HF2  ,4HF1G ,4HF1U ,4HF2G ,4HF2U ,
     24HG   ,4HH   ,4HGG  ,4HGU  ,4HHG  ,4HHU  /
      IF(IGROUP.GE.13)GO TO 2000
      GO TO (100,200,300,400,500,600,700,800,900,1000),IGROUP
C
C     -- C1 -
C
  100 S(1)=A
      GO TO 5000
C
C     -- CS -
C
  200 S(1)=AP
      S(2)=APP
      GO TO 5000
C
C     -- CI -
C
  300 S(1)=AG
      S(2)=AU
      GO TO 5000
C
C     -- CN -
C
  400 S(1)=A
      GO TO (100,420,430,440,450,460,470,480),NAXIS
  420 S(2)=B
      GO TO 5000
  430 S(2)=E
      GO TO 5000
  440 S(2)=B
      S(3)=E
      GO TO 5000
  450 S(2)=E1
      S(3)=E2
      GO TO 5000
  460 S(2)=B
      S(3)=E1
      S(4)=E2
      GO TO 5000
  470 S(2)=E1
      S(3)=E2
      S(4)=E3
      GO TO 5000
  480 S(2)=B
      S(3)=E1
      S(4)=E2
      S(5)=E3
      GO TO 5000
C
C     -- S2N -
C
  500 S(1)=A
      GO TO (300,520,530,540),NAXIS
  520 S(2)=B
      S(3)=E
      GO TO 5000
  530 S(1)=AG
      S(2)=AU
      S(3)=EG
      S(4)=EU
      GO TO 5000
  540 S(2)=B
      S(3)=E1
      S(4)=E2
      S(5)=E3
      GO TO 5000
C
C     -- CNH -
C
  600 GO TO (200,620,630,640,650,660),NAXIS
  620 S(1)=AG
      S(2)=BG
      S(3)=AU
      S(4)=BU
      GO TO 5000
  630 S(1)=AP
      S(2)=APP
      S(3)=EP
      S(4)=EPP
      GO TO 5000
  640 S(1)=AG
      S(2)=BG
      S(3)=AU
      S(4)=BU
      S(5)=EG
      S(6)=EU
      GO TO 5000
  650 S(1)=AP
      S(2)=APP
      S(3)=E1P
      S(4)=E2P
      S(5)=E1PP
      S(6)=E2PP
      GO TO 5000
  660 S(1)=AG
      S(2)=BG
      S(3)=AU
      S(4)=BU
      S(5)=E1G
      S(6)=E2G
      S(7)=E1U
      S(8)=E2U
      GO TO 5000
C
C     -- CNV -
C
  700 S(1)=A1
      S(2)=A2
      GO TO (200,720,730,740,750,760),NAXIS
  720 S(3)=B1
      S(4)=B2
      GO TO 5000
  730 S(3)=E
      GO TO 5000
  740 S(3)=B1
      S(4)=B2
      S(5)=E
      GO TO 5000
  750 S(3)=E1
      S(4)=E2
      GO TO 5000
  760 S(3)=B1
      S(4)=B2
      S(5)=E1
      S(6)=E2
      GO TO 5000
C
C     -- DN -
C
  800 S(1)=A1
      S(2)=A2
      GO TO (420,820,830,840,850,860),NAXIS
  820 S(2)=B1
      S(3)=B2
      S(4)=B3
      GO TO 5000
  830 S(3)=E
      GO TO 5000
  840 S(3)=B1
      S(4)=B2
      S(5)=E
      GO TO 5000
  850 S(3)=E1
      S(4)=E2
      GO TO 5000
  860 S(3)=B1
      S(4)=B2
      S(5)=E1
      S(6)=E2
      GO TO 5000
C
C     -- DNH -
C
  900 GO TO (720,920,930,940,950,960,970,980),NAXIS
  920 S(1)=AG
      S(2)=B1G
      S(3)=B2G
      S(4)=B3G
      S(5)=AU
      S(6)=B1U
      S(7)=B2U
      S(8)=B3U
      GO TO 5000
  930 S(1)=A1P
      S(2)=A2P
      S(3)=A1PP
      S(4)=A2PP
      S(5)=EP
      S(6)=EPP
      GO TO 5000
  940 S(1)=A1G
      S(2)=A2G
      S(3)=B1G
      S(4)=B2G
      S(5)=A1U
      S(6)=A2U
      S(7)=B1U
      S(8)=B2U
      S(9)=EG
      S(10)=EU
      GO TO 5000
  950 S(1)=A1P
      S(2)=A2P
      S(3)=A1PP
      S(4)=A2PP
      S(5)=E1P
      S(6)=E2P
      S(7)=E1PP
      S(8)=E2PP
      GO TO 5000
  960 S(1)=A1G
      S(2)=A2G
      S(3)=B1G
      S(4)=B2G
      S(5)=A1U
      S(6)=A2U
      S(7)=B1U
      S(8)=B2U
      S(9)=E1G
      S(10)=E2G
      S(11)=E1U
      S(12)=E2U
      GO TO 5000
  970 STOP
  980 S(1)=A1G
      S(2)=A2G
      S(3)=B1G
      S(4)=B2G
      S(5)=A1U
      S(6)=A2U
      S(7)=B1U
      S(8)=B2U
      S(9)=E1G
      S(10)=E2G
      S(11)=E3G
      S(12)=E1U
      S(13)=E2U
      S(14)=E3U
      GO TO 5000
C
C     -- DND -
C
 1000 GO TO (620,1020,1030,1040,1050,1060),NAXIS
 1020 S(1)=A1
      S(2)=A2
      S(3)=B1
      S(4)=B2
      S(5)=E
      GO TO 5000
 1030 S(1)=A1G
      S(2)=A2G
      S(3)=A1U
      S(4)=A2U
      S(5)=EG
      S(6)=EU
      GO TO 5000
 1040 S(1)=A1
      S(2)=A2
      S(3)=B1
      S(4)=B2
      S(5)=E1
      S(6)=E2
      S(7)=E3
      GO TO 5000
 1050 S(1)=A1G
      S(2)=A2G
      S(3)=A1U
      S(4)=A2U
      S(5)=E1G
      S(6)=E2G
      S(7)=E1U
      S(8)=E2U
      GO TO 5000
 1060 S(1)=A1
      S(2)=A2
      S(3)=B1
      S(4)=B2
      S(5)=E1
      S(6)=E2
      S(7)=E3
      S(8)=E4
      S(9)=E5
      GO TO 5000
 2000 NDEX=IGROUP-12
      GO TO (2100,2200,2300,2400,2500,2600,2700)NDEX
C
C
C   ---  TETRAHEDRAL GROUPS  ---
C
C        --- T ---
C
 2100 S(1)=A
      S(2)=E
      S(3)=F
      GO TO 5000
C
C       --- TH ---
C
 2200 S(1)=AG
      S(2)=EG
      S(3)=FG
      S(4)=AU
      S(5)=EU
      S(6)=FU
      GO TO 5000
C
C       --- TD ---
C
 2300 S(1)=A1
      S(2)=A2
      S(3)=E
      S(4)=F1
      S(5)=F2
      GO TO 5000
C
C
C   ---  OCTAHEDRAL GROUPS  ---
C
C       --- O ---
 2400 S(1)=A1
      S(2)=A2
      S(3)=E
      S(4)=F1
      S(5)=F2
      GO TO 5000
C
C       --- OH ---
C
 2500 S(1)=A1G
      S(2)=A2G
      S(3)=EG
      S(4)=F1G
      S(5)=F2G
      S(6)=A1U
      S(7)=A2U
      S(8)=EU
      S(9)=F1U
      S(10)=F2U
      GO TO 5000
C
C
C   ---  ICOSAHEDRAL GROUPS  ---
C
C       --- I ---
 2600 S(1)=A
      S(2)=F1
      S(3)=F2
      S(4)=G
      S(5)=H
      GO TO 5000
C
C       --- IH ---
 2700 S(1)=AG
      S(2)=F1G
      S(3)=F2G
      S(4)=GG
      S(5)=HG
      S(6)=AU
      S(7)=F1U
      S(8)=F2U
      S(9)=GU
      S(10)=HU
 5000 CONTINUE
      RETURN
      END
