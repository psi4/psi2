C=======================================================================
C  Modification for LINUX symstem
C  Subroutine Mwrit has been commetned out
C  on September 07, 2005 by Dr. Yukio Yamaguchi
C=======================================================================
      SUBROUTINE GETIT(KADD,LADD,IJADD,IJGRP,HQ,BUFI)
      IMPLICIT REAL*8 (A-H,O-Z)                                         DRT02210
      DIMENSION KADD(ISYMRB),LADD(ISYMRB),IJADD(NUMIJ),IJGRP(NUMIJ)     DRT02230
      DIMENSION HQ(NUMIJ),BUFI(NMAX)
      DIMENSION DUMS(56)                                                DRT02250
      LOGICAL IJEQ,IJKEQ,JKEQ                                           DRT02260
      COMMON/DIMES/ISYMRB,NUMIJ,NBF,NORBS,NMAX                          DRT02270
      COMMON/ORDCI/IDRSYM(256),IORDRT(256)                              DRT02300
      COMMON/ORDSO/IORSO(256),IDRTSO(256)                               DRT02310
      common/posit/iopci,iopmc
      COMMON/SIGNS/IOFF(256),IPRNT                                      DRT02320
      common/tapes/itap53,itap58
      DATA A00,ONE / 0.0D+00 , 1.0D+00 /                                DRT02330
      DATA D8, D4, D2 / 0.125D0, 0.25D0, 0.5D0 /                        DRT02340
    1 FORMAT(//,2X,' GET DRT MATRICES'/)                                DRT02350
    2 FORMAT(1X,7F18.8)                                                 DRT02360
    3 FORMAT(1X,'  IN GETDM, IJKLT UNKNOWN TYPE ',I5/'  ORBITALS :',4I5/DRT02370
     * '  BASIS FUNCTIONS :',4I5/'   SYMMETRIES :',4I5)                 DRT02380
    4 FORMAT(1X,'  IN GETDM, LAD TOO BIG :',I6)                         DRT02390
    5 FORMAT(2X,4I5,3X,I5,3X,3F20.10)                                   DRT02400
    6 FORMAT(//,2X,' HQ MATRIX'/)                                       DRT02410
C                                                                       DRT02420
      CALL ZERO(HQ,NUMIJ)                                               DRT02430
C                                                                       DRT02440
C   READ IN THE ONE AND TWO PARTICLE DENSITY MATRICES IN DRT ORDERING   DRT02450
c
      itape=itap53
      call rfile(itape)
  201 CONTINUE                                                          DRT02530
      NMAX2=NMAX*2                                                      DRT02540
      IBL=0                                                             DRT02550
      ISRIN=0                                                           DRT02560
      IGRP=0                                                            DRT02570
      DO 140 IOR=NORBS,1,-1                                             DRT02580
      ISM=IDRSYM(IOR)                                                   DRT02590
      ITYP=1                                                            DRT02600
      DO 130 JOR=IOR,1,-1                                               DRT02610
      IF (IGRP .NE. IJGRP(IOFF(IOR)+JOR)) THEN                          DRT02620
      CALL SREAD(ITAPE,BUFI,NMAX2)                                      DRT02630
      IF(IPRNT.GT.7)                                                    DRT02640
     *WRITE(6,2) (BUFI(IOP), IOP=1,NMAX)                                DRT02650
      ISRIN=ISRIN+NMAX/112                                              DRT02660
      IGRP=IJGRP(IOFF(IOR)+JOR)                                         DRT02670
      END IF                                                            DRT02680
      JSM=IDRSYM(JOR)                                                   DRT02690
      IJEQ=IOR.EQ.JOR                                                   DRT02700
      IJT=ITYP                                                          DRT02710
      IF (IJEQ) IJT=IJT+2                                               DRT02720
      IJSYM=IEOR((ISM-1),(JSM-1))                                       DRT02730
      KOFF=IJSYM*NORBS                                                  DRT02740
      JAD=IJADD(IOFF(IOR)+JOR)                                          DRT02750
      DO 120 KOR=1,JOR                                                  DRT02760
      KSM=IDRSYM(KOR)                                                   DRT02770
      JKEQ=JOR.EQ.KOR                                                   DRT02780
      IJKT=IJT                                                          DRT02790
      IF (JKEQ) IJKT=IJKT+1                                             DRT02800
      IJKEQ=IJEQ.AND.JKEQ                                               DRT02810
      IJKSYM=IEOR(IJSYM,(KSM-1))                                        DRT02820
      LOFF=IJKSYM*NORBS                                                 DRT02830
      KAD=JAD+KADD(KOFF+KOR)                                            DRT02840
      DO 110 LOR=1,KOR                                                  DRT02850
      IF (JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GO TO 310       DRT02860
      LSM=IDRSYM(LOR)                                                   DRT02870
      IF ((LSM-1) .NE. IJKSYM) GO TO 310                                DRT02880
      IJKLT=IJKT                                                        DRT02890
      IF (IJKEQ .AND. KOR.GT.LOR) IJKLT=IJKLT+1                         DRT02900
      IF (KOR.EQ.LOR) IJKLT=IJKLT+3                                     DRT02910
      LAD=KAD+LADD(LOFF+LOR)                                            DRT02920
      IF (IJKLT.GT.7 .OR. IJKLT.LT.1) THEN                              DRT02930
      WRITE(6,3)IJKLT,IOR,JOR,KOR,LOR,IBF,JBF,KBF,LBF,ISM,JSM,KSM,LSM   DRT02940
          STOP                                                          DRT02950
      END IF                                                            DRT02960
      IF (LAD .GT. NMAX) THEN                                           DRT02970
      WRITE(6,4) LAD                                                    DRT02980
      STOP                                                              DRT02990
      END IF                                                            DRT03000
C                                                                       DRT03010
      GOTO (310,310,310,310,305,310,307) IJKLT                          DRT03020
C                                                                       DRT03470
  305 CONTINUE                                                          DRT03480
      FAC=D2                                                            DRT03550
      VAL=BUFI(LAD+3)*FAC                                               DRT03570
      II=IDRTSO(IOR)                                                    DRT03580
      LL=IDRTSO(LOR)                                                    DRT03590
      IILL=IOFF(MAX0(II,LL))+MIN0(II,LL)                                DRT03600
      HQ(IILL)=VAL                                                      DRT03610
      GO TO 310                                                         DRT03620
C                                                                       DRT03630
  307 CONTINUE                                                          DRT03750
      VAL=BUFI(LAD+2)                                                   DRT03780
      II=IDRTSO(IOR)                                                    DRT03790
      IIII=IOFF(II+1)                                                   DRT03800
      HQ(IIII)=VAL                                                      DRT03810
C                                                                       DRT03820
  310 CONTINUE                                                          DRT03830
C*    IF(IPRNT.LE.4) GO TO 110                                          DRT03840
c*    WRITE(6,5) IOR,JOR,KOR,LOR,IJKLT,                                 DRT03850
c*   1           BUFI(LAD+1),BUFI(LAD+2),BUFI(LAD+3)                    DRT03860
C                                                                       DRT03870
  110 CONTINUE                                                          DRT03880
  120 CONTINUE                                                          DRT03890
  130 CONTINUE                                                          DRT03900
  140 CONTINUE                                                          DRT03910
C                                                                       DRT03920
      IF(IPRNT.LE.2) GO TO 202                                          DRT03930
      WRITE(6,6)                                                        DRT03940
      CALL PRINT(HQ,NUMIJ,NBF,6)                                        DRT03950
C                                                                       DRT03960
C   STORE ONE PDM ON THE MASTER FILE                                    DRT04100
  202 CONTINUE                                                          DRT04110
cyy9-07-05 call mwrit(hq,iopci)
cyy9-07-05 call mwrit(hq,iopmc)
C
C================================================
C ITAP59 has been used for pass 'one-pdm'
C================================================
      itap59 = 59
      call rfile(itap59)
      call swrit(itap59,hq,numij*2)
      call rclose(itap59,3)
C================================================
C
      RETURN                                                            DRT04180
      END                                                               DRT04190
