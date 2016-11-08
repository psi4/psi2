      SUBROUTINE GETIT(KADD,LADD,IJADD,IJGRP,HQ,BUFI,LBLO,BUFO,         DRT02190
     1                 ITAPE,NTAPE)                                     DRT02200
      IMPLICIT REAL*8 (A-H,O-Z)                                         DRT02210
      LOGICAL MOINT                                                     DRT02220
      DIMENSION KADD(ISYMRB),LADD(ISYMRB),IJADD(NUMIJ),IJGRP(NUMIJ)     DRT02230
      DIMENSION HQ(NUMIJ),BUFI(NMAX),LBLO(MAXBF4),BUFO(MAXBF2)          DRT02240
      DIMENSION DUMS(56)                                                DRT02250
      LOGICAL IJEQ,IJKEQ,JKEQ                                           DRT02260
      COMMON/DIMES/ISYMRB,NUMIJ,NBF,NORBS,NMAX                          DRT02270
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          DRT02290
      COMMON/ORDCI/IDRSYM(256),IORDRT(256)                              DRT02300
      COMMON/ORDSO/IORSO(256),IDRTSO(256)                               DRT02310
      common/posit/ioeci,ioemc,iopci,iopmc
      COMMON/SIGNS/IOFF(256),IPRNT                                      DRT02320
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
C                                                                       DRT02460
      MOINT=.FALSE.                                                     DRT02470
      CALL RFILE(ITAPE)                                                 DRT02480
      CALL RFILE(NTAPE)                                                 DRT02490
      IF(ITAPE.NE.52) GO TO 201                                         DRT02500
      MOINT=.TRUE.                                                      DRT02510
      CALL SREAD(ITAPE,DUMS,112)                                        DRT02520
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
      IF(IPRNT.GT.4)                                                    DRT02640
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
      GOTO (301,302,303,304,305,306,307) IJKLT                          DRT03020
C                                                                       DRT03030
  301 CONTINUE                                                          DRT03040
      FAC=D8                                                            DRT03050
      IF(MOINT) FAC=ONE                                                 DRT03060
      VAL=BUFI(LAD+1)*FAC                                               DRT03070
      CALL ORDER(VAL,IOR,KOR,JOR,LOR,IBL,LBLO,BUFO,NTAPE)               DRT03080
      VAL=BUFI(LAD+2)*FAC                                               DRT03090
      CALL ORDER(VAL,IOR,JOR,KOR,LOR,IBL,LBLO,BUFO,NTAPE)               DRT03100
      VAL=BUFI(LAD+3)*FAC                                               DRT03110
      CALL ORDER(VAL,IOR,LOR,JOR,KOR,IBL,LBLO,BUFO,NTAPE)               DRT03120
      GO TO 310                                                         DRT03130
C                                                                       DRT03140
  302 CONTINUE                                                          DRT03150
      FAC=D8                                                            DRT03160
      IF(MOINT) FAC=ONE                                                 DRT03170
      VAL=BUFI(LAD+1)*FAC                                               DRT03180
      CALL ORDER(VAL,IOR,JOR,JOR,LOR,IBL,LBLO,BUFO,NTAPE)               DRT03190
      FAC=D4                                                            DRT03200
      IF(MOINT) FAC=ONE                                                 DRT03210
      VAL=BUFI(LAD+2)*FAC                                               DRT03220
      CALL ORDER(VAL,IOR,LOR,JOR,JOR,IBL,LBLO,BUFO,NTAPE)               DRT03230
      GO TO 310                                                         DRT03240
C                                                                       DRT03250
  303 CONTINUE                                                          DRT03260
      FAC=D8                                                            DRT03270
      IF(MOINT) FAC=ONE                                                 DRT03280
      VAL=BUFI(LAD+1)*FAC                                               DRT03290
      CALL ORDER(VAL,IOR,KOR,IOR,LOR,IBL,LBLO,BUFO,NTAPE)               DRT03300
      FAC=D4                                                            DRT03310
      IF(MOINT) FAC=ONE                                                 DRT03320
      VAL=BUFI(LAD+2)*FAC                                               DRT03330
      CALL ORDER(VAL,IOR,IOR,KOR,LOR,IBL,LBLO,BUFO,NTAPE)               DRT03340
      GO TO 310                                                         DRT03350
C                                                                       DRT03360
  304 CONTINUE                                                          DRT03370
      FAC=D8                                                            DRT03380
      IF(MOINT) FAC=ONE                                                 DRT03390
      VAL=BUFI(LAD+1)*FAC                                               DRT03400
      CALL ORDER(VAL,IOR,LOR,JOR,LOR,IBL,LBLO,BUFO,NTAPE)               DRT03410
      FAC=D4                                                            DRT03420
      IF(MOINT) FAC=ONE                                                 DRT03430
      VAL=BUFI(LAD+2)*FAC                                               DRT03440
      CALL ORDER(VAL,IOR,JOR,LOR,LOR,IBL,LBLO,BUFO,NTAPE)               DRT03450
      GO TO 310                                                         DRT03460
C                                                                       DRT03470
  305 CONTINUE                                                          DRT03480
      FAC=D4                                                            DRT03490
      IF(MOINT) FAC=ONE                                                 DRT03500
      VAL=BUFI(LAD+1)*FAC                                               DRT03510
      CALL ORDER(VAL,IOR,IOR,IOR,LOR,IBL,LBLO,BUFO,NTAPE)               DRT03520
      VAL=BUFI(LAD+2)*FAC                                               DRT03530
      CALL ORDER(VAL,IOR,LOR,LOR,LOR,IBL,LBLO,BUFO,NTAPE)               DRT03540
      FAC=D2                                                            DRT03550
      IF(MOINT) FAC=ONE                                                 DRT03560
      VAL=BUFI(LAD+3)*FAC                                               DRT03570
      II=IDRTSO(IOR)                                                    DRT03580
      LL=IDRTSO(LOR)                                                    DRT03590
      IILL=IOFF(MAX0(II,LL))+MIN0(II,LL)                                DRT03600
      HQ(IILL)=VAL                                                      DRT03610
      GO TO 310                                                         DRT03620
C                                                                       DRT03630
  306 CONTINUE                                                          DRT03640
      FAC=D4                                                            DRT03650
      IF(MOINT) FAC=ONE                                                 DRT03660
      VAL=BUFI(LAD+1)*FAC                                               DRT03670
      CALL ORDER(VAL,IOR,LOR,IOR,LOR,IBL,LBLO,BUFO,NTAPE)               DRT03680
      FAC=D2                                                            DRT03690
      IF(MOINT) FAC=ONE                                                 DRT03700
      VAL=BUFI(LAD+2)*FAC                                               DRT03710
      CALL ORDER(VAL,IOR,IOR,LOR,LOR,IBL,LBLO,BUFO,NTAPE)               DRT03720
      GO TO 310                                                         DRT03730
C                                                                       DRT03740
  307 CONTINUE                                                          DRT03750
      VAL=BUFI(LAD+1)                                                   DRT03760
      CALL ORDER(VAL,IOR,IOR,IOR,IOR,IBL,LBLO,BUFO,NTAPE)               DRT03770
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
C   WRITE THE LAST BUFFER                                               DRT03970
  202 CONTINUE                                                          DRT03980
      IBL=IBL+1                                                         DRT03990
      LBLO(IBL+IBL-1)=0                                                 DRT04000
      LBLO(IBL+IBL)=0                                                   DRT04010
      BUFO(IBL+MAXBUF)=A00                                              DRT04020
      CALL SWRIT(NTAPE,LBLO,MAXBF4)                                     DRT04030
C                                                                       DRT04040
C   STORE ONE ELECTRON MO INTEGRALS ON THE MASTER FILE
      if(.not.moint) go to 203
      call mwrit(hq,ioeci)
      call mwrit(hq,ioemc)
      go to 205
C                                                                       DRT04090
C   STORE ONE PDM ON THE MASTER FILE                                    DRT04100
  203 CONTINUE                                                          DRT04110
      call mwrit(hq,iopci)
      call mwrit(hq,iopmc)
C                                                                       DRT04150
  205 continue
      CALL RCLOSE(ITAPE,3)                                              DRT04160
      CALL RCLOSE(NTAPE,3)                                              DRT04170
      RETURN                                                            DRT04180
      END                                                               DRT04190
