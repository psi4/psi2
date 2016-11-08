      SUBROUTINE ORDER(VAL,I,J,K,L,IBL,LBLO,BUFO,NTAPE)                 DRT04200
      IMPLICIT REAL*8 (A-H,O-Z)                                         DRT04210
      DIMENSION LBLO(MAXBF4),BUFO(MAXBF2)                               DRT04220
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          DRT04230
      COMMON/ORDSO/IORSO(256),IDRTSO(256)                               DRT04240
      COMMON/SIGNS/IOFF(256),IPRNT                                      DRT04250
    1 FORMAT(2X,4I5,5X,4I5,5X,F20.10)                                   DRT04260
C                                                                       DRT04270
      IT=IDRTSO(I)                                                      DRT04280
      JT=IDRTSO(J)                                                      DRT04290
      KT=IDRTSO(K)                                                      DRT04300
      LT=IDRTSO(L)                                                      DRT04310
      IJT=IOFF(MAX0(IT,JT))+MIN0(IT,JT)                                 DRT04320
      KLT=IOFF(MAX0(KT,LT))+MIN0(KT,LT)                                 DRT04330
      IF(IJT.LT.KLT) GO TO 201                                          DRT04340
      II=MAX0(IT,JT)                                                    DRT04350
      JJ=MIN0(IT,JT)                                                    DRT04360
      KK=MAX0(KT,LT)                                                    DRT04370
      LL=MIN0(KT,LT)                                                    DRT04380
      GO TO 202                                                         DRT04390
  201 CONTINUE                                                          DRT04400
      II=MAX0(KT,LT)                                                    DRT04410
      JJ=MIN0(KT,LT)                                                    DRT04420
      KK=MAX0(IT,JT)                                                    DRT04430
      LL=MIN0(IT,JT)                                                    DRT04440
  202 CONTINUE                                                          DRT04450
      CALL PACK(II,JJ,KK,LL,IX,JA,JAA)                                  DRT04460
      IBL=IBL+1                                                         DRT04470
      IF(IPRNT.LE.3) GO TO 203                                          DRT04480
      WRITE(6,1) I,J,K,L,II,JJ,KK,LL,VAL                                DRT04490
  203 CONTINUE                                                          DRT04500
      LBLO(IBL+IBL-1)=JA                                                DRT04510
      LBLO(IBL+IBL)=JAA                                                 DRT04520
      BUFO(IBL+MAXBUF)=VAL                                              DRT04530
      IF(IBL.LT.MAXBUF) GO TO 205                                       DRT04540
      IBL=0                                                             DRT04550
      CALL SWRIT(NTAPE,LBLO,MAXBF4)                                     DRT04560
C                                                                       DRT04570
  205 CONTINUE                                                          DRT04580
      RETURN                                                            DRT04590
      END                                                               DRT04600
