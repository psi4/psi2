      SUBROUTINE UNITRY(ESO,NIJ,RMO,R,U,W,VEC,NORD)                     SMM22380
      IMPLICIT REAL*8 (A-H,O-Z)                                         SMM22390
      DIMENSION NORD(NBASIS),VEC(NBASIS)                                SMM22400
      DIMENSION ESO(NBASIS,NBASIS),NIJ(NTRI,2),RMO(nind)
      DIMENSION R(NBASIS,NBASIS),U(NBASIS,NBASIS),W(NBASIS,NBASIS)      SMM22420
      COMMON/BASIS/NBASIS,NTRI,NTOT,NSYM
      COMMON/SIGNS/IOFF(32897),IPRNT  
      COMMON/CI102/NIND,NINTRI
      DATA A00,HALF,ONE / 0.0D+00 , 0.5D+00 , 1.0D+00 /                 SMM22470
    1 FORMAT(//,2X,' R MATRIX'/)                                        SMM22480
    2 FORMAT(//,2X,' U MATRIX'/)                                        SMM22490
    3 FORMAT(//,2X,' ORTHOGONALIZED U MATRIX'/)                         SMM22500
    4 FORMAT(//,2X,' UNITARITY OF U MATRIX'/)                           SMM22510
    5 FORMAT(//,2X,' NORD = ',(30I4))                                   SMM22520
    6 FORMAT(//,2X,' XXXXXX WARNING IN NORD XXXXXX',I4,' XXXXXX'//)     SMM22530
    7 FORMAT(//,2X,' XXXXXX PROGRAM ERROR IN NORD-----STOP ! XXXXXX'//) SMM22540
    8 FORMAT(//,2X,' REARRANGED U MATRIX'/)                             SMM22550
    9 FORMAT(//,2X,' ROTATED EIGENVECTOR (SO-MO) MATRIX'/)              SMM22560
C                                                                       SMM22570
C   FORM A UNITARY MATRIX FROM THE R MATRIX                             SMM22580
C                                                                       SMM22590
      CALL ZERO(R,NBASIS*NBASIS)                                        SMM22600
      DO 101 II=1,NIND                                                  SMM22610
      I=NIJ(II,1)                                                       SMM22620
      J=NIJ(II,2)                                                       SMM22630
      R(I,J)=RMO(II)                                             
      R(J,I)=-RMO(II)                                        
  101 CONTINUE                                                          SMM22660
      IF(IPRNT.LE.3) GO TO 201                                          SMM22670
      WRITE(6,1)                                                        SMM22680
      CALL MATOUT(R,NBASIS,NBASIS,NBASIS,NBASIS,6)                      SMM22690
C                                                                       SMM22700
C   FORM THE UNITARY MATRIX                                             SMM22710
  201 CONTINUE                                                          SMM22720
      DO 105 I=1,NBASIS                                                 SMM22730
      DO 104 J=1,NBASIS                                                 SMM22740
      VALU=A00                                                          SMM22750
      DO 103 K=1,NBASIS                                                 SMM22760
      VALU=VALU+R(I,K)*R(K,J)                                           SMM22770
  103 CONTINUE                                                          SMM22780
      U(I,J)=R(I,J)+VALU*HALF                                           SMM22790
  104 CONTINUE                                                          SMM22800
      U(I,I)=U(I,I)+ONE                                                 SMM22810
  105 CONTINUE                                                          SMM22820
      IF(IPRNT.LE.3) GO TO 202                                          SMM22830
      WRITE(6,2)                                                        SMM22840
      CALL MATOUT(U,NBASIS,NBASIS,NBASIS,NBASIS,6)                      SMM22850
C                                                                       SMM22860
C   SCHMIDT ORTHOGONALIZATION                                           SMM22870
  202 CONTINUE                                                          SMM22880
      CALL SCHMID(U,W,VEC,NBASIS,NBASIS)                                SMM22890
      IF(IPRNT.LE.3) GO TO 203                                          SMM22900
      WRITE(6,3)                                                        SMM22910
      CALL MATOUT(U,NBASIS,NBASIS,NBASIS,NBASIS,6)                      SMM22920
C                                                                       SMM22930
C   CHECK UNITARITY OF U MATRIX                                         SMM22940
  203 CONTINUE                                                          SMM22950
      DO 107 I=1,NBASIS                                                 SMM22960
      DO 107 J=1,NBASIS                                                 SMM22970
      VALU=A00                                                          SMM22980
      DO 106 K=1,NBASIS                                                 SMM22990
      VALU=VALU+U(K,I)*U(K,J)                                           SMM23000
  106 CONTINUE                                                          SMM23010
      W(I,J)=VALU                                                       SMM23020
  107 CONTINUE                                                          SMM23030
      IF(IPRNT.LE.3) GO TO 204                                          SMM23040
      WRITE(6,4)                                                        SMM23050
      CALL MATOUT(W,NBASIS,NBASIS,NBASIS,NBASIS,6)                      SMM23060
C                                                                       SMM23070
C   CHECK ORDERING OF ORBITALS                                          SMM23080
  204 CONTINUE                                                          SMM23090
      DO 111 I=1,NBASIS                                                 SMM23100
      NORD(I)=1                                                         SMM23110
      EMAX=DABS(U(1,I))                                                 SMM23120
      DO 110 J=2,NBASIS                                                 SMM23130
      IF(DABS(U(J,I)).LT.EMAX) GO TO 110                                SMM23140
      EMAX=DABS(U(J,I))                                                 SMM23150
      NORD(I)=J                                                         SMM23160
  110 CONTINUE                                                          SMM23170
  111 CONTINUE                                                          SMM23180
      IF(IPRNT.GE.3)                                                    SMM23190
     *WRITE(3,5) (NORD(I),I=1,NBASIS)                                   SMM23200
      DO 116 I=1,NBASIS                                                 SMM23210
      DO 113 J=1,NBASIS                                                 SMM23220
      IF(NORD(J).EQ.I) GO TO 116                                        SMM23230
  113 CONTINUE                                                          SMM23240
      WRITE(6,6) I                                                      SMM23250
      DO 114 K=1,NBASIS-1                                               SMM23260
      DO 114 L=K+1,NBASIS                                               SMM23270
      IF(NORD(K).EQ.NORD(L)) GO TO 115                                  SMM23280
  114 CONTINUE                                                          SMM23290
      WRITE(6,7)                                                        SMM23300
      STOP                                                              SMM23310
  115 NORD(K)=I                                                         SMM23320
  116 CONTINUE                                                          SMM23330
      IF(IPRNT.GE.3)                                                    SMM23340
     *WRITE(3,5) (NORD(I),I=1,NBASIS)                                   SMM23350
C                                                                       SMM23360
C   FORM A NEW EIGENVECTOR MATRIX                                       SMM23370
  206 CONTINUE                                                          SMM23380
      DO 121 I=1,NBASIS                                                 SMM23390
      DO 121 J=1,NBASIS                                                 SMM23400
      VALU=A00                                                          SMM23410
      DO 120 K=1,NBASIS                                                 SMM23420
      VALU=VALU+ESO(I,K)*U(K,J)                                         SMM23430
  120 CONTINUE                                                          SMM23440
      W(I,J)=VALU                                                       SMM23450
  121 CONTINUE                                                          SMM23460
      DO 122 I=1,NBASIS                                                 SMM23470
      DO 122 J=1,NBASIS                                                 SMM23480
      ESO(I,J)=W(I,J)                                                   SMM23490
  122 CONTINUE                                                          SMM23500
      IF(IPRNT.LE.3) GO TO 207                                          SMM23510
      WRITE(6,9)                                                        SMM23520
      CALL MATOUT(ESO,NBASIS,NBASIS,NBASIS,NBASIS,6)                    SMM23530
C                                                                       SMM23540
  207 CONTINUE                                                          SMM23550
      RETURN                                                            SMM23560
      END                                                               SMM23570
