      SUBROUTINE MTXMPY(A,NAD,B,NBD,C,NCD,D,NDD,N,IMPY)                 NOR08170
      IMPLICIT REAL*8 (A-H,O-Z)                                         NOR08180
      DIMENSION A(NAD,NAD),B(NBD,NBD),C(NCD,NCD),D(NDD,NDD)             NOR08190
      DATA A00 / 0.0D+00 /                                              NOR08200
C                                                                       NOR08210
      DO 101 I=1,NCD                                                    NOR08220
      DO 101 J=1,NCD                                                    NOR08230
      C(I,J)=A00                                                        NOR08240
  101 continue
      GO TO (201,202,203,204,205,206),IMPY                              NOR08250
c
C     C=A*B                                                             NOR08260
  201 continue
      DO 102 I=1,N                                                      NOR08270
      DO 102 J=1,N                                                      NOR08280
      DO 102 K=1,N                                                      NOR08290
      C(I,J)=C(I,J)+A(I,K)*B(K,J)                                       NOR08300
  102 continue
      RETURN                                                            NOR08310
c
C     C=TA*B                                                            NOR08320
  202 continue
      DO 103 I=1,N                                                      NOR08330
      DO 103 J=1,N                                                      NOR08340
      DO 103 K=1,N                                                      NOR08350
      C(I,J)=C(I,J)+A(K,I)*B(K,J)                                       NOR08360
  103 continue
      RETURN                                                            NOR08370
c
C     C=A*TB                                                            NOR08380
  203 continue
      DO 104 I=1,N                                                      NOR08390
      DO 104 J=1,N                                                      NOR08400
      DO 104 K=1,N                                                      NOR08410
      C(I,J)=C(I,J)+A(I,K)*B(J,K)                                       NOR08420
  104 continue
      RETURN                                                            NOR08430
c
C     C=TA*TB                                                           NOR08440
  204 continue
      DO 105 I=1,N                                                      NOR08450
      DO 105 J=1,N                                                      NOR08460
      DO 105 K=1,N                                                      NOR08470
      C(I,J)=C(I,J)+A(K,I)*B(J,K)                                       NOR08480
  105 continue
      RETURN                                                            NOR08490
c
C     C=TA*B*A                                                          NOR08500
  205 continue
      DO 106 I=1,N                                                      NOR08510
      DO 106 J=1,N                                                      NOR08520
      D(I,J)=A00                                                        NOR08530
      DO 106 K=1,N                                                      NOR08540
      D(I,J)=D(I,J)+A(K,I)*B(K,J)                                       NOR08550
  106 continue
      DO 107 I=1,N                                                      NOR08560
      DO 107 J=1,N                                                      NOR08570
      DO 107 K=1,N                                                      NOR08580
      C(I,J)=C(I,J)+D(I,K)*A(K,J)                                       NOR08590
  107 continue
      RETURN                                                            NOR08600
c
C     C=A*B*TA                                                          NOR08610
  206 continue
      DO 108 I=1,N                                                      NOR08620
      DO 108 J=1,N                                                      NOR08630
      D(I,J)=A00                                                        NOR08640
      DO 108 K=1,N                                                      NOR08650
      D(I,J)=D(I,J)+A(I,K)*B(K,J)                                       NOR08660
  108 continue
      DO 109 I=1,N                                                      NOR08670
      DO 109 J=1,N                                                      NOR08680
      DO 109 K=1,N                                                      NOR08690
      C(I,J)=C(I,J)+D(I,K)*A(J,K)                                       NOR08700
  109 continue
c
      RETURN                                                            NOR08710
      END                                                               NOR08720
