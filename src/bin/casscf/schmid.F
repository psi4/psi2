      SUBROUTINE SCHMID(U,E,V,NDIM,N)                                   SMM23580
      IMPLICIT REAL*8 (A-H,O-Z)                                         SMM23590
      DIMENSION U(NDIM,NDIM),E(NDIM,NDIM),V(NDIM)                       SMM23600
      DATA A00 / 0.0D+00 /                                              SMM23610
C                                                                       SMM23620
C   N=1                                                                 SMM23630
      DO 101 I=1,N                                                      SMM23640
  101 V(I)=U(I,1)                                                       SMM23650
      CALL NORMLZ(V,NDIM,N,1)                                           SMM23660
      DO 102 I=1,N                                                      SMM23670
  102 E(I,1)=V(I)                                                       SMM23680
C                                                                       SMM23690
C   N=>2                                                                SMM23700
      DO 110 II=2,N                                                     SMM23710
      DO 103 I=1,N                                                      SMM23720
  103 V(I)=U(I,II)                                                      SMM23730
      IR=II-1                                                           SMM23740
      DO 106 K=1,IR                                                     SMM23750
      EU=A00                                                            SMM23760
      DO 104 I=1,N                                                      SMM23770
  104 EU=EU+E(I,K)*U(I,II)                                              SMM23780
      DO 105 I=1,N                                                      SMM23790
  105 V(I)=V(I)-EU*E(I,K)                                               SMM23800
  106 CONTINUE                                                          SMM23810
      CALL NORMLZ(V,NDIM,N,II)                                          SMM23820
      DO 107 I=1,N                                                      SMM23830
  107 E(I,II)=V(I)                                                      SMM23840
  110 CONTINUE                                                          SMM23850
      DO 111 I=1,N                                                      SMM23860
      DO 111 J=1,N                                                      SMM23870
      U(I,J)=E(I,J)                                                     SMM23880
  111 CONTINUE                                                          SMM23890
C                                                                       SMM23900
      RETURN                                                            SMM23910
      END                                                               SMM23920
