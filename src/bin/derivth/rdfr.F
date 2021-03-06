      SUBROUTINE RDFR(T,MINK,MAXK,LKT,MINL,MAXL,LLT,NTR)
      IMPLICIT REAL*8 (A-H,O-Z)
CTPH  DIMENSION V(10),T(10,10)
      DIMENSION V(35),T(35,35)
CTPH  COMMON/GAUSC/TT(36),PTR(3,3,8),DTR(6,6,8)
      COMMON/GAUSC/TT(225),PTR(3,3,8),DTR(6,6,8),FTR(10,10,8),
     1             GTR(15,15,8)
      COMMON/PARA1/NSHELL,NT,MXT,NISO
      COMMON/TOLER/DLIM
      DATA ZERO / 0.0D+00 /
C
C   RIGHT MULTIPLY T BY R
C   RESULT BACK IN T
C
      GO TO (204,203,202,201,200),LLT
C
C     G SHELL
C
  200 CONTINUE
      DO 4 K=MINK,MAXK
      DO 2 L=21,35
      DUM=ZERO
      DO 1 N=21,35
    1 DUM=DUM+T(K,N)*GTR(N-20,L-20,NTR)
    2 V(L)=DUM
      DO 3 L=21,35
    3 T(K,L)=V(L)
    4 CONTINUE
      GO TO 204
C
C     F SHELL
C
  201 CONTINUE
      DO 104 K=MINK,MAXK
      DO 102 L=11,20
      DUM=ZERO
      DO 101 N=11,20
  101 DUM=DUM+T(K,N)*FTR(N-10,L-10,NTR)
  102 V(L)=DUM
      DO 103 L=11,20
  103 T(K,L)=V(L)
  104 CONTINUE
      GO TO 204
C
C     D SHELL
C
  202 CONTINUE
      DO 108 K=MINK,MAXK
      DO 106 L=5,10
      DUM=ZERO
      DO 105 N=5,10
  105 DUM=DUM+T(K,N)*DTR(N-4,L-4,NTR)
  106 V(L)=DUM
      DO 107 L=5,10
  107 T(K,L)=V(L)
  108 CONTINUE
      GO TO 204
C
C     P SHELL
C
  203 CONTINUE
      DO 112 K=MINK,MAXK
      DO 110 L=2,4
      DUM=ZERO
      DO 109 N=2,4
  109 DUM=DUM+T(K,N)*PTR(N-1,L-1,NTR)
  110 V(L)=DUM
      DO 111 L=2,4
  111 T(K,L)=V(L)
  112 CONTINUE
C
C   LEFT MULTIPLY T BY R
C   RESULT BACK IN T
C
  204 CONTINUE
      GO TO (209,208,207,206,205), LKT
C
C     G SHELL
C
  205 CONTINUE
      DO 216 L=MINL,MAXL
      DO 214 K=21,35
      DUM=ZERO
      DO 213 N=21,35
      DUM=DUM+T(N,L)*GTR(N-20,K-20,NTR)
  213 CONTINUE
      V(K)=DUM
  214 CONTINUE
      DO 215 K=21,35
      T(K,L)=V(K)
  215 CONTINUE
  216 CONTINUE
      GO TO 209
C
C     F SHELL
C
  206 CONTINUE
      DO 116 L=MINL,MAXL
      DO 114 K=11,20
      DUM=ZERO
      DO 113 N=11,20
      DUM=DUM+T(N,L)*FTR(N-10,K-10,NTR)
  113 CONTINUE
      V(K)=DUM
  114 CONTINUE
      DO 115 K=11,20
      T(K,L)=V(K)
  115 CONTINUE
  116 CONTINUE
      GO TO 209
C
C     D SHELL
C
  207 CONTINUE
      DO 120 L=MINL,MAXL
      DO 118 K=5,10
      DUM=ZERO
      DO 117 N=5,10
      DUM=DUM+T(N,L)*DTR(N-4,K-4,NTR)
  117 CONTINUE
      V(K)=DUM
  118 CONTINUE
      DO 119 K=5,10
      T(K,L)=V(K)
  119 CONTINUE
  120 CONTINUE
      GO TO 209
C
C     P SHELL
C
  208 CONTINUE
      DO 124 L=MINL,MAXL
      DO 122 K=2,4
      DUM=ZERO
      DO 121 N=2,4
      DUM=DUM+T(N,L)*PTR(N-1,K-1,NTR)
  121 CONTINUE
      V(K)=DUM
  122 CONTINUE
      DO 123 K=2,4
      T(K,L)=V(K)
  123 CONTINUE
  124 CONTINUE
C
  209 CONTINUE
      RETURN
      END
