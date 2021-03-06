C     /////////////////////////////////////////////////////////////////
      SUBROUTINE HIJS8(NAD,K1,K2,K3,K4,XA,H11,H21,H31,H41,
     $               H22,H32,H42,H33,H43,H44)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),E2(3),E4(3),Q1(3),Q2(3),Q3(3)
      DIMENSION E22(3,3),E44(3,3)
      DIMENSION H11(3,3),H21(3,3),H31(3,3),H41(3,3),H22(3,3),H32(3,3)
      DIMENSION H42(3,3),H33(3,3),H43(3,3),H44(3,3)
      DIMENSION Q11(3,3),Q21(3,3),Q31(3,3),Q33(3,3),Q22(3,3),Q32(3,3)
      DIMENSION Q111(3,3,3),Q222(3,3,3),Q333(3,3,3),Q444(3,3,3)
      DIMENSION Q112(3,3,3),Q223(3,3,3),Q331(3,3,3)
      DIMENSION Q221(3,3,3),Q332(3,3,3),Q113(3,3,3),Q123(3,3,3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      CALL VECT1(NAD,K2,K3,E2,XA,T32)
      CALL VECT1(NAD,K4,K3,E4,XA,T34)
      CALL VECT2(NAD,K1,K2,K3,Q1,Q2,Q3,XA,T123)
      CALL SCAPRO(E4,Q3,T)
      W=-T32*T
      CALL HIJS1(NAD,K4,K3,XA,E44)
      CALL HIJS2(NAD,K1,K2,K3,XA,Q11,Q21,Q31,Q22,Q32,Q33)
      CALL HIJS1(NAD,K2,K3,XA,E22)
      CALL HIJKS1(NAD,K4,K3,XA,Q444)
      DO 1 J=1,3
      DO 1 K=1,3
      H44(J,K)=ZERO
      DO 1 I=1,3
  1   H44(J,K)=H44(J,K)-T32*Q444(I,J,K)*Q3(I)
      CALL HIJKS2(NAD,K1,K2,K3,XA,Q111,Q112,Q113,Q123,Q221,
     $    Q222,Q223,Q331,Q332,Q333)
      DO 5 K=1,3
      DO 5 J=1,3
      H41(J,K)=ZERO
      H42(J,K)=ZERO
      H11(J,K)=ZERO
      H21(J,K)=ZERO
      H22(J,K)=W*E22(J,K)/T32
      DO 5 I=1,3
      H11(J,K)=H11(J,K)-T32*E4(I)*Q113(J,K,I)
      H21(J,K)=H21(J,K)-E4(I)*( E2(J)*Q31(I,K)+T32*Q123(K,J,I) )
      H22(J,K)=H22(J,K)-E4(I)*( E2(J)*Q32(I,K)+E2(K)*Q32(I,J)+
     $      T32*Q223(J,K,I) )
      H41(J,K)=H41(J,K)-T32*E44(I,J)*Q31(I,K)
  5   H42(J,K)=H42(J,K)-E44(I,J)*( T32*Q32(I,K)+E2(K)*Q3(I) )
      DO 3 J=1,3
      DO 3 K=1,3
      H31(J,K)=-H11(J,K)-H21(J,K)-H41(J,K)
      H32(J,K)=-H21(K,J)-H22(J,K)-H42(J,K)
C     H32(J,K)=ZERO
C     DO 2 I=1,3
C     H32(J,K)=H32(J,K)-E2(K)*( -Q3(I)*E44(I,J)+E4(I)*Q33(I,J) )
C     H32(J,K)=H32(J,K)+E4(I)*(  Q3(I)*E22(K,J)+E2(J)*Q32(I,K) )
C 2   H32(J,K)=H32(J,K)-T32*( E4(I)*Q332(I,J,K)-E44(J,I)*Q32(I,K) )
  3   H43(J,K)=-H41(J,K)-H42(J,K)-H44(J,K)
      DO 4 J=1,3
      DO 4 K=1,3
  4   H33(J,K)=-H31(J,K)-H32(J,K)-H43(K,J)
      RETURN
      END
