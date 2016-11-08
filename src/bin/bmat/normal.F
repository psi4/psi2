C     ////////////////////////////////////////////////
      SUBROUTINE NORMAL (U,V,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION U(3),V(3),W(3)
C..... W WIRD EIN SENKRECHT AUF DIE EBENE (U,V) STEHENDER EINHEITSVECTOR
      W(1)=U(2)*V(3)-U(3)*V(2)
      W(2)=U(3)*V(1)-U(1)*V(3)
      W(3) = U(1)*V(2)-U(2)*V(1)
      CALL NOM  (W)
      RETURN
      END
