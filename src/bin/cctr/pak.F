C
C                           *****************
C---------------------------*****  PAK  *****--------------------------
C                           *****************
C
      SUBROUTINE PAK
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PAKQ/ IX(2),I5(2)
C
      INTEGER*2 I5
C     LOGICAL*4 TRANS(2)
C     EQUIVALENCE (IX(1),TRANS(1))
C
      DO 1 I=1,2
         IF (IX(I) .GT. 65535) STOP ' PACKING PROBLEM '
         I5(I)=IX(I)
    1 CONTINUE
C     I5=ITEMP
C
      RETURN
C
C----UNPAK-----------------------------------------------UNPAK
C
      ENTRY UNPAK
C
C     ITEMP=I5
      DO 2 I=1,2
         IX(I)=I5(I)
    2 CONTINUE
C
      RETURN
      END
