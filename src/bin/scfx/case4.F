C     //////////////////////////////////////////////////////////////
      SUBROUTINE CASE4(V,L,D,P,NOINBF)
C
C     ----- SUBROUTINE TO FORM FOCK MATRIX FOR CLOSED SHELL CALCULATAION
C           FROM PACKED SUPERMATRIX ON TAPE302                      ----
C
      REAL*8 T,V(1),D(1),P(1),AIJKL
C
C
      DIMENSION L(1)
C     DIMENSION IJKL(2)
C     EQUIVALENCE (AIJKL,IJKL),(IJ,IJKL(1)),(KL,IJKL(2))
C
C     WRITE (3,10) NOINBF
C  10 FORMAT (' ENTERED INTO CASE4, NOINBF=',I8)
C     DO 1 I=1,NOINBF*2-1,2
C        AIJKL=V(I)
C        T=V(I+1)
CC    WRITE (3,11) IJ,KL,T
CC 11 FORMAT (' IJ:',I10,' KL:',I10,E18.8)
C        P(IJ)=P(IJ)+T*D(KL)
C        P(KL)=P(KL)+T*D(IJ)
C   1 CONTINUE
      DO 1 I=1,NOINBF
         KL=L(2*I-1)
         IJ=L(2*I  )
         T =V(I)
         P(IJ)=P(IJ) + T*D(KL)
         P(KL)=P(KL) + T*D(IJ)
    1 CONTINUE
C     WRITE (3,12)
C  12 FORMAT (' EXITING CASE4 *******************************')
C
      RETURN
      END
