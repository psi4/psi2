C/////////////////////////////////////////////////////////////////
      SUBROUTINE OVERLP(SUQ,SUQOLD,S,NQ,NMAT,KKQ,SQ1)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INTEGER NQ,NMAT,KKQ
      DOUBLE PRECISION SUQ(NQ),SUQOLD(NQ),S(NMAT)
      DOUBLE PRECISION SQ1,SQ2,VAL1,VAL2
C
      A0=0.0D+00
      SQ1=A0
      DO 100 I=1,NQ
        DO 150 J=1,NQ
          I1=MAX0(I,J)
          J1=MIN0(I,J)
          IJ1=I1*(I1-1)/2+J1+KKQ
          VAL1=SUQ(I)*SUQOLD(J)
          SQ1=SQ1+VAL1*S(IJ1)
 150    CONTINUE
 100  CONTINUE
      RETURN
      END
