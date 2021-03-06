C///////////////////////////////////////////////////////////////////////
      SUBROUTINE YIN2(NC,NS,Y,R,NB,ISCR)
      REAL*8 Y(NC,NC,NC)
      INTEGER R,NB,LR
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
    1 FORMAT(1X,'YR BLOCK',I4,' FOR R=',I4,' READ ON FILE ',I4,
     $    ' AT LOCATION ',I8)
    2 FORMAT(1X,'R = 0 ENCOUNTERED IN SUBROUTINE YIN2.')
      IF(R.LT.0) THEN
         LR=-2*(R+1)*NC+NB
      ELSE IF(R.GT.0) THEN
         LR=(2*R-1)*NC+NB
      ELSE
         WRITE(ICHECK,2)
         STOP
      END IF
      NLEN=NC**3
      MMM=(LR-1)*NLEN
      LLL=INTOWP(MMM)+1
      CALL WREADW(ISCR,Y,INTOWP(NLEN),LLL,KKK)
      IF(LPRT(4,NPRT).GE.1) RETURN
      WRITE(ICHECK,1) NB,R,ISCR,LLL
      RETURN
      END
