C
C----------------------------------------------------------------------
C
      SUBROUTINE OFFSET(INO,INV,INTRO,INTRV,ZLX,IOFF,
     .                  NBF,NO,NV,NOV,NTRO,NTRV,N2OV,NTR,JOUT)
      IMPLICIT INTEGER(A-Z)
      DIMENSION INO(NOV),INV(NOV),INTRO(N2OV),INTRV(NO*NV),ZLX(NV,NV),
     .          IOFF(NTR)
C
      IOFF(1)=0
      DO 73 I=1,NTR-1
         IOFF(I+1)=IOFF(I)+I
   73 CONTINUE
C
      DO 10 I=1,NOV
         INO(I)=(I-1)*NO
   10 CONTINUE
      DO 11 I=1,NOV
         INV(I)=(I-1)*NV
   11 CONTINUE
      DO 20 I=1,N2OV
         INTRO(I)=(I-1)*NTRO
   20 CONTINUE
      DO 30 I=1,NO*NV
         INTRV(I)=(I-1)*NTRV
   30 CONTINUE
      DO 40 A=1,NV
         DO 50 B=1,A
            ZLX(A,B)=1
   50    CONTINUE
         DO 60 B=A+1,NV
            ZLX(A,B)=2
   60    CONTINUE
   40 CONTINUE
C
      RETURN
      END
