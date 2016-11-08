C
C----------------------------------------------------------------------
C
      SUBROUTINE OFFSET(INO,INV,INTRO,INTRV,ZLX,UOFF,NBF,NO,NV,NOV,
     .                  NTRO,NTRV,N2OV,JOUT,ZLXO,UOFFV,IOFF,NTR)
      IMPLICIT INTEGER(A-Z)
      DIMENSION INO(NOV),INV(NOV),INTRO(N2OV),INTRV(NO*NV),ZLX(NV,NV),
     .          UOFF(NO,NO,2),ZLXO(NO,NO),UOFFV(NV,NV,2),IOFF(NTR)
C
      IOFF(1)=0
      DO 5 I=1,NTR-1
        IOFF(I+1)=IOFF(I)+I
    5 CONTINUE
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
      DO 70 I=1,NO
         DO 80 J=1,NO
            UOFF(I,J,1)=(I-1)*NO+J
            UOFF(I,J,2)=(J-1)*NO+I
   80    CONTINUE
   70 CONTINUE
      DO 41 I=1,NO
         DO 51 J=1,I
            ZLXO(I,J)=1
   51    CONTINUE
         DO 61 J=I+1,NO
            ZLXO(I,J)=2
   61    CONTINUE
   41 CONTINUE
      DO 71 A=1,NV
         DO 81 B=1,NV
            UOFFV(A,B,1)=(A-1)*NV+B
            UOFFV(A,B,2)=(B-1)*NV+A
   81    CONTINUE
   71 CONTINUE
C
      RETURN
      END
