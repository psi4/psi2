      SUBROUTINE PRIME(BUFFER)
C
C START THE ROUND ROBIN BY FILLING CORE
C
      IMPLICIT INTEGER (A-Z)
C
C
      COMMON /T102/ NOUT,HALF,THIRD,MAXNO,OUTBUF,KINADR
      COMMON /RROB/ BTAG,POINTR,TAGMAX,RDCNT
      DIMENSION BUFFER(1)
C     REWIND (NOUT)
      CALL SREW(NOUT)
      POINTR=1
      BTAG=1
      LIM=MIN(TAGMAX,MAXNO)
      RDCNT=LIM
      DO 1 I=1,LIM
         CALL SREAD(NOUT,BUFFER(POINTR),OUTBUF)
         POINTR=POINTR+OUTBUF
    1 CONTINUE
      POINTR=1
      RETURN
      END
