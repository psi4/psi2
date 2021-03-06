      SUBROUTINE RENUMB(LEVPT,LEVNR,ARC)
C
      IMPLICIT INTEGER (A-Z)
      INTEGER NUMINT
C
      COMMON /DIMENS/ NBF,NSYM,NORBS,NROWSP,NROW4P,NROWS,NROWS4
     #,               NLEVS,NREFS,NROWOC,NROW4O,NWKS,NWKSOC,NLEVOC
     #,               ORBFRM,SYMORB,NUMIJ,NGROUP,NUMINT,NMAX,NSPC,NVREF
     #,               NIJVIR
C
      DIMENSION LEVPT(NLEVS),LEVNR(NLEVS),ARC(NROWS4)
C
C
C     ----- CHANGE ARC ARRAY TO INCORPORATE LEVPT OFFSETS -----
C           SO DONT HAVE TO ADD LEVPT(LEVM1) ALL THE TIME
C
      DO 3 LEV=2,NLEVS
         PONTM1=LEVPT(LEV-1)
         DO 2 ROW=LEVPT(LEV)+1,LEVPT(LEV)+LEVNR(LEV)
            DO 1 CASE=(ROW-1)*4+1,(ROW-1)*4+4
               IF (ARC(CASE).GT.0) ARC(CASE)=ARC(CASE)+PONTM1
    1       CONTINUE
    2    CONTINUE
    3 CONTINUE
      RETURN
      END
