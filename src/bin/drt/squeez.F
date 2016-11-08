      SUBROUTINE SQUEEZ(BP,SP,ARCP,NLWKSP,B,S,ARC,NLWKS)
C
C**********************************************************************
C     PACK DOWN THE DRT ARRAYS NOW THAT INVALID ARCS AND POINTS
C     HAVE BEEN ELIMINATE.
C**********************************************************************
C
      IMPLICIT INTEGER (A-Z)
      INTEGER NUMINT
C
      COMMON /DIMENS/ NBF,NSYM,NORBS,NROWSP,NROW4P,NROWS,NROWS4
     #,               NLEVS,NREFS,NROWOC,NROW4O,NWKS,NWKSOC,NLEVOC
     #,               ORBFRM,SYMORB,NUMIJ,NGROUP,NUMINT,NMAX,NSPC,NVREF
     #,               NIJVIR
C
      DIMENSION BP(NROWSP),SP(NROWSP),ARCP(NROW4P),NLWKSP(NROWSP)
      DIMENSION B(NROWS),S(NROWS),ARC(NROWS4),NLWKS(NROWS)
C
      DO 1 I=1,NROWS
         B(I)=BP(I)
    1 CONTINUE
      DO 2 I=1,NROWS
         S(I)=SP(I)
    2 CONTINUE
      DO 3 I=1,NROWS4
         ARC(I)=ARCP(I)
    3 CONTINUE
      DO 4 I=1,NROWS
         NLWKS(I)=NLWKSP(I)
    4 CONTINUE
      RETURN
      END
