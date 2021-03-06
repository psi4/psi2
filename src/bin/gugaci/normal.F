      SUBROUTINE NORMAL(IWALK,WGT,VEC2,FWORD)
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IWALK(100),WGT(100),VEC2(NWKSMX)
      INTEGER FWORD,BMAX,ORBFRM
      COMMON /TAPES/ITAP52,ITAPE5,ITAPE6,ITAP58,ITAP54,ITAP99,ITAP94
     *,             ITAPE3,ITAP95,ITAP96
      COMMON /DIMS/ NBF,NSYM,NORBS,NROWS,NROWS4,NWKS,NWKS2,NLEVS
     *,             NROWOC,LEVFRM
     *,             NWKSMX,NLWKMX,NUWKMX,BMAX,NROOTS,ORBFRM
      READ(ITAPE5,5) NWALK
    5 FORMAT(I5)
      IF(NWALK.GT.100) THEN
          WRITE(ITAPE6,*) ' NWALK.GT.100'
          CALL MABORT
      END IF
      WRITE(ITAPE6,6)NWALK
    6 FORMAT(1X,I3,' ELEMENTS IN GUESS.')
      RENORM=0.0D+00
      DO 15 I=1,NWALK
         READ(ITAPE5,10) IWALK(I),WGT(I)
   10    FORMAT(I5,F10.6)
C        WRITE(ITAPE6,11) IWALK(I),WGT(I)
   11    FORMAT(1X,I5,F10.6)
   15 RENORM=RENORM+WGT(I)*WGT(I)
      RENORM=DSQRT(RENORM)
      DO 20 I=1,NWKSMX
   20 VEC2(I)=0.0D+00
      DO 25 I=1,NWALK
         WGT(I)=WGT(I)/RENORM
         IF(IWALK(I).LT.NWKSMX) VEC2(IWALK(I))=WGT(I)
         IWORD=FWORD+(IWALK(I)-1)*2
   25 CALL WWRITW(ITAP99,WGT(I),2,IWORD,JUNK)
      RETURN
      END
