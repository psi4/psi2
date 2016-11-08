      SUBROUTINE INITX2(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Z,VAL1,VAL2,VAL3
      INTEGER ARR,TR1,TR2,ASM,AOS,OS,WTW,WTX,WTY,WAB,SS,SSI,SSJ,SYMORB
      INTEGER BMAX,ORBFRM
      COMMON /IENT/ IDEN,IEXT,IDUM,IDU2,IDU3,INTOD
      COMMON /DIMS/ NBF,NSYM,NORBS,NROWS,NROWS4,NWKS,NWKS2,NLEVS
     *,             NROWOC,NROW4O,NWKSOC,NLEVOC,NORBOC,LEVFRM
     *,             NWKSMX,NLWKMX,NUWKMX,BMAX,NROOTS,ORBFRM
      COMMON /INTS/ NMAX,NMAX2,NGROUP,NBLKOC,NUMIJ,SYMORB,INTSRT
C  UNIVERSAL IDENTITY OF THE OBJECTS IN THESE COMMON
C     COMMON /ALL/ ARR,VAL1,VAL2,VAL3,ITR1,ITR2,IA,JA,ITYPE,ISEGT
C    *,LVFRM1,NLWKI,NLWKJ,IMAX,IMIN
C     COMMON /SYMM/ JSM,JFSYM,IFSYM,MAXSYM(8),MINSYM(8),ISMOFF(8)
C    #,             NUMSYM(8)
      COMMON /SYMM/ ASM,JS,IS,MX(8),MN(8),OS(8),NUMSYM(8)
      COMMON /TAPES/ITAP52,ITAPE5,ITAPE6,ITAP58,ITAP54,ITAP03,ITAP04
     *,             ITAPE3,ITAP05,ITAP06,ITAP56
      COMMON /ALL/VAL1,VAL2,VAL3,ARR,TR1,TR2,IA,JA,M,ISEG,N,N1,N2
     *,           IMAX,IMIN
C
      SQRT2=DSQRT(2.0D+00)
      SQRT3=DSQRT(3.0D+00)
      SQT1P5=DSQRT(1.5D+00)
      IEXT=0
      RETURN
c
c3-24-88  the lines below (between the return above and the entry below)
c       can not be reached and are therefore commented out.
c  10 CONTINUE
c     ITEST=0
c     RETURN
C3-24-88
C
      ENTRY EXTERN(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
CMC   WRITE(6,*) ' ENTRY EXTERNAL '
      ENTRY FOURX(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
CMC   WRITE(6,*) ' ENTRY FOURX '
      IEXT=IEXT+1
C
      AOS=OS(ASM+1)
C*********************************************************************
C     WRITE(*,*) ' AOS,ASM,M = ',AOS,ASM,M
C*********************************************************************
      GO TO (313,305,304,312,311,310,307,308,303,306,309,302,301,314
     1,315,316,317,318,319),M
C GF  GO TO (300,300,300,300,300,300,300,300,300,300,300,300,300,314
C GF 1,315,316,300,300,300),M
C
      PRINT *,' TROUBLE IN INITX2 ENTRY EXTERNAL',M
  300 PRINT *,'  M IS NOT ALL EXTERNAL:',M
      CALL DUMP
      STOP 'FUNNY M'
  301 CONTINUE
C     ----- YZ -----
      GO TO 1
  302 CONTINUE
C     ----- XZ -----
      GO TO 1
  303 CONTINUE
C     ----- XY -----
C GF  CALL XY
C GF  GO TO 9
      GOTO 1
  304 CONTINUE
C     ----- YY -----
      GO TO 1
  305 CONTINUE
C     ----- XX -----
C GF  CALL XX
C GF  GO TO 9
      GOTO 1
  306 CONTINUE
C     ----- ZY -----
      GO TO 1
  307 CONTINUE
C     ----- YX -----
C GF  CALL YX
C GF  GO TO 9
      GOTO 1
  308 CONTINUE
C     ----- YW -----
C GF  CALL YW
C GF  GO TO 9
      GOTO 1
  309 CONTINUE
C     ----- WZ -----
      GO TO 1
  310 CONTINUE
C     ----- WY -----
C GF  CALL WY
C GF  GO TO 9
      GOTO 1
  311 CONTINUE
C     ----- WX -----
C GF  CALL WX
C GF  GO TO 9
      GOTO 1
  312 CONTINUE
C     ----- XW -----
C GF  CALL XW
C GF  GO TO 9
      GOTO 1
  313 CONTINUE
C     ----- WW -----
C GF  CALL WW
C GF  GO TO 9
      GOTO 1
  314 CONTINUE
C     ----- YY 4X -----
      CALL SHAPE4(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
      GO TO 9
  315 CONTINUE
C     ----- XX 4X -----
C GF  CALL XX4X
      CALL SHAPE4(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
      GO TO 9
  316 CONTINUE
C     ----- WW 4X -----
C GF  CALL WW4X
      CALL SHAPE4(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
      GO TO 9
  317 CONTINUE
C     ----- ZZ ENTRIES ----
      GO TO 1
  318 CONTINUE
C     ----- 4 INTERNAL -----
      GO TO 1
  319 CONTINUE
C     ----- 4 INTERNAL, (3,2,1) CASE -----
      GO TO 1
C
C
C
    1 CONTINUE
      CALL SHAPES(int,c,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
    9 CONTINUE
      RETURN
      END
