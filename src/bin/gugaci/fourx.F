      SUBROUTINE fourx(INT,C,S,IJADD,KADD,LADD,WTW,WTX,WTY,WAB,SS,
     1 ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww,klxx,klww)
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 sqrt2,sqrt3,sqt1p5,asqrt2,hsqrt3
      INTEGER ARR,TR1,TR2,ASM,AOS,OS,WTW,WTX,WTY,WAB,SS,SYMORB
      INTEGER BMAX,ORBFRM
      REAL*8 INT(NMAX),C(NWKS),S(NWKS)
      real*8 ci(1),cj(1),si(1),sj(1),h1(1),h2(1)
      integer ijxx(numij),nklxx(nsym,orbfrm),ijww(numij)
      integer nklww(nsym,orbfrm),klxx(1),klww(1)
      DIMENSION KADD(SYMORB),LADD(SYMORB),IJADD(NUMIJ),WTW(ORBFRM,NSYM)
      DIMENSION WTX(ORBFRM,NSYM),WTY(ORBFRM),WAB(ORBFRM),SS(NORBS)
      COMMON /DIMS/ NBF,NSYM,NORBS,NROWS,NROWS4,NWKS,NWKS2,NLEVS
     *,             NROWOC,LEVFRM
     *,             NWKSMX,NLWKMX,NUWKMX,BMAX,NROOTS,ORBFRM
      COMMON /OPS/  IGUESS,IPRINT,IRSTRT,IROOTI,IROOTF,I34X
      COMMON /INTS/   NMAX,NMAX2,NGROUP,NBLKOC,NUMIJ,SYMORB,INTSRT
C  UNIVERSAL IDENTITY OF THE OBJECTS IN THESE COMMON
C     COMMON /ALL/ ARR,VAL1,VAL2,VAL3,ITR1,ITR2,IA,JA,ITYPE,ISEGT
C    *,LVFRM1,NLWKI,NLWKJ,IMAX,IMIN
C     COMMON /SYMM/ JSM,JFSYM,IFSYM,MAXSYM(8),MINSYM(8),ISMOFF(8)
C    #,             NUMSYM(8)
      COMMON /SYMM/ ASM,JS,IS,MX(8),MN(8),OS(8),NUMSYM(8)
      COMMON /TAPES/ITAP52,ITAPE5,ITAPE6,ITAP58,ITAP54,ITAP99,ITAP94
     *,             ITAPE3,ITAP95,ITAP96
      COMMON /ALL/VAL1,VAL2,VAL3,ARR,TR1,TR2,IA,JA,M,ISEG,N,N1,N2
     *,           IMAX,IMIN
      COMMON /COUNT/  ICOUNT,IXX4,IXX5,IXX6,IXX8,IWW4,IWW5,IWW6,IWW8
     #,               IWX7,IXW9,IXY3,IXY16,IXY18,IXY22,IWY2,IWY15,IWY17
     #,               IWY21,IYX20,IYW19
      COMMON /COUNT2/ IENTRY(20),TIME(20)
      common /sqrts/ sqrt2,sqrt3,sqt1p5,asqrt2,hsqrt3
      REAL*8 VAL1,VAL2,VAL3
cets111090      SQRT2=DSQRT(2.0D+00)
cets111090      SQRT3=DSQRT(3.0D+00)
cets111090      SQT1P5=DSQRT(1.5D+00)
cets110990      ICOUNT=0
cets110990      IXX4=0
cets110990      IXX5=0
cets110990      IXX6=0
cets110990      IXX8=0
cets110990      IWW4=0
cets110990      IWW5=0
cets110990      IWW6=0
cets110990      IWW8=0
cets110990      IWX7=0
cets110990      IXW9=0
cets110990      RETURN
C
C--------------------------------------------------------EXTERNAL
C--------------------------------------------------------FOURX---
C
cets110990      ENTRY EXTERN
cets110990      ENTRY FOURX
C
C
      IF (I34X.NE.0) THEN
         DO 11 IQ=1,NSYM
            MINS=MN(IQ)
            IF (MINS.GT.N) THEN
               NUMSYM(IQ)=0
            ELSE
               MAXS=MX(IQ)
               IF (MAXS.LE.N) THEN
                  NUMSYM(IQ)=MAXS-MINS+1
               ELSE
                  NUMSYM(IQ)=N-MINS+1
               END IF
            END IF
   11    CONTINUE
      END IF
C
C
      IENTRY(M)=IENTRY(M)+1
      ICOUNT=ICOUNT+1
C
      AOS=OS(ASM+1)
      GO TO (313,305,304,312,311,310,307,308,303,306,309,302,301,314
C     GO TO (301,302,303,304,305,306,307,308,309,310,311,312,313,314
C     GO TO ( YZ, XZ, XY, YY, XX, ZY, YX, YW, WZ, WY, WX, XW, WW,---YY
C
     1,315,316,317,318,319),M
C       XX, WW, ZZ,---4 INTERNALS),M
C
      WRITE(ITAPE6,*) ' FUNNY M'
      CALL MABORT
  301 CONTINUE
C     ----- YZ -----
      GO TO 1
  302 CONTINUE
C     ----- XZ -----
      GO TO 1
  303 CONTINUE
C     ----- XY -----
      IF (IAND(I34X,1).NE.0.AND.ISEG.EQ.3) GO TO 1
      IF (IAND(I34X,1).EQ.0.AND.N.LT.ORBFRM.AND.ISEG.EQ.18) GO TO 9
      CALL XY(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,ci,si,cj,sj,
     1 h1,h2,ijxx,nklxx,ijww,nklww,klxx,klww)
      GO TO 9
  304 CONTINUE
C     ----- YY -----
      GO TO 1
  305 CONTINUE
C     ----- XX -----
      CALL XX(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,ci,si,cj,sj,
     1 h1,h2,ijxx,nklxx,ijww,nklww,klxx,klww)
      GO TO 9
  306 CONTINUE
C     ----- ZY -----
      GO TO 1
  307 CONTINUE
C     ----- YX -----
      CALL YX(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,ci,si,cj,sj,
     1 h1,h2,ijxx,nklxx,ijww,nklww,klxx,klww)
      GO TO 9
  308 CONTINUE
C     ----- YW -----
      CALL YW(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,ci,si,cj,sj,
     1 h1,h2,ijxx,nklxx,ijww,nklww,klxx,klww)
      GO TO 9
  309 CONTINUE
C     ----- WZ -----
      GO TO 1
  310 CONTINUE
C     ----- WY -----
      IF (IAND(I34X,2).NE.0.AND.ISEG.EQ.2) GO TO 1
      IF (IAND(I34X,2).EQ.0.AND.N.LT.ORBFRM.AND.ISEG.EQ.17) GO TO 9
      CALL WY(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,ci,si,cj,sj,
     1 h1,h2,ijxx,nklxx,ijww,nklww,klxx,klww)
      GO TO 9
  311 CONTINUE
C     ----- WX -----
      CALL WX(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,ci,si,cj,sj,
     1 h1,h2,ijxx,nklxx,ijww,nklww,klxx,klww)
      GO TO 9
  312 CONTINUE
C     ----- XW -----
      CALL XW(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,ci,si,cj,sj,
     1 h1,h2,ijxx,nklxx,ijww,nklww,klxx,klww)
      GO TO 9
  313 CONTINUE
C     ----- WW -----
      CALL WW(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,ci,si,cj,sj,
     1 h1,h2,ijxx,nklxx,ijww,nklww,klxx,klww)
      GO TO 9
  314 CONTINUE
C     ----- YY 4X -----
      CALL SHAPE4(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
      GO TO 9
  315 CONTINUE
C     ----- XX 4X -----
C     IF (I34X.EQ.0) THEN
      CALL XX4XS(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,ci,si,cj,sj,
     1 h1,h2,ijxx,nklxx,ijww,nklww,klxx,klww)
C     ELSE
C        CALL SHAPE4(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
C     END IF
      GO TO 9
  316 CONTINUE
C     ----- WW 4X -----
C     IF (I34X.EQ.0) THEN
      CALL WW4X(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,ci,si,cj,sj,
     1 h1,h2,ijxx,nklxx,ijww,nklww,klxx,klww)
C     ELSE
C        CALL SHAPE4(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
C     END IF
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
    1 CONTINUE
      CALL SHAPES(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
    9 CONTINUE
      RETURN
      END
