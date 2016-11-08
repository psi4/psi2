cets110990      SUBROUTINE INITLP
      subroutine loopex(int,d,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss
     1, vector,ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww,klxx,klww)
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER ARR,REFWLK,PUWKT,SYMORB,BMAX,ORBFRM
      real*8 int(nmax),d(nwks),vector(nwks)
      real*8 ci(1),cj(1),si(1),sj(1),h1(1),h2(1)
      integer ijxx(numij),nklxx(nsym,orbfrm),ijww(numij)
      integer nklww(nsym,orbfrm),klxx(1),klww(1)
      integer wab(orbfrm),wtw(orbfrm,nsym),wty(orbfrm),wtx(orbfrm,nsym)
      integer ss(norbs)
      COMMON /LOOPS/NUWK,PUWKT,IUWK,JUWK,ITRACK,IPT1,IPT2
      COMMON /TAPES/ITAP52,ITAPE5,ITAPE6,ITAP58,ITAP54,ITAP99,ITAP94
     *,             ITAPE3,ITAP95,ITAP96
      COMMON /ALL/  VAL1,VAL2,VAL3,ARR,ITR1,ITR2,IA,JA,ITYPE,ISEG
     *,             LVFRM1,NLWKI,NLWKJ,IMAX,IMIN
      COMMON /DIMS/ NBF,NSYM,NORBS,NROWS,NROWS4,NWKS,NWKS2,NLEVS
     *,             NROWOC,LEVFRM
     *,             NWKSMX,NLWKMX,NUWKMX,BMAX,NROOTS,ORBFRM
      COMMON /INTS/   NMAX,NMAX2,NGROUP,NBLKOC,NUMIJ,SYMORB,INTSRT
      COMMON /OPS/  IGUESS,IPRINT,IRSTRT,IROOTI,IROOTF
      COMMON /SYMM/ JSM,JFSYM,IFSYM,MAXSYM(8),MINSYM(8),ISMOFF(8)
      COMMON /DIAG/ REP,FZCORE,EGUESS,ECI,REFWLK,MXITER,CNVERG,ICNVG
     *,             ITER,SQCDIF,CZERO,NROOT
      COMMON /RSTRT/  IBLCK1,INXT,LOWTRI,NDVDIT,IBLOCK
      COMMON /PAGE/ IWORD3,IWORD4,PAGEIN,NOFFI,NOFFJ,NPASS,NWORDL,IWORD6
      LOGICAL PAGEIN
      DIMENSION ITRAK(2,10),ITYPE1(4,4)
      DATA ITRAK/1,0,2,0,3,0,1,3,3,2,3,2,3,1,0,0,1,2,2,1/
      DATA ITYPE1/1,4,8,-1,5,2,7,-1,6,9,3,10,11,12,13,17/
      NWKDV2=NWKSMX/2
cets110990      RETURN
C
C-------------------------------------------------------LOOPEX
C
cets110990      ENTRY LOOPEX
C
C     WRITE(ITAPE3,*)' ENTRY LOOPEX'
      ITR1=ITRAK(1,ITRACK)
      ITR2=ITRAK(2,ITRACK)
      ITYPE=ITYPE1(IPT2,IPT1)
      IF(ITRACK.EQ.6) VAL2=1.0D+00
      IF(ITRACK.EQ.7) VAL2=1.0D+00
      IF(ITR2.EQ.0) VAL2=0.0D+00
      IF(NDVDIT.GT.0) GO TO 40
      IF(IUWK.NE.JUWK) RETURN
      JA=PUWKT+JUWK-NOFFJ
      CALL DIAGNL(int,d,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
      RETURN
C   CALCULATE NON-DIAGONAL PARTIAL EXTERNALS
   40 CONTINUE
      IF(NDVDIT.GT.1) GO TO 45
      IF(IGUESS.EQ.0.AND.ITYPE.LT.10) RETURN
   45 CONTINUE
      IA=PUWKT+IUWK-NOFFI
      JA=PUWKT+JUWK-NOFFJ
cets110990      CALL EXTERN
      call fourx(int,vector,d,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     1 ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww,klxx,klww)
      RETURN
      END
