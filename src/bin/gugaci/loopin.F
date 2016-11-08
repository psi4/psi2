cets110990      SUBROUTINE INITLP
      subroutine loopin(int,d,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     1 vector,ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww,klxx,klww)
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER PUWK,ARR,REFWLK,PUWKT,SYMORB,BMAX,ORBFRM
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
C-----------------------------------------------------------LOOPIN
C
cets110990      ENTRY LOOPIN
C
C     WRITE(ITAPE3,*) ' ENTRY LOOPIN'
      PUWK=PUWKT
      ITR1=ITRAK(1,ITRACK)
      ITR2=ITRAK(2,ITRACK)
      IF(ITR2.EQ.0) VAL2=0.0D+00
      IF(ITRACK.EQ.6) VAL2=1.0D+00
      IF(ITRACK.EQ.7) VAL2=1.0D+00
      IF(NDVDIT.GT.0) GO TO 25
      ITYPE=7
      IF(IUWK.NE.JUWK) RETURN
      JA=PUWKT+IUWK
      IF(PAGEIN) GO TO 21
      NLWKT=NLWKI
      JA=1
   20 IF(NLWKT.LE.0)RETURN
      NLWKI=MIN(NLWKT,NWKSMX)
      CALL PAGED(vector,d)
      CALL DIAGNL(int,d,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
      CALL PAGEOT(vector,d)
      NLWKT=NLWKT-NLWKI
      IUWK=IUWK+NLWKI
      GO TO 20
   21 JA=JA-NOFFJ
      CALL DIAGNL(int,d,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
      RETURN
   25 IF(IUWK.EQ.JUWK) RETURN
      ITYPE=18
      IF(ITRACK.EQ.8)ITYPE=19
      IA=PUWKT+IUWK
      JA=PUWKT+JUWK
      IF(PAGEIN) GO TO 27
      IA=1
      NLWKT=NLWKI
   26 IF(NLWKT.LE.0)RETURN
      NLWKI=MIN(NLWKT,NWKDV2)
      NLWKJ=NLWKI
      JA=NLWKI+1
      CALL PAGEO(vector,d)
      CALL fourx(int,vector,d,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     1 ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww,klxx,klww)
      CALL PAGEOT(vector,d)
      IUWK=IUWK+NLWKI
      JUWK=JUWK+NLWKI
      NLWKT=NLWKT-NLWKI
      GO TO 26
   27 JA=JA-NOFFJ
      IA=IA-NOFFI
      CALL fourx(int,vector,d,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     1 ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww,klxx,klww)
      RETURN
C
      END
