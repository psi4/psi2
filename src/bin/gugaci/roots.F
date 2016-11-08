      SUBROUTINE ROOTS(INTS,V,S,ROOT
     1,                 NABC,NLWKS,           IWGHT,NABCA,IJADD
     2,                 IJGRP,KADD,LADD,INEXT,IARC,JMNNXT,JMXNXT
     3,                 NABCS,ISYM,            ISEGM,JSEGM,IMAIN
     4,                 ISUB,IUWKMN,IUWKSB,ITRACK,LMIN,ACOEF,BCOEF
     5,                 HDWGT,HDSEG,HDPNT,IMXGRP,IMNGRP,JMXGRP,JMNGRP
     6,                 wtw,wtx,wty,wab,dvdvec,dvdmat,ndvdmx
     7,           ci,si,cj,sj,h1,h2,ijxx,nklxx, ijww,nklww,klxx,klww)

C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER REFWLK,SYMORB,BMAX,FWORD,ORBFRM
      REAL*8 INTS(NMAX),V(NWKSMX),S(NWKSMX),ROOT(NROOTS)
      real*8 dvdvec(ndvdmx,nroots),dvdmat(lowtri)
      real*8 ci(1),cj(1),si(1),sj(1),h1(1),h2(1)
      integer ijxx(numij),nklxx(nsym,orbfrm),ijww(numij)
      integer nklww(nsym,orbfrm),klxx(1),klww(1)
      integer ndvdmx
      integer nabc(nrows),nlwks(nrows),iwght(nrows4),nabca(nrows)
      integer ijadd(numij),ijgrp(numij),kadd(symorb),ladd(symorb)
      integer inext(norbs),jmnnxt(norbs),jmxnxt(norbs),nabcs(nrows)
      integer isym(norbs),isegm(nlevs),jsegm(nlevs),imain(nlevs)
      integer isub(nlevs),iuwkmn(nlevs),iuwksb(nlevs),itrack(nlevs)
      integer lmin(nlevs),hdwgt(nlevs),hdseg(nlevs),hdpnt(nlevs)
      integer imxgrp(ngroup),imngrp(ngroup),jmxgrp(ngroup)
      integer jmngrp(ngroup),iarc(nrows4)
      integer wab(orbfrm),wtw(orbfrm,nsym),wty(orbfrm),wtx(orbfrm,nsym)
      real*8 acoef(nlevs),bcoef(nlevs)
      COMMON /TAPES/ITAP52,ITAPE5,ITAPE6,ITAP58,ITAP54,ITAP99,ITAP94
     *,             ITAPE3,ITAP95,ITAP96
      COMMON /DIMS/ NBF,NSYM,NORBS,NROWS,NROWS4,NWKS,NWKS2,NLEVS
     *,             NROWOC,LEVFRM
     *,             NWKSMX,NLWKMX,NUWKMX,BMAX,NROOTS,ORBFRM
      COMMON /INTS/   NMAX,NMAX2,NGROUP,NBLKOC,NUMIJ,SYMORB,INTSRT
      COMMON /OPS/  IGUESS,IPRINT,IRSTRT,IROOTI,IROOTF
      COMMON /DIAG/ REP,FZCORE,EGUESS,ECI,REFWLK,MXITER,CNVERG,ICNVG
     *,             ITER,SQCDIF,CZERO,NROOT
      common /rstrt/ iblck1,inxt,lowtri,ndvdit,iblock
C   RESET INTEGRALS FOR NEXT ITERATION
   10 IF(NGROUP.EQ.1) GO TO 15
      CALL SREW(ITAP52)
      CALL RSETSA(ITAP52,INTSRT)
      CALL SREAD(ITAP52,INTS,NMAX2)
   15 ITER=ITER+1
      CALL LOOPY(nabc,nlwks,iwght,nabca,ijadd,ijgrp,kadd,ladd,inext
     1,          iarc,jmnnxt,jmxnxt,nabcs,isym,isegm,jsegm,imain
     2,          isub,iuwkmn,iuwksb,itrack,lmin,acoef,bcoef
     3,          hdwgt,hdseg,hdpnt,imxgrp,imngrp,jmxgrp,jmngrp
     4,          ints,s,wtw,wtx,wty,wab,v
     5,ci,si,cj,sj,h1,h2,ijxx,nklxx, ijww,nklww,klxx,klww)
CTJL MODIFIED TO GET REFERENCE COEFFICIENT BEFORE CALLING DAVIDS
CTJL           8/3/85
      FWORD=INTOWP(REFWLK-1)+IROUND(112)+1
      CALL WREADW(ITAP54,CZERO,INTOWP(1),FWORD,JUNK)
      CALL DAVIDS(v,s,root,dvdvec,dvdmat,ndvdmx)
      IF(ICNVG.EQ.0) GO TO 110
      IF(ITER.LT.MXITER) GO TO 10
      WRITE(ITAPE6,30)
   30 FORMAT(1H ,'ITERATIONS EXHAUSTED, CALCULATION DID NOT CONVERGE')
C   CALCULATE DAVIDSON CORRECTION TO ENERGY
  110 FWORD=INTOWP(REFWLK-1)+IROUND(112)+1
      CALL WREADW(ITAP54,CZERO,INTOWP(1),FWORD,JUNK)
      DO 120 I=1,NROOTS
  120 ROOT(I)=ROOT(I)+REP+FZCORE
      EDVDSN=ROOT(1)+(ROOT(1)-EGUESS)*(1.0D+00-CZERO*CZERO)
C   WRITE FINAL CALCULATION DATA
      WRITE(ITAPE6,130)
  130 FORMAT(/,1H*,8('**********')/)
      WRITE(ITAPE6,140)EGUESS
  140 FORMAT(6X,'CALCULATION RESULTS',/,1H ,'EGUESS        = ',F20.12)
      WRITE(ITAPE6,150)(I,ROOT(I),I=1,NROOTS)
  150 FORMAT(1H ,'ROOT',I3,' ECI   = ',F20.12)
      WRITE(ITAPE6,160)
  160 FORMAT(/,1H ,'DAVIDSON CORRECTION FOR LOWEST ROOT')
      WRITE(ITAPE6,170)CZERO
  170 FORMAT(1H ,3X,'REFERENCE COEFFICIENT = ',F20.12)
      WRITE(ITAPE6,180)EDVDSN
  180 FORMAT(1H ,3X,'CORRECTED ENERGY      = ',F20.12,/)
      RETURN
      END
