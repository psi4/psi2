      SUBROUTINE FDIAG(INTS,V,D,DVDMAT
     1,                 NABC,NLWKS,           IWGHT,NABCA,IJADD
     2,                 IJGRP,KADD,LADD,INEXT,IARC,JMNNXT,JMXNXT
     3,                 NABCS,ISYM,            ISEGM,JSEGM,IMAIN
     4,                 ISUB,IUWKMN,IUWKSB,ITRACK,LMIN,ACOEF,BCOEF
     5,                 HDWGT,HDSEG,HDPNT,IMXGRP,IMNGRP,JMXGRP,JMNGRP
     6,                 wtw,wtx,wty,wab,dvdvec,ndvdmx,root
     7,         ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww,klxx,klww)
C
      IMPLICIT REAL*8 (A-H,O-Z)
cets030291
#include <error.h>
      integer errcod, frdi
C
      INTEGER ARR,REFWLK,SYMORB,BMAX,FWORD,FWORD2,ORBFRM
      REAL*8 D(NWKSMX),V(NWKSMX),INTS(NMAX),DVDMAT(LOWTRI)
      real*8 dvdvec(ndvdmx,nroots),root(nroots)
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
      COMMON /DIMS/ NBF,NSYM,NORBS,NROWS,NROWS4,NWKS,NWKS2,NLEVS
     *,             NROWOC,LEVFRM
     *,             NWKSMX,NLWKMX,NUWKMX,BMAX,NROOTS,ORBFRM
      COMMON /INTS/   NMAX,NMAX2,NGROUP,NBLKOC,NUMIJ,SYMORB,INTSRT
      COMMON /OPS/  IGUESS,IPRINT,IRSTRT,IROOTI,IROOTF
      COMMON /SYMM/ JSM,JFSYM,IFSYM,MAXSYM(8),MINSYM(8),ISMOFF(8)
      COMMON /DIAG/ REP,FZCORE,EGUESS,ECI,REFWLK,MXITER,CNVERG,ICNVG
     *,             ITER,SQCDIF,CZERO,NROOT
      COMMON /TAPES/ITAP52,ITAPE5,ITAPE6,ITAP58,ITAP54,ITAP99,ITAP94
     *,             ITAPE3,ITAP95,ITAP96
      COMMON /ALL/  VAL1,VAL2,VAL3,ARR,ITR1,ITR2,IA,JA,ITYPE,ISEGT
     *,             LVFRM1,NLWKI,NLWKJ,IMAX,IMIN
      COMMON /PAGE/ IWORD3,IWORD4,PAGEIN,NOFFI,NOFFJ,NPASS,NWORDL,IWORD6
      COMMON /RSTRT/  IBLCK1,INXT,LOWTRI,NDVDIT,IBLOCK
      DIMENSION IJUNK(100)
      LOGICAL PAGEIN
clj120788 begin.
      integer gugref,ierr
clj120788 end.
      SQCDIF=1.0D+00
      EGUESS=0.0D+00
      ICNVG=-1
clj120788 begin.
      call locate(5,'# GUGAREF ',ierr)
      if (ierr .eq. 0) then
        read(5,*) gugref
        write(6,*) ' ''# GUGAREF #'' found-using alternate reference'
      else
cets030291
        gugref = 0
        errcod = frdi('REFERENCE',gugref)
        if(errcod.eq.EOK) then
          write(6,*) '  using walk number',gugref,' as reference'
        endif
      endif
      refwlk = -1
clj120788 end.
C   PREPARE FOR UNIT GUESS ON LOWEST DIAGONAL MATRIX ELEMENT
      JJ=0
      NWKMX2=NWKSMX+NWKSMX
      NWORD=NWKMX2
      FWORD=1
      DO 60 I=1,NPASS
         IF(I.EQ.NPASS) NWORD=NWORDL
         CALL WREADW(ITAP94,D,NWORD,FWORD,JUNK)
         DO 50 J=1,NWORD/2
clj120788 begin.
            if (jj+j.eq.gugref) then
              refwlk = gugref
              eguess = d(j)
              goto 40
            else if (gugref .gt. 0) then
              goto 40
              endif
clj120788 end.
            IF(D(J).GT.EGUESS) GO TO 40
            REFWLK=JJ+J
            EGUESS=D(J)
  40        D(J)=0.0D+00
  50     CONTINUE
         JJ=JJ+NWORD/2
         IF(IRSTRT.NE.0) CALL WWRITW(ITAP96,D,NWORD,FWORD,JUNK)
         FWORD=FWORD+NWORD
  60  CONTINUE
clj120788 begin.
      if (refwlk .le. 0) then
        write(6,*) 'fdiag: no reference walk was found'
        stop
        endif
clj120788 end.
      IWORD4=NWKS2+1
      IF(IRSTRT.EQ.0) IWORD6=IWORD4
      EGUESS=EGUESS+REP+FZCORE
      WRITE(ITAPE6,70)REFWLK
  70  FORMAT(1H ,'REFERENCE CONFIGURATION   =',I10)
      WRITE(ITAPE6,80)EGUESS
  80  FORMAT(1H ,'REFERENCE ENERGY          = ',G16.10)
      IF(IRSTRT.GT.0) GO TO 300
      ITER=0
      IF(IGUESS.NE.0) GO TO 200
C   USE A UNIT VECTOR AS A STARTING GUESS
      ECI=EGUESS
      NROOT=1
      NDVDIT=1
      WRITE(ITAPE6,100)
 100  FORMAT(/,1H ,'UNIT VECTOR USED AS A STARTING GUESS')
      call fioflu(6)
C   WRITE OUT VECTOR FOR DAVIDSON ALGORITHM
      FWORD=REFWLK+REFWLK-1
      ONE=1.0D+00
      CALL WWRITW(ITAP99,ONE,INTOWP(1),FWORD,JUNK)
      IF(NWKS.LE.NWKSMX) V(REFWLK)=1.0D+00
C   PREPARE FINAL VECTOR FILE
CSEL  NSECT=NROOTS*I2SEC(NWKS2)+I2SEC(112)
CSEL  CALL WFILE(ITAP54,NSECT)
CSEL  CALL ZFILE(ITAP54,NSECT)
      CALL SREW(ITAP54)
C   DO PARTIAL ITERATION FOR REFERENCE INTERACTIONS
      IF(NGROUP.EQ.1) GO TO 110
      CALL SREW(ITAP52)
      CALL RSETSA(ITAP52,INTSRT)
      CALL SREAD(ITAP52,INTS,NMAX2)
  110 CALL LOOPY(nabc,nlwks,iwght,nabca,ijadd,ijgrp,kadd,ladd,inext
     1,          iarc,jmnnxt,jmxnxt,nabcs,isym,isegm,jsegm,imain
     2,          isub,iuwkmn,iuwksb,itrack,lmin,acoef,bcoef
     3,          hdwgt,hdseg,hdpnt,imxgrp,imngrp,jmxgrp,jmngrp
     4,          ints,d,wtw,wtx,wty,wab,v
     5, ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww,klxx,klww)
      CALL DAVIDS(v,d,root,dvdvec,dvdmat,ndvdmx)
      RETURN
C   READ IN OLD VECTOR FROM ITAP54 AS A STARTING GUESS
  200 CONTINUE
      IF(IGUESS.GT.0) GO TO 205
      NVEC=-IGUESS
      WRITE(ITAPE6,201) NVEC
  201 FORMAT(1X,'READING ',I2,' APPROXIMATE VECTOR(S) FROM INPUT.')
      GO TO 210
  205 CONTINUE
      NVEC=IGUESS
      WRITE(ITAPE6,206) IGUESS
  206 FORMAT(1X,'READING ',I2,' OLD VECTOR(S) FROM TAPE 54.')
      CALL SREW(ITAP54)
      CALL DCDLBL(ITAP54,IJUNK,IJUNK,REFWLK,IJUNK,NREFS,EGOLD
     *,           ECIOLD,CZERO)
      WRITE(ITAPE6,207) REFWLK,ECIOLD
  207 FORMAT(1X,'OLD REFERENCE=',I8,/,' OLD CI RESULT=',G17.11)
  210 CONTINUE
      NROOT=IROOTI
      IF(NROOT.GT.NVEC) NROOT=NVEC
      NDVDIT=NVEC
      FWORD=IROUND(112)+1
      FWORD2=1
      DO 230 I=1,NVEC
         IF(IGUESS.GT.0) GO TO 215
         WRITE(ITAPE6,211) I
  211    FORMAT(1X,'VECTOR ',I2)
         CALL NORMAL(IJUNK,D,V,FWORD2)
         RENORM=1.0D+00
         FWORD2=FWORD2+NWKS2
         GO TO 225
  215    CONTINUE
         NWORD=NWKMX2
         RENORM=0.0D+00
         DO 220 J=1,NPASS
            IF(J.EQ.NPASS) NWORD=NWORDL
            NLOOP=NWORD/2
            DO 216 K=1,NLOOP
  216       RENORM=RENORM+V(K)*V(K)
            CALL WREADW(ITAP54,V,NWORD,FWORD,FWORD)
  220    CALL WWRITW(ITAP99,V,NWORD,FWORD2,FWORD2)
  225    IF(I.EQ.1) GO TO 230
         CALL ORTHOG(D,V,INTS,I,RENORM)
  230 CONTINUE
      IF(NVEC.EQ.1) RETURN
      DO 235 J=1,NWKSMX
  235 D(J)=0.0D+00
      FWORD2=1
      IF(NPASS.EQ.1) CALL WREADW(ITAP99,V,NWKS2,FWORD2,JUNK)
      DO 280 I=1,NVEC-1
         CALL SREW(ITAP52)
         CALL RSETSA(ITAP52,INTSRT)
         CALL SREAD(ITAP52,INTS,NMAX2)
         CALL LOOPY(nabc,nlwks,iwght,nabca,ijadd,ijgrp,kadd,ladd,inext
     1,          iarc,jmnnxt,jmxnxt,nabcs,isym,isegm,jsegm,imain
     2,          isub,iuwkmn,iuwksb,itrack,lmin,acoef,bcoef
     3,          hdwgt,hdseg,hdpnt,imxgrp,imngrp,jmxgrp,jmngrp
     4,          ints,d,wtw,wtx,wty,wab,v
     5,ci,si,cj,sj,h1,h2,ijxx,nklxx, ijww,nklww,klxx,klww)
         OFFSET=1
         NWORD=NWKMX2
         DO 270 J=1,NPASS
            IF(J.EQ.NPASS) NWORD=NWORDL
            NLOOP=NWORD/2
            FWORD=OFFSET
            KSIZE=(I-1)*I/2
            CALL WREADW(ITAP94,D,NWORD,FWORD,JUNK)
            CALL WREADW(ITAP99,V,NWORD,FWORD2,JUNK)
            DO 240 K=1,NLOOP
               V(K)=V(K)*D(K)
  240       CONTINUE
            FWORD2=FWORD2+NWKS2
            CALL WREADW(ITAP94,D,NWORD,FWORD2,JUNK)
            DO 245 K=1,NLOOP
               V(K)=V(K)+D(K)
  245       CONTINUE
            DO 260 L=1,I
               Z=0.0D+00
               KSIZE=KSIZE+1
               CALL WREADW(ITAP99,D,NWORD,FWORD,JUNK)
               DO 250 K=1,NLOOP
                  Z=Z+V(K)*D(K)
  250          CONTINUE
               DVDMAT(KSIZE)=DVDMAT(KSIZE)+Z
               FWORD=FWORD+NWKS2
  260       CONTINUE
            OFFSET=OFFSET+NWORD
  270    FWORD2=FWORD2+NWORD-NWKS2
         IWORD4=IWORD4+NWKS2
         IWORD3=IWORD3+NWKS2
         IF(NWKS.EQ.NWKSMX) CALL WREADW(ITAP99,V,NWKS2,IWORD3,JUNK)
         IF(IRSTRT.EQ.0) IWORD6=IWORD4
         DO 280 K=1,NWKSMX
            D(K)=0.0D+00
  280 CONTINUE
      RETURN
C
C   READ IN RESTART INFORMATION AND STARTED MID INTERATION
 300  CONTINUE
      IRSTRT=-1
      WRITE(ITAPE6,310) ITER
 310  FORMAT(1H ,80('*')/'0CALCULATION RESTARTED FROM ITERATION',I4)
      CALL DAVIDS(v,d,root,dvdvec,dvdmat,ndvdmx)
      DELTA=ECI-EGUESS
      WRITE(ITAPE6,320)NROOT,ITER,ECI,DELTA,SQCDIF
      WRITE(ITAPE3,320)NROOT,ITER,ECI,DELTA,SQCDIF
      call fioflu(6)
 320  FORMAT(1H ,'ROOT',I3,' ITER',I3,' ENERGY ',G17.9 ,' DELTA ',E9.3
     *,      ' SQCDIF ',E9.3)
      RETURN
      END
