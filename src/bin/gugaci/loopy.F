      SUBROUTINE LOOPY(NABC,NLWKS,           IWGHT,NABCA,IJADD
     1,                 IJGRP,KADD,LADD,INEXT,IARC,JMNNXT,JMXNXT
     2,                 NABCS,ISYM,            ISEGM,JSEGM,IMAIN
     3,                 ISUB,IUWKMN,IUWKSB,ITRACK,LMIN,ACOEF,BCOEF
     4,                 HDWGT,HDSEG,HDPNT,IMXGRP,IMNGRP,JMXGRP,JMNGRP
     5, int,sigma,wtw,wtx,wty,wab,vector,ci,si,cj,sj,h1,h2,ijxx,nklxx
     6, ijww,nklww,klxx,klww)
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER PUWKT,REFWLK,SYMORB,BMAX,ORBFRM
CTJL 1               INEXT,IARC,JMNNXT,JMXNXT,NABCS,ISYM,ISEGM,JSEGM,
CTJL 2               IMAIN,ISUB,IUWKMN,IUWKSB,ITRACK,LMIN,ACOEF,BCOEF,
CTJL 3               HDWGT,HDSEG,HDPNT,IMXGRP,IMNGRP,JMXGRP,JMNGRP
      real*8 int(nmax),sigma(nwks),vector(nwks)
      real*8 ci(1),cj(1),si(1),sj(1),h1(1),h2(1)
      integer wab(orbfrm),wtw(orbfrm,nsym),wty(orbfrm),wtx(orbfrm,nsym)
      integer ijxx(numij),nklxx(nsym,orbfrm),ijww(numij)
      integer nklww(nsym,orbfrm),klxx(1),klww(1)
cets111090      DIMENSION NLCSMN(22),LCOND(8),COEFFS(20,21),CFS(420)
      dimension nlcsmn(22),cfs(420)
      INTEGER HPT,HSEG,HSEGMX,HDWGT,HDSEG,HDPNT,SEGWGT,HARPT
      LOGICAL PAGEIN
C   LEVEL CHARACTERISTICS
      DIMENSION ISEGM(NLEVS),JSEGM(NLEVS),IMAIN(NLEVS),ISUB(NLEVS)
      DIMENSION IUWKMN(NLEVS),IUWKSB(NLEVS),ITRACK(NLEVS),ISYM(NORBS)
      DIMENSION ACOEF(NLEVS),BCOEF(NLEVS),LMIN(NLEVS)
      DIMENSION HDSEG(NLEVS),HDPNT(NLEVS),HDWGT(NLEVS)
C     DIMENSION LEVPT(NLEVS),LEVNR(NLEVS)
C   GRAPH DESCRIPTION ARRAYS
      DIMENSION NABC(NROWS ),NLWKS(NROWS)              ,NABCS(NROWS )
      DIMENSION              IWGHT(NROWS4),IARC(NROWS4),NABCA(NROWS )
C   INTEGRAL ADDRESSING ARRAYS
      DIMENSION INEXT(NORBS),JMNNXT(NORBS),JMXNXT(NORBS)
      DIMENSION IJADD(NUMIJ),IJGRP(NUMIJ),KADD(SYMORB),LADD(SYMORB)
      DIMENSION IMNGRP(NGROUP),IMXGRP(NGROUP),JMNGRP(NGROUP)
      DIMENSION JMXGRP(NGROUP)
C
      COMMON /LEVEL/ LEVI
C
      COMMON /SYMM/ JSMT,JFSYM,IFSYM,MAXSYM(8),MINSYM(8),ISMOFF(8)
     #,             NUMSYM(8)
      COMMON /DIMS/ NBF,NSYM,NORBS,NROWS,NROWS4,NWKS,NWKS2,NLEVS
     *,             NROWOC,LEVFRM
     *,             NWKSMX,NLWKMX,NUWKMX,BMAX,NROOTS,ORBFRM
      COMMON /INTS/   NMAX,NMAX2,NGROUP,NBLKOC,NUMIJ,SYMORB,INTSRT
      COMMON /DIAG/ REP,FZCORE,EGUESS,ECI,REFWLK,MXITER,CNVERG,ICNVG
     *,             ITER,SQCDIF,CZERO,NROOT
      COMMON /TAPES/ITAP52,ITAPE5,ITAPE6,ITAP58,ITAP54,ITAP99,ITAP94
     *,             ITAPE3,ITAP95,ITAP96
      COMMON /OPS/  IGUESS,IPRINT,IRSTRT,IROOTI,IROOTF
      COMMON/LOOPS/ NUWK,PUWKT,IUWK,JUWK,ITRAK,IPT1,IPT2
      COMMON/ALL/   ACF,D,CCF,LADT,ITR1,ITR2,IA,JA,ITYPE,ISEGT,LVFRM1
     *,             NLWKI,NLWKJ,IMAX,IMIN
      COMMON /PAGE/ IWORD3,IWORD4,PAGEIN,NOFFI,NOFFJ,NPASS,NWORDL,IWORD6
      COMMON /RSTRT/IBLCK1,INXT,LOWTRI,NDVDIT,IBLOCK
      COMMON /TABLES/ JSEGNR(22),JSEGPT(22),IARCMN(228),IARCSB(228)
     *,ITRK(228),JCOND(228),KCOND(228),NXTSEG(228),JSEGPX(3)
      COMMON /MINMAX/ IMING,IMAXG,JMIN,JMAX
      common /lopstf/ coeffs(20,21),lcond(8)
C
      EQUIVALENCE (COEFFS,CFS)
      DATA NLCSMN/6*1,0,1,14*0/
C
      CRITE = 0.00001D+00
      ROOT2 = DSQRT(2.0D+00)
      ROOTN2 = -ROOT2
      TOOR2 = 1.0D+00 / ROOT2
      TOORN2 = -TOOR2
      HSEGMX=4
c
cets110990      ENTRY LOOPY
cets111090      entry LOOPY(NABC,NLWKS,           IWGHT,NABCA,IJADD
cets111090     1,                 IJGRP,KADD,LADD,INEXT,IARC,JMNNXT,JMXNXT
cets111090     2,                 NABCS,ISYM,            ISEGM,JSEGM,IMAIN
cets111090     3,                 ISUB,IUWKMN,IUWKSB,ITRACK,LMIN,ACOEF,BCOEF
cets111090     4,                 HDWGT,HDSEG,HDPNT,IMXGRP,IMNGRP,JMXGRP,JMNGRP
cets111090     5, int,sigma,wtw,wtx,wty,wab,vector,ci,si,cj,sj,h1,h2,ijxx,nklxx
cets111090     6, ijww,nklww,klxx,klww)
c
      DO 500 IBLOCK=1,NGROUP
C:
         IF (IBLOCK.GT.1) CALL NXTBLK(int,nmax)
C:
         IF(NDVDIT.GT.0) GOTO 188
         IF(IPRINT.GT.0) WRITE(ITAPE6,187)IBLOCK,IMXGRP(IBLOCK),
     #   IMNGRP(IBLOCK),JMXGRP(IBLOCK),JMNGRP(IBLOCK)
 187     FORMAT(' PROCESS INTEGRALS FROM GROUP',I4,' FIRST LEVEL',I4
     #   ,      ' LAST',I4,' JMAX',I4,' JMIN',I4)
C:
 188     CONTINUE
         IMAXG=IMXGRP(IBLOCK)
         IMING=IMNGRP(IBLOCK)
         I=IMAXG
C     PRINT *,' BEFORE 190,IMAX,MIN,AND,I',IMAXG,IMING,I
C:
 190     CONTINUE
C     PRINT *,' PAST 190,I=',I,IMAXG,IMING
C:    I=INEXT(INXT)
         LEVI=I+1
         NPT=1
         PUWKT=1
C:    JMAX=JMXNXT(INXT)
C:    JMIN=JMNNXT(INXT)
         JMAX=I
         JMIN=1
         IF (I.EQ.IMAXG) JMAX=JMXGRP(IBLOCK)
         IF (I.EQ.IMING) JMIN=JMNGRP(IBLOCK)
C     PRINT *,' JMIN,MAX',JMIN,JMAX
C:
         JSM=ISYM(I)
         IAD=(I*(I-1))/2
         IJ=IAD+JMAX
         NOFFI=0
         NOFFJ=0
         LEVV=NLEVS
         HSEG=4
         LEVH=NLEVS
         PAGEIN=.FALSE.
         IF(NWKS.LE.NWKSMX) PAGEIN=.TRUE.
C:    IF(IBLOCK.NE.IJGRP(IJ)) GO TO 490
C*
C     CALL NXT4X(I)
C     IF (I.LE.ORBFRM) CALL NXT4X(I)
C*
         IF(LEVI.EQ.NLEVS) GO TO 1070
C   GENERATE HEAD SEGMENTS TO LOOPS
         LEVH=NLEVS+1
         NPT=1
         HDWGT(NLEVS)=0
         HDPNT(NLEVS)=1
 1000    LEVH=LEVH-1
         SEGWGT=0
         HSEG=0
         HPT=NPT
 1010    HSEG=HSEG+1
C     PRINT *,'AFTER 1010, HSEG',HSEG
         IF(HSEG.LE.HSEGMX) GO TO 1030
         LEVH=LEVH+1
C     PRINT *,'AFTER 1020, LEVH',LEVH,'NLEVS',NLEVS
         IF(LEVH.GT.NLEVS) GO TO 480
         PUWKT=PUWKT-HDWGT(LEVH)
         HSEG=HDSEG(LEVH)
         HPT=HDPNT(LEVH)
         GO TO 1010
 1030    CONTINUE
         HARPT=(HPT-1)*4+HSEG
         NPT=IARC(HARPT)
         IF(NPT.EQ.0) GO TO 1010
         PUWKT=PUWKT-SEGWGT
         SEGWGT=IWGHT(HARPT)
         PUWKT=PUWKT+SEGWGT
         HDWGT(LEVH)=SEGWGT
         HDSEG(LEVH)=HSEG
         LEV=LEVH-1
         HDPNT(LEV)=NPT
         IUWK=0
         JUWK=0
         NLWKI=NLWKS(NPT)
         IF(PAGEIN) GO TO 1035
         CALL PAGED(vector,sigma)
         IF(PAGEIN) LEVV=LEV
         GO TO 1040
 1035    CONTINUE
         IF(LEVV.GT.LEV) GO TO 1040
         CALL PAGEOT(vector,sigma)
         CALL PAGED(vector,sigma)
         LEVV=NLEVS
         IF(PAGEIN) LEVV=LEV
 1040    CONTINUE
         IF(LEV.GT.LEVFRM) GO TO 1050
         IF((NABCA(NPT)*2+NABC(NPT)).LE.2) GO TO 1060
 1050    IF(LEV.GT.LEVI) GO TO 1000
         GO TO 1070
 1060    CONTINUE
C   LOOPS WITH FOUR EXTERNAL INDICES
         IF(NDVDIT.EQ.1.AND.IGUESS.EQ.0) GO TO 480
         ITYPE=3*NABCA(NPT)+NABC(NPT)
         IF(ITYPE.EQ.0) GO TO 1010
         IFSYM=NABCS(NPT)
         JFSYM=IFSYM
         LVFRM1=LEV-1
         IMAX=LEVI-1
         IMIN=LEVI-1
         IF(.NOT.PAGEIN) THEN
                WRITE(ITAPE6,*) ' PAGING MAFUNCTION'
                CALL MABORT
         END IF
         CALL ALLEXT(int,sigma,ijadd,kadd,ladd,wtw,wtx,wty,wab,isym,
     1 vector,ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww,klxx,klww)
         GO TO 1010
 1070    LEV=LEVI
         LEVM=LEV-1
C     PRINT *,'AFTER 1070, LEV,LEVM:',LEV,LEVM
         ISEGM(LEV)=1
         ISEG=1
         IMN=NPT
         ISB=NPT
         KSEG=0
C     PRINT *,'KSEG=0 BEFORE 200'
         KSEGMX=JSEGNR(ISEG)
         LMIN(LEV)=LCOND(JSM+1)
         IUWKMN(LEV)=0
         IUWKSB(LEV)=0
         IMAIN(LEV)=NPT
         ISUB(LEV)=NPT
         ACOEF(LEV)=1.0D+00
C  TEST NEXT SEGMENT OF GROUP
200      KSEG=KSEG+1
C     PRINT *,' KSEG',KSEG,'LEV',LEV,'KSEGMX',KSEGMX
C     WRITE (ITAPE3,'(' KSEG:',I4,' LEV',I3)') KSEG,LEV
         IF(KSEG.GT.KSEGMX) GO TO 440
         KSB=IARCSB(KSEG)
         JARPT=4*(ISB-1)+KSB
         KSB=IARC(JARPT)
         IF(KSB.EQ.0) GO TO 200
         KMN=IARCMN(KSEG)
         IARPT=4*(IMN-1)+KMN
         KMN=IARC(IARPT)
         IF(KMN.EQ.0) GO TO 200
         IF(NDVDIT.GT.0) GO TO 210
         IF(KMN.NE.KSB) GO TO 200
 210     CONTINUE
         JSEGM(LEV)=KSEG
         IUWKMN(LEVM)=IUWKMN(LEV)+IWGHT(IARPT)
         IUWKSB(LEVM)=IUWKSB(LEV)+IWGHT(JARPT)
         LMIN(LEVM)=LMIN(LEV)
         IF(JCOND(KSEG))220,240,230
 220     CONTINUE
C
C      ----- IS THIS A BUG? PWS -----
C
C     IF(LEVM.LE.JMIN)GOTO 440
         IF (LEVM.LE.JMIN) GO TO 200
C
         GOTO 240
 230     CONTINUE
         IF (LEVM.GT.JMAX) GO TO 200
C     IF(LEVM.GT.JMAX)GOTO 420
         IJ=IAD+LEVM
         JAD=IJADD(IJ)
         KSM=IEOR(JSM,ISYM(LEVM))
         LMIN(LEVM)=LCOND(KSM+1)
         KSMPTX=ISMOFF(KSM+1)
 240     CONTINUE
         IF(KCOND(KSEG).EQ.0)GOTO 260
         KSMPT=LEVM+KSMPTX
         KAD=JAD+KADD(KSMPT)
         LSM=IEOR(KSM,ISYM(LEVM))
         LMIN(LEVM)=LCOND(LSM+1)
         LSMPTX=ISMOFF(LSM+1)
 260     CONTINUE
         IF(ITRK(KSEG))280,280,270
270      ITRACK(LEVM)=ITRK(KSEG)
         GOTO 290
280      ITRACK(LEVM)=ITRACK(LEV)
290      CONTINUE
         GO TO( 1, 1, 1, 1, 1, 1,44, 1, 3, 3, 1,45, 4, 4,50,51, 1,40, 1,
     *    1,     6, 6, 7,46, 9,54, 2,52,41, 5,47, 8, 2,53, 1, 1,42,48, 2
     *   ,52,    36,55,43, 1,11,11,12,10,13,49, 2,53, 1,77,77,77, 1,79,
     *   77, 1,    80,78, 1,71,67,68,87,75,83,69,68,76,70,82,71,71,67,68
     *   ,67,87,    75,83,69,68,83,68,76,70,69,70,82,71, 1, 6,16, 6, 6,
     *   17,16,74,     8, 1, 1,18,19,18,18,22,24,20,19,24,19,23,21,20,21
     *   , 1, 1,11,    11,27,11,28,81,27,13, 1, 1, 3, 2, 4, 2, 1, 2, 2,
     *   1,71,63,72,    84,65,85,73,29,66,64,71, 1,30,56,86, 1,57, 1, 1,
     *   56, 1,32,31,    58, 1, 1,59,33,34,61,35,88,62,60, 1, 1, 6, 1, 9
     *   , 2, 5, 2, 1,     1, 2,36,11,10, 2, 1, 6, 1, 9, 2, 5, 2, 1, 1,
     *   2,36,11,10, 2,     1, 1,37,38, 4, 2, 2, 1, 1, 2,39,37, 3, 2, 1,
     *    6, 1, 9, 2, 5,     2, 1, 1, 2,36,11,10, 2) KSEG
  1      ACOEF(LEVM)=ACOEF(LEV)
         GOTO 120
  2      ACOEF(LEVM)=-ACOEF(LEV)
         GOTO 120
   3     IA = NABC(IMN) + 2
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
   4     IA = NABC(IMN) + 83
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
   5     IA = NABC(IMN) + 82
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
   6     IA = NABC(IMN) + 261
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
   7     IA = NABC(IMN) + 1
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
   8     IA = NABC(IMN) + 102
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
   9     IA = NABC(IMN) + 362
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  10     IA = NABC(IMN) + 3
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  11     IA = NABC(IMN) + 263
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  12     IA = NABC(IMN) + 84
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  13     IA = NABC(IMN) + 23
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  16     IA = NABC(IMN) + 281
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  17     IA = NABC(IMN) + 402
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  18     IA = NABC(IMN) + 162
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  19     IA = NABC(IMN) + 222
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  20     IA = NABC(IMN) + 143
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  21     IA = NABC(IMN) + 42
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  22     IA = NABC(IMN) + 302
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  23     IA = NABC(IMN) + 303
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  24     IA = NABC(IMN) + 342
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  27     IA = NABC(IMN) + 283
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  28     IA=NABC(IMN)+404
         ACOEF(LEVM)=ACOEF(LEV)*CFS(IA)
         GO TO 120
  29     ACOEF(LEVM) = ACOEF(LEV) * ROOT2
         GO TO 120
  30     IA = NABC(IMN) + 301
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  31     IA = NABC(IMN) + 304
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  32     IA = NABC(IMN) + 244
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  33     IA = NABC(IMN) + 322
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  34     IA = NABC(IMN) + 243
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  35     IA = NABC(IMN) + 242
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  36     IA = NABC(IMN) + 384
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  37     IA = NABC(IMN) + 262
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  38     IA = NABC(IMN) + 363
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  39     IA = NABC(IMN) + 383
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  86     IA = NABC(IMN) + 241
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  40     IA = NABC(IMN) + 122
         IB=IA-61
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         BCOEF(LEVM) = ACOEF(LEV) * CFS(IB)
         GO TO 120
  41     IB = NABC(IMN) + 162
         ACOEF(LEVM) = ACOEF(LEV) * TOORN2
         BCOEF(LEVM) = ACOEF(LEV) * CFS(IB)
         GO TO 120
  42     IA = NABC(IMN) + 43
         IB = IA + 81
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         BCOEF(LEVM) = ACOEF(LEV) * CFS(IB)
         GO TO 120
  43     IB = NABC(IMN) + 222
         ACOEF(LEVM) = ACOEF(LEV) * TOORN2
         BCOEF(LEVM) = ACOEF(LEV) * CFS(IB)
         GO TO 120
  44     IB=NABC(IMN)+221
         ACOEF(LEVM) = ACOEF(LEV) * TOOR2
         BCOEF(LEVM) = ACOEF(LEV) * CFS(IB)
         GO TO 120
  45     IB = NABC(IMN) + 163
         ACOEF(LEVM) = ACOEF(LEV) * TOOR2
         BCOEF(LEVM) = ACOEF(LEV) * CFS(IB)
         GO TO 120
  46     IB = NABC(IMN) + 162
         ACOEF(LEVM) = ACOEF(LEV) * TOOR2
         BCOEF(LEVM) = ACOEF(LEV) * CFS(IB)
         GO TO 120
  47     IA = NABC(IMN) + 122
         IB = IA - 81
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         BCOEF(LEVM) = ACOEF(LEV) * CFS(IB)
         GO TO 120
  48     IB = NABC(IMN) + 222
         ACOEF(LEVM) = ACOEF(LEV) * TOOR2
         BCOEF(LEVM) = ACOEF(LEV) * CFS(IB)
         GO TO 120
  49     IA = NABC(IMN) + 43
         IB = IA + 101
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         BCOEF(LEVM) = ACOEF(LEV) * CFS(IB)
         GO TO 120
  50     ACOEF(LEVM) = ACOEF(LEV) + ACOEF(LEV)
         D=0.5D+00
         GO TO 120
  51     ACOEF(LEVM)=ACOEF(LEV)*ROOT2
         GO TO 120
  52     ACOEF(LEVM) = -ACOEF(LEV)
         D= -1.0D+00
         GO TO 120
  53     ACOEF(LEVM) = -ACOEF(LEV) - ACOEF(LEV)
         D = -0.5D+00
         GO TO 120
  54     IA=NABC(IMN)+362
         D=1.0D+00/CFS(IA)
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  55     IA = NABC(IMN) + 384
         D=1.0D+00/CFS(IA)
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
  56     ACOEF(LEVM) = ACOEF(LEV)
         D = -1.0D+00
         GO TO 120
  57     IA = NABC(IMN) + 82
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         D=-1.0D+00
         GO TO 120
  58     IA = NABC(IMN) + 3
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         D=-1.0D+00
         GO TO 120
  59     IA = NABC(IMN) + 123
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         D=-1.0D+00
         GO TO 120
  60     IA = NABC(IMN) + 222
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         D=-1.0D+00
         GO TO 120
  61     IA = NABC(IMN) + 62
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         D=-1.0D+00
         GO TO 120
  62     IA = NABC(IMN) + 162
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         D=-1.0D+00
         GO TO 120
  63     IA = NABC(IMN) + 42
         IB = IA + 81
         ACOF = ACOEF(LEV) * CFS(IA)
         BCOF = BCOEF(LEV) * CFS(IB)
         D = ACOF + BCOF
         IF(ABS(D).LT.CRITE) GO TO 110
         ACOEF(LEVM) = D
         D = (ACOF - BCOF) / D
         GO TO 120
  64     IB = NABC(IMN) + 222
         ACOF = ACOEF(LEV) * TOORN2
         BCOF = BCOEF(LEV) * CFS(IB)
         D = ACOF + BCOF
         IF(ABS(D).LT.CRITE) GO TO 110
         ACOEF(LEVM) = D
         D = (ACOF - BCOF) / D
         GO TO 120
  65     IA = NABC(IMN) + 123
         IB = IA - 61
         ACOF = ACOEF(LEV) * CFS(IA)
         BCOF = BCOEF(LEV) * CFS(IB)
         D = ACOF + BCOF
         IF(ABS(D).LT.CRITE) GO TO 110
         ACOEF(LEVM) = D
         D = (ACOF - BCOF) / D
         GO TO 120
  66     IB = NABC(IMN) + 162
         ACOF = ACOEF(LEV) * TOORN2
         BCOF = BCOEF(LEV) * CFS(IB)
         D = ACOF + BCOF
         IF(ABS(D).LT.CRITE) GO TO 110
         ACOEF(LEVM) = D
         D = (ACOF - BCOF) / D
         GO TO 120
  67     IB = NABC(IMN) + 162
         DX=ACOEF(LEV)*TOORN2
         D=DX+BCOEF(LEV)*CFS(IB)
         IF(ABS(D).LT.CRITE) GO TO 111
         ACOEF(LEVM) = D
         D=-(DX+DX)/D
         GO TO 120
  68     IB = NABC(IMN) + 222
         DX=ACOEF(LEV)*TOORN2
         D=DX+BCOEF(LEV)*CFS(IB)
         IF(ABS(D).LT.CRITE) GO TO 111
         ACOEF(LEVM) = D
         D=-(DX+DX)/D
         GO TO 120
  69     IA = NABC(IMN) + 62
         IB = IA + 81
         DX=ACOEF(LEV)*CFS(IA)
         D=DX+BCOEF(LEV)*CFS(IB)
         IF(ABS(D).LT.CRITE) GO TO 111
         ACOEF(LEVM) = D
         D=-(DX+DX)/D
         GO TO 120
  70     IA = NABC(IMN) + 143
         IB = IA - 101
         DX=ACOEF(LEV)*CFS(IA)
         D=DX+BCOEF(LEV)*CFS(IB)
         IF(ABS(D).LT.CRITE) GO TO 111
         ACOEF(LEVM) = D
         D=-(DX+DX)/D
         GO TO 120
  87     IB = NABC(IMN) + 162
         DX=ACOEF(LEV)*TOORN2
         D=DX+BCOEF(LEV)*CFS(IB)
         IF(ABS(D).LT.CRITE) GO TO 111
         ACOEF(LEVM)=D
         D=-(DX+DX)/D
         GO TO 120
  71     ACOEF(LEVM) = ACOEF(LEV)
         BCOEF(LEVM) = BCOEF(LEV)
         GO TO 120
  72     IB = NABC(IMN) + 322
         ACOEF(LEVM) = -ACOEF(LEV)
         BCOEF(LEVM) = BCOEF(LEV) * CFS(IB)
         GO TO 120
  73     IB = NABC(IMN) + 323
         ACOEF(LEVM) = -ACOEF(LEV)
         BCOEF(LEVM) = BCOEF(LEV) * CFS(IB)
         GO TO 120
  74     IA=NABC(IMN)+21
         ACOEF(LEVM)=ACOEF(LEV)*CFS(IA)
         GO TO 120
  75     IB = NABC(IMN) + 302
         ACOEF(LEVM) = ACOEF(LEV)
         BCOEF(LEVM) = BCOEF(LEV) * CFS(IB)
         GO TO 120
  76     IB = NABC(IMN) + 303
         ACOEF(LEVM) = ACOEF(LEV)
         BCOEF(LEVM) = BCOEF(LEV) * CFS(IB)
         GO TO 120
  77     ACOEF(LEVM)=ACOEF(LEV)*TOORN2
         D=-2.0D+00
         GO TO 120
  78     ACOEF(LEVM)=ACOEF(LEV)*ROOTN2
         D=-2.0D+00
         GO TO 120
  79     IA=NABC(IMN)+62
         ACOEF(LEVM)=ACOEF(LEV)*CFS(IA)
         D=-2.0D+00
         GO TO 120
  80     IA=NABC(IMN)+143
         ACOEF(LEVM)=ACOEF(LEV)*CFS(IA)
         D=-2.0D+00
         GO TO 120
  81     IA=NABC(IMN)+104
         ACOEF(LEVM)=ACOEF(LEV)*CFS(IA)
         GO TO 120
  82     ACOEF(LEVM) = ACOEF(LEV) * ROOTN2
         D=-2.0D+00
         GO TO 120
  83     IA = NABC(IMN) + 342
         ACOEF(LEVM) = BCOEF(LEV) * CFS(IA)
         GO TO 120
  84     IA = NABC(IMN) + 243
         ACOEF(LEVM) = BCOEF(LEV) * CFS(IA)
         GO TO 120
  85     IA = NABC(IMN) + 242
         ACOEF(LEVM) = BCOEF(LEV) * CFS(IA)
         GO TO 120
  88     IA = NABC(IMN) + 323
         ACOEF(LEVM) = ACOEF(LEV) * CFS(IA)
         GO TO 120
 110     ITRACK(LEVM)=3
         ACOEF(LEVM)=ACOF-BCOF
         GO TO 120
111      ITRACK(LEVM) = 2
         ACOEF(LEVM)=-(DX+DX)
  120    CONTINUE
C   CHECK TO SEE IF PAGE OF VECTOR IS CORRECT
 1100    IF(PAGEIN) GO TO 1110
         IUWK=IUWKSB(LEVM)
         JUWK=IUWKMN(LEVM)
         NLWKI=NLWKS(KSB)
         NLWKJ=NLWKS(KMN)
         CALL PAGES(vector,sigma)
         IF(PAGEIN) LEVV=LEVM
         GO TO 1130
 1110    CONTINUE
         IF(LEVM.LT.LEVV) GO TO 1130
         IUWK=IUWKSB(LEVM)
         JUWK=IUWKMN(LEVM)
         NLWKI=NLWKS(KSB)
         NLWKJ=NLWKS(KMN)
         IF(LEVM.GT.LEVV) GO TO 1120
         CALL PAGEIJ(vector,sigma)
         LEVV=NLEVS
         IF(PAGEIN) LEVV=LEVM
         GO TO 1130
 1120    CONTINUE
         CALL PAGEOT(vector,sigma)
         GO TO 1100
 1130    CONTINUE
         IF(NXTSEG(KSEG).GT.0) GO TO 400
         IF(ISYM(LEVM).NE.LSM) GO TO 200
         LSMPT=LEVM+LSMPTX
         LAD=KAD+LADD(LSMPT)
         IF(KMN-KSB.EQ.0) GO TO 380
         LEVL=LEVM
         KSEGMX=4
310      LEV=LEVM
         LEVM=LEV-1
         IF(LEVM.GT.0) GO TO 315
         WRITE(ITAPE6,313)
 313     FORMAT(' PROBLEMS WITH PARTIAL SPACE')
         CALL MABORT
 315     CONTINUE
         KSEG=0
C     PRINT *,'KSEG=0 BEFORE 320'
         IMAIN(LEV)=KMN
         IMN=KMN
         ISUB(LEV)=KSB
         ISB=KSB
320      KSEG=KSEG+1
C     PRINT *,'AFTER 320, KSEG,KSEGMX,LEV',KSEG,KSEGMX,LEV
         IF(KSEG.GT.KSEGMX)GOTO 360
         JARPT=4*(ISB-1)+KSEG
         KSB=IARC(JARPT)
         IF(KSB.LE.0)GOTO 320
         IARPT=4*(IMN-1)+KSEG
         KMN=IARC(IARPT)
         IF(KMN.LE.0)GOTO 320
         JSEGM(LEV)=KSEG
         IUWKMN(LEVM)=IUWKMN(LEV)+IWGHT(IARPT)
         IUWKSB(LEVM)=IUWKSB(LEV)+IWGHT(JARPT)
         IF(LEVM.GT.LEVFRM) GOTO 310
         IF(KMN-KSB.NE.0) THEN
           WRITE(ITAPE6,*) ' KMN-KSB.NE.0'
           CALL MABORT
         END IF
         NLWKI = NLWKS(KMN)
         IUWK = IUWKSB(LEVM)
         JUWK = IUWKMN(LEVM)
         ITRAK = ITRACK(LEVL)
         ACF = ACOEF(LEVL)
         LADT=LAD
         CALL LOOPIN(int,sigma,ijadd,kadd,ladd,wtw,wtx,wty,wab,isym,
     1    vector,ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww,klxx,klww)
         GO TO 320
  360    CONTINUE
C     PRINT *,'LEV,LEVL AT 360',LEV,LEVL
         IF(LEV.EQ.LEVL) GO TO 440
         LEVM=LEV
         LEV=LEVM+1
         IMN=IMAIN(LEV)
         ISB=ISUB(LEV)
         KSEG=JSEGM(LEV)
C     PRINT *,'KSEG SET TO BEFORE 380',KSEG
         GO TO 320
 380     NLWKI = NLWKS(KMN)
         IUWK = IUWKSB(LEVM)
         JUWK = IUWKMN(LEVM)
         ITRAK = ITRACK(LEVM)
         ACF = ACOEF(LEVM)
         LADT = LAD
         CALL LOOPIN(int,sigma,ijadd,kadd,ladd,wtw,wtx,wty,wab,isym,
     1    vector,ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww,klxx,klww)
         GO TO 200
  400    CONTINUE
         IF(LEVM.LE.LMIN(LEVM))GOTO 200
         ISEG=NXTSEG(KSEG)
         IF(NLCSMN(ISEG).GT.0.AND.(NABCA(KMN)+NABC(KMN)).EQ.0)GOTO 200
         IF(LEVM.LE.LEVFRM)GOTO 460
  410    LEV=LEVM
         LEVM=LEV-1
C     PRINT *,'AFTER 410, GOING DOWN    ,LEV=',LEV
         ISEGM(LEV)=ISEG
         KSEG=JSEGPT(ISEG)
C     PRINT *,'BEFORE 420, KSEG SET TO',KSEG
         IMN=KMN
         IMAIN(LEV)=KMN
         ISB=KSB
         ISUB(LEV)=KSB
         KSEGMX=JSEGNR(ISEG)
         GO TO 200
440      CONTINUE
C   PATHS EXHAUSTED AT THIS LEVEL, BACK UP AND TRY THE NEXT SEGMENT
C     PRINT *,' AFTER 440, LEV AND LEVI',LEV,LEVI
         IF(LEV.EQ.LEVI) GOTO 1010
C     PRINT *,' PAST TEST, BACK UP LEVEL, LEV=',LEV+1
         LEVM=LEV
         LEV=LEVM+1
         ISEG=ISEGM(LEV)
         IMN=IMAIN(LEV)
         ISB=ISUB(LEV)
         KSEG=JSEGM(LEV)
C     PRINT *,'BEFORE 460, KSEG SET TO',KSEG
         KSEGMX=JSEGNR(ISEG)
         GOTO 200
460      CONTINUE
C  FINISH LOOPS WITH 3 OR FEWER EXTERNAL INDICES
         IAMN=NABCA(KMN)
         IASB=NABCA(KSB)
         IBMN=NABC(KMN)
         IBSB=NABC(KSB)
         IPT1=IAMN+IAMN+IBMN
         IPT2=IASB+IASB+IBSB
         IF(IPT1.GT.2.OR.IPT2.GT.2) GO TO 410
         IPT1=4-IPT1-IAMN
         IPT2=4-IPT2-IASB
         LVFRM1=LEVM-1
         IUWK=IUWKSB(LEVM)
         JUWK=IUWKMN(LEVM)
         ITRAK=ITRACK(LEVM)
         JFSYM=NABCS(KMN)
         IFSYM=NABCS(KSB)
         ACF=ACOEF(LEVM)
         CCF=BCOEF(LEVM)
         ISEGT=ISEG
         IF(ISEG.GT.14) GO TO 470
         IF(ISEG.GT.3)  GO TO 465
         JSMT=JSM
         LADT=IAD
         GOTO 475
 465     CONTINUE
         JSMT=KSM
         LADT=JAD
         GO TO 475
 470     CONTINUE
         LADT=KAD
         JSMT=LSM
 475     CALL LOOPEX(int,sigma,ijadd,kadd,ladd,wtw,wtx,wty,wab,isym,
     1 vector,ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww,klxx,klww)
         GOTO 200
 480     CONTINUE
         IF(PAGEIN.AND.NWKS.GT.NWKSMX) CALL PAGEOT(vector,sigma)
C:    INXT=INXT+1
cets110990         IF(IRSTRT.NE.0) CALL SAVE(sigma)
C:    IF(INXT.LE.NORBS) GO TO 190
         I=I-1
C     PRINT *,'I AND IMING',I,IMING
         IF (I.GE.IMING) GO TO 190
C:
 500  CONTINUE
C::::
      IF (PAGEIN) CALL PAGEOT(vector,sigma)
C::::
      IBLCK1=1
      INXT=1
      RETURN
      END
