      SUBROUTINE INITPG(VECTOR)
      IMPLICIT REAL*8 (A-H,O-Z)
C      EXTENDED DUMMY VECTOR
      REAL*8 VECTOR,ACF,D,CCF
      integer nwksmx
      DIMENSION VECTOR(NWKSMX)
      INTEGER PUWKT,REFWLK,BMAX,FWORD,PUWK,ORBFRM,IJUNK(26)
      LOGICAL PAGEIN,IEQJ
C
      COMMON /DIMS/ NBF,NSYM,NORBS,NROWS,NROWS4,NWKS,NWKS2,NLEVS
     *,             NROWOC,NROW4O,NWKSOC,NLEVOC,NORBOC,LEVFRM
     *,             NWKSMX,NLWKMX,NUWKMX,BMAX,NROOTS,ORBFRM
       COMMON /DIAG/ REP,FZCORE,EGUESS,ECI,REFWLK,MXITER,CNVERG,ICNVG
     *,             ITER,SQCDIF,CZERO,NROOT
CTJL  COMMON /DIAG/ REP,FZCORE,EGUESS,ECI,SQCDIF,CNVERG,CZERO
CTJL *,             MXITER,ITER,ICNVG,NROOT,REFWLK
      COMMON /TAPES/ITAP52,ITAPE5,ITAPE6,ITAP58,ITAP54,ITAP03,ITAP04
     *,             ITAPE3,ITAP05,ITAP06,ITAP56
      COMMON /OPS/  IGUESS,IRSTRT,IROOTI,IROOTF
      COMMON /LOOPS/NUWK,PUWKT,IUWK,JUWK,ITRAK,IPT1,IPT2
      COMMON /ALL/  ACF,D,CCF,LADT,ITR1,ITR2,IA,JA,ITYPE,ISEGT,LVFRM1
     *,             NLWKI,NLWKJ,IMAX,IMIN
      COMMON /PAGE/ IWORD3,IWORD4,PAGEIN,NOFFI,NOFFJ,NPASS,NWORDL,IWORD6
      COMMON /PAGDAT/IEQJ,IUWKSV,JUWKSV,NLWKI2,NLWKJ2,JWALK,IWALK,II,JJ
CWA
      IWORD3=IROUND(112)+1+(IGUESS-1)*NWKS2
CTJL  ITAP03=12
      ITAP03=54
      NPASS=(NWKS-1)/NWKSMX+1
      NWKMX2=INTOWP(NWKSMX)
      NWORDL=NWKS2-(NPASS-1)*NWKMX2
      IEQJ=.TRUE.
      IWALK=0
      II=1
      JJ=1
C**********************************************************************
C IS THIS A BUG????????????????????????????? TJL,JER 11/17/84
C     NLWKI2=INTOWP(NWLKS2)
      NLWKI2=INTOWP(NWKS2)
C**********************************************************************
C ALSO IT APPEARS NLWKI (WHICH IS USED IN ENTRY PAGED) IS NEVER GIVEN
C A VALUE. TJL,JER 11/17/84  (THIS MAY BE WRONG)
      NLWKI=NWKS2
C**********************************************************************
C     WRITE(*,*) ' NLWKI2,NWKS2,NLWKI= ',NLWKI2,NWKS2,NLWKI
C**********************************************************************
      NWKMX2=NWKSMX
C   READ IN OLD VECTOR FROM ITAP54 IF PAGING IS NOT NEEDED
      WRITE(ITAPE3,110)
 110  FORMAT(1X,'READING IN CI VECTOR FROM TAPE 54 ')
CIBM  CALL SREW(ITAP54)
CIBM  CALL DCDLBL(ITAP54,IJUNK,IJUNK,REFWLK,IJUNK,NREFS,EGUESS
CIBM *,           ECI,CZERO)
CTJL  CALL SREW(ITAP56)
CTJL  CALL DCDLBL(ITAP56,IJUNK,IJUNK,REFWLK,IJUNK,NREFS,EGUESS
CTJL *,           ECI,CZERO)
      WRITE(ITAPE6,115) REFWLK,ECI
 115  FORMAT(1X,'REFERENCE   = ',I8,/,1X,'CI ENERGY   = ',G17.11,/)
      NWORD=NWKS2
C**********************************************************************
C     WRITE(*,*) ' READING CI VECTOR, NWKS,NWKSMX=',NWKS,NWKSMX
C     WRITE(*,*) ' NWORD,IWORD3,JUNK = ',NWORD,IWORD3,JUNK
C**********************************************************************
      IF(NWKS.LE.NWKSMX) CALL WREADW(ITAP54,VECTOR,NWORD,IWORD3
     -,JUNK)
CMC   PRINT,' MWKS,NWKSMX',NWKS,NWKSMX
CTJL  WRITE(6,*) 'CI VECTOR READ IN INITPG'
CTJL  DO 9011 IMIKE=1,NWKS
CTJL       WRITE(6,9010) IMIKE,VECTOR(IMIKE)
C9010      FORMAT(I5,E20.12)
C9011 CONTINUE
      RETURN
C
C-----------------------------------------------------------------PAGES
C
      ENTRY PAGES(VECTOR)
CMC   WRITE(6,*) ' ENTRY PAGES'
      IF(IUWK.NE.JUWK) GO TO 210
C
C-----------------------------------------------------------------PAGED
C
      ENTRY PAGED(VECTOR)
CMC   WRITE(6,*) 'ENTRY PAGED'
C**********************************************************************
C     WRITE(*,*) ' NLWKI = ',NLWKI
C**********************************************************************
  200 IF(NLWKI.GT.NWKSMX) RETURN
      IUWKSV=IUWK
      JUWKSV=JUWK
      IEQJ=.TRUE.
      PAGEIN=.TRUE.
      NLWKI2=NLWKI+NLWKI
      NLWKJ2=NLWKI2
      NOFFI=PUWKT+IUWK-1
      NOFFJ=NOFFI
      IWALK=NOFFI+NOFFI
      II=1
      JJ=1
      FWORD=IWALK+IWORD3
      CALL WREADW(ITAP03,VECTOR,NLWKI2,FWORD,JUNK)
      RETURN
  210 CONTINUE
      IF((NLWKI+NLWKJ).GT.NWKSMX) RETURN
C
C-----------------------------------------------------------------PAGEO
C
      ENTRY PAGEO(VECTOR)
CMC   WRITE(6,*) 'ENTRY PAGEO'
C
      IUWKSV=IUWK
      JUWKSV=JUWK
      PAGEIN=.TRUE.
      IEQJ=.FALSE.
      NLWKI2=NLWKI+NLWKI
      NLWKJ2=NLWKJ+NLWKJ
      NOFFI=PUWKT+IUWK-1
      NOFFJ=PUWKT+JUWK-1
      IWALK=NOFFI+NOFFI
      JWALK=NOFFJ+NOFFJ
      NOFFJ=NOFFJ-NLWKI
      II=1
      JJ=NLWKI+1
      FWORD=IWALK+IWORD3
      CALL WREADW(ITAP03,VECTOR,NLWKI2,FWORD,JUNK)
      FWORD=JWALK+IWORD3
      CALL WREADW(ITAP03,VECTOR(JJ),NLWKJ2,FWORD,JUNK)
      RETURN
C
C----------------------------------------------------------------PAGEIJ
C
      ENTRY PAGEIJ(VECTOR)
CMC   WRITE(6,*) 'ENTRY PAGEIJ'
C
C   CHECK THAT SOMETHING HAS CHANGED
C
      IF(IUWKSV.EQ.IUWK.AND.JUWKSV.EQ.JUWK) RETURN
      IF(IEQJ) GO TO 400
      IF(IUWKSV.EQ.IUWK) GO TO 310
      IF(JUWKSV.EQ.JUWK) GO TO 330
  300 PAGEIN=.FALSE.
      IF(IUWK.EQ.JUWK) GO TO 200
      GO TO 210
  310 CONTINUE
      IF((NLWKI+NLWKJ).LE.NWKSMX) GO TO 320
      PAGEIN=.FALSE.
      RETURN
  320 CONTINUE
      JJ=NLWKI+1
      JUWKSV=JUWK
      NLWKJ2=NLWKJ+NLWKJ
      NOFFJ=PUWKT+JUWK-1
      JWALK=NOFFJ+NOFFJ
      NOFFJ=NOFFJ-NLWKI
      FWORD=JWALK+IWORD3
      CALL WREADW(ITAP03,VECTOR(JJ),NLWKJ2,FWORD,JUNK)
      RETURN
  330 CONTINUE
      IF(IUWK.NE.JUWK) GO TO 340
      IEQJ=.TRUE.
      II=JJ
      IWALK=JWALK
      IUWKSV=IUWK
      NLWKI2=NLWKJ2
      RETURN
  340 CONTINUE
      IF(NLWKI.LT.JJ) GO TO 350
      PAGEIN=.FALSE.
      GO TO 210
  350 IUWKSV=IUWK
      NOFFI=PUWKT+IUWK-1
      IWALK=NOFFI+NOFFI
      NLWKI2=NLWKI+NLWKI
      FWORD=IWALK+IWORD6
      II=1
      FWORD=IWALK+IWORD3
      CALL WREADW(ITAP03,VECTOR(II),NLWKI2,FWORD,JUNK)
      RETURN
  400 CONTINUE
      IF(IUWK.NE.IUWKSV) GO TO 420
      IEQJ=.FALSE.
      IF(II.EQ.1) GO TO 405
      IF(NLWKJ.GE.II) GO TO 420
      JJ=1
      GOTO 410
  405 CONTINUE
      IF((NLWKI+NLWKJ).GT.NWKSMX) GO TO 420
      JJ=II+NLWKI
  410 JUWKSV=JUWK
      NOFFJ=PUWKT+JUWK-1
      JWALK=NOFFJ+NOFFJ
      NOFFJ=NOFFJ-NLWKI
      NLWKJ2=NLWKJ+NLWKJ
      FWORD=JWALK+IWORD3
      CALL WREADW(ITAP03,VECTOR(JJ),NLWKJ2,FWORD,JUNK)
      RETURN
  420 CONTINUE
      PAGEIN=.FALSE.
      IF(IUWK.EQ.JUWK) GO TO 200
      GO TO 210
C
      ENTRY PAGEOU(VECTOR)
CMC   WRITE(6,*) 'ENTRY PAGEOUT'
C
      PAGEIN=.FALSE.
      RETURN
      END
