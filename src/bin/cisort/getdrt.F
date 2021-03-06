      SUBROUTINE GETDRT(    KADD,LADD,IJADD,IJGRP
     #,                        NINGRP,ORBSYM,NBLKMN,NBLKMX
     #,                 IJXX,KLXX,NKLXX,IJWW,KLWW,NKLWW,END)
C
      IMPLICIT INTEGER (A-Z)
C
CTJL  EXTENDED DUMMY KADD,LADD,IJADD,IJGRP,NINGRP,ORBSYM,NBLKMN,NBLKMX
CTJL  EXTENDED DUMMY IJXX,KLXX,NKLXX,IJWW,KLWW,NKLWW
C
      REAL*8 VER4X,DRTVER
C
      COMMON /LBLS/ LBLINT(26),LBLDRT(26),VER4X,DRTVER
      COMMON /TAPES/ITAP52,ITAPE5,ITAPE6,ITAP58,ITAP12,ITAP99,ITAP04
     *,             ITAPE3,ITAP05,ITAP06
      COMMON /DIMS/ NBF,NSYM,NORBS,NROWS,NROWS4,NWKS,NWKS2,NLEVS
     *,             NROWOC,NROW4O,NWKSOC,NLEVOC,NORBOC,LEVFRM
     *,             NWKSMX,NLWKMX,NUWKMX,MAXB,NROOTS,LVFRM1,NREFS
      COMMON /INTS/ NMAX,NMAX2,NGROUP,NBLKOC,NUMIJ,SYMORB,INTSRT
      COMMON /XX4X/ NIJVIR
C
      DIMENSION IJADD(NUMIJ),KADD(SYMORB),LADD(SYMORB)
      DIMENSION IJGRP(NUMIJ)
      DIMENSION ORBSYM(NORBS)
      DIMENSION NINGRP(NGROUP),NBLKMN(NGROUP),NBLKMX(NGROUP)
      DIMENSION IJXX(NUMIJ ),KLXX(NIJVIR),NKLXX(NSYM)
      DIMENSION IJWW(NUMIJ ),KLWW(NIJVIR),NKLWW(NSYM)
C
      CALL WREADW(ITAP58,KADD,SYMORB,END,END)
      CALL WREADW(ITAP58,LADD,SYMORB,END,END)
      CALL WREADW(ITAP58,IJADD,NUMIJ,END,END)
      CALL WREADW(ITAP58,IJGRP,NUMIJ,END,END)
      END=END+NORBS*4
      CALL WREADW(ITAP58,NINGRP,NGROUP,END,END)
      END=END+NBF
      CALL WREADW(ITAP58,ORBSYM,NORBS,END,END)
      END=END+NBF+NREFS*NBF+NORBS+NLEVS*2+NROWS*12+LVFRM1*2*(NSYM+1)
      CALL WREADW(ITAP58,IJXX,NUMIJ,END,END)
      CALL WREADW(ITAP58,KLXX,NIJVIR,END,END)
      CALL WREADW(ITAP58,NKLXX,NSYM,END,END)
      END=END+NSYM*(LVFRM1-1)
      CALL WREADW(ITAP58,IJWW,NUMIJ,END,END)
      CALL WREADW(ITAP58,KLWW,NIJVIR,END,END)
      CALL WREADW(ITAP58,NKLWW,NSYM,END,END)
      END=END+NSYM*(LVFRM1-1)
C
C
C
C
      IF (IDBG9 .EQ. 0) GOTO 9
      WRITE (ITAPE6,*) ' KADD '
      WRITE (ITAPE6,710) (KADD(I),I=1,SYMORB)
      WRITE (ITAPE6,*) ' LADD '
      WRITE (ITAPE6,710) (LADD(I),I=1,SYMORB)
      WRITE (ITAPE6,*) ' IJADD'
      WRITE (ITAPE6,710) (IJADD(I),I=1,NUMIJ)
      WRITE (ITAPE6,*) ' IJGRP'
      WRITE (ITAPE6,710) (IJGRP(I),I=1,NUMIJ)
      WRITE (ITAPE6,*) ' NINGRP'
      WRITE (ITAPE6,710) (NINGRP(I),I=1,NGROUP)
      WRITE (ITAPE6,*) ' ORBSYM'
      WRITE (ITAPE6,710) (ORBSYM(I),I=1,NORBS)
      WRITE (ITAPE6,*) ' IJXX '
      WRITE (ITAPE6,710) (IJXX(I),I=1,NUMIJ)
      WRITE (ITAPE6,*) ' KLXX '
      WRITE (ITAPE6,710) (KLXX(I),I=1,NIJVIR)
      WRITE (ITAPE6,*) ' NKLXX'
      WRITE (ITAPE6,710) (NKLXX(I),I=1,NSYM)
      WRITE (ITAPE6,*) ' IJWW '
      WRITE (ITAPE6,710) (IJWW(I),I=1,NUMIJ)
      WRITE (ITAPE6,*) ' KLWW '
      WRITE (ITAPE6,710) (KLWW(I),I=1,NIJVIR)
      WRITE (ITAPE6,*) ' NKLWW'
      WRITE (ITAPE6,710) (NKLWW(I),I=1,NSYM)
C
C
C
  710 FORMAT (20I5)
    9 CONTINUE
      IBLK=1
      RETURN
      END
