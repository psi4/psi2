      SUBROUTINE x3ARR(ORBSYM,FLOV,NIRRED,XOFF,XADD,no)
      IMPLICIT INTEGER (A-Z)
      DIMENSION FLOV(NIRRED,4),ORBSYM(*)  
      INTEGER XOFF(5,5,5),XADD(9,9,9)
      write(6,*)'in x3arr'
c     call flush(6)
C
C     T3  SYMMETRY  OFFSETS
C
      NOFF=0
      DO 590 TSYM=1,NIRRED
      DO 585 USYM=1,NIRRED
      DO 580 VSYM=1,NIRRED
      DO 578 WSYM=1,NIRRED
      UVSYM=IEOR(USYM-1,VSYM-1)
      UVWSYM=IEOR(UVSYM,WSYM-1)+1
      IF(UVWSYM.NE.TSYM) GO TO 575
      DO 570 U=FLOV(USYM,1),FLOV(USYM,2)
      DO 565 V=FLOV(VSYM,1),FLOV(VSYM,2)
      DO 563 W=FLOV(WSYM,1),FLOV(WSYM,2)
      NADD=0
      XOFF(U,V,W)=NOFF
      DO 560 BESYM=1,NIRRED
      DO 555 GASYM=1,NIRRED
      DO 553 EPSYM=1,NIRRED
      BEGAS=IEOR(BESYM-1,GASYM-1)
      BEGAEP=IEOR(BEGAS,EPSYM-1)+1
      IF(BEGAEP.NE.TSYM) GO TO 550
      DO 545 BE=FLOV(BESYM,3)-NO,FLOV(BESYM,4)-NO
      DO 540 GA=FLOV(GASYM,3)-NO,FLOV(GASYM,4)-NO
      DO 538 EP=FLOV(EPSYM,3)-NO,FLOV(EPSYM,4)-NO
      NOFF=NOFF+1
      NADD=NADD+1
      XADD(BE,GA,EP)=NADD
      TOT=XOFF(U,V,W)+XADD(BE,GA,EP)
c     WRITE(6,658) U,V,W,XOFF(U,V,W),BE,GA,EP,XADD(BE,GA,EP),TOT
c     call flush(6)
  538 CONTINUE
  540 CONTINUE
  545 CONTINUE
  550 CONTINUE
  553 CONTINUE
  555 CONTINUE
  560 CONTINUE
  563 CONTINUE
  565 CONTINUE
  570 CONTINUE
  575 CONTINUE
  578 CONTINUE
  580 CONTINUE
  585 CONTINUE
  590 CONTINUE
C
      WRITE (*,*)'TOTAL IN SYMARR',TOT
C     WRITE(*,*)'IN SYMARR XADD',XADD
  657 FORMAT(' U=',I3,' V=',I3,' UOFF=',I6,' BE=',I3,' GA=',I3,
     .                        ' VADD=',I6,' TOT= ',I6)
  658 FORMAT(3I4,' UOFF=',I6,3I4,' VADD=',I6,' TOT= ',I6)
C
      END
