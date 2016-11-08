C
C----------------------------------------------------------------------
C
      SUBROUTINE SYMARR(ORBSYM,FLOV,NIRRED,NO,NV,NT,NOV,NM,
     .                  FZO,FZV,NFZO,NFZV,
     .                  UOFF,VADD,ZLX,ITR,ITV,NONO,NVNV,
     .                  AOFF,AADD,BOFF,BADD,COFF,CADD,
     .                  DOFF,DADD,EOFF,EADD,FOFF,FADD,TOFF,TADD,OPTION,
     .                  i1off,i2off,i3off)
      IMPLICIT INTEGER (A-Z)
      COMMON/INTS/NSGOO,NSGVV,NSGOV,NSHOV,NSLOV,NSLVO
      INTEGER FLOV(NIRRED,4),ORBSYM(NT),UOFF(NO,NO,2),VADD(NV,NV),
     .        ZLX(NV,NV),FZO(NO),FZV(NV),NFZO(NIRRED),NFZV(NIRRED),
     .        AOFF(*),AADD(*),BOFF(*),BADD(*),COFF(*),CADD(*),
     .        DOFF(*),DADD(*),EOFF(*),EADD(*),FOFF(*),FADD(*),
     .        TOFF(NO,NO,2,NIRRED),TADD(NV,NV,NIRRED),
     .        ITR(NOV),ITV(NV),
     .        i1off(nm),i2off(nm,nirred),i3off(nm)
      CHARACTER*4 OPTION
C
C     SET UP ARRAYS TO ADRESS T COEFF. AND 2E MO INTS.
C
      NOFF=0
      DO 124 TSYM=1,NIRRED
         DO 123 USYM=1,NIRRED
            DO 122 VSYM=1,NIRRED
               UVSYM=IEOR(USYM-1,VSYM-1)+1
               IF(UVSYM.NE.TSYM) GO TO 122
               DO 121 U=FLOV(USYM,1),FLOV(USYM,2)
                  DO 120 V=FLOV(VSYM,1),FLOV(VSYM,2)
                     NADD=0
                     UOFF(U,V,1)=NOFF
                     DO 119 BE=1,NV
                        BESYM=ORBSYM(BE+NO)
                        IF(FZV(BE).EQ.1) GO TO 119
                        DO 118 GA=1,BE
                           IF(FZV(GA).EQ.1) GO TO 118
                           GASYM=ORBSYM(GA+NO)
                           BEGSYM=IEOR(BESYM,GASYM)+1
                           IF(BEGSYM.NE.TSYM)GO TO 118
                           NOFF=NOFF+1
                           NADD=NADD+1
                           VADD(BE,GA)=NADD
C                          TOT=UOFF(U,V,1)+VADD(BE,GA)
C                  WRITE(6,657) U,V,UOFF(U,V,1),BE,GA,VADD(BE,GA),TOT
  118                   CONTINUE
  119                CONTINUE
  120             CONTINUE
  121          CONTINUE
  122       CONTINUE
  123    CONTINUE
  124 CONTINUE
C
      DO 128 U=1,NO
          DO 126 V=1,NO
             UOFF(V,U,2)=UOFF(U,V,1)
  126    CONTINUE
  128 CONTINUE
C
      DO 129 A=1,NV
         DO 127 B=1,A
            ZLX(B,A)=2
            ZLX(A,B)=1
            VADD(B,A)=VADD(A,B)
  127    CONTINUE
  129 CONTINUE
C
  657 FORMAT(' U=',I3,' V=',I3,' UOFF=',I6,' BE=',I3,' GA=',I3,
     .                        ' VADD=',I6,' TOT= ',I6)
C
C     A INTS  (OO,OO)
C
      NOFF=0
      DO 280 TSYM=1,NIRRED
      DO 270 ISYM=1,NIRRED
      DO 260 JSYM=1,ISYM
      IJSYM=IEOR(ISYM-1,JSYM-1)+1
      IF(IJSYM.NE.TSYM) GO TO 260
      FI=FLOV(ISYM,1)-NFZO(ISYM)
      DO 250 I=FI,FLOV(ISYM,2)
      JLIM=FLOV(JSYM,2)
      IF(ISYM.EQ.JSYM)JLIM=I
      FJ=FLOV(JSYM,1)-NFZO(JSYM)
      DO 240 J=FJ,JLIM
      NADD=0
      IJ=ITR(I)+J
      IF(FZO(I).EQ.1.OR.FZO(J).EQ.1)THEN
      AOFF(IJ)=NSGOO+1
      GO TO 202
      END IF
      AOFF(IJ)=NOFF
  202 CONTINUE
      DO 230 K=1,I
      KSYM=ORBSYM(K)
      LL=K
      IF(K.EQ.I)LL=J
      DO 220 L=1,LL
      LSYM=ORBSYM(L)
      KL=ITR(K)+L
      KLSYM=IEOR(KSYM,LSYM)+1
      IF(KLSYM.NE.TSYM)GO TO 210
      IF(FZO(I).EQ.1.OR.FZO(J).EQ.1)GO TO 205
      IF(FZO(K).EQ.1.OR.FZO(L).EQ.1)THEN
      AADD(KL)=NSGOO+1
      GO TO 205
      END IF
      NOFF=NOFF+1
      NADD=NADD+1
      AADD(KL)=NADD
  205 CONTINUE
C     TOT=AOFF(IJ)+AADD(KL)
C     WRITE(6,657) I,J,AOFF(IJ),K,L,AADD(KL),TOT
  210 CONTINUE
  220 CONTINUE
  230 CONTINUE
  240 CONTINUE
  250 CONTINUE
  260 CONTINUE
  270 CONTINUE
  280 CONTINUE
C
C     B INTS   (VV,VV)
C
      NOFF=0
      DO 380 TSYM=1,NIRRED
      DO 370 ISYM=1,NIRRED
      DO 360 JSYM=1,ISYM
      IJSYM=IEOR(ISYM-1,JSYM-1)+1
      IF(IJSYM.NE.TSYM) GO TO 360
      FI=FLOV(ISYM,3)-NO
      LI=FLOV(ISYM,4)+NFZV(ISYM)-NO
      DO 350 I=FI,LI
      JLIM=FLOV(JSYM,4) +NFZV(JSYM)-NO
      IF(ISYM.EQ.JSYM)JLIM=I
      DO 340 J=FLOV(JSYM,3)-NO,JLIM
      NADD=0
      IJ=ITV(I)+J
      IF(FZV(I).EQ.1.OR.FZV(J).EQ.1)THEN
      BOFF(IJ)=NSGVV+1
      GO TO 302
      END IF
      BOFF(IJ)=NOFF
  302 CONTINUE
      DO 330 K=1,I
      KSYM=ORBSYM(K+NO)
      LL=K
      IF(K.EQ.I)LL=J
      DO 320 L=1,LL
      LSYM=ORBSYM(L+NO)
      KL=ITV(K)+L
      KLSYM=IEOR(KSYM,LSYM)+1
      IF(KLSYM.NE.TSYM)GO TO 310
      IF(FZV(I).EQ.1.OR.FZV(J).EQ.1)GO TO 305
      IF(FZV(K).EQ.1.OR.FZV(L).EQ.1)THEN
      BADD(KL)=NSGVV+1
      GO TO 305
      END IF
      NOFF=NOFF+1
      NADD=NADD+1
      BADD(KL)=NADD
  305 CONTINUE
C     TOT=BOFF(IJ)+BADD(KL)
C     WRITE(6,657) I,J,BOFF(IJ),K,L,BADD(KL),TOT
  310 CONTINUE
  320 CONTINUE
  330 CONTINUE
  340 CONTINUE
  350 CONTINUE
  360 CONTINUE
  370 CONTINUE
  380 CONTINUE
C
C     C INTS   (VV,OO)
C
      NOFF=0
      DO 480 TSYM=1,NIRRED
      DO 470 ISYM=1,NIRRED
      DO 460 JSYM=1,ISYM
      IJSYM=IEOR(ISYM-1,JSYM-1)+1
      IF(IJSYM.NE.TSYM) GO TO 460
      FI=FLOV(ISYM,3)-NO
      LI=FLOV(ISYM,4)+NFZV(ISYM)-NO
      DO 450 I=FI,LI
      FJ=FLOV(JSYM,3)-NO
      LJ=FLOV(JSYM,4)+NFZV(JSYM)-NO
      IF(ISYM.EQ.JSYM)LJ=I
      DO 440 J=FJ,LJ
      NADD=0
      IJ=ITV(I)+J
      IF(FZV(I).EQ.1.OR.FZV(J).EQ.1)THEN
      COFF(IJ)=NSGOV+1
      GO TO 402
      END IF
      COFF(IJ)=NOFF
  402 CONTINUE
      DO 430 K=1,NO
      KSYM=ORBSYM(K)
      DO 420 L=1,K
      LSYM=ORBSYM(L)
      KL=ITR(K)+L
      KLSYM=IEOR(KSYM,LSYM)+1
      IF(KLSYM.NE.TSYM)GO TO 410
      IF(FZV(I).EQ.1.OR.FZV(J).EQ.1)GO TO 405
      IF(FZO(K).EQ.1.OR.FZO(L).EQ.1)THEN
      CADD(KL)=NSGOV+1
      GO TO 405
      END IF
      NOFF=NOFF+1
      NADD=NADD+1
      CADD(KL)=NADD
  405 CONTINUE
C     TOT=COFF(IJ)+CADD(KL)
C     WRITE(6,657) I,J,COFF(IJ),K,L,CADD(KL),TOT
  410 CONTINUE
  420 CONTINUE
  430 CONTINUE
  440 CONTINUE
  450 CONTINUE
  460 CONTINUE
  470 CONTINUE
  480 CONTINUE
C
C     D INTS   (VO,VO)
C
      NOFF=0
      DO 580 TSYM=1,NIRRED
      DO 570 ISYM=1,NIRRED
      DO 560 JSYM=1,NIRRED
      IJSYM=IEOR(ISYM-1,JSYM-1)+1
      IF(IJSYM.NE.TSYM) GO TO 560
      FI=FLOV(ISYM,3)-NO
      LI=FLOV(ISYM,4)+NFZV(ISYM)-NO
      DO 550 I=FI,LI
      FJ=FLOV(JSYM,1)-NFZO(JSYM)
      LJ=FLOV(JSYM,2)
      DO 540 J=FJ,LJ
      NADD=0
      IJ=ITR(I)+J
      IF(FZV(I).EQ.1.OR.FZO(J).EQ.1)THEN
      DOFF(IJ)=NSHOV+1
      GO TO 502
      END IF
      DOFF(IJ)=NOFF
  502 CONTINUE
      DO 530 K=1,I
      KSYM=ORBSYM(K+NO)
      LL=NO
      IF(K.EQ.I)LL=J
      DO 520 L=1,LL
      LSYM=ORBSYM(L)
      KL=ITR(K)+L
      KLSYM=IEOR(KSYM,LSYM)+1
      IF(KLSYM.NE.TSYM)GO TO 510
      IF(FZV(I).EQ.1.OR.FZO(J).EQ.1)GO TO 505
      IF(FZV(K).EQ.1.OR.FZO(L).EQ.1)THEN
      DADD(KL)=NSHOV+1
      GO TO 505
      END IF
      NOFF=NOFF+1
      NADD=NADD+1
      DADD(KL)=NADD
  505 CONTINUE
C     TOT=DOFF(IJ)+DADD(KL)
C     WRITE(6,657) I,J,DOFF(IJ),K,L,DADD(KL),TOT
  510 CONTINUE
  520 CONTINUE
  530 CONTINUE
  540 CONTINUE
  550 CONTINUE
  560 CONTINUE
  570 CONTINUE
  580 CONTINUE
C
C     E INTS   (VO,OO)
C
      NOFF=0
      DO 680 TSYM=1,NIRRED
      DO 670 ISYM=1,NIRRED
      DO 660 JSYM=1,NIRRED
      IJSYM=IEOR(ISYM-1,JSYM-1)+1
      IF(IJSYM.NE.TSYM) GO TO 660
      FI=FLOV(ISYM,3)-NO
      LI=FLOV(ISYM,4)+NFZV(ISYM)-NO
      DO 650 I=FI,LI
      FJ=FLOV(JSYM,1)-NFZO(JSYM)
      LJ=FLOV(JSYM,2)
      DO 640 J=FJ,LJ
      NADD=0
      IJ=ITR(I)+J
      IF(FZV(I).EQ.1.OR.FZO(J).EQ.1)THEN
      EOFF(IJ)=NSLOV+1
      GO TO 602
      END IF
      EOFF(IJ)=NOFF
  602 CONTINUE
      DO 630 K=1,NO
      KSYM=ORBSYM(K)
      DO 620 L=1,K
      LSYM=ORBSYM(L)
      KL=ITR(K)+L
      KLSYM=IEOR(KSYM,LSYM)+1
      IF(KLSYM.NE.TSYM)GO TO 610
      IF(FZV(I).EQ.1.OR.FZO(J).EQ.1)GO TO 605
      IF(FZO(K).EQ.1.OR.FZO(L).EQ.1)THEN
      EADD(KL)=NSLOV+1
      GO TO 605
      END IF
      NOFF=NOFF+1
      NADD=NADD+1
      EADD(KL)=NADD
  605 CONTINUE
C     TOT=EOFF(IJ)+EADD(KL)
C     WRITE(6,657) I,J,EOFF(IJ),K,L,EADD(KL),TOT
  610 CONTINUE
  620 CONTINUE
  630 CONTINUE
  640 CONTINUE
  650 CONTINUE
  660 CONTINUE
  670 CONTINUE
  680 CONTINUE
C
C     F INTS   (VV,VO)
C
      NOFF=0
      DO 780 TSYM=1,NIRRED
      DO 770 ISYM=1,NIRRED
      DO 760 JSYM=1,ISYM
      IJSYM=IEOR(ISYM-1,JSYM-1)+1
      IF(IJSYM.NE.TSYM) GO TO 760
      FI=FLOV(ISYM,3)-NO
      LI=FLOV(ISYM,4)+NFZV(ISYM)-NO
      DO 750 I=FI,LI
      FJ=FLOV(JSYM,3)-NO
      LJ=FLOV(JSYM,4)+NFZV(JSYM)-NO
      IF(ISYM.EQ.JSYM)LJ=I
      DO 740 J=FJ,LJ
      NADD=0
      IJ=ITV(I)+J
      IF(FZV(I).EQ.1.OR.FZV(J).EQ.1)THEN
      FOFF(IJ)=NSLVO+1
      GO TO 702
      END IF
      FOFF(IJ)=NOFF
  702 CONTINUE
      DO 730 K=1,NV
      KSYM=ORBSYM(K+NO)
      DO 720 L=1,NO
      LSYM=ORBSYM(L)
      KL=ITR(K)+L
      KLSYM=IEOR(KSYM,LSYM)+1
      IF(KLSYM.NE.TSYM)GO TO 710
      IF(FZV(I).EQ.1.OR.FZV(J).EQ.1)GO TO 705
      IF(FZV(K).EQ.1.OR.FZO(L).EQ.1)THEN
      FADD(KL)=NSLVO+1
      GO TO 705
      END IF
      NOFF=NOFF+1
      NADD=NADD+1
      FADD(KL)=NADD
  705 CONTINUE
C     TOT=FOFF(IJ)+FADD(KL)
C     WRITE(6,657) I,J,FOFF(IJ),K,L,FADD(KL),TOT
  710 CONTINUE
  720 CONTINUE
  730 CONTINUE
  740 CONTINUE
  750 CONTINUE
  760 CONTINUE
  770 CONTINUE
  780 CONTINUE
C
C     FORM ARRAYS FOR T3
C
      IF(OPTION.EQ.'SDT1'.or.option.eq.'FSDT')THEN
      DO 890 T3SYM=1,NIRRED
      NOFF=0
      DO 880 TSYM=1,NIRRED
      DO 875 USYM=1,NIRRED
      DO 870 VSYM=1,NIRRED
      UVSYM=IEOR(USYM-1,VSYM-1)
      IF(UVSYM.NE.TSYM-1)GO TO 870
      DO 865 U=FLOV(USYM,1),FLOV(USYM,2)
      DO 860 V=FLOV(VSYM,1),FLOV(VSYM,2)
      NADD=0
      TOFF(U,V,1,T3SYM)=NOFF
      DO 850 BE=1,NV
      IF(FZV(BE).EQ.1)GO TO 850
      BESYM=ORBSYM(BE+NO)
      DO 840 GA=1,BE
      IF(FZV(GA).EQ.1)GO TO 840
      GASYM=ORBSYM(GA+NO)
      BEGSYM=IEOR(BESYM,GASYM)
      XSYM=IEOR(UVSYM,BEGSYM)+1
      IF(XSYM.NE.T3SYM)GO TO 820
      NOFF=NOFF+1
      NADD=NADD+1
      TADD(BE,GA,T3SYM)=NADD
C     TOT=TOFF(U,V,1,T3SYM)+TADD(BE,GA,T3SYM)
C     WRITE(6,657) U,V,TOFF(U,V,1,T3SYM),BE,GA,TADD(BE,GA,T3SYM),TOT
  820 CONTINUE
  840 CONTINUE
  850 CONTINUE
  860 CONTINUE
  865 CONTINUE
  870 CONTINUE
  875 CONTINUE
  880 CONTINUE
      DO 885 U=1,NO
      DO 885 V=1,NO
  885 TOFF(V,U,2,T3SYM)=TOFF(U,V,1,T3SYM)
      DO 887 BE=1,NV
      DO 887 GA=1,BE
  887 TADD(GA,BE,T3SYM)=TADD(BE,GA,T3SYM)
  890 CONTINUE
C
      n1=0
      n2=0
      n3=0
      DO 1131 a=1,Nv
      IF(FZv(a).EQ.1)GO TO 1131
      aSYM=ORBSYM(a+no)
      DO 1132 i=1,no
      IF(FZO(i).EQ.1)GO TO 1132
      iSYM=ORBSYM(i)
      aiS=IEOR(aSYM,iSYM)
      ai=itr(a)+i
      i1off(ai)=n1
      n2=0
      DO 1133 b=1,Nv
      IF(FZv(b).EQ.1)GO TO 1133
      bSYM=ORBSYM(b+no)
      DO 1134 j=1,no
      IF(FZO(j).EQ.1)GO TO 1134
      jSYM=ORBSYM(j)
      bjS=IEOR(bSYM,jSYM)
      bj=itr(b)+j
      aibjs=ieor(ais,bjs)
      if(bj.gt.ai)go to 1134
      i2off(bj,aibjs+1)=n2
      n3=0
      DO 1135 c=1,Nv
      IF(FZv(c).EQ.1)GO TO 1135
      cSYM=ORBSYM(c+no)
      DO 1136 k=1,no
      IF(FZO(k).EQ.1)GO TO 1136
      kSYM=ORBSYM(k)
      ckS=IEOR(cSYM,kSYM)
      ck=itr(c)+k
      if(ck.gt.bj)go to 1136
      tsym=ieor(aibjs,cks)
      if(tsym.ne.0)go to 1136
      n1=n1+1
      n2=n2+1
      n3=n3+1
      i3off(ck)=n3
      tot=i1off(ai)+i2off(bj,aibjs+1)+i3off(ck)
c     write(6,768)a+no,i,b+no,j,c+no,k,ai,bj,ck,
c    .            i1off(ai),i2off(bj,aibjs+1),i3off(ck),tot   
  768 format(6i3,5x,3i3,9x,3i4,5x,i6)
 1136 CONTINUE
 1135 CONTINUE
 1134 CONTINUE
 1133 CONTINUE
 1132 CONTINUE
 1131 CONTINUE
      WRITE(6,*)' '
      WRITE(6,*)' in symarr , new offsets for core t3'
      WRITE(6,*)'THE NUMBER OF TRIPLE EXCITATIONS IS     ',tot       
      WRITE(6,*)' '
c     call flush(6)
C
C
      ENDIF
C
      RETURN
      END
