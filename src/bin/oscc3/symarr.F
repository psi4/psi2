      SUBROUTINE SYMARR(ORBSYM,FLOV,NIRRED,NO,NV,NT,NOV,NM,
     .                  FZO,FZV,NFZO,NFZV,
     .                  UOFF,VADD,ZLX,ITR,ITV,NONO,NVNV,
     .                  AOFF,AADD,BOFF,BADD,COFF,CADD,
     .                  DOFF,DADD,EOFF,EADD,FOFF,FADD,TOFF,TADD,OPTION,
     .                  ntmo,spin,ssov,mofz,mosym,
     .                  ioff,iadd,ipq,map,nocs,nvcs,f1,
     .                  nalp,moint)
      IMPLICIT INTEGER (A-Z)
      COMMON/INTS/NSGOO,NSGVV,NSGOV,NSHOV,NSLOV,NSLVO
      INTEGER FLOV(8,2,4),ORBSYM(NT),UOFF(no,no),VADD(nv,nv),
     .        ZLX(NV,NV),FZO(NO),FZV(NV),NFZO(NIRRED),NFZV(NIRRED),
     .        AOFF(*),AADD(*),BOFF(*),BADD(*),COFF(*),CADD(*),
     .        DOFF(*),DADD(*),EOFF(*),EADD(*),FOFF(*),FADD(*),
     .        TOFF(NO,NO,2,NIRRED),TADD(NV,NV,NIRRED),
     .        ITR(NOV),ITV(NV),ssov(2,4),spin(*),
     .        mofz(ntmo),mosym(ntmo),
     .        ioff(*),iadd(*),ipq(*),map(*),f1(nt)
      CHARACTER*4 OPTION
c
c     write(6,*)
c     write(6,*)'in symarr all ints are addressed '
c     write(6,*)'need them to construct Fock matrix'
c
c     write(6,*)'mofz'
c     write(6,*) mofz 
c
c     write(6,*)'mosym'
c     write(6,*) mosym 
C
C     SET UP ARRAYS TO ADRESS T COEFF. AND 2E MO INTS.
C
c     NOFF=0
c     DO 124 TSYM=1,NIRRED
c        DO 123 USYM=1,NIRRED
c           DO 122 VSYM=1,NIRRED
c              UVSYM=IEOR(USYM-1,VSYM-1)+1
c              IF(UVSYM.NE.TSYM) GO TO 122
c              DO 121 U=FLOV(USYM,1),FLOV(USYM,2)
c                 DO 120 V=FLOV(VSYM,1),FLOV(VSYM,2)
c                    NADD=0
c                    UOFF(U,V,1)=NOFF
c                    DO 119 BE=1,NV
c                       BESYM=ORBSYM(BE+NO)
c                       IF(FZV(BE).EQ.1) GO TO 119
c                       DO 118 GA=1,BE
c                          IF(FZV(GA).EQ.1) GO TO 118
c                          GASYM=ORBSYM(GA+NO)
c                          BEGSYM=IEOR(BESYM,GASYM)+1
c                          IF(BEGSYM.NE.TSYM)GO TO 118
c                          NOFF=NOFF+1
c                          NADD=NADD+1
c                          VADD(BE,GA)=NADD
C                          TOT=UOFF(U,V,1)+VADD(BE,GA)
C                  WRITE(6,657) U,V,UOFF(U,V,1),BE,GA,VADD(BE,GA),TOT
c 118                   CONTINUE
c 119                CONTINUE
c 120             CONTINUE
c 121          CONTINUE
c 122       CONTINUE
c 123    CONTINUE
c 124 CONTINUE
c
c     DO 128 U=1,NO
c         DO 126 V=1,NO
c            UOFF(V,U,2)=UOFF(U,V,1)
c 126    CONTINUE
c 128 CONTINUE
c
c     DO 129 A=1,NV
c        DO 127 B=1,A
c           ZLX(B,A)=2
c           ZLX(A,B)=1
c           VADD(B,A)=VADD(A,B)
c 127    CONTINUE
c 129 CONTINUE
c
  657 FORMAT(' U=',I3,' V=',I3,' UOFF=',I6,' BE=',I3,' GA=',I3,
     .                        ' VADD=',I6,' TOT= ',I6)
c
C     A INTS  (OO,OO)
C
c     NOFF=0
c     DO 280 TSYM=1,NIRRED
c     DO 270 ISYM=1,NIRRED
c     DO 260 JSYM=1,ISYM
c     IJSYM=IEOR(ISYM-1,JSYM-1)+1
c     IF(IJSYM.NE.TSYM) GO TO 260
c     FI=FLOV(ISYM,1)-NFZO(ISYM)
c     DO 250 I=FI,FLOV(ISYM,2)
c     JLIM=FLOV(JSYM,2)
c     IF(ISYM.EQ.JSYM)JLIM=I
c     FJ=FLOV(JSYM,1)-NFZO(JSYM)
c     DO 240 J=FJ,JLIM
c     NADD=0
c     IJ=ITR(I)+J
c     IF(FZO(I).EQ.1.OR.FZO(J).EQ.1)THEN
c     AOFF(IJ)=NSGOO+1
c     GO TO 202
c     END IF
c     AOFF(IJ)=NOFF
c 202 CONTINUE
c     DO 230 K=1,I
c     KSYM=ORBSYM(K)
c     LL=K
c     IF(K.EQ.I)LL=J
c     DO 220 L=1,LL
c     LSYM=ORBSYM(L)
c     KL=ITR(K)+L
c     KLSYM=IEOR(KSYM,LSYM)+1
c     IF(KLSYM.NE.TSYM)GO TO 210
c     IF(FZO(I).EQ.1.OR.FZO(J).EQ.1)GO TO 205
c     IF(FZO(K).EQ.1.OR.FZO(L).EQ.1)THEN
c     AADD(KL)=NSGOO+1
c     GO TO 205
c     END IF
c     NOFF=NOFF+1
c     NADD=NADD+1
c     AADD(KL)=NADD
c 205 CONTINUE
C     TOT=AOFF(IJ)+AADD(KL)
C     WRITE(6,657) I,J,AOFF(IJ),K,L,AADD(KL),TOT
c 210 CONTINUE
c 220 CONTINUE
c 230 CONTINUE
c 240 CONTINUE
c 250 CONTINUE
c 260 CONTINUE
c 270 CONTINUE
c 280 CONTINUE
C
C     B INTS   (VV,VV)
C
c     NOFF=0
c     DO 380 TSYM=1,NIRRED
c     DO 370 ISYM=1,NIRRED
c     DO 360 JSYM=1,ISYM
c     IJSYM=IEOR(ISYM-1,JSYM-1)+1
c     IF(IJSYM.NE.TSYM) GO TO 360
c     FI=FLOV(ISYM,3)-NO
c     LI=FLOV(ISYM,4)+NFZV(ISYM)-NO
c     DO 350 I=FI,LI
c     JLIM=FLOV(JSYM,4) +NFZV(JSYM)-NO
c     IF(ISYM.EQ.JSYM)JLIM=I
c     DO 340 J=FLOV(JSYM,3)-NO,JLIM
c     NADD=0
c     IJ=ITV(I)+J
c     IF(FZV(I).EQ.1.OR.FZV(J).EQ.1)THEN
c     BOFF(IJ)=NSGVV+1
c     GO TO 302
c     END IF
c     BOFF(IJ)=NOFF
c 302 CONTINUE
c     DO 330 K=1,I
c     KSYM=ORBSYM(K+NO)
c     LL=K
c     IF(K.EQ.I)LL=J
c     DO 320 L=1,LL
c     LSYM=ORBSYM(L+NO)
c     KL=ITV(K)+L
c     KLSYM=IEOR(KSYM,LSYM)+1
c     IF(KLSYM.NE.TSYM)GO TO 310
c     IF(FZV(I).EQ.1.OR.FZV(J).EQ.1)GO TO 305
c     IF(FZV(K).EQ.1.OR.FZV(L).EQ.1)THEN
c     BADD(KL)=NSGVV+1
c     GO TO 305
c     END IF
c     NOFF=NOFF+1
c     NADD=NADD+1
c     BADD(KL)=NADD
c 305 CONTINUE
c     TOT=BOFF(IJ)+BADD(KL)
c     WRITE(6,657) I,J,BOFF(IJ),K,L,BADD(KL),TOT
c 310 CONTINUE
c 320 CONTINUE
c 330 CONTINUE
c 340 CONTINUE
c 350 CONTINUE
c 360 CONTINUE
c 370 CONTINUE
c 380 CONTINUE
C
C     C INTS   (VV,OO)
c
c     NOFF=0
c     DO 480 TSYM=1,NIRRED
c     DO 470 ISYM=1,NIRRED
c     DO 460 JSYM=1,ISYM
c     IJSYM=IEOR(ISYM-1,JSYM-1)+1
c     IF(IJSYM.NE.TSYM) GO TO 460
c     FI=FLOV(ISYM,3)-NO
c     LI=FLOV(ISYM,4)+NFZV(ISYM)-NO
c     DO 450 I=FI,LI
c     FJ=FLOV(JSYM,3)-NO
c     LJ=FLOV(JSYM,4)+NFZV(JSYM)-NO
c     IF(ISYM.EQ.JSYM)LJ=I
c     DO 440 J=FJ,LJ
c     NADD=0
c     IJ=ITV(I)+J
c     IF(FZV(I).EQ.1.OR.FZV(J).EQ.1)THEN
c     COFF(IJ)=NSGOV+1
c     GO TO 402
c     END IF
c     COFF(IJ)=NOFF
c 402 CONTINUE
c     DO 430 K=1,NO
c     KSYM=ORBSYM(K)
c     DO 420 L=1,K
c     LSYM=ORBSYM(L)
c     KL=ITR(K)+L
c     KLSYM=IEOR(KSYM,LSYM)+1
c     IF(KLSYM.NE.TSYM)GO TO 410
c     IF(FZV(I).EQ.1.OR.FZV(J).EQ.1)GO TO 405
c     IF(FZO(K).EQ.1.OR.FZO(L).EQ.1)THEN
c     CADD(KL)=NSGOV+1
c     GO TO 405
c     END IF
c     NOFF=NOFF+1
c     NADD=NADD+1
c     CADD(KL)=NADD
c 405 CONTINUE
c     TOT=COFF(IJ)+CADD(KL)
c     WRITE(6,657) I,J,COFF(IJ),K,L,CADD(KL),TOT
c 410 CONTINUE
c 420 CONTINUE
c 430 CONTINUE
c 440 CONTINUE
c 450 CONTINUE
c 460 CONTINUE
c 470 CONTINUE
c 480 CONTINUE
c
c     D INTS   (VO,VO)
c
c     NOFF=0
c     DO 580 TSYM=1,NIRRED
c     DO 570 ISYM=1,NIRRED
c     DO 560 JSYM=1,NIRRED
c     IJSYM=IEOR(ISYM-1,JSYM-1)+1
c     IF(IJSYM.NE.TSYM) GO TO 560
c     FI=FLOV(ISYM,3)-NO
c     LI=FLOV(ISYM,4)+NFZV(ISYM)-NO
c     DO 550 I=FI,LI
c     FJ=FLOV(JSYM,1)-NFZO(JSYM)
c     LJ=FLOV(JSYM,2)
c     DO 540 J=FJ,LJ
c     NADD=0
c     IJ=ITR(I)+J
c     IF(FZV(I).EQ.1.OR.FZO(J).EQ.1)THEN
c     DOFF(IJ)=NSHOV+1
c     GO TO 502
c     END IF
c     DOFF(IJ)=NOFF
c 502 CONTINUE
c     DO 530 K=1,I
c     KSYM=ORBSYM(K+NO)
c     LL=NO
c     IF(K.EQ.I)LL=J
c     DO 520 L=1,LL
c     LSYM=ORBSYM(L)
c     KL=ITR(K)+L
c     KLSYM=IEOR(KSYM,LSYM)+1
c     IF(KLSYM.NE.TSYM)GO TO 510
c     IF(FZV(I).EQ.1.OR.FZO(J).EQ.1)GO TO 505
c     IF(FZV(K).EQ.1.OR.FZO(L).EQ.1)THEN
c     DADD(KL)=NSHOV+1
c     GO TO 505
c     END IF
c     NOFF=NOFF+1
c     NADD=NADD+1
c     DADD(KL)=NADD
c 505 CONTINUE
c     TOT=DOFF(IJ)+DADD(KL)
c     WRITE(6,657) I,J,DOFF(IJ),K,L,DADD(KL),TOT
c 510 CONTINUE
c 520 CONTINUE
c 530 CONTINUE
c 540 CONTINUE
c 550 CONTINUE
c 560 CONTINUE
c 570 CONTINUE
c 580 CONTINUE
c
c     E INTS   (VO,OO)
c
c     NOFF=0
c     DO 680 TSYM=1,NIRRED
c     DO 670 ISYM=1,NIRRED
c     DO 660 JSYM=1,NIRRED
c     IJSYM=IEOR(ISYM-1,JSYM-1)+1
c     IF(IJSYM.NE.TSYM) GO TO 660
c     FI=FLOV(ISYM,3)-NO
c     LI=FLOV(ISYM,4)+NFZV(ISYM)-NO
c     DO 650 I=FI,LI
c     FJ=FLOV(JSYM,1)-NFZO(JSYM)
c     LJ=FLOV(JSYM,2)
c     DO 640 J=FJ,LJ
c     NADD=0
c     IJ=ITR(I)+J
c     IF(FZV(I).EQ.1.OR.FZO(J).EQ.1)THEN
c     EOFF(IJ)=NSLOV+1
c     GO TO 602
c     END IF
c     EOFF(IJ)=NOFF
c 602 CONTINUE
c     DO 630 K=1,NO
c     KSYM=ORBSYM(K)
c     DO 620 L=1,K
c     LSYM=ORBSYM(L)
c     KL=ITR(K)+L
c     KLSYM=IEOR(KSYM,LSYM)+1
c     IF(KLSYM.NE.TSYM)GO TO 610
c     IF(FZV(I).EQ.1.OR.FZO(J).EQ.1)GO TO 605
c     IF(FZO(K).EQ.1.OR.FZO(L).EQ.1)THEN
c     EADD(KL)=NSLOV+1
c     GO TO 605
c     END IF
c     NOFF=NOFF+1
c     NADD=NADD+1
c     EADD(KL)=NADD
c 605 CONTINUE
c     TOT=EOFF(IJ)+EADD(KL)
c     WRITE(6,657) I,J,EOFF(IJ),K,L,EADD(KL),TOT
c 610 CONTINUE
c 620 CONTINUE
c 630 CONTINUE
c 640 CONTINUE
c 650 CONTINUE
c 660 CONTINUE
c 670 CONTINUE
c 680 CONTINUE
c
c     F INTS   (VV,VO)
c
c     NOFF=0
c     DO 780 TSYM=1,NIRRED
c     DO 770 ISYM=1,NIRRED
c     DO 760 JSYM=1,ISYM
c     IJSYM=IEOR(ISYM-1,JSYM-1)+1
c     IF(IJSYM.NE.TSYM) GO TO 760
c     FI=FLOV(ISYM,3)-NO
c     LI=FLOV(ISYM,4)+NFZV(ISYM)-NO
c     DO 750 I=FI,LI
c     FJ=FLOV(JSYM,3)-NO
c     LJ=FLOV(JSYM,4)+NFZV(JSYM)-NO
c     IF(ISYM.EQ.JSYM)LJ=I
c     DO 740 J=FJ,LJ
c     NADD=0
c     IJ=ITV(I)+J
c     IF(FZV(I).EQ.1.OR.FZV(J).EQ.1)THEN
c     FOFF(IJ)=NSLVO+1
c     GO TO 702
c     END IF
c     FOFF(IJ)=NOFF
c 702 CONTINUE
c     DO 730 K=1,NV
c     KSYM=ORBSYM(K+NO)
c     DO 720 L=1,NO
c     LSYM=ORBSYM(L)
c     KL=ITR(K)+L
c     KLSYM=IEOR(KSYM,LSYM)+1
c     IF(KLSYM.NE.TSYM)GO TO 710
c     IF(FZV(I).EQ.1.OR.FZV(J).EQ.1)GO TO 705
c     IF(FZV(K).EQ.1.OR.FZO(L).EQ.1)THEN
c     FADD(KL)=NSLVO+1
c     GO TO 705
c     END IF
c     NOFF=NOFF+1
c     NADD=NADD+1
c     FADD(KL)=NADD
c 705 CONTINUE
c     TOT=FOFF(IJ)+FADD(KL)
c     WRITE(6,657) I,J,FOFF(IJ),K,L,FADD(KL),TOT
c 710 CONTINUE
c 720 CONTINUE
c 730 CONTINUE
c 740 CONTINUE
c 750 CONTINUE
c 760 CONTINUE
c 770 CONTINUE
c 780 CONTINUE
c
C     FORM ARRAYS FOR T3
c
c     IF(OPTION.EQ.'SDT1')THEN
c     DO 890 T3SYM=1,NIRRED
c     NOFF=0
c     DO 880 TSYM=1,NIRRED
c     DO 875 USYM=1,NIRRED
c     DO 870 VSYM=1,NIRRED
c     UVSYM=IEOR(USYM-1,VSYM-1)
c     IF(UVSYM.NE.TSYM-1)GO TO 870
c     DO 865 U=FLOV(USYM,1),FLOV(USYM,2)
c     DO 860 V=FLOV(VSYM,1),FLOV(VSYM,2)
c     NADD=0
c     TOFF(U,V,1,T3SYM)=NOFF
c     DO 850 BE=1,NV
c     IF(FZV(BE).EQ.1)GO TO 850
c     BESYM=ORBSYM(BE+NO)
c     DO 840 GA=1,BE
c     IF(FZV(GA).EQ.1)GO TO 840
c     GASYM=ORBSYM(GA+NO)
c     BEGSYM=IEOR(BESYM,GASYM)
c     XSYM=IEOR(UVSYM,BEGSYM)+1
c     IF(XSYM.NE.T3SYM)GO TO 820
c     NOFF=NOFF+1
c     NADD=NADD+1
c     TADD(BE,GA,T3SYM)=NADD
c     TOT=TOFF(U,V,1,T3SYM)+TADD(BE,GA,T3SYM)
c     WRITE(6,657) U,V,TOFF(U,V,1,T3SYM),BE,GA,TADD(BE,GA,T3SYM),TOT
c 820 CONTINUE
c 840 CONTINUE
c 850 CONTINUE
c 860 CONTINUE
c 865 CONTINUE
c 870 CONTINUE
c 875 CONTINUE
c 880 CONTINUE
c     DO 885 U=1,NO
c     DO 885 V=1,NO
c 885 TOFF(V,U,2,T3SYM)=TOFF(U,V,1,T3SYM)
c     DO 887 BE=1,NV
c     DO 887 GA=1,BE
c 887 TADD(GA,BE,T3SYM)=TADD(BE,GA,T3SYM)
c 890 CONTINUE
c     ENDIF
c
c     NT1=0
c     DO 120 U=1,no
c        usym=orbsym(u)
c        uspn=spin(u)
c        fbe=flov(usym+1,uspn,3)-no
c        lbe=flov(usym+1,uspn,4)-no
c        DO 100 BE=fbe,lbe
c           NT1=NT1+1
c 100    CONTINUE
c 120 CONTINUE
C
      noff=0
      do 250 u=1,no
         if(fzo(u).eq.1)go to 250
         usym=orbsym(u)
         uspn=spin(u)
         xu=map(u)
         uoff(u,u)=1
         do 240 v=1,u-1
            if(fzo(v).eq.1)go to 240
            vsym=orbsym(v)
            uvsym=ieor(usym,vsym)
            vspn=spin(v)
         xv=map(v)
            fbe=ssov(uspn,3)-no
            lbe=ssov(uspn,4)-no
            nadd=0
            uoff(u,v)=noff
            uoff(v,u)=noff
            do 230 be=fbe,lbe
               besym=orbsym(be+no)
               xbe=map(be+no)
               gasym=ieor(uvsym,besym)+1
               fga=flov(gasym,vspn,3)-no
               lga=flov(gasym,vspn,4)-no
               if(lga.ge.be)lga=be-1
               do 220 ga=fga,lga
               xga=map(ga+no)
                  noff=noff+1
                  nadd=nadd+1
                  vadd(be,ga)=nadd
                  vadd(ga,be)=nadd
                  TOT=UOFF(U,V)+VADD(BE,GA)
c      WRITE(6,659) U,xu,uspn,V,xv,vspn,BE,xbe,uspn,GA,xga,vspn,
c    .              UOFF(U,V), VADD(BE,GA),TOT
  659 format (3i4,3x,3i4,3x,3i4,3x,3i4,3x,  3x,2i4,i6)
  220          continue
  230       continue
  240    continue
  250 continue
      do 252 a=1,nv
         vadd(a,a)=1
  252 continue
C
c
c     noff=0
c     DO 490 I=1,Ntmo
c        IF(moFZ(I).EQ.1)GO TO 490
c        ISYM= moSYM(I)
c        DO 480 J=1,I
c           IF(moFZ(J).EQ.1)GO TO 480
c           JSYM= moSYM(J)
c           IJSYM=IEOR(ISYM,JSYM)
c           nadd=0
c           ij=ipq(i)+j
c           ioff(ij)=noff
c           DO 470 K=1,I
c              IF(moFZ(K).EQ.1)GO TO 470
c              KSYM= moSYM(K)
c              LL=K
c              IF(K.EQ.I)LL=J
c              DO 460 L=1,LL
c                 IF(moFZ(L).EQ.1)GO TO 460
c                 LSYM= moSYM(L)
c                 KLSYM=IEOR(KSYM,LSYM)
c                 IF(IJSYM.NE.KLSYM) GO TO 460
c                 noff=noff+1
c                 nadd=nadd+1
c                 kl=ipq(k)+l
c                 iadd(kl)=nadd
c                 tot=ioff(IJ)+iadd(KL)
c                 write(6,658)i,j,ioff(ij),k,l,iadd(KL),tot
c 460          CONTINUE
c 470       CONTINUE
c 480    CONTINUE
c 490 CONTINUE
c
c        write(6,*)'i,mosym(i),f1(i),mosym(f1(i))'
c     do 512 i=1,ntmo
c        write(6,'(4i4)')i,mosym(i),f1(i),mosym(f1(i))
c 512 continue
c
c     write(6,*)'i,j,ioff(ij),k,l,iadd(KL),tot'
      noff=0
c     DO 590 I=1,ntmo     
      DO 590 I=nvcs+1,Ntmo
c        IF(moFZ(I).EQ.1)GO TO 590
         i1=f1(i) 
         ISYM= moSYM(I1)
         DO 580 J=1,I
c           IF(moFZ(J).EQ.1)GO TO 580
            j1=f1(j)
            JSYM= moSYM(J1)
            IJSYM=IEOR(ISYM,JSYM)
            nadd=0
            ij=ipq(i)+j
            ioff(ij)=noff
            DO 570 K=1,I
c              IF(moFZ(K).EQ.1)GO TO 570
               k1=f1(k)
               KSYM= moSYM(K1)
               LL=K
               IF(K.EQ.I)LL=J
               DO 560 L=1,LL
c                 IF(moFZ(L).EQ.1)GO TO 560
                  l1=f1(l)
                  LSYM= moSYM(L1)
                  KLSYM=IEOR(KSYM,LSYM)
                  IF(IJSYM.NE.KLSYM) GO TO 560
                  noff=noff+1
                  nadd=nadd+1
                  kl=ipq(k)+l
                  iadd(kl)=nadd
                  tot=ioff(IJ)+iadd(KL)
c                 write(6,658)i,j,ioff(ij),k,l,iadd(KL),tot
  560          CONTINUE
  570       CONTINUE
  580    CONTINUE
  590 CONTINUE
c
  658 format(3i6,5x,3i6,i12)
c
      RETURN
      END
