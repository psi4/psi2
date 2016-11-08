      SUBROUTINE NCOUNT(ORBSYM,FLOV,NIRRED,NO,NV,NT,FZO,FZV,NT1,NT2,
     .                  NSGOO,NSGVV,NSGOV,NSHOV,NSLOV,NSLVO,NT3,NTAU,
     .                  OPTION,MAXOO,MAXVV,DIMAR,XT3,MAXO,MAXV,
     .                  spin,ssov,mofz,mosym,ntmo,nsint,
     .                  maxova,maxovb,nvcs,f1)
      IMPLICIT INTEGER (A-Z)
      INTEGER FLOV(8,2,4),ORBSYM(NT),FZO(no),FZV(nv),NT3(8),
     .        mosym(*),mofz(*),spin(*),ssov(2,4),xspin(4),yspin(4),
     .        f1(ntmo)
      CHARACTER*4 OPTION
C
C     COUNT THE NON-ZERO SYMMETRY AND ACTIVE T COEFFICIENTS AND
C     2E MO INTS ONLY FOR ALLOCATION PURPOSES
C
      NT1=0
      DO 120 U=1,no
         if(fzo(u).eq.1)go to 120
         usym=orbsym(u)
         uspn=spin(u)
         fbe=flov(usym+1,uspn,3)-no
         lbe=flov(usym+1,uspn,4)-no
         DO 100 BE=fbe,lbe
            NT1=NT1+1
  100    CONTINUE
  120 CONTINUE
C
      nt2=0
      do 250 u=1,no
         if(fzo(u).eq.1)go to 250
         usym=orbsym(u)
         uspn=spin(u)
         do 240 v=1,u-1
            if(fzo(v).eq.1)go to 240
            vsym=orbsym(v)
            uvsym=ieor(usym,vsym)
            vspn=spin(v)
            fbe=ssov(uspn,3)-no
            lbe=ssov(uspn,4)-no
            do 230 be=fbe,lbe
               besym=orbsym(be+no)
               gasym=ieor(uvsym,besym)+1
               fga=flov(gasym,vspn,3)-no
               lga=flov(gasym,vspn,4)-no
               if(lga.ge.be)lga=be-1
               do 220 ga=fga,lga
                  nt2=nt2+1
  220          continue
  230       continue
  240    continue
  250 continue
C
C     DIMENSION OF AR1-6
C
      MAXOO=0
      MAXVV=0
C
      xspin(1)=1
      yspin(1)=1
      xspin(2)=1
      yspin(2)=2
      xspin(3)=2
      yspin(3)=1
      xspin(4)=2
      yspin(4)=2
c 
      do 2001 case=1,4
      spin1=xspin(case)
      spin2=yspin(case)
c     WRITE(6,*)' CASE, spins=',case,spin1,spin2
      DO 2000 TSYM=1,NIRRED
      XAB=0
      Xij=0
      fa=ssov(spin1,3)-no
      la=ssov(spin1,4)-no
      DO 19490 A=fa,la
      ASYM=ORBSYM(A+no)
      fb=ssov(spin2,3)-no
      lb=ssov(spin2,4)-no
      if(lb.ge.a)lb=a-1
      DO 19480 B=fb,lb
      BSYM=ORBSYM(B+no)
      BASYM=IEOR(BSYM,ASYM)
      IF(BASYM.NE.TSYM-1)GO TO 19480
      XIJ=0
      XAB=XAB+1
      fi=ssov(spin1,1)
      li=ssov(spin1,2)
      DO 19470 I=fi,li
      ISYM=ORBSYM(I)
      JSYM=IEOR(ISYM,BASYM)+1
      FJ=FLOV(JSYM,spin2,1)
      LJ=FLOV(JSYM,spin2,2)
      IF(LJ.Ge.I)LJ=I-1
      DO 19460 J=FJ,LJ
      XIJ=XIJ+1
19460 CONTINUE
19470 CONTINUE
19480 CONTINUE
19490 CONTINUE
      IF(XIJ.GT.MAXOO)MAXOO=XIJ
      IF(XAB.GT.MAXVV)MAXVV=XAB
c     WRITE(6,*)'TSYM=',TSYM,' XiJ=',XiJ,' XAb=',XAb,'   ',Xij*Xab
 2000 CONTINUE
 2001 CONTINUE
C
      MAXOV1=0
      MAXOV2=0
      do 2002 case=1,4
      spin1=xspin(case)
      spin2=yspin(case)
c     WRITE(6,*)' CASE, spins=',case,spin1,spin2
c
      DO 2003 TSYM=1,NIRRED
      fi=ssov(spin1,1)
      li=ssov(spin1,2)
      XAI=0
      DO 15790 I=fi,li
      ISYM=ORBSYM(I)
      fa=ssov(spin1,3)-no
      la=ssov(spin1,4)-no
      DO 15770 A=fa,la
      ASYM=ORBSYM(A+no)
      IASYM=IEOR(ISYM,ASYM)
      IF(IASYM.NE.TSYM-1)GO TO 15770
      XAI=XAI+1
      XBJ=0
      fj=ssov(spin2,1)
      lj=ssov(spin2,2)
      DO 15750 J=fj,lj
      JSYM=ORBSYM(J)
      BSYM=IEOR(JSYM,IASYM)+1
      FB=FLOV(BSYM,spin2,3)-no
      LB=FLOV(BSYM,spin2,4)-no
      DO 15730 B=FB,LB
      XBJ=XBJ+1
15730 CONTINUE
15750 CONTINUE
15770 CONTINUE
15790 CONTINUE
      IF(XAI.GT.MAXOV1)MAXOV1=XAI
      IF(Xbj.GT.MAXOV2)MAXOV2=Xbj
c     WRITE(6,*)'TSYM=',TSYM,' Xai=',Xai,' Xbj=',Xbj,'   ',XAI*XBJ
c
 2003 CONTINUE
 2002 CONTINUE
C
      xspin(1)=1
      yspin(1)=2
      xspin(2)=2
      yspin(2)=1
c
      MAXOV3=0
      MAXOV4=0
      do 2004 case=1,2
      spin1=xspin(case)
      spin2=yspin(case)
c     WRITE(6,*)' CASE, spins=',case,spin1,spin2
c
      DO 2005 TSYM=1,NIRRED
      fi=ssov(spin1,1)
      li=ssov(spin1,2)
      XAI=0
      DO 16790 I=fi,li
      ISYM=ORBSYM(I)
      fa=ssov(spin2,3)-no
      la=ssov(spin2,4)-no
      DO 16770 A=fa,la
      ASYM=ORBSYM(A+no)
      IASYM=IEOR(ISYM,ASYM)
      IF(IASYM.NE.TSYM-1)GO TO 16770
      XAI=XAI+1
      XBJ=0
      fj=ssov(spin2,1)
      lj=ssov(spin2,2)
      DO 16750 J=fj,lj
      JSYM=ORBSYM(J)
      BSYM=IEOR(JSYM,IASYM)+1
      FB=FLOV(BSYM,spin1,3)-no
      LB=FLOV(BSYM,spin1,4)-no
      DO 16730 B=FB,LB
      XBJ=XBJ+1
16730 CONTINUE
16750 CONTINUE
16770 CONTINUE
16790 CONTINUE
      IF(XAI.GT.MAXOV3)MAXOV3=XAI
      IF(Xbj.GT.MAXOV4)MAXOV4=Xbj
c     WRITE(6,*)'TSYM=',TSYM,' Xai=',Xai,' Xbj=',Xbj,'   ',XAI*XBJ
c
 2005 CONTINUE
 2004 CONTINUE
C
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)' >>>  ARRAYS  DIMENSION   <<<'
      WRITE(6,*)
      WRITE(6,*)' Maxov1=',MAXOV1,'  Maxov2=',maxov2 
      WRITE(6,*)' Maxov3=',MAXOV3,'  Maxov4=',maxov4 
      WRITE(6,*)' MAXOO =',MAXOO,'  MAXVV=',MAXVV
      maxova=max0(maxov1,maxov2)
      maxovb=max0(maxov3,maxov4)
      DIMAR=MAX0(MAXOO*MAXVV,MAXOVa*MAXOVa,maxovb*maxovb)
      WRITE(6,*)' DIMAR=', DIMAR
c
C
C     COUNT SYMMETRY NON-ZERO T3
C
c     NTAU=NT2
c     IF(OPTION.EQ.'SDT1')THEN
c     WRITE(6,*)
c     DO 440 KCSYM=1,NIRRED
c     NPFF=0
c     DO 430 U=1,NO
c     IF(FZO(U).EQ.1)GO TO 430
c     USYM=ORBSYM(U)
c     DO 420 V=1,NO
c     IF(FZO(V).EQ.1)GO TO 420
c     VSYM=ORBSYM(V)
c     UVSYM=IEOR(USYM,VSYM)
c     DO 418 BE=1,NV
c     IF(FZV(BE).EQ.1)GO TO 418
c     BESYM=ORBSYM(BE+NO)
c     DO 417 GA=1,BE
c     IF(FZV(GA).EQ.1)GO TO 417
c     GASYM=ORBSYM(GA+NO)
c     BEGSYM=IEOR(BESYM,GASYM)
c     TSYM=IEOR(UVSYM,BEGSYM)
c     IF((KCSYM-1).NE.TSYM)GO TO 415
c     NPFF=NPFF+1
CC    WRITE(*,657)U,V,BE,GA,NPFF
C 657 FORMAT(' U=',I3,' V=',I3,' BE=',I3,' GA=',I3,'NPFF= ',I6)
c 415 CONTINUE
c 417 CONTINUE
c 418 CONTINUE
c 420 CONTINUE
c 430 CONTINUE
c     NT3(KCSYM)=NPFF
c 440 CONTINUE
c     WRITE(6,7)NT3(1)
c   7 FORMAT(1X,'NT3(1) =',I6)
c     IF(NT3(1).NE.NT2)WRITE(*,*)'WARNING ...   NT3(1).NE.NT2'
c     NTAU=NT3(1)
c     DO 455 I=2,NIRRED
c     IF(NT3(I).GT.NT3(1))THEN
c     WRITE(*,*)'WARNING !!!  SOME NT3.GT.NT2'
c     NTAU=NT3(I)
c     ENDIF
c     WRITE(6,8)I,NT3(I)
c   8 FORMAT(1X,'NT3(',I1,') =',I6)
c 455 CONTINUE
C
c     DO 1008 I=1,NIRRED
c1008 MT1(I)=0
c     DO 1010 I=1,NO
c     ISYM=ORBSYM(I)
c     IF(FZO(I).EQ.1)GO TO 1010
c     DO 1009 A=1,NV
c     ASYM=ORBSYM(A+NO)
c     IF(FZV(A).EQ.1)GO TO 1009
c     IASYM=IEOR(ISYM,ASYM)+1
c     MT1(IASYM)=MT1(IASYM)+1
c1009 CONTINUE
c1010 CONTINUE
C
c     MT3=0
c     DO 1020 I=1,NIRRED
c     MT3=MT3+NT3(I)*MT1(I)
c1020 CONTINUE
c     WRITE(*,*)'THE NUMBER OF T3 AMPLITUDES NEEDED IN DISK IS ',MT3
C
c     XT3=0
c     DO 1031 I=1,NO
c     IF(FZO(I).EQ.1)GO TO 1031
c     ISYM=ORBSYM(I)
c     DO 1032 J=1,I
c     IF(FZO(J).EQ.1)GO TO 1032
c     JSYM=ORBSYM(J)
c     IJS=IEOR(ISYM,JSYM)
c     DO 1033 K=1,J
c     IF(FZO(K).EQ.1)GO TO 1033
c     KSYM=ORBSYM(K)
c     IJKS=IEOR(IJS,KSYM)
c     DO 1034 A=1,NV
c     IF(FZV(A).EQ.1)GO TO 1034
c     ASYM=ORBSYM(A+NO)
c     IJKAS=IEOR(IJKS,ASYM)
c     DO 1035 B=1,A
c     IF(FZV(B).EQ.1)GO TO 1035
c     BSYM=ORBSYM(B+NO)
c     IJKABS=IEOR(IJKAS,BSYM)
c     DO 1036 C=1,B
c     IF(FZV(C).EQ.1)GO TO 1036
c     CSYM=ORBSYM(C+NO)
c     TSYM=IEOR(IJKABS,CSYM)+1
c     IF(TSYM.EQ.1)XT3=XT3+1
c1036 CONTINUE
c1035 CONTINUE
c1034 CONTINUE
c1033 CONTINUE
c1032 CONTINUE
c1031 CONTINUE
C     WRITE(6,*)' '
C     WRITE(6,*)'THE NUMBER OF TRIPLE EXCITATIONS IS     ',XT3
C     WRITE(6,*)' '
C
c     WRITE(6,*)
c     WRITE(6,*)'VOV ARRAYS FOR CCSDT-1A'
c     MAXCNT=0
c     MAXO=0
c     MAXV=0
c     DO 2390 TSYM=1,NIRRED
c     NOSYM=FLOV(TSYM,2)-FLOV(TSYM,1)+1
c     NVSYM=FLOV(TSYM,4)-FLOV(TSYM,3)+1
c     IF(MAXO.LT.NOSYM)MAXO=NOSYM
c     IF(MAXV.LT.NVSYM)MAXV=NVSYM
c     ICNT=0
c     DO 2380 A=1,NV
c     IF(FZV(A).EQ.1)GO TO 2380
c     ASYM=ORBSYM(A+NO)
c     DO 2370 I=1,NO
c     IF(FZO(I).EQ.1)GO TO 2370
c     AISYM=IEOR(ASYM,ISYM)
c     BSYM=IEOR(TSYM-1,AISYM)+1
c     FB=FLOV(BSYM,3)-NO
c     LB=FLOV(BSYM,4)-NO
c     DO 2360 B=FB,LB
c     ICNT=ICNT+1
c2360 CONTINUE
c2370 CONTINUE
c2380 CONTINUE
c     IF(MAXCNT.LT.ICNT)MAXCNT=ICNT
c     WRITE(6,*)'TSYM=',TSYM,'  ICNT=',ICNT
c2390 CONTINUE
c     IF(MAXCNT.LT.MAXOV*MAXV)THEN
c     WRITE(6,*)'ERROR  IN INE NCOUNT, MAXOV*MAXV.GT.MAXOV*MAXOV'
c     STOP
c     ENDIF
c     WRITE(6,*)
c     WRITE(6,*)'MAXO=',MAXO,'  MAXV=',MAXV
c     WRITE(6,*)
c     ENDIF
C
      NSGOO=0
c     DO 390 I=1,NO
c     IF(FZO(I).EQ.1)GO TO 390
c     ISYM=ORBSYM(I)
c     DO 380 J=1,I
c     IF(FZO(J).EQ.1)GO TO 380
c     JSYM=ORBSYM(J)
c     IJSYM=IEOR(ISYM,JSYM)
c     DO 370 K=1,I
c     IF(FZO(K).EQ.1)GO TO 370
c     KSYM=ORBSYM(K)
c     LL=K
c     IF(K.EQ.I)LL=J
c     DO 360 L=1,LL
c     IF(FZO(L).EQ.1)GO TO 360
c     LSYM=ORBSYM(L)
c     KLSYM=IEOR(KSYM,LSYM)
c     IF(IJSYM.NE.KLSYM) GO TO 350
c     NSGOO=NSGOO+1
c 350 CONTINUE
c 360 CONTINUE
c 370 CONTINUE
c 380 CONTINUE
c 390 CONTINUE
c
      NSGVV=0
c     DO 490 I=1,NV
c     IF(FZV(I).EQ.1)GO TO 490
c     ISYM=ORBSYM(I+NO)
c     DO 480 J=1,I
c     IF(FZV(J).EQ.1)GO TO 480
c     JSYM=ORBSYM(J+NO)
c     IJSYM=IEOR(ISYM,JSYM)
c     DO 470 K=1,I
c     IF(FZV(K).EQ.1)GO TO 470
c     KSYM=ORBSYM(K+NO)
c     LL=K
c     IF(K.EQ.I)LL=J
c     DO 460 L=1,LL
c     IF(FZV(L).EQ.1)GO TO 460
c     LSYM=ORBSYM(L+NO)
c     KLSYM=IEOR(KSYM,LSYM)
c     IF(IJSYM.NE.KLSYM) GO TO 450
c     NSGVV=NSGVV+1
c 450 CONTINUE
c 460 CONTINUE
c 470 CONTINUE
c 480 CONTINUE
c 490 CONTINUE
c
      NSGOV=0
c     DO 590 I=1,NO
c     IF(FZO(I).EQ.1)GO TO 590
c     ISYM=ORBSYM(I)
c     DO 580 J=1,I
c     IF(FZO(J).EQ.1)GO TO 580
c     JSYM=ORBSYM(J)
c     IJSYM=IEOR(ISYM,JSYM)
c     DO 570 K=1,NV
c     IF(FZV(K).EQ.1)GO TO 570
c     KSYM=ORBSYM(K+NO)
c     DO 560 L=1,K
c     IF(FZV(L).EQ.1)GO TO 560
c     LSYM=ORBSYM(L+NO)
c     KLSYM=IEOR(KSYM,LSYM)
c     IF(IJSYM.NE.KLSYM) GO TO 550
c     NSGOV=NSGOV+1
c 550 CONTINUE
c 560 CONTINUE
c 570 CONTINUE
c 580 CONTINUE
c 590 CONTINUE
c
      NSHOV=0
c     DO 690 I=1,NV
c     IF(FZV(I).EQ.1)GO TO 690
c     ISYM=ORBSYM(I+NO)
c     DO 680 J=1,NO
c     IF(FZO(J).EQ.1)GO TO 680
c     JSYM=ORBSYM(J)
c     IJSYM=IEOR(ISYM,JSYM)
c     DO 670 K=1,I
c     IF(FZV(K).EQ.1)GO TO 670
c     KSYM=ORBSYM(K+NO)
c     LL=NO
c     IF(K.EQ.I)LL=J
c     DO 660 L=1,LL
c     IF(FZO(L).EQ.1)GO TO 660
c     LSYM=ORBSYM(L)
c     KLSYM=IEOR(KSYM,LSYM)
c     IF(IJSYM.NE.KLSYM) GO TO 650
c     NSHOV=NSHOV+1
c 650 CONTINUE
c 660 CONTINUE
c 670 CONTINUE
c 680 CONTINUE
c 690 CONTINUE
c
      NSLOV=0
c     DO 790 I=1,NV
c     IF(FZV(I).EQ.1)GO TO 790
c     ISYM=ORBSYM(I+NO)
c     DO 780 J=1,NO
c     IF(FZO(J).EQ.1)GO TO 780
c     JSYM=ORBSYM(J)
c     IJSYM=IEOR(ISYM,JSYM)
c     DO 770 K=1,NO
c     IF(FZO(K).EQ.1)GO TO 770
c     KSYM=ORBSYM(K)
c     DO 760 L=1,K
c     IF(FZO(L).EQ.1)GO TO 760
c     LSYM=ORBSYM(L)
c     KLSYM=IEOR(KSYM,LSYM)
c     IF(IJSYM.NE.KLSYM) GO TO 750
c     NSLOV=NSLOV+1
c 750 CONTINUE
c 760 CONTINUE
c 770 CONTINUE
c 780 CONTINUE
c 790 CONTINUE
c
      NSLVO=0
c     DO 890 I=1,NV
c     IF(FZV(I).EQ.1)GO TO 890
c     ISYM=ORBSYM(I+NO)
c     DO 880 J=1,I
c     IF(FZV(J).EQ.1)GO TO 880
c     JSYM=ORBSYM(J+NO)
c     IJSYM=IEOR(ISYM,JSYM)
c     DO 870 K=1,NV
c     IF(FZV(K).EQ.1)GO TO 870
c     KSYM=ORBSYM(K+NO)
c     DO 860 L=1,NO
c     IF(FZO(L).EQ.1)GO TO 860
c     LSYM=ORBSYM(L)
c     KLSYM=IEOR(KSYM,LSYM)
c     IF(IJSYM.NE.KLSYM) GO TO 850
c     NSLVO=NSLVO+1
c 850 CONTINUE
c 860 CONTINUE
c 870 CONTINUE
c 880 CONTINUE
c 890 CONTINUE
c
c  remove comments of if-tests for version 4
c
      NSint=0
      DO 490 I=1,Ntmo
      IF(moFZ(I).EQ.1)GO TO 490
      ISYM= moSYM(I)
      DO 480 J=1,I
      IF(moFZ(J).EQ.1)GO TO 480
      JSYM= moSYM(J)
      IJSYM=IEOR(ISYM,JSYM)
      DO 470 K=1,I
      IF(moFZ(K).EQ.1)GO TO 470
      KSYM= moSYM(K)
      LL=K
      IF(K.EQ.I)LL=J
      DO 460 L=1,LL
      IF(moFZ(L).EQ.1)GO TO 460
      LSYM= moSYM(L)
      KLSYM=IEOR(KSYM,LSYM)
      IF(IJSYM.NE.KLSYM) GO TO 460
c  exclude (vv,vv) mo ints from counting
      NSint=NSint+1
  460 CONTINUE
  470 CONTINUE
  480 CONTINUE
  490 CONTINUE
      write(6,*)
      write(6,*)'total # of ints excluding frozen orbitals', nsint    

c
      tot=0
c     DO 590 I=1,ntmo     
      DO 590 I=nvcs+1,Ntmo
         IF(moFZ(I).EQ.1)GO TO 590
         i1=f1(i) 
         ISYM= moSYM(I1)
         DO 580 J=1,I
            IF(moFZ(J).EQ.1)GO TO 580
            j1=f1(j)
            JSYM= moSYM(J1)
            IJSYM=IEOR(ISYM,JSYM)
            DO 570 K=1,I
               IF(moFZ(K).EQ.1)GO TO 570
               k1=f1(k)
               KSYM= moSYM(K1)
               LL=K
               IF(K.EQ.I)LL=J
               DO 560 L=1,LL
                  IF(moFZ(L).EQ.1)GO TO 560
                  l1=f1(l)
                  LSYM= moSYM(L1)
                  KLSYM=IEOR(KSYM,LSYM)
                  IF(IJSYM.NE.KLSYM) GO TO 560
                  tot=tot+1
  560          CONTINUE
  570       CONTINUE
  580    CONTINUE
  590 CONTINUE
      write(6,*)'minus (vv,vv)   ', tot    
c
      nsint=tot
c
      return
      END
