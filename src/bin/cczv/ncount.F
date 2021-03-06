C
C-----------------------------------------------------------------------
C
      SUBROUTINE NCOUNT(ORBSYM,FLOV,NIRRED,NO,NV,NT,FZO,FZV,NT1,NT2,
     .                  NSGOO,NSGVV,NSGOV,NSHOV,NSLOV,NSLVO,
     .                  NT3,NTAU,OPTION,DIMAR,maxo,maxv)
      IMPLICIT INTEGER (A-Z)
      COMMON/NEWPAR/NAB(8),NAI(8),NIJ(8),MAXOO,MAXVV,MAXOV
      DIMENSION FLOV(NIRRED,4),ORBSYM(NT),FZO(NO),FZV(NV),NT3(NIRRED)
      CHARACTER * 4 OPTION
C
C     COUNT THE NON-ZERO SYMMETRY AND ACTIVE T COEFFICIENTS AND
C     2E MO INTS ONLY FOR ALLOCATION PURPOSES
C
      NT1=0
      DO 120 USYM=1,NIRRED
      FU=FLOV(USYM,1)
      LU=FLOV(USYM,2)
      DO 110  U=FU,LU
      FBE=FLOV(USYM,3)-NO
      LBE=FLOV(USYM,4)-NO
      DO 100 BE=FBE,LBE
      NT1=NT1+1
  100 CONTINUE
  110 CONTINUE
  120 CONTINUE
C
      NT2=0
      DO 250 TSYM=1,NIRRED
      DO 240 USYM=1,NIRRED
      DO 230 VSYM=1,NIRRED
      UVSYM=IEOR(USYM-1,VSYM-1)+1
      IF(UVSYM.NE.TSYM) GO TO 230
      DO 225 U=FLOV(USYM,1),FLOV(USYM,2)
      DO 220 V=FLOV(VSYM,1),FLOV(VSYM,2)
      DO 218 BE=1,NV
      IF(FZV(BE).EQ.1)GO TO 218
      BESYM=ORBSYM(BE+NO)
      DO 217 GA=1,BE
      IF(FZV(GA).EQ.1)GO TO 217
      GASYM=ORBSYM(GA+NO)
      BEGSYM=IEOR(BESYM,GASYM)+1
      IF(BEGSYM.NE.TSYM)GO TO 215
      NT2=NT2+1
  215 CONTINUE
  217 CONTINUE
  218 CONTINUE
  220 CONTINUE
  225 CONTINUE
  230 CONTINUE
  240 CONTINUE
  250 CONTINUE
C
C     DIMENSION OF AR1-6
C
      MAXOO=0
      MAXVV=0
      MAXOV=0
C
      DO 20000 TSYM=1,NIRRED
      XAB=0
      DO 19490 A=1,NV
      IF(FZV(A).EQ.1)GO TO 19490
      ASYM=ORBSYM(A+NO)
      DO 19480 B=1,A
      IF(FZV(B).EQ.1)GO TO 19480
      BSYM=ORBSYM(B+NO)
      BASYM=IEOR(BSYM,ASYM)
      IF(BASYM.NE.TSYM-1)GO TO 19480
      XAB=XAB+1
      XIJ=0
      DO 19470 I=1,NO
      IF(FZO(I).EQ.1)GO TO 19470
      ISYM=ORBSYM(I)
      JSYM=IEOR(ISYM,BASYM)+1
      FJ=FLOV(JSYM,1)
      LJ=FLOV(JSYM,2)
      IF(LJ.GT.I)LJ=I
      DO 19460 J=FJ,LJ
      XIJ=XIJ+1
19460 CONTINUE
19470 CONTINUE
19480 CONTINUE
19490 CONTINUE
      IF(XIJ.GT.MAXOO)MAXOO=XIJ
      IF(XAB.GT.MAXVV)MAXVV=XAB
C     WRITE(*,*)'TSYM=',TSYM,' XIJ=',XIJ,' XAB=',XAB,'   ',XIJ*XAB
20000 CONTINUE
C
      DO 20001 TSYM=1,NIRRED
      XAI=0
      DO 15790 I=1,NO
      IF(FZO(I).EQ.1)GO TO 15790
      ISYM=ORBSYM(I)
      DO 15770 A=1,NV
      IF(FZV(A).EQ.1)GO TO 15770
      ASYM=ORBSYM(A+NO)
      IASYM=IEOR(ISYM,ASYM)
      IF(IASYM.NE.TSYM-1)GO TO 15770
      XAI=XAI+1
      XBJ=0
      DO 15750 J=1,NO
      IF(FZO(J).EQ.1)GO TO 15750
      JSYM=ORBSYM(J)
      BSYM=IEOR(JSYM,IASYM)+1
      FB=FLOV(BSYM,3)-NO
      LB=FLOV(BSYM,4)-NO
      DO 15730 B=FB,LB
      XBJ=XBJ+1
15730 CONTINUE
15750 CONTINUE
15770 CONTINUE
15790 CONTINUE
      NAI(TSYM)=XAI
      IF(XAI.GT.MAXOV)MAXOV=XAI
C     WRITE(*,*)'TSYM=',TSYM,' XBJ=',XBJ,' XAI=',XAI,'   ',XAI*XBJ
C
20001 CONTINUE
C
      WRITE(6,*)' '
      WRITE(*,*)'MAXOV=',MAXOV
      WRITE(*,*)'MAXOO=',MAXOO,' MAXVV=',MAXVV
      DIMAR=MAX0(MAXOO*MAXVV,MAXOV*MAXOV)
      WRITE(*,*)'DIMAR=', DIMAR
      WRITE(6,*)' '
 
      IF(OPTION.EQ.'SDT1')THEN
      WRITE(6,*)
      WRITE(6,*)'VOV ARRAYS FOR CCSDT-1A'
      MAXCNT=0
      MAXO=0
      MAXV=0
      DO 2390 TSYM=1,NIRRED
      NOSYM=FLOV(TSYM,2)-FLOV(TSYM,1)+1
      NVSYM=FLOV(TSYM,4)-FLOV(TSYM,3)+1
      IF(MAXO.LT.NOSYM)MAXO=NOSYM
      IF(MAXV.LT.NVSYM)MAXV=NVSYM
      ICNT=0
      DO 2380 A=1,NV
      IF(FZV(A).EQ.1)GO TO 2380
      ASYM=ORBSYM(A+NO)
      DO 2370 I=1,NO
      IF(FZO(I).EQ.1)GO TO 2370
      AISYM=IEOR(ASYM,ISYM)
      BSYM=IEOR(TSYM-1,AISYM)+1
      FB=FLOV(BSYM,3)-NO
      LB=FLOV(BSYM,4)-NO
      DO 2360 B=FB,LB
      ICNT=ICNT+1
 2360 CONTINUE
 2370 CONTINUE
 2380 CONTINUE
      IF(MAXCNT.LT.ICNT)MAXCNT=ICNT
      WRITE(6,*)'TSYM=',TSYM,'  ICNT=',ICNT
 2390 CONTINUE
      IF(MAXCNT.LT.MAXOV*MAXV)THEN
      WRITE(6,*)'ERROR  IN INE NCOUNT, MAXOV*MAXV.GT.MAXOV*MAXOV'
      STOP
      ENDIF
      WRITE(6,*)
      WRITE(6,*)'MAXO=',MAXO,'  MAXV=',MAXV
      WRITE(6,*)
      ENDIF
C
C     COUNT SYMMETRY NON-ZERO T3
C
      NTAU=NT2
      IF(OPTION.EQ.'SDT1')THEN
      WRITE(6,*)
      DO 440 KCSYM=1,NIRRED
      NPFF=0
      DO 430 U=1,NO
      IF(FZO(U).EQ.1)GO TO 430
      USYM=ORBSYM(U)
      DO 420 V=1,NO
      IF(FZO(V).EQ.1)GO TO 420
      VSYM=ORBSYM(V)
      UVSYM=IEOR(USYM,VSYM)
      DO 418 BE=1,NV
      IF(FZV(BE).EQ.1)GO TO 418
      BESYM=ORBSYM(BE+NO)
C     DO 417 GA=1,NV
      DO 417 GA=1,BE
      IF(FZV(GA).EQ.1)GO TO 417
      GASYM=ORBSYM(GA+NO)
      BEGSYM=IEOR(BESYM,GASYM)
      TSYM=IEOR(UVSYM,BEGSYM)
      IF((KCSYM-1).NE.TSYM)GO TO 415
      NPFF=NPFF+1
C     WRITE(*,657)U,V,BE,GA,NPFF
  657 FORMAT(' U=',I3,' V=',I3,' BE=',I3,' GA=',I3,'NPFF= ',I6)
  415 CONTINUE
  417 CONTINUE
  418 CONTINUE
  420 CONTINUE
  430 CONTINUE
      NT3(KCSYM)=NPFF
  440 CONTINUE
      WRITE(6,7)NT3(1)
    7 FORMAT(1X,'NT3(1) =',I6)
      IF(NT3(1).NE.NT2)WRITE(*,*)'WARNING ...   NT3(1).NE.NT2'
      NTAU=NT3(1)
      DO 455 I=2,NIRRED
      IF(NT3(I).GT.NT3(1))THEN
      WRITE(*,*)'WARNING !!!  SOME NT3.GT.NT2'
      NTAU=NT3(I)
      ENDIF
      WRITE(6,8)I,NT3(I)
    8 FORMAT(1X,'NT3(',I1,') =',I6)
  455 CONTINUE
      ENDIF
C
      NSGOO=0
      DO 390 I=1,NO
      IF(FZO(I).EQ.1)GO TO 390
      ISYM=ORBSYM(I)
      DO 380 J=1,I
      IF(FZO(J).EQ.1)GO TO 380
      JSYM=ORBSYM(J)
      IJSYM=IEOR(ISYM,JSYM)
      DO 370 K=1,I
      IF(FZO(K).EQ.1)GO TO 370
      KSYM=ORBSYM(K)
      LL=K
      IF(K.EQ.I)LL=J
      DO 360 L=1,LL
      IF(FZO(L).EQ.1)GO TO 360
      LSYM=ORBSYM(L)
      KLSYM=IEOR(KSYM,LSYM)
      IF(IJSYM.NE.KLSYM) GO TO 350
      NSGOO=NSGOO+1
  350 CONTINUE
  360 CONTINUE
  370 CONTINUE
  380 CONTINUE
  390 CONTINUE
C
      NSGVV=0
      DO 490 I=1,NV
      IF(FZV(I).EQ.1)GO TO 490
      ISYM=ORBSYM(I+NO)
      DO 480 J=1,I
      IF(FZV(J).EQ.1)GO TO 480
      JSYM=ORBSYM(J+NO)
      IJSYM=IEOR(ISYM,JSYM)
      DO 470 K=1,I
      IF(FZV(K).EQ.1)GO TO 470
      KSYM=ORBSYM(K+NO)
      LL=K
      IF(K.EQ.I)LL=J
      DO 460 L=1,LL
      IF(FZV(L).EQ.1)GO TO 460
      LSYM=ORBSYM(L+NO)
      KLSYM=IEOR(KSYM,LSYM)
      IF(IJSYM.NE.KLSYM) GO TO 450
      NSGVV=NSGVV+1
  450 CONTINUE
  460 CONTINUE
  470 CONTINUE
  480 CONTINUE
  490 CONTINUE
C
      NSGOV=0
      DO 590 I=1,NO
      IF(FZO(I).EQ.1)GO TO 590
      ISYM=ORBSYM(I)
      DO 580 J=1,I
      IF(FZO(J).EQ.1)GO TO 580
      JSYM=ORBSYM(J)
      IJSYM=IEOR(ISYM,JSYM)
      DO 570 K=1,NV
      IF(FZV(K).EQ.1)GO TO 570
      KSYM=ORBSYM(K+NO)
      DO 560 L=1,K
      IF(FZV(L).EQ.1)GO TO 560
      LSYM=ORBSYM(L+NO)
      KLSYM=IEOR(KSYM,LSYM)
      IF(IJSYM.NE.KLSYM) GO TO 550
      NSGOV=NSGOV+1
  550 CONTINUE
  560 CONTINUE
  570 CONTINUE
  580 CONTINUE
  590 CONTINUE
C
      NSHOV=0
      DO 690 I=1,NV
      IF(FZV(I).EQ.1)GO TO 690
      ISYM=ORBSYM(I+NO)
      DO 680 J=1,NO
      IF(FZO(J).EQ.1)GO TO 680
      JSYM=ORBSYM(J)
      IJSYM=IEOR(ISYM,JSYM)
      DO 670 K=1,I
      IF(FZV(K).EQ.1)GO TO 670
      KSYM=ORBSYM(K+NO)
      LL=NO
      IF(K.EQ.I)LL=J
      DO 660 L=1,LL
      IF(FZO(L).EQ.1)GO TO 660
      LSYM=ORBSYM(L)
      KLSYM=IEOR(KSYM,LSYM)
      IF(IJSYM.NE.KLSYM) GO TO 650
      NSHOV=NSHOV+1
  650 CONTINUE
  660 CONTINUE
  670 CONTINUE
  680 CONTINUE
  690 CONTINUE
C
      NSLOV=0
      DO 790 I=1,NV
      IF(FZV(I).EQ.1)GO TO 790
      ISYM=ORBSYM(I+NO)
      DO 780 J=1,NO
      IF(FZO(J).EQ.1)GO TO 780
      JSYM=ORBSYM(J)
      IJSYM=IEOR(ISYM,JSYM)
      DO 770 K=1,NO
      IF(FZO(K).EQ.1)GO TO 770
      KSYM=ORBSYM(K)
      DO 760 L=1,K
      IF(FZO(L).EQ.1)GO TO 760
      LSYM=ORBSYM(L)
      KLSYM=IEOR(KSYM,LSYM)
      IF(IJSYM.NE.KLSYM) GO TO 750
      NSLOV=NSLOV+1
  750 CONTINUE
  760 CONTINUE
  770 CONTINUE
  780 CONTINUE
  790 CONTINUE
C
      NSLVO=0
      DO 890 I=1,NV
      IF(FZV(I).EQ.1)GO TO 890
      ISYM=ORBSYM(I+NO)
      DO 880 J=1,I
      IF(FZV(J).EQ.1)GO TO 880
      JSYM=ORBSYM(J+NO)
      IJSYM=IEOR(ISYM,JSYM)
      DO 870 K=1,NV
      IF(FZV(K).EQ.1)GO TO 870
      KSYM=ORBSYM(K+NO)
      DO 860 L=1,NO
      IF(FZO(L).EQ.1)GO TO 860
      LSYM=ORBSYM(L)
      KLSYM=IEOR(KSYM,LSYM)
      IF(IJSYM.NE.KLSYM) GO TO 850
      NSLVO=NSLVO+1
  850 CONTINUE
  860 CONTINUE
  870 CONTINUE
  880 CONTINUE
  890 CONTINUE
      END
