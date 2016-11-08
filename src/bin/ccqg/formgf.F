      SUBROUTINE FORMGF(BUF,IBUF,GBUF,TAU,T2,T1,Z2,Z1,NO,VALOV,
     .                  VALOV2,NV,NO2,NV2,NTRO,NTRV,ICNT,ZLX,VALV2,
     .                  VALOV3,VALO3V,IOFF,NTR,VALV,JOUT,ITAP68,OPTION,
     .                  ORBSYM,FLOV,T3Z2,ITR,ITV,NT3,NIRRED,NTAU,NOV,
     .                  T3OFF,T3ADD,W3,Y3,NDIMT2,TOFF,TADD,FZO,FZV)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*4 OPTION
      COMMON/TAPES/ITAP81,ITAP82,ITAP83,ITAP99,itap76,i76,itap74
      COMMON/BUFFER/TOL,IBFLEN,MAXBUF,IBOFF
      REAL*8 VAL,VAL1,VAL2,VAL3,VAL4,VAL5,VAL6,VAL7,VAL8,VAL9,
     .       TOL,GFVAL,FACBC
      REAL*8 BUF(GBUF),Z2(NDIMT2),T2(NDIMT2),VALV(NV),VALV2(NTRV),
     .       T1(NO,NV),Z1(NO,NV),TAU(NDIMT2),VALO3V(NV,NO,NO,NO),
     .       VALOV(NV,NO),VALOV2(NV,NV,NO),VALOV3(NV,NV,NO),
     .       T3Z2(NO,NV,NV),W3(NTAU),Y3(NTAU)
      INTEGER IOFF(NTR),IBUF(2*GBUF),ZLX(NV,NV),
     .        ORBSYM(NO+NV),FLOV(NIRRED,4),FZO(NO),FZV(NV),
     .        ITR(NO),ITV(NOV),NT3(NIRRED),
     .        T3OFF(NO,NO,2,NIRRED),T3ADD(NV,NV,NIRRED),
     .        TOFF(NO,NO,2),TADD(NV,NV)
C
      CALL ZERO(VALV2,NTRV)
      CALL ZERO(VALV,NV)
      CALL ZERO(VALOV,NV*NO)
      CALL ZERO(VALOV2,NV*NV*NO)
      CALL ZERO(VALOV3,NV*NV*NO)
      CALL ZERO(VALO3V,NO2*NO*NV)
      DO 158 B=1,NV
         IF(FZV(B).EQ.1)GO TO 158
         BSYM=ORBSYM(B+NO)
         FC=FLOV(BSYM+1,3)-NO
         TC=FLOV(BSYM+1,4)-NO
         IF(TC.GT.B)TC=B
         DO 18 C=FC,TC
            BC=IOFF(B)+C
            VAL=0.0D+00
            DO 19 D=1,NV
               IF(FZV(D).EQ.1)GO TO 19
               DSYM=ORBSYM(D+NO)
               BDSYM=IEOR(BSYM,DSYM)
               ZLB=ZLX(B,D)
               ZLC=ZLX(C,D)
               DO 20 J=1,NO
                  IF(FZO(J).EQ.1)GO TO 20
                  JSYM=ORBSYM(J)
                  LSYM=IEOR(BDSYM,JSYM)+1
                  FL=FLOV(LSYM,1)
                  TL=FLOV(LSYM,2)
                  DO 17 L=FL,TL
                     JLCD=TOFF(J,L,ZLC)+TADD(C,D)
                     JLBD=TOFF(J,L,ZLB)+TADD(B,D)
                     VAL=VAL+T2(JLBD)*Z2(JLCD)+T2(JLCD)*Z2(JLBD)
   17             CONTINUE
   20          CONTINUE
   19       CONTINUE
            VALV2(BC)=VAL+VAL
   18    CONTINUE
  158 CONTINUE
C
      DO 35 B=1,NV
         IF(FZV(B).EQ.1)GO TO 35
         BSYM=ORBSYM(B+NO)
         DO 301 I=1,NO
            IF(FZO(I).EQ.1)GO TO 301
            ISYM=ORBSYM(I)
            FD=FLOV(ISYM+1,3)-NO
            TD=FLOV(ISYM+1,4)-NO
            BISYM=IEOR(BSYM,ISYM)
            DO 302 J=1,NO
               IF(FZO(J).EQ.1)GO TO 302
               JSYM=ORBSYM(J)
               LSYM=IEOR(BISYM,JSYM)+1
               FL=FLOV(LSYM,1)
               TL=FLOV(LSYM,2)
               DO 32 L=FL,TL
                  VAL=0.0D+00
                  DO 31 D=FD,TD
                     ZL=ZLX(B,D)
                     LJBD=TOFF(L,J,ZL)+TADD(B,D)
                     VAL=VAL-Z2(LJBD)*T1(I,D)
   31             CONTINUE
                  VALO3V(B,I,L,J)=VAL
   32          CONTINUE
  302       CONTINUE
  301    CONTINUE
   35 CONTINUE
C
      I81=1
      J81=1
      Y81=1
      i74=1
      call rfile(itap74)
      DO 10 A=1,NV
         IF(FZV(A).EQ.1)GO TO 10
         ASYM=ORBSYM(A+NO)
         FC=FLOV(ASYM+1,3)-NO
         TC=FLOV(ASYM+1,4)-NO
         CALL ZERO(T3Z2,NV*NV*NO)
c
         if(option.eq.'CCT ')then
c        write(6,*)' reading fden for mo',a+no      
         call wreadw(itap74,t3z2,intowp(nv*nv*no),i74,i74)
         do 1009 i9=1,no
         do 1008 a9=1,nv
         do 1007 b9=1,nv
            t3z2(i9,a9,b9)=-t3z2(i9,a9,b9)
 1007    continue
 1008    continue
 1009    continue
c        write(6,'(8f10.6)')t3z2
         endif
c
         DO 21 C=FC,TC
            VAL=0.0D+00
            DO 22 D=1,NV
               IF(FZV(D).EQ.1)GO TO 22
               DSYM=ORBSYM(D+NO)
               ADSYM=IEOR(ASYM,DSYM)
               ZLA=ZLX(A,D)
               ZLC=ZLX(C,D)
               DO 422 J=1,NO
                  IF(FZO(J).EQ.1)GO TO 422
                  JSYM=ORBSYM(J)
                  LSYM=IEOR(ADSYM,JSYM)+1
                  FL=FLOV(LSYM,1)
                  TL=FLOV(LSYM,2)
                  DO 423 L=FL,TL
                     JLAD=TOFF(J,L,ZLA)+TADD(A,D)
                     JLCD=TOFF(J,L,ZLC)+TADD(C,D)
                     VAL=VAL-T2(JLAD)*Z2(JLCD)
  423             CONTINUE
  422          CONTINUE
   22       CONTINUE
            VALV(C)=VAL
   21    CONTINUE
C
         FJ=FLOV(ASYM+1,1)
         TJ=FLOV(ASYM+1,2)
         DO 421 C=1,NV
            IF(FZV(C).EQ.1)GO TO 421
            CSYM=ORBSYM(C+NO)
            ACSYM=IEOR(ASYM,CSYM)
            DO 24 D=1,NV
               IF(FZV(D).EQ.1)GO TO 24
               DSYM=ORBSYM(D+NO)
               LSYM=IEOR(ACSYM,DSYM)+1
               FL=FLOV(LSYM,1)
               TL=FLOV(LSYM,2)
               ZLC=ZLX(C,D)
               DO 424 L=FL,TL
                  VAL5=0.0D+00
                  VAL6=0.0D+00
                  DO 425 J=FJ,TJ
                     JLCD=TOFF(J,L,ZLC)+TADD(C,D)
                     LJCD=TOFF(L,J,ZLC)+TADD(C,D)
                     VAL5=VAL5-Z2(JLCD)*T1(J,A)
                     VAL6=VAL6-Z2(LJCD)*T1(J,A)
  425             CONTINUE
                  VALOV2(C,D,L)=VAL5
                  VALOV3(C,D,L)=VAL6
  424          CONTINUE
   24       CONTINUE
  421    CONTINUE
C
         IF(OPTION.EQ.'SDT1')THEN
         DO 190 K=1,NO
            IF(FZO(K).EQ.1)GO TO 190
            KSYM=ORBSYM(K)
            KASYM=IEOR(ASYM,KSYM)+1
            DIM=NT3(KASYM)
            CALL WREADW(ITAP81,W3,INTOWP(DIM),I81,I81)
            DO 185 I=1,NO
               IF(FZO(I).EQ.1)GO TO 185
               ISYM=ORBSYM(I)
               AISYM=IEOR(ASYM,ISYM)
               DO 182 B=1,NV
                  IF(FZV(B).EQ.1)GO TO 182
                  BSYM=ORBSYM(B+NO)
                  CSYM=IEOR(AISYM,BSYM)
                  FC=FLOV(CSYM+1,3)-NO
                  LC=FLOV(CSYM+1,4)-NO
                  LC2=LC
                  IF(LC.GT.B)LC2=B
                  DO 180 C=FC,LC2
                     DO 175 J=1,NO
                        IF(FZO(J).EQ.1)GO TO 175
                        JSYM=ORBSYM(J)
                        JBSYM=IEOR(JSYM,BSYM)
                        DSYM=IEOR(JBSYM,KSYM)
                        FD=FLOV(DSYM+1,3)-NO
                        LD=FLOV(DSYM+1,4)-NO
                        DO 170 D=FD,LD
                           ZLBD=ZLX(B,D)
                           ZLDC=ZLX(D,C)
                           JKBD=TOFF(J,K,ZLBD)+TADD(B,D)
                           KJBD=TOFF(K,J,ZLBD)+TADD(B,D)
                           IJDC=T3OFF(I,J,ZLDC,KASYM)+T3ADD(D,C,KASYM)
                           JIDC=T3OFF(J,I,ZLDC,KASYM)+T3ADD(D,C,KASYM)
                           VAL=W3(IJDC)*Z2(JKBD)+W3(JIDC)*Z2(KJBD)
                           T3Z2(I,B,C)=T3Z2(I,B,C)+VAL
  170                   CONTINUE
                        JCSYM=IEOR(JSYM,CSYM)
                        DSYM=IEOR(JCSYM,KSYM)
                        FD=FLOV(DSYM+1,3)-NO
                        LD=FLOV(DSYM+1,4)-NO
                        DO 171 D=FD,LD
                           ZLCD=ZLX(C,D)
                           ZLDB=ZLX(D,B)
                           JKCD=TOFF(J,K,ZLCD)+TADD(C,D)
                           KJCD=TOFF(K,J,ZLCD)+TADD(C,D)
                           IJDB=T3OFF(I,J,ZLDB,KASYM)+T3ADD(D,B,KASYM)
                           JIDB=T3OFF(J,I,ZLDB,KASYM)+T3ADD(D,B,KASYM)
                           VAL=W3(IJDB)*Z2(JKCD)+W3(JIDB)*Z2(KJCD)
                           T3Z2(I,B,C)=T3Z2(I,B,C)+VAL
  171                   CONTINUE
  175                CONTINUE
  180             CONTINUE
  182          CONTINUE
  185       CONTINUE
  190    CONTINUE
      ENDIF
C
      DO 11 I=1,NO
         IF(FZO(I).EQ.1)GO TO 11
         ISYM=ORBSYM(I)
         IASYM=IEOR(ASYM,ISYM)+1
         IF(OPTION.EQ.'SDT1')THEN
            DIM=NT3(IASYM)
            CALL WREADW(ITAP81,W3,INTOWP(DIM),J81,J81)
            CALL WREADW(ITAP82,Y3  ,INTOWP(DIM),Y81,Y81)
         ENDIF
         DO 26 C=1,NV
            IF(FZV(C).EQ.1)GO TO 26
            CSYM=ORBSYM(C+NO)
            JSYM=IEOR(CSYM,IASYM-1)+1
            FJ=FLOV(JSYM,1)
            TJ=FLOV(JSYM,2)
            DO 126 J=FJ,TJ
               VAL4=0.0D+00
               DO 27 D=1,NV
                  IF(FZV(D).EQ.1)GO TO 27
                  DSYM=ORBSYM(D+NO)
                  LSYM=IEOR(DSYM,IASYM-1)+1
                  FL=FLOV(LSYM,1)
                  TL=FLOV(LSYM,2)
                  ZLA=ZLX(A,D)
                  ZLC=ZLX(C,D)
                  DO 127 L=FL,TL
                     ILAD=TOFF(I,L,ZLA)+TADD(A,D)
                     LIAD=TOFF(L,I,ZLA)+TADD(A,D)
                     JLCD=TOFF(J,L,ZLC)+TADD(C,D)
                     VAL4=VAL4+(T2(ILAD)+T2(ILAD)-T2(LIAD))*Z2(JLCD)
  127             CONTINUE
   27          CONTINUE
               VALOV(C,J)=VAL4
  126       CONTINUE
   26    CONTINUE
         DO 12 B=1,NV
            IF(FZV(B).EQ.1)GO TO 12
            BSYM=ORBSYM(B+NO)
            IBSYM=IEOR(ISYM ,BSYM)
            ABSYM=IEOR(ASYM ,BSYM)
            AB=IOFF(MAX0(A,B))+MIN0(A,B)
            ZLAB=ZLX(A,B)
            CSYM=IEOR(IASYM-1,BSYM)
            FC=FLOV(CSYM+1,3)-NO
            LC=FLOV(CSYM+1,4)-NO
            IF(LC.GT.B)LC=B
            DO 13 C=FC,LC
               BC=IOFF(B)+C
               ZLAC=ZLX(A,C)
               VAL1=0.0D+00
               VAL2=0.0D+00
               VAL3=0.0D+00
               VAL4=0.0D+00
               VAL5=0.0D+00
               VAL7=0.0D+00
               VAL8=0.0D+00
               VAL9=0.0D+00
               FJ=FLOV(BSYM+1,1)
               LJ=FLOV(BSYM+1,2)
               DO 14 J=FJ,LJ
                  IJAC=TOFF(I,J,ZLAC)+TADD(A,C)
                  JIAC=TOFF(J,I,ZLAC)+TADD(A,C)
                  VAL1=VAL1+TAU(IJAC)*Z1(J,B)
                  VAL2=VAL2+TAU(JIAC)*Z1(J,B)
                  VAL3=VAL3+T1(J,B)*Z2(IJAC)
                  VAL4=VAL4+VALOV(C,J)*T1(J,B)
   14          CONTINUE
               FJ=FLOV(CSYM+1,1)
               LJ=FLOV(CSYM+1,2)
               DO 114 J=FJ,LJ
                  IJAB=TOFF(I,J,ZLAB)+TADD(A,B)
                  JIAB=TOFF(J,I,ZLAB)+TADD(A,B)
                  VAL1=VAL1+TAU(IJAB)*Z1(J,C)
                  VAL2=VAL2+TAU(JIAB)*Z1(J,C)
                  VAL3=VAL3+T1(J,C)*Z2(IJAB)
                  VAL4=VAL4+VALOV(B,J)*T1(J,C)
  114          CONTINUE
               DO 214 J=1,NO
                  IF(FZO(J).EQ.1)GO TO 214
                  JSYM=ORBSYM(J)
                  DSYM=IEOR(JSYM,IBSYM)+1
                  FD=FLOV(DSYM,3)-NO
                  LD=FLOV(DSYM,4)-NO
                  DO 28 D=FD,LD
                     ZLB=ZLX(B,D)
                     IJBD=TOFF(I,J,ZLB)+TADD(B,D)
                     JIBD=TOFF(J,I,ZLB)+TADD(B,D)
         VAL5=VAL5+VALOV2(C,D,J)*T2(IJBD)+VALOV3(C,D,J)*T2(JIBD)
   28             CONTINUE
                  DSYM=IEOR(JSYM,ABSYM)+1
                  FD=FLOV(DSYM,3)-NO
                  LD=FLOV(DSYM,4)-NO
                  DO 128 D=FD,LD
                     ZLC=ZLX(C,D)
                     IJCD=TOFF(I,J,ZLC)+TADD(C,D)
                     JICD=TOFF(J,I,ZLC)+TADD(C,D)
         VAL5=VAL5+VALOV2(B,D,J)*T2(IJCD)+VALOV3(B,D,J)*T2(JICD)
  128             CONTINUE
                  LSYM=IEOR(JSYM,IBSYM)+1
                  FL=FLOV(LSYM,1)
                  LL=FLOV(LSYM,2)
                  DO 33 L=FL,LL
                     JLAC=TOFF(J,L,ZLAC) +TADD(A,C)
                     VAL7=VAL7+VALO3V(B,I,L,J)*TAU(JLAC)
   33             CONTINUE
                  LSYM=IEOR(JSYM,ABSYM)+1
                  FL=FLOV(LSYM,1)
                  LL=FLOV(LSYM,2)
                  DO 133  L=FL,LL
                     JLAB=TOFF(J,L,ZLAB) +TADD(A,B)
                     VAL7=VAL7+VALO3V(C,I,L,J)*TAU(JLAB)
  133             CONTINUE
  214          CONTINUE
               IF(OPTION.EQ.'SDT1')THEN
                  DO 61 J=1,NO
                     IF(FZO(J).EQ.1)GO TO 61
                     JSYM=ORBSYM(J)
                     DO 62 K=1,NO
                        IF(FZO(K).EQ.1)GO TO 62
                        KSYM=ORBSYM(K)
                        JKSYM=IEOR(JSYM,KSYM)
                        DSYM=IEOR(JKSYM,CSYM)+1
                        FD=FLOV(DSYM,3)-NO
                        LD=FLOV(DSYM,4)-NO
                        DO 63 D=FD,LD
                           ZLCD=ZLX(C,D)
                           JKCD=TOFF(J,K,ZLCD)+TADD(C,D)
                           ZLBD=ZLX(B,D)
                           JKBD=T3OFF(J,K,ZLBD,IASYM)+T3ADD(B,D,IASYM)
                           VAL8=VAL8+Y3(JKBD)*T2(JKCD)
   63                   CONTINUE
                        DSYM=IEOR(JKSYM,BSYM)+1
                        FD=FLOV(DSYM,3)-NO
                        LD=FLOV(DSYM,4)-NO
                        DO 64 D=FD,LD
                           ZLBD=ZLX(B,D)
                           ZLCD=ZLX(C,D)
                           JKBD=TOFF(J,K,ZLBD)  +TADD(B,D)
                           JKCD=T3OFF(J,K,ZLCD,IASYM)+T3ADD(C,D,IASYM)
                           VAL8=VAL8+Y3(JKCD)*T2(JKBD)
   64                   CONTINUE
   62                CONTINUE
   61             CONTINUE
C
                  DO 314 J=1,NO
                     IF(FZO(J).EQ.1)GO TO 314
                     JSYM=ORBSYM(J)
                     JBSYM=IEOR(JSYM,BSYM)
                     JCSYM=IEOR(JSYM,CSYM)
                     DO 328 D=1,NV
                        IF(FZV(D).EQ.1)GO TO 328
                        DSYM=ORBSYM(D+NO)
                        ZLB=ZLX(B,D)
                        ZLC=ZLX(C,D)
                        KSYM=IEOR(JBSYM,DSYM)+1
                        FK=FLOV(KSYM,1)
                        LK=FLOV(KSYM,2)
                        DO 288 K=FK,LK
                           JKBD=TOFF(J,K,ZLB)+TADD(B,D)
                           JKCD=T3OFF(J,K,ZLC,IASYM)+T3ADD(C,D,IASYM)
                           VAL9=VAL9+W3(JKCD)*Z2(JKBD)
  288                   CONTINUE
                        KSYM=IEOR(JCSYM,DSYM)+1
                        FK=FLOV(KSYM,1)
                        LK=FLOV(KSYM,2)
                        DO 289 K=FK,LK
                           JKCD=TOFF(J,K,ZLC)+TADD(C,D)
                           JKBD=T3OFF(J,K,ZLB,IASYM)+T3ADD(B,D,IASYM)
                           VAL9=VAL9+W3(JKBD)*Z2(JKCD)
  289                   CONTINUE
  328                CONTINUE
  314             CONTINUE
               ENDIF
               GFVAL=(VAL1+VAL1-VAL2+VAL3)+VALV2(BC)*T1(I,A)+VAL7
     .               +VALV(C)*T1(I,B)+VALV(B)*T1(I,C)+VAL4+VAL5
     .               -VAL8+VAL9+VAL9-T3Z2(I,B,C)
c              WRITE (JOUT,1) B+NO,C+NO,A+NO,I,GFVAL/8.0D+00
               IF(DABS(GFVAL).GT.TOL) THEN
                  LBL=IOR(I,ISHFT(A+NO,8))
                  LBL=IOR(B+NO,ISHFT(LBL,8))
                  LBL=IOR(C+NO,ISHFT(LBL,8))
                  IF(ICNT.GT.MAXBUF) THEN
                     CALL SWRIT(ITAP68,IBUF,INTOWP(GBUF))
                     CALL ZERO(BUF,GBUF)
                      ICNT=1
                  ENDIF
                  IBUF(1)=0
                  IBUF(2)=ICNT
                  IBUF(2+ICNT)=LBL
                  BUF(ICNT+IBOFF)=GFVAL
                  ICNT=ICNT+1
               ENDIF
   13          CONTINUE
   12       CONTINUE
C
   11    CONTINUE
   10 CONTINUE
C
    1 FORMAT ('GFMAT(',I3,',',I3,',',I3,',',I3,')=',F20.12)
      RETURN
      END
