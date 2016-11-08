C
C-----------------------------------------------------------------------
C
      SUBROUTINE BACK2 (VEC,UVEC,DMAT,IMAT,GDM,IOFF,IADD,
     .                  ORBSYM,NIRRED,FLOV,ITR,ITV,IPQ,NO,NV,NT,NT2,
     .                  NTAO,NTAO2,AOPS,IJOFF,KLOFF,LNSR,MAPIJ,
     .                  KTYPE,KLOC,NSHELL,ITAP55,AOSTRT,
     .                  CCTSCF,SCFTCC,ZETA,CORE,CORE2,
     .                  WVEC,BZET,QMAT,FZT,QCOR,QACT,QDM,NFZO,NFZV,
     .                  ZETO,ZETV,IFZO,IFZV,ITAP83,SRTOUT,NRS,AUX3,
     .                  LCOR,LEFT,ic1)
      IMPLICIT INTEGER (A-Z)
      REAL*8 VAL,Z0,FAC,VHF,VCR,VCR2
      REAL*8 VEC(NT,NT),UVEC(NTAO,NT),WVEC(NTAO,NT),QCOR(NTAO,NTAO),
     .       GDM(*),DMAT(NTAO,NTAO),IMAT(NTAO,NTAO),QACT(NTAO,NTAO),
     .       QDM(*),BUF(36*36),ZETA(*),BZET(NTAO,NTAO),
     .       QMAT(NTAO,NTAO),ZETO(*),ZETV(*),AUX3(*),LCOR(LEFT)
      INTEGER ORBSYM(*),FLOV(NIRRED,2),ITR(*),ITV(*),IPQ(*),
     .        IOFF(*),IADD(*), FZT(NT),IFZO(*),IFZV(*),
     .        AOPS(*),KTYPE(*),MAPIJ(*),
     .        KLOC(*),IJOFF(*),KLOFF(*),LNSR(*),NIJ(8),
     .        CCTSCF(*),SCFTCC(*),NFZO(*),NFZV(*),NRS(*),NCNTIN(8)
      REAL*8 CORE(NT,NT,NT,NT),CORE2(NTAO,NTAO,NT,NT)
      CHARACTER*4 SRTOUT
      PARAMETER (Z0=0.0D0)
C
      ITAP95=95
      ITAP96=96
      ITAP97=97
      ITAP98=98
C
C     WRITE(6,*)' MO TO SO SCF VECTOR'
C     CALL MATOUT ( VEC,NT,NT,NT,NT,6)
C     WRITE(6,*)' SO TO AO TRANSFORMATION MATRIX'
C     CALL MATOUT (UVEC,NTAO,NT,NTAO,NT,6)
C     WRITE(6,*)' MO TO AO SCF VECTOR'
C     CALL MATOUT (WVEC,NTAO,NT,NTAO,NT,6)
C
C     WRITE(6,*)' QDM IN BACK2'
C     CALL PRINT (QDM,NT2,NT,6)
C
C >>> BACK-TRANSFORM Z VECTOR TO AO BASIS
C
      CALL ZERO (IMAT,NTAO*NTAO)
      CALL ZERO (DMAT,NTAO*NTAO)
C
      DO 130 I=1,NO
         II=CCTSCF(I)
         ISYM=ORBSYM(II)
         DO 120 A=1,NV
            AA=CCTSCF(A+NO)
            ASYM=ORBSYM(AA)
            AI=ITR(A)+I
            IF(ISYM.NE.ASYM)GO TO 120
            DO 110 MU=1,NTAO
               IMAT(MU,II)=IMAT(MU,II)+WVEC(MU,AA)*ZETA(AI)
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
C
      DO 160 MU=1,NTAO
         DO 150 NU=1,NTAO
            DO 140 I=1,NO
               II=CCTSCF(I)
               DMAT(MU,NU)=DMAT(MU,NU)+WVEC(NU,II)*IMAT(MU,II)
  140       CONTINUE
  150    CONTINUE
  160 CONTINUE
C
      CALL ZERO (IMAT,NTAO*NTAO)
C
      DO 132 I=1,NO
         II=CCTSCF(I)
         IF(FZT(II).NE.1)GO TO 132
         ISYM=ORBSYM(II)
         DO 122 J=1,NO
            JJ=CCTSCF(J)
            IF(FZT(JJ).EQ.1)GO TO 122
            JSYM=ORBSYM(JJ)
            IF(ISYM.NE.JSYM)GO TO 122
            IJZ=ITR(I)+J
            ZIJ=IFZO(IJZ)
            DO 112 MU=1,NTAO
               IMAT(MU,II)=IMAT(MU,II)+WVEC(MU,JJ)*ZETO(ZIJ)
  112       CONTINUE
  122    CONTINUE
  132 CONTINUE
C
      DO 162 MU=1,NTAO
         DO 152 NU=1,NTAO
            DO 142 I=1,NO
               II=CCTSCF(I)
               IF(FZT(II).NE.1)GO TO 142
               DMAT(MU,NU)=DMAT(MU,NU)+WVEC(NU,II)*IMAT(MU,II)
  142       CONTINUE
  152    CONTINUE
  162 CONTINUE
C
      CALL ZERO (IMAT,NTAO*NTAO)
C
      DO 133 A=1,NV
         AA=CCTSCF(A+NO)
         IF(FZT(AA).NE.1)GO TO 133
         ASYM=ORBSYM(AA)
         DO 123 B=1,NV
            BB=CCTSCF(B+NO)
            IF(FZT(BB).EQ.1)GO TO 123
            BSYM=ORBSYM(BB)
            IF(ASYM.NE.BSYM)GO TO 123
            ABZ=ITV(A)+B
            ZAB=IFZV(ABZ)
            DO 113 MU=1,NTAO
               IMAT(MU,AA)=IMAT(MU,AA)+WVEC(MU,BB)*ZETV(ZAB)
  113       CONTINUE
  123    CONTINUE
  133 CONTINUE
C
      DO 163 MU=1,NTAO
         DO 153 NU=1,NTAO
            DO 143 A=1,NV
               AA=CCTSCF(A+NO)
               IF(FZT(AA).NE.1)GO TO 143
               DMAT(MU,NU)=DMAT(MU,NU)+WVEC(NU,AA)*IMAT(MU,AA)
  143       CONTINUE
  153    CONTINUE
  163 CONTINUE
C
      DO 180 MU=1,NTAO
         DO 170 NU=1,NTAO
            BZET(MU,NU)=0.5D0*(DMAT(MU,NU)+DMAT(NU,MU))
  170    CONTINUE
  180 CONTINUE
C
C >>> BACK-TRANSFORM ACTIVE ONE-PDM TO AO BASIS
C
      CALL ZERO (IMAT,NTAO*NTAO)
      DO 131 I=1,NT
         IF(FZT(I).EQ.1)GO TO 131
         ISYM=ORBSYM(I)
         DO 121 J=1,NT
            IF(FZT(J).EQ.1)GO TO 121
            JJ=SCFTCC(J)
            JSYM=ORBSYM(J)
            IF(ISYM.NE.JSYM)GO TO 121
            IJ=IPQ(MAX0(I,J))+MIN0(I,J)
            DO 111 MU=1,NTAO
               IMAT(MU,I)=IMAT(MU,I)+WVEC(MU,J)*QDM(IJ)
  111       CONTINUE
  121    CONTINUE
  131 CONTINUE
C
      CALL ZERO (QACT,NTAO*NTAO)
      DO 161 MU=1,NTAO
         DO 151 NU=1,NTAO
            DO 141 I=1,NT
               IF(FZT(I).EQ.1)GO TO 141
               QACT(MU,NU)=QACT(MU,NU)+WVEC(NU,I)*IMAT(MU,I)
  141       CONTINUE
  151    CONTINUE
  161 CONTINUE
C     WRITE(6,*)' BACK TRANSFORMED Q ACTIVE'
C     CALL MATOUT (QACT,NTAO,NTAO,NTAO,NTAO,6)
C
C >>> FORM AUXILIARY Q MATRIX
C
      CALL ZERO (QMAT,NTAO*NTAO)
      DO 260 MU=1,NTAO
         DO 250 NU=1,NTAO
            DO 240 I=1,NO
               II=CCTSCF(I)
COK            IF(FZT(II).EQ.1)GO TO 240
               QMAT(MU,NU)=QMAT(MU,NU)+WVEC(MU,II)*WVEC(NU,II)
  240       CONTINUE
  250    CONTINUE
  260 CONTINUE
C     WRITE(6,*)' AUXILIARY Q'
C     CALL MATOUT (QMAT,NTAO,NTAO,NTAO,NTAO,6)
C
C >>> FORM AUXILIARY QCOR FOR CORE ORBITALS
C
      CALL ZERO (QCOR,NTAO*NTAO)
      DO 261 MU=1,NTAO
         DO 251 NU=1,NTAO
            DO 241 I=1,NO
               II=CCTSCF(I)
               IF(FZT(II).NE.1)GO TO 241
               QCOR(MU,NU)=QCOR(MU,NU)+WVEC(MU,II)*WVEC(NU,II)
  241       CONTINUE
  251    CONTINUE
  261 CONTINUE
C     WRITE(6,*)' AUXILIARY Q CORE ONLY'
C     CALL MATOUT (QCOR,NTAO,NTAO,NTAO,NTAO,6)
C
C
C >>> FIRST HALF-TRANSFORMATION  FROM MO TO SO BASIS
C
      CALL RFILE(ITAP83)
      GPOINT=1
      CALL RFILE(ITAP95)
      IPT95=1
C
      DO 395 TOTSYM=1,NIRRED
      DO 390 I=1,NT
         IF(FZT(I).EQ.1)GO TO 390
         ISYM=ORBSYM(I)
         DO 380 J=1,I
            IF(FZT(J).EQ.1)GO TO 380
            JSYM=ORBSYM(J)
            IJSYM=IEOR(ISYM,JSYM)
            IF(IJSYM.NE.TOTSYM-1)GO TO 380
            IJ=IPQ(I)+J
            ICNT95=0
C
C >>> READ D(IJKL) SCF ORDERING, SYMM-PACKED TRIANGLE
C
            CALL ZERO(DMAT,NTAO*NTAO)
            IF(SRTOUT.EQ.'  NO')THEN
               DO 310 K=1,NT
                  IF(FZT(K).EQ.1)GO TO 310
                  KSYM=ORBSYM(K)
                  LSYM=IEOR(IJSYM,KSYM)+1
                  FL=FLOV(LSYM,1)
                  LL=FLOV(LSYM,2)
                  DO 300 L=FL,LL
                     KL=IPQ(MAX0(K,L))+MIN0(K,L)
                     IJKL=IOFF(MAX0(IJ,KL))+IADD(MIN0(IJ,KL))
                     DMAT(K,L)=GDM(IJKL)
  300             CONTINUE
  310          CONTINUE
            ELSE
               LRS=NRS(TOTSYM)
               CALL WREADW(ITAP83,AUX3,INTOWP(LRS),GPOINT,GPOINT)
               ICNT83=0
               DO 320 K=1,NT
                  IF(FZT(K).EQ.1)GO TO 320
                  KSYM=ORBSYM(K)
                  LSYM=IEOR(IJSYM,KSYM)+1
                  FL=FLOV(LSYM,1)
                  LL=FLOV(LSYM,2)
                  IF(LL.GT.K)LL=K
                  DO 318 L=FL,LL
                     ICNT83=ICNT83+1
                     DMAT(K,L)=AUX3(ICNT83)
                     DMAT(L,K)=AUX3(ICNT83)
  318             CONTINUE
  320          CONTINUE
            ENDIF
C
            CALL ZERO (IMAT,NTAO*NTAO)
COJO        CALL MXMB(VEC,1,NT,DMAT,NT,1,IMAT,1,NT,NT,NT,NT)
 
            DO 338 K=1,NT
               IF(FZT(K).EQ.1)GO TO 338
               KSYM=ORBSYM(K)
               LSYM=IEOR(IJSYM,KSYM)+1
               FL=FLOV(LSYM,1)
               LL=FLOV(LSYM,2)
               FS=FLOV(LSYM,1)-NFZO(LSYM)
               LS=FLOV(LSYM,2)+NFZV(LSYM)
               DO 336 S=FS,LS
                  DO 334 L=FL,LL
                     IMAT(S,K)=IMAT(S,K)+VEC(S,L)*DMAT(K,L)
C                    IMAT(S,K)=IMAT(S,K)+VEC(S,L)*DMAT(L,K)
  334             CONTINUE
  336          CONTINUE
  338       CONTINUE
C
            CALL ZERO (DMAT,NTAO*NTAO)
            DO 348 R=1,NT
               RSYM=ORBSYM(R)
               FK=FLOV(RSYM+1,1)
               LK=FLOV(RSYM+1,2)
               SSYM=IEOR(IJSYM,RSYM)+1
               FS=FLOV(SSYM,1)-NFZO(SSYM)
               LS=FLOV(SSYM,2)+NFZV(SSYM)
               IF(LS.GT.R)LS=R
               DO 346 S=FS,LS
                  ICNT95=ICNT95+1
                  DO 344 K=FK,LK
                     DMAT(R,S)=DMAT(R,S)+VEC(R,K)*IMAT(S,K)
  344             CONTINUE
CSORT             WRITE(6,'(4I6,F20.12)')I,J,R,S,DMAT(R,S)
                  AUX3(ICNT95)=DMAT(R,S)
  346          CONTINUE
  348       CONTINUE
            NCNTIN(TOTSYM)=ICNT95
CSORT       WRITE(6,*)'FILE95 WRITING PAIR OF LENGTH',I,J,ICNT95
C
C >>> SORT THE INTERMEDIATE MATRIX
C
            IF(SRTOUT.EQ.'  NO')THEN
               DO 358 R=1,NT
                  DO 356 S=1,R
                     VAL=DMAT(R,S)
                     CORE(R,S,I,J)=VAL
                     CORE(R,S,J,I)=VAL
  356             CONTINUE
  358          CONTINUE
            ELSE
               CALL WWRITW(ITAP95,AUX3,INTOWP(ICNT95),IPT95,IPT95)
            ENDIF
C
  380    CONTINUE
  390 CONTINUE
  395 CONTINUE
      CALL RCLOSE(ITAP95,3)
      CALL RCLOSE(ITAP83,4)
C
      IF(SRTOUT.EQ.' YES')THEN
         ITAPIN=ITAP95
         ITAPOU=ITAP96
         CALL SORT1(NT,ITAPIN,ITAPOU,NIRRED,FLOV,ORBSYM,NCNTIN,
     .              IPQ,IOFF,IADD,AUX3,LCOR,LEFT,FZT,NFZO,NFZV,NIJ)
      ENDIF
C
C >>> 2ND HALF-TRANS MO TO SO  &  1ST HALF-TRANS SO TO AO
C
      CALL RFILE(ITAP96)
      IPT96=1
      CALL RFILE(ITAP97)
      IPT97=1
      DO 495 TSYM=1,NIRRED
      DO 490 R=1,NT
         RSYM=ORBSYM(R)
         DO 480 S=1,R
            SSYM=ORBSYM(S)
            RSSYM=IEOR(RSYM,SSYM)
            IF(RSSYM.NE.TSYM-1)GO TO 480
C
C >>> READ RESORTED INTERMEDIATE
C
            CALL ZERO(DMAT,NTAO*NTAO)
            IF(SRTOUT.EQ.' YES')THEN
               LIJ=NIJ(TSYM)
CSORT          WRITE(6,*)'FILE96 READING PAIR OF LENGTH',R,S,LIJ
               CALL WREADW(ITAP96,AUX3,INTOWP(LIJ),IPT96,IPT96)
               ICNT96=0
            ENDIF
            DO 410 I=1,NT
               IF(FZT(I).EQ.1)GO TO 410
               ISYM=ORBSYM(I)
               JSYM=IEOR(RSSYM,ISYM)+1
               FJ=FLOV(JSYM,1)
               LJ=FLOV(JSYM,2)
               IF(LJ.GT.I)LJ=I
               DO 400 J=FJ,LJ
                  IF(FZT(J).EQ.1)GO TO 400
                  IF(SRTOUT.EQ.' YES')THEN
                     ICNT96=ICNT96+1
CSORT                WRITE(6,'(4I6,F20.12)')I,J,R,S,AUX3(ICNT96)
                     DMAT(I,J)=AUX3(ICNT96)
                     DMAT(J,I)=AUX3(ICNT96)
                  ELSE
                     DMAT(I,J)=CORE(R,S,I,J)
                     DMAT(J,I)=CORE(R,S,J,I)
                  ENDIF
  400          CONTINUE
  410       CONTINUE
C
            CALL ZERO (IMAT,NTAO*NTAO)
COJO        CALL MXMB(VEC,1,NT,DMAT,NT,1,IMAT,1,NT,NT,NT,NT)
C
            DO 438 Q=1,NT
               QSYM=ORBSYM(Q)
               FJ=FLOV(QSYM+1,1)-NFZO(QSYM+1)
               LJ=FLOV(QSYM+1,2)+NFZV(QSYM+1)
               ISYM=IEOR(RSSYM,QSYM)+1
               FI=FLOV(ISYM,1)-NFZO(ISYM)
               LI=FLOV(ISYM,2)+NFZV(ISYM)
               DO 436 I=FI,LI
                  DO 434 J=FJ,LJ
                     IMAT(Q,I)=IMAT(Q,I)+VEC(Q,J)*DMAT(I,J)
  434             CONTINUE
  436          CONTINUE
  438       CONTINUE
C
            CALL ZERO (DMAT,NTAO*NTAO)
            DO 448 P=1,NT
               PSYM=ORBSYM(P)
               FI=FLOV(PSYM+1,1)-NFZO(PSYM+1)
               LI=FLOV(PSYM+1,2)+NFZV(PSYM+1)
               QSYM=IEOR(RSSYM,PSYM)+1
               FQ=FLOV(QSYM,1)-NFZO(QSYM)
               LQ=FLOV(QSYM,2)+NFZV(QSYM)
               IF(LQ.GT.P)LQ=P
               DO 446 Q=FQ,LQ
                  DO 444 I=FI,LI
                     DMAT(P,Q)=DMAT(P,Q)+VEC(P,I)*IMAT(Q,I)
  444             CONTINUE
                     DMAT(Q,P)=DMAT(P,Q)
  446          CONTINUE
  448       CONTINUE
C
            CALL ZERO(IMAT,NTAO*NTAO)
c           if(ic1.eq.0)then
            DO 458 MU=1,NTAO
               DO 456 P=1,NT
                  VAL=UVEC(MU,P)
                  IF(VAL.EQ.Z0)GO TO 456
                  DO 454 Q=1,NT
                     IMAT(MU,Q)=IMAT(MU,Q)+VAL*DMAT(P,Q)
  454             CONTINUE
  456          CONTINUE
  458       CONTINUE
C
            CALL ZERO(DMAT,NTAO*NTAO)
            DO 468 NU=1,NTAO
               DO 466 Q=1,NT
                  VAL=UVEC(NU,Q)
                  IF(VAL.EQ.Z0)GO TO 466
                  DO 464 MU=1,NTAO
                     DMAT(MU,NU)=DMAT(MU,NU)+VAL*IMAT(MU,Q)
  464             CONTINUE
  466          CONTINUE
  468       CONTINUE
c           endif
C
C >>>  SORT THE INTERMEDIATE MATRIX
C
            ICNT97=0
            DO 478 MU=1,NTAO
               DO 474 NU=1,MU
                  IF(SRTOUT.EQ.'  NO')THEN
                     CORE2(MU,NU,R,S)=DMAT(MU,NU)
                  ELSE
                     ICNT97=ICNT97+1
                     AUX3(ICNT97)=DMAT(MU,NU)
                  ENDIF
  474          CONTINUE
  478       CONTINUE
C
            IF(SRTOUT.EQ.' YES')THEN
               CALL WWRITW(ITAP97,AUX3,INTOWP(NTAO2),IPT97,IPT97)
            ENDIF
C
  480    CONTINUE
  490 CONTINUE
  495 CONTINUE
C
      IF(SRTOUT.EQ.' YES')THEN
         CALL RCLOSE(ITAP96,4)
         CALL RCLOSE(ITAP97,3)
         ITAPIN=ITAP97
         ITAPOU=ITAP98
         CALL SORT2(NT,NTAO,ITAPIN,ITAPOU,IPQ,IOFF,IADD,NT2,NTAO2,
     .              AUX3,LCOR,LEFT,KLOC,AOPS,NSHELL,ORBSYM,NIRRED)
      ENDIF
C
C >>> 2ND HALF-TRANSFORMATION SO TO AO BASIS
C
      CALL RFILE(ITAP98)
      IPT98=1
C
      DO 695 MUSHEL=1,NSHELL
         FMU=KLOC(MUSHEL)
         LMU=FMU+AOPS(MUSHEL)-1
         DO 685 NUSHEL=1,MUSHEL
            FNU=KLOC(NUSHEL)
            LNU=FNU+AOPS(NUSHEL)-1
      DO 690 MU=FMU,LMU
         MUR=MU-FMU+1
         IF(NUSHEL.EQ.MUSHEL)LNU=MU
         DO 680 NU=FNU,LNU
            NUR=NU-FNU+1
            ACNT1=(MUR-1)*AOPS(NUSHEL)+NUR
            ACNT2=(NUR-1)*AOPS(MUSHEL)+MUR
C
C >>> READ RESORTED INTERMEDIATE
C
            IF(SRTOUT.EQ.' YES')THEN
               CALL WREADW(ITAP98,AUX3,INTOWP(NT2),IPT98,IPT98)
            ENDIF
            CALL ZERO(DMAT,NTAO*NTAO)
            ICNT98=0
            DO 610 R=1,NT
               DO 600 S=1,R
                  IF(SRTOUT.EQ.'  NO')THEN
                     DMAT(R,S)=CORE2(MU,NU,R,S)
                     DMAT(S,R)=CORE2(MU,NU,R,S)
                  ELSE
                     ICNT98=ICNT98+1
                     DMAT(R,S)=AUX3(ICNT98)
                     DMAT(S,R)=AUX3(ICNT98)
                  ENDIF
  600          CONTINUE
  610       CONTINUE
C
            CALL ZERO (IMAT,NTAO*NTAO)
COJO        CALL MXMB(VEC,1,NT,DMAT,NT,1,IMAT,1,NT,NT,NT,NT)
C
c           if (ic1.eq.0)then
            DO 638 SI=1,NTAO
               DO 636 S=1,NT
                  VAL=UVEC(SI,S)
                  IF(VAL.EQ.Z0)GO TO 636
                  DO 634 R=1,NT
                     IMAT(R,SI)=IMAT(R,SI)+VAL*DMAT(R,S)
  634             CONTINUE
  636          CONTINUE
  638       CONTINUE
C
            CALL ZERO (DMAT,NTAO*NTAO)
CBUG!       DO 648 LA=1,MU
            DO 648 LA=1,LMU
               DO 646 R=1,NT
                  VAL=UVEC(LA,R)
                  IF(VAL.EQ.Z0)GO TO 646
                  LSI=LA
                  IF(LA.EQ.LMU)LSI=NU
                  DO 644 SI=1,LSI
                     DMAT(LA,SI)=DMAT(LA,SI)+VAL*IMAT(R,SI)
  644             CONTINUE
  646          CONTINUE
  648       CONTINUE
c           endif
C
C >>> ADD HF AND Z-CPHF CONTRIBUTION TO 2-PDM IN AO BASIS
C
            DO 658 LA=1,LMU
               LSI=LA
               IF(LA.EQ.LMU)LSI=NU
               DO 654 SI=1,LSI
CWRITE            WRITE(6,6542)MU,NU,LA,SI,DMAT(LA,SI)
                  VHF = 2.0 D0 * QMAT(MU,NU)*QMAT(LA,SI)
     .                - 0.5 D0 *(QMAT(MU,SI)*QMAT(NU,LA)
     .                         + QMAT(MU,LA)*QMAT(NU,SI))
                  VAL =          BZET(MU,NU)*QMAT(LA,SI)
     .                +          BZET(LA,SI)*QMAT(MU,NU)
     .                - 0.25D0 *(BZET(MU,SI)*QMAT(NU,LA)
     .                         + BZET(MU,LA)*QMAT(NU,SI)
     .                         + BZET(NU,LA)*QMAT(MU,SI)
     .                         + BZET(NU,SI)*QMAT(MU,LA))
                  VCR =          QACT(MU,NU)*QCOR(LA,SI)
     .                +          QACT(LA,SI)*QCOR(MU,NU)
     .                - 0.25D0 *(QACT(MU,SI)*QCOR(NU,LA)
     .                         + QACT(MU,LA)*QCOR(NU,SI)
     .                         + QACT(NU,LA)*QCOR(MU,SI)
     .                         + QACT(NU,SI)*QCOR(MU,LA))
C
C                 VCR2= 2.0 D0 * QCOR(MU,NU)*QCOR(LA,SI)
C    .                - 0.5 D0 *(QCOR(MU,SI)*QCOR(NU,LA)
C    .                         + QCOR(MU,LA)*QCOR(NU,SI))
C
C                WRITE(6,6542) MU,NU,LA,SI,DMAT(LA,SI)
C                WRITE(6,6214) VHF,VAL,VCR,VCR2
CCHECK            DMAT(LA,SI)=DMAT(LA,SI)+VHF
                  DMAT(LA,SI)=DMAT(LA,SI)+VHF+VAL+VCR
CCHECK            DMAT(LA,SI)=DMAT(LA,SI)+VHF+VAL+VCR+VCR2
  654          CONTINUE
  658       CONTINUE
C
      MUNUX=IPQ(MUSHEL)+NUSHEL
      MUNU=IJOFF(MUNUX)
      DO 669 LASHEL=1,MUSHEL
         FLA=KLOC(LASHEL)
         LLA=FLA+AOPS(LASHEL)-1
         LSISHL=LASHEL
         IF(MUSHEL.EQ.LASHEL)LSISHL=NUSHEL
         DO 668 SISHEL=1,LSISHL
            FSI=KLOC(SISHEL)
            LSI=FSI+AOPS(SISHEL)-1
            LASIX=IPQ(LASHEL)+SISHEL
            LASI=KLOFF(LASIX)*LNSR(MUNUX)
            SADR=INTOWP(MUNU+LASI)+AOSTRT
            NBLA=AOPS(LASHEL)
            NBSI=AOPS(SISHEL)
            SIZE1=NBLA*NBSI
            IADR1=SADR+INTOWP((ACNT1-1)*SIZE1)
                  ICNT=0
                  DO 665 LA=FLA,LLA
                     DO 664 SI=FSI,LSI
                        ICNT=ICNT+1
                        LAM=MAX0(LA,SI)
                        SIM=MIN0(LA,SI)
                        BUF(ICNT)=DMAT(LAM,SIM)
  664                CONTINUE
  665             CONTINUE
CWRITE         WRITE(6,957)MUSHEL,NUSHEL,LASHEL,SISHEL,MU,NU,IADR1,ICNT
CWRITE         WRITE(6,*)(BUF(INDEX),INDEX=1,ICNT)
                  CALL WWRITW(ITAP55,BUF,INTOWP(ICNT),IADR1,IDUM)
C
                  IF(MUSHEL.EQ.NUSHEL.AND.MU.NE.NU)THEN
                     IADR2=SADR+INTOWP((ACNT2-1)*SIZE1)
CWRITE         WRITE(6,958)MUSHEL,NUSHEL,LASHEL,SISHEL,NU,MU,IADR2,ICNT
CWRITE         WRITE(6,*)(BUF(INDEX),INDEX=1,ICNT)
                  CALL WWRITW(ITAP55,BUF,INTOWP(ICNT),IADR2,IDUM)
                  ENDIF
 
  668    CONTINUE
  669 CONTINUE
C
  680    CONTINUE
  690 CONTINUE
C
  685    CONTINUE
  695 CONTINUE
C
      IF(SRTOUT.EQ.' YES')THEN
      CALL RCLOSE(ITAP98,4)
      ENDIF
C
  957 FORMAT(4I4,5X,2I4,5X,'IADR1=',I6,'  ICNT=',I6)
  958 FORMAT(4I4,5X,2I4,5X,'IADR2=',I6,'  ICNT=',I6)
 6214 FORMAT(4F18.6)
 6542 FORMAT(4I6,F20.12)
C
      RETURN
      END
