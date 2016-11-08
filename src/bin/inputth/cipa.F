      SUBROUTINE CIPA(PTR,DTR,KATOM,KTYPE,ISC,IPC,CC,LAB,LOC,LOC2,
     1  KLASS,LAMBDA,IPA,KMIN,SLAB,IADDR,ISO,ICT,GAM,TRACEP,TRACED,
     2  CHARR,ISGRP,FTR,GTR,TRACEF,TRACEG)
C
C  GENERATES COEFFICIENTS OF SYMMETRY-ADAPTED LINEAR COMBINATIONS
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      LOGICAL SKIP(5)
C
      COMMON /MAXDIM/ MXAT,MXTRN,MXPR,MXSH,MXSYM,MXCF,MXT,MXISO,MXPSYM
      COMMON /SCHMDT/ P(18,18),PNORM(18)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV
      COMMON /SYMTRY/ NT
      COMMON /SMORB/  GRP,IGROUP,NAXIS
      COMMON /NSHEL/  NSHELL,MSHELL
      COMMON /INFOA/  NAT,NUM,NQW,NX
      COMMON /INOUT/  IOUT,ICHECK,ITAP30
      COMMON /ISOPAC/ INDIN(120),INDOUT(30)
      COMMON /SALCS/  NIRR
      COMMON /CIPAS/  NIRRED
      COMMON /OUTPUT/ NPRINT,NGEOM,NORMF,NORMP
      COMMON /PACKQ/  I1,I2,I3,I4,I5
C
      REAL*8 PTR(3,1),DTR(6,1),FTR(10,1),GTR(15,1)
      REAL*8 CC(MXCF),GAM(MXTRN,MXTRN),CHARR(MXTRN,MXSYM)
      REAL*8 TRACEP(MXTRN),TRACED(MXTRN),TRACEF(MXTRN),TRACEG(MXTRN)
      INTEGER KATOM(MXSH),KTYPE(MXSH),ISC(MXSH,MXPSYM),LAB(MXCF)
      INTEGER IPC(MXSH,MXPSYM),LOC(MXAT,5),LOC2(MXAT,5),KLASS(MXAT)
      INTEGER LAMBDA(MXSYM),IPA(MXSYM),KMIN(MXSH),IADDR(MXSYM)
      INTEGER ISO(MXSH,MXISO),ICT(MXAT,MXTRN),SLAB(MXSYM),ISGRP(MXAT)
      DIMENSION KBB(5),KCC(5)
      DIMENSION CA(18),CHAR(5,48)
      DIMENSION IPB(70,14),ISB(70,14)
C
      DATA TOL/1.0D-10/
      DATA TOLTOL/1.0D-20/
c     DATA TOLTOL/1.0D-15/
      DATA ZERO,ONE/0.0D+00,1.0D+00/
      DATA BLOCKI/8HBLOCKI  /
cges
c
c      write(6,*)' in cipa'
c      write(6,*)'mxcf=',mxcf
C
C ZERO OUT THE P ARRAY
C
      DO 21 I=1,18
        DO 20 J=1,I
          P(I,J)=ZERO
          P(J,I)=ZERO
   20   CONTINUE
   21   CONTINUE
C
      NIRR=NIRRED
      DO 100 N=1,NSHELL
         DO 100 I=1,NIRRED
            ISB(N,I)=0
  100 IPB(N,I)=0
      DO 496 J=1,5
         DO 496 I=1,5
  496 LOC(J,I)=-1
C
C     FOR EACH CENTER FIND SUBGROUP WHICH MAPS IT INTO ITSELF
C
C***  IF(NT.GT.48) CALL ABORT(9,
C*** 168HSUBROUTINE CIPA--ORDER OF POINT GROUP TOO LARGE FOR SUBGROUP PA
C*** 2CKING)
      DO 498 N=1,NAT
         K=1
         J=0
         DO 497 IT=1,NT
            IF(ICT(N,IT).NE.N) GO TO 497
            J=J+K
  497    K=K+K
         ISGRP(N)=J
C     WRITE(IW,9987) N,J
C9987    FORMAT(I5,5X,Z8)
  498 CONTINUE
      DO 499 I=1,NAT
  499 KLASS(I)=0
      ICLASS=0
      DO 501 I=1,NIRRED
  501 IPA(I)=0
C
      IC=0
      DO 1000 NATOM=1,NAT
C
C  PICK CENTER IN P1 LIST.
C
         DO 500 IT=1,NT
           IF(ICT(NATOM,IT).GT.NATOM) GO TO 1000
  500    CONTINUE
         IF(KLASS(NATOM).NE.0) GO TO 1000
         ICLASS=ICLASS+1
         KLASS(NATOM)=ICLASS
         ISG=ISGRP(NATOM)
C
C  FIND OTHER CENTERS IN SAME CLASS
C
         MIN=NATOM+1
         IF(MIN.GT.NAT) GO TO 50
         DO 40 J=MIN,NAT
            JSG=ISGRP(J)
            IF(JSG.NE.ISG) GO TO 40
            DO 30 IT=1,NT
               IF(ICT(J,IT).GT.J) GO TO 40
   30       CONTINUE
            KLASS(J)=ICLASS
   40    CONTINUE
   50    CONTINUE
         DO 60 I=1,5
   60    SKIP(I)=.TRUE.
CJ
         DO 80 I=1,NSHELL
            J=KATOM(I)
            IF(KLASS(J).NE.ICLASS) GO TO 80
            J=KTYPE(I)
            GO TO (70,71,72,73,74),J
   70       SKIP(1)=.FALSE.
            GO TO 80
   71       SKIP(2)=.FALSE.
            IF(KMIN(I).EQ.1) SKIP(1)=.FALSE.
            GO TO 80
   72       SKIP(3)=.FALSE.
            GO TO 80
   73       SKIP(4)=.FALSE.
            GO TO 80
   74       SKIP(5)=.FALSE.
            GO TO 80
   80    CONTINUE
         DO 510 IT=1,NT
            DO 510 I=1,5
  510    CHAR(I,IT)=ZERO
C
C  FIND CHARACTERS OF REDUCIBLE REPRESENTATION.
C
         DO 550 IATOM=1,NAT
            DO 515 IT=1,NT
              IF(ICT(NATOM,IT).EQ.IATOM) GO TO 516
  515       CONTINUE
            GO TO 550
  516       CONTINUE
            DO 520 IT=1,NT
              IF(ICT(IATOM,IT).NE.IATOM) GO TO 520
C
C  EQUIVALENT CENTER IATOM WILL CONTRIBUTE ONLY TO CHARACTER
C  OF OPERATIONS WHICH MAP IT INTO ITSELF.
C
              CHAR(1,IT)=CHAR(1,IT)+ONE
              CHAR(2,IT)=CHAR(2,IT)+TRACEP(IT)
              CHAR(3,IT)=CHAR(3,IT)+TRACED(IT)
              CHAR(4,IT)=CHAR(4,IT)+TRACEF(IT)
              CHAR(5,IT)=CHAR(5,IT)+TRACEG(IT)
  520       CONTINUE
  550    CONTINUE
C     WRITE(IW,9989) ((CHARR(IT,IRP),IT=1,NT),IRP=1,NIRRED)
C9989    FORMAT(//(1X,12F10.5))
C     WRITE(IW,9989) ((CHAR(KT,IT),IT=1,NT),KT=1,3)
C
C  DECOMPOSE REDUCIBLE REPRESENTATION INTO IRREDUCIBLE
C  REPRESENTATIONS.
C
         G=DBLE(NT)
         NP=0
         NF=0
         DO 800 KT=1,5
            NF=NF+NP
            NP=NP+KT
            IF(SKIP(KT)) GO TO 800
  562       IB=0
            LOC(ICLASS,KT)=IC
            DO 700 IRP=1,NIRRED
               LL=LAMBDA(IRP)
               FAC=DBLE(LL)/G
               FNA=ZERO
               DO 570 IT=1,NT
  570          FNA=FNA+CHAR(KT,IT)*CHARR(IT,IRP)
               FNA=FNA/G
               IF(FNA.LT.TOL) GO TO 700
C
C  IRREDUCIBLE REPRESENTATION IRP OCCURS FNA TIMES.
C  FORM P-MATRIX.
C
               MAX=NP*LL
               DO 585 I=1,MAX
                  DO 585 J=1,MAX
  585          P(I,J)=ZERO
               ITB=0
               DO 620 IT=1,NT
                  IF(ICT(NATOM,IT).NE.NATOM) GO TO 620
                  LB1=0
                  DO 610 L1=1,LL
                     LB2=0
                     DO 600 L2=1,LL
                        GM=GAM(IT,IB+LL*(L1-1)+L2)*FAC
c                 write(6,*)'gm=',gm
                        DO 590 I=1,NP
                           DO 590 J=1,NP
                             DUM=GM
                             IF(KT.EQ.2) DUM=DUM*PTR(I,ITB+J)
                             IF(KT.EQ.3) DUM=DUM*DTR(I,ITB+J)
                             IF(KT.EQ.4) DUM=DUM*FTR(I,ITB+J)
                             IF(KT.EQ.5) DUM=DUM*GTR(I,ITB+J)
c                            write(6,*)'dum,kt',dum ,kt
                             P(LB1+I,LB2+J)=P(LB1+I,LB2+J)+DUM
  590                   CONTINUE
  600                LB2=LB2+NP
  610             LB1=LB1+NP
  620          ITB=ITB+NP
C
C  ORTHOGONALIZE COLUMNS OF P-MATRIX.
C
               MAX=NP*LL
               DO 2001 IXC=1,MAX
                 DO 2000 IXV=1,MAX
                   IF(DABS(P(IXC,IXV)) .LE. TOLTOL) P(IXC,IXV)=ZERO
 2000            CONTINUE
 2001          CONTINUE
 2002          FORMAT(D15.10)
               CALL SHMIDT(MAX)
               DO 630 J=1,MAX
                  IF(PNORM(J).LT.TOL) GO TO 630
c     WRITE(IW,9990) NATOM,KT,IRP,FNA,(P(I,J),I=1,MAX)
 9990             FORMAT(3I5,F10.5,6F12.6/25X,6F12.6)
  630          CONTINUE
               IPB1=0
               DO 660 J=1,MAX
                  IF(PNORM(J).LT.TOL) GO TO 660
                  IPA(IRP)=IPA(IRP)+1
                  IPB1=IPB1+1
                  LB1=0
                  DO 650 L=1,LL
                     DO 640 I=1,NP
                        IC=IC+1
C***                    IF(IC.GT.MXCF) CALL ABORT(9,68HSUBROUTINE CIPA--
C*** *NO. OF SALC COEFFICIENTS EXCEEDS CURRENT ARRAY SIZE)
                        I1=I+NF
                        I2=IPA(IRP)
                        I3=IRP
                        I4=L
                        CALL PAK
                        LAB(IC)=I5
                        CC(IC)=P(LB1+I,J)
c                       write(6,*)'ic,cc(ic)',ic,cc(ic)
  640                CONTINUE
  650             LB1=LB1+NP
  660          CONTINUE
               IPB2=-IPB1
               DO 680 I=1,NSHELL
                  IF(KLASS(KATOM(I)).NE.KLASS(NATOM)) GO TO 680
                  KKT=KTYPE(I)
                  K=0
                  IF(KKT.EQ.KT) GO TO 663
                  IF(KKT.EQ.2.AND.KT.EQ.1.AND.KMIN(I).EQ.1) GO TO 661
                  GO TO 680
  661             K=1
  663             IPB2=IPB2+IPB1
                  DO 665 J=1,15
  665             INDOUT(J)=ISO(I,J)
                  CALL ISOOUT
                  IF(K.EQ.1) GO TO 675
                  IPB(I,IRP)=IPB2
                  DO 670 IT=1,NT
                     J=INDIN(IT)
  670             IPB(J,IRP)=IPB2
                  GO TO 680
  675             ISB(I,IRP)=IPB2
                  DO 677 IT=1,NT
                     J=INDIN(IT)
  677             ISB(J,IRP)=IPB2
  680          CONTINUE
               IPA(IRP)=IPA(IRP)+IPB2
  700       IB=IB+LL*LL
            LOC2(ICLASS,KT)=IC
  800    CONTINUE
         DO 820 KT=1,5
            KBB(KT)=LOC(ICLASS,KT)
  820    KCC(KT)=LOC2(ICLASS,KT)
         ITM=NT
         JCLASS=ICLASS
C
C     ROTATE COEFFICIENTS ONTO EQUIVALENT CENTERS
C
         DO 1450 IT=1,ITM
C
C     FIND EQUIVALENT CENTER
C
            ICEN=ICT(NATOM,IT)
            IF(KLASS(ICEN).NE.0) GO TO 1450
            ICLASS=ICLASS+1
            KLASS(ICEN)=ICLASS
C
C     FIND OTHER CENTERS IN SAME CLASS
C
            DO 840 I=1,NAT
               IF(KLASS(I).NE.JCLASS) GO TO 840
               J=ICT(I,IT)
               KLASS(J)=ICLASS
  840       CONTINUE
            NP=0
            DO 1400 KT=1,5
               NP=NP+KT
               IF(SKIP(KT)) GO TO 1400
               LOC(ICLASS,KT)=IC
               KB=KBB(KT)
               KC=KCC(KT)
 1165          I5=LAB(KB+1)
               CALL UNPAK
               IRP=I3
               L=I4
               IB=0
               IF(IRP.EQ.1) GO TO 1180
               MAX=IRP-1
               DO 1170 I=1,MAX
                  LL=LAMBDA(I)
 1170          IB=IB+LL*LL
 1180          LL=LAMBDA(IRP)
               MAX=NP*LL
               DO 1190 I=1,MAX
 1190          CA(I)=CC(KB+I)
               JC=0
               DO 1240 K=1,LL
                  DO 1240 I=1,NP
                     IC=IC+1
C***                 IF(IC.GT.MXCF) CALL ABORT(9,68HSUBROUTINE CIPA--NO.
C*** * OF SALC COEFFICIENTS EXCEEDS CURRENT ARRAY SIZE)
                     DUM=ZERO
                     D=ONE
                     DO 1200 J=1,NP
                        DO 1200 M=1,LL
                          IF(KT.EQ.2) D=PTR(I,J+3*(IT-1))
                          IF(KT.EQ.3) D=DTR(I,J+6*(IT-1))
                          IF(KT.EQ.4) D=FTR(I,J+10*(IT-1))
                          IF(KT.EQ.5) D=GTR(I,J+15*(IT-1))
                          DUM=DUM+D*CA(J+NP*(M-1))*GAM(IT,IB+LL*(K-1)+M)
 1200                CONTINUE
                     CC(IC)=DUM
c                  write(6,*)'at 1200 ic,cc(ic)',ic,cc(ic)
                     JC=JC+1
                     LAB(IC)=LAB(KB+JC)
 1240          CONTINUE
               KB=KB+MAX
               IF(KB.LT.KC) GO TO 1165
               LOC2(ICLASS,KT)=IC
 1400       CONTINUE
 1450    CONTINUE
 1000 CONTINUE
C     DO 1100 I=1,IC
C1100 WRITE(IW,9988) I,LAB(I),CC(I)
 9988 FORMAT(I5,5X,Z8,1PE25.14)
C
C     DETERMINE INITIAL ADDRESS OF SYMMETRY BLOCKS
C
      J=0
      K=0
      DO 1300 I=1,NIRRED
         N=IPA(I)
         IF(N.EQ.0) GO TO 1300
         IADDR(I)=J
         J=J+N*(N+1)/2
         K=K+N*N
 1300 CONTINUE
      IF (NPRINT.NE.1) GO TO 1355
      DO 1350 I=1,NSHELL
         WRITE(IW,9986) (ISB(I,J),J=1,NIRRED)
 1350 WRITE(IW,9986) (IPB(I,J),J=1,NIRRED)
      WRITE(IW,9986) (IADDR(J),J=1,NIRRED)
      WRITE(IW,9986) (IPA(J),J=1,NIRRED),K
 1355 CONTINUE
 9986 FORMAT(14I5)
      K=NT
      NT=NIRRED
      DO 1500 I=1,NSHELL
         DO 1480 J=1,NIRRED
 1480    INDIN(J)=IPB(I,J)+1
         CALL ISOIN
         IPC(I,1)=INDOUT(1)
         IPC(I,2)=INDOUT(2)
         DO 1490 J=1,NIRRED
 1490    INDIN(J)=ISB(I,J)+1
         CALL ISOIN
         ISC(I,1)=INDOUT(1)
         ISC(I,2)=INDOUT(2)
 1500 CONTINUE
      NT=K
      MXCF=IC
      IF (NPRINT.EQ.1) CALL SALC(KLASS,KTYPE,KMIN,KATOM,LOC,LOC2,
     #LAB,CC,LAMBDA,SLAB,ISC,IPC,ICT)
C
C
C
      CALL COUNT(KLASS,KTYPE,KMIN,KATOM,LOC,LOC2,LAB,CC,LAMBDA,
     #ISC,IPC,ICT)
cges  IF(BLKTYP.EQ.BLOCKI) CALL DIJGEN(LOC,ICT,KLASS,NAT,NT)
 9000 CONTINUE
      RETURN
      END
