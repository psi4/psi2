      SUBROUTINE INTGRL(BFSYM,ORBTBF,KADD,LADD,IJADD,IJGRP,ININT
     #,                 INEXT,JMNNXT,JMXNXT,NINGRP,IJXX,KLXX,NKLXX
     #,                 IJWW,KLWW,NKLWW)
C
C***********************************************************************
C                                                                      *
C PWS 19 JULY 1982                                                     *
C     MODIFIED 19 JULY 1982 TO FORM 3 AND 4-EXTERNAL ADDRESSING        *
C     SCHEME AND TO GENERATE IJ-BLOCKS RATHER THAN ONLY I-BLOCKS       *
C                                                                      *
C     COUNT AND COMPUTE ADDRESSES FOR THE INTEGRALS INCLUDING SYMMETRY.*
C     THE ADDRESS OF AN INTEGRAL WITH INDICES I>J>K>L (OR EQUAL) IS    *
C     GIVEN BY THE FOLLOWING EXPRESSION:                               *
C                                                                      *
C IJADD(IJ)+KADD(K+SYM(I)SYM(L)*NORBS)+LADD(L+SYM(I)SYM(L)SYM(K)*NORBS)*
C                                                                      *
C     WHERE IJ=I*(I-1)/2+J AND NORBS IS NUMBER OF ORBITALS IN CI. NOTE *
C     ALSO THAT SYMMETRIES NEED DIRECT PRODUCTS. THE INTEGRALS ARE     *
C     STORED WITH THE FOLLOWING OFFSET FROM THE ADDRESS ABOVE:         *
C                                                                      *
C       TYPE            1            2          3                      *
C         1          (IK,JL)      (IJ,KL)    (IL,JK)                   *
C         2          (IJ,JL)      (IL,JJ)                              *
C         3          (IK,IL)      (II,KL)                              *
C         4          (IL,JL)      (IJ,LL)                              *
C         5          (II,IL)      (IL,LL)    <I/H/L>                   *
C         6          (IL,IL)      (II,LL)                              *
C         7          (II,II)      <I/H/I>                              *
C                                                                      *
C     THE ADDRESS OF 3- AND 4-EXTERNAL MATRIX ELEMENTS IS GIVEN BY     *
C                                                                      *
C             WY ENTRY                   XY ENTRY                      *
C        IJWW(IA)+KLWW(BC)    OR    IJXX(IA)+KLXX(BC)       (3X)       *
C                             OR                                       *
C        IJWW(AB)+KLWW(CD)    OR    IJXX(AB)+KLXX(CD)       (4X)       *
C             WW ENTRY                   XX ENTRY                      *
C                                                                      *
C     THIS SUBROUTINE ALSO DETERMINES THE INTEGRAL BLOCK SIZE (NMAX).  *
C     IJGRP(IJ) GIVES BLOCK CONTAINING ALL INTEGRALS WITH I AND J      *
C     INDICES.                                                         *
C***********************************************************************
C
C
      IMPLICIT INTEGER (A-Z)
      INTEGER NUMINT
C
      CHARACTER*4 LABEL
      COMMON /OPTION/ OPTION(16),IDENT(26),LABEL(26)
      COMMON /BLOKSZ/ BLKSIZ,ABSMAX,MAXSIZ
      COMMON /TAPES/  OUT,ERROUT,INPUT,DRTTAP
      COMMON /DIMENS/ NBF,NSYM,NORBS,NROWSP,NROW4P,NROWS,NROWS4
     #,               NLEVS,NREFS,NROWOC,NROW4O,NWKS,NWKSOC,NLEVOC
     #,               ORBFRM,SYMORB,NUMIJ,NGROUP,NUMINT,NMAX,NSPC,NVREF
     #,               NIJVIR
C
      DIMENSION BFSYM(NBF),ORBTBF(NORBS),KADD(SYMORB),LADD(SYMORB)
      DIMENSION IJADD(NUMIJ),IJGRP(NUMIJ),ININT(NORBS)
      DIMENSION INEXT(NORBS),JMNNXT(NORBS),JMXNXT(NORBS),NINGRP(1)
      DIMENSION IJXX(NUMIJ),KLXX(NUMIJ),IJWW(NUMIJ)
      DIMENSION KLWW(NUMIJ),NKLWW(NSYM,ORBFRM),NKLXX(NSYM,ORBFRM)
C
      DO 101 I=1,NUMIJ
         IJXX(I)=0
         IJWW(I)=0
  101 CONTINUE
C
C     ----- COUNT INTEGRALS AND FORM KADD AND LADD ARRAYS -----
C
      IF (OPTION(5).GT.0) BLKSIZ=OPTION(5)*100
      IF (OPTION(5).GT.0) ABSMAX=BLKSIZ
      NUMINT=0
      JADMAX=0
      IJAD=0
      DO 9 I=1,NORBS
         ISYM=BFSYM(ORBTBF(I))
         JAD=0
         IJ=I*(I-1)/2
         DO 8 J=1,I
            IJ=IJ+1
            IJSYM=IEOR(ISYM,BFSYM(ORBTBF(J)))
            KAD=0
            DO 7 K=1,J
               IJKAD=K+IJSYM*NORBS
               IF (KADD(IJKAD).EQ.0.OR.KADD(IJKAD).EQ.KAD) GO TO 2
               WRITE (OUT,1)
               WRITE (ERROUT,1)
    1          FORMAT (//,' SYMMETRY PROBLEMS WITH K IN INTGRL',//)
               CALL mabort
    2          CONTINUE
               KADD(IJKAD)=KAD
               IJKSYM=IEOR(IJSYM,BFSYM(ORBTBF(K)))
               LAD=0
               DO 6 L=1,K
                  IJKLAD=L+IJKSYM*NORBS
                  IF (LADD(IJKLAD).EQ.0.OR.LADD(IJKLAD).EQ.LAD) GO TO 4
                  WRITE (OUT,3)
                  WRITE (ERROUT,3)
    3             FORMAT (//,' SYMMETRY PROBLEMS WITH L IN INTGRL',//)
                  CALL mabort
    4             CONTINUE
                  LADD(IJKLAD)=LAD
                  IF (IJKSYM.NE.BFSYM(ORBTBF(L))) GO TO 5
                  IF (L.EQ.J.AND.L.LT.I) GO TO 5
                  LAD=LAD+3
C     IF (K.EQ.L) LAD=LAD-1
    5             CONTINUE
    6          CONTINUE
               KAD=KAD+LAD
    7       CONTINUE
            JAD=JAD+KAD
            IJADD(IJ)=KAD
    8    CONTINUE
         NUMINT=NUMINT+JAD
         ININT(I)=JAD
         IF (JAD.GT.JADMAX) JADMAX=JAD
    9 CONTINUE
C
C     ----- WORK OUT KL ADDRESSING FOR 3- AND 4- EXTERNAL INTEGRALS ----
C
      DO 34 ISYM=1,NSYM
         NUM=0
         DO 33 I=2,ORBFRM
            IS=IEOR(BFSYM(ORBTBF(I)),(ISYM-1))
            IA=I*(I-1)/2
            DO 32 J=1,I-1
               IF (BFSYM(ORBTBF(J)).NE.IS) GO TO 31
               NUM=NUM+1
               KLXX(IA+J)=NUM
   31          CONTINUE
   32       CONTINUE
            NKLXX(ISYM,I)=NUM
   33    CONTINUE
C     NKLXX(ISYM)=NUM
   34 CONTINUE
C
      DO 38 ISYM=1,NSYM
         NUM=0
         DO 37 I=1,ORBFRM
            IS=IEOR((ISYM-1),BFSYM(ORBTBF(I)))
            IA=I*(I-1)/2
            DO 36 J=1,I
               IF (BFSYM(ORBTBF(J)).NE.IS) GO TO 35
               NUM=NUM+1
               KLWW(IA+J)=NUM
   35          CONTINUE
   36       CONTINUE
            NKLWW(ISYM,I)=NUM
   37    CONTINUE
C     NKLWW(ISYM)=NUM
   38 CONTINUE
C
C     ----- NUMBER OF 3- AND 4-EXTERNAL ELEMENTS IN IJ BLOCKS -----
C
      DO 44 JUNK=ORBFRM+1,NORBS
         I=NORBS-JUNK+ORBFRM+1
         IA=I*(I-1)/2
         IS=BFSYM(ORBTBF(I))
         NUM=0
         DO 42 J=1,ORBFRM
            IJS=IEOR(IS,BFSYM(ORBTBF(J)))
            IJWW(IA+J)=NKLWW(IJS+1,ORBFRM)
            NUM=NUM+NKLWW(IJS+1,ORBFRM)
   42    CONTINUE
         DO 43 J=1,ORBFRM
            IJS=IEOR(IS,BFSYM(ORBTBF(J)))
            IJXX(IA+J)=NKLXX(IJS+1,ORBFRM)
            NUM=NUM+NKLXX(IJS+1,ORBFRM)
   43    CONTINUE
         ININT(I)=ININT(I)+NUM
         NUMINT=NUMINT+NUM
   44 CONTINUE
C
      DO 41 JUNK=1,ORBFRM
         I=ORBFRM-JUNK+1
         NUM=0
         IA=I*(I-1)/2
         IS=BFSYM(ORBTBF(I))
         DO 39 J=1,I-1
            IJS=IEOR(IS,BFSYM(ORBTBF(J)))
            IJXX(IA+J)=NKLXX(IJS+1,ORBFRM)
            NUM=NUM+NKLXX(IJS+1,ORBFRM)
   39    CONTINUE
C
         DO 40 J=1,I
            IJS=IEOR(IS,BFSYM(ORBTBF(J)))
            IJWW(IA+J)=NKLWW(IJS+1,ORBFRM)
            NUM=NUM+NKLWW(IJS+1,ORBFRM)
   40    CONTINUE
         ININT(I)=ININT(I)+NUM
         NUMINT=NUMINT+NUM
   41 CONTINUE
C
C     ----- WORK OUT A CONVENIENT OUTPUT BLOCK SIZE FOR INTEGRALS -----
C
      JADMAX=0
      DO 45 I=1,NORBS
         IF (ININT(I).GT.JADMAX) JADMAX=ININT(I)
   45 CONTINUE
      NMAX=BLKSIZ
      IF (MAXSIZ-2*NWKS.LT.BLKSIZ) NMAX=MAXSIZ-2*NWKS
      IF (NMAX.LT.JADMAX) NMAX=BLKSIZ
      IF (NMAX.GT.JADMAX) GO TO 11
      WRITE (OUT,10) JADMAX
   10 FORMAT (//,' CANNOT FIT I-BLOCK OF INTEGRALS IN REQUESTED BLOCK'
     #,          ' SIZE',/,' NEED',I7,' INTEGRALS PER BLOCK, SO WILL'
     #,          ' TRY TO HOLD LARGEST I-BLOCK',//)
      NMAX=JADMAX
   11 CONTINUE
      IF (NMAX.LE.ABSMAX) GO TO 13
      WRITE (OUT,12) NMAX
   12 FORMAT (//,' CANNOT HANDLE I-BLOCK OF INTEGRALS OF SIZE',I8
     #,        /,' SO WILL SPLIT I-BLOCKS',//)
      NMAX=ABSMAX
C***********************************************************************
      WRITE(*,*) '  BLKSIZ = ',BLKSIZ
C***********************************************************************
   13 CONTINUE
      IF (NMAX.GT.NUMINT) NMAX=NUMINT
C
C     ----- GENERATE THE IJADD AND IJGRP ARRAYS -----
C
      GROUP=0
      LEFT=0
      SUM=0
      DO 59 I=NORBS,1,-1
         IA=I*(I-1)/2
         IF (ININT(I).LE.LEFT) GO TO 50
         IF (GROUP.GT.0) NINGRP(GROUP)=SUM
         GROUP=GROUP+1
         LEFT=NMAX
         SUM=0
   50    CONTINUE
         IJMAX=IA+I
         DO 55 J=I,1,-1
            IJ=IA+J
            NIJ=IJADD(IJ)+IJWW(IJ)+IJXX(IJ)
            IF (NIJ.GT.NMAX) GO TO 910
            IF (NIJ.LE.LEFT) GO TO 54
            IJMIN=IJ+1
            DO 51 IJQ=IJMAX,IJMIN,-1
               T=IJADD(IJQ)
               IJADD(IJQ)=SUM
               SUM=SUM+T
   51       CONTINUE
            DO 52 IJQ=IJMAX,IJMIN,-1
               T=IJWW(IJQ)
               IJWW(IJQ)=SUM
               SUM=SUM+T
   52       CONTINUE
            DO 53 IJQ=IJMAX,IJMIN,-1
               T=IJXX(IJQ)
               IJXX(IJQ)=SUM
               SUM=SUM+T
   53       CONTINUE
C
            IJMAX=IJ
            NINGRP(GROUP)=SUM
            GROUP=GROUP+1
            SUM=0
            LEFT=NMAX
   54       CONTINUE
            LEFT=LEFT-NIJ
            IJGRP(IJ)=GROUP
   55    CONTINUE
C
         IJMIN=IJ
         DO 56 IJQ=IJMAX,IJMIN,-1
            T=IJADD(IJQ)
            IJADD(IJQ)=SUM
            SUM=SUM+T
   56    CONTINUE
         DO 57 IJQ=IJMAX,IJMIN,-1
            T=IJWW(IJQ)
            IJWW(IJQ)=SUM
            SUM=SUM+T
   57    CONTINUE
         DO 58 IJQ=IJMAX,IJMIN,-1
            T=IJXX(IJQ)
            IJXX(IJQ)=SUM
            SUM=SUM+T
   58    CONTINUE
C
   59 CONTINUE
      NGROUP=GROUP
      WRITE (OUT,21) NUMINT,NGROUP,NMAX
   21 FORMAT (///,T4,'THERE ARE',I9,' INTEGRALS TO BE STORED IN',I3
     #,      ' GROUPS OF',I7,' EACH')
C
C     ----- GENERATE INEXT AND JNEXT ARRAYS -----
C
      DO 20 ORB=1,NORBS
         NEXT=NORBS-ORB+1
         INEXT(NEXT)=ORB
         JMNNXT(NEXT)=1
         JMXNXT(NEXT)=ORB
   20 CONTINUE
      RETURN
C
  910 CONTINUE
      WRITE (OUT,911) NMAX,NIJ
      WRITE (ERROUT,911) NMAX,NIJ
  911 FORMAT (//,' CANNOT FIT IJ-BLOCK OF INTEGRALS. NMAX=',I7
     #,        / '                        SIZE OF IJ-BLOCK',I7,//)
      CALL mabort
      END
