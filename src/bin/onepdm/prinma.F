C///////////////////////////////////////////////////
      SUBROUTINE PRINMA(ARRAY,EIG,SYM,N,IFLG,NN)
C     EXTENDED DUMMY ARRAY,EIG
      DOUBLE PRECISION ARRAY(NN,NN),EIG(NN)
      CHARACTER SYMLBL*3,SYM*3
      LOGICAL PAGEIN
      DOUBLE PRECISION ACF,DRTVER
      INTEGER BMAX,PUWK,PRPFLG,PRINT,ROOTI,ROOTF,SYMOFF,IRSYM
      INTEGER DRTLBL,SCFLBL,SYMDIF
      COMMON /TAPES/IFILE3,IFILE5,IFILE6,ITAPE7,IFIL58,IFIL54,IFIL50
     $,             IFIL51,IDEV,IFIL30
      COMMON /DIMS/ NBF,NSYM,NORBS,NROWS,NROWS4,NWKS,NWKS2,NLEVS
     $,             NWKSMX,NWKMX2,NLWKMX,BMAX,NROOTS,NORBSQ,NBFSQ
     $,             NOCCMX,NREFS,IEND
      COMMON /INTS/ NMAX,NMAX2,NGROUP,NBLKOC,NUMIJ,NUMKL,INTSRT
      COMMON /LOOPS/ACF,IUWK,JUWK,NLWK,IORB,JORB,PUWK
      COMMON /PAGE/ IWORD3,PAGEIN,NOFFI,NOFFJ,NPASS,NWORDL
      COMMON /OPS/  PRINT,MAX,NOFLG,PRPFLG,ROOTI,ROOTF,IRSYM
      COMMON /SYMM/ SYMOFF(8),NTOTS(8,2),SYMDIF,SYMLBL(15)
      COMMON /LBLS/ DRTVER,DRTLBL(26),INTLBL(26),SCFLBL(26)
      INC=4
      IF(IDEV.EQ.0) INC=9
      MINI = 1
   10 MAXI = N
      IF((MAXI-MINI).GT.INC) MAXI = MINI + INC
      WRITE(IFILE6,20)SYM,(K,K=MINI,MAXI)
   20 FORMAT(/,1H ,'IRREP ',A3,3X,10(I3,9X))
      DO 30 I=1,N
   30 WRITE(IFILE6,31)(ARRAY(I,J),J=MINI,MAXI)
   31 FORMAT(1H ,6X,10F12.8)
      IF(IFLG.EQ.0) GO TO 40
      WRITE(IFILE6,32)(EIG(J),J=MINI,MAXI)
   32 FORMAT(1H ,'EIGENVALUES',/,1H ,6X,10F12.8)
   40 MINI = MAXI + 1
      IF(MAXI.LT.N) GO TO 10
      RETURN
      END
