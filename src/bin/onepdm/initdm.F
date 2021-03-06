C//////////////////////////////////////////////////
      SUBROUTINE INITDM(DRTDM,VECTOR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C   INITIALIZE AND PREPARE TO MAKE THE DENSITY MATRIX
      integer nbfsq,nwksmx
      DOUBLE PRECISION DRTDM(NBFSQ),VECTOR(NWKSMX),DC
      LOGICAL PAGEIN
      DOUBLE PRECISION ACF,DRTVER
      CHARACTER SYMLBL*3
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
C   ZERO THE DENSITY MATRIX
      DO 10 IK=1,NBFSQ
   10 DRTDM(IK) = 0.0D0
C   THIS SECTION CREATES THE DENSITY MATRIX FROM DATA PASSED IN BY
C   LOOPY
      RETURN
      ENTRY MAKEDM(drtdm,vector)
      NLWKT=0
      II = IUWK - NOFFI
      JJ = JUWK - NOFFJ
      IF(PAGEIN) GO TO 21
      NLWKT=NLWK
   20 IF(NLWKT.LE.0) GO TO 23
      JJ=1
      II=1
      NLWK=MIN(NLWKT,NWKSMX)
      CALL PAGED(vector)
      PAGEIN=.FALSE.
   21 DC = 0.0D0
      DO 22 I=1,NLWK
      DC = DC + VECTOR(II) * VECTOR(JJ)
      II = II + 1
   22 JJ = JJ + 1
      NLWKT=NWLKT-NLWK
      IUWK=IUWK+NLWK
      GO TO 20
   23 LAD = (IORB-1) * NBF + JORB
      MAD = (JORB-1) * NBF + IORB
C   ADD THIS CONTRIBUTION TO THE DENSITY MATRIX
      DRTDM(LAD) = DRTDM(LAD) + ACF * DC
      DRTDM(MAD) = DRTDM(LAD)
      RETURN
      END
