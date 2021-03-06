C////////////////////////////////////////////////////
      SUBROUTINE DENSMA(DRTDM,TRANS,SCFDM,SCFVEC,OCC,OCCNO,XSO,
     $      ORBTBF,BFCODE,JUNK,VECWRT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
cets030291
#include <error.h>
      integer prcntr, ip
      integer nbf,norbs
      INTEGER FWORD,FZC,ORBTBF,BFCODE,SYM,ICBN13(16),LOCVEC
      DOUBLE PRECISION SCFVEC,OCC,SUM,OCCNO,XSO,XXX
      DOUBLE PRECISION DRTDM(NBF,NBF),TRANS(NBF,NBF),SCFDM(NBF,NBF)
      CHARACTER SYMLBL*3,CHTP*3
      DIMENSION SCFVEC(NBF,NBF),ORBTBF(NORBS),OCC(NBF),JUNK(NBF)
      DIMENSION XSO(1),OCCNO(NBF),VECWRT(1)
crsg021594 use correct dimensions for BFCODE so MRCI output is correct
crsg021594      DIMENSION BFCODE(NBF),FMT(26),ILBL(26)
      DIMENSION BFCODE(nrefs,NBF),FMT(26),ILBL(26)
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
      DATA FZC/1/
C    FORM THE SCF DENSITY MATRIX
      IJL=0
      IJM=0
      IJN=0
      DO 4  I=1,NBF*(NBF+1)/2
    4 XSO(I)=0.0D0
      DO 20 I=1,NBF
      DO 10 J=1,I
      SCFDM(J,I) = 0.0D0
   10 SCFDM(I,J) = 0.0D0
crsg021594 a kludge for BFCODE dimension problem
crsg021594      IF(BFCODE(I).EQ.FZC) SCFDM(I,I) = 2.0D0
      IF(BFCODE(1,I).EQ.FZC) SCFDM(I,I) = 2.0D0
   20 CONTINUE
      DO 30 I=1,NORBS
      DO 30 J=1,I
      II = ORBTBF(I)
      JJ = ORBTBF(J)
c      drtdm(x,y) is the OPDM in DRT-ordered MO basis
      ELE = DRTDM(I,J)
c      scfdm(x,y) is the OPDM in SCF-ordered MO basis
      SCFDM(JJ,II) = ELE
   30 SCFDM(II,JJ) = ELE
      WRITE(IFILE6,31)
c      printing out diagonal elements of the OPDMMO
c      those aren't exactly good definitions occupations
c      eigenvalues of OPDMMO is a much better definition :-)
   31 FORMAT(//,1H ,32X,'POPULATIONS OF THE MOLECULAR ORBITALS'/)
      WRITE(IFILE6,32) (I,SCFDM(I,I),I=1,NBF)
   32 FORMAT(8(I3,1X,F10.8))
      IF(PRPFLG.LT.1) GO TO 50
      DO 40 I=1,NBF
      DO 40 J=1,NBF
   40 SCFVEC(I,J)=0.0D0
      REWIND IFIL51
      CALL SCFLCV(LOCVEC)
      LOCVC2=LOCVEC
   50 CONTINUE
      IORB=0
      DO 190 L=1,NSYM
      NTOT=NTOTS(L,1)
      IF (NTOT.EQ.0) GO TO 190
      SYM=NTOTS(L,2)
      CHTP=SYMLBL(SYM+SYMDIF)
      II=IORB
      DO 70 I=1,NTOT
      JJ=IORB
      II=II+1
      DO 70 J=1,NTOT
      JJ=JJ+1
      SCFVEC(I,J)=SCFDM(II,JJ)
      DRTDM(I,J)=SCFDM(II,JJ)
   70 CONTINUE
      IORB=IORB+NTOT
c      hdiag doesn't preserve DMTDM
c      DMTDM will have eigenvalues on the diagonal and zeroes elsewhere
c      trans(x,y) is just NMO matrix in MO basis
c      occ is just a buffer
      CALL HDIAG(DRTDM,TRANS,OCC,NBF,NTOT,JUNK,0)
cets030291
      ip = prcntr('IS_ON OPDMMO')
      IF(PRINT.EQ.0.and.ip.eq.0) GOTO 90
      DO 80 I=1,NTOT
c      occ() contains NMO occupancies
   80 OCC(I)=DRTDM(I,I)
      IF(L.EQ.1) WRITE(IFILE6,81)
   81 FORMAT(//,1H1,21X,'DENSITY MATRIX OVER MOLECULAR ORBITALS',/)
      CALL PRINMA(SCFVEC,OCC,CHTP,NTOT,1,NBF)
   90 CONTINUE
      IF(PRPFLG.LT.1) GO TO 190
      DO 92  IJK=1,NTOT
      IJL=IJL+1
c      occno() - occupancies of NMOs
   92 OCCNO(IJL)=OCC(IJK)
c      Now reading the SCF eigenvector
      DO 96  J=1,NTOT
      CALL WREADW(IFIL30,OCC,INTOWP(NTOT),LOCVEC,LOCVEC)
      DO 96  I=1,NTOT
   96 SCFVEC(I,J)=OCC(I)
      DO 150 I=1,NTOT
      DO 150 J=1,NTOT
      SUM=0.0D0
CCC  BEFORE 7/12/84 TRANS(J,K) APPEARED HERE INCORRECTLY.  W. ALLEN .
c      Transforming NMOs to SO basis
      DO 140 K=1,NTOT
  140 SUM=SUM+SCFVEC(I,K)*TRANS(K,J)
      SCFDM(I,J)=SUM
  150 CONTINUE
c      Thus file51 contains NMOs in terms of SOs and MOs
      DO 170 I=1,NTOT
  170 WRITE(IFIL51,172) (SCFDM(J,I),J=1,NTOT)
      DO 171 I=1,NTOT
  171 WRITE(IFIL51,172) (TRANS(J,I),J=1,NTOT)
  172 FORMAT(1X,5E20.12)
  190 CONTINUE
      IF(PRPFLG.LT.1) RETURN
c      Reading SCF eigenvalues from file30
      CALL WREADW(IFIL30,OCC,INTOWP(NBF),LOCVEC,LOCVEC)
      WRITE(IFILE6,230)
  230 FORMAT(//,1X,'MOLECULAR ORBITAL ENERGIES',/)
      WRITE(IFILE6,33) (I,OCC(I),I=1,NBF)
   33 FORMAT(8(I3,1X,F10.8))
      WRITE(IFILE6,231)
  231 FORMAT(//)
cets030291
      ip=prcntr('IS_ON NOMO')
      IF(PRINT.LT.2.and.ip.eq.0) GO TO 301

      REWIND IFIL51
      DO 300  L=1,NSYM
      NTOT=NTOTS(L,1)
      IF(NTOT.EQ.0) GO TO 300
      SYM=NTOTS(L,2)
      CHTP=SYMLBL(SYM+SYMDIF)
           DO 302  I=1,NTOT
  302      READ(IFIL51,172) (SCFDM(J,I),J=1,NTOT)
           DO 304  I=1,NTOT
  304      READ(IFIL51,172) (TRANS(J,I),J=1,NTOT)
      IF(L.EQ.1) WRITE(IFILE6,200)
  200 FORMAT(//,1H1,20X,'NATURAL ORBITALS IN TERMS OF MOLECULAR'
     $       ,' ORBITALS',/)
      CALL PRINMA(TRANS,OCC,CHTP,NTOT,0,NBF)
  300 CONTINUE
  301 REWIND IFIL51
      IJM=0
      DO 310  L=1,NSYM
      NTOT=NTOTS(L,1)
      IF(NTOT.EQ.0) GO TO 310
      SYM=NTOTS(L,2)
      CHTP=SYMLBL(SYM+SYMDIF)
           DO 312  I=1,NTOT
  312      READ(IFIL51,172) (SCFDM(J,I),J=1,NTOT)
           DO 314  I=1,NTOT
  314      READ(IFIL51,172) (TRANS(J,I),J=1,NTOT)
cets030291
      ip=prcntr('IS_ON NOSO')
      IF(PRINT.LT.3.and.ip.eq.0) GO TO 217
      IF(L.EQ.1) WRITE(IFILE6,210)
  210 FORMAT(//,1H1,20X,'NATURAL ORBITALS IN TERMS OF SYMMETRY'
     $       ,' ORBITALS',/)
      CALL PRINMA(SCFDM,OCC,CHTP,NTOT,0,NBF)
  217 CONTINUE
CWA   WRITE THE NO'S TO FILE30.
      IF(PRPFLG.EQ.2) THEN
           KK=0
           DO 316  J=1,NTOT
           DO 316  I=1,NTOT
           KK=KK+1
 316       VECWRT(KK)=SCFDM(I,J)
           CALL WWRITW(IFIL30,VECWRT,INTOWP(KK),LOCVC2,LOCVC2)
           WRITE(IFILE6,211) PRPFLG
 211  FORMAT(//,1X,'WARNING: NATURAL ORBITALS HAVE BEEN WRITTEN ',
     $  'TO FILE30.'/,1X,'PRPFLG=',I5)
      END IF
CWA   FORM THE SO 1-PDM AND WRITE IT TO FILE50.
      DO 322  I=1,NTOT
      DO 322  J=1,I
        XXX=0.0D0
        IJL=IJM
        DO 324  K=1,NTOT
        IJL=IJL+1
  324   XXX=XXX+OCCNO(IJL)*SCFDM(I,K)*SCFDM(J,K)
      IJI=I+IJM
      IJJ=IJI*(IJI-1)/2+IJM+J
  322 XSO(IJJ)=XXX
      IJM=IJM+NTOT
  310 CONTINUE
      KK=NBF*(NBF+1)/2
      CALL SWRIT(IFIL50,XSO,INTOWP(KK))
      call fioflu(ifil50)
      RETURN
C 900 WRITE(IFILE6,905) CHTP,NTOT,NTOTS(L,1)
C 905 FORMAT (1X,'VECTOR DOES NOT MATCH DRT, SYMMETRY ',A3,' SCF NTOT='
C    $,       I3,' DRT NTOT=',I3)
C     CALL MABORT
      END
