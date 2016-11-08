      SUBROUTINE xMAIN(A,IA,MAXCOR)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER IA(MAXCOR)
      REAL*8   A(1)
      COMMON /INOUT/  IOUT,ICHECK,ITAP30
      COMMON /MAXDIM/ MXAT,MXTRN,MXPR,MXSH,MXSYM,MXCF,MXT,MXISO,MXPSYM
      COMMON /PSIZE/  MSU,KAORDS,MCONSU,MGCSU,MSTU,MRU,MCSU,MCTU,MCONU,
     1                MCU,MSFU,MGU,MSFRU,MNSFU,MPRU,MCXU,MCCU,MBLU,MRCRU
C***********************************************************************
C     SET UP JOB.  THE FOLLOWING PARAMETERS ARE DEFINED AS:            *
C  MXAT  = MAXIMUM # OF ATOMS          MXTRN  =                        *
C  MXPR  = MAXIMUM # OF PRIMITIVES     MXSH   = MAXIMUM # OF SHELLS    *
C  MXSYM =                             MXCF   =                        *
C  MXISO =                             MXPSYM =                        *
C  MRCRU = MAXIMUM # OF CONTRACTIONS FOR A SET OF PRIMITIVES           *
C***********************************************************************
      ICHECK = 3
      IOUT   = 6
      ITAP30 = 30
      MXAT=40
      MXTRN=48
      MXPR=300
      MXSH=100
      MXSYM=14
      MXCF=10000
      MXT=MXTRN*9
      MXISO=MXTRN/4
      MXPSYM=(MXSYM+3)/4
c     MRCRU=6
cges dec 12 1989 increase mrcru to 9 general contractions
c     MRCRU=9
cges may 24 1989 increase mrcru to 10 general contractions
      MRCRU=10
C***********************************************************************
C     START READS IN BASIS SET AND GET INITIAL MO'S.                   *
C***********************************************************************
      CALL START(A,IA,MAXCOR)
      RETURN
C
      END
