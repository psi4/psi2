C
C                         *********************
C-------------------------******  RDTWO  ******------------------------
C                         *********************
C
      SUBROUTINE RDTWO(AOINTS,LBLI,PKI,INBUFS,IDIM,ITAP34,JOUT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION LBLI(INBUFS*2),PKI(INBUFS),AOINTS(IDIM,IDIM)
      COMMON/SIGNS/IOFF(256),IPRNT
      COMMON/SYMMS/IDEG(10),NDEG(10),LSYM(10)
      DATA A0 / 0.0D+00 /
    1 FORMAT(2X,' ILSTI = ',I5/
     1       2X,' NBUF  = ',I5/)
    2 FORMAT(2X,9I5,3F15.7)
      CALL ZERO(AOINTS,IDIM*IDIM)
C
      IBUFSZ=INBUFS*2
      IBUFQT=(INBUFS-1)/2+1
      IBL=0
      INTI=0
      INTO=0
      NBLI=0
      NBLO=0
      NINT=0
  201 NBLI=NBLI+1
      CALL SREAD(ITAP34,LBLI,IBUFSZ)
      ILSTI=LBLI(1)
      NBUF=LBLI(2)
      IF(IPRNT.GE.5)
     *WRITE(6,1) ILSTI,NBUF
C***  WRITE(8,15) NBUF
C**15 FORMAT(I5)
      DO 101 I=1,NBUF
      LOR=ISHFT(LBLI(I+I+2),-2)
      LSM=ISHFT(LOR,-8)
      KOR=ISHFT(LSM,-3)
      KSM=ISHFT(KOR,-8)+1
      KOR=IAND(KOR,255)
      LSM=IAND(LSM,7)+1
      LOR=IAND(LOR,255)
      III=IAND(LBLI(I+I+2),3)
      JSM=ISHFT(LBLI(I+I+1),-8)
      IOR=ISHFT(JSM,-3)
      ISM=ISHFT(IOR,-8)+1
      IOR=IAND(IOR,255)
      JSM=IAND(JSM,7)+1
      JOR=IAND(LBLI(I+I+1),255)
      II=IOR+IDEG(ISM)
      JJ=JOR+IDEG(JSM)
      KK=KOR+IDEG(KSM)
      LL=LOR+IDEG(LSM)
CTJL  WRITE(*,*) '  IOR,JOR,KOR,LOR',IOR,JOR,KOR,LOR
CTJL  WRITE(*,*) '  II,JJ,KK,LL ',II,JJ,KK,LL,PKI(I+IBUFQT)
      PKII=PKI(I+IBUFQT)
      IJ=IOFF(II)+JJ
      KL=IOFF(KK)+LL
      AOINTS(IJ,KL)=PKII
      AOINTS(KL,IJ)=PKII
  101 CONTINUE
      INTI=INTI+NBUF
      IF(ILSTI.EQ.0) GO TO 201
C
C     WRITE (JOUT,*) ' **********  AOINTS  *********',AOINTS
      RETURN
      END
