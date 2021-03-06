cets110289      PROGRAM READ34
      subroutine fentry(cc,ia,maxcor)
C   THIS PROGRAM READS INTEGRAL FILE34
C******************************************************
C***LAST UPDATED ON January 18, 2005 BY YUKIO YAMAGUCHI*
C******************************************************
C******************************************************
C***LAST UPDATED ON JULY 29, 1985 BY YUKIO YAMAGUCHI***
C******************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C     DIMENSION CC(1200000),IA(1)
cets110289      DIMENSION CC(360000),IA(1)
      dimension cc(maxcor),ia(1)
      INTEGER BLABEL(20)
      INTEGER ND(10),MNL(142),MS(142),MTYPE(2,20),ITYP(20)
      CHARACTER*8 TAPE
      LOGICAL TAPE36
      COMMON/BASIS/NBASIS,NTRI,NSYM,NBFAO,NBFSO,NBATRI
      COMMON/ENRGY/ESCF,ENUC
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4
      COMMON/SIGNS/IOFF(256),IPRNT
      COMMON/SYMMS/IDEG(10),NDEG(10),LSYM(10)
      COMMON/TAPES/ITAP34,ITAP36
cets110289      EQUIVALENCE (CC,IA)
    1 FORMAT(//,2X,' THE PROGRAM TO READ TAPE34'/)
    2 FORMAT(5I5)
    3 FORMAT(A8,2X)
    4 FORMAT(2X,' BLABEL = ',20A4/)
    5 FORMAT(2X,' ENUC   = ',F15.8/)
    6 FORMAT(2X,' NST    = ',I5/2X,' NDEG   =',8I5/)
    7 FORMAT(2X,' NBASIS = ',I8/
     1       2X,' NTRI   = ',I8/
     2       2X,' NSYM   = ',I8/
     3       2X,' NBFAO  = ',I8/
     4       2X,' NBFSO  = ',I8/
     5       2X,' NBATRI = ',I8/
     6       2X,' IPRNT  = ',I8/
     7       2X,' TAPE   = ',A8)
    8 FORMAT(//,2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)
    9 FORMAT(//,2X,' REQUIRED MEMORY EXCEEDS MAXCOR '/
     1          2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)
   10 FORMAT(20A4)
   11 FORMAT(F20.10)
C
C   ITAP34 = SO INTEGRALS
C   ITAP36 = TWO ELECTRON SO INTEGRALS
C
cets  open(unit=8,status='new',file='file8')
      call drum
      CALL TSTART(6)
C
      ITAPE3=3
      ITAPE6=6
      ITAPE8=8
      ITAP34=34
      ITAP36=36
      call ffile(itape8,'read34',0)
C     MAXCOR=1200000
cets110289      MAXCOR=360000
      MAXBUF=4096
      MAXBF2=MAXBUF*2
      MAXBF4=MAXBUF*4
      TAPE36=.TRUE.
      INPUT=5
      CALL LOCATE(INPUT,'# READ34 #',IERR)
C
      WRITE(6,1)
      WRITE(3,1)
C
      READ(5,2) IPRNT
c     READ(5,3) TAPE
C
      IOFF(1)=0
      DO 101 I=1,255
  101 IOFF(I+1)=IOFF(I)+I
C
C   READ IN INTEGRAL INFORMATION FROM 34 FILE
      CALL RFILE(ITAP34)
      CALL WREADW(ITAP34,NKIND,1,1,IEND)
      CALL WREADW(ITAP34,ICLOS,1,IEND,IEND)
      CALL WREADW(ITAP34,BLABEL,20,IEND,IEND)
      CALL WREADW(ITAP34,ENUC,2,IEND,IEND)
      CALL WREADW(ITAP34,NST,1,IEND,IEND)
      CALL WREADW(ITAP34,ND,NST,IEND,IEND)
      CALL WREADW(ITAP34,ITYP,NST,IEND,IEND)
      CALL WREADW(ITAP34,NDEG,NST,IEND,IEND)
      CALL WREADW(ITAP34,NS,1,IEND,IEND)
      CALL WREADW(ITAP34,MTYPE,2*NS,IEND,IEND)
      CALL WREADW(ITAP34,LAPU,1,IEND,IEND)
      CALL WREADW(ITAP34,MS,LAPU,IEND,IEND)
      CALL WREADW(ITAP34,MNL,LAPU,IEND,IPARA)
C     WRITE(6,30) IPARA
C  30 FORMAT(2X,' IPARA = ',I10/)
      ISADR=I2SEC(IPARA)+1
      CALL RSETSA(ITAP34,ISADR)
C
      WRITE(6,4) BLABEL
C     WRITE(3,4) BLABEL
C     WRITE(6,5) ENUC
      WRITE(3,5) ENUC
C     WRITE(6,6) NST,(NDEG(I),I=1,NST)
      WRITE(3,6) NST,(NDEG(I),I=1,NST)
C
      NSYM=0
      NBASIS=0
      DO 102 I=1,NST
      IF(NDEG(I).GT.0) NSYM=NSYM+1
      NBASIS=NBASIS+NDEG(I)
  102 CONTINUE
      IDEG(1)=0
      LSYM(1)=0
      IF(NSYM.EQ.1) GO TO 203
      II=0
      DO 103 I=2,NST
      IF(NDEG(I).LE.0) GO TO 202
  201 II=II+1
      NN=NDEG(II)
      IF(NN.EQ.0) GO TO 201
      IDEG(I)=IDEG(I-1)+NN
      LSYM(I)=LSYM(I-1)+IOFF(NN+1)
      GO TO 103
  202 LSYM(I)=LSYM(I-1)
      IDEG(I)=IDEG(I-1)
  103 CONTINUE
  203 NTRI=IOFF(NBASIS+1)
C
      IF(TAPE.EQ.'NOTAPE36') TAPE36=.FALSE.
      IF(TAPE.EQ.'NOFILE36') TAPE36=.FALSE.
      IF(TAPE36) TAPE='  TAPE36'
      WRITE(6,7) NBASIS,NTRI,NSYM,NBFAO,NBFSO,NBATRI,IPRNT,TAPE
C
C   FORM FILE8
      REWIND ITAPE8
      WRITE(8,10) BLABEL
      WRITE(8,11) ENUC
C
C   CREATE SO INTEGRAL FILES
      WRITE(3,21)
   21 FORMAT(//,2X,' NOW YOU ARE IN RDONE'/)
      IBUFHF=4471
      IC1=1
      IC2=IC1+NTRI
      IC3=IC2+NTRI
      IC4=IC3+NTRI
      IC5=IC4+NTRI
      IA5=IC5+IC5-1
      ICMAX=IC5+IBUFHF
      WRITE(3,8) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 399
C................S       H       T       V       LBLI    STVI....
      CALL RDONE(CC(IC1),CC(IC2),CC(IC3),CC(IC4),IA(IA5),CC(IC5),
C.......................
     1           IBUFHF)
C
      IF(.NOT.TAPE36) GO TO 400
      WRITE(3,22)
   22 FORMAT(//,2X,' NOW YOU ARE IN RDTWO'/)
      IBUFHF=4471
      IC1=1
      IA1=1
      IC2=IC1+IBUFHF
      IA2=IC2+IC2-1
C................LBLI    PKI     LBLO    BUFO...........
      CALL RDTWO(IA(IA1),CC(IC1),IA(IA2),CC(IC2),IBUFHF)
      GO TO 400
C
  399 CONTINUE
      WRITE(3,9) ICMAX,MAXCOR
      WRITE(3,9) ICMAX,MAXCOR
      STOP
C
  400 CONTINUE
      REWIND ITAPE8
      CALL RCLOSE(ITAP34,3)
      CALL TSTOP(6)
      STOP
      END
      SUBROUTINE RDONE(S,H,T,V,LBLI,STVI,INBUFS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION S(NTRI),H(NTRI),T(NTRI),V(NTRI)
      DIMENSION LBLI(INBUFS*2),STVI(INBUFS)
      COMMON/BASIS/NBASIS,NTRI,NSYM,NBFAO,NBFSO,NBATRI
      COMMON/SIGNS/IOFF(256),IPRNT
    1 FORMAT(//,2X,' S MATRIX'/)
    2 FORMAT(//,2X,' T MATRIX'/)
    3 FORMAT(//,2X,' V MATRIX'/)
    4 FORMAT(//,2X,' H MATRIX'/)
C
C   READ IN ONE ELECTRON INTEGRALS
C
C***  ITAPE8=8
C***  REWIND ITAPE8
C   S INTEGRALS
      CALL ONEI(S,LBLI,STVI,INBUFS)
      IF(IPRNT.LE.3) GO TO 201
      WRITE(6,1)
      CALL PRINT(S,NTRI,NBASIS,6)
C
C   T INTEGRALS
  201 CONTINUE
      CALL ONEI(T,LBLI,STVI,INBUFS)
      IF(IPRNT.LE.4) GO TO 202
      WRITE(6,2)
      CALL PRINT(T,NTRI,NBASIS,6)
C
C   V INTEGRALS
  202 CONTINUE
      CALL ONEI(V,LBLI,STVI,INBUFS)
      IF(IPRNT.LE.4) GO TO 203
      WRITE(6,3)
      CALL PRINT(V,NTRI,NBASIS,6)
C
C   TOTAL ONE ELECTRON INTEGRALS
  203 CONTINUE
      DO 101 I=1,NTRI
  101 H(I)=T(I)+V(I)
      IF(IPRNT.LE.3) GO TO 204
      WRITE(6,4)
      CALL PRINT(H,NTRI,NBASIS,6)
C
  204 CONTINUE
CCC   WRITE(6,*) ' ONE ELECTRON INTEGRALS---S , T, V, H'
CCC   DO 102 I=1,NTRI
CCC   WRITE(6,5) I,S(I),T(I),V(I),H(I)
CC102 CONTINUE
CCC 5 FORMAT(2X,I5,5X,4F20.10)
C
      RETURN
      END
      SUBROUTINE ONEI(X,LBLI,STVI,INBUFS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(NTRI),LBLI(INBUFS*2),STVI(INBUFS)
      COMMON/BASIS/NBASIS,NTRI,NSYM,NBFAO,NBFSO,NBATRI
      COMMON/SIGNS/IOFF(256),IPRNT
      COMMON/SYMMS/IDEG(10),NDEG(10),LSYM(10)
      COMMON/TAPES/ITAP34,ITAP36
      DATA ZERO / 0.0D+00 /
    1 FORMAT(2X,' ILSTI = ',I5/
     1       2X,' NBUF  = ',I5/)
    2 FORMAT(2X,8I5,3F15.7)
C
      DO 101 I=1,NTRI
  101 X(I)=ZERO
C
      IBUFSZ=8942
      IBUFS3=1491
      INTI=0
      NBLI=0
  201 NBLI=NBLI+1
      CALL SREAD(ITAP34,LBLI,IBUFSZ)
      ILSTI=LBLI(1)
      NBUF=LBLI(2)
      IF(IPRNT.GE.5)
     *WRITE(6,1) ILSTI,NBUF
      WRITE(8,15) NBUF
   15 FORMAT(I5)
      DO 102 I=1,NBUF
      JSM=ISHFT(LBLI(I+2),-8)
      IOR=ISHFT(JSM,-3)
      ISM=ISHFT(IOR,-8)+1
      IOR=IAND(IOR,255)
      JSM=IAND(JSM,7)+1
      JOR=IAND(LBLI(I+2),255)
      II=IOR+IDEG(ISM)
      JJ=JOR+IDEG(JSM)
      IIJJ=IOFF(II)+JJ
      STV0=STVI(I+IBUFS3)
      X(IIJJ)=STV0
      WRITE(8,10) II,JJ,STV0
      IF(IPRNT.LE.6) GO TO 102
      WRITE(6,2) I,IOR,ISM,JOR,JSM,II,JJ,IIJJ,STV0
   10 FORMAT(2I5,F20.10)
  102 CONTINUE
      INTI=INTI+NBUF
      IF(ILSTI.EQ.0) GO TO 201
C
      RETURN
      END
      SUBROUTINE RDTWO(LBLI,PKI,LBLO,BUFO,INBUFS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION LBLI(INBUFS*2),PKI(INBUFS),LBLO(MAXBF4),BUFO(MAXBF2)
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4
      COMMON/SIGNS/IOFF(256),IPRNT
      COMMON/SYMMS/IDEG(10),NDEG(10),LSYM(10)
      COMMON/TAPES/ITAP34,ITAP36
      DATA ZERO / 0.0D+00 /
    1 FORMAT(2X,' ILSTI = ',I5/
     1       2X,' NBUF  = ',I5/)
    2 FORMAT(2X,9I5,3F15.7)
C
      CALL RFILE(ITAP36)
C
      IBUFSZ=8942
      IBUFQT=2236
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
      WRITE(8,15) NBUF
   15 FORMAT(I5)
      WRITE(6,*)'   III,IOR,ISM,JOR,JSM,KOR,KSM,LOR,LSM,PKII'
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
      PKII=PKI(I+IBUFQT)
C     IF(IPRNT.LE.4) GO TO 202
C     IJ=IOFF(IOR)+JOR+LSYM(ISM)
C     IF(IOR.LT.JOR) IJ=IOFF(JOR)+IOR+LSYM(JSM)
C     KL=IOFF(KOR)+LOR+LSYM(KSM)
C     IF(KOR.LT.LOR) KL=IOFF(LOR)+KOR+LSYM(LSM)
      WRITE(6,2) III,IOR,ISM,JOR,JSM,KOR,KSM,LOR,LSM,PKII
C     WRITE(6,20) II,JJ,KK,LL,PKII
      WRITE(8,20) II,JJ,KK,LL,PKII
   20 FORMAT(4I5,F20.10)
  202 CALL PACK(II,JJ,KK,LL,IX,JA,JAA)
      IBL=IBL+1
      LBLO(IBL+IBL-1)=JA
      LBLO(IBL+IBL)=JAA
      BUFO(IBL+MAXBUF)=PKII
      IF(IBL.LT.MAXBUF) GO TO 101
      IBL=0
      NBLO=NBLO+1
      CALL SWRIT(ITAP36,LBLO,MAXBF4)
  101 CONTINUE
      INTI=INTI+NBUF
      IF(ILSTI.EQ.0) GO TO 201
C
      IBL=IBL+1
      NBLO=NBLO+1
      INTO=INTO+IBL
      LBLO(IBL+IBL-1)=0
      LBLO(IBL+IBL)=0
      BUFO(IBL+MAXBUF)=ZERO
      CALL SWRIT(ITAP36,LBLO,MAXBF4)
      CALL SREW(ITAP36)
C
      CALL RCLOSE(ITAP36,3)
      RETURN
      END

