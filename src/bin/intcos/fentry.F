      subroutine fentry(cc,ia,maxcor)
c     PROGRAM INTCOS                                                    INT00010
C  THE CARTESIAN-INTERNAL COORDINATE CONVERSION PROGRAM                 INT00020
C***********************************************************
C***LAST UPDATED ON April 19, 2003 BY YUKIO YAMAGUCHI      *
c   Adjusted for the LINUX system                          *
C***********************************************************
C***LAST UPDATED ON December 03, 1996 BY YUKIO YAMAGUCHI   *
c   Units are added to geometry and gradients              *
C***********************************************************  
C***LAST UPDATED ON JULY 16, 1992 BY YUKIO YAMAGUCHI       *   
C   EXPLICIT OUT-OF-PLANE ANGLES WERE ADDED                *    
C***********************************************************  
C   THIS PROGRAM CAN TREAT MOLECULES WITH '50' ATOMS                    INT00070
C   DIMENSION STATEMENTS                                                INT00080
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT00090
      dimension cc(maxcor),ia(maxcor*2)
      DIMENSION I30(200),COORD(3,50),GRAD(3,50)                         INT00110
C   COMMON STATEMENTS                                                   INT00120
      COMMON/VIB101/NATOM,NAT1,N3N,NINT                                 INT00130
      COMMON/VIB102/ESCF                                                INT00140
      COMMON/VIB103/NST,NBND,NLIB,NDEF,NTORS                            INT00150
      COMMON/VIB104/IFORCE,IGRAD,IATOM,ICOOD,IPRNT                      INT00160
      COMMON/VIB201/CHG(50),X(50),Y(50),Z(50)                           INT00170
      COMMON/VIB202/R(1275)                                             INT00180
      COMMON/VIB203/ENGX(50),ENGY(50),ENGZ(50)                          INT00190
      COMMON/VIB204/RINT(150),ENGR(150)                                 INT00200
      COMMON/VIB205/IOFF(150)                                           INT00210
      COMMON/VIB301/KR(150),LR(150)                                     INT00220
      COMMON/VIB302/KA(150),LA(150),MA(150),AGLB(150)                   INT00230
      COMMON/VIB303/KB(100),LB(100),MB(100),AGLL(100)                   INT00240
      COMMON/VIB304/KD(100),LD(100),MD(100),ND(100)                     INT00250
      COMMON/VIB305/KT(100),LT(100),MT(100),NT(100)                     INT00260
      COMMON/VIB306/AGLT(100),AGLO(100),AGLR(100)                       INT00270
      COMMON/VIB307/AGLD(100),AGLE(100),AGLF(100)                       INT00280
C   FORMAT STATEMENTS                                                   INT00300
    1 FORMAT(//,2X,' THE INTERNAL ENERGY GRADIENT PROGRAM'/
     *         10x,' Version of  4/19/03'//)
    2 FORMAT(10I5)                                                      INT00320
    3 FORMAT(3X,20A4)                                                   INT00330
    4 FORMAT(//,2X,' NATOM = ',I5/                                      INT00340
     1 2X,' SCF ENERGY = ',F20.10)                                      INT00350
    5 FORMAT(//,2X,' MOLECULAR GEOMETRY (IN BOHR)'/
     1 5X,4H NO.,11X,7H CHARGE,13X,2H X,18X,2H Y,18X,2H Z/)             INT00370
    6 FORMAT(2X,I7,5X,4F20.10)                                          INT00380
    7 FORMAT(//,2X,' CARTESIAN ENERGY GRADIENTS (IN HARTREE/BOHR)'/
     1 5X,4H NO.,11X,7H CHARGE,13X,2H X,18X,2H Y,18X,2H Z/)             INT00400
    8 FORMAT(I5,F20.10)                                                 INT00410
    9 FORMAT(4F20.10)                                                   INT00420
   10 FORMAT(20A4)                                                      INT00430
C                                                                       INT00440
      call psinit('APPEND')
      CALL TSTART(6)                                                    INT00450
C                                                                       INT00470
      INPUT=5                                                           INT00480
      ITAPE6=6                                                          INT00490
      ITAP11=11                                                         INT00500
      ITAP12=12                                                         INT00500
      ITAP13=13                                                         INT00510
      ITAP15=15                                                         INT00510
      ITAP30=30                                                         INT00520
C                                                                       INT00540
      call ffile(itap11,'file11',0)
      call ffile(itap12,'file12',0)
      call ffile(itap13,'file13',0)
      call ffile(itap15,'file15',0)
      CALL LOCATE(INPUT,'# INTCOS #',IERR)                              INT00550
C                                                                       INT00560
      IOFF(1)=0                                                         INT00570
      DO 101 I=1,149                                                    INT00580
      IOFF(I+1)=IOFF(I)+I                                               INT00590
  101 CONTINUE                                                          INT00600
      WRITE(6,1)                                                        INT00610
      READ(5,2) NST,NBND,NLIB,NDEF,NTORS,IFORCE,IGRAD,IATOM,ICOOD,      INT00620
     1          IPRNT                                                   INT00630
      IF(IGRAD.NE.0) GO TO 201                                          INT00640
C                                                                       INT00650
      CALL RFILE(ITAP30)                                                INT00660
      CALL WREADW(ITAP30,I30,200,101,JUNK)                              INT00670
      IEND=I30(1)                                                       INT00680
      MPOINT=I30(2)                                                     INT00690
      MCONST=I30(3)                                                     INT00700
      MCALCS=I30(4)                                                     INT00710
      NCALCS=I30(5)                                                     INT00720
      NATOM=I30(19)                                                     INT00730
      NAT1=NATOM-1                                                      INT00740
      N3N=NATOM*3                                                       INT00750
      JUNK=101+MCONST                                                   INT00760
      CALL WREADW(ITAP30,IPT,1,JUNK,JUNK)                               INT00770
      CALL WREADW(ITAP30,CHG,NATOM*2,IPT,JUNK)                          INT00780
      IGEOP=100+MCONST+MPOINT+NCALCS                                    INT00790
      CALL WREADW(ITAP30,LOCCAL,1,IGEOP,JUNK)                           INT00800
      CALL WREADW(ITAP30,I30,20,LOCCAL,JUNK)                            INT00810
      WRITE(6,3) (I30(I),I=1,20)                                        INT00820
      IGEOP=LOCCAL+60                                                   INT00830
      CALL WREADW(ITAP30,I30,20,IGEOP,IGEOP)                            INT00840
      CALL WREADW(ITAP30,COORD,N3N*2,IGEOP,JUNK)                        INT00850
      JUNK=JUNK+INTOWP(1)                                               INT00860
      CALL WREADW(ITAP30,ESCF,2,JUNK,JUNK)                              INT00870
      JUNK=I30(2)                                                       INT00880
      CALL WREADW(ITAP30,GRAD,N3N*2,JUNK,JUNK)                          INT00890
      WRITE(6,4) NATOM,ESCF                                             INT00900
      WRITE(6,5)                                                        INT00910
      DO 102 I=1,NATOM                                                  INT00920
      X(I)=COORD(1,I)                                                   INT00930
      Y(I)=COORD(2,I)                                                   INT00940
      Z(I)=COORD(3,I)                                                   INT00950
      WRITE(6,6) I,CHG(I),X(I),Y(I),Z(I)                                INT00960
  102 CONTINUE                                                          INT00970
      CALL DIST(IPRNT)                                                  INT00980
      WRITE(6,7)                                                        INT00990
      DO 103 I=1,NATOM                                                  INT01000
      ENGX(I)=GRAD(1,I)                                                 INT01010
      ENGY(I)=GRAD(2,I)                                                 INT01020
      ENGZ(I)=GRAD(3,I)                                                 INT01030
      WRITE(6,6) I,CHG(I),ENGX(I),ENGY(I),ENGZ(I)                       INT01040
  103 CONTINUE                                                          INT01050
      GO TO 203                                                         INT01060
C                                                                       INT01070
  201 CONTINUE                                                          INT01080
      NTAPE=IGRAD/100                                                   INT01090
      IF(NTAPE.EQ.13) GO TO 202                                         INT01100
      NSETS=MOD(IGRAD,100)                                              INT01110
      IF(NSETS.EQ.0) NSETS=IGRAD                                        INT01120
      REWIND ITAP11                                                     INT01130
      DO 106 II=1,NSETS                                                 INT01140
      READ(ITAP11,10) (I30(I),I=1,20)                                   INT01150
      READ(ITAP11,8) NATOM,ESCF                                         INT01160
      DO 104 I=1,NATOM                                                  INT01170
      READ(ITAP11,9) CHG(I),X(I),Y(I),Z(I)                              INT01180
  104 CONTINUE                                                          INT01190
      DO 105 I=1,NATOM                                                  INT01200
      READ(ITAP11,9) CHG(I),ENGX(I),ENGY(I),ENGZ(I)                     INT01210
  105 CONTINUE                                                          INT01220
  106 CONTINUE                                                          INT01230
      WRITE(6,3) (I30(I),I=1,20)                                        INT01240
      WRITE(6,4) NATOM,ESCF                                             INT01250
      NAT1=NATOM-1                                                      INT01260
      N3N=NATOM*3                                                       INT01270
      WRITE(6,5)                                                        INT01280
      DO 107 I=1,NATOM                                                  INT01290
      WRITE(6,6) I,CHG(I),X(I),Y(I),Z(I)                                INT01300
  107 CONTINUE                                                          INT01310
      CALL DIST(IPRNT)                                                  INT01320
      WRITE(6,7)                                                        INT01330
      DO 108 I=1,NATOM                                                  INT01340
      WRITE(6,6) I,CHG(I),ENGX(I),ENGY(I),ENGZ(I)                       INT01350
  108 CONTINUE                                                          INT01360
      GO TO 203                                                         INT01370
C                                                                       INT01380
  202 CONTINUE                                                          INT01390
      NSETS=MOD(IGRAD,100)                                              INT01400
      REWIND ITAP13                                                     INT01410
      NATOM=IATOM                                                       INT01420
      READ(ITAP13,2) NSTORE                                             INT01430
      DO 113 II=1,NSETS                                                 INT01440
      READ(ITAP13,8) NUM,ESCF                                           INT01450
      DO 110 I=1,NATOM                                                  INT01460
      READ(ITAP13,9) CHG(I),X(I),Y(I),Z(I)                              INT01470
  110 CONTINUE                                                          INT01480
      DO 111 I=1,NATOM                                                  INT01490
      READ(ITAP13,9) CHG(I),ENGX(I),ENGY(I),ENGZ(I)                     INT01500
  111 CONTINUE                                                          INT01510
      DO 112 I=1,ICOOD                                                  INT01520
      READ(ITAP13,9) RR,GR                                              INT01530
  112 CONTINUE                                                          INT01540
  113 CONTINUE                                                          INT01550
      WRITE(6,4) NATOM,ESCF                                             INT01560
      NAT1=NATOM-1                                                      INT01570
      N3N=NATOM*3                                                       INT01580
      WRITE(6,5)                                                        INT01590
      DO 114 I=1,NATOM                                                  INT01600
      WRITE(6,6) I,CHG(I),X(I),Y(I),Z(I)                                INT01610
  114 CONTINUE                                                          INT01620
      CALL DIST(IPRNT)                                                  INT01630
      WRITE(6,7)                                                        INT01640
      DO 115 I=1,NATOM                                                  INT01650
      WRITE(6,6) I,CHG(I),ENGX(I),ENGY(I),ENGZ(I)                       INT01660
  115 CONTINUE                                                          INT01670
C                                                                       INT01680
  203 CONTINUE                                                          INT01690
C   NINT IS NUMBER OF INTERNAL COORDINATES                              INT01700
      NINT=NST+NBND+NLIB+NDEF+NTORS                                     INT01710
      NMAX=MAX0(NINT,N3N)                                               INT01720
      NMAXSQ=NMAX*NMAX                                                  INT01730
      NATOM2=NATOM*NATOM                                                INT01740
      IC1=1                                                             INT01750
      IC2=IC1+NMAXSQ                                                    INT01760
      IC3=IC2+NMAXSQ                                                    INT01770
      IC4=IC3+NMAXSQ                                                    INT01780
      IC5=IC4+NMAXSQ                                                    INT01790
      IC6=IC5+NMAXSQ                                                    INT01800
      IC7=IC6+NMAX                                                      INT01810
      IA7=IC7*2-1                                                       INT01820
      IC8=IC7+NMAX                                                      INT01830
      IC9=IC8+NATOM2                                                    INT01840
      IC10=IC9+NATOM2                                                   INT01850
      ICMAX=IC10+NATOM2                                                 INT01860
C................B       BIN     AA      BB      FX      FF......       INT01870
      CALL VIBGF(CC(IC1),CC(IC2),CC(IC3),CC(IC4),CC(IC5),CC(IC5),       INT01880
C................PIVOT   INDEX   EX      EY      EZ............         INT01890
     1           CC(IC6),IA(IA7),CC(IC8),CC(IC9),CC(IC10),NMAX)         INT01900
C                                                                       INT01910
      CALL TSTOP(6)                                                     INT01920
      STOP                                                              INT01930
      END                                                               INT01940
      SUBROUTINE DIST(IPRNT)                                            INT01950
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT01960
      COMMON/VIB101/NATOM,NAT1,N3N,NINT                                 INT01970
      COMMON/VIB201/CHG(50),X(50),Y(50),Z(50)                           INT01980
      COMMON/VIB202/R(1275)                                             INT01990
      DATA A00 / 0.0D+00 /                                              INT02000
      DATA LDIM / 1275 /                                                INT02010
    1 FORMAT(//,2X,' INTERATOMIC DISTANCE MATRIX (IN BOHR)'/)
C                                                                       INT02030
      IJ=0                                                              INT02040
      DO 101 I=1,NATOM                                                  INT02050
      DO 101 J=1,I                                                      INT02060
      IJ=IJ+1                                                           INT02070
      R(IJ)=A00                                                         INT02080
      IF(I.EQ.J) GO TO 101                                              INT02090
      XD=X(I)-X(J)                                                      INT02100
      YD=Y(I)-Y(J)                                                      INT02110
      ZD=Z(I)-Z(J)                                                      INT02120
      R(IJ)=DSQRT(XD*XD+YD*YD+ZD*ZD)                                    INT02130
  101 CONTINUE                                                          INT02140
      IF(IPRNT.EQ.0) GO TO 201                                          INT02150
      WRITE(6,1)                                                        INT02160
      CALL PRINT(R,LDIM,NATOM,6)                                        INT02170
  201 RETURN                                                            INT02180
      END                                                               INT02190
      SUBROUTINE VIBGF(B,BIN,AA,BB,FX,FF,PIVOT,INDEX,EX,EY,EZ,NMAX)     INT02200
C   DIMENSION AND COMMON STATEMENTS                                     INT02210
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT02220
      DIMENSION B(NMAX,NMAX),BIN(NMAX,NMAX),AA(NMAX,NMAX),BB(NMAX,NMAX) INT02230
      DIMENSION FX(NMAX,NMAX),FF(NMAX,NMAX)                             INT02240
      DIMENSION PIVOT(NMAX),INDEX(NMAX*2)                               INT02250
      DIMENSION EX(NATOM,NATOM),EY(NATOM,NATOM),EZ(NATOM,NATOM)         INT02260
      COMMON/VIB101/NATOM,NAT1,N3N,NINT                                 INT02270
      COMMON/VIB102/ESCF                                                INT02280
      COMMON/VIB103/NST,NBND,NLIB,NDEF,NTORS                            INT02290
      COMMON/VIB104/IFORCE,IGRAD,IATOM,ICOOD,IPRNT                      INT02300
      COMMON/VIB201/CHG(50),X(50),Y(50),Z(50)                           INT02310
      COMMON/VIB202/R(1275)                                             INT02320
      COMMON/VIB203/ENGX(50),ENGY(50),ENGZ(50)                          INT02330
      COMMON/VIB204/RINT(150),ENGR(150)                                 INT02340
      COMMON/VIB205/IOFF(150)                                           INT02350
      COMMON/VIB301/KR(150),LR(150)                                     INT02360
      COMMON/VIB302/KA(150),LA(150),MA(150),AGLB(150)                   INT02370
      COMMON/VIB303/KB(100),LB(100),MB(100),AGLL(100)                   INT02380
      COMMON/VIB304/KD(100),LD(100),MD(100),ND(100)                     INT02390
      COMMON/VIB305/KT(100),LT(100),MT(100),NT(100)                     INT02400
      COMMON/VIB306/AGLT(100),AGLO(100),AGLR(100)                       INT02410
      COMMON/VIB307/AGLD(100),AGLE(100),AGLF(100)                       INT02420
      DATA A00,ONE / 0.0D+00 , 1.0D+00 /                                INT02440
      DATA PI / 3.1415926536D+00 /                                      INT02450
      data he,a0 / 4.359813653d+00 , 0.52917706d+00 /
C   FORMAT STATEMENTS                                                   INT02460
    1 FORMAT(//)                                                        INT02470
    2 FORMAT(10I5)                                                      INT02480
    3 FORMAT(2X,' PARAMETERS FOR THE CALCULATION'//                     INT02490
     * 2X,' NST    = ',I5/                                              INT02500
     * 2X,' NBND   = ',I5/                                              INT02510
     $ 2X,' NLIB   = ',I5/                                              INT02520
     * 2X,' NDEF   = ',I5/                                              INT02530
     * 2X,' NTORS  = ',I5/                                              INT02540
     * 2X,' IFORCE = ',I5/                                              INT02550
     * 2X,' IGRAD  = ',I5/                                              INT02560
     * 2X,' IATOM  = ',I5/                                              INT02570
     * 2X,' ICOOD  = ',I5/                                              INT02580
     * 2X,' IPRNT  = ',I5/                                              INT02590
     * 2X,' NINT   = ',I5/                                              INT02600
     * 2X,' NATOM  = ',I5/                                              INT02610
     * 2X,' N3N    = ',I5)                                              INT02620
    4 FORMAT(2I5)                                                       INT02630
    5 FORMAT(//,2X,21H INTERNAL COORDINATES,50X,11H COORDINATE/         INT02640
     1 73X,4H NO.)                                                      INT02650
    6 FORMAT(2X,' BOND LENGTH (IN BOHR)'/
     1 2X,20H   R   K-L    K    L/)                                     INT02670
    7 FORMAT(2X,F10.5,2I5,51X,2H (,I2,1H))                              INT02680
    8 FORMAT(3I5)                                                       INT02690
    9 FORMAT(/,2X,' BOND ANGLE (IN DEGREES)'/
     1 2X,25H ANG-KLM      K    L    M/)                                INT02710
   10 FORMAT(2X,F10.5,3I5,46X,2H (,I2,1H))                              INT02720
   11 FORMAT(/,2X,' LINEAR BOND ANGLE (IN DEGREES)'/
     1 2X,8H ANG-KLM,5X,2H K,3X,2H L,3X,2H M/)                          INT02740
   12 FORMAT(4I5)                                                       INT02750
   13 FORMAT(/,2X,' OUT OF PLANE ANGLE (IN DEGREES)'/ 
     1 2X,8H ANG-LNM,7X,8H ANG-KNM,7X,8H ANG-OOP,                       INT02770
     2 8X,2H K,3X,2H L,3X,2H M,3X,2H N/)                                INT02780
   14 FORMAT(2X,F10.5,5X,F10.5,5X,F10.5,3X,4I5,8X,2H (,I2,1H))          INT02790
   15 FORMAT(/,2X,' TORSIONAL ANGLE (IN DEGREES)'/ 
     1 2X,8H ANG-KLM,7X,8H ANG-LMN,7X,8H ANG-TOR,                       INT02810
     2 8X,2H K,3X,2H L,3X,2H M,3X,2H N/)                                INT02820
   16 FORMAT(//,2X,' EX MATRIX'/)                                       INT02830
   17 FORMAT(//,2X,' EY MATRIX'/)                                       INT02840
   18 FORMAT(//,2X,' EZ MATRIX'/)                                       INT02850
   19 FORMAT(//,2X,' B MATRIX'/)                                        INT02860
   20 FORMAT(//,2X,' B INVERSE MATRIX'/)                                INT02870
   21 FORMAT(//,2X,' DETERM OF B MATRIX = ',D15.7)                      INT02880
   22 FORMAT(//,2X,' B * B(-1) = E'/)                                   INT02890
   23 FORMAT(//,2X,' INTERNAL ENERGY GRADIENTS '/                       INT02900
     *2x,' (BOND LENGTH IN HARTREE/BOHR AND BOND ANGLE IN HARTREE/RADIAN
     *)'/5X,4H NO.,10X,7H COORD.,13X,6H GRAD./)    
   24 FORMAT(2X,I7,5X,2F20.10)                                          INT02920
   25 FORMAT(2I5)                                                       INT02930
   26 FORMAT(3F20.10)                                                   INT02940
   27 FORMAT(//,2X,' FX MATRIX (IN HARTREE/BOHR+2)'/)
   28 FORMAT(//,2X,' FX MATRIX (IN MDYN/A)'/)
   29 FORMAT(//,2X,' BIN * FX  MATRIX'/)                                INT02960
   30 FORMAT(//,2X,' F MATRIX'/)                                        INT02970
   31 FORMAT(2I5,F20.10,I5)                                             INT02980
   32 FORMAT(4F20.10)                                                   INT02990
   33 FORMAT(2F20.10)                                                   INT03000
   34 FORMAT(3F20.10)                                                   INT03010
C                                                                       INT03020
      INPUT=5                                                           INT03030
      itap12=12
      itap15=15
      CALL LOCATE(INPUT,'# INTCOS #',IERR)                              INT03040
      READ(5,2) ID1,ID2,ID3,ID4,ID5,ID6,ID7,ID8,ID9,ID10                INT03050
C                                                                       INT03060
      WRITE(6,1)                                                        INT03070
      AR=PI/180.0D+00                                                   INT03080
      funit=he/(a0*a0)
      WRITE(6,3) NST,NBND,NLIB,NDEF,NTORS,IFORCE,IGRAD,IATOM,ICOOD,     INT03090
     1           IPRNT,NINT,NATOM,N3N                                   INT03100
C                                                                       INT03110
C   THE CALCULATION OF B VECTORS                                        INT03120
C   THE CALCULATION OF UNIT VECTORS                                     INT03130
      DO 101 I=1,NATOM                                                  INT03140
      DO 101 J=1,NATOM                                                  INT03150
      EX(I,J)=A00                                                       INT03160
      EY(I,J)=A00                                                       INT03170
  101 EZ(I,J)=A00                                                       INT03180
      DO 102 I=1,NINT                                                   INT03190
      DO 102 J=1,N3N                                                    INT03200
  102 B(I,J)=A00                                                        INT03210
      DO 103 I=1,NAT1                                                   INT03220
      IR=I+1                                                            INT03230
      DO 103 J=IR,NATOM                                                 INT03240
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      INT03250
      EX(I,J)=-(X(I)-X(J))/R(IJ)                                        INT03260
      EY(I,J)=-(Y(I)-Y(J))/R(IJ)                                        INT03270
      EZ(I,J)=-(Z(I)-Z(J))/R(IJ)                                        INT03280
      EX(J,I)=-EX(I,J)                                                  INT03290
      EY(J,I)=-EY(I,J)                                                  INT03300
      EZ(J,I)=-EZ(I,J)                                                  INT03310
  103 CONTINUE                                                          INT03320
C***INPUT DATA FOR STRETCHINGS***                                       INT03330
      READ(5,4) (KR(I),LR(I),I=1,NST)                                   INT03340
      WRITE(6,5)                                                        INT03350
      WRITE(6,6)                                                        INT03360
      DO 104 I=1,NST                                                    INT03370
      K=KR(I)                                                           INT03380
      L=LR(I)                                                           INT03390
      II=I                                                              INT03400
      KL=IOFF(MAX0(K,L))+MIN0(K,L)                                      INT03410
      RINT(II)=R(KL)                                                    INT03420
      WRITE(6,7) R(KL),KR(I),LR(I),II                                   INT03430
  104 CONTINUE                                                          INT03440
      CALL STRNG(B,NMAX,EX,EY,EZ)                                       INT03450
      IF(NBND.EQ.0) GO TO 201                                           INT03460
C***INPUT DATA FOR BENDINGS***                                          INT03470
      READ(5,8) (KA(I),LA(I),MA(I),I=1,NBND)                            INT03480
      WRITE(6,9)                                                        INT03490
      DO 105 I=1,NBND                                                   INT03500
      K=KA(I)                                                           INT03510
      L=LA(I)                                                           INT03520
      M=MA(I)                                                           INT03530
      II=I+NST                                                          INT03540
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::INT03550
C   BOND ANGLE CALCULATION (WILSON-DECIOUS-CROSS  P.57, EQ. 8)          INT03560
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::INT03570
      AKLM=EX(L,K)*EX(L,M)+EY(L,K)*EY(L,M)+EZ(L,K)*EZ(L,M)              INT03580
      AGLB(I)=DACOS(AKLM)                                               INT03590
      ANG=AGLB(I)/AR                                                    INT03600
      RINT(II)=AGLB(I)                                                  INT03610
      WRITE(6,10) ANG,KA(I),LA(I),MA(I),II                              INT03620
  105 CONTINUE                                                          INT03630
      CALL BENDNG(B,NMAX,EX,EY,EZ)                                      INT03640
  201 IF(NLIB.EQ.0) GO TO 202                                           INT03650
C***INPUT DATA FOR LINEAR BENDINGS***                                   INT03660
      NLHF=NLIB/2                                                       INT03670
      READ(5,8) (KB(I),LB(I),MB(I),I=1,NLHF)                            INT03680
      WRITE(6,11)                                                       INT03690
      DO 106 I=1,NLHF                                                   INT03700
      I1=2*I-1+NST+NBND                                                 INT03710
      I2=I1+1                                                           INT03720
      K=KB(I)                                                           INT03730
      L=LB(I)                                                           INT03740
      M=MB(I)                                                           INT03750
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::INT03760
C   BOND ANGLE CALCULATION (WILSON-DECIOUS-CROSS  P.57, EQ. 8)          INT03770
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::INT03780
      AKLM=EX(L,K)*EX(L,M)+EY(L,K)*EY(L,M)+EZ(L,K)*EZ(L,M)              INT03790
      IF(DABS(AKLM).GT.ONE) AKLM=DSIGN(ONE,AKLM)                        INT03800
      AGLL(I)=DACOS(AKLM)                                               INT03810
      ANG=AGLL(I)/AR                                                    INT03820
      RINT(I1)=AGLL(I)                                                  INT03830
      RINT(I2)=AGLL(I)                                                  INT03840
      WRITE(6,10) ANG,KB(I),LB(I),MB(I),I1                              INT03850
      WRITE(6,10) ANG,KB(I),LB(I),MB(I),I2                              INT03860
  106 CONTINUE                                                          INT03870
      CALL LINBND(B,NMAX,EX,EY,EZ)                                      INT03880
  202 IF(NDEF.EQ.0) GO TO 203                                           INT03890
C***INPUT DATA FOR DEFORMATIONS***                                      INT03900
      READ(5,12) (KD(I),LD(I),MD(I),ND(I),I=1,NDEF)                     INT03910
      WRITE(6,13)                                                       INT03920
      DO 107 I=1,NDEF                                                   INT03930
      II=I+NST+NBND+NLIB                                                INT03940
      K=KD(I)                                                           INT03950
      L=LD(I)                                                           INT03960
      M=MD(I)                                                           INT03970
      N=ND(I)                                                           INT03980
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::INT03990
C   CALCULATE OUT-OF-PLANE ANGLE (WILSON-DECIOUS-CROSS   P.59 EQ. 18)   INT04000
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::INT04010
      ALNM=EX(N,L)*EX(N,M)+EY(N,L)*EY(N,M)+EZ(N,L)*EZ(N,M)              INT04020
      AGLD(I)=DACOS(ALNM)                                               INT04030
      AKNM=EX(N,K)*EX(N,M)+EY(N,K)*EY(N,M)+EZ(N,K)*EZ(N,M)              INT04040
      AGLE(I)=DACOS(AKNM)                                               INT04050
      AKNL=EX(N,K)*EX(N,L)+EY(N,K)*EY(N,L)+EZ(N,K)*EZ(N,L)              INT04060
      AGLF(I)=DACOS(AKNL)                                               INT04070
      DLM=AGLD(I)/AR                                                    INT04080
      DKM=AGLE(I)/AR                                                    INT04090
      DKL=AGLF(I)/AR                                                    INT04100
      RINT(II)=AGLD(I)                                                  INT04110
      PX=(EY(N,L)*EZ(N,M)-EZ(N,L)*EY(N,M))                              INT04120
      PY=(EZ(N,L)*EX(N,M)-EX(N,L)*EZ(N,M))                              INT04130
      PZ=(EX(N,L)*EY(N,M)-EY(N,L)*EX(N,M))                              INT04140
      PROD=PX*EX(N,K)+PY*EY(N,K)+PZ*EZ(N,K)                             INT04150
      SINP=DSQRT(ONE-ALNM*ALNM)                                         INT04160
      AKLMN=PROD/SINP                                                   INT04170
      AOP=DASIN(AKLMN)/AR                                               INT04180
      WRITE(6,14) DLM,DKM,AOP,KD(I),LD(I),MD(I),ND(I),II                INT04190
  107 CONTINUE                                                          INT04200
      CALL DEFORM(B,NMAX,EX,EY,EZ)                                      INT04210
  203 IF(NTORS.EQ.0) GO TO 204                                          INT04220
C***INPUT DATA FOR TORSIONS***                                          INT04230
      READ(5,12) (KT(I),LT(I),MT(I),NT(I),I=1,NTORS)                    INT04240
      WRITE(6,15)                                                       INT04250
      DO 108 I=1,NTORS                                                  INT04260
      II=I+NST+NBND+NLIB+NDEF                                           INT04270
      K=KT(I)                                                           INT04280
      L=LT(I)                                                           INT04290
      M=MT(I)                                                           INT04300
      N=NT(I)                                                           INT04310
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::INT04320
C   CALCULATE TORSIONAL ANGLE (WILSON-DECIOUS-CROSS   P.60 EQ. 20)      INT04330
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::INT04340
      AKLM=EX(L,K)*EX(L,M)+EY(L,K)*EY(L,M)+EZ(L,K)*EZ(L,M)              INT04350
      AGLT(I)=DACOS(AKLM)                                               INT04360
      ALMN=EX(M,L)*EX(M,N)+EY(M,L)*EY(M,N)+EZ(M,L)*EZ(M,N)              INT04370
      AGLO(I)=DACOS(ALMN)                                               INT04380
      PX=(EY(K,L)*EZ(L,M)-EZ(K,L)*EY(L,M))                              INT04390
      PY=(EZ(K,L)*EX(L,M)-EX(K,L)*EZ(L,M))                              INT04400
      PZ=(EX(K,L)*EY(L,M)-EY(K,L)*EX(L,M))                              INT04410
      QX=(EY(L,M)*EZ(M,N)-EZ(L,M)*EY(M,N))                              INT04420
      QY=(EZ(L,M)*EX(M,N)-EX(L,M)*EZ(M,N))                              INT04430
      QZ=(EX(L,M)*EY(M,N)-EY(L,M)*EX(M,N))                              INT04440
      SINPQ=DSQRT(ONE-AKLM*AKLM)*DSQRT(ONE-ALMN*ALMN)                   INT04450
      AKLMN=(PX*QX+PY*QY+PZ*QZ)/SINPQ                                   INT04460
CCC   PX=(EY(L,K)*EZ(L,M)-EZ(L,K)*EY(L,M))                              INT04470
CCC   PY=(EZ(L,K)*EX(L,M)-EX(L,K)*EZ(L,M))                              INT04480
CCC   PZ=(EX(L,K)*EY(L,M)-EY(L,K)*EX(L,M))                              INT04490
CCC   QX=(EY(M,L)*EZ(M,N)-EZ(M,L)*EY(M,N))                              INT04500
CCC   QY=(EZ(M,L)*EX(M,N)-EX(M,L)*EZ(M,N))                              INT04510
CCC   QZ=(EX(M,L)*EY(M,N)-EY(M,L)*EX(M,N))                              INT04520
CCC   SINPQ=DSQRT(ONE-AKLM*AKLM)*DSQRT(ONE-ALMN*ALMN)                   INT04530
CCC   AKLMN=(PX*QX+PY*QY+PZ*QZ)/SINPQ                                   INT04540
      IF(DABS(AKLMN).GT.ONE) AKLMN=DSIGN(ONE,AKLMN)                     INT04550
      AGLR(I)=DACOS(AKLMN)                                              INT04560
      ATL=AGLT(I)/AR                                                    INT04570
      ATM=AGLO(I)/AR                                                    INT04580
      ATT=AGLR(I)/AR                                                    INT04590
      RINT(II)=AGLR(I)                                                  INT04600
      WRITE(6,14) ATL,ATM,ATT,KT(I),LT(I),MT(I),NT(I),II                INT04610
  108 CONTINUE                                                          INT04620
      CALL TORS(B,NMAX,EX,EY,EZ)                                        INT04630
C                                                                       INT04640
C   PRINT OUT E VECTORS                                                 INT04650
  204 IF(IPRNT.LE.2) GO TO 301                                          INT04660
      WRITE(6,16)                                                       INT04670
      CALL MATOUT(EX,NATOM,NATOM,NATOM,NATOM,6)                         INT04680
      WRITE(6,17)                                                       INT04690
      CALL MATOUT(EY,NATOM,NATOM,NATOM,NATOM,6)                         INT04700
      WRITE(6,18)                                                       INT04710
      CALL MATOUT(EZ,NATOM,NATOM,NATOM,NATOM,6)                         INT04720
C                                                                       INT04730
  301 IF(IPRNT.LE.2) GO TO 302                                          INT04740
      WRITE(6,19)                                                       INT04750
      CALL MATOUT(B,NMAX,NMAX,NINT,N3N,6)                               INT04760
C   THE CALCULATION OF INVERSE B MATRIX                                 INT04770
  302 DO 111 I=1,NINT                                                   INT04780
      DO 111 J=1,NINT                                                   INT04790
      BB(I,J)=A00                                                       INT04800
      DO 111 K=1,N3N                                                    INT04810
      BB(I,J)=BB(I,J)+B(I,K)*B(J,K)                                     INT04820
  111 CONTINUE                                                          INT04830
      DO 112 I=1,NINT                                                   INT04840
      DO 112 J=1,NINT                                                   INT04850
  112 BIN(I,J)=BB(I,J)                                                  INT04860
      CALL MATINV(BIN,AA,PIVOT,INDEX,NINT,NINT,DETERM,NMAX)             INT04870
      IF(IPRNT.LE.2) GO TO 303                                          INT04880
      WRITE(6,20)                                                       INT04890
      CALL MATOUT(BIN,NMAX,NMAX,NINT,NINT,6)                            INT04900
  303 WRITE(6,21) DETERM                                                INT04910
      DO 113 I=1,NINT                                                   INT04920
      DO 113 J=1,NINT                                                   INT04930
      AA(I,J)=A00                                                       INT04940
      DO 113 K=1,NINT                                                   INT04950
      AA(I,J)=AA(I,J)+BIN(I,K)*BB(K,J)                                  INT04960
  113 CONTINUE                                                          INT04970
      IF(IPRNT.LE.2) GO TO 304                                          INT04980
      WRITE(6,22)                                                       INT04990
      CALL MATOUT(AA,NMAX,NMAX,NINT,NINT,6)                             INT05000
  304 DO 114 I=1,NINT                                                   INT05010
      DO 114 J=1,N3N                                                    INT05020
      BB(I,J)=A00                                                       INT05030
      DO 114 K=1,NINT                                                   INT05040
      BB(I,J)=BB(I,J)+BIN(I,K)*B(K,J)                                   INT05050
  114 CONTINUE                                                          INT05060
ccc   write(6,*) ' bb'
ccc   call matout(bb,nmax,nmax,nint,n3n,6)
C                                                                       INT05070
C***********************************************************            INT05080
C***TRANSFORM GRADIENTS FROM CARTESIAN TO INTERNAL SYSTEM***            INT05090
C***********************************************************            INT05100
      DO 116 I=1,NINT                                                   INT05110
      ENGR(I)=A00                                                       INT05120
      DO 115 J=1,NATOM                                                  INT05130
      JX=3*J-2                                                          INT05140
      JY=3*J-1                                                          INT05150
      JZ=3*J                                                            INT05160
      ENGR(I)=ENGR(I)+ENGX(J)*BB(I,JX)+ENGY(J)*BB(I,JY)                 INT05170
     1 +ENGZ(J)*BB(I,JZ)                                                INT05180
  115 CONTINUE                                                          INT05190
  116 CONTINUE                                                          INT05200
      WRITE(6,23)                                                       INT05210
      DO 117 I=1,NINT                                                   INT05220
      WRITE(6,24) I,RINT(I),ENGR(I)                                     INT05230
  117 CONTINUE                                                          INT05240
C                                                                       INT05250
C***********************************************************************INT05260
C***TRANSFORM FORCE CONSTANT MATRIX FROM CARTESIAN TO INTERNAL SYSTEM***INT05270
C***********************************************************************INT05280
      IF(IFORCE.EQ.0) GO TO 310                                         INT05290
C   READ IN FORCE CONSTANT MATRIX                                       INT05300
      REWIND ITAP15                                                     INT05310
      READ(ITAP15,25) NATOM,NNC                                         INT05320
      READ(ITAP15,26) ((FX(I,J),J=1,N3N),I=1,N3N)                       INT05330
      if(iprn.gt.3) then
       WRITE(6,27)                        
       CALL MATOUT(FX,NMAX,NMAX,N3N,N3N,6)     
      end if
      do 118 i=1,n3n
      do 118 j=1,n3n
      fx(i,j)=fx(i,j)*funit
  118 continue
      if(iprnt.gt.3) then
       write(6,28)
       call matout(fx,nmax,nmax,n3n,n3n,6)
      end if
c
C   TRANSFORM FX INTO F                                                 INT05370
      DO 121 I=1,NINT                                                   INT05380
      DO 121 J=1,N3N                                                    INT05390
      AA(I,J)=A00                                                       INT05400
      DO 120 K=1,N3N                                                    INT05410
      AA(I,J)=AA(I,J)+BB(I,K)*FX(K,J)                                   INT05420
  120 CONTINUE                                                          INT05430
  121 CONTINUE                                                          INT05440
      if(iprnt.gt.3) then
        WRITE(6,29)     
        CALL MATOUT(AA,NMAX,NMAX,NINT,N3N,6)    
      end if
      DO 123 I=1,NINT                                                   INT05480
      DO 123 J=1,NINT                                                   INT05490
      FF(I,J)=A00                                                       INT05500
      DO 122 K=1,N3N                                                    INT05510
      FF(I,J)=FF(I,J)+AA(I,K)*BB(J,K)                                   INT05520
  122 CONTINUE                                                          INT05530
  123 CONTINUE                                                          INT05540
      if(iprnt.gt.2) then
        WRITE(6,30)    
        CALL MATOUT(FF,NMAX,NMAX,NINT,NINT,6)    
      end if
      REWIND ITAP15                                                     INT05580
C                                                                       INT05590
C   SEND INFORMATION TO GEOMETRY OPTIMIZATION PROGRAM WITH TAPE12       INT05600
  310 CONTINUE                                                          INT05610
      REWIND ITAP12                                                     INT05620
      WRITE(ITAP12,31) NATOM,NINT,ESCF,IFORCE                           INT05630
      DO 125 I=1,NATOM                                                  INT05640
      WRITE(ITAP12,32) CHG(I),X(I),Y(I),Z(I)                            INT05650
  125 CONTINUE                                                          INT05660
      DO 126 I=1,NATOM                                                  INT05670
      WRITE(ITAP12,32) CHG(I),ENGX(I),ENGY(I),ENGZ(I)                   INT05680
  126 CONTINUE                                                          INT05690
      DO 127 I=1,NINT                                                   INT05700
      WRITE(ITAP12,33) RINT(I),ENGR(I)                                  INT05710
  127 CONTINUE                                                          INT05720
      if(iforce.ne.0) then
        WRITE(ITAP12,34) ((FF(I,J),J=1,NINT),I=1,NINT)   
      end if
      REWIND ITAP12                                                     INT05760
C                                                                       INT05770
      RETURN                                                            INT05780
      END                                                               INT05790
      SUBROUTINE STRNG(B,NMAX,EX,EY,EZ)                                 INT05800
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT05810
      DIMENSION B(NMAX,NMAX)                                            INT05820
      DIMENSION EX(NATOM,NATOM),EY(NATOM,NATOM),EZ(NATOM,NATOM)         INT05830
      COMMON/VIB101/NATOM,NAT1,N3N,NINT                                 INT05840
      COMMON/VIB103/NST,NBND,NLIB,NDEF,NTORS                            INT05850
      COMMON/VIB301/KR(150),LR(150)                                     INT05860
C                                                                       INT05870
C   THE CALCULATION OF S VECTORS FOR STRETCHINGS                        INT05880
      DO 101 I=1,NST                                                    INT05890
      K=KR(I)                                                           INT05900
      L=LR(I)                                                           INT05910
      KX=3*K-2                                                          INT05920
      KY=3*K-1                                                          INT05930
      KZ=3*K                                                            INT05940
      LX=3*L-2                                                          INT05950
      LY=3*L-1                                                          INT05960
      LZ=3*L                                                            INT05970
      B(I,KX)=EX(L,K)                                                   INT05980
      B(I,KY)=EY(L,K)                                                   INT05990
      B(I,KZ)=EZ(L,K)                                                   INT06000
      B(I,LX)=EX(K,L)                                                   INT06010
      B(I,LY)=EY(K,L)                                                   INT06020
  101 B(I,LZ)=EZ(K,L)                                                   INT06030
      RETURN                                                            INT06040
      END                                                               INT06050
      SUBROUTINE BENDNG(B,NMAX,EX,EY,EZ)                                INT06060
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT06070
      DIMENSION B(NMAX,NMAX)                                            INT06080
      DIMENSION EX(NATOM,NATOM),EY(NATOM,NATOM),EZ(NATOM,NATOM)         INT06090
      COMMON/VIB101/NATOM,NAT1,N3N,NINT                                 INT06100
      COMMON/VIB103/NST,NBND,NLIB,NDEF,NTORS                            INT06110
      COMMON/VIB202/R(1275)                                             INT06120
      COMMON/VIB205/IOFF(150)                                           INT06130
      COMMON/VIB302/KA(150),LA(150),MA(150),AGLB(150)                   INT06140
C                                                                       INT06150
C   THE CALCULATION OF S VECTORS FOR ANGLE BENDINGS                     INT06160
      DO 101 I=1,NBND                                                   INT06170
      II=I+NST                                                          INT06180
      K=KA(I)                                                           INT06190
      L=LA(I)                                                           INT06200
      M=MA(I)                                                           INT06210
      KX=3*K-2                                                          INT06220
      KY=3*K-1                                                          INT06230
      KZ=3*K                                                            INT06240
      LX=3*L-2                                                          INT06250
      LY=3*L-1                                                          INT06260
      LZ=3*L                                                            INT06270
      MX=3*M-2                                                          INT06280
      MY=3*M-1                                                          INT06290
      MZ=3*M                                                            INT06300
      AG=AGLB(I)                                                        INT06310
      LK=IOFF(MAX0(L,K))+MIN0(L,K)                                      INT06320
      LM=IOFF(MAX0(L,M))+MIN0(L,M)                                      INT06330
      B(II,KX)=(DCOS(AG)*EX(L,K)-EX(L,M))/(R(LK)*DSIN(AG))              INT06340
      B(II,KY)=(DCOS(AG)*EY(L,K)-EY(L,M))/(R(LK)*DSIN(AG))              INT06350
      B(II,KZ)=(DCOS(AG)*EZ(L,K)-EZ(L,M))/(R(LK)*DSIN(AG))              INT06360
      B(II,MX)=(DCOS(AG)*EX(L,M)-EX(L,K))/(R(LM)*DSIN(AG))              INT06370
      B(II,MY)=(DCOS(AG)*EY(L,M)-EY(L,K))/(R(LM)*DSIN(AG))              INT06380
      B(II,MZ)=(DCOS(AG)*EZ(L,M)-EZ(L,K))/(R(LM)*DSIN(AG))              INT06390
      B(II,LX)=-B(II,KX)-B(II,MX)                                       INT06400
      B(II,LY)=-B(II,KY)-B(II,MY)                                       INT06410
  101 B(II,LZ)=-B(II,KZ)-B(II,MZ)                                       INT06420
      RETURN                                                            INT06430
      END                                                               INT06440
      SUBROUTINE LINBND(B,NMAX,EX,EY,EZ)                                INT06450
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT06460
      DIMENSION B(NMAX,NMAX)                                            INT06470
      DIMENSION EX(NATOM,NATOM),EY(NATOM,NATOM),EZ(NATOM,NATOM)         INT06480
      COMMON/VIB101/NATOM,NAT1,N3N,NINT                                 INT06490
      COMMON/VIB103/NST,NBND,NLIB,NDEF,NTORS                            INT06500
      COMMON/VIB201/CHG(50),X(50),Y(50),Z(50)                           INT06510
      COMMON/VIB202/R(1275)                                             INT06520
      COMMON/VIB205/IOFF(150)                                           INT06530
      COMMON/VIB303/KB(100),LB(100),MB(100),AGLL(100)                   INT06540
      DATA ONE,TWO,THREE / 1.0D+00 , 2.0D+00 , 3.0D+00 /                INT06550
      DATA PI / 3.1415926536D+00 /                                      INT06560
C                                                                       INT06570
                                                                        INT06580
      NLHF=NLIB/2                                                       INT06590
C   N IS A DUMMY ATOM                                                   INT06600
      DO 101 I=1,NLHF                                                   INT06610
      I1=2*I-1+NST+NBND                                                 INT06620
      I2=I1+1                                                           INT06630
      K=KB(I)                                                           INT06640
      L=LB(I)                                                           INT06650
      M=MB(I)                                                           INT06660
      KX=3*K-2                                                          INT06670
      KY=3*K-1                                                          INT06680
      KZ=3*K                                                            INT06690
      LX=3*L-2                                                          INT06700
      LY=3*L-1                                                          INT06710
      LZ=3*L                                                            INT06720
      MX=3*M-2                                                          INT06730
      MY=3*M-1                                                          INT06740
      MZ=3*M                                                            INT06750
      XN=X(M)+ONE                                                       INT06760
      YN=Y(M)+ONE                                                       INT06770
      ZN=Z(M)+ONE                                                       INT06780
      RMN=DSQRT(THREE)                                                  INT06790
      ML=IOFF(MAX0(M,L))+MIN0(M,L)                                      INT06800
      RML=R(ML)                                                         INT06810
      KL=IOFF(MAX0(K,L))+MIN0(K,L)                                      INT06820
      RNL=DSQRT((XN-X(L))**2+(YN-Y(L))**2+(ZN-Z(L))**2)                 INT06830
      SQCOS=(RNL**2-RMN**2-RML**2)/(TWO*RML*RMN)                        INT06840
      AGM=PI-DACOS(SQCOS)                                               INT06850
      EXMN=ONE/RMN                                                      INT06860
      EYMN=ONE/RMN                                                      INT06870
      EZMN=ONE/RMN                                                      INT06880
      B(I1,KX)=(DCOS(AGM)*EX(M,L)-EXMN)/(R(KL)*DSIN(AGM))               INT06890
      B(I1,KY)=(DCOS(AGM)*EY(M,L)-EYMN)/(R(KL)*DSIN(AGM))               INT06900
      B(I1,KZ)=(DCOS(AGM)*EZ(M,L)-EZMN)/(R(KL)*DSIN(AGM))               INT06910
      B(I1,MX)=(DCOS(AGM)*EX(M,L)-EXMN)/(R(ML)*DSIN(AGM))               INT06920
      B(I1,MY)=(DCOS(AGM)*EY(M,L)-EYMN)/(R(ML)*DSIN(AGM))               INT06930
      B(I1,MZ)=(DCOS(AGM)*EZ(M,L)-EZMN)/(R(ML)*DSIN(AGM))               INT06940
      B(I1,LX)=-B(I1,KX)-B(I1,MX)                                       INT06950
      B(I1,LY)=-B(I1,KY)-B(I1,MY)                                       INT06960
      B(I1,LZ)=-B(I1,KZ)-B(I1,MZ)                                       INT06970
      B(I2,KX)=-(EYMN*EZ(M,L)-EZMN*EY(M,L))/(R(KL)*DSIN(AGM))           INT06980
      B(I2,KY)=-(EZMN*EX(M,L)-EXMN*EZ(M,L))/(R(KL)*DSIN(AGM))           INT06990
      B(I2,KZ)=-(EXMN*EY(M,L)-EYMN*EX(M,L))/(R(KL)*DSIN(AGM))           INT07000
      B(I2,MX)=-(EYMN*EZ(M,L)-EZMN*EY(M,L))/(R(ML)*DSIN(AGM))           INT07010
      B(I2,MY)=-(EZMN*EX(M,L)-EXMN*EZ(M,L))/(R(ML)*DSIN(AGM))           INT07020
      B(I2,MZ)=-(EXMN*EY(M,L)-EYMN*EX(M,L))/(R(ML)*DSIN(AGM))           INT07030
      B(I2,LX)=-B(I2,KX)-B(I2,MX)                                       INT07040
      B(I2,LY)=-B(I2,KY)-B(I2,MY)                                       INT07050
      B(I2,LZ)=-B(I2,KZ)-B(I2,MZ)                                       INT07060
  101 CONTINUE                                                          INT07070
      RETURN                                                            INT07080
      END                                                               INT07090
      SUBROUTINE DEFORM(B,NMAX,EX,EY,EZ)                                INT07100
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT07110
      DIMENSION B(NMAX,NMAX)                                            INT07120
      DIMENSION EX(NATOM,NATOM),EY(NATOM,NATOM),EZ(NATOM,NATOM)         INT07130
      COMMON/VIB101/NATOM,NAT1,N3N,NINT                                 INT07140
      COMMON/VIB103/NST,NBND,NLIB,NDEF,NTORS                            INT07150
      COMMON/VIB202/R(1275)                                             INT07160
      COMMON/VIB205/IOFF(150)                                           INT07170
      COMMON/VIB304/KD(100),LD(100),MD(100),ND(100)                     INT07180
C                                                                       INT07190
C   THE CALCULATION OF S VECTORS FOR OUT OF PLANE DEFORMATIONS          INT07200
      DO 101 I=1,NDEF                                                   INT07210
      II=I+NST+NBND+NLIB                                                INT07220
      K=KD(I)                                                           INT07230
      L=LD(I)                                                           INT07240
      M=MD(I)                                                           INT07250
      N=ND(I)                                                           INT07260
      KX=3*K-2                                                          INT07270
      KY=3*K-1                                                          INT07280
      KZ=3*K                                                            INT07290
      LX=3*L-2                                                          INT07300
      LY=3*L-1                                                          INT07310
      LZ=3*L                                                            INT07320
      MX=3*M-2                                                          INT07330
      MY=3*M-1                                                          INT07340
      MZ=3*M                                                            INT07350
      NX=3*N-2                                                          INT07360
      NY=3*N-1                                                          INT07370
      NZ=3*N                                                            INT07380
      NK=IOFF(MAX0(N,K))+MIN0(N,K)                                      INT07390
      NL=IOFF(MAX0(N,L))+MIN0(N,L)                                      INT07400
      NM=IOFF(MAX0(N,M))+MIN0(N,M)                                      INT07410
      CEXLM=EY(N,L)*EZ(N,M)-EZ(N,L)*EY(N,M)                             INT07420
      CEYLM=EZ(N,L)*EX(N,M)-EX(N,L)*EZ(N,M)                             INT07430
      CEZLM=EX(N,L)*EY(N,M)-EY(N,L)*EX(N,M)                             INT07440
      CEXMK=EY(N,M)*EZ(N,K)-EZ(N,M)*EY(N,K)                             INT07450
      CEYMK=EZ(N,M)*EX(N,K)-EX(N,M)*EZ(N,K)                             INT07460
      CEZMK=EX(N,M)*EY(N,K)-EY(N,M)*EX(N,K)                             INT07470
      CEXKL=EY(N,K)*EZ(N,L)-EZ(N,K)*EY(N,L)                             INT07480
      CEYKL=EZ(N,K)*EX(N,L)-EX(N,K)*EZ(N,L)                             INT07490
      CEZKL=EX(N,K)*EY(N,L)-EY(N,K)*EX(N,L)                             INT07500
      B(II,KX)=CEXLM/R(NK)                                              INT07510
      B(II,KY)=CEYLM/R(NK)                                              INT07520
      B(II,KZ)=CEZLM/R(NK)                                              INT07530
      B(II,LX)=CEXMK/R(NL)                                              INT07540
      B(II,LY)=CEYMK/R(NL)                                              INT07550
      B(II,LZ)=CEZMK/R(NL)                                              INT07560
      B(II,MX)=CEXKL/R(NM)                                              INT07570
      B(II,MY)=CEYKL/R(NM)                                              INT07580
      B(II,MZ)=CEZKL/R(NM)                                              INT07590
      B(II,NX)=-B(II,KX)-B(II,LX)-B(II,MX)                              INT07600
      B(II,NY)=-B(II,KY)-B(II,LY)-B(II,MY)                              INT07610
  101 B(II,NZ)=-B(II,KZ)-B(II,LZ)-B(II,MZ)                              INT07620
      RETURN                                                            INT07630
      END                                                               INT07640
      SUBROUTINE TORS(B,NMAX,EX,EY,EZ)                                  INT07650
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT07660
      DIMENSION B(NMAX,NMAX)                                            INT07670
      DIMENSION EX(NATOM,NATOM),EY(NATOM,NATOM),EZ(NATOM,NATOM)         INT07680
      COMMON/VIB101/NATOM,NAT1,N3N,NINT                                 INT07690
      COMMON/VIB103/NST,NBND,NLIB,NDEF,NTORS                            INT07700
      COMMON/VIB202/R(1275)                                             INT07710
      COMMON/VIB205/IOFF(150)                                           INT07720
      COMMON/VIB305/KT(100),LT(100),MT(100),NT(100)                     INT07730
      COMMON/VIB306/AGLT(100),AGLO(100),AGLR(100)                       INT07740
C                                                                       INT07750
C   THE CALCULATION OF S VECTORS FOR TORSIONS                           INT07760
      DO 101 I=1,NTORS                                                  INT07770
      II=I+NST+NBND+NLIB+NDEF                                           INT07780
      K=KT(I)                                                           INT07790
      L=LT(I)                                                           INT07800
      M=MT(I)                                                           INT07810
      N=NT(I)                                                           INT07820
      KX=3*K-2                                                          INT07830
      KY=3*K-1                                                          INT07840
      KZ=3*K                                                            INT07850
      LX=3*L-2                                                          INT07860
      LY=3*L-1                                                          INT07870
      LZ=3*L                                                            INT07880
      MX=3*M-2                                                          INT07890
      MY=3*M-1                                                          INT07900
      MZ=3*M                                                            INT07910
      NX=3*N-2                                                          INT07920
      NY=3*N-1                                                          INT07930
      NZ=3*N                                                            INT07940
      KL=IOFF(MAX0(K,L))+MIN0(K,L)                                      INT07950
      MN=IOFF(MAX0(M,N))+MIN0(M,N)                                      INT07960
      ML=IOFF(MAX0(M,L))+MIN0(M,L)                                      INT07970
      AL=AGLT(I)                                                        INT07980
      AM=AGLO(I)                                                        INT07990
CCC   CEXKM=EY(L,K)*EZ(L,M)-EZ(L,K)*EY(L,M)                             INT08000
CCC   CEYKM=EZ(L,K)*EX(L,M)-EX(L,K)*EZ(L,M)                             INT08010
CCC   CEZKM=EX(L,K)*EY(L,M)-EY(L,K)*EX(L,M)                             INT08020
CCC   CEXNL=EY(M,N)*EZ(M,L)-EZ(M,N)*EY(M,L)                             INT08030
CCC   CEYNL=EZ(M,N)*EX(M,L)-EX(M,N)*EZ(M,L)                             INT08040
CCC   CEZNL=EX(M,N)*EY(M,L)-EY(M,N)*EX(M,L)                             INT08050
      CEXKM=EY(K,L)*EZ(L,M)-EZ(K,L)*EY(L,M)                             INT08060
      CEYKM=EZ(K,L)*EX(L,M)-EX(K,L)*EZ(L,M)                             INT08070
      CEZKM=EX(K,L)*EY(L,M)-EY(K,L)*EX(L,M)                             INT08080
      CEXNL=EY(N,M)*EZ(M,L)-EZ(N,M)*EY(M,L)                             INT08090
      CEYNL=EZ(N,M)*EX(M,L)-EX(N,M)*EZ(M,L)                             INT08100
      CEZNL=EX(N,M)*EY(M,L)-EY(N,M)*EX(M,L)                             INT08110
      RS2K=R(KL)*DSIN(AL)**2                                            INT08120
      RS2M=R(ML)*DSIN(AM)**2                                            INT08130
      RS2N=R(MN)*DSIN(AM)**2                                            INT08140
      RS2L=R(ML)*DSIN(AL)**2                                            INT08150
      B(II,KX)=-CEXKM/RS2K                                              INT08160
      B(II,KY)=-CEYKM/RS2K                                              INT08170
      B(II,KZ)=-CEZKM/RS2K                                              INT08180
      B(II,LX)=(R(ML)-R(KL)*DCOS(AL))*CEXKM/(R(ML)*RS2K)+               INT08190
     1 DCOS(AM)*CEXNL/RS2M                                              INT08200
      B(II,LY)=(R(ML)-R(KL)*DCOS(AL))*CEYKM/(R(ML)*RS2K)+               INT08210
     1 DCOS(AM)*CEYNL/RS2M                                              INT08220
      B(II,LZ)=(R(ML)-R(KL)*DCOS(AL))*CEZKM/(R(ML)*RS2K)+               INT08230
     1 DCOS(AM)*CEZNL/RS2M                                              INT08240
      B(II,MX)=(R(ML)-R(MN)*DCOS(AM))*CEXNL/(R(ML)*RS2N)+               INT08250
     1 DCOS(AL)*CEXKM/RS2L                                              INT08260
      B(II,MY)=(R(ML)-R(MN)*DCOS(AM))*CEYNL/(R(ML)*RS2N)+               INT08270
     1 DCOS(AL)*CEYKM/RS2L                                              INT08280
      B(II,MZ)=(R(ML)-R(MN)*DCOS(AM))*CEZNL/(R(ML)*RS2N)+               INT08290
     1 DCOS(AL)*CEZKM/RS2L                                              INT08300
      B(II,NX)=-CEXNL/RS2N                                              INT08310
      B(II,NY)=-CEYNL/RS2N                                              INT08320
  101 B(II,NZ)=-CEZNL/RS2N                                              INT08330
      RETURN                                                            INT08340
      END                                                               INT08350
