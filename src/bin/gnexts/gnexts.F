      subroutine fentry(cc,ia,maxcor)
cc    PROGRAM GNEXTS                                                    GNE00010
C**********************************************************************
C*  Last updated on November 26, 1997 by Dr. Yukio Yamaguchi          *
C*  Modification for the Power series                                 *
C**********************************************************************
C*  Last updated on October 04, 1991 by Dr. Yukio Yamaguchi           *
C*  Modification for the UNIX system                                  *
C********************************************************************** GNE00020
C*  LAST UPDATED ON DECEMBER 05, 1990 BY YUKIO YAMAGUCHI              * GNE00030
C********************************************************************** GNE00040
C   THE SOURCE CODE OF GEOMETRY OPTIMIZATION PROGRAM                    GNE00050
C   BASED ON THE STEEPEST DESCENT METHOD                                GNE00060
C   THIS PROGRAM IS ABLE TO TREAT UP TO '50' ATOMS                      GNE00070
C   THIS PROGRAM IS ABLE TO TREAT UP TO '150' INTERNAL COORDINATES      GNE00080
      IMPLICIT REAL*8 (A-H,O-Z)                                         GNE00090
      dimension cc(maxcor),ia(maxcor*2)
cyy   DIMENSION CC(360000),IA(1)                                        GNE00100
      DIMENSION I30(200),NUNQ(50),NAT(50,50),COORD(3,50)                GNE00110
      CHARACTER*8 UPD                                                   GNE00120
      LOGICAL UPDATE                                                    GNE00130
      COMMON/COM101/NATOM,N3N,NINT                                      GNE00140
      COMMON/COM102/NSORT,NCOOD,NDIF,NPLUS,NCASE                        GNE00150
      COMMON/COM103/CHG(50),XX(50),YY(50),ZZ(50)                        GNE00160
      COMMON/COM104/EX(50),EY(50),EZ(50)                                GNE00170
      COMMON/COM105/ENGR(150),RCOOD(150)                                GNE00180
      COMMON/COM106/ESCF(150),GNORM(150),SNORM(150)                     GNE00190
      COMMON/COM107/R(1275)                                             GNE00200
      COMMON/COM108/LDIM,MDIM,NDIM,IPRNT                                GNE00210
      COMMON/COM109/XMAX,YMAX,ZMAX                                      GNE00220
      COMMON/COM110/NSET(150),NORD(150),NXVAR(150)                      GNE00230
      COMMON/COM112/XOLD(150),XNEW(150)                                 GNE00240
      COMMON/COM113/STEPMX                                              GNE00250
cyy   EQUIVALENCE (CC,IA)                                               GNE00260
      DATA A00,TEN / 0.0D+00 , 10.0D+00 /                               GNE00270
      DATA FMILI / 1.0D-3 /                                             GNE00280
    1 FORMAT(//,2X,' THE GEOMETRY OPTIMIZATION PROGRAM'/                GNE00290
     1          2X,' BASED ON THE STEEPEST DESCENT METHOD'/             GNE00300
     2          2X,' CARTESIAN COORDINATE SYSTEM'/                      GNE00310
     3          2X,' (VERSION OF 11/26/97)'/)                           GNE00320
    2 FORMAT(2I5,F20.10,I5)                                             GNE00330
    3 FORMAT(12I5)                                                      GNE00340
    4 FORMAT(14I5)                                                      GNE00350
    5 FORMAT(A5,3I5)                                                    GNE00360
    6 FORMAT(2X,' NCOOD   = ',I5/                                       GNE00370
     1 2X,' NSORT   = ',I5/                                             GNE00380
     2 2X,' IMETR   = ',I5/                                             GNE00390
     3 2X,' NVAR    = ',I5/                                             GNE00400
     4 2X,' ISTEP   = ',I5/                                             GNE00410
     5 2X,' IPRNT   = ',I5/                                             GNE00420
     6 2X,' NATOM   = ',I5/                                             GNE00430
     7 2X,' N3N     = ',I5/                                             GNE00440
     A 2X,' NINT    = ',I5/                                             GNE00450
     B 2X,' IFORCE  = ',I5/                                             GNE00460
     C 2X,' NSTORE  = ',I5/                                             GNE00470
     D 2X,' NCASE   = ',I5/                                             GNE00480
     E 2X,' NPLUS   = ',I5/                                             GNE00490
     F 2X,' NDIF    = ',I5/                                             GNE00500
     G 2X,' NCHNG   = ',I5/                                             GNE00510
     H 2X,' NCONV   = ',I5/                                             GNE00520
     I 2X,' NUNIQ   = ',I5)                                             GNE00530
    7 FORMAT(//,2X,' THERE ARE NOT ENOUGH DATA TO USE THE QUASI NEWTON-RGNE00540
     1APHSON METHOD'/)                                                  GNE00550
    8 FORMAT(//,2X,' OLD COORDINATES IN TAPE30'/                        GNE00560
     1 39X,' OLD X IN A.U.'/                                            GNE00570
     2 4X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z/)               GNE00580
    9 FORMAT(3X,I5,1X,4F16.10)                                          GNE00590
   10 FORMAT(//,2X,' NEW COORDINATES IN TAPE30'/                        GNE00600
     1 39X,' NEW X IN A.U.'/                                            GNE00610
     2 4X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z/)               GNE00620
   11 FORMAT(//,2X,' GEOMETRY IN TAPE30 IS UPDATED'/)                   GNE00630
   12 FORMAT(//,2X,':::::::::::::::::::::::::::'/                       GNE00640
     1          2X,':::GEOMETRY IS OPTIMIZED:::'/                       GNE00650
     2          2X,':::::::::::::::::::::::::::'//                      GNE00660
     3          2X,' GNORM = ',F15.10,7X,' GCONV = ',F15.10/)           GNE00670
   13 FORMAT(//,2X,' TOO BIG CHANGE IN GEOMETRY---SCRATCH YOUR HEAD TO FGNE00680
     1IND A NEW GEOMETRY'//                                             GNE00690
     2          2X,' XMAX = ',F15.10,7X,' YMAX = ',F15.10,7X,' ZMAX = ',GNE00700
     3                        F15.10,7X,' GCHNG = ',F15.10/)            GNE00710
   14 FORMAT(//,2X,' SHORT OF INFORMATION---SCRATCH YOUR HEAD TO FIND A GNE00720
     1NEW GEOMETRY'/)                                                   GNE00730
C                                                                       GNE00740
      call psinit('APPEND')
      CALL TSTART(6)                                                    GNE00750
C* 10/04/91  CALL NOUNFL                                                GNE00760
C                                                                       GNE00770
      LDIM=1275                                                         GNE00780
      MDIM=150                                                          GNE00790
      NDIM=50                                                           GNE00800
      INPUT=5                                                           GNE00810
      ITAPE6=6                                                          GNE00820
      ITAP12=12                                                         GNE00830
      ITAP13=13                                                         GNE00840
      ITAP15=15                                                         GNE00850
      ITAP30=30                                                         GNE00860
      IFLAG=0                                                           GNE00870
      MAXCOR=360000                                                     GNE00880
C                                                                       GNE00890
      WRITE(6,1)                                                        GNE00900
C*  following three lines were added for the Unix Version (10/04/91)
      call ffile(itap12,'file12',0)
      call ffile(itap13,'file13',0)
      call ffile(itap15,'file15',0)
      REWIND ITAP12                                                     GNE00910
      REWIND ITAP13                                                     GNE00920
      READ(ITAP12,2) NATOM,NINT,ESCFZ,IFORCE                            GNE00930
C                                                                       GNE00940
C   NCOOD IS NUMBER OF INTERNAL COORDINATES TO BE OPTIMIZED             GNE00950
C   NSORT IS A PARAMETER FOR SORTING                                    GNE00960
C     NSORT=0---NO SORT                                                 GNE00970
C     NSORT=1---GRADIENT SORT                                           GNE00980
C     NSORT=2---ENERGY SORT                                             GNE00990
C     NSORT=3---SELECTED GRADIENT SORT                                  GNE01000
C   IMETR IS A PARAMETER FOR THE VARIABLE METRIC METHOD                 GNE01010
C     IMETR=0---SKIP THIS METHOD                                        GNE01020
C     IMETR=1---MURTAGH-SARGENT METHOD                                  GNE01030
C     IMETR=2---FLETCHER METHOD                                         GNE01040
C     IMETR=3---DAVIDON-FLETCHER-POWELL METHOD                          GNE01050
C   NVAR IS NUMBER OF CARTESIOAN COORDINATES TO BE OPTIMIZED            GNE01060
C   ISTEP IS A PARAMETER FOR MAXIMUM STEP SIZE ALLOWED                  GNE01070
C     ISTEP=0---0.1                                                     GNE01080
C     ISTEP=N---0.01*ISTEP                                              GNE01090
C   IDUM1,IDUM2,IDUM3, AND IDUM4 ARE DUMMY PARAMETERS TO ADJUST TO      GNE01100
C     INPUT FOR OTHER OPTIMIZATION PROGRAMS                             GNE01110
C   IPRNT IS A PARAMETER FOR OUTPUT                                     GNE01120
      N3N=NATOM*3                                                       GNE01130
      CALL LOCATE(INPUT,'# OPTIM ##',IERR)                              GNE01140
      READ(5,3) NCOOD,NSORT,IMETR,NVAR,ISTEP,IDUM1,IDUM2,IDUM3,IDUM4,   GNE01150
     1          IPRNT                                                   GNE01160
      READ(5,4) (NSET(I),I=1,NCOOD)                                     GNE01170
      DO 101 I=1,N3N                                                    GNE01180
  101 NXVAR(I)=I                                                        GNE01190
      IF(ISTEP.EQ.0) ISTEP=10                                           GNE01200
      STEPMX=1.0D-2*DFLOAT(ISTEP)                                       GNE01210
      IF(NVAR.EQ.0) NVAR=N3N                                            GNE01220
      IF(NVAR.EQ.N3N) GO TO 201                                         GNE01230
      READ(5,4) (NXVAR(I),I=1,NVAR)                                     GNE01240
  201 READ(ITAP13,3) NSTORE                                             GNE01250
      NCASE=NCOOD+1                                                     GNE01260
      NPLUS=NSTORE+1                                                    GNE01270
      NDIF=NPLUS-NCASE                                                  GNE01280
      READ(5,5) UPD,NCHNG,NCONV,NUNIQ                                   GNE01290
      IF(UPD.EQ.'UP     ') UPDATE=.TRUE.                                GNE01300
      IF(NCHNG.EQ.0) NCHNG=2                                            GNE01310
      IF(NCONV.EQ.0) NCONV=7                                            GNE01320
      IF(NUNIQ.EQ.0) NUNIQ=NATOM                                        GNE01330
      WRITE(6,6) NCOOD,NSORT,IMETR,NVAR,ISTEP,IPRNT,                    GNE01340
     1 NATOM,N3N,NINT,IFORCE,NSTORE,NCASE,NPLUS,NDIF,NCHNG,NCONV,NUNIQ  GNE01350
      GCHNG=TEN**(-NCHNG)                                               GNE01360
      GCONV=TEN**(-NCONV)                                               GNE01370
      IF(.NOT.UPDATE) GO TO 202                                         GNE01380
      DO 102 I=1,NATOM                                                  GNE01390
      NUNQ(I)=1                                                         GNE01400
      NAT(I,1)=I                                                        GNE01410
  102 CONTINUE                                                          GNE01420
      IF(NUNIQ.EQ.NATOM) GO TO 202                                      GNE01430
      DO 103 I=1,NUNIQ                                                  GNE01440
      READ(5,4) NN,(NAT(I,J),J=1,NN)                                    GNE01450
      NUNQ(I)=NN                                                        GNE01460
  103 CONTINUE                                                          GNE01470
C                                                                       GNE01480
C#######################################################################GNE01490
C   DYNAMIC ALLOCATION                                                  GNE01500
C#######################################################################GNE01510
C   1 : SX(N3N,NPLUS)                                                   GNE01520
C   2 : SG(N3N,NPLUS)                                                   GNE01530
C   3 : SCOOD(NINT,NPLUS)                                               GNE01540
C   4 : SENG(NINT,NPLUS)                                                GNE01550
C#######################################################################GNE01560
  202 CONTINUE                                                          GNE01570
      IC1=1                                                             GNE01580
      IC2=IC1+N3N*NPLUS                                                 GNE01590
      IC3=IC2+N3N*NPLUS                                                 GNE01600
      IC4=IC3+NINT*NPLUS                                                GNE01610
      ICSAVE=IC4+NINT*NPLUS                                             GNE01620
C#######################################################################GNE01630
C                                                                       GNE01640
C   READ IN STORED INFORMATION FROM TAPE13                              GNE01650
      IC5=ICSAVE                                                        GNE01660
      IC6=IC5+N3N*NPLUS                                                 GNE01670
      IC7=IC6+N3N*NPLUS                                                 GNE01680
      IC8=IC7+NINT*NPLUS                                                GNE01690
      ICMAX=IC8+NINT*NPLUS                                              GNE01700
C.................X       G       COOD    ENG.....                      GNE01710
      CALL READIN(CC(IC5),CC(IC6),CC(IC7),CC(IC8))                      GNE01720
C                                                                       GNE01730
C   SORT DATA W.R.T. GRADIENT OR ENERGY                                 GNE01740
      IC5=ICSAVE                                                        GNE01750
      IC6=IC5+N3N*NPLUS                                                 GNE01760
      IC7=IC6+N3N*NPLUS                                                 GNE01770
      IC8=IC7+NINT*NPLUS                                                GNE01780
      ICMAX=IC8+NINT*NPLUS                                              GNE01790
C...............SX      SG      SCOOD   SENG    X       G.......        GNE01800
      CALL SORT(CC(IC1),CC(IC2),CC(IC3),CC(IC4),CC(IC5),CC(IC6),        GNE01810
C...............COOD    ENG.....                                        GNE01820
     1          CC(IC7),CC(IC8))                                        GNE01830
C                                                                       GNE01840
      GNORML=GNORM(NPLUS)*FMILI                                         GNE01850
      IF(NSORT.EQ.3) GNORML=SNORM(NPLUS)*FMILI                          GNE01860
      IF(GNORML.LT.GCONV) GO TO 215                                     GNE01870
C                                                                       GNE01880
C*********************************                                      GNE01890
C***QUASI NEWTON-RAPHSON METHOD***                                      GNE01900
C*********************************                                      GNE01910
C:::THIS SECTION ONLY HANDLES THE INTERNAL COORDINATE SYSTEM:::         GNE01920
C                                                                       GNE01930
      IF(NDIF.LT.0) GO TO 203                                           GNE01940
      IC5=ICSAVE                                                        GNE01950
      IC6=IC5+NCOOD*NCOOD                                               GNE01960
      IC7=IC6+NCOOD*NCOOD                                               GNE01970
      IC8=IC7+NCOOD*NCOOD                                               GNE01980
      IC9=IC8+NCOOD                                                     GNE01990
      IA9=IC9+IC9-1                                                     GNE02000
      ICMAX=IC9+NCOOD                                                   GNE02010
C................SX      SG      SCOOD   SENG    FC      FT......       GNE02020
      CALL QUASI(CC(IC1),CC(IC2),CC(IC3),CC(IC4),CC(IC5),CC(IC6),       GNE02030
C................EE      PIVOT   INDEX...                               GNE02040
     1           CC(IC7),CC(IC8),IA(IA9))                               GNE02050
C                                                                       GNE02060
C   ESTIMATE FORCE CONSTANT MATRIX                                      GNE02070
      IC5=ICSAVE                                                        GNE02080
      IC6=IC5+NCOOD*NCOOD                                               GNE02090
      IC7=IC6+NCOOD*NCOOD                                               GNE02100
      IC8=IC7+NCOOD*NCOOD                                               GNE02110
      IC9=IC8+NCOOD*NCOOD                                               GNE02120
      IC10=IC9+NCOOD*NCOOD                                              GNE02130
      IC11=IC10+NCOOD                                                   GNE02140
      IA11=IC11+IC11-1                                                  GNE02150
      ICMAX=IC11+NCOOD                                                  GNE02160
C.................SCOOD   SENG    FC      DELC    FF      DELG....      GNE02170
      CALL FESTIM(CC(IC3),CC(IC4),CC(IC5),CC(IC6),CC(IC7),CC(IC8),      GNE02180
C.................EE      PIVOT    INDEX....                            GNE02190
     1            CC(IC9),CC(IC10),IA(IA11))                            GNE02200
      IF(DABS(XMAX).GT.GCHNG) GO TO 204                                 GNE02210
      IF(DABS(YMAX).GT.GCHNG) GO TO 204                                 GNE02220
      IF(DABS(ZMAX).GT.GCHNG) GO TO 204                                 GNE02230
      GO TO 212                                                         GNE02240
  203 CONTINUE                                                          GNE02250
      WRITE(6,7)                                                        GNE02260
C                                                                       GNE02270
C********************************                                       GNE02280
C***THE VARIABLE METRIC METHOD***                                       GNE02290
C********************************                                       GNE02300
C:::THIS SECTION ONLY HANDLES THE CARTESIAN COORDINATE SYSTEM:::        GNE02310
  204 CONTINUE                                                          GNE02320
      IF(.NOT.UPDATE.AND.IMETR.EQ.0) GO TO 210                          GNE02330
      IC5=ICSAVE                                                        GNE02340
      IC6=IC5+N3N*N3N                                                   GNE02350
      IC7=IC6+NVAR*NVAR                                                 GNE02360
      IC8=IC7+NVAR*NVAR                                                 GNE02370
      IC9=IC8+NVAR*NVAR                                                 GNE02380
      IC10=IC9+NVAR                                                     GNE02390
      IA10=IC10+IC10-1                                                  GNE02400
      ICMAX=IC10+NVAR                                                   GNE02410
C.................SX      SG      FX      HH      FT      EE......      GNE02420
      CALL METRIC(CC(IC1),CC(IC2),CC(IC5),CC(IC6),CC(IC7),CC(IC8),      GNE02430
C.................PIVOT   INDEX......................                   GNE02440
     1            CC(IC9),IA(IA10),NVAR,IMETR,UPDATE)                   GNE02450
C                                                                       GNE02460
C                                                                       GNE02470
C*******************************                                        GNE02480
C***UPDATE TAPE30 IF REQUIRED***                                        GNE02490
C*******************************                                        GNE02500
C                                                                       GNE02510
  210 CONTINUE                                                          GNE02520
      IF(.NOT.UPDATE.AND.NDIF.LT.0) GO TO 218                           GNE02530
C                                                                       GNE02540
      IF(DABS(XMAX).GT.GCHNG) GO TO 217                                 GNE02550
      IF(DABS(YMAX).GT.GCHNG) GO TO 217                                 GNE02560
      IF(DABS(ZMAX).GT.GCHNG) GO TO 217                                 GNE02570
C                                                                       GNE02580
C   READ IN GEOMETRY FROM TAPE30                                        GNE02590
  212 CONTINUE                                                          GNE02600
      CALL RFILE(ITAP30)                                                GNE02610
      CALL WREADW(ITAP30,I30,200,101,JUNK)                              GNE02620
      MPOINT=I30(2)                                                     GNE02630
      MCONST=I30(3)                                                     GNE02640
      NCALCS=I30(5)                                                     GNE02650
      IGEOP=100+MCONST+MPOINT+NCALCS                                    GNE02660
      CALL WREADW(ITAP30,LOCCAL,1,IGEOP,JUNK)                           GNE02670
      IGEOP=LOCCAL+80                                                   GNE02680
      CALL WREADW(ITAP30,COORD,N3N*2,IGEOP,JUNK)                        GNE02690
      WRITE(6,8)                                                        GNE02700
      DO 104 I=1,NATOM                                                  GNE02710
      WRITE(6,9) I,CHG(I),COORD(1,I),COORD(2,I),COORD(3,I)              GNE02720
  104 CONTINUE                                                          GNE02730
      DO 107 I=1,NUNIQ                                                  GNE02740
      NN=NUNQ(I)                                                        GNE02750
      SUMX=A00                                                          GNE02760
      SUMY=A00                                                          GNE02770
      SUMZ=A00                                                          GNE02780
      DO 105 J=1,NN                                                     GNE02790
      MM=NAT(I,J)                                                       GNE02800
      SUMX=SUMX+DABS(XX(MM))                                            GNE02810
      SUMY=SUMY+DABS(YY(MM))                                            GNE02820
      SUMZ=SUMZ+DABS(ZZ(MM))                                            GNE02830
  105 CONTINUE                                                          GNE02840
      ANX=DBLE(NN)                                                      GNE02850
      AVGX=SUMX/ANX                                                     GNE02860
      AVGY=SUMY/ANX                                                     GNE02870
      AVGZ=SUMZ/ANX                                                     GNE02880
      DO 106 J=1,NN                                                     GNE02890
      MM=NAT(I,J)                                                       GNE02900
      XX(MM)=DSIGN(AVGX,XX(MM))                                         GNE02910
      YY(MM)=DSIGN(AVGY,YY(MM))                                         GNE02920
      ZZ(MM)=DSIGN(AVGZ,ZZ(MM))                                         GNE02930
  106 CONTINUE                                                          GNE02940
  107 CONTINUE                                                          GNE02950
      WRITE(6,10)                                                       GNE02960
      DO 108 I=1,NATOM                                                  GNE02970
      COORD(1,I)=XX(I)                                                  GNE02980
      COORD(2,I)=YY(I)                                                  GNE02990
      COORD(3,I)=ZZ(I)                                                  GNE03000
      WRITE(6,9) I,CHG(I),COORD(1,I),COORD(2,I),COORD(3,I)              GNE03010
  108 CONTINUE                                                          GNE03020
      CALL WWRITW(ITAP30,COORD,N3N*2,IGEOP,JUNK)                        GNE03030
      WRITE(6,11)                                                       GNE03040
      CALL RCLOSE(ITAP30,3)                                             GNE03050
      GO TO 220                                                         GNE03060
C                                                                       GNE03070
C   GEOMETRY IS OPTIMIZED                                               GNE03080
  215 CONTINUE                                                          GNE03090
      WRITE(6,12) GNORML,GCONV                                          GNE03100
      GO TO 220                                                         GNE03110
C                                                                       GNE03120
C   GEOMETRY CHANGE IS TOO LARGE                                        GNE03130
  217 CONTINUE                                                          GNE03140
      WRITE(6,13) XMAX,YMAX,ZMAX,GCHNG                                  GNE03150
      GO TO 220                                                         GNE03160
C                                                                       GNE03170
C   SHORT OF INFORMATION                                                GNE03180
  218 CONTINUE                                                          GNE03190
      WRITE(6,14)                                                       GNE03200
C                                                                       GNE03210
  220 CONTINUE                                                          GNE03220
C                                                                       GNE03230
      CALL TSTOP(6)                                                     GNE03240
      STOP                                                              GNE03250
      END                                                               GNE03260
      SUBROUTINE READIN(X,G,COOD,ENG)                                   GNE03270
      IMPLICIT REAL*8 (A-H,O-Z)                                         GNE03280
      DIMENSION X(N3N,NPLUS),G(N3N,NPLUS)                               GNE03290
      DIMENSION COOD(NINT,NPLUS),ENG(NINT,NPLUS)                        GNE03300
      COMMON/COM101/NATOM,N3N,NINT                                      GNE03310
      COMMON/COM102/NSORT,NCOOD,NDIF,NPLUS,NCASE                        GNE03320
      COMMON/COM103/CHG(50),XX(50),YY(50),ZZ(50)                        GNE03330
      COMMON/COM104/EX(50),EY(50),EZ(50)                                GNE03340
      COMMON/COM105/ENGR(150),RCOOD(150)                                GNE03350
      COMMON/COM106/ESCF(150),GNORM(150),SNORM(150)                     GNE03360
      COMMON/COM108/LDIM,MDIM,NDIM,IPRNT                                GNE03370
      COMMON/COM110/NSET(150),NORD(150),NXVAR(150)                      GNE03380
      DATA ITAP12,ITAP13 / 12 , 13 /                                    GNE03390
      DATA FKILO,FMILI / 1.0D+3 , 1.0D-3 /                              GNE03400
    1 FORMAT(2I5,F20.10,I5)                                             GNE03410
    2 FORMAT(10I5)                                                      GNE03420
    3 FORMAT(I5,F20.10)                                                 GNE03430
    4 FORMAT(////,2X,10H CASE   = ,I5/                                  GNE03440
     1 2X,10H ESCF   = ,F20.10)                                         GNE03450
    5 FORMAT(//,2X,22H CARTESIAN COORDINATES,51X,10H GRADIENTS/         GNE03460
     1 4X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z,                GNE03470
     2                       17X,2H X,14X,2H Y,14X,2H Z)                GNE03480
    6 FORMAT(4F20.10)                                                   GNE03490
    7 FORMAT(3X,I5,1X,4F16.10,3X,3F16.10)                               GNE03500
    8 FORMAT(//,2X,' INTERNAL GRADIENTS (IN 1.0D+3 A.U.)'/              GNE03510
     1 5X,4H NO.,10X,7H COORD.,13X,6H GRAD./)                           GNE03520
    9 FORMAT(2F20.10)                                                   GNE03530
   10 FORMAT(2X,I7,5X,2F20.10)                                          GNE03540
   11 FORMAT(/,2X,22H NORM FOR GRADIENTS = ,F15.7)                      GNE03550
   12 FORMAT(//,2X,' ***SELECTED VARIABLES DO NOT MATCH***'/            GNE03560
     1 2X,' NCOOD = ',I5,' II = ',I5/)                                  GNE03570
   13 FORMAT(/,2X,' NORM FOR CASE ',I5,3X,' IS ',F15.7/)                GNE03580
   14 FORMAT(///,2X,' CARTESIAN COORDINATE MATRIX'/)                    GNE03590
   15 FORMAT(//,2X,' INTERNAL COORDINATE MATRIX'/)                      GNE03600
   16 FORMAT(//,2X,' CARTESIAN GRADIENT MATRIX'/)                       GNE03610
   17 FORMAT(//,2X,' INTERNAL GRADINET MATRIX'/)                        GNE03620
   18 FORMAT(//,2X,' GRADIENTS NORMS AND ENERGIES'/                     GNE03630
     1 3X,4H NO.,8X,6H GNORM,14X,6H SNORM,14X,7H ENERGY/)               GNE03640
   19 FORMAT(2X,I4,5X,3F20.10)                                          GNE03650
C                                                                       GNE03660
      REWIND ITAP12                                                     GNE03670
      READ(ITAP12,1) NATOM,NINT,ESCFZ,IFORCE                            GNE03680
      REWIND ITAP13                                                     GNE03690
      READ(ITAP13,2) NSTORE                                             GNE03700
C                                                                       GNE03710
      DO 105 ICASE=1,NPLUS                                              GNE03720
      IF(ICASE.EQ.NPLUS) GO TO 202                                      GNE03730
      READ(ITAP13,3) JCASE,ESCF(ICASE)                                  GNE03740
      GO TO 203                                                         GNE03750
  202 ESCF(ICASE)=ESCFZ                                                 GNE03760
  203 IF(IPRNT.LE.0) GO TO 301                                          GNE03770
      WRITE(6,4) ICASE,ESCF(ICASE)                                      GNE03780
  301 CONTINUE                                                          GNE03790
      IF(IPRNT.LE.0) GO TO 302                                          GNE03800
      WRITE(6,5)                                                        GNE03810
C***READ CARTESIAN COORDINATES OF EACH CASE                             GNE03820
  302 DO 102 I=1,NATOM                                                  GNE03830
      IF(ICASE.EQ.NPLUS) GO TO 204                                      GNE03840
      READ(ITAP13,6) CHG(I),XX(I),YY(I),ZZ(I)                           GNE03850
      GO TO 102                                                         GNE03860
  204 READ(ITAP12,6) CHG(I),XX(I),YY(I),ZZ(I)                           GNE03870
  102 CONTINUE                                                          GNE03880
C***READ CARTESIAN GRADIENTS OF EACH CASE                               GNE03890
      DO 103 I=1,NATOM                                                  GNE03900
      IF(ICASE.EQ.NPLUS) GO TO 205                                      GNE03910
      READ(ITAP13,6) CHG(I),EX(I),EY(I),EZ(I)                           GNE03920
      GO TO 206                                                         GNE03930
  205 READ(ITAP12,6) CHG(I),EX(I),EY(I),EZ(I)                           GNE03940
  206 IF(IPRNT.LE.0) GO TO 303                                          GNE03950
      WRITE(6,7) I,CHG(I),XX(I),YY(I),ZZ(I),EX(I),EY(I),EZ(I)           GNE03960
  303 IX=(I-1)*3+1                                                      GNE03970
      IY=IX+1                                                           GNE03980
      IZ=IX+2                                                           GNE03990
      X(IX,ICASE)=XX(I)                                                 GNE04000
      X(IY,ICASE)=YY(I)                                                 GNE04010
      X(IZ,ICASE)=ZZ(I)                                                 GNE04020
      G(IX,ICASE)=EX(I)                                                 GNE04030
      G(IY,ICASE)=EY(I)                                                 GNE04040
      G(IZ,ICASE)=EZ(I)                                                 GNE04050
  103 CONTINUE                                                          GNE04060
      IF(IPRNT.LE.0) GO TO 304                                          GNE04070
      CALL DIST(NATOM,XX,YY,ZZ,NDIM)                                    GNE04080
      WRITE(6,8)                                                        GNE04090
C***READ COORDINATES AND GRADIENTS OF EACH CASE                         GNE04100
  304 DO 104 I=1,NINT                                                   GNE04110
      IF(ICASE.EQ.NPLUS) GO TO 207                                      GNE04120
      READ(ITAP13,9) RCOOD(I),ENGR(I)                                   GNE04130
      GO TO 208                                                         GNE04140
  207 READ(ITAP12,9) RCOOD(I),ENGR(I)                                   GNE04150
  208 ENGR(I)=ENGR(I)*FKILO                                             GNE04160
      IF(IPRNT.LE.0) GO TO 305                                          GNE04170
      WRITE(6,10) I,RCOOD(I),ENGR(I)                                    GNE04180
  305 COOD(I,ICASE)=RCOOD(I)                                            GNE04190
      ENG(I,ICASE)=ENGR(I)                                              GNE04200
  104 CONTINUE                                                          GNE04210
C   CALCULATE NORM OF GRADIENTS                                         GNE04220
      CALL ANORM(CNORM,ENGR,MDIM,NINT)                                  GNE04230
      IF(IPRNT.LE.0) GO TO 306                                          GNE04240
      WRITE(6,11) CNORM                                                 GNE04250
  306 GNORM(ICASE)=CNORM                                                GNE04260
  105 CONTINUE                                                          GNE04270
C                                                                       GNE04280
      DO 108 ICASE=1,NPLUS                                              GNE04290
      II=0                                                              GNE04300
      DO 107 I=1,NINT                                                   GNE04310
      DO 106 J=1,NCOOD                                                  GNE04320
      IF(NSET(J).NE.I) GO TO 106                                        GNE04330
      II=II+1                                                           GNE04340
      ENGR(II)=ENG(I,ICASE)                                             GNE04350
  106 CONTINUE                                                          GNE04360
  107 CONTINUE                                                          GNE04370
      IF(II.EQ.NCOOD) GO TO 307                                         GNE04380
      WRITE(6,12) NCOOD,II                                              GNE04390
      STOP                                                              GNE04400
  307 CONTINUE                                                          GNE04410
      CALL ANORM(CNORM,ENGR,MDIM,NCOOD)                                 GNE04420
      IF(IPRNT.LE.0) GO TO 308                                          GNE04430
      WRITE(6,13) ICASE,CNORM                                           GNE04440
  308 SNORM(ICASE)=CNORM                                                GNE04450
  108 CONTINUE                                                          GNE04460
C                                                                       GNE04470
C   SUMMARY OF INPUT                                                    GNE04480
      IF(IPRNT.LE.0) GO TO 309                                          GNE04490
      WRITE(6,14)                                                       GNE04500
      CALL MATOUT(X,N3N,NPLUS,N3N,NPLUS,6)                              GNE04510
      WRITE(6,15)                                                       GNE04520
      CALL MATOUT(COOD,NINT,NPLUS,NINT,NPLUS,6)                         GNE04530
      WRITE(6,16)                                                       GNE04540
      CALL MATOUT(G,N3N,NPLUS,N3N,NPLUS,6)                              GNE04550
      WRITE(6,17)                                                       GNE04560
      CALL MATOUT(ENG,NINT,NPLUS,NINT,NPLUS,6)                          GNE04570
  309 WRITE(6,18)                                                       GNE04580
      DO 110 ICASE=1,NPLUS                                              GNE04590
      WRITE(6,19) ICASE,GNORM(ICASE),SNORM(ICASE),ESCF(ICASE)           GNE04600
  110 CONTINUE                                                          GNE04610
C                                                                       GNE04620
C   STORE GEOMETRIES AND GRADIENTS IN TAPE13                            GNE04630
      REWIND ITAP13                                                     GNE04640
      WRITE(ITAP13,2) NPLUS                                             GNE04650
      DO 115 ICASE=1,NPLUS                                              GNE04660
      WRITE(ITAP13,3) ICASE,ESCF(ICASE)                                 GNE04670
      DO 111 I=1,NATOM                                                  GNE04680
      IX=(I-1)*3+1                                                      GNE04690
      IY=IX+1                                                           GNE04700
      IZ=IX+2                                                           GNE04710
      XX(I)=X(IX,ICASE)                                                 GNE04720
      YY(I)=X(IY,ICASE)                                                 GNE04730
      ZZ(I)=X(IZ,ICASE)                                                 GNE04740
      WRITE(ITAP13,6) CHG(I),XX(I),YY(I),ZZ(I)                          GNE04750
  111 CONTINUE                                                          GNE04760
      DO 112 I=1,NATOM                                                  GNE04770
      IX=(I-1)*3+1                                                      GNE04780
      IY=IX+1                                                           GNE04790
      IZ=IX+2                                                           GNE04800
      EX(I)=G(IX,ICASE)                                                 GNE04810
      EY(I)=G(IY,ICASE)                                                 GNE04820
      EZ(I)=G(IZ,ICASE)                                                 GNE04830
      WRITE(ITAP13,6) CHG(I),EX(I),EY(I),EZ(I)                          GNE04840
  112 CONTINUE                                                          GNE04850
      DO 113 I=1,NINT                                                   GNE04860
      RCOOD(I)=COOD(I,ICASE)                                            GNE04870
      ENGR(I)=ENG(I,ICASE)*FMILI                                        GNE04880
      WRITE(ITAP13,9) RCOOD(I),ENGR(I)                                  GNE04890
  113 CONTINUE                                                          GNE04900
  115 CONTINUE                                                          GNE04910
C                                                                       GNE04920
      REWIND ITAP12                                                     GNE04930
      REWIND ITAP13                                                     GNE04940
C                                                                       GNE04950
      RETURN                                                            GNE04960
      END                                                               GNE04970
      SUBROUTINE DIST(NATOM,XX,YY,ZZ,NNX)                               GNE04980
      IMPLICIT REAL*8 (A-H,O-Z)                                         GNE04990
      DIMENSION XX(NNX),YY(NNX),ZZ(NNX)                                 GNE05000
      COMMON/COM107/R(1275)                                             GNE05010
      COMMON/COM108/LDIM,MDIM,NDIM,IPRNT                                GNE05020
      DATA A00 / 0.0D+00 /                                              GNE05030
    1 FORMAT(//,2X,28H INTERATOMIC DISTANCE MATRIX/)                    GNE05040
C                                                                       GNE05050
      IJ=0                                                              GNE05060
      DO 101 I=1,NATOM                                                  GNE05070
      DO 101 J=1,I                                                      GNE05080
      IJ=IJ+1                                                           GNE05090
      R(IJ)=A00                                                         GNE05100
      IF(I.EQ.J) GO TO 101                                              GNE05110
      XD=XX(I)-XX(J)                                                    GNE05120
      YD=YY(I)-YY(J)                                                    GNE05130
      ZD=ZZ(I)-ZZ(J)                                                    GNE05140
      RR=DSQRT(XD*XD+YD*YD+ZD*ZD)                                       GNE05150
      R(IJ)=RR                                                          GNE05160
  101 CONTINUE                                                          GNE05170
      WRITE(6,1)                                                        GNE05180
      CALL PRINT(R,LDIM,NATOM,6)                                        GNE05190
      RETURN                                                            GNE05200
      END                                                               GNE05210
      SUBROUTINE ANGLE(NATOM,XX,YY,ZZ,NNX)                              GNE05220
      IMPLICIT REAL*8 (A-H,O-Z)                                         GNE05230
      DIMENSION XX(NNX),YY(NNX),ZZ(NNX)                                 GNE05240
      COMMON/COM107/R(1275)                                             GNE05250
      COMMON/COM108/LDIM,MDIM,NDIM,IPRNT                                GNE05260
      DATA PI / 3.1415926536D+00 /                                      GNE05270
    1 FORMAT(//,2X,' BOND ANGLES'/                                      GNE05280
     1 2X,' I   J   K',3X,' DEGREES   ',1X,' J   I   K',3X,' DEGREES  ',GNE05290
     2 1X,' I   K   J',3X,' DEGREES   '/)                               GNE05300
    2 FORMAT(3I4,F11.3, 2X,3I4,F11.3,1X,3I4,F11.3)                      GNE05310
    3 FORMAT(//,2X,' BOND ANGLES WILL NOT BE PRINTED---IPRNT .LT. 0 '/) GNE05320
C                                                                       GNE05330
      AR=PI/180.0D+00                                                   GNE05340
      IF(IPRNT.GE.0) WRITE(6,1)                                         GNE05350
      IF(IPRNT.LT.0) WRITE(6,3)                                         GNE05360
      DO 101 I=1,NATOM-2                                                GNE05370
      DO 101 J=I+1,NATOM-1                                              GNE05380
      IJ=I+J*(J-1)/2                                                    GNE05390
      XIJ=XX(I)-XX(J)                                                   GNE05400
      YIJ=YY(I)-YY(J)                                                   GNE05410
      ZIJ=ZZ(I)-ZZ(J)                                                   GNE05420
      EXIJ=-XIJ/R(IJ)                                                   GNE05430
      EYIJ=-YIJ/R(IJ)                                                   GNE05440
      EZIJ=-ZIJ/R(IJ)                                                   GNE05450
      EXJI=-EXIJ                                                        GNE05460
      EYJI=-EYIJ                                                        GNE05470
      EZJI=-EZIJ                                                        GNE05480
      DO 101 K=J+1,NATOM                                                GNE05490
      IK=I+K*(K-1)/2                                                    GNE05500
      XIK=XX(I)-XX(K)                                                   GNE05510
      YIK=YY(I)-YY(K)                                                   GNE05520
      ZIK=ZZ(I)-ZZ(K)                                                   GNE05530
      EXIK=-XIK/R(IK)                                                   GNE05540
      EYIK=-YIK/R(IK)                                                   GNE05550
      EZIK=-ZIK/R(IK)                                                   GNE05560
      EXKI=-EXIK                                                        GNE05570
      EYKI=-EYIK                                                        GNE05580
      EZKI=-EZIK                                                        GNE05590
      JK=J+K*(K-1)/2                                                    GNE05600
      XJK=XX(J)-XX(K)                                                   GNE05610
      YJK=YY(J)-YY(K)                                                   GNE05620
      ZJK=ZZ(J)-ZZ(K)                                                   GNE05630
      EXJK=-XJK/R(JK)                                                   GNE05640
      EYJK=-YJK/R(JK)                                                   GNE05650
      EZJK=-ZJK/R(JK)                                                   GNE05660
      EXKJ=-EXJK                                                        GNE05670
      EYKJ=-EYJK                                                        GNE05680
      EZKJ=-EZJK                                                        GNE05690
      AIJK=EXJI*EXJK+EYJI*EYJK+EZJI*EZJK                                GNE05700
      AJIK=EXIJ*EXIK+EYIJ*EYIK+EZIJ*EZIK                                GNE05710
      AIKJ=EXKI*EXKJ+EYKI*EYKJ+EZKI*EZKJ                                GNE05720
      ANGIJK=DACOS(AIJK)/AR                                             GNE05730
      ANGJIK=DACOS(AJIK)/AR                                             GNE05740
      ANGIKJ=DACOS(AIKJ)/AR                                             GNE05750
      IF(IPRNT.GE.0) WRITE(6,2) I,J,K,ANGIJK,J,I,K,ANGJIK,I,K,J,ANGIKJ  GNE05760
  101 CONTINUE                                                          GNE05770
C                                                                       GNE05780
      RETURN                                                            GNE05790
      END                                                               GNE05800
      SUBROUTINE ANORM(CNORM,C,NNC,NC)                                  GNE05810
      IMPLICIT REAL*8 (A-H,O-Z)                                         GNE05820
      DIMENSION C(NNC)                                                  GNE05830
      DATA A00 / 0.0D+00 /                                              GNE05840
C                                                                       GNE05850
      CNORM=A00                                                         GNE05860
      DO 101 I=1,NC                                                     GNE05870
      CNORM=CNORM+C(I)*C(I)                                             GNE05880
  101 CONTINUE                                                          GNE05890
      CNORM=DSQRT(CNORM)                                                GNE05900
      RETURN                                                            GNE05910
      END                                                               GNE05920
      SUBROUTINE SORT(SX,SG,SCOOD,SENG,X,G,COOD,ENG)                    GNE05930
      IMPLICIT REAL*8 (A-H,O-Z)                                         GNE05940
      DIMENSION SX(N3N,NPLUS),SG(N3N,NPLUS),X(N3N,NPLUS),G(N3N,NPLUS)   GNE05950
      DIMENSION SCOOD(NINT,NPLUS),SENG(NINT,NPLUS)                      GNE05960
      DIMENSION COOD(NINT,NPLUS),ENG(NINT,NPLUS)                        GNE05970
      COMMON/COM101/NATOM,N3N,NINT                                      GNE05980
      COMMON/COM102/NSORT,NCOOD,NDIF,NPLUS,NCASE                        GNE05990
      COMMON/COM106/ESCF(150),GNORM(150),SNORM(150)                     GNE06000
      COMMON/COM108/LDIM,MDIM,NDIM,IPRNT                                GNE06010
      COMMON/COM110/NSET(150),NORD(150),NXVAR(150)                      GNE06020
    1 FORMAT(//,2X,' SORTED CARTESIAN COORDINATE MATRIX'/)              GNE06030
    2 FORMAT(//,2X,' SORTED CARTESIAN GRADIENT MATRIX'/)                GNE06040
    3 FORMAT(//,2X,' SORTED INTERNAL COORDINATE MATRIX'/)               GNE06050
    4 FORMAT(//,2X,' SORTED INTERNAL GRADIENT MATRIX'/)                 GNE06060
C                                                                       GNE06070
      DO 101 I=1,NPLUS                                                  GNE06080
      NORD(I)=I                                                         GNE06090
  101 CONTINUE                                                          GNE06100
      IF(NPLUS.EQ.1) GO TO 205                                          GNE06110
      IF(NSORT.EQ.0) GO TO 205                                          GNE06120
      GO TO (201,202,203),NSORT                                         GNE06130
C                                                                       GNE06140
C   GRADIENT SORT                                                       GNE06150
  201 CALL ORDER(GNORM,NORD,MDIM,NPLUS)                                 GNE06160
      GO TO 205                                                         GNE06170
C                                                                       GNE06180
C   ENERGY SORT                                                         GNE06190
  202 CALL ORDER(ESCF,NORD,MDIM,NPLUS)                                  GNE06200
      GO TO 205                                                         GNE06210
C                                                                       GNE06220
C   SELECTED GRADIENT SORT                                              GNE06230
  203 CALL ORDER(SNORM,NORD,MDIM,NPLUS)                                 GNE06240
C                                                                       GNE06250
C   NO SORT AND/OR AFTER SORT                                           GNE06260
  205 II=0                                                              GNE06270
      DO 104 I=1,NPLUS                                                  GNE06280
      II=II+1                                                           GNE06290
      MM=NORD(I)                                                        GNE06300
      DO 102 J=1,NINT                                                   GNE06310
C***  DO 102 J=1,NCOOD                                                  GNE06320
C***  NN=NSET(J)                                                        GNE06330
C***  SCOOD(J,II)=COOD(NN,MM)                                           GNE06340
C***  SENG(J,II)=ENG(NN,MM)                                             GNE06350
      SCOOD(J,II)=COOD(J,MM)                                            GNE06360
      SENG(J,II)=ENG(J,MM)                                              GNE06370
  102 CONTINUE                                                          GNE06380
      DO 103 J=1,N3N                                                    GNE06390
      SX(J,II)=X(J,MM)                                                  GNE06400
      SG(J,II)=G(J,MM)                                                  GNE06410
  103 CONTINUE                                                          GNE06420
  104 CONTINUE                                                          GNE06430
      IF(IPRNT.LE.0) GO TO 206                                          GNE06440
      WRITE(6,1)                                                        GNE06450
      CALL MATOUT(SX,N3N,NPLUS,N3N,NPLUS,6)                             GNE06460
      WRITE(6,2)                                                        GNE06470
      CALL MATOUT(SG,N3N,NPLUS,N3N,NPLUS,6)                             GNE06480
  206 CONTINUE                                                          GNE06490
C***  WRITE(6,3)                                                        GNE06500
C***  CALL MATOUT(SCOOD,NINT,NPLUS,NCOOD,NPLUS,6)                       GNE06510
C***  WRITE(6,4)                                                        GNE06520
C***  CALL MATOUT(SENG,NINT,NPLUS,NCOOD,NPLUS,6)                        GNE06530
      WRITE(6,3)                                                        GNE06540
      CALL MATOUT(SCOOD,NINT,NPLUS,NINT,NPLUS,6)                        GNE06550
      WRITE(6,4)                                                        GNE06560
      CALL MATOUT(SENG,NINT,NPLUS,NINT,NPLUS,6)                         GNE06570
C                                                                       GNE06580
      RETURN                                                            GNE06590
      END                                                               GNE06600
      SUBROUTINE ORDER(A,NORD,NND,N)                                    GNE06610
      IMPLICIT REAL*8 (A-H,O-Z)                                         GNE06620
      DIMENSION A(NND),NORD(NND)                                        GNE06630
C                                                                       GNE06640
      N1=N-1                                                            GNE06650
      DO 102 I=1,N1                                                     GNE06660
      IR=I+1                                                            GNE06670
      DO 101 J=IR,N                                                     GNE06680
      IF(A(I).GT.A(J)) GO TO 101                                        GNE06690
      AA=A(I)                                                           GNE06700
      A(I)=A(J)                                                         GNE06710
      A(J)=AA                                                           GNE06720
      NN=NORD(I)                                                        GNE06730
      NORD(I)=NORD(J)                                                   GNE06740
      NORD(J)=NN                                                        GNE06750
  101 CONTINUE                                                          GNE06760
  102 CONTINUE                                                          GNE06770
      RETURN                                                            GNE06780
      END                                                               GNE06790
      SUBROUTINE QUASI(SX,SG,SCOOD,SENG,FC,FT,EE,PIVOT,INDEX)           GNE06800
      IMPLICIT REAL*8 (A-H,O-Z)                                         GNE06810
      DIMENSION SX(N3N,NPLUS),SG(N3N,NPLUS)                             GNE06820
      DIMENSION SCOOD(NINT,NPLUS),SENG(NINT,NPLUS)                      GNE06830
      DIMENSION FC(NCOOD,NCOOD),FT(NCOOD,NCOOD),EE(NCOOD,NCOOD)         GNE06840
      DIMENSION PIVOT(NCOOD),INDEX(NCOOD)                               GNE06850
      COMMON/COM101/NATOM,N3N,NINT                                      GNE06860
      COMMON/COM102/NSORT,NCOOD,NDIF,NPLUS,NCASE                        GNE06870
      COMMON/COM103/CHG(50),XX(50),YY(50),ZZ(50)                        GNE06880
      COMMON/COM104/XT(150)                                             GNE06890
      COMMON/COM105/ENGR(150),RCOOD(150)                                GNE06900
      COMMON/COM107/COEF(1275)                                          GNE06910
      COMMON/COM108/LDIM,MDIM,NDIM,IPRNT                                GNE06920
      COMMON/COM109/XMAX,YMAX,ZMAX                                      GNE06930
      COMMON/COM110/NSET(150),NORD(150),NXVAR(150)                      GNE06940
      COMMON/COM112/XOLD(150),XNEW(150)                                 GNE06950
      DATA A00,ONE / 0.0D+00 , 1.0D+00 /                                GNE06960
    1 FORMAT(//,2X,'*************************************'/             GNE06970
     1          2X,'***THE QUASI NEWTON-RAPHSON METHOD***'/             GNE06980
     2          2X,'*************************************'//)           GNE06990
    2 FORMAT(//,2X,' FC MATRIX'/)                                       GNE07000
    3 FORMAT(//,2X,' INVERSE GRADIENT MATRIX'/)                         GNE07010
    4 FORMAT(//,2X,' DETERM = ',D20.10/)                                GNE07020
    5 FORMAT(//,2X,' ENG * ENG(-1) = E '/)                              GNE07030
C                                                                       GNE07040
      WRITE(6,1)                                                        GNE07050
C                                                                       GNE07060
C   CALCULATE INVERSE FORCE CONSTANT MATRIX                             GNE07070
CCC   DO 101 I=1,NCOOD                                                  GNE07080
CCC   ENGR(I)=SENG(I,NPLUS)                                             GNE07090
CCC   DO 101 J=1,NCOOD                                                  GNE07100
CCC   JJ=NDIF+J                                                         GNE07110
CCC   FC(I,J)=SENG(I,JJ)                                                GNE07120
CC101 CONTINUE                                                          GNE07130
      DO 101 I=1,NCOOD                                                  GNE07140
      II=NSET(I)                                                        GNE07150
      ENGR(I)=SENG(II,NPLUS)                                            GNE07160
      DO 101 J=1,NCOOD                                                  GNE07170
      JJ=NDIF+J                                                         GNE07180
      FC(I,J)=SENG(II,JJ)                                               GNE07190
  101 CONTINUE                                                          GNE07200
      IF(IPRNT.LE.0) GO TO 201                                          GNE07210
      WRITE(6,2)                                                        GNE07220
      CALL MATOUT(FC,NCOOD,NCOOD,NCOOD,NCOOD,6)                         GNE07230
  201 CONTINUE                                                          GNE07240
      DO 102 I=1,NCOOD                                                  GNE07250
      DO 102 J=1,NCOOD                                                  GNE07260
  102 FT(I,J)=FC(I,J)                                                   GNE07270
      CALL MATINV(FC,EE,PIVOT,INDEX,NCOOD,NCOOD,DETERM,NCOOD)           GNE07280
      IF(IPRNT.LE.0) GO TO 202                                          GNE07290
      WRITE(6,3)                                                        GNE07300
      CALL MATOUT(FC,NCOOD,NCOOD,NCOOD,NCOOD,6)                         GNE07310
  202 WRITE(6,4) DETERM                                                 GNE07320
C                                                                       GNE07330
C   CHECK UNITARITY OF THE INVERTED MATRIX                              GNE07340
      CALL MTXMPY(FT,NCOOD,FC,NCOOD,EE,NCOOD,EE,NCOOD,NCOOD,1)          GNE07350
      IF(IPRNT.LE.0) GO TO 203                                          GNE07360
      WRITE(6,5)                                                        GNE07370
      CALL MATOUT(EE,NCOOD,NCOOD,NCOOD,NCOOD,6)                         GNE07380
C                                                                       GNE07390
C   ESTIMATE NEW GEOMETRY                                               GNE07400
  203 DO 103 I=1,NCOOD                                                  GNE07410
      COEF(I)=A00                                                       GNE07420
      DO 103 J=1,NCOOD                                                  GNE07430
      COEF(I)=COEF(I)+FC(I,J)*ENGR(J)                                   GNE07440
  103 CONTINUE                                                          GNE07450
      DO 104 I=1,N3N                                                    GNE07460
      XT(I)=A00                                                         GNE07470
      DO 104 J=1,NCOOD                                                  GNE07480
      JJ=NDIF+J                                                         GNE07490
      XT(I)=XT(I)+SX(I,JJ)*COEF(J)                                      GNE07500
  104 CONTINUE                                                          GNE07510
      SUM=A00                                                           GNE07520
      DO 105 I=1,NCOOD                                                  GNE07530
      SUM=SUM+COEF(I)                                                   GNE07540
  105 CONTINUE                                                          GNE07550
      FACT=ONE/(SUM-ONE)                                                GNE07560
C                                                                       GNE07570
      DO 106 I=1,N3N                                                    GNE07580
      XOLD(I)=SX(I,NPLUS)                                               GNE07590
      XNEW(I)=(XT(I)-XOLD(I))*FACT                                      GNE07600
  106 CONTINUE                                                          GNE07610
      IFLAG=1                                                           GNE07620
      CALL OLDNEW(XOLD,XNEW,IFLAG)                                      GNE07630
C                                                                       GNE07640
      RETURN                                                            GNE07650
      END                                                               GNE07660
      SUBROUTINE FESTIM(SCOOD,SENG,FC,DELC,FF,DELG,EE,PIVOT,INDEX)      GNE07670
      IMPLICIT REAL*8 (A-H,O-Z)                                         GNE07680
      DIMENSION SCOOD(NINT,NPLUS),SENG(NINT,NPLUS)                      GNE07690
      DIMENSION FC(NCOOD,NCOOD),DELC(NCOOD,NCOOD),FF(NCOOD,NCOOD)       GNE07700
      DIMENSION DELG(NCOOD,NCOOD),EE(NCOOD,NCOOD)                       GNE07710
      DIMENSION PIVOT(NCOOD),INDEX(NCOOD)                               GNE07720
      COMMON/COM101/NATOM,N3N,NINT                                      GNE07730
      COMMON/COM102/NSORT,NCOOD,NDIF,NPLUS,NCASE                        GNE07740
      COMMON/COM108/LDIM,MDIM,NDIM,IPRNT                                GNE07750
      COMMON/COM110/NSET(150),NORD(150),NXVAR(150)                      GNE07760
      DATA HALF / 0.5D+00 /                                             GNE07770
      DATA FKILO / 1.0D+3 /                                             GNE07780
    1 FORMAT(//,2X,' DELC MATRIX'/)                                     GNE07790
    2 FORMAT(//,2X,' DELG MATRIX'/)                                     GNE07800
    3 FORMAT(//,2X,' INVERSE DELC MATRIX'/)                             GNE07810
    4 FORMAT(//,2X,' DETERM = ',D20.10/)                                GNE07820
    5 FORMAT(//,2X,' DELC * DELC(-1) = E '/)                            GNE07830
    6 FORMAT(//,2X,' ESTIMATED FORCE CONSTANT MATRIX'/)                 GNE07840
    7 FORMAT(//,2X,' AVERAGED FORCE CONSTANT MATRIX'/)                  GNE07850
C                                                                       GNE07860
C   ESTIMATE FORCE CONSTANT MATRIX                                      GNE07870
CCC   DO 101 I=1,NCOOD                                                  GNE07880
CCC   DO 101 J=1,NCOOD                                                  GNE07890
CCC   JJ=NDIF+J                                                         GNE07900
CCC   DELC(I,J)=SCOOD(I,JJ)-SCOOD(I,NPLUS)                              GNE07910
CCC   DELC(I,J)=DELC(I,J)*FKILO                                         GNE07920
CC101 CONTINUE                                                          GNE07930
      DO 101 I=1,NCOOD                                                  GNE07940
      II=NSET(I)                                                        GNE07950
      DO 101 J=1,NCOOD                                                  GNE07960
      JJ=NDIF+J                                                         GNE07970
      DELC(I,J)=SCOOD(II,JJ)-SCOOD(II,NPLUS)                            GNE07980
      DELC(I,J)=DELC(I,J)*FKILO                                         GNE07990
  101 CONTINUE                                                          GNE08000
      IF(IPRNT.LE.0) GO TO 201                                          GNE08010
      WRITE(6,1)                                                        GNE08020
      CALL MATOUT(DELC,NCOOD,NCOOD,NCOOD,NCOOD,6)                       GNE08030
CC201 DO 102 I=1,NCOOD                                                  GNE08040
CCC   DO 102 J=1,NCOOD                                                  GNE08050
CCC   JJ=NDIF+J                                                         GNE08060
CCC   DELG(I,J)=SENG(I,JJ)-SENG(I,NPLUS)                                GNE08070
CC102 CONTINUE                                                          GNE08080
  201 DO 102 I=1,NCOOD                                                  GNE08090
      II=NSET(I)                                                        GNE08100
      DO 102 J=1,NCOOD                                                  GNE08110
      JJ=NDIF+J                                                         GNE08120
      DELG(I,J)=SENG(II,JJ)-SENG(II,NPLUS)                              GNE08130
  102 CONTINUE                                                          GNE08140
      IF(IPRNT.LE.0) GO TO 202                                          GNE08150
      WRITE(6,2)                                                        GNE08160
      CALL MATOUT(DELG,NCOOD,NCOOD,NCOOD,NCOOD,6)                       GNE08170
  202 DO 103 I=1,NCOOD                                                  GNE08180
      DO 103 J=1,NCOOD                                                  GNE08190
      FC(I,J)=DELC(I,J)                                                 GNE08200
  103 CONTINUE                                                          GNE08210
      CALL MATINV(FC,EE,PIVOT,INDEX,NCOOD,NCOOD,DETERM,NCOOD)           GNE08220
      IF(IPRNT.LE.0) GO TO 203                                          GNE08230
      WRITE(6,3)                                                        GNE08240
      CALL MATOUT(FC,NCOOD,NCOOD,NCOOD,NCOOD,6)                         GNE08250
  203 WRITE(6,4) DETERM                                                 GNE08260
      CALL MTXMPY(DELC,NCOOD,FC,NCOOD,EE,NCOOD,EE,NCOOD,NCOOD,1)        GNE08270
      IF(IPRNT.LE.0) GO TO 204                                          GNE08280
      WRITE(6,5)                                                        GNE08290
      CALL MATOUT(EE,NCOOD,NCOOD,NCOOD,NCOOD,6)                         GNE08300
  204 CALL MTXMPY(FC,NCOOD,DELG,NCOOD,FF,NCOOD,EE,NCOOD,NCOOD,4)        GNE08310
      IF(IPRNT.LE.0) GO TO 205                                          GNE08320
      WRITE(6,6)                                                        GNE08330
      CALL MATOUT(FF,NCOOD,NCOOD,NCOOD,NCOOD,6)                         GNE08340
  205 DO 104 I=1,NCOOD                                                  GNE08350
      DO 104 J=I,NCOOD                                                  GNE08360
      EE(I,J)=(FF(I,J)+FF(J,I))*HALF                                    GNE08370
      IF(I.EQ.J) GO TO 104                                              GNE08380
      EE(J,I)=EE(I,J)                                                   GNE08390
  104 CONTINUE                                                          GNE08400
      WRITE(6,7)                                                        GNE08410
      CALL MATOUT(EE,NCOOD,NCOOD,NCOOD,NCOOD,6)                         GNE08420
C                                                                       GNE08430
      RETURN                                                            GNE08440
      END                                                               GNE08450
      SUBROUTINE OLDNEW(XOLD,XNEW,IFLAG)                                GNE08460
      IMPLICIT REAL*8 (A-H,O-Z)                                         GNE08470
      DIMENSION XOLD(150),XNEW(150)                                     GNE08480
      COMMON/COM101/NATOM,N3N,NINT                                      GNE08490
      COMMON/COM103/CHG(50),XX(50),YY(50),ZZ(50)                        GNE08500
      COMMON/COM104/XT(50),YT(50),ZT(50)                                GNE08510
      COMMON/COM108/LDIM,MDIM,NDIM,IPRNT                                GNE08520
      COMMON/COM109/XMAX,YMAX,ZMAX                                      GNE08530
      DATA A00 / 0.0D+00 /                                              GNE08540
    1 FORMAT(//,2X,' OLD COORDINATES'/                                  GNE08550
     1 38X,' OLD X IN A.U.'/                                            GNE08560
     2 3X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z/)               GNE08570
    2 FORMAT(3X,I5,4F16.10)                                             GNE08580
    3 FORMAT(//,2X,' NEW COORDINATES'/                                  GNE08590
     1 38X,' NEW X IN A.U.'/                                            GNE08600
     2 3X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z/)               GNE08610
C                                                                       GNE08620
      WRITE(6,1)                                                        GNE08630
      XMAX=A00                                                          GNE08640
      YMAX=A00                                                          GNE08650
      ZMAX=A00                                                          GNE08660
      DO 101 I=1,NATOM                                                  GNE08670
      IX=(I-1)*3+1                                                      GNE08680
      IY=IX+1                                                           GNE08690
      IZ=IX+2                                                           GNE08700
      IF(IFLAG.EQ.0) GO TO 201                                          GNE08710
      XDIF=XNEW(IX)-XOLD(IX)                                            GNE08720
      YDIF=XNEW(IY)-XOLD(IY)                                            GNE08730
      ZDIF=XNEW(IZ)-XOLD(IZ)                                            GNE08740
      IF(DABS(XDIF).GT.DABS(XMAX)) XMAX=XDIF                            GNE08750
      IF(DABS(YDIF).GT.DABS(YMAX)) YMAX=YDIF                            GNE08760
      IF(DABS(ZDIF).GT.DABS(ZMAX)) ZMAX=ZDIF                            GNE08770
  201 CONTINUE                                                          GNE08780
      WRITE(6,2) I,CHG(I),XOLD(IX),XOLD(IY),XOLD(IZ)                    GNE08790
  101 CONTINUE                                                          GNE08800
      WRITE(6,3)                                                        GNE08810
      DO 102 I=1,NATOM                                                  GNE08820
      IX=(I-1)*3+1                                                      GNE08830
      IY=IX+1                                                           GNE08840
      IZ=IX+2                                                           GNE08850
      WRITE(6,2) I,CHG(I),XNEW(IX),XNEW(IY),XNEW(IZ)                    GNE08860
  102 CONTINUE                                                          GNE08870
C                                                                       GNE08880
      DO 103 I=1,NATOM                                                  GNE08890
      IX=(I-1)*3+1                                                      GNE08900
      IY=IX+1                                                           GNE08910
      IZ=IX+2                                                           GNE08920
      XT(I)=XNEW(IX)                                                    GNE08930
      YT(I)=XNEW(IY)                                                    GNE08940
      ZT(I)=XNEW(IZ)                                                    GNE08950
  103 CONTINUE                                                          GNE08960
C                                                                       GNE08970
C   CALCULATE INTERATOMIC DISTANCES FOR NEXT GEOMETRY                   GNE08980
      CALL DIST(NATOM,XT,YT,ZT,NDIM)                                    GNE08990
      IF(NATOM.LE.2) GO TO 202                                          GNE09000
      CALL ANGLE(NATOM,XT,YT,ZT,NDIM)                                   GNE09010
C                                                                       GNE09020
  202 CONTINUE                                                          GNE09030
      IF(IFLAG.EQ.0) GO TO 203                                          GNE09040
      DO 104 I=1,NATOM                                                  GNE09050
      XX(I)=XT(I)                                                       GNE09060
      YY(I)=YT(I)                                                       GNE09070
      ZZ(I)=ZT(I)                                                       GNE09080
  104 CONTINUE                                                          GNE09090
C                                                                       GNE09100
  203 CONTINUE                                                          GNE09110
      RETURN                                                            GNE09120
      END                                                               GNE09130
      SUBROUTINE METRIC(SX,SG,FX,HH,FT,EE,PIVOT,INDEX,NVAR,METHOD,      GNE09140
     1                  UPDATE)                                         GNE09150
      IMPLICIT REAL*8 (A-H,O-Z)                                         GNE09160
      LOGICAL UPDATE                                                    GNE09170
      DIMENSION SX(N3N,NPLUS),SG(N3N,NPLUS),FX(N3N,N3N)                 GNE09180
      DIMENSION HH(NVAR,NVAR),FT(NVAR,NVAR),EE(NVAR,NVAR)               GNE09190
      DIMENSION PIVOT(NVAR),INDEX(NVAR)                                 GNE09200
      COMMON/COM101/NATOM,N3N,NINT                                      GNE09210
      COMMON/COM102/NSORT,NCOOD,NDIF,NPLUS,NCASE                        GNE09220
      COMMON/COM103/CHG(50),XX(50),YY(50),ZZ(50)                        GNE09230
      COMMON/COM104/XT(150)                                             GNE09240
      COMMON/COM105/ENG0(150),ENG1(150)                                 GNE09250
      COMMON/COM106/ESCF(150),GNORM(150),SNORM(150)                     GNE09260
      COMMON/COM107/X(150),D(150),DELX(150),DELG(150),HG(150),DUMA(525) GNE09270
      COMMON/COM108/LDIM,MDIM,NDIM,IPRNT                                GNE09280
      COMMON/COM110/NSET(150),NORD(150),NXVAR(150)                      GNE09290
      COMMON/COM112/XOLD(150),XNEW(150)                                 GNE09300
      COMMON/COM113/STEPMX                                              GNE09310
      DATA ITAP15 / 15 /                                                GNE09320
      DATA A00,ONE / 0.0D+00 , 1.0D+00 /                                GNE09330
    1 FORMAT(//,2X,'********************************'/                  GNE09340
     1          2X,'***THE VARIABLE METRIC METHOD***'/                  GNE09350
     2          2X,'********************************'//)                GNE09360
    2 FORMAT(2I5)                                                       GNE09370
    3 FORMAT(3F20.10)                                                   GNE09380
    4 FORMAT(//,2X,' FX MATRIX'/)                                       GNE09390
    5 FORMAT(//,2X,' EXTRACTED FX MATRIX'/)                             GNE09400
    6 FORMAT(//,2X,' HH MATRIX'/)                                       GNE09410
    7 FORMAT(//,2X,' DETERM = ',D20.10/)                                GNE09420
    8 FORMAT(//,2X,' FX * FX(-1) = E'/)                                 GNE09430
    9 FORMAT(///,2X,' ALGOR CYCLE = ',I5)                               GNE09440
   10 FORMAT(//,2X,13H NEW X MATRIX/                                    GNE09450
     1 2X,4H NO.,5X,4H II.,6X,5H OLDX,18X,5H DELX,18X,5H NEWX/)         GNE09460
   11 FORMAT(3X,I2,4X,I5,3(3X,F20.10))                                  GNE09470
   12 FORMAT(//,2X,' GEOMETRY CHANGE IS TOO LARGE---TO BE SCALED'/      GNE09480
     1          2X,' STEPMX = ',F15.8/                                  GNE09490
     2          2X,' DD     = ',F15.8/                                  GNE09500
     3          2X,' SCALE  = ',F15.8/)                                 GNE09510
   13 FORMAT(//,2X,' NEW X MATRIX (SCLAED)'/                            GNE09520
     1 2X,4H NO.,5X,4H II.,6X,5H OLDX,18X,5H DELX,18X,5H NEWX/)         GNE09530
   14 FORMAT(//,2X,20H DELG AND HG VECTORS/)                            GNE09540
   15 FORMAT(2X,I2,2(5X,F15.8))                                         GNE09550
   16 FORMAT(//,2X,7H GHG = ,F15.8/                                     GNE09560
     1 2X,7H DG  = ,F15.8)                                              GNE09570
   17 FORMAT(//,2X,' UPDATED H MATRIX'/)                                GNE09580
   18 FORMAT(//,2X,' GEOMETRY CHANGE IS TOO LARGE---TRY AGAIN'/         GNE09590
     1          2X,' NALGOR,NMAX,NADD,XMAX = ',3I5,F10.5/)              GNE09600
   19 FORMAT(//,2X,' SUBROUTINE METRIC IS NOT ABLE TO FIND A RESONABLE  GNE09610
     1GEOMETRY GUESS'/)                                                 GNE09620
C                                                                       GNE09630
      WRITE(6,1)                                                        GNE09640
C                                                                       GNE09650
      NADD=0                                                            GNE09660
      NMAX=NPLUS                                                        GNE09670
  200 CONTINUE                                                          GNE09680
C***  IF(IHESS.EQ.0) GO TO 201                                          GNE09690
C***  REWIND ITAP15                                                     GNE09700
C***  READ(ITAP15,2) NNATOM,NNC                                         GNE09710
C***  READ(ITAP15,3) ((FX(I,J),J=1,N3N),I=1,N3N)                        GNE09720
C***  REWIND ITAP15                                                     GNE09730
C***  GO TO 202                                                         GNE09740
C                                                                       GNE09750
C   UNIT MATRIX AS AN INITIAL HESSIAN                                   GNE09760
  201 CONTINUE                                                          GNE09770
      DO 101 I=1,N3N                                                    GNE09780
      FX(I,I)=ONE                                                       GNE09790
      DO 101 J=1,N3N                                                    GNE09800
      IF(I.EQ.J) GO TO 101                                              GNE09810
      FX(I,J)=A00                                                       GNE09820
  101 CONTINUE                                                          GNE09830
C                                                                       GNE09840
  202 CONTINUE                                                          GNE09850
      IF(IPRNT.LE.0) GO TO 203                                          GNE09860
      WRITE(6,4)                                                        GNE09870
      CALL MATOUT(FX,N3N,N3N,N3N,N3N,6)                                 GNE09880
C                                                                       GNE09890
  203 CONTINUE                                                          GNE09900
      DO 102 I=1,NVAR                                                   GNE09910
      II=NXVAR(I)                                                       GNE09920
      DO 102 J=1,NVAR                                                   GNE09930
      JJ=NXVAR(J)                                                       GNE09940
      HH(I,J)=FX(II,JJ)                                                 GNE09950
      FT(I,J)=HH(I,J)                                                   GNE09960
  102 CONTINUE                                                          GNE09970
      IF(IPRNT.LE.0) GO TO 204                                          GNE09980
      WRITE(6,5)                                                        GNE09990
      CALL MATOUT(HH,NVAR,NVAR,NVAR,NVAR,6)                             GNE10000
  204 CONTINUE                                                          GNE10010
      CALL MATINV(HH,EE,PIVOT,INDEX,NVAR,NVAR,DETERM,NVAR)              GNE10020
      IF(IPRNT.LE.0) GO TO 205                                          GNE10030
      WRITE(6,6)                                                        GNE10040
      CALL MATOUT(HH,NVAR,NVAR,NVAR,NVAR,6)                             GNE10050
  205 WRITE(6,7) DETERM                                                 GNE10060
      IF(IPRNT.LE.0) GO TO 206                                          GNE10070
      CALL MTXMPY(FT,NVAR,HH,NVAR,EE,NVAR,EE,NVAR,NVAR,1)               GNE10080
      WRITE(6,8)                                                        GNE10090
      CALL MATOUT(EE,NVAR,NVAR,NVAR,NVAR,6)                             GNE10100
C                                                                       GNE10110
  206 CONTINUE                                                          GNE10120
      DO 103 I=1,N3N                                                    GNE10130
      X(I)=SX(I,NPLUS)                                                  GNE10140
      XT(I)=X(I)                                                        GNE10150
  103 CONTINUE                                                          GNE10160
C                                                                       GNE10170
      NALGOR=1                                                          GNE10180
      NPOS=NALGOR+NADD                                                  GNE10190
      DO 104 I=1,NVAR                                                   GNE10200
      II=NXVAR(I)                                                       GNE10210
      ENG0(I)=SG(II,NPOS)                                               GNE10220
  104 CONTINUE                                                          GNE10230
C                                                                       GNE10240
C   START ITERATION                                                     GNE10250
  300 CONTINUE                                                          GNE10260
      IF(IPRNT.LE.0) GO TO 301                                          GNE10270
      WRITE(6,9) NALGOR                                                 GNE10280
C                                                                       GNE10290
C   CALCULATE D VECTOR                                                  GNE10300
  301 DO 105 I=1,NVAR                                                   GNE10310
      D(I)=A00                                                          GNE10320
      DO 105 J=1,NVAR                                                   GNE10330
      D(I)=D(I)-HH(I,J)*ENG0(J)                                         GNE10340
  105 CONTINUE                                                          GNE10350
C                                                                       GNE10360
C   FIND ALPHA                                                          GNE10370
      CALL INTER(ALPHA,D,ENG0,MDIM,NVAR)                                GNE10380
C                                                                       GNE10390
C   ESTIMATE NEW COORDINATES                                            GNE10400
      IF(IPRNT.LE.0) GO TO 302                                          GNE10410
      WRITE(6,10)                                                       GNE10420
  302 CONTINUE                                                          GNE10430
      XMAX=A00                                                          GNE10440
      DO 106 I=1,NVAR                                                   GNE10450
      II=NXVAR(I)                                                       GNE10460
      X00=X(II)                                                         GNE10470
      DELX(I)=ALPHA*D(I)                                                GNE10480
      DMAX=DABS(DELX(I))                                                GNE10490
      IF(DMAX.GT.XMAX) XMAX=DMAX                                        GNE10500
      X11=X00+DELX(I)                                                   GNE10510
      XT(II)=X11                                                        GNE10520
      IF(IPRNT.LE.0) GO TO 106                                          GNE10530
      WRITE(6,11) I,II,X00,DELX(I),X11                                  GNE10540
  106 CONTINUE                                                          GNE10550
C                                                                       GNE10560
C   CHECK THE STEP SIZE                                                 GNE10570
      IF(XMAX.GT.STEPMX) THEN                                           GNE10580
        CALL ANORM(DD,DELX,150,NVAR)                                    GNE10590
        SCALE=STEPMX/DD                                                 GNE10600
        WRITE(6,12) STEPMX,DD,SCALE                                     GNE10610
        IF(IPRNT.GT.0) WRITE(6,13)                                      GNE10620
        XMAX=A00                                                        GNE10630
        DO 107 I=1,NVAR                                                 GNE10640
        II=NXVAR(I)                                                     GNE10650
        X00=X(II)                                                       GNE10660
        DELOLD=DELX(I)                                                  GNE10670
        DELX(I)=DELOLD*SCALE                                            GNE10680
        DMAX=DABS(DELX(I))                                              GNE10690
        IF(DMAX.GT.XMAX) XMAX=DMAX                                      GNE10700
        X11=X00+DELX(I)                                                 GNE10710
        XT(II)=X11                                                      GNE10720
        IF(IPRNT.LE.0) GO TO 107                                        GNE10730
        WRITE(6,11) I,II,X00,DELX(I),X11                                GNE10740
  107   CONTINUE                                                        GNE10750
      END IF                                                            GNE10760
      IF(XMAX.GT.STEPMX) GO TO 310                                      GNE10770
C                                                                       GNE10780
      NALGOR=NALGOR+1                                                   GNE10790
      NPOS=NALGOR+NADD                                                  GNE10800
      IF(NALGOR.GT.NMAX) GO TO 320                                      GNE10810
      DO 108 I=1,NVAR                                                   GNE10820
      II=NXVAR(I)                                                       GNE10830
      ENG1(I)=SG(II,NPOS)                                               GNE10840
  108 CONTINUE                                                          GNE10850
C                                                                       GNE10860
C   CALCULATION OF HH MATRIX FOR NEXT ITERATION                         GNE10870
      DO 109 I=1,NVAR                                                   GNE10880
      DELG(I)=ENG1(I)-ENG0(I)                                           GNE10890
  109 CONTINUE                                                          GNE10900
C                                                                       GNE10910
      IF(METHOD-2) 401,403,405                                          GNE10920
C                                                                       GNE10930
C::::::::::::::::::::::::::::                                           GNE10940
C:::MURTAGH-SARGENT METHOD:::                                           GNE10950
C::::::::::::::::::::::::::::                                           GNE10960
  401 DO 110 I=1,NVAR                                                   GNE10970
      HG(I)=DELX(I)                                                     GNE10980
      DO 110 J=1,NVAR                                                   GNE10990
      HG(I)=HG(I)-HH(I,J)*DELG(J)                                       GNE11000
  110 CONTINUE                                                          GNE11010
      IF(IPRNT.LE.0) GO TO 402                                          GNE11020
      WRITE(6,14)                                                       GNE11030
      DO 111 I=1,NVAR                                                   GNE11040
      WRITE(6,15) I,DELG(I),HG(I)                                       GNE11050
  111 CONTINUE                                                          GNE11060
  402 GHG=A00                                                           GNE11070
      DG=A00                                                            GNE11080
      DO 112 I=1,NVAR                                                   GNE11090
      GHG=GHG+HG(I)*DELG(I)                                             GNE11100
      DG=DG+HG(I)*HG(I)                                                 GNE11110
  112 CONTINUE                                                          GNE11120
      DO 113 I=1,NVAR                                                   GNE11130
      DO 113 J=1,NVAR                                                   GNE11140
      HH(I,J)=HH(I,J)+(HG(I)*HG(J))/GHG                                 GNE11150
  113 CONTINUE                                                          GNE11160
      IF(IPRNT.LE.0) GO TO 410                                          GNE11170
      WRITE(6,16) GHG,DG                                                GNE11180
      GO TO 410                                                         GNE11190
C                                                                       GNE11200
C:::::::::::::::::::::                                                  GNE11210
C:::FLETCHER METHOD:::                                                  GNE11220
C:::::::::::::::::::::                                                  GNE11230
  403 DO 120 I=1,NVAR                                                   GNE11240
      HG(I)=A00                                                         GNE11250
      DO 120 J=1,NVAR                                                   GNE11260
      HG(I)=HG(I)+HH(I,J)*DELG(J)                                       GNE11270
  120 CONTINUE                                                          GNE11280
      IF(IPRNT.LE.0) GO TO 404                                          GNE11290
      WRITE(6,14)                                                       GNE11300
      DO 121 I=1,NVAR                                                   GNE11310
      WRITE(6,15) I,DELG(I),HG(I)                                       GNE11320
  121 CONTINUE                                                          GNE11330
  404 GHG=A00                                                           GNE11340
      DG=A00                                                            GNE11350
      DO 122 I=1,NVAR                                                   GNE11360
      GHG=GHG+HG(I)*DELG(I)                                             GNE11370
      DG=DG+DELX(I)*DELG(I)                                             GNE11380
  122 CONTINUE                                                          GNE11390
      DO 123 I=1,NVAR                                                   GNE11400
      DO 123 J=1,NVAR                                                   GNE11410
      HH(I,J)=HH(I,J)-(HG(I)*DELX(J)+HG(J)*DELX(I))/DG                  GNE11420
     1 +(ONE+GHG/DG)*((DELX(I)*DELX(J))/DG)                             GNE11430
  123 CONTINUE                                                          GNE11440
      IF(IPRNT.LE.0) GO TO 410                                          GNE11450
      WRITE(6,16) GHG,DG                                                GNE11460
      GO TO 410                                                         GNE11470
C                                                                       GNE11480
C::::::::::::::::::::::::::::::::::::                                   GNE11490
C:::DAVIDON-FLETCHER-POWELL METHOD:::                                   GNE11500
C::::::::::::::::::::::::::::::::::::                                   GNE11510
  405 DO 130 I=1,NVAR                                                   GNE11520
      HG(I)=A00                                                         GNE11530
      DO 130 J=1,NVAR                                                   GNE11540
      HG(I)=HG(I)+HH(I,J)*DELG(J)                                       GNE11550
  130 CONTINUE                                                          GNE11560
      IF(IPRNT.LE.0) GO TO 406                                          GNE11570
      WRITE(6,14)                                                       GNE11580
      DO 131 I=1,NVAR                                                   GNE11590
      WRITE(6,15) I,DELG(I),HG(I)                                       GNE11600
  131 CONTINUE                                                          GNE11610
  406 GHG=A00                                                           GNE11620
      DG=A00                                                            GNE11630
      DO 132 I=1,NVAR                                                   GNE11640
      GHG=GHG+HG(I)*DELG(I)                                             GNE11650
      DG=DG+DELX(I)*DELG(I)                                             GNE11660
  132 CONTINUE                                                          GNE11670
      DO 133 I=1,NVAR                                                   GNE11680
      DO 133 J=1,NVAR                                                   GNE11690
      HH(I,J)=HH(I,J)-(HG(I)*HG(J))/GHG+(DELX(I)*DELX(J))/DG            GNE11700
  133 CONTINUE                                                          GNE11710
      IF(IPRNT.LE.0) GO TO 410                                          GNE11720
      WRITE(6,16) GHG,DG                                                GNE11730
C                                                                       GNE11740
C#####################                                                  GNE11750
C###UPDATE H MATRIX###                                                  GNE11760
C#####################                                                  GNE11770
  410 CONTINUE                                                          GNE11780
      IF(IPRNT.LE.0) GO TO 415                                          GNE11790
      WRITE(6,17)                                                       GNE11800
      CALL MATOUT(HH,NVAR,NVAR,NVAR,NVAR,6)                             GNE11810
C                                                                       GNE11820
  415 DO 135 I=1,NVAR                                                   GNE11830
      ENG0(I)=ENG1(I)                                                   GNE11840
  135 CONTINUE                                                          GNE11850
      GO TO 300                                                         GNE11860
C                                                                       GNE11870
C   GEOMETRY CHANGE IS TOO LARGE                                        GNE11880
  310 CONTINUE                                                          GNE11890
      WRITE(6,18) NALGOR,NMAX,NADD,XMAX                                 GNE11900
      NMAX=NMAX-1                                                       GNE11910
      NADD=NADD+1                                                       GNE11920
      IF(NMAX.LE.0) GO TO 315                                           GNE11930
      GO TO 200                                                         GNE11940
C                                                                       GNE11950
C   YOU CANNOT FIND NEXT GEOMETRY                                       GNE11960
  315 WRITE(6,19)                                                       GNE11970
      RETURN                                                            GNE11980
C                                                                       GNE11990
  320 CONTINUE                                                          GNE12000
      DO 136 I=1,N3N                                                    GNE12010
      XOLD(I)=SX(I,NPLUS)                                               GNE12020
      XNEW(I)=XT(I)                                                     GNE12030
  136 CONTINUE                                                          GNE12040
      IFLAG=0                                                           GNE12050
      IF(UPDATE) IFLAG=1                                                GNE12060
      CALL OLDNEW(XOLD,XNEW,IFLAG)                                      GNE12070
C                                                                       GNE12080
      RETURN                                                            GNE12090
      END                                                               GNE12100
      SUBROUTINE INTER(ALPHA,D,ENG0,NND,NN)                             GNE12110
      IMPLICIT REAL*8 (A-H,O-Z)                                         GNE12120
      DIMENSION D(NND),ENG0(NND)                                        GNE12130
      DATA A00,ONE,TWO / 0.0D+00 , 1.0D+00 , 2.0D+00 /                  GNE12140
      DATA FDEL / 1.0D+00 /                                             GNE12150
C                                                                       GNE12160
      GXS=A00                                                           GNE12170
      DO 101 I=1,NN                                                     GNE12180
      GXS=GXS+D(I)*ENG0(I)                                              GNE12190
  101 CONTINUE                                                          GNE12200
      FX0=-(TWO*FDEL)/GXS                                               GNE12210
      ALPHA=DMIN1(ONE,FX0)                                              GNE12220
      RETURN                                                            GNE12230
      END                                                               GNE12240
      SUBROUTINE MTXMPY(A,NAD,B,NBD,C,NCD,D,NDD,N,IMPY)                 GNE12250
      IMPLICIT REAL*8 (A-H,O-Z)                                         GNE12260
      DIMENSION A(NAD,NAD),B(NBD,NBD),C(NCD,NCD),D(NDD,NDD)             GNE12270
      DATA A00 / 0.0D+00 /                                              GNE12280
C                                                                       GNE12290
      DO 101 I=1,NCD                                                    GNE12300
      DO 101 J=1,NCD                                                    GNE12310
  101 C(I,J)=A00                                                        GNE12320
      GO TO (201,202,203,204,205,206),IMPY                              GNE12330
C     C=A*B                                                             GNE12340
  201 DO 102 I=1,N                                                      GNE12350
      DO 102 J=1,N                                                      GNE12360
      DO 102 K=1,N                                                      GNE12370
  102 C(I,J)=C(I,J)+A(I,K)*B(K,J)                                       GNE12380
      RETURN                                                            GNE12390
C     C=TA*B                                                            GNE12400
  202 DO 103 I=1,N                                                      GNE12410
      DO 103 J=1,N                                                      GNE12420
      DO 103 K=1,N                                                      GNE12430
  103 C(I,J)=C(I,J)+A(K,I)*B(K,J)                                       GNE12440
      RETURN                                                            GNE12450
C     C=A*TB                                                            GNE12460
  203 DO 104 I=1,N                                                      GNE12470
      DO 104 J=1,N                                                      GNE12480
      DO 104 K=1,N                                                      GNE12490
  104 C(I,J)=C(I,J)+A(I,K)*B(J,K)                                       GNE12500
      RETURN                                                            GNE12510
C     C=TA*TB                                                           GNE12520
  204 DO 105 I=1,N                                                      GNE12530
      DO 105 J=1,N                                                      GNE12540
      DO 105 K=1,N                                                      GNE12550
  105 C(I,J)=C(I,J)+A(K,I)*B(J,K)                                       GNE12560
      RETURN                                                            GNE12570
C     C=TA*B*A                                                          GNE12580
  205 DO 106 I=1,N                                                      GNE12590
      DO 106 J=1,N                                                      GNE12600
      D(I,J)=A00                                                        GNE12610
      DO 106 K=1,N                                                      GNE12620
  106 D(I,J)=D(I,J)+A(K,I)*B(K,J)                                       GNE12630
      DO 107 I=1,N                                                      GNE12640
      DO 107 J=1,N                                                      GNE12650
      DO 107 K=1,N                                                      GNE12660
  107 C(I,J)=C(I,J)+D(I,K)*A(K,J)                                       GNE12670
      RETURN                                                            GNE12680
C     C=A*B*TA                                                          GNE12690
  206 DO 108 I=1,N                                                      GNE12700
      DO 108 J=1,N                                                      GNE12710
      D(I,J)=A00                                                        GNE12720
      DO 108 K=1,N                                                      GNE12730
  108 D(I,J)=D(I,J)+A(I,K)*B(K,J)                                       GNE12740
      DO 109 I=1,N                                                      GNE12750
      DO 109 J=1,N                                                      GNE12760
      DO 109 K=1,N                                                      GNE12770
  109 C(I,J)=C(I,J)+D(I,K)*A(J,K)                                       GNE12780
      RETURN                                                            GNE12790
      END                                                               GNE12800
