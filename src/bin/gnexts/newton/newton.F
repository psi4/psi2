      subroutine fentry(cc,ia,maxcor)
cc    PROGRAM NEWTON                                                    NEW00010
C********************************************************************
C*  Last updated on November 26, 1997 by Dr. Yukio Yamaguchi        *
C*  Modification for Power compiler sereis                          *
C********************************************************************
C*  Last updated on October 04, 1991 by Dr. Yukio Yamaguchi         *
C*  Modification for the UNIX system                                *
C********************************************************************   NEW00020
C*  LAST UPDATED ON DECEMBER 07, 1990 BY YUKIO YAMAGUCHI            *   NEW00030
C********************************************************************   NEW00040
C   THE SOURCE CODE OF GEOMETRY OPTIMIZATION PROGRAM                    NEW00050
C   BASED ON THE NEWTON-RAPHSON METHOD                                  NEW00060
C   IN THE CARTESIAN COORDINATE                                         NEW00070
C   THIS PROGRAM IS ABLE TO TREAT UP TO '50' ATOMS                      NEW00080
C       WITH UP TO '150' INTERNAL COORDINATES                           NEW00090
      IMPLICIT REAL*8 (A-H,O-Z)                                         NEW00100
      dimension cc(maxcor),ia(maxcor*2)
cyy   DIMENSION CC(360000),IA(1)                                        NEW00110
      DIMENSION I30(200),NUNQ(50),NAT(50,50),COORD(3,50)                NEW00120
      CHARACTER*8 UPD                                                   NEW00130
      LOGICAL UPDATE                                                    NEW00140
      COMMON/COM101/NATOM,N3N,NINT                                      NEW00150
      COMMON/COM102/NSORT,NCOOD,NPLUS                                   NEW00160
      COMMON/COM103/CHG(50),XX(50),YY(50),ZZ(50)                        NEW00170
      COMMON/COM104/EX(50),EY(50),EZ(50)                                NEW00180
      COMMON/COM105/ENGR(150),RCOOD(150)                                NEW00190
      COMMON/COM106/ESCF(150),GNORM(150),SNORM(150)                     NEW00200
      COMMON/COM107/R(1275)                                             NEW00210
      COMMON/COM108/LDIM,MDIM,NDIM,IPRNT                                NEW00220
      COMMON/COM109/XMAX,YMAX,ZMAX                                      NEW00230
      COMMON/COM110/NSET(150),NORD(150),NXVAR(150)                      NEW00240
      COMMON/COM112/XOLD(150),XNEW(150)                                 NEW00250
      COMMON/COM113/STEPMX                                              NEW00260
cyy   EQUIVALENCE (CC,IA)                                               NEW00270
      DATA A00,TEN / 0.0D+00 , 10.0D+00 /                               NEW00280
      DATA FMILI / 1.0D-3 /                                             NEW00290
    1 FORMAT(//,2X,' THE GEOMETRY OPTIMIZATION PROGRAM'/                NEW00300
     1          2X,' BASED ON THE NEWTON-RAPHSON METHOD'/               NEW00310
     2          2X,' CARTESIAN COORDINATE SYSTEM'/                      NEW00320
     3          2X,' (VERSION OF 11/26/97)'//)                          NEW00330
    2 FORMAT(2I5,F20.10,I5)                                             NEW00340
    3 FORMAT(12I5)                                                      NEW00350
    4 FORMAT(14I5)                                                      NEW00360
    5 FORMAT(A5,3I5)                                                    NEW00370
    6 FORMAT(2X,' NCOOD   = ',I5/                                       NEW00380
     1 2X,' NSORT   = ',I5/                                             NEW00390
     2 2X,' IMETR   = ',I5/                                             NEW00400
     3 2X,' NVAR    = ',I5/                                             NEW00410
     4 2X,' ISTEP   = ',I5/                                             NEW00420
     5 2X,' IPRNT   = ',I5/                                             NEW00430
     6 2X,' NATOM   = ',I5/                                             NEW00440
     7 2X,' N3N     = ',I5/                                             NEW00450
     8 2X,' NINT    = ',I5/                                             NEW00460
     9 2X,' IFORCE  = ',I5/                                             NEW00470
     A 2X,' NSTORE  = ',I5/                                             NEW00480
     B 2X,' NPLUS   = ',I5/                                             NEW00490
     C 2X,' NCHNG   = ',I5/                                             NEW00500
     D 2X,' NCONV   = ',I5/                                             NEW00510
     E 2X,' NUNIQ   = ',I5)                                             NEW00520
    8 FORMAT(//,2X,' OLD COORDINATES IN TAPE30'/                        NEW00530
     1 39X,' OLD X IN A.U.'/                                            NEW00540
     2 4X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z/)               NEW00550
    9 FORMAT(3X,I5,1X,4F16.10)                                          NEW00560
   10 FORMAT(//,2X,' NEW COORDINATES IN TAPE30'/                        NEW00570
     1 39X,' NEW X IN A.U.'/                                            NEW00580
     2 4X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z/)               NEW00590
   11 FORMAT(//,2X,' GEOMETRY IN TAPE30 IS UPDATED'/)                   NEW00600
   12 FORMAT(//,2X,':::::::::::::::::::::::::::'/                       NEW00610
     1          2X,':::GEOMETRY IS OPTIMIZED:::'/                       NEW00620
     2          2X,':::::::::::::::::::::::::::'//                      NEW00630
     3          2X,' GNORM = ',F15.10,7X,' GCONV = ',F15.10/)           NEW00640
   13 FORMAT(//,2X,' TOO BIG CHANGE IN GEOMETRY---SCRATCH YOUR HEAD TO FNEW00650
     1IND A NEW GEOMETRY'//                                             NEW00660
     2          2X,' XMAX = ',F15.10,7X,' YMAX = ',F15.10,7X,' ZMAX = ',NEW00670
     3                        F15.10,7X,' GCHNG = ',F15.10/)            NEW00680
   14 FORMAT(//,2X,' SHORT OF INFORMATION---SCRATCH YOUR HEAD TO FIND A NEW00690
     1NEW GEOMETRY'/)                                                   NEW00700
C                                                                       NEW00710
      call psinit('APPEND')
      CALL TSTART(6)                                                    NEW00720
C* 10/04/91  CALL NOUNFL                                                NEW00730
C                                                                       NEW00740
      LDIM=1275                                                         NEW00750
      MDIM=150                                                          NEW00760
      NDIM=50                                                           NEW00770
      INPUT=5                                                           NEW00780
      ITAPE6=6                                                          NEW00790
      ITAP12=12                                                         NEW00800
      ITAP13=13                                                         NEW00810
      ITAP15=15                                                         NEW00820
      ITAP30=30                                                         NEW00830
      IFLAG=0                                                           NEW00840
      MAXCOR=360000                                                     NEW00850
C                                                                       NEW00860
      WRITE(6,1)                                                        NEW00870
C* Folowing three lines were added for the Unix Version (10/04/91)
      call ffile(itap12,'file12',0)
      call ffile(itap13,'file13',0)
      call ffile(itap15,'file15',0)
      REWIND ITAP12                                                     NEW00880
      REWIND ITAP13                                                     NEW00890
      READ(ITAP12,2) NATOM,NINT,ESCFZ,IFORCE                            NEW00900
C                                                                       NEW00910
C   NCOOD IS NUMBER OF INTERNAL COORDINATES TO BE OPTIMIZED             NEW00920
C   NSORT IS A PARAMETER FOR SORTING                                    NEW00930
C     NSORT=0---NO SORT                                                 NEW00940
C     NSORT=1---GRADIENT SORT                                           NEW00950
C     NSORT=2---ENERGY SORT                                             NEW00960
C     NSORT=3---SELECTED GRADIENT SORT                                  NEW00970
C   IMETR IS A PARAMETER FOR THE VARIABLE METRIC METHOD                 NEW00980
C     IMETR=0 OR 1---MURTAGH-SARGENT METHOD                             NEW00990
C     IMETR=2     ---FLETCHER METHOD                                    NEW01000
C     IMETR=3     ---DAVIDON-FLETCHER-POWELL METHOD                     NEW01010
C   NVAR IS NUMBER OF CARTESIAN COORDINATES TO BE OPTIMIZED             NEW01020
C   ISTEP IS A PARAMETER FOR MAXIMUM STEP SIZE                          NEW01030
C     ISTEP=0---0.1                                                     NEW01040
C     ISTEP=N---0.01*ISTEP                                              NEW01050
C   IDUM1,IDUM2,IDUM3, AND IDUM4 ARE DUMMY PARAMETERS TO ADJUST TO      NEW01060
C     INPUT FOR OTHER OPTIMIZATION PROGRAMS                             NEW01070
C   IPRNT IS A PARAMETER FOR OUTPUT                                     NEW01080
      N3N=NATOM*3                                                       NEW01090
      CALL LOCATE(INPUT,'# OPTIM ##',IERR)                              NEW01100
      READ(5,3) NCOOD,NSORT,IMETR,NVAR,ISTEP,IDUM1,IDUM2,IDUM3,IDUM4,   NEW01110
     1          IPRNT                                                   NEW01120
      READ(5,4) (NSET(I),I=1,NCOOD)                                     NEW01130
      DO 101 I=1,N3N                                                    NEW01140
  101 NXVAR(I)=I                                                        NEW01150
      IF(ISTEP.EQ.0) ISTEP=10                                           NEW01160
      IF(NVAR.EQ.0) NVAR=N3N                                            NEW01170
      IF(NVAR.EQ.N3N) GO TO 201                                         NEW01180
      READ(5,4) (NXVAR(I),I=1,NVAR)                                     NEW01190
  201 READ(ITAP13,3) NSTORE                                             NEW01200
      NPLUS=NSTORE+1                                                    NEW01210
      READ(5,5) UPD,NCHNG,NCONV,NUNIQ                                   NEW01220
      IF(UPD.EQ.'UP     ') UPDATE=.TRUE.                                NEW01230
      IF(NCHNG.EQ.0) NCHNG=2                                            NEW01240
      IF(NCONV.EQ.0) NCONV=7                                            NEW01250
      IF(NUNIQ.EQ.0) NUNIQ=NATOM                                        NEW01260
      WRITE(6,6) NCOOD,NSORT,IMETR,NVAR,ISTEP,IPRNT,                    NEW01270
     1 NATOM,N3N,NINT,IFORCE,NSTORE,NPLUS,NCHNG,NCONV,NUNIQ             NEW01280
      GCHNG=TEN**(-NCHNG)                                               NEW01290
      GCONV=TEN**(-NCONV)                                               NEW01300
      STEPMX=1.0D-02*DFLOAT(ISTEP)                                      NEW01310
C                                                                       NEW01320
      IF(.NOT.UPDATE) GO TO 202                                         NEW01330
      DO 102 I=1,NATOM                                                  NEW01340
      NUNQ(I)=1                                                         NEW01350
      NAT(I,1)=I                                                        NEW01360
  102 CONTINUE                                                          NEW01370
      IF(NUNIQ.EQ.NATOM) GO TO 202                                      NEW01380
      DO 103 I=1,NUNIQ                                                  NEW01390
      READ(5,4) NN,(NAT(I,J),J=1,NN)                                    NEW01400
      NUNQ(I)=NN                                                        NEW01410
  103 CONTINUE                                                          NEW01420
C                                                                       NEW01430
C***********************************************************************NEW01440
C   DYNAMIC ALLOCATION                                                  NEW01450
C***********************************************************************NEW01460
C   1 : SX(N3N,NPLUS)                                                   NEW01470
C   2 : SG(N3N,NPLUS)                                                   NEW01480
C   3 : SCOOD(NINT,NPLUS)                                               NEW01490
C   4 : SENG(NINT,NPLUS)                                                NEW01500
C***********************************************************************NEW01510
  202 CONTINUE                                                          NEW01520
      IC1=1                                                             NEW01530
      IC2=IC1+N3N*NPLUS                                                 NEW01540
      IC3=IC2+N3N*NPLUS                                                 NEW01550
      IC4=IC3+NINT*NPLUS                                                NEW01560
      ICSAVE=IC4+NINT*NPLUS                                             NEW01570
C***********************************************************************NEW01580
C                                                                       NEW01590
C   READ IN STORED INFORMATION FROM TAPE13                              NEW01600
      IC5=ICSAVE                                                        NEW01610
      IC6=IC5+N3N*NPLUS                                                 NEW01620
      IC7=IC6+N3N*NPLUS                                                 NEW01630
      IC8=IC7+NINT*NPLUS                                                NEW01640
      ICMAX=IC8+NINT*NPLUS                                              NEW01650
C.................X       G       COOD    ENG.....                      NEW01660
      CALL READIN(CC(IC5),CC(IC6),CC(IC7),CC(IC8))                      NEW01670
C                                                                       NEW01680
C   SORT DATA W.R.T. GRADIENT OR ENERGY                                 NEW01690
      IC5=ICSAVE                                                        NEW01700
      IC6=IC5+N3N*NPLUS                                                 NEW01710
      IC7=IC6+N3N*NPLUS                                                 NEW01720
      IC8=IC7+NINT*NPLUS                                                NEW01730
      ICMAX=IC8+NINT*NPLUS                                              NEW01740
C...............SX      SG      SCOOD   SENG    X       G.......        NEW01750
      CALL SORT(CC(IC1),CC(IC2),CC(IC3),CC(IC4),CC(IC5),CC(IC6),        NEW01760
C...............COOD    ENG.....                                        NEW01770
     1          CC(IC7),CC(IC8))                                        NEW01780
C                                                                       NEW01790
      GNORML=GNORM(NPLUS)*FMILI                                         NEW01800
      IF(NSORT.EQ.3) GNORML=SNORM(NPLUS)*FMILI                          NEW01810
      IF(GNORML.LT.GCONV) GO TO 215                                     NEW01820
C                                                                       NEW01830
C********************************                                       NEW01840
C***THE VARIABLE METRIC METHOD***                                       NEW01850
C********************************                                       NEW01860
C:::THIS SECTION ONLY HANDLES THE CARTESIAN COORDINATE SYSTEM:::        NEW01870
      IC5=ICSAVE                                                        NEW01880
      IC6=IC5+N3N*N3N                                                   NEW01890
      IC7=IC6+NVAR*NVAR                                                 NEW01900
      IC8=IC7+NVAR*NVAR                                                 NEW01910
      IC9=IC8+NVAR*NVAR                                                 NEW01920
      IC10=IC9+NVAR                                                     NEW01930
      IA10=IC10+IC10-1                                                  NEW01940
      ICMAX=IC10+NVAR                                                   NEW01950
C.................SX      SG      FX      HH      FT      EE......      NEW01960
      CALL METRIC(CC(IC1),CC(IC2),CC(IC5),CC(IC6),CC(IC7),CC(IC8),      NEW01970
C.................PIVOT   INDEX......................                   NEW01980
     1            CC(IC9),IA(IA10),NVAR,IMETR,UPDATE)                   NEW01990
C                                                                       NEW02000
C*******************************                                        NEW02010
C***UPDATE TAPE30 IF REQUIRED***                                        NEW02020
C*******************************                                        NEW02030
C                                                                       NEW02040
      IF(.NOT.UPDATE) GO TO 218                                         NEW02050
C                                                                       NEW02060
      IF(DABS(XMAX).GT.GCHNG) GO TO 217                                 NEW02070
      IF(DABS(YMAX).GT.GCHNG) GO TO 217                                 NEW02080
      IF(DABS(ZMAX).GT.GCHNG) GO TO 217                                 NEW02090
C                                                                       NEW02100
C   READ IN GEOMETRY FROM TAPE30                                        NEW02110
      CALL RFILE(ITAP30)                                                NEW02120
      CALL WREADW(ITAP30,I30,200,101,JUNK)                              NEW02130
      MPOINT=I30(2)                                                     NEW02140
      MCONST=I30(3)                                                     NEW02150
      NCALCS=I30(5)                                                     NEW02160
      IGEOP=100+MCONST+MPOINT+NCALCS                                    NEW02170
      CALL WREADW(ITAP30,LOCCAL,1,IGEOP,JUNK)                           NEW02180
      IGEOP=LOCCAL+80                                                   NEW02190
      CALL WREADW(ITAP30,COORD,N3N*2,IGEOP,JUNK)                        NEW02200
      WRITE(6,8)                                                        NEW02210
      DO 104 I=1,NATOM                                                  NEW02220
      WRITE(6,9) I,CHG(I),COORD(1,I),COORD(2,I),COORD(3,I)              NEW02230
  104 CONTINUE                                                          NEW02240
      DO 107 I=1,NUNIQ                                                  NEW02250
      NN=NUNQ(I)                                                        NEW02260
      SUMX=A00                                                          NEW02270
      SUMY=A00                                                          NEW02280
      SUMZ=A00                                                          NEW02290
      DO 105 J=1,NN                                                     NEW02300
      MM=NAT(I,J)                                                       NEW02310
      SUMX=SUMX+DABS(XX(MM))                                            NEW02320
      SUMY=SUMY+DABS(YY(MM))                                            NEW02330
      SUMZ=SUMZ+DABS(ZZ(MM))                                            NEW02340
  105 CONTINUE                                                          NEW02350
      ANX=DBLE(NN)                                                      NEW02360
      AVGX=SUMX/ANX                                                     NEW02370
      AVGY=SUMY/ANX                                                     NEW02380
      AVGZ=SUMZ/ANX                                                     NEW02390
      DO 106 J=1,NN                                                     NEW02400
      MM=NAT(I,J)                                                       NEW02410
      XX(MM)=DSIGN(AVGX,XX(MM))                                         NEW02420
      YY(MM)=DSIGN(AVGY,YY(MM))                                         NEW02430
      ZZ(MM)=DSIGN(AVGZ,ZZ(MM))                                         NEW02440
  106 CONTINUE                                                          NEW02450
  107 CONTINUE                                                          NEW02460
      WRITE(6,10)                                                       NEW02470
      DO 108 I=1,NATOM                                                  NEW02480
      COORD(1,I)=XX(I)                                                  NEW02490
      COORD(2,I)=YY(I)                                                  NEW02500
      COORD(3,I)=ZZ(I)                                                  NEW02510
      WRITE(6,9) I,CHG(I),COORD(1,I),COORD(2,I),COORD(3,I)              NEW02520
  108 CONTINUE                                                          NEW02530
      CALL WWRITW(ITAP30,COORD,N3N*2,IGEOP,JUNK)                        NEW02540
      WRITE(6,11)                                                       NEW02550
      CALL RCLOSE(ITAP30,3)                                             NEW02560
      GO TO 220                                                         NEW02570
C                                                                       NEW02580
C   GEOMETRY IS OPTIMIZED                                               NEW02590
  215 CONTINUE                                                          NEW02600
      WRITE(6,12) GNORML,GCONV                                          NEW02610
      GO TO 220                                                         NEW02620
C                                                                       NEW02630
C   GEOMETRY CHANGE IS TOO LARGE                                        NEW02640
  217 CONTINUE                                                          NEW02650
      WRITE(6,13) XMAX,YMAX,ZMAX,GCHNG                                  NEW02660
      GO TO 220                                                         NEW02670
C                                                                       NEW02680
C   SHORT OF INFORMATION                                                NEW02690
  218 CONTINUE                                                          NEW02700
      WRITE(6,14)                                                       NEW02710
C                                                                       NEW02720
  220 CONTINUE                                                          NEW02730
C                                                                       NEW02740
      CALL TSTOP(6)                                                     NEW02750
      STOP                                                              NEW02760
      END                                                               NEW02770
      SUBROUTINE READIN(X,G,COOD,ENG)                                   NEW02780
      IMPLICIT REAL*8 (A-H,O-Z)                                         NEW02790
      DIMENSION X(N3N,NPLUS),G(N3N,NPLUS)                               NEW02800
      DIMENSION COOD(NINT,NPLUS),ENG(NINT,NPLUS)                        NEW02810
      COMMON/COM101/NATOM,N3N,NINT                                      NEW02820
      COMMON/COM102/NSORT,NCOOD,NPLUS                                   NEW02830
      COMMON/COM103/CHG(50),XX(50),YY(50),ZZ(50)                        NEW02840
      COMMON/COM104/EX(50),EY(50),EZ(50)                                NEW02850
      COMMON/COM105/ENGR(150),RCOOD(150)                                NEW02860
      COMMON/COM106/ESCF(150),GNORM(150),SNORM(150)                     NEW02870
      COMMON/COM108/LDIM,MDIM,NDIM,IPRNT                                NEW02880
      COMMON/COM110/NSET(150),NORD(150),NXVAR(150)                      NEW02890
      DATA ITAP12,ITAP13 / 12 , 13 /                                    NEW02900
      DATA FKILO,FMILI / 1.0D+3 , 1.0D-3 /                              NEW02910
    1 FORMAT(2I5,F20.10,I5)                                             NEW02920
    2 FORMAT(10I5)                                                      NEW02930
    3 FORMAT(I5,F20.10)                                                 NEW02940
    4 FORMAT(////,2X,10H CASE   = ,I5/                                  NEW02950
     1 2X,10H ESCF   = ,F20.10)                                         NEW02960
    5 FORMAT(//,2X,22H CARTESIAN COORDINATES,51X,10H GRADIENTS/         NEW02970
     1 4X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z,                NEW02980
     2                       17X,2H X,14X,2H Y,14X,2H Z)                NEW02990
    6 FORMAT(4F20.10)                                                   NEW03000
    7 FORMAT(3X,I5,1X,4F16.10,3X,3F16.10)                               NEW03010
    8 FORMAT(//,2X,' INTERNAL GRADIENTS (IN 1.0D+3 A.U.)'/              NEW03020
     1 5X,4H NO.,10X,7H COORD.,13X,6H GRAD./)                           NEW03030
    9 FORMAT(2F20.10)                                                   NEW03040
   10 FORMAT(2X,I7,5X,2F20.10)                                          NEW03050
   11 FORMAT(/,2X,22H NORM FOR GRADIENTS = ,F15.7)                      NEW03060
   12 FORMAT(//,2X,' ***SELECTED VARIABLES DO NOT MATCH***'/            NEW03070
     1 2X,' NCOOD = ',I5,' II = ',I5/)                                  NEW03080
   13 FORMAT(/,2X,' NORM FOR CASE ',I5,3X,' IS ',F15.7/)                NEW03090
   14 FORMAT(///,2X,' CARTESIAN COORDINATE MATRIX'/)                    NEW03100
   15 FORMAT(//,2X,' INTERNAL COORDINATE MATRIX'/)                      NEW03110
   16 FORMAT(//,2X,' CARTESIAN GRADIENT MATRIX'/)                       NEW03120
   17 FORMAT(//,2X,' INTERNAL GRADINET MATRIX'/)                        NEW03130
   18 FORMAT(//,2X,' GRADIENTS NORMS AND ENERGIES'/                     NEW03140
     1 3X,4H NO.,8X,6H GNORM,14X,6H SNORM,14X,7H ENERGY/)               NEW03150
   19 FORMAT(2X,I4,5X,3F20.10)                                          NEW03160
C                                                                       NEW03170
      REWIND ITAP12                                                     NEW03180
      READ(ITAP12,1) NATOM,NINT,ESCFZ,IFORCE                            NEW03190
      REWIND ITAP13                                                     NEW03200
      READ(ITAP13,2) NSTORE                                             NEW03210
C                                                                       NEW03220
      DO 105 ICASE=1,NPLUS                                              NEW03230
      IF(ICASE.EQ.NPLUS) GO TO 202                                      NEW03240
      READ(ITAP13,3) JCASE,ESCF(ICASE)                                  NEW03250
      GO TO 203                                                         NEW03260
  202 ESCF(ICASE)=ESCFZ                                                 NEW03270
  203 IF(IPRNT.LE.0) GO TO 301                                          NEW03280
      WRITE(6,4) ICASE,ESCF(ICASE)                                      NEW03290
  301 CONTINUE                                                          NEW03300
      IF(IPRNT.LE.0) GO TO 302                                          NEW03310
      WRITE(6,5)                                                        NEW03320
C***READ CARTESIAN COORDINATES OF EACH CASE                             NEW03330
  302 DO 102 I=1,NATOM                                                  NEW03340
      IF(ICASE.EQ.NPLUS) GO TO 204                                      NEW03350
      READ(ITAP13,6) CHG(I),XX(I),YY(I),ZZ(I)                           NEW03360
      GO TO 102                                                         NEW03370
  204 READ(ITAP12,6) CHG(I),XX(I),YY(I),ZZ(I)                           NEW03380
  102 CONTINUE                                                          NEW03390
C***READ CARTESIAN GRADIENTS OF EACH CASE                               NEW03400
      DO 103 I=1,NATOM                                                  NEW03410
      IF(ICASE.EQ.NPLUS) GO TO 205                                      NEW03420
      READ(ITAP13,6) CHG(I),EX(I),EY(I),EZ(I)                           NEW03430
      GO TO 206                                                         NEW03440
  205 READ(ITAP12,6) CHG(I),EX(I),EY(I),EZ(I)                           NEW03450
  206 IF(IPRNT.LE.0) GO TO 303                                          NEW03460
      WRITE(6,7) I,CHG(I),XX(I),YY(I),ZZ(I),EX(I),EY(I),EZ(I)           NEW03470
  303 IX=(I-1)*3+1                                                      NEW03480
      IY=IX+1                                                           NEW03490
      IZ=IX+2                                                           NEW03500
      X(IX,ICASE)=XX(I)                                                 NEW03510
      X(IY,ICASE)=YY(I)                                                 NEW03520
      X(IZ,ICASE)=ZZ(I)                                                 NEW03530
      G(IX,ICASE)=EX(I)                                                 NEW03540
      G(IY,ICASE)=EY(I)                                                 NEW03550
      G(IZ,ICASE)=EZ(I)                                                 NEW03560
  103 CONTINUE                                                          NEW03570
      IF(IPRNT.LE.0) GO TO 304                                          NEW03580
      CALL DIST(NATOM,XX,YY,ZZ,NDIM)                                    NEW03590
      WRITE(6,8)                                                        NEW03600
C***READ COORDINATES AND GRADIENTS OF EACH CASE                         NEW03610
  304 DO 104 I=1,NINT                                                   NEW03620
      IF(ICASE.EQ.NPLUS) GO TO 207                                      NEW03630
      READ(ITAP13,9) RCOOD(I),ENGR(I)                                   NEW03640
      GO TO 208                                                         NEW03650
  207 READ(ITAP12,9) RCOOD(I),ENGR(I)                                   NEW03660
  208 ENGR(I)=ENGR(I)*FKILO                                             NEW03670
      IF(IPRNT.LE.0) GO TO 305                                          NEW03680
      WRITE(6,10) I,RCOOD(I),ENGR(I)                                    NEW03690
  305 COOD(I,ICASE)=RCOOD(I)                                            NEW03700
      ENG(I,ICASE)=ENGR(I)                                              NEW03710
  104 CONTINUE                                                          NEW03720
C   CALCULATE NORM OF GRADIENTS                                         NEW03730
      CALL ANORM(CNORM,ENGR,MDIM,NINT)                                  NEW03740
      IF(IPRNT.LE.0) GO TO 306                                          NEW03750
      WRITE(6,11) CNORM                                                 NEW03760
  306 GNORM(ICASE)=CNORM                                                NEW03770
  105 CONTINUE                                                          NEW03780
C                                                                       NEW03790
      DO 108 ICASE=1,NPLUS                                              NEW03800
      II=0                                                              NEW03810
      DO 107 I=1,NINT                                                   NEW03820
      DO 106 J=1,NCOOD                                                  NEW03830
      IF(NSET(J).NE.I) GO TO 106                                        NEW03840
      II=II+1                                                           NEW03850
      ENGR(II)=ENG(I,ICASE)                                             NEW03860
  106 CONTINUE                                                          NEW03870
  107 CONTINUE                                                          NEW03880
      IF(II.EQ.NCOOD) GO TO 307                                         NEW03890
      WRITE(6,12) NCOOD,II                                              NEW03900
      STOP                                                              NEW03910
  307 CONTINUE                                                          NEW03920
      CALL ANORM(CNORM,ENGR,MDIM,NCOOD)                                 NEW03930
      IF(IPRNT.LE.0) GO TO 308                                          NEW03940
      WRITE(6,13) ICASE,CNORM                                           NEW03950
  308 SNORM(ICASE)=CNORM                                                NEW03960
  108 CONTINUE                                                          NEW03970
C                                                                       NEW03980
C   SUMMARY OF INPUT                                                    NEW03990
      IF(IPRNT.LE.0) GO TO 309                                          NEW04000
      WRITE(6,14)                                                       NEW04010
      CALL MATOUT(X,N3N,NPLUS,N3N,NPLUS,6)                              NEW04020
      WRITE(6,15)                                                       NEW04030
      CALL MATOUT(COOD,NINT,NPLUS,NINT,NPLUS,6)                         NEW04040
      WRITE(6,16)                                                       NEW04050
      CALL MATOUT(G,N3N,NPLUS,N3N,NPLUS,6)                              NEW04060
      WRITE(6,17)                                                       NEW04070
      CALL MATOUT(ENG,NINT,NPLUS,NINT,NPLUS,6)                          NEW04080
  309 WRITE(6,18)                                                       NEW04090
      DO 110 ICASE=1,NPLUS                                              NEW04100
      WRITE(6,19) ICASE,GNORM(ICASE),SNORM(ICASE),ESCF(ICASE)           NEW04110
  110 CONTINUE                                                          NEW04120
C                                                                       NEW04130
C   STORE GEOMETRIES AND GRADIENTS IN TAPE13                            NEW04140
      REWIND ITAP13                                                     NEW04150
      WRITE(ITAP13,2) NPLUS                                             NEW04160
      DO 115 ICASE=1,NPLUS                                              NEW04170
      WRITE(ITAP13,3) ICASE,ESCF(ICASE)                                 NEW04180
      DO 111 I=1,NATOM                                                  NEW04190
      IX=(I-1)*3+1                                                      NEW04200
      IY=IX+1                                                           NEW04210
      IZ=IX+2                                                           NEW04220
      XX(I)=X(IX,ICASE)                                                 NEW04230
      YY(I)=X(IY,ICASE)                                                 NEW04240
      ZZ(I)=X(IZ,ICASE)                                                 NEW04250
      WRITE(ITAP13,6) CHG(I),XX(I),YY(I),ZZ(I)                          NEW04260
  111 CONTINUE                                                          NEW04270
      DO 112 I=1,NATOM                                                  NEW04280
      IX=(I-1)*3+1                                                      NEW04290
      IY=IX+1                                                           NEW04300
      IZ=IX+2                                                           NEW04310
      EX(I)=G(IX,ICASE)                                                 NEW04320
      EY(I)=G(IY,ICASE)                                                 NEW04330
      EZ(I)=G(IZ,ICASE)                                                 NEW04340
      WRITE(ITAP13,6) CHG(I),EX(I),EY(I),EZ(I)                          NEW04350
  112 CONTINUE                                                          NEW04360
      DO 113 I=1,NINT                                                   NEW04370
      RCOOD(I)=COOD(I,ICASE)                                            NEW04380
      ENGR(I)=ENG(I,ICASE)*FMILI                                        NEW04390
      WRITE(ITAP13,9) RCOOD(I),ENGR(I)                                  NEW04400
  113 CONTINUE                                                          NEW04410
  115 CONTINUE                                                          NEW04420
C                                                                       NEW04430
      REWIND ITAP12                                                     NEW04440
      REWIND ITAP13                                                     NEW04450
C                                                                       NEW04460
      RETURN                                                            NEW04470
      END                                                               NEW04480
      SUBROUTINE DIST(NATOM,XX,YY,ZZ,NNX)                               NEW04490
      IMPLICIT REAL*8 (A-H,O-Z)                                         NEW04500
      DIMENSION XX(NNX),YY(NNX),ZZ(NNX)                                 NEW04510
      COMMON/COM107/R(1275)                                             NEW04520
      COMMON/COM108/LDIM,MDIM,NDIM,IPRNT                                NEW04530
      DATA A00 / 0.0D+00 /                                              NEW04540
    1 FORMAT(//,2X,28H INTERATOMIC DISTANCE MATRIX/)                    NEW04550
C                                                                       NEW04560
      IJ=0                                                              NEW04570
      DO 101 I=1,NATOM                                                  NEW04580
      DO 101 J=1,I                                                      NEW04590
      IJ=IJ+1                                                           NEW04600
      R(IJ)=A00                                                         NEW04610
      IF(I.EQ.J) GO TO 101                                              NEW04620
      XD=XX(I)-XX(J)                                                    NEW04630
      YD=YY(I)-YY(J)                                                    NEW04640
      ZD=ZZ(I)-ZZ(J)                                                    NEW04650
      RR=DSQRT(XD*XD+YD*YD+ZD*ZD)                                       NEW04660
      R(IJ)=RR                                                          NEW04670
  101 CONTINUE                                                          NEW04680
      WRITE(6,1)                                                        NEW04690
      CALL PRINT(R,LDIM,NATOM,6)                                        NEW04700
      RETURN                                                            NEW04710
      END                                                               NEW04720
      SUBROUTINE ANGLE(NATOM,XX,YY,ZZ,NNX)                              NEW04730
      IMPLICIT REAL*8 (A-H,O-Z)                                         NEW04740
      DIMENSION XX(NNX),YY(NNX),ZZ(NNX)                                 NEW04750
      COMMON/COM107/R(1275)                                             NEW04760
      COMMON/COM108/LDIM,MDIM,NDIM,IPRNT                                NEW04770
      DATA PI / 3.1415926536D+00 /                                      NEW04780
    1 FORMAT(//,2X,' BOND ANGLES'/                                      NEW04790
     1 2X,' I   J   K',3X,' DEGREES   ',1X,' J   I   K',3X,' DEGREES  ',NEW04800
     2 1X,' I   K   J',3X,' DEGREES   '/)                               NEW04810
    2 FORMAT(3I4,F11.3, 2X,3I4,F11.3,1X,3I4,F11.3)                      NEW04820
    3 FORMAT(//,2X,' BOND ANGLES WILL NOT BE PRINTED---IPRNT .LT. 0 '/) NEW04830
C                                                                       NEW04840
      AR=PI/180.0D+00                                                   NEW04850
      IF(IPRNT.GE.0) WRITE(6,1)                                         NEW04860
      IF(IPRNT.LT.0) WRITE(6,3)                                         NEW04870
      DO 101 I=1,NATOM-2                                                NEW04880
      DO 101 J=I+1,NATOM-1                                              NEW04890
      IJ=I+J*(J-1)/2                                                    NEW04900
      XIJ=XX(I)-XX(J)                                                   NEW04910
      YIJ=YY(I)-YY(J)                                                   NEW04920
      ZIJ=ZZ(I)-ZZ(J)                                                   NEW04930
      EXIJ=-XIJ/R(IJ)                                                   NEW04940
      EYIJ=-YIJ/R(IJ)                                                   NEW04950
      EZIJ=-ZIJ/R(IJ)                                                   NEW04960
      EXJI=-EXIJ                                                        NEW04970
      EYJI=-EYIJ                                                        NEW04980
      EZJI=-EZIJ                                                        NEW04990
      DO 101 K=J+1,NATOM                                                NEW05000
      IK=I+K*(K-1)/2                                                    NEW05010
      XIK=XX(I)-XX(K)                                                   NEW05020
      YIK=YY(I)-YY(K)                                                   NEW05030
      ZIK=ZZ(I)-ZZ(K)                                                   NEW05040
      EXIK=-XIK/R(IK)                                                   NEW05050
      EYIK=-YIK/R(IK)                                                   NEW05060
      EZIK=-ZIK/R(IK)                                                   NEW05070
      EXKI=-EXIK                                                        NEW05080
      EYKI=-EYIK                                                        NEW05090
      EZKI=-EZIK                                                        NEW05100
      JK=J+K*(K-1)/2                                                    NEW05110
      XJK=XX(J)-XX(K)                                                   NEW05120
      YJK=YY(J)-YY(K)                                                   NEW05130
      ZJK=ZZ(J)-ZZ(K)                                                   NEW05140
      EXJK=-XJK/R(JK)                                                   NEW05150
      EYJK=-YJK/R(JK)                                                   NEW05160
      EZJK=-ZJK/R(JK)                                                   NEW05170
      EXKJ=-EXJK                                                        NEW05180
      EYKJ=-EYJK                                                        NEW05190
      EZKJ=-EZJK                                                        NEW05200
      AIJK=EXJI*EXJK+EYJI*EYJK+EZJI*EZJK                                NEW05210
      AJIK=EXIJ*EXIK+EYIJ*EYIK+EZIJ*EZIK                                NEW05220
      AIKJ=EXKI*EXKJ+EYKI*EYKJ+EZKI*EZKJ                                NEW05230
      ANGIJK=DACOS(AIJK)/AR                                             NEW05240
      ANGJIK=DACOS(AJIK)/AR                                             NEW05250
      ANGIKJ=DACOS(AIKJ)/AR                                             NEW05260
      IF(IPRNT.GE.0) WRITE(6,2) I,J,K,ANGIJK,J,I,K,ANGJIK,I,K,J,ANGIKJ  NEW05270
  101 CONTINUE                                                          NEW05280
C                                                                       NEW05290
      RETURN                                                            NEW05300
      END                                                               NEW05310
      SUBROUTINE ANORM(CNORM,C,NNC,NC)                                  NEW05320
      IMPLICIT REAL*8 (A-H,O-Z)                                         NEW05330
      DIMENSION C(NNC)                                                  NEW05340
      DATA A00 / 0.0D+00 /                                              NEW05350
C                                                                       NEW05360
      CNORM=A00                                                         NEW05370
      DO 101 I=1,NC                                                     NEW05380
      CNORM=CNORM+C(I)*C(I)                                             NEW05390
  101 CONTINUE                                                          NEW05400
      CNORM=DSQRT(CNORM)                                                NEW05410
      RETURN                                                            NEW05420
      END                                                               NEW05430
      SUBROUTINE SORT(SX,SG,SCOOD,SENG,X,G,COOD,ENG)                    NEW05440
      IMPLICIT REAL*8 (A-H,O-Z)                                         NEW05450
      DIMENSION SX(N3N,NPLUS),SG(N3N,NPLUS),X(N3N,NPLUS),G(N3N,NPLUS)   NEW05460
      DIMENSION SCOOD(NINT,NPLUS),SENG(NINT,NPLUS)                      NEW05470
      DIMENSION COOD(NINT,NPLUS),ENG(NINT,NPLUS)                        NEW05480
      COMMON/COM101/NATOM,N3N,NINT                                      NEW05490
      COMMON/COM102/NSORT,NCOOD,NPLUS                                   NEW05500
      COMMON/COM106/ESCF(150),GNORM(150),SNORM(150)                     NEW05510
      COMMON/COM108/LDIM,MDIM,NDIM,IPRNT                                NEW05520
      COMMON/COM110/NSET(150),NORD(150),NXVAR(150)                      NEW05530
    1 FORMAT(//,2X,' SORTED CARTESIAN COORDINATE MATRIX'/)              NEW05540
    2 FORMAT(//,2X,' SORTED CARTESIAN GRADIENT MATRIX'/)                NEW05550
    3 FORMAT(//,2X,' SORTED INTERNAL COORDINATE MATRIX'/)               NEW05560
    4 FORMAT(//,2X,' SORTED INTERNAL GRADIENT MATRIX'/)                 NEW05570
C                                                                       NEW05580
      DO 101 I=1,NPLUS                                                  NEW05590
      NORD(I)=I                                                         NEW05600
  101 CONTINUE                                                          NEW05610
      IF(NPLUS.EQ.1) GO TO 205                                          NEW05620
      IF(NSORT.EQ.0) GO TO 205                                          NEW05630
      GO TO (201,202,203),NSORT                                         NEW05640
C                                                                       NEW05650
C   GRADIENT SORT                                                       NEW05660
  201 CALL ORDER(GNORM,NORD,MDIM,NPLUS)                                 NEW05670
      GO TO 205                                                         NEW05680
C                                                                       NEW05690
C   ENERGY SORT                                                         NEW05700
  202 CALL ORDER(ESCF,NORD,MDIM,NPLUS)                                  NEW05710
      GO TO 205                                                         NEW05720
C                                                                       NEW05730
C   SELECTED GRADIENT SORT                                              NEW05740
  203 CALL ORDER(SNORM,NORD,MDIM,NPLUS)                                 NEW05750
C                                                                       NEW05760
C   NO SORT AND/OR AFTER SORT                                           NEW05770
  205 II=0                                                              NEW05780
      DO 104 I=1,NPLUS                                                  NEW05790
      II=II+1                                                           NEW05800
      MM=NORD(I)                                                        NEW05810
      DO 102 J=1,NINT                                                   NEW05820
      SCOOD(J,II)=COOD(J,MM)                                            NEW05830
      SENG(J,II)=ENG(J,MM)                                              NEW05840
  102 CONTINUE                                                          NEW05850
      DO 103 J=1,N3N                                                    NEW05860
      SX(J,II)=X(J,MM)                                                  NEW05870
      SG(J,II)=G(J,MM)                                                  NEW05880
  103 CONTINUE                                                          NEW05890
  104 CONTINUE                                                          NEW05900
      IF(IPRNT.LE.0) GO TO 206                                          NEW05910
      WRITE(6,1)                                                        NEW05920
      CALL MATOUT(SX,N3N,NPLUS,N3N,NPLUS,6)                             NEW05930
      WRITE(6,2)                                                        NEW05940
      CALL MATOUT(SG,N3N,NPLUS,N3N,NPLUS,6)                             NEW05950
  206 CONTINUE                                                          NEW05960
      WRITE(6,3)                                                        NEW05970
      CALL MATOUT(SCOOD,NINT,NPLUS,NINT,NPLUS,6)                        NEW05980
      WRITE(6,4)                                                        NEW05990
      CALL MATOUT(SENG,NINT,NPLUS,NINT,NPLUS,6)                         NEW06000
C                                                                       NEW06010
      RETURN                                                            NEW06020
      END                                                               NEW06030
      SUBROUTINE ORDER(A,NORD,NND,N)                                    NEW06040
      IMPLICIT REAL*8 (A-H,O-Z)                                         NEW06050
      DIMENSION A(NND),NORD(NND)                                        NEW06060
C                                                                       NEW06070
      N1=N-1                                                            NEW06080
      DO 102 I=1,N1                                                     NEW06090
      IR=I+1                                                            NEW06100
      DO 101 J=IR,N                                                     NEW06110
      IF(A(I).GT.A(J)) GO TO 101                                        NEW06120
      AA=A(I)                                                           NEW06130
      A(I)=A(J)                                                         NEW06140
      A(J)=AA                                                           NEW06150
      NN=NORD(I)                                                        NEW06160
      NORD(I)=NORD(J)                                                   NEW06170
      NORD(J)=NN                                                        NEW06180
  101 CONTINUE                                                          NEW06190
  102 CONTINUE                                                          NEW06200
      RETURN                                                            NEW06210
      END                                                               NEW06220
      SUBROUTINE OLDNEW(XOLD,XNEW,IFLAG)                                NEW06230
      IMPLICIT REAL*8 (A-H,O-Z)                                         NEW06240
      DIMENSION XOLD(150),XNEW(150)                                     NEW06250
      COMMON/COM101/NATOM,N3N,NINT                                      NEW06260
      COMMON/COM103/CHG(50),XX(50),YY(50),ZZ(50)                        NEW06270
      COMMON/COM104/XT(50),YT(50),ZT(50)                                NEW06280
      COMMON/COM108/LDIM,MDIM,NDIM,IPRNT                                NEW06290
      COMMON/COM109/XMAX,YMAX,ZMAX                                      NEW06300
      DATA A00 / 0.0D+00 /                                              NEW06310
    1 FORMAT(//,2X,' OLD COORDINATES'/                                  NEW06320
     1 38X,' OLD X IN A.U.'/                                            NEW06330
     2 3X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z/)               NEW06340
    2 FORMAT(3X,I5,4F16.10)                                             NEW06350
    3 FORMAT(//,2X,' NEW COORDINATES'/                                  NEW06360
     1 38X,' NEW X IN A.U.'/                                            NEW06370
     2 3X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z/)               NEW06380
C                                                                       NEW06390
      WRITE(6,1)                                                        NEW06400
      XMAX=A00                                                          NEW06410
      YMAX=A00                                                          NEW06420
      ZMAX=A00                                                          NEW06430
      DO 101 I=1,NATOM                                                  NEW06440
      IX=(I-1)*3+1                                                      NEW06450
      IY=IX+1                                                           NEW06460
      IZ=IX+2                                                           NEW06470
      IF(IFLAG.EQ.0) GO TO 201                                          NEW06480
      XDIF=XNEW(IX)-XOLD(IX)                                            NEW06490
      YDIF=XNEW(IY)-XOLD(IY)                                            NEW06500
      ZDIF=XNEW(IZ)-XOLD(IZ)                                            NEW06510
      IF(DABS(XDIF).GT.DABS(XMAX)) XMAX=XDIF                            NEW06520
      IF(DABS(YDIF).GT.DABS(YMAX)) YMAX=YDIF                            NEW06530
      IF(DABS(ZDIF).GT.DABS(ZMAX)) ZMAX=ZDIF                            NEW06540
  201 CONTINUE                                                          NEW06550
      WRITE(6,2) I,CHG(I),XOLD(IX),XOLD(IY),XOLD(IZ)                    NEW06560
  101 CONTINUE                                                          NEW06570
      WRITE(6,3)                                                        NEW06580
      DO 102 I=1,NATOM                                                  NEW06590
      IX=(I-1)*3+1                                                      NEW06600
      IY=IX+1                                                           NEW06610
      IZ=IX+2                                                           NEW06620
      WRITE(6,2) I,CHG(I),XNEW(IX),XNEW(IY),XNEW(IZ)                    NEW06630
  102 CONTINUE                                                          NEW06640
C                                                                       NEW06650
      DO 103 I=1,NATOM                                                  NEW06660
      IX=(I-1)*3+1                                                      NEW06670
      IY=IX+1                                                           NEW06680
      IZ=IX+2                                                           NEW06690
      XT(I)=XNEW(IX)                                                    NEW06700
      YT(I)=XNEW(IY)                                                    NEW06710
      ZT(I)=XNEW(IZ)                                                    NEW06720
  103 CONTINUE                                                          NEW06730
C                                                                       NEW06740
C   CALCULATE INTERATOMIC DISTANCES FOR NEXT GEOMETRY                   NEW06750
      CALL DIST(NATOM,XT,YT,ZT,NDIM)                                    NEW06760
      IF(NATOM.LE.2) GO TO 202                                          NEW06770
      CALL ANGLE(NATOM,XT,YT,ZT,NDIM)                                   NEW06780
C                                                                       NEW06790
  202 CONTINUE                                                          NEW06800
      IF(IFLAG.EQ.0) GO TO 203                                          NEW06810
      DO 104 I=1,NATOM                                                  NEW06820
      XX(I)=XT(I)                                                       NEW06830
      YY(I)=YT(I)                                                       NEW06840
      ZZ(I)=ZT(I)                                                       NEW06850
  104 CONTINUE                                                          NEW06860
C                                                                       NEW06870
  203 CONTINUE                                                          NEW06880
      RETURN                                                            NEW06890
      END                                                               NEW06900
      SUBROUTINE METRIC(SX,SG,FX,HH,FT,EE,PIVOT,INDEX,NVAR,METHOD,      NEW06910
     1                  UPDATE)                                         NEW06920
      IMPLICIT REAL*8 (A-H,O-Z)                                         NEW06930
      LOGICAL UPDATE                                                    NEW06940
      DIMENSION SX(N3N,NPLUS),SG(N3N,NPLUS),FX(N3N,N3N)                 NEW06950
      DIMENSION HH(NVAR,NVAR),FT(NVAR,NVAR),EE(NVAR,NVAR)               NEW06960
      DIMENSION PIVOT(NVAR),INDEX(NVAR)                                 NEW06970
      COMMON/COM101/NATOM,N3N,NINT                                      NEW06980
      COMMON/COM102/NSORT,NCOOD,NPLUS                                   NEW06990
      COMMON/COM103/CHG(50),XX(50),YY(50),ZZ(50)                        NEW07000
      COMMON/COM104/XT(150)                                             NEW07010
      COMMON/COM105/ENG0(150),ENG1(150)                                 NEW07020
      COMMON/COM106/ESCF(150),GNORM(150),SNORM(150)                     NEW07030
      COMMON/COM107/X(150),D(150),DELX(150),DELG(150),HG(150),DUMA(525) NEW07040
      COMMON/COM108/LDIM,MDIM,NDIM,IPRNT                                NEW07050
      COMMON/COM110/NSET(150),NORD(150),NXVAR(150)                      NEW07060
      COMMON/COM112/XOLD(150),XNEW(150)                                 NEW07070
      COMMON/COM113/STEPMX                                              NEW07080
      DATA ITAP15 / 15 /                                                NEW07090
      DATA A00,ONE / 0.0D+00 , 1.0D+00 /                                NEW07100
    1 FORMAT(//,2X,'********************************'/                  NEW07110
     1          2X,'***THE VARIABLE METRIC METHOD***'/                  NEW07120
     2          2X,'********************************'//)                NEW07130
    2 FORMAT(2I5)                                                       NEW07140
    3 FORMAT(3F20.10)                                                   NEW07150
    4 FORMAT(//,2X,' FX MATRIX'/)                                       NEW07160
    5 FORMAT(//,2X,' EXTRACTED FX MATRIX'/)                             NEW07170
    6 FORMAT(//,2X,' HH MATRIX'/)                                       NEW07180
    7 FORMAT(//,2X,' DETERM = ',D20.10/)                                NEW07190
    8 FORMAT(//,2X,' FX * FX(-1) = E'/)                                 NEW07200
    9 FORMAT(///,2X,' ALGOR CYCLE = ',I5)                               NEW07210
   10 FORMAT(//,2X,13H NEW X MATRIX/                                    NEW07220
     1 2X,4H NO.,5X,4H II.,6X,5H OLDX,18X,5H DELX,18X,5H NEWX/)         NEW07230
   11 FORMAT(3X,I2,4X,I5,3(3X,F20.10))                                  NEW07240
   12 FORMAT(//,2X,' GEOMETRY CHANGE IS TOO LARGE---TO BE SCALED'/      NEW07250
     1          2X,' STEPMX = ',F15.8/                                  NEW07260
     2          2X,' DD     = ',F15.8/                                  NEW07270
     3          2X,' SCALE  = ',F15.8/)                                 NEW07280
   13 FORMAT(//,2X,' NEW X MATRIX (SCALED)'/                            NEW07290
     1 2X,4H NO.,5X,4H II.,6X,5H OLDX,18X,5H DELX,18X,5H NEWX/)         NEW07300
   14 FORMAT(//,2X,20H DELG AND HG VECTORS/)                            NEW07310
   15 FORMAT(2X,I2,2(5X,F15.8))                                         NEW07320
   16 FORMAT(//,2X,7H GHG = ,F15.8/                                     NEW07330
     1 2X,7H DG  = ,F15.8)                                              NEW07340
   17 FORMAT(//,2X,' UPDATED HH MATRIX'/)                               NEW07350
   18 FORMAT(//,2X,' GEOMETRY CHANGE IS TOO LARGE---TRY AGAIN'/         NEW07360
     1          2X,' NALGOR,NMAX,NADD,XMAX = ',3I5,F10.5/)              NEW07370
   19 FORMAT(//,2X,' SUBROUTINE METRIC IS NOT ABLE TO FIND A RESONABLE  NEW07380
     1GEOMETRY GUESS'/)                                                 NEW07390
C                                                                       NEW07400
      WRITE(6,1)                                                        NEW07410
C                                                                       NEW07420
      NADD=0                                                            NEW07430
      NMAX=NPLUS                                                        NEW07440
  200 CONTINUE                                                          NEW07450
C                                                                       NEW07460
C   READ IN AN INITIAL HESSIAN MATRIX                                   NEW07470
      REWIND ITAP15                                                     NEW07480
      READ(ITAP15,2) NNATOM,NNC                                         NEW07490
      READ(ITAP15,3) ((FX(I,J),J=1,N3N),I=1,N3N)                        NEW07500
      REWIND ITAP15                                                     NEW07510
      IF(IPRNT.LE.0) GO TO 203                                          NEW07520
      WRITE(6,4)                                                        NEW07530
      CALL MATOUT(FX,N3N,N3N,N3N,N3N,6)                                 NEW07540
C                                                                       NEW07550
C   EXTRACT NECESSARY FORCE CONSTANTS                                   NEW07560
  203 CONTINUE                                                          NEW07570
      DO 102 I=1,NVAR                                                   NEW07580
      II=NXVAR(I)                                                       NEW07590
      DO 102 J=1,NVAR                                                   NEW07600
      JJ=NXVAR(J)                                                       NEW07610
      HH(I,J)=FX(II,JJ)                                                 NEW07620
      FT(I,J)=HH(I,J)                                                   NEW07630
  102 CONTINUE                                                          NEW07640
      IF(IPRNT.LE.0) GO TO 204                                          NEW07650
      WRITE(6,5)                                                        NEW07660
      CALL MATOUT(HH,NVAR,NVAR,NVAR,NVAR,6)                             NEW07670
C                                                                       NEW07680
C   OBTAIN THE H INVERSE MATRIX                                         NEW07690
  204 CONTINUE                                                          NEW07700
      CALL MATINV(HH,EE,PIVOT,INDEX,NVAR,NVAR,DETERM,NVAR)              NEW07710
      IF(IPRNT.LE.0) GO TO 205                                          NEW07720
      WRITE(6,6)                                                        NEW07730
      CALL MATOUT(HH,NVAR,NVAR,NVAR,NVAR,6)                             NEW07740
  205 WRITE(6,7) DETERM                                                 NEW07750
      IF(IPRNT.LE.0) GO TO 206                                          NEW07760
      CALL MTXMPY(FT,NVAR,HH,NVAR,EE,NVAR,EE,NVAR,NVAR,1)               NEW07770
      WRITE(6,8)                                                        NEW07780
      CALL MATOUT(EE,NVAR,NVAR,NVAR,NVAR,6)                             NEW07790
C                                                                       NEW07800
  206 CONTINUE                                                          NEW07810
      DO 103 I=1,N3N                                                    NEW07820
      X(I)=SX(I,NPLUS)                                                  NEW07830
      XT(I)=X(I)                                                        NEW07840
  103 CONTINUE                                                          NEW07850
C                                                                       NEW07860
      NALGOR=1                                                          NEW07870
      NPOS=NALGOR+NADD                                                  NEW07880
      DO 104 I=1,NVAR                                                   NEW07890
      II=NXVAR(I)                                                       NEW07900
      ENG0(I)=SG(II,NPOS)                                               NEW07910
  104 CONTINUE                                                          NEW07920
C                                                                       NEW07930
C   START ITERATION                                                     NEW07940
  300 CONTINUE                                                          NEW07950
C*    IF(IPRNT.LE.0) GO TO 301                                          NEW07960
      WRITE(6,9) NALGOR                                                 NEW07970
C                                                                       NEW07980
C   CALCULATE D VECTOR                                                  NEW07990
C   BASED ON THE NEWTON-RAPHSON TYPE STEP ( D = - H**(-1) X G )         NEW08000
  301 DO 105 I=1,NVAR                                                   NEW08010
      D(I)=A00                                                          NEW08020
      DO 105 J=1,NVAR                                                   NEW08030
      D(I)=D(I)-HH(I,J)*ENG0(J)                                         NEW08040
  105 CONTINUE                                                          NEW08050
C                                                                       NEW08060
C   FIND ALPHA                                                          NEW08070
      CALL INTER(ALPHA,D,ENG0,MDIM,NVAR)                                NEW08080
C                                                                       NEW08090
C   ESTIMATE NEW COORDINATES                                            NEW08100
C*    IF(IPRNT.LE.0) GO TO 302                                          NEW08110
      WRITE(6,10)                                                       NEW08120
  302 CONTINUE                                                          NEW08130
      XMAX=A00                                                          NEW08140
      DO 106 I=1,NVAR                                                   NEW08150
      II=NXVAR(I)                                                       NEW08160
      X00=X(II)                                                         NEW08170
      DELX(I)=ALPHA*D(I)                                                NEW08180
      DMAX=DABS(DELX(I))                                                NEW08190
      IF(DMAX.GT.XMAX) XMAX=DMAX                                        NEW08200
      X11=X00+DELX(I)                                                   NEW08210
      XT(II)=X11                                                        NEW08220
C*    IF(IPRNT.LE.0) GO TO 106                                          NEW08230
      WRITE(6,11) I,II,X00,DELX(I),X11                                  NEW08240
  106 CONTINUE                                                          NEW08250
C                                                                       NEW08260
C   CHECK THE STEP SIZE                                                 NEW08270
      IF(XMAX.GT.STEPMX) THEN                                           NEW08280
        CALL ANORM(DD,DELX,150,NVAR)                                    NEW08290
        SCALE=STEPMX/DD                                                 NEW08300
        WRITE(6,12) STEPMX,DD,SCALE                                     NEW08310
C*      IF(IPRNT.GT.0) WRITE(6,13)                                      NEW08320
        WRITE(6,13)                                                     NEW08330
        XMAX=A00                                                        NEW08340
        DO 107 I=1,NVAR                                                 NEW08350
        II=NXVAR(I)                                                     NEW08360
        X00=X(II)                                                       NEW08370
        DELOLD=DELX(I)                                                  NEW08380
        DELX(I)=DELOLD*SCALE                                            NEW08390
        DMAX=DABS(DELX(I))                                              NEW08400
        IF(DMAX.GT.XMAX) XMAX=DMAX                                      NEW08410
        X11=X00+DELX(I)                                                 NEW08420
        XT(II)=X11                                                      NEW08430
C*      IF(IPRNT.LE.0) GO TO 107                                        NEW08440
        WRITE(6,11) I,II,X00,DELX(I),X11                                NEW08450
  107   CONTINUE                                                        NEW08460
      END IF                                                            NEW08470
      IF(XMAX.GT.STEPMX) GO TO 310                                      NEW08480
C                                                                       NEW08490
      NALGOR=NALGOR+1                                                   NEW08500
      NPOS=NALGOR+NADD                                                  NEW08510
      IF(NALGOR.GT.NMAX) GO TO 320                                      NEW08520
      DO 108 I=1,NVAR                                                   NEW08530
      II=NXVAR(I)                                                       NEW08540
      ENG1(I)=SG(II,NPOS)                                               NEW08550
  108 CONTINUE                                                          NEW08560
C                                                                       NEW08570
C   CALCULATION OF HH MATRIX FOR NEXT ITERATION                         NEW08580
      DO 109 I=1,NVAR                                                   NEW08590
      DELG(I)=ENG1(I)-ENG0(I)                                           NEW08600
  109 CONTINUE                                                          NEW08610
C                                                                       NEW08620
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        NEW08630
C%  NOTICE : THE HH MATRIX CORESPONDS TO THE H INVERSE MATRIX  %        NEW08640
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        NEW08650
C                                                                       NEW08660
      IF(METHOD-2) 401,403,405                                          NEW08670
C                                                                       NEW08680
C::::::::::::::::::::::::::::                                           NEW08690
C:::MURTAGH-SARGENT METHOD:::                                           NEW08700
C::::::::::::::::::::::::::::                                           NEW08710
  401 DO 110 I=1,NVAR                                                   NEW08720
      HG(I)=DELX(I)                                                     NEW08730
      DO 110 J=1,NVAR                                                   NEW08740
      HG(I)=HG(I)-HH(I,J)*DELG(J)                                       NEW08750
  110 CONTINUE                                                          NEW08760
      IF(IPRNT.LE.0) GO TO 402                                          NEW08770
      WRITE(6,14)                                                       NEW08780
      DO 111 I=1,NVAR                                                   NEW08790
      WRITE(6,15) I,DELG(I),HG(I)                                       NEW08800
  111 CONTINUE                                                          NEW08810
  402 GHG=A00                                                           NEW08820
      DG=A00                                                            NEW08830
      DO 112 I=1,NVAR                                                   NEW08840
      GHG=GHG+HG(I)*DELG(I)                                             NEW08850
      DG=DG+HG(I)*HG(I)                                                 NEW08860
  112 CONTINUE                                                          NEW08870
      DO 113 I=1,NVAR                                                   NEW08880
      DO 113 J=1,NVAR                                                   NEW08890
      HH(I,J)=HH(I,J)+(HG(I)*HG(J))/GHG                                 NEW08900
  113 CONTINUE                                                          NEW08910
      IF(IPRNT.LE.0) GO TO 410                                          NEW08920
      WRITE(6,16) GHG,DG                                                NEW08930
      GO TO 410                                                         NEW08940
C                                                                       NEW08950
C:::::::::::::::::::::                                                  NEW08960
C:::FLETCHER METHOD:::                                                  NEW08970
C:::::::::::::::::::::                                                  NEW08980
  403 DO 120 I=1,NVAR                                                   NEW08990
      HG(I)=A00                                                         NEW09000
      DO 120 J=1,NVAR                                                   NEW09010
      HG(I)=HG(I)+HH(I,J)*DELG(J)                                       NEW09020
  120 CONTINUE                                                          NEW09030
      IF(IPRNT.LE.0) GO TO 404                                          NEW09040
      WRITE(6,14)                                                       NEW09050
      DO 121 I=1,NVAR                                                   NEW09060
      WRITE(6,15) I,DELG(I),HG(I)                                       NEW09070
  121 CONTINUE                                                          NEW09080
  404 GHG=A00                                                           NEW09090
      DG=A00                                                            NEW09100
      DO 122 I=1,NVAR                                                   NEW09110
      GHG=GHG+HG(I)*DELG(I)                                             NEW09120
      DG=DG+DELX(I)*DELG(I)                                             NEW09130
  122 CONTINUE                                                          NEW09140
      DO 123 I=1,NVAR                                                   NEW09150
      DO 123 J=1,NVAR                                                   NEW09160
      HH(I,J)=HH(I,J)-(HG(I)*DELX(J)+HG(J)*DELX(I))/DG                  NEW09170
     1 +(ONE+GHG/DG)*((DELX(I)*DELX(J))/DG)                             NEW09180
  123 CONTINUE                                                          NEW09190
      IF(IPRNT.LE.0) GO TO 410                                          NEW09200
      WRITE(6,16) GHG,DG                                                NEW09210
      GO TO 410                                                         NEW09220
C                                                                       NEW09230
C::::::::::::::::::::::::::::::::::::                                   NEW09240
C:::DAVIDON-FLETCHER-POWELL METHOD:::                                   NEW09250
C::::::::::::::::::::::::::::::::::::                                   NEW09260
  405 DO 130 I=1,NVAR                                                   NEW09270
      HG(I)=A00                                                         NEW09280
      DO 130 J=1,NVAR                                                   NEW09290
      HG(I)=HG(I)+HH(I,J)*DELG(J)                                       NEW09300
  130 CONTINUE                                                          NEW09310
      IF(IPRNT.LE.0) GO TO 406                                          NEW09320
      WRITE(6,14)                                                       NEW09330
      DO 131 I=1,NVAR                                                   NEW09340
      WRITE(6,15) I,DELG(I),HG(I)                                       NEW09350
  131 CONTINUE                                                          NEW09360
  406 GHG=A00                                                           NEW09370
      DG=A00                                                            NEW09380
      DO 132 I=1,NVAR                                                   NEW09390
      GHG=GHG+HG(I)*DELG(I)                                             NEW09400
      DG=DG+DELX(I)*DELG(I)                                             NEW09410
  132 CONTINUE                                                          NEW09420
      DO 133 I=1,NVAR                                                   NEW09430
      DO 133 J=1,NVAR                                                   NEW09440
      HH(I,J)=HH(I,J)-(HG(I)*HG(J))/GHG+(DELX(I)*DELX(J))/DG            NEW09450
  133 CONTINUE                                                          NEW09460
      IF(IPRNT.LE.0) GO TO 410                                          NEW09470
      WRITE(6,16) GHG,DG                                                NEW09480
C                                                                       NEW09490
C#####################                                                  NEW09500
C###UPDATE H MATRIX###                                                  NEW09510
C#####################                                                  NEW09520
  410 CONTINUE                                                          NEW09530
      IF(IPRNT.LE.0) GO TO 415                                          NEW09540
      WRITE(6,17)                                                       NEW09550
      CALL MATOUT(HH,NVAR,NVAR,NVAR,NVAR,6)                             NEW09560
C                                                                       NEW09570
  415 DO 135 I=1,NVAR                                                   NEW09580
      ENG0(I)=ENG1(I)                                                   NEW09590
  135 CONTINUE                                                          NEW09600
      GO TO 300                                                         NEW09610
C                                                                       NEW09620
C   GEOMETRY CHANGE IS TOO LARGE                                        NEW09630
  310 CONTINUE                                                          NEW09640
      WRITE(6,18) NALGOR,NMAX,NADD,XMAX                                 NEW09650
      NMAX=NMAX-1                                                       NEW09660
      NADD=NADD+1                                                       NEW09670
      IF(NMAX.LE.0) GO TO 315                                           NEW09680
      GO TO 200                                                         NEW09690
C                                                                       NEW09700
C   YOU CANNOT FIND NEXT GEOMETRY                                       NEW09710
  315 WRITE(6,19)                                                       NEW09720
      RETURN                                                            NEW09730
C                                                                       NEW09740
  320 CONTINUE                                                          NEW09750
      DO 136 I=1,N3N                                                    NEW09760
      XOLD(I)=SX(I,NPLUS)                                               NEW09770
      XNEW(I)=XT(I)                                                     NEW09780
  136 CONTINUE                                                          NEW09790
      IFLAG=0                                                           NEW09800
      IF(UPDATE) IFLAG=1                                                NEW09810
      CALL OLDNEW(XOLD,XNEW,IFLAG)                                      NEW09820
C                                                                       NEW09830
      RETURN                                                            NEW09840
      END                                                               NEW09850
      SUBROUTINE INTER(ALPHA,D,ENG0,NND,NN)                             NEW09860
      IMPLICIT REAL*8 (A-H,O-Z)                                         NEW09870
      DIMENSION D(NND),ENG0(NND)                                        NEW09880
      DATA A00,ONE,TWO / 0.0D+00 , 1.0D+00 , 2.0D+00 /                  NEW09890
      DATA FDEL / 1.0D+00 /                                             NEW09900
C                                                                       NEW09910
      GXS=A00                                                           NEW09920
      DO 101 I=1,NN                                                     NEW09930
      GXS=GXS+D(I)*ENG0(I)                                              NEW09940
  101 CONTINUE                                                          NEW09950
      FX0=-(TWO*FDEL)/GXS                                               NEW09960
      ALPHA=DMIN1(ONE,FX0)                                              NEW09970
      RETURN                                                            NEW09980
      END                                                               NEW09990
      SUBROUTINE MTXMPY(A,NAD,B,NBD,C,NCD,D,NDD,N,IMPY)                 NEW10000
      IMPLICIT REAL*8 (A-H,O-Z)                                         NEW10010
      DIMENSION A(NAD,NAD),B(NBD,NBD),C(NCD,NCD),D(NDD,NDD)             NEW10020
      DATA A00 / 0.0D+00 /                                              NEW10030
C                                                                       NEW10040
      DO 101 I=1,NCD                                                    NEW10050
      DO 101 J=1,NCD                                                    NEW10060
  101 C(I,J)=A00                                                        NEW10070
      GO TO (201,202,203,204,205,206),IMPY                              NEW10080
C     C=A*B                                                             NEW10090
  201 DO 102 I=1,N                                                      NEW10100
      DO 102 J=1,N                                                      NEW10110
      DO 102 K=1,N                                                      NEW10120
  102 C(I,J)=C(I,J)+A(I,K)*B(K,J)                                       NEW10130
      RETURN                                                            NEW10140
C     C=TA*B                                                            NEW10150
  202 DO 103 I=1,N                                                      NEW10160
      DO 103 J=1,N                                                      NEW10170
      DO 103 K=1,N                                                      NEW10180
  103 C(I,J)=C(I,J)+A(K,I)*B(K,J)                                       NEW10190
      RETURN                                                            NEW10200
C     C=A*TB                                                            NEW10210
  203 DO 104 I=1,N                                                      NEW10220
      DO 104 J=1,N                                                      NEW10230
      DO 104 K=1,N                                                      NEW10240
  104 C(I,J)=C(I,J)+A(I,K)*B(J,K)                                       NEW10250
      RETURN                                                            NEW10260
C     C=TA*TB                                                           NEW10270
  204 DO 105 I=1,N                                                      NEW10280
      DO 105 J=1,N                                                      NEW10290
      DO 105 K=1,N                                                      NEW10300
  105 C(I,J)=C(I,J)+A(K,I)*B(J,K)                                       NEW10310
      RETURN                                                            NEW10320
C     C=TA*B*A                                                          NEW10330
  205 DO 106 I=1,N                                                      NEW10340
      DO 106 J=1,N                                                      NEW10350
      D(I,J)=A00                                                        NEW10360
      DO 106 K=1,N                                                      NEW10370
  106 D(I,J)=D(I,J)+A(K,I)*B(K,J)                                       NEW10380
      DO 107 I=1,N                                                      NEW10390
      DO 107 J=1,N                                                      NEW10400
      DO 107 K=1,N                                                      NEW10410
  107 C(I,J)=C(I,J)+D(I,K)*A(K,J)                                       NEW10420
      RETURN                                                            NEW10430
C     C=A*B*TA                                                          NEW10440
  206 DO 108 I=1,N                                                      NEW10450
      DO 108 J=1,N                                                      NEW10460
      D(I,J)=A00                                                        NEW10470
      DO 108 K=1,N                                                      NEW10480
  108 D(I,J)=D(I,J)+A(I,K)*B(K,J)                                       NEW10490
      DO 109 I=1,N                                                      NEW10500
      DO 109 J=1,N                                                      NEW10510
      DO 109 K=1,N                                                      NEW10520
  109 C(I,J)=C(I,J)+D(I,K)*A(J,K)                                       NEW10530
      RETURN                                                            NEW10540
      END                                                               NEW10550
