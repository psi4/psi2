      subroutine fentry(cc,ia,maxcor)
C**********************************************************************
c*      PROGRAM MODOPTX                                               *
C**********************************************************************
C*  Last updated on December 10, 1996 by Dr. Yukio Yamaguchi          *
c   Units are added to geometry and gradients                         *
C**********************************************************************
C*  Last updated on September 11, 1993 by Dr. Yukio Yamaguchi         *
C   Further modification for the UNIX System                          *
C**********************************************************************
C*  Last updated on October 04, 1991 by Dr. Yukio Yamaguchi           *
C   Modification for the UNIX System                                  *
C**********************************************************************
C*  LAST UPDATED ON FEBRUARY 05, 1991 BY YUKIO YAMAGUCHI              *
C********************************************************************** 
C   THE SOURCE CODE OF GEOMETRY OPTIMIZATION PROGRAM                  *
C   BASED ON THE MODE FOLLWOING METHOD                                *
C   IN THE CARTESIAN COORDINATE SYSTEM                                *
C   THIS PROGRAM IS ABLE TO TREAT UP TO '50' ATOMS                    *
C   WITH UP TO '150' INTERNAL COORDINATES                             *
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)                              
      dimension cc(maxcor),ia(maxcor*2)
      DIMENSION I30(200),NUNQ(50),NAT(50,50),COORD(3,50)                MOD00120
      CHARACTER*8 UPD                                                   MOD00130
      LOGICAL UPDATE                                                    MOD00140
      COMMON/COM101/NATOM,N3N,NINT,NCOOD,NPLUS                          MOD00150
      COMMON/COM102/ISTAT,NVAR,IMODE,NVARS                              MOD00160
      COMMON/COM103/CHG(50),XX(50),YY(50),ZZ(50)                        MOD00170
      COMMON/COM104/EX(50),EY(50),EZ(50)                                MOD00180
      COMMON/COM105/ENGR(150),RCOOD(150)                                MOD00190
      COMMON/COM106/ESCF(150),GNORM(150),SNORM(150)                     MOD00200
      COMMON/COM107/R(1275)                                             MOD00210
      COMMON/COM108/IOFF(150),IPRNT                                     MOD00220
      COMMON/COM109/XMAX,YMAX,ZMAX                                      MOD00230
      COMMON/COM110/NSET(150),NORD(150),NXVAR(150)                      MOD00240
      COMMON/COM111/IHUPD                                               MOD00250
      COMMON/COM112/XOLD(150),XNEW(150)                                 MOD00260
      COMMON/COM113/ENG0(150),ENG1(150)                                 MOD00270
      COMMON/COM114/EIGMAX,EIGMIN,STEPMX                                MOD00280
      DATA A00,TEN / 0.0D+00 , 10.0D+00 /                               MOD00300
      DATA FMILI / 1.0D-3 /                                             MOD00310
    1 FORMAT(//,2X,' THE GEOMETRY OPTIMIZATION PROGRAM'/                MOD00320
     1          2X,' BASED ON THE MODE FOLLOWING METHOD'/               MOD00330
     2          2X,' CARTESIAN COORDINATE SYSTEM'/                      MOD00340
     3          2X,' (VERSION OF 12/10/96)'//)                          MOD00350
    2 FORMAT(2I5,F20.10,I5)                                             MOD00360
    3 FORMAT(14I5)                                                      MOD00370
    4 FORMAT(14I5)                                                      MOD00380
    5 FORMAT(A5,3I5)                                                    MOD00390
    6 FORMAT(2X,' NCOOD   = ',I5/                                       MOD00400
     1 2X,' NSORT   = ',I5/                                             MOD00410
     2 2X,' IHUPD   = ',I5/                                             MOD00420
     3 2X,' NVAR    = ',I5/                                             MOD00430
     4 2X,' ISTEP   = ',I5/                                             MOD00440
     5 2X,' ISTAT   = ',I5/                                             MOD00450
     6 2X,' IMODE   = ',I5/                                             MOD00460
     7 2X,' IHMAX   = ',I5/                                             MOD00470
     8 2X,' IHMIN   = ',I5/                                             MOD00480
     9 2X,' IPRNT   = ',I5)                                             MOD00490
    7 FORMAT(2X,' NATOM   = ',I5/                                       MOD00500
     1 2X,' N3N     = ',I5/                                             MOD00510
     2 2X,' NINT    = ',I5/                                             MOD00520
     3 2X,' NVARS   = ',I5/                                             MOD00530
     4 2X,' IFORCE  = ',I5/                                             MOD00540
     5 2X,' NSTORE  = ',I5/                                             MOD00550
     6 2X,' NPLUS   = ',I5/                                             MOD00560
     7 2X,' NCHNG   = ',I5/                                             MOD00570
     8 2X,' NCONV   = ',I5/                                             MOD00580
     9 2X,' NUNIQ   = ',I5)                                             MOD00590
    8 FORMAT(//,2X,' OLD COORDINATES IN TAPE30'/                        MOD00600
     1 39X,' OLD X IN BOHR'/
     2 4X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z/)               MOD00620
    9 FORMAT(3X,I5,1X,4F16.10)                                          MOD00630
   10 FORMAT(//,2X,' NEW COORDINATES IN TAPE30'/                        MOD00640
     1 39X,' NEW X IN BOHR'/
     2 4X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z/)               MOD00660
   11 FORMAT(//,2X,' GEOMETRY IN TAPE30 IS UPDATED'/)                   MOD00670
   12 FORMAT(//,2X,':::::::::::::::::::::::::::'/                       MOD00680
     1          2X,':::GEOMETRY IS OPTIMIZED:::'/                       MOD00690
     2          2X,':::::::::::::::::::::::::::'//                      MOD00700
     3          2X,' GNORM = ',F15.10,7X,' GCONV = ',F15.10/)           MOD00710
   13 FORMAT(//,2X,' TOO BIG CHANGE IN GEOMETRY---SCRATCH YOUR HEAD TO FMOD00720
     1IND A NEW GEOMETRY'//                                             MOD00730
     2          2X,' XMAX = ',F15.10,7X,' YMAX = ',F15.10,7X,' ZMAX = ',MOD00740
     3                        F15.10,7X,' GCHNG = ',F15.10/)            MOD00750
   14 FORMAT(//,2X,' GEOMETRY IN TAPE30 IS NOT UPDATED'/)               MOD00760
C                                                                       MOD00800
      INPUT=5                                                           MOD00810
      ITAPE6=6                                                          MOD00820
      ITAP12=12                                                         MOD00830
      ITAP13=13                                                         MOD00840
      ITAP15=15                                                         MOD00850
      ITAP30=30                                                         MOD00860
      mxatom=50
      IFLAG=0                                                           MOD00870
C                                                                       MOD00890
      call psinit('APPEND')
      CALL TSTART(6)                                                    MOD00780
      WRITE(6,1)                                                        MOD00900
      write(6,*) '  maxcor = ',maxcor
      write(6,*) '  '
      call ffile(itap12,'file12',0)
      call ffile(itap13,'file13',0)
      call ffile(itap15,'file15',0)
      REWIND ITAP12                                                     MOD00910
      REWIND ITAP13                                                     MOD00920
      READ(ITAP12,2) NATOM,NINT,ESCFZ,IFORCE                            MOD00930
C                                                                       MOD00940
C   NCOOD IS NUMBER OF INTERNAL COORDINATES TO BE OPTIMIZED             MOD00950
C   NSORT IS A PARAMETER FOR SORTING                                    MOD00960
C     NSORT=0---NO SORT                                                 MOD00970
C     NSORT=1---GRADIENT SORT                                           MOD00980
C     NSORT=2---ENERGY SORT                                             MOD00990
C     NSORT=3---SELECTTED GRADIENT SORT                                 MOD01000
C   IHUPD IS A PARAMETER FOR UPDATING THE HESSIAN MATRIX                MOD01010
C     IHUPD=0---DEFAULT                                                 MOD01020
C     IHUPD=1---THE BFGS UPDATE   (FOR A MINIMUM)                       MOD01030
C     IHUPD=2---THE POWELL UDPATE (FOR A TRANSITION STATE)              MOD01040
C     IHUPD=3---NO UPDATE                                               MOD01050
C   NVAR IS NUMBER OF CARTESIAN COORDINATES TO BE OPTIMIZED             MOD01060
C   ISTEP IS A PARAMETER FOR MAXIMUM STEP SIZE ALLOWED                  MOD01070
C     ISTEP= 0---0.10                                                   MOD01080
C     ISTEP= N---0.01*ISTEP                                             MOD01090
C     ISTEP=-N---0.01/ISTEP                                             MOD01100
C   ISTAT IS A PARAMETER FOR STATIONARY POINT                           MOD01110
C     ISTAT=0---MINIMUM                                                 MOD01120
C     ISTAT=1---TRANSITION STATE                                        MOD01130
C   IMODE IS A PARAMETER FOR MODE SEARCHING                             MOD01140
C     IMODE=0---FOLLOW THE LOWEST MODE                                  MOD01150
C     IMODE=N---FOLLOW THE N-TH EIGENVECTOR                             MOD01160
C   IHMAX IS A PARAMETER FOR MAXIMUM MAGNITUDE OF HESSIAN EIGENVALUES   MOD01170
C     IHMAX=0---25.0                                                    MOD01180
C     IHMAX=N---DFLOAT(N)                                               MOD01190
C   IHMIN IS A PARAMETER FOR MINIMUM MAGNITUDE OF HESSIAN EIGENVALUES   MOD01200
C     IHMIN=0---0.001
C     IHMIN=N---10.0D(-N)                                               MOD01220
C   IPRNT IS A PARAMETER FOR OUTPUT                                     MOD01230
C                                                                       MOD01240
      N3N=NATOM*3                                                       MOD01250
      CALL LOCATE(INPUT,'# OPTIM ##',IERR)                              MOD01260
      READ(5,3) NCOOD,NSORT,IHUPD,NVAR,ISTEP,ISTAT,IMODE,IHMAX,IHMIN,   MOD01270
     1          IPRNT                                                   MOD01280
      IF(IHMAX.EQ.0) IHMAX=25                                           MOD01290
      IF(IHMIN.EQ.0) IHMIN=3                                            MOD01300
      EIGMAX=DFLOAT(IHMAX)                                              MOD01310
      EIGMIN=TEN**(-IHMIN)                                              MOD01320
      IF(ISTEP.EQ.0) ISTEP=10                                           MOD01330
      IF(ISTEP.GT.0) THEN                                               MOD01340
        STEPMX=1.0D-02*DFLOAT(ISTEP)                                    MOD01350
      ELSE                                                              MOD01360
        STEPMX=1.0D-02/DFLOAT(IABS(ISTEP))                              MOD01370
      END IF                                                            MOD01380
      IF(ISTAT.EQ.0.AND.IHUPD.EQ.0) IHUPD=1                             MOD01390
      IF(ISTAT.NE.0.AND.IHUPD.EQ.0) IHUPD=2                             MOD01400
C                                                                       MOD01410
      IOFF(1)=0                                                         MOD01420
      DO 101 I=1,149                                                    MOD01430
      IOFF(I+1)=IOFF(I)+I                                               MOD01440
  101 CONTINUE                                                          MOD01450
C                                                                       MOD01460
      READ(5,4) (NSET(I),I=1,NCOOD)                                     MOD01470
      DO 102 I=1,N3N                                                    MOD01480
  102 NXVAR(I)=I                                                        MOD01490
      IF(NVAR.EQ.0) NVAR=N3N                                            MOD01500
      IF(NVAR.EQ.N3N) GO TO 201                                         MOD01510
      READ(5,4) (NXVAR(I),I=1,NVAR)                                     MOD01520
  201 READ(ITAP13,3) NSTORE                                             MOD01530
      NPLUS=NSTORE+1                                                    MOD01540
C                                                                       MOD01550
      NVARS=NVAR                                                        MOD01560
      IF(ISTAT.EQ.0) NVARS=NVAR+1                                       MOD01570
C                                                                       MOD01580
C  READ IN UPDATING INFORMATION                                         MOD01590
      READ(5,5) UPD,NCHNG,NCONV,NUNIQ                                   MOD01600
      IF(UPD.EQ.'UP     ') UPDATE=.TRUE.                                MOD01610
      IF(UPD.EQ.'UPDATE ') UPDATE=.TRUE. 
      IF(UPD.EQ.'up     ') UPDATE=.TRUE.                                MOD01610
      IF(UPD.EQ.'update ') UPDATE=.TRUE. 
      IF(NCHNG.EQ.0) NCHNG=2                                            MOD01620
      IF(NCONV.EQ.0) NCONV=7                                            MOD01630
      IF(NUNIQ.EQ.0) NUNIQ=NATOM                                        MOD01640
      WRITE(6,6) NCOOD,NSORT,IHUPD,NVAR,ISTEP,ISTAT,IMODE,IHMAX,IHMIN,  MOD01650
     *           IPRNT                                                  MOD01660
      WRITE(6,7) NATOM,N3N,NINT,NVARS,IFORCE,NSTORE,NPLUS,NCHNG,        MOD01670
     *           NCONV,NUNIQ                                            MOD01680
      GCHNG=TEN**(-NCHNG)                                               MOD01690
      GCONV=TEN**(-NCONV)                                               MOD01700
C                                                                       MOD01710
      IF(.NOT.UPDATE) GO TO 202                                         MOD01720
      DO 103 I=1,NATOM                                                  MOD01730
      NUNQ(I)=1                                                         MOD01740
      NAT(I,1)=I                                                        MOD01750
  103 CONTINUE                                                          MOD01760
      IF(NUNIQ.EQ.NATOM) GO TO 202                                      MOD01770
      DO 104 I=1,NUNIQ                                                  MOD01780
      READ(5,4) NN,(NAT(I,J),J=1,NN)                                    MOD01790
      NUNQ(I)=NN                                                        MOD01800
  104 CONTINUE                                                          MOD01810
C                                                                       MOD01820
C***********************************************************************MOD01830
C   DYNAMIC ALLOCATION                                                 *MOD01840
C***********************************************************************MOD01850
C   1 : SX(N3N,NPLUS)                                                  *MOD01860
C   2 : SG(N3N,NPLUS)                                                  *MOD01870
C   3 : SCOOD(NINT,NPLUS)                                              *MOD01880
C   4 : SENG(NINT,NPLUS)                                               *MOD01890
C***********************************************************************MOD01900
  202 CONTINUE                                                          MOD01910
      IC1=1                                                             MOD01920
      IC2=IC1+N3N*NPLUS                                                 MOD01930
      IC3=IC2+N3N*NPLUS                                                 MOD01940
      IC4=IC3+NINT*NPLUS                                                MOD01950
      ICSAVE=IC4+NINT*NPLUS                                             MOD01960
C***********************************************************************MOD01970
C                                                                       MOD01980
C   READ IN STORED INFORMATION FROM TAPE13                              MOD01990
      IC5=ICSAVE                                                        MOD02000
      IC6=IC5+N3N*NPLUS                                                 MOD02010
      IC7=IC6+N3N*NPLUS                                                 MOD02020
      IC8=IC7+NINT*NPLUS                                                MOD02030
      ICMAX=IC8+NINT*NPLUS                                              MOD02040
C.................X       G       COOD    ENG.....                      MOD02050
      CALL READIN(CC(IC5),CC(IC6),CC(IC7),CC(IC8))                      MOD02060
C                                                                       MOD02070
C   SORT DATA W.R.T. GRADIENT OR ENERGY                                 MOD02080
C...............SX      SG      SCOOD   SENG    X       G.......        MOD02090
      CALL SORT(CC(IC1),CC(IC2),CC(IC3),CC(IC4),CC(IC5),CC(IC6),        MOD02100
C...............COOD    ENG...........                                  MOD02110
     1          CC(IC7),CC(IC8),NSORT)                                  MOD02120
C                                                                       MOD02130
      GNORML=GNORM(NPLUS)*FMILI                                         MOD02140
      IF(NSORT.EQ.3) GNORML=SNORM(NPLUS)*FMILI                          MOD02150
      IF(GNORML.LT.GCONV) GO TO 215                                     MOD02160
C                                                                       MOD02170
C********************************************                           MOD02180
C***THE MODE FOLLOWING OPTIMIZATION METHOD***                           MOD02190
C********************************************                           MOD02200
C                                                                       MOD02210
      IC5=ICSAVE                                                        MOD02220
      IC6=IC5+N3N*N3N                                                   MOD02230
      IC7=IC6+NVAR*NVAR                                                 MOD02240
      IC8=IC7+NVARS*NVARS
      IC9=IC8+NVAR*NVAR                                                 MOD02260
      IC10=IC9+NVAR                                                     MOD02270
      IC11=IC10+NVARS*NVARS                                             MOD02280
      IC12=IC11+NVARS*NVARS                                             MOD02290
      IC13=IC12+NVARS                                                   MOD02300
      IC14=IC13+NVARS                                                   MOD02310
      IC15=IC14+NVARS*NVARS                                             MOD02320
      IC16=IC15+NVARS                                                   MOD02330
      ICMAX=IC16+NVARS                                                  MOD02340
C................SX      SG      FX      HH      HT      ELH.....       MOD02350
      CALL OPTIM(CC(IC1),CC(IC2),CC(IC5),CC(IC6),CC(IC7),CC(IC8),       MOD02360
C................EVALU   EE       GG       FV1      FV2......           MOD02370
     1           CC(IC9),CC(IC10),CC(IC11),CC(IC12),CC(IC13),           MOD02380
C................ELZ      EZVZL    SAVEK...........                     MOD02390
     2           CC(IC14),CC(IC15),CC(IC16),UPDATE)                     MOD02400
C                                                                       MOD02410
C*******************************                                        MOD02420
C***UPDATE TAPE30 IF REQUIRED***                                        MOD02430
C*******************************                                        MOD02440
C                                                                       MOD02450
      IF(.NOT.UPDATE) GO TO 218                                         MOD02460
C                                                                       MOD02470
      IF(DABS(XMAX).GT.GCHNG) GO TO 217                                 MOD02480
      IF(DABS(YMAX).GT.GCHNG) GO TO 217                                 MOD02490
      IF(DABS(ZMAX).GT.GCHNG) GO TO 217                                 MOD02500
C                                                                       MOD02510
C   READ IN GEOMETRY FROM TAPE30                                        MOD02520
      CALL RFILE(ITAP30)                                                MOD02530
      CALL WREADW(ITAP30,I30,200,101,JUNK)                              MOD02540
      MPOINT=I30(2)                                                     MOD02550
      MCONST=I30(3)                                                     MOD02560
      NCALCS=I30(5)                                                     MOD02570
      IGEOP=100+MCONST+MPOINT+NCALCS                                    MOD02580
      CALL WREADW(ITAP30,LOCCAL,1,IGEOP,JUNK)                           MOD02590
      IGEOP=LOCCAL+80                                                   MOD02600
C   UPDATE GEOMETRY ON TAPE30                                           MOD02610
      CALL WREADW(ITAP30,COORD,N3N*2,IGEOP,JUNK)                        MOD02620
      WRITE(6,8)                                                        MOD02630
      DO 105 I=1,NATOM                                                  MOD02640
      WRITE(6,9) I,CHG(I),COORD(1,I),COORD(2,I),COORD(3,I)              MOD02650
  105 CONTINUE                                                          MOD02660
C                                                                       MOD02670
C   AVERAGE OUT CARTESIAN COORDINATES FOR EQUIVALENT ATOMS              MOD02680
      DO 108 I=1,NUNIQ                                                  MOD02690
      NN=NUNQ(I)                                                        MOD02700
      SUMX=A00                                                          MOD02710
      SUMY=A00                                                          MOD02720
      SUMZ=A00                                                          MOD02730
      DO 106 J=1,NN                                                     MOD02740
      MM=NAT(I,J)                                                       MOD02750
      SUMX=SUMX+DABS(XX(MM))                                            MOD02760
      SUMY=SUMY+DABS(YY(MM))                                            MOD02770
      SUMZ=SUMZ+DABS(ZZ(MM))                                            MOD02780
  106 CONTINUE                                                          MOD02790
      ANX=DBLE(NN)                                                      MOD02800
      AVGX=SUMX/ANX                                                     MOD02810
      AVGY=SUMY/ANX                                                     MOD02820
      AVGZ=SUMZ/ANX                                                     MOD02830
      DO 107 J=1,NN                                                     MOD02840
      MM=NAT(I,J)                                                       MOD02850
      XX(MM)=DSIGN(AVGX,XX(MM))                                         MOD02860
      YY(MM)=DSIGN(AVGY,YY(MM))                                         MOD02870
      ZZ(MM)=DSIGN(AVGZ,ZZ(MM))                                         MOD02880
  107 CONTINUE                                                          MOD02890
  108 CONTINUE                                                          MOD02900
      WRITE(6,10)                                                       MOD02910
      DO 109 I=1,NATOM                                                  MOD02920
      COORD(1,I)=XX(I)                                                  MOD02930
      COORD(2,I)=YY(I)                                                  MOD02940
      COORD(3,I)=ZZ(I)                                                  MOD02950
      WRITE(6,9) I,CHG(I),COORD(1,I),COORD(2,I),COORD(3,I)              MOD02960
  109 CONTINUE                                                          MOD02970
      CALL WWRITW(ITAP30,COORD,N3N*2,IGEOP,JUNK)                        MOD02980
      WRITE(6,11)                                                       MOD02990
      CALL RCLOSE(ITAP30,3)                                             MOD03000
      GO TO 220                                                         MOD03010
C                                                                       MOD03020
C   GEOMETRY IS OPTIMIZED                                               MOD03030
  215 CONTINUE                                                          MOD03040
      WRITE(6,12) GNORML,GCONV                                          MOD03050
      GO TO 220                                                         MOD03060
C                                                                       MOD03070
C   GEOMETRY CHANGE IS TOO LARGE                                        MOD03080
  217 CONTINUE                                                          MOD03090
      WRITE(6,13) XMAX,YMAX,ZMAX,GCHNG                                  MOD03100
      GO TO 220                                                         MOD03110
C                                                                       MOD03120
C   GEOMETRY IN TAPE30 IS NOT UPDATED                                   MOD03130
  218 CONTINUE                                                          MOD03140
      WRITE(6,14)                                                       MOD03150
C                                                                       MOD03160
  220 CONTINUE                                                          MOD03170
C                                                                       MOD03180
      CALL TSTOP(6)                                                     MOD03190
      STOP                                                              MOD03200
      END                                                               MOD03210
      SUBROUTINE READIN(X,G,COOD,ENG)                                   MOD03220
      IMPLICIT REAL*8 (A-H,O-Z)                                         MOD03230
      DIMENSION X(N3N,NPLUS),G(N3N,NPLUS)                               MOD03240
      DIMENSION COOD(NINT,NPLUS),ENG(NINT,NPLUS)                        MOD03250
      COMMON/COM101/NATOM,N3N,NINT,NCOOD,NPLUS                          MOD03260
      COMMON/COM103/CHG(50),XX(50),YY(50),ZZ(50)                        MOD03270
      COMMON/COM104/EX(50),EY(50),EZ(50)                                MOD03280
      COMMON/COM105/ENGR(150),RCOOD(150)                                MOD03290
      COMMON/COM106/ESCF(150),GNORM(150),SNORM(150)                     MOD03300
      COMMON/COM108/IOFF(150),IPRNT                                     MOD03310
      COMMON/COM110/NSET(150),NORD(150),NXVAR(150)                      MOD03320
      DATA ITAP12,ITAP13 / 12 , 13 /                                    MOD03330
      DATA FKILO,FMILI / 1.0D+3 , 1.0D-3 /                              MOD03340
    1 FORMAT(2I5,F20.10,I5)                                             MOD03350
    2 FORMAT(10I5)                                                      MOD03360
    3 FORMAT(I5,F20.10)                                                 MOD03370
    4 FORMAT(////,2X,10H CASE   = ,I5/                                  MOD03380
     1 2X,10H ESCF   = ,F20.10)                                         MOD03390
    5 FORMAT(//,2X,' CARTESIAN COORDINATES (IN BOHR)',42X,' GRADIENTS (I 
     1N HARTREE/BOHR)'/
     2 4X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z,                MOD03410
     3                       17X,2H X,14X,2H Y,14X,2H Z)                MOD03420
    6 FORMAT(4F20.10)                                                   MOD03430
    7 FORMAT(3X,I5,1X,4F16.10,3X,3F16.10)                               MOD03440
    8 FORMAT(//,2X,' INTERNAL GRADIENTS (IN 1.0D+3 BOHR)'/              MOD03450
     *2x,' (BOND LENGTH IN HARTREE/BOHR AND BOND ANGLE IN HARTREE/RADIAN
     *)'/5X,4H NO.,10X,7H COORD.,13X,6H GRAD./)                           MOD03460
    9 FORMAT(2F20.10)                                                   MOD03470
   10 FORMAT(2X,I7,5X,2F20.10)                                          MOD03480
   11 FORMAT(/,2X,22H NORM FOR GRADIENTS = ,F15.7)                      MOD03490
   12 FORMAT(//,2X,' ***SELECTED VARIABLES DO NOT MATCH***'/            MOD03500
     1 2X,' NCOOD = ',I5,' II = ',I5/)                                  MOD03510
   13 FORMAT(/,2X,' NORM FOR CASE ',I5,3X,' IS ',F15.7/)                MOD03520
   14 FORMAT(///,2X,' CARTESIAN COORDINATE MATRIX'/)                    MOD03530
   15 FORMAT(//,2X,' INTERNAL COORDINATE MATRIX'/)                      MOD03540
   16 FORMAT(//,2X,' CARTESIAN GRADIENT MATRIX'/)                       MOD03550
   17 FORMAT(//,2X,' INTERNAL GRADINET MATRIX'/)                        MOD03560
   18 FORMAT(//,2X,' GRADIENTS NORMS AND ENERGIES'/                     MOD03570
     1 3X,4H NO.,8X,6H GNORM,14X,6H SNORM,14X,7H ENERGY/)               MOD03580
   19 FORMAT(2X,I4,5X,3F20.10)                                          MOD03590
C                                                                       MOD03600
      REWIND ITAP12                                                     MOD03610
      READ(ITAP12,1) NATOM,NINT,ESCFZ,IFORCE                            MOD03620
      REWIND ITAP13                                                     MOD03630
      READ(ITAP13,2) NSTORE                                             MOD03640
C                                                                       MOD03650
      DO 105 ICASE=1,NPLUS                                              MOD03660
      IF(ICASE.EQ.NPLUS) GO TO 202                                      MOD03670
      READ(ITAP13,3) JCASE,ESCF(ICASE)                                  MOD03680
      GO TO 203                                                         MOD03690
  202 ESCF(ICASE)=ESCFZ                                                 MOD03700
  203 IF(IPRNT.LE.0) GO TO 301                                          MOD03710
      WRITE(6,4) ICASE,ESCF(ICASE)                                      MOD03720
C                                                                       MOD03730
  301 CONTINUE                                                          MOD03740
      IF(IPRNT.LE.0) GO TO 302                                          MOD03750
      WRITE(6,5)                                                        MOD03760
C*********************************************                          MOD03770
C***READ CARTESIAN COORDINATES OF EACH CASE***                          MOD03780
C*********************************************                          MOD03790
  302 DO 102 I=1,NATOM                                                  MOD03800
      IF(ICASE.EQ.NPLUS) GO TO 204                                      MOD03810
      READ(ITAP13,6) CHG(I),XX(I),YY(I),ZZ(I)                           MOD03820
      GO TO 102                                                         MOD03830
  204 READ(ITAP12,6) CHG(I),XX(I),YY(I),ZZ(I)                           MOD03840
  102 CONTINUE                                                          MOD03850
C*******************************************                            MOD03860
C***READ CARTESIAN GRADIENTS OF EACH CASE***                            MOD03870
C*******************************************                            MOD03880
      DO 103 I=1,NATOM                                                  MOD03890
      IF(ICASE.EQ.NPLUS) GO TO 205                                      MOD03900
      READ(ITAP13,6) CHG(I),EX(I),EY(I),EZ(I)                           MOD03910
      GO TO 206                                                         MOD03920
  205 READ(ITAP12,6) CHG(I),EX(I),EY(I),EZ(I)                           MOD03930
  206 IF(IPRNT.LE.0) GO TO 303                                          MOD03940
      WRITE(6,7) I,CHG(I),XX(I),YY(I),ZZ(I),EX(I),EY(I),EZ(I)           MOD03950
  303 IX=(I-1)*3+1                                                      MOD03960
      IY=IX+1                                                           MOD03970
      IZ=IX+2                                                           MOD03980
      X(IX,ICASE)=XX(I)                                                 MOD03990
      X(IY,ICASE)=YY(I)                                                 MOD04000
      X(IZ,ICASE)=ZZ(I)                                                 MOD04010
      G(IX,ICASE)=EX(I)                                                 MOD04020
      G(IY,ICASE)=EY(I)                                                 MOD04030
      G(IZ,ICASE)=EZ(I)                                                 MOD04040
  103 CONTINUE                                                          MOD04050
      IF(IPRNT.LE.0) GO TO 304                                          MOD04060
      CALL DIST(NATOM,XX,YY,ZZ,50)                                      MOD04070
      WRITE(6,8)                                                        MOD04080
C*************************************************                      MOD04090
C***READ COORDINATES AND GRADIENTS OF EACH CASE***                      MOD04100
C*************************************************                      MOD04110
  304 DO 104 I=1,NINT                                                   MOD04120
      IF(ICASE.EQ.NPLUS) GO TO 207                                      MOD04130
      READ(ITAP13,9) RCOOD(I),ENGR(I)                                   MOD04140
      GO TO 208                                                         MOD04150
  207 READ(ITAP12,9) RCOOD(I),ENGR(I)                                   MOD04160
  208 CONTINUE                                                          MOD04170
      ENGR(I)=ENGR(I)*FKILO                                             MOD04180
      IF(IPRNT.LE.0) GO TO 305                                          MOD04190
      WRITE(6,10) I,RCOOD(I),ENGR(I)                                    MOD04200
  305 COOD(I,ICASE)=RCOOD(I)                                            MOD04210
      ENG(I,ICASE)=ENGR(I)                                              MOD04220
  104 CONTINUE                                                          MOD04230
C                                                                       MOD04240
C   CALCULATE NORM OF GRADIENTS                                         MOD04250
      CALL ANORM(CNORM,ENGR,150,NINT)                                   MOD04260
      IF(IPRNT.LE.0) GO TO 306                                          MOD04270
      WRITE(6,11) CNORM                                                 MOD04280
  306 GNORM(ICASE)=CNORM                                                MOD04290
  105 CONTINUE                                                          MOD04300
C                                                                       MOD04310
      DO 108 ICASE=1,NPLUS                                              MOD04320
      II=0                                                              MOD04330
      DO 107 I=1,NINT                                                   MOD04340
      DO 106 J=1,NCOOD                                                  MOD04350
      IF(NSET(J).NE.I) GO TO 106                                        MOD04360
      II=II+1                                                           MOD04370
      ENGR(II)=ENG(I,ICASE)                                             MOD04380
  106 CONTINUE                                                          MOD04390
  107 CONTINUE                                                          MOD04400
      IF(II.EQ.NCOOD) GO TO 307                                         MOD04410
      WRITE(6,12) NCOOD,II                                              MOD04420
      STOP                                                              MOD04430
  307 CONTINUE                                                          MOD04440
      CALL ANORM(CNORM,ENGR,150,NCOOD)                                  MOD04450
      IF(IPRNT.LE.0) GO TO 308                                          MOD04460
      WRITE(6,13) ICASE,CNORM                                           MOD04470
  308 SNORM(ICASE)=CNORM                                                MOD04480
  108 CONTINUE                                                          MOD04490
C                                                                       MOD04500
C   SUMMARY OF INPUT                                                    MOD04510
      IF(IPRNT.LE.0) GO TO 309                                          MOD04520
      WRITE(6,14)                                                       MOD04530
      CALL MATOUT(X,N3N,NPLUS,N3N,NPLUS,6)                              MOD04540
      WRITE(6,15)                                                       MOD04550
      CALL MATOUT(COOD,NINT,NPLUS,NINT,NPLUS,6)                         MOD04560
      WRITE(6,16)                                                       MOD04570
      CALL MATOUT(G,N3N,NPLUS,N3N,NPLUS,6)                              MOD04580
      WRITE(6,17)                                                       MOD04590
      CALL MATOUT(ENG,NINT,NPLUS,NINT,NPLUS,6)                          MOD04600
  309 WRITE(6,18)                                                       MOD04610
      DO 110 ICASE=1,NPLUS                                              MOD04620
      WRITE(6,19) ICASE,GNORM(ICASE),SNORM(ICASE),ESCF(ICASE)           MOD04630
  110 CONTINUE                                                          MOD04640
C                                                                       MOD04650
C   STORE GEOMETRIES AND GRADIENTS IN TAPE13                            MOD04660
      REWIND ITAP13                                                     MOD04670
      WRITE(ITAP13,2) NPLUS                                             MOD04680
      DO 115 ICASE=1,NPLUS                                              MOD04690
      WRITE(ITAP13,3) ICASE,ESCF(ICASE)                                 MOD04700
      DO 111 I=1,NATOM                                                  MOD04710
      IX=(I-1)*3+1                                                      MOD04720
      IY=IX+1                                                           MOD04730
      IZ=IX+2                                                           MOD04740
      XX(I)=X(IX,ICASE)                                                 MOD04750
      YY(I)=X(IY,ICASE)                                                 MOD04760
      ZZ(I)=X(IZ,ICASE)                                                 MOD04770
      WRITE(ITAP13,6) CHG(I),XX(I),YY(I),ZZ(I)                          MOD04780
  111 CONTINUE                                                          MOD04790
      DO 112 I=1,NATOM                                                  MOD04800
      IX=(I-1)*3+1                                                      MOD04810
      IY=IX+1                                                           MOD04820
      IZ=IX+2                                                           MOD04830
      EX(I)=G(IX,ICASE)                                                 MOD04840
      EY(I)=G(IY,ICASE)                                                 MOD04850
      EZ(I)=G(IZ,ICASE)                                                 MOD04860
      WRITE(ITAP13,6) CHG(I),EX(I),EY(I),EZ(I)                          MOD04870
  112 CONTINUE                                                          MOD04880
      DO 113 I=1,NINT                                                   MOD04890
      RCOOD(I)=COOD(I,ICASE)                                            MOD04900
      ENGR(I)=ENG(I,ICASE)*FMILI                                        MOD04910
      WRITE(ITAP13,9) RCOOD(I),ENGR(I)                                  MOD04920
  113 CONTINUE                                                          MOD04930
  115 CONTINUE                                                          MOD04940
C                                                                       MOD04950
      REWIND ITAP12                                                     MOD04960
      REWIND ITAP13                                                     MOD04970
C                                                                       MOD04980
      RETURN                                                            MOD04990
      END                                                               MOD05000
      SUBROUTINE DIST(NATOM,XX,YY,ZZ,NNX)                               MOD05010
      IMPLICIT REAL*8 (A-H,O-Z)                                         MOD05020
      DIMENSION XX(NNX),YY(NNX),ZZ(NNX)                                 MOD05030
      COMMON/COM107/R(1275)                                             MOD05040
      COMMON/COM108/IOFF(150),IPRNT                                     MOD05050
      DATA A00 / 0.0D+00 /                                              MOD05060
    1 FORMAT(//,2X,' INTERATOMIC DISTANCE MATRIX (IN BOHR)'/)
C                                                                       MOD05080
      IJ=0                                                              MOD05090
      DO 101 I=1,NATOM                                                  MOD05100
      DO 101 J=1,I                                                      MOD05110
      IJ=IJ+1                                                           MOD05120
      R(IJ)=A00                                                         MOD05130
      IF(I.EQ.J) GO TO 101                                              MOD05140
      XD=XX(I)-XX(J)                                                    MOD05150
      YD=YY(I)-YY(J)                                                    MOD05160
      ZD=ZZ(I)-ZZ(J)                                                    MOD05170
      RR=DSQRT(XD*XD+YD*YD+ZD*ZD)                                       MOD05180
      R(IJ)=RR                                                          MOD05190
  101 CONTINUE                                                          MOD05200
      WRITE(6,1)                                                        MOD05210
      CALL PRINT(R,1275,NATOM,6)                                        MOD05220
      RETURN                                                            MOD05230
      END                                                               MOD05240
      SUBROUTINE ANGLE(NATOM,XX,YY,ZZ,NNX)                              MOD05250
      IMPLICIT REAL*8 (A-H,O-Z)                                         MOD05260
      DIMENSION XX(NNX),YY(NNX),ZZ(NNX)                                 MOD05270
      COMMON/COM107/R(1275)                                             MOD05280
      COMMON/COM108/IOFF(150),IPRNT                                     MOD05290
      DATA PI / 3.1415926536D+00 /                                      MOD05300
    1 FORMAT(//,2X,' BOND ANGLES'/                                      MOD05310
     1 2X,' I   J   K',3X,' DEGREES   ',1X,' J   I   K',3X,' DEGREES  ',MOD05320
     2 1X,' I   K   J',3X,' DEGREES   '/)                               MOD05330
    2 FORMAT(3I4,F11.3, 2X,3I4,F11.3,1X,3I4,F11.3)                      MOD05340
    3 FORMAT(//,2X,' BOND ANGLES WILL NOT BE PRINTED---IPRNT .LT. 0 '/) MOD05350
C                                                                       MOD05360
      AR=PI/180.0D+00                                                   MOD05370
      IF(IPRNT.GE.0) WRITE(6,1)                                         MOD05380
      IF(IPRNT.LT.0) WRITE(6,3)                                         MOD05390
      DO 101 I=1,NATOM-2                                                MOD05400
      DO 101 J=I+1,NATOM-1                                              MOD05410
      IJ=I+IOFF(J)                                                      MOD05420
      XIJ=XX(I)-XX(J)                                                   MOD05430
      YIJ=YY(I)-YY(J)                                                   MOD05440
      ZIJ=ZZ(I)-ZZ(J)                                                   MOD05450
      EXIJ=-XIJ/R(IJ)                                                   MOD05460
      EYIJ=-YIJ/R(IJ)                                                   MOD05470
      EZIJ=-ZIJ/R(IJ)                                                   MOD05480
      EXJI=-EXIJ                                                        MOD05490
      EYJI=-EYIJ                                                        MOD05500
      EZJI=-EZIJ                                                        MOD05510
      DO 101 K=J+1,NATOM                                                MOD05520
      IK=I+IOFF(K)                                                      MOD05530
      XIK=XX(I)-XX(K)                                                   MOD05540
      YIK=YY(I)-YY(K)                                                   MOD05550
      ZIK=ZZ(I)-ZZ(K)                                                   MOD05560
      EXIK=-XIK/R(IK)                                                   MOD05570
      EYIK=-YIK/R(IK)                                                   MOD05580
      EZIK=-ZIK/R(IK)                                                   MOD05590
      EXKI=-EXIK                                                        MOD05600
      EYKI=-EYIK                                                        MOD05610
      EZKI=-EZIK                                                        MOD05620
      JK=J+IOFF(K)                                                      MOD05630
      XJK=XX(J)-XX(K)                                                   MOD05640
      YJK=YY(J)-YY(K)                                                   MOD05650
      ZJK=ZZ(J)-ZZ(K)                                                   MOD05660
      EXJK=-XJK/R(JK)                                                   MOD05670
      EYJK=-YJK/R(JK)                                                   MOD05680
      EZJK=-ZJK/R(JK)                                                   MOD05690
      EXKJ=-EXJK                                                        MOD05700
      EYKJ=-EYJK                                                        MOD05710
      EZKJ=-EZJK                                                        MOD05720
      AIJK=EXJI*EXJK+EYJI*EYJK+EZJI*EZJK                                MOD05730
      AJIK=EXIJ*EXIK+EYIJ*EYIK+EZIJ*EZIK                                MOD05740
      AIKJ=EXKI*EXKJ+EYKI*EYKJ+EZKI*EZKJ                                MOD05750
      ANGIJK=DACOS(AIJK)/AR                                             MOD05760
      ANGJIK=DACOS(AJIK)/AR                                             MOD05770
      ANGIKJ=DACOS(AIKJ)/AR                                             MOD05780
      IF(IPRNT.GE.0) WRITE(6,2) I,J,K,ANGIJK,J,I,K,ANGJIK,I,K,J,ANGIKJ  MOD05790
  101 CONTINUE                                                          MOD05800
C                                                                       MOD05810
      RETURN                                                            MOD05820
      END                                                               MOD05830
      SUBROUTINE SORT(SX,SG,SCOOD,SENG,X,G,COOD,ENG,NSORT)              MOD05840
      IMPLICIT REAL*8 (A-H,O-Z)                                         MOD05850
      DIMENSION SX(N3N,NPLUS),SG(N3N,NPLUS),X(N3N,NPLUS),G(N3N,NPLUS)   MOD05860
      DIMENSION SCOOD(NINT,NPLUS),SENG(NINT,NPLUS)                      MOD05870
      DIMENSION COOD(NINT,NPLUS),ENG(NINT,NPLUS)                        MOD05880
      COMMON/COM101/NATOM,N3N,NINT,NCOOD,NPLUS                          MOD05890
      COMMON/COM106/ESCF(150),GNORM(150),SNORM(150)                     MOD05900
      COMMON/COM108/IOFF(150),IPRNT                                     MOD05910
      COMMON/COM110/NSET(150),NORD(150),NXVAR(150)                      MOD05920
    1 FORMAT(//,2X,' SORTED CARTESIAN COORDINATE MATRIX'/)              MOD05930
    2 FORMAT(//,2X,' SORTED CARTESIAN GRADIENT MATRIX'/)                MOD05940
    3 FORMAT(//,2X,' SORTED INTERNAL COORDINATE MATRIX'/)               MOD05950
    4 FORMAT(//,2X,' SORTED INTERNAL GRADIENT MATRIX'/)                 MOD05960
C                                                                       MOD05970
      DO 101 I=1,NPLUS                                                  MOD05980
      NORD(I)=I                                                         MOD05990
  101 CONTINUE                                                          MOD06000
      IF(NPLUS.EQ.1) GO TO 205                                          MOD06010
      IF(NSORT.EQ.0) GO TO 205                                          MOD06020
      GO TO (201,202,203),NSORT                                         MOD06030
C                                                                       MOD06040
C   GRADIENT SORT                                                       MOD06050
  201 CALL ORDER(GNORM,NORD,150,NPLUS)                                  MOD06060
      GO TO 205                                                         MOD06070
C                                                                       MOD06080
C   ENERGY SORT                                                         MOD06090
  202 CALL ORDER(ESCF,NORD,150,NPLUS)                                   MOD06100
      GO TO 205                                                         MOD06110
C                                                                       MOD06120
C   SELECTED GRADIENT SORT                                              MOD06130
  203 CALL ORDER(SNORM,NORD,150,NPLUS)                                  MOD06140
C                                                                       MOD06150
C   NO SORT AND/OR AFTER SORT                                           MOD06160
  205 II=0                                                              MOD06170
      DO 104 I=1,NPLUS                                                  MOD06180
      II=II+1                                                           MOD06190
      MM=NORD(I)                                                        MOD06200
      DO 102 J=1,NINT                                                   MOD06210
      SCOOD(J,II)=COOD(J,MM)                                            MOD06220
      SENG(J,II)=ENG(J,MM)                                              MOD06230
  102 CONTINUE                                                          MOD06240
      DO 103 J=1,N3N                                                    MOD06250
      SX(J,II)=X(J,MM)                                                  MOD06260
      SG(J,II)=G(J,MM)                                                  MOD06270
  103 CONTINUE                                                          MOD06280
  104 CONTINUE                                                          MOD06290
      IF(IPRNT.LE.0) GO TO 206                                          MOD06300
      WRITE(6,1)                                                        MOD06310
      CALL MATOUT(SX,N3N,NPLUS,N3N,NPLUS,6)                             MOD06320
      WRITE(6,2)                                                        MOD06330
      CALL MATOUT(SG,N3N,NPLUS,N3N,NPLUS,6)                             MOD06340
  206 CONTINUE                                                          MOD06350
      WRITE(6,3)                                                        MOD06360
      CALL MATOUT(SCOOD,NINT,NPLUS,NINT,NPLUS,6)                        MOD06370
      WRITE(6,4)                                                        MOD06380
      CALL MATOUT(SENG,NINT,NPLUS,NINT,NPLUS,6)                         MOD06390
C                                                                       MOD06400
      RETURN                                                            MOD06410
      END                                                               MOD06420
      SUBROUTINE ORDER(A,NORD,NND,N)                                    MOD06430
      IMPLICIT REAL*8 (A-H,O-Z)                                         MOD06440
      DIMENSION A(NND),NORD(NND)                                        MOD06450
C                                                                       MOD06460
      N1=N-1                                                            MOD06470
      DO 102 I=1,N1                                                     MOD06480
      IR=I+1                                                            MOD06490
      DO 101 J=IR,N                                                     MOD06500
      IF(A(I).GT.A(J)) GO TO 101                                        MOD06510
      AA=A(I)                                                           MOD06520
      A(I)=A(J)                                                         MOD06530
      A(J)=AA                                                           MOD06540
      NN=NORD(I)                                                        MOD06550
      NORD(I)=NORD(J)                                                   MOD06560
      NORD(J)=NN                                                        MOD06570
  101 CONTINUE                                                          MOD06580
  102 CONTINUE                                                          MOD06590
      RETURN                                                            MOD06600
      END                                                               MOD06610
      SUBROUTINE OLDNEW(XOLD,XNEW,IFLAG)                                MOD06620
      IMPLICIT REAL*8 (A-H,O-Z)                                         MOD06630
      DIMENSION XOLD(150),XNEW(150)                                     MOD06640
      COMMON/COM101/NATOM,N3N,NINT,NCOOD,NPLUS                          MOD06650
      COMMON/COM103/CHG(50),XX(50),YY(50),ZZ(50)                        MOD06660
      COMMON/COM104/XT(50),YT(50),ZT(50)                                MOD06670
      COMMON/COM108/IOFF(150),IPRNT                                     MOD06680
      COMMON/COM109/XMAX,YMAX,ZMAX                                      MOD06690
      DATA A00 / 0.0D+00 /                                              MOD06700
    1 FORMAT(//,2X,' OLD COORDINATES'/                                  MOD06710
     1 38X,' OLD X IN BOHR'/
     2 3X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z/)               MOD06730
    2 FORMAT(3X,I5,4F16.10)                                             MOD06740
    3 FORMAT(//,2X,' NEW COORDINATES'/                                  MOD06750
     1 38X,' NEW X IN BOHR'/
     2 3X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z/)               MOD06770
C                                                                       MOD06780
      WRITE(6,1)                                                        MOD06790
      XMAX=A00                                                          MOD06800
      YMAX=A00                                                          MOD06810
      ZMAX=A00                                                          MOD06820
      DO 101 I=1,NATOM                                                  MOD06830
      IX=(I-1)*3+1                                                      MOD06840
      IY=IX+1                                                           MOD06850
      IZ=IX+2                                                           MOD06860
      IF(IFLAG.EQ.0) GO TO 201                                          MOD06870
      XDIF=XNEW(IX)-XOLD(IX)                                            MOD06880
      YDIF=XNEW(IY)-XOLD(IY)                                            MOD06890
      ZDIF=XNEW(IZ)-XOLD(IZ)                                            MOD06900
      IF(DABS(XDIF).GT.DABS(XMAX)) XMAX=XDIF                            MOD06910
      IF(DABS(YDIF).GT.DABS(YMAX)) YMAX=YDIF                            MOD06920
      IF(DABS(ZDIF).GT.DABS(ZMAX)) ZMAX=ZDIF                            MOD06930
  201 CONTINUE                                                          MOD06940
      WRITE(6,2) I,CHG(I),XOLD(IX),XOLD(IY),XOLD(IZ)                    MOD06950
  101 CONTINUE                                                          MOD06960
      WRITE(6,3)                                                        MOD06970
      DO 102 I=1,NATOM                                                  MOD06980
      IX=(I-1)*3+1                                                      MOD06990
      IY=IX+1                                                           MOD07000
      IZ=IX+2                                                           MOD07010
      WRITE(6,2) I,CHG(I),XNEW(IX),XNEW(IY),XNEW(IZ)                    MOD07020
  102 CONTINUE                                                          MOD07030
C                                                                       MOD07040
      DO 103 I=1,NATOM                                                  MOD07050
      IX=(I-1)*3+1                                                      MOD07060
      IY=IX+1                                                           MOD07070
      IZ=IX+2                                                           MOD07080
      XT(I)=XNEW(IX)                                                    MOD07090
      YT(I)=XNEW(IY)                                                    MOD07100
      ZT(I)=XNEW(IZ)                                                    MOD07110
  103 CONTINUE                                                          MOD07120
C                                                                       MOD07130
C   CALCULATE INTERATOMIC DISTANCES FOR NEXT GEOMETRY                   MOD07140
      CALL DIST(NATOM,XT,YT,ZT,50)                                      MOD07150
      IF(NATOM.LE.2) GO TO 202                                          MOD07160
      CALL ANGLE(NATOM,XT,YT,ZT,50)                                     MOD07170
C                                                                       MOD07180
  202 CONTINUE                                                          MOD07190
      IF(IFLAG.EQ.0) GO TO 203                                          MOD07200
      DO 104 I=1,NATOM                                                  MOD07210
      XX(I)=XT(I)                                                       MOD07220
      YY(I)=YT(I)                                                       MOD07230
      ZZ(I)=ZT(I)                                                       MOD07240
  104 CONTINUE                                                          MOD07250
C                                                                       MOD07260
  203 CONTINUE                                                          MOD07270
      RETURN                                                            MOD07280
      END                                                               MOD07290
      SUBROUTINE OPTIM(SX,SG,FX,HH,HT,ELH,EVALU,EE,GG,FV1,FV2,          MOD07300
     1                 ELZ,EZVAL,SAVEK,UPDATE)                          MOD07310
      IMPLICIT REAL*8 (A-H,O-Z)                                         MOD07320
      LOGICAL UPDATE                                                    MOD07330
      DIMENSION SX(N3N,NPLUS),SG(N3N,NPLUS)                             MOD07340
      DIMENSION FX(N3N,N3N)                                             MOD07350
      DIMENSION HH(NVAR,NVAR),HT(NVARS*NVARS)
      DIMENSION ELH(NVAR,NVAR),EVALU(NVAR)                              MOD07370
      DIMENSION EE(NVARS*NVARS),GG(NVARS*NVARS)                         MOD07380
      DIMENSION FV1(NVARS),FV2(NVARS)                                   MOD07390
      DIMENSION ELZ(NVARS,NVARS),EZVAL(NVARS)                           MOD07400
      DIMENSION SAVEK(NVARS)                                            MOD07410
      COMMON/COM101/NATOM,N3N,NINT,NCOOD,NPLUS                          MOD07420
      COMMON/COM102/ISTAT,NVAR,IMODE,NVARS                              MOD07430
      COMMON/COM103/CHG(50),XX(50),YY(50),ZZ(50)                        MOD07440
      COMMON/COM107/XT(150),DELX(150),DUMA(975)                         MOD07450
      COMMON/COM108/IOFF(150),IPRNT                                     MOD07460
      COMMON/COM110/NSET(150),NORD(150),NXVAR(150)                      MOD07470
      COMMON/COM111/IHUPD                                               MOD07480
      COMMON/COM112/XOLD(150),XNEW(150)                                 MOD07490
      COMMON/COM113/ENG0(150),ENG1(150)                                 MOD07500
      COMMON/COM114/EIGMAX,EIGMIN,STEPMX                                MOD07510
      DATA ITAP15 / 15 /                                                MOD07520
      DATA A00,HALF,ONE,FOUR / 0.0D+00 , 0.5D+00 , 1.0D+00 , 4.0D+00 /  MOD07530
      DATA STEP,BIG,TOLL / 0.05D+00 , 1.0D+3 , 1.0D-8 /                 MOD07540
      DATA MAXIT / 10 /                                                 MOD07550
    1 FORMAT(//,2X,61(1H*)/                                             MOD07560
     1          2X,'***OPTIMIZATION PROCEDURE+++CARTESIAN COORDINATE SYSMOD07570
     2TEM+++***'/                                                       MOD07580
     3          2X,61(1H*)//)                                           MOD07590
    2 FORMAT(2I5)                                                       MOD07600
    3 FORMAT(3F20.10)                                                   MOD07610
    4 FORMAT(//,2X,' INITITAL FX MATRIX (IN HARTREE/BOHR+2)'/)
    5 FORMAT(//,2X,' IOPTIM = ',I5/)                                    MOD07630
    6 FORMAT(//,2X,' EXTRACTED FX MATRIX (SQUARE FORM IN HARTREE/BOHR+2)
     *'/) 
    7 FORMAT(//,2X,' EXTRACTED FX MATRIX (TRIANGLE FORM IN HARTREE/BOHR+
     *2)'/)
    8 FORMAT(//,2X,' EIGENVALUES AND EIGENVECTORS OF THE HESSIAN MATRIX'MOD07660
     */)                                                                MOD07670
    9 FORMAT(//,2X,' ELX * ELX(-1) = E'/)                               MOD07680
   10 FORMAT(//,2X,' ADJUSTED EIGENVALUES AND EIGENVECTORS'/)           MOD07690
   11 FORMAT(//,2X,' ************************************'/             MOD07700
     *          2X,' * NOTICE : IMODE HAS BEEN SWITCHED *'/             MOD07710
     *          2X,' *          FROM',I5,' TO ',I5,'  *'/               MOD07720
     *          2X,' ************************************'/)            MOD07730
   12 FORMAT(//,2X,' NUMBER OF NEGATIVE EIGENVALUES FOR ACTIVE MODES = 'MOD07740
     *,I5/)                                                             MOD07750
   13 FORMAT(//,2X,' TRANSFORMED GRADIENTS'/)                           MOD07760
   14 FORMAT(2X,I5,F15.8,5X,F15.8)                                      MOD07770
   15 FORMAT(/,2X,' LAMBDA-P = ',F20.10)                                MOD07780
   16 FORMAT(//,2X,' EIGENVALUE EQUATION'/)                             MOD07790
   17 FORMAT(//,2X,' EIGENVALUES AND EIGENVECTORS OF THE EIGENVALUE EQUAMOD07800
     *TION'/)                                                           MOD07810
   18 FORMAT(2X,' LAMBDA-N = ',F20.10/)                                 MOD07820
   19 FORMAT(//,2X,' GEOMETRY CHANGE'/)                                 MOD07830
   20 FORMAT(2X,I5,2F15.8)                                              MOD07840
   21 FORMAT(//,2X,' STEP IS TOO LARGE---TO BE SCALED'/                 MOD07850
     1          2X,' STEPMX = ',F20.10/                                 MOD07860
     2          2X,' DD     = ',F20.10/                                 MOD07870
     3          2X,' SCALE  = ',F20.10/)                                MOD07880
   22 FORMAT(2X,I5,2F15.8)                                              MOD07890
   23 FORMAT(//,2X,' NEW X MATRIX'/                                     MOD07900
     1 2X,4H NO.,5X,4H II.,6X,5H OLDX,18X,5H DELX,18X,5H NEWX/)         MOD07910
   24 FORMAT(3X,I2,4X,I5,3(3X,F20.10))                                  MOD07920
C                                                                       MOD07930
      WRITE(6,1)                                                        MOD07940
C                                                                       MOD07950
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::MOD07960
C   READ IN AN INITIAL HESSIAN MATRIX W.R.T. CARTESIAN COORDINATE       MOD07970
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::MOD07980
      REWIND ITAP15                                                     MOD07990
      READ(ITAP15,2) NNATOM,NNC                                         MOD08000
      READ(ITAP15,3) ((FX(I,J),J=1,N3N),I=1,N3N)                        MOD08010
      REWIND ITAP15                                                     MOD08020
      IF(IPRNT.LE.0) GO TO 201                                          MOD08030
      WRITE(6,4)                                                        MOD08040
      CALL MATOUT(FX,N3N,N3N,N3N,N3N,6)                                 MOD08050
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::MOD08060
C                                                                       MOD08070
C   EXTRACT FORCE CONSTANTS                                             MOD08080
  201 CONTINUE                                                          MOD08090
      NVTRI=IOFF(NVAR+1)                                                MOD08100
      DO 101 I=1,NVAR                                                   MOD08110
      II=NXVAR(I)                                                       MOD08120
      DO 101 J=1,NVAR                                                   MOD08130
      JJ=NXVAR(J)                                                       MOD08140
      HH(I,J)=FX(II,JJ)                                                 MOD08150
  101 CONTINUE                                                          MOD08160
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::MOD08170
C                                                                       MOD08180
C   ITERATION FOR OPTIMIZATION                                          MOD08190
      IOPTIM=0                                                          MOD08200
  400 IOPTIM=IOPTIM+1                                                   MOD08210
      IF(IOPTIM.GT.NPLUS) GO TO 320                                     MOD08220
      WRITE(6,5) IOPTIM                                                 MOD08230
      PLAMDA=A00                                                        MOD08240
      TLAMDA=A00                                                        MOD08250
c++   micro=-1
c++   dlamb=eigmin
c+600 continue
c++   micro=micro+1
c++   fact=one
c++   if(micro.eq.1) fact=a00
c++   dincr=dlamb*fact
c++   write(3,*) ' micro = ',micro
c++   write(3,*) ' dlamb = ',dlamb
c++   write(3,*) ' dincr = ',dincr
c++   do 102 i=1,nvar
c++   hh(i,i)=hh(i,i)+dincr
c+102 continue
      IIJJ=0                                                            MOD08260
      DO 103 I=1,NVAR                                                   MOD08270
      DO 103 J=1,I                                                      MOD08280
      IIJJ=IIJJ+1                                                       MOD08290
      HT(IIJJ)=HH(I,J)                                                  MOD08300
  103 CONTINUE                                                          MOD08310
cc*   IF(IPRNT.LE.0) GO TO 204                                          MOD08320
      WRITE(6,6)                                                        MOD08330
      CALL MATOUT(HH,NVAR,NVAR,NVAR,NVAR,6)                             MOD08340
      WRITE(6,7)                                                        MOD08350
      CALLP RINT(HT,NVTRI,NVAR,6)                                       MOD08360
C                                                                       MOD08370
C   DIAGONALIZE THE EXTRACTED HESSIAN MATRIX                            MOD08380
  204 CONTINUE                                                          MOD08390
      CALL RSP(NVAR,NVAR,NVTRI,HT,EVALU,1,ELH,FV1,FV2)                  MOD08400
      WRITE(6,8)                                                        MOD08410
      CALL EIVOUT(ELH,EVALU,NVAR,NVAR,NVAR,NVAR,6)                      MOD08420
C                                                                       MOD08430
C   CHECK UNITARITY OF EIGENVECTORS                                     MOD08440
      IF(IPRNT.LE.0) GO TO 205                                          MOD08450
      CALL MTXMPY(ELH,NVAR,ELH,NVAR,EE,NVAR,GG,NVAR,NVAR,2)             MOD08460
      WRITE(6,9)                                                        MOD08470
      CALL MATOUT(EE,NVAR,NVAR,NVAR,NVAR,6)                             MOD08480
C                                                                       MOD08490
C   CHECK THE MAGNITUDE OF THE LOWEST EIGENVALUE                        MOD08500
  205 CONTINUE                                                          MOD08510
c++   lflag=0
      DO 107 I=1,NVAR                                                   MOD08520
      IF(EVALU(I).GT.EIGMIN) GO TO 206                                  MOD08530
      IF(DABS(EVALU(I)).LT.EIGMIN) THEN                                 MOD08540
c++     lflag=1
        IF(EVALU(I).LT.A00) EVALU(I)=-EIGMIN                            MOD08550
        IF(EVALU(I).GE.A00) EVALU(I)= EIGMIN                            MOD08560
      END IF                                                            MOD08570
  107 CONTINUE                                                          MOD08580
  206 CONTINUE                                                          MOD08590
c
c++ Return for Level shifting of the Hessian matrix
c++   if(lflag.ne.0) then
c++      go to 600
c++   end if
C                                                                       MOD08600
C   CHECK THE MAGNITUDE OF THE LARGEST EIGENVALUE                       MOD08610
      DO 108 I=1,NVAR                                                   MOD08620
      II=NVAR-I+1                                                       MOD08630
      IF(EVALU(II).LT.EIGMAX) GO TO 207                                 MOD08640
      IF(DABS(EVALU(II)).GT.EIGMAX) THEN                                MOD08650
        IF(EVALU(II).LT.A00) EVALU(II)=-EIGMAX                          MOD08660
        IF(EVALU(II).GE.A00) EVALU(II)= EIGMAX                          MOD08670
      END IF                                                            MOD08680
  108 CONTINUE                                                          MOD08690
  207 CONTINUE                                                          MOD08700
      IF(IPRNT.LE.0) GO TO 208                                          MOD08710
      WRITE(6,10)                                                       MOD08720
      CALL EIVOUT(ELH,EVALU,NVAR,NVAR,NVAR,NVAR,6)                      MOD08730
C                                                                       MOD08740
C   SAVE THE EIGENVECTOR OF THE MODE TO BE FOLLOWED                     MOD08750
C   FOR THE FIRST ITERATION                                             MOD08760
  208 CONTINUE                                                          MOD08770
      IF(IOPTIM.EQ.1.AND.IMODE.NE.0) THEN                               MOD08780
        DO 110 I=1,NVAR                                                 MOD08790
        SAVEK(I)=ELH(I,IMODE)                                           MOD08800
  110   CONTINUE                                                        MOD08810
      END IF                                                            MOD08820
C                                                                       MOD08830
C   FIND THE MAXIMUM OVERLAP FOR THE MODE                               MOD08840
C   FOR THE SECOND AND LATER ITERATIONS                                 MOD08850
      IF(IOPTIM.GE.2.AND.IMODE.NE.0) THEN                               MOD08860
        AXMAX=A00                                                       MOD08870
        NMODE=0                                                         MOD08880
        DO 111 I=1,NVAR                                                 MOD08890
        CALL DOTAB(SAVEK,ELH(1,I),OVER,NVAR)                            MOD08900
        IF(DABS(OVER).LE.AXMAX) GO TO 111                               MOD08910
        AXMAX=DABS(OVER)                                                MOD08920
        NMODE=I                                                         MOD08930
  111   CONTINUE                                                        MOD08940
        IF(NMODE.NE.IMODE) THEN                                         MOD08950
          WRITE(6,11) IMODE,NMODE                                       MOD08960
          IMODE=NMODE                                                   MOD08970
        END IF                                                          MOD08980
        DO 112 I=1,NVAR                                                 MOD08990
        SAVEK(I)=ELH(I,IMODE)                                           MOD09000
  112   CONTINUE                                                        MOD09010
      END IF                                                            MOD09020
C                                                                       MOD09030
C   FIND NUMBER OF NEGATIVE EIGENVALUES                                 MOD09040
      NEGEV=0                                                           MOD09050
      DO 113 I=1,NVAR                                                   MOD09060
      IF(EVALU(I).GT.A00) GO TO 113                                     MOD09070
      NEGEV=NEGEV+1                                                     MOD09080
  113 CONTINUE                                                          MOD09090
      WRITE(6,12) NEGEV                                                 MOD09100
C                                                                       MOD09110
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::MOD09120
C   TRANSFER GRADIENT W.R.T. CARTESIAN COORDINATE INTO ENG0 ARRAY       MOD09130
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::MOD09140
      NPOS=IOPTIM                                                       MOD09150
      DO 114 I=1,NVAR                                                   MOD09160
      II=NXVAR(I)                                                       MOD09170
      XT(I)=SX(II,NPOS)                                                 MOD09180
      ENG0(I)=SG(II,NPOS)                                               MOD09190
  114 CONTINUE                                                          MOD09200
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::MOD09210
C                                                                       MOD09220
C   TRANSFORM GRADIENT INTO THE LOCAL HESSIAN MODE                      MOD09230
      IF(IPRNT.GT.0)                                                    MOD09240
     *WRITE(6,13)                                                       MOD09250
      DO 116 I=1,NVAR                                                   MOD09260
      VALU=A00                                                          MOD09270
      DO 115 J=1,NVAR                                                   MOD09280
      VALU=VALU+ELH(J,I)*ENG0(J)                                        MOD09290
  115 CONTINUE                                                          MOD09300
      ENG1(I)=VALU                                                      MOD09310
      IF(IPRNT.LE.0) GO TO 116                                          MOD09320
      WRITE(6,14) I,ENG0(I),ENG1(I)                                     MOD09330
  116 CONTINUE                                                          MOD09340
C                                                                       MOD09350
C   DETERMIN SEARCH DIRECTION                                           MOD09360
C   FIND A MINIMUM POINT                                                MOD09370
      IF(ISTAT.EQ.0) THEN                                               MOD09380
        IMODE=0                                                         MOD09390
        KMODE=0                                                         MOD09400
        LMODE=1                                                         MOD09410
      END IF                                                            MOD09420
C   FIND A TRANSITION STATE                                             MOD09430
      IF(ISTAT.EQ.1) THEN                                               MOD09440
        KMODE=1                                                         MOD09450
        LMODE=2                                                         MOD09460
      END IF                                                            MOD09470
C   FOLLOW THE K-TH EIGENVECTOR                                         MOD09480
      IF(IMODE.NE.0) THEN                                               MOD09490
        KMODE=IMODE                                                     MOD09500
        LMODE=KMODE+1                                                   MOD09510
        IF(LMODE.GT.2) LMODE=1                                          MOD09520
      END IF                                                            MOD09530
C                                                                       MOD09540
C   FIND LAMBDA VALUE FOR KTH MODE (LAMBDA-P)                           MOD09550
      IF((KMODE).EQ.0) GO TO 209                                        MOD09560
      ROOTK=EVALU(KMODE)                                                MOD09570
      GRADK=ENG1(KMODE)                                                 MOD09580
      VALUK=DSQRT(ROOTK*ROOTK+GRADK*GRADK*FOUR)                         MOD09590
      PLAMDA=(ROOTK+VALUK)*HALF                                         MOD09600
      IF(IPRNT.GT.2) THEN                                               MOD09610
        WRITE(6,15) PLAMDA                                              MOD09620
      END IF                                                            MOD09630
C                                                                       MOD09640
C   FIND LAMBDA VALUE FOR OTHER MODES (LAMBDA-N)                        MOD09650
C   SET UP A SIMULTANEOUS EQUATION FOR THE MODE FOLLOWING STEP          MOD09660
c=======================================================
c   Note: the dimension of the HT array increased here
c=======================================================
  209 CONTINUE                                                          MOD09670
      NVMTRI=IOFF(NVARS)                                                MOD09680
      NVSTRI=IOFF(NVARS+1)                                              MOD09690
      CALL ZERO(HT,NVSTRI)                                              MOD09700
      II=0                                                              MOD09710
      DO 117 I=1,NVARS                                                  MOD09720
      IF(I.EQ.KMODE) GO TO 117                                          MOD09730
      II=II+1                                                           MOD09740
      III=II*(II+1)/2                                                   MOD09750
      IIJ=NVMTRI+II                                                     MOD09760
      IF(I.GT.NVAR) GO TO 117                                           MOD09770
      HT(III)=EVALU(I)                                                  MOD09780
      HT(IIJ)=ENG1(I)                                                   MOD09790
  117 CONTINUE                                                          MOD09800
      IF(IPRNT.GT.0) THEN                                               MOD09810
        WRITE(6,16)                                                     MOD09820
        CALL PRINT(HT,NVSTRI,NVARS,6)                                   MOD09830
      END IF                                                            MOD09840
C                                                                       MOD09850
C   DIAGONALIZE THE EIGENVALUE EQUATION                                 MOD09860
      CALL RSP(NVARS,NVARS,NVSTRI,HT,EZVAL,1,ELZ,FV1,FV2)               MOD09870
      IF(IPRNT.LE.0) GO TO 210                                          MOD09880
      WRITE(6,17)                                                       MOD09890
      CALL EIVOUT(ELZ,EZVAL,NVARS,NVARS,NVARS,NVARS,6)                  MOD09900
C                                                                       MOD09910
C   CHECK UNITARITY OF EIGENVECTORS                                     MOD09920
      CALL MTXMPY(ELZ,NVARS,ELZ,NVARS,EE,NVARS,GG,NVARS,NVARS,2)        MOD09930
      WRITE(6,9)                                                        MOD09940
      CALL MATOUT(EE,NVARS,NVARS,NVARS,NVARS,6)                         MOD09950
C                                                                       MOD09960
C   CHECK ACCEPTABILITY OF LAMBDA-N VALUE                               MOD09970
  210 CONTINUE                                                          MOD09980
      TLAMDA=EZVAL(1)                                                   MOD09990
      WRITE(6,15) PLAMDA                                                MOD10000
      WRITE(6,18) TLAMDA                                                MOD10010
      IF(TLAMDA.GT.EVALU(LMODE)) GO TO 220                              MOD10020
      IF(EVALU(LMODE).GT.A00.AND.TLAMDA.GT.A00) GO TO 230               MOD10030
C                                                                       MOD10040
C   EVALUATE CHANGES IN COORDINATES                                     MOD10050
      IF(IPRNT.GT.0)                                                    MOD10060
     *WRITE(6,19)                                                       MOD10070
      DO 121 I=1,NVAR                                                   MOD10080
      VALU=A00                                                          MOD10090
      DO 120 J=1,NVAR                                                   MOD10100
      ROOTJ=EVALU(J)                                                    MOD10110
      DEL=ONE/(ROOTJ-TLAMDA)                                            MOD10120
      IF(J.EQ.KMODE) THEN                                               MOD10130
        DEL=ONE/(ROOTJ-PLAMDA)                                          MOD10140
      END IF                                                            MOD10150
      VALU=VALU-(ELH(I,J)*ENG1(J))*DEL                                  MOD10160
  120 CONTINUE                                                          MOD10170
      DELX(I)=VALU                                                      MOD10180
      IF(IPRNT.LE.0) GO TO 121                                          MOD10190
      WRITE(6,20) I,VALU                                                MOD10200
  121 CONTINUE                                                          MOD10210
C                                                                       MOD10220
C   CHECK THE STEP SIZE                                                 MOD10230
      CALL ANORM(DD,DELX,150,NVAR)                                      MOD10240
      IF(DD.GT.STEPMX) THEN                                             MOD10250
        SCALE=STEPMX/DD                                                 MOD10260
        WRITE(6,21)  STEPMX,DD,SCALE                                    MOD10270
        DO 122 I=1,NVAR                                                 MOD10280
        DELOLD=DELX(I)                                                  MOD10290
        DELX(I)=DELOLD*SCALE                                            MOD10300
        IF(IPRNT.LE.0) GO TO 122                                        MOD10310
        WRITE(6,22) I,DELOLD,DELX(I)                                    MOD10320
  122   CONTINUE                                                        MOD10330
      END IF                                                            MOD10340
C                                                                       MOD10350
C   CHANGES IN GEOMETRICAL PARAMETERS                                   MOD10360
      WRITE(6,23)                                                       MOD10370
      XMAX=A00                                                          MOD10380
      DO 123 I=1,NVAR                                                   MOD10390
      II=NXVAR(I)                                                       MOD10400
      VALOLD=XT(I)                                                      MOD10410
      DMAX=DABS(DELX(I))                                                MOD10420
      IF(DMAX.GT.XMAX) XMAX=DMAX                                        MOD10430
      VALNEW=VALOLD+DELX(I)                                             MOD10440
      WRITE(6,24) I,II,VALOLD,DELX(I),VALNEW                            MOD10450
  123 CONTINUE                                                          MOD10460
C                                                                       MOD10470
      NPOS=NPOS+1                                                       MOD10480
      IF(NPOS.GT.NPLUS) GO TO 320                                       MOD10490
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::MOD10500
      DO 124 I=1,NVAR                                                   MOD10510
      II=NXVAR(I)                                                       MOD10520
      ENG1(I)=SG(II,NPOS)                                               MOD10530
  124 CONTINUE                                                          MOD10540
C                                                                       MOD10550
C   UPDATE THE HH MATRIX FOR NEXT ITERATION                             MOD10560
      CALL HUPDAT(HH,NVAR)                                              MOD10570
      GO TO 400                                                         MOD10580
C                                                                       MOD10590
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::MOD10600
C   EXCHANGE ARRAYS FOR CARTESIAN COORDINATES                           MOD10610
  320 CONTINUE                                                          MOD10620
      DO 135 I=1,NATOM                                                  MOD10630
      IX=(I-1)*3+1                                                      MOD10640
      IY=IX+1                                                           MOD10650
      IZ=IY+1                                                           MOD10660
      XOLD(IX)=XX(I)                                                    MOD10670
      XOLD(IY)=YY(I)                                                    MOD10680
      XOLD(IZ)=ZZ(I)                                                    MOD10690
      XNEW(IX)=XX(I)                                                    MOD10700
      XNEW(IY)=YY(I)                                                    MOD10710
      XNEW(IZ)=ZZ(I)                                                    MOD10720
  135 CONTINUE                                                          MOD10730
      DO 136 I=1,NVAR                                                   MOD10740
      II=NXVAR(I)                                                       MOD10750
      XNEW(II)=XOLD(II)+DELX(I)                                         MOD10760
  136 CONTINUE                                                          MOD10770
      IFLAG=0                                                           MOD10780
      IF(UPDATE) IFLAG=1                                                MOD10790
      CALL OLDNEW(XOLD,XNEW,IFLAG)                                      MOD10800
C                                                                       MOD10810
      GO TO 299                                                         MOD10820
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::MOD10830
C                                                                       MOD10840
C   LAMBDA-N > EVALU(LMODE)                                             MOD10850
  220 CONTINUE                                                          MOD10860
      WRITE(6,*) ' LAMBD-N > EVALU(LMODE) ---> STOP !'                  MOD10870
      STOP                                                              MOD10880
C                                                                       MOD10890
C   EVALU(LMODE) > 0.0 AND LAMBDA-N > 0.0                               MOD10900
  230 CONTINUE                                                          MOD10910
      WRITE(6,*) ' EVALU(LMODE) > 0.0 AND LAMBDA-N > 0.0 ---> STOP !'   MOD10920
      STOP                                                              MOD10930
C                                                                       MOD10940
  299 CONTINUE                                                          MOD10950
      RETURN                                                            MOD10960
      END                                                               MOD10970
      SUBROUTINE HUPDAT(H,NVAR)                                         MOD10980
      IMPLICIT REAL*8 (A-H,O-Z)                                         MOD10990
      DIMENSION H(NVAR,NVAR)                                            MOD11000
      COMMON/COM107/XT(150),DELX(150),DELG(150),HX(150),TK(150),        MOD11010
     1              DUMA(225)                                           MOD11020
      COMMON/COM108/IOFF(150),IPRNT                                     MOD11030
      COMMON/COM111/IHUPD                                               MOD11040
      COMMON/COM113/ENG0(150),ENG1(150)                                 MOD11050
      DATA A00,ONE / 0.0D+00 , 1.0D+00 /                                MOD11060
    1 FORMAT(//,2X,' PRESENT H MATRIX (IN HARTREE/BOHR+2)'/)
    2 FORMAT(//,2X,' DELX , DELG AND HX VECTORS'/)                      MOD11080
    3 FORMAT(2X,I2,3(5X,F15.8))                                         MOD11090
    4 FORMAT(//,2X,' DPX = ',F15.8/                                     MOD11100
     1 2X,' DHX = ',F15.8)                                              MOD11110
    5 FORMAT(//,2X,' DELX , DELG , HX AND TK VECTORS'/)                 MOD11120
    6 FORMAT(2X,I2,4(5X,F15.8))                                         MOD11130
    7 FORMAT(//,2X,' DTX = ',F15.8/                                     MOD11140
     1 2X,' DXX = ',F15.8)                                              MOD11150
    8 FORMAT(//,2X,' UPDATED H MATRIX (IN HARTREE/BOHR+2)'/)
C                                                                       MOD11170
      IF(IPRNT.LE.0) GO TO 201                                          MOD11180
      WRITE(6,1)                                                        MOD11190
      CALL MATOUT(H,NVAR,NVAR,NVAR,NVAR,6)                              MOD11200
C                                                                       MOD11210
C   UPDATE THE H MATRIX FOR NEXT ITERATION                              MOD11220
  201 CONTINUE                                                          MOD11230
      DO 101 I=1,NVAR                                                   MOD11240
      DELG(I)=ENG1(I)-ENG0(I)                                           MOD11250
  101 CONTINUE                                                          MOD11260
C                                                                       MOD11270
      GO TO (401,403,405),IHUPD                                         MOD11280
C                                                                       MOD11290
C:::::::::::::::::                                                      MOD11300
C:::BFGS METHOD:::                                                      MOD11310
C:::::::::::::::::                                                      MOD11320
  401 CONTINUE                                                          MOD11330
      DO 110 I=1,NVAR                                                   MOD11340
      HX(I)=A00                                                         MOD11350
      DO 110 J=1,NVAR                                                   MOD11360
      HX(I)=HX(I)+H(I,J)*DELX(J)                                        MOD11370
  110 CONTINUE                                                          MOD11380
      IF(IPRNT.LE.0) GO TO 402                                          MOD11390
      WRITE(6,2)                                                        MOD11400
      DO 111 I=1,NVAR                                                   MOD11410
      WRITE(6,3) I,DELX(I),DELG(I),HX(I)                                MOD11420
  111 CONTINUE                                                          MOD11430
  402 CONTINUE                                                          MOD11440
      DPX=A00                                                           MOD11450
      DHX=A00                                                           MOD11460
      DO 112 I=1,NVAR                                                   MOD11470
      DPX=DPX+DELG(I)*DELX(I)                                           MOD11480
      DHX=DHX+DELX(I)*HX(I)                                             MOD11490
  112 CONTINUE                                                          MOD11500
      DO 113 I=1,NVAR                                                   MOD11510
      DO 113 J=1,NVAR                                                   MOD11520
      TEMP1=(DELG(I)*DELG(J))/DPX                                       MOD11530
      TEMP2=-(HX(I)*HX(J))/DHX                                          MOD11540
      H(I,J)=H(I,J)+TEMP1+TEMP2                                         MOD11550
  113 CONTINUE                                                          MOD11560
      IF(IPRNT.LE.0) GO TO 410                                          MOD11570
      WRITE(6,4) DPX,DHX                                                MOD11580
      GO TO 410                                                         MOD11590
C                                                                       MOD11600
C:::::::::::::::::::                                                    MOD11610
C:::POWELL METHOD:::                                                    MOD11620
C:::::::::::::::::::                                                    MOD11630
  403 CONTINUE                                                          MOD11640
      DO 120 I=1,NVAR                                                   MOD11650
      HX(I)=A00                                                         MOD11660
      DO 120 J=1,NVAR                                                   MOD11670
      HX(I)=HX(I)+H(I,J)*DELX(J)                                        MOD11680
  120 CONTINUE                                                          MOD11690
      IF(IPRNT.GT.0)                                                    MOD11700
     *WRITE(6,5)                                                        MOD11710
      DO 121 I=1,NVAR                                                   MOD11720
      TK(I)=DELG(I)-HX(I)                                               MOD11730
      IF(IPRNT.LE.0) GO TO 121                                          MOD11740
      WRITE(6,6) I,DELX(I),DELG(I),HX(I),TK(I)                          MOD11750
  121 CONTINUE                                                          MOD11760
      DTX=A00                                                           MOD11770
      DXX=A00                                                           MOD11780
      DO 122 I=1,NVAR                                                   MOD11790
      DTX=DTX+TK(I)*DELX(I)                                             MOD11800
      DXX=DXX+DELX(I)*DELX(I)                                           MOD11810
  122 CONTINUE                                                          MOD11820
      FACT=DTX/DXX                                                      MOD11830
      DO 123 I=1,NVAR                                                   MOD11840
      DO 123 J=1,NVAR                                                   MOD11850
      TEMP1=(TK(I)*DELX(J)+DELX(I)*TK(J))/DXX                           MOD11860
      TEMP2=-(DELX(I)*DELX(J)*FACT)/DXX                                 MOD11870
      H(I,J)=H(I,J)+TEMP1+TEMP2                                         MOD11880
  123 CONTINUE                                                          MOD11890
      IF(IPRNT.LE.0) GO TO 410                                          MOD11900
      WRITE(6,7) DTX,DXX                                                MOD11910
      GO TO 410                                                         MOD11920
C                                                                       MOD11930
C###############                                                        MOD11940
C###NO UPDATE###                                                        MOD11950
C###############                                                        MOD11960
  405 CONTINUE                                                          MOD11970
C                                                                       MOD11980
C#######################                                                MOD11990
C###UPDATEED H MATRIX###                                                MOD12000
C#######################                                                MOD12010
  410 CONTINUE                                                          MOD12020
      IF(IPRNT.LE.0) GO TO 415                                          MOD12030
      WRITE(6,8)                                                        MOD12040
      CALL MATOUT(H,NVAR,NVAR,NVAR,NVAR,6)                              MOD12050
C                                                                       MOD12060
  415 CONTINUE                                                          MOD12070
      RETURN                                                            MOD12080
      END                                                               MOD12090
      SUBROUTINE ANORM(CNORM,C,NNC,NC)                                  MOD12100
      IMPLICIT REAL*8 (A-H,O-Z)                                         MOD12110
      DIMENSION C(NNC)                                                  MOD12120
      DATA A00 / 0.0D+00 /                                              MOD12130
C                                                                       MOD12140
      CNORM=A00                                                         MOD12150
      DO 101 I=1,NC                                                     MOD12160
      CNORM=CNORM+C(I)*C(I)                                             MOD12170
  101 CONTINUE                                                          MOD12180
      CNORM=DSQRT(CNORM)                                                MOD12190
      RETURN                                                            MOD12200
      END                                                               MOD12210
      SUBROUTINE DOTAB(A,B,CPROD,NC)                                    MOD12220
      IMPLICIT REAL*8 (A-H,O-Z)                                         MOD12230
      DIMENSION A(1),B(1)                                               MOD12240
      DATA A00 / 0.0D+00 /                                              MOD12250
C                                                                       MOD12260
      CPROD=A00                                                         MOD12270
      DO 101 I=1,NC                                                     MOD12280
      CPROD=CPROD+A(I)*B(I)                                             MOD12290
  101 CONTINUE                                                          MOD12300
      RETURN                                                            MOD12310
      END                                                               MOD12320
      SUBROUTINE MTXMPY(A,NAD,B,NBD,C,NCD,D,NDD,N,IMPY)                 MOD12330
      IMPLICIT REAL*8 (A-H,O-Z)                                         MOD12340
      DIMENSION A(NAD,NAD),B(NBD,NBD),C(NCD,NCD),D(NDD,NDD)             MOD12350
      DATA A00 / 0.0D+00 /                                              MOD12360
C                                                                       MOD12370
      DO 101 I=1,NCD                                                    MOD12380
      DO 101 J=1,NCD                                                    MOD12390
  101 C(I,J)=A00                                                        MOD12400
      GO TO (201,202,203,204,205,206),IMPY                              MOD12410
C     C=A*B                                                             MOD12420
  201 DO 102 I=1,N                                                      MOD12430
      DO 102 J=1,N                                                      MOD12440
      DO 102 K=1,N                                                      MOD12450
  102 C(I,J)=C(I,J)+A(I,K)*B(K,J)                                       MOD12460
      RETURN                                                            MOD12470
C     C=TA*B                                                            MOD12480
  202 DO 103 I=1,N                                                      MOD12490
      DO 103 J=1,N                                                      MOD12500
      DO 103 K=1,N                                                      MOD12510
  103 C(I,J)=C(I,J)+A(K,I)*B(K,J)                                       MOD12520
      RETURN                                                            MOD12530
C     C=A*TB                                                            MOD12540
  203 DO 104 I=1,N                                                      MOD12550
      DO 104 J=1,N                                                      MOD12560
      DO 104 K=1,N                                                      MOD12570
  104 C(I,J)=C(I,J)+A(I,K)*B(J,K)                                       MOD12580
      RETURN                                                            MOD12590
C     C=TA*TB                                                           MOD12600
  204 DO 105 I=1,N                                                      MOD12610
      DO 105 J=1,N                                                      MOD12620
      DO 105 K=1,N                                                      MOD12630
  105 C(I,J)=C(I,J)+A(K,I)*B(J,K)                                       MOD12640
      RETURN                                                            MOD12650
C     C=TA*B*A                                                          MOD12660
  205 DO 106 I=1,N                                                      MOD12670
      DO 106 J=1,N                                                      MOD12680
      D(I,J)=A00                                                        MOD12690
      DO 106 K=1,N                                                      MOD12700
  106 D(I,J)=D(I,J)+A(K,I)*B(K,J)                                       MOD12710
      DO 107 I=1,N                                                      MOD12720
      DO 107 J=1,N                                                      MOD12730
      DO 107 K=1,N                                                      MOD12740
  107 C(I,J)=C(I,J)+D(I,K)*A(K,J)                                       MOD12750
      RETURN                                                            MOD12760
C     C=A*B*TA                                                          MOD12770
  206 DO 108 I=1,N                                                      MOD12780
      DO 108 J=1,N                                                      MOD12790
      D(I,J)=A00                                                        MOD12800
      DO 108 K=1,N                                                      MOD12810
  108 D(I,J)=D(I,J)+A(I,K)*B(K,J)                                       MOD12820
      DO 109 I=1,N                                                      MOD12830
      DO 109 J=1,N                                                      MOD12840
      DO 109 K=1,N                                                      MOD12850
  109 C(I,J)=C(I,J)+D(I,K)*A(J,K)                                       MOD12860
      RETURN                                                            MOD12870
      END                                                               MOD12880
