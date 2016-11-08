C**********************************************************************
C     PROGRAM CPGRAO                                                  *
C                ( "B" FOR BIG NON-CORE CASES )                       *
C   COUPLED PERTURBED HARTREE-FOCK PROGRAM                            *
C   FOR THE GENERAL OPEN-SHELL SYSTEM IN AO BASIS                     *
C**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
C**********************************************************************
C   Last updated on April 17, 2007 by Dr. Yukio Yamaguchi             *
C**********************************************************************
C   Last updated on August 29, 2000 by Dr. Yukio Yamaguchi            *
c   ``call initmf'' and ``call tstart'' were moved to                 *
c   after ``call psiinit''                                            *
C**********************************************************************
C   Last updated on August 03, 1993 by Dr. Yukio Yamaguchi            *
c   ``call initmf'' and ``call tstart'' were moved to                 *
c   after ``call psiinit''                                            *
C**********************************************************************
C   Last updated on October 16, 1991 by Dr. Yukio Yamaguchi           *
C   Modification for the Unix Version                                 *
C*  MODIFICATION FOR IMS VERSION                                      *
C**********************************************************************
C*  BY:  YUKIO YAMAGUCHI                                              * CPG00090
C*  DATE:  FEBRUARY 01, 1989                                          * CPG00100
C********************************************************************** CPG00110
C*  BY:  RICHARD REMINGTON                         search:  c8-08-88  * CPG00120
C*  BY:  RICHARD REMINGTON                         search:  c3-24-88  * CPG00130
C*  DATE:  August  8,  1988                                           * CPG00140
C*  DATE:  MARCH  24,  1988                                           * CPG00150
C*  REASON: DECREASE CORE TO RUN IN 7MB ON 9370                       * CPG00160
C********************************************************************** CPG00170
C*  LAST UPDATED ON FEB 13, 1987 BY RICHARD REMINGTON                 * CPG00180
C*  CHANGE THE LARGE CASE "CPGRAOB" VERSION TO THE STD. "CPGRAO" AND  * CPG00190
C*  INCREASE MAXVEC,MAXVX FROM 500,12 TO 600,20 AND                   * CPG00200
C*  INCREASE CC(1200000) TO CC(1750000) AND MAXCOR TO 1750000         * CPG00210
C*  AND ZERO OUT CC USING VXINIT                                      * CPG00220
C********************************************************************** CPG00230
C*  LAST UPDATED ON JULY 19, 1985 BY YUKIO YAMAGUCHI                  * CPG00240
C********************************************************************** CPG00250
C                                                                       CPG00260
      subroutine fentry(cc,ia,maxcor)
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG00270
      dimension cc(maxcor),ia(maxcor*2)
      CHARACTER*80 ALABEL                                               CPG00310
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG00320
      COMMON/CALIF/LPARA(1024),APARA(1024)                              CPG00330
      COMMON/DIPOL/DIPDA(3,45),DIPDM(3,45),DIPDN(3,45),DIPDT(3,45)      CPG00340
      COMMON/ENRGY/ESCF,ENUC                                            CPG00350
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3          CPG00360
      COMMON/GVBAB/ALPA(10,10),BETA(10,10)                              CPG00370
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPG00380
      COMMON/MFSEC/MASTER,NSECT                                         CPG00390
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPG00400
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG00410
      COMMON/SIZES/MAXVEC,MAXVX                                         CPG00420
      COMMON/START/FOCC(10),NSORB(10),NSTR(10),NEND(10),MOTYP(256)      CPG00430
      COMMON/TAPES/ITAP36,ITAP42,ITAP43,ITAP44                          CPG00440
      COMMON/TOLER/ICONV                                                CPG00450
      COMMON/WORKS/JTAPE1                                               CPG00460
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPG00470
      COMMON/CI102/NIND,NDEP                                            CPG00480
      DATA ALABEL / ' SCF FIRST DERIVATIVE BROUGHT YOU BY Y & Y CO. EST.CPG00500
     1 IN MARCH , 1982             ' /                                  CPG00510
      DATA ZERO,QURT,HALF,ONE / 0.0D+00 , 0.25D+00 , 0.5D+00 , 1.0D+00 /CPG00520
    1 FORMAT(//,2X,' COUPLED PERTURBED HARTREE-FOCK PROGRAM FOR GENERAL CPG00530
     1OPEN-SHELL SCF IN AO BASIS (VERSION OF 4/17/07)'/)                CPG00540
    2 FORMAT(8I5)                                                       CPG00550
    3 FORMAT(2X,10I5)                                                   CPG00560
    4 FORMAT(//,2X,' PARAMETERS'/                                       CPG00570
     * 2X,' NBASIS = ',I8/                                              CPG00580
     * 2X,' NTRI   = ',I8/                                              CPG00590
     * 2X,' NBFAO  = ',I8/                                              CPG00600
     * 2X,' NBATRI = ',I8/                                              CPG00610
     * 2X,' NATOMS = ',I8/                                              CPG00620
     * 2X,' N3N    = ',I8/                                              CPG00630
     * 2X,' N3NP3  = ',I8/                                              CPG00640
     * 2X,' N3NX3  = ',I8/                                              CPG00650
     * 2X,' N3NXX  = ',I8/                                              CPG00660
     * 2X,' IOCC   = ',I8/                                              CPG00670
     * 2X,' JOCC   = ',I8/                                              CPG00680
     * 2X,' KOCC   = ',I8)                                              CPG00690
    5 FORMAT(2X,' NTYPES = ',I8/                                        CPG00700
     * 2X,' NTYPEP = ',I8/                                              CPG00710
     * 2X,' IPOL   = ',I8/                                              CPG00720
     * 2X,' ICONV  = ',I8/                                              CPG00730
     * 2X,' IPRNT  = ',I8/)                                             CPG00740
    6 FORMAT(//,2X,' SCF EIGENVECTOR (AO BASIS) MATRIX'/)               CPG00750
    7 FORMAT(//,2X,' ALPA MATRIX'/)                                     CPG00760
    8 FORMAT(//,2X,' BETA MATRIX'/)                                     CPG00770
    9 FORMAT(//,2X,' SCF EIGENVECTOR (SO BASIS) MATRIX'/)               CPG00780
   10 FORMAT(//,2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)                CPG00790
   11 FORMAT(//,2X,' REQUIRED MEMORY EXCEEDS MAXCOR'/                   CPG00800
     1          2X,' ICMAX = ',I7,5X,' MAXCOR = ',I7/)                  CPG00810
C                                                                       CPG00820
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::            CPG00830
C   REFERENCE :                                            :            CPG00840
C   Y.OSAMURA, Y.YAMAGUCHI, P.SAXE, D.J.FOX, M.A.VINCENT,  :            CPG00850
C   AND H.F.SCHAFER III, J.MOL.STRUC. 103,183 (1983)       :            CPG00860
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::            CPG00870
C                                                                       CPG00880
C   ITAP11 = GEOMETRY AND GRADIENT                                      CPG00890
C   ITAP15 = SCF SECOND DERIVATIVE                                      CPG00900
C   ITAP36 = TWO ELECTRON AO INTEGRALS                                  CPG00910
C   ITAP42 = FIRST AND SECOND DERIVATIVE INFORMATION                    CPG00920
C   ITAP43 = DERIVATIVE DIPOLE MOMENTS                                  CPG00930
C   ITAP44 = FIRST DERIVATIVES AND U'S                                  CPG00940
C   JTAPE1 = WORKING TAPE FOR THE CPHF PROCEDURE                        CPG00950
C                                                                       CPG00960
C::::::::::::::::::::::::::::::::::::::::                               CPG01010
C:::DEFINE MACHINE-DEPENDENT CONSTANTS:::                               CPG01020
C::::::::::::::::::::::::::::::::::::::::                               CPG01030
      MAXVEC=    600                                                    CPG01070
      MAXVX =     20                                                    CPG01080
C                                                                       CPG01090
      ITAPE3=3                                                          CPG01150
      INPUT=5                                                           CPG01160
      ITAPE6=6                                                          CPG01170
      ITAP11=11                                                         CPG01180
      ITAP15=15                                                         CPG01190
      ITAP17=17                                                         CPG01190
      ITAP36=36                                                         CPG01200
      ITAP42=42                                                         CPG01210
      ITAP43=43                                                         CPG01220
      ITAP44=44                                                         CPG01230
      JTAPE1=91                                                         CPG01240
C                                                                       CPG01250
      call psinit('APPEND')
      call ffile(itap11,'file11',0)
      call ffile(itap15,'file15',0)
      call ffile(itap17,'file17',0)
      CALL TSTART(6)                                                    CPG00970
      CALL INITMF(1)                                                    CPG00990
      CALL LOCATE(INPUT,'# CPHFAO #',IERR)                              CPG01260
C                                                                       CPG01270
C   LOCATION OF MATRICES IN ITAP44                                      CPG01280
C   (1) SA (N3NXX*NTRIL)                                                CPG01290
C   (2) HA (N3NXX*NTRIL)                                                CPG01300
C   (3) TA AND FA (N3NXX*NTRIL)                                         CPG01310
C   (4) EA (N3NXX*NBASQL)                                               CPG01320
C   (5) BA (N3NXX*NTRIL)                                                CPG01330
C   (6) UA (N3NXX*NBASQL)                                               CPG01340
C                                                                       CPG01350
      IOFF(1)=0                                                         CPG01360
      DO 101 I=1,255                                                    CPG01370
  101 IOFF(I+1)=IOFF(I)+I                                               CPG01380
C                                                                       CPG01390
      WRITE(6,1)                                                        CPG01400
      WRITE(3,1)                                                        CPG01410
C                                                                       CPG01420
C   READ IN PARAMETERS                                                  CPG01430
      READ(5,2) IPOL,ICONV,IPRNT                                        CPG01440
C                                                                       CPG01450
C   READ IN SCF INFORMATION                                             CPG01460
      NBASIS=LPARA(3)                                                   CPG01470
      NBFAO=LPARA(11)                                                   CPG01480
      NTRI=LPARA(4)                                                     CPG01490
      NTRI2=NTRI*2                                                      CPG01500
      NBATRI=IOFF(NBFAO+1)                                              CPG01510
      NBATR2=NBATRI*2                                                   CPG01520
      IOCC=LPARA(7)                                                     CPG01530
      KOCC=LPARA(8)                                                     CPG01540
      JOCC=LPARA(9)                                                     CPG01550
      NATOMS=LPARA(10)                                                  CPG01560
      N3N=LPARA(17)                                                     CPG01570
      NATRI=IOFF(N3N+1)                                                 CPG01580
      N3NP3=N3N+3                                                       CPG01590
      N3NX3=N3N*3                                                       CPG01600
      N3NXX=N3NP3                                                       CPG01610
      IF(IPOL.EQ.0) N3NXX=N3N                                           CPG01620
      NTRIL=((NTRI2-1)/NSECT+1)                                         CPG01630
      NBASQ=NBASIS*NBASIS*2                                             CPG01640
      NBASQL=((NBASQ-1)/NSECT+1)                                        CPG01650
      N3TRL=N3NXX*NTRIL                                                 CPG01660
      N3QRL=N3NXX*NBASQL                                                CPG01670
      NTYPES=LPARA(18)                                                  CPG01680
      NTYPEP=NTYPES+1                                                   CPG01690
      DO 102 I=1,N3NXX                                                  CPG01700
      ISA(I)=(I-1)*NTRIL+1                                              CPG01710
      IHA(I)=(I-1)*NTRIL+N3TRL+1                                        CPG01720
      IFA(I)=(I-1)*NTRIL+N3TRL+N3TRL+1                                  CPG01730
      IEA(I)=(I-1)*NBASQL+N3TRL+N3TRL+N3TRL+1                           CPG01740
      IBA(I)=(I-1)*NTRIL+N3TRL+N3TRL+N3TRL+N3QRL+1                      CPG01750
      IUA(I)=(I-1)*NBASQL+N3TRL+N3TRL+N3TRL+N3QRL+N3TRL+1               CPG01760
      WRITE(3,3) I,ISA(I),IHA(I),IFA(I),IEA(I),IBA(I),IUA(I)            CPG01770
  102 CONTINUE                                                          CPG01780
      ENUC=APARA(1)                                                     CPG01790
      ESCF=APARA(2)                                                     CPG01800
      DO 103 I=1,10                                                     CPG01810
      NSORB(I)=LPARA(I+40)                                              CPG01820
      FOCC(I)=APARA(I+40)                                               CPG01830
  103 CONTINUE                                                          CPG01840
      WRITE(6,4) NBASIS,NTRI,NBFAO,NBATRI,NATOMS,N3N,N3NP3,N3NX3,N3NXX, CPG01850
     1 IOCC,JOCC,KOCC                                                   CPG01860
      WRITE(6,5) NTYPES,NTYPEP,IPOL,ICONV,IPRNT                         CPG01870
C                                                                       CPG01880
C   READ IN AO-MO EIGENVECTORS                                          CPG01890
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      CPG01900
      IC1=1                                                             CPG01910
      IC2=IC1+NBASIS                                                    CPG01920
      IC3=IC2+NBASIS                                                    CPG01930
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      CPG01940
      CALL MREAD(CC(IC1),17)                                            CPG01950
      CALL MREAD(CC(IC2),18)                                            CPG01960
      CALL MREAD(CC(IC3),20)                                            CPG01970
      IF(IPRNT.LE.2) GO TO 301                                          CPG01980
      WRITE(6,6)                                                        CPG01990
      CALL EIGOUT(CC(IC3),CC(IC1),CC(IC2),NBFAO,NBASIS,NBFAO,NBASIS,6)  CPG02000
C                                                                       CPG02010
C   CALCULATE ALPA AND BETA MATRICES                                    CPG02020
  301 CONTINUE                                                          CPG02030
      DO 104 I=1,NTYPEP                                                 CPG02040
      DO 104 J=1,NTYPEP                                                 CPG02050
      ALPA(I,J)=ZERO                                                    CPG02060
      BETA(I,J)=ZERO                                                    CPG02070
  104 CONTINUE                                                          CPG02080
      IJ=0                                                              CPG02090
      DO 105 I=1,NTYPES                                                 CPG02100
      DO 105 J=1,I                                                      CPG02110
      IJ=IJ+1                                                           CPG02120
      ALPA(I,J)= (ONE-APARA(10+IJ))*FOCC(I)*FOCC(J)*HALF                CPG02130
      BETA(I,J)=-(ONE-APARA(25+IJ))*FOCC(I)*FOCC(J)*QURT                CPG02140
      IF(I.EQ.J) GO TO 105                                              CPG02150
      ALPA(J,I)=ALPA(I,J)                                               CPG02160
      BETA(J,I)=BETA(I,J)                                               CPG02170
  105 CONTINUE                                                          CPG02180
C     IF(IPRNT.LE.2) GO TO 302                                          CPG02190
      WRITE(6,7)                                                        CPG02200
      CALL MATOUT(ALPA,10,10,NTYPEP,NTYPEP,6)                           CPG02210
      WRITE(6,8)                                                        CPG02220
      CALL MATOUT(BETA,10,10,NTYPEP,NTYPEP,6)                           CPG02230
C                                                                       CPG02240
C**************************************                                 CPG02250
C***PREPARATION FOR FIRST ORDER CPHF***                                 CPG02260
C**************************************                                 CPG02270
C                                                                       CPG02280
C   FORM REGISTER MATRIX                                                CPG02290
  302 CONTINUE                                                          CPG02300
      WRITE(3,20)                                                       CPG02310
   20 FORMAT(//,2X,' NOW YOU ARE IN REGIST'/)                           CPG02320
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CPG02330
      IC4=IC3+NBFAO*NBASIS                                              CPG02340
      IA4=IC4+IC4-1                                                     CPG02350
      IC5=IC4+NTRI                                                      CPG02360
      IA5=IC5+IC5-1                                                     CPG02370
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CPG02380
C.................NIJ     KIJ.....                                      CPG02390
      CALL REGIST(IA(IA4),IA(IA5))                                      CPG02400
C                                                                       CPG02410
C   CHECK SIZE OF THE CORE MEMORY                                       CPG02420
      CALL CORSIZ(MAXCOR,N3NXX,NADD,MAXVC,IPOL)                         CPG02430
C                                                                       CPG02440
C   READ IN FIRST AND SECOND DERIVATIVE INFORMATION                     CPG02450
C   AND FORM SM, HM AND FM MATRICES                                     CPG02460
      WRITE(3,21)                                                       CPG02470
   21 FORMAT(//,2X,' NOW YOU ARE IN DERMAT'/)                           CPG02480
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CPG02490
      IC6=IC5+NTRI                                                      CPG02500
      IC7=IC6+N3N                                                       CPG02510
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CPG02520
      IC8=IC7+N3N*N3N                                                   CPG02530
      IC9=IC8+NBATRI                                                    CPG02540
      IC10=IC9+NBFAO*NBFAO                                              CPG02550
      ICMAX=IC10+NBFAO*NBFAO                                            CPG02560
      WRITE(3,10) ICMAX,MAXCOR                                          CPG02570
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPG02580
C.................EAO     E1A     E2A     SS      U       T........     CPG02590
      CALL DERMAT(CC(IC3),CC(IC6),CC(IC7),CC(IC8),CC(IC9),CC(IC10))     CPG02600
C                                                                       CPG02610
C   READ IN SO-MO EIGENVECTORS                                          CPG02620
      CALL MREAD(CC(IC3),19)                                            CPG02630
      IF(IPRNT.LE.2) GO TO 305                                          CPG02640
      WRITE(6,9)                                                        CPG02650
      CALL EIGOUT(CC(IC3),CC(IC1),CC(IC2),NBFAO,NBASIS,NBFAO,NBASIS,6)  CPG02660
C                                                                       CPG02670
C:::::::::::::::::::::::::::::::::::                                    CPG02680
C:::TEST ROUTINES FOR INPUT CHECK:::                                    CPG02690
C:::::::::::::::::::::::::::::::::::                                    CPG02700
  305 CONTINUE                                                          CPG02710
C                                                                       CPG02720
C   CALCULATE SCF ENERGY FOR A TEST                                     CPG02730
      WRITE(3,22)                                                       CPG02740
   22 FORMAT(//,2X,' NOW YOU ARE IN E0CAL'/)                            CPG02750
      IC8=IC7+N3N*N3N                                                   CPG02760
      IC9=IC8+NTRI                                                      CPG02770
      IC10=IC9+NTRI                                                     CPG02780
      IC11=IC10+NTRI                                                    CPG02790
      IC12=IC11+NBFAO*NBFAO                                             CPG02800
      IC13=IC12+NBFAO*NBFAO                                             CPG02810
      IC14=IC13+NTRI*NTYPES                                             CPG02820
      IC15=IC14+NTRI*NTYPES                                             CPG02830
      IC16=IC15+NTRI*NTYPES                                             CPG02840
      IA16=IC16+IC16-1                                                  CPG02850
      ICMAX=IC16+MAXBF2                                                 CPG02860
      WRITE(3,10) ICMAX,MAXCOR                                          CPG02870
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPG02880
C................OCC     ESO     HSO     HMO     TM       U........     CPG02890
      CALL E0CAL(CC(IC2),CC(IC3),CC(IC8),CC(IC9),CC(IC10),CC(IC11),     CPG02900
C................T        DP       DQ       DM       LBLI     BUFI..... CPG02910
     1           CC(IC12),CC(IC13),CC(IC14),CC(IC15),IA(IA16),CC(IC16)) CPG02920
C                                                                       CPG02930
C   CALCULATE SCF FIRST ENERGY DERIVATIVES                              CPG02940
      WRITE(3,23)                                                       CPG02950
   23 FORMAT(//,2X,' NOW YOU ARE IN EADER'/)                            CPG02960
      IC8=IC7+N3N*N3N                                                   CPG02970
      IC9=IC8+N3N                                                       CPG02980
      IC10=IC9+N3N                                                      CPG02990
      IC11=IC10+N3N                                                     CPG03000
      IC12=IC11+N3N                                                     CPG03010
      IC13=IC12+NTRI                                                    CPG03020
      IC14=IC13+NTRI                                                    CPG03030
      IC15=IC14+NTRI                                                    CPG03040
      ICMAX=IC15+NTRI*N3N*NTYPES                                        CPG03050
      WRITE(3,10) ICMAX,MAXCOR                                          CPG03060
C................OCC     D1O     D1T     D1W      E1T      SM.......    CPG03070
      CALL EADER(CC(IC2),CC(IC8),CC(IC9),CC(IC10),CC(IC11),CC(IC12),    CPG03080
C................HM       EPS      TM.......                            CPG03090
     1           CC(IC13),CC(IC14),CC(IC15))                            CPG03100
C                                                                       CPG03110
C   FORM DERIVATIVE LAGRANGIAN MATRIX                                   CPG03120
      WRITE(3,24)                                                       CPG03130
   24 FORMAT(//,2X,' NOW YOU ARE IN FAMAT'/)                            CPG03140
      IC8=IC7+N3N*N3N                                                   CPG03150
      IC9=IC8+NTRI                                                      CPG03160
      IC10=IC9+NTRI                                                     CPG03170
      IC11=IC10+NBASIS*NBASIS                                           CPG03180
      ICMAX=IC11+NTRI*N3N*NTYPEP                                        CPG03190
      WRITE(3,10) ICMAX,MAXCOR                                          CPG03200
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPG03210
C................HH      TT      EPA      FM.......                     CPG03220
      CALL FAMAT(CC(IC8),CC(IC9),CC(IC10),CC(IC11))                     CPG03230
C                                                                       CPG03240
C   CALCULATE BA MATRIX                                                 CPG03250
C   FOR INDEPENDENT PAIRS                                               CPG03260
      WRITE(3,25)                                                       CPG03270
   25 FORMAT(//,2X,' NOW YOU ARE IN BAIND'/)                            CPG03280
      IC8=IC7+N3N*N3N                                                   CPG03290
      IC9=IC8+NBFAO*NBFAO                                               CPG03300
      IC10=IC9+NBFAO*NBFAO                                              CPG03310
      IC11=IC10+NTRI*NTYPEP                                             CPG03320
      IC12=IC11+NTRI                                                    CPG03330
      IC13=IC12+NBASIS*NBASIS                                           CPG03340
      IC14=IC13+NTRI*N3N*NTYPEP                                         CPG03350
      IC15=IC14+NTRI*N3N*NTYPEP                                         CPG03360
      IC16=IC15+NTRI*N3N*NTYPEP                                         CPG03370
      IA16=IC16+IC16-1                                                  CPG03380
      ICMAX=IC16+MAXBF2                                                 CPG03390
      WRITE(3,10) ICMAX,MAXCOR                                          CPG03400
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPG03410
C................ESO     NIJ     U       T       ZETA     B0.......     CPG03420
      CALL BAIND(CC(IC3),IA(IA4),CC(IC8),CC(IC9),CC(IC10),CC(IC11),     CPG03430
C................EPA      DP       DQ       DM       SM.......          CPG03440
     1           CC(IC12),CC(IC13),CC(IC14),CC(IC15),CC(IC15),          CPG03450
C................LBLI     BUFI.....                                     CPG03460
     2           IA(IA16),CC(IC16))                                     CPG03470
      IF(IPRNT.LE.4) GO TO 308                                          CPG03480
      CALL WRBMAT(CC(IC15),N3N)                                         CPG03490
C                                                                       CPG03500
C   FOR DEPENDENT PAIRS                                                 CPG03510
  308 CONTINUE                                                          CPG03520
      WRITE(3,26)                                                       CPG03530
   26 FORMAT(//,2X,' NOW YOU ARE IN BADEP'/)                            CPG03540
      IC8=IC7+N3N*N3N                                                   CPG03550
      IC9=IC8+NBFAO*NBFAO                                               CPG03560
      IC10=IC9+NBFAO*NBFAO                                              CPG03570
      IC11=IC10+NTRI                                                    CPG03580
      IC12=IC11+NTRI                                                    CPG03590
      IC13=IC12+NTRI*N3N                                                CPG03600
      IC14=IC13+NTRI*N3N                                                CPG03610
      IC15=IC14+NTRI*N3N                                                CPG03620
      IA15=IC15+IC15-1                                                  CPG03630
      ICMAX=IC15+MAXBF2                                                 CPG03640
      WRITE(3,10) ICMAX,MAXCOR                                          CPG03650
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPG03660
C................EIG     OCC     ESO     KIJ     U       T.......       CPG03670
      CALL BADEP(CC(IC1),CC(IC2),CC(IC3),IA(IA5),CC(IC8),CC(IC9),       CPG03680
C................B0       FF       DP       DQ       DM.......          CPG03690
     1           CC(IC10),CC(IC11),CC(IC12),CC(IC13),CC(IC14),          CPG03700
C................LBLI     BUFI.....                                     CPG03710
     2           IA(IA15),CC(IC15))                                     CPG03720
      IF(IPRNT.LE.4) GO TO 309                                          CPG03730
      CALL WRBMAT(CC(IC13),N3N)                                         CPG03740
C                                                                       CPG03750
C   FORM BF MATRIX                                                      CPG03760
  309 CONTINUE                                                          CPG03770
      IF(IPOL.EQ.0) GO TO 310                                           CPG03780
      WRITE(3,27)                                                       CPG03790
   27 FORMAT(//,2X,' NOW YOU ARE IN BFMAT'/)                            CPG03800
      IC8=IC7+N3N*N3N                                                   CPG03810
      IC9=IC8+NTRI                                                      CPG03820
      IC10=IC9+NTRI*NTYPES*3                                            CPG03830
      IC11=IC10+NBASIS*NBASIS                                           CPG03840
      ICMAX=IC11+NTRI                                                   CPG03850
      WRITE(3,10) ICMAX,MAXCOR                                          CPG03860
C................NIJ     KIJ     HF      HM      EPF      BF.......     CPG03870
      CALL BFMAT(IA(IA4),IA(IA5),CC(IC8),CC(IC9),CC(IC10),CC(IC11))     CPG03880
      IF(IPRNT.LE.2) GO TO 310                                          CPG03890
      CALL WRBMAT(CC(IC11),N3NXX)                                       CPG03900
C                                                                       CPG03910
C********************************                                       CPG03920
C***FIRST ORDER CPHF PROCEDURE***                                       CPG03930
C********************************                                       CPG03940
C                                                                       CPG03950
  310 CONTINUE                                                          CPG03960
      WRITE(3,28)                                                       CPG03970
   28 FORMAT(//,2X,' NOW YOU ARE IN CPHF'/)                             CPG03980
      IC8=IC7+N3N*N3N                                                   CPG03990
      IC9=IC8+NBFAO*NBFAO                                               CPG04000
      IC10=IC9+NBFAO*NBFAO                                              CPG04010
      IC11=IC10+NTRI*NTYPEP                                             CPG04020
      IC12=IC11+NTRI                                                    CPG04030
      IC13=IC12+NTRI                                                    CPG04040
      IC14=IC13+NTRI                                                    CPG04050
      IC15=IC14+NTYPEP*(NTRI*NADD)                                      CPG04060
      IC16=IC15+NTYPEP*(NTRI*NADD)                                      CPG04070
      IC17=IC16+NTYPEP*(NTRI*NADD)                                      CPG04080
      IC18=IC17+MAXVC*(MAXVC+NADD)                                      CPG04090
      IC19=IC18+NIND*MAXVC*2                                            CPG04100
      IA19=IC19+IC19-1                                                  CPG04110
      ICMAX=IC19+MAXBF2                                                 CPG04120
      WRITE(3,10) ICMAX,MAXCOR                                          CPG04130
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPG04140
C...............ESO     NIJ     KIJ     U       A       T.......        CPG04150
      CALL CPHF(CC(IC3),IA(IA4),IA(IA5),CC(IC8),CC(IC8),CC(IC9),        CPG04160
C...............ZETA     B0       AA       TT       DP.......           CPG04170
     1          CC(IC10),CC(IC11),CC(IC12),CC(IC13),CC(IC14),           CPG04180
C...............DQ       DM       W        BB       LBLI     BUFI.....  CPG04190
     2          CC(IC15),CC(IC16),CC(IC17),CC(IC18),IA(IA19),CC(IC19),  CPG04200
C...........................                                            CPG04210
     3          MAXVC,NADD,N3NXX)                                       CPG04220
      IF(IPRNT.LE.4) GO TO 311                                          CPG04230
      CALL WRBMAT(CC(IC18),N3NXX)                                       CPG04240
C                                                                       CPG04250
C   CALCULATE DEPENDENT PAIRS OF U MATRIX                               CPG04260
  311 CONTINUE                                                          CPG04270
      WRITE(3,29)                                                       CPG04280
   29 FORMAT(//,2X,' NOW YOU ARE IN UCMAT'/)                            CPG04290
      IC8=IC7+N3N*N3N                                                   CPG04300
      IC9=IC8+NBFAO*NBFAO                                               CPG04310
      IC10=IC9+NBFAO*NBFAO                                              CPG04320
      IC11=IC10+NTRI                                                    CPG04330
      IC12=IC11+NTRI*N3NXX                                              CPG04340
      IC13=IC12+NTRI*N3NXX                                              CPG04350
      IC14=IC13+NTRI*N3NXX                                              CPG04360
      IC15=IC14+NTRI                                                    CPG04370
      IA15=IC15+IC15-1                                                  CPG04380
      ICMAX=IC15+MAXBF2                                                 CPG04390
      WRITE(3,10) ICMAX,MAXCOR                                          CPG04400
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPG04410
C................EIG     OCC     ESO     NIJ     KIJ     U.......       CPG04420
      CALL UCMAT(CC(IC1),CC(IC2),CC(IC3),IA(IA4),IA(IA5),CC(IC8),       CPG04430
C................T       B        DP       DQ       DM       DEL......  CPG04440
     1           CC(IC9),CC(IC10),CC(IC11),CC(IC12),CC(IC13),CC(IC14),  CPG04450
C................LBLI     BUFI...........                               CPG04460
     2           IA(IA15),CC(IC15),N3NXX)                               CPG04470
      IF(IPRNT.LE.2) GO TO 312                                          CPG04480
      CALL WRBMAT(CC(IC13),N3NXX)                                       CPG04490
C                                                                       CPG04500
C   COMPLETE U MATRIX                                                   CPG04510
  312 CONTINUE                                                          CPG04520
      WRITE(3,30)                                                       CPG04530
   30 FORMAT(//,2X,' NOW YOU ARE IN UXMAT'/)                            CPG04540
      IC8=IC7+N3N*N3N                                                   CPG04550
      IC9=IC8+NTRI*N3NXX                                                CPG04560
      IC10=IC9+NTRI*N3NXX                                               CPG04570
      ICMAX=IC10+NBASIS*NBASIS                                          CPG04580
      WRITE(3,10) ICMAX,MAXCOR                                          CPG04590
C................B       SM      U..............                        CPG04600
      CALL UXMAT(CC(IC8),CC(IC9),CC(IC10),N3NXX)                        CPG04610
C     CALL WRBMAT(CC(IC10),N3NXX)                                       CPG04620
C                                                                       CPG04630
C   CALCULATE SCF SECOND DERIVATIVES                                    CPG04640
      WRITE(3,31)                                                       CPG04650
   31 FORMAT(//,2X,' NOW YOU ARE IN WAMAT'/)                            CPG04660
      IC8=IC7+N3N*N3N                                                   CPG04670
      IC9=IC8+NBASIS*NBASIS                                             CPG04680
      IC10=IC9+NBASIS*NBASIS                                            CPG04690
      IC11=IC10+NBASIS*NBASIS                                           CPG04700
      IC12=IC11+NTRI*NTYPEP                                             CPG04710
      IC13=IC12+NBASIS*NBASIS                                           CPG04720
      IC14=IC13+NTRI*N3N*NTYPEP                                         CPG04730
      IC15=IC14+NTRI*N3N*NTYPEP                                         CPG04740
      IC16=IC15+NTRI*N3N*NTYPEP                                         CPG04750
      IA16=IC16+IC16-1                                                  CPG04760
      ICMAX=IC16+MAXBF2                                                 CPG04770
      WRITE(3,10) ICMAX,MAXCOR                                          CPG04780
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPG04790
C................ESO     U       T       EPA      ZETA     WA.......    CPG04800
      CALL WAMAT(CC(IC3),CC(IC8),CC(IC9),CC(IC10),CC(IC11),CC(IC12),    CPG04810
C................SA       DP       DQ       DM       LBLI     BUFI..... CPG04820
     1           CC(IC13),CC(IC14),CC(IC15),CC(IC13),IA(IA16),CC(IC16)) CPG04830
C                                                                       CPG04840
      WRITE(3,32)                                                       CPG04850
   32 FORMAT(//,2X,' NOW YOU ARE IN SCF2ND'/)                           CPG04860
      IC8=IC7+N3N*N3N                                                   CPG04870
      IC9=IC8+N3N*N3N                                                   CPG04880
      IC10=IC9+N3N*N3N                                                  CPG04890
      IC11=IC10+NTRI*N3N                                                CPG04900
      IC12=IC11+NBASIS*NBASIS*N3N                                       CPG04910
      IC13=IC12+NBASIS*NBASIS                                           CPG04920
      ICMAX=IC13+NBASIS*NBASIS                                          CPG04930
      WRITE(3,10) ICMAX,MAXCOR                                          CPG04940
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPG04950
C.................E2A     E2M     E2P     SA       WA       EPA......   CPG04960
      CALL SCF2ND(CC(IC7),CC(IC8),CC(IC9),CC(IC10),CC(IC11),CC(IC12),   CPG04970
C.................UA.......                                             CPG04980
     1            CC(IC13))                                             CPG04990
C                                                                       CPG05000
      IF(IPOL.EQ.0) GO TO 320                                           CPG05010
      WRITE(3,33)                                                       CPG05020
   33 FORMAT(//,2X,' NOW YOU ARE IN DIPMO'/)                            CPG05030
      IC8=IC7+N3N*N3N                                                   CPG05040
      IC9=IC8+NTRI*3                                                    CPG05050
      ICMAX=IC9+NBASIS*NBASIS*3                                         CPG05060
      WRITE(3,10) ICMAX,MAXCOR                                          CPG05070
C................OCC     DIP     UA......                               CPG05080
      CALL DIPMO(CC(IC2),CC(IC8),CC(IC9))                               CPG05090
C                                                                       CPG05100
      WRITE(3,34)                                                       CPG05110
   34 FORMAT(//,2X,' NOW YOU ARE IN POLAR'/)                            CPG05120
      IC8=IC7+N3N*N3N                                                   CPG05130
      IC9=IC8+NTRI*3                                                    CPG05140
      ICMAX=IC9+NBASIS*NBASIS                                           CPG05150
      WRITE(3,10) ICMAX,MAXCOR                                          CPG05160
C................OCC     HF      UF......                               CPG05170
      CALL POLAR(CC(IC2),CC(IC8),CC(IC9))                               CPG05180
  320 CONTINUE                                                          CPG05190
      GO TO 400                                                         CPG05200
C                                                                       CPG05210
  399 WRITE(6,11) ICMAX,MAXCOR                                          CPG05220
  400 CONTINUE                                                          CPG05230
C                                                                       CPG05240
      CALL RCLOSE(ITAP42,3)                                             CPG05250
      CALL RCLOSE(ITAP44,3)                                             CPG05270
      CALL RCLOSE(MASTER,3)                                             CPG05280
      CALL TSTOP(6)                                                     CPG05290
C                                                                       CPG05300
      STOP                                                              CPG05310
      END                                                               CPG05320
      SUBROUTINE REGIST(NIJ,KIJ)                                        CPG05330
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG05340
      DIMENSION NIJ(NTRI,2),KIJ(NTRI,2)                                 CPG05350
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG05360
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3          CPG05370
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG05380
      COMMON/START/FOCC(10),NSORB(10),NSTR(10),NEND(10),MOTYP(256)      CPG05390
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPG05400
      COMMON/CI102/NIND,NDEP                                            CPG05410
    1 FORMAT(//,2X,' REGISTERED NUMBERS'/                               CPG05420
     1 2X,' KVIR = ',I5/                                                CPG05430
     2 2X,' NIND = ',I5/                                                CPG05440
     3 2X,' NDEP = ',I5/                                                CPG05450
     4 2X,' NDIF = ',I5)                                                CPG05460
    2 FORMAT(//,2X,' NIJ MATRIX'/2X,' NO.',5X,' NIJ(1)',5X,' NIJ(2)'/)  CPG05470
    3 FORMAT(2X,I4,4X,I5,7X,I5)                                         CPG05480
    4 FORMAT(//,2X,' KIJ MATRIX'/2X,' NO.',5X,' KIJ(1)',5X,' KIJ(2)'/)  CPG05490
    5 FORMAT(//,2X,' NO.',8X,' FOCC',6X,' NSORB',7X,' NSTR',7X,' NEND'/)CPG05500
    6 FORMAT(2X,I4,3X,F10.5,7X,I5,7X,I5,7X,I5)                          CPG05510
C                                                                       CPG05520
      KVIR=NBASIS-JOCC                                                  CPG05530
C                                                                       CPG05540
      NSTR(1)=1                                                         CPG05550
      NEND(1)=NSORB(1)                                                  CPG05560
      DO 101 ITYP=2,NTYPEP                                              CPG05570
      NSTR(ITYP)=NEND(ITYP-1)+1                                         CPG05580
      NEND(ITYP)=NSTR(ITYP)+NSORB(ITYP)-1                               CPG05590
  101 CONTINUE                                                          CPG05600
      DO 103 ITYP=1,NTYPEP                                              CPG05610
      DO 102 I=NSTR(ITYP),NEND(ITYP)                                    CPG05620
  102 MOTYP(I)=ITYP                                                     CPG05630
  103 CONTINUE                                                          CPG05640
C                                                                       CPG05650
      II=0                                                              CPG05660
      DO 105 ITYP=2,NTYPEP                                              CPG05670
      DO 105 JTYP=1,ITYP-1                                              CPG05680
      DO 104 I=NSTR(ITYP),NEND(ITYP)                                    CPG05690
      DO 104 J=NSTR(JTYP),NEND(JTYP)                                    CPG05700
      II=II+1                                                           CPG05710
      NIJ(II,1)=I                                                       CPG05720
      NIJ(II,2)=J                                                       CPG05730
  104 CONTINUE                                                          CPG05740
  105 CONTINUE                                                          CPG05750
      NIND=II                                                           CPG05760
C                                                                       CPG05770
      II=0                                                              CPG05780
      DO 107 ITYP=1,NTYPEP                                              CPG05790
C*C   IF(NSORB(ITYP).EQ.1) GO TO 107                                    CPG05800
C*C   DO 106 I=NSTR(ITYP)+1,NEND(ITYP)                                  CPG05810
C*C   DO 106 J=NSTR(ITYP),I-1                                           CPG05820
      DO 106 I=NSTR(ITYP),NEND(ITYP)                                    CPG05830
      DO 106 J=NSTR(ITYP),I                                             CPG05840
      II=II+1                                                           CPG05850
      KIJ(II,1)=I                                                       CPG05860
      KIJ(II,2)=J                                                       CPG05870
  106 CONTINUE                                                          CPG05880
  107 CONTINUE                                                          CPG05890
      NDEP=II                                                           CPG05900
C                                                                       CPG05910
      NDIF=NTRI-NIND-NDEP                                               CPG05920
      IF(IPRNT.LE.2) GO TO 205                                          CPG05930
      WRITE(6,1) KVIR,NIND,NDEP,NDIF                                    CPG05940
      WRITE(6,2)                                                        CPG05950
      DO 108 I=1,NIND                                                   CPG05960
  108 WRITE(6,3) I,NIJ(I,1),NIJ(I,2)                                    CPG05970
      WRITE(6,4)                                                        CPG05980
      DO 109 I=1,NDEP                                                   CPG05990
  109 WRITE(6,3) I,KIJ(I,1),KIJ(I,2)                                    CPG06000
      WRITE(6,5)                                                        CPG06010
      DO 110 I=1,NTYPEP                                                 CPG06020
  110 WRITE(6,6) I,FOCC(I),NSORB(I),NSTR(I),NEND(I)                     CPG06030
C                                                                       CPG06040
  205 RETURN                                                            CPG06050
      END                                                               CPG06060
      SUBROUTINE CORSIZ(MAXCOR,N3NXX,NADD,MAXVC,IPOL)                   CPG06070
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG06080
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG06090
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3          CPG06100
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPG06110
      COMMON/SIZES/MAXVEC,MAXVX                                         CPG06120
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPG06130
      COMMON/CI102/NIND,NDEP                                            CPG06140
    1 FORMAT(//,2X,' %%%%%%%%%%%%%%%%%%%%%%%%%%%'/                      CPG06150
     1          2X,' %%%NOW YOU ARE IN CORSIZ%%%'/                      CPG06160
     2          2X,' %%%%%%%%%%%%%%%%%%%%%%%%%%%'/)                     CPG06170
    2 FORMAT(//,2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)                CPG06180
    3 FORMAT(//,2X,' *************************************************'/CPG06190
     1          2X,' ***REQUIRED MEMORY EXCEEDS MAXCOR-----STOP !!!***'/CPG06200
     2          2X,' *************************************************'/CPG06210
     3          2X,' IGMAX = ',I10,5X,' MAXCOR = ',I10/)                CPG06220
    4 FORMAT(//,2X,' ::::::::::::::::::::::::::::::::::'/               CPG06230
     1          2X,' :::REQUIRED MEMORY FOR THIS RUN:::'/               CPG06240
     2          2X,' ::::::::::::::::::::::::::::::::::'/               CPG06250
     3          2X,' NADD  = ',I10,5X,' MAXVC  = ',I10/                 CPG06260
     4          2X,' IGMAX = ',I10,5X,' MAXCOR = ',I10/)                CPG06270
    5 FORMAT(//,2X,' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'/                 CPG06280
     1          2X,' %%%NOW YOU ARE LEAVING CORSIZ%%%'/                 CPG06290
     2          2X,' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'/)                CPG06300
C                                                                       CPG06310
C   THIS SUBROUTINE WILL CHECK THE MAXIMUM CORE MEMORY REQUIRED         CPG06320
C   FOR THE CALCULATION                                                 CPG06330
C                                                                       CPG06340
C   DYNAMIC ALLOCATION FOR THE CORE MEMORY                              CPG06350
C   1 : EIG        (NBASIS)                                             CPG06360
C   2 : OCC        (NBASIS)                                             CPG06370
C   3 : ESO & EAO  (NBFAO*NBASIS)                                       CPG06380
C   4 : E1A        (N3N)                                                CPG06390
C   5 : E2A        (N3N*N3N)                                            CPG06400
C   6 : NIJ        (NTRI)                                               CPG06410
C   7 : KIJ        (NTRI)                                               CPG06420
C                                                                       CPG06430
      WRITE(3,1)                                                        CPG06440
      NADD=N3NXX*2                                                      CPG06450
C                                                                       CPG06460
C   READ IN AO-MO EIGENVECTORS                                          CPG06470
      IC1=1                                                             CPG06480
      IC2=IC1+NBASIS                                                    CPG06490
      IC3=IC2+NBASIS                                                    CPG06500
C                                                                       CPG06510
C   FORM REGISTER MATRIX                                                CPG06520
      WRITE(3,20)                                                       CPG06530
   20 FORMAT(//,2X,' NOW YOU ARE IN REGIST'/)                           CPG06540
      IC4=IC3+NBFAO*NBASIS                                              CPG06550
      IA4=IC4+IC4-1                                                     CPG06560
      IC5=IC4+NTRI                                                      CPG06570
      IA5=IC5+IC5-1                                                     CPG06580
C.................NIJ     KIJ.....                                      CPG06590
C***  CALL REGIST(IA(IA4),IA(IA5))                                      CPG06600
C                                                                       CPG06610
C   READ IN FIRST AND SECOND DERIVATIVE INFORMATION                     CPG06620
      WRITE(3,21)                                                       CPG06630
   21 FORMAT(//,2X,' NOW YOU ARE IN DERMAT'/)                           CPG06640
      IC6=IC5+NTRI                                                      CPG06650
      IC7=IC6+N3N                                                       CPG06660
      IC8=IC7+N3N*N3N                                                   CPG06670
      IC9=IC8+NBATRI                                                    CPG06680
      IC10=IC9+NBFAO*NBFAO                                              CPG06690
      ICMAX=IC10+NBFAO*NBFAO                                            CPG06700
      WRITE(3,2) ICMAX,MAXCOR                                           CPG06710
      IGMAX=ICMAX                                                       CPG06720
C.................EAO     E1A     E2A     SS      U       T........     CPG06730
C***  CALL DERMAT(CC(IC3),CC(IC6),CC(IC7),CC(IC8),CC(IC9),CC(IC10))     CPG06740
C                                                                       CPG06750
C   CALCULATE SCF ENERGY FOR A TEST                                     CPG06760
      WRITE(3,22)                                                       CPG06770
   22 FORMAT(//,2X,' NOW YOU ARE IN E0CAL'/)                            CPG06780
      IC8=IC7+N3N*N3N                                                   CPG06790
      IC9=IC8+NTRI                                                      CPG06800
      IC10=IC9+NTRI                                                     CPG06810
      IC11=IC10+NTRI                                                    CPG06820
      IC12=IC11+NBFAO*NBFAO                                             CPG06830
      IC13=IC12+NBFAO*NBFAO                                             CPG06840
      IC14=IC13+NTRI*NTYPES                                             CPG06850
      IC15=IC14+NTRI*NTYPES                                             CPG06860
      IC16=IC15+NTRI*NTYPES                                             CPG06870
      IA16=IC16+IC16-1                                                  CPG06880
      ICMAX=IC16+MAXBF2                                                 CPG06890
      WRITE(3,2) ICMAX,MAXCOR                                           CPG06900
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPG06910
C................OCC     ESO     HSO     HMO     TM       U........     CPG06920
C***  CALL E0CAL(CC(IC2),CC(IC3),CC(IC8),CC(IC9),CC(IC10),CC(IC11),     CPG06930
C................T        DP       DQ       DM       LBLI     BUFI..... CPG06940
C*** 1           CC(IC12),CC(IC13),CC(IC14),CC(IC15),IA(IA16),CC(IC16)) CPG06950
C                                                                       CPG06960
C   CALCULATE SCF FIRST ENERGY DERIVATIVES                              CPG06970
      WRITE(3,23)                                                       CPG06980
   23 FORMAT(//,2X,' NOW YOU ARE IN EADER'/)                            CPG06990
      IC8=IC7+N3N*N3N                                                   CPG07000
      IC9=IC8+N3N                                                       CPG07010
      IC10=IC9+N3N                                                      CPG07020
      IC11=IC10+N3N                                                     CPG07030
      IC12=IC11+N3N                                                     CPG07040
      IC13=IC12+NTRI                                                    CPG07050
      IC14=IC13+NTRI                                                    CPG07060
      IC15=IC14+NTRI                                                    CPG07070
      ICMAX=IC15+NTRI*N3N*NTYPES                                        CPG07080
      WRITE(3,2) ICMAX,MAXCOR                                           CPG07090
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPG07100
C................OCC     D1O     D1T     D1W      E1T      SM.......    CPG07110
C***  CALL EADER(CC(IC2),CC(IC8),CC(IC9),CC(IC10),CC(IC11),CC(IC12),    CPG07120
C................HM       EPS      TM.......                            CPG07130
C*** 1           CC(IC13),CC(IC14),CC(IC15))                            CPG07140
C                                                                       CPG07150
C   FORM DERIVATIVE LAGRANGIAN MATRIX                                   CPG07160
      WRITE(3,24)                                                       CPG07170
   24 FORMAT(//,2X,' NOW YOU ARE IN FAMAT'/)                            CPG07180
      IC8=IC7+N3N*N3N                                                   CPG07190
      IC9=IC8+NTRI                                                      CPG07200
      IC10=IC9+NTRI                                                     CPG07210
      IC11=IC10+NBASIS*NBASIS                                           CPG07220
      ICMAX=IC11+NTRI*N3N*NTYPEP                                        CPG07230
      WRITE(3,2) ICMAX,MAXCOR                                           CPG07240
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPG07250
C................HH      TT      EPA      FM.......                     CPG07260
C***  CALL FAMAT(CC(IC8),CC(IC9),CC(IC10),CC(IC11))                     CPG07270
C                                                                       CPG07280
C   CALCULATE BA MATRIX                                                 CPG07290
C   FOR INDEPENDENT PAIRS                                               CPG07300
      WRITE(3,25)                                                       CPG07310
   25 FORMAT(//,2X,' NOW YOU ARE IN BAIND'/)                            CPG07320
      IC8=IC7+N3N*N3N                                                   CPG07330
      IC9=IC8+NBFAO*NBFAO                                               CPG07340
      IC10=IC9+NBFAO*NBFAO                                              CPG07350
      IC11=IC10+NTRI*NTYPEP                                             CPG07360
      IC12=IC11+NTRI                                                    CPG07370
      IC13=IC12+NBASIS*NBASIS                                           CPG07380
      IC14=IC13+NTRI*N3N*NTYPEP                                         CPG07390
      IC15=IC14+NTRI*N3N*NTYPEP                                         CPG07400
      IC16=IC15+NTRI*N3N*NTYPEP                                         CPG07410
      IA16=IC16+IC16-1                                                  CPG07420
      ICMAX=IC16+MAXBF2                                                 CPG07430
      WRITE(3,2) ICMAX,MAXCOR                                           CPG07440
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPG07450
C................ESO     NIJ     U       T       ZETA     B0.......     CPG07460
C***  CALL BAIND(CC(IC3),IA(IA4),CC(IC8),CC(IC9),CC(IC10),CC(IC11),     CPG07470
C................EPA      DP       DQ       DM       SM.......          CPG07480
C*** 1           CC(IC12),CC(IC13),CC(IC14),CC(IC15),CC(IC15),          CPG07490
C................LBLI     BUFI.....                                     CPG07500
C*** 2           IA(IA16),CC(IC16))                                     CPG07510
C                                                                       CPG07520
C   FOR DEPENDENT PAIRS                                                 CPG07530
      WRITE(3,26)                                                       CPG07540
   26 FORMAT(//,2X,' NOW YOU ARE IN BADEP'/)                            CPG07550
      IC8=IC7+N3N*N3N                                                   CPG07560
      IC9=IC8+NBFAO*NBFAO                                               CPG07570
      IC10=IC9+NBFAO*NBFAO                                              CPG07580
      IC11=IC10+NTRI                                                    CPG07590
      IC12=IC11+NTRI                                                    CPG07600
      IC13=IC12+NTRI*N3N                                                CPG07610
      IC14=IC13+NTRI*N3N                                                CPG07620
      IC15=IC14+NTRI*N3N                                                CPG07630
      IA15=IC15+IC15-1                                                  CPG07640
      ICMAX=IC15+MAXBF2                                                 CPG07650
      WRITE(3,2) ICMAX,MAXCOR                                           CPG07660
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPG07670
C................EIG     OCC     ESO     KIJ     U       T.......       CPG07680
C***  CALL BADEP(CC(IC1),CC(IC2),CC(IC3),IA(IA5),CC(IC8),CC(IC9),       CPG07690
C................B0       FF       DP       DQ       DM.......          CPG07700
C*** 1           CC(IC10),CC(IC11),CC(IC12),CC(IC13),CC(IC14),          CPG07710
C................LBLI     BUFI.....                                     CPG07720
C*** 2           IA(IA15),CC(IC15))                                     CPG07730
C                                                                       CPG07740
C   FORM BF MATRIX                                                      CPG07750
      IF(IPOL.EQ.0) GO TO 310                                           CPG07760
      WRITE(3,27)                                                       CPG07770
   27 FORMAT(//,2X,' NOW YOU ARE IN BFMAT'/)                            CPG07780
      IC8=IC7+N3N*N3N                                                   CPG07790
      IC9=IC8+NTRI                                                      CPG07800
      IC10=IC9+NTRI*NTYPES*3                                            CPG07810
      IC11=IC10+NBASIS*NBASIS                                           CPG07820
      ICMAX=IC11+NTRI                                                   CPG07830
      WRITE(3,2) ICMAX,MAXCOR                                           CPG07840
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPG07850
C................NIJ     KIJ     HF      HM      EPF      BF.......     CPG07860
C***  CALL BFMAT(IA(IA4),IA(IA5),CC(IC8),CC(IC9),CC(IC10),CC(IC11))     CPG07870
C                                                                       CPG07880
  310 CONTINUE                                                          CPG07890
      WRITE(3,28)                                                       CPG07900
   28 FORMAT(//,2X,' NOW YOU ARE IN CPHF'/)                             CPG07910
  200 CONTINUE                                                          CPG07920
      NADD=NADD/2                                                       CPG07930
      MAXVC=MIN0(MAXVEC,NADD*MAXVX)                                     CPG07940
      IC8=IC7+N3N*N3N                                                   CPG07950
      IC9=IC8+NBFAO*NBFAO                                               CPG07960
      IC10=IC9+NBFAO*NBFAO                                              CPG07970
      IC11=IC10+NTRI*NTYPEP                                             CPG07980
      IC12=IC11+NTRI                                                    CPG07990
      IC13=IC12+NTRI                                                    CPG08000
      IC14=IC13+NTRI                                                    CPG08010
      IC15=IC14+NTYPEP*(NTRI*NADD)                                      CPG08020
      IC16=IC15+NTYPEP*(NTRI*NADD)                                      CPG08030
      IC17=IC16+NTYPEP*(NTRI*NADD)                                      CPG08040
      IC18=IC17+MAXVC*(MAXVC+NADD)                                      CPG08050
      IC19=IC18+NIND*MAXVC*2                                            CPG08060
      IA19=IC19+IC19-1                                                  CPG08070
      ICMAX=IC19+MAXBF2                                                 CPG08080
      WRITE(3,29) NIND,N3NXX,NADD,MAXVC,ICMAX                           CPG08090
   29 FORMAT(2X,' NIND  = ',I6,3X,' N3NXX = ',I6,3X,' NADD = ',I6/      CPG08100
     1       2X,' MAXVC = ',I6,3X,' ICMAX = ',I6/)                      CPG08110
      IF(NADD.EQ.1) GO TO 201                                           CPG08120
      IF(ICMAX.LT.MAXCOR) GO TO 201                                     CPG08130
      GO TO 200                                                         CPG08140
  201 CONTINUE                                                          CPG08150
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPG08160
C...............ESO     NIJ     KIJ     U       A       T.......        CPG08170
C***  CALL CPHF(CC(IC3),IA(IA4),IA(IA5),CC(IC8),CC(IC8),CC(IC9),        CPG08180
C...............ZETA     B0       AA       TT       DP.......           CPG08190
C*** 1          CC(IC10),CC(IC11),CC(IC12),CC(IC13),CC(IC14),           CPG08200
C...............DQ       DM       W        BB       LBLI     BUFI.....  CPG08210
C*** 2          CC(IC15),CC(IC16),CC(IC17),CC(IC18),IA(IA19),CC(IC19),  CPG08220
C...........................                                            CPG08230
C*** 3          MAXVC,N3NXX)                                            CPG08240
C                                                                       CPG08250
C   CALCULATE DEPENDENT PAIRS OF U MATRIX                               CPG08260
      WRITE(3,30)                                                       CPG08270
   30 FORMAT(//,2X,' NOW YOU ARE IN UCMAT'/)                            CPG08280
      IC8=IC7+N3N*N3N                                                   CPG08290
      IC9=IC8+NBFAO*NBFAO                                               CPG08300
      IC10=IC9+NBFAO*NBFAO                                              CPG08310
      IC11=IC10+NTRI                                                    CPG08320
      IC12=IC11+NTRI*N3NXX                                              CPG08330
      IC13=IC12+NTRI*N3NXX                                              CPG08340
      IC14=IC13+NTRI*N3NXX                                              CPG08350
      IC15=IC14+NTRI                                                    CPG08360
      IA15=IC15+IC15-1                                                  CPG08370
      ICMAX=IC15+MAXBF2                                                 CPG08380
      WRITE(3,2) ICMAX,MAXCOR                                           CPG08390
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPG08400
C................EIG     OCC     ESO     NIJ     KIJ     U.......       CPG08410
C***  CALL UCMAT(CC(IC1),CC(IC2),CC(IC3),IA(IA4),IA(IA5),CC(IC8),       CPG08420
C................T       B        DP       DQ       DM       DEL......  CPG08430
C*** 1           CC(IC9),CC(IC10),CC(IC11),CC(IC12),CC(IC13),CC(IC14),  CPG08440
C................LBLI     BUFI...........                               CPG08450
C*** 2           IA(IA15),CC(IC15),N3NXX)                               CPG08460
C                                                                       CPG08470
C   COMPLETE U MATRIX                                                   CPG08480
      WRITE(3,31)                                                       CPG08490
   31 FORMAT(//,2X,' NOW YOU ARE IN UXMAT'/)                            CPG08500
      IC8=IC7+N3N*N3N                                                   CPG08510
      IC9=IC8+NTRI*N3NXX                                                CPG08520
      IC10=IC9+NTRI*N3NXX                                               CPG08530
      ICMAX=IC10+NBASIS*NBASIS                                          CPG08540
      WRITE(3,2) ICMAX,MAXCOR                                           CPG08550
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPG08560
C................B       SM      U..............                        CPG08570
C***  CALL UXMAT(CC(IC8),CC(IC9),CC(IC10),N3NXX)                        CPG08580
C                                                                       CPG08590
C   CALCULATE SCF SECOND DERIVATIVES                                    CPG08600
      WRITE(3,32)                                                       CPG08610
   32 FORMAT(//,2X,' NOW YOU ARE IN WAMAT'/)                            CPG08620
      IC8=IC7+N3N*N3N                                                   CPG08630
      IC9=IC8+NBASIS*NBASIS                                             CPG08640
      IC10=IC9+NBASIS*NBASIS                                            CPG08650
      IC11=IC10+NBASIS*NBASIS                                           CPG08660
      IC12=IC11+NTRI*NTYPEP                                             CPG08670
      IC13=IC12+NBASIS*NBASIS                                           CPG08680
      IC14=IC13+NTRI*N3N*NTYPEP                                         CPG08690
      IC15=IC14+NTRI*N3N*NTYPEP                                         CPG08700
      IC16=IC15+NTRI*N3N*NTYPEP                                         CPG08710
      IA16=IC16+IC16-1                                                  CPG08720
      ICMAX=IC16+MAXBF2                                                 CPG08730
      WRITE(3,2) ICMAX,MAXCOR                                           CPG08740
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPG08750
C................ESO     U       T       EPA      ZETA     WA.......    CPG08760
C***  CALL WAMAT(CC(IC3),CC(IC8),CC(IC9),CC(IC10),CC(IC11),CC(IC12),    CPG08770
C................SA       DP       DQ       DM       LBLI     BUFI..... CPG08780
C*** 1           CC(IC13),CC(IC14),CC(IC15),CC(IC13),IA(IA16),CC(IC16)) CPG08790
C                                                                       CPG08800
      WRITE(3,33)                                                       CPG08810
   33 FORMAT(//,2X,' NOW YOU ARE IN SCF2ND'/)                           CPG08820
      IC8=IC7+N3N*N3N                                                   CPG08830
      IC9=IC8+N3N*N3N                                                   CPG08840
      IC10=IC9+N3N*N3N                                                  CPG08850
      IC11=IC10+NTRI*N3N                                                CPG08860
      IC12=IC11+NBASIS*NBASIS*N3N                                       CPG08870
      IC13=IC12+NBASIS*NBASIS                                           CPG08880
      ICMAX=IC13+NBASIS*NBASIS                                          CPG08890
      WRITE(3,2) ICMAX,MAXCOR                                           CPG08900
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPG08910
C.................E2A     E2M     E2P     SA       WA       EPA......   CPG08920
C***  CALL SCF2ND(CC(IC7),CC(IC8),CC(IC9),CC(IC10),CC(IC11),CC(IC12),   CPG08930
C.................UA.......                                             CPG08940
C*** 1            CC(IC13))                                             CPG08950
C                                                                       CPG08960
      IF(IPOL.EQ.0) GO TO 320                                           CPG08970
      WRITE(3,34)                                                       CPG08980
   34 FORMAT(//,2X,' NOW YOU ARE IN DIPMO'/)                            CPG08990
      IC8=IC7+N3N*N3N                                                   CPG09000
      IC9=IC8+NTRI*3                                                    CPG09010
      ICMAX=IC9+NBASIS*NBASIS*3                                         CPG09020
      WRITE(3,2) ICMAX,MAXCOR                                           CPG09030
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPG09040
C................OCC     DIP     UA......                               CPG09050
C***  CALL DIPMO(CC(IC2),CC(IC8),CC(IC9))                               CPG09060
C                                                                       CPG09070
      WRITE(3,35)                                                       CPG09080
   35 FORMAT(//,2X,' NOW YOU ARE IN POLAR'/)                            CPG09090
      IC8=IC7+N3N*N3N                                                   CPG09100
      IC9=IC8+NTRI*3                                                    CPG09110
      ICMAX=IC9+NBASIS*NBASIS                                           CPG09120
      WRITE(3,2) ICMAX,MAXCOR                                           CPG09130
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPG09140
C................OCC     HF      UF......                               CPG09150
C***  CALL POLAR(CC(IC2),CC(IC8),CC(IC9))                               CPG09160
  320 CONTINUE                                                          CPG09170
C                                                                       CPG09180
      IF(IGMAX.LT.MAXCOR) GO TO 400                                     CPG09190
C                                                                       CPG09200
C   REQUIRED MEMORY EXCEEDS MAXCOR-----STOP !!!                         CPG09210
      WRITE(6,3) IGMAX,MAXCOR                                           CPG09220
      WRITE(3,3) IGMAX,MAXCOR                                           CPG09230
      STOP                                                              CPG09240
C                                                                       CPG09250
  400 CONTINUE                                                          CPG09260
      WRITE(3,4) NADD,MAXVC,IGMAX,MAXCOR                                CPG09270
      WRITE(6,4) NADD,MAXVC,IGMAX,MAXCOR                                CPG09280
      WRITE(3,5)                                                        CPG09290
C                                                                       CPG09300
      RETURN                                                            CPG09310
      END                                                               CPG09320
      SUBROUTINE DERMAT(EAO,E1A,E2A,SS,U,T)                             CPG09330
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG09340
      DIMENSION E1A(N3N),E2A(N3N,N3N),SS(NBATRI)                        CPG09350
      DIMENSION EAO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)         CPG09360
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG09370
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3          CPG09380
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPG09390
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG09400
      COMMON/TAPES/ITAP36,ITAP42,ITAP43,ITAP44                          CPG09410
      COMMON/WORKS/JTAPE1                                               CPG09420
    1 FORMAT(//,2X,' SCF FIRST DERIVATIVE MATRIX'/)                     CPG09430
    2 FORMAT(//,2X,' SCF SECOND DERIVATIVE MATRIX'/)                    CPG09440
    3 FORMAT(//,2X,' SA (AO BASIS) MATRIX, IABC = ',I5/)                CPG09450
    4 FORMAT(//,2X,' SA (MO BASIS) MATRIX, IABC = ',I5/)                CPG09460
    5 FORMAT(//,2X,' HA (AO BASIS) MATRIX, IABC = ',I5/)                CPG09470
    6 FORMAT(//,2X,' HA (MO BASIS) MATRIX, IABC = ',I5/)                CPG09480
    7 FORMAT(//,2X,' TA (AO BASIS) MATRIX, IABC = ',I5,' ITYP = ',I5/)  CPG09490
    8 FORMAT(//,2X,' TA (MO BASIS) MATRIX, IABC = ',I5,' ITYP = ',I5/)  CPG09500
C                                                                       CPG09510
      CALL RFILE(ITAP42)                                                CPG09520
      CALL RFILE(ITAP44)                                                CPG09530
      CALL RFILE(JTAPE1)                                                CPG09540
      NBSET=NTYPES                                                      CPG09550
      NTREAD=NBSET                                                      CPG09560
C                                                                       CPG09570
C   READ IN SCF FIRST DERIVATIVES                                       CPG09580
      CALL SREAD(ITAP42,E1A,N3N*2)                                      CPG09590
      WRITE(6,1)                                                        CPG09600
      CALL MATOUT(E1A,3,NATOMS,3,NATOMS,6)                              CPG09610
C                                                                       CPG09620
C   READ IN SCF SECOND DERIVATIVES                                      CPG09630
      CALL SREAD(ITAP42,E2A,N3N*N3N*2)                                  CPG09640
      WRITE(6,2)                                                        CPG09650
      CALL MATOUT(E2A,N3N,N3N,N3N,N3N,6)                                CPG09660
C                                                                       CPG09670
C   READ IN S FIRST DERIVATIVE INTEGRALS                                CPG09680
C   AND TRANSFORM THEM FROM AO TO MO BASIS                              CPG09690
      DO 101 IABC=1,N3N                                                 CPG09700
      IS=ISA(IABC)                                                      CPG09710
      CALL SREAD(ITAP42,SS,NTRI2)                                       CPG09720
      IF(IPRNT.LE.4) GO TO 201                                          CPG09730
      WRITE(6,3) IABC                                                   CPG09740
      CALL PRINT(SS,NBATRI,NBASIS,6)                                    CPG09750
  201 CONTINUE                                                          CPG09760
      CALL MOCONV(SS,NBATRI,SS,NBATRI,EAO,U,T)                          CPG09770
      CALL RWRIT(ITAP44,SS,NTRI2,IS)                                    CPG09780
      IF(IPRNT.LE.4) GO TO 101                                          CPG09790
      WRITE(6,4) IABC                                                   CPG09800
      CALL PRINT(SS,NBATRI,NBASIS,6)                                    CPG09810
  101 CONTINUE                                                          CPG09820
C                                                                       CPG09830
C   READ IN ONE ELECTRON FIRST DERIVATIVE INTEGRALS                     CPG09840
C   AND TRANSFORM THEM FROM AO TO MO BASIS                              CPG09850
      DO 102 IABC=1,N3N                                                 CPG09860
      IH=IHA(IABC)                                                      CPG09870
      CALL SREAD(ITAP42,SS,NTRI2)                                       CPG09880
      IF(IPRNT.LE.4) GO TO 202                                          CPG09890
      WRITE(6,5) IABC                                                   CPG09900
      CALL PRINT(SS,NBATRI,NBFAO,6)                                     CPG09910
  202 CONTINUE                                                          CPG09920
      CALL MOCONV(SS,NBATRI,SS,NBATRI,EAO,U,T)                          CPG09930
      CALL RWRIT(ITAP44,SS,NTRI2,IH)                                    CPG09940
      IF(IPRNT.LE.4) GO TO 102                                          CPG09950
      WRITE(6,6) IABC                                                   CPG09960
      CALL PRINT(SS,NBATRI,NBASIS,6)                                    CPG09970
  102 CONTINUE                                                          CPG09980
C                                                                       CPG09990
C   READ IN TWO ELECTRON FIRST DERIVATIVE INTEGRALS                     CPG10000
C   AND TRANSFORM THEM FROM AO TO MO BASIS                              CPG10010
      DO 104 ITYP=1,NBSET                                               CPG10020
      DO 103 IABC=1,N3N                                                 CPG10030
      CALL SREAD(ITAP42,SS,NTRI2)                                       CPG10040
      IF(IPRNT.LE.4) GO TO 203                                          CPG10050
      WRITE(6,7) ITYP,IABC                                              CPG10060
      CALL PRINT(SS,NBATRI,NBFAO,6)                                     CPG10070
  203 CALL MOCONV(SS,NBATRI,SS,NBATRI,EAO,U,T)                          CPG10080
      CALL SWRIT(JTAPE1,SS,NTRI2)                                       CPG10090
      IF(IPRNT.LE.4) GO TO 103                                          CPG10100
      WRITE(6,8) ITYP,IABC                                              CPG10110
      CALL PRINT(SS,NBATRI,NBASIS,6)                                    CPG10120
  103 CONTINUE                                                          CPG10130
  104 CONTINUE                                                          CPG10140
C                                                                       CPG10150
      CALL SREW(ITAP42)                                                 CPG10160
      CALL SREW(ITAP44)                                                 CPG10170
      CALL SREW(JTAPE1)                                                 CPG10180
      RETURN                                                            CPG10190
      END                                                               CPG10200
      SUBROUTINE MOCONV(SA,NNA,SM,NNM,EAO,U,T)                          CPG10210
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG10220
      DIMENSION SA(NNA),SM(NNM)                                         CPG10230
      DIMENSION EAO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)         CPG10240
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG10250
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG10260
      DATA ZERO / 0.0D+00 /                                             CPG10270
C                                                                       CPG10280
C   TRANSFORM INTEGRALS FROM AO TO MO BASIS                             CPG10290
      DO 101 II=1,NBFAO                                                 CPG10300
      DO 101 JJ=1,NBFAO                                                 CPG10310
      IIJJ=IOFF(MAX0(II,JJ))+MIN0(II,JJ)                                CPG10320
      SAX=SA(IIJJ)                                                      CPG10330
      T(II,JJ)=SAX                                                      CPG10340
  101 CONTINUE                                                          CPG10350
      DO 103 II=1,NBFAO                                                 CPG10360
      DO 103 J=1,NBASIS                                                 CPG10370
      SUM=ZERO                                                          CPG10380
      DO 102 JJ=1,NBFAO                                                 CPG10390
      SUM=SUM+EAO(JJ,J)*T(II,JJ)                                        CPG10400
  102 CONTINUE                                                          CPG10410
      U(II,J)=SUM                                                       CPG10420
  103 CONTINUE                                                          CPG10430
      DO 105 I=1,NBASIS                                                 CPG10440
      DO 105 J=1,I                                                      CPG10450
      IJ=IOFF(I)+J                                                      CPG10460
      SUM=ZERO                                                          CPG10470
      DO 104 II=1,NBFAO                                                 CPG10480
      SUM=SUM+EAO(II,I)*U(II,J)                                         CPG10490
  104 CONTINUE                                                          CPG10500
      SM(IJ)=SUM                                                        CPG10510
  105 CONTINUE                                                          CPG10520
C                                                                       CPG10530
      RETURN                                                            CPG10540
      END                                                               CPG10550
      SUBROUTINE MOCONS(SA,NNA,SM,NNM,ESO,U,T)                          CPG10560
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG10570
      DIMENSION SA(NNA),SM(NNM)                                         CPG10580
      DIMENSION ESO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)         CPG10590
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG10600
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG10610
      DATA ZERO / 0.0D+00 /                                             CPG10620
C                                                                       CPG10630
C   TRANSFORM INTEGRALS FROM SO TO MO BASIS                             CPG10640
      DO 101 II=1,NBASIS                                                CPG10650
      DO 101 JJ=1,II                                                    CPG10660
      IIJJ=IOFF(II)+JJ                                                  CPG10670
      SAX=SA(IIJJ)                                                      CPG10680
      T(II,JJ)=SAX                                                      CPG10690
      IF(II.EQ.JJ) GO TO 101                                            CPG10700
      T(JJ,II)=SAX                                                      CPG10710
  101 CONTINUE                                                          CPG10720
      DO 103 II=1,NBASIS                                                CPG10730
      DO 103 J=1,NBASIS                                                 CPG10740
      SUM=ZERO                                                          CPG10750
      DO 102 JJ=1,NBASIS                                                CPG10760
      SUM=SUM+ESO(JJ,J)*T(II,JJ)                                        CPG10770
  102 CONTINUE                                                          CPG10780
      U(II,J)=SUM                                                       CPG10790
  103 CONTINUE                                                          CPG10800
      DO 105 I=1,NBASIS                                                 CPG10810
      DO 105 J=1,I                                                      CPG10820
      IJ=IOFF(I)+J                                                      CPG10830
      SUM=ZERO                                                          CPG10840
      DO 104 II=1,NBASIS                                                CPG10850
      SUM=SUM+ESO(II,I)*U(II,J)                                         CPG10860
  104 CONTINUE                                                          CPG10870
      SM(IJ)=SUM                                                        CPG10880
  105 CONTINUE                                                          CPG10890
C                                                                       CPG10900
      RETURN                                                            CPG10910
      END                                                               CPG10920
      SUBROUTINE E0CAL(OCC,ESO,HSO,HMO,TM,U,T,DP,DQ,DM,LBLI,BUFI)       CPG10930
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG10940
      DIMENSION OCC(NBASIS),HSO(NTRI),HMO(NTRI),TM(NTRI)                CPG10950
      DIMENSION ESO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)         CPG10960
      DIMENSION DP(NTRI,NTYPES),DQ(NTRI,NTYPES),DM(NTRI,NTYPES)         CPG10970
      DIMENSION LBLI(MAXBF4),BUFI(MAXBF2)                               CPG10980
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG10990
      COMMON/ENRGY/ESCF,ENUC                                            CPG11000
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3          CPG11010
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPG11020
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG11030
      COMMON/START/FOCC(10),NSORB(10),NSTR(10),NEND(10),MOTYP(256)      CPG11040
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPG11050
      DATA ZERO / 0.0D+00 /                                             CPG11060
    1 FORMAT(//,2X,' H (SO BASIS) MATRIX'/)                             CPG11070
    2 FORMAT(//,2X,' H (MO BASIS) MATRIX'/)                             CPG11080
    3 FORMAT(//,2X,' SCF ENERGIES'/                                     CPG11090
     * 2X,' ESCF       = ',F20.10/                                      CPG11100
     * 2X,' ENUC       = ',F20.10/                                      CPG11110
     * 2X,' ELEC-ONE   = ',F20.10/                                      CPG11120
     * 2X,' ELEC-TWO   = ',F20.10/                                      CPG11130
     * 2X,' ELEC-TOTAL = ',F20.10/                                      CPG11140
     * 2X,' ETOT       = ',F20.10)                                      CPG11150
C                                                                       CPG11160
C   CALCULATE SCF ENERGY FOR A TEST                                     CPG11170
C   SEE EQ.(1)                                                          CPG11180
C                                                                       CPG11190
C   READ IN ONE ELCTRON INTEGRALS (SO BASIS)                            CPG11200
      CALL MREAD(HSO,15)                                                CPG11210
      CALL MOCONS(HSO,NTRI,HMO,NTRI,ESO,U,T)                            CPG11220
      IF(IPRNT.LE.3) GO TO 201                                          CPG11230
      WRITE(6,1)                                                        CPG11240
      CALL PRINT(HSO,NTRI,NBASIS,6)                                     CPG11250
      WRITE(6,2)                                                        CPG11260
      CALL PRINT(HMO,NTRI,NBASIS,6)                                     CPG11270
C                                                                       CPG11280
C   CALCULATE ONE-ELCTRON ENERGY                                        CPG11290
  201 CONTINUE                                                          CPG11300
      ELEC1=ZERO                                                        CPG11310
      DO 101 I=1,JOCC                                                   CPG11320
      II=IOFF(I+1)                                                      CPG11330
      ELEC1=ELEC1+HMO(II)*OCC(I)                                        CPG11340
  101 CONTINUE                                                          CPG11350
C                                                                       CPG11360
C   FORM DENSITY MATRIX IN SO BASIS                                     CPG11370
      CALL PQAMAT(DP,DQ,ESO,NTYPES)                                     CPG11380
C                                                                       CPG11390
C   FORM HALF-TRANSFORMED TWO ELECTRON MATRIX                           CPG11400
      CALL PQMAT(DP,DQ,DM,NTYPES,LBLI,BUFI)                             CPG11410
C                                                                       CPG11420
C   CALCULATE TWO-ELECTRON ENERGY                                       CPG11430
      ELEC2=ZERO                                                        CPG11440
      DO 103 ITYP=1,NTYPES                                              CPG11450
      CALL MOCONS(DM(1,ITYP),NTRI,TM,NTRI,ESO,U,T)                      CPG11460
      DO 102 I=NSTR(ITYP),NEND(ITYP)                                    CPG11470
      II=IOFF(I+1)                                                      CPG11480
      ELEC2=ELEC2+TM(II)                                                CPG11490
  102 CONTINUE                                                          CPG11500
  103 CONTINUE                                                          CPG11510
C                                                                       CPG11520
      ELECT=ELEC1+ELEC2                                                 CPG11530
      ETOT=ENUC+ELECT                                                   CPG11540
      WRITE(6,3) ESCF,ENUC,ELEC1,ELEC2,ELECT,ETOT                       CPG11550
C                                                                       CPG11560
      RETURN                                                            CPG11570
      END                                                               CPG11580
      SUBROUTINE PQAMAT(DP,DQ,ESO,NABC)                                 CPG11590
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG11600
      DIMENSION DP(NTRI,NABC),DQ(NTRI,NABC),ESO(NBFAO,NBASIS)           CPG11610
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG11620
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3          CPG11630
      COMMON/GVBAB/A(10,10),B(10,10)                                    CPG11640
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG11650
      COMMON/START/FOCC(10),NSORB(10),NSTR(10),NEND(10),MOTYP(256)      CPG11660
      DATA ZERO / 0.0D+00 /                                             CPG11670
    1 FORMAT(//,2X,' DP MATRIX'/)                                       CPG11680
    2 FORMAT(//,2X,' DQ MATRIX'/)                                       CPG11690
C                                                                       CPG11700
      IJ=0                                                              CPG11710
      DO 105 I=1,NBASIS                                                 CPG11720
      DO 105 J=1,I                                                      CPG11730
      IJ=IJ+1                                                           CPG11740
      DO 104 ITYP=1,NTYPES                                              CPG11750
      SUMP=ZERO                                                         CPG11760
      SUMQ=ZERO                                                         CPG11770
      DO 103 JTYP=1,NTYPES                                              CPG11780
      DO 102 K=NSTR(JTYP),NEND(JTYP)                                    CPG11790
      FAC=ESO(I,K)*ESO(J,K)                                             CPG11800
      SUMP=SUMP+A(ITYP,JTYP)*FAC                                        CPG11810
      SUMQ=SUMQ+B(ITYP,JTYP)*FAC                                        CPG11820
  102 CONTINUE                                                          CPG11830
  103 CONTINUE                                                          CPG11840
      DP(IJ,ITYP)=SUMP                                                  CPG11850
      DQ(IJ,ITYP)=SUMQ                                                  CPG11860
  104 CONTINUE                                                          CPG11870
  105 CONTINUE                                                          CPG11880
C                                                                       CPG11890
      IF(IPRNT.LE.4) GO TO 201                                          CPG11900
      WRITE(6,1)                                                        CPG11910
      CALL MATOUT(DP,NTRI,NTYPES,NTRI,NTYPES,6)                         CPG11920
      WRITE(6,2)                                                        CPG11930
      CALL MATOUT(DQ,NTRI,NTYPES,NTRI,NTYPES,6)                         CPG11940
C                                                                       CPG11950
  201 CONTINUE                                                          CPG11960
      RETURN                                                            CPG11970
      END                                                               CPG11980
      SUBROUTINE PQMAT(DP,DQ,DM,NABC,LBLI,BUFI)                         CPG11990
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG12000
      DIMENSION DP(NTRI,NABC),DQ(NTRI,NABC),DM(NTRI,NABC)               CPG12010
      DIMENSION LBLI(MAXBF4),BUFI(MAXBF2)                               CPG12020
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG12030
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPG12040
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG12050
      COMMON/TAPES/ITAP36,ITAP42,ITAP43,ITAP44                          CPG12060
      DATA ZERO,TWO / 0.0D+00 , 2.0D+00 /                               CPG12070
    1 FORMAT(2X,6I5,F15.8)                                              CPG12080
    2 FORMAT(//,2X,' DM MATRIX'/)                                       CPG12090
C                                                                       CPG12100
C   DP :  COULOMB DENSITY TYPE SO MATRIX                                CPG12110
C   DQ :  EXCHANGE DENSITY TYPE SO MATRIX                               CPG12120
C   DM :  HALF TRANSFORMED SO MATRIX                                    CPG12130
C                                                                       CPG12140
C   FORM DM MATRICES                                                    CPG12150
      DO 101 IABC=1,NABC                                                CPG12160
      DO 101 I=1,NTRI                                                   CPG12170
      DM(I,IABC)=ZERO                                                   CPG12180
  101 CONTINUE                                                          CPG12190
C                                                                       CPG12200
      CALL RFILE(ITAP36)                                                CPG12210
  500 CALL SREAD(ITAP36,LBLI,MAXBF4)                                    CPG12220
      DO 160 I=1,MAXBUF                                                 CPG12230
      IBA=LBLI(I+I-1)                                                   CPG12240
      IBB=LBLI(I+I)                                                     CPG12250
      IF(IBA.EQ.0) GO TO 501                                            CPG12260
      VALU=BUFI(MAXBUF+I)                                               CPG12270
      CALL UNPACK(II,JJ,KK,LL,IX,IBA,IBB)                               CPG12280
C     IF(IPRNT.LE.9) GO TO 301                                          CPG12290
C     WRITE(6,1) I,II,JJ,KK,LL,IX,VALU                                  CPG12300
  301 GO TO (201,202,203,204,205,206,207,208,209,210,211,212,212,212),IXCPG12310
C                                                                       CPG12320
C (II/II)   1=(11/11)                                                   CPG12330
C                                                                       CPG12340
  201 IIII=IOFF(II+1)                                                   CPG12350
      DO 112 IABC=1,NABC                                                CPG12360
      DM(IIII,IABC)=DM(IIII,IABC)+(DP(IIII,IABC)+DQ(IIII,IABC))*VALU    CPG12370
  112 CONTINUE                                                          CPG12380
      GO TO 160                                                         CPG12390
C                                                                       CPG12400
C  (II/KK)  2=(22/11)                                                   CPG12410
C                                                                       CPG12420
  202 IIII=IOFF(II+1)                                                   CPG12430
      KKKK=IOFF(KK+1)                                                   CPG12440
      IIKK=IOFF(II)+KK                                                  CPG12450
      DO 116 IABC=1,NABC                                                CPG12460
      DM(IIII,IABC)=DM(IIII,IABC)+DP(KKKK,IABC)*VALU                    CPG12470
      DM(KKKK,IABC)=DM(KKKK,IABC)+DP(IIII,IABC)*VALU                    CPG12480
      DM(IIKK,IABC)=DM(IIKK,IABC)+DQ(IIKK,IABC)*VALU                    CPG12490
  116 CONTINUE                                                          CPG12500
      GO TO 160                                                         CPG12510
C                                                                       CPG12520
C  (II/IL)  3=(22/21)                                                   CPG12530
C                                                                       CPG12540
  203 IIII=IOFF(II+1)                                                   CPG12550
      IILL=IOFF(II)+LL                                                  CPG12560
      VAL2=VALU+VALU                                                    CPG12570
      DO 120 IABC=1,NABC                                                CPG12580
      DM(IIII,IABC)=DM(IIII,IABC)+(DP(IILL,IABC)+DQ(IILL,IABC))*VAL2    CPG12590
      DM(IILL,IABC)=DM(IILL,IABC)+(DP(IIII,IABC)+DQ(IIII,IABC))*VALU    CPG12600
  120 CONTINUE                                                          CPG12610
      GO TO 160                                                         CPG12620
C                                                                       CPG12630
C  (IJ/JJ)  4=(21/11)                                                   CPG12640
C                                                                       CPG12650
  204 IIJJ=IOFF(II)+JJ                                                  CPG12660
      JJJJ=IOFF(JJ+1)                                                   CPG12670
      VAL2=VALU+VALU                                                    CPG12680
      DO 124 IABC=1,NABC                                                CPG12690
      DM(IIJJ,IABC)=DM(IIJJ,IABC)+(DP(JJJJ,IABC)+DQ(JJJJ,IABC))*VALU    CPG12700
      DM(JJJJ,IABC)=DM(JJJJ,IABC)+(DP(IIJJ,IABC)+DQ(IIJJ,IABC))*VAL2    CPG12710
  124 CONTINUE                                                          CPG12720
      GO TO 160                                                         CPG12730
C                                                                       CPG12740
C  (IJ/IJ)  5=(21/21)                                                   CPG12750
C                                                                       CPG12760
  205 IIJJ=IOFF(II)+JJ                                                  CPG12770
      IIII=IOFF(II+1)                                                   CPG12780
      JJJJ=IOFF(JJ+1)                                                   CPG12790
      DO 128 IABC=1,NABC                                                CPG12800
      DM(IIJJ,IABC)=DM(IIJJ,IABC)+(DP(IIJJ,IABC)*TWO+DQ(IIJJ,IABC))*VALUCPG12810
      DM(IIII,IABC)=DM(IIII,IABC)+DQ(JJJJ,IABC)*VALU                    CPG12820
      DM(JJJJ,IABC)=DM(JJJJ,IABC)+DQ(IIII,IABC)*VALU                    CPG12830
  128 CONTINUE                                                          CPG12840
      GO TO 160                                                         CPG12850
C                                                                       CPG12860
C  (II/KL)  6=(33/21)                                                   CPG12870
C                                                                       CPG12880
  206 IIII=IOFF(II+1)                                                   CPG12890
      KKLL=IOFF(KK)+LL                                                  CPG12900
      IIKK=IOFF(II)+KK                                                  CPG12910
      IILL=IOFF(II)+LL                                                  CPG12920
      VAL2=VALU+VALU                                                    CPG12930
      DO 132 IABC=1,NABC                                                CPG12940
      DM(IIII,IABC)=DM(IIII,IABC)+DP(KKLL,IABC)*VAL2                    CPG12950
      DM(KKLL,IABC)=DM(KKLL,IABC)+DP(IIII,IABC)*VALU                    CPG12960
      DM(IIKK,IABC)=DM(IIKK,IABC)+DQ(IILL,IABC)*VALU                    CPG12970
      DM(IILL,IABC)=DM(IILL,IABC)+DQ(IIKK,IABC)*VALU                    CPG12980
  132 CONTINUE                                                          CPG12990
      GO TO 160                                                         CPG13000
C                                                                       CPG13010
C  (IJ/KK)  7=(32/11)                                                   CPG13020
C                                                                       CPG13030
  207 IIJJ=IOFF(II)+JJ                                                  CPG13040
      KKKK=IOFF(KK+1)                                                   CPG13050
      IIKK=IOFF(II)+KK                                                  CPG13060
      JJKK=IOFF(JJ)+KK                                                  CPG13070
      VAL2=VALU+VALU                                                    CPG13080
      DO 136 IABC=1,NABC                                                CPG13090
      DM(IIJJ,IABC)=DM(IIJJ,IABC)+DP(KKKK,IABC)*VALU                    CPG13100
      DM(KKKK,IABC)=DM(KKKK,IABC)+DP(IIJJ,IABC)*VAL2                    CPG13110
      DM(IIKK,IABC)=DM(IIKK,IABC)+DQ(JJKK,IABC)*VALU                    CPG13120
      DM(JJKK,IABC)=DM(JJKK,IABC)+DQ(IIKK,IABC)*VALU                    CPG13130
  136 CONTINUE                                                          CPG13140
      GO TO 160                                                         CPG13150
C                                                                       CPG13160
C  (IJ/KK)  8=(31/22)                                                   CPG13170
C                                                                       CPG13180
  208 IIJJ=IOFF(II)+JJ                                                  CPG13190
      KKKK=IOFF(KK+1)                                                   CPG13200
      IIKK=IOFF(II)+KK                                                  CPG13210
      JJKK=JJ+IOFF(KK)                                                  CPG13220
      VAL2=VALU+VALU                                                    CPG13230
      DO 140 IABC=1,NABC                                                CPG13240
      DM(IIJJ,IABC)=DM(IIJJ,IABC)+DP(KKKK,IABC)*VALU                    CPG13250
      DM(KKKK,IABC)=DM(KKKK,IABC)+DP(IIJJ,IABC)*VAL2                    CPG13260
      DM(IIKK,IABC)=DM(IIKK,IABC)+DQ(JJKK,IABC)*VALU                    CPG13270
      DM(JJKK,IABC)=DM(JJKK,IABC)+DQ(IIKK,IABC)*VALU                    CPG13280
  140 CONTINUE                                                          CPG13290
      GO TO 160                                                         CPG13300
C                                                                       CPG13310
C (IJ/IL)  9=(32/31)                                                    CPG13320
C                                                                       CPG13330
  209 IIJJ=IOFF(II)+JJ                                                  CPG13340
      IILL=IOFF(II)+LL                                                  CPG13350
      IIII=IOFF(II+1)                                                   CPG13360
      JJLL=IOFF(JJ)+LL                                                  CPG13370
      VAL2=VALU+VALU                                                    CPG13380
      DO 144 IABC=1,NABC                                                CPG13390
      DM(IIJJ,IABC)=DM(IIJJ,IABC)+(DP(IILL,IABC)*TWO+DQ(IILL,IABC))*VALUCPG13400
      DM(IILL,IABC)=DM(IILL,IABC)+(DP(IIJJ,IABC)*TWO+DQ(IIJJ,IABC))*VALUCPG13410
      DM(IIII,IABC)=DM(IIII,IABC)+DQ(JJLL,IABC)*VAL2                    CPG13420
      DM(JJLL,IABC)=DM(JJLL,IABC)+DQ(IIII,IABC)*VALU                    CPG13430
  144 CONTINUE                                                          CPG13440
      GO TO 160                                                         CPG13450
C                                                                       CPG13460
C  (IJ/JL)  10=(32/21)                                                  CPG13470
C                                                                       CPG13480
  210 IIJJ=IOFF(II)+JJ                                                  CPG13490
      JJLL=IOFF(JJ)+LL                                                  CPG13500
      IILL=IOFF(II)+LL                                                  CPG13510
      JJJJ=IOFF(JJ+1)                                                   CPG13520
      VAL2=VALU+VALU                                                    CPG13530
      DO 148 IABC=1,NABC                                                CPG13540
      DM(IIJJ,IABC)=DM(IIJJ,IABC)+(DP(JJLL,IABC)*TWO+DQ(JJLL,IABC))*VALUCPG13550
      DM(JJLL,IABC)=DM(JJLL,IABC)+(DP(IIJJ,IABC)*TWO+DQ(IIJJ,IABC))*VALUCPG13560
      DM(IILL,IABC)=DM(IILL,IABC)+DQ(JJJJ,IABC)*VALU                    CPG13570
      DM(JJJJ,IABC)=DM(JJJJ,IABC)+DQ(IILL,IABC)*VAL2                    CPG13580
  148 CONTINUE                                                          CPG13590
      GO TO 160                                                         CPG13600
C                                                                       CPG13610
C  (IJ/KJ)  11=(31/21)                                                  CPG13620
C                                                                       CPG13630
  211 IIJJ=IOFF(II)+JJ                                                  CPG13640
      JJKK=JJ+IOFF(KK)                                                  CPG13650
      IIKK=IOFF(II)+KK                                                  CPG13660
      JJJJ=IOFF(JJ+1)                                                   CPG13670
      VAL2=VALU+VALU                                                    CPG13680
      DO 152 IABC=1,NABC                                                CPG13690
      DM(IIJJ,IABC)=DM(IIJJ,IABC)+(DP(JJKK,IABC)*TWO+DQ(JJKK,IABC))*VALUCPG13700
      DM(JJKK,IABC)=DM(JJKK,IABC)+(DP(IIJJ,IABC)*TWO+DQ(IIJJ,IABC))*VALUCPG13710
      DM(IIKK,IABC)=DM(IIKK,IABC)+DQ(JJJJ,IABC)*VALU                    CPG13720
      DM(JJJJ,IABC)=DM(JJJJ,IABC)+DQ(IIKK,IABC)*VAL2                    CPG13730
  152 CONTINUE                                                          CPG13740
      GO TO 160                                                         CPG13750
C                                                                       CPG13760
C  (IJ/KL)  12=(43/21)                                                  CPG13770
C  (IJ/KL)  13=(42/31)                                                  CPG13780
C  (IJ/KL)  14=(41/32)                                                  CPG13790
C                                                                       CPG13800
  212 IIJJ=IOFF(II)+JJ                                                  CPG13810
      KKLL=IOFF(KK)+LL                                                  CPG13820
      IIKK=IOFF(II)+KK                                                  CPG13830
      IILL=IOFF(II)+LL                                                  CPG13840
      JJLL=IOFF(MAX0(JJ,LL))+MIN0(JJ,LL)                                CPG13850
      JJKK=IOFF(MAX0(JJ,KK))+MIN0(JJ,KK)                                CPG13860
      VAL2=VALU+VALU                                                    CPG13870
      DO 156 IABC=1,NABC                                                CPG13880
      DM(IIJJ,IABC)=DM(IIJJ,IABC)+DP(KKLL,IABC)*VAL2                    CPG13890
      DM(KKLL,IABC)=DM(KKLL,IABC)+DP(IIJJ,IABC)*VAL2                    CPG13900
      DM(IIKK,IABC)=DM(IIKK,IABC)+DQ(JJLL,IABC)*VALU                    CPG13910
      DM(JJLL,IABC)=DM(JJLL,IABC)+DQ(IIKK,IABC)*VALU                    CPG13920
      DM(IILL,IABC)=DM(IILL,IABC)+DQ(JJKK,IABC)*VALU                    CPG13930
      DM(JJKK,IABC)=DM(JJKK,IABC)+DQ(IILL,IABC)*VALU                    CPG13940
  156 CONTINUE                                                          CPG13950
C                                                                       CPG13960
  160 CONTINUE                                                          CPG13970
      GO TO 500                                                         CPG13980
C                                                                       CPG13990
  501 CALL RCLOSE(ITAP36,3)                                             CPG14000
C                                                                       CPG14010
      IF (IPRNT.LE.7) GO TO 502                                         CPG14020
      WRITE(6,2)                                                        CPG14030
      CALL MATOUT(DM,NTRI,NABC,NTRI,NABC,6)                             CPG14040
C                                                                       CPG14050
  502 CONTINUE                                                          CPG14060
      RETURN                                                            CPG14070
      END                                                               CPG14080
      SUBROUTINE EADER(OCC,D1O,D1T,D1W,E1T,SM,HM,EPS,TM)                CPG14090
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG14100
      DIMENSION OCC(NBASIS),SM(NTRI),HM(NTRI),EPS(NTRI)                 CPG14110
      DIMENSION TM(NTRI,N3N,NTYPES)                                     CPG14120
      DIMENSION D1O(N3N),D1T(N3N),D1W(N3N),E1T(N3N)                     CPG14130
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG14140
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3          CPG14150
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPG14160
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG14170
      COMMON/START/FOCC(10),NSORB(10),NSTR(10),NEND(10),MOTYP(256)      CPG14180
      COMMON/TAPES/ITAP36,ITAP42,ITAP43,ITAP44                          CPG14190
      COMMON/WORKS/JTAPE1                                               CPG14200
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPG14210
      DATA ZERO,TWO / 0.0D+00 , 2.0D+00 /                               CPG14220
    1 FORMAT(//,2X,' LAGRANGIAN MATRIX'/)                               CPG14230
    2 FORMAT(//,2X,' NUCLEAR REPULSION FIRST DERIVATIVE MATRIX'/)       CPG14240
    3 FORMAT(//,2X,' SCF ONE-ELECTRON FIRST DERIVATIVE MATRIX'/)        CPG14250
    4 FORMAT(//,2X,' SCF TWO-ELECTRON FIRST DERIVATIVE MATRIX'/)        CPG14260
    5 FORMAT(//,2X,' SCF OVERLAP CONTRIBUTION MATRIX'/)                 CPG14270
    6 FORMAT(//,2X,' TOTAL SCF FIRST DERIVATIVE MATRIX'/)               CPG14280
C                                                                       CPG14290
C   CALCULATE SCF FIRST DERIVATIVES FOR A TEST                          CPG14300
C   SEE EQ.(2)                                                          CPG14310
C                                                                       CPG14320
      CALL MREAD(EPS,25)                                                CPG14330
      IF(IPRNT.LE.3) GO TO 201                                          CPG14340
      WRITE(6,1)                                                        CPG14350
      CALL PRINT(EPS,NTRI,NBASIS,6)                                     CPG14360
C                                                                       CPG14370
C   CALL BACK DERIVATIVE TM MATRICES                                    CPG14380
  201 CONTINUE                                                          CPG14390
      CALL SREW(JTAPE1)                                                 CPG14400
      DO 102 ITYP=1,NTYPES                                              CPG14410
      DO 101 IABC=1,N3N                                                 CPG14420
      CALL SREAD(JTAPE1,TM(1,IABC,ITYP),NTRI2)                          CPG14430
  101 CONTINUE                                                          CPG14440
  102 CONTINUE                                                          CPG14450
C                                                                       CPG14460
      CALL SREW(ITAP44)                                                 CPG14470
      DO 110 IABC=1,N3N                                                 CPG14480
      IS=ISA(IABC)                                                      CPG14490
      IH=IHA(IABC)                                                      CPG14500
      CALL RREAD(ITAP44,SM,NTRI2,IS)                                    CPG14510
      CALL RREAD(ITAP44,HM,NTRI2,IH)                                    CPG14520
C  ONE-ELECTRON CONTRIBUTION                                            CPG14530
      VALU=ZERO                                                         CPG14540
      DO 106 I=1,JOCC                                                   CPG14550
      II=IOFF(I+1)                                                      CPG14560
      VALU=VALU+HM(II)*OCC(I)                                           CPG14570
  106 CONTINUE                                                          CPG14580
      D1O(IABC)=VALU                                                    CPG14590
C  TWO-ELECTRON CONTRIBUTION                                            CPG14600
      VALU=ZERO                                                         CPG14610
      DO 108 ITYP=1,NTYPES                                              CPG14620
      DO 107 I=NSTR(ITYP),NEND(ITYP)                                    CPG14630
      II=IOFF(I+1)                                                      CPG14640
      VALU=VALU+TM(II,IABC,ITYP)                                        CPG14650
  107 CONTINUE                                                          CPG14660
      D1T(IABC)=VALU                                                    CPG14670
  108 CONTINUE                                                          CPG14680
C   OVERLAP CONTRIBUTION                                                CPG14690
      VALU=ZERO                                                         CPG14700
      DO 109 I=1,JOCC                                                   CPG14710
      DO 109 J=1,JOCC                                                   CPG14720
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPG14730
      VALU=VALU-SM(IJ)*EPS(IJ)                                          CPG14740
  109 CONTINUE                                                          CPG14750
      D1W(IABC)=VALU*TWO                                                CPG14760
C                                                                       CPG14770
C   SUM UP ALL CONTRIBUTIONS                                            CPG14780
      E1T(IABC)=D1O(IABC)+D1T(IABC)+D1W(IABC)                           CPG14790
C                                                                       CPG14800
  110 CONTINUE                                                          CPG14810
C                                                                       CPG14820
C     WRITE(6,2)                                                        CPG14830
C     CALL MATOUT(E1N,3,NATOMS,3,NATOMS,6)                              CPG14840
      WRITE(6,3)                                                        CPG14850
      CALL MATOUT(D1O,3,NATOMS,3,NATOMS,6)                              CPG14860
      WRITE(6,4)                                                        CPG14870
      CALL MATOUT(D1T,3,NATOMS,3,NATOMS,6)                              CPG14880
      WRITE(6,5)                                                        CPG14890
      CALL MATOUT(D1W,3,NATOMS,3,NATOMS,6)                              CPG14900
      WRITE(6,6)                                                        CPG14910
      CALL MATOUT(E1T,3,NATOMS,3,NATOMS,6)                              CPG14920
C                                                                       CPG14930
      CALL SREW(ITAP44)                                                 CPG14940
      CALL SREW(JTAPE1)                                                 CPG14950
      RETURN                                                            CPG14960
      END                                                               CPG14970
      SUBROUTINE FAMAT(HH,TT,EPA,FM)                                    CPG14980
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG14990
      DIMENSION HH(NTRI),TT(NTRI),FM(NTRI,N3N,NTYPEP)                   CPG15000
      DIMENSION EPA(NBASIS,NBASIS)                                      CPG15010
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG15020
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3          CPG15030
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPG15040
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG15050
      COMMON/START/FOCC(10),NSORB(10),NSTR(10),NEND(10),MOTYP(256)      CPG15060
      COMMON/TAPES/ITAP36,ITAP42,ITAP43,ITAP44                          CPG15070
      COMMON/WORKS/JTAPE1                                               CPG15080
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPG15090
      DATA ZERO,HALF / 0.0D+00 , 0.5D+00 /                              CPG15100
    1 FORMAT(//,2X,' EPA MATRIX, IABC = ',I5/)                          CPG15110
C                                                                       CPG15120
C   CALL BACK HM AND TM MATRICES AND FORM FA MATRICES                   CPG15130
C   SEE EQ.(7)                                                          CPG15140
C                                                                       CPG15150
      CALL SREW(ITAP44)                                                 CPG15160
      CALL SREW(JTAPE1)                                                 CPG15170
      DO 103 IABC=1,N3N                                                 CPG15180
      IH=IHA(IABC)                                                      CPG15190
      CALL RREAD(ITAP44,HH,NTRI2,IH)                                    CPG15200
      DO 102 ITYP=1,NTREAD                                              CPG15210
      FAC=FOCC(ITYP)*HALF                                               CPG15220
      DO 101 I=1,NTRI                                                   CPG15230
  101 FM(I,IABC,ITYP)=HH(I)*FAC                                         CPG15240
  102 CONTINUE                                                          CPG15250
  103 CONTINUE                                                          CPG15260
C                                                                       CPG15270
      DO 106 ITYP=1,NTREAD                                              CPG15280
      DO 105 IABC=1,N3N                                                 CPG15290
      CALL SREAD(JTAPE1,TT,NTRI2)                                       CPG15300
      DO 104 I=1,NTRI                                                   CPG15310
      FM(I,IABC,ITYP)=FM(I,IABC,ITYP)+TT(I)                             CPG15320
  104 CONTINUE                                                          CPG15330
  105 CONTINUE                                                          CPG15340
  106 CONTINUE                                                          CPG15350
C                                                                       CPG15360
C   STORE DERIVATIVE LAGRANGIAN MATRICES                                CPG15370
      DO 120 IABC=1,N3N                                                 CPG15380
      IE=IEA(IABC)                                                      CPG15390
      DO 110 I=1,NBASIS                                                 CPG15400
      DO 110 J=1,NBASIS                                                 CPG15410
  110 EPA(I,J)=ZERO                                                     CPG15420
C                                                                       CPG15430
      DO 112 ITYP=1,NTYPES                                              CPG15440
      DO 111 I=NSTR(ITYP),NEND(ITYP)                                    CPG15450
      DO 111 J=1,NBASIS                                                 CPG15460
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPG15470
      EPA(I,J)=FM(IJ,IABC,ITYP)                                         CPG15480
  111 CONTINUE                                                          CPG15490
  112 CONTINUE                                                          CPG15500
C                                                                       CPG15510
      CALL RWRIT(ITAP44,EPA,NBASQ,IE)                                   CPG15520
      IF(IPRNT.LE.4) GO TO 120                                          CPG15530
      WRITE(6,1) IABC                                                   CPG15540
      CALL MATOUT(EPA,NBASIS,NBASIS,NBASIS,NBASIS,6)                    CPG15550
  120 CONTINUE                                                          CPG15560
C                                                                       CPG15570
C   STORE DERIVATIVE AVERAGED FOCK                                      CPG15580
CC*   IF(ICORE.NE.0) GO TO 201                                          CPG15590
CC*   DO 121 I=1,NTRI                                                   CPG15600
CC121 TT(I)=ZERO                                                        CPG15610
CC*   DO 122 I=1,N3N                                                    CPG15620
CC*   IT=IFA(IABC)                                                      CPG15630
CC*   CALL RWRIT(ITAP44,TT,NTRI2,IT)                                    CPG15640
CC122 CONTINUE                                                          CPG15650
CC*   GO TO 205                                                         CPG15660
CC201 CONTINUE                                                          CPG15670
      DO 125 IABC=1,N3N                                                 CPG15680
      IT=IFA(IABC)                                                      CPG15690
      CALL RWRIT(ITAP44,FM(1,IABC,NTREAD),NTRI2,IT)                     CPG15700
  125 CONTINUE                                                          CPG15710
  205 CONTINUE                                                          CPG15720
C                                                                       CPG15730
      CALL SREW(ITAP44)                                                 CPG15740
      CALL SREW(JTAPE1)                                                 CPG15750
      RETURN                                                            CPG15760
      END                                                               CPG15770
      SUBROUTINE BAIND(ESO,NIJ,U,T,ZETA,B0,EPA,DP,DQ,DM,SM,LBLI,BUFI)   CPG15780
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG15790
      DIMENSION ESO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)         CPG15800
      DIMENSION NIJ(NTRI,2),ZETA(NTRI,NTYPEP)                           CPG15810
      DIMENSION B0(NTRI),EPA(NBASIS,NBASIS)                             CPG15820
      DIMENSION DP(NTRI,N3N,NTYPEP),DQ(NTRI,N3N,NTYPEP)                 CPG15830
      DIMENSION DM(NTRI,N3N,NTYPEP),SM(NTRI,N3N*NTYPEP)                 CPG15840
      DIMENSION LBLI(MAXBF4),BUFI(MAXBF2)                               CPG15850
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG15860
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3          CPG15870
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPG15880
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPG15890
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG15900
      COMMON/START/FOCC(10),NSORB(10),NSTR(10),NEND(10),MOTYP(256)      CPG15910
      COMMON/TAPES/ITAP36,ITAP42,ITAP43,ITAP44                          CPG15920
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPG15930
      COMMON/CI102/NIND,NDEP                                            CPG15940
      DATA ZERO / 0.0D+00 /                                             CPG15950
    1 FORMAT(//,2X,' ZETA MATRIX, ITYP = ',I5/)                         CPG15960
    2 FORMAT(//,2X,' FM MATRIX, IABC = ',I5/)                           CPG15970
    3 FORMAT(2X,4I5,3F15.8)                                             CPG15980
    4 FORMAT(//,2X,' B0 MATRIX IN BAIND'/)                              CPG15990
    5 FORMAT(//,2X,' DP MATRIX IN BAIND'/)                              CPG16000
C                                                                       CPG16010
C   FORM B0 MATRIX                                                      CPG16020
C   SEE EQ.(13)                                                         CPG16030
C                                                                       CPG16040
      DO 101 ITYP=1,NTYPES                                              CPG16050
      CALL MREAD(ZETA(1,ITYP),26+ITYP)                                  CPG16060
      IF(IPRNT.LE.2) GO TO 101                                          CPG16070
      WRITE(6,1) ITYP                                                   CPG16080
      CALL PRINT(ZETA(1,ITYP),NTRI,NBASIS,6)                            CPG16090
  101 CONTINUE                                                          CPG16100
      DO 102 I=1,NTRI                                                   CPG16110
  102 ZETA(I,NTYPEP)=ZERO                                               CPG16120
c
cyy   WRITE(6,*) 'II,I,J,K,L,SVALU,FAC2,VALU1,VALU2'
      CALL SREW(ITAP44)                                                 CPG16140
      DO 120 IABC=1,N3N                                                 CPG16150
      IS=ISA(IABC)                                                      CPG16160
      IE=IEA(IABC)                                                      CPG16170
      IB=IBA(IABC)                                                      CPG16180
      CALL RREAD(ITAP44,SM(1,IABC),NTRI2,IS)                            CPG16190
      CALL RREAD(ITAP44,EPA,NBASQ,IE)                                   CPG16200
      DO 103 I=1,NTRI                                                   CPG16210
  103 B0(I)=ZERO                                                        CPG16220
C                                                                       CPG16230
C   ELEMENTS FOR INDEPENDENT PAIRS                                      CPG16240
c     IF(IPRNT.LE.4) GO TO 301                                          CPG16250
c     WRITE(6,2) IABC                                                   CPG16260
C     CALL MATOUT(EPA,NBASIS,NBASIS,NBASIS,NBASIS,6)                    CPG16270
  301 DO 110 II=1,NIND                                                  CPG16280
      I=NIJ(II,1)                                                       CPG16290
      J=NIJ(II,2)                                                       CPG16300
      ITYP=MOTYP(I)                                                     CPG16310
      JTYP=MOTYP(J)                                                     CPG16320
      IJ=IOFF(I)+J                                                      CPG16330
      VALU1=EPA(I,J)-EPA(J,I)                                           CPG16340
      VALU2=ZERO                                                        CPG16350
C                                                                       CPG16360
      DO 107 K=2,NBASIS                                                 CPG16370
      LIM=MIN0(JOCC,K-1)                                                CPG16380
      DO 107 L=1,LIM                                                    CPG16390
      KL=IOFF(K)+L                                                      CPG16400
      FAC2=ZERO                                                         CPG16410
      IF(K.NE.J) GO TO 205                                              CPG16420
      IL=IOFF(MAX0(I,L))+MIN0(I,L)                                      CPG16430
      FAC2=FAC2+ZETA(IL,ITYP)-ZETA(IL,JTYP)                             CPG16440
  205 IF(K.NE.I) GO TO 206                                              CPG16450
      JL=IOFF(MAX0(J,L))+MIN0(J,L)                                      CPG16460
      FAC2=FAC2-ZETA(JL,JTYP)+ZETA(JL,ITYP)                             CPG16470
  206 continue
      VALU2=VALU2-SM(KL,IABC)*FAC2                                      CPG16480
  107 CONTINUE                                                          CPG16490
      VALUT=VALU1+VALU2                                                 CPG16500
      B0(IJ)=VALUT                                                      CPG16510
C     IF(IPRNT.LE.4) GO TO 110                                          CPG16520
c     WRITE(6,3) II,I,J,IJ,VALU1,VALU2,VALUT
  110 CONTINUE                                                          CPG16540
      CALL RWRIT(ITAP44,B0,NTRI2,IB)                                    CPG16550
c     IF(IPRNT.LE.4) GO TO 120                                          CPG16560
c     WRITE(6,4)                                                        CPG16570
c     CALL PRINT(B0,NTRI,NBASIS,6)                                      CPG16580
  120 CONTINUE                                                          CPG16590
C                                                                       CPG16600
C   FORM DENSITY LIKE MATRIX IN SO BASIS                                CPG16610
      CALL PQBMAT(DP,DQ,SM,ESO,N3N)                                     CPG16620
C                                                                       CPG16630
C   FORM HALF-TRANSFORMED DENSITY MATRIX                                CPG16640
      CALL PQMAT(DP,DQ,DM,N3N*NTYPES,LBLI,BUFI)                         CPG16650
C                                                                       CPG16660
C   FORM DENSITY LIKE MATRIX IN MO BASIS                                CPG16670
      DO 125 IABC=1,N3N                                                 CPG16680
      DO 125 ITYP=1,NTYPES                                              CPG16690
      CALL MOCONS(DM(1,IABC,ITYP),NTRI,DP(1,IABC,ITYP),NTRI,ESO,U,T)    CPG16700
  125 CONTINUE                                                          CPG16710
C                                                                       CPG16720
c     IF(IPRNT.LE.4) GO TO 303                                          CPG16730
c     WRITE(6,5)                                                        CPG16740
c     CALL MATOUT(DP,NTRI,NTYPEP*N3N,NTRI,N3N*NTYPEP,6)                 CPG16750
                                                                        CPG16760
C   COMPLETE B0 MATRIX                                                  CPG16770
  303 CONTINUE                                                          CPG16780
      DO 130 IABC=1,N3N                                                 CPG16790
      IB=IBA(IABC)                                                      CPG16800
      CALL RREAD(ITAP44,B0,NTRI2,IB)                                    CPG16810
      DO 129 II=1,NIND                                                  CPG16820
      I=NIJ(II,1)                                                       CPG16830
      ITYP=MOTYP(I)                                                     CPG16840
      J=NIJ(II,2)                                                       CPG16850
      JTYP=MOTYP(J)                                                     CPG16860
      IJ=IOFF(I)+J                                                      CPG16870
      B0(IJ)=B0(IJ)+DP(IJ,IABC,ITYP)-DP(IJ,IABC,JTYP)                   CPG16880
  129 CONTINUE                                                          CPG16890
      CALL RWRIT(ITAP44,B0,NTRI2,IB)                                    CPG16900
  130 CONTINUE                                                          CPG16910
C                                                                       CPG16920
      CALL SREW(ITAP44)                                                 CPG16930
      RETURN                                                            CPG16940
      END                                                               CPG16950
      SUBROUTINE PQBMAT(DP,DQ,SM,ESO,NABC)                              CPG16960
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG16970
      DIMENSION DP(NTRI,NABC,NTYPEP),DQ(NTRI,NABC,NTYPEP)               CPG16980
      DIMENSION SM(NTRI,NABC*NTYPEP)                                    CPG16990
      DIMENSION ESO(NBFAO,NBASIS)                                       CPG17000
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG17010
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3          CPG17020
      COMMON/GVBAB/A(10,10),B(10,10)                                    CPG17030
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG17040
      COMMON/START/FOCC(10),NSORB(10),NSTR(10),NEND(10),MOTYP(256)      CPG17050
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPG17060
      DATA ZERO / 0.0D+00 /                                             CPG17070
    1 FORMAT(//,2X,' DP MATRIX'/)                                       CPG17080
    2 FORMAT(//,2X,' DQ MATRIX'/)                                       CPG17090
C                                                                       CPG17100
C   FORM DENSITY LIKE MATRIX                                            CPG17110
C   SEE EQ.(A-2)                                                        CPG17120
C                                                                       CPG17130
      DO 101 IABC=1,NABC                                                CPG17140
      DO 101 ITYP=1,NTYPEP                                              CPG17150
      DO 101 I=1,NTRI                                                   CPG17160
      DP(I,IABC,ITYP)=ZERO                                              CPG17170
      DQ(I,IABC,ITYP)=ZERO                                              CPG17180
  101 CONTINUE                                                          CPG17190
      IJ=0                                                              CPG17200
      DO 110 I=1,NBASIS                                                 CPG17210
      DO 110 J=1,I                                                      CPG17220
      IJ=IJ+1                                                           CPG17230
      DO 106 K=1,JOCC                                                   CPG17240
      KTYP=MOTYP(K)                                                     CPG17250
      DO 106 L=1,K                                                      CPG17260
      KL=IOFF(K)+L                                                      CPG17270
      FAC=ESO(I,K)*ESO(J,L)+ESO(J,K)*ESO(I,L)                           CPG17280
      IF(K.EQ.L) FAC=ESO(I,K)*ESO(J,K)                                  CPG17290
      DO 105 ITYP=1,NTYPES                                              CPG17300
      AVAL=A(ITYP,KTYP)                                                 CPG17310
      BVAL=B(ITYP,KTYP)                                                 CPG17320
      FACA=AVAL*FAC                                                     CPG17330
      FACB=BVAL*FAC                                                     CPG17340
      DO 104 IABC=1,NABC                                                CPG17350
      DP(IJ,IABC,ITYP)=DP(IJ,IABC,ITYP)-SM(KL,IABC)*FACA                CPG17360
      DQ(IJ,IABC,ITYP)=DQ(IJ,IABC,ITYP)-SM(KL,IABC)*FACB                CPG17370
  104 CONTINUE                                                          CPG17380
  105 CONTINUE                                                          CPG17390
  106 CONTINUE                                                          CPG17400
  110 CONTINUE                                                          CPG17410
C                                                                       CPG17420
      IF(IPRNT.LE.4) GO TO 201                                          CPG17430
      WRITE(6,1)                                                        CPG17440
      CALL MATOUT(DP,NTRI,NTYPES*NABC,NTRI,NTYPES*NABC,6)               CPG17450
      WRITE(6,2)                                                        CPG17460
      CALL MATOUT(DP,NTRI,NTYPES*NABC,NTRI,NTYPES*NABC,6)               CPG17470
C                                                                       CPG17480
  201 CONTINUE                                                          CPG17490
      RETURN                                                            CPG17500
      END                                                               CPG17510
      SUBROUTINE BADEP(EIG,OCC,ESO,KIJ,U,T,B0,FF,DP,DQ,DM,LBLI,BUFI)    CPG17520
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG17530
      DIMENSION EIG(NBASIS),OCC(NBASIS),KIJ(NTRI,2)                     CPG17540
      DIMENSION ESO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)         CPG17550
      DIMENSION B0(NTRI),FF(NTRI)                                       CPG17560
      DIMENSION DP(NTRI,N3N),DQ(NTRI,N3N),DM(NTRI,N3N)                  CPG17570
      DIMENSION LBLI(MAXBF4),BUFI(MAXBF2)                               CPG17580
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG17590
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3          CPG17600
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPG17610
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPG17620
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG17630
      COMMON/TAPES/ITAP36,ITAP42,ITAP43,ITAP44                          CPG17640
      COMMON/CI102/NIND,NDEP                                            CPG17650
    1 FORMAT(//,2X,' B0 MATRIX IN BADEP'/)                              CPG17660
    2 FORMAT(//,2X,' DP MATRIX IN BADEP'/)                              CPG17670
C                                                                       CPG17680
      CALL SREW(ITAP44)                                                 CPG17690
C   DM ARRAY STORES DERIVATIVE OVERLAP MATRIX                           CPG17700
      DO 102 IABC=1,N3N                                                 CPG17710
      IS=ISA(IABC)                                                      CPG17720
      IT=IFA(IABC)                                                      CPG17730
      IB=IBA(IABC)                                                      CPG17740
      CALL RREAD(ITAP44,DM(1,IABC),NTRI2,IS)                            CPG17750
      CALL RREAD(ITAP44,FF,NTRI2,IT)                                    CPG17760
      CALL RREAD(ITAP44,B0,NTRI2,IB)                                    CPG17770
      IJ=0                                                              CPG17780
      DO 101 II=1,NDEP                                                  CPG17790
      I=KIJ(II,1)                                                       CPG17800
      J=KIJ(II,2)                                                       CPG17810
      IJ=IOFF(I)+J                                                      CPG17820
      B0(IJ)=FF(IJ)-EIG(J)*DM(IJ,IABC)                                  CPG17830
  101 CONTINUE                                                          CPG17840
      CALL RWRIT(ITAP44,B0,NTRI2,IB)                                    CPG17850
      IF(IPRNT.LE.4) GO TO 102                                          CPG17860
      WRITE(6,1)                                                        CPG17870
      CALL PRINT(B0,NTRI,NBASIS,6)                                      CPG17880
  102 CONTINUE                                                          CPG17890
C                                                                       CPG17900
C   FORM DENSITY LIKE MATRIX                                            CPG17910
C   DM ARRAY CONTAINES DERIVATIVE OVERLAP MATRIX                        CPG17920
C   DP ARRAY STORES DENSITY LIKE MATRIX IN SO BASIS                     CPG17930
  301 CONTINUE                                                          CPG17940
      CALL PQCMAT(DP,DQ,DM,OCC,ESO,N3N)                                 CPG17950
C                                                                       CPG17960
C   FORM HALF-TRANSFORMED DENSITY LIKE MATRIX                           CPG17970
C   AND TRANSFORM THEM INTO MO BASIS                                    CPG17980
C   DP ARRAY CONTAINS DENSITY LIKE MATRIX IN SO BASIS                   CPG17990
C   DM ARRAY STORES HALF-TRANSFORMED DENSITY LIKE MATRIX                CPG18000
      CALL PQMAT(DP,DQ,DM,N3N,LBLI,BUFI)                                CPG18010
      IF(IPRNT.LE.4) GO TO 302                                          CPG18020
      WRITE(6,2)                                                        CPG18030
      CALL MATOUT(DP,NTRI,N3N,NTRI,N3N,6)                               CPG18040
C                                                                       CPG18050
C   COMPLETE B0 MATRICES                                                CPG18060
  302 CONTINUE                                                          CPG18070
      CALL SREW(ITAP44)                                                 CPG18080
      DO 105 IABC=1,N3N                                                 CPG18090
      IB=IBA(IABC)                                                      CPG18100
      CALL RREAD(ITAP44,B0,NTRI2,IB)                                    CPG18110
      CALL MOCONS(DM(1,IABC),NTRI,DP(1,IABC),NTRI,ESO,U,T)              CPG18120
      DO 104 II=1,NDEP                                                  CPG18130
      I=KIJ(II,1)                                                       CPG18140
      J=KIJ(II,2)                                                       CPG18150
      IJ=IOFF(I)+J                                                      CPG18160
      B0(IJ)=B0(IJ)+DP(IJ,IABC)                                         CPG18170
  104 CONTINUE                                                          CPG18180
      CALL RWRIT(ITAP44,B0,NTRI2,IB)                                    CPG18190
      IF(IPRNT.LE.4) GO TO 105                                          CPG18200
      WRITE(6,1)                                                        CPG18210
      CALL PRINT(B0,NTRI,N3N,6)   
  105 CONTINUE                                                          CPG18230
C                                                                       CPG18240
  303 CONTINUE                                                          CPG18250
      CALL SREW(ITAP44)                                                 CPG18260
      RETURN                                                            CPG18270
      END                                                               CPG18280
      SUBROUTINE PQCMAT(DP,DQ,SM,OCC,ESO,NABC)                          CPG18290
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG18300
      DIMENSION DP(NTRI,NABC),DQ(NTRI,NABC),SM(NTRI,NABC)               CPG18310
      DIMENSION OCC(NBASIS),ESO(NBFAO,NBASIS)                           CPG18320
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG18330
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG18340
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPG18350
      DATA ZERO,HALF / 0.0D+00 , 0.5D+00 /                              CPG18360
    1 FORMAT(//,2X,' DP MATRIX'/)                                       CPG18370
    2 FORMAT(//,2X,' DQ MATRIX'/)                                       CPG18380
C                                                                       CPG18390
C   FORM DENSITY LIKE MATRICES                                          CPG18400
C                                                                       CPG18410
      DO 101 IABC=1,NABC                                                CPG18420
      DO 101 I=1,NTRI                                                   CPG18430
      DP(I,IABC)=ZERO                                                   CPG18440
      DQ(I,IABC)=ZERO                                                   CPG18450
  101 CONTINUE                                                          CPG18460
C                                                                       CPG18470
      IJ=0                                                              CPG18480
      DO 107 I=1,NBASIS                                                 CPG18490
      DO 107 J=1,I                                                      CPG18500
      IJ=IJ+1                                                           CPG18510
      DO 106 K=1,JOCC                                                   CPG18520
      FACO=OCC(K)*HALF                                                  CPG18530
      DO 106 L=1,K                                                      CPG18540
      KL=IOFF(MAX0(K,L))+MIN0(K,L)                                      CPG18550
      FACV=ESO(I,K)*ESO(J,L)+ESO(J,K)*ESO(I,L)                          CPG18560
      IF(K.EQ.L) FACV=ESO(I,K)*ESO(J,K)                                 CPG18570
      FAC=FACO*FACV                                                     CPG18580
      DO 105 IABC=1,NABC                                                CPG18590
      FACT=-FAC*SM(KL,IABC)                                             CPG18600
      FAC2=FACT+FACT                                                    CPG18610
      DP(IJ,IABC)=DP(IJ,IABC)+FAC2                                      CPG18620
      DQ(IJ,IABC)=DQ(IJ,IABC)-FACT                                      CPG18630
  105 CONTINUE                                                          CPG18640
  106 CONTINUE                                                          CPG18650
  107 CONTINUE                                                          CPG18660
C                                                                       CPG18670
      IF(IPRNT.LE.4) GO TO 201                                          CPG18680
      WRITE(6,1)                                                        CPG18690
      CALL MATOUT(DP,NTRI,NABC,NTRI,NABC,6)                             CPG18700
      WRITE(6,2)                                                        CPG18710
      CALL MATOUT(DQ,NTRI,NABC,NTRI,NABC,6)                             CPG18720
C                                                                       CPG18730
  201 CONTINUE                                                          CPG18740
      RETURN                                                            CPG18750
      END                                                               CPG18760
      SUBROUTINE BFMAT(NIJ,KIJ,HF,HM,EPF,BF)                            CPG18770
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG18780
      DIMENSION NIJ(NTRI,2),KIJ(NTRI,2)                                 CPG18790
      DIMENSION HF(NTRI),HM(NTRI,3,NTYPES),EPF(NBASIS,NBASIS),BF(NTRI)  CPG18800
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG18810
      COMMON/DIPOL/DIPDA(135),DIPDM(135),DIPDN(135),DIPDT(135)          CPG18820
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3          CPG18830
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPG18840
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG18850
      COMMON/START/FOCC(10),NSORB(10),NSTR(10),NEND(10),MOTYP(256)      CPG18860
      COMMON/TAPES/ITAP36,ITAP42,ITAP43,ITAP44                          CPG18870
      COMMON/CI102/NIND,NDEP                                            CPG18880
      DATA ZERO,HALF / 0.0D+00 , 0.5D+00 /                              CPG18890
    1 FORMAT(//,2X,' DIPDA MATRIX'/)                                    CPG18900
    2 FORMAT(//,2X,' DIPDN MATRIX'/)                                    CPG18910
    3 FORMAT(//,2X,' HF MATRIX, IXYZ = ',I5/)                           CPG18920
    4 FORMAT(//,2X,' EPF MATRIX, IXYZ = ',I5/)                          CPG18930
C                                                                       CPG18940
      CALL RFILE(ITAP43)                                                CPG18950
      CALL SREAD(ITAP43,DIPDA,270)                                      CPG18960
      CALL SREAD(ITAP43,DIPDN,270)                                      CPG18970
      IF(IPRNT.LE.4) GO TO 201                                          CPG18980
      WRITE(6,1)                                                        CPG18990
      CALL MATOUT(DIPDA,3,45,3,N3N,6)                                   CPG19000
      WRITE(6,2)                                                        CPG19010
      CALL MATOUT(DIPDN,3,45,3,N3N,6)                                   CPG19020
C                                                                       CPG19030
  201 CONTINUE                                                          CPG19040
      CALL SREW(ITAP44)                                                 CPG19050
      DO 103 IXYZ=1,3                                                   CPG19060
      IH=IHA(IXYZ+N3N)                                                  CPG19070
      CALL SREAD(ITAP43,HF,NTRI2)                                       CPG19080
      CALL RWRIT(ITAP44,HF,NTRI2,IH)                                    CPG19090
      DO 102 ITYP=1,NTYPES                                              CPG19100
      FAC=FOCC(ITYP)*HALF                                               CPG19110
      DO 101 I=1,NTRI                                                   CPG19120
  101 HM(I,IXYZ,ITYP)=HF(I)*FAC                                         CPG19130
  102 CONTINUE                                                          CPG19140
      IF(IPRNT.LE.4) GO TO 103                                          CPG19150
      WRITE(6,3) IXYZ                                                   CPG19160
      CALL PRINT(HF,NTRI,NBASIS,6)                                      CPG19170
  103 CONTINUE                                                          CPG19180
C                                                                       CPG19190
C   STORE DERIVATIVE LAGRANGIAN MATRICES                                CPG19200
      DO 110 IXYZ=1,3                                                   CPG19210
      IH=IHA(IXYZ+N3N)                                                  CPG19220
      IE=IEA(IXYZ+N3N)                                                  CPG19230
      IB=IBA(IXYZ+N3N)                                                  CPG19240
      CALL RREAD(ITAP44,HF,NTRI2,IH)                                    CPG19250
      DO 104 I=1,NTRI                                                   CPG19260
  104 BF(I)=ZERO                                                        CPG19270
      DO 105 I=1,NBASIS                                                 CPG19280
      DO 105 J=1,NBASIS                                                 CPG19290
  105 EPF(I,J)=ZERO                                                     CPG19300
C                                                                       CPG19310
      DO 107 ITYP=1,NTYPES                                              CPG19320
      DO 106 I=NSTR(ITYP),NEND(ITYP)                                    CPG19330
      DO 106 J=1,NBASIS                                                 CPG19340
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPG19350
      EPF(I,J)=HM(IJ,IXYZ,ITYP)                                         CPG19360
  106 CONTINUE                                                          CPG19370
  107 CONTINUE                                                          CPG19380
C                                                                       CPG19390
C   ELEMENTS FOR INDEPENDENT PAIRS                                      CPG19400
      DO 108 II=1,NIND                                                  CPG19410
      I=NIJ(II,1)                                                       CPG19420
      J=NIJ(II,2)                                                       CPG19430
      IJ=IOFF(I)+J                                                      CPG19440
      BF(IJ)=EPF(I,J)-EPF(J,I)                                          CPG19450
  108 CONTINUE                                                          CPG19460
C                                                                       CPG19470
C   ELEMENTS FOR DEPENDENT PAIRS                                        CPG19480
      DO 109 II=1,NDEP                                                  CPG19490
      I=KIJ(II,1)                                                       CPG19500
      J=KIJ(II,2)                                                       CPG19510
      IJ=IOFF(I)+J                                                      CPG19520
      BF(IJ)=HF(IJ)                                                     CPG19530
  109 CONTINUE                                                          CPG19540
C                                                                       CPG19550
  202 CONTINUE                                                          CPG19560
      CALL RWRIT(ITAP44,BF,NTRI2,IB)                                    CPG19570
      CALL RWRIT(ITAP44,EPF,NBASQ,IE)                                   CPG19580
      IF(IPRNT.LE.4) GO TO 110                                          CPG19590
      WRITE(6,4) IXYZ                                                   CPG19600
      CALL MATOUT(EPF,NBASIS,NBASIS,NBASIS,NBASIS,6)                    CPG19610
  110 CONTINUE                                                          CPG19620
      CALL SREW(ITAP43)                                                 CPG19630
      CALL SREW(ITAP44)                                                 CPG19640
      call rclose(itap43,3)
C                                                                       CPG19650
      RETURN                                                            CPG19660
      END                                                               CPG19670
      SUBROUTINE CPHF(ESO,NIJ,KIJ,U,A,T,ZETA,B0,AA,TT,DP,DQ,DM,W,BB,    CPG19680
     1                LBLI,BUFI,NMAX,NADD,N3NXX)                        CPG19690
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG19700
cyy   REAL*16 DLIM,DETERM                                               CPG19710
      REAL*8  DLIM,DETERM                                               CPG19710
      DIMENSION ESO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)         CPG19720
      DIMENSION A(NBFAO,NBFAO),ZETA(NTRI,NTYPEP)                        CPG19730
      DIMENSION BB(NIND,NMAX*2),W(NMAX,NMAX+NADD)                       CPG19740
      DIMENSION DP(NTRI,NADD,NTYPEP),DQ(NTRI,NADD,NTYPEP)               CPG19750
      DIMENSION DM(NTRI,NADD,NTYPEP)                                    CPG19760
      DIMENSION AA(NTRI),TT(NTRI),B0(NTRI),NIJ(NTRI,2),KIJ(NTRI,2)      CPG19770
      DIMENSION LBLI(MAXBF4),BUFI(MAXBF2)                               CPG19780
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG19790
      COMMON/CALIF/LPARA(1024),BNORM(1024)                              CPG19800
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3          CPG19810
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPG19820
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPG19830
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG19840
      COMMON/START/FOCC(10),NSORB(10),NSTR(10),NEND(10),MOTYP(256)      CPG19850
      COMMON/TAPES/ITAP36,ITAP42,ITAP43,ITAP44                          CPG19860
      COMMON/TOLER/ICONV                                                CPG19870
      COMMON/WORKS/JTAPE1                                               CPG19880
      COMMON/CI102/NIND,NDEP                                            CPG19890
      DATA ZERO,ONE,TEN / 0.0D+00 , 1.0D+00 , 10.0D+00 /                CPG19900
    1 FORMAT(//,2X,' A MATRIX'/)                                        CPG19910
    2 FORMAT(//,2X,' IGROUP = ',I5,5X,' IFIRST = ',I5,5X,               CPG19920
     1             ' NABC   = ',I5/)                                    CPG19930
    3 FORMAT(//,2X,' BB MATRIX'/)                                       CPG19940
    4 FORMAT(//,2X,' SCALED BB MATRIX'/)                                CPG19950
    5 FORMAT(//,2X,' ORTHONORMALIZED BB MATRIX'/)                       CPG19960
    6 FORMAT(//,2X,' ORTHONORMALITY OF BB MATRIX'/)                     CPG19970
    7 FORMAT(//,2X,' ITERATION NO. = ',I5/)                             CPG19980
    8 FORMAT(//,2X,' (1-A)*B MATRIX'/)                                  CPG19990
    9 FORMAT(//,2X,' KVEC = ',I5,' K = ',I5,' TT = ',(10F10.5))         CPG20000
   10 FORMAT(/,2X,' NITER  = ',I5,5X,' NVEC   = ',I5,5X,' NADDV  = ',   CPG20010
     1         I5/)                                                     CPG20020
   11 FORMAT(//,2X,' NORM OF ADDED VECTOR',(5(I5,F15.8)))               CPG20030
   12 FORMAT(//,2X,' BB MATRIX'/)                                       CPG20040
   13 FORMAT(/,2X,' END OF ITERATION    NVTOT  = ',I5/)                 CPG20050
   14 FORMAT(//,2X,'***NVTOT EXCEEDS NMAX***',5X,' NMAX = ',I5/)        CPG20060
   15 FORMAT(//,2X,' EXPANDED VECTORS BB AND A*BB MATRICES'/)           CPG20070
   16 FORMAT(//,2X,' W MATRIX'/)                                        CPG20080
   17 FORMAT(//,2X,' SCALED W MATRIX'/)                                 CPG20090
   18 FORMAT(//,2X,'***SIMULTANEOUS EQUATION CANNOT BE SOLVED***',      CPG20100
     1 5X,' DETERM = ',D20.10/)                                         CPG20110
cyy  1 5X,' DETERM = ',Q20.10/)                                         CPG20110
   19 FORMAT(/,2X,' LINEAR EQUATION HAS BEEN SOLVED',                   CPG20120
     1 5X,' DETERM = ',D20.10/)                                         CPG20130
cyy  1 5X,' DETERM = ',Q20.10/)                                         CPG20130
   20 FORMAT(//,2X,' SOLUTION MATRIX OF SCALED W MATRIX'/)              CPG20140
   21 FORMAT(//,2X,' SOLUTION MATRIX'/)                                 CPG20150
   22 FORMAT(//,2X,' B0 MATRIX WITH INDEP. OF UA, IABC = ',I5/)         CPG20160
C                                                                       CPG20170
      CRIT=1.0D-10                                                      CPG20180
      IF(ICONV.NE.0) CRIT=ONE/(TEN**ICONV)                              CPG20190
      DLIM=1.0D-30                                                      CPG20200
cyy   DLIM=1.0Q-76                                                      CPG20200
      NIND2=NIND*2                                                      CPG20210
C                                                                       CPG20220
      DO 101 ITYP=1,NTYPES                                              CPG20230
      CALL MREAD(ZETA(1,ITYP),26+ITYP)                                  CPG20240
C     CALL PRINT(ZETA(1,ITYP),NTRI,NBASIS,6)                            CPG20250
  101 CONTINUE                                                          CPG20260
      DO 102 I=1,NTRI                                                   CPG20270
  102 ZETA(I,NTYPEP)=ZERO                                               CPG20280
C                                                                       CPG20290
C   FORM A MATRIX                                                       CPG20300
C   SEE EQ.(A-1)                                                        CPG20310
C                                                                       CPG20320
      DO 103 I=1,NBASIS                                                 CPG20330
      DO 103 J=1,NBASIS                                                 CPG20340
  103 A(I,J)=ZERO                                                       CPG20350
      DO 104 I=1,NBASIS                                                 CPG20360
      DO 104 J=1,NBASIS                                                 CPG20370
      ITYP=MOTYP(I)                                                     CPG20380
      JTYP=MOTYP(J)                                                     CPG20390
      II=IOFF(I+1)                                                      CPG20400
      JJ=IOFF(J+1)                                                      CPG20410
      A(I,J)=A(I,J)-ZETA(II,ITYP)+ZETA(II,JTYP)                         CPG20420
     1             -ZETA(JJ,JTYP)+ZETA(JJ,ITYP)                         CPG20430
  104 CONTINUE                                                          CPG20440
      IF(IPRNT.LE.4) GO TO 301                                          CPG20450
      WRITE(6,1)                                                        CPG20460
      CALL MATOUT(A,NBFAO,NBFAO,NBASIS,NBASIS,6)                        CPG20470
C                                                                       CPG20480
C   SCALING MATRIX                                                      CPG20490
C   SEE EQS.(24) AND (25)                                               CPG20500
C                                                                       CPG20510
  301 CONTINUE                                                          CPG20520
      DO 105 II=1,NIND                                                  CPG20530
      I=NIJ(II,1)                                                       CPG20540
      J=NIJ(II,2)                                                       CPG20550
      AVAL=A(I,J)                                                       CPG20560
      AVAL=DABS(AVAL)                                                   CPG20570
      AA(II)=ONE/DSQRT(AVAL)                                            CPG20580
  105 CONTINUE                                                          CPG20590
C                                                                       CPG20600
C   START GROUPING THE B VECTORS                                        CPG20610
      IGROUP=0                                                          CPG20620
      ILAST=0                                                           CPG20630
  500 CONTINUE                                                          CPG20640
      IGROUP=IGROUP+1                                                   CPG20650
      IFIRST=ILAST+1                                                    CPG20660
      ILAST=IFIRST+NADD-1                                               CPG20670
      IF(ILAST.GT.N3NXX) ILAST=N3NXX                                    CPG20680
      NABC=ILAST-IFIRST+1                                               CPG20690
      WRITE(6,2) IGROUP,IFIRST,NABC                                     CPG20700
      CALL SREW(JTAPE1)                                                 CPG20710
C                                                                       CPG20720
C   CALL BACK B0 MATRICES                                               CPG20730
      CALL SREW(ITAP44)                                                 CPG20740
      DO 107 IABC=1,NABC                                                CPG20750
      III=IABC+IFIRST-1                                                 CPG20760
      IB=IBA(III)                                                       CPG20770
      CALL RREAD(ITAP44,B0,NTRI2,IB)                                    CPG20780
      DO 106 II=1,NIND                                                  CPG20790
      I=NIJ(II,1)                                                       CPG20800
      J=NIJ(II,2)                                                       CPG20810
      IJ=IOFF(I)+J                                                      CPG20820
      BB(II,IABC)=B0(IJ)                                                CPG20830
  106 CONTINUE                                                          CPG20840
  107 CONTINUE                                                          CPG20850
      IF(IPRNT.LE.4) GO TO 302                                          CPG20860
      WRITE(6,3)                                                        CPG20870
      CALL MATOUT(BB,NIND,NABC,NIND,NABC,6)                             CPG20880
C                                                                       CPG20890
C   SCALE B0 MATRICES                                                   CPG20900
  302 CONTINUE                                                          CPG20910
      DO 109 II=1,NIND                                                  CPG20920
      DO 108 IABC=1,NABC                                                CPG20930
  108 BB(II,IABC)=BB(II,IABC)*AA(II)                                    CPG20940
  109 CONTINUE                                                          CPG20950
      IF(IPRNT.LE.4) GO TO 303                                          CPG20960
      WRITE(6,4)                                                        CPG20970
      CALL MATOUT(BB,NIND,NABC,NIND,NABC,6)                             CPG20980
C                                                                       CPG20990
C   ORTHONORMALIZE B0 MATRIX                                            CPG21000
C                                                                       CPG21010
  303 CONTINUE                                                          CPG21020
      CALL ORTHOG(BB,BB(1,NABC+1),NIND,NABC,NIND,NABC,NVEC)             CPG21030
      IF(IPRNT.LE.4) GO TO 304                                          CPG21040
      WRITE(6,5)                                                        CPG21050
      CALL MATOUT(BB,NIND,NVEC,NIND,NVEC,6)                             CPG21060
  304 CONTINUE                                                          CPG21070
      DO 110 IABC=1,NVEC                                                CPG21080
  110 BNORM(IABC)=ONE                                                   CPG21090
C                                                                       CPG21100
C   CHECK ORTHONORMALITY                                                CPG21110
      IF(IPRNT.LE.4) GO TO 305                                          CPG21120
      DO 112 I=1,NVEC                                                   CPG21130
      DO 112 J=1,NVEC                                                   CPG21140
      VALU=ZERO                                                         CPG21150
      DO 111 K=1,NIND                                                   CPG21160
  111 VALU=VALU+BB(K,I)*BB(K,J)                                         CPG21170
      W(I,J)=VALU                                                       CPG21180
  112 CONTINUE                                                          CPG21190
      WRITE(6,6)                                                        CPG21200
      CALL MATOUT(W,NMAX,NMAX+NABC,NVEC,NVEC,6)                         CPG21210
C                                                                       CPG21220
  305 NPVEC=0                                                           CPG21230
      NITER=0                                                           CPG21240
C:::::::::::::::::::::                                                  CPG21250
C:::START ITERATION:::                                                  CPG21260
C:::::::::::::::::::::                                                  CPG21270
  200 NITER=NITER+1                                                     CPG21280
      IF(IPRNT.GT.2)                                                    CPG21290
     *WRITE(6,7) NITER                                                  CPG21300
C                                                                       CPG21310
      DO 114 II=1,NIND                                                  CPG21320
      DO 113 IABC=1,NVEC                                                CPG21330
  113 BB(II,NPVEC+IABC)=BB(II,NPVEC+IABC)*AA(II)                        CPG21340
  114 CONTINUE                                                          CPG21350
C                                                                       CPG21360
C   FORM A*B MATRIX                                                     CPG21370
      CALL MAKEAB(BB(1,NPVEC+1),DP,DQ,DM,BB(1,NPVEC+NVEC+1),NVEC,       CPG21380
     1 ZETA,NIJ,ESO,U,T,LBLI,BUFI)                                      CPG21390
C                                                                       CPG21400
      DO 116 II=1,NIND                                                  CPG21410
      DO 115 IABC=1,NVEC                                                CPG21420
      BB(II,NPVEC+IABC)=BB(II,NPVEC+IABC)/AA(II)                        CPG21430
      BB(II,NPVEC+NVEC+IABC)=BB(II,NPVEC+NVEC+IABC)*AA(II)              CPG21440
  115 CONTINUE                                                          CPG21450
  116 CONTINUE                                                          CPG21460
cyy     write(6,*) ' updated BB matrix'
cyy     CALL MATOUT(BB,NIND,NPVEC+NVEC*2,NIND,NPVEC+NVEC*2,6) 
C                                                                       CPG21470
C   STORE A*B MATRIX                                                    CPG21480
      DO 117 IABC=NPVEC+NVEC+1,NPVEC+NVEC*2                             CPG21490
cyy   write(6,*)  ' iabc = ', iabc
cyy   call matout(BB(1,IABC),nind, 1, nind, 1, 6)
      CALL SWRIT(JTAPE1,BB(1,IABC),NIND2)                               CPG21500
  117 CONTINUE                                                          CPG21510
C                                                                       CPG21520
C   FORM (1-A)*B                                                        CPG21530
      DO 120 IABC=NPVEC+1,NPVEC+NVEC                                    CPG21540
      DO 119 II=1,NIND                                                  CPG21550
C*119 BB(II,NVEC+IABC)=-BB(II,IABC)-BB(II,NVEC+IABC)                    CPG21560
  119 BB(II,NVEC+IABC)=BB(II,IABC)-BB(II,NVEC+IABC)                     CPG21570
  120 CONTINUE                                                          CPG21580
      IF(IPRNT.LE.4) GO TO 306                                          CPG21590
      WRITE(6,8)                                                        CPG21600
      CALL MATOUT(BB(1,NPVEC+NVEC+1),NIND,NVEC,NIND,NVEC,6)             CPG21610
cyy   write(6,*) '  (1-a) * b '
cyy   CALL MATOUT(BB,NIND,NPVEC+NVEC*2,NIND,NPVEC+NVEC*2,6)  
C                                                                       CPG21620
C   FORM ADDITIONAL EXPANSION VECTORS                                   CPG21630
  306 NADDV=0                                                           CPG21640
      K=NPVEC+NVEC                                                      CPG21650
      DO 130 IABC=NPVEC+NVEC+1,NPVEC+NVEC*2                             CPG21660
C                                                                       CPG21670
      DO 121 II=1,NIND                                                  CPG21680
  121 TT(II)=BB(II,IABC)                                                CPG21690
      DO 124 L=1,K                                                      CPG21700
      CALL NORM(XBB,BB(1,L),BB(1,IABC),NIND)                            CPG21710
      FAC=XBB/BNORM(L)                                                  CPG21720
      DO 123 II=1,NIND                                                  CPG21730
  123 TT(II)=TT(II)-BB(II,L)*FAC                                        CPG21740
  124 CONTINUE                                                          CPG21750
      IF(IPRNT.LE.4) GO TO 307                                          CPG21760
      WRITE(6,9) IABC,K,(TT(I),I=1,NIND)                                CPG21770
  307 CONTINUE                                                          CPG21780
      CALL NORM(XTT,TT,TT,NIND)                                         CPG21790
      IF(DABS(XTT).LT.CRIT) GO TO 130                                   CPG21800
      NADDV=NADDV+1                                                     CPG21810
      K=K+1                                                             CPG21820
      BNORM(K)=XTT                                                      CPG21830
      DO 126 II=1,NIND                                                  CPG21840
C*126 BB(II,K)=-TT(II)                                                  CPG21850
  126 BB(II,K)=TT(II)                                                   CPG21860
C                                                                       CPG21870
  130 CONTINUE                                                          CPG21880
C                                                                       CPG21890
      WRITE(6,10) NITER,NVEC,NADDV                                      CPG21900
      WRITE(3,10) NITER,NVEC,NADDV                                      CPG21910
      NPVEC=NPVEC+NVEC                                                  CPG21920
      NVEC=NADDV                                                        CPG21930
      IF(NADDV.EQ.0) GO TO 310                                          CPG21940
      IF(IPRNT.LE.4) GO TO 308                                          CPG21950
      IST=NPVEC+1                                                       CPG21960
      IED=NPVEC+NVEC                                                    CPG21970
      WRITE(6,11) (I,BNORM(I),I=IST,IED)                                CPG21980
      WRITE(6,12)                                                       CPG21990
      CALL MATOUT(BB(1,NPVEC+1),NIND,NVEC,NIND,NVEC,6)                  CPG22000
  308 GO TO 200                                                         CPG22010
C                                                                       CPG22020
  310 NVTOT=NPVEC                                                       CPG22030
      WRITE(6,13) NVTOT                                                 CPG22040
      WRITE(3,13) NVTOT                                                 CPG22050
      IF(NVTOT.LE.NMAX) GO TO 400                                       CPG22060
      WRITE(6,14) NMAX                                                  CPG22070
      STOP                                                              CPG22080
C                                                                       CPG22090
C   CALL BACK A*B MATRICES                                              CPG22100
  400 CONTINUE                                                          CPG22110
      CALL SREW(JTAPE1)                                                 CPG22120
      DO 131 IABC=NVTOT+1,NVTOT*2                                       CPG22130
  131 CALL SREAD(JTAPE1,BB(1,IABC),NIND2)                               CPG22140
      IF(IPRNT.LE.4) GO TO 311                                          CPG22150
      WRITE(6,15)                                                       CPG22160
      CALL MATOUT(BB,NIND,NVTOT*2,NIND,NVTOT*2,6)                       CPG22170
C                                                                       CPG22180
  311 CONTINUE                                                          CPG22190
      DO 133 I=1,NVTOT                                                  CPG22200
      DO 133 J=1,NVTOT                                                  CPG22210
      SUM=ZERO                                                          CPG22220
      DO 132 II=1,NIND                                                  CPG22230
C*132 SUM=SUM-BB(II,I)*BB(II,NVTOT+J)                                   CPG22240
  132 SUM=SUM+BB(II,I)*BB(II,NVTOT+J)                                   CPG22250
      W(I,J)=SUM                                                        CPG22260
  133 CONTINUE                                                          CPG22270
C                                                                       CPG22280
C   READ BACK ORIGINAL B0 MATRICES                                      CPG22290
      CALL SREW(ITAP44)                                                 CPG22300
      DO 135 IABC=1,NABC                                                CPG22310
      III=IABC+IFIRST-1                                                 CPG22320
      IB=IBA(III)                                                       CPG22330
      CALL RREAD(ITAP44,B0,NTRI2,IB)                                    CPG22340
      DO 134 II=1,NIND                                                  CPG22350
      I=NIJ(II,1)                                                       CPG22360
      J=NIJ(II,2)                                                       CPG22370
      IJ=IOFF(I)+J                                                      CPG22380
      BB(II,NVTOT+IABC)=B0(IJ)*AA(II)                                   CPG22390
  134 CONTINUE                                                          CPG22400
  135 CONTINUE                                                          CPG22410
      IF(IPRNT.LE.4) GO TO 312                                          CPG22420
      WRITE(6,4)                                                        CPG22430
      CALL MATOUT(BB(1,NVTOT+1),NIND,NABC,NIND,NABC,6)                  CPG22440
C                                                                       CPG22450
C   SET B0*B ON W                                                       CPG22460
  312 CONTINUE                                                          CPG22470
      DO 137 I=1,NVTOT                                                  CPG22480
      DO 137 J=1,NABC                                                   CPG22490
      SUM=ZERO                                                          CPG22500
      DO 136 II=1,NIND                                                  CPG22510
  136 SUM=SUM+BB(II,I)*BB(II,NVTOT+J)                                   CPG22520
      W(I,NVTOT+J)=SUM                                                  CPG22530
  137 CONTINUE                                                          CPG22540
C                                                                       CPG22550
      IF(IPRNT.LE.3) GO TO 313                                          CPG22560
      WRITE(6,16)                                                       CPG22570
      CALL MATOUT(W,NMAX,NMAX+NABC,NVTOT,NVTOT+NABC,6)                  CPG22580
  313 CONTINUE                                                          CPG22590
      DO 138 I=1,NVTOT                                                  CPG22600
      AVAL=ONE/DSQRT(BNORM(I))                                          CPG22610
      DO 138 J=1,NVTOT                                                  CPG22620
      BVAL=ONE/DSQRT(BNORM(J))                                          CPG22630
      W(I,J)=W(I,J)*AVAL*BVAL                                           CPG22640
  138 CONTINUE                                                          CPG22650
      DO 139 I=1,NVTOT                                                  CPG22660
      AVAL=ONE/DSQRT(BNORM(I))                                          CPG22670
      DO 139 J=1,NABC                                                   CPG22680
      W(I,NVTOT+J)=W(I,NVTOT+J)*AVAL                                    CPG22690
  139 CONTINUE                                                          CPG22700
      IF(IPRNT.LE.3) GO TO 315                                          CPG22710
      WRITE(6,17)                                                       CPG22720
      CALL MATOUT(W,NMAX,NMAX+NABC,NVTOT,NVTOT+NABC,6)                  CPG22730
C                                                                       CPG22740
C   SOLVE SMALL SIMULTANEOUS EQUATION                                   CPG22750
  315 CONTINUE                                                          CPG22760
      CALL FLINQ(W,NMAX,NVTOT,NABC,DETERM)                              CPG22770
      IF(DABS(DETERM).GT.DLIM) GO TO 320                                CPG22780
cyy   IF(QABS(DETERM).GT.DLIM) GO TO 320                                CPG22780
C                                                                       CPG22790
C   SIMULTANEOUS EQUATION CANNOT BE SOLVED---STOP!!!                    CPG22800
      WRITE(6,18) DETERM                                                CPG22810
      STOP                                                              CPG22820
C                                                                       CPG22830
C::::::::::::::::::::::::::::::::                                       CPG22840
C:::LINEAR EQUATION WAS SOLVED:::                                       CPG22850
C::::::::::::::::::::::::::::::::                                       CPG22860
  320 WRITE(6,19) DETERM                                                CPG22870
      IF(IPRNT.LE.4) GO TO 322                                          CPG22880
      WRITE(6,20)                                                       CPG22890
      CALL MATOUT(W(1,NVTOT+1),NMAX,NMAX+NABC,NVTOT,NABC,6)             CPG22900
  322 CONTINUE                                                          CPG22910
      DO 151 I=1,NVTOT                                                  CPG22920
      AVAL=ONE/DSQRT(BNORM(I))                                          CPG22930
      DO 151 J=1,NABC                                                   CPG22940
      W(I,NVTOT+J)=W(I,NVTOT+J)*AVAL                                    CPG22950
  151 CONTINUE                                                          CPG22960
      IF(IPRNT.LE.4) GO TO 323                                          CPG22970
      WRITE(6,21)                                                       CPG22980
      CALL MATOUT(W(1,NVTOT+1),NMAX,NMAX+NABC,NVTOT,NABC,6)             CPG22990
C                                                                       CPG23000
C   FORM INDEPENDENT ELEMENTS OF UA MATRIX                              CPG23010
  323 CONTINUE                                                          CPG23020
      DO 153 II=1,NIND                                                  CPG23030
      DO 153 IABC=1,NABC                                                CPG23040
      SUM=ZERO                                                          CPG23050
      DO 152 I=1,NVTOT                                                  CPG23060
C*152 SUM=SUM-BB(II,I)*W(I,NVTOT+IABC)                                  CPG23070
  152 SUM=SUM+BB(II,I)*W(I,NVTOT+IABC)                                  CPG23080
      BB(II,NVTOT+IABC)=SUM*AA(II)                                      CPG23090
  153 CONTINUE                                                          CPG23100
C                                                                       CPG23110
      CALL SREW(ITAP44)                                                 CPG23120
      DO 155 IABC=1,NABC                                                CPG23130
      III=IABC+IFIRST-1                                                 CPG23140
      IB=IBA(III)                                                       CPG23150
      CALL RREAD(ITAP44,B0,NTRI2,IB)                                    CPG23160
      DO 154 II=1,NIND                                                  CPG23170
      I=NIJ(II,1)                                                       CPG23180
      J=NIJ(II,2)                                                       CPG23190
      IJ=IOFF(I)+J                                                      CPG23200
      B0(IJ)=BB(II,NVTOT+IABC)                                          CPG23210
  154 CONTINUE                                                          CPG23220
      CALL RWRIT(ITAP44,B0,NTRI2,IB)                                    CPG23230
      IF(IPRNT.LE.3) GO TO 155                                          CPG23240
      WRITE(6,22) IABC                                                  CPG23250
      CALL PRINT(B0,NTRI,NBASIS,6)                                      CPG23260
  155 CONTINUE                                                          CPG23270
C                                                                       CPG23280
      IF(ILAST.GE.N3NXX) GO TO 501                                      CPG23290
      GO TO 500                                                         CPG23300
C                                                                       CPG23310
  501 CONTINUE                                                          CPG23320
      CALL SREW(ITAP44)                                                 CPG23330
      CALL SREW(JTAPE1)                                                 CPG23340
      RETURN                                                            CPG23350
      END                                                               CPG23360
      SUBROUTINE MAKEAB(B0,DP,DQ,DM,B1,NABC,ZETA,NIJ,ESO,U,T,LBLI,BUFI) CPG23370
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG23380
      DIMENSION B0(NIND,NABC),B1(NIND,NABC)                             CPG23390
      DIMENSION ZETA(NTRI,NTYPEP),NIJ(NTRI,2)                           CPG23400
      DIMENSION DP(NTRI,NABC,NTYPEP),DQ(NTRI,NABC,NTYPEP)               CPG23410
      DIMENSION DM(NTRI,NABC,NTYPEP)                                    CPG23420
      DIMENSION ESO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)         CPG23430
      DIMENSION LBLI(MAXBF4),BUFI(MAXBF2)                               CPG23440
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG23450
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3          CPG23460
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPG23470
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG23480
      COMMON/START/FOCC(10),NSORB(10),NSTR(10),NEND(10),MOTYP(256)      CPG23490
      COMMON/CI102/NIND,NDEP                                            CPG23500
      DATA ZERO / 0.0D+00 /                                             CPG23510
    1 FORMAT(//,2X,' B1 MATRIX'/)                                       CPG23520
    2 FORMAT(//,2X,' FINAL B1 MATRIX'/)                                 CPG23530
C                                                                       CPG23540
C   SEE EQ.(21)                                                         CPG23550
C   FORM DENSITY LIKE MATRIX IN SO BASIS                                CPG23560
C                                                                       CPG23570
      CALL PQDMAT(DP,DQ,B0,NIJ,ESO,NABC)                                CPG23580
C                                                                       CPG23590
C   FORM HALF-TRANSFORMED DENSITY LIKE MATRIX                           CPG23600
      CALL PQMAT(DP,DQ,DM,NTYPES*NABC,LBLI,BUFI)                        CPG23610
C                                                                       CPG23620
C   COMPLETE TRANSFORMATION                                             CPG23630
      DO 101 IABC=1,NABC                                                CPG23640
      DO 101 ITYP=1,NTYPES                                              CPG23650
      CALL MOCONS(DM(1,IABC,ITYP),NTRI,DP(1,IABC,ITYP),NTRI,ESO,U,T)    CPG23660
  101 CONTINUE                                                          CPG23670
C                                                                       CPG23680
      DO 103 II=1,NIND                                                  CPG23690
      I=NIJ(II,1)                                                       CPG23700
      J=NIJ(II,2)                                                       CPG23710
      ITYP=MOTYP(I)                                                     CPG23720
      JTYP=MOTYP(J)                                                     CPG23730
      IJ=IOFF(I)+J                                                      CPG23740
      DO 102 IABC=1,NABC                                                CPG23750
  102 B1(II,IABC)=DP(IJ,IABC,ITYP)-DP(IJ,IABC,JTYP)                     CPG23760
  103 CONTINUE                                                          CPG23770
C                                                                       CPG23780
      IF(IPRNT.LE.4) GO TO 301                                          CPG23790
      WRITE(6,1)                                                        CPG23800
      CALL MATOUT(B1,NIND,NABC,NIND,NABC,6)                             CPG23810
C                                                                       CPG23820
C   SEE EQ.(12)                                                         CPG23830
C   REST OF EQ.(11)                                                     CPG23840
C                                                                       CPG23850
  301 CONTINUE                                                          CPG23860
      DO 106 II=1,NIND                                                  CPG23870
      I=NIJ(II,1)                                                       CPG23880
      J=NIJ(II,2)                                                       CPG23890
      ITYP=MOTYP(I)                                                     CPG23900
      JTYP=MOTYP(J)                                                     CPG23910
      DO 105 JJ=1,NIND                                                  CPG23920
      K=NIJ(JJ,1)                                                       CPG23930
      L=NIJ(JJ,2)                                                       CPG23940
      KTYP=MOTYP(K)                                                     CPG23950
      LTYP=MOTYP(L)                                                     CPG23960
      ZVAL1=ZERO                                                        CPG23970
      ZVAL2=ZERO                                                        CPG23980
      ZVAL3=ZERO                                                        CPG23990
      ZVAL4=ZERO                                                        CPG24000
      IF(J.NE.K) GO TO 201                                              CPG24010
      IL=IOFF(MAX0(I,L))+MIN0(I,L)                                      CPG24020
      ZVAL1=ZETA(IL,ITYP)-ZETA(IL,JTYP)                                 CPG24030
  201 IF(I.NE.K) GO TO 202                                              CPG24040
      JL=IOFF(MAX0(J,L))+MIN0(J,L)                                      CPG24050
      ZVAL2=-ZETA(JL,JTYP)+ZETA(JL,ITYP)                                CPG24060
  202 IF(J.NE.L) GO TO 203                                              CPG24070
      IK=IOFF(MAX0(I,K))+MIN0(I,K)                                      CPG24080
      ZVAL3=-ZETA(IK,ITYP)+ZETA(IK,JTYP)                                CPG24090
  203 IF(I.NE.L) GO TO 204                                              CPG24100
      JK=IOFF(MAX0(J,K))+MIN0(J,K)                                      CPG24110
      ZVAL4=ZETA(JK,JTYP)-ZETA(JK,ITYP)                                 CPG24120
  204 ZVALT=ZVAL1+ZVAL2+ZVAL3+ZVAL4                                     CPG24130
C                                                                       CPG24140
      DO 104 IABC=1,NABC                                                CPG24150
      B1(II,IABC)=B1(II,IABC)+B0(JJ,IABC)*ZVALT                         CPG24160
  104 CONTINUE                                                          CPG24170
  105 CONTINUE                                                          CPG24180
  106 CONTINUE                                                          CPG24190
C                                                                       CPG24200
      IF(IPRNT.LE.4) GO TO 302                                          CPG24210
      WRITE(6,2)                                                        CPG24220
      CALL MATOUT(B1,NIND,NABC,NIND,NABC,6)                             CPG24230
C                                                                       CPG24240
  302 CONTINUE                                                          CPG24250
      RETURN                                                            CPG24260
      END                                                               CPG24270
      SUBROUTINE PQDMAT(DP,DQ,B0,NIJ,ESO,NABC)                          CPG24280
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG24290
      DIMENSION DP(NTRI,NABC,NTYPEP),DQ(NTRI,NABC,NTYPEP)               CPG24300
      DIMENSION B0(NIND,NABC),NIJ(NTRI,2),ESO(NBFAO,NBASIS)             CPG24310
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG24320
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3          CPG24330
      COMMON/GVBAB/A(10,10),B(10,10)                                    CPG24340
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG24350
      COMMON/START/FOCC(10),NSORB(10),NSTR(10),NEND(10),MOTYP(256)      CPG24360
      COMMON/CI102/NIND,NDEP                                            CPG24370
      DATA ZERO / 0.0D+00 /                                             CPG24380
    1 FORMAT(//,2X,' DP MATRIX'/)                                       CPG24390
    2 FORMAT(//,2X,' DQ MATRIX'/)                                       CPG24400
C                                                                       CPG24410
C   SEE EQS.(22) AND (23)                                               CPG24420
C                                                                       CPG24430
cyy   write(6,*) '  nabc in pqdmat = ', nabc
cyy   write(6,*) '  B0 matrix in pqdmat'
cyy   call matout(B0,nind,nabc,nind,nabc,6)
      DO 101 IABC=1,NABC                                                CPG24440
      DO 101 ITYP=1,NTYPEP                                              CPG24450
      DO 101 I=1,NTRI                                                   CPG24460
      DP(I,IABC,ITYP)=ZERO                                              CPG24470
      DQ(I,IABC,ITYP)=ZERO                                              CPG24480
  101 CONTINUE                                                          CPG24490
C                                                                       CPG24500
      IJ=0                                                              CPG24510
      DO 106 I=1,NBASIS                                                 CPG24520
      DO 106 J=1,I                                                      CPG24530
      IJ=IJ+1                                                           CPG24540
      DO 105 II=1,NIND                                                  CPG24550
      K=NIJ(II,1)                                                       CPG24560
      L=NIJ(II,2)                                                       CPG24570
      KTYP=MOTYP(K)                                                     CPG24580
      LTYP=MOTYP(L)                                                     CPG24590
      FAC=ESO(I,K)*ESO(J,L)+ESO(I,L)*ESO(J,K)                           CPG24600
      DO 104 ITYP=1,NTYPES                                              CPG24610
      AVAL=A(ITYP,KTYP)-A(ITYP,LTYP)                                    CPG24620
      BVAL=B(ITYP,KTYP)-B(ITYP,LTYP)                                    CPG24630
      FACA=AVAL*FAC                                                     CPG24640
      FACB=BVAL*FAC                                                     CPG24650
      DO 103 IABC=1,NABC                                                CPG24660
      DP(IJ,IABC,ITYP)=DP(IJ,IABC,ITYP)+B0(II,IABC)*FACA                CPG24670
      DQ(IJ,IABC,ITYP)=DQ(IJ,IABC,ITYP)+B0(II,IABC)*FACB                CPG24680
  103 CONTINUE                                                          CPG24690
  104 CONTINUE                                                          CPG24700
  105 CONTINUE                                                          CPG24710
  106 CONTINUE                                                          CPG24720
  110 CONTINUE                                                          CPG24730
      IF(IPRNT.LE.7) GO TO 201                                          CPG24740
      WRITE(6,1)                                                        CPG24750
      CALL MATOUT(DP,NTRI,NTYPES*NABC,NTRI,NTYPES*NABC,6)               CPG24760
      WRITE(6,2)                                                        CPG24770
      CALL MATOUT(DQ,NTRI,NTYPES*NABC,NTRI,NTYPES*NABC,6)               CPG24780
C                                                                       CPG24790
  201 CONTINUE                                                          CPG24800
      RETURN                                                            CPG24810
      END                                                               CPG24820
      SUBROUTINE UCMAT(EIG,OCC,ESO,NIJ,KIJ,U,T,B,DP,DQ,DM,DEL,          CPG24830
     1                 LBLI,BUFI,NABC)                                  CPG24840
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG24850
      DIMENSION EIG(NBASIS),OCC(NBASIS)                                 CPG24860
      DIMENSION ESO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)         CPG24870
      DIMENSION NIJ(NTRI,2),KIJ(NTRI,2),B(NTRI)                         CPG24880
      DIMENSION DP(NTRI,NABC),DQ(NTRI,NABC),DM(NTRI,NABC)               CPG24890
      DIMENSION DEL(NTRI),LBLI(MAXBF4),BUFI(MAXBF2)                     CPG24900
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG24910
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3          CPG24920
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPG24930
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPG24940
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG24950
      COMMON/TAPES/ITAP36,ITAP42,ITAP43,ITAP44                          CPG24960
      COMMON/CI102/NIND,NDEP                                            CPG24970
      DATA DLIMIT / 1.0D-8 /                                            CPG24980
      DATA ZERO,ONE / 0.0D+00 , 1.0D+00 /                               CPG24990
    1 FORMAT(//,2X,' B0 MATRIX IN UCMAT'/)
    2 FORMAT(//,2X,' DP MATRIX (SO BASIS) IN UCMAT'/) 
    3 FORMAT(//,2X,' DQ MATRIX (SO BASIS) IN UCMAT'/) 
    4 FORMAT(//,2X,' DM MATRIX (HALF-MO BASIS) IN UCMAT'/) 
    5 FORMAT(//,2X,' DP MATRIX (MO BASIS) IN UCMAT'/) 
C                                                                       CPG25030
C   CALL BACK B0 MATRICES                                               CPG25040
      CALL SREW(ITAP44)                                                 CPG25050
      DO 101 IABC=1,NABC                                                CPG25060
      IB=IBA(IABC)                                                      CPG25070
      CALL RREAD(ITAP44,DM(1,IABC),NTRI2,IB)                            CPG25080
  101 CONTINUE                                                          CPG25090
      IF(IPRNT.LE.4) GO TO 301                                          CPG25100
      WRITE(6,1)                                                        CPG25110
      CALL MATOUT(B,NTRI,NABC,NTRI,NABC,6)                              CPG25120
C                                                                       CPG25130
C   FORM DENSITY LIKE MATRIX IN SO BASIS                                CPG25140
  301 CONTINUE                                                          CPG25150
      CALL PQEMAT(DP,DQ,DM,NIJ,OCC,ESO,NABC)                            CPG25160
      IF(IPRNT.LE.4) GO TO 302                                          CPG25170
      WRITE(6,2)                                                        CPG25180
      CALL MATOUT(DP,NTRI,NABC,NTRI,NABC,6)                             CPG25190
      WRITE(6,3)                                                        CPG25180
      CALL MATOUT(DQ,NTRI,NABC,NTRI,NABC,6)                             CPG25190
C                                                                       CPG25200
C   FORM HALF-TRANSFORMED DENSITY LIKE MATRIX                           CPG25210
C   AND TRANSFORM THEM INTO MO BASIS                                    CPG25220
  302 CONTINUE                                                          CPG25230
      CALL PQMAT(DP,DQ,DM,NABC,LBLI,BUFI)                               CPG25240
      IF(IPRNT.LE.4) GO TO 303                                          CPG25280
      WRITE(6,4)                                                        CPG25290
      CALL MATOUT(DM,NTRI,NABC,NTRI,NABC,6)                             CPG25300
  303 continue
      DO 102 IABC=1,NABC                                                CPG25250
      CALL MOCONS(DM(1,IABC),NTRI,DP(1,IABC),NTRI,ESO,U,T)              CPG25260
  102 CONTINUE                                                          CPG25270
      IF(IPRNT.LE.4) GO TO 304                                          CPG25280
      WRITE(6,5)                                                        CPG25290
      CALL MATOUT(DP,NTRI,NABC,NTRI,NABC,6)                             CPG25300
C                                                                       CPG25310
  304 CONTINUE                                                          CPG25320
      DO 103 II=1,NDEP                                                  CPG25330
      I=KIJ(II,1)                                                       CPG25340
      J=KIJ(II,2)                                                       CPG25350
      DEL(II)=ZERO                                                      CPG25360
      DELJI=EIG(J)-EIG(I)                                               CPG25370
      IF(I.EQ.J) DELJI=ONE                                              CPG25380
      IF(DABS(DELJI).LE.DLIMIT) GO TO 103                               CPG25390
      DEL(II)=ONE/DELJI                                                 CPG25400
  103 CONTINUE                                                          CPG25410
C                                                                       CPG25420
      DO 105 IABC=1,NABC                                                CPG25430
      IB=IBA(IABC)                                                      CPG25440
      CALL RREAD(ITAP44,B,NTRI2,IB)                                     CPG25450
      DO 104 II=1,NDEP                                                  CPG25460
      I=KIJ(II,1)                                                       CPG25470
      J=KIJ(II,2)                                                       CPG25480
      IJ=IOFF(I)+J                                                      CPG25490
      TDEL=DEL(II)                                                      CPG25500
      IF(TDEL.EQ.ZERO) GO TO 104
      B(IJ)=(B(IJ)+DP(IJ,IABC))*TDEL                                    CPG25520
  104 CONTINUE                                                          CPG25530
      call rwrit(ITAP44,B,NTRI2,IB)  
cyy   write(6,*) ' B matrix in ucmat, iabc = ',iabc
cyy   call print(b,ntri,nbasis,6)
  105 CONTINUE                                                          CPG25540
C                                                                       CPG25550
      CALL SREW(ITAP44)                                                 CPG25560
      RETURN                                                            CPG25570
      END                                                               CPG25580
      SUBROUTINE PQEMAT(DP,DQ,B,NIJ,OCC,ESO,NABC)                       CPG25590
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG25600
      DIMENSION DP(NTRI,NABC),DQ(NTRI,NABC),B(NTRI,NABC)                CPG25610
      DIMENSION OCC(NBASIS),ESO(NBFAO,NBASIS),NIJ(NTRI,2)               CPG25620
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG25630
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3          CPG25640
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG25650
      COMMON/CI102/NIND,NDEP                                            CPG25660
      DATA ZERO,HALF / 0.0D+00 , 0.5D+00 /                              CPG25670
    1 FORMAT(//,2X,' DP MATRIX IN PQEMAT'/)
    2 FORMAT(//,2X,' DQ MATRIX IN PQEMAT'/)
C                                                                       CPG25700
C   FORM DENSITY MATRICES                                               CPG25710
C                                                                       CPG25720
C   FORM DP MATRIX                                                      CPG25730
      DO 101 IABC=1,NABC                                                CPG25740
      DO 101 I=1,NTRI                                                   CPG25750
      DP(I,IABC)=ZERO                                                   CPG25760
      DQ(I,IABC)=ZERO                                                   CPG25770
  101 CONTINUE                                                          CPG25780
C                                                                       CPG25790
      IJ=0                                                              CPG25800
      DO 107 I=1,NBASIS                                                 CPG25810
      DO 107 J=1,I                                                      CPG25820
      IJ=IJ+1                                                           CPG25830
      DO 106 II=1,NIND                                                  CPG25840
      K=NIJ(II,1)                                                       CPG25850
      L=NIJ(II,2)                                                       CPG25860
      KL=IOFF(K)+L                                                      CPG25870
      FAC=(ESO(I,K)*ESO(J,L)+ESO(I,L)*ESO(J,K))*(OCC(L)-OCC(K))*HALF    CPG25880
      DO 105 IABC=1,NABC                                                CPG25890
      FACT=B(KL,IABC)*FAC                                               CPG25900
      FAC2=FACT+FACT                                                    CPG25910
      DP(IJ,IABC)=DP(IJ,IABC)+FAC2                                      CPG25920
      DQ(IJ,IABC)=DQ(IJ,IABC)-FACT                                      CPG25930
  105 CONTINUE                                                          CPG25940
  106 CONTINUE                                                          CPG25950
  107 CONTINUE                                                          CPG25960
C                                                                       CPG25970
      IF(IPRNT.LE.4) GO TO 201                                          CPG25980
      WRITE(6,1)                                                        CPG25990
      CALL MATOUT(DP,NTRI,NABC,NTRI,NABC,6)                             CPG26000
      WRITE(6,2)                                                        CPG26010
      CALL MATOUT(DQ,NTRI,NABC,NTRI,NABC,6)                             CPG26020
C                                                                       CPG26030
  201 CONTINUE                                                          CPG26040
      RETURN                                                            CPG26050
      END                                                               CPG26060
      SUBROUTINE ORTHOG(U,W,NDIM,MDIM,N,NTEMP,NVEC)                     CPG26070
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG26080
      DIMENSION NDUM(360)                                               CPG26090
      DIMENSION U(NDIM,MDIM),W(NDIM)                                    CPG26100
      DATA ZERO / 0.0D+00 /                                             CPG26110
C                                                                       CPG26120
      DO 101 I=1,NTEMP                                                  CPG26130
  101 NDUM(I)=0                                                         CPG26140
C   N=1                                                                 CPG26150
      DO 102 I=1,N                                                      CPG26160
  102 W(I)=U(I,1)                                                       CPG26170
      CALL NORMLZ(W,NDIM,N,IFLAG)                                       CPG26180
      DO 103 I=1,N                                                      CPG26190
  103 U(I,1)=W(I)                                                       CPG26200
C                                                                       CPG26210
C   N=>2                                                                CPG26220
C                                                                       CPG26230
      NVEC=1                                                            CPG26240
      IF(NTEMP.EQ.1) RETURN                                             CPG26250
      DO 110 II=2,NTEMP                                                 CPG26260
      DO 104 I=1,N                                                      CPG26270
  104 W(I)=U(I,II)                                                      CPG26280
      DO 107 K=1,II-1                                                   CPG26290
      IF(NDUM(K).NE.0) GO TO 107                                        CPG26300
      EU=ZERO                                                           CPG26310
      DO 105 I=1,N                                                      CPG26320
  105 EU=EU+U(I,K)*U(I,II)                                              CPG26330
      DO 106 I=1,N                                                      CPG26340
  106 W(I)=W(I)-EU*U(I,K)                                               CPG26350
  107 CONTINUE                                                          CPG26360
      CALL NORMLZ(W,NDIM,N,IFLAG)                                       CPG26370
      IF(IFLAG.NE.0) NDUM(II)=1                                         CPG26380
      IF(IFLAG.NE.0) GO TO 110                                          CPG26390
      NVEC=NVEC+1                                                       CPG26400
      DO 108 I=1,N                                                      CPG26410
  108 U(I,II)=W(I)                                                      CPG26420
  110 CONTINUE                                                          CPG26430
C                                                                       CPG26440
      IX=0                                                              CPG26450
      DO 120 II=1,NTEMP                                                 CPG26460
      IF(NDUM(II).NE.0) GO TO 120                                       CPG26470
      IX=IX+1                                                           CPG26480
      DO 115 I=1,N                                                      CPG26490
  115 W(I)=U(I,II)                                                      CPG26500
      DO 116 I=1,N                                                      CPG26510
  116 U(I,IX)=W(I)                                                      CPG26520
  120 CONTINUE                                                          CPG26530
c************
c     write(6,*) '  nvec in orthog= ',nvec
c     write(6,*) '  U matrix after 120 in orthog'
c     CALL MATOUT(U,ndim,mdim,n,ntemp,6)
c************
C                                                                       CPG26540
      RETURN                                                            CPG26550
      END                                                               CPG26560
      SUBROUTINE NORMLZ(V,NDIM,N,IFLAG)                                 CPG26570
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG26580
      DIMENSION V(NDIM)                                                 CPG26590
    1 FORMAT(//,2X,'***ERROR IN NORMALIZATION***'/)                     CPG26600
      DATA RLIM / 1.0D-38 /                                             CPG26610
      DATA VLIM / 1.0D-6 /                                              CPG26620
      DATA ZERO / 0.0D+00 /                                             CPG26630
C                                                                       CPG26640
      IFLAG=0                                                           CPG26650
      SQL=ZERO                                                          CPG26660
      DO 101 I=1,N                                                      CPG26670
  101 SQL=SQL+V(I)*V(I)                                                 CPG26680
      R=DSQRT(SQL)                                                      CPG26690
      IF(R.LE.RLIM) GO TO 202                                           CPG26700
      IF(R.LE.VLIM) GO TO 204                                           CPG26710
C                                                                       CPG26720
  201 DO 102 I=1,N                                                      CPG26730
  102 V(I)=V(I)/R                                                       CPG26740
      GO TO 205                                                         CPG26750
  202 WRITE(6,1)                                                        CPG26760
      STOP                                                              CPG26770
C                                                                       CPG26780
  204 IFLAG=1                                                           CPG26790
  205 RETURN                                                            CPG26800
      END                                                               CPG26810
      SUBROUTINE NORM(XAB,A,B,N)                                        CPG26820
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG26830
      DIMENSION A(1),B(1)                                               CPG26840
      DATA ZERO / 0.0D+00 /                                             CPG26850
C                                                                       CPG26860
      XAB=ZERO                                                          CPG26870
      DO 101 I=1,N                                                      CPG26880
  101 XAB=XAB+A(I)*B(I)                                                 CPG26890
      RETURN                                                            CPG26900
      END                                                               CPG26910
      SUBROUTINE UXMAT(B,SM,U,NABC)                                     CPG26920
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG26930
      DIMENSION B(NTRI,NABC),SM(NTRI,NABC),U(NBASIS,NBASIS)             CPG26940
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG26950
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3          CPG26960
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPG26970
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG26980
      COMMON/TAPES/ITAP36,ITAP42,ITAP43,ITAP44                          CPG26990
      DATA ZERO,HALF / 0.0D+00 , 0.5D+00 /                              CPG27000
    1 FORMAT(//,2X,' UA MATRIX IN UXMAT, IABC= ',I5/)
C                                                                       CPG27020
      CALL SREW(ITAP44)                                                 CPG27030
      DO 101 IABC=1,NABC                                                CPG27040
      IB=IBA(IABC)                                                      CPG27050
      CALL RREAD(ITAP44,B(1,IABC),NTRI2,IB)                             CPG27060
      IF(IABC.GT.N3N) GO TO 101                                         CPG27070
      IS=ISA(IABC)                                                      CPG27080
      CALL RREAD(ITAP44,SM(1,IABC),NTRI2,IS)                            CPG27090
  101 CONTINUE                                                          CPG27100
C                                                                       CPG27110
C   FILL UP THE REST OF U MATRIX                                        CPG27120
      DO 105 IABC=1,NABC                                                CPG27130
      IU=IUA(IABC)                                                      CPG27140
      IF(IABC.GT.N3N) GO TO 301                                         CPG27150
C                                                                       CPG27160
      IJ=0                                                              CPG27170
      DO 102 I=1,NBASIS                                                 CPG27180
      DO 102 J=1,I                                                      CPG27190
      IJ=IJ+1                                                           CPG27200
      IF(I.EQ.J) GO TO 201                                              CPG27210
      U(I,J)=B(IJ,IABC)                                                 CPG27220
      U(J,I)=-U(I,J)-SM(IJ,IABC)                                        CPG27230
      GO TO 102                                                         CPG27240
  201 U(I,I)=-SM(IJ,IABC)*HALF                                          CPG27250
  102 CONTINUE                                                          CPG27260
      GO TO 302                                                         CPG27270
C                                                                       CPG27280
  301 CONTINUE                                                          CPG27290
      IJ=0                                                              CPG27300
      DO 103 I=1,NBASIS                                                 CPG27310
      DO 103 J=1,I                                                      CPG27320
      IJ=IJ+1                                                           CPG27330
      IF(I.EQ.J) GO TO 202                                              CPG27340
      U(I,J)=B(IJ,IABC)                                                 CPG27350
      U(J,I)=-U(I,J)                                                    CPG27360
      GO TO 103                                                         CPG27370
  202 U(I,I)=ZERO                                                       CPG27380
  103 CONTINUE                                                          CPG27390
C                                                                       CPG27400
  302 CONTINUE                                                          CPG27410
      CALL RWRIT(ITAP44,U,NBASQ,IU)                                     CPG27420
      IF(IPRNT.LE.2) GO TO 105                                          CPG27430
      WRITE(6,1) IABC                                                   CPG27440
      CALL MATOUT(U,NBASIS,NBASIS,NBASIS,NBASIS,6)                      CPG27450
C                                                                       CPG27460
  105 CONTINUE                                                          CPG27470
C                                                                       CPG27480
      CALL SREW(ITAP44)                                                 CPG27490
      RETURN                                                            CPG27500
      END                                                               CPG27510
      SUBROUTINE WAMAT(ESO,U,T,EPA,ZETA,WA,SA,DP,DQ,DM,LBLI,BUFI)       CPG27520
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG27530
      DIMENSION ESO(NBFAO,NBASIS),U(NBASIS,NBASIS),T(NBASIS,NBASIS)     CPG27540
      DIMENSION EPA(NBASIS,NBASIS)                                      CPG27550
      DIMENSION ZETA(NTRI,NTYPEP),WA(NBASIS,NBASIS),SA(NTRI,N3N*NTYPEP) CPG27560
      DIMENSION DP(NTRI,N3N,NTYPEP),DQ(NTRI,N3N,NTYPEP)                 CPG27570
      DIMENSION DM(NTRI,N3N,NTYPEP)                                     CPG27580
      DIMENSION LBLI(MAXBF4),BUFI(MAXBF2)                               CPG27590
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG27600
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3          CPG27610
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPG27620
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPG27630
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG27640
      COMMON/START/FOCC(10),NSORB(10),NSTR(10),NEND(10),MOTYP(256)      CPG27650
      COMMON/TAPES/ITAP36,ITAP42,ITAP43,ITAP44                          CPG27660
      COMMON/WORKS/JTAPE1                                               CPG27670
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPG27680
      DATA ZERO,TWO / 0.0D+00 , 2.0D+00 /                               CPG27690
    1 FORMAT(//,2X,' DP MATRIX (SO BASIS) IN WAMAT'/)
    2 FORMAT(//,2X,' DQ MATRIX (SO BASIS) IN WAMAT'/)
    3 FORMAT(//,2X,' DM MATRIX (HALF-MO BASIS) IN WAMAT'/)
    4 FORMAT(//,2X,' DP MATRIX (MO BASIS) IN WAMAT'/)
    5 FORMAT(//,2X,' WA MATRIX IN WAMAT, IABC = ',I5/)

C                                                                       CPG27710
C   FORM WA MATRICES                                                    CPG27720
C   SEE EQ.(6)                                                          CPG27730
C                                                                       CPG27740
C   CALL BACK LAGRANGIAN MATRICES                                       CPG27750
      DO 101 ITYP=1,NTYPES                                              CPG27760
      CALL MREAD(ZETA(1,ITYP),26+ITYP)                                  CPG27770
  101 CONTINUE                                                          CPG27780
      DO 102 I=1,NTRI                                                   CPG27790
  102 ZETA(I,NTYPEP)=ZERO                                               CPG27800
C                                                                       CPG27810
C   CALL BACK DERIVATIVE OVERLAP INTEGRALS                              CPG27820
      CALL SREW(ITAP44)                                                 CPG27830
      DO 103 IABC=1,N3N                                                 CPG27840
      IS=ISA(IABC)                                                      CPG27850
      CALL RREAD(ITAP44,SA(1,IABC),NTRI2,IS)                            CPG27860
  103 CONTINUE                                                          CPG27870
C                                                                       CPG27880
C   FORM DENSITY LIKE MATRIX                                            CPG27890
      CALL PQFMAT(DP,DQ,SA,ESO,N3N)                                     CPG27900
      IF(IPRNT.LE.4) GO TO 301                                          CPG25170
      WRITE(6,1)                                                        CPG25180
      CALL MATOUT(DP,NTRI,N3N*NTYPES,NTRI,N3N*NTYPES,6) 
      WRITE(6,2)                                                        CPG25180
      CALL MATOUT(DQ,NTRI,N3N*NTYPES,NTRI,N3N*NTYPES,6) 
C                                                                       CPG27910
C   FORM HALF-TRANSFORMED DENSITY LIKE MATRIX                           CPG27920
  301 continue
      CALL PQMAT(DP,DQ,DM,NTYPES*N3N,LBLI,BUFI)                         CPG27930
      IF(IPRNT.LE.4) GO TO 302                                          CPG25170
      WRITE(6,3)                                                        CPG25180
      CALL MATOUT(DM,NTRI,N3N*NTYPES,NTRI,N3N*NTYPES,6) 
C                                                                       CPG27940
  302 continue
      DO 104 IABC=1,N3N                                                 CPG27950
      DO 104 ITYP=1,NTYPES                                              CPG27960
      CALL MOCONS(DM(1,IABC,ITYP),NTRI,DP(1,IABC,ITYP),NTRI,ESO,U,T)    CPG27970
  104 CONTINUE                                                          CPG27980
cyy   WRITE(6,4)                                                        CPG25180
cyy   CALL MATOUT(DP,NTRI,N3N*NTYPES,NTRI,N3N*NTYPES,6) 
C                                                                       CPG27990
C   COMPLETE WA MATRICES                                                CPG28000
      CALL SREW(JTAPE1)                                                 CPG28010
      DO 110 IABC=1,N3N                                                 CPG28020
      IS=ISA(IABC)                                                      CPG28030
      IE=IEA(IABC)                                                      CPG28040
      CALL RREAD(ITAP44,SA(1,IABC),NTRI2,IS)                            CPG28050
      CALL RREAD(ITAP44,EPA,NBASQ,IE)                                   CPG28060
C                                                                       CPG28070
      DO 107 I=1,NBASIS                                                 CPG28080
      ITYP=MOTYP(I)                                                     CPG28090
      DO 107 J=1,NBASIS                                                 CPG28100
      JTYP=MOTYP(J)                                                     CPG28110
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPG28120
      VALU1=EPA(J,I)+EPA(J,I)                                           CPG28130
      VALU2=DP(IJ,IABC,JTYP)                                            CPG28140
      VALU3=ZERO                                                        CPG28150
      DO 106 K=1,JOCC                                                   CPG28160
      IK=IOFF(MAX0(I,K))+MIN0(I,K)                                      CPG28170
      JK=IOFF(MAX0(J,K))+MIN0(J,K)                                      CPG28180
      FAC1=ZETA(IK,ITYP)+ZETA(IK,JTYP)                                  CPG28190
      FAC2=ZETA(JK,JTYP)+ZETA(JK,JTYP)                                  CPG28200
      VALU3=VALU3+SA(JK,IABC)*FAC1+SA(IK,IABC)*FAC2                     CPG28210
  106 CONTINUE                                                          CPG28220
      WA(I,J)=VALU1-VALU2-VALU3                                         CPG28230
  107 CONTINUE                                                          CPG28240
      CALL SWRIT(JTAPE1,WA,NBASQ)                                       CPG28250
      IF(IPRNT.LE.2) GO TO 110                                          CPG28260
      WRITE(6,5) IABC                                                   CPG28270
      CALL MATOUT(WA,NBASIS,NBASIS,NBASIS,NBASIS,6)                     CPG28280
  110 CONTINUE                                                          CPG28290
C                                                                       CPG28300
      CALL SREW(ITAP44)                                                 CPG28310
      CALL SREW(JTAPE1)                                                 CPG28320
      RETURN                                                            CPG28330
      END                                                               CPG28340
      SUBROUTINE PQFMAT(DP,DQ,SA,ESO,NABC)                              CPG28350
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG28360
      DIMENSION DP(NTRI,NABC,NTYPEP),DQ(NTRI,NABC,NTYPEP)               CPG28370
      DIMENSION ESO(NBFAO,NBASIS),SA(NTRI,NABC*NTYPEP)                  CPG28380
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG28390
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3          CPG28400
      COMMON/GVBAB/A(10,10),B(10,10)                                    CPG28410
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG28420
      COMMON/START/FOCC(10),NSORB(10),NSTR(10),NEND(10),MOTYP(256)      CPG28430
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPG28440
      DATA ZERO / 0.0D+00 /                                             CPG28450
    1 FORMAT(//,2X,' DP MATRIX IN PQFMAT'/)
    2 FORMAT(//,2X,' DQ MATRIX IN PQFMAT'/)
C                                                                       CPG28460
      DO 101 IABC=1,NABC                                                CPG28470
      DO 101 ITYP=1,NTYPEP                                              CPG28480
      DO 101 I=1,NTRI                                                   CPG28490
      DP(I,IABC,ITYP)=ZERO                                              CPG28500
      DQ(I,IABC,ITYP)=ZERO                                              CPG28510
  101 CONTINUE                                                          CPG28520
C                                                                       CPG28530
      IJ=0                                                              CPG28540
      DO 105 I=1,NBASIS                                                 CPG28550
      DO 105 J=1,I                                                      CPG28560
      IJ=IJ+1                                                           CPG28570
      DO 104 K=1,JOCC                                                   CPG28580
      KTYP=MOTYP(K)                                                     CPG28590
      DO 104 L=1,JOCC                                                   CPG28600
      LTYP=MOTYP(L)                                                     CPG28610
      KL=IOFF(MAX0(K,L))+MIN0(K,L)                                      CPG28620
      FAC=ESO(I,K)*ESO(J,L)+ESO(J,K)*ESO(I,L)                           CPG28630
      DO 103 ITYP=1,NTYPES                                              CPG28640
      AVAL=A(ITYP,LTYP)                                                 CPG28650
      BVAL=B(ITYP,LTYP)                                                 CPG28660
      FACA=AVAL*FAC                                                     CPG28670
      FACB=BVAL*FAC                                                     CPG28680
      DO 102 IABC=1,NABC                                                CPG28690
      DP(IJ,IABC,ITYP)=DP(IJ,IABC,ITYP)+SA(KL,IABC)*FACA                CPG28700
      DQ(IJ,IABC,ITYP)=DQ(IJ,IABC,ITYP)+SA(KL,IABC)*FACB                CPG28710
  102 CONTINUE                                                          CPG28720
  103 CONTINUE                                                          CPG28730
  104 CONTINUE                                                          CPG28740
  105 CONTINUE                                                          CPG28750
C                                                                       CPG28760
      IF(IPRNT.LE.4) GO TO 201                                          CPG25980
      WRITE(6,1) 
      CALL MATOUT(DP,NTRI,NABC*NYYPES,NTRI,NABC*NTYPES,6) 
      WRITE(6,2)   
      CALL MATOUT(DQ,NTRI,NABC*NYYPES,NTRI,NABC*NTYPES,6) 
c
  201 continue
      RETURN                                                            CPG28770
      END                                                               CPG28780
      SUBROUTINE SCF2ND(E2A,E2M,E2P,SA,WA,EPA,UA)                       CPG28790
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG28800
      DIMENSION E2A(N3N,N3N),E2M(N3N,N3N),E2P(N3N,N3N)                  CPG28810
      DIMENSION SA(NTRI,N3N),WA(NBASIS,NBASIS,N3N)                      CPG28820
      DIMENSION EPA(NBASIS,NBASIS),UA(NBASIS,NBASIS)                    CPG28830
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG28840
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3          CPG28850
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPG28860
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG28870
      COMMON/START/FOCC(10),NSORB(10),NSTR(10),NEND(10),MOTYP(256)      CPG28880
      COMMON/TAPES/ITAP36,ITAP42,ITAP43,ITAP44                          CPG28890
      COMMON/WORKS/JTAPE1                                               CPG28900
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPG28910
      DATA ZERO,TWO / 0.0D+00 , 2.0D+00 /                               CPG28920
    1 FORMAT(//,2X,' WA MATRIX, IABC = ',I5/)                           CPG28930
    2 FORMAT(//,2X,' E2P MATRIX'/)                                      CPG28940
    3 FORMAT(//,2X,' E2M MATRIX'/)                                      CPG28950
    4 FORMAT(///                                                        CPG28960
     1 2X,':::::::::::::::::::::::::::::::::::::::::::::'/              CPG28970
     2 2X,':::SUMMARY OF FIRST ORDER CPHF CALCULATION:::'/              CPG28980
     3 2X,':::::::::::::::::::::::::::::::::::::::::::::'/)             CPG28990
    5 FORMAT(//,2X,' E2A MATRIX'/)                                      CPG29000
    6 FORMAT(//,2X,' E2M MATRIX'/)                                      CPG29010
    7 FORMAT(//,2X,' SCF SECOND DERIVATIVE MATRIX'/)                    CPG29020
    8 FORMAT(2I5)                                                       CPG29030
    9 FORMAT(3F20.10)                                                   CPG29040
C                                                                       CPG29050
C   CALCULATE SCF SECOND DERIVATIVES                                    CPG29060
C   SEE EQ.(5)                                                          CPG29070
C                                                                       CPG29080
C   CALL BACK SA MATRICES                                               CPG29090
      CALL SREW(ITAP44)                                                 CPG29100
      CALL SREW(JTAPE1)                                                 CPG29110
      DO 101 IABC=1,N3N                                                 CPG29120
      IS=ISA(IABC)                                                      CPG29130
      CALL RREAD(ITAP44,SA(1,IABC),NTRI2,IS)                            CPG29140
      CALL SREAD(JTAPE1,WA(1,1,IABC),NBASQ)                             CPG29150
CCC   WRITE(6,1) IABC                                                   CPG29160
CCC   CALL MATOUT(WA(1,1,IABC),NBASIS,NBASIS,NBASIS,NBASIS,6)           CPG29170
  101 CONTINUE                                                          CPG29180
C                                                                       CPG29190
C   (B COORDINATE)                                                      CPG29200
      DO 110 IABC=1,N3N                                                 CPG29210
      IE=IEA(IABC)                                                      CPG29220
      IU=IUA(IABC)                                                      CPG29230
      CALL RREAD(ITAP44,EPA,NBASQ,IE)                                   CPG29240
      CALL RREAD(ITAP44,UA,NBASQ,IU)                                    CPG29250
C                                                                       CPG29260
C   (A COORDINATE)                                                      CPG29270
      DO 105 JABC=1,N3N                                                 CPG29280
      VALUP=ZERO                                                        CPG29290
      VALUM=ZERO                                                        CPG29300
C                                                                       CPG29310
      DO 104 I=1,JOCC                                                   CPG29320
      DO 102 J=1,JOCC                                                   CPG29330
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPG29340
      VALUP=VALUP-SA(IJ,JABC)*EPA(I,J)                                  CPG29350
  102 CONTINUE                                                          CPG29360
      DO 103 J=1,NBASIS                                                 CPG29370
      VALUM=VALUM+WA(J,I,JABC)*UA(J,I)                                  CPG29380
  103 CONTINUE                                                          CPG29390
  104 CONTINUE                                                          CPG29400
      E2M(IABC,JABC)=VALUM*TWO                                          CPG29410
      E2P(IABC,JABC)=VALUP*TWO                                          CPG29420
  105 CONTINUE                                                          CPG29430
  110 CONTINUE                                                          CPG29440
      IF(IPRNT.LE.3) GO TO 201                                          CPG29450
      WRITE(6,2)                                                        CPG29460
      CALL MATOUT(E2P,N3N,N3N,N3N,N3N,6)                                CPG29470
      WRITE(6,3)                                                        CPG29480
      CALL MATOUT(E2M,N3N,N3N,N3N,N3N,6)                                CPG29490
C                                                                       CPG29500
  201 CONTINUE                                                          CPG29510
      DO 115 I=1,N3N                                                    CPG29520
      DO 115 J=1,N3N                                                    CPG29530
  115 E2M(I,J)=E2M(I,J)+E2P(I,J)                                        CPG29540
      WRITE(6,4)                                                        CPG29550
      WRITE(6,5)                                                        CPG29560
      CALL MATOUT(E2A,N3N,N3N,N3N,N3N,6)                                CPG29570
      WRITE(6,6)                                                        CPG29580
      CALL MATOUT(E2M,N3N,N3N,N3N,N3N,6)                                CPG29590
      DO 116 I=1,N3N                                                    CPG29600
      DO 116 J=1,N3N                                                    CPG29610
  116 E2M(I,J)=E2A(I,J)+E2M(I,J)                                        CPG29620
      WRITE(6,7)                                                        CPG29630
      CALL MATOUT(E2M,N3N,N3N,N3N,N3N,6)                                CPG29640
C                                                                       CPG29650
      ITAP15=15                                                         CPG29660
      REWIND ITAP15                                                     CPG29670
      N6N=N3N*2                                                         CPG29680
      WRITE(ITAP15,8) NATOMS,N6N                                        CPG29690
      WRITE(ITAP15,9) ((E2M(I,J),J=1,N3N),I=1,N3N)                      CPG29700
      REWIND ITAP15                                                     CPG29710
C                                                                       CPG29720
      CALL SREW(ITAP44)                                                 CPG29730
      CALL SREW(JTAPE1)                                                 CPG29740
      RETURN                                                            CPG29750
      END                                                               CPG29760
      SUBROUTINE DIPMO(OCC,DIP,UA)                                      CPG29770
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG29780
      DIMENSION OCC(NBASIS),DIP(NTRI,3),UA(NBASIS,NBASIS,3)             CPG29790
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG29800
      COMMON/CALIF/LPARA(1024),DIPDE(3,45),DUM(889)                     CPG29810
      COMMON/DIPOL/DIPDA(3,45),DIPDM(3,45),DIPDN(3,45),DIPDT(3,45)      CPG29820
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3          CPG29830
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPG29840
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG29850
      COMMON/TAPES/ITAP36,ITAP42,ITAP43,ITAP44                          CPG29860
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPG29870
      DATA ZERO,HALF,ONE,FOUR / 0.0D+00 , 0.5D+00 , 1.0D+00 , 4.0D+00 / CPG29880
      DATA DEBYE,BOHR / 2.541765480D+00 , 0.52917706D+00 /              CPG29890
    1 FORMAT(//,2X,' FIRST DERIVATIVES OF DIPOLE MOMENTS'/)             CPG29900
    2 FORMAT(//,2X,' UA MATRIX, IABC = ',I5,' IXYZ = ',I5/)             CPG29910
    3 FORMAT(//,2X,' DIPDA MATRIX'/)                                    CPG29920
    4 FORMAT(//,2X,' DIPDM MATRIX'/)                                    CPG29930
    5 FORMAT(//,2X,' DIPDE MATRIX'/)                                    CPG29940
    6 FORMAT(//,2X,' DIPDN MATRIX'/)                                    CPG29950
    7 FORMAT(//,2X,' FIRST DERIVATIVES OF DIPOLE MOMENTS (IN DEBYE/BOHR)CPG29960
     *'/)                                                               CPG29970
    8 FORMAT(//,2X,' FIRST DERIVATIVES OF DIPOLE MOMENTS (IN DEBYE/A)'/)CPG29980
    9 FORMAT(2I5)                                                       CPG29990
   10 FORMAT(3F20.10)                                                   CPG30000
C                                                                       CPG30010
      CALL SREW(ITAP44)                                                 CPG30020
      DEB4=DEBYE*FOUR                                                   CPG30030
      RBOHR=ONE/BOHR                                                    CPG30040
C     WRITE(6,1)                                                        CPG30050
C                                                                       CPG30060
      DO 101 I=1,3                                                      CPG30070
      DO 101 J=1,N3N                                                    CPG30080
      DIPDM(I,J)=ZERO                                                   CPG30090
  101 CONTINUE                                                          CPG30100
C                                                                       CPG30110
      DO 102 IXYZ=1,3                                                   CPG30120
      IH=IHA(N3N+IXYZ)                                                  CPG30130
      CALL RREAD(ITAP44,DIP(1,IXYZ),NTRI2,IH)                           CPG30140
  102 CONTINUE                                                          CPG30150
C                                                                       CPG30160
      KABC=0                                                            CPG30170
      DO 105 IABC=1,NATOMS                                              CPG30180
      I00=(IABC-1)*3                                                    CPG30190
      IXX=I00+1                                                         CPG30200
      IYY=I00+2                                                         CPG30210
      IZZ=I00+3                                                         CPG30220
      DO 103 IXYZ=1,3                                                   CPG30230
      KABC=KABC+1                                                       CPG30240
      IU=IUA(KABC)                                                      CPG30250
      CALL RREAD(ITAP44,UA(1,1,IXYZ),NBASQ,IU)                          CPG30260
      IF(IPRNT.LE.3) GO TO 103                                          CPG30270
      WRITE(6,2) IABC,IXYZ                                              CPG30280
      CALL MATOUT(UA(1,1,IXYZ),NBASIS,NBASIS,NBASIS,NBASIS,6)           CPG30290
  103 CONTINUE                                                          CPG30300
C                                                                       CPG30310
      DO 104 I=1,JOCC                                                   CPG30320
      FAC=OCC(I)*HALF                                                   CPG30330
      DO 104 J=1,NBASIS                                                 CPG30340
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPG30350
      DIPDM(1,IXX)=DIPDM(1,IXX)+UA(J,I,1)*DIP(IJ,1)*FAC                 CPG30360
      DIPDM(1,IYY)=DIPDM(1,IYY)+UA(J,I,2)*DIP(IJ,1)*FAC                 CPG30370
      DIPDM(1,IZZ)=DIPDM(1,IZZ)+UA(J,I,3)*DIP(IJ,1)*FAC                 CPG30380
      DIPDM(2,IXX)=DIPDM(2,IXX)+UA(J,I,1)*DIP(IJ,2)*FAC                 CPG30390
      DIPDM(2,IYY)=DIPDM(2,IYY)+UA(J,I,2)*DIP(IJ,2)*FAC                 CPG30400
      DIPDM(2,IZZ)=DIPDM(2,IZZ)+UA(J,I,3)*DIP(IJ,2)*FAC                 CPG30410
      DIPDM(3,IXX)=DIPDM(3,IXX)+UA(J,I,1)*DIP(IJ,3)*FAC                 CPG30420
      DIPDM(3,IYY)=DIPDM(3,IYY)+UA(J,I,2)*DIP(IJ,3)*FAC                 CPG30430
      DIPDM(3,IZZ)=DIPDM(3,IZZ)+UA(J,I,3)*DIP(IJ,3)*FAC                 CPG30440
  104 CONTINUE                                                          CPG30450
  105 CONTINUE                                                          CPG30460
C                                                                       CPG30470
      DO 106 I=1,3                                                      CPG30480
      DO 106 J=1,N3N                                                    CPG30490
      DIPDM(I,J)=-DIPDM(I,J)*DEB4                                       CPG30500
  106 CONTINUE                                                          CPG30510
C                                                                       CPG30520
C   SUM UP CONTRIBUTIONS FROM MO CHANGE                                 CPG30530
      DO 107 I=1,3                                                      CPG30540
      DO 107 J=1,N3N                                                    CPG30550
      DIPDE(I,J)=DIPDA(I,J)+DIPDM(I,J)                                  CPG30560
  107 CONTINUE                                                          CPG30570
      IF(IPRNT.LE.2) GO TO 201                                          CPG30580
      WRITE(6,3)                                                        CPG30590
      CALL MATOUT(DIPDA,3,45,3,N3N,6)                                   CPG30600
      WRITE(6,4)                                                        CPG30610
      CALL MATOUT(DIPDM,3,45,3,N3N,6)                                   CPG30620
C                                                                       CPG30630
C   SUM UP ALL CONTRIBUTIONS                                            CPG30640
  201 CONTINUE                                                          CPG30650
      DO 108 I=1,3                                                      CPG30660
      DO 108 J=1,N3N                                                    CPG30670
      DIPDT(I,J)=DIPDE(I,J)+DIPDN(I,J)                                  CPG30680
  108 CONTINUE                                                          CPG30690
C                                                                       CPG30700
      IF(IPRNT.LE.2) GO TO 202                                          CPG30710
      WRITE(6,5)                                                        CPG30720
      CALL MATOUT(DIPDE,3,45,3,N3N,6)                                   CPG30730
      WRITE(6,6)                                                        CPG30740
      CALL MATOUT(DIPDN,3,45,3,N3N,6)                                   CPG30750
  202 CONTINUE                                                          CPG30760
      WRITE(6,7)                                                        CPG30770
      CALL MATOUT(DIPDT,3,45,3,N3N,6)                                   CPG30780
C                                                                       CPG30790
      DO 109 I=1,3                                                      CPG30800
      DO 109 J=1,N3N                                                    CPG30810
      DIPDT(I,J)=DIPDT(I,J)*RBOHR                                       CPG30820
  109 CONTINUE                                                          CPG30830
      WRITE(6,8)                                                        CPG30840
      CALL MATOUT(DIPDT,3,45,3,N3N,6)                                   CPG30850
C                                                                       CPG30860
      ITAP17=17                                                         CPG30870
      REWIND ITAP17                                                     CPG30880
      WRITE(ITAP17,9) NATOMS,N3N                                        CPG30890
      WRITE(ITAP17,10) ((DIPDT(I,J),J=1,N3N),I=1,3)                     CPG30900
      REWIND ITAP17                                                     CPG30910
C                                                                       CPG30920
      CALL SREW(ITAP44)                                                 CPG30930
      RETURN                                                            CPG30940
      END                                                               CPG30950
      SUBROUTINE POLAR(OCC,HF,UF)                                       CPG30960
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG30970
      DIMENSION OCC(NBASIS),HF(NTRI,3),UF(NBASIS,NBASIS)                CPG30980
      DIMENSION POL(3,3),POLC(3,3),EIG(3),EIV(3,3),A(6),FV1(3),FV2(3)   CPG30990
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG31000
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOMS,N3N,N3NP3,N3NX3          CPG31010
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPG31020
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPG31030
      COMMON/TAPES/ITAP36,ITAP42,ITAP43,ITAP44                          CPG31040
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPG31050
      DATA ZERO,HALF,FOUR / 0.0D+00 , 0.5D+00 , 4.0D+00 /               CPG31060
      DATA BOHR / 0.52917706D+00 /                                      CPG31070
    1 FORMAT(//,2X,' HF MATRIX, IXYZ = ',I5/)                           CPG31080
    2 FORMAT(//,2X,' UF MATRIX, IXYZ = ',I5/)                           CPG31090
    3 FORMAT(//,2X,' POLARIZABILITY TENSOR MATRIX (IN A.U.)'/)          CPG31100
    4 FORMAT(//,2X,' POLARIZABILITY TENSOR MATRIX (IN A**3)'/)          CPG31110
    5 FORMAT(//,2X,' A MATRIX'/)                                        CPG31120
    6 FORMAT(//,2X,' PRINCIPAL POLARIZABILITY TENSOR MATRIX (IN A.U.)'/)CPG31130
    7 FORMAT(//,2X,' PRINCIPAL POLARIZABILITY TENSOR MATRIX (IN A**3)'/)CPG31140
    8 FORMAT(//,2X,' EIGENVECTOR MATRIX'/)                              CPG31150
C                                                                       CPG31160
      CALL SREW(ITAP44)                                                 CPG31170
      A33=BOHR*BOHR*BOHR                                                CPG31180
C                                                                       CPG31190
C   CALL BACK HF MATRICES                                               CPG31200
      DO 101 IXYZ=1,3                                                   CPG31210
      IH=IHA(N3N+IXYZ)                                                  CPG31220
      CALL RREAD(ITAP44,HF(1,IXYZ),NTRI2,IH)                            CPG31230
      IF(IPRNT.LE.3) GO TO 101                                          CPG31240
      WRITE(6,1) IXYZ                                                   CPG31250
      CALL PRINT(HF(1,IXYZ),NTRI,NBASIS,6)                              CPG31260
  101 CONTINUE                                                          CPG31270
C                                                                       CPG31280
C:::B COORDINATE:::                                                     CPG31290
      DO 105 IXYZ=1,3                                                   CPG31300
      IU=IUA(N3N+IXYZ)                                                  CPG31310
      CALL RREAD(ITAP44,UF,NBASQ,IU)                                    CPG31320
      IF(IPRNT.LE.3) GO TO 301                                          CPG31330
      WRITE(6,2) IXYZ                                                   CPG31340
      CALL MATOUT(UF,NBASIS,NBASIS,NBASIS,NBASIS,6)                     CPG31350
C                                                                       CPG31360
C:::A COORDINATES:::                                                    CPG31370
  301 CONTINUE                                                          CPG31380
      DO 103 JXYZ=1,3                                                   CPG31390
      VALUP=ZERO                                                        CPG31400
      DO 102 I=1,JOCC                                                   CPG31410
      FAC=OCC(I)*HALF                                                   CPG31420
      DO 102 J=1,NBASIS                                                 CPG31430
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPG31440
      VALUP=VALUP-UF(J,I)*HF(IJ,JXYZ)*FAC                               CPG31450
  102 CONTINUE                                                          CPG31460
      POL(IXYZ,JXYZ)=VALUP*FOUR                                         CPG31470
  103 CONTINUE                                                          CPG31480
  105 CONTINUE                                                          CPG31490
      WRITE(6,3)                                                        CPG31500
      CALL MATOUT(POL,3,3,3,3,6)                                        CPG31510
      DO 106 I=1,3                                                      CPG31520
      DO 106 J=1,3                                                      CPG31530
      POLC(I,J)=POL(I,J)*A33                                            CPG31540
  106 CONTINUE                                                          CPG31550
      WRITE(6,4)                                                        CPG31560
      CALL MATOUT(POLC,3,3,3,3,6)                                       CPG31570
C                                                                       CPG31580
      IJ=0                                                              CPG31590
      DO 107 I=1,3                                                      CPG31600
      DO 107 J=1,I                                                      CPG31610
      IJ=IJ+1                                                           CPG31620
      A(IJ)=POL(I,J)                                                    CPG31630
  107 CONTINUE                                                          CPG31640
      WRITE(6,5)                                                        CPG31650
      CALL PRINT(A,6,3,6)                                               CPG31660
      CALL RSP(3,3,6,A,EIG,1,EIV,FV1,FV2)                               CPG31670
      DO 108 I=1,3                                                      CPG31680
      DO 108 J=1,3                                                      CPG31690
      POL(I,J)=ZERO                                                     CPG31700
  108 CONTINUE                                                          CPG31710
      DO 109 I=1,3                                                      CPG31720
  109 POL(I,I)=EIG(I)                                                   CPG31730
      WRITE(6,6)                                                        CPG31740
      CALL MATOUT(POL,3,3,3,3,6)                                        CPG31750
      DO 110 I=1,3                                                      CPG31760
      DO 110 J=1,3                                                      CPG31770
      POLC(I,J)=POL(I,J)*A33                                            CPG31780
  110 CONTINUE                                                          CPG31790
      WRITE(6,7)                                                        CPG31800
      CALL MATOUT(POLC,3,3,3,3,6)                                       CPG31810
      WRITE(6,8)                                                        CPG31820
      CALL MATOUT(EIV,3,3,3,3,6)                                        CPG31830
C                                                                       CPG31840
      CALL SREW(ITAP44)                                                 CPG31850
      RETURN                                                            CPG31860
      END                                                               CPG31870
      SUBROUTINE WRBMAT(BB,NABC)                                        CPG31880
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPG31890
      DIMENSION BB(NTRI,NABC)                                           CPG31900
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPG31910
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPG31920
      COMMON/TAPES/ITAP36,ITAP42,ITAP43,ITAP44                          CPG31930
    1 FORMAT(//,2X,' B0 MATRIX'/)                                       CPG31940
C                                                                       CPG31950
      CALL SREW(ITAP44)                                                 CPG31960
      DO 101 IABC=1,NABC                                                CPG31970
      IB=IBA(IABC)                                                      CPG31980
      CALL RREAD(ITAP44,BB(1,IABC),NTRI2,IB)                            CPG31990
  101 CONTINUE                                                          CPG32000
      WRITE(6,1)                                                        CPG32010
      CALL MATOUT(BB,NTRI,NABC,NTRI,NABC,6)                             CPG32020
      CALL SREW(ITAP44)                                                 CPG32030
C                                                                       CPG32040
      RETURN                                                            CPG32050
      END                                                               CPG32060
