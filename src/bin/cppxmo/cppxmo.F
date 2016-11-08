      PROGRAM CPPXMO                                                    CPP00010
C   COUPLED PERTURBED HARTREE-FOCK PROGRAM                              CPP00020
C   FOR THE CLOSED SHELL GVBSCF AND PAIR EXCITATION MCSCF WAVEFUNCTIONS CPP00030
C   IN THE MO BASIS                                                     CPP00040
C*******************************************************                CPP00050
C***LAST UPDATED ON MARCH 20, 1987 BY YUKIO YAMAGUCHI***                CPP00060
C*******************************************************                CPP00070
C                                                                       CPP00080
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP00090
      DIMENSION CC(1200000),IA(1),I30(200)                              CPP00100
      CHARACTER*80 ALABEL                                               CPP00110
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP00120
      COMMON/CALIF/LPARA(1024),APARA(1024)                              CPP00130
      COMMON/CIOFF/IIUNQ,IJUNQ,ICIUNQ,IPX(5050),JPX(5050)               CPP00140
      COMMON/CIPEX/PEXC(100),INDPEX(20,100),MOPEX(20),ICIMO(20)         CPP00150
      COMMON/CIPOS/ICITYP(5050),ICPOS(200,100)                          CPP00160
      COMMON/CIPMO/IIMO(5050),JJMO(5050)                                CPP00170
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP00180
      COMMON/CITYP/IJGRP,IBM(200),JBM(200),IBNUM(200)                   CPP00190
      COMMON/DENST/ALPC(465),BETC(465)                                  CPP00200
      COMMON/DIPOL/DIPDA(3,45),DIPDM(3,45),DIPDN(3,45),DIPDT(3,45)      CPP00210
      COMMON/ENRGY/ESCF,ENUC,ELEC                                       CPP00220
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP00230
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPP00240
      COMMON/MFSEC/MASTER,NSECT                                         CPP00250
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP00260
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP00270
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP00280
      COMMON/TAPES/ITAP35,ITAP42,ITAP43,ITAP44,ITAP47                   CPP00290
      COMMON/WORKS/JTAPE1,JTAPE2                                        CPP00300
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP00310
      COMMON/CI102/NIND,NDEP                                            CPP00320
      COMMON/CI103/CI(100)                                              CPP00330
      COMMON/CI104/IEIJ(5050)                                           CPP00340
      EQUIVALENCE (CC,IA)                                               CPP00350
      DATA ALABEL / ' SCF FIRST DERIVATIVE BROUGHT YOU BY Y & Y CO. EST.CPP00360
     1 IN MARCH , 1982             ' /                                  CPP00370
    1 FORMAT(//,2X,' COUPLED PERTURBED HARTREE-FOCK PROGRAM FOR GVB-SCF CPP00380
     1IN THE MO BASIS'/)                                                CPP00390
    2 FORMAT(10I5)                                                      CPP00400
    3 FORMAT(2X,10I5)                                                   CPP00410
    4 FORMAT(//,2X,' PARAMETERS'/                                       CPP00420
     1 2X,' NBASIS = ',I8/                                              CPP00430
     2 2X,' NTRI   = ',I8/                                              CPP00440
     3 2X,' NBFAO  = ',I8/                                              CPP00450
     4 2X,' NBATRI = ',I8/                                              CPP00460
     5 2X,' NATOMS = ',I8/                                              CPP00470
     6 2X,' N3N    = ',I8/                                              CPP00480
     7 2X,' N3NP3  = ',I8/                                              CPP00490
     8 2X,' N3NXX  = ',I8/                                              CPP00500
     9 2X,' IOCC   = ',I8/                                              CPP00510
     X 2X,' JOCC   = ',I8/                                              CPP00520
     1 2X,' KOCC   = ',I8/                                              CPP00530
     2 2X,' NTYPES = ',I8/                                              CPP00540
     3 2X,' NTYPEP = ',I8/                                              CPP00550
     4 2X,' NTREAD = ',I8/                                              CPP00560
     5 2X,' NTYTRI = ',I8)                                              CPP00570
    5 FORMAT(2X,' NCPEX  = ',I8/                                        CPP00580
     1 2X,' NORPEX = ',I8/                                              CPP00590
     2 2X,' NCTRI  = ',I8/                                              CPP00600
     3 2X,' LVPEX  = ',I8/                                              CPP00610
     4 2X,' ITEST  = ',I8/                                              CPP00620
     5 2X,' IORB   = ',I8/                                              CPP00630
     6 2X,' IPOL   = ',I8/                                              CPP00640
     7 2X,' ICORE  = ',I8/                                              CPP00650
     8 2X,' IPRNT  = ',I8)                                              CPP00660
    6 FORMAT(//,2X,' ALPC MATRIX'/)                                     CPP00670
    7 FORMAT(//,2X,' BETC MATRIX'/)                                     CPP00680
    8 FORMAT(//,2X,' SCF EIGENVECTOR (AO BASIS) MATRIX'/)               CPP00690
    9 FORMAT(//,2X,' SCF EIGENVECTOR (SO BASIS) MATRIX'/)               CPP00700
   10 FORMAT(//,2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)                CPP00710
   11 FORMAT(//,2X,' REQUIRED MEMORY EXCEEDS MAXCOR'/                   CPP00720
     1          2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)                CPP00730
C                                                                       CPP00740
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::             CPP00750
C   REFERENCES FOR TCSCF & GRSCF :                        :             CPP00760
C   Y.YAMAGUCHI, Y.OSAMURA, AND H.F.SCHAEFER III,         :             CPP00770
C   J.AM.CHEM.SOC. 105,7506(1983)                         :             CPP00780
C                 &                                       :             CPP00790
C   Y.OSAMURA, Y.YAMAGUCHI, P.SAXE, D.J.FOX, M.A.VINCENT, :             CPP00800
C   AND H.F.SCHAEFER III, J.MOL.STRUC. 103,183 (1983)     :             CPP00810
C                 &                                       :             CPP00820
C   Y.OSAMURA, Y.YAMAGUCHI, P.SAXE, M.A.VINCENT, J.F.GAW, :             CPP00830
C   AND H.F.SCHAEFER III, CHEM.PHYS. 72,131(1982)         :             CPP00840
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::             CPP00850
C                                                                       CPP00860
C   ITAP11 = GEOMETRY AND GRADIENTS                                     CPP00870
C   ITAP15 = SCF SECOND DERIVATIVES                                     CPP00880
C   ITAP17 = DIPOLE DERIVATIVES                                         CPP00890
C   ITAP25 = FIRST DERIVATIVES OF ORBITAL ENERGIES                      CPP00900
C   ITAP35 = TWO ELECTRON MO INTEGRALS                                  CPP00910
C   ITAP42 = FIRST AND SECOND DERIVATIVE INFORMATION                    CPP00920
C   ITAP43 = FIRST DERIVATIVES OF DIPOLE MOMENTS                        CPP00930
C   ITAP44 = FIRST DERIVATIVE MATRICES AND U'S                          CPP00940
C   ITAP47 = THE "BARE" GENERALIZED LAGRANGIAN MATRICES                 CPP00950
C   JTAPE1 = WORKING TAPE FOR CPHF                                      CPP00960
C   JTAPE2 = WORKING TAPE FOR E MATRICES                                CPP00970
C                                                                       CPP00980
      CALL TSTART(6)                                                    CPP00990
      CALL NOUNFL                                                       CPP01000
      CALL INITMF(1)                                                    CPP01010
C                                                                       CPP01020
      ITAPE3=3                                                          CPP01030
      ITAPE4=4                                                          CPP01040
      INPUT=5                                                           CPP01050
      ITAPE6=6                                                          CPP01060
      ITAP11=11                                                         CPP01070
      ITAP15=15                                                         CPP01080
      ITAP17=17                                                         CPP01090
      ITAP25=25                                                         CPP01100
      ITAP30=30                                                         CPP01110
      ITAP35=35                                                         CPP01120
      ITAP42=42                                                         CPP01130
      ITAP43=43                                                         CPP01140
      ITAP44=44                                                         CPP01150
      ITAP47=47                                                         CPP01160
      JTAPE1=91                                                         CPP01170
      JTAPE2=92                                                         CPP01180
      MAXCOR=1200000                                                    CPP01190
C                                                                       CPP01200
      IOFF(1)=0                                                         CPP01210
      DO 101 I=1,1378                                                   CPP01220
  101 IOFF(I+1)=IOFF(I)+I                                               CPP01230
C                                                                       CPP01240
      CALL LOCATE(INPUT,'# CPHFMO #',IERR)                              CPP01250
      WRITE(6,1)                                                        CPP01260
      WRITE(3,1)                                                        CPP01270
      WRITE(4,1)                                                        CPP01280
C                                                                       CPP01290
C   READIN PARAMETERS                                                   CPP01300
      READ(5,2) ITEST,IORB,IPOL,ICORE,IPRNT                             CPP01310
C                                                                       CPP01320
C   READIN SCF INFORMATION                                              CPP01330
      NBASIS=LPARA(3)                                                   CPP01340
      NBFAO=LPARA(11)                                                   CPP01350
      NTRI=LPARA(4)                                                     CPP01360
      NTRI2=NTRI*2                                                      CPP01370
      NTOT=IOFF(NTRI+1)                                                 CPP01380
      NBASQ=NBASIS*NBASIS*2                                             CPP01390
      NBATRI=IOFF(NBFAO+1)                                              CPP01400
      NBATR2=NBATRI*2                                                   CPP01410
      IOCC=LPARA(7)                                                     CPP01420
      KOCC=LPARA(8)                                                     CPP01430
      JOCC=LPARA(9)                                                     CPP01440
      NATOMS=LPARA(10)                                                  CPP01450
      N3N=LPARA(17)                                                     CPP01460
      NATRI=IOFF(N3N+1)                                                 CPP01470
      N3NP3=N3N+3                                                       CPP01480
      N3NXX=N3NP3                                                       CPP01490
      IF(IPOL.EQ.0) N3NXX=N3N                                           CPP01500
      NTYPES=LPARA(18)                                                  CPP01510
      NTYPEP=NTYPES+1                                                   CPP01520
      NTREAD=NTYPES*2                                                   CPP01530
      NTRIL=((NTRI2-1)/NSECT+1)                                         CPP01540
      NBASQL=((NBASQ-1)/NSECT+1)                                        CPP01550
      N3TRL=N3NXX*NTRIL                                                 CPP01560
      N3QRL=N3NXX*NBASQL                                                CPP01570
      ENUC=APARA(1)                                                     CPP01580
      ESCF=APARA(2)                                                     CPP01590
      DO 102 I=1,NTYPEP                                                 CPP01600
      NSORB(I)=LPARA(I+40)                                              CPP01610
      FOCC(I)=APARA(I+10)                                               CPP01620
  102 CONTINUE                                                          CPP01630
      NTYTRI=IOFF(NTYPEP+1)                                             CPP01640
      DO 103 I=1,NTYTRI                                                 CPP01650
      ALPC(I)=APARA(I+40)                                               CPP01660
      BETC(I)=APARA(I+505)                                              CPP01670
  103 CONTINUE                                                          CPP01680
      DO 104 I=1,N3NXX                                                  CPP01690
      ISA(I)=(I-1)*NTRIL+1                                              CPP01700
      IHA(I)=(I-1)*NTRIL+N3TRL+1                                        CPP01710
      IFA(I)=(I-1)*NTRIL+N3TRL+N3TRL+1                                  CPP01720
      IEA(I)=(I-1)*NBASQL+N3TRL+N3TRL+N3TRL+1                           CPP01730
      IBA(I)=(I-1)*NTRIL+N3TRL+N3TRL+N3TRL+N3QRL+1                      CPP01740
      IUA(I)=(I-1)*NBASQL+N3TRL+N3TRL+N3TRL+N3QRL+N3TRL+1               CPP01750
      WRITE(3,3) I,ISA(I),IHA(I),IFA(I),IEA(I),IBA(I),IUA(I)            CPP01760
  104 CONTINUE                                                          CPP01770
      NSECMX=IUA(N3NXX)+NBASQL                                          CPP01780
C                                                                       CPP01790
C   READ IN CI INFORMATION                                              CPP01800
      CALL RFILE(ITAP30)                                                CPP01810
      CALL WREADW(ITAP30,I30,200,101,JUNK)                              CPP01820
      NCPEX=I30(46)                                                     CPP01830
      NORPEX=I30(47)                                                    CPP01840
      LVPEX=I30(48)                                                     CPP01850
      NCTRI=IOFF(NCPEX+1)                                               CPP01860
      CALL WREADW(ITAP30,PEXC,2*NCPEX,LVPEX,LOCVEC)                     CPP01870
      CALL WREADW(ITAP30,INDPEX,20*NCPEX,LOCVEC,LOCVEC)                 CPP01880
      CALL WREADW(ITAP30,MOPEX,20,LOCVEC,LOCVEC)                        CPP01890
C                                                                       CPP01900
      WRITE(6,4) NBASIS,NTRI,NBFAO,NBATRI,NATOMS,N3N,N3NP3,N3NXX,       CPP01910
     1           IOCC,JOCC,KOCC,NTYPES,NTYPEP,NTREAD,NTYTRI             CPP01920
      WRITE(6,5) NCPEX,NORPEX,NCTRI,LVPEX,ITEST,IORB,IPOL,ICORE,        CPP01930
     1           IPRNT                                                  CPP01940
      WRITE(6,6)                                                        CPP01950
      CALL PRINT(ALPC,465,NTYPEP,6)                                     CPP01960
      WRITE(6,7)                                                        CPP01970
      CALL PRINT(BETC,465,NTYPEP,6)                                     CPP01980
C                                                                       CPP01990
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::CPP02000
C                                                                       CPP02010
C  DYNAMIC ALLOCATION                                                   CPP02020
C     IC1  = EIG(NBASIS)                                                CPP02030
C     IC2  = OCC(NBASIS)                                                CPP02040
C     IC3  = EAO OR ESO(NBFAO,NBASIS)                                   CPP02050
C     IC4  = NIJ(NTRI,2)                                                CPP02060
C     IC5  = KIJ(NTRI,2)                                                CPP02070
C     IC6  = ALPA(NTRI)                                                 CPP02080
C     IC7  = BETA(NTRI)                                                 CPP02090
C     IC8  = FIJ(NTYPEP,ICIUNQ)                                         CPP02100
C     IC9  = AIJ(NTYTRI,ICIUNQ)                                         CPP02110
C     IC10 = BIJ(NTYTRI,ICIUNQ)                                         CPP02120
C     IC11 = E1A(N3N)                                                   CPP02130
C     IC12 = E2A(N3N,N3N)                                               CPP02140
C     IC13 = D1W(N3N)                                                   CPP02150
C     IC14 = E1T(N3NXX)                                                 CPP02160
C     IC15 = EPS(NTRI)                                                  CPP02170
C     IC16 = ZETA(NTRI,NTYPEP)                                          CPP02180
C     IC17 = HCIA(NCTRI,N3NXX)                                          CPP02190
C     IC18 = HCI(NCTRI)                                                 CPP02200
C     IC19 = BAF2(NCPEX,N3NXX)                                          CPP02210
C     IC20 = CIA(NCPEX,N3NXX)                                           CPP02220
C     IC21 = A22(NCPEX,NCPEX)                                           CPP02230
C     IC22 = A12(NIND,NCPEX)                                            CPP02240
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::CPP02250
C                                                                       CPP02260
      IC1=1                                                             CPP02270
      IC2=IC1+NBASIS                                                    CPP02280
      IC3=IC2+NBASIS                                                    CPP02290
      IC4=IC3+NBFAO*NBASIS                                              CPP02300
      IA4=IC4+IC4-1                                                     CPP02310
      IC5=IC4+NTRI                                                      CPP02320
      IA5=IC5+IC5-1                                                     CPP02330
      IC6=IC5+NTRI                                                      CPP02340
      IC7=IC6+NTRI                                                      CPP02350
      IC8=IC7+NTRI                                                      CPP02360
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::CPP02370
C                                                                       CPP02380
C   READ IN AO-MO EIGENVECTORS                                          CPP02390
      CALL MREAD(CC(IC1),18)                                            CPP02400
      CALL MREAD(CC(IC2),19)                                            CPP02410
      CALL MREAD(CC(IC3),22)                                            CPP02420
      IF(IPRNT.LE.2) GO TO 301                                          CPP02430
      WRITE(6,8)                                                        CPP02440
      CALL EIGOUT(CC(IC3),CC(IC1),CC(IC2),NBFAO,NBASIS,NBFAO,NBASIS,6)  CPP02450
C                                                                       CPP02460
C**************************                                             CPP02470
C***PREPARATION FOR CPHF***                                             CPP02480
C**************************                                             CPP02490
C                                                                       CPP02500
C   FORM REGISTER MATRICES                                              CPP02510
  301 CONTINUE                                                          CPP02520
      WRITE(3,20)                                                       CPP02530
      WRITE(4,20)                                                       CPP02540
   20 FORMAT(//,2X,' NOW YOU ARE IN REGIST'/)                           CPP02550
C.................NIJ     KIJ.....                                      CPP02560
      CALL REGIST(IA(IA4),IA(IA5))                                      CPP02570
C                                                                       CPP02580
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::CPP02590
C  REST OF DYNAMIC ALLOCATION                                           CPP02600
      IC9=IC8+NTYPEP*ICIUNQ                                             CPP02610
      IC10=IC9+NTYTRI*ICIUNQ                                            CPP02620
      IC11=IC10+NTYTRI*ICIUNQ                                           CPP02630
      IC12=IC11+N3N                                                     CPP02640
      IC13=IC12+N3N*N3N                                                 CPP02650
      IC14=IC13+N3N                                                     CPP02660
      IC15=IC14+N3NXX                                                   CPP02670
      IC16=IC15+NTRI                                                    CPP02680
      IC17=IC16+NTRI*NTYPEP                                             CPP02690
      IC18=IC17+NCTRI*N3NXX                                             CPP02700
      IC19=IC18+NCTRI                                                   CPP02710
      IC20=IC19+NCPEX*N3NXX                                             CPP02720
      IC21=IC20+NCPEX*N3NXX                                             CPP02730
      IC22=IC21+NCPEX*NCPEX                                             CPP02740
      ICKEEP=IC22+NIND*NCPEX                                            CPP02750
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::CPP02760
C                                                                       CPP02770
      DO 105 I=1,IIUNQ+IJGRP                                            CPP02780
      IEIJ(I)=(I-1)*NBASQL+1                                            CPP02790
      IF(IPRNT.LE.2) GO TO 105                                          CPP02800
      WRITE(3,3) I,IEIJ(I)                                              CPP02810
  105 CONTINUE                                                          CPP02820
C                                                                       CPP02830
C   FORM ALPA AND BETA MATRICES                                         CPP02840
      WRITE(3,21)                                                       CPP02850
      WRITE(4,21)                                                       CPP02860
   21 FORMAT(//,2X,' NOW YOU ARE IN ABMAT'/)                            CPP02870
C................ALPA    BETA....                                       CPP02880
      CALL ABMAT(CC(IC6),CC(IC7))                                       CPP02890
C                                                                       CPP02900
C   FORM THE COUPLING CONSTANT MATRICES                                 CPP02910
      WRITE(3,22)                                                       CPP02920
      WRITE(4,22)                                                       CPP02930
   22 FORMAT(//,2X,' NOW YOU ARE IN BARES'/)                            CPP02940
C................FIJ     AIJ     BIJ......                              CPP02950
      CALL BARES(CC(IC8),CC(IC9),CC(IC10))                              CPP02960
C                                                                       CPP02970
C   READ IN FIRST AND SECOND DERIVATIVE INFORMATION                     CPP02980
C   AND FORM SM, HM AND FM MATRICES                                     CPP02990
      WRITE(3,23)                                                       CPP03000
      WRITE(4,23)                                                       CPP03010
   23 FORMAT(//,2X,' NOW YOU ARE IN DERMAT'/)                           CPP03020
      IC23=ICKEEP                                                       CPP03030
      IC24=IC23+NBATRI                                                  CPP03040
      IC25=IC24+NBFAO*NBFAO                                             CPP03050
      ICMAX=IC25+NBFAO*NBFAO                                            CPP03060
      WRITE(3,10) ICMAX,MAXCOR                                          CPP03070
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP03080
C.................EAO     E1A      E2A      SS       U        T........ CPP03090
      CALL DERMAT(CC(IC3),CC(IC11),CC(IC12),CC(IC23),CC(IC24),CC(IC25)) CPP03100
C                                                                       CPP03110
C   READ IN SO-MO EIGENVECTORS                                          CPP03120
      CALL MREAD(CC(IC3),21)                                            CPP03130
      IF(IPRNT.LE.2) GO TO 302                                          CPP03140
      WRITE(6,9)                                                        CPP03150
      CALL EIGOUT(CC(IC3),CC(IC1),CC(IC2),NBFAO,NBASIS,NBFAO,NBASIS,6)  CPP03160
C                                                                       CPP03170
C:::::::::::::::::::::::::::::::::::                                    CPP03180
C:::TEST ROUTINES FOR INPUT CHECK:::                                    CPP03190
C:::::::::::::::::::::::::::::::::::                                    CPP03200
  302 CONTINUE                                                          CPP03210
      IF(ITEST.EQ.0) GO TO 303                                          CPP03220
      WRITE(3,24)                                                       CPP03230
      WRITE(4,24)                                                       CPP03240
   24 FORMAT(//,2X,' NOW YOU ARE IN E0CAL'/)                            CPP03250
C   CALCULATE SCF ENERGY FOR A TEST                                     CPP03260
      IC23=ICKEEP                                                       CPP03270
      IC24=IC23+NTRI                                                    CPP03280
      IC25=IC24+NTOT                                                    CPP03290
      IA25=IC25+IC25-1                                                  CPP03300
      ICMAX=IC25+MAXBF2                                                 CPP03310
      WRITE(3,10) ICMAX,MAXCOR                                          CPP03320
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP03330
C................ALPA    BETA    H        CC       LBLI     BUFI.....   CPP03340
      CALL E0CAL(CC(IC6),CC(IC7),CC(IC23),CC(IC24),IA(IA25),CC(IC25),   CPP03350
C.....................                                                  CPP03360
     1           NTOT)                                                  CPP03370
C                                                                       CPP03380
C   CALCULATE LAGRANGIAN AND ZETA MATRICES                              CPP03390
  303 CONTINUE                                                          CPP03400
      WRITE(3,25)                                                       CPP03410
      WRITE(4,25)                                                       CPP03420
   25 FORMAT(//,2X,' NOW YOU ARE IN E0LAG'/)                            CPP03430
      IC23=ICKEEP                                                       CPP03440
      IC24=IC23+NTRI                                                    CPP03450
      IC25=IC24+NTOT                                                    CPP03460
      IA25=IC25+IC25-1                                                  CPP03470
      ICMAX=IC25+MAXBF2                                                 CPP03480
      WRITE(3,10) ICMAX,MAXCOR                                          CPP03490
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP03500
C................ALPA    BETA    EPS      ZETA     H        CC.......   CPP03510
      CALL E0LAG(CC(IC6),CC(IC7),CC(IC15),CC(IC16),CC(IC23),CC(IC24),   CPP03520
C................LBLI     BUFI..........                                CPP03530
     1           IA(IA25),CC(IC25),NTOT)                                CPP03540
C                                                                       CPP03550
C   CALCULATE SCF FIRST ENERGY DERIVATIVES                              CPP03560
      IF(ITEST.EQ.0) GO TO 304                                          CPP03570
      WRITE(3,26)                                                       CPP03580
      WRITE(4,26)                                                       CPP03590
   26 FORMAT(//,2X,' NOW YOU ARE IN MODER'/)                            CPP03600
      IC23=ICKEEP                                                       CPP03610
      IC24=IC23+N3N                                                     CPP03620
      IC25=IC24+N3N                                                     CPP03630
      IC26=IC25+NTRI*N3N                                                CPP03640
      IC27=IC26+NTRI*N3N                                                CPP03650
      ICMAX=IC27+NTRI*N3N*NTREAD                                        CPP03660
      WRITE(3,10) ICMAX,MAXCOR                                          CPP03670
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP03680
C................ALPA    BETA    D1W      E1T      EPS......            CPP03690
      CALL MODER(CC(IC6),CC(IC7),CC(IC13),CC(IC14),CC(IC15),            CPP03700
C................D1O,     D1T      SM       HM       TM.......          CPP03710
     1           CC(IC23),CC(IC24),CC(IC25),CC(IC26),CC(IC27))          CPP03720
C                                                                       CPP03730
C   FORM HA MATRICES                                                    CPP03740
  304 CONTINUE                                                          CPP03750
      WRITE(3,27)                                                       CPP03760
      WRITE(4,27)                                                       CPP03770
   27 FORMAT(//,2X,' NOW YOU ARE IN HAMAT'/)                            CPP03780
      IC23=ICKEEP                                                       CPP03790
      IC24=IC23+N3N                                                     CPP03800
      IC25=IC24+NTRI*N3N                                                CPP03810
      ICMAX=IC25+NTRI*N3N*NTREAD                                        CPP03820
      WRITE(3,10) ICMAX,MAXCOR                                          CPP03830
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP03840
C................FIJ     AIJ     BIJ      D1W      E1T      HCIA.....   CPP03850
      CALL HAMAT(CC(IC8),CC(IC9),CC(IC10),CC(IC13),CC(IC14),CC(IC17),   CPP03860
C................D1E      HM       TM........................           CPP03870
     1           CC(IC23),CC(IC24),CC(IC25),N3N*NTREAD,N3NXX)           CPP03880
C                                                                       CPP03890
C   FORM DERIVATIVE LAGRANGIAN MATRICES                                 CPP03900
      WRITE(3,28)                                                       CPP03910
      WRITE(4,28)                                                       CPP03920
   28 FORMAT(//,2X,' NOW YOU ARE IN FAMAT'/)                            CPP03930
      IC23=ICKEEP                                                       CPP03940
      IC24=IC23+NTRI                                                    CPP03950
      IC25=IC24+NBASIS*NBASIS                                           CPP03960
      IC26=IC25+NTRI*N3N                                                CPP03970
      ICMAX=IC26+NTRI*N3N*NTREAD                                        CPP03980
      WRITE(3,10) ICMAX,MAXCOR                                          CPP03990
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP04000
C................ALPA    BETA    FF       EPA      HM       TM.......   CPP04010
      CALL FAMAT(CC(IC6),CC(IC7),CC(IC23),CC(IC24),CC(IC25),CC(IC26),   CPP04020
C...........................                                            CPP04030
     1           N3N*NTREAD)                                            CPP04040
C                                                                       CPP04050
C   CALCULATE BA MATRICES                                               CPP04060
      WRITE(3,29)                                                       CPP04070
      WRITE(4,29)                                                       CPP04080
   29 FORMAT(//,2X,' NOW YOU ARE IN BAMAT'/)                            CPP04090
      IC23=ICKEEP                                                       CPP04100
      IC24=IC23+NTRI*N3N                                                CPP04110
      IC25=IC24+NTRI                                                    CPP04120
      IC26=IC25+NBASIS*NBASIS                                           CPP04130
      IC27=IC26+NTRI                                                    CPP04140
      IC28=IC27+NTOT                                                    CPP04150
      IA28=IC28+IC28-1                                                  CPP04160
      ICMAX=IC28+MAXBF2                                                 CPP04170
      WRITE(3,10) ICMAX,MAXCOR                                          CPP04180
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP04190
C................EIG     NIJ     KIJ     ALPA    BETA    ZETA.....      CPP04200
      CALL BAMAT(CC(IC1),IA(IA4),IA(IA5),CC(IC6),CC(IC7),CC(IC16),      CPP04210
C................SM       B0       EPA      FA       CC.......          CPP04220
     1           CC(IC23),CC(IC24),CC(IC25),CC(IC26),CC(IC27),          CPP04230
C................LBLI     BUFI..........                                CPP04240
     2           IA(IA28),CC(IC28),NTOT)                                CPP04250
C     WRITE(6,*) ' AFTER BAMAT'                                         CPP04260
C     CALL WRBMAT(CC(IC23),N3NXX)                                       CPP04270
C                                                                       CPP04280
C   CALCULATE CI MATRIX AND BARE LAGRANGIAN MATRICES                    CPP04290
      WRITE(3,30)                                                       CPP04300
      WRITE(4,30)                                                       CPP04310
   30 FORMAT(//,2X,' NOW YOU ARE IN HEMAT'/)                            CPP04320
      IC23=ICKEEP                                                       CPP04330
      IC24=IC23+NBASIS*NBASIS                                           CPP04340
      IC25=IC24+NTRI                                                    CPP04350
      IC26=IC25+NBASIS*NBASIS                                           CPP04360
      IC27=IC26+NTRI                                                    CPP04370
      IC28=IC27+NTOT                                                    CPP04380
      IA28=IC28+IC28-1                                                  CPP04390
      ICMAX=IC28+MAXBF2                                                 CPP04400
      WRITE(3,10) ICMAX,MAXCOR                                          CPP04410
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP04420
C................FIJ     AIJ     BIJ      HCI      EIJ      ZIJ......   CPP04430
      CALL HEMAT(CC(IC8),CC(IC9),CC(IC10),CC(IC18),CC(IC23),CC(IC24),   CPP04440
C................EE       H        CC       LBLI     BUFI..........     CPP04450
     1           CC(IC25),CC(IC26),CC(IC27),IA(IA28),CC(IC28),NTOT)     CPP04460
C                                                                       CPP04470
C   CALCULATE A MATRIX                                                  CPP04480
      WRITE(3,31)                                                       CPP04490
      WRITE(4,31)                                                       CPP04500
   31 FORMAT(//,2X,' NOW YOU ARE IN AMAT'/)                             CPP04510
      IC23=ICKEEP                                                       CPP04520
      IC24=IC23+NTOT                                                    CPP04530
      IA24=IC24+IC24-1                                                  CPP04540
      ICMAX=IC24+MAXBF2                                                 CPP04550
      WRITE(3,10) ICMAX,MAXCOR                                          CPP04560
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP04570
C...............NIJ     KIJ     ALPA    BETA    ZETA     CC.......      CPP04580
      CALL AMAT(IA(IA4),IA(IA5),CC(IC6),CC(IC7),CC(IC16),CC(IC23),      CPP04590
C...............LBLI     BUFI..........                                 CPP04600
     1          IA(IA24),CC(IC24),NTOT)                                 CPP04610
C                                                                       CPP04620
C   FORM BF MATRICES                                                    CPP04630
      WRITE(3,32)                                                       CPP04640
      WRITE(4,32)                                                       CPP04650
   32 FORMAT(//,2X,' NOW YOU ARE IN BFMAT'/)                            CPP04660
      IC23=ICKEEP                                                       CPP04670
      IC24=IC23+NTRI*3                                                  CPP04680
      IC25=IC24+NTRI*NTYPES*3                                           CPP04690
      IC26=IC25+NBASIS*NBASIS                                           CPP04700
      ICMAX=IC26+NTRI                                                   CPP04710
      WRITE(3,10) ICMAX,MAXCOR                                          CPP04720
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP04730
C................NIJ     KIJ     FIJ     E1T      HCIA     HF.......    CPP04740
      CALL BFMAT(IA(IA4),IA(IA5),CC(IC8),CC(IC14),CC(IC17),CC(IC23),    CPP04750
C................HM       EPF      BF.............                      CPP04760
     1           CC(IC24),CC(IC25),CC(IC26),N3NXX)                      CPP04770
C*    WRITE(6,*) ' AFTER BFMAT'                                         CPP04780
C*    CALL WRBMAT(CC(IC23),N3NXX)                                       CPP04790
C                                                                       CPP04800
C   FORM BAF MATRICES                                                   CPP04810
      WRITE(3,33)                                                       CPP04820
      WRITE(4,33)                                                       CPP04830
   33 FORMAT(//,2X,' NOW YOU ARE IN AB2CI'/)                            CPP04840
      IC23=ICKEEP                                                       CPP04850
      IC24=IC23+NTRI*N3N                                                CPP04860
      ICMAX=IC24+NBASIS*NBASIS                                          CPP04870
      WRITE(3,10) ICMAX,MAXCOR                                          CPP04880
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP04890
C................NIJ     E1T      HCIA     HCI      BAF2     A22......  CPP04900
      CALL AB2CI(IA(IA4),CC(IC14),CC(IC17),CC(IC18),CC(IC19),CC(IC21),  CPP04910
C................A12      SA       EIJ............                      CPP04920
     1           CC(IC22),CC(IC23),CC(IC24),N3NXX)                      CPP04930
C                                                                       CPP04940
C********************                                                   CPP04950
C***CPHF PROCEDURE***                                                   CPP04960
C********************                                                   CPP04970
C                                                                       CPP04980
C   SOLVE SIMULTANEOUS EQUATIONS DIRECTLY                               CPP04990
      WRITE(3,34)                                                       CPP05000
      WRITE(4,34)                                                       CPP05010
   34 FORMAT(//,2X,' NOW YOU ARE IN CPHF'/)                             CPP05020
      NINDPX=NIND+NCPEX                                                 CPP05030
      IC23=ICKEEP                                                       CPP05040
      IC24=IC23+NTRI                                                    CPP05050
      IC25=IC24+NBASIS*NBASIS                                           CPP05060
      IC26=IC25+NTRI*N3N                                                CPP05070
      IC27=IC26+NTRI*N3NXX                                              CPP05080
      IC28=IC27+NINDPX*(NINDPX+N3NXX)                                   CPP05090
      IA28=IC28+IC28-1                                                  CPP05100
      ICMAX=IC28+MAXBF2                                                 CPP05110
      WRITE(3,10) ICMAX,MAXCOR                                          CPP05120
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP05130
C................EIG     NIJ     KIJ     BAF2     CIA......             CPP05140
       CALL CPHF(CC(IC1),IA(IA4),IA(IA5),CC(IC19),CC(IC20),             CPP05150
C................A22      A12      B0       UA       SA       BA....... CPP05160
     1           CC(IC21),CC(IC22),CC(IC23),CC(IC24),CC(IC25),CC(IC26), CPP05170
C................C        LBLI     BUFI..................               CPP05180
     2           CC(IC27),IA(IA28),CC(IC28),NINDPX,N3NXX)               CPP05190
C                                                                       CPP05200
C   FORM THE WA MATRICES                                                CPP05210
      WRITE(3,35)                                                       CPP05220
      WRITE(4,35)                                                       CPP05230
   35 FORMAT(//,2X,' NOW YOU ARE IN WAMAT'/)                            CPP05240
      IC23=ICKEEP                                                       CPP05250
      IC24=IC23+N3N*N3N                                                 CPP05260
      IC25=IC24+N3N*N3N                                                 CPP05270
      IC26=IC25+N3N*N3N                                                 CPP05280
      IC27=IC26+NBASIS*NBASIS                                           CPP05290
      IC28=IC27+NBASIS*NBASIS                                           CPP05300
      IC29=IC28+NTRI*N3N                                                CPP05310
      IC30=IC29+NBASIS*NBASIS*N3N                                       CPP05320
      IC31=IC30+NTOT                                                    CPP05330
      IA31=IC31+IC31-1                                                  CPP05340
      ICMAX=IC31+MAXBF2                                                 CPP05350
      WRITE(3,10) ICMAX,MAXCOR                                          CPP05360
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP05370
C................ALPA    BETA    EPS      ZETA     E2M      E2P......   CPP05380
      CALL WAMAT(CC(IC6),CC(IC7),CC(IC15),CC(IC16),CC(IC23),CC(IC24),   CPP05390
C................E2Q      EPA      UA       SA       WA       CC....... CPP05400
     1           CC(IC25),CC(IC26),CC(IC27),CC(IC28),CC(IC29),CC(IC30), CPP05410
C................LBLI     BUFI..........                                CPP05420
     2           IA(IA31),CC(IC31),NTOT)                                CPP05430
C                                                                       CPP05440
C   CALCULATE SCF SECOND DERIVATIVES                                    CPP05450
      WRITE(3,36)                                                       CPP05460
      WRITE(4,36)                                                       CPP05470
   36 FORMAT(//,2X,' NOW YOU ARE IN SCF2ND'/)                           CPP05480
      IC23=ICKEEP                                                       CPP05490
      IC24=IC23+N3N*N3N                                                 CPP05500
      IC25=IC24+N3N*N3N                                                 CPP05510
      IC26=IC25+N3N*N3N                                                 CPP05520
      IC27=IC26+NTRI*N3N                                                CPP05530
      ICMAX=IC27+NBASIS*NBASIS                                          CPP05540
      WRITE(3,10) ICMAX,MAXCOR                                          CPP05550
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP05560
C.................E2A      HCIA     CIA      E2M      E2C      E2T......CPP05570
      CALL SCF2ND(CC(IC12),CC(IC17),CC(IC20),CC(IC23),CC(IC24),CC(IC25),CPP05580
C.................SA       EIJ......                                    CPP05590
     1            CC(IC26),CC(IC27))                                    CPP05600
C                                                                       CPP05610
C   CALCULATE DIPOLE MOMENT DERIVATIVES                                 CPP05620
      WRITE(3,37)                                                       CPP05630
      WRITE(4,37)                                                       CPP05640
   37 FORMAT(//,2X,' NOW YOU ARE IN DIPMO'/)                            CPP05650
      IC23=ICKEEP                                                       CPP05660
      IC24=IC23+NTRI*3                                                  CPP05670
      ICMAX=IC24+NBASIS*NBASIS*N3N                                      CPP05680
      WRITE(3,10) ICMAX,MAXCOR                                          CPP05690
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP05700
C................HCIA     CIA      HF       UA.......                   CPP05710
      CALL DIPMO(CC(IC17),CC(IC20),CC(IC23),CC(IC24))                   CPP05720
C                                                                       CPP05730
C   FIND ORBITAL ENERGY DERIVATIVES                                     CPP05740
      IF(IORB.EQ.0) GO TO 305                                           CPP05750
      WRITE(3,38)                                                       CPP05760
      WRITE(4,38)                                                       CPP05770
   38 FORMAT(//,2X,' NOW YOU ARE IN ORB1ST'/)                           CPP05780
      IC23=ICKEEP                                                       CPP05790
      IC24=IC23+NTRI                                                    CPP05800
      ICMAX=IC24+NBASIS*N3N                                             CPP05810
      WRITE(3,10) ICMAX,MAXCOR                                          CPP05820
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP05830
C.................B0       CC.........                                  CPP05840
      CALL ORB1ST(CC(IC23),CC(IC24),1)                                  CPP05850
      IF(IPOL.EQ.0) GO TO 305                                           CPP05860
C.................B0       CC.........                                  CPP05870
      CALL ORB1ST(CC(IC23),CC(IC24),2)                                  CPP05880
C                                                                       CPP05890
C   CALCULATE ELECTRIC POLARIZABILITIES                                 CPP05900
  305 CONTINUE                                                          CPP05910
      IF(IPOL.EQ.0) GO TO 310                                           CPP05920
      WRITE(3,39)                                                       CPP05930
      WRITE(4,39)                                                       CPP05940
   39 FORMAT(//,2X,' NOW YOU ARE IN POLAR'/)                            CPP05950
      IC23=ICKEEP                                                       CPP05960
      IC24=IC23+NTRI*3                                                  CPP05970
      ICMAX=IC24+NBASIS*NBASIS*3                                        CPP05980
      WRITE(3,10) ICMAX,MAXCOR                                          CPP05990
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP06000
C................HCIA     CIA      HF       UF.......                   CPP06010
      CALL POLAR(CC(IC17),CC(IC20),CC(IC23),CC(IC24))                   CPP06020
C                                                                       CPP06030
C   STORE CONTRIBUTION TO THIRD DERIVATIVES                             CPP06040
  310 CONTINUE                                                          CPP06050
      WRITE(3,40)                                                       CPP06060
      WRITE(4,40)                                                       CPP06070
   40 FORMAT(//,2X,' NOW YOU ARE IN THIRD'/)                            CPP06080
C................E1T      HCIA     CIA.............                     CPP06090
      CALL THIRD(CC(IC14),CC(IC17),CC(IC20),NSECMX)                     CPP06100
C                                                                       CPP06110
      GO TO 400                                                         CPP06120
  399 CONTINUE                                                          CPP06130
      WRITE(6,11) ICMAX,MAXCOR                                          CPP06140
      WRITE(3,11) ICMAX,MAXCOR                                          CPP06150
  400 CONTINUE                                                          CPP06160
      CALL RCLOSE(MASTER,3)                                             CPP06170
      CALL RCLOSE(ITAP30,3)                                             CPP06180
      CALL RCLOSE(ITAP42,3)                                             CPP06190
      CALL RCLOSE(ITAP43,3)                                             CPP06200
      CALL RCLOSE(ITAP44,3)                                             CPP06210
      CALL RCLOSE(ITAP47,3)                                             CPP06220
      CALL RCLOSE(JTAPE1,3)                                             CPP06230
      CALL RCLOSE(JTAPE2,3)                                             CPP06240
C                                                                       CPP06250
      CALL TSTOP(6)                                                     CPP06260
      STOP                                                              CPP06270
      END                                                               CPP06280
      SUBROUTINE REGIST(NIJ,KIJ)                                        CPP06290
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP06300
      DIMENSION NIJ(NTRI,2),KIJ(NTRI,2)                                 CPP06310
      DIMENSION LOTYP(256)                                              CPP06320
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP06330
      COMMON/CIOFF/IIUNQ,IJUNQ,ICIUNQ,IPX(5050),JPX(5050)               CPP06340
      COMMON/CIPEX/PEXC(100),INDPEX(20,100),MOPEX(20),ICIMO(20)         CPP06350
      COMMON/CIPMO/IIMO(5050),JJMO(5050)                                CPP06360
      COMMON/CIPOS/ICITYP(5050),ICPOS(200,100)                          CPP06370
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP06380
      COMMON/CITYP/IJGRP,IBM(200),JBM(200),IBNUM(200)                   CPP06390
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP06400
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP06410
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP06420
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP06430
      COMMON/CI102/NIND,NDEP                                            CPP06440
      DATA A0,ONE,TWO / 0.0D+00 , 1.0D+00 , 2.0D+00 /                   CPP06450
    1 FORMAT(//,2X,' REGISTERED NUMBERS'/                               CPP06460
     1 2X,' KVIR = ',I5/                                                CPP06470
     2 2X,' NIND = ',I5/                                                CPP06480
     3 2X,' NDEP = ',I5/                                                CPP06490
     4 2X,' NDIF = ',I5)                                                CPP06500
    2 FORMAT(//,2X,' NIJ MATRIX'/2X,' NO.',5X,' NIJ(1)',5X,' NIJ(2)'/)  CPP06510
    3 FORMAT(2X,I4,4X,I5,7X,I5)                                         CPP06520
    4 FORMAT(//,2X,' KIJ MATRIX'/2X,' NO.',5X,' KIJ(1)',5X,' KIJ(2)'/)  CPP06530
    5 FORMAT(//,2X,' NO.',8X,' FOCC',6X,' NSORB',7X,' NSTR',7X,' NEND'/)CPP06540
    6 FORMAT(2X,I4,3X,F10.5,7X,I5,7X,I5,7X,I5)                          CPP06550
    7 FORMAT(//,2X,' MOPEX VECTOR '//16I5)                              CPP06560
    8 FORMAT(//,2X,' LOTYP VECTOR '//16I5)                              CPP06570
    9 FORMAT(//,2X,' ICIMO VECTOR '//16I5)                              CPP06580
   10 FORMAT(//,2X,' INDPEX VECTOR'/)                                   CPP06590
   11 FORMAT(2X,I5,5X,10I5)                                             CPP06600
   12 FORMAT(//,2X,' IIUNQ = ',I6,4X,' IJUNQ = ',I6,4X,' ICIUNQ = ',I6) CPP06610
   13 FORMAT(//,4X,' IJCI    IEX    JEX    III    JJJ   IIMO   JJMO',   CPP06620
     1          2X,' ICMO   JCMO  IJGRP'/)                              CPP06630
   14 FORMAT(2X,10I7)                                                   CPP06640
   15 FORMAT(//,2X,' IJGRP = ',I6)                                      CPP06650
   16 FORMAT(//,4X,' IBCI    ICM    JCM    NUM'/)                       CPP06660
C                                                                       CPP06670
      KVIR=NBASIS-JOCC                                                  CPP06680
C                                                                       CPP06690
      NSTR(1)=1                                                         CPP06700
      NEND(1)=NSORB(1)                                                  CPP06710
      DO 101 ITYP=2,NTYPEP                                              CPP06720
      NSTR(ITYP)=NEND(ITYP-1)+1                                         CPP06730
      NEND(ITYP)=NSTR(ITYP)+NSORB(ITYP)-1                               CPP06740
  101 CONTINUE                                                          CPP06750
      DO 103 ITYP=1,NTYPEP                                              CPP06760
      DO 102 I=NSTR(ITYP),NEND(ITYP)                                    CPP06770
  102 MOTYP(I)=ITYP                                                     CPP06780
  103 CONTINUE                                                          CPP06790
C                                                                       CPP06800
      II=0                                                              CPP06810
      DO 105 ITYP=2,NTYPEP                                              CPP06820
      DO 105 JTYP=1,ITYP-1                                              CPP06830
      DO 104 I=NSTR(ITYP),NEND(ITYP)                                    CPP06840
      DO 104 J=NSTR(JTYP),NEND(JTYP)                                    CPP06850
      II=II+1                                                           CPP06860
      NIJ(II,1)=I                                                       CPP06870
      NIJ(II,2)=J                                                       CPP06880
  104 CONTINUE                                                          CPP06890
  105 CONTINUE                                                          CPP06900
      NIND=II                                                           CPP06910
C                                                                       CPP06920
      II=0                                                              CPP06930
      DO 107 ITYP=1,NTYPEP                                              CPP06940
C*C   IF(NSORB(ITYP).EQ.1) GO TO 107                                    CPP06950
C*C   DO 106 I=NSTR(ITYP)+1,NEND(ITYP)                                  CPP06960
C*C   DO 106 J=NSTR(ITYP),I-1                                           CPP06970
      DO 106 I=NSTR(ITYP),NEND(ITYP)                                    CPP06980
      DO 106 J=NSTR(ITYP),I                                             CPP06990
      II=II+1                                                           CPP07000
      KIJ(II,1)=I                                                       CPP07010
      KIJ(II,2)=J                                                       CPP07020
  106 CONTINUE                                                          CPP07030
  107 CONTINUE                                                          CPP07040
      NDEP=II                                                           CPP07050
C                                                                       CPP07060
      NDIF=NTRI-NIND-NDEP                                               CPP07070
      IF(IPRNT.LE.2) GO TO 201                                          CPP07080
      WRITE(6,1) KVIR,NIND,NDEP,NDIF                                    CPP07090
      WRITE(6,2)                                                        CPP07100
      DO 108 I=1,NIND                                                   CPP07110
  108 WRITE(6,3) I,NIJ(I,1),NIJ(I,2)                                    CPP07120
      WRITE(6,4)                                                        CPP07130
      DO 109 I=1,NDEP                                                   CPP07140
  109 WRITE(6,3) I,KIJ(I,1),KIJ(I,2)                                    CPP07150
      WRITE(6,5)                                                        CPP07160
      DO 110 I=1,NTYPEP                                                 CPP07170
  110 WRITE(6,6) I,FOCC(I),NSORB(I),NSTR(I),NEND(I)                     CPP07180
C                                                                       CPP07190
C   INTERRELATE THE MO'S IN THE ACTIVE SPACE                            CPP07200
  201 CONTINUE                                                          CPP07210
      CALL MREAD(LOTYP,20)                                              CPP07220
      IF(IPRNT.GT.2)                                                    CPP07230
     *WRITE(6,7) (MOPEX(I),I=1,NORPEX)                                  CPP07240
      IF(IPRNT.GT.2)                                                    CPP07250
     *WRITE(6,8) (LOTYP(I),I=1,NBASIS)                                  CPP07260
      DO 113 I=1,NORPEX                                                 CPP07270
      MM=MOPEX(I)                                                       CPP07280
      DO 112 J=1,NBASIS                                                 CPP07290
      IF(LOTYP(J).NE.MM) GO TO 112                                      CPP07300
      MTEMP=J                                                           CPP07310
  112 CONTINUE                                                          CPP07320
      ICIMO(I)=MTEMP                                                    CPP07330
  113 CONTINUE                                                          CPP07340
      IF(IPRNT.GT.2)                                                    CPP07350
     *WRITE(6,9) (ICIMO(I),I=1,NORPEX)                                  CPP07360
C                                                                       CPP07370
C   GENERATE DIFFERENCE VECTOR                                          CPP07380
      IF(IPRNT.LE.2) GO TO 205                                          CPP07390
      WRITE(6,10)                                                       CPP07400
      DO 115 I=1,NCPEX                                                  CPP07410
      WRITE(6,11) I,(INDPEX(J,I),J=1,NORPEX)                            CPP07420
  115 CONTINUE                                                          CPP07430
  205 CONTINUE                                                          CPP07440
      IIUNQ=NCPEX                                                       CPP07450
      IJUNQ=0                                                           CPP07460
      DO 117 I=1,NCPEX                                                  CPP07470
      IPX(I)=I                                                          CPP07480
      JPX(I)=I                                                          CPP07490
      ICITYP(I)=I                                                       CPP07500
      DO 117 J=1,I                                                      CPP07510
      NN=0                                                              CPP07520
      DO 116 K=1,NORPEX                                                 CPP07530
  116 IF(INDPEX(K,I).NE.INDPEX(K,J)) NN=NN+1                            CPP07540
      NN=NN/2                                                           CPP07550
      IF(NN.NE.1) GO TO 117                                             CPP07560
      IJUNQ=IJUNQ+1                                                     CPP07570
      IPX(IIUNQ+IJUNQ)=I                                                CPP07580
      JPX(IIUNQ+IJUNQ)=J                                                CPP07590
  117 CONTINUE                                                          CPP07600
      ICIUNQ=IIUNQ+IJUNQ                                                CPP07610
      WRITE(6,12) IIUNQ,IJUNQ,ICIUNQ                                    CPP07620
      IF(IPRNT.GT.2)                                                    CPP07630
     *WRITE(6,13)                                                       CPP07640
C                                                                       CPP07650
      IJGRP=0                                                           CPP07660
      CALL IZERO(IBNUM,200)                                             CPP07670
      DO 120 IJCI=1,IJUNQ                                               CPP07680
      IEX=IPX(IIUNQ+IJCI)                                               CPP07690
      JEX=JPX(IIUNQ+IJCI)                                               CPP07700
      DO 118 K=1,NORPEX                                                 CPP07710
      IF(INDPEX(K,IEX).EQ.INDPEX(K,JEX)) GO TO 118                      CPP07720
      IF(INDPEX(K,IEX).NE.0) III=K                                      CPP07730
      IF(INDPEX(K,JEX).NE.0) JJJ=K                                      CPP07740
  118 CONTINUE                                                          CPP07750
      IIMO(IJCI)=MOTYP(ICIMO(III))                                      CPP07760
      JJMO(IJCI)=MOTYP(ICIMO(JJJ))                                      CPP07770
      ICM=ICIMO(III)                                                    CPP07780
      JCM=ICIMO(JJJ)                                                    CPP07790
      IF(IJCI.EQ.1) GO TO 206                                           CPP07800
      ITEST=0                                                           CPP07810
      DO 119 IBCI=1,IJGRP                                               CPP07820
      IF(ICM.EQ.IBM(IBCI).AND.JCM.EQ.JBM(IBCI)) THEN                    CPP07830
        IBNUM(IBCI)=IBNUM(IBCI)+1                                       CPP07840
        INUM=IBNUM(IBCI)                                                CPP07850
        ITEST=ITEST+1                                                   CPP07860
        ICITYP(IIUNQ+IJCI)=IIUNQ+IBCI                                   CPP07870
        ICPOS(IBCI,INUM)=IJCI                                           CPP07880
      END IF                                                            CPP07890
  119 CONTINUE                                                          CPP07900
      IF(ITEST.NE.0) GO TO 207                                          CPP07910
  206 CONTINUE                                                          CPP07920
      IJGRP=IJGRP+1                                                     CPP07930
      IBM(IJGRP)=ICM                                                    CPP07940
      JBM(IJGRP)=JCM                                                    CPP07950
      IBNUM(IJGRP)=1                                                    CPP07960
      ICITYP(IIUNQ+IJCI)=IIUNQ+IJGRP                                    CPP07970
      ICPOS(IJGRP,1)=IJCI                                               CPP07980
C                                                                       CPP07990
  207 CONTINUE                                                          CPP08000
      IF(IPRNT.LE.2) GO TO 120                                          CPP08010
      WRITE(6,14) IJCI,IEX,JEX,III,JJJ,IIMO(IJCI),JJMO(IJCI),           CPP08020
     *            ICIMO(III),ICIMO(JJJ),ICITYP(IIUNQ+IJCI)              CPP08030
  120 CONTINUE                                                          CPP08040
C                                                                       CPP08050
      WRITE(6,15) IJGRP                                                 CPP08060
      IF(IPRNT.LE.2) GO TO 210                                          CPP08070
      WRITE(6,16)                                                       CPP08080
      DO 122 IBCI=1,IJGRP                                               CPP08090
      ICM=IBM(IBCI)                                                     CPP08100
      JCM=JBM(IBCI)                                                     CPP08110
      NUM=IBNUM(IBCI)                                                   CPP08120
      WRITE(6,14) IBCI,ICM,JCM,NUM                                      CPP08130
      DO 121 INUM=1,NUM                                                 CPP08140
      WRITE(6,14) ICPOS(IBCI,INUM)                                      CPP08150
  121 CONTINUE                                                          CPP08160
  122 CONTINUE                                                          CPP08170
C                                                                       CPP08180
  210 CONTINUE                                                          CPP08190
      RETURN                                                            CPP08200
      END                                                               CPP08210
      SUBROUTINE BARES(FIJ,AIJ,BIJ)                                     CPP08220
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP08230
      DIMENSION FIJ(NTYPEP,ICIUNQ),AIJ(NTYTRI,ICIUNQ),BIJ(NTYTRI,ICIUNQ)CPP08240
      DIMENSION F(20),A(20,20),B(20,20)                                 CPP08250
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP08260
      COMMON/CIOFF/IIUNQ,IJUNQ,ICIUNQ,IPX(5050),JPX(5050)               CPP08270
      COMMON/CIPEX/PEXC(100),INDPEX(20,100),MOPEX(20),ICIMO(20)         CPP08280
      COMMON/CIPMO/IIMO(5050),JJMO(5050)                                CPP08290
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP08300
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP08310
      COMMON/DENST/ALPC(465),BETC(465)                                  CPP08320
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP08330
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP08340
      COMMON/CI103/CI(100)                                              CPP08350
      DATA A00,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /   CPP08360
    1 FORMAT(//,2X,' THE CI COEFFICIENTS'/)                             CPP08370
    2 FORMAT(2X,I5,F15.8)                                               CPP08380
    3 FORMAT(//,2X,' THE NORMALIZATION OF CI COEFFICIENTS',5X,F15.8/)   CPP08390
    4 FORMAT(//,2X,' FIJ VECTOR, ICI = ',I5/)                           CPP08400
    5 FORMAT(8F10.5)                                                    CPP08410
    6 FORMAT(//,2X,' AIJ MATRIX, ICI = ',I5/)                           CPP08420
    7 FORMAT(//,2X,' BIJ MATRIX, ICI = ',I5/)                           CPP08430
    8 FORMAT(//,2X,' FOCC VECTOR FROM FI  (NON-ZERO CI''S)'/5X,10F10.5/)CPP08440
    9 FORMAT(//,2X,' ALPA MATRIX FROM AIJ (NON-ZERO CI''S)'/)           CPP08450
   10 FORMAT(//,2X,' BETA MATRIX FROM BIJ (NON-ZERO CI''S)'/)           CPP08460
C                                                                       CPP08470
C   TRANSFER THE CI COEFFICEINTS                                        CPP08480
C                                                                       CPP08490
      IF(IPRNT.GT.2)                                                    CPP08500
     *WRITE(6,1)                                                        CPP08510
      CTOT=A00                                                          CPP08520
      DO 101 I=1,NCPEX                                                  CPP08530
      CI(I)=PEXC(I)                                                     CPP08540
      CTOT=CTOT+CI(I)*CI(I)                                             CPP08550
      IF(IPRNT.LT.2) GO TO 101                                          CPP08560
      WRITE(6,2) I,CI(I)                                                CPP08570
  101 CONTINUE                                                          CPP08580
      IF(IPRNT.GT.2)                                                    CPP08590
     *WRITE(6,3) CTOT                                                   CPP08600
C                                                                       CPP08610
C   FORM THE "BARE" COUPLING CONSTANT MATRICES                          CPP08620
      CALL ZERO(FIJ,NTYPEP*ICIUNQ)                                      CPP08630
      CALL ZERO(AIJ,NTYTRI*ICIUNQ)                                      CPP08640
      CALL ZERO(BIJ,NTYTRI*ICIUNQ)                                      CPP08650
C                                                                       CPP08660
C   ONE-ELECTRON COUPLING CONSTANTS                                     CPP08670
      DO 105 ICI=1,NCPEX                                                CPP08680
      DO 102 ITYP=1,NTYPES                                              CPP08690
      IF(FOCC(ITYP).NE.TWO.AND.FOCC(ITYP).NE.ONE) GO TO 102             CPP08700
      FIJ(ITYP,ICI)=FOCC(ITYP)*HALF                                     CPP08710
  102 CONTINUE                                                          CPP08720
      DO 103 I=1,NORPEX                                                 CPP08730
      IF(INDPEX(I,ICI).EQ.0) GO TO 103                                  CPP08740
      MTYP=MOTYP(ICIMO(I))                                              CPP08750
      FIJ(MTYP,ICI)=ONE                                                 CPP08760
  103 CONTINUE                                                          CPP08770
      IF(IPRNT.LE.2) GO TO 105                                          CPP08780
      WRITE(6,4) ICI                                                    CPP08790
      WRITE(6,5) (FIJ(ITYP,ICI),ITYP=1,NTYPEP)                          CPP08800
  105 CONTINUE                                                          CPP08810
C                                                                       CPP08820
C   TWO-ELECTRON COUPLING CONSTANTS                                     CPP08830
C   FOR THE CI DIAGONAL ELEMENTS                                        CPP08840
      DO 115 ICI=1,NCPEX                                                CPP08850
C   FOR CLOSED---OPEN-SHELL                                             CPP08860
      DO 106 ITYP=1,NTYPES                                              CPP08870
      DO 106 JTYP=1,ITYP                                                CPP08880
      IJTYP=IOFF(MAX0(ITYP,JTYP))+MIN0(ITYP,JTYP)                       CPP08890
      IF(FOCC(ITYP).NE.TWO.AND.FOCC(ITYP).NE.ONE) GO TO 106             CPP08900
      IF(FOCC(JTYP).NE.TWO.AND.FOCC(JTYP).NE.ONE) GO TO 106             CPP08910
      AIJ(IJTYP,ICI)=ALPC(IJTYP)                                        CPP08920
      BIJ(IJTYP,ICI)=BETC(IJTYP)                                        CPP08930
  106 CONTINUE                                                          CPP08940
C   FOR CLOSED AND OPEN-SHELLS---PEX SHELLS                             CPP08950
      DO 108 ITYP=1,NTYPES                                              CPP08960
      IF(FOCC(ITYP).NE.TWO.AND.FOCC(ITYP).NE.ONE) GO TO 108             CPP08970
      DO 107 I=1,NORPEX                                                 CPP08980
      IF(INDPEX(I,ICI).EQ.0) GO TO 107                                  CPP08990
      MTYP=MOTYP(ICIMO(I))                                              CPP09000
      IMTYP=IOFF(MAX0(ITYP,MTYP))+MIN0(ITYP,MTYP)                       CPP09010
      AIJ(IMTYP,ICI)=FOCC(ITYP)                                         CPP09020
      BIJ(IMTYP,ICI)=-FOCC(ITYP)*HALF                                   CPP09030
  107 CONTINUE                                                          CPP09040
  108 CONTINUE                                                          CPP09050
C   FOR PEX-PEX SHELLS                                                  CPP09060
      DO 110 I=1,NORPEX                                                 CPP09070
      IF(INDPEX(I,ICI).EQ.0) GO TO 110                                  CPP09080
      MTYP=MOTYP(ICIMO(I))                                              CPP09090
      DO 109 J=1,I                                                      CPP09100
      IF(INDPEX(J,ICI).EQ.0) GO TO 109                                  CPP09110
      NTYP=MOTYP(ICIMO(J))                                              CPP09120
      MNTYP=IOFF(MAX0(MTYP,NTYP))+MIN0(MTYP,NTYP)                       CPP09130
      AIJ(MNTYP,ICI)=TWO                                                CPP09140
      BIJ(MNTYP,ICI)=-ONE                                               CPP09150
  109 CONTINUE                                                          CPP09160
  110 CONTINUE                                                          CPP09170
      IF(IPRNT.LE.2) GO TO 115                                          CPP09180
      WRITE(6,6) ICI                                                    CPP09190
      CALL PRINT(AIJ(1,ICI),NCTRI,NTYPEP,6)                             CPP09200
      WRITE(6,7) ICI                                                    CPP09210
      CALL PRINT(BIJ(1,ICI),NCTRI,NTYPEP,6)                             CPP09220
  115 CONTINUE                                                          CPP09230
C                                                                       CPP09240
C   TWO-ELECTRON COUPLING CONSTANTS                                     CPP09250
C   FOR THE CI OFF-DIAGONAL ELEMENTS (BIJ ONLY)                         CPP09260
      DO 116 IJCI=1,IJUNQ                                               CPP09270
      II=IIMO(IJCI)                                                     CPP09280
      JJ=JJMO(IJCI)                                                     CPP09290
      IIJJ=IOFF(MAX0(II,JJ))+MIN0(II,JJ)                                CPP09300
      IJEX=IIUNQ+IJCI                                                   CPP09310
      BIJ(IIJJ,IJEX)=HALF                                               CPP09320
      IF(IPRNT.LE.2) GO TO 116                                          CPP09330
      WRITE(6,7) IJCI                                                   CPP09340
      CALL PRINT(BIJ(1,IJEX),NCTRI,NTYPEP,6)                            CPP09350
  116 CONTINUE                                                          CPP09360
C                                                                       CPP09370
C   TEST THE "BARE" COUPLING CONSTANTS USING NON-ZERO CI'S              CPP09380
C*    IF(IPRNT.LE.2) GO TO 201                                          CPP09390
      CALL ZERO(F,20)                                                   CPP09400
      CALL ZERO(A,400)                                                  CPP09410
      CALL ZERO(B,400)                                                  CPP09420
      DO 120 IJCI=1,ICIUNQ                                              CPP09430
      IEX=IPX(IJCI)                                                     CPP09440
      JEX=JPX(IJCI)                                                     CPP09450
      FAC=TWO                                                           CPP09460
      IF(IEX.EQ.JEX) FAC=ONE                                            CPP09470
      DO 117 ITYP=1,NTYPES                                              CPP09480
      F(ITYP)=F(ITYP)+CI(IEX)*CI(JEX)*FIJ(ITYP,IJCI)*FAC                CPP09490
  117 CONTINUE                                                          CPP09500
      DO 118 ITYP=1,NTYPES                                              CPP09510
      DO 118 JTYP=1,NTYPES                                              CPP09520
      IJTYP=IOFF(MAX0(ITYP,JTYP))+MIN0(ITYP,JTYP)                       CPP09530
      A(ITYP,JTYP)=A(ITYP,JTYP)+CI(IEX)*CI(JEX)*AIJ(IJTYP,IJCI)*FAC     CPP09540
      B(ITYP,JTYP)=B(ITYP,JTYP)+CI(IEX)*CI(JEX)*BIJ(IJTYP,IJCI)*FAC     CPP09550
  118 CONTINUE                                                          CPP09560
  120 CONTINUE                                                          CPP09570
      WRITE(6,8) (F(I),I=1,NTYPEP)                                      CPP09580
      WRITE(6,9)                                                        CPP09590
      CALL MATOUT(A,20,20,NTYPEP,NTYPEP,6)                              CPP09600
      WRITE(6,10)                                                       CPP09610
      CALL MATOUT(B,20,20,NTYPEP,NTYPEP,6)                              CPP09620
C                                                                       CPP09630
  201 CONTINUE                                                          CPP09640
      RETURN                                                            CPP09650
      END                                                               CPP09660
      SUBROUTINE DERMAT(EAO,E1A,E2A,SS,U,T)                             CPP09670
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP09680
      DIMENSION E1A(N3N),E2A(N3N,N3N),SS(NBATRI)                        CPP09690
      DIMENSION EAO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)         CPP09700
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP09710
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP09720
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP09730
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP09740
      COMMON/TAPES/ITAP35,ITAP42,ITAP43,ITAP44,ITAP47                   CPP09750
      COMMON/WORKS/JTAPE1,JTAPE2                                        CPP09760
    1 FORMAT(//,2X,' SCF FIRST DERIVATIVE MATRIX'/)                     CPP09770
    2 FORMAT(//,2X,' SCF SECOND DERIVATIVE MATRIX'/)                    CPP09780
    3 FORMAT(//,2X,' SA (AO BASIS) MATRIX, IABC = ',I5/)                CPP09790
    4 FORMAT(//,2X,' SA (MO BASIS) MATRIX, IABC = ',I5/)                CPP09800
    5 FORMAT(//,2X,' HA (AO BASIS) MATRIX, IABC = ',I5/)                CPP09810
    6 FORMAT(//,2X,' HA (MO BASIS) MATRIX, IABC = ',I5/)                CPP09820
    7 FORMAT(//,2X,' TA (AO BASIS) MATRIX, IABC = ',I5,' ITYP = ',I5/)  CPP09830
    8 FORMAT(//,2X,' TA (MO BASIS) MATRIX, IABC = ',I5,' ITYP = ',I5/)  CPP09840
C                                                                       CPP09850
      CALL RFILE(ITAP42)                                                CPP09860
      CALL RFILE(ITAP44)                                                CPP09870
      CALL RFILE(JTAPE1)                                                CPP09880
      NBSET=NTREAD                                                      CPP09890
C                                                                       CPP09900
C   READ IN SCF FIRST DERIVATIVES                                       CPP09910
      CALL SREAD(ITAP42,E1A,N3N*2)                                      CPP09920
      WRITE(6,1)                                                        CPP09930
      CALL MATOUT(E1A,3,NATOMS,3,NATOMS,6)                              CPP09940
C                                                                       CPP09950
C   READ IN SCF SECOND DERIVATIVES                                      CPP09960
      CALL SREAD(ITAP42,E2A,N3N*N3N*2)                                  CPP09970
      WRITE(6,2)                                                        CPP09980
      CALL MATOUT(E2A,N3N,N3N,N3N,N3N,6)                                CPP09990
C                                                                       CPP10000
C   READ IN S FIRST DERIVATIVE INTEGRALS                                CPP10010
C   AND TRANSFORM THEM FROM AO TO MO BASIS                              CPP10020
      DO 101 IABC=1,N3N                                                 CPP10030
      IS=ISA(IABC)                                                      CPP10040
      CALL SREAD(ITAP42,SS,NTRI2)                                       CPP10050
      IF(IPRNT.LE.4) GO TO 201                                          CPP10060
      WRITE(6,3) IABC                                                   CPP10070
      CALL PRINT(SS,NBATRI,NBASIS,6)                                    CPP10080
  201 CONTINUE                                                          CPP10090
      CALL MOCONV(SS,NBATRI,SS,NBATRI,EAO,U,T)                          CPP10100
      CALL RWRIT(ITAP44,SS,NTRI2,IS)                                    CPP10110
      IF(IPRNT.LE.4) GO TO 101                                          CPP10120
      WRITE(6,4) IABC                                                   CPP10130
      CALL PRINT(SS,NBATRI,NBASIS,6)                                    CPP10140
  101 CONTINUE                                                          CPP10150
C                                                                       CPP10160
C   READ IN ONE ELECTRON FIRST DERIVATIVE INTEGRALS                     CPP10170
C   AND TRANSFORM THEM FROM AO TO MO BASIS                              CPP10180
      DO 102 IABC=1,N3N                                                 CPP10190
      IH=IHA(IABC)                                                      CPP10200
      CALL SREAD(ITAP42,SS,NTRI2)                                       CPP10210
      IF(IPRNT.LE.4) GO TO 202                                          CPP10220
      WRITE(6,5) IABC                                                   CPP10230
      CALL PRINT(SS,NBATRI,NBFAO,6)                                     CPP10240
  202 CONTINUE                                                          CPP10250
      CALL MOCONV(SS,NBATRI,SS,NBATRI,EAO,U,T)                          CPP10260
      CALL RWRIT(ITAP44,SS,NTRI2,IH)                                    CPP10270
      IF(IPRNT.LE.4) GO TO 102                                          CPP10280
      WRITE(6,6) IABC                                                   CPP10290
      CALL PRINT(SS,NBATRI,NBASIS,6)                                    CPP10300
  102 CONTINUE                                                          CPP10310
C                                                                       CPP10320
C   READ IN TWO ELECTRON FIRST DERIVATIVE INTEGRALS                     CPP10330
C   AND TRANSFORM THEM FROM AO TO MO BASIS                              CPP10340
      DO 104 ITYP=1,NBSET                                               CPP10350
      DO 103 IABC=1,N3N                                                 CPP10360
      CALL SREAD(ITAP42,SS,NTRI2)                                       CPP10370
      IF(IPRNT.LE.4) GO TO 203                                          CPP10380
      WRITE(6,7) ITYP,IABC                                              CPP10390
      CALL PRINT(SS,NBATRI,NBFAO,6)                                     CPP10400
  203 CALL MOCONV(SS,NBATRI,SS,NBATRI,EAO,U,T)                          CPP10410
      CALL SWRIT(JTAPE1,SS,NTRI2)                                       CPP10420
      IF(IPRNT.LE.4) GO TO 103                                          CPP10430
      WRITE(6,8) ITYP,IABC                                              CPP10440
      CALL PRINT(SS,NBATRI,NBASIS,6)                                    CPP10450
  103 CONTINUE                                                          CPP10460
  104 CONTINUE                                                          CPP10470
C                                                                       CPP10480
      CALL SREW(ITAP42)                                                 CPP10490
      CALL SREW(ITAP44)                                                 CPP10500
      CALL SREW(JTAPE1)                                                 CPP10510
      RETURN                                                            CPP10520
      END                                                               CPP10530
      SUBROUTINE MOCONV(SA,NNA,SM,NNM,EAO,U,T)                          CPP10540
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP10550
      DIMENSION SA(NNA),SM(NNM)                                         CPP10560
      DIMENSION EAO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)         CPP10570
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP10580
C                                                                       CPP10590
C   TRANSFORM INTEGRALS FROM AO TO MO BASIS                             CPP10600
      NBFASQ=NBFAO*NBFAO                                                CPP10610
      IIJJ=0                                                            CPP10620
      DO 101 II=1,NBFAO                                                 CPP10630
      DO 101 JJ=1,II                                                    CPP10640
      IIJJ=IIJJ+1                                                       CPP10650
      T(II,JJ)=SA(IIJJ)                                                 CPP10660
      T(JJ,II)=SA(IIJJ)                                                 CPP10670
  101 CONTINUE                                                          CPP10680
      CALL ZERO(U,NBFASQ)                                               CPP10690
      CALL MXMB(EAO,NBFAO,1,T,1,NBFAO,U,1,NBFAO,NBFAO,NBFAO,NBFAO)      CPP10700
      CALL ZERO(T,NBFASQ)                                               CPP10710
      CALL MXMB(U,1,NBFAO,EAO,1,NBFAO,T,1,NBFAO,NBFAO,NBFAO,NBFAO)      CPP10720
      IJ=0                                                              CPP10730
      DO 102 I=1,NBASIS                                                 CPP10740
      DO 102 J=1,I                                                      CPP10750
      IJ=IJ+1                                                           CPP10760
      SM(IJ)=T(I,J)                                                     CPP10770
  102 CONTINUE                                                          CPP10780
C                                                                       CPP10790
      RETURN                                                            CPP10800
      END                                                               CPP10810
      SUBROUTINE ABMAT(ALPA,BETA)                                       CPP10820
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP10830
      DIMENSION ALPA(NTRI),BETA(NTRI)                                   CPP10840
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP10850
      COMMON/DENST/ALPC(465),BETC(465)                                  CPP10860
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP10870
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP10880
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP10890
    1 FORMAT(//,2X,' ALPA MATRIX'/)                                     CPP10900
    2 FORMAT(//,2X,' BETA MATRIX'/)                                     CPP10910
C                                                                       CPP10920
C   SET UP ALPHA AND BETA MATRICES                                      CPP10930
      CALL ZERO(ALPA,NTRI)                                              CPP10940
      CALL ZERO(BETA,NTRI)                                              CPP10950
      NSTRI=0                                                           CPP10960
      NENDI=0                                                           CPP10970
      DO 103 ITYP=1,NTYPES                                              CPP10980
      MM=NSORB(ITYP)                                                    CPP10990
      IF(MM.EQ.0) GO TO 103                                             CPP11000
      NSTRI=NENDI+1                                                     CPP11010
      NENDI=NSTRI+MM-1                                                  CPP11020
      NSTRJ=0                                                           CPP11030
      NENDJ=0                                                           CPP11040
      DO 102 JTYP=1,ITYP                                                CPP11050
      NN=NSORB(JTYP)                                                    CPP11060
      IF(NN.EQ.0) GO TO 102                                             CPP11070
      NSTRJ=NENDJ+1                                                     CPP11080
      NENDJ=NSTRJ+NN-1                                                  CPP11090
      DO 101 I=NSTRI,NENDI                                              CPP11100
      DO 101 J=NSTRJ,NENDJ                                              CPP11110
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPP11120
      IIJJ=IOFF(MAX0(ITYP,JTYP))+MIN0(ITYP,JTYP)                        CPP11130
      ALPA(IJ)=ALPC(IIJJ)                                               CPP11140
      BETA(IJ)=BETC(IIJJ)                                               CPP11150
  101 CONTINUE                                                          CPP11160
  102 CONTINUE                                                          CPP11170
  103 CONTINUE                                                          CPP11180
C                                                                       CPP11190
      IF(IPRNT.LE.2) GO TO 201                                          CPP11200
      WRITE(6,1)                                                        CPP11210
      CALL PRINT(ALPA,NTRI,NBASIS,6)                                    CPP11220
      WRITE(6,2)                                                        CPP11230
      CALL PRINT(BETA,NTRI,NBASIS,6)                                    CPP11240
C                                                                       CPP11250
  201 RETURN                                                            CPP11260
      END                                                               CPP11270
      SUBROUTINE E0CAL(ALPA,BETA,H,CC,LBLI,BUFI,NNC)                    CPP11280
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP11290
      DIMENSION ALPA(NTRI),BETA(NTRI),H(NTRI),CC(NNC)                   CPP11300
      DIMENSION LBLI(MAXBF4),BUFI(MAXBF2)                               CPP11310
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP11320
      COMMON/ENRGY/ESCF,ENUC,ELEC                                       CPP11330
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPP11340
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP11350
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP11360
      COMMON/TAPES/ITAP35,ITAP42,ITAP43,ITAP44,ITAP47                   CPP11370
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP11380
      DATA A00,ONE,TWO / 0.0D+00 , 1.0D+00 , 2.0D+00 /                  CPP11390
    1 FORMAT(//,2X,' H MATRIX (MO BASIS)'/)                             CPP11400
    2 FORMAT(//,2X,' SCF ENERGIES'/                                     CPP11410
     * 2X,' ESCF       = ',F20.10/                                      CPP11420
     * 2X,' ENUC       = ',F20.10/                                      CPP11430
     * 2X,' ELEC-ONE   = ',F20.10/                                      CPP11440
     * 2X,' ELEC-TWO   = ',F20.10/                                      CPP11450
     * 2X,' ELEC-TOTAL = ',F20.10/                                      CPP11460
     * 2X,' ETOT       = ',F20.10)                                      CPP11470
C                                                                       CPP11480
C   CALCULATE SCF ENERGY FOR A TEST                                     CPP11490
C   SEE EQ.(4)                                                          CPP11500
C                                                                       CPP11510
      CALL MREAD(H,17)                                                  CPP11520
      IF(IPRNT.LE.3) GO TO 201                                          CPP11530
      WRITE(6,1)                                                        CPP11540
      CALL PRINT(H,NTRI,NBASIS,6)                                       CPP11550
C                                                                       CPP11560
C   READ IN TWO ELECTRON MO INTEGRALS                                   CPP11570
  201 CONTINUE                                                          CPP11580
      CALL RDINT2(CC,NNC,ITAP35,NNC,LBLI,BUFI,IPRNT)                    CPP11590
C                                                                       CPP11600
C   CALCULATE ONE-ELECTRON ENERGY                                       CPP11610
      ELEC1=A00                                                         CPP11620
      DO 101 I=1,JOCC                                                   CPP11630
      II=IOFF(I+1)                                                      CPP11640
      ITYP=MOTYP(I)                                                     CPP11650
      ELEC1=ELEC1+H(II)*FOCC(ITYP)                                      CPP11660
  101 CONTINUE                                                          CPP11670
C   CALCULATE TWO-ELECTRON ENRGY                                        CPP11680
      ELEC2=A00                                                         CPP11690
      DO 102 I=1,JOCC                                                   CPP11700
      DO 102 J=1,I                                                      CPP11710
      IJ=IOFF(I)+J                                                      CPP11720
      II=IOFF(I+1)                                                      CPP11730
      JJ=IOFF(J+1)                                                      CPP11740
      IIJJ=IOFF(II)+JJ                                                  CPP11750
      IJIJ=IOFF(IJ+1)                                                   CPP11760
      FAC=TWO                                                           CPP11770
      IF(I.EQ.J) FAC=ONE                                                CPP11780
      ELEC2=ELEC2+(ALPA(IJ)*CC(IIJJ)+BETA(IJ)*CC(IJIJ))*FAC             CPP11790
  102 CONTINUE                                                          CPP11800
C                                                                       CPP11810
      ELECT=ELEC1+ELEC2                                                 CPP11820
      ETOT=ENUC+ELECT                                                   CPP11830
      WRITE(6,2) ESCF,ENUC,ELEC1,ELEC2,ELECT,ETOT                       CPP11840
C                                                                       CPP11850
      RETURN                                                            CPP11860
      END                                                               CPP11870
      SUBROUTINE E0LAG(ALPA,BETA,EPS,ZETA,H,CC,LBLI,BUFI,NNC)           CPP11880
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP11890
      DIMENSION ALPA(NTRI),BETA(NTRI),H(NTRI),CC(NNC)                   CPP11900
      DIMENSION EPS(NTRI),ZETA(NTRI,NTYPEP)                             CPP11910
      DIMENSION LBLI(MAXBF4),BUFI(MAXBF2)                               CPP11920
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP11930
      COMMON/ENRGY/ESCF,ENUC,ELEC                                       CPP11940
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP11950
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPP11960
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP11970
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP11980
      COMMON/TAPES/ITAP35,ITAP42,ITAP43,ITAP44,ITAP47                   CPP11990
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP12000
      DATA A00,HALF,TWO / 0.0D+00 , 0.5D+00 , 2.0D+00 /                 CPP12010
    1 FORMAT(//,2X,' LAGRANGIAN MATRIX'/)                               CPP12020
    2 FORMAT(//,2X,' SCF ENERGIES'/                                     CPP12030
     * 2X,' ESCF       = ',F20.10/                                      CPP12040
     * 2X,' ENUC       = ',F20.10/                                      CPP12050
     * 2X,' ELEC-ONE   = ',F20.10/                                      CPP12060
     * 2X,' ELEC-TWO   = ',F20.10/                                      CPP12070
     * 2X,' ELEC-TOTAL = ',F20.10/                                      CPP12080
     * 2X,' ETOT       = ',F20.10)                                      CPP12090
    3 FORMAT(//,2X,' ZETA MATRIX  , ITYP = ',I5/)                       CPP12100
C                                                                       CPP12110
C   FORM LAGRANGIAN MATRIX AND CALCULATE SCF ENERGY FOR A TEST          CPP12120
C   SEE EQ.(11)                                                         CPP12130
C                                                                       CPP12140
      CALL MREAD(H,17)                                                  CPP12150
C   READ IN TWO ELECTRON MO INTEGRALS                                   CPP12160
      CALL RDINT2(CC,NNC,ITAP35,NNC,LBLI,BUFI,IPRNT)                    CPP12170
C                                                                       CPP12180
C   ONE-ELECTRON TERMS                                                  CPP12190
      DO 101 I=1,NBASIS                                                 CPP12200
      ITYP=MOTYP(I)                                                     CPP12210
      FAC=FOCC(ITYP)*HALF                                               CPP12220
      DO 101 J=1,I                                                      CPP12230
      IJ=IOFF(I)+J                                                      CPP12240
      EPS(IJ)=H(IJ)*FAC                                                 CPP12250
  101 CONTINUE                                                          CPP12260
      ELEC1=A00                                                         CPP12270
      DO 102 I=1,JOCC                                                   CPP12280
      II=IOFF(I+1)                                                      CPP12290
  102 ELEC1=ELEC1+EPS(II)*TWO                                           CPP12300
C                                                                       CPP12310
C   TWO-ELECTRON TERMS                                                  CPP12320
      DO 105 I=1,NBASIS                                                 CPP12330
      DO 105 J=1,I                                                      CPP12340
      IJ=IOFF(I)+J                                                      CPP12350
      VALU=A00                                                          CPP12360
      DO 104 K=1,JOCC                                                   CPP12370
      KK=IOFF(K+1)                                                      CPP12380
      IK=IOFF(MAX0(I,K))+MIN0(I,K)                                      CPP12390
      JK=IOFF(MAX0(J,K))+MIN0(J,K)                                      CPP12400
      IJKK=IOFF(MAX0(IJ,KK))+MIN0(IJ,KK)                                CPP12410
      IKJK=IOFF(MAX0(IK,JK))+MIN0(IK,JK)                                CPP12420
      VALU=VALU+ALPA(IK)*CC(IJKK)+BETA(IK)*CC(IKJK)                     CPP12430
  104 CONTINUE                                                          CPP12440
      EPS(IJ)=EPS(IJ)+VALU                                              CPP12450
  105 CONTINUE                                                          CPP12460
      IF(IPRNT.LE.2) GO TO 201                                          CPP12470
      WRITE(6,1)                                                        CPP12480
      CALL PRINT(EPS,NTRI,NBASIS,6)                                     CPP12490
  201 ELEC2=A00                                                         CPP12500
      DO 106 I=1,JOCC                                                   CPP12510
      II=IOFF(I+1)                                                      CPP12520
  106 ELECT=ELECT+EPS(II)                                               CPP12530
C                                                                       CPP12540
      ELEC2=ELECT-ELEC1*HALF                                            CPP12550
      ELECT=ELEC1+ELEC2                                                 CPP12560
      ETOT=ENUC+ELECT                                                   CPP12570
      WRITE(6,2) ESCF,ENUC,ELEC1,ELEC2,ELECT,ETOT                       CPP12580
C                                                                       CPP12590
C   CALCULATE ZETA MATRICES                                             CPP12600
C   SEE EQ.(16)                                                         CPP12610
C                                                                       CPP12620
      CALL ZERO(ZETA,NTRI*NTYPEP)                                       CPP12630
      DO 110 ITYP=1,NTYPES                                              CPP12640
      L=NSTR(ITYP)                                                      CPP12650
      FAC=FOCC(ITYP)*HALF                                               CPP12660
      DO 109 I=1,NBASIS                                                 CPP12670
      DO 109 J=1,I                                                      CPP12680
      IJ=IOFF(I)+J                                                      CPP12690
      VALU=H(IJ)*FAC                                                    CPP12700
      DO 108 K=1,JOCC                                                   CPP12710
      KK=IOFF(K+1)                                                      CPP12720
      IK=IOFF(MAX0(I,K))+MIN0(I,K)                                      CPP12730
      JK=IOFF(MAX0(J,K))+MIN0(J,K)                                      CPP12740
      IJKK=IOFF(MAX0(IJ,KK))+MIN0(IJ,KK)                                CPP12750
      IKJK=IOFF(MAX0(IK,JK))+MIN0(IK,JK)                                CPP12760
      LK=IOFF(MAX0(L,K))+MIN0(L,K)                                      CPP12770
      VALU=VALU+ALPA(LK)*CC(IJKK)+BETA(LK)*CC(IKJK)                     CPP12780
  108 CONTINUE                                                          CPP12790
      ZETA(IJ,ITYP)=VALU                                                CPP12800
  109 CONTINUE                                                          CPP12810
      IF(IPRNT.LE.2) GO TO 110                                          CPP12820
      WRITE(6,3) ITYP                                                   CPP12830
      CALL PRINT(ZETA(1,ITYP),NTRI,NBASIS,6)                            CPP12840
  110 CONTINUE                                                          CPP12850
C                                                                       CPP12860
      RETURN                                                            CPP12870
      END                                                               CPP12880
      SUBROUTINE MODER(ALPA,BETA,D1W,E1T,EPS,D1O,D1T,SM,HM,TM)          CPP12890
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP12900
      DIMENSION ALPA(NTRI),BETA(NTRI),EPS(NTRI)                         CPP12910
      DIMENSION D1W(N3N),E1T(N3N)                                       CPP12920
      DIMENSION D1O(N3N),D1T(N3N)                                       CPP12930
      DIMENSION SM(NTRI,N3N),HM(NTRI,N3N),TM(NTRI,N3N,NTREAD)           CPP12940
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP12950
      COMMON/DENST/ALPC(465),BETC(465)                                  CPP12960
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP12970
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP12980
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP12990
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP13000
      COMMON/TAPES/ITAP35,ITAP42,ITAP43,ITAP44,ITAP47                   CPP13010
      COMMON/WORKS/JTAPE1,JTAPE2                                        CPP13020
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP13030
      DATA A00,ONE,TWO / 0.0D+00 , 1.0D+00 , 2.0D+00 /                  CPP13040
    1 FORMAT(//,2X,' D1N MATRIX IN MODER'/)                             CPP13050
    2 FORMAT(//,2X,' D1O MATRIX IN MODER'/)                             CPP13060
    3 FORMAT(//,2X,' D1T MATRIX IN MODER'/)                             CPP13070
    4 FORMAT(//,2X,' D1W MATRIX IN MODER'/)                             CPP13080
    5 FORMAT(//,2X,' EGRAD MATRIX IN MODER'/)                           CPP13090
C                                                                       CPP13100
C   CALCULATE SCF FIRST DERIVATIVES FOR A TEST                          CPP13110
C   SEE EQ.(9)                                                          CPP13120
C                                                                       CPP13130
C   READ BACK DERIVATIVE SM AND HM MATRICES                             CPP13140
      CALL SREW(ITAP44)                                                 CPP13150
      DO 101 IABC=1,N3N                                                 CPP13160
      IS=ISA(IABC)                                                      CPP13170
      IH=IHA(IABC)                                                      CPP13180
      CALL RREAD(ITAP44,SM(1,IABC),NTRI2,IS)                            CPP13190
      CALL RREAD(ITAP44,HM(1,IABC),NTRI2,IH)                            CPP13200
  101 CONTINUE                                                          CPP13210
C                                                                       CPP13220
C   READ BACK DERIVATIVE TM MATRICES                                    CPP13230
      CALL SREW(JTAPE1)                                                 CPP13240
      DO 103 ITYP=1,NTREAD                                              CPP13250
      DO 102 IABC=1,N3N                                                 CPP13260
      CALL SREAD(JTAPE1,TM(1,IABC,ITYP),NTRI2)                          CPP13270
  102 CONTINUE                                                          CPP13280
  103 CONTINUE                                                          CPP13290
C                                                                       CPP13300
      DO 110 IABC=1,N3N                                                 CPP13310
C  ONE-ELECTRON CONTRIBUTION                                            CPP13320
      VALU=A00                                                          CPP13330
      DO 106 I=1,JOCC                                                   CPP13340
      II=IOFF(I+1)                                                      CPP13350
      ITYP=MOTYP(I)                                                     CPP13360
      VALU=VALU+HM(II,IABC)*FOCC(ITYP)                                  CPP13370
  106 CONTINUE                                                          CPP13380
      D1O(IABC)=VALU                                                    CPP13390
C  TWO-ELECTRON CONTRIBUTION                                            CPP13400
      VALU=A00                                                          CPP13410
      DO 108 ITYP=1,NTYPES                                              CPP13420
      I=NSTR(ITYP)                                                      CPP13430
      DO 108 JTYP=1,NTYPES                                              CPP13440
      DO 107 J=NSTR(JTYP),NEND(JTYP)                                    CPP13450
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPP13460
      JJ=IOFF(J+1)                                                      CPP13470
      VALU=VALU+ALPA(IJ)*TM(JJ,IABC,ITYP*2-1)                           CPP13480
     1         +BETA(IJ)*TM(JJ,IABC,ITYP*2)                             CPP13490
  107 CONTINUE                                                          CPP13500
  108 CONTINUE                                                          CPP13510
      D1T(IABC)=VALU                                                    CPP13520
C   OVERLAP CONTRIBUTION                                                CPP13530
      VALU=A00                                                          CPP13540
      IJ=0                                                              CPP13550
      DO 109 I=1,JOCC                                                   CPP13560
      DO 109 J=1,I                                                      CPP13570
      IJ=IJ+1                                                           CPP13580
      FAC=TWO                                                           CPP13590
      IF(I.EQ.J) FAC=ONE                                                CPP13600
      VALU=VALU-SM(IJ,IABC)*EPS(IJ)*FAC                                 CPP13610
  109 CONTINUE                                                          CPP13620
      D1W(IABC)=VALU*TWO                                                CPP13630
C                                                                       CPP13640
C   SUM UP ALL CONTRIBUTIONS                                            CPP13650
      E1T(IABC)=D1O(IABC)+D1T(IABC)+D1W(IABC)                           CPP13660
C                                                                       CPP13670
  110 CONTINUE                                                          CPP13680
C                                                                       CPP13690
CCC   WRITE(6,1)                                                        CPP13700
CCC   CALL MATOUT(E1N,3,NATOMS,3,NATOMS,6)                              CPP13710
      WRITE(6,2)                                                        CPP13720
      CALL MATOUT(D1O,3,NATOMS,3,NATOMS,6)                              CPP13730
      WRITE(6,3)                                                        CPP13740
      CALL MATOUT(D1T,3,NATOMS,3,NATOMS,6)                              CPP13750
      WRITE(6,4)                                                        CPP13760
      CALL MATOUT(D1W,3,NATOMS,3,NATOMS,6)                              CPP13770
      WRITE(6,5)                                                        CPP13780
      CALL MATOUT(E1T,3,NATOMS,3,NATOMS,6)                              CPP13790
C                                                                       CPP13800
      CALL SREW(ITAP44)                                                 CPP13810
      CALL SREW(JTAPE1)                                                 CPP13820
      RETURN                                                            CPP13830
      END                                                               CPP13840
      SUBROUTINE HAMAT(FIJ,AIJ,BIJ,D1W,E1T,HCIA,D1E,HM,TM,NABC,N3NXX)   CPP13850
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP13860
      DIMENSION FIJ(NTYPEP,ICIUNQ),AIJ(NTYTRI,ICIUNQ),BIJ(NTYTRI,ICIUNQ)CPP13870
      DIMENSION HCIA(NCTRI,N3NXX)                                       CPP13880
      DIMENSION D1W(N3N),E1T(N3N),D1E(N3N)                              CPP13890
      DIMENSION HM(NTRI,N3N),TM(NTRI,NABC)                              CPP13900
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP13910
      COMMON/CIOFF/IIUNQ,IJUNQ,ICIUNQ,IPX(5050),JPX(5050)               CPP13920
      COMMON/CIPOS/ICITYP(5050),ICPOS(200,100)                          CPP13930
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP13940
      COMMON/CITYP/IJGRP,IBM(200),JBM(200),IBNUM(200)                   CPP13950
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP13960
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP13970
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP13980
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP13990
      COMMON/TAPES/ITAP35,ITAP42,ITAP43,ITAP44,ITAP47                   CPP14000
      COMMON/WORKS/JTAPE1,JTAPE2                                        CPP14010
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP14020
      COMMON/CI103/CI(100)                                              CPP14030
      DATA A00,ONE,TWO / 0.0D+00 , 1.0D+00 , 2.0D+00 /                  CPP14040
    1 FORMAT(//,2X,' HA MATRIX'/)                                       CPP14050
    2 FORMAT(2X,' IJCI,IABC,ELEC1,ELEC2 = ',2I5,2F15.8)                 CPP14060
    3 FORMAT(//,2X,' HCIA MATRIX'/)                                     CPP14070
    4 FORMAT(//,2X,' D1E MATRIX IN HAMAT'/)                             CPP14080
    5 FORMAT(//,2X,' D1W MATRIX IN HAMAT'/)                             CPP14090
    6 FORMAT(//,2X,' EGRAD MATRIX IN HAMAT'/)                           CPP14100
C                                                                       CPP14110
C   CALL BACK HM AND TM MATRICES AND FORM HA MATRICES                   CPP14120
C   SEE EQS.(33)-(35)                                                   CPP14130
C                                                                       CPP14140
      CALL SREW(ITAP44)                                                 CPP14150
      DO 101 IABC=1,N3N                                                 CPP14160
      IH=IHA(IABC)                                                      CPP14170
      CALL RREAD(ITAP44,HM(1,IABC),NTRI2,IH)                            CPP14180
  101 CONTINUE                                                          CPP14190
C                                                                       CPP14200
      CALL SREW(JTAPE1)                                                 CPP14210
      DO 102 ITYPX=1,NABC                                               CPP14220
      CALL SREAD(JTAPE1,TM(1,ITYPX),NTRI2)                              CPP14230
  102 CONTINUE                                                          CPP14240
C                                                                       CPP14250
      CALL ZERO(HCIA,NCTRI*N3NXX)                                       CPP14260
C                                                                       CPP14270
C     IF(IPRNT.GT.2)                                                    CPP14280
C     WRITE(6,1)                                                        CPP14290
C   FOR CI DIAGONAL ELEMENTS                                            CPP14300
      DO 110 IJCI=1,IIUNQ                                               CPP14310
      IEX=IPX(IJCI)                                                     CPP14320
      JEX=JPX(IJCI)                                                     CPP14330
      IJEX=IOFF(MAX0(IEX,JEX))+MIN(IEX,JEX)                             CPP14340
      DO 109 IABC=1,N3N                                                 CPP14350
C   ONE-ELECTRON CONTRIBUTION                                           CPP14360
      ELEC1=A00                                                         CPP14370
      DO 103 I=1,JOCC                                                   CPP14380
      II=IOFF(I+1)                                                      CPP14390
      ITYP=MOTYP(I)                                                     CPP14400
      ELEC1=ELEC1+FIJ(ITYP,IJCI)*HM(II,IABC)*TWO                        CPP14410
  103 CONTINUE                                                          CPP14420
C   TWO-ELECTRON CONTRIBUTION                                           CPP14430
      ELEC2=A00                                                         CPP14440
      DO 105 ITYP=1,NTYPES                                              CPP14450
      JAPT=(ITYP*2-2)*N3N+IABC                                          CPP14460
      KAPT=(ITYP*2-1)*N3N+IABC                                          CPP14470
      DO 105 JTYP=1,NTYPES                                              CPP14480
      IJTYP=IOFF(MAX0(ITYP,JTYP))+MIN0(ITYP,JTYP)                       CPP14490
      DO 104 J=NSTR(JTYP),NEND(JTYP)                                    CPP14500
      JJ=IOFF(J+1)                                                      CPP14510
      ELEC2=ELEC2+AIJ(IJTYP,IJCI)*TM(JJ,JAPT)                           CPP14520
     1           +BIJ(IJTYP,IJCI)*TM(JJ,KAPT)                           CPP14530
  104 CONTINUE                                                          CPP14540
  105 CONTINUE                                                          CPP14550
C                                                                       CPP14560
      HCIA(IJEX,IABC)=ELEC1+ELEC2                                       CPP14570
      IF(IPRNT.GT.3)                                                    CPP14580
     *WRITE(6,2) IJCI,IABC,ELEC1,ELEC2                                  CPP14590
  109 CONTINUE                                                          CPP14600
  110 CONTINUE                                                          CPP14610
C                                                                       CPP14620
C   FOR CI OFF-DIAGONAL ELEMENTS                                        CPP14630
      DO 115 IBCI=1,IJGRP                                               CPP14640
      ICM=IBM(IBCI)                                                     CPP14650
      JCM=JBM(IBCI)                                                     CPP14660
      NUM=IBNUM(IBCI)                                                   CPP14670
      ITYP=MOTYP(ICM)                                                   CPP14680
      JJ=IOFF(JCM+1)                                                    CPP14690
      DO 112 IABC=1,N3N                                                 CPP14700
      KAPT=(ITYP*2-1)*N3N+IABC                                          CPP14710
      ELEC2=TM(JJ,KAPT)                                                 CPP14720
      DO 111 INUM=1,NUM                                                 CPP14730
      IJCI=ICPOS(IBCI,INUM)                                             CPP14740
      IEX=IPX(IIUNQ+IJCI)                                               CPP14750
      JEX=JPX(IIUNQ+IJCI)                                               CPP14760
      IJEX=IOFF(MAX0(IEX,JEX))+MIN0(IEX,JEX)                            CPP14770
      HCIA(IJEX,IABC)=ELEC2                                             CPP14780
  111 CONTINUE                                                          CPP14790
  112 CONTINUE                                                          CPP14800
  115 CONTINUE                                                          CPP14810
      IF(IPRNT.LE.2) GO TO 201                                          CPP14820
      WRITE(6,3)                                                        CPP14830
      CALL MATOUT(HCIA,NCTRI,N3NXX,NCTRI,N3N,6)                         CPP14840
C                                                                       CPP14850
C   CALCULATE SCF DERIVATIVES FOR A TEST                                CPP14860
  201 CONTINUE                                                          CPP14870
      CALL ZERO(D1E,N3N)                                                CPP14880
      DO 120 IABC=1,N3N                                                 CPP14890
      DO 119 IJCI=1,ICIUNQ                                              CPP14900
      IEX=IPX(IJCI)                                                     CPP14910
      JEX=JPX(IJCI)                                                     CPP14920
      IJEX=IOFF(MAX0(IEX,JEX))+MIN(IEX,JEX)                             CPP14930
      FAC=TWO                                                           CPP14940
      IF(IEX.EQ.JEX) FAC=ONE                                            CPP14950
      CIFAC=CI(IEX)*CI(JEX)*FAC                                         CPP14960
      D1E(IABC)=D1E(IABC)+HCIA(IJEX,IABC)*CIFAC                         CPP14970
  119 CONTINUE                                                          CPP14980
      E1T(IABC)=D1E(IABC)+D1W(IABC)                                     CPP14990
  120 CONTINUE                                                          CPP15000
C***  IF(IPRNT.LE.2) GO TO 202                                          CPP15010
      WRITE(6,4)                                                        CPP15020
      CALL MATOUT(D1E,3,NATOMS,3,NATOMS,6)                              CPP15030
      WRITE(6,5)                                                        CPP15040
      CALL MATOUT(D1W,3,NATOMS,3,NATOMS,6)                              CPP15050
      WRITE(6,6)                                                        CPP15060
      CALL MATOUT(E1T,3,NATOMS,3,NATOMS,6)                              CPP15070
C                                                                       CPP15080
  202 CONTINUE                                                          CPP15090
      CALL SREW(ITAP44)                                                 CPP15100
      CALL SREW(JTAPE1)                                                 CPP15110
      RETURN                                                            CPP15120
      END                                                               CPP15130
      SUBROUTINE FAMAT(ALPA,BETA,FF,EPA,HM,TM,NABC)                     CPP15140
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP15150
      DIMENSION ALPA(NTRI),BETA(NTRI)                                   CPP15160
      DIMENSION FF(NTRI),HM(NTRI,N3N),TM(NTRI,NABC)                     CPP15170
      DIMENSION EPA(NBASIS,NBASIS)                                      CPP15180
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP15190
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP15200
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP15210
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP15220
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP15230
      COMMON/TAPES/ITAP35,ITAP42,ITAP43,ITAP44,ITAP47                   CPP15240
      COMMON/WORKS/JTAPE1,JTAPE2                                        CPP15250
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP15260
      DATA HALF / 0.5D+00 /                                             CPP15270
    1 FORMAT(//,2X,' EPA MATRIX, IABC = ',I5/)                          CPP15280
C                                                                       CPP15290
C   CALL BACK HM AND TM MATRICES AND FORM FA MATRICES                   CPP15300
C   SEE EQS.(15) AND (32)                                               CPP15310
C                                                                       CPP15320
      CALL SREW(ITAP44)                                                 CPP15330
      CALL SREW(JTAPE1)                                                 CPP15340
C                                                                       CPP15350
      DO 101 IABC=1,N3N                                                 CPP15360
      IH=IHA(IABC)                                                      CPP15370
      CALL RREAD(ITAP44,HM(1,IABC),NTRI2,IH)                            CPP15380
  101 CONTINUE                                                          CPP15390
C                                                                       CPP15400
      DO 102 ITYPX=1,NABC                                               CPP15410
      CALL SREAD(JTAPE1,TM(1,ITYPX),NTRI2)                              CPP15420
  102 CONTINUE                                                          CPP15430
C                                                                       CPP15440
C   STORE DERIVATIVE LAGRANGIAN MATRICES                                CPP15450
      DO 110 IABC=1,N3N                                                 CPP15460
      CALL ZERO(EPA,NBASIS*NBASIS)                                      CPP15470
      IE=IEA(IABC)                                                      CPP15480
      DO 105 I=1,JOCC                                                   CPP15490
      ITYP=MOTYP(I)                                                     CPP15500
      FAC=FOCC(ITYP)*HALF                                               CPP15510
      DO 105 J=1,NBASIS                                                 CPP15520
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPP15530
      VALU=HM(IJ,IABC)*FAC                                              CPP15540
      DO 104 KTYP=1,NTYPES                                              CPP15550
      K=NSTR(KTYP)                                                      CPP15560
      JAPT=(KTYP*2-2)*N3N+IABC                                          CPP15570
      KAPT=(KTYP*2-1)*N3N+IABC                                          CPP15580
      IK=IOFF(MAX0(I,K))+MIN0(I,K)                                      CPP15590
      VALU=VALU+ALPA(IK)*TM(IJ,JAPT)+BETA(IK)*TM(IJ,KAPT)               CPP15600
  104 CONTINUE                                                          CPP15610
      EPA(I,J)=VALU                                                     CPP15620
  105 CONTINUE                                                          CPP15630
C                                                                       CPP15640
      CALL RWRIT(ITAP44,EPA,NBASQ,IE)                                   CPP15650
      IF(IPRNT.LE.4) GO TO 110                                          CPP15660
      WRITE(6,1) IABC                                                   CPP15670
      CALL MATOUT(EPA,NBASIS,NBASIS,NBASIS,NBASIS,6)                    CPP15680
  110 CONTINUE                                                          CPP15690
C                                                                       CPP15700
C   STORE DERIVATIVE AVERAGED FOCK                                      CPP15710
      DO 115 IABC=1,N3N                                                 CPP15720
      IT=IFA(IABC)                                                      CPP15730
      DO 114 I=1,NTRI                                                   CPP15740
  114 FF(I)=HM(I,IABC)+TM(I,IABC)                                       CPP15750
      CALL RWRIT(ITAP44,FF,NTRI2,IT)                                    CPP15760
  115 CONTINUE                                                          CPP15770
C                                                                       CPP15780
      CALL SREW(ITAP44)                                                 CPP15790
      CALL SREW(JTAPE1)                                                 CPP15800
      RETURN                                                            CPP15810
      END                                                               CPP15820
      FUNCTION GSAI(I,J,K,L,ALPA,BETA,CC,NNC)                           CPP15830
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP15840
      DIMENSION ALPA(NTRI),BETA(NTRI),CC(NNC)                           CPP15850
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP15860
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP15870
      DATA TWO / 2.0D+00 /                                              CPP15880
C                                                                       CPP15890
C FORM GSAI MATRIX                                                      CPP15900
C                                                                       CPP15910
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPP15920
      IL=IOFF(MAX0(I,L))+MIN0(I,L)                                      CPP15930
      JL=IOFF(MAX0(J,L))+MIN0(J,L)                                      CPP15940
      KL=IOFF(MAX0(K,L))+MIN0(K,L)                                      CPP15950
      IK=IOFF(MAX0(I,K))+MIN0(I,K)                                      CPP15960
      JK=IOFF(MAX0(J,K))+MIN0(J,K)                                      CPP15970
      IJKL=IOFF(MAX0(IJ,KL))+MIN0(IJ,KL)                                CPP15980
      IKJL=IOFF(MAX0(IK,JL))+MIN0(IK,JL)                                CPP15990
      ILJK=IOFF(MAX0(IL,JK))+MIN0(IL,JK)                                CPP16000
      GSAI=(ALPA(IK)-ALPA(JK))*CC(IJKL)*TWO                             CPP16010
     1    +(BETA(IK)-BETA(JK))*(CC(IKJL)+CC(ILJK))                      CPP16020
C                                                                       CPP16030
      RETURN                                                            CPP16040
      END                                                               CPP16050
      SUBROUTINE BAMAT(EIG,NIJ,KIJ,ALPA,BETA,ZETA,SM,B0,EPA,FA,CC,      CPP16060
     1                 LBLI,BUFI,NNC)                                   CPP16070
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP16080
      DIMENSION EIG(NBASIS),NIJ(NTRI,2),KIJ(NTRI,2)                     CPP16090
      DIMENSION ALPA(NTRI),BETA(NTRI),ZETA(NTRI,NTYPEP)                 CPP16100
C     DIMENSION SM(NTRI,N3N),B0(NTRI,N3N),EPA(NBASIS,NBASIS),FA(NTRI)   CPP16110
      DIMENSION SM(NTRI,N3N),B0(NTRI),EPA(NBASIS,NBASIS),FA(NTRI)       CPP16120
      DIMENSION CC(NNC),LBLI(MAXBF4),BUFI(MAXBF2)                       CPP16130
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP16140
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP16150
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPP16160
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP16170
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP16180
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP16190
      COMMON/TAPES/ITAP35,ITAP42,ITAP43,ITAP44,ITAP47                   CPP16200
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP16210
      COMMON/CI102/NIND,NDEP                                            CPP16220
      DATA A00,HALF / 0.0D+00 , 0.5D+00 /                               CPP16230
    1 FORMAT(//,2X,' FM MATRIX, IABC = ',I5/)                           CPP16240
    2 FORMAT(//,2X,' B0 MATRIX IN BAMAT'/)                              CPP16250
    3 FORMAT(2X,3I5,4F15.8)                                             CPP16260
C                                                                       CPP16270
C   FORM B0 MATRICES                                                    CPP16280
C   SEE EQ.(26)                                                         CPP16290
C                                                                       CPP16300
C   READ IN TWO ELECTRON MO INTEGRALS AGAIN                             CPP16310
      CALL RDINT2(CC,NNC,ITAP35,NNC,LBLI,BUFI,IPRNT)                    CPP16320
C                                                                       CPP16330
C   CALL BACK SM MATRICES                                               CPP16340
      CALL SREW(ITAP44)                                                 CPP16350
      DO 102 IABC=1,N3N                                                 CPP16360
      IS=ISA(IABC)                                                      CPP16370
      CALL RREAD(ITAP44,SM(1,IABC),NTRI2,IS)                            CPP16380
  102 CONTINUE                                                          CPP16390
C     WRITE(6,*) ' SM MATRIX IN BAMAT'                                  CPP16400
C     CALL MATOUT(SM,NTRI,N3N,NTRI,N3N,6)                               CPP16410
C     CALL ZERO(B0,NTRI*N3N)                                            CPP16420
C                                                                       CPP16430
C   ELEMENTS FOR INDEPENDENT PAIRS                                      CPP16440
      DO 120 IABC=1,N3N                                                 CPP16450
      CALL ZERO(B0,NTRI)                                                CPP16460
      IT=IFA(IABC)                                                      CPP16470
      IE=IEA(IABC)                                                      CPP16480
      IB=IBA(IABC)                                                      CPP16490
      CALL RREAD(ITAP44,FA,NTRI2,IT)                                    CPP16500
      CALL RREAD(ITAP44,EPA,NBASQ,IE)                                   CPP16510
      IF(IPRNT.LE.4) GO TO 301                                          CPP16520
      WRITE(6,1) IABC                                                   CPP16530
      CALL MATOUT(EPA,NBASIS,NBASIS,NBASIS,NBASIS,6)                    CPP16540
      WRITE(6,2) IABC                                                   CPP16550
  301 CONTINUE                                                          CPP16560
      DO 110 II=1,NIND                                                  CPP16570
      I=NIJ(II,1)                                                       CPP16580
      J=NIJ(II,2)                                                       CPP16590
      ITYP=MOTYP(I)                                                     CPP16600
      JTYP=MOTYP(J)                                                     CPP16610
      IJ=IOFF(I)+J                                                      CPP16620
      VALU1=EPA(I,J)-EPA(J,I)                                           CPP16630
      VALU2=A00                                                         CPP16640
      VALU3=A00                                                         CPP16650
C                                                                       CPP16660
      DO 107 K=1,NBASIS                                                 CPP16670
      LIM=MIN0(JOCC,K)                                                  CPP16680
      DO 107 L=1,LIM                                                    CPP16690
      KL=IOFF(K)+L                                                      CPP16700
      FAC2=GSAI(I,J,K,L,ALPA,BETA,CC,NNC)                               CPP16710
      IF(K.EQ.L) FAC2=FAC2*HALF                                         CPP16720
      VALU2=VALU2-SM(KL,IABC)*FAC2                                      CPP16730
      FAC3=A00                                                          CPP16740
      IF(K.NE.J) GO TO 205                                              CPP16750
      IL=IOFF(MAX0(I,L))+MIN0(I,L)                                      CPP16760
      FAC3=FAC3+ZETA(IL,ITYP)-ZETA(IL,JTYP)                             CPP16770
  205 IF(K.NE.I) GO TO 206                                              CPP16780
      JL=IOFF(MAX0(J,L))+MIN0(J,L)                                      CPP16790
      FAC3=FAC3-ZETA(JL,JTYP)+ZETA(JL,ITYP)                             CPP16800
  206 VALU3=VALU3-SM(KL,IABC)*FAC3                                      CPP16810
  107 CONTINUE                                                          CPP16820
      VALUT=VALU1+VALU2+VALU3                                           CPP16830
C     B0(IJ,IABC)=VALUT                                                 CPP16840
      B0(IJ)=VALUT                                                      CPP16850
      IF(IPRNT.LE.4) GO TO 110                                          CPP16860
      WRITE(6,3) II,I,J,VALU1,VALU2,VALU3,VALUT                         CPP16870
  110 CONTINUE                                                          CPP16880
C                                                                       CPP16890
C   ELEMENTS FOR DEPENDENT PAIRS                                        CPP16900
      DO 115 II=1,NDEP                                                  CPP16910
      I=KIJ(II,1)                                                       CPP16920
      J=KIJ(II,2)                                                       CPP16930
      IJ=IOFF(I)+J                                                      CPP16940
      VALU1=FA(IJ)-EIG(J)*SM(IJ,IABC)                                   CPP16950
      VALU2=A00                                                         CPP16960
      DO 114 K=1,JOCC                                                   CPP16970
      KTYP=MOTYP(K)                                                     CPP16980
      DO 114 L=1,K                                                      CPP16990
      KL=IOFF(K)+L                                                      CPP17000
      FAC=ZIJ(I,J,K,L,CC,NNC)*FOCC(KTYP)*HALF                           CPP17010
      IF(K.EQ.L) FAC=FAC*HALF                                           CPP17020
      VALU2=VALU2-SM(KL,IABC)*FAC                                       CPP17030
  114 CONTINUE                                                          CPP17040
C     B0(IJ,IABC)=VALU1+VALU2                                           CPP17050
      B0(IJ)=VALU1+VALU2                                                CPP17060
  115 CONTINUE                                                          CPP17070
C     CALL RWRIT(ITAP44,B0(1,IABC),NTRI2,IB)                            CPP17080
      CALL RWRIT(ITAP44,B0,NTRI2,IB)                                    CPP17090
  120 CONTINUE                                                          CPP17100
C                                                                       CPP17110
C*    IF(IPRNT.LE.3) GO TO 302                                          CPP17120
C     WRITE(6,2)                                                        CPP17130
C     CALL MATOUT(B0,NTRI,N3N,NTRI,N3N,6)                               CPP17140
C                                                                       CPP17150
  302 CONTINUE                                                          CPP17160
      CALL SREW(ITAP44)                                                 CPP17170
      RETURN                                                            CPP17180
      END                                                               CPP17190
      SUBROUTINE HEMAT(FIJ,AIJ,BIJ,HCI,EIJ,ZIJ,EE,H,CC,LBLI,BUFI,NNC)   CPP17200
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP17210
      DIMENSION FIJ(NTYPEP,ICIUNQ),AIJ(NTYTRI,ICIUNQ),BIJ(NTYTRI,ICIUNQ)CPP17220
      DIMENSION HCI(NCTRI)                                              CPP17230
      DIMENSION EIJ(NBASIS,NBASIS),ZIJ(NTRI),EE(NBASIS,NBASIS)          CPP17240
      DIMENSION H(NTRI),CC(NNC),LBLI(MAXBF4),BUFI(MAXBF2)               CPP17250
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP17260
      COMMON/CIOFF/IIUNQ,IJUNQ,ICIUNQ,IPX(5050),JPX(5050)               CPP17270
      COMMON/CIPOS/ICITYP(5050),ICPOS(200,100)                          CPP17280
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP17290
      COMMON/CITYP/IJGRP,IBM(200),JBM(200),IBNUM(200)                   CPP17300
      COMMON/ENRGY/ESCF,ENUC,ELEC                                       CPP17310
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP17320
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPP17330
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP17340
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP17350
      COMMON/TAPES/ITAP35,ITAP42,ITAP43,ITAP44,ITAP47                   CPP17360
      COMMON/WORKS/JTAPE1,JTAPE2                                        CPP17370
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP17380
      COMMON/CI103/CI(100)                                              CPP17390
      COMMON/CI104/IEIJ(5050)                                           CPP17400
      DATA A00,HALF,ONE,TWO,FOUR / 0.0D+00 , 0.5D+00 , 1.0D+00 ,        CPP17410
     *                             2.0D+00 , 4.0D+00 /                  CPP17420
    1 FORMAT(//,2X,' CI MATRIX'/)                                       CPP17430
    2 FORMAT(//,2X,' EIJ MATRIX, IJCI = ',I5/)                          CPP17440
    3 FORMAT(//,2X,' ZIJ MATRIX, IJCI = ',I5,' ITYP= ',I5/)             CPP17450
    4 FORMAT(//,2X,' EIJ MATRIX, IBCI = ',I5/)                          CPP17460
    5 FORMAT(//,2X,' ZIJ MATRIX, IBCI = ',I5,' ITYP= ',I5/)             CPP17470
    6 FORMAT(//,2X,' CI ENERGIES'/                                      CPP17480
     * 2X,' ESCF       = ',F20.10/                                      CPP17490
     * 2X,' ENUC       = ',F20.10/                                      CPP17500
     * 2X,' ELEC       = ',F20.10/                                      CPP17510
     * 2X,' ECITOT     = ',F20.10/)                                     CPP17520
    7 FORMAT(//,2X,' THE LAGRANGIAN MATRIX'/)                           CPP17530
C                                                                       CPP17540
C   CALCULATE H MATRIX                                                  CPP17550
C   SEE EQS.(5)-(7)                                                     CPP17560
C                                                                       CPP17570
      CALL MREAD(H,17)                                                  CPP17580
      CALL RDINT2(CC,NNC,ITAP35,NNC,LBLI,BUFI,IPRNT)                    CPP17590
C                                                                       CPP17600
C::::::::::::::::::::::::::::::::::::                                   CPP17610
C:::FORM THE CI HAMILTONIAN MATRIX:::                                   CPP17620
C::::::::::::::::::::::::::::::::::::                                   CPP17630
      CALL ZERO(HCI,NCTRI)                                              CPP17640
C                                                                       CPP17650
C   FOR CI DIAGONAL ELEMENTS                                            CPP17660
      DO 103 IJCI=1,IIUNQ                                               CPP17670
      IEX=IPX(IJCI)                                                     CPP17680
      JEX=JPX(IJCI)                                                     CPP17690
      IJEX=IOFF(MAX0(IEX,JEX))+MIN0(IEX,JEX)                            CPP17700
C                                                                       CPP17710
C   CALCULATE ONE-ELECTRON ENERGY                                       CPP17720
      ELEC1=A00                                                         CPP17730
      DO 101 I=1,JOCC                                                   CPP17740
      ITYP=MOTYP(I)                                                     CPP17750
      II=IOFF(I+1)                                                      CPP17760
      ELEC1=ELEC1+FIJ(ITYP,IJCI)*H(II)*TWO                              CPP17770
  101 CONTINUE                                                          CPP17780
C   CALCULATE TWO-ELECTRON ENRGY                                        CPP17790
      ELEC2=A00                                                         CPP17800
      DO 102 I=1,JOCC                                                   CPP17810
      ITYP=MOTYP(I)                                                     CPP17820
      II=IOFF(I+1)                                                      CPP17830
      DO 102 J=1,I                                                      CPP17840
      JTYP=MOTYP(J)                                                     CPP17850
      IJTYP=IOFF(MAX0(ITYP,JTYP))+MIN0(ITYP,JTYP)                       CPP17860
      IJ=IOFF(I)+J                                                      CPP17870
      JJ=IOFF(J+1)                                                      CPP17880
      IIJJ=IOFF(II)+JJ                                                  CPP17890
      IJIJ=IOFF(IJ+1)                                                   CPP17900
      FAC=TWO                                                           CPP17910
      IF(I.EQ.J) FAC=ONE                                                CPP17920
      ELEC2=ELEC2+(AIJ(IJTYP,IJCI)*CC(IIJJ)+BIJ(IJTYP,IJCI)*CC(IJIJ))   CPP17930
     1           *FAC                                                   CPP17940
  102 CONTINUE                                                          CPP17950
      HCI(IJEX)=ELEC1+ELEC2                                             CPP17960
  103 CONTINUE                                                          CPP17970
C                                                                       CPP17980
C   FOR CI OFF-DIAGONAL ELEMENTS                                        CPP17990
      DO 105 IBCI=1,IJGRP                                               CPP18000
      ICM=IBM(IBCI)                                                     CPP18010
      JCM=JBM(IBCI)                                                     CPP18020
      NUM=IBNUM(IBCI)                                                   CPP18030
      IJCM=IOFF(MAX0(ICM,JCM))+MIN0(ICM,JCM)                            CPP18040
      IJIJ=IOFF(IJCM+1)                                                 CPP18050
      ELEC2=CC(IJIJ)                                                    CPP18060
      DO 104 INUM=1,NUM                                                 CPP18070
      IJCI=ICPOS(IBCI,INUM)                                             CPP18080
      IEX=IPX(IIUNQ+IJCI)                                               CPP18090
      JEX=JPX(IIUNQ+IJCI)                                               CPP18100
      IJEX=IOFF(MAX0(IEX,JEX))+MIN0(IEX,JEX)                            CPP18110
      HCI(IJEX)=ELEC2                                                   CPP18120
  104 CONTINUE                                                          CPP18130
  105 CONTINUE                                                          CPP18140
C                                                                       CPP18150
      WRITE(6,1)                                                        CPP18160
      CALL PRINT(HCI,NCTRI,NCPEX,6)                                     CPP18170
C                                                                       CPP18180
C::::::::::::::::::::::::::::::::::::::::::::::::::::                   CPP18190
C:::FORM THE "BARE" LAGRANGIAN MATRICES           :::                   CPP18200
C:::AND THE GENERALIZED "BARE" LAGRANGIAN MATRICES:::                   CPP18210
C::::::::::::::::::::::::::::::::::::::::::::::::::::                   CPP18220
      CALL RFILE(JTAPE2)                                                CPP18230
      CALL RFILE(ITAP47)                                                CPP18240
C                                                                       CPP18250
C   FOR CI DIAGONAL ELEMENTS                                            CPP18260
      DO 115 IJCI=1,IIUNQ                                               CPP18270
      CALL ZERO(EIJ,NBASIS*NBASIS)                                      CPP18280
C                                                                       CPP18290
C////////////////////////////////////                                   CPP18300
C///THE "BARE" LAGRANGIAN MATRICES///                                   CPP18310
C////////////////////////////////////                                   CPP18320
      DO 107 I=1,JOCC                                                   CPP18330
      ITYP=MOTYP(I)                                                     CPP18340
      DO 107 J=1,NBASIS                                                 CPP18350
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPP18360
C   ONE-ELECTRON PART                                                   CPP18370
      ELEC1=FIJ(ITYP,IJCI)*H(IJ)                                        CPP18380
C   TWO-ELECTRON PART                                                   CPP18390
      ELEC2=A00                                                         CPP18400
      DO 106 K=1,JOCC                                                   CPP18410
      KTYP=MOTYP(K)                                                     CPP18420
      IKTYP=IOFF(MAX0(ITYP,KTYP))+MIN0(ITYP,KTYP)                       CPP18430
      IK=IOFF(MAX0(I,K))+MIN0(I,K)                                      CPP18440
      JK=IOFF(MAX0(J,K))+MIN0(J,K)                                      CPP18450
      KK=IOFF(K+1)                                                      CPP18460
      IJKK=IOFF(MAX0(IJ,KK))+MIN0(IJ,KK)                                CPP18470
      IKJK=IOFF(MAX0(IK,JK))+MIN0(IK,JK)                                CPP18480
      ELEC2=ELEC2+AIJ(IKTYP,IJCI)*CC(IJKK)+BIJ(IKTYP,IJCI)*CC(IKJK)     CPP18490
  106 CONTINUE                                                          CPP18500
      EIJ(I,J)=ELEC1+ELEC2                                              CPP18510
  107 CONTINUE                                                          CPP18520
C                                                                       CPP18530
      IPT=IEIJ(IJCI)                                                    CPP18540
      CALL RWRIT(JTAPE2,EIJ,NBASQ,IPT)                                  CPP18550
      IF(IPRNT.LE.2) GO TO 201                                          CPP18560
      WRITE(6,2) IJCI                                                   CPP18570
      CALL MATOUT(EIJ,NBASIS,NBASIS,NBASIS,NBASIS,6)                    CPP18580
C                                                                       CPP18590
C////////////////////////////////////////////////                       CPP18600
C///THE GENERALIZED "BARE" LAGRANIGAN MATRICES///                       CPP18610
C////////////////////////////////////////////////                       CPP18620
  201 CONTINUE                                                          CPP18630
      DO 110 ITYP=1,NTYPES                                              CPP18640
      CALL ZERO(ZIJ,NTRI)                                               CPP18650
      IJ=0                                                              CPP18660
      DO 109 I=1,NBASIS                                                 CPP18670
      DO 109 J=1,I                                                      CPP18680
      IJ=IJ+1                                                           CPP18690
C   ONE-ELECTRON PART                                                   CPP18700
      ELEC1=FIJ(ITYP,IJCI)*H(IJ)                                        CPP18710
C   TWO-ELECTRON PART                                                   CPP18720
      ELEC2=A00                                                         CPP18730
      DO 108 K=1,JOCC                                                   CPP18740
      KTYP=MOTYP(K)                                                     CPP18750
      IKTYP=IOFF(MAX0(ITYP,KTYP))+MIN0(ITYP,KTYP)                       CPP18760
      IK=IOFF(MAX0(I,K))+MIN0(I,K)                                      CPP18770
      JK=IOFF(MAX0(J,K))+MIN0(J,K)                                      CPP18780
      KK=IOFF(K+1)                                                      CPP18790
      IJKK=IOFF(MAX0(IJ,KK))+MIN0(IJ,KK)                                CPP18800
      IKJK=IOFF(MAX0(IK,JK))+MIN0(IK,JK)                                CPP18810
      ELEC2=ELEC2+AIJ(IKTYP,IJCI)*CC(IJKK)+BIJ(IKTYP,IJCI)*CC(IKJK)     CPP18820
  108 CONTINUE                                                          CPP18830
      ZIJ(IJ)=ELEC1+ELEC2                                               CPP18840
  109 CONTINUE                                                          CPP18850
C                                                                       CPP18860
      CALL SWRIT(ITAP47,ZIJ,NTRI2)                                      CPP18870
      IF(IPRNT.LE.2) GO TO 110                                          CPP18880
      WRITE(6,3) IJCI,ITYP                                              CPP18890
      CALL PRINT(ZIJ,NTRI,NBASIS,6)                                     CPP18900
  110 CONTINUE                                                          CPP18910
  115 CONTINUE                                                          CPP18920
C                                                                       CPP18930
C   FOR CI OFF-DIAGONAL ELEMENTS                                        CPP18940
      DO 120 IBCI=1,IJGRP                                               CPP18950
      IP=IBM(IBCI)                                                      CPP18960
      IQ=JBM(IBCI)                                                      CPP18970
      NUM=IBNUM(IBCI)                                                   CPP18980
      IIPQ=IOFF(MAX0(IP,IQ))+MIN0(IP,IQ)                                CPP18990
C                                                                       CPP19000
C////////////////////////////////////                                   CPP19010
C///THE "BARE" LAGRANGIAN MATRICES///                                   CPP19020
C////////////////////////////////////                                   CPP19030
      CALL ZERO(EIJ,NBASIS*NBASIS)                                      CPP19040
      DO 116 J=1,NBASIS                                                 CPP19050
      JJQ=IOFF(MAX0(J,IQ))+MIN0(J,IQ)                                   CPP19060
      JJP=IOFF(MAX0(J,IP))+MIN0(J,IP)                                   CPP19070
      NPQJQ=IOFF(MAX0(IIPQ,JJQ))+MIN0(IIPQ,JJQ)                         CPP19080
      NQPJP=IOFF(MAX0(IIPQ,JJP))+MIN0(IIPQ,JJP)                         CPP19090
      EIJ(IP,J)=CC(NPQJQ)*HALF                                          CPP19100
      EIJ(IQ,J)=CC(NQPJP)*HALF                                          CPP19110
  116 CONTINUE                                                          CPP19120
      IPT=IEIJ(IIUNQ+IBCI)                                              CPP19130
      CALL RWRIT(JTAPE2,EIJ,NBASQ,IPT)                                  CPP19140
      IF(IPRNT.LE.2) GO TO 202                                          CPP19150
      WRITE(6,4) IBCI                                                   CPP19160
      CALL MATOUT(EIJ,NBASIS,NBASIS,NBASIS,NBASIS,6)                    CPP19170
C                                                                       CPP19180
C////////////////////////////////////////////////                       CPP19190
C///THE GENERALIZED "BARE" LAGRANIGAN MATRICES///                       CPP19200
C////////////////////////////////////////////////                       CPP19210
  202 CONTINUE                                                          CPP19220
      DO 119 ITYP=1,NTYPES                                              CPP19230
      CALL ZERO(ZIJ,NTRI)                                               CPP19240
      IF(ITYP.NE.MOTYP(IP).AND.ITYP.NE.MOTYP(IQ)) GO TO 203             CPP19250
      IF(ITYP.EQ.MOTYP(IP)) IPQ=IQ                                      CPP19260
      IF(ITYP.EQ.MOTYP(IQ)) IPQ=IP                                      CPP19270
C                                                                       CPP19280
      IJ=0                                                              CPP19290
      DO 118 I=1,NBASIS                                                 CPP19300
      DO 118 J=1,I                                                      CPP19310
      IJ=IJ+1                                                           CPP19320
      IIPQ=IOFF(MAX0(I,IPQ))+MIN0(I,IPQ)                                CPP19330
      JIPQ=IOFF(MAX0(J,IPQ))+MIN0(J,IPQ)                                CPP19340
      IPQJPQ=IOFF(MAX0(IIPQ,JIPQ))+MIN0(IIPQ,JIPQ)                      CPP19350
      ZIJ(IJ)=CC(IPQJPQ)*HALF                                           CPP19360
  118 CONTINUE                                                          CPP19370
C                                                                       CPP19380
  203 CONTINUE                                                          CPP19390
      CALL SWRIT(ITAP47,ZIJ,NTRI2)                                      CPP19400
      IF(IPRNT.LE.2) GO TO 119                                          CPP19410
      WRITE(6,5) IBCI,ITYP                                              CPP19420
      CALL PRINT(ZIJ,NTRI,NBASIS,6)                                     CPP19430
  119 CONTINUE                                                          CPP19440
C                                                                       CPP19450
  120 CONTINUE                                                          CPP19460
C                                                                       CPP19470
C   CALCULATE SCF ENERGY FOR A TEST                                     CPP19480
      CALL SREW(JTAPE2)                                                 CPP19490
      CALL ZERO(EE,NBASIS*NBASIS)                                       CPP19500
      ELEC=A00                                                          CPP19510
      DO 122 IJCI=1,ICIUNQ                                              CPP19520
      IEX=IPX(IJCI)                                                     CPP19530
      JEX=JPX(IJCI)                                                     CPP19540
      IJEX=IOFF(MAX0(IEX,JEX))+MIN0(IEX,JEX)                            CPP19550
      IBCI=ICITYP(IJCI)                                                 CPP19560
      IPT=IEIJ(IBCI)                                                    CPP19570
      CALL RREAD(JTAPE2,EIJ,NBASQ,IPT)                                  CPP19580
      CIFAC=CI(IEX)*CI(JEX)*TWO                                         CPP19590
      IF(IEX.EQ.JEX) CIFAC=CIFAC*HALF                                   CPP19600
      ELEC=ELEC+HCI(IJEX)*CIFAC                                         CPP19610
      DO 121 I=1,NBASIS                                                 CPP19620
      DO 121 J=1,NBASIS                                                 CPP19630
      EE(I,J)=EE(I,J)+EIJ(I,J)*CIFAC                                    CPP19640
  121 CONTINUE                                                          CPP19650
  122 CONTINUE                                                          CPP19660
      ECITOT=ENUC+ELEC                                                  CPP19670
      WRITE(6,6) ESCF,ENUC,ELEC,ECITOT                                  CPP19680
      IF(IPRNT.LE.2) GO TO 204                                          CPP19690
      WRITE(6,7)                                                        CPP19700
      CALL MATOUT(EE,NBASIS,NBASIS,NBASIS,NBASIS,6)                     CPP19710
C                                                                       CPP19720
  204 CONTINUE                                                          CPP19730
      CALL SREW(JTAPE2)                                                 CPP19740
      RETURN                                                            CPP19750
      END                                                               CPP19760
      SUBROUTINE AMAT(NIJ,KIJ,ALPA,BETA,ZETA,CC,LBLO,BUFO,NNC)          CPP19770
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP19780
      DIMENSION NIJ(NTRI,2),KIJ(NTRI,2)                                 CPP19790
      DIMENSION ALPA(NTRI),BETA(NTRI),ZETA(NTRI,NTYPEP)                 CPP19800
      DIMENSION CC(NNC)                                                 CPP19810
      DIMENSION LBLO(MAXBF4),BUFO(MAXBF2)                               CPP19820
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP19830
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP19840
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPP19850
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP19860
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP19870
      COMMON/TAPES/ITAP35,ITAP42,ITAP43,ITAP44,ITAP47                   CPP19880
      COMMON/WORKS/JTAPE1,JTAPE2                                        CPP19890
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP19900
      COMMON/CI102/NIND,NDEP                                            CPP19910
      DATA A00,HALF / 0.0D+00 , 0.5D+00 /                               CPP19920
    1 FORMAT(//,2X,' A MATRIX'/)                                        CPP19930
    2 FORMAT(2X,6I5,3F15.8)                                             CPP19940
    3 FORMAT(2X,6I5,F15.8)                                              CPP19950
C                                                                       CPP19960
C   FORM A MATRIX                                                       CPP19970
C   SEE EQ.(23)                                                         CPP19980
C                                                                       CPP19990
      CALL RDINT2(CC,NNC,ITAP35,NNC,LBLO,BUFO,IPRNT)                    CPP20000
      CALL SREW(JTAPE1)                                                 CPP20010
      IF(IPRNT.LE.4) GO TO 301                                          CPP20020
      WRITE(6,1)                                                        CPP20030
C                                                                       CPP20040
C   ELEMENTS FOR INDEPENDENT PAIRS                                      CPP20050
  301 IBL=0                                                             CPP20060
      DO 102 II=1,NIND                                                  CPP20070
      I=NIJ(II,1)                                                       CPP20080
      J=NIJ(II,2)                                                       CPP20090
      ITYP=MOTYP(I)                                                     CPP20100
      JTYP=MOTYP(J)                                                     CPP20110
      IJ=IOFF(I)+J                                                      CPP20120
C                                                                       CPP20130
      DO 101 JJ=1,II                                                    CPP20140
      K=NIJ(JJ,1)                                                       CPP20150
      L=NIJ(JJ,2)                                                       CPP20160
      KL=IOFF(K)+L                                                      CPP20170
      VALU1=A00                                                         CPP20180
      VALU2=A00                                                         CPP20190
      ZVAL1=A00                                                         CPP20200
      ZVAL2=A00                                                         CPP20210
      ZVAL3=A00                                                         CPP20220
      ZVAL4=A00                                                         CPP20230
      VALU1=VALU1+GSAI(I,J,K,L,ALPA,BETA,CC,NNC)                        CPP20240
      IF(L.GT.JOCC) GO TO 201                                           CPP20250
      VALU1=VALU1-GSAI(I,J,L,K,ALPA,BETA,CC,NNC)                        CPP20260
  201 IF(J.NE.K) GO TO 202                                              CPP20270
      IL=IOFF(MAX0(I,L))+MIN0(I,L)                                      CPP20280
      ZVAL1=ZETA(IL,ITYP)-ZETA(IL,JTYP)                                 CPP20290
  202 IF(I.NE.K) GO TO 203                                              CPP20300
      JL=IOFF(MAX0(J,L))+MIN0(J,L)                                      CPP20310
      ZVAL2=-ZETA(JL,JTYP)+ZETA(JL,ITYP)                                CPP20320
  203 IF(J.NE.L) GO TO 204                                              CPP20330
      IK=IOFF(MAX0(I,K))+MIN0(I,K)                                      CPP20340
      ZVAL3=-ZETA(IK,ITYP)+ZETA(IK,JTYP)                                CPP20350
  204 IF(I.NE.L) GO TO 205                                              CPP20360
      JK=IOFF(MAX0(J,K))+MIN0(J,K)                                      CPP20370
      ZVAL4=ZETA(JK,JTYP)-ZETA(JK,ITYP)                                 CPP20380
  205 VALU2=ZVAL1+ZVAL2+ZVAL3+ZVAL4                                     CPP20390
      VALUT=VALU1+VALU2                                                 CPP20400
C     IF(IPRNT.LE.4) GO TO 206                                          CPP20410
C     WRITE(6,2) II,JJ,I,J,K,L,VALU1,VALU2,VALUT                        CPP20420
  206 IBL=IBL+1                                                         CPP20430
      LBLO(IBL+IBL-1)=II                                                CPP20440
      LBLO(IBL+IBL)=JJ                                                  CPP20450
      BUFO(MAXBUF+IBL)=VALUT                                            CPP20460
      IF(IBL.LT.MAXBUF) GO TO 101                                       CPP20470
      IBL=0                                                             CPP20480
      CALL SWRIT(JTAPE1,LBLO,MAXBF4)                                    CPP20490
  101 CONTINUE                                                          CPP20500
  102 CONTINUE                                                          CPP20510
C                                                                       CPP20520
C   STORE THE LAST BUFFER                                               CPP20530
      IBL=IBL+1                                                         CPP20540
      LBLO(IBL+IBL-1)=0                                                 CPP20550
      LBLO(IBL+IBL)=0                                                   CPP20560
      BUFO(MAXBUF+IBL)=A00                                              CPP20570
      CALL SWRIT(JTAPE1,LBLO,MAXBF4)                                    CPP20580
C                                                                       CPP20590
C   ELEMENTS FOR DEPENDENT PAIRS                                        CPP20600
CCCCC IF(ICI.EQ.0) GO TO 210                                            CPP20610
      IBL=0                                                             CPP20620
      DO 105 II=1,NDEP                                                  CPP20630
      I=KIJ(II,1)                                                       CPP20640
      J=KIJ(II,2)                                                       CPP20650
      DO 104 JJ=1,NIND                                                  CPP20660
      K=NIJ(JJ,1)                                                       CPP20670
      L=NIJ(JJ,2)                                                       CPP20680
      KTYP=MOTYP(K)                                                     CPP20690
      LTYP=MOTYP(L)                                                     CPP20700
      A=ZIJ(I,J,K,L,CC,NNC)*(FOCC(LTYP)-FOCC(KTYP))*HALF                CPP20710
      IBL=IBL+1                                                         CPP20720
      LBLO(IBL+IBL-1)=256*I+J                                           CPP20730
      LBLO(IBL+IBL)=256*K+L                                             CPP20740
      BUFO(MAXBUF+IBL)=A                                                CPP20750
C     IF(IPRNT.LE.4) GO TO 207                                          CPP20760
C     WRITE(6,3) II,JJ,I,J,K,L,A                                        CPP20770
  207 CONTINUE                                                          CPP20780
      IF(IBL.LT.MAXBUF) GO TO 104                                       CPP20790
      IBL=0                                                             CPP20800
      CALL SWRIT(JTAPE1,LBLO,MAXBF4)                                    CPP20810
  104 CONTINUE                                                          CPP20820
  105 CONTINUE                                                          CPP20830
C                                                                       CPP20840
C   STORE THE LAST BUFFER                                               CPP20850
      IBL=IBL+1                                                         CPP20860
      LBLO(IBL+IBL-1)=0                                                 CPP20870
      LBLO(IBL+IBL)=0                                                   CPP20880
      BUFO(MAXBUF+IBL)=A00                                              CPP20890
      CALL SWRIT(JTAPE1,LBLO,MAXBF4)                                    CPP20900
C                                                                       CPP20910
  210 CALL SREW(JTAPE1)                                                 CPP20920
      RETURN                                                            CPP20930
      END                                                               CPP20940
      FUNCTION ZIJ(I,J,K,L,CC,NNC)                                      CPP20950
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP20960
      DIMENSION CC(NNC)                                                 CPP20970
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP20980
      DATA FOUR / 4.0D+00 /                                             CPP20990
C                                                                       CPP21000
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPP21010
      KL=IOFF(MAX0(K,L))+MIN0(K,L)                                      CPP21020
      IL=IOFF(MAX0(I,L))+MIN0(I,L)                                      CPP21030
      JK=IOFF(MAX0(J,K))+MIN0(J,K)                                      CPP21040
      IK=IOFF(MAX0(I,K))+MIN0(I,K)                                      CPP21050
      JL=IOFF(MAX0(J,L))+MIN0(J,L)                                      CPP21060
      IJKL=IOFF(MAX0(IJ,KL))+MIN0(IJ,KL)                                CPP21070
      ILJK=IOFF(MAX0(IL,JK))+MIN0(IL,JK)                                CPP21080
      IKJL=IOFF(MAX0(IK,JL))+MIN0(IK,JL)                                CPP21090
      ZIJ=CC(IJKL)*FOUR-CC(ILJK)-CC(IKJL)                               CPP21100
      RETURN                                                            CPP21110
      END                                                               CPP21120
      SUBROUTINE BFMAT(NIJ,KIJ,FIJ,E1T,HCIA,HF,HM,EPF,BF,NABC)          CPP21130
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP21140
      DIMENSION NIJ(NTRI,2),KIJ(NTRI,2)                                 CPP21150
      DIMENSION FIJ(NTYPEP,ICIUNQ),E1T(NABC),HCIA(NCTRI,NABC)           CPP21160
      DIMENSION HF(NTRI,3),HM(NTRI,3,NTYPES)                            CPP21170
      DIMENSION EPF(NBASIS,NBASIS),BF(NTRI)                             CPP21180
      DIMENSION DIP(3)                                                  CPP21190
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP21200
      COMMON/CIOFF/IIUNQ,IJUNQ,ICIUNQ,IPX(5050),JPX(5050)               CPP21210
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP21220
      COMMON/DIPOL/DIPDA(135),DIPDM(135),DIPDN(135),DIPDT(135)          CPP21230
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP21240
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP21250
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP21260
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP21270
      COMMON/TAPES/ITAP35,ITAP42,ITAP43,ITAP44,ITAP47                   CPP21280
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP21290
      COMMON/CI102/NIND,NDEP                                            CPP21300
      COMMON/CI103/CI(100)                                              CPP21310
      DATA A00,HALF,TWO / 0.0D+00 , 0.5D+00 , 2.0D+00 /                 CPP21320
      DATA DEBYE / 2.541765480D+00 /                                    CPP21330
    1 FORMAT(//,2X,' DIPDA MATRIX'/)                                    CPP21340
    2 FORMAT(//,2X,' DIPDN MATRIX'/)                                    CPP21350
    3 FORMAT(//,2X,' HF MATRIX, IXYZ = ',I5/)                           CPP21360
    4 FORMAT(//,2X,' EPF MATRIX, IXYZ = ',I5/)                          CPP21370
    5 FORMAT(//,2X,' HCIA MATRIX, IJCI = ',I5/)                         CPP21380
    6 FORMAT(//,2X,' DIPOLE MOMENT (ELECTRONIC PART IN DEBYE)'/)        CPP21390
C                                                                       CPP21400
      CALL RFILE(ITAP43)                                                CPP21410
      CALL SREAD(ITAP43,DIPDA,270)                                      CPP21420
      CALL SREAD(ITAP43,DIPDN,270)                                      CPP21430
      IF(IPRNT.LE.4) GO TO 201                                          CPP21440
      WRITE(6,1)                                                        CPP21450
      CALL MATOUT(DIPDA,3,45,3,N3N,6)                                   CPP21460
      WRITE(6,2)                                                        CPP21470
      CALL MATOUT(DIPDN,3,45,3,N3N,6)                                   CPP21480
C                                                                       CPP21490
  201 CONTINUE                                                          CPP21500
      CALL SREW(ITAP44)                                                 CPP21510
      DO 103 IXYZ=1,3                                                   CPP21520
      IH=IHA(IXYZ+N3N)                                                  CPP21530
      CALL SREAD(ITAP43,HF(1,IXYZ),NTRI2)                               CPP21540
      CALL RWRIT(ITAP44,HF(1,IXYZ),NTRI2,IH)                            CPP21550
      DO 102 ITYP=1,NTYPES                                              CPP21560
      FAC=FOCC(ITYP)*HALF                                               CPP21570
      DO 101 I=1,NTRI                                                   CPP21580
  101 HM(I,IXYZ,ITYP)=HF(I,IXYZ)*FAC                                    CPP21590
  102 CONTINUE                                                          CPP21600
  103 CONTINUE                                                          CPP21610
C                                                                       CPP21620
C   STORE DERIVATIVE LAGRANGIAN MATRICES                                CPP21630
      DO 110 IXYZ=1,3                                                   CPP21640
      IE=IEA(IXYZ+N3N)                                                  CPP21650
      IB=IBA(IXYZ+N3N)                                                  CPP21660
      CALL ZERO(EPF,NBASIS*NBASIS)                                      CPP21670
C                                                                       CPP21680
      DO 107 ITYP=1,NTYPES                                              CPP21690
      DO 106 I=NSTR(ITYP),NEND(ITYP)                                    CPP21700
      DO 106 J=1,NBASIS                                                 CPP21710
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPP21720
      EPF(I,J)=HM(IJ,IXYZ,ITYP)                                         CPP21730
  106 CONTINUE                                                          CPP21740
  107 CONTINUE                                                          CPP21750
C                                                                       CPP21760
C   ELEMENTS FOR INDEPENDENT PAIRS                                      CPP21770
      DO 108 II=1,NIND                                                  CPP21780
      I=NIJ(II,1)                                                       CPP21790
      J=NIJ(II,2)                                                       CPP21800
      IJ=IOFF(I)+J                                                      CPP21810
      BF(IJ)=EPF(I,J)-EPF(J,I)                                          CPP21820
  108 CONTINUE                                                          CPP21830
C                                                                       CPP21840
C   ELEMENTS FOR DEPENDENT PAIRS                                        CPP21850
      DO 109 II=1,NDEP                                                  CPP21860
      I=KIJ(II,1)                                                       CPP21870
      J=KIJ(II,2)                                                       CPP21880
      IJ=IOFF(I)+J                                                      CPP21890
      BF(IJ)=HF(IJ,IXYZ)                                                CPP21900
  109 CONTINUE                                                          CPP21910
C                                                                       CPP21920
      CALL RWRIT(ITAP44,BF,NTRI2,IB)                                    CPP21930
      CALL RWRIT(ITAP44,EPF,NBASQ,IE)                                   CPP21940
      IF(IPRNT.LE.4) GO TO 110                                          CPP21950
      WRITE(6,4) IXYZ                                                   CPP21960
      CALL MATOUT(EPF,NBASIS,NBASIS,NBASIS,NBASIS,6)                    CPP21970
  110 CONTINUE                                                          CPP21980
C                                                                       CPP21990
C   FORM THE HIJF MATRIX                                                CPP22000
C   NOTE THE SUMMATION LIMIT                                            CPP22010
      DO 115 IJCI=1,IIUNQ                                               CPP22020
      IEX=IPX(IJCI)                                                     CPP22030
      JEX=JPX(IJCI)                                                     CPP22040
      IJEX=IOFF(MAX0(IEX,JEX))+MIN0(IEX,JEX)                            CPP22050
      DO 113 IXYZ=1,3                                                   CPP22060
C                                                                       CPP22070
      ELECF=A00                                                         CPP22080
      DO 112 I=1,JOCC                                                   CPP22090
      II=IOFF(I+1)                                                      CPP22100
      ITYP=MOTYP(I)                                                     CPP22110
      ELECF=ELECF+FIJ(ITYP,IJCI)*HF(II,IXYZ)*TWO                        CPP22120
  112 CONTINUE                                                          CPP22130
      HCIA(IJEX,IXYZ+N3N)=ELECF                                         CPP22140
  113 CONTINUE                                                          CPP22150
      IF(IPRNT.LE.2) GO TO 115                                          CPP22160
      WRITE(6,5) IJEX                                                   CPP22170
      CALL PRINT(HCIA(IJEX,IXYZ+N3N),NCTRI,NCPEX,6)                     CPP22180
  115 CONTINUE                                                          CPP22190
C                                                                       CPP22200
C   CALCULATE DIPOLE MOMENTS FOR A TEST                                 CPP22210
      DO 120 IXYZ=1,3                                                   CPP22220
      VALU=A00                                                          CPP22230
      DO 118 IJCI=1,ICIUNQ                                              CPP22240
      IEX=IPX(IJCI)                                                     CPP22250
      JEX=JPX(IJCI)                                                     CPP22260
      IJEX=IOFF(MAX0(IEX,JEX))+MIN0(IEX,JEX)                            CPP22270
      FAC=CI(IEX)*CI(JEX)*TWO                                           CPP22280
      IF(IEX.EQ.JEX) FAC=FAC*HALF                                       CPP22290
      VALU=VALU+HCIA(IJEX,IXYZ+N3N)*FAC                                 CPP22300
  118 CONTINUE                                                          CPP22310
      DIP(IXYZ)=-VALU*DEBYE                                             CPP22320
      E1T(IXYZ+N3N)=VALU                                                CPP22330
  120 CONTINUE                                                          CPP22340
      WRITE(6,6)                                                        CPP22350
      CALL MATOUT(DIP,3,1,3,1,6)                                        CPP22360
C                                                                       CPP22370
      CALL SREW(ITAP43)                                                 CPP22380
      CALL SREW(ITAP44)                                                 CPP22390
C                                                                       CPP22400
      RETURN                                                            CPP22410
      END                                                               CPP22420
      SUBROUTINE AB2CI(NIJ,E1T,HCIA,HCI,BAF2,A22,A12,SA,EIJ,NABC)       CPP22430
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP22440
      DIMENSION NIJ(NTRI,2),SA(NTRI,N3N)                                CPP22450
      DIMENSION E1T(NABC),HCIA(NCTRI,NABC),HCI(NCTRI)                   CPP22460
      DIMENSION BAF2(NCPEX,NABC),A22(NCPEX,NCPEX),A12(NIND,NCPEX)       CPP22470
      DIMENSION EIJ(NBASIS,NBASIS)                                      CPP22480
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP22490
      COMMON/CIOFF/IIUNQ,IJUNQ,ICIUNQ,IPX(5050),JPX(5050)               CPP22500
      COMMON/CIPOS/ICITYP(5050),ICPOS(200,100)                          CPP22510
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP22520
      COMMON/ENRGY/ESCF,ENUC,ELEC                                       CPP22530
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP22540
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP22550
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP22560
      COMMON/TAPES/ITAP35,ITAP42,ITAP43,ITAP44,ITAP47                   CPP22570
      COMMON/WORKS/JTAPE1,JTAPE2                                        CPP22580
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP22590
      COMMON/CI102/NIND,NDEP                                            CPP22600
      COMMON/CI103/CI(100)                                              CPP22610
      COMMON/CI104/IEIJ(5050)                                           CPP22620
      DATA A00,HALF,TWO,FOUR /  0.0D+00 , 0.5D+00 , 2.0D+00 , 4.0D+00 / CPP22630
    1 FORMAT(//,2X,' BAF2 MATRIX'/)                                     CPP22640
    2 FORMAT(//,2X,' A12 MATRIX'/)                                      CPP22650
    3 FORMAT(//,2X,' A22 MATRIX'/)                                      CPP22660
C                                                                       CPP22670
C   CALL BACK SA MATRICES                                               CPP22680
      CALL SREW(ITAP44)                                                 CPP22690
      DO 101 IABC=1,N3N                                                 CPP22700
      IS=ISA(IABC)                                                      CPP22710
      CALL RREAD(ITAP44,SA(1,IABC),NTRI2,IS)                            CPP22720
  101 CONTINUE                                                          CPP22730
C                                                                       CPP22740
C   FORM THE BAF AND A12 MATRICES                                       CPP22750
      CALL ZERO(BAF2,NCPEX*NABC)                                        CPP22760
      CALL ZERO(A12,NIND*NCPEX)                                         CPP22770
      CALL SREW(JTAPE2)                                                 CPP22780
C                                                                       CPP22790
      DO 110 ICI=1,NCPEX                                                CPP22800
      DO 102 IABC=1,NABC                                                CPP22810
      BAF2(ICI,IABC)=E1T(IABC)*CI(ICI)*HALF                             CPP22820
  102 CONTINUE                                                          CPP22830
      DO 109 JCI=1,NCPEX                                                CPP22840
      IJCI=IOFF(MAX0(ICI,JCI))+MIN0(ICI,JCI)                            CPP22850
      DO 103 III=1,ICIUNQ                                               CPP22860
      IEX=IPX(III)                                                      CPP22870
      JEX=JPX(III)                                                      CPP22880
      IJEX=IOFF(MAX0(IEX,JEX))+MIN0(IEX,JEX)                            CPP22890
      IF(IJEX.EQ.IJCI) GO TO 301                                        CPP22900
  103 CONTINUE                                                          CPP22910
      GO TO 109                                                         CPP22920
C                                                                       CPP22930
  301 CONTINUE                                                          CPP22940
      IBCI=ICITYP(III)                                                  CPP22950
      IPT=IEIJ(IBCI)                                                    CPP22960
      CALL RREAD(JTAPE2,EIJ,NBASQ,IPT)                                  CPP22970
C                                                                       CPP22980
      DO 105 IABC=1,NABC                                                CPP22990
      VALU0=A00                                                         CPP23000
      IF(IABC.GT.N3N) GO TO 201                                         CPP23010
      DO 104 I=1,JOCC                                                   CPP23020
      DO 104 J=1,I                                                      CPP23030
      IJ=IOFF(I)+J                                                      CPP23040
      FAC4=FOUR                                                         CPP23050
      IF(I.EQ.J) FAC4=TWO                                               CPP23060
      VALU0=VALU0+SA(IJ,IABC)*EIJ(I,J)*FAC4                             CPP23070
  104 CONTINUE                                                          CPP23080
  201 CONTINUE                                                          CPP23090
      VALU2=(-HCIA(IJCI,IABC)+VALU0)*CI(JCI)                            CPP23100
      BAF2(ICI,IABC)=BAF2(ICI,IABC)+VALU2*HALF                          CPP23110
  105 CONTINUE                                                          CPP23120
C                                                                       CPP23130
      DO 107 II=1,NIND                                                  CPP23140
      I=NIJ(II,1)                                                       CPP23150
      J=NIJ(II,2)                                                       CPP23160
      VALU=(EIJ(I,J)-EIJ(J,I))*CI(JCI)                                  CPP23170
      A12(II,ICI)=A12(II,ICI)-VALU*TWO                                  CPP23180
  107 CONTINUE                                                          CPP23190
C                                                                       CPP23200
  109 CONTINUE                                                          CPP23210
  110 CONTINUE                                                          CPP23220
      IF(IPRNT.LE.2) GO TO 202                                          CPP23230
      WRITE(6,1)                                                        CPP23240
      CALL MATOUT(BAF2,NCPEX,NABC,NCPEX,NABC,6)                         CPP23250
      WRITE(6,2)                                                        CPP23260
      CALL MATOUT(A12,NIND,NCPEX,NIND,NCPEX,6)                          CPP23270
C                                                                       CPP23280
C   FORM THE A22 MATRIX                                                 CPP23290
  202 CONTINUE                                                          CPP23300
      CALL ZERO(A22,NCPEX*NCPEX)                                        CPP23310
      IJ=0                                                              CPP23320
      DO 115 I=1,NCPEX                                                  CPP23330
      DO 115 J=1,I                                                      CPP23340
      IJ=IJ+1                                                           CPP23350
      VALU=HCI(IJ)+CI(I)*CI(J)                                          CPP23360
      IF(I.EQ.J) VALU=VALU-ELEC                                         CPP23370
      A22(I,J)=VALU*HALF                                                CPP23380
      A22(J,I)=A22(I,J)                                                 CPP23390
  115 CONTINUE                                                          CPP23400
      IF(IPRNT.LE.2) GO TO 203                                          CPP23410
      WRITE(6,3)                                                        CPP23420
      CALL MATOUT(A22,NCPEX,NCPEX,NCPEX,NCPEX,6)                        CPP23430
C                                                                       CPP23440
  203 CONTINUE                                                          CPP23450
      CALL SREW(JTAPE2)                                                 CPP23460
      CALL SREW(ITAP44)                                                 CPP23470
      RETURN                                                            CPP23480
      END                                                               CPP23490
      SUBROUTINE CPHF(EIG,NIJ,KIJ,BAF2,CIA,A22,A12,B0,UA,SA,BA,C,       CPP23500
     1                LBLI,BUFI,NINDPX,NABC)                            CPP23510
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP23520
      REAL*16 DLIM,DETERM                                               CPP23530
      DIMENSION EIG(NBASIS),B0(NTRI),UA(NBASIS,NBASIS)                  CPP23540
      DIMENSION BAF2(NCPEX,NABC),CIA(NCPEX,NABC)                        CPP23550
      DIMENSION A22(NCPEX,NCPEX),A12(NIND,NCPEX)                        CPP23560
      DIMENSION NIJ(NTRI,2),KIJ(NTRI,2)                                 CPP23570
      DIMENSION SA(NTRI,N3N),BA(NTRI,NABC),C(NINDPX,NINDPX+NABC)        CPP23580
      DIMENSION LBLI(MAXBF4),BUFI(MAXBF2)                               CPP23590
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP23600
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP23610
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP23620
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPP23630
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP23640
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP23650
      COMMON/TAPES/ITAP35,ITAP42,ITAP43,ITAP44,ITAP47                   CPP23660
      COMMON/WORKS/JTAPE1,JTAPE2                                        CPP23670
      COMMON/CI102/NIND,NDEP                                            CPP23680
      DATA CRIT / 1.0D-6 /                                              CPP23690
      DATA DLIM / 1.0Q-76 /                                             CPP23700
      DATA A00,HALF,ONE,TWO,FOUR / 0.0D+00 , 0.5D+00 , 1.0D+00 ,        CPP23710
     1                             2.0D+00 , 4.0D+00 /                  CPP23720
    1 FORMAT(//,2X,' C MATRIX'/)                                        CPP23730
    2 FORMAT(//,2X,' DETERM = ',Q20.10/)                                CPP23740
    3 FORMAT(//,2X,'***THE VALUE OF DETERMINANT IS TOO SMALL***',       CPP23750
     1 5X,' DETERM = ',Q20.10/)                                         CPP23760
    4 FORMAT(2X,2I7,5X,F15.8)                                           CPP23770
    5 FORMAT(//,2X,' B MATRIX, IABC = ',I5/)                            CPP23780
    6 FORMAT(//,2X,' U MATRIX, IABC = ',I5/)                            CPP23790
C                                                                       CPP23800
C***********************************                                    CPP23810
C***SET UP SIMULTANEOUS EQUATIONS***                                    CPP23820
C***********************************                                    CPP23830
C                                                                       CPP23840
      CALL SREW(ITAP44)                                                 CPP23850
      DO 101 IABC=1,NABC                                                CPP23860
      IB=IBA(IABC)                                                      CPP23870
      CALL RREAD(ITAP44,BA(1,IABC),NTRI2,IB)                            CPP23880
      IF(IABC.GT.N3N) GO TO 101                                         CPP23890
      IS=ISA(IABC)                                                      CPP23900
      CALL RREAD(ITAP44,SA(1,IABC),NTRI2,IS)                            CPP23910
  101 CONTINUE                                                          CPP23920
C                                                                       CPP23930
C"""""""""""""""""""                                                    CPP23940
C"""FORM A MATRIX"""                                                    CPP23950
C"""""""""""""""""""                                                    CPP23960
C                                                                       CPP23970
C:::ORBITAL SPACE:::                                                    CPP23980
      CALL SREW(JTAPE1)                                                 CPP23990
  500 CALL SREAD(JTAPE1,LBLI,MAXBF4)                                    CPP24000
      DO 102 M=1,MAXBUF                                                 CPP24010
      II=LBLI(M+M-1)                                                    CPP24020
      IF(II.EQ.0) GO TO 600                                             CPP24030
      JJ=LBLI(M+M)                                                      CPP24040
      C(II,JJ)=BUFI(MAXBUF+M)                                           CPP24050
      C(JJ,II)=BUFI(MAXBUF+M)                                           CPP24060
  102 CONTINUE                                                          CPP24070
      GO TO 500                                                         CPP24080
  600 CONTINUE                                                          CPP24090
C                                                                       CPP24100
C:::CI SPACE:::                                                         CPP24110
C   SEE EQS.(24) AND (25)                                               CPP24120
      DO 104 II=1,NIND                                                  CPP24130
      DO 103 ICI=1,NCPEX                                                CPP24140
      C(II,NIND+ICI)=A12(II,ICI)                                        CPP24150
      C(NIND+ICI,II)=A12(II,ICI)                                        CPP24160
  103 CONTINUE                                                          CPP24170
  104 CONTINUE                                                          CPP24180
      DO 105 ICI=1,NCPEX                                                CPP24190
      DO 105 JCI=1,NCPEX                                                CPP24200
      C(NIND+ICI,NIND+JCI)=A22(ICI,JCI)                                 CPP24210
      C(NIND+JCI,NIND+ICI)=A22(ICI,JCI)                                 CPP24220
  105 CONTINUE                                                          CPP24230
C                                                                       CPP24240
C"""""""""""""""""""                                                    CPP24250
C"""FORM B MATRIX"""                                                    CPP24260
C"""""""""""""""""""                                                    CPP24270
C                                                                       CPP24280
C:::ORBITAL SPACE:::                                                    CPP24290
      DO 107 IABC=1,NABC                                                CPP24300
      DO 106 II=1,NIND                                                  CPP24310
      I=NIJ(II,1)                                                       CPP24320
      J=NIJ(II,2)                                                       CPP24330
      IJ=IOFF(I)+J                                                      CPP24340
      C(II,NINDPX+IABC)=BA(IJ,IABC)                                     CPP24350
  106 CONTINUE                                                          CPP24360
  107 CONTINUE                                                          CPP24370
C                                                                       CPP24380
C:::CI SPACE:::                                                         CPP24390
C   SEE EQ.(27)                                                         CPP24400
      DO 109 IABC=1,NABC                                                CPP24410
      DO 108 ICI=1,NCPEX                                                CPP24420
      C(NIND+ICI,NINDPX+IABC)=BAF2(ICI,IABC)                            CPP24430
  108 CONTINUE                                                          CPP24440
  109 CONTINUE                                                          CPP24450
C                                                                       CPP24460
      IF(IPRNT.LE.3) GO TO 202                                          CPP24470
      WRITE(6,1)                                                        CPP24480
      CALL MATOUT(C,NINDPX,NINDPX+NABC,NINDPX,NINDPX+NABC,6)            CPP24490
C                                                                       CPP24500
C""""""""""""""""""""""""""""""""""""""                                 CPP24510
C"""SOLVE THE SIMULTANEOUS EQUATIONS"""                                 CPP24520
C""""""""""""""""""""""""""""""""""""""                                 CPP24530
  202 CALL FLINQ(C,NINDPX,NINDPX,NABC,DETERM)                           CPP24540
      WRITE(6,2) DETERM                                                 CPP24550
      IF(IPRNT.GT.3)                                                    CPP24560
     * CALL MATOUT(C,NINDPX,NINDPX+NABC,NINDPX,NINDPX+NABC,6)           CPP24570
      IF(QABS(DETERM).LE.DLIM) GO TO 203                                CPP24580
      GO TO 205                                                         CPP24590
  203 WRITE(6,3) DETERM                                                 CPP24600
      STOP                                                              CPP24610
C                                                                       CPP24620
C"""""""""""""""""""""""""""""""""""""""""""""""                        CPP24630
C"""OBTAIN INDEPENDENT ELEMENTS OF U MATRICES"""                        CPP24640
C"""""""""""""""""""""""""""""""""""""""""""""""                        CPP24650
C                                                                       CPP24660
C:::ORBITAL SPACE:::                                                    CPP24670
  205 CONTINUE                                                          CPP24680
      DO 113 IABC=1,NABC                                                CPP24690
      DO 111 II=1,NIND                                                  CPP24700
      I=NIJ(II,1)                                                       CPP24710
      J=NIJ(II,2)                                                       CPP24720
      IJ=IOFF(I)+J                                                      CPP24730
      BA(IJ,IABC)=C(II,NINDPX+IABC)                                     CPP24740
  111 CONTINUE                                                          CPP24750
C                                                                       CPP24760
C:::CI SPACE:::                                                         CPP24770
      DO 112 ICI=1,NCPEX                                                CPP24780
      CIA(ICI,IABC)=C(NIND+ICI,NINDPX+IABC)                             CPP24790
      VALUC=CIA(ICI,IABC)                                               CPP24800
      IF(IPRNT.LE.2) GO TO 112                                          CPP24810
      WRITE(6,4) IABC,ICI,VALUC                                         CPP24820
  112 CONTINUE                                                          CPP24830
C                                                                       CPP24840
      IF(IPRNT.LE.3) GO TO 113                                          CPP24850
      WRITE(6,5)                                                        CPP24860
      CALL PRINT(BA(1,IABC),NTRI,NBASIS,6)                              CPP24870
  113 CONTINUE                                                          CPP24880
C                                                                       CPP24890
C   CALCULATE DEPENDENT ELEMENTS OF U MATRICES                          CPP24900
CCCCC IF(ICI.EQ.0) GO TO 210                                            CPP24910
  700 CONTINUE                                                          CPP24920
      CALL SREAD(JTAPE1,LBLI,MAXBF4)                                    CPP24930
      DO 115 M=1,MAXBUF                                                 CPP24940
      IB1=LBLI(M+M-1)                                                   CPP24950
      IF(IB1.EQ.0) GO TO 800                                            CPP24960
      IB2=LBLI(M+M)                                                     CPP24970
      J=MOD(IB1,256)                                                    CPP24980
      I=IB1/256                                                         CPP24990
      L=MOD(IB2,256)                                                    CPP25000
      K=IB2/256                                                         CPP25010
      IJ=IOFF(I)+J                                                      CPP25020
      KL=IOFF(K)+L                                                      CPP25030
      VALU=BUFI(MAXBUF+M)                                               CPP25040
      DO 114 IABC=1,NABC                                                CPP25050
      BA(IJ,IABC)=BA(IJ,IABC)+VALU*BA(KL,IABC)                          CPP25060
  114 CONTINUE                                                          CPP25070
  115 CONTINUE                                                          CPP25080
      GO TO 700                                                         CPP25090
  800 CALL SREW(JTAPE1)                                                 CPP25100
C                                                                       CPP25110
      DO 125 IABC=1,NABC                                                CPP25120
      IB=IBA(IABC)                                                      CPP25130
      IU=IUA(IABC)                                                      CPP25140
C                                                                       CPP25150
      DO 116 II=1,NDEP                                                  CPP25160
      I=KIJ(II,1)                                                       CPP25170
      J=KIJ(II,2)                                                       CPP25180
      IJ=IOFF(I)+J                                                      CPP25190
      DELJI=EIG(J)-EIG(I)                                               CPP25200
      IF(I.EQ.J) DELJI=ONE                                              CPP25210
      IF(DABS(DELJI).LT.CRIT) GO TO 206                                 CPP25220
      BA(IJ,IABC)=BA(IJ,IABC)/DELJI                                     CPP25230
      GO TO 116                                                         CPP25240
  206 BA(IJ,IABC)=A00                                                   CPP25250
  116 CONTINUE                                                          CPP25260
C                                                                       CPP25270
      IF(IABC.GT.N3N) GO TO 208                                         CPP25280
      DO 118 I=1,NBASIS                                                 CPP25290
      DO 118 J=1,I                                                      CPP25300
      IJ=IOFF(I)+J                                                      CPP25310
      B0(IJ)=BA(IJ,IABC)                                                CPP25320
      IF(I.EQ.J) GO TO 207                                              CPP25330
      UA(I,J)=BA(IJ,IABC)                                               CPP25340
      UA(J,I)=-UA(I,J)-SA(IJ,IABC)                                      CPP25350
      GO TO 118                                                         CPP25360
  207 UA(I,I)=-SA(IJ,IABC)*HALF                                         CPP25370
  118 CONTINUE                                                          CPP25380
      GO TO 210                                                         CPP25390
C                                                                       CPP25400
  208 CONTINUE                                                          CPP25410
      DO 119 I=1,NBASIS                                                 CPP25420
      DO 119 J=1,I                                                      CPP25430
      IJ=IOFF(I)+J                                                      CPP25440
      B0(IJ)=BA(IJ,IABC)                                                CPP25450
      IF(I.EQ.J) GO TO 209                                              CPP25460
      UA(I,J)=BA(IJ,IABC)                                               CPP25470
      UA(J,I)=-UA(I,J)                                                  CPP25480
      GO TO 119                                                         CPP25490
  209 UA(I,I)=A00                                                       CPP25500
  119 CONTINUE                                                          CPP25510
C                                                                       CPP25520
  210 CONTINUE                                                          CPP25530
      CALL RWRIT(ITAP44,B0,NTRI2,IB)                                    CPP25540
      CALL RWRIT(ITAP44,UA,NBASQ,IU)                                    CPP25550
      IF(IPRNT.LE.3) GO TO 125                                          CPP25560
      WRITE(6,6) IABC                                                   CPP25570
      CALL MATOUT(UA,NBASIS,NBASIS,NBASIS,NBASIS,6)                     CPP25580
C                                                                       CPP25590
  125 CONTINUE                                                          CPP25600
C                                                                       CPP25610
      CALL SREW(ITAP44)                                                 CPP25620
      CALL SREW(JTAPE1)                                                 CPP25630
      RETURN                                                            CPP25640
      END                                                               CPP25650
      SUBROUTINE WAMAT(ALPA,BETA,EPS,ZETA,E2M,E2P,E2Q,EPA,UA,SA,WA,     CPP25660
     1                 CC,LBLI,BUFI,NNC)                                CPP25670
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP25680
      DIMENSION ALPA(NTRI),BETA(NTRI),EPS(NTRI),ZETA(NTRI,NTYPEP)       CPP25690
      DIMENSION E2M(N3N,N3N),E2P(N3N,N3N),E2Q(N3N,N3N)                  CPP25700
      DIMENSION EPA(NBASIS,NBASIS),UA(NBASIS,NBASIS)                    CPP25710
      DIMENSION SA(NTRI,N3N),WA(NBASIS,NBASIS,N3N),CC(NNC)              CPP25720
      DIMENSION LBLI(MAXBF4),BUFI(MAXBF2)                               CPP25730
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP25740
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP25750
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPP25760
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP25770
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP25780
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP25790
      COMMON/TAPES/ITAP35,ITAP42,ITAP43,ITAP44,ITAP47                   CPP25800
      COMMON/WORKS/JTAPE1,JTAPE2                                        CPP25810
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP25820
      DATA A00,TWO,FOUR / 0.0D+00 , 2.0D+00 , 4.0D+00 /                 CPP25830
    1 FORMAT(//,2X,' WA MATRIX, IABC = ',I5/)                           CPP25840
    2 FORMAT(//,2X,' E2P MATRIX IN WAMAT'/)                             CPP25850
    3 FORMAT(//,2X,' E2Q MATRIX IN WAMAT'/)                             CPP25860
    4 FORMAT(//,2X,' E2M MATRIX IN WAMAT'/)                             CPP25870
C                                                                       CPP25880
C   FORM THE WA MATRICES                                                CPP25890
C   SEE EQ.(13)                                                         CPP25900
C                                                                       CPP25910
C   READ IN TWO ELECTRON MO INTEGRALS                                   CPP25920
      CALL RDINT2(CC,NNC,ITAP35,NNC,LBLI,BUFI,IPRNT)                    CPP25930
C                                                                       CPP25940
      CALL SREW(ITAP44)                                                 CPP25950
      DO 102 IABC=1,N3N                                                 CPP25960
      IS=ISA(IABC)                                                      CPP25970
      IE=IEA(IABC)                                                      CPP25980
      CALL RREAD(ITAP44,SA(1,IABC),NTRI2,IS)                            CPP25990
      CALL RREAD(ITAP44,EPA,NBASQ,IE)                                   CPP26000
      DO 101 I=1,NBASIS                                                 CPP26010
      DO 101 J=1,NBASIS                                                 CPP26020
      WA(I,J,IABC)=EPA(J,I)+EPA(J,I)                                    CPP26030
  101 CONTINUE                                                          CPP26040
  102 CONTINUE                                                          CPP26050
C                                                                       CPP26060
      DO 105 I=1,NBASIS                                                 CPP26070
      DO 105 J=1,NBASIS                                                 CPP26080
      NJ=MOTYP(J)                                                       CPP26090
      DO 104 K=1,JOCC                                                   CPP26100
      IK=IOFF(MAX0(I,K))+MIN0(I,K)                                      CPP26110
      JK=IOFF(MAX0(J,K))+MIN0(J,K)                                      CPP26120
      FAC1=EPS(IK)+ZETA(IK,NJ)                                          CPP26130
      FAC2=EPS(JK)+EPS(JK)                                              CPP26140
      DO 103 IABC=1,N3N                                                 CPP26150
      WA(I,J,IABC)=WA(I,J,IABC)-SA(JK,IABC)*FAC1-SA(IK,IABC)*FAC2       CPP26160
  103 CONTINUE                                                          CPP26170
  104 CONTINUE                                                          CPP26180
  105 CONTINUE                                                          CPP26190
      DO 108 I=1,NBASIS                                                 CPP26200
      DO 108 J=1,NBASIS                                                 CPP26210
      DO 107 K=1,JOCC                                                   CPP26220
      DO 107 L=1,JOCC                                                   CPP26230
      KL=IOFF(MAX0(K,L))+MIN0(K,L)                                      CPP26240
      FAC3=ABIJ(I,J,K,L,ALPA,BETA,CC,NNC)                               CPP26250
      DO 106 IABC=1,N3N                                                 CPP26260
      WA(I,J,IABC)=WA(I,J,IABC)-SA(KL,IABC)*FAC3                        CPP26270
  106 CONTINUE                                                          CPP26280
  107 CONTINUE                                                          CPP26290
  108 CONTINUE                                                          CPP26300
      IF(IPRNT.LE.3) GO TO 201                                          CPP26310
      DO 110 IABC=1,N3N                                                 CPP26320
      WRITE(6,1) IABC                                                   CPP26330
      CALL MATOUT(WA(1,1,IABC),NBASIS,NBASIS,NBASIS,NBASIS,6)           CPP26340
  110 CONTINUE                                                          CPP26350
C                                                                       CPP26360
C   CALCULATE SCF SECOND DERIVATIVES                                    CPP26370
C   SEE EQ.(12)                                                         CPP26380
C                                                                       CPP26390
C   (B COORDINATE)                                                      CPP26400
  201 CONTINUE                                                          CPP26410
      DO 115 IABC=1,N3N                                                 CPP26420
      IE=IEA(IABC)                                                      CPP26430
      IU=IUA(IABC)                                                      CPP26440
      CALL RREAD(ITAP44,EPA,NBASQ,IE)                                   CPP26450
      CALL RREAD(ITAP44,UA,NBASQ,IU)                                    CPP26460
C                                                                       CPP26470
C   (A COORDINATE)                                                      CPP26480
      DO 114 JABC=1,N3N                                                 CPP26490
      VALUP=A00                                                         CPP26500
      VALUQ=A00                                                         CPP26510
      DO 113 I=1,JOCC                                                   CPP26520
      DO 111 J=1,JOCC                                                   CPP26530
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPP26540
      VALUP=VALUP-SA(IJ,JABC)*EPA(I,J)                                  CPP26550
  111 CONTINUE                                                          CPP26560
      DO 112 J=1,NBASIS                                                 CPP26570
      VALUQ=VALUQ+WA(J,I,JABC)*UA(J,I)                                  CPP26580
  112 CONTINUE                                                          CPP26590
  113 CONTINUE                                                          CPP26600
      E2P(IABC,JABC)=VALUP*TWO                                          CPP26610
      E2Q(IABC,JABC)=VALUQ*TWO                                          CPP26620
  114 CONTINUE                                                          CPP26630
  115 CONTINUE                                                          CPP26640
      IF(IPRNT.LE.3) GO TO 202                                          CPP26650
      WRITE(6,2)                                                        CPP26660
      CALL MATOUT(E2P,N3N,N3N,N3N,N3N,6)                                CPP26670
      WRITE(6,3)                                                        CPP26680
      CALL MATOUT(E2Q,N3N,N3N,N3N,N3N,6)                                CPP26690
C                                                                       CPP26700
C*******************************************                            CPP26710
C***SUM UP CONTRIBUTIONS DUE TO MO CHANGE***                            CPP26720
C*******************************************                            CPP26730
  202 CONTINUE                                                          CPP26740
      DO 120 I=1,N3N                                                    CPP26750
      DO 120 J=1,N3N                                                    CPP26760
      E2M(I,J)=E2P(I,J)+E2Q(I,J)                                        CPP26770
  120 CONTINUE                                                          CPP26780
C                                                                       CPP26790
      IF(IPRNT.LE.2) GO TO 203                                          CPP26800
      WRITE(6,4)                                                        CPP26810
      CALL MATOUT(E2M,N3N,N3N,N3N,N3N,6)                                CPP26820
C                                                                       CPP26830
  203 CONTINUE                                                          CPP26840
      CALL SREW(ITAP44)                                                 CPP26850
      RETURN                                                            CPP26860
      END                                                               CPP26870
      FUNCTION ABIJ(I,J,K,L,ALPA,BETA,CC,NNC)                           CPP26880
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP26890
      DIMENSION ALPA(NTRI),BETA(NTRI),CC(NNC)                           CPP26900
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP26910
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP26920
      DATA TWO / 2.0D+00 /                                              CPP26930
C                                                                       CPP26940
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPP26950
      KL=IOFF(MAX0(K,L))+MIN0(K,L)                                      CPP26960
      IL=IOFF(MAX0(I,L))+MIN0(I,L)                                      CPP26970
      JK=IOFF(MAX0(J,K))+MIN0(J,K)                                      CPP26980
      IK=IOFF(MAX0(I,K))+MIN0(I,K)                                      CPP26990
      JL=IOFF(MAX0(J,L))+MIN0(J,L)                                      CPP27000
      IJKL=IOFF(MAX0(IJ,KL))+MIN0(IJ,KL)                                CPP27010
      ILJK=IOFF(MAX0(IL,JK))+MIN0(IL,JK)                                CPP27020
      IKJL=IOFF(MAX0(IK,JL))+MIN0(IK,JL)                                CPP27030
      ABIJ=ALPA(JL)*CC(IJKL)*TWO+BETA(JL)*(CC(IKJL)+CC(ILJK))           CPP27040
      RETURN                                                            CPP27050
      END                                                               CPP27060
      SUBROUTINE SCF2ND(E2A,HCIA,CIA,E2M,E2C,E2T,SA,EIJ)                CPP27070
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP27080
      DIMENSION HCIA(NCTRI,N3N),CIA(NCPEX,N3N)                          CPP27090
      DIMENSION E2A(N3N,N3N),E2M(N3N,N3N),E2C(N3N,N3N),E2T(N3N,N3N)     CPP27100
      DIMENSION SA(NTRI,N3N),EIJ(NBASIS,NBASIS)                         CPP27110
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP27120
      COMMON/CIOFF/IIUNQ,IJUNQ,ICIUNQ,IPX(5050),JPX(5050)               CPP27130
      COMMON/CIPOS/ICITYP(5050),ICPOS(200,100)                          CPP27140
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP27150
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP27160
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP27170
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP27180
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP27190
      COMMON/TAPES/ITAP35,ITAP42,ITAP43,ITAP44,ITAP47                   CPP27200
      COMMON/WORKS/JTAPE1,JTAPE2                                        CPP27210
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP27220
      COMMON/CI103/CI(100)                                              CPP27230
      COMMON/CI104/IEIJ(5050)                                           CPP27240
      DATA A00,TWO,FOUR / 0.0D+00 , 2.0D+00 , 4.0D+00 /                 CPP27250
    1 FORMAT(///                                                        CPP27260
     1 2X,'::::::::::::::::::::::::::::'/                               CPP27270
     2 2X,':::SUMMARY OF CALCULATION:::'/                               CPP27280
     3 2X,'::::::::::::::::::::::::::::'/)                              CPP27290
    2 FORMAT(//,2X,' E2A MATRIX'/)                                      CPP27300
    3 FORMAT(//,2X,' E2M MATRIX'/)                                      CPP27310
    4 FORMAT(//,2X,' E2C MATRIX'/)                                      CPP27320
    5 FORMAT(//,2X,' SCF SECOND DERIVATIVE MATRIX'/)                    CPP27330
    6 FORMAT(2I5)                                                       CPP27340
    7 FORMAT(3F20.10)                                                   CPP27350
C                                                                       CPP27360
C   CALCULATE SCF SECOND DERIVATIVES                                    CPP27370
C   SEE EQ.(12)                                                         CPP27380
C                                                                       CPP27390
C***********************************************                        CPP27400
C***CONTRIBUTION DUE TO CI COEFFICIENT CHANGE***                        CPP27410
C***********************************************                        CPP27420
C                                                                       CPP27430
      CALL SREW(ITAP44)                                                 CPP27440
      DO 101 IABC=1,N3N                                                 CPP27450
      IS=ISA(IABC)                                                      CPP27460
      CALL RREAD(ITAP44,SA(1,IABC),NTRI2,IS)                            CPP27470
  101 CONTINUE                                                          CPP27480
C                                                                       CPP27490
      CALL SREW(JTAPE2)                                                 CPP27500
      CALL ZERO(E2C,N3N*N3N)                                            CPP27510
      DO 105 ICI=1,NCPEX                                                CPP27520
      DO 105 JCI=1,NCPEX                                                CPP27530
      IJCI=IOFF(MAX0(ICI,JCI))+MIN0(ICI,JCI)                            CPP27540
      DO 102 III=1,ICIUNQ                                               CPP27550
      IEX=IPX(III)                                                      CPP27560
      JEX=JPX(III)                                                      CPP27570
      IJEX=IOFF(MAX0(IEX,JEX))+MIN0(IEX,JEX)                            CPP27580
      IF(IJEX.EQ.IJCI) GO TO 301                                        CPP27590
  102 CONTINUE                                                          CPP27600
      GO TO 105                                                         CPP27610
C                                                                       CPP27620
  301 CONTINUE                                                          CPP27630
      IBCI=ICITYP(III)                                                  CPP27640
      IPT=IEIJ(IBCI)                                                    CPP27650
      CALL RREAD(JTAPE2,EIJ,NBASQ,IPT)                                  CPP27660
C                                                                       CPP27670
      DO 104 IABC=1,N3N                                                 CPP27680
      DO 104 JABC=1,N3N                                                 CPP27690
      VALU0=A00                                                         CPP27700
      DO 103 I=1,JOCC                                                   CPP27710
      DO 103 J=1,JOCC                                                   CPP27720
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPP27730
C*    VALU0=VALU0-SA(IJ,JABC)*EIJ(I,J)                                  CPP27740
      VALU0=VALU0-SA(IJ,JABC)*EIJ(J,I)                                  CPP27750
  103 CONTINUE                                                          CPP27760
      VALUC=CIA(ICI,IABC)*(HCIA(IJCI,JABC)+VALU0*TWO)*CI(JCI)           CPP27770
      E2C(IABC,JABC)=E2C(IABC,JABC)+VALUC*TWO                           CPP27780
  104 CONTINUE                                                          CPP27790
  105 CONTINUE                                                          CPP27800
C                                                                       CPP27810
      WRITE(6,1)                                                        CPP27820
      WRITE(6,2)                                                        CPP27830
      CALL MATOUT(E2A,N3N,N3N,N3N,N3N,6)                                CPP27840
      WRITE(6,3)                                                        CPP27850
      CALL MATOUT(E2M,N3N,N3N,N3N,N3N,6)                                CPP27860
      WRITE(6,4)                                                        CPP27870
      CALL MATOUT(E2C,N3N,N3N,N3N,N3N,6)                                CPP27880
      DO 106 I=1,N3N                                                    CPP27890
      DO 106 J=1,N3N                                                    CPP27900
  106 E2T(I,J)=E2A(I,J)+E2M(I,J)+E2C(I,J)                               CPP27910
      WRITE(6,5)                                                        CPP27920
      CALL MATOUT(E2T,N3N,N3N,N3N,N3N,6)                                CPP27930
C                                                                       CPP27940
      ITAP15=15                                                         CPP27950
      REWIND ITAP15                                                     CPP27960
      N6N=N3N*2                                                         CPP27970
      WRITE(ITAP15,6) NATOMS,N6N                                        CPP27980
      WRITE(ITAP15,7) ((E2T(I,J),J=1,N3N),I=1,N3N)                      CPP27990
      REWIND ITAP15                                                     CPP28000
C                                                                       CPP28010
      CALL SREW(ITAP44)                                                 CPP28020
      RETURN                                                            CPP28030
      END                                                               CPP28040
      SUBROUTINE CIGRAD(E1A,E1T,D1W,U,W)                                CPP28050
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP28060
      DIMENSION E1A(N3N),E1T(N3N),D1W(N3N)                              CPP28070
      DIMENSION U(NBASIS,NBASIS),W(NBASIS,NBASIS)                       CPP28080
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP28090
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP28100
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP28110
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP28120
      COMMON/TAPES/ITAP35,ITAP42,ITAP43,ITAP44,ITAP47                   CPP28130
      DATA A00 / 0.0D+00 /                                              CPP28140
    1 FORMAT(//,2X,' LAGRANGIAN MATRIX'/)                               CPP28150
    2 FORMAT(//,2X,' U MATRIX, IABC = ',I5/)                            CPP28160
    3 FORMAT(//,2X,' EA1 MATRIX'/)                                      CPP28170
    4 FORMAT(//,2X,' EA2 MATRIX'/)                                      CPP28180
    5 FORMAT(//,2X,' E1T MATRIX'/)                                      CPP28190
C                                                                       CPP28200
      CALL MREAD(W,32)                                                  CPP28210
      IF(IPRNT.LE.2) GO TO 301                                          CPP28220
      WRITE(6,1)                                                        CPP28230
      CALL MATOUT(W,NBASIS,NBASIS,NBASIS,NBASIS,6)                      CPP28240
C                                                                       CPP28250
  301 CONTINUE                                                          CPP28260
      CALL SREW(ITAP44)                                                 CPP28270
      DO 105 IABC=1,N3N                                                 CPP28280
      IU=IUA(IABC)                                                      CPP28290
      CALL RREAD(ITAP44,U,NBASQ,IU)                                     CPP28300
      IF(IPRNT.LE.2) GO TO 302                                          CPP28310
      WRITE(6,2) IABC                                                   CPP28320
      CALL MATOUT(U,NBASIS,NBASIS,NBASIS,NBASIS,6)                      CPP28330
  302 VALU=A00                                                          CPP28340
      DO 103 I=1,NBASIS                                                 CPP28350
      DO 103 J=1,NBASIS                                                 CPP28360
  103 VALU=VALU+U(J,I)*W(I,J)                                           CPP28370
      D1W(IABC)=VALU                                                    CPP28380
  105 CONTINUE                                                          CPP28390
      WRITE(6,3)                                                        CPP28400
      CALL MATOUT(E1A,3,NATOMS,3,NATOMS,6)                              CPP28410
      WRITE(6,4)                                                        CPP28420
      CALL MATOUT(D1W,3,NATOMS,3,NATOMS,6)                              CPP28430
C                                                                       CPP28440
C   SUM UP ALL CONTRIBUTIONS                                            CPP28450
      DO 107 IABC=1,N3N                                                 CPP28460
  107 E1T(IABC)=E1A(IABC)+D1W(IABC)                                     CPP28470
      WRITE(6,5)                                                        CPP28480
      CALL MATOUT(E1T,3,NATOMS,3,NATOMS,6)                              CPP28490
C                                                                       CPP28500
      CALL SREW(ITAP44)                                                 CPP28510
      RETURN                                                            CPP28520
      END                                                               CPP28530
      SUBROUTINE DIPMO(HCIA,CIA,HF,UA)                                  CPP28540
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP28550
      DIMENSION HF(NTRI,3),UA(NBASIS,NBASIS,N3N)                        CPP28560
      DIMENSION HCIA(NCTRI,N3NP3),CIA(NCPEX,N3NP3)                      CPP28570
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP28580
      COMMON/CALIF/LPARA(1024),DIPDE(3,45),DIPDC(3,45),DUM(754)         CPP28590
      COMMON/CIOFF/IIUNQ,IJUNQ,ICIUNQ,IPX(5050),JPX(5050)               CPP28600
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP28610
      COMMON/DIPOL/DIPDA(3,45),DIPDM(3,45),DIPDN(3,45),DIPDT(3,45)      CPP28620
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP28630
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP28640
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP28650
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP28660
      COMMON/TAPES/ITAP35,ITAP42,ITAP43,ITAP44,ITAP47                   CPP28670
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP28680
      COMMON/CI103/CI(100)                                              CPP28690
      DATA A00,HALF,ONE,TWO,FOUR /  0.0D+00 , 0.5D+00 , 1.0D+00 ,       CPP28700
     1                              2.0D+00 , 4.0D+00 /                 CPP28710
      DATA DEBYE,BOHR / 2.541765480D+00 , 0.52917706D+00 /              CPP28720
    1 FORMAT(//,2X,' FIRST DERIVATIVES OF DIPOLE MOMENTS'/)             CPP28730
    2 FORMAT(//,2X,' UA MATRIX, IABC = ',I5,' IXYZ = ',I5/)             CPP28740
    3 FORMAT(//,2X,' DIPDA MATRIX'/)                                    CPP28750
    4 FORMAT(//,2X,' DIPDM MATRIX'/)                                    CPP28760
    5 FORMAT(//,2X,' DIPDC MATRIX'/)                                    CPP28770
    6 FORMAT(//,2X,' DIPDE MATRIX'/)                                    CPP28780
    7 FORMAT(//,2X,' DIPDN MATRIX'/)                                    CPP28790
    8 FORMAT(//,2X,' FIRST DERIVATIVES OF DIPOLE MOMENTS (IN DEBYE/BOHR)CPP28800
     *'/)                                                               CPP28810
    9 FORMAT(//,2X,' FIRST DERIVATIVES OF DIPOLE MOMENTS (IN DEBYE/A)'/)CPP28820
   10 FORMAT(2I5)                                                       CPP28830
   11 FORMAT(3F20.10)                                                   CPP28840
C                                                                       CPP28850
      CALL SREW(ITAP44)                                                 CPP28860
      DEB4=DEBYE*FOUR                                                   CPP28870
      DEB2=DEBYE*TWO                                                    CPP28880
      RBOHR=ONE/BOHR                                                    CPP28890
C     WRITE(6,1)                                                        CPP28900
C                                                                       CPP28910
      CALL ZERO(DIPDM,N3N*3)                                            CPP28920
      CALL ZERO(DIPDC,N3N*3)                                            CPP28930
C                                                                       CPP28940
      DO 101 IXYZ=1,3                                                   CPP28950
      IH=IHA(IXYZ+N3N)                                                  CPP28960
      CALL RREAD(ITAP44,HF(1,IXYZ),NTRI2,IH)                            CPP28970
  101 CONTINUE                                                          CPP28980
      DO 102 IABC=1,N3N                                                 CPP28990
      IU=IUA(IABC)                                                      CPP29000
      CALL RREAD(ITAP44,UA(1,1,IABC),NBASQ,IU)                          CPP29010
      IF(IPRNT.LE.3) GO TO 102                                          CPP29020
      WRITE(6,2) IABC                                                   CPP29030
      CALL MATOUT(UA(1,1,IABC),NBASIS,NBASIS,NBASIS,NBASIS,6)           CPP29040
  102 CONTINUE                                                          CPP29050
C                                                                       CPP29060
      DO 106 IABC=1,N3N                                                 CPP29070
      DO 105 IXYZ=1,3                                                   CPP29080
      VALUM=A00                                                         CPP29090
      DO 103 I=1,JOCC                                                   CPP29100
      ITYP=MOTYP(I)                                                     CPP29110
      FAC=FOCC(ITYP)*HALF                                               CPP29120
      DO 103 J=1,NBASIS                                                 CPP29130
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPP29140
      VALUM=VALUM+UA(J,I,IABC)*HF(IJ,IXYZ)*FAC                          CPP29150
  103 CONTINUE                                                          CPP29160
      VALUC=A00                                                         CPP29170
      DO 104 IJCI=1,ICIUNQ                                              CPP29180
      IEX=IPX(IJCI)                                                     CPP29190
      JEX=JPX(IJCI)                                                     CPP29200
      FAC=TWO                                                           CPP29210
      IF(IEX.EQ.JEX) FAC=ONE                                            CPP29220
      IJEX=IOFF(MAX0(IEX,JEX))+MIN0(IEX,JEX)                            CPP29230
      VALUC=VALUC+CIA(IEX,IABC)*CI(JEX)*HCIA(IJEX,IXYZ+N3N)*FAC         CPP29240
  104 CONTINUE                                                          CPP29250
      DIPDM(IXYZ,IABC)=VALUM                                            CPP29260
      DIPDC(IXYZ,IABC)=VALUC                                            CPP29270
  105 CONTINUE                                                          CPP29280
  106 CONTINUE                                                          CPP29290
C                                                                       CPP29300
      DO 107 I=1,3                                                      CPP29310
      DO 107 J=1,N3N                                                    CPP29320
      DIPDM(I,J)=-DIPDM(I,J)*DEB4                                       CPP29330
      DIPDC(I,J)=-DIPDC(I,J)*DEB2                                       CPP29340
  107 CONTINUE                                                          CPP29350
C                                                                       CPP29360
C   SUM UP CONTRIBUTIONS FROM MO CHANGE                                 CPP29370
      DO 108 I=1,3                                                      CPP29380
      DO 108 J=1,N3N                                                    CPP29390
      DIPDE(I,J)=DIPDA(I,J)+DIPDM(I,J)+DIPDC(I,J)                       CPP29400
  108 CONTINUE                                                          CPP29410
C*    IF(IPRNT.LE.2) GO TO 201                                          CPP29420
      WRITE(6,3)                                                        CPP29430
      CALL MATOUT(DIPDA,3,45,3,N3N,6)                                   CPP29440
      WRITE(6,4)                                                        CPP29450
      CALL MATOUT(DIPDM,3,45,3,N3N,6)                                   CPP29460
      WRITE(6,5)                                                        CPP29470
      CALL MATOUT(DIPDC,3,45,3,N3N,6)                                   CPP29480
C                                                                       CPP29490
C   SUM UP ALL CONTRIBUTIONS                                            CPP29500
  201 CONTINUE                                                          CPP29510
      DO 109 I=1,3                                                      CPP29520
      DO 109 J=1,N3N                                                    CPP29530
      DIPDT(I,J)=DIPDE(I,J)+DIPDN(I,J)                                  CPP29540
  109 CONTINUE                                                          CPP29550
C                                                                       CPP29560
C*    IF(IPRNT.LE.2) GO TO 202                                          CPP29570
      WRITE(6,6)                                                        CPP29580
      CALL MATOUT(DIPDE,3,45,3,N3N,6)                                   CPP29590
      WRITE(6,7)                                                        CPP29600
      CALL MATOUT(DIPDN,3,45,3,N3N,6)                                   CPP29610
  202 CONTINUE                                                          CPP29620
      WRITE(6,8)                                                        CPP29630
      CALL MATOUT(DIPDT,3,45,3,N3N,6)                                   CPP29640
C                                                                       CPP29650
      DO 110 I=1,3                                                      CPP29660
      DO 110 J=1,N3N                                                    CPP29670
      DIPDT(I,J)=DIPDT(I,J)*RBOHR                                       CPP29680
  110 CONTINUE                                                          CPP29690
      WRITE(6,9)                                                        CPP29700
      CALL MATOUT(DIPDT,3,45,3,N3N,6)                                   CPP29710
C                                                                       CPP29720
      ITAP17=17                                                         CPP29730
      REWIND ITAP17                                                     CPP29740
      WRITE(ITAP17,10) NATOMS,N3N                                       CPP29750
      WRITE(ITAP17,11) ((DIPDT(I,J),J=1,N3N),I=1,3)                     CPP29760
      REWIND ITAP17                                                     CPP29770
C                                                                       CPP29780
      CALL SREW(ITAP44)                                                 CPP29790
      RETURN                                                            CPP29800
      END                                                               CPP29810
      SUBROUTINE POLAR(HCIA,CIA,HF,UF)                                  CPP29820
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP29830
      DIMENSION HF(NTRI,3),UF(NBASIS,NBASIS,3)                          CPP29840
      DIMENSION HCIA(NCTRI,N3NP3),CIA(NCPEX,N3NP3)                      CPP29850
      DIMENSION POM(3,3),POC(3,3)                                       CPP29860
      DIMENSION POL(3,3),POLC(3,3),EIG(3),EIV(3,3),A(6),FV1(3),FV2(3)   CPP29870
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP29880
      COMMON/CIOFF/IIUNQ,IJUNQ,ICIUNQ,IPX(5050),JPX(5050)               CPP29890
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP29900
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP29910
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP29920
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP29930
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP29940
      COMMON/TAPES/ITAP35,ITAP42,ITAP43,ITAP44,ITAP47                   CPP29950
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP29960
      COMMON/CI103/CI(100)                                              CPP29970
      DATA A00,HALF,ONE,TWO,FOUR / 0.0D+00 , 0.5D+00 , 1.0D+00 ,        CPP29980
     1                             2.0D+00 , 4.0D+00 /                  CPP29990
      DATA BOHR / 0.52917706D+00 /                                      CPP30000
    1 FORMAT(//,2X,' HF MATRIX, IXYZ = ',I5/)                           CPP30010
    2 FORMAT(//,2X,' UF MATRIX, IXYZ = ',I5/)                           CPP30020
    3 FORMAT(//,2X,' POM MATRIX'/)                                      CPP30030
    4 FORMAT(//,2X,' POC MATRIX'/)                                      CPP30040
    5 FORMAT(//,2X,' POLARIZABILITY TENSOR MATRIX (IN A.U.)'/)          CPP30050
    6 FORMAT(//,2X,' POLARIZABILITY TENSOR MATRIX (IN A**3)'/)          CPP30060
    7 FORMAT(//,2X,' A MATRIX'/)                                        CPP30070
    8 FORMAT(//,2X,' PRINCIPAL POLARIZABILITY TENSOR MATRIX (IN A.U.)'/)CPP30080
    9 FORMAT(//,2X,' PRINCIPAL POLARIZABILITY TENSOR MATRIX (IN A**3)'/)CPP30090
   10 FORMAT(//,2X,' EIGENVECTOR MATRIX'/)                              CPP30100
C                                                                       CPP30110
      CALL SREW(ITAP44)                                                 CPP30120
      A33=BOHR*BOHR*BOHR                                                CPP30130
C                                                                       CPP30140
C   CALL BACK HF MATRICES                                               CPP30150
      DO 101 IXYZ=1,3                                                   CPP30160
      IH=IHA(IXYZ+N3N)                                                  CPP30170
      IU=IUA(IXYZ+N3N)                                                  CPP30180
      CALL RREAD(ITAP44,HF(1,IXYZ),NTRI2,IH)                            CPP30190
      CALL RREAD(ITAP44,UF(1,1,IXYZ),NBASQ,IU)                          CPP30200
      IF(IPRNT.LE.3) GO TO 101                                          CPP30210
      WRITE(6,1) IXYZ                                                   CPP30220
      CALL PRINT(HF(1,IXYZ),NTRI,NBASIS,6)                              CPP30230
      WRITE(6,2) IXYZ                                                   CPP30240
      CALL MATOUT(UF(1,1,IXYZ),NBASIS,NBASIS,NBASIS,NBASIS,6)           CPP30250
  101 CONTINUE                                                          CPP30260
C                                                                       CPP30270
C:::B (G) COORDINATE:::                                                 CPP30280
      DO 105 IXYZ=1,3                                                   CPP30290
C                                                                       CPP30300
C:::A (F) COORDINATES:::                                                CPP30310
  301 CONTINUE                                                          CPP30320
      DO 104 JXYZ=1,3                                                   CPP30330
      VALUP=A00                                                         CPP30340
      DO 102 I=1,JOCC                                                   CPP30350
      ITYP=MOTYP(I)                                                     CPP30360
      FAC=FOCC(ITYP)*HALF                                               CPP30370
      DO 102 J=1,NBASIS                                                 CPP30380
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPP30390
      VALUP=VALUP-UF(J,I,IXYZ)*HF(IJ,JXYZ)*FAC                          CPP30400
  102 CONTINUE                                                          CPP30410
      VALUC=A00                                                         CPP30420
      DO 103 IJCI=1,NCPEX                                               CPP30430
      IEX=IPX(IJCI)                                                     CPP30440
      JEX=JPX(IJCI)                                                     CPP30450
      FAC=TWO                                                           CPP30460
      IF(IEX.EQ.JEX) FAC=ONE                                            CPP30470
      IJEX=IOFF(MAX0(IEX,JEX))+MIN0(IEX,JEX)                            CPP30480
      VALUC=VALUC-CIA(IEX,IXYZ+N3N)*CI(JEX)*HCIA(IJEX,JXYZ+N3N)*FAC     CPP30490
  103 CONTINUE                                                          CPP30500
      POM(IXYZ,JXYZ)=VALUP*FOUR                                         CPP30510
      POC(IXYZ,JXYZ)=VALUC*TWO                                          CPP30520
  104 CONTINUE                                                          CPP30530
  105 CONTINUE                                                          CPP30540
C                                                                       CPP30550
      DO 106 IXYZ=1,3                                                   CPP30560
      DO 106 JXYZ=1,3                                                   CPP30570
      POL(IXYZ,JXYZ)=POM(IXYZ,JXYZ)+POC(IXYZ,JXYZ)                      CPP30580
  106 CONTINUE                                                          CPP30590
      WRITE(6,3)                                                        CPP30600
      CALL MATOUT(POM,3,3,3,3,6)                                        CPP30610
      WRITE(6,4)                                                        CPP30620
      CALL MATOUT(POC,3,3,3,3,6)                                        CPP30630
      WRITE(6,5)                                                        CPP30640
      CALL MATOUT(POL,3,3,3,3,6)                                        CPP30650
      DO 107 I=1,3                                                      CPP30660
      DO 107 J=1,3                                                      CPP30670
      POLC(I,J)=POL(I,J)*A33                                            CPP30680
  107 CONTINUE                                                          CPP30690
      WRITE(6,6)                                                        CPP30700
      CALL MATOUT(POLC,3,3,3,3,6)                                       CPP30710
C                                                                       CPP30720
      IJ=0                                                              CPP30730
      DO 108 I=1,3                                                      CPP30740
      DO 108 J=1,I                                                      CPP30750
      IJ=IJ+1                                                           CPP30760
      A(IJ)=POL(I,J)                                                    CPP30770
  108 CONTINUE                                                          CPP30780
      WRITE(6,7)                                                        CPP30790
      CALL PRINT(A,6,3,6)                                               CPP30800
      CALL RSP(3,3,6,A,EIG,1,EIV,FV1,FV2)                               CPP30810
C                                                                       CPP30820
      CALL ZERO(POL,9)                                                  CPP30830
      DO 110 I=1,3                                                      CPP30840
  110 POL(I,I)=EIG(I)                                                   CPP30850
      WRITE(6,8)                                                        CPP30860
      CALL MATOUT(POL,3,3,3,3,6)                                        CPP30870
      DO 111 I=1,3                                                      CPP30880
      DO 111 J=1,3                                                      CPP30890
      POLC(I,J)=POL(I,J)*A33                                            CPP30900
  111 CONTINUE                                                          CPP30910
      WRITE(6,9)                                                        CPP30920
      CALL MATOUT(POLC,3,3,3,3,6)                                       CPP30930
      WRITE(6,10)                                                       CPP30940
      CALL MATOUT(EIV,3,3,3,3,6)                                        CPP30950
C                                                                       CPP30960
      CALL SREW(ITAP44)                                                 CPP30970
      RETURN                                                            CPP30980
      END                                                               CPP30990
      SUBROUTINE ORB1ST(B0,CC,IFLAG)                                    CPP31000
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP31010
      DIMENSION B0(NTRI),CC(NBASIS,1)                                   CPP31020
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP31030
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP31040
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP31050
      COMMON/SIGNS/IOFF(1379),IPRNT                                     CPP31060
      COMMON/TAPES/ITAP35,ITAP42,ITAP43,ITAP44,ITAP47                   CPP31070
    1 FORMAT(//,2X,' FIRST DERIVATIVES OF ORBITAL ENERGIES W.R.T. NUCLEACPP31080
     1R COORDINATE'/)                                                   CPP31090
    2 FORMAT(//,2X,' FIRST DERIVATIVES OF ORBITAL ENERGIES W.R.T. ELECTRCPP31100
     1IC FIELD'/)                                                       CPP31110
C                                                                       CPP31120
C   FOR NUCLEAR COORDINATE PERTURBATION                                 CPP31130
      IF(IFLAG.EQ.2) GO TO 201                                          CPP31140
      CALL SREW(ITAP44)                                                 CPP31150
      DO 102 IABC=1,N3N                                                 CPP31160
      IB=IBA(IABC)                                                      CPP31170
      CALL RREAD(ITAP44,B0,NTRI2,IB)                                    CPP31180
      DO 101 I=1,NBASIS                                                 CPP31190
      II=IOFF(I+1)                                                      CPP31200
      CC(I,IABC)=B0(II)                                                 CPP31210
  101 CONTINUE                                                          CPP31220
  102 CONTINUE                                                          CPP31230
      WRITE(6,1)                                                        CPP31240
      CALL MATOUT(CC,NBASIS,N3N,NBASIS,N3N,6)                           CPP31250
      CALL SREW(ITAP44)                                                 CPP31260
      RETURN                                                            CPP31270
C                                                                       CPP31280
C   FOR ELECTRIC FIELD PERTURBATION                                     CPP31290
  201 CONTINUE                                                          CPP31300
      CALL SREW(ITAP44)                                                 CPP31310
      DO 104 IXYZ=1,3                                                   CPP31320
      IB=IBA(IXYZ+N3N)                                                  CPP31330
      CALL RREAD(ITAP44,B0,NTRI2,IB)                                    CPP31340
      DO 103 I=1,NBASIS                                                 CPP31350
      II=IOFF(I+1)                                                      CPP31360
      CC(I,IXYZ)=B0(II)                                                 CPP31370
  103 CONTINUE                                                          CPP31380
  104 CONTINUE                                                          CPP31390
      WRITE(6,2)                                                        CPP31400
      CALL MATOUT(CC,NBASIS,3,NBASIS,3,6)                               CPP31410
      CALL SREW(ITAP44)                                                 CPP31420
C                                                                       CPP31430
      RETURN                                                            CPP31440
      END                                                               CPP31450
      SUBROUTINE WRBMAT(BB,NABC)                                        CPP31460
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP31470
      DIMENSION BB(NTRI,NABC)                                           CPP31480
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP31490
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP31500
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP31510
      COMMON/TAPES/ITAP35,ITAP42,ITAP43,ITAP44,ITAP47                   CPP31520
    1 FORMAT(//,2X,' B0 MATRIX'/)                                       CPP31530
C                                                                       CPP31540
      CALL SREW(ITAP44)                                                 CPP31550
      DO 101 IABC=1,NABC                                                CPP31560
      IB=IBA(IABC)                                                      CPP31570
      CALL RREAD(ITAP44,BB(1,IABC),NTRI2,IB)                            CPP31580
  101 CONTINUE                                                          CPP31590
      WRITE(6,1)                                                        CPP31600
      CALL MATOUT(BB,NTRI,NABC,NTRI,NABC,6)                             CPP31610
      CALL SREW(ITAP44)                                                 CPP31620
C                                                                       CPP31630
      RETURN                                                            CPP31640
      END                                                               CPP31650
      SUBROUTINE THIRD(E1T,HCIA,CIA,NSECMX)                             CPP31660
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP31670
      DIMENSION E1T(N3N),HCIA(NCTRI*N3N),CIA(NCPEX*N3N)                 CPP31680
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP31690
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP31700
      COMMON/MFSEC/MASTER,NSECT                                         CPP31710
      COMMON/TAPES/ITAP35,ITAP42,ITAP43,ITAP44,ITAP47                   CPP31720
    1 FORMAT(//,2X,' E1T MATRIX IN THIRD'/)                             CPP31730
    2 FORMAT(//,2X,' HCIA MATRIX IN THIRD'/)                            CPP31740
    3 FORMAT(//,2X,' CIA MATRIX IN THIRD'/)                             CPP31750
C                                                                       CPP31760
      IF(IPRNT.LE.2) GO TO 201                                          CPP31770
      WRITE(6,1)                                                        CPP31780
      CALL MATOUT(E1T,3,NATOMS,3,NATOMS,6)                              CPP31790
      WRITE(6,2)                                                        CPP31800
      CALL MATOUT(HCIA,NCTRI,N3N,NCTRI,N3N,6)                           CPP31810
      WRITE(6,3)                                                        CPP31820
      CALL MATOUT(CIA,NCPEX,N3N,NCPEX,N3N,6)                            CPP31830
C                                                                       CPP31840
  201 CONTINUE                                                          CPP31850
      CALL SREW(ITAP44)                                                 CPP31860
      N3N2=N3N*2                                                        CPP31870
      NPOINT=NSECMX                                                     CPP31880
      CALL RWRIT(ITAP44,E1T,N3N2,NPOINT)                                CPP31890
      NSEC3L=((N3N2-1)/NSECT+1)                                         CPP31900
      NPOINT=NPOINT+NSEC3L                                              CPP31910
      CALL RWRIT(ITAP44,HCIA,N3N2*NCTRI,NPOINT)                         CPP31920
      NSECNC=((N3N2*NCTRI-1)/NSECT+1)                                   CPP31930
      NPOINT=NPOINT+NSECNC                                              CPP31940
      CALL RWRIT(ITAP44,CIA,N3N2*NCPEX,NPOINT)                          CPP31950
C                                                                       CPP31960
      CALL SREW(ITAP44)                                                 CPP31970
      RETURN                                                            CPP31980
      END                                                               CPP31990
      SUBROUTINE RDINT2(BUFO,NDIMB,ITAPE,MAXINT,LBLI,BUFI,IPRNT)        CPP32000
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP32010
      DIMENSION BUFO(NDIMB),LBLI(MAXBF4),BUFI(MAXBF2)                   CPP32020
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPP32030
      COMMON/SIGNS/IOFF(1379),JPRNT                                     CPP32040
    1 FORMAT(2X,4I6,5X,4I6,5X,F15.8)                                    CPP32050
C                                                                       CPP32060
C   READ IN TWO ELECTRON INTEGRALS                                      CPP32070
      CALL RFILE(ITAPE)                                                 CPP32080
      CALL ZERO(BUFO,MAXINT)                                            CPP32090
      INTI=0                                                            CPP32100
      NBLI=0                                                            CPP32110
  200 IBL=0                                                             CPP32120
      NBLI=NBLI+1                                                       CPP32130
      CALL SREAD(ITAPE,LBLI,MAXBF4)                                     CPP32140
  201 IBL=IBL+1                                                         CPP32150
      JA=LBLI(IBL+IBL-1)                                                CPP32160
      IF(JA.EQ.0) GO TO 205                                             CPP32170
      JAA=LBLI(IBL+IBL)                                                 CPP32180
      CALL UNPACK(I,J,K,L,IX,JA,JAA)                                    CPP32190
      IJ=IOFF(I)+J                                                      CPP32200
      KL=IOFF(K)+L                                                      CPP32210
C     IJKL=IJ*(IJ-1)/2+KL                                               CPP32220
      IJKL=IOFF(IJ)+KL                                                  CPP32230
      BUFO(IJKL)=BUFI(MAXBUF+IBL)                                       CPP32240
C***  IF(IPRNT.LE.5) GO TO 202                                          CPP32250
C***  WRITE(6,1) I,J,K,L,IJ,KL,IX,IJKL,BUFO(IJKL)                       CPP32260
  202 IF(IBL.EQ.MAXBUF) GO TO 203                                       CPP32270
      GO TO 201                                                         CPP32280
  203 INTI=INTI+MAXBUF                                                  CPP32290
      GO TO 200                                                         CPP32300
  205 INTI=INTI+IBL                                                     CPP32310
C                                                                       CPP32320
      CALL RCLOSE(ITAPE,3)                                              CPP32330
      RETURN                                                            CPP32340
      END                                                               CPP32350
