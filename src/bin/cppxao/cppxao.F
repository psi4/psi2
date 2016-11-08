      PROGRAM CPPXAO                                                    CPP00010
C   COUPLED PERTURBED HARTREE-FOCK PROGRAM IN AO BASIS                  CPP00020
C   FOR THE GVB AND PAIR EXCITATION SELF-CONSISTENT-FILED WAVEFUNCTIONS CPP00030
C   THIS PROGRAM EMPLOYS I/O WITH THE SUPERMATRICES                     CPP00040
C***********************************************************************CPP00050
C*  MODIFICATION FOR IMS VERSION                                       *CPP00060
C*  BY YUKIO YAMAGUCHI                                                 *CPP00070
C*  DATE FEBRUARY 01, 1989                                             *CPP00080
C***********************************************************************CPP00090
C***LAST UPDATED ON NOVEMBER 17, 1987 BY YUKIO YAMAGUCHI               *CPP00100
C***********************************************************************CPP00110
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP00120
CCC   DIMENSION CC(1650000),IA(1),I30(200),NADD(10)                     CPP00130
      DIMENSION CC(650000),IA(1),I30(200),NADD(10)                      CPP00140
      CHARACTER*80 ALABEL                                               CPP00150
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP00160
      COMMON/CALIF/LPARA(1024),APARA(1024)                              CPP00170
      COMMON/CIOFF/IIUNQ,IJUNQ,ICIUNQ,IPX(5050),JPX(5050)               CPP00180
      COMMON/CIPEX/PEXC(100),INDPEX(20,100),MOPEX(20),ICIMO(20)         CPP00190
      COMMON/CIPMO/IIMO(5050),JJMO(5050)                                CPP00200
      COMMON/CIPOS/ICITYP(5050),ICPOS(200,100)                          CPP00210
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP00220
      COMMON/CITYP/IJGRP,IBM(200),JBM(200),IBNUM(200)                   CPP00230
      COMMON/DENST/ALPC(465),BETC(465)                                  CPP00240
      COMMON/DIPOL/DIPDA(3,45),DIPDM(3,45),DIPDN(3,45),DIPDT(3,45)      CPP00250
      COMMON/ENRGY/ESCF,ENUC,ELEC                                       CPP00260
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP00270
      COMMON/GVBAB/ALPA(30,30),BETA(30,30)                              CPP00280
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPP00290
      COMMON/MFSEC/MASTER,NSECT                                         CPP00300
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP00310
      COMMON/RSTRT/IGSTOP,ILSTOP                                        CPP00320
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP00330
      COMMON/SIZES/MAXVEC,MAXVX                                         CPP00340
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP00350
      COMMON/SUPER/MAXIJ                                                CPP00360
      COMMON/TAPES/ITAP37,ITAP42,ITAP43,ITAP44,ITAP47                   CPP00370
      COMMON/TOLER/CRIT,TOLR                                            CPP00380
      COMMON/WORKS/JTAPE1,JTAPE2                                        CPP00390
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP00400
      COMMON/CI102/NIND,NDEP                                            CPP00410
      COMMON/CI103/CI(100)                                              CPP00420
      COMMON/CI104/IEIJ(5050)                                           CPP00430
      EQUIVALENCE (CC,IA)                                               CPP00440
      DATA ALABEL / ' SCF FIRST DERIVATIVE BROUGHT YOU BY Y & Y CO. EST.CPP00450
     1 IN MARCH , 1982             ' /                                  CPP00460
    1 FORMAT(//,2X,' COUPLED PERTURBED HARTREE-FOCK PROGRAM FOR GVB-SCF CPP00470
     1IN THE AO BASIS (VERSION OF 2/01/89)'/)                           CPP00480
    2 FORMAT(7I5)                                                       CPP00490
    3 FORMAT(A8)                                                        CPP00500
    4 FORMAT(2X,10I5)                                                   CPP00510
    5 FORMAT(//,2X,' PARAMETERS'/                                       CPP00520
     * 2X,' NBASIS = ',I8/                                              CPP00530
     * 2X,' NTRI   = ',I8/                                              CPP00540
     * 2X,' NBFAO  = ',I8/                                              CPP00550
     * 2X,' NBATRI = ',I8/                                              CPP00560
     * 2X,' NATOMS = ',I8/                                              CPP00570
     * 2X,' N3N    = ',I8/                                              CPP00580
     * 2X,' N3NP3  = ',I8/                                              CPP00590
     * 2X,' N3NXX  = ',I8/                                              CPP00600
     * 2X,' IOCC   = ',I8/                                              CPP00610
     * 2X,' JOCC   = ',I8/                                              CPP00620
     * 2X,' KOCC   = ',I8/                                              CPP00630
     * 2X,' NTYPES = ',I8/                                              CPP00640
     * 2X,' NTYPEP = ',I8/                                              CPP00650
     * 2X,' NTREAD = ',I8/                                              CPP00660
     * 2X,' NTYTRI = ',I8)                                              CPP00670
    6 FORMAT(2X,' NCPEX  = ',I8/                                        CPP00680
     * 2X,' NORPEX = ',I8/                                              CPP00690
     * 2X,' NCTRI  = ',I8/                                              CPP00700
     * 2X,' LVPEX  = ',I8/                                              CPP00710
     * 2X,' IPOL   = ',I8/                                              CPP00720
     * 2X,' ICONV  = ',I8/                                              CPP00730
     * 2X,' IPRNT  = ',I8)                                              CPP00740
    7 FORMAT(2X,' ISTEP  = ',I8/                                        CPP00750
     * 2X,' IGSTOP = ',I8/                                              CPP00760
     * 2X,' ILSTOP = ',I8/)                                             CPP00770
    8 FORMAT(//,2X,' SCF EIGENVECTOR (AO BASIS) MATRIX'/)               CPP00780
    9 FORMAT(//,2X,' ALPA MATRIX'/)                                     CPP00790
   10 FORMAT(//,2X,' BETA MATRIX'/)                                     CPP00800
   11 FORMAT(//,2X,' SCF EIGENVECTOR (SO BASIS) MATRIX'/)               CPP00810
   12 FORMAT(//,2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)                CPP00820
   13 FORMAT(//,2X,' REQUIRED MEMORY EXCEEDS MAXCOR'/                   CPP00830
     1          2X,' ICMAX = ',I7,5X,' MAXCOR = ',I7/)                  CPP00840
   14 FORMAT(//,2X,' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'/   CPP00850
     1          2X,' %%%THE CPHF STEP HAS ALREADY BEEN COMPLETED%%%'/   CPP00860
     2          2X,' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'/)  CPP00870
C                                                                       CPP00880
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::            CPP00890
C   REFERENCES FOR TCSCF AND GRSCF :                       :            CPP00900
C   Y.YAMAGUCHI, Y.OSAMURA, AND H.F.SCHAEFER III,          :            CPP00910
C   J.AM.CHEM.SOC. 105,7506(1983)                          :            CPP00920
C                     &                                    :            CPP00930
C   Y.OSAMURA, Y.YAMAGUCHI, P.SAXE, D.J.FOX, M.A.VINCENT,  :            CPP00940
C   AND H.F.SCHAEFER III, J.MOL.STRUC. 103,183 (1983)      :            CPP00950
C                     &                                    :            CPP00960
C   Y.OSAMURA, Y.YAMAGUCHI, P.SAXE, M.A.VINCENT, J.F.GAW,  :            CPP00970
C   AND H.F.SCHAEFER III, CHEM.PHYS. 72,131 (1982)         :            CPP00980
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::            CPP00990
C                                                                       CPP01000
C   ITAP11 = GEOMETRY AND GRADIENTS                                     CPP01010
C   ITAP15 = SCF SECOND DERIVATIVES                                     CPP01020
C   ITAP17 = DIPOLE MOMENT DERIVATIVES                                  CPP01030
C   ITAP37 = TWO ELECTRON AO INTEGRALS (PK-FILE)                        CPP01040
C   ITAP42 = FIRST AND SECOND DERIVATIVE INFORMATION                    CPP01050
C   ITAP43 = DERIVATIVE DIPOLE MOMENTS                                  CPP01060
C   ITAP44 = FIRST DERIVATIVE                                           CPP01070
C   ITAP47 = THE "BARE" GENERALIZED LAGRANGIAN MATRICES                 CPP01080
C   JTAPE1 = WORKING TAPE FOR CPHF                                      CPP01090
C   JTAPE2 = WORKING TAPE FOR E MATRICES                                CPP01100
C                                                                       CPP01110
      CALL TSTART(6)                                                    CPP01120
      CALL NOUNFL                                                       CPP01130
      CALL INITMF(1)                                                    CPP01140
      CALL RCLOSE(MASTER,3)                                             CPP01150
C                                                                       CPP01160
C::::::::::::::::::::::::::::::::::::::::                               CPP01170
C:::DEFINE MACHINE-DEPENDENT CONSTANTS:::                               CPP01180
C::::::::::::::::::::::::::::::::::::::::                               CPP01190
CCC   MAXCOR=1650000                                                    CPP01200
      MAXCOR=650000                                                     CPP01210
      MAXVEC=600                                                        CPP01220
      MAXVX=30                                                          CPP01230
C                                                                       CPP01240
      ITAPE3=3                                                          CPP01250
      ITAPE4=4                                                          CPP01260
      INPUT=5                                                           CPP01270
      ITAPE6=6                                                          CPP01280
      ITAP11=11                                                         CPP01290
      ITAP15=15                                                         CPP01300
      ITAP17=17                                                         CPP01310
      ITAP30=30                                                         CPP01320
      ITAP37=37                                                         CPP01330
      ITAP42=42                                                         CPP01340
      ITAP43=43                                                         CPP01350
      ITAP44=44                                                         CPP01360
      ITAP47=47                                                         CPP01370
      JTAPE1=91                                                         CPP01380
      JTAPE2=92                                                         CPP01390
C                                                                       CPP01400
C   LOCATION OF MATRICES IN ITAP44                                      CPP01410
C   (1) SA (N3NXX*NTRIL)                                                CPP01420
C   (2) HA (N3NXX*NTRIL)                                                CPP01430
C   (3) TA AND FA (N3NXX*NTRIL)                                         CPP01440
C   (4) EA (N3NXX*NBASQL)                                               CPP01450
C   (5) BA (N3NXX*NTRIL)                                                CPP01460
C   (6) UA (N3NXX*NBASQL)                                               CPP01470
C                                                                       CPP01480
      IOFF(1)=0                                                         CPP01490
      DO 101 I=1,255                                                    CPP01500
  101 IOFF(I+1)=IOFF(I)+I                                               CPP01510
C                                                                       CPP01520
      CALL LOCATE(INPUT,'# CPHFAO #',IERR)                              CPP01530
      WRITE(6,1)                                                        CPP01540
      WRITE(3,1)                                                        CPP01550
      DO 102 I=1,10                                                     CPP01560
  102 WRITE(4,*) ' '                                                    CPP01570
      WRITE(4,*) '  *************************************************'  CPP01580
      WRITE(4,*) '  ***THE COUPLED PERTUBED HATREE-FOCK PROGRAM******'  CPP01590
      WRITE(4,*) '  ***FOR THE GVB AND PEX SCF WAVE FUNCTIONS********'  CPP01600
      WRITE(4,*) '  *************************************************'  CPP01610
      DO 103 I=1,2                                                      CPP01620
  103 WRITE(4,*) ' '                                                    CPP01630
C                                                                       CPP01640
C   READ IN PARAMETERS                                                  CPP01650
      READ(5,2) IPOL,ICONV,IPRNT                                        CPP01660
C#######################################################################CPP01670
C   DEFINE CPHF CONVERGENCE TOLERANCE                                   CPP01680
      CRIT=1.0D-12                                                      CPP01690
      TOLR=1.0D-12                                                      CPP01700
      ICRIT=MOD(ICONV,100)                                              CPP01710
      IF(ICRIT.NE.0) CRIT=1.0D+00/(10.0D+00**ICRIT)                     CPP01720
      ITOLR=ICONV/100                                                   CPP01730
      IF(ITOLR.NE.0) TOLR=1.0D+00/(10.0D+00**ITOLR)                     CPP01740
      WRITE(3,16) CRIT,TOLR                                             CPP01750
   16 FORMAT(//,2X,' CRIT = ',D20.10/                                   CPP01760
     1          2X,' TOLR = ',D20.10/)                                  CPP01770
C#######################################################################CPP01780
C                                                                       CPP01790
C   READ IN SCF INFORMATION                                             CPP01800
      NBASIS=LPARA(3)                                                   CPP01810
      NBFAO=LPARA(11)                                                   CPP01820
      NTRI=LPARA(4)                                                     CPP01830
      NTRI2=NTRI*2                                                      CPP01840
      NBATRI=IOFF(NBFAO+1)                                              CPP01850
      NBATR2=NBATRI*2                                                   CPP01860
      IOCC=LPARA(7)                                                     CPP01870
      KOCC=LPARA(8)                                                     CPP01880
      JOCC=LPARA(9)                                                     CPP01890
      NATOMS=LPARA(10)                                                  CPP01900
      N3N=LPARA(17)                                                     CPP01910
      NATRI=IOFF(N3N+1)                                                 CPP01920
      N3NP3=N3N+3                                                       CPP01930
      N3NXX=N3NP3                                                       CPP01940
      IF(IPOL.NE.0) N3NXX=N3N                                           CPP01950
      NTRIL=((NTRI2-1)/NSECT+1)                                         CPP01960
      NBASQ=NBASIS*NBASIS*2                                             CPP01970
      NBASQL=((NBASQ-1)/NSECT+1)                                        CPP01980
      NSEC3L=((N3N*2-1)/NSECT+1)                                        CPP01990
      N3TRL=N3NXX*NTRIL                                                 CPP02000
      N3QRL=N3NXX*NBASQL                                                CPP02010
      NTYPES=LPARA(18)                                                  CPP02020
      NTYPEP=NTYPES+1                                                   CPP02030
      NTREAD=NTYPES*2                                                   CPP02040
      MAXIJ=(MAXBUF*16)/20                                              CPP02050
C                                                                       CPP02060
      ISTEP=0                                                           CPP02070
      IGSTOP=0                                                          CPP02080
      ILSTOP=0                                                          CPP02090
C                                                                       CPP02100
      DO 104 I=1,N3NXX                                                  CPP02110
      ISA(I)=(I-1)*NTRIL+1                                              CPP02120
      IHA(I)=(I-1)*NTRIL+N3TRL+1                                        CPP02130
      IFA(I)=(I-1)*NTRIL+N3TRL+N3TRL+1                                  CPP02140
      IEA(I)=(I-1)*NBASQL+N3TRL+N3TRL+N3TRL+1                           CPP02150
      IBA(I)=(I-1)*NTRIL+N3TRL+N3TRL+N3TRL+N3QRL+1                      CPP02160
      IUA(I)=(I-1)*NBASQL+N3TRL+N3TRL+N3TRL+N3QRL+N3TRL+1               CPP02170
      WRITE(3,4) I,ISA(I),IHA(I),IFA(I),IEA(I),IBA(I),IUA(I)            CPP02180
  104 CONTINUE                                                          CPP02190
      ENUC=APARA(1)                                                     CPP02200
      ESCF=APARA(2)                                                     CPP02210
      NSECMX=IUA(N3NXX)+NBASQL                                          CPP02220
      DO 105 I=1,NTYPEP                                                 CPP02230
      NSORB(I)=LPARA(I+40)                                              CPP02240
      FOCC(I)=APARA(I+10)                                               CPP02250
  105 CONTINUE                                                          CPP02260
      NTYTRI=IOFF(NTYPEP+1)                                             CPP02270
C                                                                       CPP02280
C   READ IN CI INFORMATION                                              CPP02290
      CALL RFILE(ITAP30)                                                CPP02300
      CALL WREADW(ITAP30,I30,200,101,JUNK)                              CPP02310
      NCPEX=I30(46)                                                     CPP02320
      NORPEX=I30(47)                                                    CPP02330
      LVPEX=I30(48)                                                     CPP02340
      NCTRI=IOFF(NCPEX+1)                                               CPP02350
      CALL WREADW(ITAP30,PEXC,2*NCPEX,LVPEX,LOCVEC)                     CPP02360
      CALL WREADW(ITAP30,INDPEX,20*NCPEX,LOCVEC,LOCVEC)                 CPP02370
      CALL WREADW(ITAP30,MOPEX,20,LOCVEC,LOCVEC)                        CPP02380
      CALL RCLOSE(ITAP30,3)                                             CPP02390
C                                                                       CPP02400
      WRITE(6,5) NBASIS,NTRI,NBFAO,NBATRI,NATOMS,N3N,N3NP3,N3NXX,       CPP02410
     1           IOCC,JOCC,KOCC,NTYPES,NTYPEP,NTREAD,NTYTRI             CPP02420
      WRITE(6,6) NCPEX,NORPEX,NCTRI,LVPEX,IPOL,ICONV,IPRNT              CPP02430
C                                                                       CPP02440
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::CPP02450
C                                                                       CPP02460
C   DYNAMIC ALLOCATION                                                  CPP02470
C     IC1  = EIG(NBASIS)                                                CPP02480
C     IC2  = OCC(NBASIS)                                                CPP02490
C     IC3  = EAO OR ESO(NBFAO,NBASIS)                                   CPP02500
C     IC4  = NIJ(NTRI,2)                                                CPP02510
C     IC5  = KIJ(NTRI,2)                                                CPP02520
C     IC6  = FIJ(NTYPEP,ICIUNQ)                                         CPP02530
C     IC7  = AIJ(NTYTRI,ICIUNQ)                                         CPP02540
C     IC8  = BIJ(NTYTRI,ICIUNQ)                                         CPP02550
C     IC9  = E1A(N3N)                                                   CPP02560
C     IC10 = E2A(N3N,N3N)                                               CPP02570
C     IC11 = D1W(N3N)                                                   CPP02580
C     IC12 = E1T(N3NXX)                                                 CPP02590
C     IC13 = EPS(NTRI)                                                  CPP02600
C     IC14 = ZETA(NTRI,NTYPEP)                                          CPP02610
C     IC15 = HCIA(NCTRI,N3NXX)                                          CPP02620
C     IC16 = HCI(NCTRI)                                                 CPP02630
C     IC17 = BAF2(NCPEX,N3NXX)                                          CPP02640
C     IC18 = CIA(NCPEX,N3NXX)                                           CPP02650
C     IC19 = A22(NCPEX,NCPEX)                                           CPP02660
C     IC20 = A12(NIND,NCPEX)                                            CPP02670
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::CPP02680
C                                                                       CPP02690
      IC1=1                                                             CPP02700
      IC2=IC1+NBASIS                                                    CPP02710
      IC3=IC2+NBASIS                                                    CPP02720
      IC4=IC3+NBFAO*NBASIS                                              CPP02730
      IA4=IC4+IC4-1                                                     CPP02740
      IC5=IC4+NTRI                                                      CPP02750
      IA5=IC5+IC5-1                                                     CPP02760
      IC6=IC5+NTRI                                                      CPP02770
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::CPP02780
C                                                                       CPP02790
C   READ IN AO-MO EIGENVECTORS                                          CPP02800
      CALL RMREAD(CC(IC1),18)                                           CPP02810
      CALL RMREAD(CC(IC2),19)                                           CPP02820
      CALL RMREAD(CC(IC3),22)                                           CPP02830
      IF(IPRNT.LE.2) GO TO 301                                          CPP02840
      WRITE(6,8)                                                        CPP02850
      CALL EIGOUT(CC(IC3),CC(IC1),CC(IC2),NBFAO,NBASIS,NBFAO,NBASIS,6)  CPP02860
C                                                                       CPP02870
C   CALCULATE ALPA AND BETA MATRICES                                    CPP02880
  301 CONTINUE                                                          CPP02890
      NTYTRI=IOFF(NTYPEP+1)                                             CPP02900
      DO 106 I=1,NTYTRI                                                 CPP02910
      ALPC(I)=APARA(I+40)                                               CPP02920
      BETC(I)=APARA(I+505)                                              CPP02930
  106 CONTINUE                                                          CPP02940
      CALL ZERO(ALPA,NTYPEP*NTYPEP)                                     CPP02950
      CALL ZERO(BETA,NTYPEP*NTYPEP)                                     CPP02960
      IJ=0                                                              CPP02970
      DO 107 I=1,NTYPES                                                 CPP02980
      DO 107 J=1,I                                                      CPP02990
      IJ=IJ+1                                                           CPP03000
      ALPA(I,J)=ALPC(IJ)                                                CPP03010
      BETA(I,J)=BETC(IJ)                                                CPP03020
      IF(I.EQ.J) GO TO 107                                              CPP03030
      ALPA(J,I)=ALPA(I,J)                                               CPP03040
      BETA(J,I)=BETA(I,J)                                               CPP03050
  107 CONTINUE                                                          CPP03060
C*    IF(IPRNT.LE.2) GO TO 302                                          CPP03070
      WRITE(6,9)                                                        CPP03080
      CALL MATOUT(ALPA,30,30,NTYPEP,NTYPEP,6)                           CPP03090
      WRITE(6,10)                                                       CPP03100
      CALL MATOUT(BETA,30,30,NTYPEP,NTYPEP,6)                           CPP03110
C                                                                       CPP03120
C**************************************                                 CPP03130
C***PREPARATION FOR FIRST ORDER CPHF***                                 CPP03140
C**************************************                                 CPP03150
C                                                                       CPP03160
C   FORM REGISTER MATRICES                                              CPP03170
  302 CONTINUE                                                          CPP03180
      WRITE(3,20)                                                       CPP03190
      WRITE(4,20)                                                       CPP03200
   20 FORMAT(//,2X,' NOW YOU ARE IN REGIST'/)                           CPP03210
C.................NIJ     KIJ.....                                      CPP03220
      CALL REGIST(IA(IA4),IA(IA5))                                      CPP03230
C                                                                       CPP03240
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::CPP03250
C   REST OF DYNAMIC ALLOCATION                                          CPP03260
      IC7=IC6+NTYPEP*ICIUNQ                                             CPP03270
      IC8=IC7+NTYTRI*ICIUNQ                                             CPP03280
      IC9=IC8+NTYTRI*ICIUNQ                                             CPP03290
      IC10=IC9+N3N                                                      CPP03300
      IC11=IC10+N3N*N3N                                                 CPP03310
      IC12=IC11+N3N                                                     CPP03320
      IC13=IC12+N3NXX                                                   CPP03330
      IC14=IC13+NTRI                                                    CPP03340
      IC15=IC14+NTRI*NTYPEP                                             CPP03350
      IC16=IC15+NCTRI*N3NXX                                             CPP03360
      IC17=IC16+NCTRI                                                   CPP03370
      IC18=IC17+NCPEX*N3NXX                                             CPP03380
      IC19=IC18+NCPEX*N3NXX                                             CPP03390
      IC20=IC19+NCPEX*NCPEX                                             CPP03400
      ICKEEP=IC20+NIND*NCPEX                                            CPP03410
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::CPP03420
C                                                                       CPP03430
      DO 108 I=1,IIUNQ+IJGRP                                            CPP03440
      IEIJ(I)=(I-1)*NBASQL+1                                            CPP03450
      IF(IPRNT.LE.2) GO TO 108                                          CPP03460
      WRITE(3,4) I,IEIJ(I)                                              CPP03470
  108 CONTINUE                                                          CPP03480
C                                                                       CPP03490
C   CALCULATE REQUIRED CORE MEMORY FOR THE CALCULATION                  CPP03500
      CALL CORSIZ(MAXCOR,N3NXX,NADD,MAXVC,IPOL)                         CPP03510
C                                                                       CPP03520
C   FORM THE COUPLING CONSTANT MATRICES                                 CPP03530
      WRITE(3,21)                                                       CPP03540
      WRITE(4,21)                                                       CPP03550
   21 FORMAT(//,2X,' NOW YOU ARE IN BARES'/)                            CPP03560
C................FIJ     AIJ     BIJ.....                               CPP03570
      CALL BARES(CC(IC6),CC(IC7),CC(IC8))                               CPP03580
C                                                                       CPP03590
C#########################                                              CPP03600
C###BEGINNING OF STEP 1###                                              CPP03610
C#########################                                              CPP03620
C                                                                       CPP03630
C   READ IN FIRST AND SECOND DERIVATIVE INFORMATION                     CPP03640
C   AND FORM SM, HM AND FM MATRICES                                     CPP03650
C**   IF(ISTEP.GE.1) GO TO 303                                          CPP03660
      WRITE(3,22)                                                       CPP03670
      WRITE(4,22)                                                       CPP03680
   22 FORMAT(//,2X,' NOW YOU ARE IN DERMAT'/)                           CPP03690
      IC21=ICKEEP                                                       CPP03700
      IC22=IC21+NTRI                                                    CPP03710
      IC23=IC22+NBFAO*NBFAO                                             CPP03720
      ICMAX=IC23+NBFAO*NBFAO                                            CPP03730
      WRITE(3,12) ICMAX,MAXCOR                                          CPP03740
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP03750
C.................EAO     E1A     E2A      SS       U        T........  CPP03760
      CALL DERMAT(CC(IC3),CC(IC9),CC(IC10),CC(IC21),CC(IC22),CC(IC23),  CPP03770
C..............................                                         CPP03780
     1            NSECMX,ISTEP)                                         CPP03790
C                                                                       CPP03800
      LPARA(101)=1                                                      CPP03810
      CALL RMWRIT(LPARA,2)                                              CPP03820
C                                                                       CPP03830
C###################                                                    CPP03840
C###END OF STEP 1###                                                    CPP03850
C###################                                                    CPP03860
C                                                                       CPP03870
C   READ IN SO-MO EIGENVECTORS                                          CPP03880
  303 CONTINUE                                                          CPP03890
      CALL RMREAD(CC(IC3),21)                                           CPP03900
      IF(IPRNT.LE.2) GO TO 304                                          CPP03910
      WRITE(6,11)                                                       CPP03920
      CALL EIGOUT(CC(IC3),CC(IC1),CC(IC2),NBFAO,NBASIS,NBFAO,NBASIS,6)  CPP03930
C                                                                       CPP03940
C#########################                                              CPP03950
C###BEGINNING OF STEP 2###                                              CPP03960
C#########################                                              CPP03970
C                                                                       CPP03980
  304 CONTINUE                                                          CPP03990
      IF(ISTEP.GE.2) GO TO 305                                          CPP04000
C                                                                       CPP04010
C   CALCULATE SCF ENERGY FOR A TEST                                     CPP04020
      WRITE(3,23)                                                       CPP04030
      WRITE(4,23)                                                       CPP04040
   23 FORMAT(//,2X,' NOW YOU ARE IN E0CAL'/)                            CPP04050
      IC21=ICKEEP                                                       CPP04060
      IC22=IC21+NTRI                                                    CPP04070
      IC23=IC22+NTRI                                                    CPP04080
      IC24=IC23+NTRI                                                    CPP04090
      IC25=IC24+NBFAO*NBFAO                                             CPP04100
      IC26=IC25+NBFAO*NBFAO                                             CPP04110
      IC27=IC26+NTRI*NTYPES                                             CPP04120
      IC28=IC27+NTRI*NTYPES                                             CPP04130
      IC29=IC28+NTRI*NTYPES                                             CPP04140
      IA29=IC29+IC29-1                                                  CPP04150
      ICMAX=IC29+MAXBF2                                                 CPP04160
      WRITE(3,12) ICMAX,MAXCOR                                          CPP04170
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP04180
C................ESO     HSO      HMO      TM       U        T........  CPP04190
      CALL E0CAL(CC(IC3),CC(IC21),CC(IC22),CC(IC23),CC(IC24),CC(IC25),  CPP04200
C................DP       DQ       DM       MIJ      PQ.............    CPP04210
     1           CC(IC26),CC(IC27),CC(IC28),IA(IA29),IA(IA29+MAXIJ))    CPP04220
C                                                                       CPP04230
      LPARA(101)=2                                                      CPP04240
      CALL RMWRIT(LPARA,2)                                              CPP04250
C                                                                       CPP04260
C###################                                                    CPP04270
C###END OF STEP 2###                                                    CPP04280
C###################                                                    CPP04290
C                                                                       CPP04300
C#########################                                              CPP04310
C###BEGINNING OF STEP 3###                                              CPP04320
C#########################                                              CPP04330
C                                                                       CPP04340
C   CALCULATE SCF FIRST ENERGY DERIVATIVES                              CPP04350
  305 CONTINUE                                                          CPP04360
      IF(ISTEP.GE.3) GO TO 306                                          CPP04370
      WRITE(3,24)                                                       CPP04380
      WRITE(4,24)                                                       CPP04390
   24 FORMAT(//,2X,' NOW YOU ARE IN EADER'/)                            CPP04400
      IC21=ICKEEP                                                       CPP04410
      IC22=IC21+N3N                                                     CPP04420
      IC23=IC22+N3N                                                     CPP04430
      IC24=IC23+NTRI                                                    CPP04440
      IC25=IC24+NTRI                                                    CPP04450
      ICMAX=IC25+NTRI*N3N*NTREAD                                        CPP04460
      WRITE(3,12) ICMAX,MAXCOR                                          CPP04470
C................D1W      E1T      EPS      D1O      D1T......          CPP04480
      CALL EADER(CC(IC11),CC(IC12),CC(IC13),CC(IC21),CC(IC22),          CPP04490
C................SM       HM       TM..................                 CPP04500
     1           CC(IC23),CC(IC24),CC(IC25),N3N*NTREAD)                 CPP04510
C                                                                       CPP04520
      LPARA(101)=3                                                      CPP04530
      CALL RMWRIT(LPARA,2)                                              CPP04540
C                                                                       CPP04550
C###################                                                    CPP04560
C###END OF STEP 3###                                                    CPP04570
C###################                                                    CPP04580
C                                                                       CPP04590
C#########################                                              CPP04600
C###BEGINNING OF STEP 4###                                              CPP04610
C#########################                                              CPP04620
C                                                                       CPP04630
C   FORM HA MATRICES                                                    CPP04640
  306 CONTINUE                                                          CPP04650
      IF(ISTEP.GE.4) GO TO 307                                          CPP04660
      WRITE(3,25)                                                       CPP04670
      WRITE(4,25)                                                       CPP04680
   25 FORMAT(//,2X,' NOW YOU ARE IN HAMAT'/)                            CPP04690
      IC21=ICKEEP                                                       CPP04700
      IC22=IC21+N3N                                                     CPP04710
      IC23=IC22+NTRI*N3N                                                CPP04720
      ICMAX=IC23+NTRI*N3N*NTREAD                                        CPP04730
      WRITE(3,12) ICMAX,MAXCOR                                          CPP04740
C................FIJ     AIJ     BIJ     D1W      E1T      HCIA.....    CPP04750
      CALL HAMAT(CC(IC6),CC(IC7),CC(IC8),CC(IC11),CC(IC12),CC(IC15),    CPP04760
C................D1E      HM       TM........................           CPP04770
     1           CC(IC21),CC(IC22),CC(IC23),N3N*NTREAD,N3NXX)           CPP04780
C                                                                       CPP04790
      LPARA(101)=4                                                      CPP04800
      CALL RMWRIT(LPARA,2)                                              CPP04810
C                                                                       CPP04820
C###################                                                    CPP04830
C###END OF STEP 4###                                                    CPP04840
C###################                                                    CPP04850
C                                                                       CPP04860
C#########################                                              CPP04870
C###BEGINNING OF STEP 5###                                              CPP04880
C#########################                                              CPP04890
C                                                                       CPP04900
C   FORM DERIVATIVE LAGRANGIAN MATRICES                                 CPP04910
  307 CONTINUE                                                          CPP04920
      IF(ISTEP.GE.5) GO TO 308                                          CPP04930
      WRITE(3,26)                                                       CPP04940
      WRITE(4,26)                                                       CPP04950
   26 FORMAT(//,2X,' NOW YOU ARE IN FAMAT'/)                            CPP04960
      IC21=ICKEEP                                                       CPP04970
      IC22=IC21+NTRI                                                    CPP04980
      IC23=IC22+NBASIS*NBASIS                                           CPP04990
      IC24=IC23+NTRI*N3N                                                CPP05000
      ICMAX=IC24+NTRI*N3N*NTREAD                                        CPP05010
      WRITE(3,12) ICMAX,MAXCOR                                          CPP05020
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP05030
C................FF       EPA      HM       TM..................        CPP05040
      CALL FAMAT(CC(IC21),CC(IC22),CC(IC23),CC(IC24),N3N*NTREAD)        CPP05050
C                                                                       CPP05060
      LPARA(101)=5                                                      CPP05070
      CALL RMWRIT(LPARA,2)                                              CPP05080
C                                                                       CPP05090
C###################                                                    CPP05100
C###END OF STEP 5###                                                    CPP05110
C###################                                                    CPP05120
C                                                                       CPP05130
C#########################                                              CPP05140
C###BEGINNING OF STEP 6###                                              CPP05150
C#########################                                              CPP05160
C                                                                       CPP05170
C   CALCULATE BA MATRICES                                               CPP05180
C   FOR INDEPENDENT PAIRS                                               CPP05190
  308 CONTINUE                                                          CPP05200
      IF(ISTEP.GE.6) GO TO 309                                          CPP05210
      NGROUP=(N3N-1)/NADD(1)+1                                          CPP05220
      WRITE(3,27) NGROUP                                                CPP05230
      WRITE(4,27) NGROUP                                                CPP05240
   27 FORMAT(//,2X,' NOW YOU ARE IN BAIND (# OF GROUPS = ',I3,')'/)     CPP05250
      IC21=ICKEEP                                                       CPP05260
      IC22=IC21+NBFAO*NBFAO                                             CPP05270
      IC23=IC22+NBFAO*NBFAO                                             CPP05280
      IC24=IC23+NTRI                                                    CPP05290
      IC25=IC24+NBASIS*NBASIS                                           CPP05300
      IC26=IC25+NTRI*NTYPEP*NADD(1)                                     CPP05310
      IC27=IC26+NTRI*NTYPEP*NADD(1)                                     CPP05320
      IC28=IC27+NTRI*NTYPEP*NADD(1)                                     CPP05330
      IA28=IC28+IC28-1                                                  CPP05340
      ICMAX=IC28+MAXBF2                                                 CPP05350
      WRITE(3,12) ICMAX,MAXCOR                                          CPP05360
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP05370
C................ESO     NIJ     ZETA     U        T        B0.......   CPP05380
      CALL BAIND(CC(IC3),IA(IA4),CC(IC14),CC(IC21),CC(IC22),CC(IC23),   CPP05390
C................EPA      DP       DQ       DM       SM.......          CPP05400
     1           CC(IC24),CC(IC25),CC(IC26),CC(IC27),CC(IC27),          CPP05410
C................MIJ      PQ.....................                       CPP05420
     2           IA(IA28),IA(IA28+MAXIJ),NADD(1))                       CPP05430
C                                                                       CPP05440
      LPARA(101)=6                                                      CPP05450
      CALL RMWRIT(LPARA,2)                                              CPP05460
C*    IF(IPRNT.LE.2) GO TO 309                                          CPP05470
C*    CALL WRBMAT(CC(IC21),N3N)                                         CPP05480
C                                                                       CPP05490
C###################                                                    CPP05500
C###END OF STEP 6###                                                    CPP05510
C###################                                                    CPP05520
C                                                                       CPP05530
C#########################                                              CPP05540
C###BEGINNING OF STEP 7###                                              CPP05550
C#########################                                              CPP05560
C                                                                       CPP05570
C   FOR DEPENDENT PAIRS                                                 CPP05580
  309 CONTINUE                                                          CPP05590
      IF(ISTEP.GE.7) GO TO 310                                          CPP05600
      NGROUP=(N3N-1)/NADD(2)+1                                          CPP05610
      WRITE(3,28) NGROUP                                                CPP05620
      WRITE(4,28) NGROUP                                                CPP05630
   28 FORMAT(//,2X,' NOW YOU ARE IN BADEP (# OF GROUPS = ',I3,')'/)     CPP05640
      IC21=ICKEEP                                                       CPP05650
      IC22=IC21+NBFAO*NBFAO                                             CPP05660
      IC23=IC22+NBFAO*NBFAO                                             CPP05670
      IC24=IC23+NTRI                                                    CPP05680
      IC25=IC24+NTRI                                                    CPP05690
      IC26=IC25+NTRI*NADD(2)                                            CPP05700
      IC27=IC26+NTRI*NADD(2)                                            CPP05710
      IC28=IC27+NTRI*NADD(2)                                            CPP05720
      IA28=IC28+IC28-1                                                  CPP05730
      ICMAX=IC28+MAXBF2                                                 CPP05740
      WRITE(3,12) ICMAX,MAXCOR                                          CPP05750
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP05760
C................EIG     ESO     KIJ     U        T.......              CPP05770
      CALL BADEP(CC(IC1),CC(IC3),IA(IA5),CC(IC21),CC(IC22),             CPP05780
C................B0       FF       DP       DQ       DM.......          CPP05790
     1           CC(IC23),CC(IC24),CC(IC25),CC(IC26),CC(IC27),          CPP05800
C................MIJ      PQ.....................                       CPP05810
     2           IA(IA28),IA(IA28+MAXIJ),NADD(2))                       CPP05820
C                                                                       CPP05830
      LPARA(101)=7                                                      CPP05840
      CALL RMWRIT(LPARA,2)                                              CPP05850
C*    IF(IPRNT.LE.2) GO TO 310                                          CPP05860
C*    CALL WRBMAT(CC(IC21),N3N)                                         CPP05870
C                                                                       CPP05880
C###################                                                    CPP05890
C###END OF STEP 7###                                                    CPP05900
C###################                                                    CPP05910
C                                                                       CPP05920
C#########################                                              CPP05930
C###BEGINNING OF STEP 8###                                              CPP05940
C#########################                                              CPP05950
C                                                                       CPP05960
C   CALCULATE H AND E MATRICES                                          CPP05970
  310 CONTINUE                                                          CPP05980
      IF(ISTEP.GE.8) GO TO 311                                          CPP05990
      WRITE(3,29)                                                       CPP06000
      WRITE(4,29)                                                       CPP06010
   29 FORMAT(//,2X,' NOW YOU ARE IN HEMAT'/)                            CPP06020
      IC21=ICKEEP                                                       CPP06030
      IC22=IC21+NTRI                                                    CPP06040
      IC23=IC22+NTRI                                                    CPP06050
      IC24=IC23+NTRI*NTYPEP                                             CPP06060
      IC25=IC24+NTRI*NTYPEP                                             CPP06070
      IC26=IC25+NTRI*NTYPEP                                             CPP06080
      IC27=IC26+NBFAO*NBFAO                                             CPP06090
      IC28=IC27+NTRI                                                    CPP06100
      IC29=IC28+NBFAO*NBFAO                                             CPP06110
      IC30=IC29+NBFAO*NBFAO                                             CPP06120
      IC31=IC30+NBFAO*NBFAO                                             CPP06130
      IA31=IC31+IC31-1                                                  CPP06140
      ICMAX=IC31+MAXBF2                                                 CPP06150
      WRITE(3,12) ICMAX,MAXCOR                                          CPP06160
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP06170
C................ESO     FIJ     AIJ     BIJ     HCI......              CPP06180
      CALL HEMAT(CC(IC3),CC(IC6),CC(IC7),CC(IC8),CC(IC16),              CPP06190
C................HSO      HMO      DP       DQ       DM       EIJ...... CPP06200
     1           CC(IC21),CC(IC22),CC(IC23),CC(IC24),CC(IC25),CC(IC26), CPP06210
C................ZIJ      EE       U        T........                   CPP06220
     2           CC(IC27),CC(IC28),CC(IC29),CC(IC30),                   CPP06230
C................MIJ      PQ.............                               CPP06240
     3           IA(IA31),IA(IA31+MAXIJ))                               CPP06250
C                                                                       CPP06260
      LPARA(101)=8                                                      CPP06270
      CALL RMWRIT(LPARA,2)                                              CPP06280
C*    IF(IPRNT.LE.2) GO TO 311                                          CPP06290
C*    CALL WRBMAT(CC(IC21),N3N)                                         CPP06300
C                                                                       CPP06310
C###################                                                    CPP06320
C###END OF STEP 8###                                                    CPP06330
C###################                                                    CPP06340
C                                                                       CPP06350
C#########################                                              CPP06360
C###BEGINNING OF STEP 9###                                              CPP06370
C#########################                                              CPP06380
C                                                                       CPP06390
C   FORM BF MATRICES                                                    CPP06400
  311 CONTINUE                                                          CPP06410
      IF(IPOL.NE.0) GO TO 312                                           CPP06420
C**   IF(ISTEP.GE.9) GO TO 312                                          CPP06430
      WRITE(3,30)                                                       CPP06440
      WRITE(4,30)                                                       CPP06450
   30 FORMAT(//,2X,' NOW YOU ARE IN BFMAT'/)                            CPP06460
      IC21=ICKEEP                                                       CPP06470
      IC22=IC21+NTRI*3                                                  CPP06480
      IC23=IC22+NTRI*NTYPES*3                                           CPP06490
      IC24=IC23+NBASIS*NBASIS                                           CPP06500
      ICMAX=IC24+NTRI                                                   CPP06510
      WRITE(3,12) ICMAX,MAXCOR                                          CPP06520
C................NIJ     KIJ     FIJ     E1T      HCIA     HF.......    CPP06530
      CALL BFMAT(IA(IA4),IA(IA5),CC(IC6),CC(IC12),CC(IC15),CC(IC21),    CPP06540
C................HM       EPF      BF...................                CPP06550
     1           CC(IC22),CC(IC23),CC(IC24),N3NXX,ISTEP)                CPP06560
C                                                                       CPP06570
      LPARA(101)=9                                                      CPP06580
      CALL RMWRIT(LPARA,2)                                              CPP06590
      IF(IPRNT.LE.2) GO TO 312                                          CPP06600
      WRITE(6,*) ' AFTER BFMAT'                                         CPP06610
      CALL WRBMAT(CC(IC21),N3NXX)                                       CPP06620
C                                                                       CPP06630
C###################                                                    CPP06640
C###END OF STEP 9###                                                    CPP06650
C###################                                                    CPP06660
C                                                                       CPP06670
C##########################                                             CPP06680
C###BEGINNING OF STEP 10###                                             CPP06690
C##########################                                             CPP06700
C                                                                       CPP06710
C   FORM BAF MATRICES                                                   CPP06720
  312 CONTINUE                                                          CPP06730
      IF(IGSTEP.GE.10) GO TO 313                                        CPP06740
      WRITE(3,31)                                                       CPP06750
      WRITE(4,31)                                                       CPP06760
   31 FORMAT(//,2X,' NOW YOU ARE IN AB2CI'/)                            CPP06770
      IC21=ICKEEP                                                       CPP06780
      IC22=IC21+NTRI*N3N                                                CPP06790
      ICMAX=IC22+NBASIS*NBASIS                                          CPP06800
      WRITE(3,12) ICMAX,MAXCOR                                          CPP06810
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP06820
C................NIJ     E1T      HCIA     HCI      BAF2     A22......  CPP06830
      CALL AB2CI(IA(IA4),CC(IC12),CC(IC15),CC(IC16),CC(IC17),CC(IC19),  CPP06840
C................A12      SA       EIJ............                      CPP06850
     1           CC(IC20),CC(IC21),CC(IC22),N3NXX)                      CPP06860
C                                                                       CPP06870
      LPARA(101)=10                                                     CPP06880
      CALL RMWRIT(LPARA,2)                                              CPP06890
C*    IF(IPRNT.LE.2) GO TO 313                                          CPP06900
C*    CALL WRBMAT(CC(IC21),N3N)                                         CPP06910
C                                                                       CPP06920
C####################                                                   CPP06930
C###END OF STEP 10###                                                   CPP06940
C####################                                                   CPP06950
C                                                                       CPP06960
C********************************                                       CPP06970
C***FIRST ORDER CPHF PROCEDURE***                                       CPP06980
C********************************                                       CPP06990
C                                                                       CPP07000
C##########################                                             CPP07010
C###BEGINNING OF STEP 11###                                             CPP07020
C##########################                                             CPP07030
C                                                                       CPP07040
C   SOLVE SIMULTANEOUS EQUATIONS ITERATIVELY                            CPP07050
  313 CONTINUE                                                          CPP07060
      IF(ISTEP.GE.11) GO TO 314                                         CPP07070
      NGROUP=(N3NXX-1)/NADD(3)+1                                        CPP07080
      WRITE(3,32) NGROUP                                                CPP07090
      WRITE(4,32) NGROUP                                                CPP07100
   32 FORMAT(//,2X,' NOW YOU ARE IN CPHF (# OF GROUPS = ',I3,')'/)      CPP07110
      NINDPX=NIND+NCPEX                                                 CPP07120
      IC21=ICKEEP                                                       CPP07130
      IC22=IC21+NBFAO*NBFAO                                             CPP07140
      IC23=IC22+NBFAO*NBFAO                                             CPP07150
      IC24=IC23+NBFAO*NBFAO                                             CPP07160
      IC25=IC24+NTRI                                                    CPP07170
      IC26=IC25+NINDPX                                                  CPP07180
      IC27=IC26+NINDPX                                                  CPP07190
      IC28=IC27+NTYPEP*(NTRI*NADD(3))                                   CPP07200
      IC29=IC28+NTYPEP*(NTRI*NADD(3))                                   CPP07210
      IC30=IC29+NTYPEP*(NTRI*NADD(3))                                   CPP07220
      IC31=IC30+MAXVC*(MAXVC+NADD(3))                                   CPP07230
      IC32=IC31+NINDPX*MAXVC*2                                          CPP07240
      IA32=IC32+IC32-1                                                  CPP07250
      ICMAX=IC32+MAXBF2                                                 CPP07260
      WRITE(3,12) ICMAX,MAXCOR                                          CPP07270
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP07280
C...............ESO     NIJ     KIJ     ZETA     BAF2     CIA......     CPP07290
      CALL CPHF(CC(IC3),IA(IA4),IA(IA5),CC(IC14),CC(IC17),CC(IC18),     CPP07300
C...............A22      A12      A........                             CPP07310
     1          CC(IC19),CC(IC20),CC(IC21),                             CPP07320
C...............U        T        B0       AA       TT       DP.......  CPP07330
     2          CC(IC22),CC(IC23),CC(IC24),CC(IC25),CC(IC26),CC(IC27),  CPP07340
C...............DQ       DM       W        BB       MIJ......           CPP07350
     3          CC(IC28),CC(IC29),CC(IC30),CC(IC31),IA(IA32),           CPP07360
C...............PQ........................................              CPP07370
     4          IA(IA32+MAXIJ),NINDPX,MAXVC,NADD(3),N3NXX)              CPP07380
C                                                                       CPP07390
      LPARA(101)=11                                                     CPP07400
      CALL RMWRIT(LPARA,2)                                              CPP07410
      IF(IPRNT.LE.2) GO TO 314                                          CPP07420
      CALL WRBMAT(CC(IC21),N3N)                                         CPP07430
C                                                                       CPP07440
C####################                                                   CPP07450
C###END OF STEP 11###                                                   CPP07460
C####################                                                   CPP07470
C                                                                       CPP07480
C##########################                                             CPP07490
C###BEGINNING OF STEP 12###                                             CPP07500
C##########################                                             CPP07510
C                                                                       CPP07520
C   CALCULATE DEPENDENT PAIRS OF U MATRICES                             CPP07530
  314 CONTINUE                                                          CPP07540
      IF(ISTEP.GE.12) GO TO 315                                         CPP07550
      NGROUP=(N3NXX-1)/NADD(4)+1                                        CPP07560
      WRITE(3,33) NGROUP                                                CPP07570
      WRITE(4,33) NGROUP                                                CPP07580
   33 FORMAT(//,2X,' NOW YOU ARE IN UCMAT (# OF GROUPS = ',I3,')'/)     CPP07590
      IC21=ICKEEP                                                       CPP07600
      IC22=IC21+NBFAO*NBFAO                                             CPP07610
      IC23=IC22+NBFAO*NBFAO                                             CPP07620
      IC24=IC23+NTRI                                                    CPP07630
      IC25=IC24+NTRI*NADD(4)                                            CPP07640
      IC26=IC25+NTRI*NADD(4)                                            CPP07650
      IC27=IC26+NTRI*NADD(4)                                            CPP07660
      IC28=IC27+NTRI                                                    CPP07670
      IA28=IC28+IC28-1                                                  CPP07680
      ICMAX=IC28+MAXBF2                                                 CPP07690
      WRITE(3,12) ICMAX,MAXCOR                                          CPP07700
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP07710
C................EIG     ESO     NIJ     KIJ     U........              CPP07720
      CALL UCMAT(CC(IC1),CC(IC3),IA(IA4),IA(IA5),CC(IC21),              CPP07730
C................T        B        DP       DQ       DM       DEL...... CPP07740
     1           CC(IC22),CC(IC23),CC(IC24),CC(IC25),CC(IC26),CC(IC27), CPP07750
C................MIJ      PQ...........................                 CPP07760
     2           IA(IA28),IA(IA28+MAXIJ),NADD(4),N3NXX)                 CPP07770
C                                                                       CPP07780
      LPARA(101)=12                                                     CPP07790
      CALL RMWRIT(LPARA,2)                                              CPP07800
      IF(IPRNT.LE.2) GO TO 315                                          CPP07810
      CALL WRBMAT(CC(IC21),N3N)                                         CPP07820
C                                                                       CPP07830
C####################                                                   CPP07840
C###END OF STEP 12###                                                   CPP07850
C####################                                                   CPP07860
C                                                                       CPP07870
C##########################                                             CPP07880
C###BEGINNING OF STEP 13###                                             CPP07890
C##########################                                             CPP07900
C                                                                       CPP07910
C   COMPLETE U MATRICES                                                 CPP07920
  315 CONTINUE                                                          CPP07930
      IF(ISTEP.GE.13) GO TO 316                                         CPP07940
      WRITE(3,34)                                                       CPP07950
      WRITE(4,34)                                                       CPP07960
   34 FORMAT(//,2X,' NOW YOU ARE IN UXMAT'/)                            CPP07970
      IC21=ICKEEP                                                       CPP07980
      IC22=IC21+NTRI*N3NXX                                              CPP07990
      IC23=IC22+NTRI*N3NXX                                              CPP08000
      ICMAX=IC23+NBASIS*NBASIS                                          CPP08010
      WRITE(3,12) ICMAX,MAXCOR                                          CPP08020
C................KIJ     B        SM       U..............              CPP08030
      CALL UXMAT(IA(IA5),CC(IC21),CC(IC22),CC(IC23),N3NXX)              CPP08040
C                                                                       CPP08050
      LPARA(101)=13                                                     CPP08060
      CALL RMWRIT(LPARA,2)                                              CPP08070
      IF(IPRNT.LE.2) GO TO 316                                          CPP08080
      CALL WRBMAT(CC(IC21),N3N)                                         CPP08090
C                                                                       CPP08100
C####################                                                   CPP08110
C###END OF STEP 13###                                                   CPP08120
C####################                                                   CPP08130
C                                                                       CPP08140
C##########################                                             CPP08150
C###BEGINNING OF STEP 14###                                             CPP08160
C##########################                                             CPP08170
C                                                                       CPP08180
C   CALCULATE WA MATRICES                                               CPP08190
  316 CONTINUE                                                          CPP08200
      IF(ISTEP.GE.14) GO TO 317                                         CPP08210
      NGROUP=(N3N-1)/NADD(5)+1                                          CPP08220
      WRITE(3,35) NGROUP                                                CPP08230
      WRITE(4,35) NGROUP                                                CPP08240
   35 FORMAT(//,2X,' NOW YOU ARE IN WAMAT (# OF GROUPS = ',I3,')'/)     CPP08250
      IC21=ICKEEP                                                       CPP08260
      IC22=IC21+NBASIS*NBASIS                                           CPP08270
      IC23=IC22+NBASIS*NBASIS                                           CPP08280
      IC24=IC23+NBASIS*NBASIS                                           CPP08290
      IC25=IC24+NBASIS*NBASIS                                           CPP08300
      IC26=IC25+NTRI*NTYPEP*NADD(5)                                     CPP08310
      IC27=IC26+NTRI*NTYPEP*NADD(5)                                     CPP08320
      IC28=IC27+NTRI*NTYPEP*NADD(5)                                     CPP08330
      IA28=IC28+IC28-1                                                  CPP08340
      ICMAX=IC28+MAXBF2                                                 CPP08350
      WRITE(3,12) ICMAX,MAXCOR                                          CPP08360
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP08370
C................ESO     ZETA     U        T        EPA      WA.......  CPP08380
      CALL WAMAT(CC(IC3),CC(IC14),CC(IC21),CC(IC22),CC(IC23),CC(IC24),  CPP08390
C................SA       DP       DQ       DM       MIJ......          CPP08400
     1           CC(IC25),CC(IC26),CC(IC27),CC(IC25),IA(IA28),          CPP08410
C................PQ.....................                                CPP08420
     2           IA(IA28+MAXIJ),NADD(5))                                CPP08430
C                                                                       CPP08440
      LPARA(101)=14                                                     CPP08450
      CALL RMWRIT(LPARA,2)                                              CPP08460
C*    IF(IPRNT.LE.2) GO TO 317                                          CPP08470
C*    CALL WRBMAT(CC(IC21),N3N)                                         CPP08480
C                                                                       CPP08490
C####################                                                   CPP08500
C###END OF STEP 14###                                                   CPP08510
C####################                                                   CPP08520
C                                                                       CPP08530
C##########################                                             CPP08540
C###BEGINNING OF STEP 15###                                             CPP08550
C##########################                                             CPP08560
C                                                                       CPP08570
C   CALCULATE SCF SECOND DERIVATIVES                                    CPP08580
  317 CONTINUE                                                          CPP08590
      IF(ISTEP.GE.15) GO TO 318                                         CPP08600
      WRITE(3,36)                                                       CPP08610
      WRITE(4,36)                                                       CPP08620
   36 FORMAT(//,2X,' NOW YOU ARE IN SCF2W'/)                            CPP08630
      IC21=ICKEEP                                                       CPP08640
      IC22=IC21+N3N*N3N                                                 CPP08650
      IC23=IC22+N3N*N3N                                                 CPP08660
      IC24=IC23+N3N*N3N                                                 CPP08670
      IC25=IC24+NBASIS*NBASIS                                           CPP08680
      IC26=IC25+NBASIS*NBASIS                                           CPP08690
      IC27=IC26+NTRI*N3N                                                CPP08700
      ICMAX=IC27+NBASIS*NBASIS*N3N                                      CPP08710
      WRITE(3,12) ICMAX,MAXCOR                                          CPP08720
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP08730
C................E2M      E2P      E2Q      EPA      UA       SA....... CPP08740
      CALL SCF2W(CC(IC21),CC(IC22),CC(IC23),CC(IC24),CC(IC25),CC(IC26), CPP08750
C................WA.......                                              CPP08760
     1           CC(IC27))                                              CPP08770
C                                                                       CPP08780
      WRITE(3,37)                                                       CPP08790
      WRITE(4,37)                                                       CPP08800
   37 FORMAT(//,2X,' NOW YOU ARE IN SCF2ND'/)                           CPP08810
      IC21=ICKEEP                                                       CPP08820
      IC22=IC21+N3N*N3N                                                 CPP08830
      IC23=IC22+N3N*N3N                                                 CPP08840
      IC24=IC23+N3N*N3N                                                 CPP08850
      IC25=IC24+NTRI*N3N                                                CPP08860
      ICMAX=IC25+NBASIS*NBASIS                                          CPP08870
      WRITE(3,12) ICMAX,MAXCOR                                          CPP08880
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     CPP08890
C.................E2A      HCIA     CIA      E2M      E2C      E2T      CPP08900
      CALL SCF2ND(CC(IC10),CC(IC15),CC(IC18),CC(IC21),CC(IC22),CC(IC23),CPP08910
C.................SA       EIJ                                          CPP08920
     1            CC(IC24),CC(IC25))                                    CPP08930
C                                                                       CPP08940
      LPARA(101)=15                                                     CPP08950
      CALL RMWRIT(LPARA,2)                                              CPP08960
C*    IF(IPRNT.LE.2) GO TO 318                                          CPP08970
C*    CALL WRBMAT(CC(IC21),N3N)                                         CPP08980
C                                                                       CPP08990
C####################                                                   CPP09000
C###END OF STEP 15###                                                   CPP09010
C####################                                                   CPP09020
C                                                                       CPP09030
C##########################                                             CPP09040
C###BEGINNING OF STEP 16###                                             CPP09050
C##########################                                             CPP09060
C                                                                       CPP09070
C   CALCULATE DIPOLE MOMENT DERIVATIVES                                 CPP09080
  318 CONTINUE                                                          CPP09090
      IF(IPOL.NE.0) GO TO 325                                           CPP09100
      IF(ISTEP.GE.16) GO TO 319                                         CPP09110
      WRITE(3,38)                                                       CPP09120
      WRITE(4,38)                                                       CPP09130
   38 FORMAT(//,2X,' NOW YOU ARE IN DIPMO'/)                            CPP09140
      IC21=ICKEEP                                                       CPP09150
      IC22=IC21+NTRI*3                                                  CPP09160
      ICMAX=IC22+NBASIS*NBASIS*N3N                                      CPP09170
      WRITE(3,12) ICMAX,MAXCOR                                          CPP09180
C................HCIA     CIA      HF       UA.......                   CPP09190
      CALL DIPMO(CC(IC15),CC(IC18),CC(IC21),CC(IC22))                   CPP09200
C                                                                       CPP09210
      LPARA(101)=16                                                     CPP09220
      CALL RMWRIT(LPARA,2)                                              CPP09230
C*    IF(IPRNT.LE.2) GO TO 319                                          CPP09240
C*    CALL WRBMAT(CC(IC21),N3N)                                         CPP09250
C                                                                       CPP09260
C####################                                                   CPP09270
C###END OF STEP 16###                                                   CPP09280
C####################                                                   CPP09290
C                                                                       CPP09300
C##########################                                             CPP09310
C###BEGINNING OF STEP 17###                                             CPP09320
C##########################                                             CPP09330
C                                                                       CPP09340
  319 CONTINUE                                                          CPP09350
      IF(ISTEP.GE.17) GO TO 325                                         CPP09360
      WRITE(3,39)                                                       CPP09370
      WRITE(4,39)                                                       CPP09380
   39 FORMAT(//,2X,' NOW YOU ARE IN POLAR'/)                            CPP09390
      IC21=ICKEEP                                                       CPP09400
      IC22=IC21+NTRI*3                                                  CPP09410
      ICMAX=IC22+NBASIS*NBASIS*3                                        CPP09420
      WRITE(3,12) ICMAX,MAXCOR                                          CPP09430
C................HCIA     CIA      HF       UF.......                   CPP09440
      CALL POLAR(CC(IC15),CC(IC18),CC(IC21),CC(IC22))                   CPP09450
C                                                                       CPP09460
      LPARA(101)=17                                                     CPP09470
      CALL RMWRIT(LPARA,2)                                              CPP09480
C*    IF(IPRNT.LE.2) GO TO 325                                          CPP09490
C*    CALL WRBMAT(CC(IC21),N3N)                                         CPP09500
C                                                                       CPP09510
C####################                                                   CPP09520
C###END OF STEP 17###                                                   CPP09530
C####################                                                   CPP09540
C                                                                       CPP09550
C   STORE CONTRIBUTION TO THIRD DERIVATIVES                             CPP09560
  325 CONTINUE                                                          CPP09570
      WRITE(3,40)                                                       CPP09580
      WRITE(4,40)                                                       CPP09590
   40 FORMAT(//,2X,' NOW YOU ARE IN THIRD'/)                            CPP09600
C................E1T      HCIA     CIA.............                     CPP09610
      CALL THIRD(CC(IC12),CC(IC15),CC(IC18),NSECMX)                     CPP09620
C                                                                       CPP09630
      GO TO 400                                                         CPP09640
C                                                                       CPP09650
  399 WRITE(6,13) ICMAX,MAXCOR                                          CPP09660
      WRITE(3,13) ICMAX,MAXCOR                                          CPP09670
  400 CONTINUE                                                          CPP09680
      IF(ISTEP.GE.17) THEN                                              CPP09690
        WRITE(6,14)                                                     CPP09700
        WRITE(4,14)                                                     CPP09710
      END IF                                                            CPP09720
      DO 109 I=1,2                                                      CPP09730
  109 WRITE(4,*) ' '                                                    CPP09740
      WRITE(4,*) '  *******************************'                    CPP09750
      WRITE(4,*) '  ***END OF THE CPHF PROCEDURE***'                    CPP09760
      WRITE(4,*) '  *******************************'                    CPP09770
      DO 110 I=1,2                                                      CPP09780
  110 WRITE(4,*) ' '                                                    CPP09790
C                                                                       CPP09800
      CALL TSTOP(6)                                                     CPP09810
      STOP                                                              CPP09820
      END                                                               CPP09830
      SUBROUTINE REGIST(NIJ,KIJ)                                        CPP09840
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP09850
      DIMENSION NIJ(NTRI,2),KIJ(NTRI,2)                                 CPP09860
      DIMENSION LOTYP(256)                                              CPP09870
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP09880
      COMMON/CIOFF/IIUNQ,IJUNQ,ICIUNQ,IPX(5050),JPX(5050)               CPP09890
      COMMON/CIPEX/PEXC(100),INDPEX(20,100),MOPEX(20),ICIMO(20)         CPP09900
      COMMON/CIPMO/IIMO(5050),JJMO(5050)                                CPP09910
      COMMON/CIPOS/ICITYP(5050),ICPOS(200,100)                          CPP09920
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP09930
      COMMON/CITYP/IJGRP,IBM(200),JBM(200),IBNUM(200)                   CPP09940
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP09950
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP09960
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP09970
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP09980
      COMMON/CI102/NIND,NDEP                                            CPP09990
      DATA A00,ONE,TWO / 0.0D+00 , 1.0D+00 , 2.0D+00 /                  CPP10000
    1 FORMAT(//,2X,' REGISTERED NUMBERS'/                               CPP10010
     1 2X,' KVIR = ',I5/                                                CPP10020
     2 2X,' NIND = ',I5/                                                CPP10030
     3 2X,' NDEP = ',I5/                                                CPP10040
     4 2X,' NDIF = ',I5)                                                CPP10050
    2 FORMAT(//,2X,' NIJ MATRIX'/2X,' NO.',5X,' NIJ(1)',5X,' NIJ(2)'/)  CPP10060
    3 FORMAT(2X,I4,4X,I5,7X,I5)                                         CPP10070
    4 FORMAT(//,2X,' KIJ MATRIX'/2X,' NO.',5X,' KIJ(1)',5X,' KIJ(2)'/)  CPP10080
    5 FORMAT(//,2X,' NO.',8X,' FOCC',6X,' NSORB',7X,' NSTR',7X,' NEND'/)CPP10090
    6 FORMAT(2X,I4,3X,F10.5,7X,I5,7X,I5,7X,I5)                          CPP10100
    7 FORMAT(//,2X,' MOPEX VECTOR '//16I5)                              CPP10110
    8 FORMAT(//,2X,' LOTYP VECTOR '//16I5)                              CPP10120
    9 FORMAT(//,2X,' ICIMO VECTOR '//16I5)                              CPP10130
   10 FORMAT(//,2X,' INDPEX VECTOR'/)                                   CPP10140
   11 FORMAT(2X,I5,5X,10I5)                                             CPP10150
   12 FORMAT(//,2X,' IIUNQ = ',I6,4X,' IJUNQ = ',I6,4X,' ICIUNQ = ',I6) CPP10160
   13 FORMAT(//,4X,' IJCI    IEX    JEX    III    JJJ   IIMO   JJMO',   CPP10170
     1          2X,' ICMO   JCMO  IJGRP'/)                              CPP10180
   14 FORMAT(2X,10I7)                                                   CPP10190
   15 FORMAT(//,2X,' IJGRP = ',I6)                                      CPP10200
   16 FORMAT(//,4X,' IBCI    ICM    JCM    NUM'/)                       CPP10210
C                                                                       CPP10220
      KVIR=NBASIS-JOCC                                                  CPP10230
C                                                                       CPP10240
      NSTR(1)=1                                                         CPP10250
      NEND(1)=NSORB(1)                                                  CPP10260
      DO 101 ITYP=2,NTYPEP                                              CPP10270
      NSTR(ITYP)=NEND(ITYP-1)+1                                         CPP10280
      NEND(ITYP)=NSTR(ITYP)+NSORB(ITYP)-1                               CPP10290
  101 CONTINUE                                                          CPP10300
      DO 103 ITYP=1,NTYPEP                                              CPP10310
      DO 102 I=NSTR(ITYP),NEND(ITYP)                                    CPP10320
  102 MOTYP(I)=ITYP                                                     CPP10330
  103 CONTINUE                                                          CPP10340
C                                                                       CPP10350
      II=0                                                              CPP10360
      DO 105 ITYP=2,NTYPEP                                              CPP10370
      DO 105 JTYP=1,ITYP-1                                              CPP10380
      DO 104 I=NSTR(ITYP),NEND(ITYP)                                    CPP10390
      DO 104 J=NSTR(JTYP),NEND(JTYP)                                    CPP10400
      II=II+1                                                           CPP10410
      NIJ(II,1)=I                                                       CPP10420
      NIJ(II,2)=J                                                       CPP10430
  104 CONTINUE                                                          CPP10440
  105 CONTINUE                                                          CPP10450
      NIND=II                                                           CPP10460
C                                                                       CPP10470
      II=0                                                              CPP10480
      DO 107 ITYP=1,NTYPEP                                              CPP10490
C*C   IF(NSORB(ITYP).EQ.1) GO TO 107                                    CPP10500
C*C   DO 106 I=NSTR(ITYP)+1,NEND(ITYP)                                  CPP10510
C*C   DO 106 J=NSTR(ITYP),I-1                                           CPP10520
      DO 106 I=NSTR(ITYP),NEND(ITYP)                                    CPP10530
      DO 106 J=NSTR(ITYP),I                                             CPP10540
      II=II+1                                                           CPP10550
      KIJ(II,1)=I                                                       CPP10560
      KIJ(II,2)=J                                                       CPP10570
  106 CONTINUE                                                          CPP10580
  107 CONTINUE                                                          CPP10590
      NDEP=II                                                           CPP10600
C                                                                       CPP10610
      NDIF=NTRI-NIND-NDEP                                               CPP10620
      IF(IPRNT.LE.2) GO TO 201                                          CPP10630
      WRITE(6,1) KVIR,NIND,NDEP,NDIF                                    CPP10640
      WRITE(6,2)                                                        CPP10650
      DO 108 I=1,NIND                                                   CPP10660
  108 WRITE(6,3) I,NIJ(I,1),NIJ(I,2)                                    CPP10670
      WRITE(6,4)                                                        CPP10680
      DO 109 I=1,NDEP                                                   CPP10690
  109 WRITE(6,3) I,KIJ(I,1),KIJ(I,2)                                    CPP10700
      WRITE(6,5)                                                        CPP10710
      DO 110 I=1,NTYPEP                                                 CPP10720
  110 WRITE(6,6) I,FOCC(I),NSORB(I),NSTR(I),NEND(I)                     CPP10730
C                                                                       CPP10740
C   INTERRELATE THE MO'S IN THE ACTIVE SPACE                            CPP10750
  201 CONTINUE                                                          CPP10760
      CALL RMREAD(LOTYP,20)                                             CPP10770
      IF(IPRNT.GT.2)                                                    CPP10780
     *WRITE(6,7) (MOPEX(I),I=1,NORPEX)                                  CPP10790
      IF(IPRNT.GT.2)                                                    CPP10800
     *WRITE(6,8) (LOTYP(I),I=1,NORPEX)                                  CPP10810
      DO 113 I=1,NORPEX                                                 CPP10820
      MM=MOPEX(I)                                                       CPP10830
      DO 112 J=1,NBASIS                                                 CPP10840
      IF(LOTYP(J).NE.MM) GO TO 112                                      CPP10850
      MTEMP=J                                                           CPP10860
  112 CONTINUE                                                          CPP10870
      ICIMO(I)=MTEMP                                                    CPP10880
  113 CONTINUE                                                          CPP10890
      IF(IPRNT.GT.2)                                                    CPP10900
     *WRITE(6,9) (ICIMO(I),I=1,NORPEX)                                  CPP10910
C                                                                       CPP10920
C   GENERATE DIFFERENCE VECTOR                                          CPP10930
      IF(IPRNT.LE.2) GO TO 205                                          CPP10940
      WRITE(6,10)                                                       CPP10950
      DO 115 I=1,NCPEX                                                  CPP10960
      WRITE(6,11) I,(INDPEX(J,I),J=1,NORPEX)                            CPP10970
  115 CONTINUE                                                          CPP10980
  205 CONTINUE                                                          CPP10990
      IIUNQ=NCPEX                                                       CPP11000
      IJUNQ=0                                                           CPP11010
      DO 117 I=1,NCPEX                                                  CPP11020
      IPX(I)=I                                                          CPP11030
      JPX(I)=I                                                          CPP11040
      ICITYP(I)=I                                                       CPP11050
      DO 117 J=1,I                                                      CPP11060
      NN=0                                                              CPP11070
      DO 116 K=1,NORPEX                                                 CPP11080
  116 IF(INDPEX(K,I).NE.INDPEX(K,J)) NN=NN+1                            CPP11090
      NN=NN/2                                                           CPP11100
      IF(NN.NE.1) GO TO 117                                             CPP11110
      IJUNQ=IJUNQ+1                                                     CPP11120
      IPX(IIUNQ+IJUNQ)=I                                                CPP11130
      JPX(IIUNQ+IJUNQ)=J                                                CPP11140
  117 CONTINUE                                                          CPP11150
      ICIUNQ=IIUNQ+IJUNQ                                                CPP11160
      WRITE(6,12) IIUNQ,IJUNQ,ICIUNQ                                    CPP11170
      IF(IPRNT.GT.2)                                                    CPP11180
     *WRITE(6,13)                                                       CPP11190
C                                                                       CPP11200
      IJGRP=0                                                           CPP11210
      CALL IZERO(IBNUM,200)                                             CPP11220
      DO 120 IJCI=1,IJUNQ                                               CPP11230
      IEX=IPX(IIUNQ+IJCI)                                               CPP11240
      JEX=JPX(IIUNQ+IJCI)                                               CPP11250
      DO 118 K=1,NORPEX                                                 CPP11260
      IF(INDPEX(K,IEX).EQ.INDPEX(K,JEX)) GO TO 118                      CPP11270
      IF(INDPEX(K,IEX).NE.0) III=K                                      CPP11280
      IF(INDPEX(K,JEX).NE.0) JJJ=K                                      CPP11290
  118 CONTINUE                                                          CPP11300
      IIMO(IJCI)=MOTYP(ICIMO(III))                                      CPP11310
      JJMO(IJCI)=MOTYP(ICIMO(JJJ))                                      CPP11320
      ICM=ICIMO(III)                                                    CPP11330
      JCM=ICIMO(JJJ)                                                    CPP11340
      IF(IJCI.EQ.1) GO TO 206                                           CPP11350
      ITEST=0                                                           CPP11360
      DO 119 IBCI=1,IJGRP                                               CPP11370
      IF(ICM.EQ.IBM(IBCI).AND.JCM.EQ.JBM(IBCI)) THEN                    CPP11380
        IBNUM(IBCI)=IBNUM(IBCI)+1                                       CPP11390
        INUM=IBNUM(IBCI)                                                CPP11400
        ITEST=ITEST+1                                                   CPP11410
        ICITYP(IIUNQ+IJCI)=IIUNQ+IBCI                                   CPP11420
        ICPOS(IBCI,INUM)=IJCI                                           CPP11430
      END IF                                                            CPP11440
  119 CONTINUE                                                          CPP11450
      IF(ITEST.NE.0) GO TO 207                                          CPP11460
  206 CONTINUE                                                          CPP11470
      IJGRP=IJGRP+1                                                     CPP11480
      IBM(IJGRP)=ICM                                                    CPP11490
      JBM(IJGRP)=JCM                                                    CPP11500
      IBNUM(IJGRP)=1                                                    CPP11510
      ICITYP(IIUNQ+IJCI)=IIUNQ+IJGRP                                    CPP11520
      ICPOS(IJGRP,1)=IJCI                                               CPP11530
C                                                                       CPP11540
  207 CONTINUE                                                          CPP11550
      IF(IPRNT.LE.2) GO TO 120                                          CPP11560
      WRITE(6,14) IJCI,IEX,JEX,III,JJJ,IIMO(IJCI),JJMO(IJCI),           CPP11570
     *            ICIMO(III),ICIMO(JJJ),ICITYP(IIUNQ+IJCI)              CPP11580
  120 CONTINUE                                                          CPP11590
C                                                                       CPP11600
      WRITE(6,15) IJGRP                                                 CPP11610
      IF(IPRNT.LE.2) GO TO 210                                          CPP11620
      WRITE(6,16)                                                       CPP11630
      DO 122 IBCI=1,IJGRP                                               CPP11640
      ICM=IBM(IBCI)                                                     CPP11650
      JCM=JBM(IBCI)                                                     CPP11660
      NUM=IBNUM(IBCI)                                                   CPP11670
      WRITE(6,14) IBCI,ICM,JCM,NUM                                      CPP11680
      DO 121 INUM=1,NUM                                                 CPP11690
      WRITE(6,14) ICPOS(IBCI,INUM)                                      CPP11700
  121 CONTINUE                                                          CPP11710
  122 CONTINUE                                                          CPP11720
C                                                                       CPP11730
  210 CONTINUE                                                          CPP11740
      RETURN                                                            CPP11750
      END                                                               CPP11760
      SUBROUTINE CORSIZ(MAXCOR,N3NXX,NADD,MAXVC,IPOL)                   CPP11770
      DIMENSION NADD(10)                                                CPP11780
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP11790
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP11800
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP11810
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPP11820
      COMMON/SIZES/MAXVEC,MAXVX                                         CPP11830
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP11840
      COMMON/CI102/NIND,NDEP                                            CPP11850
    1 FORMAT(//,2X,' %%%%%%%%%%%%%%%%%%%%%%%%%%%'/                      CPP11860
     1          2X,' %%%NOW YOU ARE IN CORSIZ%%%'/                      CPP11870
     2          2X,' %%%%%%%%%%%%%%%%%%%%%%%%%%%'/)                     CPP11880
    2 FORMAT(//,2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)                CPP11890
    3 FORMAT(//,2X,' *************************************************'/CPP11900
     1          2X,' ***REQUIRED MEMORY EXCEEDS MAXCOR-----STOP !!!***'/CPP11910
     2          2X,' *************************************************'/CPP11920
     3          2X,' IGMAX = ',I10,5X,' MAXCOR = ',I10/)                CPP11930
    4 FORMAT(//,2X,' ::::::::::::::::::::::::::::::::::'/               CPP11940
     1          2X,' :::REQUIRED MEMORY FOR THIS RUN:::'/               CPP11950
     2          2X,' ::::::::::::::::::::::::::::::::::'/               CPP11960
     3          2X,' NADD  = ',I10,5X,' MAXVC  = ',I10/                 CPP11970
     4          2X,' IGMAX = ',I10,5X,' MAXCOR = ',I10/)                CPP11980
    5 FORMAT(//,2X,' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'/                 CPP11990
     1          2X,' %%%NOW YOU ARE LEAVING CORSIZ%%%'/                 CPP12000
     2          2X,' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'/)                CPP12010
C                                                                       CPP12020
C   THIS SUBROUTINE WILL CHECK THE MAXIMUM CORE MEMORY REQUIRED         CPP12030
C   FOR THE CALCULATION                                                 CPP12040
C                                                                       CPP12050
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::CPP12060
C     IC1  = EIG(NBASIS)                                                CPP12070
C     IC2  = OCC(NBASIS)                                                CPP12080
C     IC3  = EAO OR ESO(NBFAO,NBASIS)                                   CPP12090
C     IC4  = NIJ(NTRI,2)                                                CPP12100
C     IC5  = KIJ(NTRI,2)                                                CPP12110
C     IC6  = FIJ(NTYPEP,NCTRI)                                          CPP12120
C     IC7  = AIJ(NTYTRI,NCTRI)                                          CPP12130
C     IC8  = BIJ(NTYTRI,NCTRI)                                          CPP12140
C     IC9  = E1A(N3N)                                                   CPP12150
C     IC10 = E2A(N3N,N3N)                                               CPP12160
C     IC11 = D1W(N3N)                                                   CPP12170
C     IC12 = E1T(N3NXX)                                                 CPP12180
C     IC13 = EPS(NTRI)                                                  CPP12190
C     IC14 = ZETA(NTRI,NTYPEP)                                          CPP12200
C     IC15 = HCIA(NCTRI,N3NXX)                                          CPP12210
C     IC16 = HCI(NCTRI)                                                 CPP12220
C     IC17 = BAF2(NCPEX,N3NXX)                                          CPP12230
C     IC18 = CIA(NCPEX,N3NXX)                                           CPP12240
C     IC19 = A22(NCPEX,NCPEX)                                           CPP12250
C     IC20 = A12(NIND,NCPEX)                                            CPP12260
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::CPP12270
C                                                                       CPP12280
      IC1=1                                                             CPP12290
      IC2=IC1+NBASIS                                                    CPP12300
      IC3=IC2+NBASIS                                                    CPP12310
      IC4=IC3+NBFAO*NBASIS                                              CPP12320
      IA4=IC4+IC4-1                                                     CPP12330
      IC5=IC4+NTRI                                                      CPP12340
      IA5=IC5+IC5-1                                                     CPP12350
      IC6=IC5+NTRI                                                      CPP12360
      IC7=IC6+NTYPEP*NCTRI                                              CPP12370
      IC8=IC7+NTYTRI*NCTRI                                              CPP12380
      IC9=IC8+NTYTRI*NCTRI                                              CPP12390
      IC10=IC9+N3N                                                      CPP12400
      IC11=IC10+N3N*N3N                                                 CPP12410
      IC12=IC11+N3N                                                     CPP12420
      IC13=IC12+N3NXX                                                   CPP12430
      IC14=IC13+NTRI                                                    CPP12440
      IC15=IC14+NTRI*NTYPEP                                             CPP12450
      IC16=IC15+NCTRI*N3NXX                                             CPP12460
      IC17=IC16+NCTRI                                                   CPP12470
      IC18=IC17+NCPEX*N3NXX                                             CPP12480
      IC19=IC18+NCPEX*N3NXX                                             CPP12490
      IC20=IC19+NCPEX*NCPEX                                             CPP12500
      ICKEEP=IC20+NIND*NCPEX                                            CPP12510
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::CPP12520
C                                                                       CPP12530
C                                                                       CPP12540
      WRITE(3,1)                                                        CPP12550
      NADD(1)=N3N*2                                                     CPP12560
      NADD(2)=N3N*2                                                     CPP12570
      NADD(3)=N3NXX*2                                                   CPP12580
      NADD(4)=N3NXX*2                                                   CPP12590
      NADD(5)=N3N*2                                                     CPP12600
C                                                                       CPP12610
C   READ IN AO-MO EIGENVECTORS                                          CPP12620
C                                                                       CPP12630
C   FORM REGISTER MATRICES                                              CPP12640
      WRITE(3,20)                                                       CPP12650
   20 FORMAT(//,2X,' NOW YOU ARE IN REGIST'/)                           CPP12660
C.................NIJ     KIJ.....                                      CPP12670
C***  CALL REGIST(IA(IA4),IA(IA5))                                      CPP12680
C                                                                       CPP12690
C   READ IN FIRST AND SECOND DERIVATIVE INFORMATION                     CPP12700
      WRITE(3,21)                                                       CPP12710
   21 FORMAT(//,2X,' NOW YOU ARE IN DERMAT'/)                           CPP12720
      IC21=ICKEEP                                                       CPP12730
      IC22=IC21+NBATRI                                                  CPP12740
      IC23=IC22+NBFAO*NBFAO                                             CPP12750
      ICMAX=IC23+NBFAO*NBFAO                                            CPP12760
      WRITE(3,2) ICMAX,MAXCOR                                           CPP12770
      IGMAX=ICMAX                                                       CPP12780
C.................EAO     E1A     E2A      SS       U        T........  CPP12790
C***  CALL DERMAT(CC(IC3),CC(IC9),CC(IC10),CC(IC21),CC(IC22),CC(IC23),  CPP12800
C..............................                                         CPP12810
C*** 1            NSECMX,ISTEP)                                         CPP12820
C                                                                       CPP12830
C   CALCULATE SCF ENERGY FOR A TEST                                     CPP12840
      WRITE(3,22)                                                       CPP12850
   22 FORMAT(//,2X,' NOW YOU ARE IN E0CAL'/)                            CPP12860
      IC21=ICKEEP                                                       CPP12870
      IC22=IC21+NTRI                                                    CPP12880
      IC23=IC22+NTRI                                                    CPP12890
      IC24=IC23+NTRI                                                    CPP12900
      IC25=IC24+NBFAO*NBFAO                                             CPP12910
      IC26=IC25+NBFAO*NBFAO                                             CPP12920
      IC27=IC26+NTRI*NTYPES                                             CPP12930
      IC28=IC27+NTRI*NTYPES                                             CPP12940
      IC29=IC28+NTRI*NTYPES                                             CPP12950
      IA29=IC29+IC29-1                                                  CPP12960
      ICMAX=IC29+MAXBF2                                                 CPP12970
      WRITE(3,2) ICMAX,MAXCOR                                           CPP12980
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPP12990
C................ESO     HSO      HMO      TM       U        T........  CPP13000
C***  CALL E0CAL(CC(IC3),CC(IC21),CC(IC22),CC(IC23),CC(IC24),CC(IC25),  CPP13010
C................DP       DQ       DM       MIJ      PQ.............    CPP13020
C*** 1           CC(IC26),CC(IC27),CC(IC28),IA(IA29),IA(IA29+MAXIJ))    CPP13030
C                                                                       CPP13040
C   CALCULATE SCF FIRST ENERGY DERIVATIVES                              CPP13050
      WRITE(3,23)                                                       CPP13060
   23 FORMAT(//,2X,' NOW YOU ARE IN EADER'/)                            CPP13070
      IC21=ICKEEP                                                       CPP13080
      IC22=IC21+N3N                                                     CPP13090
      IC23=IC22+N3N                                                     CPP13100
      IC24=IC23+NTRI                                                    CPP13110
      IC25=IC24+NTRI                                                    CPP13120
      ICMAX=IC25+NTRI*N3N*NTREAD                                        CPP13130
      WRITE(3,2) ICMAX,MAXCOR                                           CPP13140
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPP13150
C................D1W      E1T      EPS      D1O      D1T......          CPP13160
C***  CALL EADER(CC(IC11),CC(IC12),CC(IC13),CC(IC21),CC(IC22),          CPP13170
C................SM       HM       TM..................                 CPP13180
C*** 1           CC(IC23),CC(IC24),CC(IC25),N3N*NTREAD)                 CPP13190
C                                                                       CPP13200
C   FORM HA MATRICES                                                    CPP13210
      WRITE(3,24)                                                       CPP13220
   24 FORMAT(//,2X,' NOW YOU ARE IN HAMAT'/)                            CPP13230
      IC21=ICKEEP                                                       CPP13240
      IC22=IC21+N3N                                                     CPP13250
      IC23=IC22+NTRI*N3N                                                CPP13260
      ICMAX=IC23+NTRI*N3N*NTREAD                                        CPP13270
      WRITE(3,2) ICMAX,MAXCOR                                           CPP13280
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPP13290
C................FIJ     AIJ     BIJ     D1W      E1T      HCIA.....    CPP13300
C***  CALL HAMAT(CC(IC6),CC(IC7),CC(IC8),CC(IC11),CC(IC12),CC(IC15),    CPP13310
C................D1E      HM       TM........................           CPP13320
C*** 1           CC(IC21),CC(IC22),CC(IC23),N3N*NTREAD,N3NXX)           CPP13330
C                                                                       CPP13340
C   FORM DERIVATIVE LAGRANGIAN MATRICES                                 CPP13350
      WRITE(3,25)                                                       CPP13360
   25 FORMAT(//,2X,' NOW YOU ARE IN FAMAT'/)                            CPP13370
      IC21=ICKEEP                                                       CPP13380
      IC22=IC21+NTRI                                                    CPP13390
      IC23=IC22+NBASIS*NBASIS                                           CPP13400
      IC24=IC23+NTRI*N3N                                                CPP13410
      ICMAX=IC24+NTRI*N3N*NTREAD                                        CPP13420
      WRITE(3,2) ICMAX,MAXCOR                                           CPP13430
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPP13440
C................FF       EPA      HM       TM...................       CPP13450
C***  CALL FAMAT(CC(IC21),CC(IC22),CC(IC23),CC(IC24),N3N*NTREAD))       CPP13460
C                                                                       CPP13470
C   CALCULATE BA MATRICES                                               CPP13480
C   FOR INDEPENDENT PAIRS                                               CPP13490
      WRITE(3,26)                                                       CPP13500
   26 FORMAT(//,2X,' NOW YOU ARE IN BAIND'/)                            CPP13510
  210 CONTINUE                                                          CPP13520
      NADD(1)=NADD(1)/2                                                 CPP13530
      IC21=ICKEEP                                                       CPP13540
      IC22=IC21+NBFAO*NBFAO                                             CPP13550
      IC23=IC22+NBFAO*NBFAO                                             CPP13560
      IC24=IC23+NTRI                                                    CPP13570
      IC25=IC24+NBASIS*NBASIS                                           CPP13580
      IC26=IC25+NTRI*NTYPEP*NADD(1)                                     CPP13590
      IC27=IC26+NTRI*NTYPEP*NADD(1)                                     CPP13600
      IC28=IC27+NTRI*NTYPEP*NADD(1)                                     CPP13610
      IA28=IC28+IC28-1                                                  CPP13620
      ICMAX=IC28+MAXBF2                                                 CPP13630
      WRITE(3,27) N3N,NADD(1),ICMAX                                     CPP13640
   27 FORMAT(2X,' N3N   = ',I8,4X,' NADD  = ',I8,4X,' ICMAX = ',I8/)    CPP13650
      IF(NADD(1).EQ.1) GO TO 211                                        CPP13660
      IF(ICMAX.LT.MAXCOR) GO TO 211                                     CPP13670
      GO TO 210                                                         CPP13680
  211 CONTINUE                                                          CPP13690
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPP13700
C................ESO     NIJ     ZETA     U        T        B0.......   CPP13710
C***  CALL BAIND(CC(IC3),IA(IA4),CC(IC14),CC(IC21),CC(IC22),CC(IC23),   CPP13720
C................EPA      DP       DQ       DM       SM.......          CPP13730
C*** 1           CC(IC24),CC(IC25),CC(IC26),CC(IC27),CC(IC26),          CPP13740
C................MIJ      PQ.....................                       CPP13750
C*** 2           IA(IA28),IA(IA28+MAXIJ),NADD(1))                       CPP13760
C                                                                       CPP13770
C   FOR DEPENDENT PAIRS                                                 CPP13780
      WRITE(3,28)                                                       CPP13790
   28 FORMAT(//,2X,' NOW YOU ARE IN BADEP'/)                            CPP13800
  220 CONTINUE                                                          CPP13810
      NADD(2)=NADD(2)/2                                                 CPP13820
      IC21=ICKEEP                                                       CPP13830
      IC22=IC21+NBFAO*NBFAO                                             CPP13840
      IC23=IC22+NBFAO*NBFAO                                             CPP13850
      IC24=IC23+NTRI                                                    CPP13860
      IC25=IC24+NTRI                                                    CPP13870
      IC26=IC25+NTRI*NADD(2)                                            CPP13880
      IC27=IC26+NTRI*NADD(2)                                            CPP13890
      IC28=IC27+NTRI*NADD(2)                                            CPP13900
      IA28=IC28+IC28-1                                                  CPP13910
      ICMAX=IC28+MAXBF2                                                 CPP13920
      WRITE(3,29) N3N,NADD(2),ICMAX                                     CPP13930
   29 FORMAT(2X,' N3N   = ',I8,4X,' NADD  = ',I8,4X,' ICMAX = ',I8/)    CPP13940
      IF(NADD(2).EQ.1) GO TO 221                                        CPP13950
      IF(ICMAX.LT.MAXCOR) GO TO 221                                     CPP13960
      GO TO 220                                                         CPP13970
  221 CONTINUE                                                          CPP13980
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPP13990
C................EIG     ESO     KIJ     U        T.......              CPP14000
C***  CALL BADEP(CC(IC1),CC(IC3),IA(IA5),CC(IC21),CC(IC22),             CPP14010
C................B0       FF       DP       DQ       DM.......          CPP14020
C*** 1           CC(IC23),CC(IC24),CC(IC25),CC(IC26),CC(IC27),          CPP14030
C................MIJ      PQ.....................                       CPP14040
C*** 2           IA(IA28),IA(IA28+MAXIJ),NADD(2))                       CPP14050
C                                                                       CPP14060
C   CALCULATE H AND E MATRICES                                          CPP14070
      WRITE(3,30)                                                       CPP14080
   30 FORMAT(//,2X,' NOW YOU ARE IN HEMAT'/)                            CPP14090
      IC21=ICKEEP                                                       CPP14100
      IC22=IC21+NTRI                                                    CPP14110
      IC23=IC22+NTRI                                                    CPP14120
      IC24=IC23+NTRI*NTYPEP                                             CPP14130
      IC25=IC24+NTRI*NTYPEP                                             CPP14140
      IC26=IC25+NTRI*NTYPEP                                             CPP14150
      IC27=IC26+NBASIS*NBASIS                                           CPP14160
      IC28=IC27+NTRI                                                    CPP14170
      IC29=IC28+NBASIS*NBASIS                                           CPP14180
      IC30=IC29+NBASIS*NBASIS                                           CPP14190
      IC31=IC30+NBASIS*NBASIS                                           CPP14200
      IA31=IC31+IC31-1                                                  CPP14210
      ICMAX=IC31+MAXBF2                                                 CPP14220
      WRITE(3,2) ICMAX,MAXCOR                                           CPP14230
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPP14240
C................ESO     FIJ     AIJ     BIJ     HCI......              CPP14250
C***  CALL HEMAT(CC(IC3),CC(IC6),CC(IC7),CC(IC8),CC(IC16),              CPP14260
C................HSO      HMO      DP       DQ       DM       EIJ...... CPP14270
C*** 1           CC(IC21),CC(IC22),CC(IC23),CC(IC24),CC(IC25),CC(IC26), CPP14280
C................ZIJ      EE       U        T........                   CPP14290
C*** 2           CC(IC27),CC(IC28),CC(IC29),CC(IC30),                   CPP14300
C................MIJ      PQ.............                               CPP14310
C*** 3           IA(IA31),IA(IA30+MAXIJ))                               CPP14320
C                                                                       CPP14330
C   FORM BF MATRICES                                                    CPP14340
      WRITE(3,31)                                                       CPP14350
   31 FORMAT(//,2X,' NOW YOU ARE IN BFMAT'/)                            CPP14360
      IC21=ICKEEP                                                       CPP14370
      IC22=IC21+NTRI*3                                                  CPP14380
      IC23=IC22+NTRI*NTYPES*3                                           CPP14390
      IC24=IC23+NBASIS*NBASIS                                           CPP14400
      ICMAX=IC24+NTRI                                                   CPP14410
      WRITE(3,2) ICMAX,MAXCOR                                           CPP14420
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPP14430
C................NIJ     KIJ     FIJ     E1T      HCIA     HF.......    CPP14440
C***  CALL BFMAT(IA(IA4),IA(IA5),CC(IC6),CC(IC12),CC(IC15),CC(IC21),    CPP14450
C................HM       EPF      BF...................                CPP14460
C*** 1           CC(IC22),CC(IC23),CC(IC24),N3NXX,ISTEP)                CPP14470
C                                                                       CPP14480
C   FORM BAF MATRICES                                                   CPP14490
      WRITE(3,32)                                                       CPP14500
   32 FORMAT(//,2X,' NOW YOU ARE IN AB2CI'/)                            CPP14510
      IC21=ICKEEP                                                       CPP14520
      IC22=IC21+NTRI*N3N                                                CPP14530
      ICMAX=IC22+NBASIS*NBASIS                                          CPP14540
      WRITE(3,2) ICMAX,MAXCOR                                           CPP14550
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPP14560
C................NIJ     E1T      HCIA     HCI      BAF2     A22......  CPP14570
C***  CALL AB2CI(IA(IA4),CC(IC12),CC(IC15),CC(IC16),CC(IC17),CC(IC19),  CPP14580
C................A12      SA       EIJ............                      CPP14590
C*** 1           CC(IC20),CC(IC21),CC(IC22),N3NXX)                      CPP14600
C                                                                       CPP14610
C   SOLVE SIMULTANEOUS EQUATIONS ITERATIVELY                            CPP14620
      WRITE(3,33)                                                       CPP14630
   33 FORMAT(//,2X,' NOW YOU ARE IN CPHF'/)                             CPP14640
  230 CONTINUE                                                          CPP14650
      NADD(3)=NADD(3)/2                                                 CPP14660
      MAXVC=MIN0(MAXVEC,NADD(3)*MAXVX)                                  CPP14670
      NINDPX=NIND+NCPEX                                                 CPP14680
      IC21=ICKEEP                                                       CPP14690
      IC22=IC21+NBFAO*NBFAO                                             CPP14700
      IC23=IC22+NBFAO*NBFAO                                             CPP14710
      IC24=IC23+NBFAO*NBFAO                                             CPP14720
      IC25=IC24+NTRI                                                    CPP14730
      IC26=IC25+NINDPX                                                  CPP14740
      IC27=IC26+NINDPX                                                  CPP14750
      IC28=IC27+NTYPEP*(NTRI*NADD(3))                                   CPP14760
      IC29=IC28+NTYPEP*(NTRI*NADD(3))                                   CPP14770
      IC30=IC29+NTYPEP*(NTRI*NADD(3))                                   CPP14780
      IC31=IC30+MAXVC*(MAXVC+NADD(3))                                   CPP14790
      IC32=IC31+NINDPX*MAXVC*2                                          CPP14800
      IA32=IC32+IC32-1                                                  CPP14810
      ICMAX=IC32+MAXBF2                                                 CPP14820
      WRITE(3,34) NIND,N3NXX,NADD(3),MAXVC,ICMAX                        CPP14830
   34 FORMAT(2X,' NIND  = ',I8,4X,' N3NXX = ',I8,4X,' NADD = ',I8/      CPP14840
     1       2X,' MAXVC = ',I8,4X,' ICMAX = ',I8/)                      CPP14850
      IF(NADD(3).EQ.1) GO TO 231                                        CPP14860
      IF(ICMAX.LT.MAXCOR) GO TO 231                                     CPP14870
      GO TO 230                                                         CPP14880
  231 CONTINUE                                                          CPP14890
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPP14900
C...............ESO     NIJ     KIJ     ZETA     BAF2     CIA......     CPP14910
C***  CALL CPHF(CC(IC3),IA(IA4),IA(IA5),CC(IC14),CC(IC17),CC(IC18),     CPP14920
C...............A22      A12      A........                             CPP14930
C*** 1          CC(IC19),CC(IC20),CC(IC21),                             CPP14940
C...............U        T        B0       AA       TT       DP.......  CPP14950
C*** 2          CC(IC22),CC(IC23),CC(IC24),CC(IC25),CC(IC26),CC(IC27),  CPP14960
C...............DQ       DM       W        BB       MIJ......           CPP14970
C*** 3          CC(IC28),CC(IC29),CC(IC30),CC(IC31),IA(IA32),           CPP14980
C...............PQ........................................              CPP14990
C*** 4          IA(IA32+MAXIJ),NINDPX,MAXVC,NADD(3),N3NXX)              CPP15000
C                                                                       CPP15010
C   CALCULATE DEPENDENT PAIRS OF U MATRICES                             CPP15020
      WRITE(3,35)                                                       CPP15030
   35 FORMAT(//,2X,' NOW YOU ARE IN UCMAT'/)                            CPP15040
  240 CONTINUE                                                          CPP15050
      NADD(4)=NADD(4)/2                                                 CPP15060
      IC21=ICKEEP                                                       CPP15070
      IC22=IC21+NBFAO*NBFAO                                             CPP15080
      IC23=IC22+NBFAO*NBFAO                                             CPP15090
      IC24=IC23+NTRI                                                    CPP15100
      IC25=IC24+NTRI*NADD(4)                                            CPP15110
      IC26=IC25+NTRI*NADD(4)                                            CPP15120
      IC27=IC26+NTRI*NADD(4)                                            CPP15130
      IC28=IC27+NTRI                                                    CPP15140
      IA28=IC28+IC28-1                                                  CPP15150
      ICMAX=IC28+MAXBF2                                                 CPP15160
      WRITE(3,36) N3NXX,NADD(4),ICMAX                                   CPP15170
   36 FORMAT(2X,' N3NXX  = ',I8,4X,' NADD = ',I8,4X,' ICMAX = ',I8/)    CPP15180
      IF(NADD(4).EQ.1) GO TO 241                                        CPP15190
      IF(ICMAX.LT.MAXCOR) GO TO 241                                     CPP15200
      GO TO 240                                                         CPP15210
  241 CONTINUE                                                          CPP15220
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPP15230
C................EIG     ESO     NIJ     KIJ     U........              CPP15240
C***  CALL UCMAT(CC(IC1),CC(IC3),IA(IA4),IA(IA5),CC(IC21),              CPP15250
C................T        B        DP       DQ       DM       DEL...... CPP15260
C*** 1           CC(IC22),CC(IC23),CC(IC24),CC(IC25),CC(IC26),CC(IC27), CPP15270
C................MIJ      PQ...........................                 CPP15280
C*** 2           IA(IA28),IA(IA28+MAXIJ),NADD(4),N3NXX)                 CPP15290
C                                                                       CPP15300
C   COMPLETE U MATRICES                                                 CPP15310
      WRITE(3,37)                                                       CPP15320
   37 FORMAT(//,2X,' NOW YOU ARE IN UXMAT'/)                            CPP15330
      IC21=ICKEEP                                                       CPP15340
      IC22=IC21+NTRI*N3NXX                                              CPP15350
      IC23=IC22+NTRI*N3NXX                                              CPP15360
      ICMAX=IC23+NBASIS*NBASIS                                          CPP15370
      WRITE(3,2) ICMAX,MAXCOR                                           CPP15380
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPP15390
C................KIJ     B        SM       U..............              CPP15400
C***  CALL UXMAT(IA(IA5),CC(IC21),CC(IC22),CC(IC23),N3NXX)              CPP15410
C                                                                       CPP15420
C   CALCULATE WA MATRICES                                               CPP15430
      WRITE(3,38)                                                       CPP15440
   38 FORMAT(//,2X,' NOW YOU ARE IN WAMAT'/)                            CPP15450
  250 CONTINUE                                                          CPP15460
      NADD(5)=NADD(5)/2                                                 CPP15470
      IC21=ICKEEP                                                       CPP15480
      IC22=IC21+NBASIS*NBASIS                                           CPP15490
      IC23=IC22+NBASIS*NBASIS                                           CPP15500
      IC24=IC23+NBASIS*NBASIS                                           CPP15510
      IC25=IC24+NBASIS*NBASIS                                           CPP15520
      IC26=IC25+NTRI*NTYPEP*NADD(5)                                     CPP15530
      IC27=IC26+NTRI*NTYPEP*NADD(5)                                     CPP15540
      IC28=IC27+NTRI*NTYPEP*NADD(5)                                     CPP15550
      IA28=IC28+IC28-1                                                  CPP15560
      ICMAX=IC28+MAXBF2                                                 CPP15570
      WRITE(3,39) N3N,NADD(5),ICMAX                                     CPP15580
   39 FORMAT(2X,' N3N   = ',I8,4X,' NADD  = ',I8,4X,' ICMAX = ',I8/)    CPP15590
      IF(NADD(5).EQ.1) GO TO 251                                        CPP15600
      IF(ICMAX.LT.MAXCOR) GO TO 251                                     CPP15610
      GO TO 250                                                         CPP15620
  251 CONTINUE                                                          CPP15630
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPP15640
C................ESO     ZETA     U        T        EPA      WA.......  CPP15650
C***  CALL WAMAT(CC(IC3),CC(IC14),CC(IC21),CC(IC22),CC(IC23),CC(IC24),  CPP15660
C................SA       DP       DQ       DM       MIJ......          CPP15670
C*** 1           CC(IC25),CC(IC26),CC(IC27),CC(IC25),IA(IA28),          CPP15680
C................PQ.....................                                CPP15690
C    2           IA(IA28+MAXIJ),NADD(5))                                CPP15700
C                                                                       CPP15710
C   CALCULATE SCF SECOND DERIVATIVES                                    CPP15720
      WRITE(3,40)                                                       CPP15730
   40 FORMAT(//,2X,' NOW YOU ARE IN SCF2W'/)                            CPP15740
      IC21=ICKEEP                                                       CPP15750
      IC22=IC21+N3N*N3N                                                 CPP15760
      IC23=IC22+N3N*N3N                                                 CPP15770
      IC24=IC23+N3N*N3N                                                 CPP15780
      IC25=IC24+NBASIS*NBASIS                                           CPP15790
      IC26=IC25+NBASIS*NBASIS                                           CPP15800
      IC27=IC26+NTRI*N3N                                                CPP15810
      ICMAX=IC27+NBASIS*NBASIS*N3N                                      CPP15820
      WRITE(3,2) ICMAX,MAXCOR                                           CPP15830
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPP15840
C................E2M      E2P      E2Q      EPA      UA       SA....... CPP15850
C***  CALL SCF2W(CC(IC21),CC(IC22),CC(IC23),CC(IC24),CC(IC25),CC(IC26), CPP15860
C................WA.......                                              CPP15870
C*** 1           CC(IC27))                                              CPP15880
C                                                                       CPP15890
      WRITE(3,41)                                                       CPP15900
   41 FORMAT(//,2X,' NOW YOU ARE IN SCF2ND'/)                           CPP15910
      IC21=ICKEEP                                                       CPP15920
      IC22=IC21+N3N*N3N                                                 CPP15930
      IC23=IC22+N3N*N3N                                                 CPP15940
      IC24=IC23+N3N*N3N                                                 CPP15950
      IC25=IC24+NTRI*N3N                                                CPP15960
      ICMAX=IC25+NBASIS*NBASIS                                          CPP15970
      WRITE(3,2) ICMAX,MAXCOR                                           CPP15980
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPP15990
C.................E2A      HCIA     CIA      E2M      E2C      E2T......CPP16000
C***  CALL SCF2ND(CC(IC10),CC(IC15),CC(IC18),CC(IC21),CC(IC22),CC(IC23),CPP16010
C.................SA       EIJ......                                    CPP16020
C*** 1            CC(IC24),CC(IC25))                                    CPP16030
C                                                                       CPP16040
C   CALCULATE DIPOLE MOMENT DERIVATIVES                                 CPP16050
      IF(IPOL.NE.0) GO TO 320                                           CPP16060
      WRITE(3,42)                                                       CPP16070
   42 FORMAT(//,2X,' NOW YOU ARE IN DIPMO'/)                            CPP16080
      IC21=ICKEEP                                                       CPP16090
      IC22=IC21+NTRI*3                                                  CPP16100
      ICMAX=IC22+NBASIS*NBASIS*N3N                                      CPP16110
      WRITE(3,2) ICMAX,MAXCOR                                           CPP16120
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPP16130
C................HCIA     CIA      HF       UA.......                   CPP16140
C***  CALL DIPMO(CC(IC15),CC(IC18),CC(IC21),CC(IC22))                   CPP16150
C                                                                       CPP16160
C   CALCULATE ELECTRIC POLARIZABILITIES                                 CPP16170
      WRITE(3,43)                                                       CPP16180
   43 FORMAT(//,2X,' NOW YOU ARE IN POLAR'/)                            CPP16190
      IC21=ICKEEP                                                       CPP16200
      IC22=IC21+NTRI*3                                                  CPP16210
      ICMAX=IC22+NBASIS*NBASIS*3                                        CPP16220
      WRITE(3,2) ICMAX,MAXCOR                                           CPP16230
      IGMAX=MAX0(ICMAX,IGMAX)                                           CPP16240
C................HCIA     CIA      HF       UF.......                   CPP16250
C***  CALL POLAR(CC(IC15),CC(IC18),CC(IC21),CC(IC22))                   CPP16260
  320 CONTINUE                                                          CPP16270
C                                                                       CPP16280
      IF(IGMAX.LT.MAXCOR) GO TO 400                                     CPP16290
C                                                                       CPP16300
C   REQUIRED MEMORY EXCEEDS MAXCOR-----STOP !!!                         CPP16310
      WRITE(6,3) ICMAX,MAXCOR                                           CPP16320
      WRITE(3,3) ICMAX,MAXCOR                                           CPP16330
      STOP                                                              CPP16340
C                                                                       CPP16350
  400 CONTINUE                                                          CPP16360
      WRITE(3,4) NADD(3),MAXVC,IGMAX,MAXCOR                             CPP16370
      WRITE(6,4) NADD(3),MAXVC,IGMAX,MAXCOR                             CPP16380
C                                                                       CPP16390
      WRITE(3,5)                                                        CPP16400
      RETURN                                                            CPP16410
      END                                                               CPP16420
      SUBROUTINE BARES(FIJ,AIJ,BIJ)                                     CPP16430
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP16440
      DIMENSION FIJ(NTYPEP,ICIUNQ),AIJ(NTYTRI,ICIUNQ),BIJ(NTYTRI,ICIUNQ)CPP16450
      DIMENSION F(20),A(20,20),B(20,20)                                 CPP16460
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP16470
      COMMON/CIOFF/IIUNQ,IJUNQ,ICIUNQ,IPX(5050),JPX(5050)               CPP16480
      COMMON/CIPEX/PEXC(100),INDPEX(20,100),MOPEX(20),ICIMO(20)         CPP16490
      COMMON/CIPMO/IIMO(5050),JJMO(5050)                                CPP16500
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP16510
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP16520
      COMMON/DENST/ALPC(465),BETC(465)                                  CPP16530
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP16540
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP16550
      COMMON/CI103/CI(100)                                              CPP16560
      DATA A00,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /   CPP16570
    1 FORMAT(//,2X,' THE CI COEFFICIENTS'/)                             CPP16580
    2 FORMAT(2X,I5,F15.8)                                               CPP16590
    3 FORMAT(//,2X,' THE NORMALIZATION OF CI COEFFICIENTS',5X,F15.8/)   CPP16600
    4 FORMAT(//,2X,' FIJ VECTOR, ICI = ',I5/)                           CPP16610
    5 FORMAT(8F10.5)                                                    CPP16620
    6 FORMAT(//,2X,' AIJ MATRIX, ICI = ',I5/)                           CPP16630
    7 FORMAT(//,2X,' BIJ MATRIX, ICI = ',I5/)                           CPP16640
    8 FORMAT(//,2X,' FOCC VECTOR FROM FIJ (NON-ZERO CI''S)'/5X,10F10.5/)CPP16650
    9 FORMAT(//,2X,' ALPA MATRIX FROM AIJ (NON-ZERO CI''S)'/)           CPP16660
   10 FORMAT(//,2X,' BETA MATRIX FROM BIJ (NON-ZERO CI''S)'/)           CPP16670
C                                                                       CPP16680
C   TRANSFER THE CI COEFFICEINTS                                        CPP16690
C                                                                       CPP16700
      IF(IPRNT.GT.2)                                                    CPP16710
     *WRITE(6,1)                                                        CPP16720
      CTOT=A00                                                          CPP16730
      DO 101 I=1,NCPEX                                                  CPP16740
      CI(I)=PEXC(I)                                                     CPP16750
      CTOT=CTOT+CI(I)*CI(I)                                             CPP16760
      IF(IPRNT.LT.2) GO TO 101                                          CPP16770
      WRITE(6,2) I,CI(I)                                                CPP16780
  101 CONTINUE                                                          CPP16790
      IF(IPRNT.GT.2)                                                    CPP16800
     *WRITE(6,3) CTOT                                                   CPP16810
C                                                                       CPP16820
C   FORM THE "BARE" COUPLING CONSTANT MATRICES                          CPP16830
      CALL ZERO(FIJ,NTYPEP*ICIUNQ)                                      CPP16840
      CALL ZERO(AIJ,NTYTRI*ICIUNQ)                                      CPP16850
      CALL ZERO(BIJ,NTYTRI*ICIUNQ)                                      CPP16860
C                                                                       CPP16870
C   ONE-ELECTRON COUPLING CONSTANTS                                     CPP16880
      DO 105 ICI=1,NCPEX                                                CPP16890
      DO 102 ITYP=1,NTYPES                                              CPP16900
      IF(FOCC(ITYP).NE.TWO.AND.FOCC(ITYP).NE.ONE) GO TO 102             CPP16910
      FIJ(ITYP,ICI)=FOCC(ITYP)*HALF                                     CPP16920
  102 CONTINUE                                                          CPP16930
      DO 103 I=1,NORPEX                                                 CPP16940
      IF(INDPEX(I,ICI).EQ.0) GO TO 103                                  CPP16950
      MTYP=MOTYP(ICIMO(I))                                              CPP16960
      FIJ(MTYP,ICI)=ONE                                                 CPP16970
  103 CONTINUE                                                          CPP16980
      IF(IPRNT.LE.2) GO TO 105                                          CPP16990
      WRITE(6,4) ICI                                                    CPP17000
      WRITE(6,5) (FIJ(ITYP,ICI),ITYP=1,NTYPEP)                          CPP17010
  105 CONTINUE                                                          CPP17020
C                                                                       CPP17030
C   TWO-ELECTRON COUPLING CONSTANTS                                     CPP17040
C   FOR THE CI DIAGONAL ELEMENTS                                        CPP17050
      DO 115 ICI=1,NCPEX                                                CPP17060
C   FOR CLOSED---OPEN-SHELL                                             CPP17070
      DO 106 ITYP=1,NTYPES                                              CPP17080
      DO 106 JTYP=1,ITYP                                                CPP17090
      IJTYP=IOFF(MAX0(ITYP,JTYP))+MIN0(ITYP,JTYP)                       CPP17100
      IF(FOCC(ITYP).NE.TWO.AND.FOCC(ITYP).NE.ONE) GO TO 106             CPP17110
      IF(FOCC(JTYP).NE.TWO.AND.FOCC(JTYP).NE.ONE) GO TO 106             CPP17120
      AIJ(IJTYP,ICI)=ALPC(IJTYP)                                        CPP17130
      BIJ(IJTYP,ICI)=BETC(IJTYP)                                        CPP17140
  106 CONTINUE                                                          CPP17150
C   FOR CLOSED AND OPEN-SHELLS---PEX SHELLS                             CPP17160
      DO 108 ITYP=1,NTYPES                                              CPP17170
      IF(FOCC(ITYP).NE.TWO.AND.FOCC(ITYP).NE.ONE) GO TO 108             CPP17180
      DO 107 I=1,NORPEX                                                 CPP17190
      IF(INDPEX(I,ICI).EQ.0) GO TO 107                                  CPP17200
      MTYP=MOTYP(ICIMO(I))                                              CPP17210
      IMTYP=IOFF(MAX0(ITYP,MTYP))+MIN0(ITYP,MTYP)                       CPP17220
      AIJ(IMTYP,ICI)=FOCC(ITYP)                                         CPP17230
      BIJ(IMTYP,ICI)=-FOCC(ITYP)*HALF                                   CPP17240
  107 CONTINUE                                                          CPP17250
  108 CONTINUE                                                          CPP17260
C   FOR PEX-PEX SHELLS                                                  CPP17270
      DO 110 I=1,NORPEX                                                 CPP17280
      IF(INDPEX(I,ICI).EQ.0) GO TO 110                                  CPP17290
      MTYP=MOTYP(ICIMO(I))                                              CPP17300
      DO 109 J=1,I                                                      CPP17310
      IF(INDPEX(J,ICI).EQ.0) GO TO 109                                  CPP17320
      NTYP=MOTYP(ICIMO(J))                                              CPP17330
      MNTYP=IOFF(MAX0(MTYP,NTYP))+MIN0(MTYP,NTYP)                       CPP17340
      AIJ(MNTYP,ICI)=TWO                                                CPP17350
      BIJ(MNTYP,ICI)=-ONE                                               CPP17360
  109 CONTINUE                                                          CPP17370
  110 CONTINUE                                                          CPP17380
      IF(IPRNT.LE.2) GO TO 115                                          CPP17390
      WRITE(6,6) ICI                                                    CPP17400
      CALL PRINT(AIJ(1,ICI),NCTRI,NTYPEP,6)                             CPP17410
      WRITE(6,7) ICI                                                    CPP17420
      CALL PRINT(BIJ(1,ICI),NCTRI,NTYPEP,6)                             CPP17430
  115 CONTINUE                                                          CPP17440
C                                                                       CPP17450
C   TWO-ELECTRON COUPLING CONSTANTS                                     CPP17460
C   FOR THE CI OFF-DIAGONAL ELEMENTS (BIJ ONLY)                         CPP17470
      DO 116 IJCI=1,IJUNQ                                               CPP17480
      II=IIMO(IJCI)                                                     CPP17490
      JJ=JJMO(IJCI)                                                     CPP17500
      IIJJ=IOFF(MAX0(II,JJ))+MIN0(II,JJ)                                CPP17510
      IJEX=IIUNQ+IJCI                                                   CPP17520
      BIJ(IIJJ,IJEX)=HALF                                               CPP17530
      IF(IPRNT.LE.2) GO TO 116                                          CPP17540
      WRITE(6,7) IJCI                                                   CPP17550
      CALL PRINT(BIJ(1,IJEX),NCTRI,NTYPEP,6)                            CPP17560
  116 CONTINUE                                                          CPP17570
C                                                                       CPP17580
C   TEST THE "BARE" COUPLING CONSTANTS USING NON-ZERO CI'S              CPP17590
      CALL ZERO(F,20)                                                   CPP17600
      CALL ZERO(A,400)                                                  CPP17610
      CALL ZERO(B,400)                                                  CPP17620
      IJCI=0                                                            CPP17630
      DO 120 IJCI=1,ICIUNQ                                              CPP17640
      IEX=IPX(IJCI)                                                     CPP17650
      JEX=JPX(IJCI)                                                     CPP17660
      FAC=TWO                                                           CPP17670
      IF(IEX.EQ.JEX) FAC=ONE                                            CPP17680
      DO 117 ITYP=1,NTYPES                                              CPP17690
      F(ITYP)=F(ITYP)+CI(IEX)*CI(JEX)*FIJ(ITYP,IJCI)*FAC                CPP17700
  117 CONTINUE                                                          CPP17710
      DO 118 ITYP=1,NTYPES                                              CPP17720
      DO 118 JTYP=1,NTYPES                                              CPP17730
      IJTYP=IOFF(MAX0(ITYP,JTYP))+MIN0(ITYP,JTYP)                       CPP17740
      A(ITYP,JTYP)=A(ITYP,JTYP)+CI(IEX)*CI(JEX)*AIJ(IJTYP,IJCI)*FAC     CPP17750
      B(ITYP,JTYP)=B(ITYP,JTYP)+CI(IEX)*CI(JEX)*BIJ(IJTYP,IJCI)*FAC     CPP17760
  118 CONTINUE                                                          CPP17770
  120 CONTINUE                                                          CPP17780
      WRITE(6,8) (F(I),I=1,NTYPEP)                                      CPP17790
      WRITE(6,9)                                                        CPP17800
      CALL MATOUT(A,20,20,NTYPEP,NTYPEP,6)                              CPP17810
      WRITE(6,10)                                                       CPP17820
      CALL MATOUT(B,20,20,NTYPEP,NTYPEP,6)                              CPP17830
C                                                                       CPP17840
      RETURN                                                            CPP17850
      END                                                               CPP17860
      SUBROUTINE DERMAT(EAO,E1A,E2A,SS,U,T,NSECMX,ISTEP)                CPP17870
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP17880
      DIMENSION E1A(N3N),E2A(N3N,N3N),SS(NBATRI)                        CPP17890
      DIMENSION EAO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)         CPP17900
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP17910
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP17920
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP17930
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP17940
      COMMON/TAPES/ITAP37,ITAP42,ITAP43,ITAP44,ITAP47                   CPP17950
      COMMON/WORKS/JTAPE1,JTAPE2                                        CPP17960
    1 FORMAT(//,2X,' SCF FIRST DERIVATIVE MATRIX'/)                     CPP17970
    2 FORMAT(//,2X,' SCF SECOND DERIVATIVE MATRIX'/)                    CPP17980
    3 FORMAT(//,2X,' SA (AO BASIS) MATRIX, IABC = ',I5/)                CPP17990
    4 FORMAT(//,2X,' SA (MO BASIS) MATRIX, IABC = ',I5/)                CPP18000
    5 FORMAT(//,2X,' HA (AO BASIS) MATRIX, IABC = ',I5/)                CPP18010
    6 FORMAT(//,2X,' HA (MO BASIS) MATRIX, IABC = ',I5/)                CPP18020
    7 FORMAT(//,2X,' TA (AO BASIS) MATRIX, IABC = ',I5,' ITYP = ',I5/)  CPP18030
    8 FORMAT(//,2X,' TA (MO BASIS) MATRIX, IABC = ',I5,' ITYP = ',I5/)  CPP18040
C                                                                       CPP18050
      CALL RFILE(ITAP42)                                                CPP18060
      NBSET=NTREAD                                                      CPP18070
C                                                                       CPP18080
C   READ IN SCF FIRST DERIVATIVES                                       CPP18090
      CALL SREAD(ITAP42,E1A,N3N*2)                                      CPP18100
      WRITE(6,1)                                                        CPP18110
      CALL MATOUT(E1A,3,NATOMS,3,NATOMS,6)                              CPP18120
C                                                                       CPP18130
C   READ IN SCF SECOND DERIVATIVES                                      CPP18140
      CALL SREAD(ITAP42,E2A,N3N*N3N*2)                                  CPP18150
      WRITE(6,2)                                                        CPP18160
      CALL MATOUT(E2A,N3N,N3N,N3N,N3N,6)                                CPP18170
      IF(ISTEP.GE.5) GO TO 205                                          CPP18180
C                                                                       CPP18190
C   READ IN S FIRST DERIVATIVE INTEGRALS                                CPP18200
C   AND TRANSFORM THEM FROM AO TO MO BASIS                              CPP18210
      CALL RFILE(ITAP44)                                                CPP18220
      CALL RFILE(JTAPE1)                                                CPP18230
      CALL ZFILE(ITAP44,NSECMX)                                         CPP18240
      DO 101 IABC=1,N3N                                                 CPP18250
      IS=ISA(IABC)                                                      CPP18260
      CALL SREAD(ITAP42,SS,NTRI2)                                       CPP18270
      IF(IPRNT.LE.4) GO TO 201                                          CPP18280
      WRITE(6,3) IABC                                                   CPP18290
      CALL PRINT(SS,NBATRI,NBASIS,6)                                    CPP18300
  201 CONTINUE                                                          CPP18310
      CALL MOCONV(SS,NBATRI,SS,NBATRI,EAO,U,T)                          CPP18320
      CALL RWRIT(ITAP44,SS,NTRI2,IS)                                    CPP18330
      IF(IPRNT.LE.4) GO TO 101                                          CPP18340
      WRITE(6,4) IABC                                                   CPP18350
      CALL PRINT(SS,NBATRI,NBASIS,6)                                    CPP18360
  101 CONTINUE                                                          CPP18370
C                                                                       CPP18380
C   READ IN ONE ELECTRON FIRST DERIVATIVE INTEGRALS                     CPP18390
C   AND TRANSFORM THEM FROM AO TO MO BASIS                              CPP18400
      DO 102 IABC=1,N3N                                                 CPP18410
      IH=IHA(IABC)                                                      CPP18420
      CALL SREAD(ITAP42,SS,NTRI2)                                       CPP18430
      IF(IPRNT.LE.4) GO TO 202                                          CPP18440
      WRITE(6,5) IABC                                                   CPP18450
      CALL PRINT(SS,NBATRI,NBFAO,6)                                     CPP18460
  202 CONTINUE                                                          CPP18470
      CALL MOCONV(SS,NBATRI,SS,NBATRI,EAO,U,T)                          CPP18480
      CALL RWRIT(ITAP44,SS,NTRI2,IH)                                    CPP18490
      IF(IPRNT.LE.4) GO TO 102                                          CPP18500
      WRITE(6,6) IABC                                                   CPP18510
      CALL PRINT(SS,NBATRI,NBASIS,6)                                    CPP18520
  102 CONTINUE                                                          CPP18530
C                                                                       CPP18540
C   READ IN TWO ELECTRON FIRST DERIVATIVE INTEGRALS                     CPP18550
C   AND TRANSFORM THEM FROM AO TO MO BASIS                              CPP18560
      DO 104 ITYP=1,NBSET                                               CPP18570
      DO 103 IABC=1,N3N                                                 CPP18580
      CALL SREAD(ITAP42,SS,NTRI2)                                       CPP18590
      IF(IPRNT.LE.4) GO TO 203                                          CPP18600
      WRITE(6,7) ITYP,IABC                                              CPP18610
      CALL PRINT(SS,NBATRI,NBFAO,6)                                     CPP18620
  203 CALL MOCONV(SS,NBATRI,SS,NBATRI,EAO,U,T)                          CPP18630
      CALL SWRIT(JTAPE1,SS,NTRI2)                                       CPP18640
      IF(IPRNT.LE.4) GO TO 103                                          CPP18650
      WRITE(6,8) ITYP,IABC                                              CPP18660
      CALL PRINT(SS,NBATRI,NBASIS,6)                                    CPP18670
  103 CONTINUE                                                          CPP18680
  104 CONTINUE                                                          CPP18690
      CALL RCLOSE(ITAP44,3)                                             CPP18700
      CALL RCLOSE(JTAPE1,3)                                             CPP18710
C                                                                       CPP18720
  205 CONTINUE                                                          CPP18730
      CALL RCLOSE(ITAP42,3)                                             CPP18740
      RETURN                                                            CPP18750
      END                                                               CPP18760
      SUBROUTINE MOCONV(SA,NNA,SM,NNM,EAO,U,T)                          CPP18770
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP18780
      DIMENSION SA(NNA),SM(NNM)                                         CPP18790
      DIMENSION EAO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)         CPP18800
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP18810
C                                                                       CPP18820
C   TRANSFORM INTEGRALS FROM AO TO MO BASIS                             CPP18830
      NBFASQ=NBFAO*NBFAO                                                CPP18840
      IIJJ=0                                                            CPP18850
      DO 101 II=1,NBFAO                                                 CPP18860
      DO 101 JJ=1,II                                                    CPP18870
      IIJJ=IIJJ+1                                                       CPP18880
      T(II,JJ)=SA(IIJJ)                                                 CPP18890
      T(JJ,II)=SA(IIJJ)                                                 CPP18900
  101 CONTINUE                                                          CPP18910
      CALL ZERO(U,NBFASQ)                                               CPP18920
      CALL MXMB(EAO,NBFAO,1,T,1,NBFAO,U,1,NBFAO,NBFAO,NBFAO,NBFAO)      CPP18930
      CALL ZERO(T,NBFASQ)                                               CPP18940
      CALL MXMB(U,1,NBFAO,EAO,1,NBFAO,T,1,NBFAO,NBFAO,NBFAO,NBFAO)      CPP18950
      IJ=0                                                              CPP18960
      DO 102 I=1,NBASIS                                                 CPP18970
      DO 102 J=1,I                                                      CPP18980
      IJ=IJ+1                                                           CPP18990
      SM(IJ)=T(I,J)                                                     CPP19000
  102 CONTINUE                                                          CPP19010
C                                                                       CPP19020
      RETURN                                                            CPP19030
      END                                                               CPP19040
      SUBROUTINE MOCONS(SA,NNA,SM,NNM,ESO,U,T)                          CPP19050
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP19060
      DIMENSION SA(NNA),SM(NNM)                                         CPP19070
      DIMENSION ESO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)         CPP19080
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP19090
C                                                                       CPP19100
C   TRANSFORM INTEGRALS FROM SO TO MO BASIS                             CPP19110
      NBFSSQ=NBASIS*NBASIS                                              CPP19120
      IIJJ=0                                                            CPP19130
      DO 101 II=1,NBASIS                                                CPP19140
      DO 101 JJ=1,II                                                    CPP19150
      IIJJ=IIJJ+1                                                       CPP19160
      T(II,JJ)=SA(IIJJ)                                                 CPP19170
      IF(II.EQ.JJ) GO TO 101                                            CPP19180
      T(JJ,II)=SA(IIJJ)                                                 CPP19190
  101 CONTINUE                                                          CPP19200
      CALL ZERO(U,NBFSSQ)                                               CPP19210
      CALL MXMB(ESO,NBFAO,1,T,1,NBFAO,U,1,NBFAO,NBFAO,NBFAO,NBFAO)      CPP19220
      CALL ZERO(T,NBFSSQ)                                               CPP19230
      CALL MXMB(U,1,NBFAO,ESO,1,NBFAO,T,1,NBFAO,NBFAO,NBFAO,NBFAO)      CPP19240
      IJ=0                                                              CPP19250
      DO 102 I=1,NBASIS                                                 CPP19260
      DO 102 J=1,I                                                      CPP19270
      IJ=IJ+1                                                           CPP19280
      SM(IJ)=T(I,J)                                                     CPP19290
  102 CONTINUE                                                          CPP19300
C                                                                       CPP19310
      RETURN                                                            CPP19320
      END                                                               CPP19330
      SUBROUTINE E0CAL(ESO,HSO,HMO,TM,U,T,DP,DQ,DM,MIJ,PQ)              CPP19340
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP19350
      INTEGER*2 MIJ(*)                                                  CPP19360
      DIMENSION HSO(NTRI),HMO(NTRI),TM(NTRI)                            CPP19370
      DIMENSION ESO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)         CPP19380
      DIMENSION DP(NTRI,NTYPES),DQ(NTRI,NTYPES),DM(NTRI,NTYPES)         CPP19390
      DIMENSION PQ(*)                                                   CPP19400
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP19410
      COMMON/ENRGY/ESCF,ENUC,ELEC                                       CPP19420
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP19430
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP19440
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP19450
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP19460
      DATA A00 / 0.0D+00 /                                              CPP19470
    1 FORMAT(//,2X,' H (SO BASIS) MATRIX'/)                             CPP19480
    2 FORMAT(//,2X,' H (MO BASIS) MATRIX'/)                             CPP19490
    3 FORMAT(//,2X,' SCF ENERGIES'/                                     CPP19500
     * 2X,' ESCF       = ',F20.10/                                      CPP19510
     * 2X,' ENUC       = ',F20.10/                                      CPP19520
     * 2X,' ELEC-ONE   = ',F20.10/                                      CPP19530
     * 2X,' ELEC-TWO   = ',F20.10/                                      CPP19540
     * 2X,' ELEC-TOTAL = ',F20.10/                                      CPP19550
     * 2X,' ETOT       = ',F20.10)                                      CPP19560
C                                                                       CPP19570
C   CALCULATE SCF ENERGY FOR A TEST                                     CPP19580
C   SEE EQ.(4)                                                          CPP19590
C                                                                       CPP19600
C   READ IN ONE ELECTRON INTEGRALS (SO BASIS)                           CPP19610
      CALL RMREAD(HSO,16)                                               CPP19620
      CALL MOCONS(HSO,NTRI,HMO,NTRI,ESO,U,T)                            CPP19630
      IF(IPRNT.LE.3) GO TO 201                                          CPP19640
      WRITE(6,1)                                                        CPP19650
      CALL PRINT(HSO,NTRI,NBASIS,6)                                     CPP19660
      WRITE(6,2)                                                        CPP19670
      CALL PRINT(HMO,NTRI,NBASIS,6)                                     CPP19680
C                                                                       CPP19690
C   CALCULATE ONE-ELCTRON ENERGY                                        CPP19700
  201 CONTINUE                                                          CPP19710
      ELEC1=A00                                                         CPP19720
      DO 101 I=1,JOCC                                                   CPP19730
      II=IOFF(I+1)                                                      CPP19740
      ITYP=MOTYP(I)                                                     CPP19750
      ELEC1=ELEC1+HMO(II)*FOCC(ITYP)                                    CPP19760
  101 CONTINUE                                                          CPP19770
C                                                                       CPP19780
C   FORM DENSITY MATRIX IN SO BASIS                                     CPP19790
      CALL PQAMAT(DP,DQ,ESO,NTYPES)                                     CPP19800
C                                                                       CPP19810
C   FORM HALF-TRANSFORMED TWO ELECTRON MATRIX                           CPP19820
      CALL PQMAT(DP,DQ,DM,NTYPES,MIJ,PQ)                                CPP19830
C                                                                       CPP19840
C   CALCULATE TWO-ELECTRON ENERGY                                       CPP19850
      ELEC2=A00                                                         CPP19860
      DO 103 ITYP=1,NTYPES                                              CPP19870
      CALL MOCONS(DM(1,ITYP),NTRI,TM,NTRI,ESO,U,T)                      CPP19880
      DO 102 I=NSTR(ITYP),NEND(ITYP)                                    CPP19890
      II=IOFF(I+1)                                                      CPP19900
      ELEC2=ELEC2+TM(II)                                                CPP19910
  102 CONTINUE                                                          CPP19920
  103 CONTINUE                                                          CPP19930
C                                                                       CPP19940
      ELECT=ELEC1+ELEC2                                                 CPP19950
      ETOT=ENUC+ELECT                                                   CPP19960
      WRITE(6,3) ESCF,ENUC,ELEC1,ELEC2,ELECT,ETOT                       CPP19970
C                                                                       CPP19980
      RETURN                                                            CPP19990
      END                                                               CPP20000
      SUBROUTINE PQAMAT(DP,DQ,ESO,NABC)                                 CPP20010
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP20020
      DIMENSION DP(NTRI,NABC),DQ(NTRI,NABC),ESO(NBFAO,NBASIS)           CPP20030
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP20040
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP20050
      COMMON/GVBAB/A(30,30),B(30,30)                                    CPP20060
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP20070
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP20080
      COMMON/TOLER/CRIT,TOLR                                            CPP20090
      DATA A00,TWO / 0.0D+00 , 2.0D+00 /                                CPP20100
    1 FORMAT(//,2X,' DP MATRIX'/)                                       CPP20110
    2 FORMAT(//,2X,' DQ MATRIX'/)                                       CPP20120
C                                                                       CPP20130
      IJ=0                                                              CPP20140
      DO 105 I=1,NBASIS                                                 CPP20150
      DO 105 J=1,I                                                      CPP20160
      IJ=IJ+1                                                           CPP20170
      DO 104 ITYP=1,NTYPES                                              CPP20180
      SUMP=A00                                                          CPP20190
      SUMQ=A00                                                          CPP20200
      DO 103 JTYP=1,NTYPES                                              CPP20210
      DO 102 K=NSTR(JTYP),NEND(JTYP)                                    CPP20220
      FAC=ESO(I,K)*ESO(J,K)                                             CPP20230
      IF(DABS(FAC).LT.TOLR) GO TO 102                                   CPP20240
      SUMP=SUMP+A(ITYP,JTYP)*FAC                                        CPP20250
      SUMQ=SUMQ+B(ITYP,JTYP)*FAC                                        CPP20260
  102 CONTINUE                                                          CPP20270
  103 CONTINUE                                                          CPP20280
      DP(IJ,ITYP)=SUMP                                                  CPP20290
      DQ(IJ,ITYP)=SUMQ                                                  CPP20300
  104 CONTINUE                                                          CPP20310
  105 CONTINUE                                                          CPP20320
C                                                                       CPP20330
      IF(IPRNT.LE.4) GO TO 201                                          CPP20340
      WRITE(6,1)                                                        CPP20350
      CALL MATOUT(DP,NTRI,NTYPES,NTRI,NTYPES,6)                         CPP20360
      WRITE(6,2)                                                        CPP20370
      CALL MATOUT(DQ,NTRI,NTYPES,NTRI,NTYPES,6)                         CPP20380
C                                                                       CPP20390
  201 CONTINUE                                                          CPP20400
      DO 107 I=1,NBASIS                                                 CPP20410
      II=IOFF(I+1)                                                      CPP20420
      DO 106 ITYP=1,NTYPES                                              CPP20430
      DP(II,ITYP)=DP(II,ITYP)/TWO                                       CPP20440
      DQ(II,ITYP)=DQ(II,ITYP)/TWO                                       CPP20450
  106 CONTINUE                                                          CPP20460
  107 CONTINUE                                                          CPP20470
      CALL SCLMPY(DP,NTRI*NTYPES,TWO)                                   CPP20480
      CALL SCLMPY(DQ,NTRI*NTYPES,TWO)                                   CPP20490
C                                                                       CPP20500
      RETURN                                                            CPP20510
      END                                                               CPP20520
      SUBROUTINE PQMAT(DP,DQ,DM,NABC,MIJ,PQ)                            CPP20530
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP20540
      INTEGER*2 MIJ(*)                                                  CPP20550
      DIMENSION DP(NTRI,NABC),DQ(NTRI,NABC),DM(NTRI,NABC)               CPP20560
      DIMENSION PQ(*)                                                   CPP20570
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP20580
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPP20590
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP20600
      COMMON/SUPER/MAXIJ                                                CPP20610
      COMMON/TAPES/ITAP37,ITAP42,ITAP43,ITAP44,ITAP47                   CPP20620
    1 FORMAT(2X,3I5,4F15.8)                                             CPP20630
    2 FORMAT(//,2X,' DM MATRIX'/)                                       CPP20640
C                                                                       CPP20650
C   DP :  COULOMB DENSITY TYPE SO MATRIX                                CPP20660
C   DQ :  EXCHANGE DENSITY TYPE SO MATRIX                               CPP20670
C   DM :  HALF TRANSFORMED SO MATRIX                                    CPP20680
C                                                                       CPP20690
C   FORM DM MATRICES                                                    CPP20700
      CALL ZERO(DM,NTRI*NABC)                                           CPP20710
C                                                                       CPP20720
      CALL RFILE(ITAP37)                                                CPP20730
  200 CALL SREAD(ITAP37,MIJ,MAXBF4)                                     CPP20740
      LIMINT=MIJ(10*MAXIJ+1)                                            CPP20750
      DO 102 M=1,LIMINT                                                 CPP20760
      IJ=MIJ(M)                                                         CPP20770
      KL=MIJ(M+MAXIJ)                                                   CPP20780
      P=PQ(M)                                                           CPP20790
      Q=PQ(M+MAXIJ)                                                     CPP20800
      CJ=P+Q                                                            CPP20810
      EK=Q+Q                                                            CPP20820
C     IF(IPRNT.LE.5) GO TO 201                                          CPP20830
C     WRITE(6,1) M,IJ,KL,P,Q,CJ,EK                                      CPP20840
  201 CONTINUE                                                          CPP20850
      DO 101 IABC=1,NABC                                                CPP20860
      DM(IJ,IABC)=DM(IJ,IABC)+DP(KL,IABC)*CJ+DQ(KL,IABC)*EK             CPP20870
      DM(KL,IABC)=DM(KL,IABC)+DP(IJ,IABC)*CJ+DQ(IJ,IABC)*EK             CPP20880
  101 CONTINUE                                                          CPP20890
  102 CONTINUE                                                          CPP20900
      IF(LIMINT.NE.MAXIJ) GO TO 202                                     CPP20910
      GO TO 200                                                         CPP20920
C                                                                       CPP20930
  202 CONTINUE                                                          CPP20940
      CALL SREW(ITAP37)                                                 CPP20950
      IF (IPRNT.LE.4) GO TO 203                                         CPP20960
      WRITE(6,2)                                                        CPP20970
      CALL MATOUT(DM,NTRI,NABC,NTRI,NABC,6)                             CPP20980
C                                                                       CPP20990
  203 CONTINUE                                                          CPP21000
      CALL RCLOSE(ITAP37,3)                                             CPP21010
      RETURN                                                            CPP21020
      END                                                               CPP21030
      SUBROUTINE PQQMAT(DQ,DM,NABC,MIJ,PQ)                              CPP21040
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP21050
      INTEGER*2 MIJ(*)                                                  CPP21060
      DIMENSION DQ(NTRI,NABC),DM(NTRI,NABC)                             CPP21070
      DIMENSION PQ(*)                                                   CPP21080
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP21090
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPP21100
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP21110
      COMMON/SUPER/MAXIJ                                                CPP21120
      COMMON/TAPES/ITAP37,ITAP42,ITAP43,ITAP44,ITAP47                   CPP21130
    1 FORMAT(2X,3I5,4F15.8)                                             CPP21140
    2 FORMAT(//,2X,' DM MATRIX'/)                                       CPP21150
C                                                                       CPP21160
C   DQ :  EXCHANGE DENSITY TYPE SO MATRIX                               CPP21170
C   DM :  HALF TRANSFORMED SO MATRIX                                    CPP21180
C                                                                       CPP21190
C   FORM DM MATRICES                                                    CPP21200
      CALL ZERO(DM,NTRI*NABC)                                           CPP21210
C                                                                       CPP21220
      CALL RFILE(ITAP37)                                                CPP21230
  200 CALL SREAD(ITAP37,MIJ,MAXBF4)                                     CPP21240
      LIMINT=MIJ(10*MAXIJ+1)                                            CPP21250
      DO 102 M=1,LIMINT                                                 CPP21260
      IJ=MIJ(M)                                                         CPP21270
      KL=MIJ(M+MAXIJ)                                                   CPP21280
C*    P=PQ(M)                                                           CPP21290
      Q=PQ(M+MAXIJ)                                                     CPP21300
C*    CJ=P+Q                                                            CPP21310
      EK=Q+Q                                                            CPP21320
C     IF(IPRNT.LE.5) GO TO 201                                          CPP21330
C     WRITE(6,1) M,IJ,KL,P,Q,CJ,EK                                      CPP21340
  201 CONTINUE                                                          CPP21350
      DO 101 IABC=1,NABC                                                CPP21360
      DM(IJ,IABC)=DM(IJ,IABC)+DQ(KL,IABC)*EK                            CPP21370
      DM(KL,IABC)=DM(KL,IABC)+DQ(IJ,IABC)*EK                            CPP21380
  101 CONTINUE                                                          CPP21390
  102 CONTINUE                                                          CPP21400
      IF(LIMINT.NE.MAXIJ) GO TO 202                                     CPP21410
      GO TO 200                                                         CPP21420
C                                                                       CPP21430
  202 CONTINUE                                                          CPP21440
      CALL SREW(ITAP37)                                                 CPP21450
      IF (IPRNT.LE.4) GO TO 203                                         CPP21460
      WRITE(6,2)                                                        CPP21470
      CALL MATOUT(DM,NTRI,NABC,NTRI,NABC,6)                             CPP21480
C                                                                       CPP21490
  203 CONTINUE                                                          CPP21500
      CALL RCLOSE(ITAP37,3)                                             CPP21510
      RETURN                                                            CPP21520
      END                                                               CPP21530
      SUBROUTINE EADER(D1W,E1T,EPS,D1O,D1T,SM,HM,TM,NABC)               CPP21540
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP21550
      DIMENSION D1W(N3N),E1T(N3N),EPS(NTRI)                             CPP21560
      DIMENSION D1O(N3N),D1T(N3N)                                       CPP21570
      DIMENSION SM(NTRI),HM(NTRI),TM(NTRI,NABC)                         CPP21580
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP21590
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP21600
      COMMON/GVBAB/ALPA(30,30),BETA(30,30)                              CPP21610
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP21620
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP21630
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP21640
      COMMON/TAPES/ITAP37,ITAP42,ITAP43,ITAP44,ITAP47                   CPP21650
      COMMON/WORKS/JTAPE1,JTAPE2                                        CPP21660
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP21670
      DATA A00,ONE,TWO / 0.0D+00 , 1.0D+00 , 2.0D+00 /                  CPP21680
    1 FORMAT(//,2X,' LAGRANGIAN MATRIX'/)                               CPP21690
    2 FORMAT(//,2X,' D1O MATRIX IN EADER'/)                             CPP21700
    3 FORMAT(//,2X,' D1T MATRIX IN EADER'/)                             CPP21710
    4 FORMAT(//,2X,' D1W MATRIX IN EADER'/)                             CPP21720
    5 FORMAT(//,2X,' EGRAD MATRIX IN EADER'/)                           CPP21730
C                                                                       CPP21740
C   CALCULATE SCF FIRST DERIVATIVES FOR A TEST                          CPP21750
C   SEE EQ.(9)                                                          CPP21760
C                                                                       CPP21770
      CALL RMREAD(EPS,27)                                               CPP21780
      IF(IPRNT.LE.2) GO TO 201                                          CPP21790
      WRITE(6,1)                                                        CPP21800
      CALL PRINT(EPS,NTRI,NBASIS,6)                                     CPP21810
C                                                                       CPP21820
C   CALL BACK DERIVATIVE TM MATRICES                                    CPP21830
  201 CONTINUE                                                          CPP21840
      CALL RFILE(JTAPE1)                                                CPP21850
      DO 101 ITYPX=1,NABC                                               CPP21860
      CALL SREAD(JTAPE1,TM(1,ITYPX),NTRI2)                              CPP21870
  101 CONTINUE                                                          CPP21880
C                                                                       CPP21890
      CALL RFILE(ITAP44)                                                CPP21900
      DO 110 IABC=1,N3N                                                 CPP21910
      IS=ISA(IABC)                                                      CPP21920
      IH=IHA(IABC)                                                      CPP21930
      CALL RREAD(ITAP44,SM,NTRI2,IS)                                    CPP21940
      CALL RREAD(ITAP44,HM,NTRI2,IH)                                    CPP21950
C  ONE-ELECTRON CONTRIBUTION                                            CPP21960
      VALU=A00                                                          CPP21970
      DO 106 I=1,JOCC                                                   CPP21980
      II=IOFF(I+1)                                                      CPP21990
      ITYP=MOTYP(I)                                                     CPP22000
      VALU=VALU+HM(II)*FOCC(ITYP)                                       CPP22010
  106 CONTINUE                                                          CPP22020
      D1O(IABC)=VALU                                                    CPP22030
C  TWO-ELECTRON CONTRIBUTION                                            CPP22040
      VALU=A00                                                          CPP22050
      DO 108 ITYP=1,NTYPES                                              CPP22060
      JAPT=(ITYP*2-2)*N3N+IABC                                          CPP22070
      KAPT=(ITYP*2-1)*N3N+IABC                                          CPP22080
      DO 108 JTYP=1,NTYPES                                              CPP22090
      DO 107 J=NSTR(JTYP),NEND(JTYP)                                    CPP22100
      JJ=IOFF(J+1)                                                      CPP22110
      VALU=VALU+ALPA(ITYP,JTYP)*TM(JJ,JAPT)                             CPP22120
     1         +BETA(ITYP,JTYP)*TM(JJ,KAPT)                             CPP22130
  107 CONTINUE                                                          CPP22140
  108 CONTINUE                                                          CPP22150
      D1T(IABC)=VALU                                                    CPP22160
C   OVERLAP CONTRIBUTION                                                CPP22170
      VALU=A00                                                          CPP22180
      IJ=0                                                              CPP22190
      DO 109 I=1,JOCC                                                   CPP22200
      DO 109 J=1,I                                                      CPP22210
      IJ=IJ+1                                                           CPP22220
      FAC=TWO                                                           CPP22230
      IF(I.EQ.J) FAC=ONE                                                CPP22240
      VALU=VALU-SM(IJ)*EPS(IJ)*FAC                                      CPP22250
  109 CONTINUE                                                          CPP22260
      D1W(IABC)=VALU*TWO                                                CPP22270
C                                                                       CPP22280
C   SUM UP ALL CONTRIBUTIONS                                            CPP22290
      E1T(IABC)=D1O(IABC)+D1T(IABC)+D1W(IABC)                           CPP22300
C                                                                       CPP22310
  110 CONTINUE                                                          CPP22320
C                                                                       CPP22330
      WRITE(6,2)                                                        CPP22340
      CALL MATOUT(D1O,3,NATOMS,3,NATOMS,6)                              CPP22350
      WRITE(6,3)                                                        CPP22360
      CALL MATOUT(D1T,3,NATOMS,3,NATOMS,6)                              CPP22370
      WRITE(6,4)                                                        CPP22380
      CALL MATOUT(D1W,3,NATOMS,3,NATOMS,6)                              CPP22390
      WRITE(6,5)                                                        CPP22400
      CALL MATOUT(E1T,3,NATOMS,3,NATOMS,6)                              CPP22410
C                                                                       CPP22420
      CALL RCLOSE(ITAP44,3)                                             CPP22430
      CALL RCLOSE(JTAPE1,3)                                             CPP22440
      RETURN                                                            CPP22450
      END                                                               CPP22460
      SUBROUTINE HAMAT(FIJ,AIJ,BIJ,D1W,E1T,HCIA,D1E,HM,TM,NABC,N3NXX)   CPP22470
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP22480
      DIMENSION FIJ(NTYPEP,ICIUNQ),AIJ(NTYTRI,ICIUNQ),BIJ(NTYTRI,ICIUNQ)CPP22490
      DIMENSION HCIA(NCTRI,N3NXX)                                       CPP22500
      DIMENSION D1W(N3N),E1T(N3N),D1E(N3N)                              CPP22510
      DIMENSION HM(NTRI,N3N),TM(NTRI,NABC)                              CPP22520
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP22530
      COMMON/CIOFF/IIUNQ,IJUNQ,ICIUNQ,IPX(5050),JPX(5050)               CPP22540
      COMMON/CIPOS/ICITYP(5050),ICPOS(200,100)                          CPP22550
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP22560
      COMMON/CITYP/IJGRP,IBM(200),JBM(200),IBNUM(200)                   CPP22570
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP22580
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP22590
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP22600
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP22610
      COMMON/TAPES/ITAP37,ITAP42,ITAP43,ITAP44,ITAP47                   CPP22620
      COMMON/WORKS/JTAPE1,JTAPE2                                        CPP22630
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP22640
      COMMON/CI103/CI(100)                                              CPP22650
      DATA A00,HALF,ONE,TWO,FOUR / 0.0D+00 , 0.5D+00 , 1.0D+00 ,        CPP22660
     1                             2.0D+00 , 4.0D+00 /                  CPP22670
    1 FORMAT(//,2X,' HA MATRIX'/)                                       CPP22680
    2 FORMAT(2X,' IJCI,IABC,ELEC1,ELEC2 = ',2I5,2F15.8)                 CPP22690
    3 FORMAT(//,2X,' HCIA MATRIX'/)                                     CPP22700
    4 FORMAT(//,2X,' D1E MATRIX IN HAMAT'/)                             CPP22710
    5 FORMAT(//,2X,' D1W MATRIX IN HAMAT'/)                             CPP22720
    6 FORMAT(//,2X,' EGRAD MATRIX IN HAMAT'/)                           CPP22730
C                                                                       CPP22740
C   CALL BACK HM AND TM MATRICES AND FORM HA MATRICES                   CPP22750
C   SEE EQS.(33)-(35)                                                   CPP22760
C                                                                       CPP22770
      CALL RFILE(ITAP44)                                                CPP22780
      DO 101 IABC=1,N3N                                                 CPP22790
      IH=IHA(IABC)                                                      CPP22800
      CALL RREAD(ITAP44,HM(1,IABC),NTRI2,IH)                            CPP22810
  101 CONTINUE                                                          CPP22820
C                                                                       CPP22830
      CALL RFILE(JTAPE1)                                                CPP22840
      DO 102 ITYPX=1,NABC                                               CPP22850
      CALL SREAD(JTAPE1,TM(1,ITYPX),NTRI2)                              CPP22860
  102 CONTINUE                                                          CPP22870
C                                                                       CPP22880
      CALL ZERO(HCIA,NCTRI*N3NXX)                                       CPP22890
C                                                                       CPP22900
C     IF(IPRNT.GT.2)                                                    CPP22910
C     WRITE(6,1)                                                        CPP22920
C   FOR CI DIAGONAL ELEMENTS                                            CPP22930
      DO 110 IJCI=1,IIUNQ                                               CPP22940
      IEX=IPX(IJCI)                                                     CPP22950
      JEX=JPX(IJCI)                                                     CPP22960
      IJEX=IOFF(MAX0(IEX,JEX))+MIN0(IEX,JEX)                            CPP22970
      DO 109 IABC=1,N3N                                                 CPP22980
C   ONE-ELECTRON CONTRIBUTION                                           CPP22990
      ELEC1=A00                                                         CPP23000
      DO 103 I=1,JOCC                                                   CPP23010
      II=IOFF(I+1)                                                      CPP23020
      ITYP=MOTYP(I)                                                     CPP23030
      ELEC1=ELEC1+FIJ(ITYP,IJCI)*HM(II,IABC)*TWO                        CPP23040
  103 CONTINUE                                                          CPP23050
C   TWO-ELECTRON CONTRIBUTION                                           CPP23060
      ELEC2=A00                                                         CPP23070
      DO 105 ITYP=1,NTYPES                                              CPP23080
      JAPT=(ITYP*2-2)*N3N+IABC                                          CPP23090
      KAPT=(ITYP*2-1)*N3N+IABC                                          CPP23100
      DO 105 JTYP=1,NTYPES                                              CPP23110
      IJTYP=IOFF(MAX0(ITYP,JTYP))+MIN0(ITYP,JTYP)                       CPP23120
      DO 104 J=NSTR(JTYP),NEND(JTYP)                                    CPP23130
      JJ=IOFF(J+1)                                                      CPP23140
      ELEC2=ELEC2+AIJ(IJTYP,IJCI)*TM(JJ,JAPT)                           CPP23150
     1           +BIJ(IJTYP,IJCI)*TM(JJ,KAPT)                           CPP23160
  104 CONTINUE                                                          CPP23170
  105 CONTINUE                                                          CPP23180
C                                                                       CPP23190
      HCIA(IJEX,IABC)=ELEC1+ELEC2                                       CPP23200
      IF(IPRNT.GT.2)                                                    CPP23210
     *WRITE(6,2) IJCI,IABC,ELEC1,ELEC2                                  CPP23220
  109 CONTINUE                                                          CPP23230
  110 CONTINUE                                                          CPP23240
C                                                                       CPP23250
C   FOR CI OFF-DIAGONAL ELEMENTS                                        CPP23260
      DO 115 IBCI=1,IJGRP                                               CPP23270
      ICM=IBM(IBCI)                                                     CPP23280
      JCM=JBM(IBCI)                                                     CPP23290
      NUM=IBNUM(IBCI)                                                   CPP23300
      ITYP=MOTYP(ICM)                                                   CPP23310
      JJ=IOFF(JCM+1)                                                    CPP23320
      DO 112 IABC=1,N3N                                                 CPP23330
      KAPT=(ITYP*2-1)*N3N+IABC                                          CPP23340
      ELEC2=TM(JJ,KAPT)                                                 CPP23350
      DO 111 INUM=1,NUM                                                 CPP23360
      IJCI=ICPOS(IBCI,INUM)                                             CPP23370
      IEX=IPX(IIUNQ+IJCI)                                               CPP23380
      JEX=JPX(IIUNQ+IJCI)                                               CPP23390
      IJEX=IOFF(MAX0(IEX,JEX))+MIN0(IEX,JEX)                            CPP23400
      HCIA(IJEX,IABC)=ELEC2                                             CPP23410
  111 CONTINUE                                                          CPP23420
  112 CONTINUE                                                          CPP23430
  115 CONTINUE                                                          CPP23440
      IF(IPRNT.LE.2) GO TO 201                                          CPP23450
      WRITE(6,3)                                                        CPP23460
      CALL MATOUT(HCIA,NCTRI,N3NXX,NCTRI,N3N,6)                         CPP23470
C                                                                       CPP23480
C   CALCULATE SCF DERIVATIVES FOR A TEST                                CPP23490
  201 CONTINUE                                                          CPP23500
      CALL ZERO(D1E,N3N)                                                CPP23510
      DO 120 IABC=1,N3N                                                 CPP23520
      DO 119 IJCI=1,ICIUNQ                                              CPP23530
      IEX=IPX(IJCI)                                                     CPP23540
      JEX=JPX(IJCI)                                                     CPP23550
      IJEX=IOFF(MAX0(IEX,JEX))+MIN0(IEX,JEX)                            CPP23560
      FAC=TWO                                                           CPP23570
      IF(IEX.EQ.JEX) FAC=ONE                                            CPP23580
      CIFAC=CI(IEX)*CI(JEX)*FAC                                         CPP23590
      D1E(IABC)=D1E(IABC)+HCIA(IJEX,IABC)*CIFAC                         CPP23600
  119 CONTINUE                                                          CPP23610
      E1T(IABC)=D1E(IABC)+D1W(IABC)                                     CPP23620
  120 CONTINUE                                                          CPP23630
C***  IF(IPRNT.LE.2) GO TO 202                                          CPP23640
      WRITE(6,4)                                                        CPP23650
      CALL MATOUT(D1E,3,NATOMS,3,NATOMS,6)                              CPP23660
      WRITE(6,5)                                                        CPP23670
      CALL MATOUT(D1W,3,NATOMS,3,NATOMS,6)                              CPP23680
      WRITE(6,6)                                                        CPP23690
      CALL MATOUT(E1T,3,NATOMS,3,NATOMS,6)                              CPP23700
C                                                                       CPP23710
  202 CONTINUE                                                          CPP23720
      CALL RCLOSE(ITAP44,3)                                             CPP23730
      CALL RCLOSE(JTAPE1,3)                                             CPP23740
      RETURN                                                            CPP23750
      END                                                               CPP23760
      SUBROUTINE FAMAT(FF,EPA,HM,TM,NABC)                               CPP23770
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP23780
      DIMENSION FF(NTRI),HM(NTRI,N3N),TM(NTRI,NABC)                     CPP23790
      DIMENSION EPA(NBASIS,NBASIS)                                      CPP23800
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP23810
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP23820
      COMMON/GVBAB/ALPA(30,30),BETA(30,30)                              CPP23830
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP23840
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP23850
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP23860
      COMMON/TAPES/ITAP37,ITAP42,ITAP43,ITAP44,ITAP47                   CPP23870
      COMMON/WORKS/JTAPE1,JTAPE2                                        CPP23880
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP23890
      DATA HALF,TWO / 0.5D+00  , 2.0D+00 /                              CPP23900
    1 FORMAT(//,2X,' EPA MATRIX, IABC = ',I5/)                          CPP23910
C                                                                       CPP23920
C   CALL BACK HM AND TM MATRICES AND FORM FA MATRICES                   CPP23930
C   SEE EQS.(15) AND (32)                                               CPP23940
C                                                                       CPP23950
      CALL RFILE(ITAP44)                                                CPP23960
      CALL RFILE(JTAPE1)                                                CPP23970
C                                                                       CPP23980
      DO 101 IABC=1,N3N                                                 CPP23990
      IH=IHA(IABC)                                                      CPP24000
      CALL RREAD(ITAP44,HM(1,IABC),NTRI2,IH)                            CPP24010
  101 CONTINUE                                                          CPP24020
C                                                                       CPP24030
      DO 102 ITYPX=1,NABC                                               CPP24040
      CALL SREAD(JTAPE1,TM(1,ITYPX),NTRI2)                              CPP24050
  102 CONTINUE                                                          CPP24060
C                                                                       CPP24070
C   STORE DERIVATIVE LAGRANGIAN MATRICES                                CPP24080
      DO 110 IABC=1,N3N                                                 CPP24090
      CALL ZERO(EPA,NBASIS*NBASIS)                                      CPP24100
      IE=IEA(IABC)                                                      CPP24110
      DO 105 I=1,JOCC                                                   CPP24120
      ITYP=MOTYP(I)                                                     CPP24130
      FAC=FOCC(ITYP)*HALF                                               CPP24140
      DO 105 J=1,NBASIS                                                 CPP24150
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPP24160
      VALU=HM(IJ,IABC)*FAC                                              CPP24170
      DO 104 KTYP=1,NTYPES                                              CPP24180
      JAPT=(KTYP*2-2)*N3N+IABC                                          CPP24190
      KAPT=(KTYP*2-1)*N3N+IABC                                          CPP24200
      VALU=VALU+ALPA(ITYP,KTYP)*TM(IJ,JAPT)                             CPP24210
     1         +BETA(ITYP,KTYP)*TM(IJ,KAPT)                             CPP24220
  104 CONTINUE                                                          CPP24230
      EPA(I,J)=VALU                                                     CPP24240
  105 CONTINUE                                                          CPP24250
C                                                                       CPP24260
      CALL RWRIT(ITAP44,EPA,NBASQ,IE)                                   CPP24270
      IF(IPRNT.LE.4) GO TO 110                                          CPP24280
      WRITE(6,1) IABC                                                   CPP24290
      CALL MATOUT(EPA,NBASIS,NBASIS,NBASIS,NBASIS,6)                    CPP24300
  110 CONTINUE                                                          CPP24310
C                                                                       CPP24320
C   STORE DERIVATIVE AVERAGED FOCKS                                     CPP24330
      DO 115 IABC=1,N3N                                                 CPP24340
      IT=IFA(IABC)                                                      CPP24350
      DO 114 I=1,NTRI                                                   CPP24360
  114 FF(I)=HM(I,IABC)+TM(I,IABC)                                       CPP24370
      CALL RWRIT(ITAP44,FF,NTRI2,IT)                                    CPP24380
  115 CONTINUE                                                          CPP24390
C                                                                       CPP24400
      CALL RCLOSE(ITAP44,3)                                             CPP24410
      CALL RCLOSE(JTAPE1,3)                                             CPP24420
      RETURN                                                            CPP24430
      END                                                               CPP24440
      SUBROUTINE BAIND(ESO,NIJ,ZETA,U,T,B0,EPA,DP,DQ,DM,SM,MIJ,PQ,      CPP24450
     1                 NADD)                                            CPP24460
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP24470
      INTEGER*2 MIJ(*)                                                  CPP24480
      DIMENSION ESO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)         CPP24490
      DIMENSION NIJ(NTRI,2),ZETA(NTRI,NTYPEP)                           CPP24500
      DIMENSION B0(NTRI),EPA(NBASIS,NBASIS)                             CPP24510
      DIMENSION DP(NTRI,NADD*NTYPEP),DQ(NTRI,NADD*NTYPEP)               CPP24520
      DIMENSION DM(NTRI,NADD*NTYPEP),SM(NTRI,NADD*NTYPEP)               CPP24530
      DIMENSION PQ(*)                                                   CPP24540
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP24550
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP24560
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP24570
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP24580
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP24590
      COMMON/TAPES/ITAP37,ITAP42,ITAP43,ITAP44,ITAP47                   CPP24600
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP24610
      COMMON/CI102/NIND,NDEP                                            CPP24620
      DATA A00 / 0.0D+00 /                                              CPP24630
    1 FORMAT(//,2X,' ZETA MATRIX, ITYP = ',I5/)                         CPP24640
    2 FORMAT(//,2X,' IGROUP = ',I5,5X,' IFIRST = ',I5,5X,               CPP24650
     1             ' NABC   = ',I5/)                                    CPP24660
    3 FORMAT(//,2X,' FM MATRIX, IABC = ',I5/)                           CPP24670
    4 FORMAT(2X,3I5,4F15.8)                                             CPP24680
    5 FORMAT(//,2X,' B0 MATRIX IN BAIND'/)                              CPP24690
    6 FORMAT(//,2X,' DP MATRIX IN BAIND'/)                              CPP24700
C                                                                       CPP24710
C   FORM B0 MATRICES                                                    CPP24720
C   SEE EQ.(26)                                                         CPP24730
C                                                                       CPP24740
      DO 101 ITYP=1,NTYPES                                              CPP24750
      CALL RMREAD(ZETA(1,ITYP),35+ITYP)                                 CPP24760
      IF(IPRNT.LE.2) GO TO 101                                          CPP24770
      WRITE(6,1) ITYP                                                   CPP24780
      CALL PRINT(ZETA(1,ITYP),NTRI,NBASIS,6)                            CPP24790
  101 CONTINUE                                                          CPP24800
      CALL ZERO(ZETA(1,NTYPEP),NTRI)                                    CPP24810
C                                                                       CPP24820
      CALL RFILE(ITAP44)                                                CPP24830
C                                                                       CPP24840
C   START GROUPING THE BA MATRICES                                      CPP24850
      IGROUP=0                                                          CPP24860
      ILAST=0                                                           CPP24870
  500 CONTINUE                                                          CPP24880
      IGROUP=IGROUP+1                                                   CPP24890
      IFIRST=ILAST+1                                                    CPP24900
      ILAST=IFIRST+NADD-1                                               CPP24910
      NABC=ILAST-IFIRST+1                                               CPP24920
      WRITE(3,2) IGROUP,IFIRST,NABC                                     CPP24930
      WRITE(4,2) IGROUP,IFIRST,NABC                                     CPP24940
C                                                                       CPP24950
      DO 120 IABC=1,NABC                                                CPP24960
      III=IABC+IFIRST-1                                                 CPP24970
      IS=ISA(III)                                                       CPP24980
      IE=IEA(III)                                                       CPP24990
      IB=IBA(III)                                                       CPP25000
      CALL RREAD(ITAP44,SM(1,IABC),NTRI2,IS)                            CPP25010
      CALL RREAD(ITAP44,EPA,NBASQ,IE)                                   CPP25020
      CALL ZERO(B0,NTRI)                                                CPP25030
C                                                                       CPP25040
C   ELEMENTS FOR INDEPENDENT PAIRS                                      CPP25050
C     IF(IPRNT.LE.4) GO TO 301                                          CPP25060
C     WRITE(6,3) IABC                                                   CPP25070
C     CALL MATOUT(EPA,NBASIS,NBASIS,NBASIS,NBASIS,6)                    CPP25080
  301 DO 110 II=1,NIND                                                  CPP25090
      I=NIJ(II,1)                                                       CPP25100
      J=NIJ(II,2)                                                       CPP25110
      ITYP=MOTYP(I)                                                     CPP25120
      JTYP=MOTYP(J)                                                     CPP25130
      IJ=IOFF(I)+J                                                      CPP25140
      VALU1=EPA(I,J)-EPA(J,I)                                           CPP25150
      VALU2=A00                                                         CPP25160
C                                                                       CPP25170
      DO 107 K=2,NBASIS                                                 CPP25180
      LIM=MIN0(JOCC,K-1)                                                CPP25190
      DO 107 L=1,LIM                                                    CPP25200
      KL=IOFF(K)+L                                                      CPP25210
      FAC2=A00                                                          CPP25220
      IF(K.NE.J) GO TO 205                                              CPP25230
      IL=IOFF(MAX0(I,L))+MIN0(I,L)                                      CPP25240
      FAC2=FAC2+ZETA(IL,ITYP)-ZETA(IL,JTYP)                             CPP25250
  205 IF(K.NE.I) GO TO 206                                              CPP25260
      JL=IOFF(MAX0(J,L))+MIN0(J,L)                                      CPP25270
      FAC2=FAC2-ZETA(JL,JTYP)+ZETA(JL,ITYP)                             CPP25280
  206 VALU2=VALU2-SM(KL,IABC)*FAC2                                      CPP25290
  107 CONTINUE                                                          CPP25300
      VALUT=VALU1+VALU2                                                 CPP25310
      B0(IJ)=VALUT                                                      CPP25320
C     IF(IPRNT.LE.4) GO TO 110                                          CPP25330
C     WRITE(6,4) II,I,J,VALU1,VALU2,VALUT                               CPP25340
  110 CONTINUE                                                          CPP25350
      CALL RWRIT(ITAP44,B0,NTRI2,IB)                                    CPP25360
      IF(IPRNT.LE.4) GO TO 120                                          CPP25370
      WRITE(6,5)                                                        CPP25380
      CALL PRINT(B0,NTRI,NBASIS,6)                                      CPP25390
  120 CONTINUE                                                          CPP25400
C                                                                       CPP25410
C   FORM DENSITY LIKE MATRIX IN SO BASIS                                CPP25420
  302 CONTINUE                                                          CPP25430
      CALL PQBMAT(DP,DQ,SM,ESO,NABC)                                    CPP25440
C                                                                       CPP25450
C   FORM HALF-TRANSFORMED DENSITY MATRIX                                CPP25460
      NTYABC=NTYPES*NABC                                                CPP25470
      CALL PQMAT(DP,DQ,DM,NTYABC,MIJ,PQ)                                CPP25480
C                                                                       CPP25490
C   FORM DENSITY LIKE MATRIX IN MO BASIS                                CPP25500
      DO 125 ITYPX=1,NTYABC                                             CPP25510
      CALL MOCONS(DM(1,ITYPX),NTRI,DP(1,ITYPX),NTRI,ESO,U,T)            CPP25520
  125 CONTINUE                                                          CPP25530
      IF(IPRNT.LE.4) GO TO 303                                          CPP25540
      WRITE(6,6)                                                        CPP25550
      CALL MATOUT(DP,NTRI,NTYABC,NTRI,NTYABC,6)                         CPP25560
C                                                                       CPP25570
C   COMPLETE B0 MATRICES                                                CPP25580
  303 CONTINUE                                                          CPP25590
      DO 130 IABC=1,NABC                                                CPP25600
      III=IABC+IFIRST-1                                                 CPP25610
      IB=IBA(III)                                                       CPP25620
      CALL RREAD(ITAP44,B0,NTRI2,IB)                                    CPP25630
      DO 129 II=1,NIND                                                  CPP25640
      I=NIJ(II,1)                                                       CPP25650
      ITYP=MOTYP(I)                                                     CPP25660
      IAPT=(ITYP-1)*NABC+IABC                                           CPP25670
      J=NIJ(II,2)                                                       CPP25680
      JTYP=MOTYP(J)                                                     CPP25690
      JAPT=(JTYP-1)*NABC+IABC                                           CPP25700
      IJ=IOFF(I)+J                                                      CPP25710
      B0(IJ)=B0(IJ)+DP(IJ,IAPT)-DP(IJ,JAPT)                             CPP25720
  129 CONTINUE                                                          CPP25730
      CALL RWRIT(ITAP44,B0,NTRI2,IB)                                    CPP25740
  130 CONTINUE                                                          CPP25750
      IF(ILAST.GE.N3N) GO TO 501                                        CPP25760
      GO TO 500                                                         CPP25770
C                                                                       CPP25780
  501 CONTINUE                                                          CPP25790
      CALL RCLOSE(ITAP44,3)                                             CPP25800
      RETURN                                                            CPP25810
      END                                                               CPP25820
      SUBROUTINE PQBMAT(DP,DQ,SM,ESO,NABC)                              CPP25830
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP25840
      DIMENSION DP(NTRI,NABC*NTYPEP),DQ(NTRI,NABC*NTYPEP)               CPP25850
      DIMENSION SM(NTRI,NABC*NTYPEP)                                    CPP25860
      DIMENSION ESO(NBFAO,NBASIS)                                       CPP25870
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP25880
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP25890
      COMMON/GVBAB/A(30,30),B(30,30)                                    CPP25900
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP25910
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP25920
      COMMON/TOLER/CRIT,TOLR                                            CPP25930
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP25940
      DATA TWO / 2.0D+00 /                                              CPP25950
    1 FORMAT(//,2X,' DP MATRIX'/)                                       CPP25960
    2 FORMAT(//,2X,' DQ MATRIX'/)                                       CPP25970
C                                                                       CPP25980
C   FORM DENSITY LIKE MATRIX                                            CPP25990
C   SEE EQ.(A-2)                                                        CPP26000
C                                                                       CPP26010
      CALL ZERO(DP,NTRI*NTYPEP*NABC)                                    CPP26020
      CALL ZERO(DQ,NTRI*NTYPEP*NABC)                                    CPP26030
      IJ=0                                                              CPP26040
      DO 105 I=1,NBASIS                                                 CPP26050
      DO 105 J=1,I                                                      CPP26060
      IJ=IJ+1                                                           CPP26070
      DO 103 K=1,JOCC                                                   CPP26080
      KTYP=MOTYP(K)                                                     CPP26090
      DO 103 L=1,K                                                      CPP26100
      KL=IOFF(K)+L                                                      CPP26110
      FAC=ESO(I,K)*ESO(J,L)+ESO(J,K)*ESO(I,L)                           CPP26120
      IF(K.EQ.L) FAC=ESO(I,K)*ESO(J,K)                                  CPP26130
      IF(DABS(FAC).LT.TOLR) GO TO 103                                   CPP26140
      DO 102 ITYP=1,NTYPES                                              CPP26150
      AVAL=A(ITYP,KTYP)                                                 CPP26160
      BVAL=B(ITYP,KTYP)                                                 CPP26170
      FACA=AVAL*FAC                                                     CPP26180
      FACB=BVAL*FAC                                                     CPP26190
      IADD=(ITYP-1)*NABC                                                CPP26200
      DO 101 IABC=1,NABC                                                CPP26210
      IAPT=IADD+IABC                                                    CPP26220
      DP(IJ,IAPT)=DP(IJ,IAPT)-SM(KL,IABC)*FACA                          CPP26230
      DQ(IJ,IAPT)=DQ(IJ,IAPT)-SM(KL,IABC)*FACB                          CPP26240
  101 CONTINUE                                                          CPP26250
  102 CONTINUE                                                          CPP26260
  103 CONTINUE                                                          CPP26270
  105 CONTINUE                                                          CPP26280
C                                                                       CPP26290
      NTYABC=NTYPES*NABC                                                CPP26300
      IF(IPRNT.LE.4) GO TO 201                                          CPP26310
      WRITE(6,1)                                                        CPP26320
      CALL MATOUT(DP,NTRI,NTYABC,NTRI,NTYABC,6)                         CPP26330
      WRITE(6,2)                                                        CPP26340
      CALL MATOUT(DQ,NTRI,NTYABC,NTRI,NTYABC,6)                         CPP26350
C                                                                       CPP26360
  201 CONTINUE                                                          CPP26370
      DO 107 I=1,NBASIS                                                 CPP26380
      II=IOFF(I+1)                                                      CPP26390
      DO 106 ITYPX=1,NTYABC                                             CPP26400
      DP(II,ITYPX)=DP(II,ITYPX)/TWO                                     CPP26410
      DQ(II,ITYPX)=DQ(II,ITYPX)/TWO                                     CPP26420
  106 CONTINUE                                                          CPP26430
  107 CONTINUE                                                          CPP26440
      CALL SCLMPY(DP,NTRI*NTYABC,TWO)                                   CPP26450
      CALL SCLMPY(DQ,NTRI*NTYABC,TWO)                                   CPP26460
C                                                                       CPP26470
      RETURN                                                            CPP26480
      END                                                               CPP26490
      SUBROUTINE BADEP(EIG,ESO,KIJ,U,T,B0,FF,DP,DQ,DM,MIJ,PQ,           CPP26500
     1                 NADD)                                            CPP26510
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP26520
      INTEGER*2 MIJ(*)                                                  CPP26530
      DIMENSION EIG(NBASIS),KIJ(NTRI,2)                                 CPP26540
      DIMENSION ESO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)         CPP26550
      DIMENSION B0(NTRI),FF(NTRI)                                       CPP26560
      DIMENSION DP(NTRI,NADD),DQ(NTRI,NADD),DM(NTRI,NADD)               CPP26570
      DIMENSION PQ(*)                                                   CPP26580
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP26590
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP26600
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          CPP26610
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP26620
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP26630
      COMMON/TAPES/ITAP37,ITAP42,ITAP43,ITAP44,ITAP47                   CPP26640
      COMMON/CI102/NIND,NDEP                                            CPP26650
    1 FORMAT(//,2X,' IGROUP = ',I5,5X,' IFIRST = ',I5,5X,               CPP26660
     1             ' NABC   = ',I5/)                                    CPP26670
    2 FORMAT(//,2X,' B0 MATRIX IN BADEP'/)                              CPP26680
    3 FORMAT(//,2X,' DM MATRIX IN BADEP'/)                              CPP26690
C                                                                       CPP26700
      CALL RFILE(ITAP44)                                                CPP26710
C                                                                       CPP26720
C   START GROUPING BA MATRICES                                          CPP26730
      IGROUP=0                                                          CPP26740
      ILAST=0                                                           CPP26750
  500 CONTINUE                                                          CPP26760
      IGROUP=IGROUP+1                                                   CPP26770
      IFIRST=ILAST+1                                                    CPP26780
      ILAST=IFIRST+NADD-1                                               CPP26790
      IF(ILAST.GT.N3N) ILAST=N3N                                        CPP26800
      NABC=ILAST-IFIRST+1                                               CPP26810
      WRITE(3,1) IGROUP,IFIRST,NABC                                     CPP26820
      WRITE(4,1) IGROUP,IFIRST,NABC                                     CPP26830
C                                                                       CPP26840
C   DM ARRAY STORES DERIVATIVE OVERLAP MATRIX                           CPP26850
      DO 102 IABC=1,NABC                                                CPP26860
      III=IABC+IFIRST-1                                                 CPP26870
      IS=ISA(III)                                                       CPP26880
      IT=IFA(III)                                                       CPP26890
      IB=IBA(III)                                                       CPP26900
      CALL RREAD(ITAP44,DM(1,IABC),NTRI2,IS)                            CPP26910
      CALL RREAD(ITAP44,FF,NTRI2,IT)                                    CPP26920
      CALL RREAD(ITAP44,B0,NTRI2,IB)                                    CPP26930
      IJ=0                                                              CPP26940
      DO 101 II=1,NDEP                                                  CPP26950
      I=KIJ(II,1)                                                       CPP26960
      J=KIJ(II,2)                                                       CPP26970
      IJ=IOFF(I)+J                                                      CPP26980
      B0(IJ)=FF(IJ)-EIG(J)*DM(IJ,IABC)                                  CPP26990
  101 CONTINUE                                                          CPP27000
      CALL RWRIT(ITAP44,B0,NTRI2,IB)                                    CPP27010
      IF(IPRNT.LE.4) GO TO 102                                          CPP27020
      WRITE(6,2)                                                        CPP27030
      CALL PRINT(B0,NTRI,NBASIS,6)                                      CPP27040
  102 CONTINUE                                                          CPP27050
C                                                                       CPP27060
C   FORM DENSITY LIKE MATRIX                                            CPP27070
C   DM ARRAY CONTAINES DERIVATIVE OVERLAP MATRIX                        CPP27080
C   DP ARRAY STORES DENSITY LIKE MATRIX IN SO BASIS                     CPP27090
      CALL PQCMAT(DP,DQ,DM,ESO,NABC)                                    CPP27100
C                                                                       CPP27110
C   FORM HALF-TRANSFORMED DENSITY LIKE MATRIX                           CPP27120
C   AND TRANSFORM THEM INTO MO BASIS                                    CPP27130
C   DP ARRAY CONTAINS DENSITY LIKE MATRIX IN SO BASIS                   CPP27140
C   DM ARRAY STORES HALF-TRANSFORMED DENSITY LIKE MATRIX                CPP27150
      CALL PQMAT(DP,DQ,DM,NABC,MIJ,PQ)                                  CPP27160
      IF(IPRNT.LE.4) GO TO 302                                          CPP27170
      WRITE(6,3)                                                        CPP27180
      CALL MATOUT(DM,NTRI,NABC,NTRI,NABC,6)                             CPP27190
C                                                                       CPP27200
C   COMPLETE B0 MATRICES                                                CPP27210
  302 CONTINUE                                                          CPP27220
      CALL SREW(ITAP44)                                                 CPP27230
      DO 105 IABC=1,NABC                                                CPP27240
      III=IABC+IFIRST-1                                                 CPP27250
      IB=IBA(III)                                                       CPP27260
      CALL RREAD(ITAP44,B0,NTRI2,IB)                                    CPP27270
      CALL MOCONS(DM(1,IABC),NTRI,DP(1,IABC),NTRI,ESO,U,T)              CPP27280
      DO 104 II=1,NDEP                                                  CPP27290
      I=KIJ(II,1)                                                       CPP27300
      J=KIJ(II,2)                                                       CPP27310
      IJ=IOFF(I)+J                                                      CPP27320
      B0(IJ)=B0(IJ)+DP(IJ,IABC)                                         CPP27330
  104 CONTINUE                                                          CPP27340
      CALL RWRIT(ITAP44,B0,NTRI2,IB)                                    CPP27350
  105 CONTINUE                                                          CPP27360
      IF(ILAST.GE.N3N) GO TO 501                                        CPP27370
      GO TO 500                                                         CPP27380
C                                                                       CPP27390
  501 CONTINUE                                                          CPP27400
      CALL RCLOSE(ITAP44,3)                                             CPP27410
      RETURN                                                            CPP27420
      END                                                               CPP27430
      SUBROUTINE PQCMAT(DP,DQ,SM,ESO,NABC)                              CPP27440
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP27450
      DIMENSION DP(NTRI,NABC),DQ(NTRI,NABC),SM(NTRI,NABC)               CPP27460
      DIMENSION ESO(NBFAO,NBASIS)                                       CPP27470
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP27480
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP27490
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP27500
      COMMON/TOLER/CRIT,TOLR                                            CPP27510
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP27520
      DATA HALF,TWO / 0.5D+00 , 2.0D+00 /                               CPP27530
    1 FORMAT(//,2X,' DP MATRIX IN PQCMAT'/)                             CPP27540
    2 FORMAT(//,2X,' DQ MATRIX IN PQCMAT'/)                             CPP27550
C                                                                       CPP27560
C   FORM DENSITY LIKE MATRICES                                          CPP27570
C                                                                       CPP27580
      CALL ZERO(DP,NTRI*NABC)                                           CPP27590
      CALL ZERO(DQ,NTRI*NABC)                                           CPP27600
C                                                                       CPP27610
      IJ=0                                                              CPP27620
      DO 103 I=1,NBASIS                                                 CPP27630
      DO 103 J=1,I                                                      CPP27640
      IJ=IJ+1                                                           CPP27650
      DO 102 K=1,JOCC                                                   CPP27660
      KTYP=MOTYP(K)                                                     CPP27670
      FACO=FOCC(KTYP)*HALF                                              CPP27680
      DO 102 L=1,K                                                      CPP27690
      KL=IOFF(MAX0(K,L))+MIN0(K,L)                                      CPP27700
      FACV=ESO(I,K)*ESO(J,L)+ESO(J,K)*ESO(I,L)                          CPP27710
      IF(K.EQ.L) FACV=ESO(I,K)*ESO(J,K)                                 CPP27720
      FAC=FACO*FACV                                                     CPP27730
      IF(DABS(FAC).LT.TOLR) GO TO 102                                   CPP27740
      DO 101 IABC=1,NABC                                                CPP27750
      FACT=-SM(KL,IABC)*FAC                                             CPP27760
      DP(IJ,IABC)=DP(IJ,IABC)+FACT*TWO                                  CPP27770
      DQ(IJ,IABC)=DQ(IJ,IABC)-FACT                                      CPP27780
  101 CONTINUE                                                          CPP27790
  102 CONTINUE                                                          CPP27800
  103 CONTINUE                                                          CPP27810
C                                                                       CPP27820
      IF(IPRNT.LE.4) GO TO 201                                          CPP27830
      WRITE(6,1)                                                        CPP27840
      CALL MATOUT(DP,NTRI,NABC,NTRI,NABC,6)                             CPP27850
      WRITE(6,2)                                                        CPP27860
      CALL MATOUT(DQ,NTRI,NABC,NTRI,NABC,6)                             CPP27870
  201 CONTINUE                                                          CPP27880
      DO 105 I=1,NBASIS                                                 CPP27890
      II=IOFF(I+1)                                                      CPP27900
      DO 104 IABC=1,NABC                                                CPP27910
      DP(II,IABC)=DP(II,IABC)/TWO                                       CPP27920
      DQ(II,IABC)=DQ(II,IABC)/TWO                                       CPP27930
  104 CONTINUE                                                          CPP27940
  105 CONTINUE                                                          CPP27950
      CALL SCLMPY(DP,NTRI*NABC,TWO)                                     CPP27960
      CALL SCLMPY(DQ,NTRI*NABC,TWO)                                     CPP27970
C                                                                       CPP27980
      RETURN                                                            CPP27990
      END                                                               CPP28000
      SUBROUTINE HEMAT(ESO,FIJ,AIJ,BIJ,HCI,HSO,HMO,DP,DQ,DM,EIJ,ZIJ,EE, CPP28010
     1                 U,T,MIJ,PQ)                                      CPP28020
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP28030
      INTEGER*2 MIJ(*)                                                  CPP28040
      DIMENSION FIJ(NTYPEP,ICIUNQ),AIJ(NTYTRI,ICIUNQ),BIJ(NTYTRI,ICIUNQ)CPP28050
      DIMENSION ESO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)         CPP28060
      DIMENSION HSO(NTRI),HMO(NTRI),HCI(NCTRI),ZIJ(NTRI)                CPP28070
      DIMENSION EIJ(NBASIS,NBASIS),EE(NBASIS,NBASIS)                    CPP28080
      DIMENSION DP(NTRI,NTYPEP),DQ(NTRI,NTYPEP),DM(NTRI,NTYPEP)         CPP28090
      DIMENSION PQ(*)                                                   CPP28100
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP28110
      COMMON/CIOFF/IIUNQ,IJUNQ,ICIUNQ,IPX(5050),JPX(5050)               CPP28120
      COMMON/CIPOS/ICITYP(5050),ICPOS(200,100)                          CPP28130
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP28140
      COMMON/CITYP/IJGRP,IBM(200),JBM(200),IBNUM(200)                   CPP28150
      COMMON/ENRGY/ESCF,ENUC,ELEC                                       CPP28160
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP28170
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP28180
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP28190
      COMMON/TAPES/ITAP37,ITAP42,ITAP43,ITAP44,ITAP47                   CPP28200
      COMMON/WORKS/JTAPE1,JTAPE2                                        CPP28210
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP28220
      COMMON/CI103/CI(100)                                              CPP28230
      COMMON/CI104/IEIJ(5050)                                           CPP28240
      DATA A00,HALF,TWO / 0.0D+00 , 0.5D+00 , 2.0D+00 /                 CPP28250
    1 FORMAT(//,2X,' H (SO BASIS) MATRIX'/)                             CPP28260
    2 FORMAT(//,2X,' H (MO BASIS) MATRIX'/)                             CPP28270
    3 FORMAT(//,2X,' DM MATRIX'/)                                       CPP28280
    4 FORMAT(//,2X,' DP MATRIX'/)                                       CPP28290
    5 FORMAT(2X,' ICI , JCI , IJCI , ELEC1 , ELEC2 ',3I5,2F15.8)        CPP28300
    6 FORMAT(//,2X,' EIJ MATRIX, IJCI = ',I5,' IJEX = ',I5/)            CPP28310
    7 FORMAT(//,2X,' ZIJ MATRIX, IJCI = ',I5,' IJEX = ',I5,' ITYP = ',  CPP28320
     * I5/)                                                             CPP28330
    8 FORMAT(//,2X,' EIJ MATRIX, IBCI = ',I5/)                          CPP28340
    9 FORMAT(//,2X,' ZIJ MATRIX, IBCI = ',I5,' ITYP = ',I5/)            CPP28350
   10 FORMAT(//,2X,' CI MATRIX'/)                                       CPP28360
   11 FORMAT(//,2X,' CI ENERGIES'/                                      CPP28370
     * 2X,' ESCF       = ',F20.10/                                      CPP28380
     * 2X,' ENUC       = ',F20.10/                                      CPP28390
     * 2X,' ELEC       = ',F20.10/                                      CPP28400
     * 2X,' ECITOT     = ',F20.10/)                                     CPP28410
   12 FORMAT(//,2X,' THE LAGRANGIAN MATRIX'/)                           CPP28420
C                                                                       CPP28430
C   CALCULATE H MATRIX                                                  CPP28440
C   SEE EQS.(5)-(7)                                                     CPP28450
C                                                                       CPP28460
C   READ IN ONE ELECTRON INTEGRALS (SO BASIS)                           CPP28470
      CALL RMREAD(HSO,16)                                               CPP28480
      CALL MOCONS(HSO,NTRI,HMO,NTRI,ESO,U,T)                            CPP28490
      IF(IPRNT.LE.4) GO TO 201                                          CPP28500
      WRITE(6,1)                                                        CPP28510
      CALL PRINT(HSO,NTRI,NBASIS,6)                                     CPP28520
      WRITE(6,2)                                                        CPP28530
      CALL PRINT(HMO,NTRI,NBASIS,6)                                     CPP28540
C                                                                       CPP28550
C   FORM HIJ AND EIJ MATRICES                                           CPP28560
  201 CONTINUE                                                          CPP28570
      CALL RFILE(ITAP47)                                                CPP28580
      CALL RFILE(JTAPE2)                                                CPP28590
      CALL ZERO(HCI,NCTRI)                                              CPP28600
C                                                                       CPP28610
C   FOR CI DIAGONAL ELEMENTS                                            CPP28620
      DO 110 IJCI=1,IIUNQ                                               CPP28630
      IEX=IPX(IJCI)                                                     CPP28640
      JEX=JPX(IJCI)                                                     CPP28650
      IJEX=IOFF(MAX0(IEX,JEX))+MIN0(IEX,JEX)                            CPP28660
C                                                                       CPP28670
C   FORM DENSITY MATRIX IN SO BASIS                                     CPP28680
C   AND FORM HALF-TRANSFORMED TWO ELECTRON MATRIX                       CPP28690
      CALL PQDMAT(DP,DQ,ESO,AIJ,BIJ,IJCI)                               CPP28700
      CALL PQMAT(DP,DQ,DM,NTYPES,MIJ,PQ)                                CPP28710
C                                                                       CPP28720
C   FINAL HALF TRANSFORMATION                                           CPP28730
      DO 102 ITYP=1,NTYPES                                              CPP28740
      CALL MOCONS(DM(1,ITYP),NTRI,DP(1,ITYP),NTRI,ESO,U,T)              CPP28750
  102 CONTINUE                                                          CPP28760
      IF(IPRNT.LE.4) GO TO 202                                          CPP28770
      WRITE(6,3)                                                        CPP28780
      CALL MATOUT(DM,NTRI,NTYPEP,NTRI,NTYPEP,6)                         CPP28790
      WRITE(6,4)                                                        CPP28800
      CALL MATOUT(DP,NTRI,NTYPEP,NTRI,NTYPEP,6)                         CPP28810
C                                                                       CPP28820
C   FORM THE CI HAMILTONIAN MATRIX                                      CPP28830
  202 CONTINUE                                                          CPP28840
      ELEC1=A00                                                         CPP28850
      ELEC2=A00                                                         CPP28860
      DO 103 I=1,JOCC                                                   CPP28870
      II=IOFF(I+1)                                                      CPP28880
      ITYP=MOTYP(I)                                                     CPP28890
      ELEC1=ELEC1+HMO(II)*FIJ(ITYP,IJCI)*TWO                            CPP28900
      ELEC2=ELEC2+DP(II,ITYP)                                           CPP28910
  103 CONTINUE                                                          CPP28920
      HCI(IJEX)=ELEC1+ELEC2                                             CPP28930
      IF(IPRNT.LE.2) GO TO 203                                          CPP28940
      WRITE(6,5) IEX,JEX,IJEX,ELEC1,ELEC2                               CPP28950
C                                                                       CPP28960
C   FORM THE "BARE" LAGRANGIAN MATRICES                                 CPP28970
  203 CONTINUE                                                          CPP28980
      CALL ZERO(EIJ,NBASIS*NBASIS)                                      CPP28990
      DO 104 I=1,JOCC                                                   CPP29000
      ITYP=MOTYP(I)                                                     CPP29010
      FAC=FIJ(ITYP,IJCI)                                                CPP29020
      DO 104 J=1,NBASIS                                                 CPP29030
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPP29040
      ELEC1=HMO(IJ)*FAC                                                 CPP29050
      EIJ(I,J)=HMO(IJ)*FAC+DP(IJ,ITYP)                                  CPP29060
  104 CONTINUE                                                          CPP29070
      IPT=IEIJ(IJCI)                                                    CPP29080
      CALL RWRIT(JTAPE2,EIJ,NBASQ,IPT)                                  CPP29090
      IF(IPRNT.LE.2) GO TO 204                                          CPP29100
      WRITE(6,6) IJCI,IJEX                                              CPP29110
      CALL MATOUT(EIJ,NBASIS,NBASIS,NBASIS,NBASIS,6)                    CPP29120
C                                                                       CPP29130
  204 CONTINUE                                                          CPP29140
C   FORM THE "BARE" GENERALIZED LAGRANGIAN MATRICES                     CPP29150
      DO 106 ITYP=1,NTYPES                                              CPP29160
      CALL ZERO(ZIJ,NTRI)                                               CPP29170
      FAC=FIJ(ITYP,IJCI)                                                CPP29180
      IJ=0                                                              CPP29190
      DO 105 I=1,NBASIS                                                 CPP29200
      DO 105 J=1,I                                                      CPP29210
      IJ=IJ+1                                                           CPP29220
      ELEC1=HMO(IJ)*FAC                                                 CPP29230
      ZIJ(IJ)=HMO(IJ)*FAC+DP(IJ,ITYP)                                   CPP29240
  105 CONTINUE                                                          CPP29250
      CALL SWRIT(ITAP47,ZIJ,NTRI2)                                      CPP29260
      IF(IPRNT.LE.2) GO TO 106                                          CPP29270
      WRITE(6,7) IJCI,IJEX,ITYP                                         CPP29280
      CALL PRINT(ZIJ,NTRI,NBASIS,6)                                     CPP29290
  106 CONTINUE                                                          CPP29300
C                                                                       CPP29310
  110 CONTINUE                                                          CPP29320
C                                                                       CPP29330
C::::::::::::::::::::::::::::::::::                                     CPP29340
C:::FOR CI OFF-DIAGONAL ELEMENTS:::                                     CPP29350
C::::::::::::::::::::::::::::::::::                                     CPP29360
      DO 120 IBCI=1,IJGRP                                               CPP29370
      IP=IBM(IBCI)                                                      CPP29380
      IQ=JBM(IBCI)                                                      CPP29390
      NUM=IBNUM(IBCI)                                                   CPP29400
C                                                                       CPP29410
C   FORM DENSITY MATRIX IN SO BASIS                                     CPP29420
C   AND FORM HALF-TRANSFORMED TWO ELECTRON MATRIX                       CPP29430
      CALL PQDMTB(DQ,ESO,IP,IQ)                                         CPP29440
      CALL PQQMAT(DQ,DM,NTYPES,MIJ,PQ)                                  CPP29450
C                                                                       CPP29460
C   FINAL HALF TRANSFORMATION                                           CPP29470
      DO 111 ITYP=1,NTYPES                                              CPP29480
      CALL MOCONS(DM(1,ITYP),NTRI,DP(1,ITYP),NTRI,ESO,U,T)              CPP29490
  111 CONTINUE                                                          CPP29500
      IF(IPRNT.LE.4) GO TO 205                                          CPP29510
      WRITE(6,3) IBCI                                                   CPP29520
      CALL MATOUT(DM,NTRI,NTYPEP,NTRI,NTYPEP,6)                         CPP29530
      WRITE(6,4)                                                        CPP29540
      CALL MATOUT(DP,NTRI,NTYPEP,NTRI,NTYPEP,6)                         CPP29550
C                                                                       CPP29560
C   FORM THE CI HAMILTONIAN MATRIX                                      CPP29570
  205 CONTINUE                                                          CPP29580
      ITYP=MOTYP(IP)                                                    CPP29590
      II=IOFF(IP+1)                                                     CPP29600
      ELEC2=DP(II,ITYP)                                                 CPP29610
C***********************************************************************CPP29620
C*  FOLLOWING CODES GIVE THE SAME RESULT                               *CPP29630
C*    ELEC2=A00                                                        *CPP29640
C*    DO 112 I=1,JOCC                                                  *CPP29650
C*    II=IOFF(I+1)                                                     *CPP29660
C*    ITYP=MOTYP(I)                                                    *CPP29670
C*    ELEC2=ELEC2+DP(II,ITYP)*HALF                                     *CPP29680
C*112 CONTINUE                                                         *CPP29690
C***********************************************************************CPP29700
      DO 113 INUM=1,NUM                                                 CPP29710
      IJCI=ICPOS(IBCI,INUM)                                             CPP29720
      IEX=IPX(IIUNQ+IJCI)                                               CPP29730
      JEX=JPX(IIUNQ+IJCI)                                               CPP29740
      IJEX=IOFF(MAX0(IEX,JEX))+MIN0(IEX,JEX)                            CPP29750
      HCI(IJEX)=ELEC2                                                   CPP29760
  113 CONTINUE                                                          CPP29770
C                                                                       CPP29780
C   FORM THE "BARE" LAGRANGIAN MATRICES                                 CPP29790
      CALL ZERO(EIJ,NBASIS*NBASIS)                                      CPP29800
      DO 114 I=1,JOCC                                                   CPP29810
      ITYP=MOTYP(I)                                                     CPP29820
      DO 114 J=1,NBASIS                                                 CPP29830
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPP29840
      EIJ(I,J)=DP(IJ,ITYP)*HALF                                         CPP29850
  114 CONTINUE                                                          CPP29860
      IPT=IEIJ(IIUNQ+IBCI)                                              CPP29870
      CALL RWRIT(JTAPE2,EIJ,NBASQ,IPT)                                  CPP29880
      IF(IPRNT.LE.2) GO TO 206                                          CPP29890
      WRITE(6,8) IBCI                                                   CPP29900
      CALL MATOUT(EIJ,NBASIS,NBASIS,NBASIS,NBASIS,6)                    CPP29910
C                                                                       CPP29920
  206 CONTINUE                                                          CPP29930
C   FORM THE "BARE" GENERALIZED LAGRANGIAN MATRICES                     CPP29940
      DO 116 ITYP=1,NTYPES                                              CPP29950
      CALL ZERO(ZIJ,NTRI)                                               CPP29960
      IF(ITYP.NE.MOTYP(IP).AND.ITYP.NE.MOTYP(IQ)) GO TO 207             CPP29970
      IJ=0                                                              CPP29980
      DO 115 I=1,NBASIS                                                 CPP29990
      DO 115 J=1,I                                                      CPP30000
      IJ=IJ+1                                                           CPP30010
      ZIJ(IJ)=DP(IJ,ITYP)*HALF                                          CPP30020
  115 CONTINUE                                                          CPP30030
  207 CONTINUE                                                          CPP30040
      CALL SWRIT(ITAP47,ZIJ,NTRI2)                                      CPP30050
      IF(IPRNT.LE.2) GO TO 116                                          CPP30060
      WRITE(6,9) IBCI,ITYP                                              CPP30070
      CALL PRINT(ZIJ,NTRI,NBASIS,6)                                     CPP30080
  116 CONTINUE                                                          CPP30090
C                                                                       CPP30100
  120 CONTINUE                                                          CPP30110
      WRITE(6,10)                                                       CPP30120
      CALL PRINT(HCI,NCTRI,NCPEX,6)                                     CPP30130
C                                                                       CPP30140
C   CALCULATE SCF ENERGY FOR A TEST                                     CPP30150
C   CALCULATE THE LAGRANGIAN MATRIX FOR A TEST                          CPP30160
      CALL SREW(JTAPE2)                                                 CPP30170
      CALL ZERO(EE,NBASIS*NBASIS)                                       CPP30180
      ELEC=A00                                                          CPP30190
      DO 122 IJCI=1,ICIUNQ                                              CPP30200
      IEX=IPX(IJCI)                                                     CPP30210
      JEX=JPX(IJCI)                                                     CPP30220
      IJEX=IOFF(MAX0(IEX,JEX))+MIN0(IEX,JEX)                            CPP30230
      IBCI=ICITYP(IJCI)                                                 CPP30240
      IPT=IEIJ(IBCI)                                                    CPP30250
      CALL RREAD(JTAPE2,EIJ,NBASQ,IPT)                                  CPP30260
      CIFAC=CI(IEX)*CI(JEX)*TWO                                         CPP30270
      IF(IEX.EQ.JEX) CIFAC=CIFAC*HALF                                   CPP30280
      ELEC=ELEC+HCI(IJEX)*CIFAC                                         CPP30290
      DO 121 I=1,NBASIS                                                 CPP30300
      DO 121 J=1,NBASIS                                                 CPP30310
      EE(I,J)=EE(I,J)+EIJ(I,J)*CIFAC                                    CPP30320
  121 CONTINUE                                                          CPP30330
  122 CONTINUE                                                          CPP30340
      ECITOT=ENUC+ELEC                                                  CPP30350
      WRITE(6,11) ESCF,ENUC,ELEC,ECITOT                                 CPP30360
      IF(IPRNT.LE.2) GO TO 210                                          CPP30370
      WRITE(6,12)                                                       CPP30380
      CALL MATOUT(EE,NBASIS,NBASIS,NBASIS,NBASIS,6)                     CPP30390
C                                                                       CPP30400
  210 CONTINUE                                                          CPP30410
      CALL RCLOSE(ITAP47,3)                                             CPP30420
      CALL RCLOSE(JTAPE2,3)                                             CPP30430
      RETURN                                                            CPP30440
      END                                                               CPP30450
      SUBROUTINE PQDMAT(DP,DQ,ESO,AIJ,BIJ,IJCI)                         CPP30460
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP30470
      DIMENSION ESO(NBFAO,NBASIS)                                       CPP30480
      DIMENSION DP(NTRI,NTYPEP),DQ(NTRI,NTYPEP)                         CPP30490
      DIMENSION AIJ(NTYTRI,ICIUNQ),BIJ(NTYTRI,ICIUNQ)                   CPP30500
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP30510
      COMMON/CIOFF/IIUNQ,IJUNQ,ICIUNQ,IPX(5050),JPX(5050)               CPP30520
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP30530
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP30540
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP30550
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP30560
      COMMON/TOLER/CRIT,TOLR                                            CPP30570
      DATA A00,TWO / 0.0D+00 , 2.0D+00 /                                CPP30580
    1 FORMAT(//,2X,' DP MATRIX'/)                                       CPP30590
    2 FORMAT(//,2X,' DQ MATRIX'/)                                       CPP30600
C                                                                       CPP30610
      CALL ZERO(DP,NTRI*NTYPEP)                                         CPP30620
      CALL ZERO(DQ,NTRI*NTYPEP)                                         CPP30630
C                                                                       CPP30640
      IJ=0                                                              CPP30650
      DO 104 I=1,NBASIS                                                 CPP30660
      DO 104 J=1,I                                                      CPP30670
      IJ=IJ+1                                                           CPP30680
      DO 103 ITYP=1,NTYPES                                              CPP30690
      SUMP=A00                                                          CPP30700
      SUMQ=A00                                                          CPP30710
      DO 102 JTYP=1,NTYPES                                              CPP30720
      IJTYP=IOFF(MAX0(ITYP,JTYP))+MIN0(ITYP,JTYP)                       CPP30730
      DO 101 K=NSTR(JTYP),NEND(JTYP)                                    CPP30740
      FAC=ESO(I,K)*ESO(J,K)                                             CPP30750
      IF(DABS(FAC).LT.TOLR) GO TO 101                                   CPP30760
      SUMP=SUMP+AIJ(IJTYP,IJCI)*FAC                                     CPP30770
      SUMQ=SUMQ+BIJ(IJTYP,IJCI)*FAC                                     CPP30780
  101 CONTINUE                                                          CPP30790
  102 CONTINUE                                                          CPP30800
      DP(IJ,ITYP)=SUMP                                                  CPP30810
      DQ(IJ,ITYP)=SUMQ                                                  CPP30820
  103 CONTINUE                                                          CPP30830
  104 CONTINUE                                                          CPP30840
      IF(IPRNT.LE.4) GO TO 201                                          CPP30850
      WRITE(6,1)                                                        CPP30860
      CALL MATOUT(DP,NTRI,NTYPEP,NTRI,NTYPEP,6)                         CPP30870
      WRITE(6,2)                                                        CPP30880
      CALL MATOUT(DQ,NTRI,NTYPEP,NTRI,NTYPEP,6)                         CPP30890
                                                                        CPP30900
  201 CONTINUE                                                          CPP30910
      DO 108 I=1,NBASIS                                                 CPP30920
      II=IOFF(I+1)                                                      CPP30930
      DO 106 ITYP=1,NTYPES                                              CPP30940
      DP(II,ITYP)=DP(II,ITYP)/TWO                                       CPP30950
      DQ(II,ITYP)=DQ(II,ITYP)/TWO                                       CPP30960
  106 CONTINUE                                                          CPP30970
  108 CONTINUE                                                          CPP30980
      CALL SCLMPY(DP,NTRI*NTYPES,TWO)                                   CPP30990
      CALL SCLMPY(DQ,NTRI*NTYPES,TWO)                                   CPP31000
C                                                                       CPP31010
      RETURN                                                            CPP31020
      END                                                               CPP31030
      SUBROUTINE PQDMTB(DQ,ESO,IP,IQ)                                   CPP31040
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP31050
      DIMENSION ESO(NBFAO,NBASIS)                                       CPP31060
      DIMENSION DQ(NTRI,NTYPEP)                                         CPP31070
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP31080
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP31090
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP31100
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP31110
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP31120
      COMMON/TOLER/CRIT,TOLR                                            CPP31130
      DATA TWO / 2.0D+00 /                                              CPP31140
    1 FORMAT(//,2X,' DQ MATRIX'/)                                       CPP31150
C                                                                       CPP31160
      CALL ZERO(DQ,NTRI*NTYPEP)                                         CPP31170
C                                                                       CPP31180
      IJ=0                                                              CPP31190
      DO 103 I=1,NBASIS                                                 CPP31200
      DO 103 J=1,I                                                      CPP31210
      IJ=IJ+1                                                           CPP31220
      DO 102 K=1,NBASIS                                                 CPP31230
      KTYP=MOTYP(K)                                                     CPP31240
      IF(KTYP.NE.MOTYP(IP).AND.KTYP.NE.MOTYP(IQ)) GO TO 102             CPP31250
      FAC=ESO(I,K)*ESO(J,K)                                             CPP31260
      IF(DABS(FAC).LT.TOLR) GO TO 102                                   CPP31270
      DO 101 ITYP=1,NTYPES                                              CPP31280
      IF(KTYP.EQ.MOTYP(IP).AND.ITYP.EQ.MOTYP(IQ)) GO TO 201             CPP31290
      IF(KTYP.EQ.MOTYP(IQ).AND.ITYP.EQ.MOTYP(IP)) GO TO 201             CPP31300
      GO TO 101                                                         CPP31310
  201 CONTINUE                                                          CPP31320
      DQ(IJ,ITYP)=DQ(IJ,ITYP)+FAC                                       CPP31330
  101 CONTINUE                                                          CPP31340
  102 CONTINUE                                                          CPP31350
  103 CONTINUE                                                          CPP31360
      IF(IPRNT.LE.4) GO TO 202                                          CPP31370
      WRITE(6,1)                                                        CPP31380
      CALL MATOUT(DQ,NTRI,NTYPEP,NTRI,NTYPEP,6)                         CPP31390
                                                                        CPP31400
  202 CONTINUE                                                          CPP31410
      DO 105 I=1,NBASIS                                                 CPP31420
      II=IOFF(I+1)                                                      CPP31430
      DO 104 ITYP=1,NTYPES                                              CPP31440
      DQ(II,ITYP)=DQ(II,ITYP)/TWO                                       CPP31450
  104 CONTINUE                                                          CPP31460
  105 CONTINUE                                                          CPP31470
      CALL SCLMPY(DQ,NTRI*NTYPES,TWO)                                   CPP31480
C                                                                       CPP31490
      RETURN                                                            CPP31500
      END                                                               CPP31510
      SUBROUTINE BFMAT(NIJ,KIJ,FIJ,E1T,HCIA,HF,HM,EPF,BF,NABC,ISTEP)    CPP31520
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP31530
      DIMENSION NIJ(NTRI,2),KIJ(NTRI,2)                                 CPP31540
      DIMENSION FIJ(NTYPEP,ICIUNQ),E1T(NABC),HCIA(NCTRI,NABC)           CPP31550
      DIMENSION HF(NTRI,3),HM(NTRI,3,NTYPES)                            CPP31560
      DIMENSION EPF(NBASIS,NBASIS),BF(NTRI)                             CPP31570
      DIMENSION DIP(3)                                                  CPP31580
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP31590
      COMMON/CIOFF/IIUNQ,IJUNQ,ICIUNQ,IPX(5050),JPX(5050)               CPP31600
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP31610
      COMMON/DIPOL/DIPDA(135),DIPDM(135),DIPDN(135),DIPDT(135)          CPP31620
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP31630
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP31640
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP31650
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP31660
      COMMON/TAPES/ITAP37,ITAP42,ITAP43,ITAP44,ITAP47                   CPP31670
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP31680
      COMMON/CI102/NIND,NDEP                                            CPP31690
      COMMON/CI103/CI(100)                                              CPP31700
      DATA A00,HALF,TWO / 0.0D+00 , 0.5D+00 , 2.0D+00 /                 CPP31710
      DATA DEBYE / 2.541765480D+00 /                                    CPP31720
    1 FORMAT(//,2X,' DIPDA MATRIX'/)                                    CPP31730
    2 FORMAT(//,2X,' DIPDN MATRIX'/)                                    CPP31740
    3 FORMAT(//,2X,' HF MATRIX, IXYZ = ',I5/)                           CPP31750
    4 FORMAT(//,2X,' EPF MATRIX, IXYZ = ',I5/)                          CPP31760
    5 FORMAT(//,2X,' HCIA MATRIX, IJCI = ',I5/)                         CPP31770
    6 FORMAT(//,2X,' DIPOLE MOMENT (ELECTRONIC PART IN DEBYE)'/)        CPP31780
C                                                                       CPP31790
      CALL RFILE(ITAP43)                                                CPP31800
      CALL SREAD(ITAP43,DIPDA,270)                                      CPP31810
      CALL SREAD(ITAP43,DIPDN,270)                                      CPP31820
      IF(IPRNT.LE.4) GO TO 201                                          CPP31830
      WRITE(6,1)                                                        CPP31840
      CALL MATOUT(DIPDA,3,45,3,N3N,6)                                   CPP31850
      WRITE(6,2)                                                        CPP31860
      CALL MATOUT(DIPDN,3,45,3,N3N,6)                                   CPP31870
C                                                                       CPP31880
  201 CONTINUE                                                          CPP31890
      CALL RFILE(ITAP44)                                                CPP31900
      DO 103 IXYZ=1,3                                                   CPP31910
      IH=IHA(IXYZ+N3N)                                                  CPP31920
      CALL SREAD(ITAP43,HF(1,IXYZ),NTRI2)                               CPP31930
      CALL RWRIT(ITAP44,HF(1,IXYZ),NTRI2,IH)                            CPP31940
      DO 102 ITYP=1,NTYPES                                              CPP31950
      FAC=FOCC(ITYP)*HALF                                               CPP31960
      DO 101 I=1,NTRI                                                   CPP31970
  101 HM(I,IXYZ,ITYP)=HF(I,IXYZ)*FAC                                    CPP31980
  102 CONTINUE                                                          CPP31990
  103 CONTINUE                                                          CPP32000
C                                                                       CPP32010
C   STORE DERIVATIVE LAGRANGIAN MATRICES                                CPP32020
      DO 110 IXYZ=1,3                                                   CPP32030
      IE=IEA(IXYZ+N3N)                                                  CPP32040
      IB=IBA(IXYZ+N3N)                                                  CPP32050
      CALL ZERO(EPF,NBASIS*NBASIS)                                      CPP32060
C                                                                       CPP32070
      DO 107 ITYP=1,NTYPES                                              CPP32080
      DO 106 I=NSTR(ITYP),NEND(ITYP)                                    CPP32090
      DO 106 J=1,NBASIS                                                 CPP32100
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPP32110
      EPF(I,J)=HM(IJ,IXYZ,ITYP)                                         CPP32120
  106 CONTINUE                                                          CPP32130
  107 CONTINUE                                                          CPP32140
C                                                                       CPP32150
C   ELEMENTS FOR INDEPENDENT PAIRS                                      CPP32160
      DO 108 II=1,NIND                                                  CPP32170
      I=NIJ(II,1)                                                       CPP32180
      J=NIJ(II,2)                                                       CPP32190
      IJ=IOFF(I)+J                                                      CPP32200
      BF(IJ)=EPF(I,J)-EPF(J,I)                                          CPP32210
  108 CONTINUE                                                          CPP32220
C                                                                       CPP32230
C   ELEMENTS FOR DEPENDENT PAIRS                                        CPP32240
      DO 109 II=1,NDEP                                                  CPP32250
      I=KIJ(II,1)                                                       CPP32260
      J=KIJ(II,2)                                                       CPP32270
      IJ=IOFF(I)+J                                                      CPP32280
      BF(IJ)=HF(IJ,IXYZ)                                                CPP32290
  109 CONTINUE                                                          CPP32300
C                                                                       CPP32310
      CALL RWRIT(ITAP44,BF,NTRI2,IB)                                    CPP32320
      CALL RWRIT(ITAP44,EPF,NBASQ,IE)                                   CPP32330
      IF(IPRNT.LE.4) GO TO 110                                          CPP32340
      WRITE(6,4) IXYZ                                                   CPP32350
      CALL MATOUT(EPF,NBASIS,NBASIS,NBASIS,NBASIS,6)                    CPP32360
  110 CONTINUE                                                          CPP32370
C                                                                       CPP32380
C   FORM THE HIJF MATRIX                                                CPP32390
      DO 115 IJCI=1,IIUNQ                                               CPP32400
      IEX=IPX(IJCI)                                                     CPP32410
      JEX=JPX(IJCI)                                                     CPP32420
      IJEX=IOFF(MAX0(IEX,JEX))+MIN0(IEX,JEX)                            CPP32430
      DO 113 IXYZ=1,3                                                   CPP32440
C                                                                       CPP32450
      ELECF=A00                                                         CPP32460
      DO 112 I=1,JOCC                                                   CPP32470
      II=IOFF(I+1)                                                      CPP32480
      ITYP=MOTYP(I)                                                     CPP32490
      ELECF=ELECF+FIJ(ITYP,IJCI)*HF(II,IXYZ)*TWO                        CPP32500
  112 CONTINUE                                                          CPP32510
      HCIA(IJEX,IXYZ+N3N)=ELECF                                         CPP32520
  113 CONTINUE                                                          CPP32530
      IF(IPRNT.LE.2) GO TO 115                                          CPP32540
      WRITE(6,5) IJEX                                                   CPP32550
      CALL PRINT(HCIA(IJEX,IXYZ+N3N),NCTRI,NCPEX,6)                     CPP32560
  115 CONTINUE                                                          CPP32570
C                                                                       CPP32580
C   CALCULATE DIPOLE MOMENTS FOR A TEST                                 CPP32590
      DO 120 IXYZ=1,3                                                   CPP32600
      VALU=A00                                                          CPP32610
      DO 118 IJCI=1,ICIUNQ                                              CPP32620
      IEX=IPX(IJCI)                                                     CPP32630
      JEX=JPX(IJCI)                                                     CPP32640
      IJEX=IOFF(MAX0(IEX,JEX))+MIN0(IEX,JEX)                            CPP32650
      FAC=CI(IEX)*CI(JEX)*TWO                                           CPP32660
      IF(IEX.EQ.JEX) FAC=FAC*HALF                                       CPP32670
      VALU=VALU+HCIA(IJEX,IXYZ+N3N)*FAC                                 CPP32680
  118 CONTINUE                                                          CPP32690
      DIP(IXYZ)=-VALU*DEBYE                                             CPP32700
      E1T(IXYZ+N3N)=VALU                                                CPP32710
  120 CONTINUE                                                          CPP32720
      WRITE(6,6)                                                        CPP32730
      CALL MATOUT(DIP,3,1,3,1,6)                                        CPP32740
C                                                                       CPP32750
      CALL RCLOSE(ITAP43,3)                                             CPP32760
      CALL RCLOSE(ITAP44,3)                                             CPP32770
C                                                                       CPP32780
      RETURN                                                            CPP32790
      END                                                               CPP32800
      SUBROUTINE AB2CI(NIJ,E1T,HCIA,HCI,BAF2,A22,A12,SA,EIJ,NABC)       CPP32810
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP32820
      DIMENSION NIJ(NTRI,2),SA(NTRI,N3N)                                CPP32830
      DIMENSION E1T(NABC),HCIA(NCTRI,NABC),HCI(NCTRI)                   CPP32840
      DIMENSION BAF2(NCPEX,NABC),A22(NCPEX,NCPEX),A12(NIND,NCPEX)       CPP32850
      DIMENSION EIJ(NBASIS,NBASIS)                                      CPP32860
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP32870
      COMMON/CIOFF/IIUNQ,IJUNQ,ICIUNQ,IPX(5050),JPX(5050)               CPP32880
      COMMON/CIPOS/ICITYP(5050),ICPOS(200,100)                          CPP32890
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP32900
      COMMON/ENRGY/ESCF,ENUC,ELEC                                       CPP32910
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP32920
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP32930
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP32940
      COMMON/TAPES/ITAP37,ITAP42,ITAP43,ITAP44,ITAP47                   CPP32950
      COMMON/WORKS/JTAPE1,JTAPE2                                        CPP32960
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP32970
      COMMON/CI102/NIND,NDEP                                            CPP32980
      COMMON/CI103/CI(100)                                              CPP32990
      COMMON/CI104/IEIJ(5050)                                           CPP33000
      DATA A00,HALF,TWO,FOUR /  0.0D+00 , 0.5D+00 , 2.0D+00 , 4.0D+00 / CPP33010
    1 FORMAT(//,2X,' BAF2 MATRIX'/)                                     CPP33020
    2 FORMAT(//,2X,' A12 MATRIX'/)                                      CPP33030
    3 FORMAT(//,2X,' A22 MATRIX'/)                                      CPP33040
C                                                                       CPP33050
C   CALL BACK SA MATRICES                                               CPP33060
      CALL RFILE(ITAP44)                                                CPP33070
      DO 101 IABC=1,N3N                                                 CPP33080
      IS=ISA(IABC)                                                      CPP33090
      CALL RREAD(ITAP44,SA(1,IABC),NTRI2,IS)                            CPP33100
  101 CONTINUE                                                          CPP33110
C                                                                       CPP33120
C   FORM THE BAF AND A12 MATRICES                                       CPP33130
      CALL ZERO(BAF2,NCPEX*NABC)                                        CPP33140
      CALL ZERO(A12,NIND*NCPEX)                                         CPP33150
      CALL RFILE(JTAPE2)                                                CPP33160
C                                                                       CPP33170
      DO 110 ICI=1,NCPEX                                                CPP33180
      DO 102 IABC=1,NABC                                                CPP33190
      BAF2(ICI,IABC)=E1T(IABC)*CI(ICI)*HALF                             CPP33200
  102 CONTINUE                                                          CPP33210
      DO 109 JCI=1,NCPEX                                                CPP33220
      IJCI=IOFF(MAX0(ICI,JCI))+MIN0(ICI,JCI)                            CPP33230
      DO 103 III=1,ICIUNQ                                               CPP33240
      IEX=IPX(III)                                                      CPP33250
      JEX=JPX(III)                                                      CPP33260
      IJEX=IOFF(MAX0(IEX,JEX))+MIN0(IEX,JEX)                            CPP33270
      IF(IJCI.EQ.IJEX) GO TO 301                                        CPP33280
  103 CONTINUE                                                          CPP33290
      GO TO 109                                                         CPP33300
C                                                                       CPP33310
C   CUMMULATE ONLY NON-ZERO CI'S                                        CPP33320
  301 CONTINUE                                                          CPP33330
      IBCI=ICITYP(III)                                                  CPP33340
      IPT=IEIJ(IBCI)                                                    CPP33350
      CALL RREAD(JTAPE2,EIJ,NBASQ,IPT)                                  CPP33360
C                                                                       CPP33370
      DO 105 IABC=1,NABC                                                CPP33380
      VALU0=A00                                                         CPP33390
      IF(IABC.GT.N3N) GO TO 201                                         CPP33400
      DO 104 I=1,JOCC                                                   CPP33410
      DO 104 J=1,I                                                      CPP33420
      IJ=IOFF(I)+J                                                      CPP33430
      FAC4=FOUR                                                         CPP33440
      IF(I.EQ.J) FAC4=TWO                                               CPP33450
      VALU0=VALU0+SA(IJ,IABC)*EIJ(I,J)*FAC4                             CPP33460
  104 CONTINUE                                                          CPP33470
  201 CONTINUE                                                          CPP33480
      VALU2=(-HCIA(IJCI,IABC)+VALU0)*CI(JCI)                            CPP33490
      BAF2(ICI,IABC)=BAF2(ICI,IABC)+VALU2*HALF                          CPP33500
  105 CONTINUE                                                          CPP33510
C                                                                       CPP33520
      DO 107 II=1,NIND                                                  CPP33530
      I=NIJ(II,1)                                                       CPP33540
      J=NIJ(II,2)                                                       CPP33550
      VALU=(EIJ(I,J)-EIJ(J,I))*CI(JCI)                                  CPP33560
      A12(II,ICI)=A12(II,ICI)-VALU*TWO                                  CPP33570
  107 CONTINUE                                                          CPP33580
C                                                                       CPP33590
  109 CONTINUE                                                          CPP33600
  110 CONTINUE                                                          CPP33610
      IF(IPRNT.LE.2) GO TO 202                                          CPP33620
      WRITE(6,1)                                                        CPP33630
      CALL MATOUT(BAF2,NCPEX,NABC,NCPEX,NABC,6)                         CPP33640
      WRITE(6,2)                                                        CPP33650
      CALL MATOUT(A12,NIND,NCPEX,NIND,NCPEX,6)                          CPP33660
C                                                                       CPP33670
C   FORM THE A22 MATRIX                                                 CPP33680
  202 CONTINUE                                                          CPP33690
      CALL ZERO(A22,NCPEX*NCPEX)                                        CPP33700
      IJ=0                                                              CPP33710
      DO 115 I=1,NCPEX                                                  CPP33720
      DO 115 J=1,I                                                      CPP33730
      IJ=IJ+1                                                           CPP33740
      VALU=HCI(IJ)+CI(I)*CI(J)                                          CPP33750
      IF(I.EQ.J) VALU=VALU-ELEC                                         CPP33760
      A22(I,J)=VALU*HALF                                                CPP33770
      A22(J,I)=A22(I,J)                                                 CPP33780
  115 CONTINUE                                                          CPP33790
      IF(IPRNT.LE.2) GO TO 203                                          CPP33800
      WRITE(6,3)                                                        CPP33810
      CALL MATOUT(A22,NCPEX,NCPEX,NCPEX,NCPEX,6)                        CPP33820
C                                                                       CPP33830
  203 CONTINUE                                                          CPP33840
      CALL RCLOSE(JTAPE2,3)                                             CPP33850
      CALL RCLOSE(ITAP44,3)                                             CPP33860
      RETURN                                                            CPP33870
      END                                                               CPP33880
      SUBROUTINE CPHF(ESO,NIJ,KIJ,ZETA,BAF2,CIA,A22,A12,A,U,T,B0,AA,TT, CPP33890
     1                DP,DQ,DM,W,BB,MIJ,PQ,NINDPX,NMAX,NADD,N3NXX)      CPP33900
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP33910
      REAL*16 DLIM,DETERM                                               CPP33920
      INTEGER*2 MIJ                                                     CPP33930
      DIMENSION ESO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)         CPP33940
      DIMENSION A(NBFAO,NBFAO),ZETA(NTRI,NTYPEP)                        CPP33950
      DIMENSION BB(NINDPX,NMAX*2),W(NMAX,NMAX+NADD)                     CPP33960
      DIMENSION DP(NTRI,NADD*NTYPEP),DQ(NTRI,NADD*NTYPEP)               CPP33970
      DIMENSION DM(NTRI,NADD*NTYPEP)                                    CPP33980
      DIMENSION AA(NINDPX),TT(NINDPX),B0(NTRI),NIJ(NTRI,2),KIJ(NTRI,2)  CPP33990
      DIMENSION BAF2(NCPEX,N3NXX),CIA(NCPEX,N3NXX)                      CPP34000
      DIMENSION A22(NCPEX,NCPEX),A12(NIND,NCPEX)                        CPP34010
      DIMENSION PQ(*)                                                   CPP34020
      DIMENSION NDUM(1024),BNORM(1024)                                  CPP34030
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP34040
      COMMON/CALIF/LPARA(1024),ANORM(1024)                              CPP34050
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP34060
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP34070
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP34080
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP34090
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP34100
      COMMON/TAPES/ITAP37,ITAP42,ITAP43,ITAP44,ITAP47                   CPP34110
      COMMON/TOLER/CRIT,TOLR                                            CPP34120
      COMMON/WORKS/JTAPE1,JTAPE2                                        CPP34130
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP34140
      COMMON/CI102/NIND,NDEP                                            CPP34150
      DATA A00,HALF,ONE,TWO,FOUR,TEN / 0.0D+00 , 0.5D+00 , 1.0D+00 ,    CPP34160
     1                                 2.0D+00 , 4.0D+00 , 10.0D+00 /   CPP34170
    1 FORMAT(//,2X,' A MATRIX'/)                                        CPP34180
    2 FORMAT(//,2X,' IGROUP = ',I5,5X,' IFIRST = ',I5,5X,               CPP34190
     1             ' NABC   = ',I5/)                                    CPP34200
    3 FORMAT(//,2X,' BB MATRIX'/)                                       CPP34210
    4 FORMAT(//,2X,' SCALED BB MATRIX'/)                                CPP34220
    5 FORMAT(//,2X,' ORTHONORMALIZED BB MATRIX'/)                       CPP34230
    6 FORMAT(//,2X,' ORTHONORMALITY OF BB MATRIX'/)                     CPP34240
    7 FORMAT(//,2X,' ITERATION NO. = ',I5/)                             CPP34250
    8 FORMAT(//,2X,' (1-A)*B MATRIX'/)                                  CPP34260
    9 FORMAT(//,2X,' KVEC = ',I5,' K = ',I5,' TT = ',(10F10.5))         CPP34270
   10 FORMAT(/,2X,' NITER  = ',I5,5X,' NVEC   = ',I5,5X,' NADDV  = ',I5/CPP34280
     1)                                                                 CPP34290
   11 FORMAT(//,2X,' NORM OF ADDED VECTOR',(5(I5,F15.8)))               CPP34300
   12 FORMAT(//,2X,' BB MATRIX'/)                                       CPP34310
   13 FORMAT(/,2X,' END OF ITERATION    NVTOT  = ',I5/)                 CPP34320
   14 FORMAT(//,2X,'***NVTOT EXCEEDS NMAX***',5X,' NMAX = ',I5/)        CPP34330
   15 FORMAT(//,2X,' EXPANDED VECTORS BB AND A*BB MATRICES'/)           CPP34340
   16 FORMAT(//,2X,' W MATRIX'/)                                        CPP34350
   17 FORMAT(//,2X,' SCALED W MATRIX'/)                                 CPP34360
   18 FORMAT(//,2X,'***SIMULTANEOUS EQUATION CANNOT BE SOLVED***',      CPP34370
     1 5X,' DETERM = ',Q20.10/)                                         CPP34380
   19 FORMAT(/,2X,' LINEAR EQUATION HAS BEEN SOLVED',                   CPP34390
     1 5X,' DETERM = ',Q20.10/)                                         CPP34400
   20 FORMAT(//,2X,' SOLUTION MATRIX OF SCALED W MATRIX'/)              CPP34410
   21 FORMAT(//,2X,' SOLUTION MATRIX'/)                                 CPP34420
   22 FORMAT(2X,2I7,F15.8)                                              CPP34430
C                                                                       CPP34440
      DLIM=1.0Q-76                                                      CPP34450
      CALL RFILE(JTAPE1)                                                CPP34460
      CALL RFILE(ITAP44)                                                CPP34470
C                                                                       CPP34480
      NINDPX=NIND+NCPEX                                                 CPP34490
      NINP2=NINDPX*2                                                    CPP34500
C                                                                       CPP34510
C   CALL BACK LAGRANGIAN AND ZETA MATRICES                              CPP34520
      DO 101 ITYP=1,NTYPES                                              CPP34530
      CALL RMREAD(ZETA(1,ITYP),35+ITYP)                                 CPP34540
C     CALL PRINT(ZETA(1,ITYP),NTRI,NBASIS,6)                            CPP34550
  101 CONTINUE                                                          CPP34560
      CALL ZERO(ZETA(1,NTYPEP),NTRI)                                    CPP34570
C                                                                       CPP34580
C   FORM A MATRIX                                                       CPP34590
C   SEE EQ.(A-1)                                                        CPP34600
C                                                                       CPP34610
      CALL ZERO(A,NBASIS*NBASIS)                                        CPP34620
      DO 102 I=1,NBASIS                                                 CPP34630
      DO 102 J=1,NBASIS                                                 CPP34640
      ITYP=MOTYP(I)                                                     CPP34650
      JTYP=MOTYP(J)                                                     CPP34660
      II=IOFF(I+1)                                                      CPP34670
      JJ=IOFF(J+1)                                                      CPP34680
      A(I,J)=A(I,J)-ZETA(II,ITYP)+ZETA(II,JTYP)                         CPP34690
     1             -ZETA(JJ,JTYP)+ZETA(JJ,ITYP)                         CPP34700
  102 CONTINUE                                                          CPP34710
      IF(IPRNT.LE.2) GO TO 301                                          CPP34720
      WRITE(6,1)                                                        CPP34730
      CALL MATOUT(A,NBFAO,NBFAO,NBASIS,NBASIS,6)                        CPP34740
C                                                                       CPP34750
C   FORM SCALING VECTOR                                                 CPP34760
C   SEE EQS.(24) AND (25)                                               CPP34770
C                                                                       CPP34780
  301 CONTINUE                                                          CPP34790
C                                                                       CPP34800
C:::ORBITAL SPACE:::                                                    CPP34810
      DO 103 II=1,NIND                                                  CPP34820
      I=NIJ(II,1)                                                       CPP34830
      J=NIJ(II,2)                                                       CPP34840
      AA(II)=ONE/DSQRT(DABS(A(I,J)))                                    CPP34850
  103 CONTINUE                                                          CPP34860
C                                                                       CPP34870
C:::CISPACE:::                                                          CPP34880
      DO 104 ICI=1,NCPEX                                                CPP34890
      AA(NIND+ICI)=ONE/DSQRT(DABS(A22(ICI,ICI)))                        CPP34900
  104 CONTINUE                                                          CPP34910
C                                                                       CPP34920
C   START GROUPING THE B VECTORS                                        CPP34930
      IGROUP=0                                                          CPP34940
      ILAST=0                                                           CPP34950
  500 CONTINUE                                                          CPP34960
      IGROUP=IGROUP+1                                                   CPP34970
      IFIRST=ILAST+1                                                    CPP34980
      ILAST=IFIRST+NADD-1                                               CPP34990
      IF(ILAST.GT.N3NXX) ILAST=N3NXX                                    CPP35000
      NABC=ILAST-IFIRST+1                                               CPP35010
      WRITE(3,2) IGROUP,IFIRST,NABC                                     CPP35020
      WRITE(4,2) IGROUP,IFIRST,NABC                                     CPP35030
      WRITE(6,2) IGROUP,IFIRST,NABC                                     CPP35040
      CALL SREW(JTAPE1)                                                 CPP35050
C                                                                       CPP35060
C   CALL BACK B0 MATRICES                                               CPP35070
      CALL SREW(ITAP44)                                                 CPP35080
      DO 107 IABC=1,NABC                                                CPP35090
      III=IABC+IFIRST-1                                                 CPP35100
C                                                                       CPP35110
C:::ORBITAL SPACE:::                                                    CPP35120
      IB=IBA(III)                                                       CPP35130
      CALL RREAD(ITAP44,B0,NTRI2,IB)                                    CPP35140
      DO 105 II=1,NIND                                                  CPP35150
      I=NIJ(II,1)                                                       CPP35160
      J=NIJ(II,2)                                                       CPP35170
      IJ=IOFF(I)+J                                                      CPP35180
      BB(II,IABC)=B0(IJ)                                                CPP35190
  105 CONTINUE                                                          CPP35200
C                                                                       CPP35210
C:::CI SPACE:::                                                         CPP35220
      DO 106 ICI=1,NCPEX                                                CPP35230
      BB(NIND+ICI,IABC)=BAF2(ICI,III)                                   CPP35240
  106 CONTINUE                                                          CPP35250
C                                                                       CPP35260
  107 CONTINUE                                                          CPP35270
      IF(IPRNT.LE.2) GO TO 302                                          CPP35280
      WRITE(6,3)                                                        CPP35290
      CALL MATOUT(BB,NINDPX,NABC,NINDPX,NABC,6)                         CPP35300
  302 CONTINUE                                                          CPP35310
C                                                                       CPP35320
      DO 109 II=1,NINDPX                                                CPP35330
      DO 108 IABC=1,NABC                                                CPP35340
  108 BB(II,IABC)=BB(II,IABC)*AA(II)                                    CPP35350
  109 CONTINUE                                                          CPP35360
      IF(IPRNT.LE.2) GO TO 303                                          CPP35370
      WRITE(6,4)                                                        CPP35380
      CALL MATOUT(BB,NINDPX,NABC,NINDPX,NABC,6)                         CPP35390
C                                                                       CPP35400
C   ORTHONORMALIZE B0 MATRIX                                            CPP35410
C                                                                       CPP35420
  303 CONTINUE                                                          CPP35430
      CALL ORTHOG(BB,BB(1,NABC+1),NINDPX,NABC,NINDPX,NDUM,NABC,NVEC)    CPP35440
      IF(IPRNT.LE.2) GO TO 304                                          CPP35450
      WRITE(6,5)                                                        CPP35460
      CALL MATOUT(BB,NINDPX,NVEC,NINDPX,NVEC,6)                         CPP35470
  304 CONTINUE                                                          CPP35480
      DO 110 IABC=1,NVEC                                                CPP35490
  110 BNORM(IABC)=ONE                                                   CPP35500
C                                                                       CPP35510
C   CHECK ORTHONORMALITY                                                CPP35520
      IF(IPRNT.LE.2) GO TO 305                                          CPP35530
      DO 112 I=1,NVEC                                                   CPP35540
      DO 112 J=1,NVEC                                                   CPP35550
      VALU=A00                                                          CPP35560
      DO 111 K=1,NINDPX                                                 CPP35570
  111 VALU=VALU+BB(K,I)*BB(K,J)                                         CPP35580
      W(I,J)=VALU                                                       CPP35590
  112 CONTINUE                                                          CPP35600
      WRITE(6,6)                                                        CPP35610
      CALL MATOUT(W,NMAX,NMAX+NABC,NVEC,NVEC,6)                         CPP35620
C                                                                       CPP35630
  305 NPVEC=0                                                           CPP35640
      NITER=0                                                           CPP35650
C:::::::::::::::::::::                                                  CPP35660
C:::START ITERATION:::                                                  CPP35670
C:::::::::::::::::::::                                                  CPP35680
  200 NITER=NITER+1                                                     CPP35690
      IF(IPRNT.GT.2)                                                    CPP35700
     *WRITE(6,7) NITER                                                  CPP35710
C                                                                       CPP35720
      DO 114 II=1,NINDPX                                                CPP35730
      DO 113 IABC=1,NVEC                                                CPP35740
  113 BB(II,NPVEC+IABC)=BB(II,NPVEC+IABC)*AA(II)                        CPP35750
  114 CONTINUE                                                          CPP35760
C                                                                       CPP35770
C   FORM A*B MATRIX                                                     CPP35780
      CALL MAKEAB(BB(1,NPVEC+1),DP,DQ,DM,BB(1,NPVEC+NVEC+1),ZETA,       CPP35790
     1            NIJ,A22,A12,ESO,U,T,MIJ,PQ,NINDPX,NVEC)               CPP35800
C                                                                       CPP35810
      DO 116 II=1,NINDPX                                                CPP35820
      DO 115 IABC=1,NVEC                                                CPP35830
      BB(II,NPVEC+IABC)=BB(II,NPVEC+IABC)/AA(II)                        CPP35840
      BB(II,NPVEC+NVEC+IABC)=BB(II,NPVEC+NVEC+IABC)*AA(II)              CPP35850
  115 CONTINUE                                                          CPP35860
  116 CONTINUE                                                          CPP35870
C                                                                       CPP35880
C   STORE A*B MATRIX                                                    CPP35890
      DO 117 IABC=NPVEC+NVEC+1,NPVEC+NVEC*2                             CPP35900
      CALL SWRIT(JTAPE1,BB(1,IABC),NINP2)                               CPP35910
  117 CONTINUE                                                          CPP35920
C                                                                       CPP35930
C   FORM (1-A)*B                                                        CPP35940
      DO 120 IABC=NPVEC+1,NPVEC+NVEC                                    CPP35950
      DO 119 II=1,NINDPX                                                CPP35960
C*119 BB(II,NVEC+IABC)=-BB(II,IABC)-BB(II,NVEC+IABC)                    CPP35970
  119 BB(II,NVEC+IABC)=BB(II,IABC)-BB(II,NVEC+IABC)                     CPP35980
  120 CONTINUE                                                          CPP35990
      IF(IPRNT.LE.2) GO TO 306                                          CPP36000
      WRITE(6,8)                                                        CPP36010
      CALL MATOUT(BB(1,NPVEC+NVEC+1),NINDPX,NVEC,NINDPX,NVEC,6)         CPP36020
C                                                                       CPP36030
C   FORM ADDITIONAL EXPANSION VECTORS                                   CPP36040
  306 NADDV=0                                                           CPP36050
      K=NPVEC+NVEC                                                      CPP36060
      DO 130 IABC=NPVEC+NVEC+1,NPVEC+NVEC*2                             CPP36070
C                                                                       CPP36080
      DO 121 II=1,NINDPX                                                CPP36090
  121 TT(II)=BB(II,IABC)                                                CPP36100
      DO 124 L=1,K                                                      CPP36110
      CALL NORM(XBB,BB(1,L),BB(1,IABC),NINDPX)                          CPP36120
      FAC=XBB/BNORM(L)                                                  CPP36130
      DO 123 II=1,NINDPX                                                CPP36140
  123 TT(II)=TT(II)-BB(II,L)*FAC                                        CPP36150
  124 CONTINUE                                                          CPP36160
      IF(IPRNT.LE.2) GO TO 307                                          CPP36170
      WRITE(6,9) IABC,K,(TT(I),I=1,NINDPX)                              CPP36180
  307 CONTINUE                                                          CPP36190
      CALL NORM(XTT,TT,TT,NINDPX)                                       CPP36200
      IF(DABS(XTT).LT.CRIT) GO TO 130                                   CPP36210
      NADDV=NADDV+1                                                     CPP36220
      K=K+1                                                             CPP36230
      BNORM(K)=XTT                                                      CPP36240
      DO 126 II=1,NINDPX                                                CPP36250
C*126 BB(II,K)=-TT(II)                                                  CPP36260
  126 BB(II,K)=TT(II)                                                   CPP36270
C                                                                       CPP36280
  130 CONTINUE                                                          CPP36290
C                                                                       CPP36300
      WRITE(6,10) NITER,NVEC,NADDV                                      CPP36310
      WRITE(3,10) NITER,NVEC,NADDV                                      CPP36320
      WRITE(4,10) NITER,NVEC,NADDV                                      CPP36330
      NPVEC=NPVEC+NVEC                                                  CPP36340
      NVEC=NADDV                                                        CPP36350
      IF(NADDV.EQ.0) GO TO 310                                          CPP36360
      IF(IPRNT.LE.2) GO TO 308                                          CPP36370
      IST=NPVEC+1                                                       CPP36380
      IED=NPVEC+NVEC                                                    CPP36390
      WRITE(6,11) (I,BNORM(I),I=IST,IED)                                CPP36400
      WRITE(6,12)                                                       CPP36410
      CALL MATOUT(BB(1,NPVEC+1),NINDPX,NVEC,NINDPX,NVEC,6)              CPP36420
  308 GO TO 200                                                         CPP36430
C                                                                       CPP36440
  310 NVTOT=NPVEC                                                       CPP36450
      WRITE(6,13) NVTOT                                                 CPP36460
      WRITE(3,13) NVTOT                                                 CPP36470
      WRITE(4,13) NVTOT                                                 CPP36480
      IF(NVTOT.LE.NMAX) GO TO 400                                       CPP36490
      WRITE(6,14) NMAX                                                  CPP36500
      STOP                                                              CPP36510
C                                                                       CPP36520
C   CALL BACK A*B MATRICES                                              CPP36530
  400 CONTINUE                                                          CPP36540
      CALL SREW(JTAPE1)                                                 CPP36550
      DO 131 IABC=NVTOT+1,NVTOT*2                                       CPP36560
  131 CALL SREAD(JTAPE1,BB(1,IABC),NINP2)                               CPP36570
      IF(IPRNT.LE.2) GO TO 311                                          CPP36580
      WRITE(6,15)                                                       CPP36590
      CALL MATOUT(BB,NINDPX,NVTOT*2,NINDPX,NVTOT*2,6)                   CPP36600
C                                                                       CPP36610
  311 CONTINUE                                                          CPP36620
      DO 133 I=1,NVTOT                                                  CPP36630
      DO 133 J=1,NVTOT                                                  CPP36640
      SUM=A00                                                           CPP36650
      DO 132 II=1,NINDPX                                                CPP36660
C*132 SUM=SUM-BB(II,I)*BB(II,NVTOT+J)                                   CPP36670
  132 SUM=SUM+BB(II,I)*BB(II,NVTOT+J)                                   CPP36680
      W(I,J)=SUM                                                        CPP36690
  133 CONTINUE                                                          CPP36700
C                                                                       CPP36710
C   READ BACK ORIGINAL B0 MATRICES                                      CPP36720
      CALL SREW(ITAP44)                                                 CPP36730
      DO 136 IABC=1,NABC                                                CPP36740
      III=IABC+IFIRST-1                                                 CPP36750
      IB=IBA(III)                                                       CPP36760
      CALL RREAD(ITAP44,B0,NTRI2,IB)                                    CPP36770
      DO 134 II=1,NIND                                                  CPP36780
      I=NIJ(II,1)                                                       CPP36790
      J=NIJ(II,2)                                                       CPP36800
      IJ=IOFF(I)+J                                                      CPP36810
      BB(II,NVTOT+IABC)=B0(IJ)*AA(II)                                   CPP36820
  134 CONTINUE                                                          CPP36830
      DO 135 ICI=1,NCPEX                                                CPP36840
      BB(NIND+ICI,NVTOT+IABC)=BAF2(ICI,III)*AA(NIND+ICI)                CPP36850
  135 CONTINUE                                                          CPP36860
  136 CONTINUE                                                          CPP36870
      IF(IPRNT.LE.2) GO TO 312                                          CPP36880
      WRITE(6,4)                                                        CPP36890
      CALL MATOUT(BB(1,NVTOT+1),NINDPX,NABC,NINDPX,NABC,6)              CPP36900
C                                                                       CPP36910
C   SET B0*B ON W                                                       CPP36920
  312 CONTINUE                                                          CPP36930
      DO 138 I=1,NVTOT                                                  CPP36940
      DO 138 J=1,NABC                                                   CPP36950
      SUM=A00                                                           CPP36960
      DO 137 II=1,NINDPX                                                CPP36970
  137 SUM=SUM+BB(II,I)*BB(II,NVTOT+J)                                   CPP36980
      W(I,NVTOT+J)=SUM                                                  CPP36990
  138 CONTINUE                                                          CPP37000
C                                                                       CPP37010
      IF(IPRNT.LE.2) GO TO 313                                          CPP37020
      WRITE(6,16)                                                       CPP37030
      CALL MATOUT(W,NMAX,NMAX+NABC,NVTOT,NVTOT+NABC,6)                  CPP37040
  313 CONTINUE                                                          CPP37050
      DO 139 I=1,NVTOT                                                  CPP37060
      AVAL=ONE/DSQRT(BNORM(I))                                          CPP37070
      DO 139 J=1,NVTOT                                                  CPP37080
      BVAL=ONE/DSQRT(BNORM(J))                                          CPP37090
      W(I,J)=W(I,J)*AVAL*BVAL                                           CPP37100
  139 CONTINUE                                                          CPP37110
      DO 140 I=1,NVTOT                                                  CPP37120
      AVAL=ONE/DSQRT(BNORM(I))                                          CPP37130
      DO 140 J=1,NABC                                                   CPP37140
      W(I,NVTOT+J)=W(I,NVTOT+J)*AVAL                                    CPP37150
  140 CONTINUE                                                          CPP37160
      IF(IPRNT.LE.3) GO TO 315                                          CPP37170
      WRITE(6,17)                                                       CPP37180
      CALL MATOUT(W,NMAX,NMAX+NABC,NVTOT,NVTOT+NABC,6)                  CPP37190
C                                                                       CPP37200
C   SOLVE SMALL SIMULTANEOUS EQUATION                                   CPP37210
  315 CONTINUE                                                          CPP37220
      CALL FLINQ(W,NMAX,NVTOT,NABC,DETERM)                              CPP37230
      IF(QABS(DETERM).GT.DLIM) GO TO 320                                CPP37240
C                                                                       CPP37250
C   SIMULTANEOUS EQUATION CANNOT BE SOLVED---STOP!!!                    CPP37260
      WRITE(6,18) DETERM                                                CPP37270
      STOP                                                              CPP37280
C                                                                       CPP37290
C::::::::::::::::::::::::::::::::                                       CPP37300
C:::LINEAR EQUATION WAS SOLVED:::                                       CPP37310
C::::::::::::::::::::::::::::::::                                       CPP37320
  320 CONTINUE                                                          CPP37330
      WRITE(6,19) DETERM                                                CPP37340
      WRITE(3,19) DETERM                                                CPP37350
      WRITE(4,19) DETERM                                                CPP37360
      IF(IPRNT.LE.2) GO TO 322                                          CPP37370
      WRITE(6,20)                                                       CPP37380
      CALL MATOUT(W(1,NVTOT+1),NMAX,NMAX+NABC,NVTOT,NABC,6)             CPP37390
  322 CONTINUE                                                          CPP37400
      DO 151 I=1,NVTOT                                                  CPP37410
      AVAL=ONE/DSQRT(BNORM(I))                                          CPP37420
      DO 151 J=1,NABC                                                   CPP37430
      W(I,NVTOT+J)=W(I,NVTOT+J)*AVAL                                    CPP37440
  151 CONTINUE                                                          CPP37450
      IF(IPRNT.LE.4) GO TO 323                                          CPP37460
      WRITE(6,21)                                                       CPP37470
      CALL MATOUT(W(1,NVTOT+1),NMAX,NMAX+NABC,NVTOT,NABC,6)             CPP37480
C                                                                       CPP37490
C   OBTAIN INDEPENDENT ELEMENTS OF UA MATRICES                          CPP37500
  323 CONTINUE                                                          CPP37510
      DO 153 II=1,NINDPX                                                CPP37520
      DO 153 IABC=1,NABC                                                CPP37530
      SUM=A00                                                           CPP37540
      DO 152 I=1,NVTOT                                                  CPP37550
C*152 SUM=SUM-BB(II,I)*W(I,NVTOT+IABC)                                  CPP37560
  152 SUM=SUM+BB(II,I)*W(I,NVTOT+IABC)                                  CPP37570
      BB(II,NVTOT+IABC)=SUM*AA(II)                                      CPP37580
  153 CONTINUE                                                          CPP37590
C                                                                       CPP37600
      CALL SREW(ITAP44)                                                 CPP37610
      DO 156 IABC=1,NABC                                                CPP37620
      III=IABC+IFIRST-1                                                 CPP37630
C                                                                       CPP37640
C:::ORBITAL SPACE:::                                                    CPP37650
      IB=IBA(III)                                                       CPP37660
      CALL RREAD(ITAP44,B0,NTRI2,IB)                                    CPP37670
      DO 154 II=1,NIND                                                  CPP37680
      I=NIJ(II,1)                                                       CPP37690
      J=NIJ(II,2)                                                       CPP37700
      IJ=IOFF(I)+J                                                      CPP37710
      B0(IJ)=BB(II,NVTOT+IABC)                                          CPP37720
  154 CONTINUE                                                          CPP37730
      CALL RWRIT(ITAP44,B0,NTRI2,IB)                                    CPP37740
C                                                                       CPP37750
C:::CI SPACE:::                                                         CPP37760
      DO 155 ICI=1,NCPEX                                                CPP37770
      CIA(ICI,III)=BB(NIND+ICI,NVTOT+IABC)                              CPP37780
      IF(IPRNT.LE.2) GO TO 155                                          CPP37790
      WRITE(6,22) ICI,III,CIA(ICI,III)                                  CPP37800
  155 CONTINUE                                                          CPP37810
  156 CONTINUE                                                          CPP37820
C                                                                       CPP37830
      LPARA(101)=10                                                     CPP37840
      LPARA(102)=IGROUP                                                 CPP37850
      LPARA(103)=ILAST                                                  CPP37860
      CALL RMWRIT(LPARA,2)                                              CPP37870
      IF(ILAST.GE.N3NXX) GO TO 501                                      CPP37880
      GO TO 500                                                         CPP37890
C                                                                       CPP37900
  501 CONTINUE                                                          CPP37910
      CALL RCLOSE(ITAP44,3)                                             CPP37920
      CALL RCLOSE(JTAPE1,3)                                             CPP37930
      RETURN                                                            CPP37940
      END                                                               CPP37950
      SUBROUTINE MAKEAB(B0,DP,DQ,DM,B1,ZETA,NIJ,A22,A12,ESO,U,T,MIJ,PQ, CPP37960
     1                  NINDPX,NABC)                                    CPP37970
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP37980
      INTEGER*2 MIJ(*)                                                  CPP37990
      DIMENSION B0(NINDPX,NABC),B1(NINDPX,NABC)                         CPP38000
      DIMENSION ZETA(NTRI,NTYPEP),NIJ(NTRI,2)                           CPP38010
      DIMENSION A22(NCPEX,NCPEX),A12(NIND,NCPEX)                        CPP38020
      DIMENSION DP(NTRI,NABC*NTYPEP),DQ(NTRI,NABC*NTYPEP)               CPP38030
      DIMENSION DM(NTRI,NABC*NTYPEP)                                    CPP38040
      DIMENSION ESO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)         CPP38050
      DIMENSION PQ(*)                                                   CPP38060
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP38070
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP38080
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP38090
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP38100
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP38110
      COMMON/TOLER/CRIT,TOLR                                            CPP38120
      COMMON/CI102/NIND,NDEP                                            CPP38130
      DATA A00 / 0.0D+00 /                                              CPP38140
    1 FORMAT(//,2X,' B1 MATRIX'/)                                       CPP38150
    2 FORMAT(//,2X,' FINAL B1 MATRIX'/)                                 CPP38160
C                                                                       CPP38170
C   SEE EQ.(21)                                                         CPP38180
C   FORM DENSITY LIKE MATRIX IN SO BASIS                                CPP38190
C                                                                       CPP38200
      CALL PQEMAT(DP,DQ,B0,NIJ,ESO,NINDPX,NABC)                         CPP38210
C                                                                       CPP38220
C   FORM HALF-TRANSFORMED DENSITY LIKE MATRIX                           CPP38230
      NTYABC=NTYPES*NABC                                                CPP38240
      CALL PQMAT(DP,DQ,DM,NTYABC,MIJ,PQ)                                CPP38250
C                                                                       CPP38260
C   COMPLETE TRANSFORMATIONTION                                         CPP38270
      DO 101 ITYPX=1,NTYABC                                             CPP38280
      CALL MOCONS(DM(1,ITYPX),NTRI,DP(1,ITYPX),NTRI,ESO,U,T)            CPP38290
  101 CONTINUE                                                          CPP38300
C                                                                       CPP38310
      CALL ZERO(B1,NINDPX*NABC)                                         CPP38320
      DO 104 II=1,NIND                                                  CPP38330
      I=NIJ(II,1)                                                       CPP38340
      J=NIJ(II,2)                                                       CPP38350
      ITYP=MOTYP(I)                                                     CPP38360
      JTYP=MOTYP(J)                                                     CPP38370
      IJ=IOFF(I)+J                                                      CPP38380
      IADD=(ITYP-1)*NABC                                                CPP38390
      JADD=(JTYP-1)*NABC                                                CPP38400
      DO 103 IABC=1,NABC                                                CPP38410
      IAPT=IADD+IABC                                                    CPP38420
      JAPT=JADD+IABC                                                    CPP38430
  103 B1(II,IABC)=DP(IJ,IAPT)-DP(IJ,JAPT)                               CPP38440
  104 CONTINUE                                                          CPP38450
      IF(IPRNT.LE.3) GO TO 301                                          CPP38460
      WRITE(6,1)                                                        CPP38470
      CALL MATOUT(B1,NINDPX,NABC,NINDPX,NABC,6)                         CPP38480
C                                                                       CPP38490
C   SEE EQ.(12)                                                         CPP38500
C   REST OF EQ.(11)                                                     CPP38510
C                                                                       CPP38520
  301 CONTINUE                                                          CPP38530
      DO 110 II=1,NIND                                                  CPP38540
      I=NIJ(II,1)                                                       CPP38550
      J=NIJ(II,2)                                                       CPP38560
      ITYP=MOTYP(I)                                                     CPP38570
      JTYP=MOTYP(J)                                                     CPP38580
      DO 106 JJ=1,NIND                                                  CPP38590
      K=NIJ(JJ,1)                                                       CPP38600
      L=NIJ(JJ,2)                                                       CPP38610
      KTYP=MOTYP(K)                                                     CPP38620
      LTYP=MOTYP(L)                                                     CPP38630
      ZVAL1=A00                                                         CPP38640
      ZVAL2=A00                                                         CPP38650
      ZVAL3=A00                                                         CPP38660
      ZVAL4=A00                                                         CPP38670
      IF(J.NE.K) GO TO 201                                              CPP38680
      IL=IOFF(MAX0(I,L))+MIN0(I,L)                                      CPP38690
      ZVAL1=ZETA(IL,ITYP)-ZETA(IL,JTYP)                                 CPP38700
  201 IF(I.NE.K) GO TO 202                                              CPP38710
      JL=IOFF(MAX0(J,L))+MIN0(J,L)                                      CPP38720
      ZVAL2=-ZETA(JL,JTYP)+ZETA(JL,ITYP)                                CPP38730
  202 IF(J.NE.L) GO TO 203                                              CPP38740
      IK=IOFF(MAX0(I,K))+MIN0(I,K)                                      CPP38750
      ZVAL3=-ZETA(IK,ITYP)+ZETA(IK,JTYP)                                CPP38760
  203 IF(I.NE.L) GO TO 204                                              CPP38770
      JK=IOFF(MAX0(J,K))+MIN0(J,K)                                      CPP38780
      ZVAL4=ZETA(JK,JTYP)-ZETA(JK,ITYP)                                 CPP38790
  204 ZVALT=ZVAL1+ZVAL2+ZVAL3+ZVAL4                                     CPP38800
C                                                                       CPP38810
      IF(DABS(ZVALT).LT.TOLR) GO TO 106                                 CPP38820
      DO 105 IABC=1,NABC                                                CPP38830
      B1(II,IABC)=B1(II,IABC)+B0(JJ,IABC)*ZVALT                         CPP38840
  105 CONTINUE                                                          CPP38850
  106 CONTINUE                                                          CPP38860
C                                                                       CPP38870
      DO 108 JJ=1,NCPEX                                                 CPP38880
      DO 107 IABC=1,NABC                                                CPP38890
      B1(II,IABC)=B1(II,IABC)+B0(JJ+NIND,IABC)*A12(II,JJ)               CPP38900
C*    B1(II,IABC)=B1(II,IABC)-B0(JJ+NIND,IABC)*A12(II,JJ)               CPP38910
  107 CONTINUE                                                          CPP38920
  108 CONTINUE                                                          CPP38930
C                                                                       CPP38940
  110 CONTINUE                                                          CPP38950
C                                                                       CPP38960
      DO 120 II=1,NCPEX                                                 CPP38970
      DO 115 JJ=1,NIND                                                  CPP38980
      DO 114 IABC=1,NABC                                                CPP38990
      B1(II+NIND,IABC)=B1(II+NIND,IABC)+B0(JJ,IABC)*A12(JJ,II)          CPP39000
C*    B1(II+NIND,IABC)=B1(II+NIND,IABC)-B0(JJ,IABC)*A12(JJ,II)          CPP39010
  114 CONTINUE                                                          CPP39020
  115 CONTINUE                                                          CPP39030
C                                                                       CPP39040
      DO 117 JJ=1,NCPEX                                                 CPP39050
      DO 116 IABC=1,NABC                                                CPP39060
      B1(II+NIND,IABC)=B1(II+NIND,IABC)+B0(JJ+NIND,IABC)*A22(II,JJ)     CPP39070
  116 CONTINUE                                                          CPP39080
  117 CONTINUE                                                          CPP39090
C                                                                       CPP39100
  120 CONTINUE                                                          CPP39110
      IF(IPRNT.LE.4) GO TO 302                                          CPP39120
      WRITE(6,2)                                                        CPP39130
      CALL MATOUT(B1,NINDPX,NABC,NINDPX,NABC,6)                         CPP39140
C                                                                       CPP39150
  302 CONTINUE                                                          CPP39160
      RETURN                                                            CPP39170
      END                                                               CPP39180
      SUBROUTINE PQEMAT(DP,DQ,B0,NIJ,ESO,NINDPX,NABC)                   CPP39190
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP39200
      DIMENSION DP(NTRI,NABC*NTYPEP),DQ(NTRI,NABC*NTYPEP)               CPP39210
      DIMENSION B0(NINDPX,NABC),NIJ(NTRI,2),ESO(NBFAO,NBASIS)           CPP39220
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP39230
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP39240
      COMMON/GVBAB/A(30,30),B(30,30)                                    CPP39250
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP39260
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP39270
      COMMON/TOLER/CRIT,TOLR                                            CPP39280
      COMMON/CI102/NIND,NDEP                                            CPP39290
      DATA TWO / 2.0D+00 /                                              CPP39300
    1 FORMAT(//,2X,' DP MATRIX IN PQDMAT'/)                             CPP39310
    2 FORMAT(//,2X,' DQ MATRIX IN PQDMAT'/)                             CPP39320
C                                                                       CPP39330
C   SEE EQS.(22) AND (23)                                               CPP39340
C                                                                       CPP39350
      CALL ZERO(DP,NTRI*NTYPEP*NABC)                                    CPP39360
      CALL ZERO(DQ,NTRI*NTYPEP*NABC)                                    CPP39370
C                                                                       CPP39380
      IJ=0                                                              CPP39390
      DO 105 I=1,NBASIS                                                 CPP39400
      DO 105 J=1,I                                                      CPP39410
      IJ=IJ+1                                                           CPP39420
      DO 103 ITYP=1,NTYPES                                              CPP39430
      DO 102 II=1,NIND                                                  CPP39440
      K=NIJ(II,1)                                                       CPP39450
      L=NIJ(II,2)                                                       CPP39460
      KTYP=MOTYP(K)                                                     CPP39470
      LTYP=MOTYP(L)                                                     CPP39480
      AVAL=A(ITYP,KTYP)-A(ITYP,LTYP)                                    CPP39490
      BVAL=B(ITYP,KTYP)-B(ITYP,LTYP)                                    CPP39500
      FAC=ESO(I,K)*ESO(J,L)+ESO(I,L)*ESO(J,K)                           CPP39510
      IF(DABS(FAC).LT.TOLR) GO TO 102                                   CPP39520
      FACA=AVAL*FAC                                                     CPP39530
      FACB=BVAL*FAC                                                     CPP39540
      IADD=(ITYP-1)*NABC                                                CPP39550
      DO 101 IABC=1,NABC                                                CPP39560
      IAPT=IADD+IABC                                                    CPP39570
      DP(IJ,IAPT)=DP(IJ,IAPT)+B0(II,IABC)*FACA                          CPP39580
      DQ(IJ,IAPT)=DQ(IJ,IAPT)+B0(II,IABC)*FACB                          CPP39590
  101 CONTINUE                                                          CPP39600
  102 CONTINUE                                                          CPP39610
  103 CONTINUE                                                          CPP39620
  105 CONTINUE                                                          CPP39630
      NTYABC=NTYPES*NABC                                                CPP39640
      IF(IPRNT.LE.4) GO TO 201                                          CPP39650
      WRITE(6,1)                                                        CPP39660
      CALL MATOUT(DP,NTRI,NTYABC,NTRI,NTYABC,6)                         CPP39670
      WRITE(6,2)                                                        CPP39680
      CALL MATOUT(DQ,NTRI,NTYABC,NTRI,NTYABC,6)                         CPP39690
C                                                                       CPP39700
  201 CONTINUE                                                          CPP39710
      DO 107 I=1,NBASIS                                                 CPP39720
      II=IOFF(I+1)                                                      CPP39730
      DO 106 ITYPX=1,NTYABC                                             CPP39740
      DP(II,ITYPX)=DP(II,ITYPX)/TWO                                     CPP39750
      DQ(II,ITYPX)=DQ(II,ITYPX)/TWO                                     CPP39760
  106 CONTINUE                                                          CPP39770
  107 CONTINUE                                                          CPP39780
      CALL SCLMPY(DP,NTRI*NTYABC,TWO)                                   CPP39790
      CALL SCLMPY(DQ,NTRI*NTYABC,TWO)                                   CPP39800
C                                                                       CPP39810
      RETURN                                                            CPP39820
      END                                                               CPP39830
      SUBROUTINE UCMAT(EIG,ESO,NIJ,KIJ,U,T,B,DP,DQ,DM,DEL,MIJ,PQ,       CPP39840
     1                 NADD,N3NXX)                                      CPP39850
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP39860
      INTEGER*2 MIJ(*)                                                  CPP39870
      DIMENSION EIG(NBASIS)                                             CPP39880
      DIMENSION ESO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)         CPP39890
      DIMENSION NIJ(NTRI,2),KIJ(NTRI,2),B(NTRI),DEL(NTRI)               CPP39900
      DIMENSION DP(NTRI,NADD),DQ(NTRI,NADD),DM(NTRI,NADD)               CPP39910
      DIMENSION PQ(*)                                                   CPP39920
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP39930
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP39940
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP39950
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP39960
      COMMON/TAPES/ITAP37,ITAP42,ITAP43,ITAP44,ITAP47                   CPP39970
      COMMON/CI102/NIND,NDEP                                            CPP39980
      DATA DLIMIT / 1.0D-8 /                                            CPP39990
      DATA A00,ONE / 0.0D+00 , 1.0D+00 /                                CPP40000
    1 FORMAT(//,2X,' IGROUP = ',I5,5X,' IFIRST = ',I5,5X,               CPP40010
     1             ' NABC   = ',I5/)                                    CPP40020
    2 FORMAT(//,2X,' B0 MATRIX'/)                                       CPP40030
    3 FORMAT(//,2X,' DP MATRIX'/)                                       CPP40040
    4 FORMAT(//,2X,' DM MATRIX'/)                                       CPP40050
C                                                                       CPP40060
C   START GROUPING B0 MATRICES                                          CPP40070
      CALL RFILE(ITAP44)                                                CPP40080
      IGROUP=0                                                          CPP40090
      ILAST=0                                                           CPP40100
  500 CONTINUE                                                          CPP40110
      IGROUP=IGROUP+1                                                   CPP40120
      IFIRST=ILAST+1                                                    CPP40130
      ILAST=IFIRST+NADD-1                                               CPP40140
      IF(ILAST.GT.N3NXX) ILAST=N3NXX                                    CPP40150
      NABC=ILAST-IFIRST+1                                               CPP40160
      WRITE(3,1) IGROUP,IFIRST,NABC                                     CPP40170
      WRITE(4,1) IGROUP,IFIRST,NABC                                     CPP40180
C                                                                       CPP40190
C   CALL BACK B0 MATRICES                                               CPP40200
      DO 101 IABC=1,NABC                                                CPP40210
      III=IABC+IFIRST-1                                                 CPP40220
      IB=IBA(III)                                                       CPP40230
      CALL RREAD(ITAP44,DM(1,IABC),NTRI2,IB)                            CPP40240
  101 CONTINUE                                                          CPP40250
      IF(IPRNT.LE.4) GO TO 301                                          CPP40260
      WRITE(6,2)                                                        CPP40270
      CALL MATOUT(DM,NTRI,NABC,NTRI,NABC,6)                             CPP40280
C                                                                       CPP40290
C   FORM DENSITY LIKE MATRIX IN SO BASIS                                CPP40300
  301 CONTINUE                                                          CPP40310
      CALL PQFMAT(DP,DQ,DM,NIJ,ESO,NABC)                                CPP40320
      IF(IPRNT.LE.4) GO TO 302                                          CPP40330
      WRITE(6,2)                                                        CPP40340
      CALL MATOUT(DP,NTRI,NABC,NTRI,NABC,6)                             CPP40350
C                                                                       CPP40360
C   FORM HALF-TRANSFORMED DENSITY LIKE MATRIX                           CPP40370
C   AND TRANSFORM THEM INTO MO BASIS                                    CPP40380
  302 CONTINUE                                                          CPP40390
      CALL PQMAT(DP,DQ,DM,NABC,MIJ,PQ)                                  CPP40400
      DO 102 IABC=1,NABC                                                CPP40410
      CALL MOCONS(DM(1,IABC),NTRI,DP(1,IABC),NTRI,ESO,U,T)              CPP40420
  102 CONTINUE                                                          CPP40430
      IF(IPRNT.LE.4) GO TO 303                                          CPP40440
      WRITE(6,3)                                                        CPP40450
      CALL MATOUT(DP,NTRI,NABC,NTRI,NABC,6)                             CPP40460
C                                                                       CPP40470
  303 CONTINUE                                                          CPP40480
      DO 103 II=1,NDEP                                                  CPP40490
      I=KIJ(II,1)                                                       CPP40500
      J=KIJ(II,2)                                                       CPP40510
      DEL(II)=A00                                                       CPP40520
      DELJI=EIG(J)-EIG(I)                                               CPP40530
      IF(I.EQ.J) DELJI=ONE                                              CPP40540
      IF(DABS(DELJI).LE.DLIMIT) GO TO 103                               CPP40550
      DEL(II)=ONE/DELJI                                                 CPP40560
  103 CONTINUE                                                          CPP40570
C                                                                       CPP40580
      DO 105 IABC=1,NABC                                                CPP40590
      III=IABC+IFIRST-1                                                 CPP40600
      IB=IBA(III)                                                       CPP40610
      CALL RREAD(ITAP44,B,NTRI2,IB)                                     CPP40620
      DO 104 II=1,NDEP                                                  CPP40630
      I=KIJ(II,1)                                                       CPP40640
      J=KIJ(II,2)                                                       CPP40650
      IJ=IOFF(I)+J                                                      CPP40660
      TDEL=DEL(II)                                                      CPP40670
      B(IJ)=(B(IJ)+DP(IJ,IABC))*TDEL                                    CPP40680
  104 CONTINUE                                                          CPP40690
      CALL RWRIT(ITAP44,B,NTRI2,IB)                                     CPP40700
  105 CONTINUE                                                          CPP40710
      IF(ILAST.GE.N3NXX) GO TO 501                                      CPP40720
      GO TO 500                                                         CPP40730
C                                                                       CPP40740
  501 CONTINUE                                                          CPP40750
      CALL RCLOSE(ITAP44,3)                                             CPP40760
      RETURN                                                            CPP40770
      END                                                               CPP40780
      SUBROUTINE PQFMAT(DP,DQ,B,NIJ,ESO,NABC)                           CPP40790
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP40800
      DIMENSION DP(NTRI,NABC),DQ(NTRI,NABC),B(NTRI,NABC)                CPP40810
      DIMENSION ESO(NBFAO,NBASIS),NIJ(NTRI,2)                           CPP40820
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP40830
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP40840
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP40850
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP40860
      COMMON/TOLER/CRIT,TOLR                                            CPP40870
      COMMON/CI102/NIND,NDEP                                            CPP40880
      DATA A00,HALF,TWO / 0.0D+00 , 0.5D+00 , 2.0D+00 /                 CPP40890
    1 FORMAT(//,2X,' DP MATRIX'/)                                       CPP40900
    2 FORMAT(//,2X,' DQ MATRIX'/)                                       CPP40910
C                                                                       CPP40920
C   FORM DENSITY MATRICES                                               CPP40930
C                                                                       CPP40940
C   FORM DP MATRIX                                                      CPP40950
      CALL ZERO(DP,NTRI*NABC)                                           CPP40960
      CALL ZERO(DQ,NTRI*NABC)                                           CPP40970
C                                                                       CPP40980
      IJ=0                                                              CPP40990
      DO 103 I=1,NBASIS                                                 CPP41000
      DO 103 J=1,I                                                      CPP41010
      IJ=IJ+1                                                           CPP41020
      DO 102 II=1,NIND                                                  CPP41030
      K=NIJ(II,1)                                                       CPP41040
      L=NIJ(II,2)                                                       CPP41050
      KTYP=MOTYP(K)                                                     CPP41060
      LTYP=MOTYP(L)                                                     CPP41070
      KL=IOFF(K)+L                                                      CPP41080
      FAC=(ESO(I,K)*ESO(J,L)+ESO(I,L)*ESO(J,K))                         CPP41090
     1   *(FOCC(LTYP)-FOCC(KTYP))*HALF                                  CPP41100
      IF(DABS(FAC).LT.TOLR) GO TO 102                                   CPP41110
      DO 101 IABC=1,NABC                                                CPP41120
      FACT=B(KL,IABC)*FAC                                               CPP41130
      DP(IJ,IABC)=DP(IJ,IABC)+FACT*TWO                                  CPP41140
      DQ(IJ,IABC)=DQ(IJ,IABC)-FACT                                      CPP41150
  101 CONTINUE                                                          CPP41160
  102 CONTINUE                                                          CPP41170
  103 CONTINUE                                                          CPP41180
C                                                                       CPP41190
      IF(IPRNT.LE.4) GO TO 201                                          CPP41200
      WRITE(6,1)                                                        CPP41210
      CALL MATOUT(DP,NTRI,NABC,NTRI,NABC,6)                             CPP41220
      WRITE(6,2)                                                        CPP41230
      CALL MATOUT(DQ,NTRI,NABC,NTRI,NABC,6)                             CPP41240
C                                                                       CPP41250
  201 CONTINUE                                                          CPP41260
      DO 105 I=1,NBASIS                                                 CPP41270
      II=IOFF(I+1)                                                      CPP41280
      DO 104 ITYPX=1,NABC                                               CPP41290
      DP(II,ITYPX)=DP(II,ITYPX)/TWO                                     CPP41300
      DQ(II,ITYPX)=DQ(II,ITYPX)/TWO                                     CPP41310
  104 CONTINUE                                                          CPP41320
  105 CONTINUE                                                          CPP41330
      CALL SCLMPY(DP,NTRI*NABC,TWO)                                     CPP41340
      CALL SCLMPY(DQ,NTRI*NABC,TWO)                                     CPP41350
C                                                                       CPP41360
      RETURN                                                            CPP41370
      END                                                               CPP41380
      SUBROUTINE ORTHOG(U,W,NDIM,MDIM,N,NDUM,NTEMP,NVEC)                CPP41390
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP41400
      DIMENSION U(NDIM,MDIM),W(NDIM),NDUM(1)                            CPP41410
      DATA A00 / 0.0D+00 /                                              CPP41420
C                                                                       CPP41430
      DO 101 I=1,NTEMP                                                  CPP41440
  101 NDUM(I)=0                                                         CPP41450
C   N=1                                                                 CPP41460
      DO 102 I=1,N                                                      CPP41470
  102 W(I)=U(I,1)                                                       CPP41480
      CALL NORMLZ(W,NDIM,N,IFLAG)                                       CPP41490
      DO 103 I=1,N                                                      CPP41500
  103 U(I,1)=W(I)                                                       CPP41510
C                                                                       CPP41520
C   N=>2                                                                CPP41530
C                                                                       CPP41540
      NVEC=1                                                            CPP41550
      IF(NTEMP.EQ.1) RETURN                                             CPP41560
      DO 110 II=2,NTEMP                                                 CPP41570
      DO 104 I=1,N                                                      CPP41580
  104 W(I)=U(I,II)                                                      CPP41590
      DO 107 K=1,II-1                                                   CPP41600
      IF(NDUM(K).NE.0) GO TO 107                                        CPP41610
      EU=A00                                                            CPP41620
      DO 105 I=1,N                                                      CPP41630
  105 EU=EU+U(I,K)*U(I,II)                                              CPP41640
      DO 106 I=1,N                                                      CPP41650
  106 W(I)=W(I)-EU*U(I,K)                                               CPP41660
  107 CONTINUE                                                          CPP41670
      CALL NORMLZ(W,NDIM,N,IFLAG)                                       CPP41680
      IF(IFLAG.NE.0) NDUM(II)=1                                         CPP41690
      IF(IFLAG.NE.0) GO TO 110                                          CPP41700
      NVEC=NVEC+1                                                       CPP41710
      DO 108 I=1,N                                                      CPP41720
  108 U(I,II)=W(I)                                                      CPP41730
  110 CONTINUE                                                          CPP41740
C                                                                       CPP41750
      IX=0                                                              CPP41760
      DO 120 II=1,NTEMP                                                 CPP41770
      IF(NDUM(II).NE.0) GO TO 120                                       CPP41780
      IX=IX+1                                                           CPP41790
      DO 115 I=1,N                                                      CPP41800
  115 W(I)=U(I,II)                                                      CPP41810
      DO 116 I=1,N                                                      CPP41820
  116 U(I,IX)=W(I)                                                      CPP41830
  120 CONTINUE                                                          CPP41840
C                                                                       CPP41850
      RETURN                                                            CPP41860
      END                                                               CPP41870
      SUBROUTINE NORMLZ(V,NDIM,N,IFLAG)                                 CPP41880
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP41890
      DIMENSION V(NDIM)                                                 CPP41900
    1 FORMAT(//,2X,'***ERROR IN NORMALIZATION***'/)                     CPP41910
      DATA RLIM / 1.0D-38 /                                             CPP41920
      DATA VLIM / 1.0D-06 /                                             CPP41930
      DATA A00 / 0.0D+00 /                                              CPP41940
C                                                                       CPP41950
      IFLAG=0                                                           CPP41960
      SQL=A00                                                           CPP41970
      DO 101 I=1,N                                                      CPP41980
  101 SQL=SQL+V(I)*V(I)                                                 CPP41990
      R=DSQRT(SQL)                                                      CPP42000
      IF(R.LE.RLIM) GO TO 202                                           CPP42010
      IF(R.LE.VLIM) GO TO 204                                           CPP42020
C                                                                       CPP42030
  201 DO 102 I=1,N                                                      CPP42040
  102 V(I)=V(I)/R                                                       CPP42050
      GO TO 205                                                         CPP42060
  202 WRITE(6,1)                                                        CPP42070
      STOP                                                              CPP42080
C                                                                       CPP42090
  204 IFLAG=1                                                           CPP42100
  205 RETURN                                                            CPP42110
      END                                                               CPP42120
      SUBROUTINE NORM(XAB,A,B,N)                                        CPP42130
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP42140
      DIMENSION A(1),B(1)                                               CPP42150
      DATA A00 / 0.0D+00 /                                              CPP42160
C                                                                       CPP42170
      XAB=A00                                                           CPP42180
      DO 101 I=1,N                                                      CPP42190
  101 XAB=XAB+A(I)*B(I)                                                 CPP42200
      RETURN                                                            CPP42210
      END                                                               CPP42220
      SUBROUTINE UXMAT(KIJ,B,SM,U,NABC)                                 CPP42230
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP42240
      DIMENSION KIJ(NTRI,2),B(NTRI,NABC),SM(NTRI,NABC),U(NBASIS,NBASIS) CPP42250
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP42260
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP42270
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP42280
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP42290
      COMMON/TAPES/ITAP37,ITAP42,ITAP43,ITAP44,ITAP47                   CPP42300
      COMMON/CI102/NIND,NDEP                                            CPP42310
      DATA A00,HALF / 0.0D+00 , 0.5D+00 /                               CPP42320
    1 FORMAT(//,2X,' U MATRIX, IABC= ',I5/)                             CPP42330
C                                                                       CPP42340
      CALL RFILE(ITAP44)                                                CPP42350
      DO 101 IABC=1,NABC                                                CPP42360
      IB=IBA(IABC)                                                      CPP42370
      CALL RREAD(ITAP44,B(1,IABC),NTRI2,IB)                             CPP42380
      IF(IABC.GT.N3N) GO TO 101                                         CPP42390
      IS=ISA(IABC)                                                      CPP42400
      CALL RREAD(ITAP44,SM(1,IABC),NTRI2,IS)                            CPP42410
  101 CONTINUE                                                          CPP42420
C                                                                       CPP42430
C   FILL UP THE REST OF U MATRIX                                        CPP42440
      DO 105 IABC=1,NABC                                                CPP42450
      IU=IUA(IABC)                                                      CPP42460
      IF(IABC.GT.N3N) GO TO 301                                         CPP42470
C                                                                       CPP42480
      IJ=0                                                              CPP42490
      DO 102 I=1,NBASIS                                                 CPP42500
      DO 102 J=1,I                                                      CPP42510
      IJ=IJ+1                                                           CPP42520
      IF(I.EQ.J) GO TO 201                                              CPP42530
      U(I,J)=B(IJ,IABC)                                                 CPP42540
      U(J,I)=-U(I,J)-SM(IJ,IABC)                                        CPP42550
      GO TO 102                                                         CPP42560
  201 U(I,I)=-SM(IJ,IABC)*HALF                                          CPP42570
  102 CONTINUE                                                          CPP42580
C   SPECIFIY NON-INDEPENDENT VARIABLES                                  CPP42590
      DO 103 II=1,NDEP                                                  CPP42600
      I=KIJ(II,1)                                                       CPP42610
      J=KIJ(II,2)                                                       CPP42620
      IJ=IOFF(I)+J                                                      CPP42630
      U(I,J)=-SM(IJ,IABC)*HALF                                          CPP42640
      U(J,I)=U(I,J)                                                     CPP42650
  103 CONTINUE                                                          CPP42660
      GO TO 302                                                         CPP42670
C                                                                       CPP42680
  301 CONTINUE                                                          CPP42690
      IJ=0                                                              CPP42700
      DO 104 I=1,NBASIS                                                 CPP42710
      DO 104 J=1,I                                                      CPP42720
      IJ=IJ+1                                                           CPP42730
      IF(I.EQ.J) GO TO 202                                              CPP42740
      U(I,J)=B(IJ,IABC)                                                 CPP42750
      U(J,I)=-U(I,J)                                                    CPP42760
      GO TO 104                                                         CPP42770
  202 U(I,I)=A00                                                        CPP42780
  104 CONTINUE                                                          CPP42790
C                                                                       CPP42800
  302 CONTINUE                                                          CPP42810
      CALL RWRIT(ITAP44,U,NBASQ,IU)                                     CPP42820
      IF(IPRNT.LE.3) GO TO 105                                          CPP42830
      WRITE(6,1) IABC                                                   CPP42840
      CALL MATOUT(U,NBASIS,NBASIS,NBASIS,NBASIS,6)                      CPP42850
C                                                                       CPP42860
  105 CONTINUE                                                          CPP42870
C                                                                       CPP42880
      CALL RCLOSE(ITAP44,3)                                             CPP42890
      RETURN                                                            CPP42900
      END                                                               CPP42910
      SUBROUTINE WAMAT(ESO,ZETA,U,T,EPA,WA,SA,DP,DQ,DM,MIJ,PQ,NADD)     CPP42920
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP42930
      INTEGER*2 MIJ(*)                                                  CPP42940
      DIMENSION ESO(NBFAO,NBASIS),U(NBASIS,NBASIS),T(NBASIS,NBASIS)     CPP42950
      DIMENSION EPA(NBASIS,NBASIS)                                      CPP42960
      DIMENSION ZETA(NTRI,NTYPEP),WA(NBASIS,NBASIS),SA(NTRI,NADD*NTYPEP)CPP42970
      DIMENSION DP(NTRI,NADD*NTYPEP),DQ(NTRI,NADD*NTYPEP)               CPP42980
      DIMENSION DM(NTRI,NADD*NTYPEP)                                    CPP42990
      DIMENSION PQ(*)                                                   CPP43000
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP43010
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP43020
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP43030
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP43040
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP43050
      COMMON/TAPES/ITAP37,ITAP42,ITAP43,ITAP44,ITAP47                   CPP43060
      COMMON/WORKS/JTAPE1,JTAPE2                                        CPP43070
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP43080
      DATA A00 / 0.0D+00 /                                              CPP43090
    1 FORMAT(//,2X,' IGROUP = ',I5,5X,' IFIRST = ',I5,5X,               CPP43100
     1             ' NABC   = ',I5/)                                    CPP43110
    2 FORMAT(//,2X,' WA MATRIX, IABC = ',I5/)                           CPP43120
C                                                                       CPP43130
C   FORM WA MATRICES                                                    CPP43140
C   SEE EQ.(13)                                                         CPP43150
C                                                                       CPP43160
      CALL RFILE(ITAP44)                                                CPP43170
      CALL RFILE(JTAPE1)                                                CPP43180
C                                                                       CPP43190
C   CALL BACK LAGRANGIAN MATRICES                                       CPP43200
      DO 101 ITYP=1,NTYPES                                              CPP43210
      CALL RMREAD(ZETA(1,ITYP),35+ITYP)                                 CPP43220
  101 CONTINUE                                                          CPP43230
      CALL ZERO(ZETA(1,NTYPEP),NTRI)                                    CPP43240
C                                                                       CPP43250
C   START GROUPING WA MATRICES                                          CPP43260
      IGROUP=0                                                          CPP43270
      ILAST=0                                                           CPP43280
  500 CONTINUE                                                          CPP43290
      IGROUP=IGROUP+1                                                   CPP43300
      IFIRST=ILAST+1                                                    CPP43310
      ILAST=IFIRST+NADD-1                                               CPP43320
      IF(ILAST.GT.N3N) ILAST=N3N                                        CPP43330
      NABC=ILAST-IFIRST+1                                               CPP43340
      WRITE(3,1) IGROUP,IFIRST,NABC                                     CPP43350
      WRITE(4,1) IGROUP,IFIRST,NABC                                     CPP43360
C                                                                       CPP43370
C   CALL BACK DERIVATIVE OVERLAP INTEGRALS                              CPP43380
      DO 103 IABC=1,NABC                                                CPP43390
      III=IABC+IFIRST-1                                                 CPP43400
      IS=ISA(III)                                                       CPP43410
      CALL RREAD(ITAP44,SA(1,IABC),NTRI2,IS)                            CPP43420
  103 CONTINUE                                                          CPP43430
C                                                                       CPP43440
C   FORM DENSITY LIKE MATRIX                                            CPP43450
      CALL PQGMAT(DP,DQ,SA,ESO,NABC)                                    CPP43460
C                                                                       CPP43470
C   FORM HALF-TRANSFORMED DENSITY LIKE MATRIX                           CPP43480
      NTYABC=NTYPES*NABC                                                CPP43490
      CALL PQMAT(DP,DQ,DM,NTYABC,MIJ,PQ)                                CPP43500
C                                                                       CPP43510
C   COMPLETE TRANSFORMATION                                             CPP43520
      DO 105 ITYPX=1,NTYABC                                             CPP43530
      CALL MOCONS(DM(1,ITYPX),NTRI,DP(1,ITYPX),NTRI,ESO,U,T)            CPP43540
  105 CONTINUE                                                          CPP43550
C                                                                       CPP43560
      DO 110 IABC=1,NABC                                                CPP43570
      III=IABC+IFIRST-1                                                 CPP43580
      IS=ISA(III)                                                       CPP43590
      IE=IEA(III)                                                       CPP43600
      CALL RREAD(ITAP44,SA(1,IABC),NTRI2,IS)                            CPP43610
      CALL RREAD(ITAP44,EPA,NBASQ,IE)                                   CPP43620
      DO 107 I=1,NBASIS                                                 CPP43630
      ITYP=MOTYP(I)                                                     CPP43640
      DO 107 J=1,NBASIS                                                 CPP43650
      JTYP=MOTYP(J)                                                     CPP43660
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPP43670
      JAPT=(JTYP-1)*NABC+IABC                                           CPP43680
      VALU1=EPA(J,I)+EPA(J,I)                                           CPP43690
      VALU2=DP(IJ,JAPT)                                                 CPP43700
      VALU3=A00                                                         CPP43710
      DO 106 K=1,JOCC                                                   CPP43720
      IK=IOFF(MAX0(I,K))+MIN0(I,K)                                      CPP43730
      JK=IOFF(MAX0(J,K))+MIN0(J,K)                                      CPP43740
      FAC1=ZETA(IK,ITYP)+ZETA(IK,JTYP)                                  CPP43750
      FAC2=ZETA(JK,JTYP)+ZETA(JK,JTYP)                                  CPP43760
      VALU3=VALU3+SA(JK,IABC)*FAC1+SA(IK,IABC)*FAC2                     CPP43770
  106 CONTINUE                                                          CPP43780
      WA(I,J)=VALU1-VALU2-VALU3                                         CPP43790
  107 CONTINUE                                                          CPP43800
      CALL SWRIT(JTAPE1,WA,NBASQ)                                       CPP43810
      IF(IPRNT.LE.3) GO TO 110                                          CPP43820
      WRITE(6,2) IABC                                                   CPP43830
      CALL MATOUT(WA,NBASIS,NBASIS,NBASIS,NBASIS,6)                     CPP43840
  110 CONTINUE                                                          CPP43850
      IF(ILAST.GE.N3N) GO TO 501                                        CPP43860
      GO TO 500                                                         CPP43870
C                                                                       CPP43880
  501 CONTINUE                                                          CPP43890
      CALL RCLOSE(ITAP44,3)                                             CPP43900
      CALL RCLOSE(JTAPE1,3)                                             CPP43910
      RETURN                                                            CPP43920
      END                                                               CPP43930
      SUBROUTINE PQGMAT(DP,DQ,SA,ESO,NABC)                              CPP43940
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP43950
      DIMENSION DP(NTRI,NABC*NTYPEP),DQ(NTRI,NABC*NTYPEP)               CPP43960
      DIMENSION ESO(NBFAO,NBASIS),SA(NTRI,NABC*NTYPEP)                  CPP43970
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP43980
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP43990
      COMMON/GVBAB/A(30,30),B(30,30)                                    CPP44000
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP44010
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP44020
      COMMON/TOLER/CRIT,TOLR                                            CPP44030
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP44040
      DATA A00,TWO / 0.0D+00 , 2.0D+00 /                                CPP44050
    1 FORMAT(//,2X,' DP MATRIX IN PQGMAT'/)                             CPP44060
    2 FORMAT(//,2X,' DQ MATRIX IN PQGMAT'/)                             CPP44070
C                                                                       CPP44080
      CALL ZERO(DP,NTRI*NTYPEP*NABC)                                    CPP44090
      CALL ZERO(DQ,NTRI*NTYPEP*NABC)                                    CPP44100
C                                                                       CPP44110
      IJ=0                                                              CPP44120
      DO 105 I=1,NBASIS                                                 CPP44130
      DO 105 J=1,I                                                      CPP44140
      IJ=IJ+1                                                           CPP44150
      DO 104 K=1,JOCC                                                   CPP44160
      KTYP=MOTYP(K)                                                     CPP44170
      DO 104 L=1,JOCC                                                   CPP44180
      LTYP=MOTYP(L)                                                     CPP44190
      KL=IOFF(MAX0(K,L))+MIN0(K,L)                                      CPP44200
      FAC=ESO(I,K)*ESO(J,L)+ESO(J,K)*ESO(I,L)                           CPP44210
      IF(DABS(FAC).LT.TOLR) GO TO 104                                   CPP44220
      DO 103 ITYP=1,NTYPES                                              CPP44230
      AVAL=A(ITYP,LTYP)                                                 CPP44240
      BVAL=B(ITYP,LTYP)                                                 CPP44250
      FACA=AVAL*FAC                                                     CPP44260
      FACB=BVAL*FAC                                                     CPP44270
      IADD=(ITYP-1)*NABC                                                CPP44280
      DO 102 IABC=1,NABC                                                CPP44290
      IAPT=IADD+IABC                                                    CPP44300
      DP(IJ,IAPT)=DP(IJ,IAPT)+SA(KL,IABC)*FACA                          CPP44310
      DQ(IJ,IAPT)=DQ(IJ,IAPT)+SA(KL,IABC)*FACB                          CPP44320
  102 CONTINUE                                                          CPP44330
  103 CONTINUE                                                          CPP44340
  104 CONTINUE                                                          CPP44350
  105 CONTINUE                                                          CPP44360
      NTYABC=NTYPES*NABC                                                CPP44370
      IF(IPRNT.LE.3) GO TO 201                                          CPP44380
      WRITE(6,1)                                                        CPP44390
      CALL MATOUT(DP,NTRI,NTYABC,NTRI,NTYABC,6)                         CPP44400
      WRITE(6,2)                                                        CPP44410
      CALL MATOUT(DQ,NTRI,NTYABC,NTRI,NTYABC,6)                         CPP44420
C                                                                       CPP44430
  201 CONTINUE                                                          CPP44440
      DO 107 I=1,NBASIS                                                 CPP44450
      II=IOFF(I+1)                                                      CPP44460
      DO 106 ITYPX=1,NTYABC                                             CPP44470
      DP(II,ITYPX)=DP(II,ITYPX)/TWO                                     CPP44480
      DQ(II,ITYPX)=DQ(II,ITYPX)/TWO                                     CPP44490
  106 CONTINUE                                                          CPP44500
  107 CONTINUE                                                          CPP44510
      CALL SCLMPY(DP,NTRI*NTYABC,TWO)                                   CPP44520
      CALL SCLMPY(DQ,NTRI*NTYABC,TWO)                                   CPP44530
C                                                                       CPP44540
      RETURN                                                            CPP44550
      END                                                               CPP44560
      SUBROUTINE SCF2W(E2M,E2P,E2Q,EPA,UA,SA,WA)                        CPP44570
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP44580
      DIMENSION E2M(N3N,N3N),E2P(N3N,N3N),E2Q(N3N,N3N)                  CPP44590
      DIMENSION EPA(NBASIS,NBASIS),UA(NBASIS,NBASIS)                    CPP44600
      DIMENSION SA(NTRI,N3N),WA(NBASIS,NBASIS,N3N)                      CPP44610
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP44620
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP44630
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP44640
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP44650
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP44660
      COMMON/TAPES/ITAP37,ITAP42,ITAP43,ITAP44,ITAP47                   CPP44670
      COMMON/WORKS/JTAPE1,JTAPE2                                        CPP44680
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP44690
      DATA A00,TWO,FOUR / 0.0D+00 , 2.0D+00 , 4.0D+00 /                 CPP44700
    1 FORMAT(//,2X,' SA MATRIX, IABC = ',I5/)                           CPP44710
    2 FORMAT(//,2X,' WA MATRIX, IABC = ',I5/)                           CPP44720
    3 FORMAT(//,2X,' E2P MATRIX IN SCF2W'/)                             CPP44730
    4 FORMAT(//,2X,' E2Q MATRIX IN SCF2W'/)                             CPP44740
    5 FORMAT(//,2X,' E2M MATRIX IN SCF2W'/)                             CPP44750
C                                                                       CPP44760
C   CALL BACK SA , UA  AND WA MATRICES                                  CPP44770
      CALL RFILE(ITAP44)                                                CPP44780
      CALL RFILE(JTAPE1)                                                CPP44790
      DO 101 IABC=1,N3N                                                 CPP44800
      IS=ISA(IABC)                                                      CPP44810
      CALL RREAD(ITAP44,SA(1,IABC),NTRI2,IS)                            CPP44820
      CALL SREAD(JTAPE1,WA(1,1,IABC),NBASQ)                             CPP44830
      IF(IPRNT.LE.3) GO TO 101                                          CPP44840
      WRITE(6,1) IABC                                                   CPP44850
      CALL PRINT(SA(1,IABC),NTRI,NBASIS,6)                              CPP44860
      WRITE(6,2) IABC                                                   CPP44870
      CALL MATOUT(WA(1,1,IABC),NBASIS,NBASIS,NBASIS,NBASIS,6)           CPP44880
  101 CONTINUE                                                          CPP44890
C                                                                       CPP44900
C   CALCULATE SCF SECOND DERIVATIVES                                    CPP44910
C   SEE EQ.(12)                                                         CPP44920
C                                                                       CPP44930
C   (B COORDINATE)                                                      CPP44940
      DO 106 IABC=1,N3N                                                 CPP44950
      IE=IEA(IABC)                                                      CPP44960
      IU=IUA(IABC)                                                      CPP44970
      CALL RREAD(ITAP44,EPA,NBASQ,IE)                                   CPP44980
      CALL RREAD(ITAP44,UA,NBASQ,IU)                                    CPP44990
C                                                                       CPP45000
C   (A COORDINATE)                                                      CPP45010
      DO 105 JABC=1,N3N                                                 CPP45020
      VALUP=A00                                                         CPP45030
      VALUQ=A00                                                         CPP45040
      DO 104 I=1,JOCC                                                   CPP45050
      DO 102 J=1,JOCC                                                   CPP45060
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPP45070
      VALUP=VALUP-SA(IJ,JABC)*EPA(I,J)                                  CPP45080
  102 CONTINUE                                                          CPP45090
      DO 103 J=1,NBASIS                                                 CPP45100
      VALUQ=VALUQ+WA(J,I,JABC)*UA(J,I)                                  CPP45110
  103 CONTINUE                                                          CPP45120
  104 CONTINUE                                                          CPP45130
      E2P(IABC,JABC)=VALUP*TWO                                          CPP45140
      E2Q(IABC,JABC)=VALUQ*TWO                                          CPP45150
  105 CONTINUE                                                          CPP45160
  106 CONTINUE                                                          CPP45170
      IF(IPRNT.LE.3) GO TO 201                                          CPP45180
      WRITE(6,3)                                                        CPP45190
      CALL MATOUT(E2P,N3N,N3N,N3N,N3N,6)                                CPP45200
      WRITE(6,4)                                                        CPP45210
      CALL MATOUT(E2Q,N3N,N3N,N3N,N3N,6)                                CPP45220
C                                                                       CPP45230
C*******************************************                            CPP45240
C***SUM UP CONTRIBUTIONS DUE TO MO CHANGE***                            CPP45250
C*******************************************                            CPP45260
  201 CONTINUE                                                          CPP45270
      DO 110 I=1,N3N                                                    CPP45280
      DO 110 J=1,N3N                                                    CPP45290
      E2M(I,J)=E2P(I,J)+E2Q(I,J)                                        CPP45300
  110 CONTINUE                                                          CPP45310
C                                                                       CPP45320
      IF(IPRNT.LE.2) GO TO 202                                          CPP45330
      WRITE(6,5)                                                        CPP45340
      CALL MATOUT(E2M,N3N,N3N,N3N,N3N,6)                                CPP45350
C                                                                       CPP45360
  202 CONTINUE                                                          CPP45370
      CALL RCLOSE(ITAP44,3)                                             CPP45380
      CALL RCLOSE(JTAPE1,3)                                             CPP45390
      RETURN                                                            CPP45400
      END                                                               CPP45410
      SUBROUTINE SCF2ND(E2A,HCIA,CIA,E2M,E2C,E2T,SA,EIJ)                CPP45420
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP45430
      DIMENSION HCIA(NCTRI,N3N),CIA(NCPEX,N3N)                          CPP45440
      DIMENSION E2A(N3N,N3N),E2M(N3N,N3N),E2C(N3N,N3N),E2T(N3N,N3N)     CPP45450
      DIMENSION SA(NTRI,N3N),EIJ(NBASIS,NBASIS)                         CPP45460
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP45470
      COMMON/CIOFF/IIUNQ,IJUNQ,ICIUNQ,IPX(5050),JPX(5050)               CPP45480
      COMMON/CIPOS/ICITYP(5050),ICPOS(200,100)                          CPP45490
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP45500
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP45510
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP45520
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP45530
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP45540
      COMMON/TAPES/ITAP35,ITAP42,ITAP43,ITAP44,ITAP47                   CPP45550
      COMMON/WORKS/JTAPE1,JTAPE2                                        CPP45560
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP45570
      COMMON/CI103/CI(100)                                              CPP45580
      COMMON/CI104/IEIJ(5050)                                           CPP45590
      DATA A00,HALF,TWO / 0.0D+00 , 0.5D+00 , 2.0D+00 /                 CPP45600
    1 FORMAT(///                                                        CPP45610
     1 2X,'::::::::::::::::::::::::::::'/                               CPP45620
     2 2X,':::SUMMARY OF CALCULATION:::'/                               CPP45630
     3 2X,'::::::::::::::::::::::::::::'/)                              CPP45640
    2 FORMAT(//,2X,' E2A MATRIX'/)                                      CPP45650
    3 FORMAT(//,2X,' E2M MATRIX'/)                                      CPP45660
    4 FORMAT(//,2X,' E2C MATRIX'/)                                      CPP45670
    5 FORMAT(//,2X,' SCF SECOND DERIVATIVE MATRIX'/)                    CPP45680
    6 FORMAT(2I5)                                                       CPP45690
    7 FORMAT(3F20.10)                                                   CPP45700
C                                                                       CPP45710
C   CALCULATE SCF SECOND DERIVATIVES                                    CPP45720
C   SEE EQ.(12)                                                         CPP45730
C                                                                       CPP45740
C***********************************************                        CPP45750
C***CONTRIBUTION DUE TO CI COEFFICIENT CHANGE***                        CPP45760
C***********************************************                        CPP45770
C                                                                       CPP45780
      CALL RFILE(ITAP44,3)                                              CPP45790
      DO 101 IABC=1,N3N                                                 CPP45800
      IS=ISA(IABC)                                                      CPP45810
      CALL RREAD(ITAP44,SA(1,IABC),NTRI2,IS)                            CPP45820
  101 CONTINUE                                                          CPP45830
C                                                                       CPP45840
      CALL RFILE(JTAPE2)                                                CPP45850
      CALL ZERO(E2C,N3N*N3N)                                            CPP45860
      DO 105 ICI=1,NCPEX                                                CPP45870
      DO 105 JCI=1,NCPEX                                                CPP45880
      IJCI=IOFF(MAX0(ICI,JCI))+MIN0(ICI,JCI)                            CPP45890
      DO 102 III=1,ICIUNQ                                               CPP45900
      IEX=IPX(III)                                                      CPP45910
      JEX=JPX(III)                                                      CPP45920
      IJEX=IOFF(MAX0(IEX,JEX))+MIN0(IEX,JEX)                            CPP45930
      IF(IJCI.EQ.IJEX) GO TO 301                                        CPP45940
  102 CONTINUE                                                          CPP45950
      GO TO 105                                                         CPP45960
C                                                                       CPP45970
C   CUMMULATE ONLY NON-ZERO CI'S                                        CPP45980
  301 CONTINUE                                                          CPP45990
      IBCI=ICITYP(III)                                                  CPP46000
      IPT=IEIJ(IBCI)                                                    CPP46010
      CALL RREAD(JTAPE2,EIJ,NBASQ,IPT)                                  CPP46020
C                                                                       CPP46030
      DO 104 IABC=1,N3N                                                 CPP46040
      DO 104 JABC=1,N3N                                                 CPP46050
      VALU0=A00                                                         CPP46060
      DO 103 I=1,JOCC                                                   CPP46070
      DO 103 J=1,JOCC                                                   CPP46080
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPP46090
C*    VALU0=VALU0-SA(IJ,JABC)*EIJ(I,J)                                  CPP46100
      VALU0=VALU0-SA(IJ,JABC)*EIJ(J,I)                                  CPP46110
  103 CONTINUE                                                          CPP46120
      VALUC=CIA(ICI,IABC)*(HCIA(IJCI,JABC)+VALU0*TWO)*CI(JCI)           CPP46130
      E2C(IABC,JABC)=E2C(IABC,JABC)+VALUC*TWO                           CPP46140
  104 CONTINUE                                                          CPP46150
  105 CONTINUE                                                          CPP46160
C                                                                       CPP46170
      WRITE(6,1)                                                        CPP46180
      WRITE(6,2)                                                        CPP46190
      CALL MATOUT(E2A,N3N,N3N,N3N,N3N,6)                                CPP46200
      WRITE(6,3)                                                        CPP46210
      CALL MATOUT(E2M,N3N,N3N,N3N,N3N,6)                                CPP46220
      WRITE(6,4)                                                        CPP46230
      CALL MATOUT(E2C,N3N,N3N,N3N,N3N,6)                                CPP46240
      DO 106 I=1,N3N                                                    CPP46250
      DO 106 J=1,N3N                                                    CPP46260
  106 E2T(I,J)=E2A(I,J)+E2M(I,J)+E2C(I,J)                               CPP46270
      WRITE(6,5)                                                        CPP46280
      CALL MATOUT(E2T,N3N,N3N,N3N,N3N,6)                                CPP46290
C                                                                       CPP46300
      ITAP15=15                                                         CPP46310
      REWIND ITAP15                                                     CPP46320
      N6N=N3N*2                                                         CPP46330
      WRITE(ITAP15,6) NATOMS,N6N                                        CPP46340
      WRITE(ITAP15,7) ((E2T(I,J),J=1,N3N),I=1,N3N)                      CPP46350
      REWIND ITAP15                                                     CPP46360
C                                                                       CPP46370
      CALL RCLOSE(ITAP44,3)                                             CPP46380
      CALL RCLOSE(JTAPE2,3)                                             CPP46390
      RETURN                                                            CPP46400
      END                                                               CPP46410
      SUBROUTINE DIPMO(HCIA,CIA,HF,UA)                                  CPP46420
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP46430
      DIMENSION HF(NTRI,3),UA(NBASIS,NBASIS,N3N)                        CPP46440
      DIMENSION HCIA(NCTRI,N3NP3),CIA(NCPEX,N3NP3)                      CPP46450
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP46460
      COMMON/CALIF/LPARA(1024),DIPDE(3,45),DIPDC(3,45),DUM(754)         CPP46470
      COMMON/CIOFF/IIUNQ,IJUNQ,ICIUNQ,IPX(5050),JPX(5050)               CPP46480
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP46490
      COMMON/DIPOL/DIPDA(3,45),DIPDM(3,45),DIPDN(3,45),DIPDT(3,45)      CPP46500
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP46510
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP46520
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP46530
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP46540
      COMMON/TAPES/ITAP37,ITAP42,ITAP43,ITAP44,ITAP47                   CPP46550
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP46560
      COMMON/CI103/CI(100)                                              CPP46570
      DATA A00,HALF,ONE,TWO,FOUR / 0.0D+00 , 0.5D+00 , 1.0D+00 ,        CPP46580
     1                             2.0D+00 , 4.0D+00 /                  CPP46590
      DATA DEBYE,BOHR / 2.541765480D+00 , 0.52917706D+00 /              CPP46600
    1 FORMAT(//,2X,' FIRST DERIVATIVES OF DIPOLE MOMENTS'/)             CPP46610
    2 FORMAT(//,2X,' UA MATRIX, IABC = ',I5,' IXYZ = ',I5/)             CPP46620
    3 FORMAT(//,2X,' DIPDA MATRIX'/)                                    CPP46630
    4 FORMAT(//,2X,' DIPDM MATRIX'/)                                    CPP46640
    5 FORMAT(//,2X,' DIPDC MATRIX'/)                                    CPP46650
    6 FORMAT(//,2X,' DIPDE MATRIX'/)                                    CPP46660
    7 FORMAT(//,2X,' DIPDN MATRIX'/)                                    CPP46670
    8 FORMAT(//,2X,' FIRST DERIVATIVES OF DIPOLE MOMENTS (IN DEBYE/BOHR)CPP46680
     *'/)                                                               CPP46690
    9 FORMAT(//,2X,' FIRST DERIVATIVES OF DIPOLE MOMENTS (IN DEBYE/A)'/)CPP46700
   10 FORMAT(2I5)                                                       CPP46710
   11 FORMAT(3F20.10)                                                   CPP46720
C                                                                       CPP46730
      CALL RFILE(ITAP44)                                                CPP46740
      DEB4=DEBYE*FOUR                                                   CPP46750
      DEB2=DEBYE*TWO                                                    CPP46760
      RBOHR=ONE/BOHR                                                    CPP46770
C     WRITE(6,1)                                                        CPP46780
C                                                                       CPP46790
      CALL ZERO(DIPDM,N3N*3)                                            CPP46800
      CALL ZERO(DIPDC,N3N*3)                                            CPP46810
C                                                                       CPP46820
      DO 101 IXYZ=1,3                                                   CPP46830
      IH=IHA(IXYZ+N3N)                                                  CPP46840
      CALL RREAD(ITAP44,HF(1,IXYZ),NTRI2,IH)                            CPP46850
  101 CONTINUE                                                          CPP46860
      DO 102 IABC=1,N3N                                                 CPP46870
      IU=IUA(IABC)                                                      CPP46880
      CALL RREAD(ITAP44,UA(1,1,IABC),NBASQ,IU)                          CPP46890
      IF(IPRNT.LE.3) GO TO 102                                          CPP46900
      WRITE(6,2) IABC                                                   CPP46910
      CALL MATOUT(UA(1,1,IABC),NBASIS,NBASIS,NBASIS,NBASIS,6)           CPP46920
  102 CONTINUE                                                          CPP46930
C                                                                       CPP46940
      DO 106 IABC=1,N3N                                                 CPP46950
      DO 105 IXYZ=1,3                                                   CPP46960
      VALUM=A00                                                         CPP46970
      DO 103 I=1,JOCC                                                   CPP46980
      ITYP=MOTYP(I)                                                     CPP46990
      FAC=FOCC(ITYP)*HALF                                               CPP47000
      DO 103 J=1,NBASIS                                                 CPP47010
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPP47020
      VALUM=VALUM+UA(J,I,IABC)*HF(IJ,IXYZ)*FAC                          CPP47030
  103 CONTINUE                                                          CPP47040
      VALUC=A00                                                         CPP47050
      DO 104 IJCI=1,ICIUNQ                                              CPP47060
      IEX=IPX(IJCI)                                                     CPP47070
      JEX=JPX(IJCI)                                                     CPP47080
      FAC=TWO                                                           CPP47090
      IF(IEX.EQ.JEX) FAC=ONE                                            CPP47100
      IJEX=IOFF(MAX0(IEX,JEX))+MIN0(IEX,JEX)                            CPP47110
      VALUC=VALUC+CIA(IEX,IABC)*CI(JEX)*HCIA(IJEX,IXYZ+N3N)*FAC         CPP47120
  104 CONTINUE                                                          CPP47130
      DIPDM(IXYZ,IABC)=VALUM                                            CPP47140
      DIPDC(IXYZ,IABC)=VALUC                                            CPP47150
  105 CONTINUE                                                          CPP47160
  106 CONTINUE                                                          CPP47170
C                                                                       CPP47180
      DO 107 I=1,3                                                      CPP47190
      DO 107 J=1,N3N                                                    CPP47200
      DIPDM(I,J)=-DIPDM(I,J)*DEB4                                       CPP47210
      DIPDC(I,J)=-DIPDC(I,J)*DEB2                                       CPP47220
  107 CONTINUE                                                          CPP47230
C                                                                       CPP47240
C   SUM UP CONTRIBUTIONS FROM MO CHANGE                                 CPP47250
      DO 108 I=1,3                                                      CPP47260
      DO 108 J=1,N3N                                                    CPP47270
      DIPDE(I,J)=DIPDA(I,J)+DIPDM(I,J)+DIPDC(I,J)                       CPP47280
  108 CONTINUE                                                          CPP47290
C*    IF(IPRNT.LE.2) GO TO 201                                          CPP47300
      WRITE(6,3)                                                        CPP47310
      CALL MATOUT(DIPDA,3,45,3,N3N,6)                                   CPP47320
      WRITE(6,4)                                                        CPP47330
      CALL MATOUT(DIPDM,3,45,3,N3N,6)                                   CPP47340
      WRITE(6,5)                                                        CPP47350
      CALL MATOUT(DIPDC,3,45,3,N3N,6)                                   CPP47360
C                                                                       CPP47370
C   SUM UP ALL CONTRIBUTIONS                                            CPP47380
  201 CONTINUE                                                          CPP47390
      DO 109 I=1,3                                                      CPP47400
      DO 109 J=1,N3N                                                    CPP47410
      DIPDT(I,J)=DIPDE(I,J)+DIPDN(I,J)                                  CPP47420
  109 CONTINUE                                                          CPP47430
C                                                                       CPP47440
C*    IF(IPRNT.LE.2) GO TO 202                                          CPP47450
      WRITE(6,6)                                                        CPP47460
      CALL MATOUT(DIPDE,3,45,3,N3N,6)                                   CPP47470
      WRITE(6,7)                                                        CPP47480
      CALL MATOUT(DIPDN,3,45,3,N3N,6)                                   CPP47490
  202 CONTINUE                                                          CPP47500
      WRITE(6,8)                                                        CPP47510
      CALL MATOUT(DIPDT,3,45,3,N3N,6)                                   CPP47520
C                                                                       CPP47530
      DO 110 I=1,3                                                      CPP47540
      DO 110 J=1,N3N                                                    CPP47550
      DIPDT(I,J)=DIPDT(I,J)*RBOHR                                       CPP47560
  110 CONTINUE                                                          CPP47570
      WRITE(6,9)                                                        CPP47580
      CALL MATOUT(DIPDT,3,45,3,N3N,6)                                   CPP47590
C                                                                       CPP47600
      ITAP17=17                                                         CPP47610
      REWIND ITAP17                                                     CPP47620
      WRITE(ITAP17,10) NATOMS,N3N                                       CPP47630
      WRITE(ITAP17,11) ((DIPDT(I,J),J=1,N3N),I=1,3)                     CPP47640
      REWIND ITAP17                                                     CPP47650
C                                                                       CPP47660
      CALL RCLOSE(ITAP44,3)                                             CPP47670
      RETURN                                                            CPP47680
      END                                                               CPP47690
      SUBROUTINE POLAR(HCIA,CIA,HF,UF)                                  CPP47700
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP47710
      DIMENSION HF(NTRI,3),UF(NBASIS,NBASIS,3)                          CPP47720
      DIMENSION HCIA(NCTRI,N3NP3),CIA(NCPEX,N3NP3)                      CPP47730
      DIMENSION POM(3,3),POC(3,3)                                       CPP47740
      DIMENSION POL(3,3),POLC(3,3),EIG(3),EIV(3,3),A(6),FV1(3),FV2(3)   CPP47750
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP47760
      COMMON/CIOFF/IIUNQ,IJUNQ,ICIUNQ,IPX(5050),JPX(5050)               CPP47770
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP47780
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP47790
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP47800
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP47810
      COMMON/START/FOCC(30),NSORB(30),NSTR(30),NEND(30),MOTYP(256)      CPP47820
      COMMON/TAPES/ITAP37,ITAP42,ITAP43,ITAP44,ITAP47                   CPP47830
      COMMON/CI101/IOCC,JOCC,KOCC                                       CPP47840
      COMMON/CI103/CI(100)                                              CPP47850
      DATA A00,HALF,ONE,TWO,FOUR / 0.0D+00 , 0.5D+00 , 1.0D+00 ,        CPP47860
     1                             2.0D+00 , 4.0D+00 /                  CPP47870
      DATA BOHR / 0.52917706D+00 /                                      CPP47880
    1 FORMAT(//,2X,' HF MATRIX, IXYZ = ',I5/)                           CPP47890
    2 FORMAT(//,2X,' UF MATRIX, IXYZ = ',I5/)                           CPP47900
    3 FORMAT(//,2X,' POM MATRIX'/)                                      CPP47910
    4 FORMAT(//,2X,' POC MATRIX'/)                                      CPP47920
    5 FORMAT(//,2X,' POLARIZABILITY TENSOR MATRIX (IN A.U.)'/)          CPP47930
    6 FORMAT(//,2X,' POLARIZABILITY TENSOR MATRIX (IN A**3)'/)          CPP47940
    7 FORMAT(//,2X,' A MATRIX'/)                                        CPP47950
    8 FORMAT(//,2X,' PRINCIPAL POLARIZABILITY TENSOR MATRIX (IN A.U.)'/)CPP47960
    9 FORMAT(//,2X,' PRINCIPAL POLARIZABILITY TENSOR MATRIX (IN A**3)'/)CPP47970
   10 FORMAT(//,2X,' EIGENVECTOR MATRIX'/)                              CPP47980
C                                                                       CPP47990
      CALL RFILE(ITAP44)                                                CPP48000
      A33=BOHR*BOHR*BOHR                                                CPP48010
C                                                                       CPP48020
C   CALL BACK HF MATRICES                                               CPP48030
      DO 101 IXYZ=1,3                                                   CPP48040
      IH=IHA(IXYZ+N3N)                                                  CPP48050
      IU=IUA(IXYZ+N3N)                                                  CPP48060
      CALL RREAD(ITAP44,HF(1,IXYZ),NTRI2,IH)                            CPP48070
      CALL RREAD(ITAP44,UF(1,1,IXYZ),NBASQ,IU)                          CPP48080
      IF(IPRNT.LE.3) GO TO 101                                          CPP48090
      WRITE(6,1) IXYZ                                                   CPP48100
      CALL PRINT(HF(1,IXYZ),NTRI,NBASIS,6)                              CPP48110
      WRITE(6,2) IXYZ                                                   CPP48120
      CALL MATOUT(UF(1,1,IXYZ),NBASIS,NBASIS,NBASIS,NBASIS,6)           CPP48130
  101 CONTINUE                                                          CPP48140
C                                                                       CPP48150
C:::B (G) COORDINATE:::                                                 CPP48160
      DO 105 IXYZ=1,3                                                   CPP48170
C                                                                       CPP48180
C:::A (F) COORDINATES:::                                                CPP48190
  301 CONTINUE                                                          CPP48200
      DO 104 JXYZ=1,3                                                   CPP48210
      VALUP=A00                                                         CPP48220
      DO 102 I=1,JOCC                                                   CPP48230
      ITYP=MOTYP(I)                                                     CPP48240
      FAC=FOCC(ITYP)*HALF                                               CPP48250
      DO 102 J=1,NBASIS                                                 CPP48260
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)                                      CPP48270
      VALUP=VALUP-UF(J,I,IXYZ)*HF(IJ,JXYZ)*FAC                          CPP48280
  102 CONTINUE                                                          CPP48290
      VALUC=A00                                                         CPP48300
      DO 103 IJCI=1,ICIUNQ                                              CPP48310
      IEX=IPX(IJCI)                                                     CPP48320
      JEX=JPX(IJCI)                                                     CPP48330
      FAC=TWO                                                           CPP48340
      IF(IEX.EQ.JEX) FAC=ONE                                            CPP48350
      IJEX=IOFF(MAX0(IEX,JEX))+MIN0(IEX,JEX)                            CPP48360
      VALUC=VALUC-CIA(IEX,IXYZ+N3N)*CI(JEX)*HCIA(IJEX,JXYZ+N3N)*FAC     CPP48370
  103 CONTINUE                                                          CPP48380
      POM(IXYZ,JXYZ)=VALUP*FOUR                                         CPP48390
      POC(IXYZ,JXYZ)=VALUC*TWO                                          CPP48400
  104 CONTINUE                                                          CPP48410
  105 CONTINUE                                                          CPP48420
C                                                                       CPP48430
      DO 106 IXYZ=1,3                                                   CPP48440
      DO 106 JXYZ=1,3                                                   CPP48450
      POL(IXYZ,JXYZ)=POM(IXYZ,JXYZ)+POC(IXYZ,JXYZ)                      CPP48460
  106 CONTINUE                                                          CPP48470
      WRITE(6,3)                                                        CPP48480
      CALL MATOUT(POM,3,3,3,3,6)                                        CPP48490
      WRITE(6,4)                                                        CPP48500
      CALL MATOUT(POC,3,3,3,3,6)                                        CPP48510
      WRITE(6,5)                                                        CPP48520
      CALL MATOUT(POL,3,3,3,3,6)                                        CPP48530
      DO 107 I=1,3                                                      CPP48540
      DO 107 J=1,3                                                      CPP48550
      POLC(I,J)=POL(I,J)*A33                                            CPP48560
  107 CONTINUE                                                          CPP48570
      WRITE(6,6)                                                        CPP48580
      CALL MATOUT(POLC,3,3,3,3,6)                                       CPP48590
C                                                                       CPP48600
      IJ=0                                                              CPP48610
      DO 108 I=1,3                                                      CPP48620
      DO 108 J=1,I                                                      CPP48630
      IJ=IJ+1                                                           CPP48640
      A(IJ)=POL(I,J)                                                    CPP48650
  108 CONTINUE                                                          CPP48660
      WRITE(6,7)                                                        CPP48670
      CALL PRINT(A,6,3,6)                                               CPP48680
      CALL RSP(3,3,6,A,EIG,1,EIV,FV1,FV2)                               CPP48690
C                                                                       CPP48700
      CALL ZERO(POL,9)                                                  CPP48710
      DO 110 I=1,3                                                      CPP48720
  110 POL(I,I)=EIG(I)                                                   CPP48730
      WRITE(6,8)                                                        CPP48740
      CALL MATOUT(POL,3,3,3,3,6)                                        CPP48750
      DO 111 I=1,3                                                      CPP48760
      DO 111 J=1,3                                                      CPP48770
      POLC(I,J)=POL(I,J)*A33                                            CPP48780
  111 CONTINUE                                                          CPP48790
      WRITE(6,9)                                                        CPP48800
      CALL MATOUT(POLC,3,3,3,3,6)                                       CPP48810
      WRITE(6,10)                                                       CPP48820
      CALL MATOUT(EIV,3,3,3,3,6)                                        CPP48830
C                                                                       CPP48840
      CALL RCLOSE(ITAP44,3)                                             CPP48850
      RETURN                                                            CPP48860
      END                                                               CPP48870
      SUBROUTINE WRBMAT(BB,NABC)                                        CPP48880
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP48890
      DIMENSION BB(NTRI,NABC)                                           CPP48900
      COMMON/BASIS/NBASIS,NTRI,NTRI2,NBFAO,NBATRI,NBASQ                 CPP48910
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP48920
      COMMON/POSIT/ISA(153),IHA(153),IFA(153),IEA(153),IBA(153),IUA(153)CPP48930
      COMMON/TAPES/ITAP37,ITAP42,ITAP43,ITAP44,ITAP47                   CPP48940
    1 FORMAT(//,2X,' B0 MATRIX'/)                                       CPP48950
C                                                                       CPP48960
      CALL RFILE(ITAP44)                                                CPP48970
      DO 101 IABC=1,NABC                                                CPP48980
      IB=IBA(IABC)                                                      CPP48990
      CALL RREAD(ITAP44,BB(1,IABC),NTRI2,IB)                            CPP49000
  101 CONTINUE                                                          CPP49010
      WRITE(6,1)                                                        CPP49020
      CALL MATOUT(BB,NTRI,NABC,NTRI,NABC,6)                             CPP49030
      CALL RCLOSE(ITAP44,3)                                             CPP49040
C                                                                       CPP49050
      RETURN                                                            CPP49060
      END                                                               CPP49070
      SUBROUTINE THIRD(E1T,HCIA,CIA,NSECMX)                             CPP49080
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP49090
      DIMENSION E1T(N3N),HCIA(NCTRI*N3N),CIA(NCPEX*N3N)                 CPP49100
      COMMON/CISET/NCPEX,NORPEX,NCTRI                                   CPP49110
      COMMON/FUNCS/NTYPES,NTYPEP,NTYTRI,NTREAD,NATOMS,N3N,N3NP3         CPP49120
      COMMON/MFSEC/MASTER,NSECT                                         CPP49130
      COMMON/SIGNS/IOFF(256),IPRNT                                      CPP49140
      COMMON/TAPES/ITAP37,ITAP42,ITAP43,ITAP44,ITAP47                   CPP49150
    1 FORMAT(//,2X,' E1T MATRIX IN THIRD'/)                             CPP49160
    2 FORMAT(//,2X,' HCIA MATRIX IN THIRD'/)                            CPP49170
    3 FORMAT(//,2X,' CIA MATRIX IN THIRD'/)                             CPP49180
C                                                                       CPP49190
      IF(IPRNT.LE.2) GO TO 201                                          CPP49200
      WRITE(6,1)                                                        CPP49210
      CALL MATOUT(E1T,3,NATOMS,3,NATOMS,6)                              CPP49220
      WRITE(6,2)                                                        CPP49230
      CALL MATOUT(HCIA,NCTRI,N3N,NCTRI,N3N,6)                           CPP49240
      WRITE(6,3)                                                        CPP49250
      CALL MATOUT(CIA,NCPEX,N3N,NCPEX,N3N,6)                            CPP49260
C                                                                       CPP49270
  201 CONTINUE                                                          CPP49280
      CALL RFILE(ITAP44)                                                CPP49290
      N3N2=N3N*2                                                        CPP49300
      NPOINT=NSECMX                                                     CPP49310
      CALL RWRIT(ITAP44,E1T,N3N2,NPOINT)                                CPP49320
      NSEC3L=((N3N2-1)/NSECT+1)                                         CPP49330
      NPOINT=NPOINT+NSEC3L                                              CPP49340
      CALL RWRIT(ITAP44,HCIA,N3N2*NCTRI,NPOINT)                         CPP49350
      NSECNC=((N3N2*NCTRI-1)/NSECT+1)                                   CPP49360
      NPOINT=NPOINT+NSECNC                                              CPP49370
      CALL RWRIT(ITAP44,CIA,N3N2*NCPEX,NPOINT)                          CPP49380
C                                                                       CPP49390
      CALL RCLOSE(ITAP44,3)                                             CPP49400
      RETURN                                                            CPP49410
      END                                                               CPP49420
      SUBROUTINE RMREAD(IA,IADD)                                        CPP49430
      IMPLICIT INTEGER (A-Z)                                            CPP49440
      DIMENSION IA(1)                                                   CPP49450
      COMMON/MFSEC/MASTER,NSECT                                         CPP49460
C                                                                       CPP49470
      CALL RFILE(MASTER)                                                CPP49480
      CALL MREAD(IA,IADD)                                               CPP49490
      CALL RCLOSE(MASTER,3)                                             CPP49500
      RETURN                                                            CPP49510
C                                                                       CPP49520
      ENTRY RMWRIT(IA,IADD)                                             CPP49530
      CALL RFILE(MASTER)                                                CPP49540
      CALL MWRIT(IA,IADD)                                               CPP49550
      CALL RCLOSE(MASTER,3)                                             CPP49560
      RETURN                                                            CPP49570
      END                                                               CPP49580
      SUBROUTINE SCLMPY(A,NDIM,SCALE)                                   CPP49590
      IMPLICIT REAL*8 (A-H,O-Z)                                         CPP49600
      DIMENSION A(1)                                                    CPP49610
C                                                                       CPP49620
      DO 101 I=1,NDIM                                                   CPP49630
      A(I)=A(I)*SCALE                                                   CPP49640
  101 CONTINUE                                                          CPP49650
C                                                                       CPP49660
      RETURN                                                            CPP49670
      END                                                               CPP49680
