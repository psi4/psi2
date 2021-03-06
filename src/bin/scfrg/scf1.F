      SUBROUTINE SCF1(CM,NLAMDA,NV,SU,ITYP,C1,C2,NX,NORBS,ALPHA,BETA
     1,               NC,NO,S,H,P1,P2,D1,D2,SD,ALL,D1OLD
     2,          EIG,SF1,SF2,SF11,SF12,SF21,SF22,FC,FO,IBUFFR,BUFFER
     3,          ICM,ISUQ,IALPHA,IBETA,IEIG,UNIT
     4,          FSAVE,ERROR,NFOCK,NDIIS,SMHALF,U,EIGVAL,T1,T2,TRIANG
     5,          BMAT,
     6           wmot,sd2,iwmsiz)
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
cets022791
#include <error.h>
C
CTJL 1,               NC,NO,S,H,P1,P2,D1,D2,SD,ALL,D1OLD,UNIT
C
      integer mxcoef,nsfmax,msfru,nsym,nmat,ntotal
      INTEGER ALABEL,BLABEL
      INTEGER ICM(1),ISUQ(1),IALPHA(1),IBETA(1),IEIG(1)
      DIMENSION SD(1),BUFFER(1),IBUFFR(1),UNIT(NMAT)
      DIMENSION ALPHA(1),BETA(1),ALL(1)
      DIMENSION SU(MXCOEF),NX(NSFMAX),CM(2,MSFRU)
      DIMENSION NC(NSYM),NO(NSYM),NORBS(NSYM)
      DIMENSION NLAMDA(NSYM),NV(2,NSYM),C1(NSYM),C2(NSYM),ITYP(NSYM)
      DIMENSION D1(NMAT),D2(NMAT),D1OLD(NMAT),S(NMAT),H(NMAT),P1(NMAT)
      DIMENSION P2(NMAT),SF1(NMAT),SF2(NMAT)
      DIMENSION SF11(NMAT),SF12(NMAT),SF21(NMAT),SF22(NMAT)
      DIMENSION FSAVE(NMAT,NFOCK,NDIIS),ERROR(NMAT,NDIIS),SMHALF(NMAT)
      DIMENSION U(MXCOEF),EIGVAL(NTOTAL),T1(*),T2(*),TRIANG(NMAT)
      DIMENSION BMAT(*)
      dimension wmot(*),sd2(msfru)
      DIMENSION FC(2),FO(2),EIG(NTOTAL)
      INTEGER PRINT,NAME(6),FCB49(16)
      INTEGER START, END
      LOGICAL TWOCON,IPOPLE,EXTEST,FIRST
      CHARACTER*6 CALTYP,DERTYP,CITYP
      COMMON /DERT/ CALTYP,DERTYP,CITYP
      COMMON /POP/ IPOPLE
      COMMON /PRINTQ/ PRINT
      COMMON /TAPES/ NIN,IR,IOUT,ICHECK,IFLAG3,ITAPE,I30
      COMMON /GVB/ ALPHA1,ALPHA2,ALPHA3,IJUP,TWOCON
      COMMON /FLAGS/ IDEBUG,ISTO,IPUNCH,ISAVE,IRORDR
      COMMON /IFDMP/ DAMPSV,INFLG,MICMX,IDFIX
      COMMON /SCFIN/ REPNUC,IOPEN,ISCF,IPRCT,IVECT,NBFLD,ICUT
      COMMON /INPT/  NSYM,IPMAX,NMAT,ITMAX,NSFMAX,NTOTAL,MXCOEF,MSFRU
      COMMON /DIM21/ IPQ(256)
      COMMON /NEWRUN/ ITPREP,ITCON
      COMMON /T102/ NOUT,IHALF,ITHIRD,MAXNO,IOUTBF,KINADR
      COMMON /EXTR/ NMAT2,NMAT3,NMAT4,NMAT5,NMAT6,NMAT7
      COMMON /NEXTFG/ NXTFLG
      COMMON/TIMX/ITIME,ICONT0,JCOUNT,ITIME2,ICONT2,ITIME3,ICONT3,
     1ITIME4,ICONT4,ITIME1,ICONT1
      COMMON /HOLL/ BLABEL(20),ALABEL(26),FMT(10),JL
cets022791
      integer frdboo,frdc,errcod
      common /damp/ dam1,damp
      logical oldinp
      common /oldin/ oldinp
cets
CTJL
      COMMON /TJL/ ASHIFO,DIISCT,AERRO,AFAC
CTJL
      DATA  EPS2,   A0, APT1, APT2, A1S2,   A1,A2PT2,   A10
     1/1.0D-6,0.0D0,0.1D0,0.2D0,0.5D0,1.0D0,2.2D0,10.0D0/
      DATA AERRC,DIISER,OLDER,ASHIFC /1.0D-01,0.0D+00,1.0D+20,0.05D+00/
      DATA A2,A3S2,ITJL /2.0D+00,1.5D+00,0/
CTJL  DATA AERRO /1.0D-10/
CTJL  DATA ASHIFO /2.0D+00,1.5D+00,0,0.05D+00/
C
      EPSLON=1.0D-15
C
      IPRINT=-1
C FORM SMHALF MATRIX
      CALL SINV(S,SMHALF,U,EIGVAL,T1,T2,NSYM,NMAT,NLAMDA,IPRINT
     1,MXCOEF,NTOTAL,TRIANG,IOUT)
C
C     ----- CORE GUESS IF INFLG.LT.0 -----
C
      IF (INFLG.LT.0) THEN
      CALL NEWVC(SU,H,SMHALF,T1,T2,U,NLAMDA,NSYM,NMAT,MXCOEF,
     #            EIGVAL,NTOTAL,IOUT,TRIANG)
CP    WRITE(IOUT,11115)
CP115 FORMAT(' SCF VECTOR AFTER NEWVC ')
CP    CALL VECOUT(SU,EIGVAL,MXCOEF,NTOTAL,NLAMDA,NSYM,IOUT)
      INFLG=0
      END IF
C
      EXTEST=.FALSE.
      CUT=FLOAT(ICUT)
      IF (CUT .LT. 0.0) CUT=ABS(1/CUT)
      IF (INFLG .LT. 0) FIRST=.TRUE.
      SQDIF=1.0D+00
      DO 2000 I=1,NMAT
        D1OLD(I)=0.0D+00
        D1(I)=0.0D+00
        D2(I)=0.0D+00
 2000 CONTINUE
      THRSH=A10**(-ISTO)
      MXITR=ITMAX-1
      EPS=A10**(-ISCF)
      IF (TWOCON) THEN
        ISCFT=1
        EPS=0.1
        IVECTS=IVECT
        IVECT=-1
        KLI=0
        KL=0
        DO 2002 L=1,NSYM
          NN=NLAMDA(L)
          DO 2001 I=1,NN
            KL=KL+1
            IF (CM(2,KL).NE.0.0) THEN
              IF (KLI.NE.0) THEN
                KLII=KL
                LII=L
                GO TO 1040
              ELSE
                KLI=KL
                LI=L
              END IF
            END IF
 2001     CONTINUE
 2002   CONTINUE
      END IF
 1040 NOC=1
      IF (IOPEN.NE.0) NOC=2
      DO 1 I=1,NSYM
 1    NV(2,I)=0
      NEXT=0
      ELEC=A0
      KL=0
      DO 3 I=1,NSYM
        C1(I)=A0
        C2(I)=A0
        NN=NLAMDA(I)
        DO 2 J=1,NN
          KL=KL+1
          IF(CM(1,KL).NE.A0) C1(I)=CM(1,KL)
          IF(CM(2,KL).NE.A0) C2(I)=CM(2,KL)
 2      CONTINUE
        IF (TWOCON) THEN
          IF (I.EQ.LI) C1I=C1(I)
          IF (I.EQ.LII) C1II=C1(I)
        END IF
        DEN=C1(I)-C2(I)
        IF(DEN.EQ.A0) GO TO 3
        C1(I)=C1(I)/DEN
        C2(I)=C2(I)/DEN
 3    CONTINUE
      DO 5 M=1,NOC
        KL=0
        DO 5 I=1,NSYM
          NZCOUN=0
          NN=NLAMDA(I)
          DO 4 J=1,NN
            KL=KL+1
            IF(CM(M,KL).EQ.A0) GO TO 4
            ELEC=ELEC+CM(M,KL)
            NZCOUN=NZCOUN+1
 4        CONTINUE
 5    NV(M,I)=NZCOUN
C
      CALL SCHMIT(NSYM,NV,0,THRSH,NLAMDA,SU,S,SD)
C
CX    KK=1
CX    DO 15 L=1,NSYM
CX      NN=NLAMDA(L)
CX      NNV=NORBS(L)
CX      DO 6 I=1,NN
C6      NX(I)=I
CX      IF(NNV.EQ.0) GO TO 9
CX      DO 8 I=1,NNV
CX        XMAX=A0
CX        DO 7 J=1,NN
CX          IF(XMAX.GT.ABS(SU(KK))) GO TO 7
CX          IF(NX(J).EQ.0) GO TO 7
CX          XMAX=ABS(SU(KK))
CX          JMAX=J
C7        KK=KK+1
CX        NX(JMAX)=0
C8      CONTINUE
C9      IU=NN-NNV
CX      IF(IU.EQ.0) GO TO 14
CX      JMAX=1
CX      DO 13 I=1,IU
C10       IF(NX(JMAX).NE.0) GO TO 11
CX        JMAX=JMAX+1
CX        GO TO 10
C11       DO 12 J=1,NN
CX          SU(KK)=A0
CX          IF(J.EQ.JMAX) SU(KK)=A1
C12       KK=KK+1
CX        JMAX=JMAX+1
C13     CONTINUE
C14     CONTINUE
C15   CONTINUE
C
CX    CALL SCHMIT(NSYM,NV,1,THRSH,NLAMDA,SU,S,SD)
C
CTEST RETURN
C
C-----------------------------------------------------------SCF2--------
C
CTEST ENTRY SCF2(EIG,SF1,SF2,SF11,SF12,SF21,SF22,FC,FO,IBUFFR,BUFFER)
C
      IF(IDEBUG.NE.0.OR.IVECT.GT.0)CALL CPRINT(NEXT,NV,NLAMDA,ITYP,
     1SU,EIG)
      IF (TWOCON) THEN
        E1CIP=0.0D+00
        E2CIP=0.0D+00
        KEYCIF=0
      END IF
      E1=A0
      E2=A0
      ELAST=A0
      DE=A0
      DEP=A0
      DEAVG=A0
      INEW=0
      DE1=A0
      KEYF=0
cets022791
      if(oldinp) READ (IR,16) DAMP, DAM1
   16 FORMAT(2F20.10)
C***************** TJL; CHANGED TO ACCOMODATE THE DIIS ROUTINE FOR
C   WHICH LEVEL SHIFTING IS NOT DESIRABLE.
      IF(IOPEN.EQ.0.AND.DABS(DAM1).GT.A0) DAM1 = APT1
C*****************
      WRITE(IOUT,240) DAM1,DAMP
      IF (IFLAG3.GT.0) WRITE (ICHECK,240) DAM1,DAMP
240   FORMAT(/,' LEVEL SHIFT  = ',F11.6,/,' DAMPING      = ',F11.6,/)
      IF(ITCON.EQ.0) GO TO 18
      READ(ITCON)E1,E2,DE,DE1,KEYF,DEAVG,INEW,DEP,ITMAXO,
     1(SF11(I),I=1,NMAT),(SF12(I),I=1,NMAT),(SF21(I),I=1,NMAT),
     2(SF22(I),I=1,NMAT),DAMPSV
C     REWIND ITCON
      CALL RFILE(ITCON)
CTJL  CALL SREW(ITCON)   THIS IS DONE IN RFILE
      IF(IDFIX.NE.0) DAMPSV=DAMP
      WRITE (ICHECK,17)
 17   FORMAT (38H ITERATIONS CONTINUED FROM EARLIER RUN)
      MXITR=ITMAXO+MXITR
      NEXT=NEXT+ITMAXO
 18   WRITE (IOUT,19)
      IF(IFLAG3.GT.0) WRITE(ICHECK,19)
C19   FORMAT (10H ITERATION,6X,6HENERGY,14X,2HDE,12X,5HSQDIF,12X,4HDAMP,
c19   FORMAT(10H ITERATION,6X,6HENERGY,14X,2HDE,12X,5HSQDIF,12X,6HDIISER
c    *,8X,19HONE-ELECTRON ENERGY,3X,19HTWO-ELECTRON ENERGY/)
C
C
 1210 IF (TWOCON) THEN
        CM(1,KLII)=2.0D+00
        CALL DMAT (NLAMDA,CM,SU,1,D1)
        CM(1,KLII)=0.0D+00
        CM(1,KLI)=2.0D+00
        CALL DMAT (NLAMDA,CM,SU,1,D2)
        DO 2003 I=1,NMAT
          AMAX=D1(I)
          TERM=D2(I)
          D1(I)=AMAX+TERM
          D2(I)=AMAX-TERM
 2003   CONTINUE
        CALL CONTR2 (BUFFER,IBUFFR,D1,D2,P1,P2,NMAT,.FALSE.)
        ECI1=0.0D+00
        ECI2=0.0D+00
        ECI3=0.0D+00
        DO 2004 I=1,NMAT
          ECI1=ECI1+D2(I)*H(I)
          ECI2=ECI2+D2(I)*P1(I)
          ECI3=ECI3+D2(I)*P2(I)
 2004   CONTINUE
        TERM=0.5*(ECI1+0.5*ECI2)
        HOFF=-0.25*ECI3
        AMAX=DSQRT(TERM*TERM+HOFF*HOFF)
        IF (TERM.LT.0.0D+00) AMAX=-AMAX
        CN=DSQRT((AMAX+TERM)/(AMAX+AMAX))
        SN=HOFF/(CN*(AMAX+AMAX))
        IF (TERM.GE.0.0D+00) THEN
          CI1=CN
          CI2=-SN
        ELSE
          CI1=SN
          CI2=CN
        END IF
C
        WRITE (IOUT,900) CI1,CI2
        IF (IFLAG3.GT.0) WRITE (IFLAG3,900) CI1,CI2
  900   FORMAT (' CI COEFFICIENTS ',2E18.8)
C
        CM(1,KLI)=0.0D+00
        CM(2,KLI)=2.0D+00*CI1**2
        CM(2,KLII)=2.0D+00*CI2**2
        DEN=C1I-CM(2,KLI)
        C1(LI)=C1I/DEN
        C2(LI)=CM(2,KLI)/DEN
        DEN=C1II-CM(2,KLII)
        C1(LII)=C1II/DEN
        C2(LII)=CM(2,KLII)/DEN
        ALPHA1=1.0D+00-2.0D+00/CM(2,KLI)
        ALPHA2=-1.0D+00/(CI1*CI2)
        ALPHA3=1.0D+00-2.0D+00/CM(2,KLII)
        KEYF=0
      END IF
C
C
   20 CONTINUE
      IF (FIRST) THEN
        FIRST=.FALSE.
      ELSE
CSEL       CALL ETIME2(START)
        CALL DMAT1(NLAMDA,CM,SU,1,D1,D2,IOPEN,SQDIF,NEXT,ALL,ICOUNT)
CSEL       CALL ETIME2(END)
        ITIME2=ITIME2+END-START
        IF(IOPEN.EQ.0)GO TO 39
CSEL       CALL ETIME2(START)
        CALL DMAT1(NLAMDA,CM,SU,2,D1,D2,IOPEN,SQDIF,NEXT,ALL,ICOUNT)
CSEL       CALL ETIME2(END)
        ITIME2=ITIME2+END-START
        DO 40 I=1,NMAT
          D1(I)=D1(I)+D2(I)
  40    CONTINUE
        IF(IDEBUG.NE.0)PRINT 75,(D1(I),I=1,NMAT)
  39    CALL CONTR2 (BUFFER,IBUFFR,D1,D2,P1,P2,NMAT,TWOCON)
      END IF
      IF (IDEBUG.NE.0) THEN
        WRITE (IOUT,901) P1
  901   FORMAT (/,' CLOSED SHELL TWO-ELECTRON CONTRIBUTION TO FOCK '
     #  ,            'MATRIX',//,(1X,7E18.8))
        IF (IOPEN.NE.0) WRITE (IOUT,902) P2
  902   FORMAT (//,' OPEN-SHELL TWO-ELECTRON CONTRIBUTION TO FOCK '
     #  ,            'MATRIX',//,(1X,7E18.8))
      END IF
      KEY=0
      E1P=E1
      E2P=E2
      ELASTP=ELAST
      DE1P=DE1
      DEP=DE
      E1=A0
      E2=A0
      IF(NEXT.EQ.1) DEAVG=ABS(DE)
      IF(IOPEN.NE.0) GO TO 22
      DO 21 I=1,NMAT
        E1=E1+D1(I)*H(I)
  21  E2=E2+D1(I)*P1(I)
      GO TO 24
 22   DO 23 I=1,NMAT
        E1=E1+D1(I)*H(I)
  23  E2=E2+D1(I)*P1(I)-D2(I)*P2(I)
 24   E2=E2*A1S2
      ELAST=E1+E2+REPNUC
      IF(E2.NE.A0) DE1=(E2-E2P)/ABS(E2)
      IF(INEW.EQ.0) GO TO 26
      IF(ABS(ELAST-ELASTP).GT.EPS) GO TO 25
      DEAVG=(ABS(DE1)+ABS(DE1P)+APT2*DEAVG)/A2PT2
      DE=DE1
      DEP=DE1P
      GO TO 27
 25   INEW=0
 26   DE=ELAST-ELASTP
      IF(NEXT.GE.2) DEAVG=(ABS(DE)+ABS(DEP)+APT2*DEAVG)/A2PT2
      IF(NEXT.GE.2.AND.DEAVG.LT.EPS2) INEW=1
      IF(NEXT.EQ.0) DE=A0
 27   NEXT=NEXT+1
CTJL  WRITE(IOUT,72) NEXT,ELAST,DE,SQDIF,DAMP,E1,E2
c
clj ------------------------------------------ old output with E1 and E2
c     WRITE(IOUT,72) NEXT,ELAST,DE,SQDIF,DIISER,E1,E2
clj This is the iteration header
c19   FORMAT('ITER','        ENERGY     ',
c    &       '       DE  ','     SQDIF ','    DIISER ',
c    &       ' ONE-ELECTRON ENERGY',' TWO-ELECTRON ENERGY',/)
c
clj and this is for each iteration
c72   FORMAT (I4,1X,G19.11,3(1X,E10.3),2(1X,G19.11))
clj ------------------------------------------   new output without E1 and E2
      WRITE(IOUT,72) NEXT,ELAST,DE,SQDIF,DIISER
clj This is the iteration header
 19   FORMAT('  ITER','        ENERGY     ',
     &       '       DE  ','     SQDIF ','    DIISER ')
clj and this is for each iteration
 72   FORMAT (1X,I5,1X,G19.11,3(1X,E10.3),2(1X,G19.11))
c
C01-16-89  BELOW ADDED TO FORCE WRITE TO FILE6 AT EACH SCF ITERATION
      call fioflu(iout)
      CALL RCLOSE(I30,3)
      CALL RFILE(I30)
C01-16-89  ABOVE ADDED TO FORCE WRITE TO FILE6 AT EACH SCF ITERATION
      IF(IFLAG3.GT.0)WRITE (ICHECK,72) NEXT,ELAST,DE,SQDIF,DAMP,E1,E2
      IF(ABS(E1-E1P).LT.ABS(E1*EPS).AND.ABS(E2-E2P).LT.ABS(E2*EPS))
     1KEY=1
      IF(IDFIX.NE.0.OR.NEXT.LT.3) GO TO 28
CTJL  CALL CNTRL (DE,DEP,DEAVG,DAMP,EPS,INEW)
 28   IF(KEYCNV.EQ.1) GO TO 31
      IFEXT=0
      CALL EXTRL(NEXT,NMAT,P1,SF12,SF22,IPRCT,DAMP,IFEXT)
      IF(IOPEN.EQ.0) GO TO 29
      CALL EXTRL(NEXT,NMAT,P2,SF11,SF21,IPRCT,DAMP,IFEXT)
 29   IF(IFEXT.EQ.1) WRITE(IOUT,30)
      IF(IFEXT.EQ.1.AND.IFLAG3.GT.0) PRINT 30
 30   FORMAT (28X,10(2H* ),21HEXTRAPOLATE P MATRIX ,10(2H* ))
 31   CONTINUE
      IF(IOPEN.NE.0) GO TO 33
      DO 32 I=1,NMAT
  32  SF2(I)=H(I)+P1(I)
      GO TO 35
 33   DO 34 I=1,NMAT
        SF2(I)=H(I)+P1(I)
  34  SF1(I)=H(I)+P1(I)-P2(I)
 35   IF(IDEBUG.EQ.0) GO TO 36
      PRINT 73,(SF2(I),I=1,NMAT)
      IF(IOPEN.NE.0) PRINT 74,(SF1(I),I=1,NMAT)
 36   KK=0
C
      IF (IOPEN.EQ.0) THEN
      IF (NEXT.GT.NXTFLG) THEN
C***************** TJL; CHANGED TO ACCOMODATE THE DIIS ROUTINE FOR
C   WHICH LEVEL SHIFTING IS NOT DESIRABLE.
CTJL  IF(NEXT.EQ.1) DAM1 = A0
      CALL DIIS(SF2,SF1,SU,FSAVE,ERROR,SMHALF,U,T1,T2,BMAT
     1,NLAMDA,CM,C1,C2
     2,DIISER
     3,NMAT,MXCOEF,NSYM,IOPEN,NDIIS,NFOCK,IPRINT,IOUT)
CTJL  WRITE(IOUT,11111) DIISER
C1111 FORMAT(' DIISER AFTER CALL TO DIIS =',E20.12)
CTJL  IF(DIISER.LT.AERR) DAM1 = APT1
      IF(IOPEN.NE.0) THEN
      IF(SQDIF.LT.AERRO) THEN
      WRITE(*,*) ' LEVEL SHIFT = ',DAM1
CTJL  IF(MOD(ITJL,5).EQ.0) DAM1 = DAM1/A2
      IF(ITJL.EQ.0) DAM1 = DAM1/AFAC
      ITJL = ITJL + 1
      ELSE
      IF(DIISER.GT.OLDER.AND.NEXT.GT.1) DAM1 = DAM1 + ASHIFO
      END IF
      ELSE
      IF(DIISER.GT.OLDER) DAM1 = DAM1 + ASHIFC
      IF(DIISER.LT.AERRC) DAM1 = A0
      END IF
      OLDER = DIISER
      END IF
      END IF
C
      L0=0
      ISU=0
      DO 59 L=1,NSYM
        NN=NLAMDA(L)
        IF (NN.EQ.0) GOTO 59
        NNV=NV(1,L)+NV(2,L)
        K0P=L0
        IF(NNV.NE.0.OR.KEYCNV.EQ.1.OR.NEXT.GT.MXITR) GO TO 38
CWA     DO 37 I=1,NN
C         K0P=K0P+1
C37     EIG(K0P)=A0
        GO TO 58
 38     MIC=MAX0(1,MICMX)
        IF(NV(2,L).EQ.0.OR.NV(1,L).EQ.0.OR.KEYCNV.EQ.1)MIC=1
        DO 57 MICRO=1,MIC
          K0P=L0
C
C
          CALL TRITRN(NN,ISU,KK,FC,SF2,SU,SD)
C
          IF(NV(2,L).EQ.0) GO TO 210
C
          CALL TRITRN(NN,ISU,KK,FO,SF1,SU,SD)
          IJ=1
          IF ((FLOAT(ISCFT)/CUT.GE.FLOAT(ISCF)) .AND. TWOCON) GOTO
     *    120
          IF(KEYCNV.NE.1.AND. NEXT .LE. MXITR) GOTO 77
          IF (.NOT.TWOCON) GOTO 120
          IF ((FLOAT(ISCFT)/CUT.LT.FLOAT(ISCF)).AND.NEXT.LE.MXITR)
     *    GOTO 77
  120     DO 80 I=1,NN
            DO 90 J=1,I
              IF(CM(1,L0+I).EQ.A0)GO TO 100
              IF(CM(2,L0+J).EQ.A0)GO TO 110
 140          FC(IJ)=C1(L)*FC(IJ)-C2(L)*FO(IJ)
              GO TO 110
 100          IF(CM(1,L0+J).EQ.A0) FC(IJ)=FO(IJ)
              IF(CM(2,L0+I).NE.A0)GO TO 140
 110          IJ=IJ+1
  90        CONTINUE
  80      CONTINUE
          GO TO 210
  77      DO 200 I=1,NN
            DO 200 J=1,I
              CIP=CM(1,L0+I) + CM(2,L0+I)
              CJP=CM(1,L0+J) + CM(2,L0+J)
              IF(CIP.EQ.A0 .OR. CJP.EQ.A0) GO TO 150
 145          FC(IJ)=C1(L)*FC(IJ)-C2(L)*FO(IJ)
              GO TO 190
 150          IF( (CIP+CJP) .EQ. A0) GO TO 145
              IF( (CM(2,L0+I)+CM(2,L0+J)).NE. A0) FC(IJ)=FO(
     *        IJ)
 190          IJ=IJ+1
 200      CONTINUE
  210     DAM2=DAM1
          II=0
          CIP=CM(1,L0+1)+CM(2,L0+1)
          DO 51 I=1,NNV
            II=II+I
            CJP=CM(1,L0+I)+CM(2,L0+I)
            IF(CIP.EQ.CJP) GO TO 51
            CIP=CJP
            DAM2=A1S2*DAM2
   51     IF(CIP.NE.A0) FC(II)=FC(II)-DAM2
CWA
          CALL ERDUW(FC,SU(ISU+1),NN,EPSLON,CM(1,L0+1),IRORDR)
C     IF(KEYCNV.EQ.1.OR.NEXT.GT.MXITR) GO TO 55
          DAM2=DAM1
          II=0
          CIP=CM(1,L0+1)+CM(2,L0+1)
          DO 53 I=1,NNV
            II=II+I
            CJP=CM(1,L0+I)+CM(2,L0+I)
            IF(CIP.EQ.CJP) GO TO 53
            CIP=CJP
            DAM2=A1S2*DAM2
   53     IF(CIP.NE.A0) FC(II)=FC(II)+DAM2
 55       DO 56 J=1,NN
            K0P=K0P+1
            EIG(K0P)=FC(IPQ(J+1))
 56       CONTINUE
 57     CONTINUE
 58     L0=L0+NN
        KK=KK+IPQ(NN+1)
        ISU=ISU+NN*NN
 59   CONTINUE
      IF(NEXT.GT.MXITR) THEN
        IVECT=IVECTS
        GOTO 60
      END IF
      IF(KEYCNV.EQ.1) GO TO 660
      SUM=0.0D0
      DO 2005 I=1,NMAT
        SUM=SUM+(D1(I)-D1OLD(I))**2
        D1OLD(I)=D1(I)
 2005 CONTINUE
      SQDIF=DSQRT(SUM)/NMAT
      IF (.NOT. IPOPLE) EXTEST=.TRUE.
      IF (IPOPLE .AND. ICOUNT .GT. 1) EXTEST=.TRUE.
      IF (INFLG .LT. 0 .AND. NEXT .LT. 3) EXTEST=.FALSE.
      IF (SQDIF .LE. 10.**(-ISCF) .AND. EXTEST) KEYCNV=1
      IF ((SQDIF .LE. 1.0D-05 .AND. NEXT .GT.2).OR. IPOPLE) THEN
        DAMP=0.0D+00
        DAMPSV=0.0D+00
      END IF
      IF ((SQDIF.LE.EPS).AND.TWOCON) GOTO 660
      CALL SCHMIT(NSYM,NV,1,THRSH,NLAMDA,SU,S,SD)
      IF(IVECT.GT.0)CALL CPRINT(NEXT,NV,NLAMDA,ITYP,SU,EIG)
      GO TO 20
  660 IF (TWOCON) THEN
C     PRINT, ' IN TWOCON LOOP. KEYCIF=',KEYCIF,'  ISCFT=',ISCFT,' EPS=',
C    1 EPS,' KEYCNV=',KEYCNV
        KEYCIF=KEYCIF+1
        IF (FLOAT(ISCFT)/CUT.GT.FLOAT(ISCF)) THEN
          IVECT=IVECTS
          GO TO 61
        END IF
        IF (KEYCIF.GT.1) THEN
          KEYCIF=0
          ISCFT=ISCFT+1
          EPS=10.0**(-FLOAT(ISCFT)/CUT)
          E1CIP=E1
          E2CIP=E2
          CALL SCHMIT (NSYM,NV,1,THRSH,NLAMDA,SU,S,SD)
          KEYCNV=0
          GO TO 1210
        ELSE
          GOTO 20
        END IF
      END IF
      GO TO 61
 60   WRITE (ICHECK,68)
   61 CONTINUE
CTJL  WRITE(ICHECK,*) '  CALLING OUTPTX'
      IF(JL.EQ.0) JL = 1
      CALL OUTPTX(ISUQ, ALPHA, BETA,NLAMDA,NORBS,NC,NO,ITYP,CM,IEIG
     1,           ELAST,TWOCON,KLI,KLII,CI1,CI2,IALPHA,IBETA,
     2            sf1,sf2,fc,fo,wmot,sd2,iwmsiz)
C
      IF(ITPREP.EQ.0) GO TO 63
      WRITE(ITPREP)E1,E2,ELAST,DE,DE1,KEYF,DEAVG,INEW,DEP,NEXT,
     1(SF11(I),I=1,NMAT),(SF12(I),I=1,NMAT),(SF21(I),I=1,NMAT),
     2(SF22(I),I=1,NMAT),DAMPSV
C     REWIND ITPREP
      CALL RFILE(ITPREP)
CTJL  CALL SREW(ITPREP)      THIS IS DONE IN RFILE
      WRITE (IOUT,62) ITPREP
      IF(IFLAG3.GT.0) WRITE(ICHECK,62) ITPREP
 62   FORMAT (53H ITERATION CONTINUATION DATA HAS BEEN WRITTEN ON UNIT,I
     13)
 63   IF(IVECT.GT.0) THEN
        PRINT 64
 64     FORMAT (1H )
        CALL CPRINT(NEXT,NV,NLAMDA,ITYP,SU,EIG)
      END IF
      IF(IVECT.GE.0) THEN
        CALL VIRTUL(NLAMDA,ITYP,SU,CM,EIG)
      ELSE IF(IVECT.GT.-2) THEN
        CALL CPRINT(NEXT,NV,NLAMDA,ITYP,SU,EIG)
      END IF
CTJL  CALL RSETSA(NOUT,KINADR)
C
C  IN CASE OF ORBITAL ROTATION FOR CI, NEED TO SAVE H
C
      DO 297 II=1,NMAT
        UNIT(II)=H(II)
  297 CONTINUE
C
CTJL  CALL SREAD(NOUT,H,INTOWP(NMAT))
CTJL  CALL SREW(NOUT)
      OVLP=A0
      EKIN=A0
      DO 66 I=1,NMAT
        DT=D1(I)
        OVLP=OVLP+DT*S(I)
 66   EKIN=EKIN+DT*H(I)
      OVLP=OVLP/ELEC
      EPOT=ELAST-EKIN
      VIRIAL=EPOT/ELAST
C
      IF (TWOCON) WRITE (IOUT,904) CI1,CI2,2*CI1**2,2*CI2**2
  904 FORMAT (//,' CI COEFFICIENTS:',21X,2F18.8,/,' ORBITAL OCCUPANCIES'
     #,       ' FOR GVB ORBITALS:',2F18.8,//)
C
      WRITE(IOUT,71) ELAST,EKIN,EPOT,VIRIAL,OVLP
      IF(KEYCNV.NE.1) WRITE (IOUT,68)
      IF(IFLAG3.GT.0) THEN
        IF (TWOCON) WRITE (ICHECK,904) CI1,CI2,2*CI1**2,2*CI2**2
        WRITE(ICHECK,71) ELAST,EKIN,EPOT,VIRIAL,OVLP
        IF(KEYCNV.NE.1) WRITE (ICHECK,68)
      END IF
C
C *********  THIS DOES THE ORBITAL ROTATION FOR A CI GRADIENT *******
C
      CALL LOCATE(5,'# TFOCK ##',IERROR)
C*********************************************************************
      irot=1
      IF (IERROR.EQ.0) then
        READ(5,922) CALTYP,CITYP,DERTYP

        if(cityp.eq.'SCF   '.or.dertyp.eq.'SECOND') irot=0
      else
        dertyp = 'FIRST '
        cityp = 'SCF   '

        errcod = frdc('WFN',cityp)
        errcod = frdc('DERTYPE',dertyp)

        if(dertyp.ne.'FIRST '.and.dertyp.ne.'NONE  ') irot=0
        if(cityp.eq.'SCF   ') irot=0

        errcod = frdboo('ROTATE',irot)
      endif
 922  FORMAT(A6,4X,A6,4X,A6)
CTJL  WRITE(*,*) ' THE MO FOCK MATRIX '
CTJL    CALL PRINT(FC,NN*(NN+1)/2,NN,6)
CTJL  CITYP = ' '
CTJL  WRITE(*,*) ' CITYPE ',CITYP
cets030291 do not rotate if scf, rotate for cc and ci
c      IF (CITYP.NE.'CI    '.AND.CITYP.NE.'GVBCI ') GOTO 1010
      if(cityp.eq.'SCF   '.or.irot.eq.0) goto 1010
      IF (IOPEN .EQ. 0) GOTO 1010
      KK=0
      ISU=0
      L0=0
      CALL DMAT1(NLAMDA,CM,SU,1,D1,D2,IOPEN,SQDIF,NEXT,ALL,ICOUNT)
      CALL DMAT1(NLAMDA,CM,SU,2,D1,D2,IOPEN,SQDIF,NEXT,ALL,ICOUNT)
      DO 1069 II=1,NMAT
 1069 D1(II)=D1(II)+D2(II)
cets071291      CALL CONTR2(BUFFER,IBUFFR,D1,D2,P1,P2,NMAT,TWOCON)
      call contr2g(buffer,ibuffr,d1,d2,p1,p2,nmat,twocon)
      DO 1023 II=1,NMAT
 1023 SF2(II)=UNIT(II)+P1(II)
C
      do 1000 i=1,nmat
         d1(i)=p1(i)-p2(i)+unit(i)
         d2(i)=p1(i)+p2(i)+unit(i)
 1000 continue

      itap79=79
      call rfile(itap79)
      ntmo2=ntotal*(ntotal+1)/2
      call zero(sf11,ntmo2)
cges leave space for one-e ham ints in mo basis
      call wwritw(itap79,sf11,intowp(ntmo2),1,i79)
      call wwritw(itap79,d1,intowp(nmat),i79,i79)
      call wwritw(itap79,d2,intowp(nmat),i79,i79)
      call rclose(itap79,3)

      DO 1059 L=1,NSYM
        NN=NLAMDA(L)
        IF (NN .EQ. 0) GOTO 1059
        CALL TRITRN(NN,ISU,KK,FC,SF2,SU,SD)
        IJ=0
        DO 1060 IOR=1,NN
          CIP=CM(1,L0+IOR)+CM(2,L0+IOR)
          DO 1060 JOR=1,IOR
            IJ=IJ+1
            IF (CIP .NE. CM(1,L0+JOR)+CM(2,L0+JOR)) FC(IJ)=0.0D+00
 1060   CONTINUE
C
C       WRITE(*,*) '  FC BEFORE:'
C       CALL PRINT(FC,NN*(NN+1)/2,NN,6)
CWA
        CALL ERDUW(FC,SU(ISU+1),NN,EPSLON,CM(1,L0+1),IRORDR)
C       CALL PRINT(FC,NN*(NN+1)/2,NN,6)
        DO 1062 II=1,NN
          IIJJ=II*(II+1)/2
 1062   EIG(II+L0)=FC(IIJJ)
        L0=L0+NN
        KK=KK+IPQ(NN+1)
        ISU=ISU+NN*NN
 1059 CONTINUE
      WRITE(3,1061)
      WRITE(ICHECK,1061)
 1061 FORMAT(/'  THE VECTOR HAS BEEN ROTATED FOR THE CI GRADIENT.'/)
      CALL OUTPTX(ISUQ,ALPHA,BETA,NLAMDA,NORBS,NC,NO,ITYP,CM,IEIG,
     1            ELAST,TWOCON,KLI,KLII,CI1,CI2,IALPHA,IBETA,
     2            sf1,sf2,fc,fo,wmot,sd2,iwmsiz)
      write(6,*)
      write(6,*)'MO rotation has been done'
      write(6,*)'and fock matrices have been written to file79 '
      WRITE(IOUT,71) ELAST,EKIN,EPOT,VIRIAL,OVLP

      CALL VIRTUL(NLAMDA,ITYP,SU,CM,EIG)
1010  CONTINUE
C
C     ***** THESE ADDED STATEMENTS ARE FOR OPEN-SHELL DERIVATIVE PROGRAM
C     *****                            Y. OSAMURA
C
      ITAP49=49
CCRBR CALL RFILE(ITAP49)    THE RFILE IS MOVED TO FRONT OF PROGRAM SO
CTJLRRCALL SREW(ITAP49)    THIS IS DONE IN RFILE       USE THIS SREW
      CALL SREW(ITAP49)
      CALL SWRIT(ITAP49,SU,INTOWP(MXCOEF))
      CALL SWRIT(ITAP49,EIG,INTOWP(NTOTAL))
      WRITE(ICHECK,*)'MXCOEF,NTOTAL,IOPEN,NSYM',MXCOEF,NTOTAL,IOPEN,NSYM
      IF (IOPEN .EQ. 0)  GO TO 600
      KK=0
      ISU=0
      DO 559 L=1,NSYM
        NN=NLAMDA(L)
        CALL TRITRN(NN,ISU,KK,FC,SF2,SU,SD)
        CALL TRITRN(NN,ISU,KK,FO,SF1,SU,SD)
        NDET=NN*(NN+1)/2
C
CC    WRITE(ITAP49) (FC(J),J=1,NDET)
CC    WRITE(ITAP49) (FO(J),J=1,NDET)
        CALL SWRIT(ITAP49,FC,INTOWP(NDET))
        CALL SWRIT(ITAP49,FO,INTOWP(NDET))
      WRITE(ICHECK,*) ' NSYM,NN,NDET,L ',NSYM,NN,NDET,L
        IF (IAND(PRINT,1024) .NE. 0) THEN
          WRITE(ICHECK,555) L,(FC(J),J=1,NDET)
          WRITE(ICHECK,554) L,(FO(J),J=1,NDET)
          PRINT 75,(D1(J), J=1,NMAT)
          PRINT 76,(D2(J), J=1,NMAT)
  555     FORMAT(' FC(',I3,')', (T11,10F12.6 ) )
  554     FORMAT(' FO(',I3,')', (T11,10F12.6 ) )
        END IF
C
        KK=KK+NDET
        ISU=ISU+NN*NN
  559 CONTINUE
C
  600 CONTINUE
      CALL SREW(ITAP49)
C
C     ... END OF INSERT CARDS FOR OPEN-SHELL DERIVATIVE ...
C
      IF (IAND(PRINT,512) .NE. 0) THEN
        PRINT 69,(P1(I),I=1,NMAT)
        PRINT 73,(SF2(I),I=1,NMAT)
        PRINT 75,(D1(I),I=1,NMAT)
        IF (IOPEN.EQ.0) GO TO 67
        PRINT 70,(P2(I),I=1,NMAT)
        PRINT 74,(SF1(I),I=1,NMAT)
        PRINT 76,(D2(I),I=1,NMAT)
      END IF
 67   RETURN
 68   FORMAT (30H CALCULATION HAS NOT CONVERGED)
 69   FORMAT (9H P-MATRIX/(7E18.8))
 70   FORMAT (9H Q-MATRIX/(7E18.8))
 71   FORMAT (//,'     TOTAL ENERGY       =',F18.10,/,
     1'     KINETIC ENERGY     =',F18.10,/,
     2'     POTENTIAL ENERGY   =',F18.10,/,
     3'     VIRIAL THEOREM     =',F18.10,/,
     4'     WAVEFUNCTION NORM  =',F18.10,/)
c72   FORMAT (1X,I5,5X,G19.11,3(2X,E14.5),2(3X,G19.11))
 73   FORMAT (25H CLOSED-SHELL FOCK MATRIX/(7E18.8))
 74   FORMAT (25H   OPEN-SHELL FOCK MATRIX/(7E18.8))
 75   FORMAT (19H      TOTAL DENSITY/(7E18.8))
 76   FORMAT (19H OPEN-SHELL DENSITY/(7E18.8))
      END
