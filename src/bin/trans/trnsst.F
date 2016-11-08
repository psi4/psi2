      SUBROUTINE TRNSST(A,B,C,U,IOFF,NSOFF,TTT,NSYMHF,NSO,CMAX,BMAX,
     1AMAX,MOINT,NSOMO,MOOFF,LENGTH,MOARR,MOBKT,MOLBL,MOSRT,BFDRT,
     2ORBSYM,IJGRP,IJADD,KADD,LADD,LENMO,ACRCY,LCIBUF,ITTT,IMOINT)
C
C  TRANSFORM INTEGRALS OF TYPE (AB,AB)
C
C
C  MOARR : MOBKT : MOLBL : MOSRT ARE IMPLICITLY EQUIVALENCED
C
      IMPLICIT REAL*8 (A-H,O,T-Z)
      IMPLICIT INTEGER (I-N,P,Q,R,S)
C
C
      integer norbs,nsymhf,numij,lenr,lenmo
      INTEGER TIME,ENDD,STARTT
      INTEGER H,CMAX,BMAX,AMAX,T,MOARR(LCIBUF),MOBKT(LENMO),MOLBL(LENMO)
      INTEGER BFDRT(NORBS),NSOMO(NSYMHF),MOOFF(NSYMHF),SYMORB
      INTEGER ORBSYM(NORBS),IJGRP(NUMIJ),IJADD(NUMIJ)
      INTEGER KADD(SYMORB),LADD(SYMORB),MULTS(8),LKUPSM(64)
      INTEGER AB,ASM,BSM,MINI,ITTT(1),IMOINT(1)
      REAL*8 SEC
      REAL*8 MOINT(LENR),MOSRT(LENMO),ACRCY
C
      COMMON /DIMS/   SYMORB,NUMIJ,NIJ,NBF,NORBS,NMAX,NBFORB,NSYM0
      COMMON /TAPES/  IT34,IT58,IT47,IT52,IT91,IT93,IT95
      COMMON /OUT/    IOUT6,ICHECK,IDBG,ICNT
      COMMON /BDIM/   IMAX,JMAX,KMAX,LENR,LENAOR,IBIG,IBIGR,LENPWS
C
      DIMENSION A(AMAX),B(BMAX),C(CMAX),IOFF(NIJ),NSOFF(NSYMHF)
      DIMENSION U(NBFORB),TTT(LENR),NSO(NSYMHF)
C
      DATA MULTS/0,8,16,24,32,40,48,56/
      DATA LKUPSM/1,2,3,4,5,6,7,8,2,1,4,3,6,5,8,7,3,4,1,2,7,8,5,6,4,3,
     A2,1,8,7,6,5,5,6,7,8,1,2,3,4,6,5,8,7,2,1,4,3,7,8,5,6,3,4,1,2,
     B8,7,6,5,4,3,2,1/
C
C
C
      WRITE(IOUT6,791)
  791 FORMAT (/,' TRANSFORMING INTEGRALS OF SYMMETRIES [AB;AB]:')
C
CDB   IDBG32=0
CDB   IF((IDBG.AND.32).NE.0) IDBG32=1
      LENOUT=I2SEC(LENMO*4+2)
      SRIN=0
      SWOUT=0
      AB=1
      DO 900 ASM=2,NSYMHF
           IF (NSO(ASM).EQ.0) THEN
                AB=AB+ASM
                GOTO 99
           END IF
           NA=NSOFF(ASM)
           NAMO=MOOFF(ASM)
           DO 800 BSM=1,ASM-1
                N105=0
                AB=AB+1
                IF (NSO(BSM).EQ.0) GOTO 88
C     PRINT,'  AB=',AB,'  ASM=',ASM,'  BSM=',BSM
                NB=NSOFF(BSM)
                NBMO=MOOFF(BSM)
                NSOMOB=NSOMO(BSM)*NSOMO(BSM)
                NSOB2=NSO(BSM)*NSO(BSM)
                NSOAB=NSOMO(ASM)*NSOMO(BSM)
                IPNT=0
                IF (NSO(ASM) .EQ. 0) THEN
                     MAXSZ=0
                ELSE
                     MAXSZ=(IOFF(NSO(ASM))+NSO(ASM))*NSOB2*2
                END IF
                NWRT=MAXSZ/LENGTH
                IF (NWRT*LENGTH .NE. MAXSZ) NWRT=NWRT+1
                NRED=1
C     PRINT,'  NWRT=',NWRT,'  MAXSZ=',MAXSZ
                CALL SREW(IT95)
                DO 1 II=1,LENGTH/2
                     MOINT(II)=0.0D+00
    1           CONTINUE
                IWRD=LENGTH
                IF (NWRT .EQ. 1) IWRD=MAXSZ
CDB   CALL RGETSA(IT93,IX)
                CALL SREAD(IT93,ITTT,IWRD)
                SRIN=SRIN+I2SEC(IWRD)
CDB   PRINT,'  READ IT93 IADR=',IX,'  WORDS=',IWRD
C     WRITE(IOUT6,'(7G18.8)') (TTT(IOP), IOP=1,MAXSZ/2)
                IF (NSOMO(ASM)*NSOMO(BSM) .EQ. 0) GOTO 88
                PQRS=0
                DO 200 P=1,NSO(ASM)
                     PR1=(IOFF(P)-1)*NSOB2
                     PNA=(P+NA-1)*NORBS+NAMO
                     E=1.0D+00
                     DO 2 II=1,CMAX
                          C(II)=0.0D+00
    2                CONTINUE
                     DO 13 Q=1,NSO(BSM)
                          QS1=(Q-1)*NSO(BSM)
                          QNB=(Q+NB-1)*NORBS+NBMO
                          DO 3 II=1,BMAX
                               B(II)=0.0D+00
    3                     CONTINUE
                          PR=PR1
                          DO 9 R=1,P
                               PR=PR+NSOB2
                               RNA=(R+NA-1)*NORBS+NAMO
                               DO 4 II=1,AMAX
                                    A(II)=0.0D+00
    4                          CONTINUE
                               PRQS=PR+QS1
                               T=NSO(BSM)
C     IF (P.EQ.R) T=Q
                               DO 6 S=1,T
C     PRQS=PRQS+1
                                    PQRS=PQRS+1
                                    IPNT=IOFF(P)*NSO(BSM)*NSO(BSM)+(Q-1)
     *                              *P*NSO(BSM)+(R-1)*NSO(BSM)+S-(NRED-1
     *                              )*LENGTH/2
                                    SNB=(S+NB-1)*NORBS+NBMO
                                    EE=E
                                    IF ((S.EQ.Q).AND.(P.EQ.R)) EE=0.5D+
     *                              00
                                    F=TTT(IPNT)*EE
                                    IF (F .EQ. 0.0D+00) GOTO 15
                                    DO 5 L=1,NSOMO(BSM)
                                         A(L)=A(L)+U(SNB+L)*F
    5                               CONTINUE
 15                                 CONTINUE
                                    IF ((IPNT .GE. LENGTH/2) .AND. (NRED
     *                               .LT. NWRT))THEN
                                         NRED=NRED+1
C     PRINT,'  IPNT=',IPNT,' NRED',NRED,' LENGTH/2',LENGTH/2
C     PRINT,'  P=',P,'  Q=',Q,'  R=',R,'  S=',S
                                         IWRD=LENGTH
C     IF (NRED*LENGTH .GT. MAXSZ) IWRD=MAXSZ-(NRED-1)*LENGTH
                                         IPNT=0
CDB   CALL RGETSA(IT93,IX)
                                         CALL SREAD(IT93,ITTT,IWRD)
                                         SRIN=SRIN+I2SEC(IWRD)
CDB   PRINT,'  READ IT93, IADR=',IX,'  WORDS=',IWRD
                                    END IF
    6                          CONTINUE
                               DO 8 K=1,NSOMO(ASM)
                                    URNAK=U(RNA+K)
                                    KL=(K-1)*NSOMO(BSM)
                                    DO 7 L=1,NSOMO(BSM)
                                         KL=KL+1
                                         B(KL)=B(KL)+URNAK*A(L)
    7                               CONTINUE
    8                          CONTINUE
    9                     CONTINUE
                          DO 12 J=1,NSOMO(BSM)
                               UQNBJ=U(QNB+J)
                               IF (UQNBJ .EQ. 0.0D+00) GOTO 25
                               JKL=(J-1)*NSOAB
                               DO 11 K=1,NSOMO(ASM)
                                    KL=(K-1)*NSOMO(BSM)
                                    DO 10 L=1,NSOMO(BSM)
                                         KL=KL+1
                                         JKL=JKL+1
                                         C(JKL)=C(JKL)+UQNBJ*B(KL)
   10                               CONTINUE
   11                          CONTINUE
 25                            CONTINUE
   12                     CONTINUE
   13                CONTINUE
                     IK=0
                     IJ=0
                     IPOUT=0
                     IF (N105 .NE. 0) THEN
                          CALL SREW(IT95)
                          CALL RGETSA(IT95,IADR)
                          CALL SREAD(IT95,IMOINT,LENGTH)
CDB   PRINT,'  READ IT95 FROM ADR',IADR
                          CALL  RSETSA(IT95,IADR)
                     END IF
                     DO 100 I=1,NSOMO(ASM)
                          UPNAI=U(PNA+I)
                          IJG=IOFF(I+NAMO)+NBMO
                          IJ=(I-1)*NSOMO(BSM)
                          IK1=(IOFF(I)-1)*NSOMOB
                          DO 90 J=1,NSOMO(BSM)
                               IJG=IJG+1
                               IJKG=IOFF(IJG)
                               IJ=IJ+1
                               IK=IK1
                               JL1=(J-1)*NSOMO(BSM)
                               JKL=(J-1)*NSOAB
                               IF (IPOUT+NSOAB+1 .GT. LENGTH/2) THEN
                                    N105=1
CDB   CALL RGETSA(IT95,IX)
CDB   PRINT,'  WRITE TO IT95,'
CDB   PRINT,'  IPOUT+NSOAB+1=',IPOUT+NSOAB+1,'  IADR=',IX
                                    CALL SWRIT(IT95,IMOINT,LENGTH)
                                    IPOUT=0
                                    IF (P .NE. 1) THEN
                                         CALL RGETSA(IT95,IADR)
                                         CALL SREAD(IT95,IMOINT,
     *                                   LENGTH)
                                         CALL  RSETSA(IT95,IADR)
                                    ELSE
                                         DO 60 II=1,LENGTH/2
                                              MOINT(II)=0.0D+00
   60                                    CONTINUE
                                    END IF
                               END IF
                               DO 80 K=1,I
                                    UPNAK=U(PNA+K)
                                    IJKLG=IJKG+IOFF(K+NAMO)+NBMO
                                    IK=IK+NSOMOB
                                    IJKL=IK+JL1
                                    H=NSOMO(BSM)
                                    IF (I.EQ.K) H=J
                                    DO 70 L=1,H
                                         IPOUT=IPOUT+1
                                         IJKLG=IJKLG+1
                                         LIJ=(L-1)*NSOAB+IJ
                                         IJKL=IJKL+1
                                         JKL=JKL+1
                                         MOINT(IPOUT)=MOINT(IPOUT)+UPNAI
     *                                   *C(JKL)+UPNAK*C(LIJ)
   70                               CONTINUE
   80                          CONTINUE
   90                     CONTINUE
  100                CONTINUE
                     IF ((IPOUT .NE. 0) .AND. (N105 .NE. 0))THEN
CDB   CALL RGETSA(IT95,IX)
CDB   PRINT,' LAST WRITE TO IT95. IADR=',IX
                          CALL SWRIT(IT95,IMOINT,LENGTH)
                     END IF
  200           CONTINUE
                IJ=0
C     ICNT=0
                IRIN=0
                IF (N105 .NE. 0) THEN
                     CALL SREW(IT95)
C     PRINT,' REWOUND IT95'
                     CALL SREAD(IT95,IMOINT,LENGTH)
CDB   PRINT,'  READ FIRST MO BLOCK FROM IT95 FOR FINAL SORT'
                END IF
                DO 500 IORB=1,NSOMO(ASM)
                     DO 480 JORB=1,NSOMO(BSM)
                          IJ=IJ+1
                          KL=0
                          IF (IRIN+NSOAB+1 .GT. LENGTH/2) THEN
CDB   CALL RGETSA(IT95,IX)
CDB   PRINT,' IRIN+NSOAB+1=',NSOAB+1+IRIN,'  IADR=',IX
                               CALL SREAD(IT95,IMOINT,LENGTH)
                               IRIN=0
                          END IF
                          DO 400 KORB=1,IORB
                               LMX=NSOMO(BSM)
                               IF (IORB .EQ. KORB) LMX=JORB
                               DO 390 LORB=1,LMX
                                    II=BFDRT(IORB+NAMO)
                                    JJ=BFDRT(JORB+NBMO)
                                    K=BFDRT(KORB+NAMO)
                                    L=BFDRT(LORB+NBMO)
                                    KL=KL+1
                                    IJKL=IJKL+1
                                    IRIN=IRIN+1
                                    VAL=MOINT(IRIN)
                                    IF (ABS(VAL) .LT. ACRCY) GOTO 385
C  PUT CANNONICAL MO INT INTO CANNONICAL DRT INT
                                    ICNT=ICNT+1
                                    IF (II .LT. JJ) THEN
                                         ITMP=II
                                         II=JJ
                                         JJ=ITMP
                                    END IF
                                    IJ1=IOFF(II)+JJ
                                    IF (K .LT. L) THEN
                                         ITMP=K
                                         K=L
                                         L=ITMP
                                    END IF
                                    KL1=IOFF(K)+L
                                    IF (IJ1 .LT. KL1) THEN
                                         ITMP=II
                                         II=K
                                         K=ITMP
                                         ITMP=JJ
                                         JJ=L
                                         L=ITMP
                                    END IF
                                    IF(JJ-K) 310,330,360
 310                                JA=K
                                    KA=JJ
                                    IA=II
                                    LA=L
                                    IF(KA.GE.LA) GO TO 350
                                    LA=KA
                                    KA=L
                                    IF(JA.GT.KA) GO TO 370
                                    MB=2
                                    GO TO 380
 330                                IA=II
                                    JA=JJ
                                    KA=K
                                    LA=L
                                    IF(JA.GT.LA) GO TO 350
                                    IF(IA.EQ.JA) GO TO 350
                                    JA=IA
                                    KA=IA
                                    MB=2
                                    GO TO 380
 360                                IA=II
                                    JA=JJ
                                    KA=K
                                    LA=L
                                    MB=2
                                    GO TO 380
 370                                MB=3
                                    GO TO 380
 350                                MB=1
 380                                IJA=IOFF(IA)+JA
                                    JSM=ORBSYM(IA)
                                    LKUP=MULTS(JSM)+ORBSYM(JA)
                                    KSM=LKUPSM(LKUP)
                                    LKUP=MULTS(KSM)+ORBSYM(KA)
                                    LSM=LKUPSM(LKUP)
                                    IF(LSM.NE.ORBSYM(LA)) THEN
                                         WRITE(IOUT6,701)
  701                                    FORMAT (//,' ***** SYMMETRY ERR
     *OR IN TWO-',                            'ELECTRON INTEGRALS *****'
     *                                   ,//)
                                         CALL MABORT
                                    END IF
                                    KSMPT=KA+(KSM-1)*NORBS
                                    LSMPT=LA+(LSM-1)*NORBS
                                    MOLBL(ICNT)=IJADD(IJA)+KADD(KSMPT)+
     *                              LADD(LSMPT)+MB
                                    MOBKT(ICNT)=IJGRP(IJA)
                                    MOSRT(ICNT)=VAL
                                    IF (ICNT .GE. LENMO) THEN
                                         MOARR(1)=1
                                         MOARR(2)=ICNT
CDB   CALL RGETSA(IT91,IX)
CDB   PRINT,'WRITING TO 101:  NUMINT=',MOARR(2),'  IADR=',IX
C     CALL BUBY(MOLBL,MOBKT,MOSRT,ICNT)
                                         CALL SWRIT(IT91,MOARR,LENMO*
     *                                   4+2)
                                         SWOUT=SWOUT+LENOUT
                                         DO 382 IOP=1,LENMO*4+2
                                              MOARR(IOP)=0
  382                                    CONTINUE
                                         ICNT=0
                                    END IF
C     WRITE(ICHECK,'(G18.8,I5)') VAL,LAD
 385                                CONTINUE
 390                           CONTINUE
 400                      CONTINUE
 480                 CONTINUE
 500            CONTINUE
C     IADR=1
C     MOARR(1)=IADR
C     MOARR(2)=ICNT
C     PRINT,'  TO IT91'
C     CALL RGETSA(IT91,IX)
C     PRINT,'  NUMINT=',MOARR(2), '  IADR=',IX
C     PRINT,'  BUBY FOR FINAL LIST:'
C     CALL BUBY(MOLBL,MOBKT,MOSRT,ICNT)
C     CALL SWRIT(IT91,MOARR,LENMO*4+2)
C     PRINT,'  BWROTE FINAL LIST TO IT91'
C     ICNT=0
C     FOR IOP=1,LENMO*4+2
C     MOARR(IOP)=0
C     END FOR
C     PRINT,'  ZEROED MOARR'
 88             CONTINUE
  800      CONTINUE
           AB=AB+1
 99        CONTINUE
  900 CONTINUE
      IF (NSYM0 .LT. 4) THEN
C     PRINT,'  NSYM0=',NSYM0
           MOARR(1)=-1
           MOARR(2)=ICNT
C     CALL RGETSA(IT91,IX)
C     PRINT, '  LAST BUFFER AT ADR ',IX
           CALL SWRIT(IT91,MOARR,LENMO*4+2)
           SWOUT=SWOUT+LENOUT
      END IF
C
      WRITE(IOUT6,702) SRIN,IT93,SWOUT,IT91
  702 FORMAT (5X,I8,' SECTORS  READ FROM TAPE ',I3,/
     #,       5X,I8,' SECTORS WRITTEN TO TAPE ',I3)
C
      RETURN
      END
