      SUBROUTINE TRNTTT(A,B,C,U,IOFF,NSOFF,TTT,NSYMHF,NSO,CMAX,BMAX,
     1AMAX,MOINT,NSOMO,MOOFF,LENGTH,MOARR,MOBKT,MOLBL,MOSRT,BFDRT,
     2ORBSYM,IJGRP,IJADD,KADD,LADD,LENMO,ACRCY,LCIBUF,ITTT,IMOINT)
C
C  TRANSFORM INTEGRALS OF TYPE (AA,AA)
C
      IMPLICIT REAL*8 (A-H,O,T-Z)
      IMPLICIT INTEGER (I-N,P,Q,R,S)
C
C
      integer lcibuf,lenmo,norbs,nsymhf,numij,lenr
      INTEGER TIME,STARTT,ENDD
      INTEGER H,CMAX,BMAX,AMAX,T,MOARR(LCIBUF)
      INTEGER MOBKT(LENMO),MOLBL(LENMO),BFDRT(NORBS),MINI
      INTEGER NSOMO(NSYMHF),MOOFF(NSYMHF),SYMORB,ITTT(1),IMOINT(1)
      INTEGER ORBSYM(NORBS),IJGRP(NUMIJ),IJADD(NUMIJ)
      INTEGER KADD(SYMORB),LADD(SYMORB),MULTS(8),LKUPSM(64)
      REAL*8 SEC
CDB   REAL*8 SUM
      REAL*8 MOINT(LENR),MOSRT(LENMO),ACRCY
C
      COMMON /DIMS/   SYMORB,NUMIJ,NIJ,NBF,NORBS,NMAX,NBFORB,NSYM0
      COMMON /BDIM/   IMAX2,JMAX2,KMAX2,LENR,LENAOR,IBIG,IBIGR,LENPWS
      COMMON /TAPES/  IT34,IT58,IT47,IT52,IT91,IT93,IT95
      COMMON /OUT/    IOUT6,ICHECK,IDBG,ICNT
      COMMON /TEST/ ITEST
C
      DIMENSION A(AMAX),B(BMAX),C(CMAX),IOFF(NIJ),NSOFF(NSYMHF)
      DIMENSION U(NBFORB),NSO(NSYMHF),TTT(LENR)
C
      DATA MULTS/0,8,16,24,32,40,48,56/
      DATA LKUPSM/1,2,3,4,5,6,7,8,2,1,4,3,6,5,8,7,3,4,1,2,7,8,5,6,4,3,
     A2,1,8,7,6,5,5,6,7,8,1,2,3,4,6,5,8,7,2,1,4,3,7,8,5,6,3,4,1,2,
     B8,7,6,5,4,3,2,1/
C
      WRITE(IOUT6,791)
  791 FORMAT (/,' TRANSFORMING INTEGRALS OF SYMMETRIES [AA;AA]:')
CDB   IDBG32=0
CDB   IF((IDBG.AND.32).NE.0) IDBG32=1
CDB   PRINT,'  LENGTH=',LENGTH,'  LENMO=',LENMO
      LENOUT=I2SEC(LENMO*4+2)
      SWOUT=0
      SRIN=0
      DO 900 ISM=1,NSYMHF
           N105=0
C     PRINT *,' ISM=',ISM
           IF (NSO(ISM).EQ.0) GOTO 99
           CALL SREW(IT95)
           IF (NSO(ISM) .EQ. 0) THEN
                MAXSZ=0
           ELSE
                ITMP=IOFF(NSO(ISM))+NSO(ISM)
                MAXSZ=(IOFF(ITMP)+ITMP)*2
           END IF
           NWRT=MAXSZ/LENGTH
           IF (LENGTH*NWRT .NE. MAXSZ) NWRT=NWRT+1
           NRED=1
C     PRINT *,'  NWRT=',NWRT,'  MAXSZ=',MAXSZ
           DO 1 II=1,LENGTH/2
                MOINT(II)=0.0D+00
    1      CONTINUE
           IWRD=LENGTH
           IF (NWRT .EQ. 1) IWRD=MAXSZ
           CALL RGETSA(IT93,IX)
           CALL SREAD(IT93,ITTT,IWRD)
           SRIN=SRIN+I2SEC(IWRD)
C     PRINT *, '  READ FROM IT93, IADR=',IX,' WORDS=',IWRD
C     ITOPT=IWRD/INTOWP(1)
C     ITOP9=ITOPT-10
C     WRITE(IOUT6,801) (TTT(IOP),IOP=1,10)
C     WRITE(IOUT6,801) (TTT(IOP),IOP=ITOP9,ITOPT)
  801 FORMAT (' READING IN TTT',7G18.8)
           IF (NSOMO(ISM).EQ.0) GOTO 99
           IPNT=0
           PQ=0
           NN=NSOFF(ISM)
           NNMO=MOOFF(ISM)
           IF (NSOMO(ISM) .EQ. 0) THEN
                NNNMO=NSOMO(ISM)
           ELSE
                NNNMO=IOFF(NSOMO(ISM))+NSOMO(ISM)
           END IF
C
CDB   TUM=0.0D+00
CDB   SUM=0.0D+00
C
           DO 200 P=1,NSO(ISM)
                PNN=(P+NN-1)*NORBS+NNMO
                DO 2 II=1,CMAX
                     C(II)=0.0D+00
    2           CONTINUE
                E=1.0D+00
                DO 13 Q=1,P
                     QNN=(Q+NN-1)*NORBS+NNMO
                     PQ=PQ+1
                     DO 3 II=1,BMAX
                          B(II)=0.0D+00
    3                CONTINUE
                     IF (P.EQ.Q) E=0.5D+00
                     RS=0
                     DO 9 R=1,P
                          RNN=(R+NN-1)*NORBS+NNMO
                          DO 4 II=1,AMAX
                               A(II)=0.0D+00
    4                     CONTINUE
                          EE=E
                          T=R
                          IF (R.EQ.P) T=Q
                          DO 6 S=1,T
                               SNN=(S+NN-1)*NORBS+NNMO
                               RS=RS+1
                               IPNT=IPNT+1
C
CDB   TUM=TUM+ABS(TTT(IPNT))
C
                               IF (R.EQ.S) EE=EE*0.5D+00
                               IF (RS.EQ.PQ) EE=EE*0.5D+00
                               F=TTT(IPNT)*EE
                               IF (F .EQ. 0.0D+00) GOTO 15
                               DO 5 L=1,NSOMO(ISM)
                                    A(L)=A(L)+U(SNN+L)*F
    5                          CONTINUE
 15                            CONTINUE
                               IF ((IPNT .GE. LENGTH/2) .AND. (NRED .LT.
     *                          NWRT))THEN
                                    NRED=NRED+1
                                    IWRD=LENGTH
C     IF (NRED*LENGTH .GT. MAXSZ) IWRD=MAXSZ-(NRED-1)*LENGTH
                                    IPNT=0
CDB   CALL RGETSA(IT93,IX)
CDB   PRINT *,' SHOULDNT SEE THIS'
                                    CALL SREAD(IT93,ITTT,IWRD)
                                    SRIN=SRIN+I2SEC(IWRD)
CDB   PRINT,'  READ FROM IT93, IADR=',IX,'  WORDS=',IWRD
C     PRINT,'  NRED=',NRED
                               END IF
    6                     CONTINUE
                          KL=0
                          DO 8 K=1,NSOMO(ISM)
                               URNNK=U(RNN+K)
                               AK=A(K)
                               DO 7 L=1,K
                                    KL=KL+1
                                    B(KL)=B(KL)+URNNK*A(L)+U(RNN+L)*AK
    7                          CONTINUE
    8                     CONTINUE
    9                CONTINUE
                     JKL=0
                     DO 12 J=1,NSOMO(ISM)
                          UQNNJ=U(QNN+J)
                          IF (UQNNJ .EQ. 0.0D+00) THEN
                               JKL=JKL+IOFF(NSOMO(ISM)+1)
                               GOTO 20
                          END IF
                          KL=0
                          DO 11 K=1,NSOMO(ISM)
                               DO 10 L=1,K
                                    JKL=JKL+1
                                    KL=KL+1
                                    C(JKL)=C(JKL)+UQNNJ*B(KL)
   10                          CONTINUE
   11                     CONTINUE
 20                       CONTINUE
   12                CONTINUE
   13           CONTINUE
                IJ=0
                IPOUT=0
                IF (N105.NE.0) THEN
                     CALL SREW(IT95)
                     CALL RGETSA(IT95,IADR)
CDB   PRINT,'  READING IT95 AT ADR ',IADR
C******** TJL
cets  ITEST = 1
                     CALL SREAD(IT95,IMOINT,LENGTH)
                     CALL  RSETSA(IT95,IADR)
                END IF
                DO 100 I=1,NSOMO(ISM)
                     UPNNI=U(PNN+I)
                     IK=(I-1)*NNNMO
                     DO 90 J=1,I
                          UPNNJ=U(PNN+J)
                          IJ=IJ+1
                          IJKL=IOFF(IJ)
                          IKL=IK
                          JKL=(J-1)*NNNMO
                          IF (IPOUT+IJ+1 .GT. LENGTH/2) THEN
                               N105=1
CDB   CALL RGETSA(IT95,IX)
CDB   PRINT *,' WRITING TO IT95'
CDB   PRINT,'  LENGTH/2=',LENGTH/2,'  IPOUT+IJ+1=',IPOUT+IJ+1,' IADR=',
CDB  1               IX
                               CALL SWRIT(IT95,IMOINT,LENGTH)
                               IPOUT=0
                               IF (P .NE. 1) THEN
                                    CALL RGETSA(IT95,IADR)
CDB   PRINT,'  READING FROM IT95 AT ADR ',IADR
                                    CALL SREAD(IT95,IMOINT,LENGTH)
                                    CALL  RSETSA(IT95,IADR)
                               ELSE
                                    DO 60 II=1,LENGTH/2
                                         MOINT(II)=0.0D+00
   60                               CONTINUE
                               END IF
                          END IF
                          DO 80 K=1,I
                               UPNNK=U(PNN+K)
                               KIJ=(K-1)*NNNMO+IJ
                               H=K
                               IF (K.EQ.I) H=J
                               DO 70 L=1,H
                                    IPOUT=IPOUT+1
                                    IJKL=IJKL+1
                                    IKL=IKL+1
                                    JKL=JKL+1
                                    LIJ=(L-1)*NNNMO+IJ
                                    MOINT(IPOUT)=MOINT(IPOUT)+UPNNI*C(
     *                              JKL)+UPNNJ*C(IKL)+UPNNK*C(LIJ)+U(PNN
     *                              +L)*C(KIJ)
CDB   SUM=SUM+ABS(MOINT(IPOUT))
   70                          CONTINUE
   80                     CONTINUE
   90                CONTINUE
  100           CONTINUE
                IF ((IPOUT .NE. 0) .AND. (N105 .NE. 0))THEN
CDB   CALL RGETSA(IT95,IX)
CDB   WRITE(IOUT6,'('  WRITING',I8,' WORDS TO IT95 AT SEC',I8)')LENGTH,I
                     CALL SWRIT(IT95,IMOINT,LENGTH)
                END IF
  200      CONTINUE
C
CDB   WRITE(IOUT6,'(' --IN TTT, MOINT SUM IS',G20.10)') SUM
CDB   WRITE(IOUT6,'('           TTT   SUM IS',G20.10)') TUM
C
           IJ=0
C     ICNT=0
           IRIN=0
           IF (N105 .NE. 0) THEN
                CALL SREW(IT95)
CDB   CALL RGETSA(IT95,IX)
CDB   PRINT,'  READING IT95 FOR RESORT. IADR=',IX
                CALL SREAD(IT95,IMOINT,LENGTH)
           END IF
C     PRINT *,'  FINISHED TRANSFORM. DOING RE-SORT. '
C     PRINT *,'  LENMO=',MAXSZ,'  MAXSZ=',MAXSZ
           DO 500 IORB=1,NSOMO(ISM)
                DO 480 JORB=1,IORB
                     IJ=IJ+1
                     KL=0
                     IF (IRIN+IJ+1 .GT. LENGTH/2) THEN
CDB   CALL RGETSA(IT95,IX)
CDB   PRINT,'  IRIN+IJ+1=',IRIN+IJ+1, ' IADR=',IX
                          CALL SREAD(IT95,IMOINT,LENGTH)
                          IRIN=0
                     END IF
                     DO 400 KORB=1,IORB
                          LMX=KORB
                          IF (KORB .EQ. IORB) LMX=JORB
                          DO 390 LORB=1,LMX
                               II=BFDRT(IORB+NNMO)
                               JJ=BFDRT(JORB+NNMO)
                               K=BFDRT(KORB+NNMO)
                               L=BFDRT(LORB+NNMO)
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
 310                           JA=K
                               KA=JJ
                               IA=II
                               LA=L
                               IF(KA.GE.LA) GO TO 350
                               LA=KA
                               KA=L
                               IF(JA.GT.KA) GO TO 370
                               MB=2
                               GO TO 380
 330                           IA=II
                               JA=JJ
                               KA=K
                               LA=L
                               IF(JA.GT.LA) GO TO 350
                               IF(IA.EQ.JA) GO TO 350
                               JA=IA
                               KA=IA
                               MB=2
                               GO TO 380
 360                           IA=II
                               JA=JJ
                               KA=K
                               LA=L
                               MB=2
                               GO TO 380
 370                           MB=3
                               GO TO 380
 350                           MB=1
 380                           IJA=IOFF(IA)+JA
                               JSM=ORBSYM(IA)
                               LKUP=MULTS(JSM)+ORBSYM(JA)
                               KSM=LKUPSM(LKUP)
                               LKUP=MULTS(KSM)+ORBSYM(KA)
                               LSM=LKUPSM(LKUP)
                               IF(LSM.NE.ORBSYM(LA)) THEN
                                    WRITE(IOUT6,701)
  701                               FORMAT (//,' ***** SYMMETRY ERROR IN
     * TWO-',                         'ELECTRON INTEGRALS *****',//)
                                    CALL MABORT
                               END IF
                               KSMPT=KA+(KSM-1)*NORBS
                               LSMPT=LA+(LSM-1)*NORBS
                               MOLBL(ICNT)=IJADD(IJA)+KADD(KSMPT)+LADD(
     *                         LSMPT)+MB
                               MOBKT(ICNT)=IJGRP(IJA)
                               MOSRT(ICNT)=VAL
C     IF ((ICNT.LT.10) .OR. (ICNT.GT.(MAXSZ/2-10))) THEN
C     WRITE(IOUT6,802) ICNT,MOLBL(ICNT),MOBKT(ICNT),MOSRT(ICNT)
C 802 FORMAT (1X,3I5,E15.5)
C     END IF
                               IF (ICNT .GE. LENMO) THEN
                                    MOARR(1)=1
                                    MOARR(2)=ICNT
C     CALL BUBY(MOLBL,MOBKT,MOSRT,ICNT)
      CALL RGETSA(IT91,IX)
C     WRITE(IOUT6,1101) IX,MOARR(2)
C1101 FORMAT('  WRITING TO IT91. IADR=',I5,'  NUM INTS=',I9)
CDB   PRINT,' WRITING TO IT91. IADR=',IX,' NUM INT=',MOARR(2)
                                    CALL SWRIT(IT91,MOARR,LENMO*4+2)
                                    SWOUT=SWOUT+LENOUT
                                    DO 382 IOP=1,LENMO*4+2
                                         MOARR(IOP)=0
  382                               CONTINUE
                                    ICNT=0
                               END IF
C     WRITE(ICHECK,'(G18.8,I5)') VAL,LAD
 385                           CONTINUE
 390                      CONTINUE
 400                 CONTINUE
 480            CONTINUE
 500       CONTINUE
 99        CONTINUE
  900 CONTINUE
      IF (NSYM0 .EQ. 1) THEN
           MOARR(1)=-1
           MOARR(2)=ICNT
      CALL RGETSA(IT91,IX)
C     WRITE(IOUT6,1101)IX,MOARR(2)
C     PRINT,' WRITNG TO IT91. IADR=',IX,' NUM INT=',MOARR(2),' SIZE=',
C    1     LENMO*4+2
C     CALL BUBY(MOLBL,MOBKT,MOSRT,ICNT)
           CALL SWRIT(IT91,MOARR,LENMO*4+2)
           SWOUT=SWOUT+LENOUT
           ICNT=0
           DO 910 IOP=1,LENMO*4+2
                MOARR(IOP)=0
  910      CONTINUE
      END IF
C
      WRITE(IOUT6,702) SRIN,IT93,SWOUT,IT91
  702 FORMAT (5X,I8,' SECTORS  READ FROM TAPE ',I3,/
     #,       5X,I8,' SECTORS WRITTEN TO TAPE ',I3)
C
      RETURN
      END
