      SUBROUTINE TRNSSS(A,B,C,U,IOFF,NSOFF,TTT,NSYMHF,NSO,CMAX,BMAX,
     1AMAX,MOINT,NSOMO,MOOFF,LENGTH,MOARR,MOBKT,MOLBL,MOSRT,BFDRT,
     2ORBSYM,IJGRP,IJADD,KADD,LADD,LENMO,ACRCY,LCIBUF,ITTT,IMOINT)
C
C  TRANSFORM INTEGRALS OF TYPE (AB,CD)
C
C
C  MOARR : MOBKT : MOLBL : MOSRT ARE IMPLICITLY EQUIVALENCED
C
      IMPLICIT REAL*8 (A-H,O,T-Z)
      IMPLICIT INTEGER (I-N,P,Q,R,S)
C
C
      integer lcibuf,lenmo,norbs,nsymhf,numij,lenr
      INTEGER TIME,ENDD,STARTT
      INTEGER H,CMAX,BMAX,AMAX,T,MOARR(LCIBUF),MOBKT(LENMO),MOLBL(LENMO)
      INTEGER BFDRT(NORBS),NSOMO(NSYMHF),MOOFF(NSYMHF),SYMORB
      INTEGER ORBSYM(NORBS),IJGRP(NUMIJ),IJADD(NUMIJ),MINI
      INTEGER KADD(SYMORB),LADD(SYMORB),MULTS(8),LKUPSM(64)
      INTEGER ID(42),IC(42),IB(42),IAA(42),ITTT(1),IMOINT(1)
      REAL*8 SEC
      REAL*8 MOINT(LENR),MOSRT(LENMO),ACRCY
C
      COMMON /D2H/    MAXD2H(42)
      COMMON /BDIM/   IMAX,JMAX,KMAX,LENR,LENAOR,IBIG,IBIGR,LENPWS
      COMMON /DIMS/   SYMORB,NUMIJ,NIJ,NBF,NORBS,NMAX,NBFORB,NSYM0
      COMMON /TAPES/  IT34,IT58,IT47,IT52,IT91,IT93,IT95
      COMMON /OUT/    IOUT6,ICHECK,IDBG,ICNT
C
      DIMENSION A(AMAX),B(BMAX),C(CMAX),IOFF(NIJ),NSOFF(NSYMHF)
      DIMENSION U(NBFORB),TTT(LENR),NSO(NSYMHF)
C
      DATA MULTS/0,8,16,24,32,40,48,56/
      DATA LKUPSM/1,2,3,4,5,6,7,8,2,1,4,3,6,5,8,7,3,4,1,2,7,8,5,6,4,3,
     A2,1,8,7,6,5,5,6,7,8,1,2,3,4,6,5,8,7,2,1,4,3,7,8,5,6,3,4,1,2,
     B8,7,6,5,4,3,2,1/
      DATA IAA /2,1,1,2,1,4,3,1,3,3,4,4,3,1,2,2,1,1,2,2,1,4,3,2,3,4,1,
     12,1,4,1,2,3,2,1,6,1,2,5,1,3,5/
      DATA IB /3,3,2,5,5,5,5,2,4,5,6,5,6,5,6,5,6,3,4,3,4,5,6,7,5,6,7,
     15,6,7,5,6,7,3,4,7,3,4,7,2,4,6/
      DATA IC /1,2,3,1,2,3,4,5,5,1,1,2,2,3,3,4,4,5,5,6,6,1,1,1,2,2,2,
     13,3,3,4,4,4,5,5,5,6,6,6,7,7,7/
      DATA ID /4,4,4,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,
     18,8,8,8,8,8,8,8,8,8,8,8,8,8,8/
C
C
C
      WRITE(IOUT6,791)
  791 FORMAT (/,' TRANSFORMING INTEGRALS OF SYMMETRIES [AB;CD]:')
C
CDB   IDBG32=0
CDB   IF((IDBG.AND.32).NE.0) IDBG32=1
      LENOUT=I2SEC(LENMO*4+2)
      SWOUT=0
      SRIN=0
      ISTOP=3
      IF (NSYMHF .EQ. 8) ISTOP=42
      DO 900 JJJ=1,ISTOP
           ND=NSO(ID(JJJ))
           NDMO=NSOMO(ID(JJJ))
           NDO=NSOFF(ID(JJJ))
           NDOMO=MOOFF(ID(JJJ))
           N105=0
           NRED=1
CDB   PRINT,' JJJ=',JJJ
C     WRITE(IOUT6,'('  ISM=',I4,' JSM=',I4,' KSM=',I4,' LSM=',I4)')
C    1 ID(JJJ),IC(JJJ),IB(JJJ),IAA(JJJ)
           NA=NSO(IAA(JJJ))
           NAMO=NSOMO(IAA(JJJ))
           NAO=NSOFF(IAA(JJJ))
           NAOMO=MOOFF(IAA(JJJ))
           NB=NSO(IB(JJJ))
           NBMO=NSOMO(IB(JJJ))
           NBO=NSOFF(IB(JJJ))
           NBOMO=MOOFF(IB(JJJ))
           NC=NSO(IC(JJJ))
           NCMO=NSOMO(IC(JJJ))
           NCO=NSOFF(IC(JJJ))
           NCOMO=MOOFF(IC(JJJ))
           NABC=NA*NB*NC
           NAB=NA*NB
           NABCMO=NAMO*NBMO*NCMO
           NABMO=NAMO*NBMO
           MAXSZ=MAXD2H(JJJ)*2
C     PRINT,'  MAXSZ=',MAXSZ
           IF (MAXSZ .EQ. 0) GOTO 99
           NWRT=MAXSZ/LENGTH
           IF (LENGTH*NWRT .NE. MAXSZ) NWRT=NWRT+1
CDB   PRINT,'  NWRT=',NWRT,'  MAXSZ=',MAXSZ
           IPNT=0
           CALL SREW(IT95)
           DO 1 II=1,LENGTH/2
                MOINT(II)=0.0D+00
    1      CONTINUE
           IWRD=LENGTH
           IF (NWRT .EQ. 1) IWRD=MAXSZ
CDB   CALL RGETSA(IT93,IX)
           CALL SREAD(IT93,ITTT,IWRD)
           SRIN=SRIN+I2SEC(IWRD)
CDB   PRINT,'  READ IT93,IADR=',IX,'  WORDS=',IWRD
C     WRITE(IOUT6,'(7G18.8)') (TTT(IOP), IOP=1,MAXSZ/2)
           IF (NAMO*NBMO*NCMO*NDMO .EQ. 0) GOTO 99
           DO 200 P=1,ND
                PND=(P+NDO-1)*NORBS+NDOMO
C     PX=(P-1)*NABC-NAB-NA+ABCD(JJJ)
                DO 2 II=1,CMAX
                     C(II)=0.0D+00
    2           CONTINUE
                DO 13 Q=1,NC
                     QNC=(Q+NCO-1)*NORBS+NCOMO
                     PX=PX+NAB
                     DO 3 II=1,BMAX
                          B(II)=0.0D+00
    3                CONTINUE
                     PQ=PX
                     DO 9 R=1,NB
                          RNB=(R+NBO-1)*NORBS+NBOMO
                          PQ=PQ+NA
                          DO 4 II=1,AMAX
                               A(II)=0.0D+00
    4                     CONTINUE
                          PQRS=PQ
                          DO 6 S=1,NA
                               SNA=(S+NAO-1)*NORBS+NAOMO
                               PQRS=PQRS+1
                               IPNT=IPNT+1
                               F=TTT(IPNT)
C     WRITE(ICHECK,'(G18.8)') TTT(IPNT)
                               IF (F .EQ. 0.0D+00) GOTO 15
                               DO 5 L=1,NAMO
                                    A(L)=A(L)+U(SNA+L)*F
    5                          CONTINUE
 15                            CONTINUE
                               IF ((IPNT .GE. LENGTH/2) .AND. (NRED .LT.
     *                          NWRT))THEN
                                    NRED=NRED+1
C     PRINT,'  NRED=',NRED
                                    IWRD=LENGTH
C     IF (NRED*LENGTH .GT. MAXSZ) IWRD=MAXSZ-(NRED-1)*LENGTH
C     PRINT,'  IWRD=',IWRD
                                    IPNT=0
CDB   CALL RGETSA(IT93,IX)
                                    CALL SREAD(IT93,ITTT,IWRD)
                                    SRIN=SRIN+I2SEC(IWRD)
CDB   PRINT,' READ IT93, IADR=',IX,'  WORDS=',IWRD
                               END IF
    6                     CONTINUE
                          KL=0
                          DO 8 K=1,NBMO
                               URNBK=U(RNB+K)
                               IF (URNBK .EQ. 0.0D+00) THEN
                                    KL=KL+NAMO
                                    GOTO 20
                               END IF
                               DO 7 L=1,NAMO
                                    KL=KL+1
                                    B(KL)=B(KL)+URNBK*A(L)
    7                          CONTINUE
 20                            CONTINUE
    8                     CONTINUE
    9                CONTINUE
                     DO 12 J=1,NCMO
                          UQNCJ=U(QNC+J)
                          IF (UQNCJ .EQ. 0.0D+00) GOTO 25
                          JKL=(J-1)*NABMO
                          KL=0
                          DO 11 K=1,NBMO
                               DO 10 JUNKN=1,NAMO
                                    KL=KL+1
                                    JKL=JKL+1
                                    C(JKL)=C(JKL)+UQNCJ*B(KL)
   10                          CONTINUE
   11                     CONTINUE
 25                       CONTINUE
   12                CONTINUE
   13           CONTINUE
                IPOUT=0
                IF (N105 .NE. 0) THEN
                     CALL SREW(IT95)
                     CALL RGETSA(IT95,IADR)
                     CALL SREAD(IT95,IMOINT,LENGTH)
CDB   PRINT, ' READ IT95 AT ADR',IADR
                     CALL  RSETSA(IT95,IADR)
                END IF
                DO 100 I=1,NDMO
                     UPNDI=U(PND+I)
C     IF (UPNDI .EQ. 0.0D+00) THEN
C     IPOUT=IPOUT+NAMO*NBMO*NCMO
C     GOTO 30
C     END IF
                     IJG=IOFF(I+NDOMO)+NCOMO
                     IJ=(I-1)*NABCMO
                     DO 90 J=1,NCMO
                          IJG=IJG+1
                          IJKG=IOFF(IJG)+NAOMO
                          JK=(J-1)*NABMO
                          IF (IPOUT+NABMO+1 .GT. LENGTH/2) THEN
                               N105=1
CDB   CALL RGETSA(IT95,IX)
CDB   PRINT,'  IPOUT+NABMO+1=',IPOUT+NABMO+1,' IADR=',IX
                               CALL SWRIT(IT95,IMOINT,LENGTH)
                               IPOUT=0
                               IF (P .NE. 1) THEN
                                    CALL RGETSA(IT95,IADR)
                                    CALL SREAD(IT95,IMOINT,LENGTH)
                                    CALL  RSETSA(IT95,IADR)
                               ELSE
                                    DO 60 II=1,LENGTH/2
                                         MOINT(II)=0.0D+00
   60                               CONTINUE
                               END IF
                          END IF
                          IF (UPNDI .EQ. 0.0D+00) THEN
                               IPOUT=IPOUT+NAMO*NBMO
                               GOTO 30
                          END IF
                          DO 80 K=1,NBMO
                               IJKLG=IJKG+IOFF(K+NBOMO)
                               JKL=JK+(K-1)*NAMO
                               DO 70 JUNKN=1,NAMO
                                    IJKLG=IJKLG+1
                                    JKL=JKL+1
                                    IJKL=IJ+JKL
                                    IPOUT=IPOUT+1
                                    MOINT(IPOUT)=MOINT(IPOUT)+UPNDI*C(
     *                              JKL)
   70                          CONTINUE
   80                     CONTINUE
 30                       CONTINUE
   90                CONTINUE
  100           CONTINUE
                IF ((IPOUT .NE. 0) .AND. (N105 .NE. 0))THEN
CDB   CALL RGETSA(IT95,IX)
CDB   PRINT,' LAST WRITE TO IT95. IADR=',IX
                     CALL SWRIT(IT95,IMOINT,LENGTH)
                END IF
  200      CONTINUE
C     ICNT=0
           IRIN=0
           IF (N105 .NE. 0) THEN
                CALL SREW(IT95)
                CALL SREAD(IT95,IMOINT,LENGTH)
CDB   PRINT,'  READ IT95 FOR RESRT'
           END IF
           DO 500 IORB=1,NDMO
                DO 480 JORB=1,NCMO
                     IF (IRIN+NABMO+1 .GT. LENGTH/2) THEN
CDB   CALL RGETSA(IT95,IX)
CDB   PRINT,'  IRIN+NABMO+1=',IRIN+NABMO+1,'  IADR=',IX
                          CALL SREAD(IT95,IMOINT,LENGTH)
                          IRIN=0
                     END IF
                     DO 400 KORB=1,NBMO
                          DO 390 LORB=1,NAMO
                               II=BFDRT(IORB+NDOMO)
                               JJ=BFDRT(JORB+NCOMO)
                               K=BFDRT(KORB+NBOMO)
                               L=BFDRT(LORB+NAOMO)
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
                               IF (ICNT .GE. LENMO) THEN
                                    MOARR(1)=1
                                    MOARR(2)=ICNT
C     CALL BUBY(MOLBL,MOBKT,MOSRT,ICNT)
CDB   CALL RGETSA(IT91,IX)
CDB   PRINT,'  TO IT91,  NUMINT=',MOARR(2),'  IADR=',IX
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
C     IF (ICNT .GT. 0) THEN
C     IADR=1
C     IF (JJJ .EQ. 3) IADR=0
C     MOARR(1)=IADR
C     MOARR(2)=ICNT
C     CALL RGETSA(IT91,IX)
C     PRINT,'  TO  IT91,  NUMINT=',MOARR(2),'  IADR=',IX
C     CALL BUBY(MOLBL,MOBKT,MOSRT,ICNT)
C     CALL SWRIT(IT91,MOARR,LENMO*4+2)
C     ICNT=0
C     END IF
C     FOR IOP=1,LENMO*4+2
C     MOARR(IOP)=0
C     END FOR
 99        CONTINUE
  900 CONTINUE
CDB   CALL RGETSA(IT91,IX)
CDB   PRINT,'  FINAL BUFFER TO IT91 AT ADR',IX,' NUMINT=',ICNT
      MOARR(1)=-1
      MOARR(2)=ICNT
      CALL SWRIT(IT91,MOARR,LENMO*4+2)
      SWOUT=SWOUT+LENOUT
C
      WRITE(IOUT6,702) SRIN,IT93,SWOUT,IT91
  702 FORMAT (5X,I8,' SECTORS  READ FROM TAPE ',I3,/
     #5X,I8,' SECTORS WRITTEN TO TAPE ',I3)
C
      RETURN
      END
