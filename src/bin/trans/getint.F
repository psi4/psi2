      SUBROUTINE GETINT(IOFF,NSOFF,OFFA1,NINTS,PKI,SRTI,SRTAR,RESRT,F,
     *P,H,S,U,UMO,LENGTH,NSO,IOUT,IFZ,NBKT,BKTSZ,BKTPNT,START,ILAST,
     *INTOFF,IJBKT,GFMAX,LBLS,INTS,LEN2,NSYMHF,BKTR,ILAST2,SPACE
     #,                 IRESRT,IFF)
C
C
C  THE FOLLOWING DUMMY ARGUMENTS ARE IMPLICITLY EQUIVALENCED
C  REAL/INTEGER PAIRS:
C     PKI : NINTS           SRTAR : SRTI        LBLS : INTS
C
C
      IMPLICIT INTEGER (A-Z)
C
C
      integer nij,lenaor,nbf,norbs,ibigr,lenr,bktr
      INTEGER ID(42),IC(42),IB(42),IAA(42),IRESRT(1),IFF(1)
      REAL*8 DSQRT
      REAL*8 X,ACRCY,FZC,REPNRG
      REAL*8 S(NIJ),H(NIJ),PKI(LENAOR),VAL,U(NBF,NBF),P(NIJ)
      REAL*8 UMO(NORBS,NBF),RESRT(LENR),INTS(BKTR),F(NIJ),SRTAR(IBIGR)
      integer errcod,iacrcy
C
      COMMON /BDIM/   IMAX2,JMAX2,KMAX2,LENR,LENAOR,IBIG,IBIGR,LENPWS
      COMMON /DIMS/   SYMORB,NUMIJ,NIJ,NBF,NORBS,NMAX,NBFORB,NSYM0
      COMMON /OUT/    IOUT6,ICHECK,IDBG,ICNT
      COMMON /COR/    FZC,REPNRG
      COMMON /SRT/    VAL,IADR,IBKT,IADBK,MNBK,MXBK
      COMMON /D2H/    MAXD2H(42)
      COMMON /TAPES/  IT34,IT58,IT47,IT52,IT91,IT93,IT95
C
      DIMENSION NSO(NSYMHF),IOFF(NIJ),NSOFF(NSYMHF),OFFA1(NSYMHF),
     *NINTS(LENPWS),IVI(3),IFZ(NBF),IJBKT(36),ILAST(NBKT),INTOFF(NBKT)
      DIMENSION ABCDBK(615),LBLS(BKTSZ),IOUT(NBF)
      DIMENSION SRTI(IBIG),BKTPNT(NBKT),START(NBKT),ILAST2(NBKT)
C
      DATA IAA /2,1,1,2,1,4,3,1,3,3,4,4,3,1,2,2,1,1,2,2,1,4,3,2,3,4,1,
     12,1,4,1,2,3,2,1,6,1,2,5,1,3,5/
      DATA IB /3,3,2,5,5,5,5,2,4,5,6,5,6,5,6,5,6,3,4,3,4,5,6,7,5,6,7,
     15,6,7,5,6,7,3,4,7,3,4,7,2,4,6/
      DATA IC /1,2,3,1,2,3,4,5,5,1,1,2,2,3,3,4,4,5,5,6,6,1,1,1,2,2,2,
     13,3,3,4,4,4,5,5,5,6,6,6,7,7,7/
      DATA ID /4,4,4,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,
     18,8,8,8,8,8,8,8,8,8,8,8,8,8,8/
C
CDB   IDBG4=0
CDB   IF((IDBG.AND.4).NE.0) IDBG4=1
C     PRINT,'  IN GETINT '
C     PRINT,'  NSYM0=',NSYM0
      SWOUT=0
      SRIN=0
 1000 FORMAT(7G18.8)
 1001 FORMAT(G18.8,2I6)
      CALL SETD2H(ABCDBK,NSO,NSYMHF)
      LENSEC=I2SEC(LENGTH)
      errcod = frdi('INT_ACCURACY',IACRCY)
      ACRCY = 10.0D+00**(-IACRCY)
      if (iacrcy.eq.0) ACRCY=1.0D-12
      IADBK=NBKT
      MNBK=1
      MXBK=NBKT
C
C     -----  FOR 1 ELECTRON INTS -----
C
      CALL STUFIT(IT34,NINTS,PKI,LENGTH,S,NSOFF,IOFF,NSYMHF,NIJ,SRIN)
      IMO=0
      DO 600 IBF=1,NBF
           IF(IOUT(IBF).EQ.0) GOTO 600
           IJ=0
           X=0.0D+00
           DO 540 I=1,NBF
                DO 540 J=1,I
                     IJ=IJ+1
                     VAL=U(I,IBF)*U(J,IBF)*S(IJ)
                     X=X+VAL
                     IF(I.NE.J) X=X+VAL
 540       CONTINUE
           IF(ABS(1.0D+00-X).GT.ACRCY*2.0D+00) WRITE(IOUT6,541) IBF,X
 541       FORMAT(10H *********,' WARNING , SCF ORBITAL',I4,' IS NOT
     *     NORMALIZED. INNER PRODUCT IS',F14.8,10H *********)
           IORB=IOUT(IBF)
           IF(IORB.LT.0) GOTO 560
           IMO=IMO+1
           IFZ(IBF)=IMO
           X=1.0D+00/DSQRT(X)
           DO 550 I=1,NBF
                UMO(IMO,I)=U(I,IBF)*X
 550       CONTINUE
           GO TO 600
 560       X=1.0D+00/X
           IJ=0
           DO 580 I=1,NBF
                DO 580 J=1,I
                     IJ=IJ+1
                     P(IJ)=P(IJ)+U(I,IBF)*U(J,IBF)*X
 580       CONTINUE
 600  CONTINUE
CDB   PRINT,' NSOMO= ',NSOMO
CDB   PRINT,' MOOFF= ',MOOFF
CDB   PRINT,' IFZ= ',IFZ
      CALL STUFIT(IT34,NINTS,PKI,LENGTH,H,NSOFF,IOFF,NSYMHF,NIJ,SRIN)
      CALL STUFIT(IT34,NINTS,PKI,LENGTH,H,NSOFF,IOFF,NSYMHF,NIJ,SRIN)
      DO 610 II=1,NIJ
           F(II)=H(II)
  610 CONTINUE
CDB   IF (IDBG4 .NE. 0) WRITE(IOUT6,'(' H INTS:'/)')
CDB   IF (IDBG4 .NE. 0) WRITE(IOUT6,'(7G18.8)') (H(II),II=1,NIJ)
C     WRITE(IOUT6,9001) (F(II),II=1,NIJ)
C9001 FORMAT (/,(5E15.5))
C  OUTPUT DATA TO TAPE 2
CTJL  IVI(1)=NBF
CTJL  IVI(2)=NORBS
CTJL  IVI(3)=3HINT
      NINC=NBF*NORBS*2
      NBFSQ=NBF*NBF*2
      NIJR=NIJ*2
C1    CALL SWRITW (IT52,IVI,3)
C1    CALL SWRIT(IT52,IOUT,128)
C1    CALL SWRIT(IT52,ISYM,128)
C1    CALL SWRIT(IT52,U,NINC)
C1    CALL SWRIT(IT52,U,NBFSQ)
C1    CALL SWRIT(IT52,S,NIJR)
C1    CALL SWRIT(IT52,S,NIJR)
C  CHANGE FIRST U TO C, SECOND S TO P,(LATER)
CDB   IF (IDBG4.NE.0) THEN
CDB   WRITE(IOUT6,'('    IOR   JOR   KOR   LOR   ISM   JSM   KSM   LSM
CDB  *VAL ')')
CDB   END IF
      CALL RGETSA(IT34,IXXX)
      WRITE(IOUT6,*) 'IT34 STOPPED AT ADR ',IXXX
      NPASS=4
      IPCNT=0
      DO 900 IPASS=1,NPASS
C     PRINT, '  IPASS=',IPASS
           IPCNT=IPCNT+1
           CALL  RSETSA(IT34,IXXX)
C     PRINT,'  IT34 RESET TO ADR',IXXX
           CALL SRTONE(BKTPNT,START,ILAST,BKTSZ,NBKT,SRTAR,SRTI,INTOFF,
     *     IT91)
 10        CONTINUE
           CALL SREAD(IT34,NINTS,LENGTH)
           SRIN=SRIN+LENSEC
           ILSTI=NINTS(1)
           NBUF=NINTS(2)
           DO 200 IINT=1,NBUF
C     LOR=NINTS(2*IINT+2).SHIFT.-2
C     LSM=LOR.SHIFT.-8
C     KOR=LSM.SHIFT.-3
C     KSM=(KOR.SHIFT.-8)+1
C     KOR=KOR.AND.'377
C     LSM=(LSM.AND.'7)+1
C     LOR=LOR.AND.'377
C     III=NINTS(2*IINT+2).AND.'3
C     JSM=NINTS(2*IINT+1).SHIFT.-8
C     IOR=JSM.SHIFT.-3
C     ISM=(IOR.SHIFT.-8)+1
C     IOR=IOR.AND.'377
C     JSM=(JSM.AND.'7)+1
C     JOR=NINTS(2*IINT+1).AND.'377
                LOR=ISHFT(NINTS(2*IINT+2),-2)
                LSM=ISHFT(LOR,-8)
                KOR=ISHFT(LSM,-3)
                KSM=ISHFT(KOR,-8)+1
                KOR=IAND(KOR,255)
                LSM=IAND(LSM,7)+1
                LOR=IAND(LOR,255)
                III=IAND(NINTS(2*IINT+2),3)
                JSM=ISHFT(NINTS(2*IINT+1),-8)
                IOR=ISHFT(JSM,-3)
                ISM=ISHFT(IOR,-8)+1
                IOR=IAND(IOR,255)
                JSM=IAND(JSM,7)+1
                JOR=IAND(NINTS(2*IINT+1),255)
                VAL=PKI(IINT+2236)
C     IF (IPASS.EQ.1)WRITE(IOUT6,'(8I6,G18.8)')IOR,JOR,KOR,LOR,ISM,JSM,
C    1KSM,LSM,VAL
C
C     ----- CONSTRUCT FROZEN CORE OPERATOR -----
C
                IF (IPASS.EQ.1) THEN
                     I=IOR+NSOFF(ISM)
                     J=JOR+NSOFF(JSM)
                     K=KOR+NSOFF(KSM)
                     L=LOR+NSOFF(LSM)
                     IJ=IOFF(I)+J
                     KL=IOFF(K)+L
                     IK=IOFF(I)+K
                     IL=IOFF(I)+L
                     IF(IJ.EQ.KL) GOTO 105
                     IF(J-K) 60,85,55
  55                 CONTINUE
                     IF(I.EQ.J) GO TO 90
                     IF(K.EQ.L) GO TO 70
C
C     ----- FOR INTEGRALS LIKE (43,21) -----
C
                     JL=IOFF(J)+L
                     JK=IOFF(J)+K
                     GO TO 67
  60                 CONTINUE
                     IF(J-L) 65,75,61
  61                 CONTINUE
                     IF(I.EQ.K) GO TO 95
C
C     ----- FOR INTEGRALS LIKE (42,31) -----
C
                     JL=IOFF(J)+L
                     JK=IOFF(K)+J
                     GO TO 67
  65                 CONTINUE
                     IF(K.EQ.L) GO TO 80
C
C     ----- FOR INTEGRALS LIKE (41,32) -----
C
                     JL=IOFF(L)+J
                     JK=IOFF(K)+J
  67                 CONTINUE
                     F(IJ)=F(IJ)+4.0D+00*P(KL)*VAL
                     F(KL)=F(KL)+4.0D+00*P(IJ)*VAL
                     F(IL)=F(IL)-P(JK)*VAL
                     F(IK)=F(IK)-P(JL)*VAL
                     F(JL)=F(JL)-P(IK)*VAL
                     F(JK)=F(JK)-P(IL)*VAL
                     GO TO 130
  70                 CONTINUE
C
C     ----- FOR INTEGRALS LIKE (32,11) -----
C
                     JK=IOFF(J)+K
  72                 CONTINUE
                     F(IJ)=F(IJ)+2.0D+00*P(KL)*VAL
                     F(KL)=F(KL)+4.0D+00*P(IJ)*VAL
                     F(IL)=F(IL)-P(JK)*VAL
                     F(JK)=F(JK)-P(IL)*VAL
                     GO TO 130
  75                 CONTINUE
C
C     ----- FOR INTEGRALS LIKE (31,21) -----
C
                     JK=IOFF(K)+J
                     JL=IOFF(J)+L
                     F(JL)=F(JL)-P(IK)*VAL
                     GO TO 67
  80                 CONTINUE
C
C     ----- FOR INTEGRALS LIKE (31,22) -----
C
                     JK=IOFF(K)+J
                     GO TO 72
  85                 CONTINUE
                     IF(I.EQ.J.OR.K.EQ.L) GO TO 110
C
C     ----- FOR INTEGRALS LIKE (32,21) -----
C
                     JK=IOFF(J)+K
                     JL=IOFF(J)+L
                     F(JK)=F(JK)-P(IL)*VAL
                     GO TO 67
  90                 CONTINUE
                     IF(K.EQ.L) GO TO 100
C
C     ----- FOR INTEGRALS LIKE (33,21) -----
C
                     JK=IOFF(J)+K
                     JL=IJ
                     IJ=KL
                     KL=JL
                     GO TO 72
  95                 CONTINUE
C
C     ----- FOR INTEGRALS LIKE (32,31) -----
C
                     JK=IOFF(K)+J
                     JL=IOFF(J)+L
                     F(IK)=F(IK)-P(JL)*VAL
                     GO TO 67
 100                 CONTINUE
C
C     ----- FOR INTEGRALS LIKE (22,11) -----
C
                     F(IJ)=F(IJ)+2.0D+00*P(KL)*VAL
                     F(KL)=F(KL)+2.0D+00*P(IJ)*VAL
                     F(IL)=F(IL)-P(IL)*VAL
                     GO TO 130
 105                 CONTINUE
                     IF(I.EQ.J) GO TO 120
C
C     ----- FOR INTEGRALS LIKE (21,21) -----
C
                     JL=IOFF(J)+L
                     F(IJ)=F(IJ)+3.0D+00*P(IJ)*VAL
                     F(IK)=F(IK)-P(JL)*VAL
                     F(JL)=F(JL)-P(IK)*VAL
                     GO TO 130
 110                 CONTINUE
C
C     ----- FOR INTEGRALS LIKE (21,11) -----
C
                     JK=IOFF(J)+K
                     F(IL)=F(IL)+P(JK)*VAL
                     F(JK)=F(JK)+2.0D+00*P(IL)*VAL
                     GO TO 130
 120                 CONTINUE
C
C     ----- FOR INTEGRALS LIKE (11,11) -----
C
                     F(IJ)=F(IJ)+P(IJ)*VAL
 130                 CONTINUE
                END IF
C
C     ----- END CONSTRUCTION OF THE FROZEN CORE OPERATOR -----
C
                IF (IPASS .EQ. 1) THEN
C
C     ----- SORT INTS OF TYPE (AA,AA) -----
C
                     IF ((ISM.EQ.JSM) .AND. (ISM.EQ.KSM)) THEN
                          IJ=IOFF(IOR)+JOR
                          KL=IOFF(KOR)+LOR
                          IADR=IOFF(IJ)+KL
                          IBKT=ISM
C     WRITE(IOUT6,1001)VAL,IADR,IBKT
                          CALL SRTOUT(BKTPNT,START,ILAST,BKTSZ,NBKT,
     *                    SRTAR,SRTI,INTOFF,IT91)
                     END IF
                ELSE IF (IPASS .EQ. 2) THEN
C
C     ----- SORT INTS OF TYPE (AA,BB) -----
C
                     IF ((ISM.EQ.JSM) .AND. (KSM.EQ.LSM) .AND. (ISM.NE.
     *               KSM))THEN
                          IJ=IOFF(IOR)+JOR
                          KL=IOFF(KOR)+LOR
                          IADR=(IJ-1)*(IOFF(NSO(KSM))+NSO(KSM))+KL
                          IBKT=IJBKT(IOFF(ISM)+KSM)
C     WRITE(IOUT6,'('  ISM=',I4,' JSM=',I4,' KSM=',I4,' LSM=',I4,
C    1' IBKT=',I5,' IADR=',I9)')ISM,JSM,KSM,LSM,IBKT,IADR
C     WRITE(IOUT6,1001)VAL,IADR,IBKT
                          IF (IBKT .EQ. 0) THEN
                               WRITE(IOUT6,701)
  701                          FORMAT (//,' ***** IN GETINT, IBKT=0 ')
                               CALL MABORT
                          END IF
                          CALL SRTOUT(BKTPNT,START,ILAST,BKTSZ,NBKT,
     *                    SRTAR,SRTI,INTOFF,IT91)
                     END IF
                ELSE IF (IPASS .EQ. 3) THEN
C
C     ----- SORT INTS OF TYPE (AB,AB) -----
C
                     IF ((ISM.EQ.KSM) .AND. (JSM.EQ.LSM) .AND. (ISM.NE.
     *               JSM))THEN
C     IK=(IOFF(IOR)+KOR-1)*NSO(JSM)*NSO(JSM)
C     JL=(JOR-1)*NSO(JSM)+LOR
C     IADR=IK+JL
                          IADR=IOFF(IOR)*NSO(JSM)*NSO(JSM)+(JOR-1)*IOR*
     *                    NSO(JSM)+(KOR-1)*NSO(JSM)+LOR
                          IBKT=IJBKT(IOFF(ISM)+JSM)
C     WRITE(IOUT6,1001)VAL,IADR,IBKT
C     WRITE(ICHECK,'('  VAL=',G18.8,'  IK=',I5,'  JL=',I5,' IADR=',
C    1 I7,' IBKT=',I5)')VAL,IK,JL,IADR,IBKT
                          IF (IBKT.EQ.0) THEN
                               WRITE(IOUT6,701)
                               CALL MABORT
                          END IF
                          CALL SRTOUT(BKTPNT,START,ILAST,BKTSZ,NBKT,
     *                    SRTAR,SRTI,INTOFF,IT91)
                     END IF
                ELSE IF (IPASS .EQ.4) THEN
C
C     ----- SORT INT OF TYPE (AB,CD) -----
C
                     IF ((ISM.NE.JSM) .AND. (ISM.NE.KSM)) THEN
                          I=(IOR-1)*NSO(JSM)*NSO(KSM)*NSO(LSM)
                          J=(JOR-1)*NSO(KSM)*NSO(LSM)
                          K=(KOR-1)*NSO(LSM)
                          IADR=I+J+K+LOR
                          IAB=IOFF(ISM)+JSM
                          ISYM=IOFF(IAB)+IOFF(KSM)+LSM
                          IBKT=ABCDBK(ISYM)
C     WRITE(IOUT6,1001)VAL,IADR,IBKT
C     WRITE(IOUT6,'(9I6,G18.8)')IOR,JOR,KOR,LOR,ISM,JSM,KSM,
C    1 LSM,IBKT,VAL
                          IF ((IBKT .LT. 1) .OR. (IBKT .GT. 42)) THEN
                               WRITE(IOUT6,702) ISM,JSM,KSM,LSM,IBKT
  702                          FORMAT (' IN GETINT, ISM=',I4,' JSM=',I4,
     *' KSM=',                      I4,' LSM=',I4,' IBKT=',I8,//)
                               CALL MABORT
                          END IF
                          CALL SRTOUT(BKTPNT,START,ILAST,BKTSZ,NBKT,
     *                    SRTAR,SRTI,INTOFF,IT91)
                     END IF
                END IF
  200      CONTINUE
           IF (ILSTI.EQ.0) GOTO 10
           CALL SRTLST(BKTPNT,START,ILAST,BKTSZ,NBKT,SRTAR,SRTI,INTOFF,
     *     IT91)
C     PRINT *,'  COMPLETED SORT, READYING FOR RESORT'
C     WRITE(ICHECK,'('  READYING FOR RESORT')')
           IF (IPASS.EQ.1) THEN
                IJ=0
                FZC=0.0D+00
                DO 220 I=1,NBF
                     DO 210 J=1,I
                          IJ=IJ+1
                          VAL=P(IJ)*(F(IJ)+H(IJ))
                          FZC=FZC+VAL
                          IF (I.NE.J) FZC=FZC+VAL
  210                CONTINUE
  220           CONTINUE
C     WRITE(IOUT6,9001) (F(II),II=1,NIJ)
C     WRITE(IOUT6,1000) (F(IOP),IOP=1,NIJ)
CDB   CALL RGETSA(IT93,IX)
                CALL SWRIT(IT93,IFF,INTOWP(NIJ))
CDB   PRINT,'  ONE E INTS TO IT93 AT IADR=',IX, '  WORDS=',NIJ*2
                SWOUT=SWOUT+I2SEC(INTOWP(NIJ))
           END IF
 301       FORMAT(/2X,'INTEGRALS OF SYMMETRY TYPE: ',I1,1X,I1,1H,,
     *      I1,1X,I1,'  REQUIRE AN ADDITONAL SORT')
           INTF=INTOFF(1)
           DO 310 III=1,LENR
                RESRT(III)=0.0D+00
  310      CONTINUE
C     PRINT,'  INTF=',INTF,'  ILAST=', ILAST
           IF (IPASS.EQ.1) THEN
                DO 315 II=1,NSYMHF
C     PRINT *,' ISM=',II
C     WRITE(ICHECK,'('  SYMM=',I4)')II
                     IADR=ILAST(II)
C     PRINT,' IADR=',IADR
                     IF (IADR .EQ. -1) GOTO 20
                     DO 314 JJ=1,LENR
                          RESRT(JJ)=0.0D+00
  314                CONTINUE
C     PRINT,'  ZEROED RESRT'
                     MAXSZ=(IOFF(IOFF(NSO(II)+1))+IOFF(NSO(II)+1))*2
                     NWRT=MAXSZ/LEN2
                     IF (NWRT*LEN2 .NE. MAXSZ) NWRT=NWRT+1
C     PRINT *,'  NWRT=',NWRT
C     PRINT *,'  MAXSZ=',MAXSZ
C     PRINT *,'  LEN2=',LEN2
                     IF (NWRT .EQ. 1) THEN
 30                       CALL RREAD(IT91,LBLS,BKTSZ,IADR)
                          IADR=LBLS(1)
C     PRINT *,' IADR=',IADR,'  NUM INTS=',LBLS(2)
                          DO 313 JJ=1,LBLS(2)
C     IF (LBLS(JJ+2) .GT. GFMAX) THEN
C     PRINT,'  LBLS(JJ+2)=',LBLS(JJ+2),'  VAL=',INTS(JJ+INTF)
C     CALL MABORT
C     END IF
                               RESRT(LBLS(JJ+2))=INTS(JJ+INTF)
  313                     CONTINUE
                          IF (IADR .NE. -1) GOTO 30
C     WRITE(IOUT6,1100)
 1100 FORMAT('  THE RESRT ARRAY')
C     WRITE(IOUT6,'(7G18.8)') (RESRT(IOP), IOP=1,10)
C     MTOP=MAXSZ/2
C     M1=MTOP-10
C     WRITE(IOUT6,'(7G18.8)') (RESRT(IOP), IOP=M1,MTOP)
                          IWRD=MAXSZ
                          IST=1
CDB   CALL RGETSA(IT93,IX)
                          CALL SWRIT(IT93,IRESRT,MAXSZ)
                          SWOUT=SWOUT+I2SEC(MAXSZ)
CDB   PRINT,'  WRITE TO IT93, IADR=',IX,'  WORDS=',MAXSZ
                     ELSE
                          WRITE(IOUT6,301)II,II,II,II
C     WRITE(ICHECK,301)II,II,II,II
C     NBKT2=MAXSZ/LEN2
C     IF (LEN2*NBKT2 .NE. MAXSZ) NBKT2=NBKT2+1
                          NBKT2=NWRT
                          CALL SORTWO(BKTPNT,START,ILAST2,BKTSZ,SRTAR,
     *                    INTOFF,LEN2,NWRT,MAXSZ,SPACE,RESRT,IADR,INTF,
     *                    NBKT,LBLS,INTS,LEN2/2,BKTR,NBKT2,SWOUT,SRTI)
                     END IF
 20                  CONTINUE
  315           CONTINUE
                IF (NSYM0 .EQ. 1) GOTO 250
           ELSE IF (IPASS .EQ. 2) THEN
                II=0
                DO 500 ISM=1,NSYMHF
                     DO 501 JSM=1,ISM
                          IF (ISM .EQ. JSM) GOTO 21
                          II=II+1
                          IADR=ILAST(II)
                          IF (IADR .EQ. -1) GOTO 21
C     PRINT,'  ASM=',ISM,' BSM=',JSM,'  IBKT=',II
C     WRITE(ICHECK,'('  ASM=',I4,'  BSM=',I4,'  IBKT=',I4)')ISM,JSM,II
                          DO 502 JJ=1,LENR
                               RESRT(JJ)=0.0D+00
  502                     CONTINUE
                          MAXSZ=(IOFF(NSO(ISM))+NSO(ISM))*(IOFF(NSO(JSM)
     *                    )+NSO(JSM))*2
C     WRITE(IOUT6,1000) (RESRT(IOP),IOP=1,MAXSZ/2)
                          NWRT=MAXSZ/LEN2
                          IF (NWRT*LEN2 .NE. MAXSZ) NWRT=NWRT+1
C     PRINT,'  NWRT=',NWRT,'  MAXSZ=',MAXSZ
                          IF (NWRT .EQ. 1) THEN
 31                            CALL RREAD(IT91,LBLS,BKTSZ,IADR)
                               IADR=LBLS(1)
C     PRINT,'  NEXT ADR=',IADR
                               DO 503 JJ=1,LBLS(2)
                                    RESRT(LBLS(JJ+2))=INTS(JJ+INTF)
  503                          CONTINUE
                               IF (IADR .NE. -1) GOTO 31
                               IWRD=MAXSZ
                               IST=1
CDB   CALL RGETSA(IT93,IX)
                               CALL SWRIT(IT93,IRESRT,MAXSZ)
                               SWOUT=SWOUT+I2SEC(MAXSZ)
CDB   PRINT,'  WRITE TO IT93, IADR=',IX, '  WORD=',MAXSZ
                          ELSE
                               WRITE(IOUT6,301) ISM,ISM,JSM,JSM
C     WRITE(ICHECK,301) ISM,ISM,JSM,JSM
                               NBKT2=NWRT
                               CALL SORTWO(BKTPNT,START,ILAST2,BKTSZ,
     *                         SRTAR,INTOFF,LEN2,NWRT,MAXSZ,SPACE,RESRT,
     *                         IADR,INTF,NBKT,LBLS,INTS,LEN2/2,BKTR,
     *                         NBKT2,SWOUT,SRTI)
                          END IF
 21                       CONTINUE
  501                CONTINUE
  500           CONTINUE
           ELSE IF (IPASS .EQ. 3) THEN
                II=0
                DO 510 ISM=1,NSYMHF
                     DO 511 JSM=1,ISM
                          IF (ISM .EQ. JSM) GOTO 22
                          II=II+1
                          IADR=ILAST(II)
                          IF (IADR .EQ. -1) GOTO 22
C     WRITE(IOUT6,'('  ASM=',I4,'  BSM=',I4,'  IBKT=',I4)')ISM,JSM,II
                          DO 512 JJ=1,LENR
                               RESRT(JJ)=0.0D+00
  512                     CONTINUE
                          MAXSZ=(IOFF(NSO(ISM))+NSO(ISM))*NSO(JSM)*NSO(
     *                    JSM)*2
C     WRITE(IOUT6,1000) (RESRT(IOP) ,IOP=1,MAXSZ/2)
                          NWRT=MAXSZ/LEN2
                          IF (NWRT*LEN2 .NE. MAXSZ) NWRT=NWRT+1
C     PRINT,'  NWRT=',NWRT,' MAXSZ=',MAXSZ
                          IF (NWRT .EQ. 1) THEN
 32                            CALL RREAD(IT91,LBLS,BKTSZ,IADR)
C     PRINT,BKTSZ,' WRDS READ FROM IT91 AT ADR ',IADR
                               IADR=LBLS(1)
C     PRINT,'  INXT=',IADR,'  NUMINTS=',LBLS(2)
                               DO 513 JJ=1,LBLS(2)
                                    RESRT(LBLS(JJ+2))=INTS(JJ+INTF)
C     WRITE(IOUT6,'(' RESRT(',I6,')=',G18.8)')LBLS(JJ+2),INTS(JJ+INTF)
  513                          CONTINUE
                               IF (IADR .NE. -1) GOTO 32
                               IWRD=MAXSZ
                               IST=1
CDB   CALL RGETSA(IT93,IX)
                               CALL SWRIT(IT93,IRESRT,MAXSZ)
                               SWOUT=SWOUT+I2SEC(MAXSZ)
C     PRINT,'  WRITE TO IT93, IADR=',IX, '  WORDS=',MAXSZ
                          ELSE
                               WRITE(IOUT6,301)ISM,JSM,ISM,JSM
C     WRITE(ICHECK,301)ISM,JSM,ISM,JSM
                               NBKT2=NWRT
                               CALL SORTWO(BKTPNT,START,ILAST2,BKTSZ,
     *                         SRTAR,INTOFF,LEN2,NWRT,MAXSZ,SPACE,RESRT,
     *                         IADR,INTF,NBKT,LBLS,INTS,LEN2/2,BKTR,
     *                         NBKT2,SWOUT,SRTI)
                          END IF
 22                       CONTINUE
  511                CONTINUE
  510           CONTINUE
                IF (NSYM0 .LT. 4) GOTO 250
           ELSE IF (IPASS .EQ. 4) THEN
                ISTOP=3
                IF (NSYMHF .EQ. 8) ISTOP=42
                DO 520 II=1,ISTOP
C     PRINT,' II=',II
                     IADR=ILAST(II)
                     IF (IADR .EQ. -1) GOTO 23
                     DO 521 JJ=1,LENR
                          RESRT(JJ)=0.0D+00
  521                CONTINUE
                     MAXSZ=MAXD2H(II)*2
C     PRINT,'  MAXSZ=',MAXSZ
C     WRITE(ICHECK,'('  MAXSZ='I6)')MAXSZ
                     NWRT=MAXSZ/LEN2
                     IF (NWRT*LEN2 .NE. MAXSZ) NWRT=NWRT+1
C     PRINT,'  NWRT=',NWRT, '  MAXSZ=',MAXSZ
C     WRITE(IOUT6,'(7G18.8)') (RESRT(IOP),IOP=1,MAXSZ/2)
                     IF (NWRT .EQ. 1) THEN
 33                       CALL RREAD(IT91,LBLS,BKTSZ,IADR)
                          IADR=LBLS(1)
                          DO 522 JJ=1,LBLS(2)
                               RESRT(LBLS(JJ+2))=INTS(JJ+INTF)
  522                     CONTINUE
                          IF (IADR .NE. -1) GOTO 33
                          IWRD=MAXSZ
                          IST=1
CDB   CALL RGETSA(IT93,IX)
                          CALL SWRIT(IT93,IRESRT,MAXSZ)
                          SWOUT=SWOUT+I2SEC(MAXSZ)
CDB   PRINT,'  WRITE TO IT93, IADR=',IX,'  WORDS=',MAXSZ
                     ELSE
                          WRITE(IOUT6,301)IAA(II),IB(II),IC(II),ID(II)
C     WRITE(ICHECK,301)IAA(II),IB(II),IC(II),ID(II)
                          NBKT2=NWRT
                          CALL SORTWO(BKTPNT,START,ILAST2,BKTSZ,SRTAR,
     *                    INTOFF,LEN2,NWRT,MAXSZ,SPACE,RESRT,IADR,INTF,
     *                    NBKT,LBLS,INTS,LEN2/2,BKTR,NBKT2,SWOUT,SRTI)
                     END IF
 23                  CONTINUE
  520           CONTINUE
           END IF
  900 CONTINUE
 250  CONTINUE
C
C
C
      WRITE(IOUT6,703) SRIN,IT34,IPCNT,SWOUT,IT93,FZC
  703 FORMAT (/,5X,I8,' SECTORS  READ FROM TAPE ',I4,' IN ',I3,' PASSES'
     #,       /,5X,I8,' SECTORS WRITTEN TO TAPE ',I4,//,' FROZEN CORE'
     #,       ' ENERGY = ',G19.9)
C
C50   FORMAT(7G18.8)
      RETURN
      END
