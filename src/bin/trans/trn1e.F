      SUBROUTINE TRN1E(NSOFF,H,U,OFFA1,P1,NSO,MAX,NSYMHF,MOINT1,
     1IOFF,NSOMO,MOOFF,IJGRP,IJADD,KADD,LADD,MOARR,MOBKT,MOLBL,MOSRT,
     2LENMO,BFDRT,ORBSYM,ACRCY,LCIBUF,IH)
C
C  TRANSFORM ONE ELECTRON INTEGRALS
C
C  MOARR : MOBKT : MOLBL : MOSRT ARE ALL IMPLICITLY EQUIVALENCED
C
C
      IMPLICIT REAL*8 (A-H,O,R-Z)
      IMPLICIT INTEGER (I-N,P,Q)
C
C
      integer norbs,nsymhf,nij,numij,lenmo,lcibuf,max
      INTEGER STARTT,ENDD,TIME,IH(1)
      INTEGER SYMORB,ORBSYM(NORBS),MULTS(8),LKUPSM(64),SRIN,SWOUT
      INTEGER NSOFF(NSYMHF),NSO(NSYMHF),OFFA1(NSYMHF),IOFF(NIJ),
     1NSOMO(NSYMHF),MOOFF(NSYMHF)
      INTEGER IJGRP(1),IJADD(NUMIJ),KADD(SYMORB),LADD(SYMORB)
      INTEGER MOLBL(LENMO),MOARR(LCIBUF),MOBKT(LENMO),BFDRT(NORBS)
      REAL*8 H(NIJ),P1(MAX),MOINT1(NUMIJ),P1I,MOSRT(LENMO),
     1ACRCY
C
      COMMON /TAPES/  IT34,IT58,IT47,IT52,IT91,IT93,IT95
      COMMON /OUT/    IOUT6,ICHECK,IDBG,ICNT
      COMMON /DIMS/   SYMORB,NUMIJ,NIJ,NBF,NORBS,NMAX,NBFORB,NSYM0
C
      DIMENSION U(NBFORB)
C
      DATA MULTS/0,8,16,24,32,40,48,56/
      DATA LKUPSM/1,2,3,4,5,6,7,8,2,1,4,3,6,5,8,7,3,4,1,2,7,8,5,6,4,3,
     A2,1,8,7,6,5,5,6,7,8,1,2,3,4,6,5,8,7,2,1,4,3,7,8,5,6,3,4,1,2,
     B8,7,6,5,4,3,2,1/

clj
c     write(6,*) 'trn1e: u:'
c     call matout(u,norbs,nbf,norbs,nbf,6)
C
CDB   IDBG16=0
CDB   IF((IDBG.AND.16).NE.0) IDBG16=1
      SRIN=0
      SWOUT=0
      WRITE(IOUT6,701)
  701 FORMAT (/,' TRANSFORMING ONE-ELECTRON INTEGRALS:')
      CALL SREW(IT91)
      CALL SREW(IT93)
      CALL RGETSA(IT93,IX)
      CALL SREAD(IT93,IH,INTOWP(NIJ))
      SRIN=SRIN+I2SEC(INTOWP(NIJ))
clj
c     write(6,*) 'trn1e: h:'
c     call print(h,nij,nbf,6)
c     write(6,*) 'nsomo', (nsomo(i),i=1,nsymhf)
c     write(6,*) 'nsoff', (nsoff(i),i=1,nsymhf)
c     write(6,*) 'mooff', (mooff(i),i=1,nsymhf)
C     PRINT,'  READ 1E INT IADR=',IX,'  WORD=',NIJ*2
C     WRITE(IOUT6,703) (H(IOP),IOP=1,NIJ)
C 703 FORMAT (/,(5E15.5))
      DO 100 ISM=1,NSYMHF
           IF (NSOMO(ISM).EQ.0) GOTO 99
           NN=NSOFF(ISM)
           NNMO=MOOFF(ISM)
           DO 40 P=1,NSO(ISM)
                PNN=(P+NN-1)*NORBS+NNMO
                PQ=IOFF(P+NN)+NN
                DO 5 II=1,MAX
                     P1(II)=0.0D+00
    5           CONTINUE
                E=1.0D+00
                DO 20 Q=1,P
                     QNN=(Q+NN-1)*NORBS+NNMO
                     IF (P.EQ.Q) E=0.5D+00
                     PQ=PQ+1
                     F=H(PQ)*E
                     IF (F .EQ. 0.0D+00) GOTO 15
                     DO 10 J=1,NSOMO(ISM)
                          P1(J)=P1(J)+U(QNN+J)*F
   10                CONTINUE
 15                  CONTINUE
   20           CONTINUE
                IJ=0
                DO 30 I=1,NSOMO(ISM)
                     UPNNI=U(PNN+I)
                     P1I=P1(I)
                     IJG=IOFF(I+NNMO)+NNMO
                     DO 25 J=1,I
                          IJG=IJG+1
                          IJ=IJ+1
                          MOINT1(IJG)=MOINT1(IJG)+P1I*U(PNN+J)+P1(J)*
     *                    UPNNI
   25                CONTINUE
   30           CONTINUE
   40      CONTINUE
 99        CONTINUE
  100 CONTINUE
clj
c     write(6,*) 'trn1e: moint1:'
c     call print(moint1,numij,norbs,6)
C     PRINT, '  MOINT1:'
C     WRITE(IOUT6,'(7G18.8)') (MOINT1(IOP), IOP=1,NUMIJ)
      ICNT=0
      KL=0
      DO 220 KORB=1,NORBS
           DO 215 LORB=1,KORB
                KL=KL+1
                VAL=MOINT1(KL)
                IF (ABS(VAL) .LT. ACRCY) GOTO 210
                K=BFDRT(KORB)
                L=BFDRT(LORB)
                IF (K .LT. L) THEN
                     ITMP=K
                     K=L
                     L=ITMP
                END IF
                MB=3
                IF (K .EQ. L) MB=2
                IJA=IOFF(K+1)
                LSM=ORBSYM(K)
                IF (LSM .NE. ORBSYM(L)) GOTO 210
                ICNT=ICNT+1
                LSMPT=L+(LSM-1)*NORBS
                MOLBL(ICNT)=IJADD(IJA)+KADD(K)+LADD(LSMPT)+MB
                MOBKT(ICNT)=IJGRP(IJA)
                MOSRT(ICNT)=VAL
C     WRITE(19,'(G18.8,I5)') VAL,LAD
 210            CONTINUE
  215      CONTINUE
  220 CONTINUE
      MOARR(1)=1
      MOARR(2)=ICNT
CDB   CALL RGETSA(IT91,IX)
CDB   PRINT,'  WRITTEN ',ICNT,' WORDS TO IT91 AT ADR',IX
      CALL SWRIT(IT91,MOARR,LENMO*4+2)
      SWOUT=SWOUT+I2SEC(LENMO*4+2)
C     WRITE(IOUT6,870) (IQQ,MOBKT(IQQ),MOLBL(IQQ),MOSRT(IQQ),IQQ=1,ICNT)
C 870 FORMAT (1X,3I4,E15.5)
C     PRINT,' ICNT=',ICNT
C     PRINT,' MOSRT'
C     WRITE(IOUT6,'(7G18.8)') (MOSRT(IOP), IOP=1,ICNT)
C     PRINT,' MOLBL'
C     WRITE(IOUT6,'(7I8)') (MOLBL(IOP),IOP=1,ICNT)
C     PRINT,'  MOBKT'
C     WRITE(IOUT6,'(7I8)') (MOBKT(IOP),IOP=1,ICNT)
      ICNT=0
      WRITE(IOUT6,702) SRIN,IT93,SWOUT,IT91
  702 FORMAT (5X,I8,' SECTORS  READ FROM TAPE ',I3,/
     #,       5X,I8,' SECTORS WRITTEN TO TAPE ',I3)
      RETURN
      END
