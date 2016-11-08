C
C                           *****************
C---------------------------***** SORT2 *****--------------------------
C                           *****************
C
      SUBROUTINE SORT2(ITAP77,ITAP78,BUF,IBUF,LENBUF,BKT,IBKT,IBKTSP,
     .                 ITRI,MCHAIN,MAXBKT,JOUT,NBUCK,NIOBF,LENGTH,
     .                 JBUF,MAXVAL)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION BUF(LENBUF),IBUF(LENBUF*2),MCHAIN(MAXBKT)
      DIMENSION BKT(IBKTSP),IBKT(IBKTSP*2)
      COMMON/PAKQ/I1,I2,I3
C
   11 FORMAT ('                                                ')
   12 FORMAT ('************************************************')
      MINSEC=10
C     WRITE(*,*) '   IN SORT2 '
      CALL ZERO(BKT,IBKTSP)
C
CC    IF (NIOBF.LT.MINSEC) THEN
CC       WRITE (JOUT,13)
CC 13    FORMAT (' INSUFFICIENT CC TO SORT 2 ELECTRON INTEGRALS')
CC       CALL MABORT
CC    END IF
C
      INTLEN=(INTOWP(LENGTH)-2)/INTOWP(1)
      MAXVAL=INTOWP(INTLEN)/(1+INTOWP(1))
CTJL  IVOFF=(MAXVAL+2)/INTOWP(1)
      IVOFF=(MAXVAL+3)/INTOWP(1)
C
      IBFLEN=(INTOWP(LENBUF)-2)/INTOWP(1)
      MAXBUF=INTOWP(IBFLEN)/(1+INTOWP(1))
CTJL  IBOFF=(MAXBUF+2)/INTOWP(1)
      IBOFF=(MAXBUF+3)/INTOWP(1)
C     WRITE (JOUT,*) ' LENBUF=',LENBUF
C     WRITE (JOUT,*) ' MAXBUF=',MAXBUF
      CALL SREW(ITAP78)
C
      DO 212 IBFCNT=1,JBUF
CTJL     CALL SREAD(ITAP78,IBUF,INTOWP(LENBUF))
         CALL SREAD(ITAP78,BUF,INTOWP(LENBUF))
C
         CALL SETMBF(MBUF,IBUF(2))
C        WRITE (JOUT,*) ' MBUF=',MBUF
         DO 101 I=1,MBUF
                CALL SETMBF(I3,IBUF(2+I))
CTJL            I3=IBUF(I+2)
            CALL UNPAK
            KL=I1
C           IF (KL.EQ.0) GOTO 101
            MUNU=I2
            RINT=BUF(IBOFF+I)
C           WRITE (JOUT,11)
C           WRITE (JOUT,12)
C           WRITE (JOUT,*) ' *********  INTS AS READ INTO SORT2 *******'
C           WRITE (JOUT,11)
C           WRITE (JOUT,12)
C           WRITE (JOUT,*) KL,MUNU,BUF(IBOFF+I)
            IBUCK1=((KL-1)/ITRI)+1
            NOFF1=(IBUCK1-1)*LENGTH
            IFILL1=IBKT(INTOWP(NOFF1)+2)+1
            IF (IFILL1.GT.MAXVAL) THEN
C              WRITE (JOUT,11)
C              WRITE (JOUT,12)
C              WRITE (JOUT,*) ' IBUCK1=',IBUCK1
C              WRITE (JOUT,11)
C              WRITE (JOUT,12)
C              DO 14 M=1,MAXVAL
C                 I3=IBKT(INTOWP(NOFF1)+2+M)
C                 CALL UNPAK
C                 IJ=I1
C                 KL=I2
C                 WRITE (JOUT,*) I1,I2,BKT(NOFF1+IVOFF+M)
C  14          CONTINUE
               CALL BKTDMP(BKT(NOFF1+1),LENGTH,ITAP77,ICHAN1)
C              WRITE (JOUT,*) ICHAN1
               MCHAIN(IBUCK1)=ICHAN1
               CALL ZERO(BKT(NOFF1+1),LENGTH)
               IBKT(INTOWP(NOFF1)+1)=ICHAN1
               IFILL1=1
            END IF
CTJL
C     IF(IVOFF+IFILL1.GE.1.OR.IVOFF+IFILL1.LE.10) THEN
C     WRITE(*,*) ' IN SORT2:KL,MN,VAL=',KL,MUNU,RINT
C     END IF
CTJL
            IBKT(INTOWP(NOFF1)+2+IFILL1)=I3
            BKT(IVOFF+NOFF1+IFILL1)=RINT
            IBKT(INTOWP(NOFF1)+2)=IFILL1
  101    CONTINUE
  212 CONTINUE
C
      DO 22 KBKT=1,NBUCK
         NOFF=(KBKT-1)*LENGTH
C        WRITE (JOUT,11)
C        WRITE (JOUT,12)
C        WRITE (JOUT,*) ' KBKT=',KBKT
C        WRITE (JOUT,11)
C        WRITE (JOUT,12)
C        DO 24 M=1,MAXVAL
C           I3=IBKT(INTOWP(NOFF)+2+M)
C           CALL UNPAK
C           IJ=I1
C           KL=I2
C           WRITE (JOUT,*) IJ,KL,BKT(NOFF+IVOFF+M)
C  24    CONTINUE
         CALL BKTDMP(BKT(NOFF+1),LENGTH,ITAP77,ICHAN)
         MCHAIN(KBKT)=ICHAN
C        WRITE(*,*) '  IN SORT2 WROTE LAST BUFFER AT SECTOR',ICHAN
   22 CONTINUE
C
      CALL ZERO(BKT,IBKTSP)
      CALL ZERO(BUF,LENBUF)
C
      RETURN
      END
