C
C                           ******************
C---------------------------***** BCKCHN *****-------------------------
C                           ******************
C
      SUBROUTINE BCKCHN(MCHAIN,IBKT,ITAP77,BUF,IBUF,INTBUF,MAXBKT,
     .                  MAXVAL,RINTS,LENINT,NTRI,JOUT,ITRI)
 
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER IBUF(INTBUF*2),MCHAIN(MAXBKT)
      REAL*8 BUF(INTBUF),RINTS(LENINT)
      COMMON /PAKQ/ I1,I2,I3
C
      ISCTR=MCHAIN(IBKT)
C     WRITE(*,*) '  IN BCKCH READING FROM SECTOR ',ISCTR
 1111 CALL RREAD(ITAP77,BUF,INTOWP(INTBUF),ISCTR)
C     WRITE (JOUT,*) ' ISCTR=',ISCTR
C     WRITE(*,*) ' IBUF ',IBUF(1),IBUF(2)
C     WRITE(*,*) ' BUF ',BUF(1)
C     WRITE (JOUT,*) ' READ TAPE 77'
CTJL  ISCTR=IBUF(1)
      CALL SETMBF(ISCTR,IBUF(1))
C     WRITE (JOUT,*) ' ISCTR=',ISCTR
      CALL SETMBF(MBUF,IBUF(2))
C     WRITE (JOUT,*) ' MBUF=',MBUF,IBUF(2)
C     WRITE (JOUT,*) ' MAXVAL=',MAXVAL
CTJL  IVOFF=(MAXVAL+2)/INTOWP(1)
      IVOFF=(MAXVAL+3)/INTOWP(1)
C     IBOFF=ITRI*(IBKT-1)
C     WRITE (JOUT,*) ' IVOFF=',IVOFF
      DO 101 I=1,MBUF
C        WRITE (JOUT,*) ' I=',I
         CALL SETMBF(I3,IBUF(2+I))
CTJL     I3=IBUF(I+2)
         CALL UNPAK
         IJ=I1
         KL=I2
CTJL     IND=NTRI*(IJ-1)+KL
         IND=NTRI*(IJ-1-ITRI*(IBKT-1))+KL
C        WRITE (JOUT,*) ' IND=',IND
         RINTS(IND)=BUF(IVOFF+I)
CTJL
C     IF(IND.GE.1.OR.IND.LE.10) THEN
C     WRITE(*,*) ' IN BCKCH:IJ,KL,VAL=',IJ,KL,RINTS(IND)
C     END IF
CTJL
C        WRITE (JOUT,*) I1,I2,IND,IND-IBOFF,RINTS(IND-IBOFF)
  101 CONTINUE
      CALL ZERO(BUF,INTBUF)
      IF (ISCTR.NE.0) GOTO 1111
C
C     WRITE (JOUT,*) ' EXITING BACKCHAIN'
      RETURN
      END
