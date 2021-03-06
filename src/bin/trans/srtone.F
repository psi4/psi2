      SUBROUTINE SRTONE(BKTPNT,START,ILAST,BKTSZ,NBKT,TTT,ITTT,
     1INTOFF,ITXXX)
C
C
C  TTT : ITTT ARE AN IMPLICITLY EQUIVALENCED REAL : INTEGER PAIR
C
C
      IMPLICIT INTEGER (A-Z)
C
C
      integer ibigr
      REAL*8 TTT(IBIGR),VAL
C
      COMMON /OUT/    IOUT6,ICHECK,IDBG,ICNT
      COMMON /BDIM/   IMAX2,JMAX2,KMAX2,LENR,LENAOR,IBIG,IBIGR,LENPWS
      COMMON /OFF/    MAXINT
      COMMON /SRT/    VAL,IADR,IDUM,IADBK,MNBK,MXBK
C
      DIMENSION ITTT(IBIG),BKTPNT(NBKT),START(NBKT)
      DIMENSION ILAST(NBKT),INTOFF(NBKT)
C
C     PRINT, '  IN SRTONE'
      CALL SREW(ITXXX)
      NUM=NUM+1
C     SWOUT=0
      BKTPNT(1)=1
      START(1)=0
      ILAST(1)=-1
      DO 1 II=2,NBKT
           ILAST(II)=-1
           START(II)=START(II-1)+BKTSZ
           BKTPNT(II)=1
    1 CONTINUE
C
C POINTR SHOULD BE SET EQUAL TO THE NUMBER OF THE FIRST SECTOR
C
      POINTR=1
      MAXINT=((BKTSZ-2)/6)*2
      IF(MAXINT.LE.0) THEN
C     WRITE(IOUT6,'(/' IN SRTONE, MAXINT=',I5)')MAXINT
           CALL MABORT
      END IF
      DO 2 II=1,NBKT
           INTOFF(II)=(START(II)+MAXINT)/2+1
    2 CONTINUE
C     PRINT,' NPASS= ',NPASS
C     PRINT,' IADBK= ',IADBK
C     PRINT,' NBKT= ',NBKT
C     PRINT,' BKTSZ= ',BKTSZ
C     PRINT,' MAXINT= ',MAXINT
C     PRINT,' BKTPNT= ',BKTPNT
C     PRINT,' START= ',START
C     PRINT,' ILAST= ',ILAST
C     PRINT,' INTOFF= ',INTOFF
C     PRINT,' DONE SRTONE  NUM=',NUM
      RETURN
C
C--------------------------------------------------------SRTOUT---
C
      ENTRY SRTOUT(BKTPNT,START,ILAST,BKTSZ,NBKT,TTT,ITTT,
     *INTOFF,ITXXX)
C
C     WRITE(16,'(' BKTPNT=',10I7)') (BKTPNT(II),II=1,NBKT)
C     WRITE(16,'(' START=',10I7)') (START(II),II=1,NBKT)
C     WRITE(16,'(' ILAST=',10I7)') (ILAST(II),II=1,NBKT)
C     IDUM=IJBKT(IJ)
      STRT=START(IDUM)
      PNT=BKTPNT(IDUM)
      OFSET=INTOFF(IDUM)
      ADRL=STRT+PNT+2
      ADRV=OFSET+PNT
      ITTT(ADRL)=IADR
      TTT(ADRV)=VAL
C     WRITE(IOUT6,'(' TTT(',I6,')= ',G18.8,'  ITTT(',I6,')= '
C    1,I6,' IDUM= ',I3)')ADRV,VAL,ADRL,IADR,IDUM
      BKTPNT(IDUM)=BKTPNT(IDUM)+1
      IF (ADRL .LT. (OFSET*2)) RETURN
      ITTT(STRT+1)=ILAST(IDUM)
      ILAST(IDUM)=POINTR
      POINTR=POINTR+I2SEC(BKTSZ)
      ITTT(STRT+2)=MAXINT
CDB   CALL RGETSA(ITXXX,IX)
      CALL SWRIT(ITXXX,ITTT(STRT+1),BKTSZ)
C     SWOUT=SWOUT+I2SEC(BKTSZ)
CDB   WRITE(IOUT6,'(/2X,I6,' WORDS WRITTEN TO ITXXX AT SEC ',I9,
CDB  1 ' IJBKT IS ',I5/)')BKTSZ,IX,IDUM
CDB   WRITE(IOUT6,'(' WRITTEN TO 101'/' ITTT(',I6,')= ',I6,' ITTT(',
CDB  1I6,')= ',I6,' POINTER= ',I6,' IDUM= ',I4)')STRT+1,ITTT(STRT+1),
CDB  1STRT+2,ITTT(STRT+2),POINTR,IDUM
C     WRITE(IOUT6,'(' ITTT:')')
C     WRITE(IOUT6,'(7I10)') (ITTT(II), II=STRT+1,STRT+MAXINT+2)
C     WRITE(IOUT6,'(' TTT')')
C     WRITE(IOUT6,'(7G18.8)') (TTT(II),II=OFSET+1,OFSET+MAXINT)
      BKTPNT(IDUM)=1
      RETURN
C
C----------------------------------------------------------SRTLST-----
C
      ENTRY SRTLST(BKTPNT,START,ILAST,BKTSZ,NBKT,TTT,ITTT,
     *INTOFF,ITXXX)
C
C     WRITE(IOUT6,'(' IN SRTLST')')
C     WRITE(IOUT6,'(' ILAST= ',5I7)') (ILAST(KKK), KKK=1,NBKT)
      DO 3 II=1,NBKT
           STRT=START(II)
C     PRINT, 'STRT=',STRT
           IF (BKTPNT(II) .EQ. 1) GOTO 201
           ITTT(STRT+1)=ILAST(II)
           ILAST(II)=POINTR
           POINTR=POINTR+I2SEC(BKTSZ)
           ITTT(STRT+2)=BKTPNT(II)-1
C     IF (NUM .EQ. 3)
C     WRITE(IOUT6,'(' ITTT:')')
C     WRITE(IOUT6,'(7I10)') (ITTT(JJ), JJ=STRT+1,BKTPNT(II)+STRT+1)
C     WRITE(IOUT6,'(' TTT')')
C     WRITE(IOUT6,'(7G18.8)') (TTT(JJ), JJ=INTOFF(II)+1,
C    1 INTOFF(II)+BKTPNT(II)-1)
C     END IF
C     PRINT,' CALLING RGETDA'
CDB   CALL RGETSA(ITXXX,IX)
           CALL SWRIT(ITXXX,ITTT(STRT+1),BKTSZ)
CDB   WRITE(IOUT6,'(/2X,I6,' WORDS WRITTEN TO ITXXX AT SEC ',I9,
CDB  1 ' IJBKT IS ',I5/)')BKTSZ,IX,II
CDB   WRITE(IOUT6,'(' WRITTEN TO 101'/' ITTT(',I6,')= ',I6,' ITTT(',
CDB  1 I6,')= ',I6,' POINTER= ',I6,' IDUM= ',I4)')STRT+1,ITTT(STRT+1),
CDB  1   STRT+2,ITTT(STRT+2),POINTR,II
C     SWOUT=SWOUT+I2SEC(BKTSZ)
 201       CONTINUE
    3 CONTINUE
C     PRINT, ' ILAST= ',ILAST
C     WRITE(IOUT6,'(/2X,I6,' SECTORS WRITTEN TO TAPE ',I4,
C    1 ' AFTER SORT')') SWOUT,ITXXX
      RETURN
      END
