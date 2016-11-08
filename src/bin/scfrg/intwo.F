      SUBROUTINE INTWO (ALPHA,BETA,IJL,IJU,VAL1,VAL2,OBUF,BUF,INBUF
     1,                 PA,PB,PC,LBLINT,VAL0)
C
C
C
      IMPLICIT INTEGER(A -Z)
C
CTJL #,              LBLINT,VAL0
C
      INTEGER START,END,FCB102(16)
      REAL*8 BUF(1),VAL0(2),VAL1(1),VAL2(1),ALPHA(1),BETA(1)
      REAL*8 AIJKL,PA(1),PB(1),PC(1),HKT,PT,PKII
      DIMENSION OBUF(1),IJL(1),IJU(1),IJKL(2),LBLINT(2)
      EQUIVALENCE (AIJKL,IJKL)
C     EQUIVALENCE (VAL0(1),LBLINT(1))
      COMMON /POINT/ POINTR(4)
      COMMON /KIND/ NKIND,ICLOS
      COMMON /STIMEX/ KTIMES,ITIME5,ITIME6,ITIME7,ITIME8,ICONT5,ICONT6
     1,              ICONT7,ICONT8,ITIMED,ICONTD,ITEMPT
      COMMON /TAPES/ NIN,IR,IOUT,ICHECK,IFLAG3,ITAPE,I30
      COMMON /DIM21/ IPQ(256)
      COMMON /T102/ NOUT,HALF,THIRD,MAXNO,OUTBUF,KINADR
      COMMON /PRINTQ/ PRINT
      COMMON /ABOPEN/ NOST
      COMMON /SYMMS/ND(10),LLSYM(10)
      COMMON /PASS1/IBL,NINT,INTI,NBLI,INTO,NBLO,NEWBUF
      COMMON /PASS2/ PKII,III,ISM,JSM,KSM,LSM,IOR,JOR,KOR,LOR
      COMMON /SIGNS/ IOFF(200),MASK,IPRNT
      DATA NUMP /0/, NUMQ /0/
C
      OBUFLD=0
      IBUFLD=0
      OUTBUF=(OUTBUF-2)/12*6
      HALF=OUTBUF/2
      THIRD=OUTBUF/3
      OUTBUF=OUTBUF*2+2
      IF (NOST.EQ.0) THEN
         POINTR(1)=OUTBUF*2+1
         POINTR(2)=1
         POINTR(3)=OUTBUF*2+1
         POINTR(4)=OUTBUF+1
      ELSE
         DO 1 I=1,4
            POINTR(I)=(I-1)*OUTBUF+1
    1    CONTINUE
      END IF
      NEXT=OUTBUF+POINTR(4)
      N1I=POINTR(1)+1
      N2I=POINTR(2)+1
      N3I=POINTR(3)+1
      N4I=POINTR(4)+1
      N1IMAX=N1I+THIRD
      N2IMAX=N2I+HALF
C     N3RAMX=N3I/2+THIRD
      N3IMAX=N3I+THIRD*2
C     N4RAMX=N4I/2+HALF*2-1
cges  N4IMAX=N4I+HALF*2
      N4IMAX=N4I+third*2
      N1RA=N1I/2+THIRD
      N1RB=N1RA+THIRD
      N2R=N2I/2+HALF
C     N3RA=N3I/2
C     N3RB=N3RA+THIRD
      N3RB=N3I/2+THIRD
      N3RC=N3RB+THIRD
C     N4RA=N4I/2-1
c     N4R=N4I/2+HALF
      N4Rb=N4I/2+third
      N4Rc=N4rb+third
      OBUF(POINTR(1))=1
      OBUF(POINTR(1)+1)=THIRD
      OBUF(POINTR(2))=2
      OBUF(POINTR(2)+1)=HALF
      OBUF(POINTR(3))=3
      OBUF(POINTR(3)+1)=THIRD
      OBUF(POINTR(4))=4
C     OBUF(POINTR(4)+1)=HALF*4
cges  OBUF(POINTR(4)+1)=HALF
      OBUF(POINTR(4)+1)=third
C     CALL RFILE(NOUT)
      CALL SREW(NOUT)
      RETURN
C
      ENTRY HUSTON(ALPHA,BETA,IJL,IJU,VAL1,VAL2,OBUF,BUF,INBUF
     1,                 PA,PB,PC,LBLINT,VAL0)
CSEL  CALL ETIME2(START)
      ICONT7=ICONT7+1
      IBUFLD=IBUFLD+1
      IF (IAND(PRINT,32).NE.0) WRITE (ICHECK,900) IBUFLD,NINT,NINT
  900 FORMAT (/,' INPUT BLOCK ',I5,I15,' P AND ',I15,' K/2 INTEGRALS',/)
      DO 20 I=1,NINT
         IJ=LBLINT(I*2)
         KL=LBLINT(I*2-1)
         IF(IAND(PRINT,32).NE.0) PRINT '(8X,I6,I5,1P,2E14.6)',IJ,KL,
     *   VAL1(I),VAL2(I)
         IF (IJ.EQ.KL) THEN
            IF (NOST.NE.0) THEN
               DO 2 IOST=1,NOST
                  IF (IJ.LT.IJL(IOST)) GO TO 3
                  IF (IJ.LE.IJU(IOST)) THEN
                     N1I=N1I+1
                     N1RA=N1RA+1
                     N1RB=N1RB+1
                     OBUF(N1I)=IJ
                     BUF(N1RA)=VAL1(I)
                     IJOST=IPQ(IOST+1)
                     BUF(N1RB)=ALPHA(IJOST)*VAL1(I)+BETA(IJOST
     *               )*VAL2(I)
                     IF (N1I.GE.N1IMAX) THEN
                        OBUFLD=OBUFLD+1
CTJL                    CALL AWAIT(NOUT)
                        CALL SWRIT(NOUT,OBUF(POINTR(1)
     *                  ),                            OUTBUF
     *                  )
                        NUMP=NUMP+THIRD
                        NUMQ=NUMQ+THIRD
                        IF (IAND(PRINT,64).NE.0) THEN
                           N1I=POINTR(1)+1
                           N1RA=N1I/2+THIRD
                           N1RB=N1RA+THIRD
                           WRITE (ICHECK,901) OBUFLD,THIRD
     #                     ,
     *                     THIRD
                           DO 100 J=1,THIRD
                              PRINT '(I6,1P,2E14.6)',
     *                        OBUF(J+N1I),BUF(J+N1RA),
     *                        BUF(J+N1RB)
  100                      CONTINUE
                        END IF
                        SAVE=POINTR(1)
                        POINTR(1)=NEXT
                        NEXT=SAVE
                        N1I=POINTR(1)+1
                        N1IMAX=N1I+THIRD
                        N1RA=N1I/2+THIRD
                        N1RB=N1RA+THIRD
                        OBUF(POINTR(1))=1
                        OBUF(POINTR(1)+1)=THIRD
                     END IF
                     GO TO 20
                  END IF
    2          CONTINUE
            END IF
    3       CONTINUE
            N2I=N2I+1
            N2R=N2R+1
            OBUF(N2I)=IJ
            BUF(N2R)=VAL1(I)
            IF (N2I.GE.N2IMAX) THEN
               OBUFLD=OBUFLD+1
CTJL           CALL AWAIT(NOUT)
               CALL SWRIT(NOUT,OBUF(POINTR(2)),OUTBUF)
               NUMP=NUMP+HALF
               IF (IAND(PRINT,64).NE.0) THEN
                  N2I=POINTR(2)+1
                  N2R=N2I/2+HALF
                  WRITE (ICHECK,902) OBUFLD,HALF
                  DO 101 J=1,HALF
                     PRINT '(I14,1PE14.6)', OBUF(J+N2I),BUF(J+
     *               N2R)
  101             CONTINUE
               END IF
               SAVE=POINTR(2)
               POINTR(2)=NEXT
               NEXT=SAVE
               N2I=POINTR(2)+1
               N2R=N2I/2+HALF
               N2IMAX=N2I+HALF
               OBUF(POINTR(2))=2
               OBUF(POINTR(2)+1)=HALF
            END IF
         ELSE
            IF (NOST.NE.0) THEN
               DO 5 IOST=1,NOST
                  IF (IJ.LT.IJL(IOST)) GO TO 6
                  IF (IJ.LE.IJU(IOST)) THEN
                     DO 4 JOST=1,IOST
                        IF (KL.LT.IJL(JOST)) GO TO 10
                        IF (KL.LE.IJU(JOST)) THEN
C     N3RA=N3RA+1
                           N3I=N3I+2
                           N3RB=N3RB+1
                           N3RC=N3RC+1
C     BUF(N3RA)=VAL0(I)
                           OBUF(N3I-1)=KL
                           OBUF(N3I  )=IJ
                           BUF(N3RB)=VAL1(I)
                           IJOST=IPQ(IOST)+JOST
                           BUF(N3RC)=ALPHA(IJOST)*VAL1(I)+
     *                     BETA(IJOST)*VAL2(I)
C     IF (N3RA.GE.N3RAMX) THEN
                           IF (N3I.GE.N3IMAX) THEN
                              OBUFLD=OBUFLD+1
CTJL                          CALL AWAIT(NOUT)
                              CALL SWRIT(NOUT,OBUF(
     *                        POINTR(3)),OUTBUF)
                              NUMP=NUMP+THIRD
                              NUMQ=NUMQ+THIRD
                              IF (IAND(PRINT,64).NE.0)
     *                        THEN
                                 WRITE (ICHECK,901) OBUFLD,
     #                           THIRD,THIRD
                                 N3I=POINTR(3)
                                 N3RB=N3I/2+THIRD
                                 N3RC=N3RB+THIRD
                                 DO 102 J=1,THIRD
                                    PRINT '(I14,I5,1P,2E14.6)', OBUF(J*2
     *                              +1+N3I),OBUF(J*2+2+N3I),BUF(J+N3RB),
     *                              BUF(J+N3RC)
  102                            CONTINUE
                              END IF
                              SAVE=POINTR(3)
                              POINTR(3)=NEXT
                              NEXT=SAVE
                              N3I=POINTR(3)+1
C     N3RA=N3I/2
C     N3RB=N3RA+THIRD
                              N3RB=N3I/2+THIRD
                              N3RC=N3RB+THIRD
C     N3RAMX=N3RA+THIRD
                              N3IMAX=N3I+THIRD*2
                              OBUF(POINTR(3))=3
                              OBUF(POINTR(3)+1)=THIRD
                           END IF
                           GO TO 20
                        END IF
    4                CONTINUE
                  END IF
    5          CONTINUE
            END IF
    6       CONTINUE
C  10 N4RA=N4RA+2
C     IJKL(1)=KL
C     IJKL(2)=IJ
C     PRINT '(' IJ=',I5,' KL=',I5,'  VAL=',E18.8)',IJ,KL,VAL1(I)
C***C
C     BUF(N4RA)=AIJKL
C     BUF(N4RA+1)=VAL1(I)
C     IF (N4RA.GE.N4RAMX) THEN
   10       CONTINUE
            N4I=N4I+2
cges        N4R=N4R+1
            N4Rb=N4Rb+1
            N4Rc=N4Rc+1
            OBUF(N4I-1)=KL
            OBUF(N4I  )=IJ
cges        BUF(N4R  )=VAL1(I)
            BUF(N4Rb )=VAL1(I)
            BUF(N4Rc )=VAL2(I)
            IF (N4I.GE.N4IMAX) THEN
               OBUFLD=OBUFLD+1
               IF(OBUF(POINTR(4)).NE.4) STOP
CTJL           CALL AWAIT(NOUT)
               CALL SWRIT(NOUT,OBUF(POINTR(4)),OUTBUF)
               NUMP=NUMP+HALF
c              IF (IAND(PRINT,64).NE.0) THEN
c                 WRITE (ICHECK,902) OBUFLD,HALF
c                 N4I=POINTR(4)
c                 N4RB=N4I/2+HALF
c                 DO 103 J=1,HALF
c                    PRINT '(I14,I5,1PE14.6)', OBUF(J*2+1+N4I)
c    *               ,OBUF(J*2+2+N4I),BUF(J+N4RB)
c 103             CONTINUE
c              END IF
               SAVE=POINTR(4)
               POINTR(4)=NEXT
               NEXT=SAVE
               N4I=POINTR(4)+1
C     N4RA=N4I/2-1
C     N4RAMX=N4RA+HALF*2
cges           N4R=N4I/2+HALF
               N4Rb=N4I/2+third
               N4Rc=N4rb+third
               N4IMAX=N4I+third*2
               OBUF(POINTR(4))=4
C     OBUF(POINTR(4)+1)=HALF*4
               OBUF(POINTR(4)+1)=third
            END IF
         END IF
   20 CONTINUE
CSEL  CALL ETIME2(END)
      ITIME7=ITIME7+END-START
      RETURN
C
      ENTRY CLNUP(ALPHA,BETA,IJL,IJU,VAL1,VAL2,OBUF,BUF,INBUF
     1,                 PA,PB,PC,LBLINT,VAL0)
      IF (N1I.NE.POINTR(1)+1) THEN
         OBUFLD=OBUFLD+1
         NOINBF=N1I-POINTR(1)-1
         OBUF(POINTR(1)+1)=NOINBF
CTJL     CALL AWAIT(NOUT)
         CALL SWRIT(NOUT,OBUF(POINTR(1)),OUTBUF)
         NUMP=NUMP+NOINBF
         NUMQ=NUMQ+NOINBF
         IF (IAND(PRINT,64).NE.0) THEN
            N1I=POINTR(1)+1
            N1RA=N1I/2+THIRD
            N1RB=N1RA+THIRD
            WRITE (ICHECK,901) OBUFLD,NOINBF,NOINBF
            DO 104 J=1,NOINBF
               PRINT '(I6,1P,2E14.6)', OBUF(J+N1I),BUF(J+N1RA),BUF(
     *         J+N1RB)
  104       CONTINUE
         END IF
      END IF
      IF (N2I.NE.POINTR(2)+1) THEN
         OBUFLD=OBUFLD+1
         NOINBF=N2I-POINTR(2)-1
         OBUF(POINTR(2)+1)=NOINBF
CTJL     CALL AWAIT(NOUT)
         CALL SWRIT(NOUT,OBUF(POINTR(2)),OUTBUF)
         NUMP=NUMP+NOINBF
         IF (IAND(PRINT,64).NE.0) THEN
            N2I=POINTR(2)+1
            N2R=N2I/2+HALF
            WRITE (ICHECK,902) OBUFLD,NOINBF
            DO 105 J=1,NOINBF
               PRINT '(I14,1PE14.6)', OBUF(J+N2I),BUF(J+N2R)
  105       CONTINUE
         END IF
      END IF
C     IF (N3RA.GT.(POINTR(3)+1)/2) THEN
      IF (N3I.GT.(POINTR(3)+1)) THEN
C     NOINBF=N3RA-(POINTR(3)+1)/2
         NOINBF=(N3I-(POINTR(3)+1))/2
         OBUFLD=OBUFLD+1
         OBUF(POINTR(3)+1)=NOINBF
CTJL     CALL AWAIT(NOUT)
         CALL SWRIT(NOUT,OBUF(POINTR(3)),OUTBUF)
         NUMP=NUMP+NOINBF
         NUMQ=NUMQ+NOINBF
         IF (IAND(PRINT,64).NE.0) THEN
            WRITE (ICHECK,901) OBUFLD,NOINBF,NOINBF
            N3I=POINTR(3)+1
            N3RB=N3I/2+THIRD
            N3RC=N3RB+THIRD
            DO 106 J=1,NOINBF
               PRINT '(I14,I5,1P,2E14.6)', OBUF(J*2-1+N3I),OBUF(J*2
     *         +N3I),BUF(J+N3RB),BUF(J+N3RC)
  106       CONTINUE
         END IF
      END IF
C     IF (N4RA.GT.(POINTR(4)+1)/2-1) THEN
C     NOINBF=(N4RA-(POINTR(4)+1)/2+1)*2
      IF (N4I.GT.(POINTR(4)+1)) THEN
         NOINBF=(N4I-(POINTR(4)+1))/2
         OBUFLD=OBUFLD+1
         OBUF(POINTR(4)+1)=NOINBF
CTJL     CALL AWAIT(NOUT)
         CALL SWRIT(NOUT,OBUF(POINTR(4)),OUTBUF)
C     NUMP=NUMP+NOINBF/4
         NUMP=NUMP+NOINBF
         IF (IAND(PRINT,64).NE.0) THEN
            WRITE (ICHECK,902) OBUFLD,NOINBF
            N4I=POINTR(4)+1
            N4RB=N4I/2+HALF
            DO 107 J=1,NOINBF
               PRINT '(I14,I5,1PE14.6)', OBUF(J*2-1+N4I),OBUF(J*2+
     *         N4I),BUF(J+N4RB)
  107       CONTINUE
         END IF
      END IF
      MAXNO=OBUFLD
      WRITE (IOUT,903) NUMP,NUMQ,MAXNO,MAXNO*I2SEC(OUTBUF)
      IF (IFLAG3.GT.0) WRITE (IFLAG3,903) NUMP,NUMQ,MAXNO,MAXNO*
     #I2SEC(OUTBUF)
C
CTJL  CALL AWAIT(NOUT)
C
      RETURN
C
  901 FORMAT (/,' OUTPUT BLOCK ',I5,I15,' P AND ',I6,' Q INTEGRALS',/)
  902 FORMAT (/,' OUTPUT BLOCK ',I5,I15,' P INTEGRALS',/)
  903 FORMAT (/,1X,I12,' P AND ',I12,' Q INTEGRALS STORED ON FILE92 IN '
     #,       I4,' BUFFERS',/,T35,'TOTAL SPACE USED IS ',I6,' SECTORS')
C
      END
