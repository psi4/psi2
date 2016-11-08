C
C----------------------------------------------------------------------
C
      SUBROUTINE SORT (BUF,IBUF,LENGTH,CCTSCF,SCFTCC,NPQ,NRS,
     .                 NO,NV,NT,ITAPIN,ITAPOU,JOUT,NIRRED,FLOV,ORBSYM,
     .                 IPQ,ITV,IOFF,IADD,GBUCK,CORE,LEFT,FZT,XAR1,
     .                 I1OFF,I1ADD,X1INT,C2OFF,C2ADD,C2INT,
     .                 D2OFF,D2ADD,D2INT,E3OFF,E3ADD,E3INT)
      IMPLICIT INTEGER (A-Z)
      REAL*8 BUF(LENGTH),RINT,CORE(LEFT),GBUCK(NT*NT),VAL,XAR1(*),
     .       X1INT(*),C2INT(*),D2INT(*),E3INT(*)
      INTEGER IBUF(LENGTH*2),FZT(*),E3OFF(*),E3ADD(*),
     .        ORBSYM(NT),FLOV(NIRRED,2),IOFF(*),IADD(*),
     .        ITV(*),IPQ(*),CCTSCF(NT),SCFTCC(NT),NPQ(*),NRS(*),
     .        I1OFF(*),C2OFF(*),I1ADD(*),C2ADD(*),D2OFF(*),D2ADD(*)
C
C >>> THIS SORT LIMITED TO 500 LOOPS
C
      INTEGER FPAIR(500),LPAIR(500),FINT(500),LINT(500)
      DATA ITEMP /255/
C
      ITAP68=68
      ITAP78=78
      CALL RFILE(ITAPOU)
C
      WRITE(6,*)
      WRITE(6,*) '>>>>>>   IN SORT <<<<<<<<<'
      WRITE(6,*)
C
C >>> CONSTRUCT SYMMETRY-PACKING ARRAYS FOR EXIT FILE
C
      TOFF=0
      SOFF=0
      PQPAIR=0
      LOOP=1
      FPAIR(1)=1
      LPAIR(1)=NT*(NT+1)/2
      FINT(1)=1
      LINT(1)=0
      DO 190 TSYM=1,NIRRED
         NOFF=0
         NPQ(TSYM)=0
         DO 180 P=1,NT
            IF(FZT(P).EQ.1) GO TO 180
            PSYM=ORBSYM(P)
            DO 170 Q=1,P
               IF(FZT(Q).EQ.1) GO TO 170
               QSYM=ORBSYM(Q)
               PQSYM=IEOR(PSYM,QSYM)
               IF(PQSYM.NE.TSYM-1)GO TO 170
               PQ=IPQ(P)+Q
               IOFF(PQ)=TOFF
               NPQ(TSYM)=NPQ(TSYM)+1
               PQPAIR=PQPAIR+1
               NRS(TSYM)=0
               NADD=0
               DO 150 R=1,NT
                  IF(FZT(R).EQ.1) GO TO 150
                  RSYM=ORBSYM(R)
                  SSYM=IEOR(RSYM,PQSYM)+1
                  FS=FLOV(SSYM,1)
                  LS=FLOV(SSYM,2)
                  IF(LS.GT.R)LS=R
                  DO 140 S=FS,LS
                     RS=IPQ(R)+S
                     NOFF=NOFF+1
                     NADD=NADD+1
                     IADD(RS)=NADD
                     TOFF=TOFF+1
                     SOFF=SOFF+1
                     TOT=IOFF(PQ)+IADD(RS)
                     NRS(TSYM)=NRS(TSYM)+1
C     WRITE(6,634)P,Q,PQ,IOFF(PQ),R,S,RS,IADD(RS),TOT,TOFF
  634 FORMAT(10(1X,I6))
  140             CONTINUE
  150          CONTINUE
               IF(SOFF.LT.LEFT)THEN
                  LINT(LOOP)=TOFF
                  LPAIR(LOOP)=PQPAIR
               ELSE
                  LOOP=LOOP+1
                  FINT(LOOP)=LINT(LOOP-1)+1
                  FPAIR(LOOP)=LPAIR(LOOP-1)+1
                  LINT(LOOP)=TOFF
                  LPAIR(LOOP)=PQPAIR
                  SOFF=NADD
               ENDIF
  170       CONTINUE
  180    CONTINUE
  190 CONTINUE
C
      WRITE(6,*)
      WRITE(6,*)LOOP,' LOOPS ARE NEEDED TO RESORT FILE',ITAPIN
      WRITE(6,*)
      IF(LOOP.GT.500)THEN
         WRITE(6,*)' TOO MANY LOOPS !!   GET LOST! '
         STOP ' TOO MANY LOOPS IN SORT '
      ENDIF
      WRITE(6,*  )'    I   FPAIR(I)    LPAIR(I)    FINT(I)  LINT(I)'
      DO 231 I=1,LOOP
      WRITE(6,654)I,FPAIR(I),LPAIR(I),FINT(I),LINT(I)
  654 FORMAT(5(2X,I8))
  231 CONTINUE
      WRITE(6,*)
      WRITE(6,*  )'   ISYM        NPQ(ISYM)      NRS(ISYM)'
      TPAR=0
      SIZEOU=0
      DO 232 I=1,NIRRED
      WRITE(6,652)I,NPQ(I),NRS(I)
      TPAR=TPAR+NPQ(I)
      SIZEOU=SIZEOU+NPQ(I)*NRS(I)
  652 FORMAT(3(2X,I8))
  232 CONTINUE
      WRITE(6,*)
      WRITE(6,*)'NUMBER OF ELEMENTS IN DISK IS =',SIZEOU
      WRITE(6,*)
C
      INTLEN=(INTOWP(LENGTH)-2)/INTOWP(1)
      MAXVAL=INTOWP(INTLEN)/(1+INTOWP(1))
      IVOFF=(MAXVAL+3)/INTOWP(1)
      WRITE (JOUT,*) ' LENGTH',LENGTH
      WRITE (JOUT,*) ' INTLEN',INTLEN
      WRITE (JOUT,*) ' MAXVAL',MAXVAL
      WRITE (JOUT,*) ' IVOFF',IVOFF
C
      POINT=1
      DO 1000 ILOOP=1,LOOP
      CALL ZERO(GBUCK,NT*NT)
      CALL ZERO(CORE ,LEFT )
      FFINT =FINT (ILOOP)
      LLINT =LINT (ILOOP)
      FFPAIR=FPAIR(ILOOP)
      LLPAIR=LPAIR(ILOOP)
      SHIFT=FFINT-1
      WRITE(6,*)'LOOP=',ILOOP,' SHIFT=',SHIFT
      CALL SREW(ITAPIN)
      IF(ITAPIN.EQ.ITAP78)CALL SREAD(ITAP78,JUNK,1)
  111 CALL ZERO(BUF,LENGTH)
      CALL SREAD(ITAPIN,IBUF,INTOWP(LENGTH))
C
      IFLG=IBUF(1)
C     WRITE (JOUT,*) ' IFLG2=',IFLG
      MBUF=IBUF(2)
C     WRITE (JOUT,*) ' MBUF2=',MBUF
      DO 101 II=1,MBUF
         ABCD=IBUF(2+II)
         A=ISHFT(ABCD,-24)
         B=IAND(ITEMP,ISHFT(ABCD,-16))
         C=IAND(ITEMP,ISHFT(ABCD,-8))
         D=IAND(ITEMP,ABCD)
          A=CCTSCF(A)
          B=CCTSCF(B)
          C=CCTSCF(C)
          D=CCTSCF(D)
         RINT=BUF(IVOFF+II)
C
         FREEZE=FZT(A)+FZT(B)+FZT(C)+FZT(D)
C
C >>> PUT IN CANONICAL ORDER
C
      IF (A.LT.B) THEN
         ISWOP=B
         B=A
         A=ISWOP
      END IF
      IF (C.LT.D) THEN
         ISWOP=D
         D=C
         C=ISWOP
      END IF
      IF (IPQ(A)+B.LT.IPQ(C)+D) THEN
         ISWOP=C
         C=A
         A=ISWOP
         ISWOP=D
         D=B
         B=ISWOP
      END IF
      AB=IPQ(A)+B
      CD=IPQ(C)+D
C
         IF(ITAPIN.EQ.ITAP68)THEN
            RINT=RINT*0.125D0
         ENDIF
C
   19 FORMAT (4I3,F20.12,3X,'* CANONICAL SCF ORDER *')
C
      IF(FREEZE.EQ.0)THEN
C
         ABCD=IOFF(AB)+IADD(CD)
         IF(ABCD.GE.FFINT.AND.ABCD.LE.LLINT) THEN
            ABCD=ABCD-SHIFT
            CORE(ABCD)=RINT
         ENDIF
C        WRITE (JOUT,619) FREEZE,A,B,C,D,ABCD,RINT
C
         CDAB=IOFF(CD)+IADD(AB)
         IF(CDAB.GE.FFINT.AND.CDAB.LE.LLINT) THEN
            CDAB=CDAB-SHIFT
            CORE(CDAB)=RINT
         ENDIF
C
      ELSE IF(ILOOP.EQ.LOOP.AND.FREEZE.EQ.1.AND.ITAPIN.EQ.ITAP78)THEN
         IF(FZT(A).EQ.1.OR.FZT(B).EQ.1)THEN
            ABCD=I1OFF(AB)+I1ADD(CD)
         ELSE
            ABCD=I1OFF(CD)+I1ADD(AB)
         ENDIF
         X1INT(ABCD)=RINT
C        WRITE (JOUT,619) FREEZE,A,B,C,D,ABCD,RINT
      ELSE IF(ILOOP.EQ.LOOP.AND.FREEZE.EQ.2.AND.ITAPIN.EQ.ITAP78)THEN
         IF(FZT(A).EQ.1.AND.FZT(B).EQ.1)THEN
            ABCD=C2OFF(AB)+C2ADD(CD)
            C2INT(ABCD)=RINT
C        WRITE (JOUT,618) FREEZE,A,B,C,D,ABCD,RINT
         ELSE IF(FZT(C).EQ.1.AND.FZT(D).EQ.1)THEN
            ABCD=C2OFF(CD)+C2ADD(AB)
            C2INT(ABCD)=RINT
C        WRITE (JOUT,618) FREEZE,A,B,C,D,ABCD,RINT
crsg     ELSE
crsg        ABCD=D2OFF(AB)+D2ADD(CD)
crsg        D2INT(ABCD)=RINT
crsg     ENDIF
crsg non-canonical addressing scheme of fz-act-fz-act integrals
crsg in symcnt requires this additional check
c  THIS BLOCK OF CODE WORKS BUT HAS MANY UNNEESSARY CHECKS
c        else if(fzt(b).eq.1.and.fzt(d).eq.1) then
c  I think the second half of this test i always false (c>a)
c           if(b.lt.d.and.a.gt.c.and.a.gt.d.or
c    1                 .d.lt.b.and.c.gt.a.and.c.gt.b) then
c               abcd = d2off(cd)+d2add(ab)
c           else
c               abcd = d2off(ab)+d2add(cd)
c           endif
c           d2int(abcd)=rint
c        else if(fzt(a).eq.1.and.fzt(d).eq.1) then
c  I think this test is always false (a<d or c>a)
c           if(a.lt.d.and.b.gt.c.and.b.gt.d.or
c    1                 .d.lt.a.and.c.gt.a.and.c.gt.b) then
c               abcd = d2off(cd)+d2add(ab)
c           else
c               abcd = d2off(ab)+d2add(cd)
c           endif
c           d2int(abcd)=rint
c        else if(fzt(b).eq.1.and.fzt(c).eq.1) then
c  I think the second half of this test is always false (d>a)
c           if(b.lt.c.and.a.gt.d.and.a.gt.c.or
c    1                 .c.lt.b.and.d.gt.a.and.d.gt.b) then
c               abcd = d2off(cd)+d2add(ab)
c           else
c               abcd = d2off(ab)+d2add(cd)
c           endif
c           d2int(abcd)=rint
c        else if(fzt(a).eq.1.and.fzt(c).eq.1) then
c  I think this test is always false (a<c or d>a)
c           if(a.lt.c.and.b.gt.c.and.b.gt.d.or
c    1                 .c.lt.a.and.d.gt.a.and.d.gt.b) then
c               abcd = d2off(cd)+d2add(ab)
c           else
c               abcd = d2off(ab)+d2add(cd)
c           endif
c           d2int(abcd)=rint
c        WRITE (JOUT,617) FREEZE,A,B,C,D,ABCD,RINT
c        endif
         else if(fzt(b).eq.1.and.fzt(d).eq.1) then
            if(b.lt.d.and.a.gt.c.and.a.gt.d) then
                abcd = d2off(cd)+d2add(ab)
            else
                abcd = d2off(ab)+d2add(cd)
            endif
            d2int(abcd)=rint
         else if(fzt(b).eq.1.and.fzt(c).eq.1) then
            if(b.lt.c.and.a.gt.d.and.a.gt.c) then
                abcd = d2off(cd)+d2add(ab)
            else
                abcd = d2off(ab)+d2add(cd)
            endif
            d2int(abcd)=rint
         else 
            abcd = d2off(ab)+d2add(cd)
            d2int(abcd)=rint
C        WRITE (JOUT,617) FREEZE,A,B,C,D,ABCD,RINT
         endif
crsg end changes
      ELSE IF(ILOOP.EQ.LOOP.AND.FREEZE.EQ.3.AND.ITAPIN.EQ.ITAP78)THEN
         IF(FZT(A).EQ.0.OR.FZT(B).EQ.0)THEN
            ABCD=E3OFF(AB)+E3ADD(CD)
         ELSE
            ABCD=E3OFF(CD)+E3ADD(AB)
         ENDIF
            E3INT(ABCD)=RINT
C        WRITE (JOUT,616) FREEZE,A,B,C,D,ABCD,RINT
      ENDIF
  616 FORMAT(' E3INT =',I4,4(2X,I4),I6,5X,F12.8)
  617 FORMAT(' D2INT =',I4,4(2X,I4),I6,5X,F12.8)
  618 FORMAT(' C2INT =',I4,4(2X,I4),I6,5X,F12.8)
  619 FORMAT(' FREEZE=',I4,4(2X,I4),I6,5X,F12.8)
C
  101 CONTINUE
      IF (IFLG.EQ.0) GO TO 111
C
      TOFF=0
      XPQ=0
      DO 390 TSYM=1,NIRRED
      DO 380 P=1,NT
      IF(FZT(P).EQ.1)GO TO 380
      PSYM=ORBSYM(P)
      DO 370 Q=1,P
      IF(FZT(Q).EQ.1)GO TO 370
      QSYM=ORBSYM(Q)
      PQSYM=IEOR(PSYM,QSYM)
      IF(PQSYM.NE.TSYM-1)GO TO 370
      XPQ=XPQ+1
      IF(XPQ.LT.FFPAIR.OR.XPQ.GT.LLPAIR)GO TO 365
      NADD=0
      DO 350 R=1,NT
      IF(FZT(R).EQ.1)GO TO 350
      RSYM=ORBSYM(R)
      SSYM=IEOR(RSYM,PQSYM)+1
      FS=FLOV(SSYM,1)
      LS=FLOV(SSYM,2)
      IF(LS.GT.R)LS=R
      DO 340 S=FS,LS
      NADD=NADD+1
      TOFF=TOFF+1
      GBUCK(NADD)=CORE(TOFF)
  340 CONTINUE
  350 CONTINUE
      LRS=NRS(TSYM)
C     WRITE(6,651)XPQ,P,Q,NRS(TSYM)
  651 FORMAT(' WRITING PAIR',I4,' P=',I3,' P=',I3,' OF LENGTH ',I6)
      CALL WWRITW(ITAPOU,GBUCK,INTOWP(LRS),POINT,POINT)
  365 CONTINUE
  370 CONTINUE
  380 CONTINUE
  390 CONTINUE
C
 1000 CONTINUE
      CALL RCLOSE(ITAPOU,3)
C
C >>> READ ONE-ELECTRON MATRIX
C
 1111 CONTINUE
      CALL SREAD(ITAPIN,IBUF,INTOWP(LENGTH))
      IFLG = IBUF(1)
      MBUF = IBUF(2)
      WRITE (6,*) ' IFLG1 =',IFLG
      WRITE (6,*) ' MBUF1 =',MBUF
      DO 1011 II = 1,MBUF
         IJKL   = IBUF(2+II)
         I      = ISHFT(IJKL,-8)
         L      = IAND(ITEMP,IJKL)
         VAL    = BUF(IVOFF+II)
C
         IF(ITAPIN.EQ.ITAP68)VAL=VAL*0.5D0
         I      = CCTSCF(I)
         L      = CCTSCF(L)
C
COUT     WRITE (6,*) 'I,L, ONEM(I,L)= ',I,L,VAL
C
         IF(I.LT.L) THEN
            ISWOP = L
            L     = I
            I     = ISWOP
         END IF
         IL = IPQ(I)+L
         XAR1(IL) =  VAL
 1011 CONTINUE
      IF ( IFLG .EQ. 0 ) GOTO 1111
C
C
C     READ ITAPOU TO CHECK
C
C     CALL RFILE(ITAPOU)
C     POINT=1
C     DO 590 TSYM=1,NIRRED
C     DO 580 P=1,NT
C     IF(FZT(P).EQ.1)GO TO 580
C     PSYM=ORBSYM(P)
C     DO 570 Q=1,P
C     IF(FZT(Q).EQ.1)GO TO 570
C     QSYM=ORBSYM(Q)
C     PQSYM=IEOR(PSYM,QSYM)
C     IF(PQSYM.NE.TSYM-1)GO TO 570
C     LRS=NRS(TSYM)
C     CALL WREADW(ITAPOU,GBUCK,INTOWP(LRS),POINT,POINT)
C     NADD=0
C     DO 550 R=1,NT
C     IF(FZT(R).EQ.1)GO TO 550
C     RSYM=ORBSYM(R)
C     SSYM=IEOR(RSYM,PQSYM)+1
C     FS=FLOV(SSYM,1)
C     LS=FLOV(SSYM,2)
C     IF(LS.GT.R)LS=R
C     DO 540 S=FS,LS
C     NADD=NADD+1
C     WRITE(6,656)P  ,Q  ,R  ,S  ,GBUCK(NADD)
C 656 FORMAT (4I3,F20.12,3X,'OU')
C 540 CONTINUE
C 550 CONTINUE
C 570 CONTINUE
C 580 CONTINUE
C 590 CONTINUE
C     CALL RCLOSE(ITAPOU,3)
C
      RETURN
      END
