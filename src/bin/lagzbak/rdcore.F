C
C-----------------------------------------------------------------------
C
      SUBROUTINE RDCORE(ITAPNN,XAR2,LENAR2,JOUT,BUF,IBUF,BUFSZ,
     .                  IPQ,NT2,CCTSCF,NO,NT,XAR1,IOFF,IADD,FZT,
     .                  I1OFF,I1ADD,X1INT,C2OFF,C2ADD,C2INT,
     .                  D2OFF,D2ADD,D2INT,E3OFF,E3ADD,E3INT)
      IMPLICIT INTEGER (A-Z)
      REAL*8 VAL,FAC1,FAC2
      REAL*8 XAR2(LENAR2),XAR1(NT2),BUF(BUFSZ),
     .       X1INT(*),C2INT(*),D2INT(*),E3INT(*)
      INTEGER IBUF(*),IPQ(*),CCTSCF(*),IOFF(*),IADD(*),FZT(*),
     .        I1OFF(*),C2OFF(*),I1ADD(*),C2ADD(*),D2OFF(*),D2ADD(*),
     .        E3OFF(*),E3ADD(*)
C
      PARAMETER ( ITEMP = 255, FAC1 = 0.5D0, FAC2 = 0.125D0 )
      LOGICAL    ALLOCC
C
      IND(I,J) = I*(I-1)/2 + J
C
      ITAP68=68
      ITAP78=78
      IBFLEN = (INTOWP(BUFSZ)-2) / INTOWP(1)
      MAXBUF = INTOWP(IBFLEN)    / (1+INTOWP(1))
      IBOFF  = (MAXBUF+3)        / INTOWP(1)
C
C     WRITE(6,*)' '
C     WRITE(6,*)'IN RDCORE READING TAPE',ITAPNN
C     WRITE(6,*)' '
C
C     WRITE (6,*) ' MAXBUF=',MAXBUF
C     WRITE (6,*) ' IBFLEN=',IBFLEN
C     WRITE (6,*) ' IBOFF =',IBOFF
C
C >>> TWO-ELECTRON MATRIX
C
      CALL SREW(ITAPNN)
      IF(ITAPNN.EQ.ITAP78)CALL SREAD(ITAP78,JUNK,1)
 1111 CONTINUE
      CALL SREAD(ITAPNN,IBUF,INTOWP(BUFSZ))
      IFLG = IBUF(1)
      MBUF = IBUF(2)
C     WRITE (6,*) ' IFLG =',IFLG
C     WRITE (6,*) ' MBUF =',MBUF
      DO 1101 II = 1,MBUF
         IJKL   = IBUF(2+II)
         I1     = ISHFT(IJKL,-24)
         J1     = IAND(ITEMP,ISHFT(IJKL,-16))
         K1     = IAND(ITEMP,ISHFT(IJKL,-8))
         L1     = IAND(ITEMP,IJKL)
         VAL    = BUF(IBOFF+II)
          I      = CCTSCF(I1)
          J      = CCTSCF(J1)
          K      = CCTSCF(K1)
          L      = CCTSCF(L1)
C
         FREEZE=FZT(I)+FZT(J)+FZT(K)+FZT(L)
C
C >>> PUT IN CANONICAL ORDER
C
      IF (I.LT.J) THEN
         ISWOP=J
         J=I
         I=ISWOP
      END IF
      IF (K.LT.L) THEN
         ISWOP=L
         L=K
         K=ISWOP
      END IF
      IF (IND(I,J).LT.IND(K,L)) THEN
         ISWOP=K
         K=I
         I=ISWOP
         ISWOP=L
         L=J
         J=ISWOP
      END IF
C
         IF(ITAPNN.EQ.ITAP68)THEN
            VAL=VAL*FAC2
         ENDIF
C
COUT     WRITE (6,613) I ,J ,K ,L ,VAL
         IJ   = IND(I,J)
         KL   = IND(K,L)
C
      IF(FREEZE.EQ.0)THEN
         IJKL = IOFF(IJ)+IADD(KL)
         XAR2(IJKL) =  VAL
C        WRITE (6,619) FREEZE,I ,J ,K ,L ,IJKL,VAL
      ELSE IF(FREEZE.EQ.1.AND.ITAPNN.EQ.ITAP78) THEN
         IF(FZT(I).EQ.1.OR.FZT(J).EQ.1)THEN
            IJKL = I1OFF(IJ)+I1ADD(KL)
         ELSE
            IJKL = I1OFF(KL)+I1ADD(IJ)
         ENDIF
         X1INT(IJKL)=VAL
C        WRITE (6,619) FREEZE,I ,J ,K ,L ,IJKL,VAL
      ELSE IF(FREEZE.EQ.2.AND.ITAPNN.EQ.ITAP78) THEN
         IF(FZT(I).EQ.1.AND.FZT(J).EQ.1) THEN
            IJKL = C2OFF(IJ)+C2ADD(KL)
            C2INT(IJKL)=VAL
         ELSE IF(FZT(K).EQ.1.AND.FZT(L).EQ.1) THEN
            IJKL = C2OFF(KL)+C2ADD(IJ)
            C2INT(IJKL)=VAL
crsg     ELSE
crsg        IJKL = D2OFF(IJ)+D2ADD(KL)
crsg        D2INT(IJKL)=VAL
crsg     ENDIF
crsg non-canonical addressing scheme of fz-act-fz-act integrals
crsg in symcnt requires replacing the simple check above as follows
c   NOTE: THIS BLOCK OF CODE WORKS, BUT MAY HAVE UNNECESSARY CHECKS
c        else if(fzt(j).eq.1.and.fzt(l).eq.1) then
c   I think the second half of this test is always false (k.gt.i)
c           if(j.lt.l.and.i.gt.k.and.i.gt.l.or
c    1                 .l.lt.j.and.k.gt.i.and.k.gt.j) then
c               ijkl = d2off(kl)+d2add(ij)
c           else
c               ijkl = d2off(ij)+d2add(kl)
c           endif
c           d2int(ijkl)=val
c        else if(fzt(i).eq.1.and.fzt(l).eq.1) then
c   I think the test below is always false (i<l or k>i)
c           if(i.lt.l.and.j.gt.k.and.j.gt.l.or
c    1                 .l.lt.i.and.k.gt.i.and.k.gt.j) then
c               ijkl = d2off(kl)+d2add(ij)
c           else
c               ijkl = d2off(ij)+d2add(kl)
c           endif
c           d2int(ijkl)=val
c        else if(fzt(j).eq.1.and.fzt(k).eq.1) then
c   I think the second half of the test below is always false (l>i)
c           if(j.lt.k.and.i.gt.l.and.i.gt.k.or
c    1                 .k.lt.j.and.l.gt.i.and.l.gt.j) then
c               ijkl = d2off(kl)+d2add(ij)
c           else
c               ijkl = d2off(ij)+d2add(kl)
c           endif
c           d2int(ijkl)=val
c        else if(fzt(i).eq.1.and.fzt(k).eq.1) then
c   I think the test below is always false (i<k or l>i)
c           if(i.lt.k.and.j.gt.l.and.j.gt.k.or
c    1                 .k.lt.i.and.l.gt.i.and.l.gt.j) then
c               ijkl = d2off(kl)+d2add(ij)
c           else
c               ijkl = d2off(ij)+d2add(kl)
c           endif
c           d2int(ijkl)=val
c        endif
         else if(fzt(j).eq.1.and.fzt(l).eq.1) then
            if(j.lt.l.and.i.gt.k.and.i.gt.l) then
                ijkl = d2off(kl)+d2add(ij)
            else
                ijkl = d2off(ij)+d2add(kl)
            endif
            d2int(ijkl)=val
         else if(fzt(j).eq.1.and.fzt(k).eq.1) then
            if(j.lt.k.and.i.gt.l.and.i.gt.k) then
                ijkl = d2off(kl)+d2add(ij)
            else
                ijkl = d2off(ij)+d2add(kl)
            endif
            d2int(ijkl)=val
         else
            ijkl = d2off(ij)+d2add(kl)
            d2int(ijkl)=val
         endif
C        WRITE (6,619) FREEZE,I ,J ,K ,L ,IJKL,VAL
      ELSE IF(FREEZE.EQ.3.AND.ITAPNN.EQ.ITAP78)THEN
         IF(FZT(I).EQ.0.OR.FZT(J).EQ.0)THEN
            IJKL=E3OFF(IJ)+E3ADD(KL)
         ELSE
            IJKL=E3OFF(KL)+E3ADD(IJ)
         ENDIF
            E3INT(IJKL)=VAL
C        WRITE (6,619) FREEZE,I ,J ,K ,L ,IJKL,VAL
      ENDIF
  619 FORMAT(' FREEZE=',I4,4(2X,I3),5X,I6,3X,F12.8)
C
 1101 CONTINUE
      IF (IFLG.EQ.0) GOTO 1111
  613 FORMAT(4I4,I6,F20.12)
C
C >>> ONE-ELECTRON MATRIX
C
  111 CONTINUE
      CALL SREAD(ITAPNN,IBUF,INTOWP(BUFSZ))
      IFLG = IBUF(1)
      MBUF = IBUF(2)
C     WRITE (6,*) ' IFLG =',IFLG
C     WRITE (6,*) ' MBUF =',MBUF
      DO 101 II = 1,MBUF
         IJKL   = IBUF(2+II)
         I      = ISHFT(IJKL,-8)
         L      = IAND(ITEMP,IJKL)
         VAL    = BUF(IBOFF+II)
C
C        WRITE (6,*) 'UNNORMAL CC I,L,1M(I,L)= ',I,L,VAL
C
         IF(ITAPNN.EQ.ITAP68)VAL=VAL*FAC1
         I      = CCTSCF(I)
         L      = CCTSCF(L)
C
C        WRITE (6,*) 'CANONICAL SCF I,L,1M(I,L)= ',I,L,VAL
C
         IF(I.LT.L) THEN
            ISWOP = L
            L     = I
            I     = ISWOP
         END IF
         IL = IND(I,L)
         XAR1(IL) =  VAL
  101 CONTINUE
      IF ( IFLG .EQ. 0 ) GOTO 111
C
      RETURN
      END
