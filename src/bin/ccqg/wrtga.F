C
C----------------------------------------------------------------------
C
      SUBROUTINE WRTGA(BUF,IBUF,GBUF,NBF,NO,NV,ITAP68,GMAT,NDIMG,
     .                 IOFF,NTR,JOUT,ICNT,AOFF,AADD,ITR,
     .                 FLOV,ORBSYM,NIRRED,FZO,FZV)
      IMPLICIT INTEGER (A-Z)
      REAL*8 BUF(GBUF),GMAT(NDIMG),TOL
      INTEGER IBUF(GBUF*2),IOFF(NTR),FZO(NO),FZV(NV),
     .        ITR(NO),AOFF(NO*NO),AADD(NO*NO),ORBSYM(NBF),FLOV(NIRRED,4)
      COMMON/BUFFER/TOL,IBFLEN,MAXBUF,IBOFF
C
      DO 223 I=1,NO
         IF(FZO(I).EQ.1)GO TO 223
         ISYM=ORBSYM(I)
         DO 222 J=1,I
            IF(FZO(J).EQ.1)GO TO 222
            JSYM=ORBSYM(J)
            IJSYM=IEOR(ISYM,JSYM)
            IJ=ITR(I)+J
            DO 221 K=1,I
               IF(FZO(K).EQ.1)GO TO 221
               KSYM=ORBSYM(K)
               LSYM=IEOR(IJSYM,KSYM)+1
               FL=FLOV(LSYM,1)
               TL=FLOV(LSYM,2)
               TL2=TL
               IF(TL.GT.K)TL2=K
               IF(K.EQ.I)TL2=J
               DO 220 L=FL,TL2
                  KL=ITR(K)+L
                  IJKL=AOFF(IJ)+AADD(KL)
c                 WRITE (JOUT,1) I,J,K,L,GMAT(IJKL)/8.0D+00
                  IF(DABS(GMAT(IJKL)).GT.TOL) THEN
                     LBL=IOR(J,ISHFT(I,8))
                     LBL=IOR(K,ISHFT(LBL,8))
                     LBL=IOR(L,ISHFT(LBL,8))
                     IF(ICNT.GT.MAXBUF) THEN
                        CALL SWRIT(ITAP68,IBUF,INTOWP(GBUF))
                        CALL ZERO(BUF,GBUF)
                        ICNT=1
                     ENDIF
                     IBUF(1)=0
                     IBUF(2)=ICNT
                     IBUF(2+ICNT)=LBL
                     BUF(ICNT+IBOFF)=GMAT(IJKL)
                     ICNT=ICNT+1
                  ENDIF
  220          CONTINUE
  221       CONTINUE
  222    CONTINUE
  223 CONTINUE
C
    1 FORMAT ('GMAT(',I3,',',I3,',',I3,',',I3,')=',F20.12)
      RETURN
      END
