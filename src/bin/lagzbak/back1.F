C
C-----------------------------------------------------------------------
C
      SUBROUTINE BACK1 (QDMT,ZETA,XLAG,VEC,IMAT,DMAT,ITAP79,
     .                  ORBSYM,NIRRED,FLOV,IPQ,NO,NV,NT,NT2,
     .                  QDMAO,ITR,CCTSCF,NTAO,NTAO2)
      IMPLICIT INTEGER (A-Z)
      REAL*8 VEC(NTAO,NT),QDMT(*),IMAT(NTAO,NT),DMAT(NTAO,NTAO),
     .       XLAG(NT,NT),ZETA(*),QDMAO(NTAO2),VAL,VAL1,VAL2,SUM
      INTEGER ORBSYM(*),FLOV(NIRRED,2),IPQ(*),ITR(*),CCTSCF(*)
C
C >>> WRITE Q,X (MO BASIS) TO TAPE79
C
      START=1+INTOWP(NT2+NT2+1)
      I79=START
      CALL WWRITW(ITAP79,QDMT,INTOWP(NT2),I79,I79)
C
c     WRITE(6,*)
c     WRITE(6,*)'TOTAL ONEPDM IN mO BASIS in  back1 as written to 79'
c     CALL PRINT(QDMt,NT2,NT,6)
C
      CALL WWRITW(ITAP79,XLAG,INTOWP(NT*NT),I79,I79)
C
c     WRITE(6,*)
c     WRITE(6,*)'TOTAL lagran IN mO BASIS in  back1 as written to 79'
c     CALL matout(xlag,NT,nt,nt,nt,6)
C
C
C >>> BACK TRANSFORM ONE-PARTICLE OPERATORS
C
      CALL ZERO (IMAT,NTAO*NT)
      DO 330 I=1,NT
         DO 320 J=1,NT
            JI=IPQ(MAX0(I,J))+MIN0(I,J)
            DO 310 S=1,NTAO
               IMAT(S,I)=IMAT(S,I)+VEC(S,J)*QDMT(JI)
  310       CONTINUE
  320    CONTINUE
  330 CONTINUE
C
      CALL ZERO (DMAT,NTAO*NTAO)
      DO 430 S=1,NTAO
         DO 420 R=1,NTAO
            DO 410 I=1,NT
               DMAT(S,R)=DMAT(S,R)+VEC(R,I)*IMAT(S,I)
  410       CONTINUE
  420    CONTINUE
  430 CONTINUE
C
      IJ=0
      DO 780 I=1,NTAO
      DO 770 J=1,I
      IJ=IJ+1
      QDMAO(IJ)=DMAT(I,J)
  770 CONTINUE
  780 CONTINUE
C
c     WRITE(6,*)
c     WRITE(6,*)' TOTAL ONEPDM IN AO BASIS '
c     CALL PRINT(QDMAO,NTAO2,NTAO,6)
C
      CALL WWRITW(ITAP79,QDMAO,INTOWP(NTAO2),I79,I79)
C
      CALL ZERO (IMAT,NTAO*NT)
      DO 530 I=1,NT
         DO 520 J=1,NT
            DO 510 S=1,NTAO
               IMAT(S,I)=IMAT(S,I)+VEC(S,J)*XLAG(J,I)
  510       CONTINUE
  520    CONTINUE
  530 CONTINUE
C
      CALL ZERO (DMAT,NTAO*NTAO)
      DO 630 S=1,NTAO
         DO 620 R=1,NTAO
            DO 610 I=1,NT
               DMAT(S,R)=DMAT(S,R)+VEC(R,I)*IMAT(S,I)
  610       CONTINUE
  620    CONTINUE
  630 CONTINUE
C
      DO 730 I=1,NTAO
         DO 720 J=1,I
             IJ=IPQ(I)+J
             QDMAO(IJ)=DMAT(I,J)
  720    CONTINUE
  730 CONTINUE
C
c     WRITE(6,*)
c     WRITE(6,*)' TOTAL LAGRANGIAN IN AO BASIS'
c     CALL PRINT(QDMAO,NTAO2,NTAO,6)
C
      CALL WWRITW(ITAP79,QDMAO,INTOWP(NTAO2),I79,I79)
      CALL WWRITW(ITAP79,ZETA,INTOWP(NO*NV),I79,I79)
C
      itap80=80
      call rfile(itap80)
      call wwritw(itap80,zeta,intowp(no*nv),1,iend)
      call rclose(itap80,3)
c     write(6,*)
c     write(6,*)'z vector as written to file80'
c     call matout(zeta,no,nv,no,nv,6)
c
      RETURN
      END
