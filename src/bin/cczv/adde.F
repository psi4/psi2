C
C-----------------------------------------------------------------------
C
      SUBROUTINE ADDE(EMAT,I1MAT,JMAT,T1,T2,IOFF,NO,NV,NTRO,NTRV,
     .                NO2,NV2,NTR,JOUT,INO,INTRO,INTRV,NOV,INV,N2OV,ZLX,
     .                ZLXO,UOFF,FLOV,ORBSYM,NIRRED,
     .                TOFF,TADD,NDIMT2,EMOFF,EMADD,LEMAT,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      REAL*8 T2(NDIMT2),T1(NO,NV),JMAT(NV*NO),EMAT(LEMAT),
     .       I1MAT(NO,NV)
      INTEGER IOFF(NTR),EMOFF(NO*NV),EMADD(NO*NO),FZO(*),FZV(*)
      INTEGER INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV)
      INTEGER TOFF(NO,NO,2),TADD(NV,NV)
      INTEGER ZLX(NV,NV),UOFF(NO,NO,2),ZLXO(NO,NO)
      DIMENSION FLOV(NIRRED,4),ORBSYM(NO+NV)
      COMMON/TAPES/TPBMAT,TPFMAT,TPHMAT,TPRMAT,TPQMAT,TPJMAT,ITAP71,
     .             TPAMAT
C
      CALL RFILE(TPJMAT)
      POINT=1
      DO 10 U=1,NO
         IF(FZO(U).EQ.1)GO TO 10
         USYM=ORBSYM(U)
         DO 720 IRI=1,NIRRED
         ISYM=IRI-1
         UISYM=IEOR(USYM,ISYM)
         FI=FLOV(IRI,1)
         LI=FLOV(IRI,2)
         DO 20 I=FI,LI
            CALL WREADW(TPJMAT,IDIM,1,POINT,POINT)
            CALL WREADW(TPJMAT,JMAT,INTOWP(IDIM),POINT,POINT)
            ICNTJ=0
            UIX=INO(U)+I
            DO 730 IRJ=1,NIRRED
               JSYM=IRJ-1
               FJ=FLOV(IRJ,1)
               LJ=FLOV(IRJ,2)
               ASYM=IEOR(UISYM,JSYM)+1
               FA=FLOV(ASYM,3)-NO
               LA=FLOV(ASYM,4)-NO
               FK=FLOV(ASYM,1)
               LK=FLOV(ASYM,2)
C
            DO 30 J=FJ,LJ
               DO 40 A=FA,LA
                  FB=FLOV(USYM+1,3)-NO
                  LB=FLOV(USYM+1,4)-NO
                  JAX=INO(A)+J
                  UIJA=EMOFF(JAX)+EMADD(UIX)
       IF(UIJA.GT.LEMAT)WRITE(6,*)'ERROR UIJA ADDE'
                  DO 70 B=FB,LB
                     ZL=ZLX(A,B)
                     JIAB=TOFF(J,I,ZL)+TADD(A,B)
                     EMAT(UIJA)=EMAT(UIJA)+T2(JIAB)*I1MAT(U,B)
   70             CONTINUE
   40          CONTINUE
   30       CONTINUE
  730       CONTINUE
C
            DO 90 K=1,NO
            IF(FZO(K).EQ.1)GO TO 90
            KSYM=ORBSYM(K)
            BSYM=IEOR(UISYM,KSYM)+1
            FB=FLOV(BSYM,3)-NO
            LB=FLOV(BSYM,4)-NO
               DO 100 B=FB,LB
      ICNTJ=ICNTJ+1
                  DO 110 J=1,NO
                  IF(FZO(J).EQ.1)GO TO 110
                  JSYM=ORBSYM(J)
                  ASYM=IEOR(UISYM,JSYM)+1
                  FA=FLOV(ASYM,3)-NO
                  LA=FLOV(ASYM,4)-NO
                     DO 120 A=FA,LA
                     ZL=ZLX(A,B)
                     JKAB=TOFF(J,K,ZL)+TADD(A,B)
CJ                 EMAT(U,I,J,A)=EMAT(U,I,J,A)+JMAT(B,K,U,I)*T2(JKAB)
                   JAX=INO(A)+J
                   UIJA=EMOFF(JAX)+EMADD(UIX)
                   EMAT(UIJA)=EMAT(UIJA)+JMAT(ICNTJ)*T2(JKAB)
  120          CONTINUE
  110       CONTINUE
  100                CONTINUE
   90             CONTINUE
   20    CONTINUE
  720    CONTINUE
   10 CONTINUE
      CALL RCLOSE(TPJMAT,3)
      RETURN
      END
