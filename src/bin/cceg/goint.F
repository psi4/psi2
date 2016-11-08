C
C ----------------------------------------------------------------------
C
      SUBROUTINE GOINT (T1O,T2O,EIN,FOCK,HO,GO,
     .                  ORBSYM,FLOV,EOFF,EADD,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL
      REAL*8  T1O(NO,NV),T2O(*),EIN(*),FOCK(*),HO(NO,NO),GO(NO,NO)
      INTEGER ORBSYM(*),FLOV(NIRRED,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO,2),TADD(NV,NV),ZLX(NV,NV),
     .        EOFF(*),EADD(*)
C
      CALL ZERO(GO,NO*NO)
C
      DO 190 U=1,NO
         IF(FZO(U).EQ.1)GO TO 190
         USYM=ORBSYM(U)
         DO 180 I=1,NO
            IF(FZO(I).EQ.1)GO TO 180
            ISYM=ORBSYM(I)
            IF(USYM.NE.ISYM)GO TO 180
            VAL=0.0D0
            IU=ITR(MAX0(I,U))+MIN0(I,U)
            VAL=VAL+HO(U,I)
            FA=FLOV(ISYM+1,3)-NO
            LA=FLOV(ISYM+1,4)-NO
            DO 120 A=FA,LA
               FIA=IPQ(A+NO)+I
               VAL=VAL+FOCK(FIA)*T1O(U,A)
  120       CONTINUE
            DO 170 A=1,NV
               IF(FZV(A).EQ.1)GO TO 170
               ASYM=ORBSYM(A+NO)
               IA=ITR(A)+I
               FJ=FLOV(ASYM+1,1)
               LJ=FLOV(ASYM+1,2)
               DO 160 J=FJ,LJ
                  JA=ITR(A)+J
                  IUJA=EOFF(JA)+EADD(IU)
                  JU=ITR(MAX0(J,U))+MIN0(J,U)
                  JUIA=EOFF(IA)+EADD(JU)
                  VAL=VAL+(EIN(IUJA)+EIN(IUJA)-EIN(JUIA))*T1O(J,A)
  160          CONTINUE
  170       CONTINUE
            GO(U,I)=VAL
  180    CONTINUE
  190 CONTINUE
C
      RETURN
      END
