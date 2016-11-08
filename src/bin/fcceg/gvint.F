C
C ----------------------------------------------------------------------
C
      SUBROUTINE GVINT (T1O,T2O,FIN,FOCK,HV,GV,
     .                  ORBSYM,FLOV,FOFF,FADD,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL
      REAL*8  T1O(NO,NV),T2O(*),FIN(*),FOCK(*),HV(NV,NV),GV(NV,NV)
      INTEGER ORBSYM(*),FLOV(NIRRED,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO,2),TADD(NV,NV),ZLX(NV,NV),
     .        FOFF(*),FADD(*)
C
      CALL ZERO(GV,NV*NV)
C
      DO 190 A=1,NV
         IF(FZV(A).EQ.1)GO TO 190
         ASYM=ORBSYM(A+NO)
         DO 180 BE=1,NV
            IF(FZV(BE).EQ.1)GO TO 180
            BESYM=ORBSYM(BE+NO)
            IF(BESYM.NE.ASYM)GO TO 180
            VAL=0.0D0
            BEA=ITV(MAX0(A,BE))+MIN0(A,BE)
            VAL=VAL+HV(BE,A)
            FI=FLOV(ASYM+1,1)
            LI=FLOV(ASYM+1,2)
            DO 110 I=FI,LI
               FIA=IPQ(A+NO)+I
               VAL=VAL-FOCK(FIA)*T1O(I,BE)
110         CONTINUE
            DO 170 I=1,NO
               IF(FZO(I).EQ.1)GO TO 170
               ISYM=ORBSYM(I)
               IA=ITR(A)+I
               FB=FLOV(ISYM+1,3)-NO
               LB=FLOV(ISYM+1,4)-NO
               DO 160 B=FB,LB
                  IB=ITR(B)+I
                  BEAIB=FOFF(BEA)+FADD(IB)
                  BEB=ITV(MAX0(B,BE))+MIN0(B,BE)
                  BEBIA=FOFF(BEB)+FADD(IA)
                  VAL=VAL+(FIN(BEAIB)+FIN(BEAIB)-FIN(BEBIA))*T1O(I,B)
  160          CONTINUE
  170       CONTINUE
            GV(BE,A)=VAL
  180    CONTINUE
  190 CONTINUE
C
      RETURN
      END
