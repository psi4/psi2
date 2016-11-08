C
C ----------------------------------------------------------------------
C
      SUBROUTINE HT1   (T1,T1O,T2O,CIN,DIN,FOCK,HO,HV,HIA,
     .                  ORBSYM,FLOV,COFF,CADD,DOFF,DADD,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL,VAL1,VAL2,FVAL
      REAL*8  T1(NO,NV),T1O(NO,NV),T2O(*),CIN(*),DIN(*),FOCK(*),
     .        HO(NO,NO),HV(NV,NV),HIA(NO,NV)
      INTEGER ORBSYM(*),FLOV(NIRRED,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO,2),TADD(NV,NV),ZLX(NV,NV),
     .        COFF(*),CADD(*),DOFF(*),DADD(*)
C
      DO 191 U=1,NO
         IF(FZO(U).EQ.1)GO TO 191
         USYM=ORBSYM(U)
         FBE=FLOV(USYM+1,3)-NO
         LBE=FLOV(USYM+1,4)-NO
         FI =FLOV(USYM+1,1)
         LI =FLOV(USYM+1,2)
         DO 171 BE=FBE,LBE
            BEU=ITR(BE)+U
            FBEU=IPQ(BE+NO)+U
            VAL=0.0D0
            DO 135 I=FI,LI
               VAL=VAL-HO(U,I)*T1O(I,BE)
  135       CONTINUE
            DO 138 A=FBE,LBE
               VAL=VAL+HV(BE,A)*T1O(U,A)
  138       CONTINUE
            DO 148 I=1,NO
               IF(FZO(I).EQ.1)GO TO 148
               IU=ITR(MAX0(I,U))+MIN0(I,U)
               ISYM=ORBSYM(I)
               FA=FLOV(ISYM+1,3)-NO
               LA=FLOV(ISYM+1,4)-NO
               DO 145 A=FA,LA
                  ZL=ZLX(A,BE)
                  IA=ITR(A)+I
                  BEUIA=DOFF(MAX0(BEU,IA))+DADD(MIN0(BEU,IA))
                  BEA=ITV(MAX0(BE,A))+MIN0(BE,A)
                  IUBEA=COFF(BEA)+CADD(IU)
                  VAL=VAL+(DIN(BEUIA)+DIN(BEUIA)-CIN(IUBEA))*T1O(I,A)
                  IUABE=TOFF(I,U,ZL)+TADD(A,BE)
                  UIABE=TOFF(U,I,ZL)+TADD(A,BE)
                  VAL1=T2O(IUABE)
                  VAL2=T2O(UIABE)-T1O(U,A)*T1O(I,BE)
                  VAL=VAL+HIA(I,A)*(VAL1+VAL1-VAL2)
                  FIA=IPQ(A+NO)+I
                  FVAL=FOCK(FIA)*T1O(U,A)*T1O(I,BE)
                  VAL=VAL-FVAL-FVAL
  145          CONTINUE
  148       CONTINUE
            VAL=VAL+FOCK(FBEU)
            T1(U,BE)=T1(U,BE)+VAL
  171    CONTINUE
  191 CONTINUE
C
      RETURN
      END
