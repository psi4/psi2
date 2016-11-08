      SUBROUTINE fterm (T1O,T2O,T2,FIN,AR3,AR4,
     .                  SUX1,SUX2,SUX3,SUX4,
     .                  ORBSYM,FLOV,nFOFF,nFADD,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED,MAXOO,MAXVV,nfints)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL1,VAL2,TMP,TMV
      REAL*8  T1O(NO,NV),T2O(*),T2(*),FIN(*),
     .        SUX1(*),SUX2(*),SUX3(*),SUX4(*),
     .        AR3(MAXVV,MAXOO),AR4(MAXVV,MAXOO)
      INTEGER ORBSYM(*),FLOV(NIRRED,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO,2),TADD(NV,NV),ZLX(NV,NV),
     .        nFOFF(*),nFADD(*),nfints(*)
C
      itap75=75
      CALL RFILE(ITAP75)
C
      DO 992 TSYM=1,NIRRED
         XAB=0
         DO 290 A=1,NV
            IF(FZV(A).EQ.1)GO TO 290
            ASYM=ORBSYM(A+NO)
            DO 280 B=1,A
               IF(FZV(B).EQ.1)GO TO 280
               BSYM=ORBSYM(B+NO)
               BASYM=IEOR(BSYM,ASYM)
               IF(BASYM.NE.TSYM-1)GO TO 280
               XAB=XAB+1
               XIJ=0
               DO 270 I=1,NO
                  IF(FZO(I).EQ.1)GO TO 270
                  ISYM=ORBSYM(I)
                  JSYM=IEOR(ISYM,BASYM)+1
                  FJ=FLOV(JSYM,1)
                  LJ=FLOV(JSYM,2)
                  IF(LJ.GT.I)LJ=I
                  DO 260 J=FJ,LJ
                     XIJ=XIJ+1
                     IJAB=TOFF(I,J,1)+TADD(A,B)
                     JIAB=TOFF(J,I,1)+TADD(A,B)
                     VAL1=T2O(IJAB)+T1O(I,A)*T1O(J,B)
                     VAL2=T2O(JIAB)+T1O(I,B)*T1O(J,A)
                     AR3(XAB,XIJ)=VAL1+VAL2
                     AR4(XAB,XIJ)=VAL1-VAL2
  260             CONTINUE
  270          CONTINUE
  280       CONTINUE
  290    CONTINUE
         LIMO=XIJ
         LIMV=XAB
C
         pon75=1
         shift=0
         DO 991 BE=1,NV
            IF(FZV(BE).EQ.1)GO TO 991
            dim75=nfints(be)
            BESYM=ORBSYM(BE+NO)
            CALL WREADW(ITAP75,FIN,INTOWP(dim75),PON75,PON75)
            DO 981 i=1,no
               IF(FZo(i).EQ.1)GO TO 981
               iSYM=ORBSYM(i)
               fga=flov(isym+1,3)-no
               lga=flov(isym+1,4)-no
               lga2=lga
               if(lga.gt.be)lga2=be
               fga2=fga
               if(fga.lt.be)fga2=be
               BEisym=IEOR(iSYM,BESYM)
               IF(BEisym.NE.TSYM-1)GO TO 981
C
               XBA=0
               DO 425 A=1,NV
                  IF(FZV(A).EQ.1)GO TO 425
                  ASYM=ORBSYM(A+NO)
                  BEA=ITV(BE)+A
                  ia=itr(a)+i
                  BSYM=IEOR(ASYM,BEisym)
                  FB=FLOV(BSYM+1,3)-NO
                  LB=FLOV(BSYM+1,4)-NO
                  IF(LB.GT.A)LB=A
                  DO 420 B=FB,LB
                     XBA=XBA+1
                     BEB=ITV(BE)+B
                     ib=itr(b)+i
                     beaib=nfoff(bea)+nfadd(ib)-shift
                     bebia=nfoff(beb)+nfadd(ia)-shift
                     TMP=fIN(beaib)+fIN(bebia)
                     TMV=fIN(beaib)-fIN(bebia)
                     SUX3(XBA)=TMP*0.5D0
                     SUX4(XBA)=TMV*0.5D0
  420             CONTINUE
                  IF(ASYM.EQ.BSYM)SUX3(XBA)=SUX3(XBA)*0.5D0
  425          CONTINUE
C
C              GO FOR THE N**6 LOOP
C
               DO 577 UV=1,LIMO
                  TMP = 0.0D0
                  TMV = 0.0D0
                  DO 529 AB=1,LIMV
                     TMP=TMP+SUX3(AB)* AR3(AB,UV)
                     TMV=TMV+SUX4(AB)* AR4(AB,UV)
  529             CONTINUE
                  SUX1(UV)=-TMP
                  SUX2(UV)=-TMV
  577          CONTINUE
C
               UV=0
               DO 668 U=1,NO
                  IF(FZO(U).EQ.1)GO TO 668
                  USYM=ORBSYM(U)
                  VSYM=IEOR(USYM,BEisym)+1
                  SV=FLOV(VSYM,1)
                  LV=FLOV(VSYM,2)
                  IF(LV.GT.U )LV=U
                  DO 663 V=SV,LV
                     UV=UV+1
                     do 661 ga=fga,lga2
                     zl=zlx(be,ga)
                     UVBEGA=TOFF(U,V,zl)+TADD(BE,GA)
                     T2(UVBEGA)=T2(UVBEGA)+(SUX1(UV)+SUX2(UV))*t1o(i,ga)
                 IF(U.NE.V)THEN
                    VUBEGA=TOFF(V,U,zl)+TADD(BE,GA)
                    T2(VUBEGA)=T2(VUBEGA)+(SUX1(UV)-SUX2(UV))*t1o(i,ga)
                 ENDIF
  661                CONTINUE
                     do 662 ga=fga2,lga 
                     zl=zlx(ga,be)
                     vugabe=TOFF(v,u,zl)+TADD(ga,be)
                     T2(vugabe)=T2(vugabe)+(SUX1(UV)+SUX2(UV))*t1o(i,ga)
                 IF(U.NE.V)THEN
                    uvgabe=TOFF(u,v,zl)+TADD(ga,be)
                    T2(uvgabe)=T2(uvgabe)+(SUX1(UV)-SUX2(UV))*t1o(i,ga)
                 ENDIF
  662                CONTINUE
  663             CONTINUE
  668          CONTINUE
C
  981       CONTINUE
            shift=shift+dim75
  991    CONTINUE
  992 CONTINUE
      CALL RCLOSE(ITAP75,3)
C
      RETURN
      END
