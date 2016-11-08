      SUBROUTINE INITZ(Z1,Z2,I1MAT,MOINTS,LENINT,EIGVAL,IOFF,NO,NV,
     .              NTRO,NTRV,NO2,NV2,NTR,NBF,JOUT,IFLG,INO,INTRO,INTRV,
     .              NOV,INV,N2OV,FLOV,ORBSYM,NIRRED,
     .              DOFF,NOFF,DADD,NADD,ZLX,TOFF,TADD,NDIMT2,FZO,FZV,
     .              option)
      IMPLICIT INTEGER(A-Z)
      character*4 option
      REAL*8 Z1(NO,NV),Z2(NDIMT2),I1MAT(NO,NV),MOINTS(LENINT),
     .       EIGVAL(NBF),VAL2N,VAL2D,VAL1N,VAL1D
      INTEGER IOFF(NTR),FZO(*),FZV(*)
      INTEGER INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV)
      INTEGER DOFF(NOFF),DADD(NADD),TOFF(NO,NO,2),TADD(NV,NV),ZLX(NV,NV)
      DIMENSION FLOV(NIRRED,4),ORBSYM(NO+NV)
C
      CALL ZERO(Z1,NO*NV)
      CALL ZERO(Z2,NDIMT2)
c
      if (option.eq.'CCT ')then
          itap70=70
          call rfile(itap70)
          call wreadw(itap70,z1,intowp(no*nv),1,i70)
          call wreadw(itap70,z2,intowp(ndimt2),i70,i70)
          call rclose(itap70,3)
c         if(iflg.eq.0)then
c         write(6,*)'j1 is ',no*nv,'  long'
c         write(6,'(8f10.6)')z1
c         write(6,*)'j2 is ',ndimt2,'  long'
c         write(6,'(8f10.6)')(z2(i),i=1,ndimt2)
c         endif
      endif
C
      DO 10 U=1,NO
         IF(FZO(U).EQ.1)GO TO 10
         USYM=ORBSYM(U)
         DO 20 BE=1,NV
            IF(FZV(BE).EQ.1)GO TO 20
            BESYM=ORBSYM(BE+NO)
            UBESYM=IEOR(USYM,BESYM)
            VAL1N=I1MAT(U,BE)+I1MAT(U,BE)
            VAL1D=1.0D+00
            BEBE=IOFF(BE+NO)+BE+NO
            UU  =IOFF(U)+U
            IF (IFLG.EQ.0) VAL1D=EIGVAL(BEBE)-EIGVAL(UU)
            Z1(U,BE)=(-z1(u,be)-VAL1N)/VAL1D
            BEU=INO(BE)+U
            DO 30 GA=1,BE
               IF(FZV(GA).EQ.1)GO TO 30
               GAU=INO(GA)+U
               GASYM=ORBSYM(GA+NO)
               VSYM=IEOR(UBESYM,GASYM)+1
               FV=FLOV(VSYM,1)
               LV=FLOV(VSYM,2)
               DO 40 V=FV,LV
                  BEV=INO(BE)+V
                  GAV=INO(GA)+V
                  BEUGAV=DOFF(MAX0(BEU,GAV))+DADD(MIN0(BEU,GAV))
                  BEVGAU=DOFF(MAX0(BEV,GAU))+DADD(MIN0(BEV,GAU))
                  VAL2N=MOINTS(BEUGAV)+MOINTS(BEUGAV)+MOINTS(BEUGAV)+
     .                MOINTS(BEUGAV)-MOINTS(BEVGAU)-MOINTS(BEVGAU)
                  VAL2D=1.0D+00
            BEBE=IOFF(BE+NO)+BE+NO
            UU  =IOFF(U)+U
            GAGA=IOFF(GA+NO)+GA+NO
            VV  =IOFF(V)+V
                  IF (IFLG.EQ.0)
     .  VAL2D=EIGVAL(BEBE)+EIGVAL(GAGA)-EIGVAL(UU)-EIGVAL(VV)
                  UVBEGA=TOFF(U,V,1)+TADD(BE,GA)
                  Z2(UVBEGA)=(-z2(uvbega)-VAL2N)/VAL2D
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
      RETURN
      END
