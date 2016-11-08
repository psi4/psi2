c
      subroutine fint(t1,t1o,t2,t2o,gv,hv,fock,fbuff,aux4,
     .                toff,tadd,zlx,nfints,nfoff,nfadd,
     .                fzo,fzv,orbsym,flov,ipq,itr,itv,
     .                no,nv,nirred)
      implicit integer (a-z)
      real*8 t2(*),t2o(*),fock(*),fbuff(*),gv(nv,nv),
     .       hv(nv,nv),t1o(no,nv),t1(no,nv),aux4(no,nv)
      integer toff(no,no,2),tadd(nv,nv),
     .        zlx(nv,nv),fzo(no),fzv(nv),orbsym(*),flov(nirred,4),
     .        ipq(*),itr(*),itv(*),nfoff(*),nfadd(*),nfints(*)
      real*8 val,val1,val2,tmp,tmv,delta,valgv,val3
c
C     F INTS CONTRIBUTIONS:   F * Tau TO T1 EQN
C                             F * T1 TO GV
C
      itap75=75
      CALL RFILE(ITAP75)
      PON75=1
      shift=0
      DO 13799 BE=1,NV
      if(fzv(be).eq.1)go to 13799
      dim75=nfints(be)
c     size75=size75+dim75
c     write(6,*)'be ,nfints(be)',be,nfints(be)
      CALL WREADW(ITAP75,FBUFF,INTOWP(DIM75),PON75,PON75)
      BESYM=ORBSYM(BE+NO)
      DO 13791 A=1,NV
      IF(FZV(A).EQ.1)GO TO 13791
      ASYM=ORBSYM(A+NO)
      BEASYM=IEOR(BESYM,ASYM)
      BEA=ITV(BE)+A
      FBEA=NFOFF(BEA)-shift
C
      VALGV=hv(be,a)
      FU=FLOV(BESYM+1,1)
      LU=FLOV(BESYM+1,2)
      DO 13781 I=1,NO
      IF(FZO(I).EQ.1)GO TO 13781
      ISYM=ORBSYM(I)
      IA=ITR(A)+I
      FFIA=NFADD(IA)-shift
CFOK  VALGV=VALGV-FOCK(FIA)*T1O(I,BE)
      BSYM=IEOR(ISYM,BEASYM)+1
      FB=FLOV(BSYM,3)-NO
      LB=FLOV(BSYM,4)-NO
      DO 13771 B=FB,LB
      ZL=ZLX(A,B)
      IB=ITR(B)+I
      BEB=ITV(BE)+B
      BEAIB=FBEA+NFADD(IB)
      BEBIA=NFOFF(BEB)+FFIA
      VAL1=FBUFF(BEAIB)
      VAL2=FBUFF(BEBIA)
      TMP=VAL1+VAL1-VAL2
      VALGV=VALGV+TMP*T1O(I,B)
      DO 13761 U=FU,LU
      UIAB=TOFF(U,I,ZL)+TADD(A,B)
      VAL3=T2O(UIAB)+T1O(U,A)*T1O(I,B)
      aux4(u,be)=+aux4(u,be)+TMP*VAL3
13761 CONTINUE
13771 CONTINUE
13781 CONTINUE
      GV(BE,A)=+VALGV
C
13791 CONTINUE
      shift=shift+dim75
13799 CONTINUE
      CALL RCLOSE(ITAP75,3)
C
C     ADD INTO T1
C
      DO 291 U=1,NO
         IF(FZO(U).EQ.1)GO TO 291
         USYM=ORBSYM(U)
         FBE=FLOV(USYM+1,3)-NO
         LBE=FLOV(USYM+1,4)-NO
         DO 271 BE=FBE,LBE
            T1(U,BE)=T1(U,BE)+AUX4(U,BE)
  271    CONTINUE
  291 CONTINUE
C
c     write(6,*)' size75=',size75
      return
      end
