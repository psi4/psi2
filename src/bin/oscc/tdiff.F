      SUBROUTINE TDIFF(T2,T1,T2O,T1O,T1NORM,DELT,NO,NV,NDIMT2,spin,ssov,
     .                 ORBSYM,FLOV,NIRRED,TOFF,TADD,FZO,FZV,XE,XEO,
     .                 ncl,nop,t1diag)
      IMPLICIT INTEGER(A-Z)
      REAL*8 XE,XEO,DELE,T1NORM,DELT
      REAL*8 T2(NDIMT2),T1(NO,NV),T2O(NDIMT2),T1O(NO,NV),t1sum(2)
      INTEGER ORBSYM(NO+NV),FLOV(8,2,4),TOFF(NO,NO),TADD(NV,NV),
     .        FZO(NO),FZV(NV),spin(*),ssov(2,4)
ctdc.1.30.95
      integer y,i,a,gamma,closed,opn,max,fa,la,ncl(8),nop(8)
      real*8 sum, tia_a(65025), tia_b(65025), t1diag
C
      T1NORM=0.0D+00
      t1sum(1)=0.0d0
      t1sum(2)=0.0d0
      DELT=0.0D+00

      do i=1,65025
            tia_a(i) = 0.0d0
            tia_b(i) = 0.0d0
      enddo
C
      iocc=0
      DO 10 U=1,NO
         IF(FZO(U).EQ.1)GO TO 10
         iocc=iocc+1
         USYM=ORBSYM(U)
         uspn=spin(u)
         fbe=ssov(uspn,3)-no
         lbe=ssov(uspn,4)-no
         DO 20 BE=fbe,lbe
            BESYM=ORBSYM(BE+NO)
            IF(USYM.EQ.BESYM)THEN
               DELT=DELT+(T1(U,BE)-T1O(U,BE))*(T1(U,BE)-T1O(U,BE))
               T1sum(uspn) =T1sum(uspn) +T1(U,BE)*T1(U,BE)
            ENDIF
            UBESYM=IEOR(USYM,BESYM)
            DO 30 V=1,U-1
               IF(FZO(V).EQ.1)GO TO 30
               vspn=spin(v)
               VSYM=ORBSYM(V)
               GASYM=IEOR(UBESYM,VSYM)+1
               FGA=FLOV(GASYM,vspn,3)-NO
               LGA=FLOV(GASYM,vspn,4)-NO
               if(lga.ge.be)lga=be-1
               DO 40 GA=FGA,LGA
                  UVBEGA=TOFF(U,V)+TADD(BE,GA)
                  DELT=DELT
     .     +      (T2O(UVBEGA)-T2(UVBEGA))*(T2O(UVBEGA)-T2(UVBEGA))
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
C
      DELT=DSQRT(DELT)
      T1NORM=0.5d0*(DSQRT(T1sum(1))+dsqrt(t1sum(2)))/dsqrt(1.0d0*iocc)
C
Ctdc.1.30.95  Start new T1 Diagnostic section
Ctdc  This diagnostic excludes certain T1 amplitudes in order to maintain
Ctdc  consistency with the closed-shell diagnostic.  Note that MOs(spin-
Ctdc  orbitals) are ordered by alpha-occupied/irrep, then
Ctdc  beta-occupied/irrep, then alpha-unoccupied/irrep, and then
Ctdc  beta-unoccupied/irrep. tia_a() and tia_b are alpha and beta spin t1s
Ctdc  with notation the same as the original paper (see fentry()).  They are
Ctdc  dimensioned to 65025 such that the largest possible case (with 255 max
Ctdc  basis functions) will fit.  If the basis set limit is expanded, simply
Ctdc  enlarge this array to max(no*nv), where no == #occupied sos and nv =
Ctdc  #unoccupied sos.
      y = 1
      i = 1
      sum = 0.0
      do 100 gamma = 1, nirred
         closed = ncl(gamma)
         opn = nop(gamma)
         max = i + closed - 1
         do 110 i = i, max
            if(fzo(i).eq.1) go to 110
            fa = flov(gamma,1,3) - no
            la = flov(gamma,1,4) - no
            do 120 a = fa,la
               if(orbsym(a+no)+1.eq.gamma) then
                  tia_a(y) = t1(i,a)
                  y = y + 1
               endif
 120        continue
 110     continue
         max = i + opn - 1
         do 130 i = i,max
            fa = flov(gamma,1,3) - no
            la = flov(gamma,1,4) - no
            do 140 a = fa,la
               sum = sum + t1(i,a) * t1(i,a)
 140        continue
 130     continue
 100  continue

      y = 1
      do 150 gamma = 1, nirred
         closed = ncl(gamma)
         opn = nop(gamma)
         max = i + closed - 1
         do 160 i = i,max
            if(fzo(i).eq.1) go to 160
            fa = flov(gamma,2,3) - no + opn
            la = flov(gamma,2,4) - no
            do 170 a = fa,la
               if(orbsym(a+no)+1.eq.gamma) then
                  tia_b(y) = t1(i,a)
                  y = y + 1
               endif
 170        continue
            fa = flov(gamma,2,3) - no
            la = fa + opn - 1
            do 180 a = fa,la
               sum = sum + t1(i,a) * t1(i,a)
 180        continue
 160     continue
 150  continue

      sum = sum*2.0d0

      do 200 i = 1,(y-1)
            sum = sum + (tia_a(i) + tia_b(i))*
     .                  (tia_a(i) + tia_b(i))
 200  continue

      t1diag = 0.5d0 * (dsqrt(sum)) / dsqrt(1.0d0*iocc)
ctdc.1.30.95
          
      RETURN
      END
