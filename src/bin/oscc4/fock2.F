      SUBROUTINE FOCK2 (NT,HH,FOCK,EREF,ENUC,EIG,
     .                  IPQ,ITR,ITV,NO,NV,xIN,
     .                  iOFF,iADD,
     .                  ORBSYM,FLOV,NIRRED,I79,FZO,FZV,
     .                  ntmo,spin,map,ctscf2,nfzo,nfzv,omap,
     .                  dvec,xvec,wvec,ss,gg,xtri,gaux,
     .                  xarr,gamma,beta,ctscf1,
     .                  nlamda,nsymhf,scftc1,scftc2)
      IMPLICIT INTEGER (A-Z)
      REAL*8 HH(*),FOCK(*),xIN(*),EIG(*)
      REAL*8 EREF,ENUC
      INTEGER IPQ(*),ITR(*),ITV(*),FZO(no),FZV(nv),
     .        iOFF(*),iADD(*),nfzo(*),nfzv(*),
     .        ORBSYM(nt),FLOV(8,2,4),scftc1(*),scftc2(*),
     .        spin(nt),map(nt),ctscf2(*),omap(nt),ctscf1(*),
     .        nlamda(*)
      REAL * 8 DVEC(ntmo,ntmo),XVEC(ntmo,ntmo),WVEC(ntmo*ntmo),
     .         SS(*),gg(*),XTRI(*),gaux(*),
     .         XARR(ntmo),GAMMA(ntmo),BETA(ntmo)
      common /occup/noa,nob
      it30=30

C
      NT2=NT*(NT+1)/2
      ntmo2=ntmo*(ntmo+1)/2
c
      write(6,*)
      write(6,*)'reading sym. orb. alpha & beta fock matrices from 79'
c 
      nmat=0
      do 127 i=1,nsymhf
         nl=nlamda(i)
         nl2=nl*(nl+1)/2
         nmat=nmat+nl2
  127 continue
c     write(6,*)'nmat=',nmat
C
C >>> GET OVERLAP MATRIX IN SO BASIS FROM TAPE 67
C
c     ITAP67=67
c     CALL RFILE(ITAP67)
c     CALL SREAD(ITAP67,SS,INTOWP(ntmo2))
c     WRITE(6,*)
c     WRITE(6,*)' SO OVERLAP MATRIX'
c     CALL PRINT(SS,ntmo2,ntmo,6)
C
C     DVEC IS THE SO TO MO EIGENVECTOR
C     WVEC IS THE AO TO MO EIGENVECTOR
C     XVEC IS THE SO TO AO TRANSFORMATION MATRIX
C
      CALL VECS(DVEC,XVEC,WVEC,'READ',ntmo,ntmo,it30)
C
c     Ct * S * C = 1
c
c     CALL ZERO (XVEC,ntmo*ntmo)
c     DO 148  I=1,ntmo
c        DO 146  A=1,ntmo
c           DO 142 MU=1,ntmo
c              DO 140 NU=1,ntmo
c                 MA=MAX0(MU,NU)
c                 MI=MIN0(MU,NU)
c                 MUNU=(MA*(MA-1))/2+MI
c                 XVEC(I,A)=XVEC(I,A)+DVEC(MU,I)*DVEC(NU,A)*SS(MUNU)
c 140          CONTINUE
c 142       CONTINUE
c 146    CONTINUE
c 148 CONTINUE
c     WRITE(6,*)
c     WRITE(6,*)' METRIC OF THE BASIS IN FILE30'
c     CALL MATOUT(XVEC,ntmo,ntmo,ntmo,ntmo,6)
c
      CALL ZERO (FOCK,NT2)
C
C >>> READ ONE-E MO INTS FROM FILE79
C
      ITAP79=79
COK   CALL RFILE(ITAP79)
c    this is no longer necessary except for the ref energy
      CALL WREADW(ITAP79,HH,INTOWP(NTmo2),1,ipn79)
c     WRITE(6,*)'ONE E HAMILTONIAN MO BASIS cc-order'
c     CALL PRINT(HH,NTmo2,NTmo,6)
C
c  read alpha fock matrix from 79
c
      ipn79=1+intowp(ntmo2)
      call wreadw(itap79,gg,intowp(nmat),ipn79,ipn79)
c
c     expand into full-matrix
c
      call zero(gaux,ntmo2)
      ij=0
      iof=0
      jof=0
      do 297 is=1,nsymhf
         nis=nlamda(is)
         do 296 i=1,nis
            do 295 j=1,i
               ij=ij+1
               it=i+iof
               jt=j+jof
               itjt=ipq(it)+jt
               gaux(itjt)=gg(ij)
  295       continue
  296    continue
         iof=iof+nis
         jof=jof+nis
  297 continue
c
c     write(6,*)
c     write(6,*)'alpha-fock matrix in symorb basis from 79'
c     call print(gaux,ntmo2,ntmo,6)
c
c     transform to mo basis and add one-el ham
c
      call zero(wvec,ntmo2)
      call zero(  gg,ntmo2)
      ij=0
      do 378 i=1,ntmo
         do 376 j=1,i
            ij=ij+1
            do 374 k=1,ntmo
            do 372 l=1,ntmo
            ma=max0(k,l)
            mi=min0(k,l)
            kl=(ma*(ma-1))/2+mi
            wvec(ij)=wvec(ij)+dvec(k,i)*gaux(kl)*dvec(l,j)
  372       continue
  374       continue
  376    continue
  378 continue
      do 379 i=1,ntmo2
         gg(i)=wvec(i)
  379 continue
c     write(6,*)
c     write(6,*)'alpha-fock in mo basis scf-order'
c     call print(gg,ntmo2,ntmo,6)
c
c     sort to cc order
c
      call zero(gaux,ntmo2)
      do 381 i=1,ntmo
         do 380 j=1,i
            ic=scftc1(i)
            jc=scftc1(j)
            ij=ipq(i)+j
            icjc=ipq(max0(ic,jc))+min0(ic,jc)
            gaux(icjc)=gg(ij)
  380    continue
  381 continue
c     write(6,*)
c     write(6,*)'alpha-fock in mo basis cc-order'
c     call print(gaux,ntmo2,ntmo,6)
c
c     put in fock 
c
      do 385 i=1,nt
         ispn=spin(i)
         if(ispn.ne.1)go to 385
c        xi=omap(i)
         xi=i
         if(i.gt.no)xi=i-nob
         do 383 j=1,i
            jspn=spin(j)
            if(jspn.ne.ispn)go to 383
c           xj=omap(j)
            xj=j
            if(j.gt.no)xj=j-nob
            ij=ipq(i)+j
            xij=ipq(max0(xi,xj))+min0(xi,xj)
            fock(ij)=gaux(xij)
  383    continue
  385 continue
c     write(6,*)
c     write(6,*)'fock-matrix (alpha only) in so basis '
c     call print(fock,nt2,nt,6)
 
C
c  read beta fock matrix from 79
c
      call wreadw(itap79,gg,intowp(nmat),ipn79,ipn79)
c
c     expand into full-matrix
c
      call zero(gaux,ntmo2)
      ij=0
      iof=0
      jof=0
      do 597 is=1,nsymhf
         nis=nlamda(is)
         do 596 i=1,nis
            do 595 j=1,i
               ij=ij+1
               it=i+iof
               jt=j+jof
               itjt=ipq(it)+jt
               gaux(itjt)=gg(ij)
  595       continue
  596    continue
         iof=iof+nis
         jof=jof+nis
  597 continue
c
c     write(6,*)
c     write(6,*)'beta-fock matrix in symorb basis from 79'
c     call print(gaux,ntmo2,ntmo,6)
c
c     transform to mo basis 
c
      call zero(wvec,ntmo2)
      call zero(  gg,ntmo2)
      ij=0
      do 578 i=1,ntmo
         do 576 j=1,i
            ij=ij+1
            do 574 k=1,ntmo
            do 572 l=1,ntmo
            ma=max0(k,l)
            mi=min0(k,l)
            kl=(ma*(ma-1))/2+mi
            wvec(ij)=wvec(ij)+dvec(k,i)*gaux(kl)*dvec(l,j)
  572       continue
  574       continue
  576    continue
  578 continue
      do 579 i=1,ntmo2
         gg(i)=wvec(i)
  579 continue
c     write(6,*)
c     write(6,*)'beta-fock in mo basis scf-order'
c     call print(gg,ntmo2,ntmo,6)
c
c     sort to cc order
c
      call zero(gaux,ntmo2)
      do 581 i=1,ntmo
         do 580 j=1,i
            ic=scftc2(i)
            jc=scftc2(j)
            ij=ipq(i)+j
            icjc=ipq(max0(ic,jc))+min0(ic,jc)
            gaux(icjc)=gg(ij)
  580    continue
  581 continue
c     write(6,*)
c     write(6,*)'beta-fock in mo basis cc-order'
c     call print(gaux,ntmo2,ntmo,6)
c
c     put in fock 
c
      do 585 i=1,nt
         ispn=spin(i)
         if(ispn.ne.2)go to 585
c        xi=omap(i)
         xi=i-noa
         if(i.gt.no)xi=i-ntmo
         do 583 j=1,i
            jspn=spin(j)
            if(jspn.ne.ispn)go to 583
c           xj=omap(j)
            xj=j-noa
            if(j.gt.no)xj=j-ntmo
            ij=ipq(i)+j
            xij=ipq(max0(xi,xj))+min0(xi,xj)
            fock(ij)=gaux(xij)
  583    continue
  585 continue
c     write(6,*)
c     write(6,*)'fock-matrix (total) in so basis '
c     call print(fock,nt2,nt,6)
c
      EREF=0.0D0
      DO 700 I=1,no
         xi=omap(i)
         xii=ipq(xi)+xi
         FII=IPQ(I)+I
         EREF=EREF+FOCK(FII)+hh(xii)
  700 CONTINUE
      EREF=EREF*0.5d0+ENUC
      WRITE(6,*)' '
      WRITE(6,638)EREF
  638 FORMAT(/,' REFERENCE ENERGY  =  ',F20.12)
c
      RETURN
      END
