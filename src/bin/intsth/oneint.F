      subroutine oneint(a,chg,cx,eta,ica,icb,icxst,il,
     1   iprst,lmnp1,mcons,nc,ncon,nf,nfct,
     2   ngw,npair,nprir,nrcr,nt,ntl,
     3   ntu,x,y,z,zet)
C
C  CALCULATE ONE-ELECTRON INTEGRALS (OVERLAP, KE, OR PE)
C
*mdc*on vax ibm cray
      implicit real*8 (a-h,o-z)
*mdc*off
      character*2 nam
      logical es, esf, esfc, igueq1, jgueq1
      COMMON /OUT/    IOUT,ICHECK,IFLAG3,NSECT,NINTS
      common /parmr/ cutoff, tol
      common /parmi/ mblu, mccu, mconu, mcu, mcxu,
     1   mrcru, mstu, msu, nblu, ng, ns, nst
      common /ntgr/ kaords, mconsu, mgcsu, mru, mcsu, mctu, mcru,
     1   msfu, mgu, msfru, mnru, mnsfu, mpru, ngcs, nu, mcxu2,
     2   lxyzir(3), inam(5), nnam, mdum(30)
      common /one/dxij,dyij,dzij,fnfct,rr,xij,xijm,yij,yijm,zij,zijm,
     1  ibl1,ibl2,icxi1,icxi2,ij,ijsf,ijx(225),ijy(225),ijz(225),ic,
     2  icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl,
     3  jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
      common /dim21/ ipq(256)
      COMMON /TAPE/   NTAPE
      COMMON /OUTBUF/ IBLK, IBUF, LBL(2980), SPI(2980)
c
      dimension a(*), chg(*), cx(*), eta(mrcru,mconu,*),
     1   ica(mcu,msu,*), icb(4,24,*), icxst(2,*), il(*),
     2   iprst(*), lmnp1(*), mcons(*),
     3   nc(*), ncon(*), nf(*), nfct(*), ngw(*),
     5   npair(2,*), nprir(2,mstu,*), nrcr(*), nt(*), ntl(*),
     6   ntu(*), x(mcu,*), y(mcu,*), z(mcu,*), zet(mconu,*)
*mdc*on ibm cray vax
*      parameter ( maskij=2**30-1)
*mdc*off
      dimension nam(5)
      data nam/' s',' t',' v',' c',' z'/
      DATA IDTR/2980/
c
      iblk=0
      IBUF = 0
      NSCTLD=0
      NBFLD=0
      ijsf = 0
      isf = 0
      do 712 is = 1,ns
      ifu = nf(is)
      if(ifu.eq.0) go to 712
      icu=nc(is)
      isfis=isf
      jsf = 0
      do 708 js = 1,is
      jsfjs=jsf
      es=is.eq.js
      ic=icu
      ncc = 0
      if(.not.es) go to 200
c
c  center combinations for is=js
c
      do 140 jc = 1, ic
      igwu=0
      do 132 ig = 1, ng
      ic1 = ica(ic,is,ig)
      jc1 = ica(jc,is,ig)
      if(ic1.ne.ic.and.jc1.ne.ic) go to 132
      ic2=max(ic1,jc1)
      jc2=min(ic1,jc1)
c
c  do 112 icc=ncc,1,-1
c
      icc = ncc
  104 if (icc .eq. 0) go to 116
      if(jc2-icb(2,1,icc)) 108,136,116
  108 icc = icc - 1
      go to 104
c
c 112  end icc loop
c
  116 if(jc2.ne.jc) go to 132
c
c  do 124 igw=igwu,1,-1
c
      igw=igwu
  120 if(igw.eq.0) go to 128
      if(ic1.eq.icb(1,igw,ncc+1).and.jc1.eq.icb(2,igw,ncc+1)) go to 132
      igw=igw-1
      go to 120
c
c 124  end igw loop
c
  128 igwu=igwu+1
      icb(1,igwu,ncc+1)=ic1
      icb(2,igwu,ncc+1)=jc1
  132 continue
      ncc = ncc + 1
      nfct(ncc) = ic
      ngw(ncc)=igwu
      go to 140
  136 nfct(icc) = nfct(icc) + ic
  140 continue
      iccu = ncc - 1
      if (iccu .eq. 0 ) go to 300
      do 144 icc = 1, iccu
  144 nfct(icc) = nfct(icc)/2
      go to 300
c
c  center combinations for is>js
c
  200 jfu = nf(js)
      if(jfu.eq.0) go to 708
      jcu = nc(js)
      do 240 jc = 1, jcu
      do 232 ig = 1, ng
      if (ica(ic,is,ig) .ne. ic) go to 232
      jc1 = ica (jc,js,ig)
c
c  do 212 icc=ncc,1,-1
c
      icc = ncc
  204 if (icc .eq. 0) go to 232
      if(jc1-icb(2,1,icc)) 208,236,232
  208 icc = icc - 1
      go to 204
c
c 212  end icc loop
c
  232 continue
      ncc = ncc + 1
      icb(1,1,ncc)=ic
      icb(2,1,ncc)=jc
      nfct(ncc) = ic
      go to 240
  236 nfct(icc) = nfct(icc) + ic
  240 continue
c
c  evaluate integrals
c
  300 IF(NCC.GT.MCCU) THEN
        WRITE(IOUT,988) NCC
  988   FORMAT(' MCCU TOO SMALL',I4)
        call qabort
      END IF
      isf=isfis
      do 704 if = 1, ifu
      isf = isf + 1
      icons=mcons(isf)
      ircru=nrcr(icons)
      lit=lmnp1(icons)
      igu=ncon(icons)
      itl=ntl(isf)
      itu=ntu(isf)
      jsf=jsfjs
      if(es) jfu=if
      do 700 jf = 1, jfu
      jsf = jsf + 1
      esf=isf.eq.jsf
      jcons=mcons(jsf)
      jrcru=nrcr(jcons)
      ijsf = ijsf + 1
      ibl1=0
      ibl2=0
c     do 308 ist=1,nst
c     iop=nopir(ist)
c     if(iop.eq.0) go to 308
c     if(esf) ibl1=ibl1+iop*nprir(1,ist,ijsf)
      if(esf) ibl1=ibl1+nprir(1,1,ijsf)
c     ibl2=ibl2+iop*nprir(2,ist,ijsf)
      ibl2=ibl2+nprir(2,1,ijsf)
c 308 continue
      ntij2=(nt(isf)*nt(jsf))*npair(2,ijsf)
      icxi2=icxst(2,ijsf)
      if(.not.esf) go to 312
      nblt1=npair(1,ijsf)
      ntij1=(nt(isf)*nt(jsf))*nblt1
      nblu=ircru*ibl1+ipq(ircru)*ibl2
      icxi1=icxst(1,ijsf)
      go to 316
  312 nblu=ircru*jrcru*ibl2
  316 continue
      if(nblu.eq.0) go to 700
      ipr=iprst(ijsf)
cges 16jan90  this should be nblu?  zet found it.
c     IF(NBL.GT.MBLU) THEN
c       WRITE(ICHECK,989) NBL
      IF(NBLu.GT.MBLU) THEN
        WRITE(ICHECK,989) NBLu
  989   FORMAT(' MBLU TOO SMALL',I8)
        call qabort
      END IF
      mblu=mblu-nblu
      ljt=lmnp1(jcons)
      jgu=ncon(jcons)
      jtl=ntl(jsf)
      jtu=ntu(jsf)
      call wzero(nblu,a,1)
      if(es.and..not.esf) go to 500
c
      do 424 icc = 1, ncc
      ic=icb(1,1,icc)
      jc=icb(2,1,icc)
      fnfct=nfct(icc)
  424 call stvcz(a,chg,cx,eta,a(nblu+1),nc,npair,
     1   nprir,nt,x,y,z,zet)
      go to 600
c
  500 do 524 icc = 1, ncc
      igwu=ngw(icc)
      itst=0
      do 520 iprm=1,2
      if(iprm.ne.1) go to 504
      inx=1
      jnx=2
      go to 508
  504 inx=2
      jnx=1
  508 jc=icb(jnx,1,icc)
      if(jc.eq.itst) go to 520
      ic=icb(inx,1,icc)
      nwt=1
      if(igwu.eq.1) go to 516
      do 512 igw=2,igwu
      if(icb(jnx,igw,icc).lt.jc) go to 520
  512 nwt=nwt+1
  516 fnfct=(nwt*nfct(icc))
      call stvcz(a,chg,cx,eta,a(nblu+1),nc,npair,
     1   nprir,nt,x,y,z,zet)
      itst=jc
  520 continue
  524 continue
c
c  put integrals and labels into buffers
c
  600 mblu=mblu+nblu
      ibl=0
      do 640 ircr=1,ircru
      if(esf) jrcru=ircr
      do 640 jrcr=1,jrcru
      if(esf.and.jrcr.eq.ircr) then
        npairi=nprir(1,1,ijsf)
        npairt=npair(1,ijsf)
      else
        npairi=nprir(2,1,ijsf)
        npairt=npair(2,ijsf)
      endif
      if(npairi.eq.0) go to 636
      do 624 ipair=1,npairi
      ibl=ibl+1
      if(abs(a(ibl)).le.cutoff) go to 624
c
c  write integral record
c
      IF(IBUF.EQ.IDTR) THEN
        NBFLD=NBFLD+1
        NSCTLD=NSCTLD+I2SEC(8924)
        CALL SWRIT(NTAPE,IBLK,8942)
C       WRITE (NTAPE) IBLK, IBUF, LBL, SPI
        IBUF=0
      END IF
      IBUF = IBUF + 1
      LBL(IBUF)=IL(IPR+IPAIR)
      spi(ibuf)=a(ibl)
  620 continue
  624 continue
  632 ipr=ipr+npairt
  636 continue
  640 continue
  700 continue
  704 continue
  708 continue
  712 continue
c
c  write out the last record of integrals
c
C     WRITE(ICHECK,920) (SPI(IQ),IQ=1,IBUF)
C 920 FORMAT (/,(1X,4F10.1))
C
      IBLK=1
      NSCTLD=NSCTLD+I2SEC(8924)
      CALL SWRIT(NTAPE,IBLK,8942)
C     WRITE (NTAPE) IBLK, IBUF, LBL, SPI
      NSECT=NSECT+NSCTLD
      IBUF=IDTR*NBFLD+IBUF
      NBFLD=NBFLD+1
      IF (IFLAG3.GT.0) WRITE(ICHECK,900) IBUF,NAM(NNAM),NSCTLD
      WRITE (IOUT,900) IBUF, NAM(NNAM), NSCTLD
      NINTS=NINTS+IBUF
      RETURN
  900 FORMAT(I19,A2,' INTEGRALS IN ',I9,' SECTORS')
      end
