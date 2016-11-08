      subroutine jandk(eta,g4,is,isf,ic,js,jsf,jc,ks,ksf,kc,ls,lsf,
     1   lc,lmnp1,mcons,ncon,nrcr,nt,ntl,ntu,x,y,z,zet)
*mdc*on vax ibm cray
      implicit real*8 (a-h,o-z)
*mdc*off
      logical esfc, esfcij, esfckl, igueq1, jgueq1, kgueq1, lgueq1
      COMMON /OUT/    IOUT,ICHECK,IFLAG3,NSECT,NINTS
      common /parmi/ mblu, mccu, mconu, mcu, mcxu,
     1   mrcru, mstu, msu, nblu, ng, ns, nst
ctph  common/root/xx,u(9),w(9),nroots
      common/root/xx,u(10),w(10),nroots
      common /shlnos/lit,ljt,lkt,llt,ij,ijkl
      common /index/ ijgt(225), ijx(225), ijy(225), ijz(225), ik(225),
     1               klgt(225), klx(225), kly(225), klz(225)
      common /shlinf/ xi,yi,zi,rij,xk,yk,zk,rkl,igu,jgu,kgu,lgu,nc4,
     1   nc3,nc2,nc1,igueq1,jgueq1,kgueq1,lgueq1
      common/misc/xij,xijm,yij,yijm,zij,zijm,xkl,xklm,ykl,yklm,zkl,zklm,
     1icons,ircru,jcons,jrcru,kcons,krcru,lcons,lrcru,esfc,esfcij,esfckl
      common /setint/ dxij,dyij,dzij,dxkl,dykl,dzkl,bp01,b00,b10,xcp00,
ctph 1  xc00,ycp00,yc00,zcp00,zc00,f00,in(9),kn(9),ni,nj,nk,nl,nmax,mmax
     1  xc00,ycp00,yc00,zcp00,zc00,f00,in(10),kn(10),ni,nj,nk,nl,nmax,
     2  mmax
      common /dim21/ ipq(256)
      common / maxval / mxblu
c
      parameter ( a1s2=0.5d0 )
      dimension eta(mrcru,mconu,*), g4(*), lmnp1(*), mcons(*),
     1   ncon(*), nrcr(*), nt(*), ntl(*), ntu(*), x(mcu,*),
     2   y(mcu,*), z(mcu,*), zet(mconu,*)
      dimension ix(35), iy(35), iz(35), jx(35), jy(35), jz(35),
     1          kx(35), ky(35), kz(35), lx(35), ly(35), lz(35)
      data lx /   0,  1,  0,  0,  2,  0,  0,  1,  1,  0,
     1            3,  0,  0,  2,  2,  1,  0,  1,  0,  1,
     2            4,  0,  0,  3,  3,  1,  0,  1,  0,  2,
     3            2,  0,  2,  1,  1/
      data kx /   0,  5,  0,  0, 10,  0,  0,  5,  5,  0,
     1           15,  0,  0, 10, 10,  5,  0,  5,  0,  5,
     2           20,  0,  0, 15, 15,  5,  0,  5,  0, 10,
     3           10,  0, 10,  5,  5/
      data jx /   0, 25,  0,  0, 50,  0,  0, 25, 25,  0,
     1           75,  0,  0, 50, 50, 25,  0, 25,  0, 25,
     2          100,  0,  0, 75, 75, 25,  0, 25,  0, 50,
     3           50,  0, 50, 25, 25/
      data ix /   1,126,  1,  1,251,  1,  1,126,126,  1,
     1          376,  1,  1,251,251,126,  1,126,  1,126,
     2          501,  1,  1,376,376,126,  1,126,  1,251,
     3          251,  1,251,126,126/
      data ly /   0,  0,  1,  0,  0,  2,  0,  1,  0,  1,
     1            0,  3,  0,  1,  0,  2,  2,  0,  1,  1,
     2            0,  4,  0,  1,  0,  3,  3,  0,  1,  2,
     3            0,  2,  1,  2,  1/
      data ky /   0,  0,  5,  0,  0, 10,  0,  5,  0,  5,
     1            0, 15,  0,  5,  0, 10, 10,  0,  5,  5,
     2            0, 20,  0,  5,  0, 15, 15,  0,  5, 10,
     3            0, 10,  5, 10,  5/
      data jy /   0,  0, 25,  0,  0, 50,  0, 25,  0, 25,
     1            0, 75,  0, 25,  0, 50, 50,  0, 25, 25,
     2            0,100,  0, 25,  0, 75, 75,  0, 25, 50,
     3            0, 50, 25, 50, 25/
      data iy /   1,  1,126,  1,  1,251,  1,126,  1,126,
     1            1,376,  1,126,  1,251,251,  1,126,126,
     2            1,501,  1,126,  1,376,376,  1,126,251,
     3            1,251,126,251,126/
      data lz /   0,  0,  0,  1,  0,  0,  2,  0,  1,  1,
     1            0,  0,  3,  0,  1,  0,  1,  2,  2,  1,
     2            0,  0,  4,  0,  1,  0,  1,  3,  3,  0,
     3            2,  2,  1,  1,  2/
      data kz /   0,  0,  0,  5,  0,  0, 10,  0,  5,  5,
     1            0,  0, 15,  0,  5,  0,  5, 10, 10,  5,
     2            0,  0, 20,  0,  5,  0,  5, 15, 15,  0,
     3           10, 10,  5,  5, 10/
      data jz /   0,  0,  0, 25,  0,  0, 50,  0, 25, 25,
     1            0,  0, 75,  0, 25,  0, 25, 50, 50, 25,
     2            0,  0,100,  0, 25,  0, 25, 75, 75,  0,
     3           50, 50, 25, 25, 50/
      data iz /   1,  1,  1,126,  1,  1,251,  1,126,126,
     1            1,  1,376,  1,126,  1,126,251,251,126,
     2            1,  1,501,  1,126,  1,126,376,376,  1,
     3          251,251,126,126,251/
c
      ngtl=1
      ngtk=nt(lsf)
      ngtj=ngtk*nt(ksf)
      ngti=ngtj*nt(jsf)
c     ----- ishell -----
      xi=x(ic,is)
      yi=y(ic,is)
      zi=z(ic,is)
      icons=mcons(isf)
      igu=ncon(icons)
      lit=lmnp1(icons)
      mini=ntl(isf)
      maxi=ntu(isf)
      ircru=nrcr(icons)
c     ----- jshell -----
      xj=x(jc,js)
      yj=y(jc,js)
      zj=z(jc,js)
      jcons=mcons(jsf)
      jgu=ncon(jcons)
      ljt=lmnp1(jcons)
      minj=ntl(jsf)
      maxj=ntu(jsf)
      jrcru=nrcr(jcons)
c     ----- prepare items for pairs of (i,j) functions -----
      esfcij=isf.eq.jsf.and.ic.eq.jc
      xij=a1s2*(xi+xj)
      dxij=xi-xj
      xijm=a1s2*dxij
      yij=a1s2*(yi+yj)
      dyij=yi-yj
      yijm=a1s2*dyij
      zij=a1s2*(zi+zj)
      dzij=zi-zj
      zijm=a1s2*dzij
      rij=dxij*dxij+dyij*dyij+dzij*dzij
      ij=0
      do 150 i=mini,maxi
      nx=ix(i)
      ny=iy(i)
      nz=iz(i)
c     if(esfcij) maxj=i
      do 150 j=minj,maxj
      ij=ij+1
      ijx(ij)=nx+jx(j)
      ijy(ij)=ny+jy(j)
      ijz(ij)=nz+jz(j)
  150 ijgt(ij)=ngti*(i-mini)+ngtj*(j-minj)+1
c     ----- kshell -----
      xk=x(kc,ks)
      yk=y(kc,ks)
      zk=z(kc,ks)
      kcons=mcons(ksf)
      kgu=ncon(kcons)
      lkt=lmnp1(kcons)
      mink=ntl(ksf)
      maxk=ntu(ksf)
      krcru=nrcr(kcons)
c     ----- lshell -----
      xl=x(lc,ls)
      yl=y(lc,ls)
      zl=z(lc,ls)
      lcons=mcons(lsf)
      lgu=ncon(lcons)
      llt=lmnp1(lcons)
      minl=ntl(lsf)
      maxl=ntu(lsf)
      lrcru=nrcr(lcons)
c     ----- prepare items for pairs of (k,l) functions -----
      esfckl=ksf.eq.lsf.and.kc.eq.lc
      xkl=a1s2*(xk+xl)
      dxkl=xk-xl
      xklm=a1s2*dxkl
      ykl=a1s2*(yk+yl)
      dykl=yk-yl
      yklm=a1s2*dykl
      zkl=a1s2*(zk+zl)
      dzkl=zk-zl
      zklm=a1s2*dzkl
      rkl=dxkl*dxkl+dykl*dykl+dzkl*dzkl
      kl=0
      do 350 k=mink,maxk
      nx=kx(k)
      ny=ky(k)
      nz=kz(k)
c     if(esfckl) maxl=k
      do 350 l=minl,maxl
      kl=kl+1
      klx(kl)=nx+lx(l)
      kly(kl)=ny+ly(l)
      klz(kl)=nz+lz(l)
  350 klgt(kl)=ngtk*(k-mink)+ngtl*(l-minl)
c     ----- prepare miscellaneous items -----
      nroots=(lit+ljt+lkt+llt-2)/2
      esfc=isf.eq.ksf.and.ic.eq.kc.and.jsf.eq.lsf.and.jc.eq.lc
      do 500 i=1,ij
  500 ik(i)=kl
      ijkl=ij*kl
      ngout=ngti*nt(isf)
c     ----- allocate space for transforming ao integrals -----
      igueq1=igu.eq.1
      jgueq1=jgu.eq.1
      kgueq1=kgu.eq.1
      lgueq1=lgu.eq.1
      nc1=lrcru*ngout
      if(.not.esfckl) go to 600
      nc2=ipq(krcru+1)*ngout
      go to 610
  600 nc2=krcru*nc1
      if((.not.igueq1).or.(.not.esfc).or.jgueq1) go to 610
      nc3=ipq(jrcru+1)*ngout
      nc4=nc3
      go to 660
  610 nc3=jrcru*nc2
      if(.not.esfc) go to 640
      if(.not.esfcij) go to 620
      nc4=ipq(ircru+1)
      go to 630
  620 nc4=ircru*jrcru
  630 nc4=ipq(nc4+1)*ngout
      go to 660
  640 if(.not.esfcij) go to 650
      nc4=ipq(ircru+1)*nc2
      go to 660
  650 nc4=ircru*nc3
  660 if(.not.igueq1) go to 670
      ic3=0
      go to 680
  670 ic3=nc4
  680 if(.not.jgueq1) go to 690
      ic2=ic3
      go to 700
  690 ic2=ic3+nc3
  700 if(.not.kgueq1) go to 710
      ic1=ic2
      go to 720
  710 ic1=ic2+nc2
  720 if(.not.lgueq1) go to 730
      ic0=ic1
      go to 740
  730 ic0=ic1+nc1
  740 mgout=ic0+ngout
      mxblu = max( mxblu, mcxu+nblu+mgout )
      if(mgout.gt.mblu)then
        mgout=mcxu+nblu+mgout
        write (iout,989) mgout
  989   format(' jandk: mblu too small',i8)
        call qabort
      endif
c     ----- compute two-electron integrals ----
      call wzero(nc4,g4,1)
      if(ijkl.gt.1) go to 760
      call s0000(eta,g4,g4(ic3+1),g4(ic2+1),g4(ic1+1),g4(ic0+1),zet)
      return
  760 call genral(eta,g4,g4(ic3+1),g4(ic2+1),g4(ic1+1),g4(ic0+1),zet)
      return
      end
