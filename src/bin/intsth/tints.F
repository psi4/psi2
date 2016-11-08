      subroutine tints(eta,g2,g1,gout,zet)
*mdc*on vax ibm cray
      implicit real*8 (a-h,o-z)
*mdc*off
      logical esf, esfc, igueq1, jgueq1
      common /parmr/ cutoff, tol
      common /parmi/ mblu, mccu, mconu, mcu, mcxu,
     1   mrcru, mstu, msu, nblu, ng, ns, nst
      common /one/dxij,dyij,dzij,fnfct,rr,xij,xijm,yij,yijm,zij,zijm,
     1  ibl1,ibl2,icxi1,icxi2,ij,ijsf,ijx(225),ijy(225),ijz(225),ic,
     2  icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl,
     3  jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
      common /dim21/ ipq(256)
      common/stv/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      dimension eta(mrcru,mconu,*), g2(*), g1(*), gout(*), zet(mconu,*)
      dimension  xin(50), yin(50), zin(50)
      parameter ( a0=0.0d0, a1s2=0.5d0, a1=1.0d0, a2=2.0d0 )
      if(.not.igueq1) go to 130
      fctr1=eta(1,1,icons)
      go to 140
  130 fctr1=a1
  140 if(jgueq1) fctr1=fctr1*eta(1,1,jcons)
c
c  i primitive
c
      do 490 ig=1,igu
      if(.not.igueq1) call wzero(nc1,g1,1)
      ai=zet(ig,icons)
      arri=ai*rr
c
c  j primitive
c
      if(esfc) jgu=ig
      do 390 jg=1,jgu
      if(.not.jgueq1) call wzero(ij,gout,1)
      aj=zet(jg,jcons)
      aa=ai+aj
      aaa=(ai-aj)/aa
      dum=aj*arri/aa
      if(dum.gt.tol) go to 390
c
c  form kinetic energy integrals
c
      x0=xij+aaa*xijm
      y0=yij+aaa*yijm
      z0=zij+aaa*zijm
      t=sqrt(aa)
      fctr2=fctr1*exp(-dum)/(t*aa**3)
      if((esfc.and..not.igueq1).and.ig.eq.jg) fctr2=a1s2*fctr2
      dum=ai*aj
      tx=dum*dxij
      ty=dum*dyij
      tz=dum*dzij
      dum=dum*aa
      txx=dum-a2*tx*tx
      tyy=dum-a2*ty*ty
      tzz=dum-a2*tz*tz
      nij=0
      tii=a0
      ti=a0
      do 330 in=1,lit
      ni=in
      tjj=a0
      tj=a0
      do 320 jn=1,ljt
      nj=jn
      nij=nij+1
      call stvint
      xin(nij)=xint
      xin(nij+25)=txx*xint
      yin(nij)=yint
      yin(nij+25)=tyy*yint
      zin(nij)=zint
      zin(nij+25)=tzz*zint
      if(jn.eq.1) go to 210
      dum=tj+tj
      xin(nij+25)=xin(nij+25)+dum*tx*xin(nij-1)
      yin(nij+25)=yin(nij+25)+dum*ty*yin(nij-1)
      zin(nij+25)=zin(nij+25)+dum*tz*zin(nij-1)
      if(jn.eq.2) go to 210
      xin(nij+25)=xin(nij+25)+tjj*xin(nij-2)
      yin(nij+25)=yin(nij+25)+tjj*yin(nij-2)
      zin(nij+25)=zin(nij+25)+tjj*zin(nij-2)
  210 if(in.eq.1) go to 310
      dum=ti+ti
      xin(nij+25)=xin(nij+25)+dum*tx*xin(nij-5)
      yin(nij+25)=yin(nij+25)+dum*ty*yin(nij-5)
      zin(nij+25)=zin(nij+25)+dum*tz*zin(nij-5)
      if(jn.eq.1) go to 220
      dum=-tj*ti
      xin(nij+25)=xin(nij+25)+dum*xin(nij-6)
      yin(nij+25)=yin(nij+25)+dum*yin(nij-6)
      zin(nij+25)=zin(nij+25)+dum*zin(nij-6)
  220 if(in.eq.2) go to 310
      xin(nij+25)=xin(nij+25)+tii*xin(nij-10)
      yin(nij+25)=yin(nij+25)+tii*yin(nij-10)
      zin(nij+25)=zin(nij+25)+tii*zin(nij-10)
  310 tjj=tjj-tj*ai
  320 tj=tj+ai
      tii=tii+ti*aj
      ti=ti-aj
  330 nij=nij+(5-ljt)
      do 350 i=1,ij
      nx=ijx(i)
      ny=ijy(i)
      nz=ijz(i)
  350 gout(i)=gout(i)+fctr2*(xin(nx+25)*yin(ny)*zin(nz)+xin(nx)*
     1   yin(ny+25)*zin(nz)+xin(nx)*yin(ny)*zin(nz+25))
c
c  j transformation
c
      if(jgueq1) go to 400
      j1=0
      do 380 jrcr=1,jrcru
      do 370 i=1,ij
  370 g1(j1+i)=g1(j1+i)+gout(i)*eta(jrcr,jg,jcons)
  380 j1=j1+ij
  390 continue
c
c  i transformation
c
  400 if(igueq1) return
      if(esfc) i1=0
      ij2=0
      do 470 ircr=1,ircru
      j1=0
      if(.not.esfc) go to 440
      do 430 jrcr=1,ircr
      do 410 i=1,ij
  410 g2(ij2+i)=g2(ij2+i)+g1(j1+i)*eta(ircr,ig,icons)+
     1                    g1(i1+i)*eta(jrcr,ig,icons)
      ij2=ij2+ij
  430 j1=j1+ij
      i1=i1+ij
      go to 470
  440 do 460 jrcr=1,jrcru
      do 450 i=1,ij
  450 g2(ij2+i)=g2(ij2+i)+g1(j1+i)*eta(ircr,ig,icons)
      ij2=ij2+ij
  460 j1=j1+ij
  470 continue
  490 continue
      return
      end
