      subroutine vints(chg,eta,g2,g1,gout,nc,x,y,z,zet)
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
ctph  common/root/xx,u(9),w(9),nroots
      common/root/xx,u(10),w(10),nroots
      common/stv/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      dimension chg(*), eta(mrcru,mconu,*), g2(*), g1(*), gout(*),
     1   nc(*), x(mcu,*), y(mcu,*), z(mcu,*), zet(mconu,*)
      dimension xin(125), yin(125), zin(125)
      parameter ( a1s2=0.5d0, asrtpi=1.1283791670955126d0 )
      nroots=(lit+ljt)/2
      mmu=25*nroots
      if(.not.igueq1) go to 130
      fctr1=asrtpi*eta(1,1,icons)
      go to 140
  130 fctr1=asrtpi
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
c  form nuclear attraction integrals
c
      ax=xij+aaa*xijm
      ay=yij+aaa*yijm
      az=zij+aaa*zijm
      fctr2=fctr1*exp(-dum)/aa
      if((esfc.and..not.igueq1).and.ig.eq.jg) fctr2=a1s2*fctr2
      aax=aa*ax
      aay=aa*ay
      aaz=aa*az
      do 350 ks=1,ns
      fctr3=-chg(ks)*fctr2
      kcu=nc(ks)
      do 350 kc=1,kcu
      cx=x(kc,ks)
      cy=y(kc,ks)
      cz=z(kc,ks)
      xx=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
      if(nroots-4) 230,220,210
  210 call root5
      go to 290
  220 call root4
      go to 290
  230 call rt123
  290 mm=0
      do 340 k=1,nroots
      uu=aa*u(k)
      ww=w(k)*fctr3
      tt=aa+uu
      t=sqrt(tt)
      x0=(aax+uu*cx)/tt
      y0=(aay+uu*cy)/tt
      z0=(aaz+uu*cz)/tt
      do 330 in=1,lit
      ni=in
      do 320 jn=1,ljt
      nj=jn
      mm=mm+1
      call stvint
      xin(mm)=xint
      yin(mm)=yint
      zin(mm)=zint*ww
  320 continue
  330 mm=mm+(5-ljt)
  340 mm=mm+(5*(5-lit))
      do 350 i=1,ij
      nx=ijx(i)-25
      ny=ijy(i)-25
      nz=ijz(i)-25
      do 350 mm=25,mmu,25
  350 gout(i)=gout(i)+xin(nx+mm)*yin(ny+mm)*zin(nz+mm)
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
