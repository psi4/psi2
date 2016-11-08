      subroutine genral(eta,g4,g3,g2,g1,gout,zet)
*mdc*on vax ibm cray
      implicit real*8 (a-h,o-z)
*mdc*off
      logical esfc, esfcij, esfckl, igueq1, jgueq1, kgueq1, lgueq1
      common /parmr/ cutoff, tol
      common /parmi/ mblu, mccu, mconu, mcu, mcxu,
     1   mrcru, mstu, msu, nblu, ng, ns, nst
      common /shlinf/ xi,yi,zi,rij,xk,yk,zk,rkl,igu,jgu,kgu,lgu,nc4,
     1   nc3,nc2,nc1,igueq1,jgueq1,kgueq1,lgueq1
      common /shlnos/lit,ljt,lkt,llt,ij,ijkl
ctph  common/root/xx,u(9),w(9),nroots
      common/root/xx,u(10),w(10),nroots
      common/misc/xij,xijm,yij,yijm,zij,zijm,xkl,xklm,ykl,yklm,zkl,zklm,
     1icons,ircru,jcons,jrcru,kcons,krcru,lcons,lrcru,esfc,esfcij,esfckl
      common /setint/ dxij,dyij,dzij,dxkl,dykl,dzkl,bp01,b00,b10,xcp00,
ctph 1  xc00,ycp00,yc00,zcp00,zc00,f00,in(9),kn(9),ni,nj,nk,nl,nmax,mmax
     1  xc00,ycp00,yc00,zcp00,zc00,f00,in(10),kn(10),ni,nj,nk,nl,nmax,
     2  mmax
      common/xyz/xint(5625),yint(5625),zint(5625)
      dimension eta(mrcru,mconu,*), g4(*), g3(*), g2(*), g1(*),
     1   gout(*), zet(mconu,*)
ctph  dimension in1(9)
      dimension in1(10)
      parameter ( a1s2=0.5d0, asrtpi=1.1283791670955126d0 )
      ni=lit-1
      nj=ljt-1
      nk=lkt-1
      nl=llt-1
      nmax=ni+nj
      mmax=nk+nl
      ihi=nmax+1
      do 12 i=1,ihi
      if(i.gt.lit) go to 10
      in1(i)=125*i-124
      go to 12
   10 in1(i)=(100*lit-124)+25*i
   12 continue
      ihi=mmax+1
      do 22 i=1,ihi
      if(i.gt.lkt) go to 20
      kn(i)=5*i-5
      go to 22
   20 kn(i)=(4*lkt-5)+i
   22 continue
      if(.not.igueq1) go to 30
      fctr1=asrtpi*eta(1,1,icons)
      go to 32
   30 fctr1=asrtpi
   32 if(jgueq1) fctr1=fctr1*eta(1,1,jcons)
      if(kgueq1) fctr1=fctr1*eta(1,1,kcons)
      if(lgueq1) fctr1=fctr1*eta(1,1,lcons)
      kgl=1
c     ----- i primitive -----
      do 550 ig=1,igu
      if(.not.igueq1) call wzero(nc3,g3,1)
      fctr2=fctr1
      ai=zet(ig,icons)
      arij=ai*rij
      if(esfc) kgl=ig
c     ----- j primitive -----
      if(esfcij) jgu=ig
      do 350 jg=1,jgu
      if(.not.jgueq1) call wzero(nc2,g2,1)
      if((esfcij.and..not.igueq1).and.jg.eq.ig) fctr2=a1s2*fctr2
      aj=zet(jg,jcons)
      aa=ai+aj
      aarij=aj*arij/aa
      if(aarij.gt.tol) go to 350
      fctr3=fctr2/aa
      aaa=(ai-aj)/aa
      xa=xij+aaa*xijm
      ya=yij+aaa*yijm
      za=zij+aaa*zijm
      axai=aa*(xa-xi)
      ayai=aa*(ya-yi)
      azai=aa*(za-zi)
      axak=aa*(xa-xk)
      ayak=aa*(ya-yk)
      azak=aa*(za-zk)
c     ----- k primitive -----
      do 250 kg=kgl,kgu
      if(.not.kgueq1) call wzero(nc1,g1,1)
      bk=zet(kg,kcons)
      brkl=bk*rkl
c     ----- l primitive -----
      if(.not.(esfc.and.(kg.eq.ig))) go to 36
      lgl=jg
      go to 37
   36 lgl=1
   37 if(esfckl) lgu=kg
      do 150 lg=lgl,lgu
      if(.not.lgueq1) call wzero(ijkl,gout,1)
      bl=zet(lg,lcons)
      bb=bk+bl
      dum=aarij+bl*brkl/bb
      if(dum.gt.tol) go to 150
      apb=aa+bb
      expe=fctr3*exp(-dum)/(bb*sqrt(apb))
      if(((esfc.and..not.(igueq1.and.jgueq1)).and.(kg.eq.ig)).and.
     1   lg.eq.jg) expe=a1s2*expe
      if((esfckl.and..not.kgueq1).and.lg.eq.kg) expe=a1s2*expe
      bbb=(bk-bl)/bb
      xb=xkl+bbb*xklm
      yb=ykl+bbb*yklm
      zb=zkl+bbb*zklm
      atb=aa*bb
      rho=atb/apb
      xx=rho*((xa-xb)**2+(ya-yb)**2+(za-zb)**2)
      bxbi=bb*(xb-xi)
      bybi=bb*(yb-yi)
      bzbi=bb*(zb-zi)
      bxbk=bb*(xb-xk)
      bybk=bb*(yb-yk)
      bzbk=bb*(zb-zk)
      c1x=bxbk+axak
      c2x=aa*bxbk
      c3x=bxbi+axai
      c4x=bb*axai
      c1y=bybk+ayak
      c2y=aa*bybk
      c3y=bybi+ayai
      c4y=bb*ayai
      c1z=bzbk+azak
      c2z=aa*bzbk
      c3z=bzbi+azai
      c4z=bb*azai
c     ----- roots and weights for quadrature -----
      if(nroots-4) 88,86,82
   82 if(nroots.eq.5) go to 84
      call droot
      go to 90
   84 call root5
      go to 90
   86 call root4
      go to 90
   88 call rt123
c     ----- compute two-electron  integrals for each root -----
   90 mm=0
      ihi=nmax+1
      do 120 m=1,nroots
      u2=u(m)*rho
      f00=expe*w(m)
      do 100 i=1,ihi
  100 in(i)=in1(i)+mm
      dum=atb+u2*apb
      dum2=dum+dum
      bp01=(aa+u2)/dum2
      b00=u2/dum2
      b10=(bb+u2)/dum2
      xcp00=(u2*c1x+c2x)/dum
      xc00 =(u2*c3x+c4x)/dum
      ycp00=(u2*c1y+c2y)/dum
      yc00 =(u2*c3y+c4y)/dum
      zcp00=(u2*c1z+c2z)/dum
      zc00 =(u2*c3z+c4z)/dum
      call xyzint
  120 mm=mm+625
c     ----- form (i,j//k,l) integrals over functions -----
      call forms(gout)
c     ----- l transformation -----
      if(lgueq1) go to 200
      l1=0
      do 140 lrcr=1,lrcru
      do 135 n=1,ijkl
  135 g1(l1+n)=g1(l1+n)+gout(n)*eta(lrcr,lg,lcons)
  140 l1=l1+ijkl
  150 continue
c     ----- k transformation -----
  200 if(kgueq1) go to 280
      if(esfckl) k1=0
      kl2=0
      do 240 krcr=1,krcru
      l1=0
      if(.not.esfckl) go to 220
      do 210 lrcr=1,krcr
      do 205 n=1,ijkl
  205 g2(kl2+n)=g2(kl2+n)+g1(l1+n)*eta(krcr,kg,kcons)+
     1                    g1(k1+n)*eta(lrcr,kg,kcons)
      kl2=kl2+ijkl
  210 l1=l1+ijkl
      k1=k1+ijkl
      go to 240
  220 do 230 lrcr=1,lrcru
      do 225 n=1,ijkl
  225 g2(kl2+n)=g2(kl2+n)+g1(l1+n)*eta(krcr,kg,kcons)
      kl2=kl2+ijkl
  230 l1=l1+ijkl
  240 continue
  250 continue
c     ----- j transformation -----
  280 if(jgueq1) go to 410
      jkl3=0
      if(.not.(esfc.and.igueq1)) go to 320
      ij2=0
      do 310 jrcr=1,jrcru
      kl2=0
      do 300 lrcr=1,jrcr
      do 300 nij=1,ij
      ij2n=ij2+nij
      do 290 nkl=1,ij
      g3(jkl3+nkl)=g3(jkl3+nkl)+g2(kl2+nkl)*eta(jrcr,jg,jcons)+
     1                          g2(ij2n)*eta(lrcr,jg,jcons)
  290 ij2n=ij2n+ij
      jkl3=jkl3+ij
  300 kl2=kl2+ij
  310 ij2=ij2+ijkl
      go to 350
  320 do 340 jrcr=1,jrcru
      kl2=0
      do 340 krcr=1,krcru
      if(esfckl) lrcru=krcr
      do 340 lrcr=1,lrcru
      do 325 n=1,ijkl
  325 g3(jkl3+n)=g3(jkl3+n)+g2(kl2+n)*eta(jrcr,jg,jcons)
      jkl3=jkl3+ijkl
  340 kl2=kl2+ijkl
  350 continue
c     ----- i transformation -----
  410 if(igueq1) return
      if(esfc) ij3=0
      if(esfcij) ikl3i=0
      ijkl4=0
      do 550 ircr=1,ircru
      jkl3j=0
      if(esfcij) jrcru=ircr
      do 540 jrcr=1,jrcru
      if(esfc) jkl3k=jkl3j
      if(esfcij) ikl3k=ikl3i
      if(esfc.and.esfcij) kij3=ij3
      if(esfc) krcru=ircr
      do 520 krcr=1,krcru
      if(.not.esfc) go to 480
      jkl3=jkl3k
      lij3=ij3
      if(krcr.ne.ircr) go to 442
      lrcru=jrcr
      if(.not.esfckl) go to 460
      go to 444
  442 if(.not.esfckl) go to 455
      lrcru=krcr
  444 ikl3=ikl3k
      do 450 lrcr=1,lrcru
      do 445 n=1,ijkl
  445 g4(ijkl4+n)=g4(ijkl4+n)+g3(jkl3+n)*eta(ircr,ig,icons)+
     1                        g3(ikl3+n)*eta(jrcr,ig,icons)+
     2                        g3(lij3+n)*eta(krcr,ig,icons)+
     3                        g3(kij3+n)*eta(lrcr,ig,icons)
      ijkl4=ijkl4+ijkl
      jkl3=jkl3+ijkl
      ikl3=ikl3+ijkl
  450 lij3=lij3+nc2
      kij3=kij3+nc2
      jkl3k=jkl3k+(krcr*ijkl)
      ikl3k=ikl3k+(krcr*ijkl)
      go to 520
  455 lrcru=jrcru
  460 do 470 lrcr=1,lrcru
      do 465 nij=1,ij
      lij3n=lij3+nij
      do 463 nkl=1,ij
      g4(ijkl4+nkl)=g4(ijkl4+nkl)+g3(jkl3+nkl)*eta(ircr,ig,icons)+
     1                            g3(lij3n)*eta(krcr,ig,icons)
  463 lij3n=lij3n+ij
      ijkl4=ijkl4+ij
  465 jkl3=jkl3+ij
  470 lij3=lij3+nc2
      jkl3k=jkl3k+nc1
      go to 520
  480 if(esfckl) lrcru=krcr
      if(.not.esfcij) go to 500
      do 490 lrcr=1,lrcru
      do 485 n=1,ijkl
  485 g4(ijkl4+n)=g4(ijkl4+n)+g3(jkl3j+n)*eta(ircr,ig,icons)+
     1                        g3(ikl3k+n)*eta(jrcr,ig,icons)
      ijkl4=ijkl4+ijkl
      jkl3j=jkl3j+ijkl
  490 ikl3k=ikl3k+ijkl
      go to 520
  500 do 510 lrcr=1,lrcru
      do 505 n=1,ijkl
  505 g4(ijkl4+n)=g4(ijkl4+n)+g3(jkl3j+n)*eta(ircr,ig,icons)
      ijkl4=ijkl4+ijkl
  510 jkl3j=jkl3j+ijkl
  520 continue
      if(esfc) jkl3j=jkl3j+nc2
  540 if(esfc) ij3=ij3+ijkl
  550 if(esfcij) ikl3i=ikl3i+nc2
      return
      end
