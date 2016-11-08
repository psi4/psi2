      subroutine s0000(eta,g4,g3,g2,g1,gout,zet)
*mdc*on vax ibm cray
      implicit real*8 (a-h,o-z)
*mdc*off
      logical esfc, esfcij, esfckl, igueq1, jgueq1, kgueq1, lgueq1
      common /parmr/ cutoff, tol
      common /parmi/ mblu, mccu, mconu, mcu, mcxu,
     1   mrcru, mstu, msu, nblu, ng, ns, nst
      common /shlinf/ xi,yi,zi,rij,xk,yk,zk,rkl,igu,jgu,kgu,lgu,nc4,
     1   nc3,nc2,nc1,igueq1,jgueq1,kgueq1,lgueq1
      common/misc/xij,xijm,yij,yijm,zij,zijm,xkl,xklm,ykl,yklm,zkl,zklm,
     1icons,ircru,jcons,jrcru,kcons,krcru,lcons,lrcru,esfc,esfcij,esfckl
      dimension eta(mrcru,mconu,*), g4(*), g3(*), g2(*), g1(*),
     1   gout(*), zet(mconu,*)
      parameter ( a1s2=0.5d0, asrtpi=1.1283791670955126d0 )
      parameter ( pie4=7.85398163397448d-01 )
      abx=xij-xkl
      aby=yij-ykl
      abz=zij-zkl
      if(.not.igueq1) go to 10
      fctr1=asrtpi*eta(1,1,icons)
      go to 12
   10 fctr1=asrtpi
   12 if(jgueq1) fctr1=fctr1*eta(1,1,jcons)
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
      px=abx+aaa*xijm
      py=aby+aaa*yijm
      pz=abz+aaa*zijm
c     ----- k primitive -----
      do 250 kg=kgl,kgu
      if(.not.kgueq1) call wzero(nc1,g1,1)
      bk=zet(kg,kcons)
      brkl=bk*rkl
c     ----- l primitive -----
      if(.not.(esfc.and.(kg.eq.ig))) go to 16
      lgl=jg
      go to 17
   16 lgl=1
   17 if(esfckl) lgu=kg
      do 150 lg=lgl,lgu
      bl=zet(lg,lcons)
      bb=bk+bl
      dum=aarij+bl*brkl/bb
      if(dum.gt.tol) go to 150
      expe=fctr3*exp(-dum)/bb
      if(((esfc.and..not.(igueq1.and.jgueq1)).and.(kg.eq.ig)).and.
     1   lg.eq.jg) expe=a1s2*expe
      if((esfckl.and..not.kgueq1).and.lg.eq.kg) expe=a1s2*expe
      bbb=(bk-bl)/bb
      apb=aa+bb
      x=((px-bbb*xklm)**2+(py-bbb*yklm)**2+(pz-bbb*zklm)**2)*aa*bb/apb
      if(x.gt.5.0d+00) go to 50
      if(x.gt.1.0d+00) go to 30
      if(x.gt.3.0d-7) go to 20
      ww1=1.0d+00-x/3.0d+00
      go to 100
   20 f1=          ((((((((-8.36313918003957d-08*x+1.21222603512827d-06
     1)*x-1.15662609053481d-05 )*x+9.25197374512647d-05
     2)*x-6.40994113129432d-04 )*x+3.78787044215009d-03
     3)*x-1.85185172458485d-02 )*x+7.14285713298222d-02
     4)*x-1.99999999997023d-01 )*x+3.33333333333318d-01
      ww1=(x+x)*f1+exp(-x)
      go to 100
   30 if(x.gt.3.0d+00) go to 40
      y=x-2.0d+00
      f1=        ((((((((((-1.61702782425558d-10*y+1.96215250865776d-09
     1)*y-2.14234468198419d-08 )*y+2.17216556336318d-07
     2)*y-1.98850171329371d-06 )*y+1.62429321438911d-05
     3)*y-1.16740298039895d-04 )*y+7.24888732052332d-04
     4)*y-3.79490003707156d-03 )*y+1.61723488664661d-02
     5)*y-5.29428148329736d-02 )*y+1.15702180856167d-01
      ww1=(x+x)*f1+exp(-x)
      go to 100
   40 y=x-4.0d+00
      f1=        ((((((((((-2.62453564772299d-11*y+3.24031041623823d-10
     1)*y-3.614965656163d-09)*y+3.760256799971d-08)*y-3.553558319675d-07
     2)*y+3.022556449731d-06)*y-2.290098979647d-05)*y+1.526537461148d-04
     3)*y-8.81947375894379d-04 )*y+4.33207949514611d-03
     4)*y-1.75257821619926d-02 )*y+5.28406320615584d-02
      ww1=(x+x)*f1+exp(-x)
      go to 100
   50 if(x.gt.15.0d+00) go to 70
      e=exp(-x)
      if(x.gt.10.0d+00) go to 60
      ww1=    (((((( 4.6897511375022d-01/x-6.9955602298985d-01)/x
     1+5.3689283271887d-01)/x-3.2883030418398d-01)/x
     2+2.4645596956002d-01)/x-4.9984072848436d-01)/x
     3-3.1501078774085d-06)*e + sqrt(pie4/x)
      go to 100
   60 ww1=       (((-1.8784686463512d-01/x+2.2991849164985d-01)/x
     1-4.9893752514047d-01)/x-2.1916512131607d-05)*e + sqrt(pie4/x)
      go to 100
   70 if(x.gt.33.0d+00) go to 90
      e=exp(-x)
      ww1=        (( 1.9623264149430d-01/x-4.9695241464490d-01)/x
     1-6.0156581186481d-05)*e + sqrt(pie4/x)
      go to 100
   90 ww1=sqrt(pie4/x)
  100 gout(1)=ww1*expe/sqrt(apb)
c     ----- l transformation -----
      if(lgueq1) go to 200
      do 140 lrcr=1,lrcru
  140 g1(lrcr)=g1(lrcr)+gout(1)*eta(lrcr,lg,lcons)
  150 continue
c     ----- k transformation -----
  200 if(kgueq1) go to 280
      kl2=0
      do 240 krcr=1,krcru
      if(.not.esfckl) go to 220
      do 210 lrcr=1,krcr
  210 g2(kl2+lrcr)=g2(kl2+lrcr)+g1(lrcr)*eta(krcr,kg,kcons)+
     1                          g1(krcr)*eta(lrcr,kg,kcons)
      kl2=kl2+krcr
      go to 240
  220 do 230 lrcr=1,lrcru
  230 g2(kl2+lrcr)=g2(kl2+lrcr)+g1(lrcr)*eta(krcr,kg,kcons)
      kl2=kl2+lrcru
  240 continue
  250 continue
c     ----- j transformation -----
  280 if(jgueq1) go to 410
      jkl3=0
      if(.not.(esfc.and.igueq1)) go to 320
      do 310 jrcr=1,jrcru
      do 300 lrcr=1,jrcr
  300 g3(jkl3+lrcr)=g3(jkl3+lrcr)+g2(lrcr)*eta(jrcr,jg,jcons)+
     1                            g2(jrcr)*eta(lrcr,jg,jcons)
  310 jkl3=jkl3+jrcr
      go to 350
  320 do 340 jrcr=1,jrcru
      kl2=0
      do 340 krcr=1,krcru
      if(esfckl) lrcru=krcr
      do 330 lrcr=1,lrcru
  330 g3(jkl3+lrcr)=g3(jkl3+lrcr)+g2(kl2+lrcr)*eta(jrcr,jg,jcons)
      jkl3=jkl3+lrcru
  340 kl2=kl2+lrcru
  350 continue
c     ----- i transformation -----
  410 if(igueq1) return
      if(esfc) ij3=-nc2+1
      if(esfcij) ikl3i=0
      ijkl4=0
      do 550 ircr=1,ircru
      jkl3j=0
      if(esfcij) jrcru=ircr
      do 540 jrcr=1,jrcru
      jkl3=jkl3j
      if(esfcij) ikl3=ikl3i
      if(esfc) krcru=ircr
      do 520 krcr=1,krcru
      if(.not.esfc) go to 480
      if(krcr.ne.ircr) go to 442
      lrcru=jrcr
      if(.not.esfckl) go to 460
      go to 444
  442 if(.not.esfckl) go to 455
      lrcru=krcr
  444 do 450 lrcr=1,lrcru
  450 g4(ijkl4+lrcr)=g4(ijkl4+lrcr)+g3(jkl3+lrcr)*eta(ircr,ig,icons)+
     1                              g3(ikl3+lrcr)*eta(jrcr,ig,icons)+
     2                              g3(ij3+nc2*lrcr)*eta(krcr,ig,icons)+
     3                              g3(ij3+nc2*krcr)*eta(lrcr,ig,icons)
      ijkl4=ijkl4+lrcru
      jkl3=jkl3+krcr
      ikl3=ikl3+krcr
      go to 520
  455 lrcru=jrcru
  460 do 470 lrcr=1,lrcru
  470 g4(ijkl4+lrcr)=g4(ijkl4+lrcr)+g3(jkl3+lrcr)*eta(ircr,ig,icons)+
     1                              g3(ij3+nc2*lrcr)*eta(krcr,ig,icons)
      ijkl4=ijkl4+lrcru
      jkl3=jkl3+jrcru
      go to 520
  480 if(esfckl) lrcru=krcr
      if(.not.esfcij) go to 500
      do 490 lrcr=1,lrcru
  490 g4(ijkl4+lrcr)=g4(ijkl4+lrcr)+g3(jkl3+lrcr)*eta(ircr,ig,icons)+
     1                              g3(ikl3+lrcr)*eta(jrcr,ig,icons)
      ijkl4=ijkl4+lrcru
      jkl3=jkl3+lrcru
      ikl3=ikl3+lrcru
      go to 520
  500 do 510 lrcr=1,lrcru
  510 g4(ijkl4+lrcr)=g4(ijkl4+lrcr)+g3(jkl3+lrcr)*eta(ircr,ig,icons)
      ijkl4=ijkl4+lrcru
      jkl3=jkl3+lrcru
  520 continue
      jkl3j=jkl3j+nc2
  540 if(esfc) ij3=ij3+1
  550 if(esfcij) ikl3i=ikl3i+nc2
      return
      end
