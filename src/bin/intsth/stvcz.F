      subroutine stvcz(a,chg,cx,eta,g2,nc,npair,
     1   nprir,nt,x,y,z,zet)
*mdc*on vax ibm cray
      implicit real*8 (a-h,o-z)
*mdc*off
      logical ec, esf, esfb, esfc, esfbc, igueq1, jgueq1, msfbc
cges 16jan90  this common was missing. zet found it.
      common/parmr/cutoff,tol
      COMMON /OUT/    IOUT,ICHECK,IFLAG3,NSECT,NINTS
      common /parmi/ mblu, mccu, mconu, mcu, mcxu,
     1   mrcru, mstu, msu, nblu, ng, ns, nst
      common /ntgr/ kaords, mconsu, mgcsu, mru, mcsu, mctu, mcru,
     1   msfu, mgu, msfru, mnru, mnsfu, mpru, ngcs, nu, mcxu2,
     2   lxyzir(3), inam(5), nnam, mdum(30)
      common /one/dxij,dyij,dzij,fnfct,rr,xij,xijm,yij,yijm,zij,zijm,
     1  ibl1,ibl2,icxi1,icxi2,ij,ijsf,ijx(225),ijy(225),ijz(225),ic,
     2  icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl,
     3  jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
      common/stv/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      common /dim21/ ipq(256)
      dimension  a(*),chg(*),cx(*),eta(mrcru,mconu,*),g2(*),
     1   nc(*),npair(2,*),
     2   nprir(2,mstu,*),
     3   nt(*),x(mcu,*),y(mcu,*),z(mcu,*),zet(mconu,*)
      dimension ix(35), iy(35), iz(35), jx(35), jy(35), jz(35)
c
      parameter ( a1s2=0.5d0 )
      data jx / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0,
     1          3, 0, 0, 2, 2, 1, 0, 1, 0, 1,
     2          4, 0, 0, 3, 3, 1, 0, 1, 0, 2,
     3          2, 0, 2, 1, 1/
      data ix / 1, 6, 1, 1,11, 1, 1, 6, 6, 1,
     1         16, 1, 1,11,11, 6, 1, 6, 1, 6,
     2         21, 1, 1,16,16, 6, 1, 6, 1,11,
     3         11, 1,11, 6, 6/
      data jy / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1,
     1          0, 3, 0, 1, 0, 2, 2, 0, 1, 1,
     2          0, 4, 0, 1, 0, 3, 3, 0, 1, 2,
     3          0, 2, 1, 2, 1/
      data iy / 1, 1, 6, 1, 1,11, 1, 6, 1, 6,
     1          1,16, 1, 6, 1,11,11, 1, 6, 6,
     2          1,21, 1, 6, 1,16,16, 1, 6,11,
     3          1,11, 6,11, 6/
      data jz / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1,
     1          0, 0, 3, 0, 1, 0, 1, 2, 2, 1,
     2          0, 0, 4, 0, 1, 0, 1, 3, 3, 0,
     3          2, 2, 1, 1, 2/
      data iz / 1, 1, 1, 6, 1, 1,11, 1, 6, 6,
     1          1, 1,16, 1, 6, 1, 6,11,11, 6,
     2          1, 1,21, 1, 6, 1, 6,16,16, 1,
     3         11,11, 6, 6,11/
      ec=ic.eq.jc
      esfc=esf.and.ec
      icx2=icxi2+((ic-1)*nc(js)+jc-1)*ntij2
      if(esfc.or.ircru.eq.1.or..not.esf) go to 12
      ibld=ibl2-ibl1
      icx3=icx2-(ic-jc)*(nc(is)-1)*ntij2
   12 if(esf) icx1=icxi1+ipq((ic-1)*nt(isf)+1)*nblt1+(jc-1)*ntij1
c
c  ishell
c
      xi=x(ic,is)
      yi=y(ic,is)
      zi=z(ic,is)
c
c  jshell
c
      xj=x(jc,js)
      yj=y(jc,js)
      zj=z(jc,js)
      xij=a1s2*(xi+xj)
      dxij=xi-xj
      xijm=a1s2*dxij
      yij=a1s2*(yi+yj)
      dyij=yi-yj
      yijm=a1s2*dyij
      zij=a1s2*(zi+zj)
      dzij=zi-zj
      zijm=a1s2*dzij
      rr=dxij*dxij+dyij*dyij+dzij*dzij
c
c  prepare items for pairs of (i,j) functions
c
      ij=0
      do 16 it=itl,itu
      nx=ix(it)
      ny=iy(it)
      nz=iz(it)
      if(esfc) jtu=it
      do 16 jt=jtl,jtu
      ij=ij+1
      ijx(ij)=nx+jx(jt)
      ijy(ij)=ny+jy(jt)
   16 ijz(ij)=nz+jz(jt)
      nc1=jrcru*ij
      if(.not.esfc) go to 18
      nc2=ipq(ircru+1)*ij
      go to 20
   18 nc2=ircru*nc1
   20 igueq1=igu.eq.1
      jgueq1=jgu.eq.1
      if(.not.igueq1) go to 22
      ic1=0
      go to 24
   22 ic1=nc2
   24 if(.not.jgueq1) go to 26
      ic0=ic1
      go to 28
   26 ic0=ic1+nc1
   28 mgout=ic0+ij
      if(mgout.gt.mblu) then
        mgout=mcxu+nblu+mgout
        write (iout,989) mgout
  989   format(' stvcz: mblu too small',i8)
        stop
      endif
      call wzero(nc2,g2,1)
      go to (36,38,40), nnam
   36 call sints(eta,g2,g2(ic1+1),g2(ic0+1),zet)
      go to 46
   38 call tints(eta,g2,g2(ic1+1),g2(ic0+1),zet)
      go to 46
   40 call vints(chg,eta,g2,g2(ic1+1),g2(ic0+1),nc,x,y,z,zet)
   46 ndx=0
      ndxa=0
      do 88 ircr=1,ircru
      if(esfc) jrcru=ircr
      do 88 jrcr=1,jrcru
      esfb=esf.and.ircr.eq.jrcr
      esfbc=esfb.and.ec
      msfbc=(.not.esfc).or.esfbc
      if(.not.esfb) go to 50
      icx=icx1
      go to 56
   50 if((.not.esf).or.esfc) go to 54
      if(ircr.gt.jrcr) go to 52
      ndxa=(ipq(jrcr)+(ircr-1))*ibl2-(jrcr-1)*ibld
      icx=icx3
      ngti=1
      ngtj=nt(isf)
      go to 60
   52 if(jrcr.eq.1) ndxa=(ipq(ircr)+(jrcr-1))*ibl2-(ircr-1)*ibld
   54 icx=icx2
   56 ngti=nt(jsf)
      ngtj=1
      if(.not.msfbc) icnt=0
   60 ndxi=ndx
      do 86 it=itl,itu
      ndxj=ndxi
      if(.not.msfbc) go to 64
      if(.not.esfbc) go to 62
      jtu=it
      go to 66
   62 ndxi=ndxi+ngti
      go to 68
   64 isw=icnt
      icnt=icnt+1
      istp=icnt
   66 ndxi=ndxi+(it-(jtl-1))
   68 do 86 jt=jtl,jtu
      if(msfbc.or.it.ge.jt) go to 70
      index=ndxj+isw
      istp=istp+1
      isw=isw+istp
      go to 72
   70 index=ndxj
      ndxj=ndxj+ngtj
   72 ndxb=ndxa
      if(esfb) then
        npr=nprir(1,1,ijsf)
        nprt=npair(1,ijsf)
      else
        npr=nprir(2,1,ijsf)
        nprt=npair(2,ijsf)
      endif
      if(npr.eq.0) go to 84
      val=g2(index+1)
      if(abs(val).le.cutoff) go to 80
      val=fnfct*val
      do 78 ibl=1,npr
   78 a(ndxb+ibl)=a(ndxb+ibl)+cx(icx+ibl)*val
   80 ndxb=ndxb+npr
   84 icx=icx+nprt
   86 continue
      ndxa=ndxb
   88 ndx=ndx+ij
      return
      end
