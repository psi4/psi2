      subroutine sf3eq7(a,cx,eta,g4,icb,icc,iscmci,is,isf,js,jsf,ks,
     1   ksf,lmnp1,mcons,nc,ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
*mdc*on vax ibm cray
      implicit real*8 (a-h,o-z)
*mdc*off
      logical kceqlc
      common /parmr/ cutoff, tol
      common /parmi/ mblu, mccu, mconu, mcu, mcxu,
     1   mrcru, mstu, msu, nblu, ng, ns, nst
      common /ntgr/itl,itu,jtl,jtu,ktl,ktu,ltl,ltu,jcu,kcu,lcu,inx,
     1jnx,knx,lnx,nwt,nwt1(3),ijsf(3),klsf(3),icxs(2,3),kcxs(2,3),
     2npri(2,3),nprk(2,3),iesfb,kesfb,ircru,jrcru,krcru,lrcru
      common /cxindx/ ntij,ntik,ntil,ntkl,ntjl,ntjk,ntij1,ntjk1,ntkl1,
     1   ijcxst(2),ikcxst(2),ilcxst(2),klcxst(2),jlcxst(2),jkcxst(2)
      common /dim21/ ipq(256)
      common /ikpr/ npriri(2,3,8), nprirk(2,3,8)
      common /sf3ao/ ibl
      dimension a(*),cx(*),eta(mrcru,mconu,*),g4(*),icb(4,24,*),
     1   lmnp1(*),mcons(*),nc(*),ncon(*),nfct(*),ngw(*),nrcr(*),
     2   nt(*),ntl(*),ntu(*),x(mcu,*),y(mcu,*),z(mcu,*),zet(mconu,*)
      dimension icxbg(3), jbl(2), kcxbg(3,3)
      ic=icb(inx,1,icc)
      jc=icb(jnx,1,icc)
      kc=icb(knx,1,icc)
      lc=icb(lnx,1,icc)
      kceqlc=kc.eq.lc
      do 10 jesfb=iscmci,2
      jbl(jesfb)=0
      do 10 ist=1,nst
   10 jbl(jesfb)=jbl(jesfb)+npriri(2,iscmci,ist)*
     1   nprirk(jesfb,iscmci,ist)
      if(iscmci.eq.2) go to 44
      jbld=jbl(2)-jbl(1)
      jblkl=ipq(krcru)*jbl(2)+krcru*jbl(1)
      fnfct=(nwt*nfct(icc))
      icxbg(1)=ijcxst(2)+((ic-1)*nc(js)+(jc-1))*ntij
      kcxbg(1,1)=klcxst(1)+ipq((kc-1)*nt(ksf)+1)*nprk(1,1)+(lc-1)*ntkl1
      kcxbg(2,1)=klcxst(2)+((kc-1)*nc(ks)+(lc-1))*ntkl
      kcxbg(3,1)=kcxbg(2,1)-(kc-lc)*(nc(ks)-1)*ntkl
      call jandk(eta,g4,is,isf,ic,js,jsf,jc,ks,ksf,kc,ks,ksf,lc,
     1   lmnp1,mcons,ncon,nrcr,nt,ntl,ntu,x,y,z,zet)
      ndxj=1
      ngtj=nt(ksf)**2
      iblbij=0
      do 40 ircr=1,ircru
      do 40 jrcr=1,jrcru
      do 36 krcr=1,krcru
      if(kceqlc) lrcru=krcr
      do 36 lrcr=1,lrcru
      if(krcr-lrcr) 12,13,14
   12 iblb=iblbij+(ipq(lrcr)+(krcr-1))*jbl(2)-(lrcr-1)*jbld
      kesfb=2
      kcxbg1=kcxbg(3,1)
      ngtk=1
      ngtl=nt(ksf)
      go to 16
   13 iblb=iblbij+ipq(krcr)*jbl(2)+(krcr-1)*jbl(1)
      kesfb=1
      go to 15
   14 iblb=iblbij+(ipq(krcr)+(lrcr-1))*jbl(2)-(krcr-1)*jbld
      kesfb=2
   15 kcxbg1=kcxbg(kesfb,1)
      ngtk=nt(ksf)
      ngtl=1
   16 icxs(2,1)=icxbg(1)
      do 28 it=itl,itu
      do 28 jt=jtl,jtu
      kcxs(kesfb,1)=kcxbg1
      ndxk=ndxj
      do 24 kt=ktl,ktu
      if(krcr.eq.lrcr.and.kceqlc) ltu=kt
      ndxl=ndxk
      do 20 lt=ltl,ltu
      val=g4(ndxl)
      if(abs(val).le.cutoff) go to 18
      val=fnfct*val
      ibl=iblb
      call aoso2(a,cx,1,val)
   18 ndxl=ndxl+ngtl
   20 kcxs(kesfb,1)=kcxs(kesfb,1)+nprk(kesfb,1)
   24 ndxk=ndxk+ngtk
      ndxj=ndxj+ngtj
   28 icxs(2,1)=icxs(2,1)+npri(2,1)
   36 continue
   40 iblbij=iblbij+jblkl
      return
   44 icxbg(2)=ikcxst(2)+((ic-1)*nc(ks)+(kc-1))*ntik
      icxbg(3)=icxbg(2)+(lc-kc)*ntik
      kcxbg(2,3)=jkcxst(2)+((jc-1)*nc(ks)+(kc-1))*ntjk
      kcxbg(2,2)=kcxbg(2,3)+(lc-kc)*ntjk
      if(.not.kceqlc) go to 48
      fnfct=(nwt*nfct(icc))
      go to 64
   48 if(ngw(icc).eq.1) go to 56
      igwu = ngw(icc)
      do 52 igw=2,igwu
   52 if(icb(inx,igw,icc).eq.ic.and.icb(jnx,igw,icc).eq.jc.and.
     1   icb(knx,igw,icc).eq.lc) go to 60
   56 iscmu=3
      fnfct=(nwt*nfct(icc))
      go to 68
   60 fnfct=(2*nwt*nfct(icc))
   64 iscmu=2
   68 do 98 iscm=2,iscmu
      if(iscm.eq.3) go to 72
      call jandk(eta,g4,is,isf,ic,ks,ksf,kc,js,jsf,jc,ks,ksf,lc,
     1   lmnp1,mcons,ncon,nrcr,nt,ntl,ntu,x,y,z,zet)
      go to 76
   72 call jandk(eta,g4,is,isf,ic,ks,ksf,lc,js,jsf,jc,ks,ksf,kc,
     1   lmnp1,mcons,ncon,nrcr,nt,ntl,ntu,x,y,z,zet)
   76 ndx=0
      iblb=0
      do 98 ircr=1,ircru
      do 98 krcr=1,krcru
      do 98 jrcr=1,jrcru
      do 98 lrcr=1,lrcru
      icxs(2,iscm)=icxbg(iscm)
      do 88 it=itl,itu
      do 88 kt=ktl,ktu
      kcxs(2,iscm)=kcxbg(2,iscm)
      do 84 jt=jtl,jtu
      do 80 lt=ltl,ltu
      ndx=ndx+1
      val=g4(ndx)
      if(abs(val).le.cutoff) go to 80
      val=fnfct*val
      ibl=iblb
      call aoso2(a,cx,iscm,val)
   80 kcxs(2,iscm)=kcxs(2,iscm)+nprk(2,2)
   84 continue
   88 icxs(2,iscm)=icxs(2,iscm)+npri(2,2)
   98 iblb=iblb+jbl(2)
      return
      end
