      subroutine sf3eq6(a,cx,eta,g4,icb,icc,iscmci,is,isf,js,jsf,ls,
     1   lsf,lmnp1,mcons,nc,ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
*mdc*on vax ibm cray
      implicit real*8 (a-h,o-z)
*mdc*off
      logical jceqkc
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
      dimension icxbg(3),jbl(2), kcxbg(3,3)
      ic=icb(inx,1,icc)
      jc=icb(jnx,1,icc)
      kc=icb(knx,1,icc)
      lc=icb(lnx,1,icc)
      jceqkc=jc.eq.kc
      jesfbl=(5-iscmci)/2
      do 10 jesfb=jesfbl,2
      jbl(jesfb)=0
      do 10 ist=1,nst
   10 jbl(jesfb)=jbl(jesfb)+npriri(2,iscmci,ist)*
     1   nprirk(jesfb,iscmci,ist)
      if(iscmci.eq.3) go to 64
      icxbg(1)=ijcxst(2)+((ic-1)*nc(js)+(jc-1))*ntij
      icxbg(2)=icxbg(1)+(kc-jc)*ntij
      kcxbg(2,1)=klcxst(2)+((kc-1)*nc(ls)+(lc-1))*ntkl
      kcxbg(2,2)=kcxbg(2,1)+(jc-kc)*nc(ls)*ntkl
      if(.not.jceqkc) go to 12
      fnfct=(nwt*nfct(icc))
      go to 28
   12 if(ngw(icc).eq.1) go to 20
      igwu = ngw(icc)
      do 16 igw=2,igwu
   16 if(icb(inx,igw,icc).eq.ic.and.icb(jnx,igw,icc).eq.kc.and.
     1   icb(lnx,igw,icc).eq.lc) go to 24
   20 iscml=1
      fnfct=(nwt*nfct(icc))
      go to 32
   24 fnfct=(2*nwt*nfct(icc))
   28 iscml=2
   32 do 60 iscm=iscml,2
      if(iscm.eq.2) go to 36
      call jandk(eta,g4,is,isf,ic,js,jsf,jc,js,jsf,kc,ls,lsf,lc,
     1   lmnp1,mcons,ncon,nrcr,nt,ntl,ntu,x,y,z,zet)
      go to 40
   36 call jandk(eta,g4,is,isf,ic,js,jsf,kc,js,jsf,jc,ls,lsf,lc,
     1   lmnp1,mcons,ncon,nrcr,nt,ntl,ntu,x,y,z,zet)
   40 ndx=0
      iblb=0
      do 60 ircr=1,ircru
      do 60 jrcr=1,jrcru
      do 60 krcr=1,krcru
      do 60 lrcr=1,lrcru
      icxs(2,iscm)=icxbg(iscm)
      do 48 it=itl,itu
      do 48 jt=jtl,jtu
      kcxs(2,iscm)=kcxbg(2,iscm)
      do 44 kt=ktl,ktu
      do 44 lt=ltl,ltu
      ndx=ndx+1
      val=g4(ndx)
      if(abs(val).le.cutoff) go to 44
      val=fnfct*val
      ibl=iblb
      call aoso2(a,cx,iscm,val)
   44 kcxs(2,iscm)=kcxs(2,iscm)+nprk(2,1)
   48 icxs(2,iscm)=icxs(2,iscm)+npri(2,1)
      iblb=iblb+jbl(2)
   60 continue
      return
   64 fnfct=(nwt*nfct(icc))
      jbld=jbl(2)-jbl(1)
      jbljk=ipq(jrcru)*jbl(2)+jrcru*jbl(1)
      icxbg(3)=ilcxst(2)+((ic-1)*nc(ls)+(lc-1))*ntil
      kcxbg(1,3)=jkcxst(1)+ipq((jc-1)*nt(jsf)+1)*nprk(1,3)+(kc-1)*ntjk1
      kcxbg(2,3)=jkcxst(2)+((jc-1)*nc(js)+(kc-1))*ntjk
      kcxbg(3,3)=kcxbg(2,3)-(jc-kc)*(nc(js)-1)*ntjk
      call jandk(eta,g4,is,isf,ic,ls,lsf,lc,js,jsf,jc,js,jsf,kc,
     1   lmnp1,mcons,ncon,nrcr,nt,ntl,ntu,x,y,z,zet)
      ndxl=1
      ngtl=nt(jsf)**2
      iblbil=0
      do 96 ircr=1,ircru
      do 96 lrcr=1,lrcru
      do 92 jrcr=1,jrcru
      if(jceqkc) krcru=jrcr
      do 92 krcr=1,krcru
      if(jrcr-krcr) 67,69,71
   67 iblb=iblbil+(ipq(krcr)+(jrcr-1))*jbl(2)-(krcr-1)*jbld
      kesfb=2
      kcxbg1=kcxbg(3,3)
      ngtj=1
      ngtk=nt(jsf)
      go to 73
   69 iblb=iblbil+ipq(krcr)*jbl(2)+(krcr-1)*jbl(1)
      kesfb=1
      go to 72
   71 iblb=iblbil+(ipq(jrcr)+(krcr-1))*jbl(2)-(jrcr-1)*jbld
      kesfb=2
   72 kcxbg1=kcxbg(kesfb,3)
      ngtj=nt(jsf)
      ngtk=1
   73 icxs(2,3)=icxbg(3)
      do 84 it=itl,itu
      do 84 lt=ltl,ltu
      kcxs(kesfb,3)=kcxbg1
      ndxj=ndxl
      do 80 jt=jtl,jtu
      if(jrcr.eq.krcr.and.jceqkc) ktu=jt
      ndxk=ndxj
      do 76 kt=ktl,ktu
      val=g4(ndxk)
      if(abs(val).le.cutoff) go to 74
      val=fnfct*val
      ibl=iblb
      call aoso2(a,cx,3,val)
   74 ndxk=ndxk+ngtk
   76 kcxs(kesfb,3)=kcxs(kesfb,3)+nprk(kesfb,3)
   80 ndxj=ndxj+ngtj
      ndxl=ndxl+ngtl
   84 icxs(2,3)=icxs(2,3)+npri(2,3)
   92 continue
   96 iblbil=iblbil+jbljk
      return
      end
