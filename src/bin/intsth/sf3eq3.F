      subroutine sf3eq3(a,cx,eta,g4,icb,icc,iscmci,is,isf,ks,ksf,
     1   lmnp1,mcons,nc,ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
*mdc*on vax ibm cray
      implicit real*8 (a-h,o-z)
*mdc*off
      logical iceqjc, kceqlc
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
      dimension icxbg(3,3), jbl(2,2), jbld(2), kcxbg(3,3)
      ic=icb(inx,1,icc)
      jc=icb(jnx,1,icc)
      kc=icb(knx,1,icc)
      lc=icb(lnx,1,icc)
      iceqjc=ic.eq.jc
      kceqlc=kc.eq.lc
      ngt=(nt(isf)*nt(ksf))**2
      if(iscmci.eq.2) go to 50
      do 14 jesfb=1,2
      do 14 lesfb=1,2
      jbl(jesfb,lesfb)=0
      do 12 ist=1,nst
   12 jbl(jesfb,lesfb)=jbl(jesfb,lesfb)+npriri(jesfb,1,ist)*
     1   nprirk(lesfb,1,ist)
   14 jbld(jesfb)=jbl(jesfb,2)-jbl(jesfb,1)
      jblii=ipq(krcru)*jbl(1,2)+krcru*jbl(1,1)
      jblij=ipq(krcru)*jbl(2,2)+krcru*jbl(2,1)
      jblid=jblij-jblii
      fnfct=(nwt*nfct(icc))
      icxbg(1,1)=ijcxst(1)+ipq((ic-1)*nt(isf)+1)*npri(1,1)+(jc-1)*ntij1
      icxbg(2,1)=ijcxst(2)+((ic-1)*nc(is)+(jc-1))*ntij
      icxbg(3,1)=icxbg(2,1)-(ic-jc)*(nc(is)-1)*ntij
      kcxbg(1,1)=klcxst(1)+ipq((kc-1)*nt(ksf)+1)*nprk(1,1)+(lc-1)*ntkl1
      kcxbg(2,1)=klcxst(2)+((kc-1)*nc(ks)+(lc-1))*ntkl
      kcxbg(3,1)=kcxbg(2,1)-(kc-lc)*(nc(ks)-1)*ntkl
      call jandk(eta,g4,is,isf,ic,is,isf,jc,ks,ksf,kc,ks,ksf,lc,
     1   lmnp1,mcons,ncon,nrcr,nt,ntl,ntu,x,y,z,zet)
      ndx=1
      do 46 ircr=1,ircru
      if(iceqjc) jrcru=ircr
      do 46 jrcr=1,jrcru
      if(ircr-jrcr) 16,18,20
   16 iblkl=(ipq(jrcr)+(ircr-1))*jblij-(jrcr-1)*jblid
      iesfb=2
      icxbg1=icxbg(3,1)
      ngti=(nt(ksf)*nt(ksf))
      ngtj=(nt(isf)*(nt(ksf)*nt(ksf)))
      go to 24
   18 iblkl=ipq(ircr)*jblij+(ircr-1)*jblii
      iesfb=1
      go to 22
   20 iblkl=(ipq(ircr)+(jrcr-1))*jblij-(ircr-1)*jblid
      iesfb=2
   22 icxbg1=icxbg(iesfb,1)
      ngti=(nt(isf)*(nt(ksf)*nt(ksf)))
      ngtj=(nt(ksf)*nt(ksf))
   24 do 46 krcr=1,krcru
      if(kceqlc) lrcru=krcr
      do 46 lrcr=1,lrcru
      if(krcr-lrcr) 26,28,30
   26 iblb=iblkl+(ipq(lrcr)+(krcr-1))*jbl(iesfb,2)-(lrcr-1)*jbld(iesfb)
      kesfb=2
      kcxbg1=kcxbg(3,1)
      ngtk=1
      ngtl=nt(ksf)
      go to 34
   28 iblb=iblkl+ipq(krcr)*jbl(iesfb,2)+(krcr-1)*jbl(iesfb,1)
      kesfb=1
      go to 32
   30 iblb=iblkl+(ipq(krcr)+(lrcr-1))*jbl(iesfb,2)-(krcr-1)*jbld(iesfb)
      kesfb=2
   32 kcxbg1=kcxbg(kesfb,1)
      ngtk=nt(ksf)
      ngtl=1
   34 icxs(iesfb,1)=icxbg1
      ndxi=ndx
      do 44 it=itl,itu
      if(iceqjc.and.ircr.eq.jrcr) jtu=it
      ndxj=ndxi
      do 42 jt=itl,jtu
      kcxs(kesfb,1)=kcxbg1
      ndxk=ndxj
      do 40 kt=ktl,ktu
      if(kceqlc.and.krcr.eq.lrcr) ltu=kt
      ndxl=ndxk
      do 38 lt=ktl,ltu
      val=g4(ndxl)
      if(abs(val).le.cutoff) go to 36
      val=fnfct*val
      ibl=iblb
      call aoso2(a,cx,1,val)
   36 ndxl=ndxl+ngtl
   38 kcxs(kesfb,1)=kcxs(kesfb,1)+nprk(kesfb,1)
   40 ndxk=ndxk+ngtk
      ndxj=ndxj+ngtj
   42 icxs(iesfb,1)=icxs(iesfb,1)+npri(iesfb,1)
   44 ndxi=ndxi+ngti
   46 ndx=ndx+ngt
      return
   50 icxbg(2,2)=ikcxst(2)+((ic-1)*nc(ks)+(kc-1))*ntik
      icxbg(2,3)=icxbg(2,2)+(lc-kc)*ntik
      kcxbg(2,2)=icxbg(2,3)+(jc-ic)*nc(ks)*ntik
      kcxbg(2,3)=kcxbg(2,2)+(kc-lc)*ntik
      jbl1=0
      jbl2=0
      do 52 ist=1,nst
      jbl1=jbl1+ipq(npriri(2,2,ist)+1)
   52 jbl2=jbl2+npriri(2,2,ist)**2
      jbldi=jbl2-jbl1
      if(.not.iceqjc.and..not.kceqlc) go to 54
      fnfct=(nwt*nfct(icc))
      go to 62
   54 if(ngw(icc).eq.1) go to 58
      igwu = ngw(icc)
      do 56 igw=2,igwu
   56 if((icb(inx,igw,icc).eq.ic.and.icb(jnx,igw,icc).eq.jc).or.
     1   (icb(knx,igw,icc).eq.kc.and.icb(lnx,igw,icc).eq.lc)) go to 60
   58 iscmu=3
      fnfct=(nwt*nfct(icc))
      go to 64
   60 fnfct=(2*nwt*nfct(icc))
   62 iscmu=2
   64 do 96 iscm=2,iscmu
      if(iscm.eq.3) go to 66
      call jandk(eta,g4,is,isf,ic,ks,ksf,kc,is,isf,jc,ks,ksf,lc,
     1   lmnp1,mcons,ncon,nrcr,nt,ntl,ntu,x,y,z,zet)
      go to 68
   66 call jandk(eta,g4,is,isf,ic,ks,ksf,lc,is,isf,jc,ks,ksf,kc,
     1   lmnp1,mcons,ncon,nrcr,nt,ntl,ntu,x,y,z,zet)
   68 ndx=1
      ikrcr=0
      do 96 ircr=1,ircru
      if((iceqjc.and.kceqlc)) jrcru=ircr
      do 96 krcr=1,krcru
      ikrcr=ikrcr+1
      jlrcr=0
      do 96 jrcr=1,jrcru
      if(.not.(iceqjc.and.kceqlc)) go to 72
      if(ircr.ne.jrcr) go to 70
      lrcru=krcr
      go to 72
   70 lrcru=krcru
   72 do 96 lrcr=1,lrcru
      jlrcr=jlrcr+1
      if(ikrcr-jlrcr) 74,76,78
   74 iblb=(ipq(jlrcr)+(ikrcr-1))*jbl2-(jlrcr-1)*jbldi
      icxs(2,iscm)=kcxbg(2,iscm)
      kcxbg1=icxbg(2,iscm)
      ngtk=1
      ngtl=(nt(isf)*nt(ksf))
      go to 82
   76 iblb=(ipq(ikrcr+1)-1)*jbl2-(ikrcr-1)*jbldi
      go to 80
   78 iblb=(ipq(ikrcr)+(jlrcr-1))*jbl2-(ikrcr-1)*jbldi
   80 icxs(2,iscm)=icxbg(2,iscm)
      kcxbg1=kcxbg(2,iscm)
      ngtk=(nt(isf)*nt(ksf))
      ngtl=1
   82 ndxk=ndx
      do 94 it=itl,itu
      if(((iceqjc.and.kceqlc).and.ikrcr.eq.jlrcr)) jtu=it
      do 94 kt=ktl,ktu
      kcxs(2,iscm)=kcxbg1
      ndxl=ndxk
      do 92 jt=itl,jtu
      if(.not.((iceqjc.and.kceqlc).and.ikrcr.eq.jlrcr)) go to 86
      if(it.ne.jt) go to 84
      ltu=kt
      go to 86
   84 ltu=ktu
   86 do 92 lt=ktl,ltu
      val=g4(ndxl)
      if(abs(val).le.cutoff) go to 90
      val=fnfct*val
      ibl=iblb
      if(ikrcr.ne.jlrcr) go to 88
      call aoso2e(a,cx,iscm,val)
      go to 90
   88 call aoso2(a,cx,iscm,val)
   90 ndxl=ndxl+ngtl
   92 kcxs(2,iscm)=kcxs(2,iscm)+npri(2,2)
      ndxk=ndxk+ngtk
   94 icxs(2,iscm)=icxs(2,iscm)+npri(2,2)
   96 ndx=ndx+ngt
      return
      end
