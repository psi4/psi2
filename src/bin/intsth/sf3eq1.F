      subroutine sf3eq1(a,cx,eta,g4,icb,icc,is,isf,lmnp1,mcons,nc,
     1   ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
*mdc*on vax ibm cray
      implicit real*8 (a-h,o-z)
*mdc*off
      logical eij, eij1, ekl, esfc, esfcij, esfckl
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
      dimension esfc(3),esfcij(3),esfckl(3),icxbg(3,3),jbl(6),kcxbg(3,3)
      ic=icb(1,1,icc)
      jc=icb(2,1,icc)
      kc=icb(3,1,icc)
      lc=icb(4,1,icc)
      icxbg(1,1)=ijcxst(1)+ipq((ic-1)*nt(isf)+1)*npri(1,1)+(jc-1)*ntij1
      icxbg(1,2)=icxbg(1,1)+(kc-jc)*ntij1
      icxbg(1,3)=icxbg(1,2)+(lc-kc)*ntij1
      kcxbg(1,1)=ijcxst(1)+ipq((kc-1)*nt(isf)+1)*npri(1,1)+(lc-1)*ntij1
      kcxbg(1,2)=ijcxst(1)+ipq((jc-1)*nt(isf)+1)*npri(1,1)+(lc-1)*ntij1
      kcxbg(1,3)=kcxbg(1,2)+(kc-lc)*ntij1
      icxbg(2,1)=ijcxst(2)+((ic-1)*nc(is)+(jc-1))*ntij
      icxbg(3,1)=icxbg(2,1)-(ic-jc)*((nc(is)-1)*ntij)
      icxbg(2,2)=icxbg(2,1)+(kc-jc)*ntij
      icxbg(3,2)=icxbg(2,2)-(ic-kc)*((nc(is)-1)*ntij)
      icxbg(2,3)=icxbg(2,2)+(lc-kc)*ntij
      icxbg(3,3)=icxbg(2,3)-(ic-lc)*((nc(is)-1)*ntij)
      kcxbg(2,1)=icxbg(2,3)+(kc-ic)*(nc(is)*ntij)
      kcxbg(3,1)=kcxbg(2,1)-(kc-lc)*((nc(is)-1)*ntij)
      kcxbg(2,2)=kcxbg(2,1)+(jc-kc)*(nc(is)*ntij)
      kcxbg(3,2)=kcxbg(2,2)-(jc-lc)*((nc(is)-1)*ntij)
      kcxbg(2,3)=kcxbg(2,2)+(kc-lc)*ntij
      kcxbg(3,3)=kcxbg(2,3)-(jc-kc)*((nc(is)-1)*ntij)
      esfcij(1)=ic.eq.jc
      esfckl(3)=jc.eq.kc
      esfckl(1)=kc.eq.lc
      esfcij(2)=ic.eq.kc
      esfckl(2)=jc.eq.lc
      esfcij(3)=ic.eq.lc
      esfc(1)=esfcij(2).and.esfckl(2)
      esfc(2)=esfcij(1).and.esfckl(1)
      esfc(3)=esfc(2)
      CALL IZERO(JBL,6)
      do 10 ist=1,nst
      jbl(1)=jbl(1)+ipq(npriri(1,1,ist)+1)
      jbl(3)=jbl(3)+ipq(npriri(2,1,ist)+1)
      jbl(4)=jbl(4)+npriri(1,1,ist)**2
      jbl(5)=jbl(5)+npriri(2,1,ist)*npriri(1,1,ist)
   10 jbl(6)=jbl(6)+npriri(2,1,ist)**2
      if(.not.esfcij(1).and..not.esfckl(1)) go to 14
      nwt1(1)=1
      if(.not.esfckl(3)) go to 12
      iscmu=1
      go to 30
   12 iscmu=2
      go to 28
   14 if(.not.esfckl(3)) go to 16
      iscml=2
      iscmu=3
      nwt1(2)=1
      go to 32
   16 if(ngw(icc).eq.1) go to 26
      isw2 = 2
      isw3 = 3
      igwu = ngw(icc)
      do 18 igw=2,igwu
      icb1 = icb(inx,igw,icc)
      icb2 = icb(knx,igw,icc)
      itst1 = max(icb1, icb2)
      itst2 = min(icb1, icb2)
      if((itst1.eq.ic.and.itst2.eq.jc).or.(itst1.eq.kc.and.itst2.eq.lc))
     1   isw2 = 1
      icb2 = icb(lnx,igw,icc)
      itst1 = max(icb1, icb2)
      itst2 = min(icb1, icb2)
      if((itst1.eq.ic.and.itst2.eq.kc).or.(itst1.eq.jc.and.itst2.eq.lc))
     1   isw3 = min(isw3,2)
      if((itst1.eq.ic.and.itst2.eq.jc).or.(itst1.eq.kc.and.itst2.eq.lc))
     1   isw3 = 1
   18 continue
      if(isw3.ne.1) go to 22
      if(isw2.ne.1) go to 20
      iscmu = 1
      nwt1 (1) = 3
      go to 30
   20 iscmu=2
      nwt1(1) = 2
      go to 28
   22 if(isw3.ne.2) go to 24
      iscmu = 2
      nwt1(1)=1
      nwt1(2) = 2
      go to 30
   24 if(isw2.ne.1) go to 26
      iscml = 2
      iscmu = 3
      nwt1(2) = 2
      go to 32
   26 iscmu=3
      nwt1(1)=1
   28 nwt1(2)=1
   30 iscml = 1
   32 ngt=nt(isf)*(nt(isf)*(nt(isf)*nt(isf)))
      do 99 iscm=iscml,iscmu
      fnfct=nwt1(iscm)*nfct(icc)
      if(iscm-2) 34,35,36
   34 call jandk(eta,g4,is,isf,ic,is,isf,jc,is,isf,kc,is,isf,lc,
     1   lmnp1,mcons,ncon,nrcr,nt,ntl,ntu,x,y,z,zet)
      go to 38
   35 call jandk(eta,g4,is,isf,ic,is,isf,kc,is,isf,jc,is,isf,lc,
     1   lmnp1,mcons,ncon,nrcr,nt,ntl,ntu,x,y,z,zet)
      go to 38
   36 call jandk(eta,g4,is,isf,ic,is,isf,lc,is,isf,jc,is,isf,kc,
     1   lmnp1,mcons,ncon,nrcr,nt,ntl,ntu,x,y,z,zet)
   38 ndx=1
      do 99 ircr=1,ircru
      if(esfc(iscm)) krcru=ircr
      if(esfcij(iscm)) jrcru=ircr
      do 99 jrcr=1,jrcru
      if(ircr-jrcr) 40,41,42
   40 ijrcr=ipq(jrcr)+ircr
      ije=jrcr-1
      iesfb1=2
      icxbg1=icxbg(3,iscm)
      ngti1=(nt(isf)*nt(isf))
      ngtj1=(nt(isf)*(nt(isf)*nt(isf)))
      go to 46
   41 ijrcr=ipq(ircr+1)
      ije=ircr
      iesfb1=1
      go to 44
   42 ijrcr=ipq(ircr)+jrcr
      ije=ircr-1
      iesfb1=2
   44 icxbg1=icxbg(iesfb1,iscm)
      ngti1=(nt(isf)*(nt(isf)*nt(isf)))
      ngtj1=(nt(isf)*nt(isf))
   46 ijn=ijrcr-ije
      eij1=esfcij(iscm).and.ircr.eq.jrcr
      do 99 krcr=1,krcru
      if(.not.esfc(iscm)) go to 52
      if(ircr.ne.krcr) go to 48
      lrcru=jrcr
      go to 54
   48 if(.not.esfckl(iscm)) go to 50
      lrcru=krcr
      go to 54
   50 lrcru=ircru
      go to 54
   52 if(esfckl(iscm)) lrcru=krcr
   54 do 99 lrcr=1,lrcru
      if(krcr-lrcr) 56,57,58
   56 klrcr=ipq(lrcr)+krcr
      kle=lrcr-1
      kesfb=2
      kcxbg1=kcxbg(3,iscm)
      ngtk1=1
      ngtl1=nt(isf)
      go to 62
   57 klrcr=ipq(krcr+1)
      kle=krcr
      kesfb=1
      go to 60
   58 klrcr=ipq(krcr)+lrcr
      kle=krcr-1
      kesfb=2
   60 kcxbg1=kcxbg(kesfb,iscm)
      ngtk1=nt(isf)
      ngtl1=1
   62 kln=klrcr-kle
      ekl=esfckl(iscm).and.krcr.eq.lrcr
      if(ijrcr-klrcr) 64,70,72
   64 if(kesfb.eq.2) go to 66
      iblb=(kle-1)*jbl(1)+kln*jbl(3)+(ipq(kle-1)+ije)*jbl(4)+(kln*
     1   (kle-1)+ijn)*jbl(5)+ipq(kln)*jbl(6)-jbl(iesfb1+3)
      go to 68
   66 iblb=kle*jbl(1)+(kln-1)*jbl(3)+ipq(kle)*jbl(4)+((kln-1)*kle+
     1   ije)*jbl(5)+(ipq(kln)-kln+1+ijn)*jbl(6)-jbl(iesfb1+4)
   68 iesfb=kesfb
      kesfb=iesfb1
      icxs(iesfb,iscm)=kcxbg1
      kcxbg1=icxbg1
      eij=ekl
      ekl=eij1
      ngti=ngtk1
      ngtj=ngtl1
      ngtk=ngti1
      ngtl=ngtj1
      go to 78
   70 iblb=ije*jbl(1)+ijn*jbl(3)+ipq(ije)*jbl(4)+ijn*ije*jbl(5)+
     1   (ipq(ijn+1)-ijn)*jbl(6)-jbl(iesfb1+kesfb-1)
      go to 76
   72 if(iesfb1.eq.2) go to 74
      iblb=(ije-1)*jbl(1)+ijn*jbl(3)+(ipq(ije-1)+kle)*jbl(4)+(ijn*
     1   (ije-1)+kln)*jbl(5)+ipq(ijn)*jbl(6)-jbl(kesfb+3)
      go to 76
   74 iblb=ije*jbl(1)+(ijn-1)*jbl(3)+ipq(ije)*jbl(4)+((ijn-1)*ije+
     1   kle)*jbl(5)+(ipq(ijn)-ijn+1+kln)*jbl(6)-jbl(kesfb+4)
   76 iesfb=iesfb1
      icxs(iesfb,iscm)=icxbg1
      eij=eij1
      ngti=ngti1
      ngtj=ngtj1
      ngtk=ngtk1
      ngtl=ngtl1
   78 ndxi=ndx
      do 98 it=itl,itu
      if(((esfc(iscm).and.ircr.eq.krcr).and.jrcr.eq.lrcr)) ktu=it
      ndxj=ndxi
      if(eij) jtu=it
      do 96 jt=itl,jtu
      ndxk=ndxj
      kcxs(kesfb,iscm)=kcxbg1
      do 94 kt=itl,ktu
      ndxl=ndxk
      if(.not.((esfc(iscm).and.ircr.eq.krcr).and.jrcr.eq.lrcr)) go to 84
      if(it.ne.kt) go to 80
      ltu=jt
      go to 86
   80 if(.not.ekl) go to 82
      ltu=kt
      go to 86
   82 ltu=itu
      go to 86
   84 if(ekl) ltu=kt
   86 do 92 lt=itl,ltu
      val=g4(ndxl)
      if(abs(val).le.cutoff) go to 90
      val=fnfct*val
      ibl=iblb
      if(ijrcr.ne.klrcr) go to 88
      call aoso2e(a,cx,iscm,val)
      go to 90
   88 call aoso2(a,cx,iscm,val)
   90 ndxl=ndxl+ngtl
   92 kcxs(kesfb,iscm)=kcxs(kesfb,iscm)+npri(kesfb,1)
   94 ndxk=ndxk+ngtk
      ndxj=ndxj+ngtj
   96 icxs(iesfb,iscm)=icxs(iesfb,iscm)+npri(iesfb,1)
   98 ndxi=ndxi+ngti
   99 ndx=ndx+ngt
      return
      end
