      subroutine sf3eq2(a,cx,eta,g4,icb,icc,is,isf,ls,lsf,lmnp1,
     1   mcons,nc,ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
*mdc*on vax ibm cray
      implicit real*8 (a-h,o-z)
*mdc*off
      logical ec
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
      dimension ec(3), icxbg(3,3), jbl(2), kcxbg(3)
      ic=icb(inx,1,icc)
      jc=icb(jnx,1,icc)
      kc=icb(knx,1,icc)
      lc=icb(lnx,1,icc)
      icxbg(1,1)=ijcxst(1)+ipq((ic-1)*nt(isf)+1)*npri(1,1)+(jc-1)*ntij1
      icxbg(1,2)=icxbg(1,1)+(kc-jc)*ntij1
      kcxbg(3)=ilcxst(2)+((ic-1)*nc(ls)+(lc-1))*ntil
      kcxbg(1)=kcxbg(3)+(kc-ic)*(nc(ls)*ntil)
      kcxbg(2)=kcxbg(1)+(jc-kc)*(nc(ls)*ntil)
      icxbg(1,3)=ijcxst(1)+ipq((jc-1)*nt(isf)+1)*npri(1,1)+(kc-1)*ntij1
      icxbg(2,1)=ijcxst(2)+((ic-1)*nc(is)+(jc-1))*ntij
      icxbg(3,1)=icxbg(2,1)-(ic-jc)*((nc(is)-1)*ntij)
      icxbg(2,2)=icxbg(2,1)+(kc-jc)*ntij
      icxbg(3,2)=icxbg(2,2)-(ic-kc)*((nc(is)-1)*ntij)
      icxbg(2,3)=icxbg(2,2)+(jc-ic)*nc(is)*ntij
      icxbg(3,3)=icxbg(2,3)-(jc-kc)*((nc(is)-1)*ntij)
      ec(1)=ic.eq.jc
      ec(3)=jc.eq.kc
      ec(2)=ic.eq.kc
      do 10 jesfb=1,2
      jbl(jesfb)=0
      do 10 ist=1,nst
   10 jbl(jesfb)=jbl(jesfb)+npriri(jesfb,1,ist)*nprirk(2,1,ist)
      jbl1=(ircru*lrcru)*jbl(1)
      jbl2=(ircru*lrcru)*jbl(2)
      jbld=jbl2-jbl1
      if(.not.ec(1)) go to 16
      nwt1(1)=1
      if(.not.ec(3)) go to 12
      iscmu=1
      go to 48
   12 iscmu=2
      go to 44
   16 if(.not.ec(3)) go to 20
      iscml=2
      iscmu=3
      nwt1(2)=1
      go to 52
   20 if(ngw(icc).eq.1) go to 40
      isw2 = 2
      isw3 = 3
      igwu = ngw(icc)
      do 24 igw=2,igwu
      if(icb(lnx,igw,icc).ne.lc) go to 24
      if(icb(jnx,igw,icc).eq.kc) isw2=1
      if(icb(inx,igw,icc).eq.jc) isw3=min(isw3,2)
      if(icb(inx,igw,icc).eq.kc) isw3=1
   24 continue
      if(isw3.ne.1) go to 32
      if(isw2.ne.1) go to 28
      iscmu = 1
      nwt1 (1) = 3
      go to 48
   28 iscmu=2
      nwt1(1) = 2
      go to 44
   32 if(isw3.ne.2) go to 36
      iscmu = 2
      nwt1(1)=1
      nwt1(2) = 2
      go to 48
   36 if(isw2.ne.1) go to 40
      iscml = 2
      iscmu = 3
      nwt1(2) = 2
      go to 52
   40 iscmu=3
      nwt1(1)=1
   44 nwt1(2)=1
   48 iscml = 1
   52 nwt2 = nwt * nfct(icc)
      ngt=nt(isf)*(nt(isf)*(nt(isf)*nt(lsf)))
      do 96 iscm=iscml,iscmu
      fnfct=nwt1(iscm)*nwt2
      if(iscm-2) 56,60,64
   56 call jandk(eta,g4,is,isf,ic,is,isf,jc,is,isf,kc,ls,lsf,lc,
     1   lmnp1,mcons,ncon,nrcr,nt,ntl,ntu,x,y,z,zet)
      go to 68
   60 call jandk(eta,g4,is,isf,ic,is,isf,kc,is,isf,jc,ls,lsf,lc,
     1   lmnp1,mcons,ncon,nrcr,nt,ntl,ntu,x,y,z,zet)
      go to 68
   64 call jandk(eta,g4,is,isf,jc,is,isf,kc,is,isf,ic,ls,lsf,lc,
     1   lmnp1,mcons,ncon,nrcr,nt,ntl,ntu,x,y,z,zet)
   68 ndx=1
      do 96 ircr=1,ircru
      if(ec(iscm)) jrcru=ircr
      do 96 jrcr=1,jrcru
      if(ircr-jrcr) 74,76,78
   74 iblb=(ipq(jrcr)+(ircr-1))*jbl2-(jrcr-1)*jbld
      iesfb=2
      icxbg1=icxbg(3,iscm)
      ngti=(nt(isf)*nt(lsf))
      ngtj=(nt(isf)*(nt(isf)*nt(lsf)))
      go to 82
   76 iblb=ipq(ircr)*jbl2+(ircr-1)*jbl1
      iesfb=1
      go to 80
   78 iblb=(ipq(ircr)+(jrcr-1))*jbl2-(ircr-1)*jbld
      iesfb=2
   80 icxbg1=icxbg(iesfb,iscm)
      ngti=(nt(isf)*(nt(isf)*nt(lsf)))
      ngtj=(nt(isf)*nt(lsf))
   82 do 96 krcr=1,ircru
      do 96 lrcr=1,lrcru
      icxs(iesfb,iscm)=icxbg1
      ndxi=ndx
      do 92 it=itl,itu
      if(ircr.eq.jrcr.and.ec(iscm)) jtu=it
      ndxj=ndxi
      do 90 jt=itl,jtu
      kcxs(2,iscm)=kcxbg(iscm)
      ndxl=ndxj
      do 88 kt=itl,itu
      do 88 lt=ltl,ltu
      val=g4(ndxl)
      if(abs(val).le.cutoff) go to 86
      val=fnfct*val
      ibl=iblb
      call aoso2(a,cx,iscm,val)
   86 ndxl=ndxl+1
   88 kcxs(2,iscm)=kcxs(2,iscm)+nprk(2,1)
      ndxj=ndxj+ngtj
   90 icxs(iesfb,iscm)=icxs(iesfb,iscm)+npri(iesfb,1)
   92 ndxi=ndxi+ngti
      iblb=iblb+jbl(iesfb)
   96 ndx=ndx+ngt
      return
      end
