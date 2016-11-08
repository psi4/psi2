      subroutine sf3eq4(a,cx,eta,g4,icb,icc,iscmci,is,isf,ks,ksf,ls,
     1   lsf,lmnp1,mcons,nc,ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
*mdc*on vax ibm cray
      implicit real*8 (a-h,o-z)
*mdc*off
      logical iceqjc
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
      dimension icxbg(3,3), jbl(2), kcxbg(3)
      ic=icb(inx,1,icc)
      jc=icb(jnx,1,icc)
      kc=icb(knx,1,icc)
      lc=icb(lnx,1,icc)
      iceqjc=ic.eq.jc
      do 10 jesfb=iscmci,2
      jbl(jesfb)=0
      do 10 ist=1,nst
   10 jbl(jesfb)=jbl(jesfb)+npriri(jesfb,iscmci,ist)*
     1   nprirk(2,iscmci,ist)
      if(iscmci.eq.2) go to 44
      jbl1=(krcru*lrcru)*jbl(1)
      jbl2=(krcru*lrcru)*jbl(2)
      jbld=jbl2-jbl1
      fnfct=(nwt*nfct(icc))
      icxbg(1,1)=ijcxst(1)+ipq((ic-1)*nt(isf)+1)*npri(1,1)+(jc-1)*ntij1
      icxbg(2,1)=ijcxst(2)+((ic-1)*nc(is)+(jc-1))*ntij
      icxbg(3,1)=icxbg(2,1)-(ic-jc)*(nc(is)-1)*ntij
      kcxbg(1)=klcxst(2)+((kc-1)*nc(ls)+(lc-1))*ntkl
      call jandk(eta,g4,is,isf,ic,is,isf,jc,ks,ksf,kc,ls,lsf,lc,
     1   lmnp1,mcons,ncon,nrcr,nt,ntl,ntu,x,y,z,zet)
      ndx=1
      ngt=nt(isf)*(nt(isf)*(nt(ksf)*nt(lsf)))
      do 40 ircr=1,ircru
      if(iceqjc) jrcru=ircr
      do 40 jrcr=1,jrcru
      if(ircr-jrcr) 12,13,14
   12 iblb=(ipq(jrcr)+(ircr-1))*jbl2-(jrcr-1)*jbld
      iesfb=2
      icxbg1=icxbg(3,1)
      ngti=(nt(ksf)*nt(lsf))
      ngtj=(nt(isf)*(nt(ksf)*nt(lsf)))
      go to 16
   13 iblb=ipq(ircr)*jbl2+(ircr-1)*jbl1
      iesfb=1
      go to 15
   14 iblb=(ipq(ircr)+(jrcr-1))*jbl2-(ircr-1)*jbld
      iesfb=2
   15 icxbg1=icxbg(iesfb,1)
      ngti=(nt(isf)*(nt(ksf)*nt(lsf)))
      ngtj=(nt(ksf)*nt(lsf))
   16 do 40 krcr=1,krcru
      do 40 lrcr=1,lrcru
      icxs(iesfb,1)=icxbg1
      ndxi=ndx
      do 28 it=itl,itu
      if(iceqjc.and.ircr.eq.jrcr) jtu=it
      ndxj=ndxi
      do 24 jt=itl,jtu
      kcxs(2,1)=kcxbg(1)
      ndxl=ndxj
      do 20 kt=ktl,ktu
      do 20 lt=ltl,ltu
      val=g4(ndxl)
      if(abs(val).le.cutoff) go to 18
      val=fnfct*val
      ibl=iblb
      call aoso2(a,cx,1,val)
   18 ndxl=ndxl+1
   20 kcxs(2,1)=kcxs(2,1)+nprk(2,1)
      ndxj=ndxj+ngtj
   24 icxs(iesfb,1)=icxs(iesfb,1)+npri(iesfb,1)
   28 ndxi=ndxi+ngti
      iblb=iblb+jbl(iesfb)
   40 ndx=ndx+ngt
      return
   44 icxbg(2,2)=ikcxst(2)+((ic-1)*nc(ks)+(kc-1))*ntik
      icxbg(2,3)=icxbg(2,2)+(jc-ic)*nc(ks)*ntik
      kcxbg(3)=ilcxst(2)+((ic-1)*nc(ls)+(lc-1))*ntil
      kcxbg(2)=kcxbg(3)+(jc-ic)*nc(ls)*ntil
      if(.not.iceqjc) go to 48
      fnfct=(nwt*nfct(icc))
      go to 64
   48 if(ngw(icc).eq.1) go to 56
      igwu = ngw(icc)
      do 52 igw=2,igwu
   52 if(icb(knx,igw,icc).eq.kc.and.icb(lnx,igw,icc).eq.lc.and.
     1   icb(inx,igw,icc).eq.jc) go to 60
   56 iscmu=3
      fnfct=(nwt*nfct(icc))
      go to 68
   60 fnfct=(2*nwt*nfct(icc))
   64 iscmu=2
   68 do 96 iscm=2,iscmu
      if(iscm.eq.3) go to 72
      call jandk(eta,g4,is,isf,ic,ks,ksf,kc,is,isf,jc,ls,lsf,lc,
     1   lmnp1,mcons,ncon,nrcr,nt,ntl,ntu,x,y,z,zet)
      go to 76
   72 call jandk(eta,g4,is,isf,jc,ks,ksf,kc,is,isf,ic,ls,lsf,lc,
     1   lmnp1,mcons,ncon,nrcr,nt,ntl,ntu,x,y,z,zet)
   76 ndx=0
      iblb=0
      do 96 ircr=1,ircru
      do 96 krcr=1,krcru
      do 96 jrcr=1,ircru
      do 96 lrcr=1,lrcru
      icxs(2,iscm)=icxbg(2,iscm)
      do 84 it=itl,itu
      do 84 kt=ktl,ktu
      kcxs(2,iscm)=kcxbg(iscm)
      do 80 jt=itl,itu
      do 80 lt=ltl,ltu
      ndx=ndx+1
      val=g4(ndx)
      if(abs(val).le.cutoff) go to 80
      val=fnfct*val
      ibl=iblb
      call aoso2(a,cx,iscm,val)
   80 kcxs(2,iscm)=kcxs(2,iscm)+nprk(2,2)
   84 icxs(2,iscm)=icxs(2,iscm)+npri(2,2)
      iblb=iblb+jbl(2)
   96 continue
      return
      end
