      subroutine sf3eq8(a,cx,eta,g4,icb,icc,iscm,is,isf,js,jsf,ks,ksf,
     1   ls,lsf,lmnp1,mcons,nc,ncon,nfct,nrcr,nt,ntl,ntu,x,y,z,zet)
*mdc*on vax ibm cray
      implicit real*8 (a-h,o-z)
*mdc*off
      common /parmr/ cutoff, tol
      common /parmi/ mblu, mccu, mconu, mcu, mcxu,
     1   mrcru, mstu, msu, nblu, ng, ns, nst
      common /ntgr/itl,itu,jtl,jtu,ktl,ktu,ltl,ltu,jcu,kcu,lcu,inx,
     1jnx,knx,lnx,nwt,nwt1(3),ijsf(3),klsf(3),icxs(2,3),kcxs(2,3),
     2npri(2,3),nprk(2,3),iesfb,kesfb,ircru,jrcru,krcru,lrcru
      common /cxindx/ ntij,ntik,ntil,ntkl,ntjl,ntjk,ntij1,ntjk1,ntkl1,
     1   ijcxst(2),ikcxst(2),ilcxst(2),klcxst(2),jlcxst(2),jkcxst(2)
      common /ikpr/ npriri(2,3,8), nprirk(2,3,8)
      common /sf3ao/ ibl
      dimension a(*), cx(*), eta(mrcru,mconu,*), g4(*), icb(4,24,*),
     1   lmnp1(*), mcons(*), nc(*), ncon(*), nfct(*), nrcr(*), nt(*),
     2   ntl(*), ntu(*), x(mcu,*), y(mcu,*), z(mcu,*), zet(mconu,*)
      fnfct=nwt*nfct(icc)
      ic=icb(inx,1,icc)
      jc=icb(jnx,1,icc)
      kc=icb(knx,1,icc)
      lc=icb(lnx,1,icc)
      jbl2=0
      do 10 ist=1,nst
   10 jbl2=jbl2+npriri(2,iscm,ist)*nprirk(2,iscm,ist)
      ndx=0
      iblb=0
      if(iscm-2) 12,40,68
   12 icxbg=ijcxst(2)+((ic-1)*nc(js)+(jc-1))*ntij
      kcxbg=klcxst(2)+((kc-1)*nc(ls)+(lc-1))*ntkl
      call jandk(eta,g4,is,isf,ic,js,jsf,jc,ks,ksf,kc,ls,lsf,lc,
     1   lmnp1,mcons,ncon,nrcr,nt,ntl,ntu,x,y,z,zet)
      do 36 ircr=1,ircru
      do 36 jrcr=1,jrcru
      do 36 krcr=1,krcru
      do 36 lrcr=1,lrcru
      icxs(2,1)=icxbg
      do 24 it=itl,itu
      do 24 jt=jtl,jtu
      kcxs(2,1)=kcxbg
      do 20 kt=ktl,ktu
      do 16 lt=ltl,ltu
      ndx=ndx+1
      val=g4(ndx)
      if(abs(val).le.cutoff) go to 16
      val=fnfct*val
      ibl=iblb
      call aoso2(a,cx,1,val)
   16 kcxs(2,1)=kcxs(2,1)+nprk(2,1)
   20 continue
   24 icxs(2,1)=icxs(2,1)+npri(2,1)
   36 iblb=iblb+jbl2
      return
   40 icxbg=ikcxst(2)+((ic-1)*nc(ks)+(kc-1))*ntik
      kcxbg=jlcxst(2)+((jc-1)*nc(ls)+(lc-1))*ntjl
      call jandk(eta,g4,is,isf,ic,ks,ksf,kc,js,jsf,jc,ls,lsf,lc,
     1   lmnp1,mcons,ncon,nrcr,nt,ntl,ntu,x,y,z,zet)
      do 64 ircr=1,ircru
      do 64 krcr=1,krcru
      do 64 jrcr=1,jrcru
      do 64 lrcr=1,lrcru
      icxs(2,2)=icxbg
      do 52 it=itl,itu
      do 52 kt=ktl,ktu
      kcxs(2,2)=kcxbg
      do 48 jt=jtl,jtu
      do 44 lt=ltl,ltu
      ndx=ndx+1
      val=g4(ndx)
      if(abs(val).le.cutoff) go to 44
      val=fnfct*val
      ibl=iblb
      call aoso2(a,cx,2,val)
   44 kcxs(2,2)=kcxs(2,2)+nprk(2,2)
   48 continue
   52 icxs(2,2)=icxs(2,2)+npri(2,2)
   64 iblb=iblb+jbl2
      return
   68 icxbg=ilcxst(2)+((ic-1)*nc(ls)+(lc-1))*ntil
      kcxbg=jkcxst(2)+((jc-1)*nc(ks)+(kc-1))*ntjk
      call jandk(eta,g4,is,isf,ic,ls,lsf,lc,js,jsf,jc,ks,ksf,kc,
     1   lmnp1,mcons,ncon,nrcr,nt,ntl,ntu,x,y,z,zet)
      do 92 ircr=1,ircru
      do 92 lrcr=1,lrcru
      do 92 jrcr=1,jrcru
      do 92 krcr=1,krcru
      icxs(2,3)=icxbg
      do 80 it=itl,itu
      do 80 lt=ltl,ltu
      kcxs(2,3)=kcxbg
      do 76 jt=jtl,jtu
      do 72 kt=ktl,ktu
      ndx=ndx+1
      val=g4(ndx)
      if(abs(val).le.cutoff) go to 72
      val=fnfct*val
      ibl=iblb
      call aoso2(a,cx,3,val)
   72 kcxs(2,3)=kcxs(2,3)+nprk(2,3)
   76 continue
   80 icxs(2,3)=icxs(2,3)+npri(2,3)
   92 iblb=iblb+jbl2
      return
      end
