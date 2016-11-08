      subroutine twoint(a,cx,eta,ica,icb,icxst,il,iprst,lmnp1,mcons,
     1   nc,ncon,nf,nfct,ngw,npair,nprir,nrcr,nt,ntl,ntu,x,y,z,zet)
*mdc*on vax cray ibm
      implicit real*8 (a-h,o-z)
*mdc*off
      logical eijs, ejks, ekls, eprsf, esfi, esfk
      COMMON /OUT/    IOUT,ICHECK,IFLAG3,NSECT,NINTS
      COMMON /TAPE/ NTAPE
      common /parmi/ mblu, mccu, mconu, mcu, mcxu,
     1   mrcru, mstu, msu, nblu, ng, ns, nst
      common /lgcl/ eprsf(3), esfi(3), esfk(3)
      common /ntgr/itl,itu,jtl,jtu,ktl,ktu,ltl,ltu,jcu,kcu,lcu,inx,
     1jnx,knx,lnx,nwt,nwt1(3),ijsf(3),klsf(3),icxs(2,3),kcxs(2,3),
     2npri(2,3),nprk(2,3),iesfb,kesfb,ircru,jrcru,krcru,lrcru
      common /cxindx/ ntij,ntik,ntil,ntkl,ntjl,ntjk,ntij1,ntjk1,ntkl1,
     1   ijcxst(2),ikcxst(2),ilcxst(2),klcxst(2),jlcxst(2),jkcxst(2)
      COMMON /OUTBUF/ IBLK, IBUF, ALBL(2235), AI(2235)
      common /dim21/ ipq(256)
      common /ikpr/ npriri(2,3,8), nprirk(2,3,8)
      common / maxval / mxblu
c
      dimension a(*),cx(*),eta(mrcru,mconu,*),ica(mcu,msu,*),
     1  icb(4,24,*),icxst(2,*),il(*),iprst(*),lmnp1(*),mcons(*),nc(*),
     2  ncon(*),nf(*),nfct(*),ngw(*),npair(2,*),nprir(2,mstu,*),nrcr(*),
     3  nt(*),ntl(*),ntu(*),x(mcu,*),y(mcu,*),z(mcu,*),zet(mconu,*)
      DIMENSION LBL(4470)
      EQUIVALENCE (ALBL,LBL)
      dimension ijb(2,3), klb(2,3), niqd(3), nqd(3)
      dimension iperm(4,24),ip21(4),ip31(6),ip41(12),ip42(3),ip43(2),
     1   ip51(4),ip61(12),ip62(3),ip65(3),ip71(12),ip73(2),ip75(3),
     2   ip82(6),ip83(4),ip84(2),ip85(6),ip86(2),ip87(2)
c
      data iperm/1,2,3,4, 2,1,3,4, 1,3,2,4, 3,1,2,4, 2,3,1,4, 3,2,1,4,
     +           1,2,4,3, 2,1,4,3, 1,4,2,3, 4,1,2,3, 2,4,1,3, 4,2,1,3,
     +           1,3,4,2, 3,1,4,2, 1,4,3,2, 4,1,3,2, 3,4,1,2, 4,3,1,2,
     +           2,3,4,1, 3,2,4,1, 2,4,3,1, 4,2,3,1, 3,4,2,1, 4,3,2,1/
      data ip21/1,7,13,19/
      data ip31/1,3,9,5,11,17/
      data ip41/1,3,5,7,9,11,13,15,17,19,21,23/
      data ip42/1,3,5/
      data ip43/1,7/
      data ip51/10,4,2,1/
      data ip61/4,2,1,10,8,7,16,14,13,22,20,19/
      data ip62/4,2,1/
      data ip65/1,7,13/
      data ip71/17,11,9,18,5,3,12,6,1,10,4,2/
      data ip73/1,2/
      data ip75/9,3,1/
      data ip82/1,2,3,4,5,6/
      data ip83/1,2,7,8/
      data ip84/1,2/
      data ip85/1,3,7,9,13,15/
      data ip86/1,3/
      data ip87/1,7/
C
      DATA IDTR/2235/
c
      mxblu = -1
c
      nwt1(3)=1
      iblk=0
      IBUF = 0
      NSCTLD=0
      IBFLD = 0
      npkflg=0
      isf = 0
c     ----- ijsf(1) is ijsf -----
      ijsf(1) = 0
c     ----- ijsf(2) is iksf -----
      ijsf(2) = 0
c     ----- ijsf(3) is ilsf -----
      ijsf(3) = 0
      do 9984 is = 1,ns
      ifu = nf(is)
      if(ifu.eq.0) go to 9984
      isfis= isf
      iksfis = ijsf(2)
      ilsfis = ijsf(3)
      ic = nc(is)
      icu = nc(is)
      ieqs1 = 8
      jsf = 0
c     ----- klsf(2) is jlsf -----
      klsf(2) = 0
c     ----- klsf(3) is jksf -----
      klsf(3) = 0
      do 9983 js = 1,is
      eijs=is.eq.js
      if(.not.eijs) go to 10
      ieqs1=4
      go to 12
   10 jfu = nf(js)
      if(jfu.eq.0) go to 9983
   12 ieqs2 = ieqs1
      jcu = nc(js)
      jsfjs = jsf
      ijsfjs = ijsf(1)
      jlsfjs = klsf(2)
      ksf = 0
      ijsf(2) = iksfis
c     ----- klsf(1) is klsf -----
      klsf(1) = 0
      do 9982 ks = 1,js
      ejks=js.eq.ks
      if(.not.ejks) go to 14
      ieqs2=ieqs1-2
      go to 16
   14 kfu = nf(ks)
      if(kfu.eq.0) go to 9982
   16 ieqs3 = ieqs2
      kcu = nc(ks)
      ksfks= ksf
      iksfks = ijsf(2)
      jksfks = klsf(3)
      lsf = 0
      ijsf(3) = ilsfis
      klsf(2) = jlsfjs
      do 9981 ls = 1,ks
      ekls=ks.eq.ls
      if(.not.ekls) go to 18
      ieqs3=ieqs2-1
      go to 20
   18 lfu = nf(ls)
      if(lfu.eq.0) go to 9981
   20 ncc = 0
      lcu = nc(ls)
      lsfls = lsf
      klsfls = klsf(1)
      jlsfls = klsf(2)
      call cecob(ic,ica,icb,ieqs3,is,js,ks,ls,ncc,nfct,ngw)
      isf = isfis
      ijsf(1) = ijsfjs
      ijsf(2) = iksfks
      do 998 if = 1, ifu
      isf = isf + 1
      ircru=nrcr(mcons(isf))
      itl = ntl(isf)
      itu = ntu(isf)
      iksfif = ijsf(2)
      ilsfif = ijsf(3)
      ieqsf1 = 8
      jsf = jsfjs
      klsf(2) = jlsfls
      klsf(3) = jksfks
      if(eijs) jfu=if
      do 998 jf = 1, jfu
      jlsfjf = klsf(2)
      jsf = jsf + 1
      jrcru=nrcr(mcons(jsf))
      jtl = ntl(jsf)
      jtu = ntu(jsf)
      ijsf(1) = ijsf(1) + 1
      ntij=nt(isf)*nt(jsf)*npair(2,ijsf(1))
      esfi(1)=isf.eq.jsf
      if(.not.esfi(1)) go to 22
      ieqsf1=4
      ntij1=nt(isf)**2*npair(1,ijsf(1))
      ijb(1,1)=ircru
   22 ijb(2,1)=ircru*jrcru
      ieqsf2 = ieqsf1
      ksf = ksfks
      ijsf(2) = iksfif
      klsf(1) = klsfls
      if(ejks) kfu=jf
      do 998 kf = 1, kfu
      ksf = ksf + 1
      krcru=nrcr(mcons(ksf))
      ktl = ntl(ksf)
      ktu = ntu(ksf)
      ijsf(2) = ijsf(2) + 1
      ntik=nt(isf)*nt(ksf)*npair(2,ijsf(2))
      klsf(3) = klsf(3) + 1
      ntjk=nt(jsf)*nt(ksf)*npair(2,klsf(3))
      esfi(2)=isf.eq.ksf
      if(esfi(2)) ijb(1,2)=ircru
      ijb(2,2)=ircru*krcru
      esfk(3)=jsf.eq.ksf
      if(.not.esfk(3)) go to 24
      ieqsf2=ieqsf1-2
      ntjk1=nt(jsf)**2*npair(1,klsf(3))
      klb(1,3)=jrcru
   24 klb(2,3)=jrcru*krcru
      ieqsf3 = ieqsf2
      lsf = lsfls
      ijsf(3) = ilsfif
      klsf(2) = jlsfjf
      if(ekls) lfu=kf
      do 997 lf = 1, lfu
      lsf = lsf + 1
      lrcru=nrcr(mcons(lsf))
      ltl = ntl(lsf)
      ltu = ntu(lsf)
      ijsf(3) = ijsf(3) + 1
      ntil=nt(isf)*nt(lsf)*npair(2,ijsf(3))
      klsf(1) = klsf(1) + 1
      ntkl=nt(ksf)*nt(lsf)*npair(2,klsf(1))
      klsf(2) = klsf(2) + 1
      ntjl=nt(jsf)*nt(lsf)*npair(2,klsf(2))
      esfi(3)=isf.eq.lsf
      if(esfi(3)) ijb(1,3)=ircru
      ijb(2,3)=ircru*lrcru
      esfk(2)=jsf.eq.lsf
      if(esfk(2)) klb(1,2)=jrcru
      klb(2,2)=jrcru*lrcru
      esfk(1)=ksf.eq.lsf
      if(.not.esfk(1)) go to 26
      ieqsf3=ieqsf2-1
      ntkl1=nt(ksf)**2*npair(1,klsf(1))
      klb(1,1)=krcru
   26 klb(2,1)=krcru*lrcru
c     ----- set up arrays by permutation -----
      ibl = 0
c     ----- do 460 iscm=1,3 -----
      iscm = 1
      go to 280
  210 if(max(ijsf(iscm),klsf(iscm)).ne.ijsf(iscm-1)) go to 280
c     ----- same block of integrals -----
      nqd(iscm) = nqd(iscm-1)
      if(.not.esfi(iscm-1)) go to 230
      npri(1,iscm)=npri(1,iscm-1)
      do 220 ist=1,nst
  220 npriri(1,iscm,ist)=npriri(1,iscm-1,ist)
  230 npri(2,iscm)=npri(2,iscm-1)
      do 240 ist=1,nst
  240 npriri(2,iscm,ist)=npriri(2,iscm-1,ist)
      if(.not.esfk(iscm-1)) go to 260
      nprk(1,iscm)=nprk(1,iscm-1)
      do 250 ist=1,nst
  250 nprirk(1,iscm,ist)=nprirk(1,iscm-1,ist)
  260 nprk(2,iscm)=nprk(2,iscm-1)
      do 270 ist=1,nst
  270 nprirk(2,iscm,ist)=nprirk(2,iscm-1,ist)
      niqd(iscm) = 0
      go to 450
c     ----- new block of integrals -----
  280 if(.not.esfi(iscm)) go to 300
      npri(1,iscm)=npair(1,ijsf(iscm))
      do 290 ist=1,nst
  290 npriri(1,iscm,ist)=nprir(1,ist,ijsf(iscm))
  300 npri(2,iscm)=npair(2,ijsf(iscm))
      do 310 ist=1,nst
  310 npriri(2,iscm,ist)=nprir(2,ist,ijsf(iscm))
      if(.not.esfk(iscm)) go to 330
      nprk(1,iscm)=npair(1,klsf(iscm))
      do 320 ist=1,nst
  320 nprirk(1,iscm,ist)=nprir(1,ist,klsf(iscm))
  330 nprk(2,iscm)=npair(2,klsf(iscm))
      do 340 ist=1,nst
  340 nprirk(2,iscm,ist)=nprir(2,ist,klsf(iscm))
      nblu=0
      eprsf(iscm)=ijsf(iscm).eq.klsf(iscm)
      do 440 ist=1,nst
      if(.not.esfi(iscm)) go to 380
      nbli=ijb(1,iscm)*npriri(1,iscm,ist)+ipq(ijb(1,iscm))*
     1   npriri(2,iscm,ist)
      go to 390
  380 nbli=ijb(2,iscm)*npriri(2,iscm,ist)
  390 if(.not.eprsf(iscm)) go to 400
      nblu=nblu+ipq(nbli+1)
      go to 440
  400 if(.not.esfk(iscm)) go to 420
      nblk=klb(1,iscm)*nprirk(1,iscm,ist)+ipq(klb(1,iscm))*
     1   nprirk(2,iscm,ist)
      go to 430
  420 nblk=klb(2,iscm)*nprirk(2,iscm,ist)
  430 nblu=nblu+nbli*nblk
  440 continue
      nqd(iscm)=nblu
      niqd(iscm)=nblu
      if(nblu.gt.mblu) then
        nblu=mcxu+nblu
        write (iout,989) nblu
  989   format(' twoint: mblu too small',i8)
      endif
  450 if(iscm.eq.3) go to 600
      iscm=iscm+1
      go to 210
c 460 ----- continue -----
  600 do 610 jesfb=1,2
      ijcxst(jesfb)=icxst(jesfb,ijsf(1))
      ikcxst(jesfb)=icxst(jesfb,ijsf(2))
      ilcxst(jesfb)=icxst(jesfb,ijsf(3))
      klcxst(jesfb)=icxst(jesfb,klsf(1))
      jlcxst(jesfb)=icxst(jesfb,klsf(2))
  610 jkcxst(jesfb)=icxst(jesfb,klsf(3))
      ibufp=ibuf
      ibfldp=ibfld
      do 995 iscm=1,3
      if(niqd(iscm).eq.0) go to 995
      iesfb=2
      kesfb=2
      nblu=nqd(iscm)
      mblu=mblu-nblu
      call wzero(nblu,a,1)
      nwt=1
      inx = 1
      jnx = 2
      knx = 3
      lnx = 4
      go to (1000,2000,3000,4000,5000,6000,7000,8000), ieqsf3
c        isf = jsf = ksf = lsf section
 1000 do 1100 icc = 1, ncc
 1100 call sf3eq1(a,cx,eta,a(nblu+1),icb,icc,is,isf,lmnp1,mcons,nc,ncon,
     1   nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
      go to 996
c        isf = jsf = ksf > lsf section
 2000 do 2990 icc=1,ncc
      if(ieqs3.ne.1) go to 2200
      itst = 0
      igwu = ngw(icc)
      do 2190 iprm = 1, 4
      inx=iperm(1,ip21(iprm))
      jnx=iperm(2,ip21(iprm))
      knx=iperm(3,ip21(iprm))
      lnx=iperm(4,ip21(iprm))
      lc1 = icb(lnx,1,icc)
      if(lc1.eq.itst) go to 2190
      nwt = 1
      if(igwu.eq.1) go to 2180
      ngw1 = 1
      do 2170 igw=2,igwu
      if(icb(lnx,igw,icc)-lc1) 2185,2160,2170
 2160 ngw1 = ngw1 + 1
 2170 nwt = nwt + 1
      nwt = nwt/ngw1
 2180 call sf3eq2(a,cx,eta,a(nblu+1),icb,icc,is,isf,ls,lsf,lmnp1,
     1   mcons,nc,ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
 2185 itst= lc1
 2190 continue
      go to 2990
 2200 call sf3eq2(a,cx,eta,a(nblu+1),icb,icc,is,isf,ls,lsf,lmnp1,
     1   mcons,nc,ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
 2990 continue
      go to 996
c        isf = jsf > ksf = lsf section
 3000 do 3990 icc = 1, ncc
      if(ieqs3.ne.1) go to 3300
      itst = 0
      igwu = ngw(icc)
      do 3119 iprm = 1, 6
      inx=iperm(1,ip31(iprm))
      jnx=iperm(2,ip31(iprm))
      knx=iperm(3,ip31(iprm))
      lnx=iperm(4,ip31(iprm))
      kc1 = icb(knx,1,icc)
      lc1 = icb(lnx,1,icc)
      itst1=ipq(kc1)+lc1
      if (itst1 .le. itst) go to 3119
      nwt = 1
      if (igwu .eq. 1) go to 3220
      ngw1 = 1
      do 3118 igw = 2, igwu
      kc2 = icb(knx,igw,icc)
      lc2 = icb(lnx,igw,icc)
      if(lc2.lt.kc2) go to 3223
      if(ipq(lc2)+kc2-itst1) 3121,3117,3118
3223  if(ipq(kc2)+lc2-itst1) 3121,3117,3118
3117  ngw1 = ngw1 + 1
3118  nwt = nwt + 1
      nwt = nwt/ngw1
 3220 call sf3eq3(a,cx,eta,a(nblu+1),icb,icc,iscm,is,isf,ks,ksf,
     1   lmnp1,mcons,nc,ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
3121  itst= itst1
3119  continue
      go to 3990
 3300 call sf3eq3(a,cx,eta,a(nblu+1),icb,icc,iscm,is,isf,ks,ksf,
     1   lmnp1,mcons,nc,ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
 3990 continue
      go to 996
c        isf = jsf > ksf > lsf section
 4000 do 4990 icc = 1, ncc
      go to (4100,4200,4300,4400), ieqs3
 4100 itst = 0
      igwu = ngw(icc)
      do 4119 iprm = 1, 12
      inx=iperm(1,ip41(iprm))
      jnx=iperm(2,ip41(iprm))
      knx=iperm(3,ip41(iprm))
      lnx=iperm(4,ip41(iprm))
      lc1 = icb(lnx,1,icc)
      kc1 = icb(knx,1,icc)
      itst1 = kc1 + lc1 * icu
      if (itst1 .le. itst) go to 4119
      nwt = 1
      if (igwu .eq. 1) go to 4118
      ngw1 = 1
      do 4117 igw = 2, igwu
      if(icb(knx,igw,icc)+icb(lnx,igw,icc)*icu-itst1) 4120,4116,4117
4116  ngw1 = ngw1 + 1
4117  nwt = nwt + 1
      nwt = nwt/ngw1
 4118 call sf3eq4(a,cx,eta,a(nblu+1),icb,icc,iscm,is,isf,ks,ksf,ls,lsf,
     1   lmnp1,mcons,nc,ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
4120  itst= itst1
4119  continue
      go to 4990
 4200 itst = 0
      igwu = ngw(icc)
      do 4290 iprm = 1, 3
      inx=iperm(1,ip42(iprm))
      jnx=iperm(2,ip42(iprm))
      knx=iperm(3,ip42(iprm))
      lnx=iperm(4,ip42(iprm))
      kc1 = icb(knx,1,icc)
      if(kc1.eq.itst) go to 4290
      nwt = 1
      if(igwu.eq.1) go to 4280
      ngw1 = 1
      do 4270 igw = 2, igwu
      if(icb(knx,igw,icc)-kc1) 4285,4260,4270
 4260 ngw1 = ngw1 + 1
 4270 nwt = nwt + 1
      nwt = nwt /ngw1
 4280 call sf3eq4(a,cx,eta,a(nblu+1),icb,icc,iscm,is,isf,ks,ksf,ls,lsf,
     1   lmnp1,mcons,nc,ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
 4285 itst = kc1
 4290 continue
      go to 4990
 4300 itst = 0
      igwu = ngw(icc)
      do 4390 iprm = 1, 2
      inx=iperm(1,ip43(iprm))
      jnx=iperm(2,ip43(iprm))
      knx=iperm(3,ip43(iprm))
      lnx=iperm(4,ip43(iprm))
      lc1 = icb (lnx,1,icc)
      if(lc1.eq.itst) go to 4390
      nwt = 1
      if (igwu .eq. 1) go to 4380
      ngw1 = 1
      do 4370 igw = 2, igwu
      if(icb(lnx,igw,icc)-lc1) 4385,4360,4370
 4360 ngw1 = ngw1 + 1
 4370 nwt = nwt + 1
      nwt = nwt /ngw1
 4380 call sf3eq4(a,cx,eta,a(nblu+1),icb,icc,iscm,is,isf,ks,ksf,ls,lsf,
     1   lmnp1,mcons,nc,ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
 4385 itst = lc1
 4390 continue
      go to 4990
 4400 call sf3eq4(a,cx,eta,a(nblu+1),icb,icc,iscm,is,isf,ks,ksf,ls,lsf,
     1   lmnp1,mcons,nc,ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
 4990 continue
      go to 996
c        isf > jsf = ksf = lsf section
 5000 do 5990 icc = 1, ncc
      if(ieqs3.ne.1) go to 5500
      itst = 0
      igwu = ngw(icc)
      do 519 iprm = 1, 4
      inx=iperm(1,ip51(iprm))
      jnx=iperm(2,ip51(iprm))
      knx=iperm(3,ip51(iprm))
      lnx=iperm(4,ip51(iprm))
      ic1 = icb(inx,1,icc)
      if (ic1 .eq. itst) go to 519
      nwt = 1
      if (igwu .eq. 1) go to 518
      ngw1 = 1
      do 517 igw = 2, igwu
      if (icb(inx,igw,icc) - ic1) 520,516,517
516   ngw1 = ngw1 + 1
517   nwt = nwt + 1
      nwt = nwt/ngw1
  518 call sf3eq5(a,cx,eta,a(nblu+1),icb,icc,is,isf,js,jsf,lmnp1,
     1   mcons,nc,ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
520   itst = ic1
519   continue
      go to 5990
 5500 call sf3eq5(a,cx,eta,a(nblu+1),icb,icc,is,isf,js,jsf,lmnp1,
     1   mcons,nc,ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
 5990 continue
      go to 996
c        isf > jsf = ksf > lsf section
 6000 do 6900 icc=1,ncc
      go to (6100,6200,6900,6900,6500,6600), ieqs3
 6100 itst = 0
      igwu = ngw(icc)
      do 6119 iprm = 1, 12
      inx=iperm(1,ip61(iprm))
      jnx=iperm(2,ip61(iprm))
      knx=iperm(3,ip61(iprm))
      lnx=iperm(4,ip61(iprm))
      lc1 = icb(lnx,1,icc)
      ic1 = icb(inx,1,icc)
      itst1 = ic1 + lc1 * icu
      if (itst1 .le. itst) go to 6119
      nwt = 1
      if (igwu .eq. 1) go to 6118
      ngw1 = 1
      do 6117 igw = 2, igwu
      if (icb(inx,igw,icc)+icb(lnx,igw,icc)*icu-itst1) 6120,6116,6117
6116  ngw1 = ngw1+ 1
6117  nwt = nwt + 1
      nwt = nwt/ngw1
 6118 call sf3eq6(a,cx,eta,a(nblu+1),icb,icc,iscm,is,isf,js,jsf,ls,lsf,
     1   lmnp1,mcons,nc,ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
6120  itst = itst1
6119  continue
      go to 6900
 6200 itst = 0
      igwu = ngw(icc)
      do 629 iprm = 1, 3
      inx=iperm(1,ip62(iprm))
      jnx=iperm(2,ip62(iprm))
      knx=iperm(3,ip62(iprm))
      lnx=iperm(4,ip62(iprm))
      ic1 = icb(inx,1,icc)
      if (ic1 .eq. itst) go to 629
      nwt = 1
      if (igwu .eq. 1) go to 628
      ngw1 = 1
      do 627 igw = 2, igwu
      if (icb(inx,igw,icc) - ic1) 630,626,627
626   ngw1 = ngw1 + 1
627   nwt = nwt + 1
      nwt = nwt/ngw1
  628 call sf3eq6(a,cx,eta,a(nblu+1),icb,icc,iscm,is,isf,js,jsf,ls,lsf,
     1   lmnp1,mcons,nc,ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
630   itst = ic1
629   continue
      go to 6900
 6500 itst = 0
      igwu = ngw(icc)
      do 659 iprm = 1, 3
      inx=iperm(1,ip65(iprm))
      jnx=iperm(2,ip65(iprm))
      knx=iperm(3,ip65(iprm))
      lnx=iperm(4,ip65(iprm))
      lc1 = icb(lnx,1,icc)
      if (lc1 .eq. itst) go to 659
      nwt = 1
      if (igwu .eq. 1) go to 658
      ngw1 = 1
      do 657 igw = 2, igwu
      if (icb(lnx,igw,icc) - lc1) 660,656,657
656   ngw1 = ngw1 + 1
657   nwt = nwt + 1
      nwt = nwt/ngw1
  658 call sf3eq6(a,cx,eta,a(nblu+1),icb,icc,iscm,is,isf,js,jsf,ls,lsf,
     1   lmnp1,mcons,nc,ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
660   itst = lc1
659   continue
      go to 6900
 6600 call sf3eq6(a,cx,eta,a(nblu+1),icb,icc,iscm,is,isf,js,jsf,ls,lsf,
     1   lmnp1,mcons,nc,ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
 6900 continue
      go to 996
c        isf > jsf > ksf = lsf section
 7000 do 7900 icc=1,ncc
      go to (7100,7900,7300,7900,7500,7900,7700), ieqs3
 7100 itst = 0
      igwu = ngw(icc)
      do 7119 iprm = 1, 12
      inx=iperm(1,ip71(iprm))
      jnx=iperm(2,ip71(iprm))
      knx=iperm(3,ip71(iprm))
      lnx=iperm(4,ip71(iprm))
      jc1 = icb(jnx,1,icc)
      ic1 = icb(inx,1,icc)
      itst1 = ic1 + jc1 * icu
      if (itst1 .le. itst) go to 7119
      nwt = 1
      if (igwu .eq. 1) go to 7118
      ngw1 = 1
      do 7117 igw = 2, igwu
      if(icb(inx,igw,icc)+icb(jnx,igw,icc)*icu-itst1) 7120,7116,7117
7116  ngw1 = ngw1 + 1
7117  nwt = nwt + 1
      nwt = nwt/ngw1
 7118 call sf3eq7(a,cx,eta,a(nblu+1),icb,icc,iscm,is,isf,js,jsf,ks,ksf,
     1   lmnp1,mcons,nc,ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
7120  itst= itst1
7119  continue
      go to 7900
 7300 itst = 0
      igwu = ngw(icc)
      do 7848 iprm = 1, 2
      inx=iperm(1,ip73(iprm))
      jnx=iperm(2,ip73(iprm))
      knx=iperm(3,ip73(iprm))
      lnx=iperm(4,ip73(iprm))
      jc1 = icb(jnx,1,icc)
      if (jc1 .eq. itst) go to 7848
      nwt = 1
      if (igwu .eq. 1) go to 7845
      ngw1 = 1
      do 7844 igw = 2, igwu
      if(icb(jnx,igw,icc) - jc1) 7849,7846,7844
7846  ngw1 = ngw1 + 1
7844  nwt = nwt + 1
      nwt = nwt /ngw1
 7845 call sf3eq7(a,cx,eta,a(nblu+1),icb,icc,iscm,is,isf,js,jsf,ks,ksf,
     1   lmnp1,mcons,nc,ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
7849  itst = jc1
7848  continue
      go to 7900
 7500 itst = 0
      igwu = ngw(icc)
      do 759 iprm = 1, 3
      inx=iperm(1,ip75(iprm))
      jnx=iperm(2,ip75(iprm))
      knx=iperm(3,ip75(iprm))
      lnx=iperm(4,ip75(iprm))
      jc1 = icb(jnx,1,icc)
      if (jc1 .eq. itst) go to 759
      nwt = 1
      if (igwu .eq. 1) go to 758
      ngw1 = 1
      do 757  igw = 2, igwu
      if (icb(jnx,igw,icc) - jc1) 760,756,757
756   ngw1 = ngw1 + 1
757   nwt = nwt + 1
      nwt = nwt/ngw1
  758 call sf3eq7(a,cx,eta,a(nblu+1),icb,icc,iscm,is,isf,js,jsf,ks,ksf,
     1   lmnp1,mcons,nc,ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
760   itst = jc1
759   continue
      go to 7900
 7700 call sf3eq7(a,cx,eta,a(nblu+1),icb,icc,iscm,is,isf,js,jsf,ks,ksf,
     1   lmnp1,mcons,nc,ncon,nfct,ngw,nrcr,nt,ntl,ntu,x,y,z,zet)
 7900 continue
      go to 996
c        isf > jsf > ksf > lsf section
 8000 do 8900 icc=1,ncc
      go to (8100,8200,8300,8400,8500,8600,8700,8800), ieqs3
 8100 itst = 0
      igwu = ngw(icc)
      do 8129 iprm = 1, 24
      inx=iperm(1,iprm)
      jnx=iperm(2,iprm)
      knx=iperm(3,iprm)
      lnx=iperm(4,iprm)
      lc1 = icb(lnx,1,icc)
      kc1 = icb(knx,1,icc)
      jc1 = icb(jnx,1,icc)
      itst1 = jc1 + (kc1 + lc1 * icu) * icu
      if (itst1 .le. itst) go to 8129
      nwt = 1
      if (igwu .eq. 1) go to 8128
      do 8127 igw = 2, igwu
      if(icb(jnx,igw,icc)+(icb(knx,igw,icc)+icb(lnx,igw,icc)*icu)*icu
     1   .lt.itst1) go to 8130
8127  nwt = nwt + 1
 8128 call sf3eq8(a,cx,eta,a(nblu+1),icb,icc,iscm,is,isf,js,jsf,ks,ksf,
     1   ls,lsf,lmnp1,mcons,nc,ncon,nfct,nrcr,nt,ntl,ntu,x,y,z,zet)
8130  itst = itst1
8129  continue
      go to 8900
 8200 itst = 0
      igwu = ngw(icc)
      do 828 iprm = 1, 6
      inx=iperm(1,ip82(iprm))
      jnx=iperm(2,ip82(iprm))
      knx=iperm(3,ip82(iprm))
      lnx=iperm(4,ip82(iprm))
      kc1 = icb(knx,1,icc)
      jc1 = icb(jnx,1,icc)
      itst1 = jc1 + kc1 * icu
      if (itst1 .le. itst) go to 828
      nwt = 1
      if (igwu .eq. 1) go to 829
      do 8221 igw = 2, igwu
      if(icb(jnx,igw,icc)+icb(knx,igw,icc)*icu.lt.itst1) go to 8285
8221  nwt = nwt + 1
  829 call sf3eq8(a,cx,eta,a(nblu+1),icb,icc,iscm,is,isf,js,jsf,ks,ksf,
     1   ls,lsf,lmnp1,mcons,nc,ncon,nfct,nrcr,nt,ntl,ntu,x,y,z,zet)
 8285 itst = itst1
828   continue
      go to 8900
 8300 itst = 0
      igwu = ngw(icc)
      do 830 iprm = 1, 4
      inx=iperm(1,ip83(iprm))
      jnx=iperm(2,ip83(iprm))
      knx=iperm(3,ip83(iprm))
      lnx=iperm(4,ip83(iprm))
      lc1 = icb(lnx,1,icc)
      jc1 = icb(jnx,1,icc)
      itst1 = jc1 + lc1 * icu
      if (itst1 .le. itst) go to 830
      nwt = 1
      if (igwu .eq. 1) go to 836
      do 838 igw = 2, igwu
      if(icb(jnx,igw,icc)+icb(lnx,igw,icc)*icu.lt.itst1) go to 837
838   nwt = nwt + 1
  836 call sf3eq8(a,cx,eta,a(nblu+1),icb,icc,iscm,is,isf,js,jsf,ks,ksf,
     1   ls,lsf,lmnp1,mcons,nc,ncon,nfct,nrcr,nt,ntl,ntu,x,y,z,zet)
837   itst = itst1
830   continue
      go to 8900
 8400 itst = 0
      igwu = ngw(icc)
      do 848 iprm = 1, 2
      inx=iperm(1,ip84(iprm))
      jnx=iperm(2,ip84(iprm))
      knx=iperm(3,ip84(iprm))
      lnx=iperm(4,ip84(iprm))
      jc1 = icb(jnx,1,icc)
      if (jc1 .eq. itst) go to 848
      nwt = 1
      if (igwu .eq. 1) go to 845
      if(icb(jnx,2,icc).lt.jc1) go to 848
      nwt = 2
845   itst = jc1
      call sf3eq8(a,cx,eta,a(nblu+1),icb,icc,iscm,is,isf,js,jsf,ks,ksf,
     1   ls,lsf,lmnp1,mcons,nc,ncon,nfct,nrcr,nt,ntl,ntu,x,y,z,zet)
848   continue
      go to 8900
 8500 itst = 0
      igwu = ngw(icc)
      do 858 iprm = 1, 6
      inx=iperm(1,ip85(iprm))
      jnx=iperm(2,ip85(iprm))
      knx=iperm(3,ip85(iprm))
      lnx=iperm(4,ip85(iprm))
      lc1 = icb(lnx,1,icc)
      kc1 = icb(knx,1,icc)
      itst1 = kc1 + lc1 * jcu
      if (itst1 .le. itst) go to 858
      nwt = 1
      if (igwu .eq. 1) go to 859
      do 8551 igw = 2, igwu
      if(icb(knx,igw,icc)+icb(lnx,igw,icc)*jcu.lt.itst1) go to 8581
8551  nwt = nwt + 1
  859 call sf3eq8(a,cx,eta,a(nblu+1),icb,icc,iscm,is,isf,js,jsf,ks,ksf,
     1   ls,lsf,lmnp1,mcons,nc,ncon,nfct,nrcr,nt,ntl,ntu,x,y,z,zet)
8581  itst = itst1
858   continue
      go to 8900
 8600 itst = 0
      igwu = ngw(icc)
      do 866 iprm = 1, 2
      inx=iperm(1,ip86(iprm))
      jnx=iperm(2,ip86(iprm))
      knx=iperm(3,ip86(iprm))
      lnx=iperm(4,ip86(iprm))
      kc1 = icb(knx,1,icc)
      if (kc1 .eq. itst) go to 866
      nwt = 1
      if (igwu .eq. 1) go to 865
      if(icb(knx,2,icc).lt.kc1) go to 866
      nwt = 2
865   itst = kc1
      call sf3eq8(a,cx,eta,a(nblu+1),icb,icc,iscm,is,isf,js,jsf,ks,ksf,
     1   ls,lsf,lmnp1,mcons,nc,ncon,nfct,nrcr,nt,ntl,ntu,x,y,z,zet)
866   continue
      go to 8900
 8700 itst = 0
      igwu = ngw(icc)
      do 877 iprm = 1, 2
      inx=iperm(1,ip87(iprm))
      jnx=iperm(2,ip87(iprm))
      knx=iperm(3,ip87(iprm))
      lnx=iperm(4,ip87(iprm))
      lc1 = icb (lnx,1,icc)
      if (lc1 .eq. itst) go to 877
      nwt = 1
      if (igwu .eq. 1) go to 875
      if(icb(lnx,2,icc).lt.lc1) go to 877
      nwt = 2
875   itst = lc1
      call sf3eq8(a,cx,eta,a(nblu+1),icb,icc,iscm,is,isf,js,jsf,ks,ksf,
     1   ls,lsf,lmnp1,mcons,nc,ncon,nfct,nrcr,nt,ntl,ntu,x,y,z,zet)
877   continue
      go to 8900
 8800 call sf3eq8(a,cx,eta,a(nblu+1),icb,icc,iscm,is,isf,js,jsf,ks,ksf,
     1   ls,lsf,lmnp1,mcons,nc,ncon,nfct,nrcr,nt,ntl,ntu,x,y,z,zet)
 8900 continue
  996 CALL WTINT2(A,IL,IPRST,ISCM,IBFLD,NSCTLD)
      mblu=mblu+nblu
  995 continue
c     ----- complete set of integrals has been calculated for the -----
c     ----- construction of p and k elements.  set the pk flag.   -----
      if(ibufp.ne.ibuf .or. ibfldp.ne.ibfld) then
        LBL(IBUF*2)=IOR(LBL(IBUF*2),1)
        npkflg=npkflg+1
      endif
  997 continue
  998 continue
 9981 continue
 9982 continue
 9983 continue
 9984 continue
c     ----- write out the last 2-e integral record -----
      iblk=1
      NBUF=IDTR*IBFLD+IBUF
      IBFLD = IBFLD + 1
      IBLK=1
      NSCTLD=NSCTLD+I2SEC(8924)
      CALL SWRIT(NTAPE,IBLK,8942)
C     WRITE (NTAPE) IBLK, IBUF, ALBL, AI
      NSECT=NSECT+NSCTLD
      NSECT=NSECT+1
      WRITE (IOUT,9999)  NBUF, NPKFLG, NSCTLD
      IF (IFLAG3.GT.0) WRITE(ICHECK,9999) NBUF,NPKFLG,NSCTLD
      NINTS=NINTS+NBUF
 9999 format(i20,' integrals and',i9,' pk blocks in',i6,' records')
      write (icheck, * ) 'twoint:  maximum mblu needed ', mxblu
      return
      end
