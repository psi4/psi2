      subroutine socfpd(c,cx,icxst,icxsv1,icxsv2,idp,il,
     1   iprst,la,lb,maords,mau,mcons,mgcs,mics,mjcs,nblpr,nc,nct,
     2   nd,nf,nir,npair,nprir,nrcr,nt,ntl,ntu)
*mdc*on vax ibm cray
      implicit real*8 (a-h,o-z)
*mdc*off
      logical es, esf, esfb, esfbc, msfbct
      COMMON /OUT/    IOUT,ICHECK,IFLAG3,NSECT,NINTS
      common /parmi/ mblu, mccu, mconu, mcu, mcxu,
     1   mrcru, mstu, msu, nblu, ng, ns, nst
      common /ntgr/ kaords, mconsu, mgcsu, mru, mcsu, mctu, mcru,
     1   msfu, mgu, msfru, mnru, mnsfu, mpru, ngcs, nu, mcxu2,
     2   lxyzir(3), inam(5), nnam, mdum(30)
      common /dim21/ ipq(256)
      dimension c(mctu,mcsu,*),cx(*),icxst(2,*),icxsv1(*),
     1   icxsv2(mgcsu,*),idp(mstu,mstu,*),il(*),iprst(*),la(mru,*),
     2   lb(*),maords(*),mau(*),mcons(*),mgcs(*),mics(*),mjcs(*),
     3   nblpr(*),nc(*),nct(*),nd(*),nf(*),nir(*),npair(2,*),
     4   nprir(2,mstu,*),nrcr(*),nt(*),ntl(*),ntu(*)
      parameter ( a0=0.0d0 )
c
      icx=0
      ipair=0
c
c  ijsf = one dimensionalization of is, js, if, jf
c  isfr = one dimensionalization of is, if, ir
c
      ijsf=0
      isf=0
      isfr=0
      do 10 igcs=1,ngcs
      icxsv1(igcs)=0
      do 10 jgcs=1,ngcs
   10 icxsv2(jgcs,igcs)=0
      do 340 is=1,ns
      ifu=nf(is)
      if(ifu.eq.0) go to 340
      icu=nc(is)
      isfis=isf
      isfris=isfr
      jsf=0
      jsfr=0
      do 330 js=1,is
      jsfjs=jsf
      jsfrjs=jsfr
      es=is.eq.js
      if(es) go to 15
      jfu=nf(js)
      if(jfu.eq.0) go to 330
      jcu=nc(js)
   15 isf=isfis
      isfr=isfris
      do 320 if=1,ifu
      isf=isf+1
      isfrif=isfr
      igcs=mgcs(isf)
      iaords=maords(igcs)
      iru=nir(iaords)
      ircru=nrcr(mcons(isf))
      itl=ntl(isf)
      itu=ntu(isf)
      jsf=jsfjs
      jsfr=jsfrjs
      if(es) jfu=if
      do 320 jf=1,jfu
      jsf=jsf+1
      jgcs=mgcs(jsf)
      jaords=maords(jgcs)
      jtl=ntl(jsf)
      jsfrjf=jsfr
      ijsf=ijsf+1
      esf=isf.eq.jsf
      if(.not.esf) go to 25
      if(ircru.ge.2) go to 30
      npair(2,ijsf)=0
      icxst(2,ijsf)=0
      do 20 ist=1,nst
   20 nprir(2,ist,ijsf)=0
      go to 35
   25 jrcru=nrcr(mcons(jsf))
   30 jtu=ntu(jsf)
      jru=nir(jaords)
   35 iprst(ijsf)=ipair
      isfr=isfrif
      do 320 ircr=1,ircru
      isfrib=isfr
      jsfr=jsfrjf
      if(esf) jrcru=ircr
      do 310 jrcr=1,jrcru
      jsfrjb=jsfr
      esfb=esf.and.ircr.eq.jrcr
      if(.not.esfb) go to 40
      iesfb=1
      go to 45
   40 iesfb=2
   45 npr=0
      do 75 ist=1,nst
      ipr=npr
      ics=0
      isfr=isfrib
      do 70 ir=1,iru
      isfr=isfr+1
      lai=la(ir,iaords)
      iau=nd(lai)
*mdc*on ibm vax
      ILPR=IOR(ISHFT(LAI-1,8),LB(ISFR))
*mdc*off
      if(esfb) jru=ir
      jsfr=jsfrjb
      jcs=0
      do 65 jr=1,jru
      jsfr=jsfr+1
      IF (JSFR.GT.MSFRU) THEN
        WRITE(IOUT,901) MSFRU
  901   FORMAT (//' *** MSFRU TOO SMALL IN SOCFPD:',I7//)
        call qabort
      END IF
      laj=la(jr,jaords)
      if(idp(lai,laj,ist)) 50,65,50
   50 ipair =ipair+1
      if(lai.lt.laj) then
*mdc*on ibm vax
*      il(ipair)=ior(ishft(ior(ishft(laj,11),lb(jsfr)),15),
*     +                    ior(ishft(lai,11),lb(isfr)))
      IL(IPAIR)=IOR(IOR(ISHFT(LAJ-1,19),ISHFT(LB(JSFR),11)),ILPR)
*mdc*off
*mdc*on harris
*      il(ipair)=((((laj-1).shift.8).or.lb(jsfr)).shift.11).or.
*     +           (((lai-1).shift.8).or.lb(isfr))
*mdc*off
*mdc*on cray
*      il(ipair)=or(shiftl(or(shiftl(laj,11),lb(jsfr)),15),
*     +                    or(shiftl(lai,11),lb(isfr)))
*mdc*off
c
      else
c
*mdc*on ibm vax
*      il(ipair)=ior(ishft(ior(ishft(lai,11),lb(isfr)),15),
*     +                    ior(ishft(laj,11),lb(jsfr)))
      IL(IPAIR)=IOR(ISHFT(ILPR,11),IOR(ISHFT(LAJ-1,8),LB(JSFR)))
*mdc*off
*mdc*on harris
*      il(ipair)=((((lai-1).shift.8).or.lb(isfr)).shift.11).or.
*     +           (((laj-1).shift.8).or.lb(jsfr))
*mdc*off
*mdc*on crayosu
*      il(ipair)=or(shiftl(or(shiftl(lai,11),lb(isfr)),15),
*     +                    or(shiftl(laj,11),lb(jsfr)))
*mdc*off
c
      endif
c
      npr=npr+1
      IF (NPR.GT.MNRU) THEN
        WRITE(IOUT,903) NPR
  903   FORMAT (//' *** NPR TOO SMALL IN SOCFPD',I7//)
        call qabort
      END IF
      mau(npr)=iau
      mics(npr)=ics
      mjcs(npr)=jcs
   65 jcs=jcs+nd(laj)
   70 ics=ics+iau
   75 nprir(iesfb,ist,ijsf)=npr-ipr
      npair(iesfb,ijsf)=npr
      if(npr.eq.0) go to 310
      if(.not.esfb) go to 85
      isv1=icxsv1(igcs)
      if(isv1.eq.0) go to 80
      icxst(1,ijsf)=icxst(1,isv1)
      go to 310
   80 icxsv1(igcs)=ijsf
      icxst(1,ijsf)=icx
      icxu=icx+(npr*ipq(icu*nt(isf)+1))
      go to 105
   85 isv2=icxsv2(jgcs,igcs)
      if(isv2.eq.0) go to 90
      icxst(2,ijsf)=icxst(2,isv2)
      go to 310
   90 icxsv2(jgcs,igcs)=ijsf
      icxst(2,ijsf)=icx
      icxu=icx+npr*icu*nt(isf)*jcu*nt(jsf)
  105 IF(ICXU.GT.MCXU) THEN
        WRITE(IOUT,985) ICXU,MCXU
  985   FORMAT(' MCXU TOO SMALL',2I9)
        call qabort
      END IF
c
c  compute new cx block
c
      ict=0
      do 290 ic=1,icu
      ictic=ict
      if(esfb) jcu=ic
      jct=0
      do 290 jc=1,jcu
      jctjc=jct
      esfbc=esfb.and.ic.eq.jc
      ict=ictic
      do 290 it=itl,itu
      ict=ict+1
      if(esfbc) jtu=it
      jct=jctjc
      do 290 jt=jtl,jtu
      jct=jct+1
      msfbct=.not.esfb.or.(ic.eq.jc.and.it.eq.jt)
      do 290 ipr=1,npr
      icx=icx+1
      cx(icx)= a0
      iau=mau(ipr)
      ics=mics(ipr)
      jcs=mjcs(ipr)
      do 280 ia=1,iau
      ics=ics+1
      jcs=jcs+1
      cx(icx)=cx(icx)+(c(ict,ics,igcs)*c(jct,jcs,jgcs))
      if(msfbct) go to 280
      cx(icx)=cx(icx)+(c(jct,ics,igcs)*c(ict,jcs,jgcs))
  280 continue
  290 continue
  310 continue
  320 continue
  330 continue
  340 continue
      IF(IJSF.GT.MNSFU) THEN
        WRITE(IOUT,983) IJSF
  983   FORMAT(' MNSFU TOO SMALL',I5)
        call qabort
      END IF
      IF (IPAIR.GT.MPRU) THEN
        WRITE(IOUT,902) MPRU
  902   FORMAT (//' *** MPRU TOO SMALL IN SOCFPD:',I7//)
        call qabort
      END IF
      WRITE(ICHECK,900) ICX, IJSF, IPAIR
  900 FORMAT(' MCXU=',I8,' MNSFU=',I8,' MPRU=',I7)
C
      MNSFU=IJSF
      MPRU=IPAIR
      MCXU=ICX
C
      return
      end
