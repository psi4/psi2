      subroutine cecob(ic,ica,icb,ieqs3,is,js,ks,ls,ncc,nfct,ngw)
      COMMON /OUT/    IOUT,ICHECK,IFLAG3,NSECT,NINTS
      common /parmi/ mblu, mccu, mconu, mcu, mcxu,
     1   mrcru, mstu, msu, nblu, ng, ns, nst
      common /ntgr/itl,itu,jtl,jtu,ktl,ktu,ltl,ltu,jcu,kcu,lcu,inx,
     1jnx,knx,lnx,nwt,nwt1(3),ijsf(3),klsf(3),icxs(2,3),kcxs(2,3),
     2npri(2,3),nprk(2,3),iesfb,kesfb,ircru,jrcru,krcru,lrcru
      dimension ica(mcu,msu,*), icb(4,24,*), nfct(*), ngw(*)
c
      go to (22,24,26,28,30,32,34,36), ieqs3
c        center combinations for is=js=ks=ls
22    do 100 jc = 1, ic
      do 100 kc = 1, jc
      do 100 lc = 1, kc
      igw = 0
      do 70 ig = 1, ng
      ic1 = ica(ic,is,ig)
      jc1 = ica(jc,is,ig)
      kc1 = ica(kc,is,ig)
      lc1 = ica(lc,is,ig)
      if(ic1.ne.ic.and.jc1.ne.ic.and.kc1.ne.ic.and.lc1.ne.ic) go to 70
      ic2 = max(ic1,jc1)
      jc2 = min(ic1,jc1)
      kc2 = max(kc1,lc1)
      lc2 = min(kc1,lc1)
      ic3 = max(ic2,kc2)
      kc2 = min(ic2,kc2)
      lc3 = min(jc2,lc2)
      jc2 = max(jc2,lc2)
      jc3 = max(jc2,kc2)
      kc3 = min(jc2,kc2)
      icc = ncc
50    if(icc .eq. 0) go to 55
      if(jc3 - icb(2,1,icc)) 51, 52, 55
52    if(kc3 - icb(3,1,icc)) 51, 53, 55
53    if(lc3 - icb(4,1,icc)) 51, 54, 55
51    icc = icc - 1
      go to 50
55    if(jc3 .ne. jc .or. kc3 .ne. kc .or. lc3 .ne. lc) go to 70
      igw1 = igw
59    if (igw1 .eq. 0) go to 40
      if (ic1.eq.icb(1,igw1,ncc+1).and.jc1.eq.icb(2,igw1,ncc+1).and.kc1
     x.eq.icb(3,igw1,ncc+1) .and. lc1.eq.icb(4,igw1,ncc+1)) go to 70
      igw1 = igw1 - 1
      go to 59
40    igw = igw + 1
      icb(1,igw,ncc+1) = ic1
      icb(2,igw,ncc+1) = jc1
      icb(3,igw,ncc+1) = kc1
      icb(4,igw,ncc+1) = lc1
70    continue
      ncc = ncc + 1
      nfct(ncc) = ic
      ngw(ncc) = igw
      go to 80
54    nfct(icc) = nfct(icc) + ic
80    continue
100   continue
      do 10  icc = 1, ncc
      jc1  = icb(2,1,icc)
      kc1  = icb(3,1,icc)
      lc1  = icb(4,1,icc)
      idc = 1
      if (jc1 .ne. ic) idc= 2
      if (kc1 .ne. jc1) idc = idc + 1
      if (lc1 .ne. kc1) idc = idc + 1
10    nfct(icc) = nfct(icc)/idc
      go to 900
c        center combinations for is=js=ks>ls
24    do 29 jc=1,ic
      do 29 kc=1,jc
      do 29 lc=1,lcu
      igw = 0
      do 27 ig = 1, ng
      ic1 = ica(ic,is,ig)
      jc1 = ica(jc,is,ig)
      kc1 = ica(kc,is,ig)
      if(ic1.ne.ic.and.jc1.ne.ic.and.kc1.ne.ic) go to 27
      lc1 = ica(lc,ls,ig)
      ic2 = max(ic1,jc1)
      jc2 = min(ic1,jc1)
      ic3 = max(ic2,kc1)
      kc2 = min(ic2,kc1)
      jc3 = max(jc2,kc2)
      kc3 = min(jc2,kc2)
      icc = ncc
213   if(icc .eq. 0) go to 212
      if(jc3 - icb(2,1,icc)) 214,233,212
233   if(kc3 - icb(3,1,icc)) 214,235,212
235   if(lc1 - icb(4,1,icc)) 214,215,212
214   icc = icc - 1
      go to 213
212   if(jc3 .ne. jc .or. kc3 .ne. kc .or. lc1 .ne. lc) go to 27
      igw1 = igw
237   if (igw1 .eq. 0) go to 239
      if (ic1.eq.icb(1,igw1,ncc+1).and.jc1.eq.icb(2,igw1,ncc+1).and.kc1
     x.eq.icb(3,igw1,ncc+1)) go to 27
      igw1 = igw1 - 1
      go to 237
239   igw = igw + 1
      icb(1,igw,ncc+1) = ic1
      icb(2,igw,ncc+1) = jc1
      icb(3,igw,ncc+1) = kc1
      icb(4,igw,ncc+1) = lc1
27    continue
      ncc = ncc + 1
      nfct(ncc) = ic
      ngw(ncc) = igw
      go to 29
215   nfct(icc) = nfct(icc) + ic
29    continue
      do 217 icc = 1, ncc
      jc1  = icb(2,1,icc)
      kc1  = icb(3,1,icc)
      idc = 1
      if ( jc1 .ne. ic) idc = 2
      if (kc1 .ne. jc1) idc = idc + 1
217   nfct(icc) = nfct(icc)/idc
      go to 900
c        center combinations for is=js>ks=ls
26    do 379 jc = 1, ic
      do 378 kc = 1, kcu
      do 377 lc = 1, kc
      igw = 0
      do 370 ig = 1, ng
      ic1 = ica(ic,is,ig)
      jc1 = ica(jc,is,ig)
      if(ic1.ne.ic.and.jc1.ne.ic) go to 370
      kc1 = ica(kc,ks,ig)
      lc1 = ica(lc,ks,ig)
      ic2 = max(ic1,jc1)
      jc2 = min(ic1,jc1)
      kc2 = max(kc1,lc1)
      lc2 = min(kc1,lc1)
      icc = ncc
372   if (icc .eq. 0) go to 376
      if (jc2 - icb(2,1,icc)) 375,381,376
381   if (kc2 - icb(3,1,icc)) 375,382,376
382   if (lc2 - icb(4,1,icc)) 375,385,376
375   icc = icc -1
      go to 372
376   if (jc2 .ne. jc .or. kc2 .ne. kc .or. lc2 .ne. lc) go to 370
      igw1 = igw
374   if (igw1 .eq. 0) go to 373
      if (ic1 .eq.icb(1,igw1,ncc+1).and.jc1.eq.icb(2,igw1,ncc+1).and.
     xkc1.eq.icb(3,igw1,ncc+1).and.lc1.eq.icb(4,igw1,ncc+1)) go to 370
      igw1 = igw1 - 1
      go to 374
373   igw = igw + 1
      icb (1,igw,ncc+1) = ic1
      icb (2,igw,ncc+1) = jc1
      icb (3,igw,ncc+1) = kc1
      icb (4,igw,ncc+1) = lc1
370   continue
      ncc = ncc + 1
      nfct(ncc) =  ic
      ngw(ncc) = igw
      go to 377
385   nfct(icc) = nfct(icc) + ic
377   continue
378   continue
379   continue
      go to 384
c        center combinations for is=js>ks>ls
28    do 410 jc = 1, ic
      do 410 kc = 1, kcu
      do 410 lc = 1, lcu
      igw = 0
      do 400 ig = 1, ng
      ic1 = ica(ic,is,ig)
      jc1 = ica(jc,is,ig)
      if(ic1.ne.ic.and.jc1.ne.ic) go to 400
      kc1 = ica(kc,ks,ig)
      lc1 = ica(lc,ls,ig)
      ic2 = max(ic1,jc1)
      jc2 = min(ic1,jc1)
      icc = ncc
402   if (icc .eq. 0) go to 406
      if (jc2 - icb(2,1,icc)) 405,403,406
403   if (kc1 - icb(3,1,icc)) 405,404,406
404   if (lc1 - icb(4,1,icc)) 405,413,406
405   icc = icc - 1
      go to 402
406   if(jc2 .ne. jc .or. kc1 .ne. kc .or. lc1 .ne. lc) go to 400
      igw1 = igw
411   if (igw1 .eq. 0) go to 412
      if (ic1 .eq. icb(1,igw1,ncc+1).and.jc1.eq.icb(2,igw1,ncc+1)) go
     xto 400
      igw1 = igw1- 1
      go to 411
412   igw = igw + 1
      icb(1,igw,ncc+1) = ic1
      icb(2,igw,ncc+1) = jc1
      icb(3,igw,ncc+1) = kc1
      icb(4,igw,ncc+1) = lc1
400   continue
      ncc = ncc + 1
      nfct(ncc) = ic
      ngw(ncc) = igw
      go to 410
413   nfct(icc) = nfct(icc) + ic
410   continue
384   do 414 icc = 1, ncc
 414  if(icb(2,1,icc).ne.ic) nfct(icc)=nfct(icc)/2
      go to 900
c        center combinations for is>js=ks=ls
30    do 520 jc = 1, jcu
      do 520 kc = 1, jc
      do 520 lc = 1, kc
      igw =0
      do 515 ig = 1, ng
      if (ica(ic,is,ig) .ne. ic) go to 515
      jc1 = ica(jc,js,ig)
      kc1 = ica(kc,js,ig)
      lc1 = ica(lc,js,ig)
      jc2 = max (jc1,kc1)
      kc2 = min (jc1,kc1)
      jc3 = max (jc2,lc1)
      lc2 = min (jc2,lc1)
      kc3 = max (kc2,lc2)
      lc3 = min (kc2,lc2)
c     do 509 icc = ncc, 1, -1
      icc = ncc
506   if (icc .eq. 0) go to 512
      if (jc3 - icb(2,1,icc)) 509, 507, 512
507   if (kc3 - icb(3,1,icc)) 509, 508, 512
508   if (lc3 - icb(4,1,icc)) 509, 516, 512
509   icc = icc - 1
      go to 506
512   if (jc3.ne.jc.or.kc3.ne.kc.or.lc3.ne.lc) go to 515
c     do --- igw1 = igw, 1, -1
      igw1 = igw
511   if (igw1 .eq. 0) go to 513
      if (jc1.eq.icb(2,igw1,ncc+1).and.kc1.eq.icb(3,igw1,ncc+1).and.lc1
     x.eq.icb(4,igw1,ncc+1)) go to 515
      igw1 = igw1 - 1
      go to 511
513   igw = igw + 1
      icb (1,igw,ncc+1) = ic
      icb (2,igw,ncc+1) = jc1
      icb (3,igw,ncc+1) = kc1
      icb (4,igw,ncc+1) = lc1
515   continue
      ncc = ncc + 1
      nfct(ncc) = ic
      ngw(ncc) = igw
      go to 520
516   nfct(icc) = nfct(icc) + ic
520   continue
      go to 900
c        center combinations for is>js=ks>ls
32    do 220 jc = 1, jcu
      do 220 kc = 1, jc
      do 220 lc = 1, lcu
      igw = 0
      do 210 ig = 1, ng
      if (ica(ic,is,ig) .ne. ic) go to 210
      jc1 = ica(jc,js,ig)
      kc1 = ica(kc,js,ig)
      lc1 = ica(lc,ls,ig)
      jc2 = max (jc1, kc1)
      kc2 = min (jc1, kc1)
c     do 204 icc= ncc, 1, -1
      icc = ncc
201   if (icc .eq. 0) go to 206
      if (jc2 - icb(2,1,icc)) 204, 202, 206
202   if (kc2 - icb(3,1,icc)) 204, 203, 206
203   if (lc1 - icb(4,1,icc)) 204, 211, 206
204   icc = icc - 1
      go to 201
206   if (jc2.ne.jc.or.kc2.ne.kc.or.lc1.ne.lc) go to 210
c     do --- igw1 = igw, 1, -1
      igw1 = igw
205   if (igw1 .eq. 0) go to 207
      if (jc1.eq.icb(2,igw1,ncc+1).and.kc1.eq.icb(3,igw1,ncc+1)) go
     cto 210
      igw1 = igw1 - 1
      go to 205
207   igw = igw + 1
      icb (1,igw,ncc+1) = ic
      icb (2,igw,ncc+1) = jc1
      icb (3,igw,ncc+1) = kc1
      icb (4,igw,ncc+1) = lc1
210   continue
      ncc = ncc + 1
      nfct(ncc) = ic
      ngw (ncc) = igw
      go to 220
211   nfct(icc) = nfct(icc) + ic
220   continue
      go to 900
c        center combinations for is>js>ks=ls
34    do 720 jc = 1, jcu
      do 720 kc = 1, kcu
      do 720 lc = 1, kc
      igw = 0
      do 710 ig = 1, ng
      if (ica(ic,is,ig) .ne. ic) go to 710
      jc1 = ica(jc,js,ig)
      kc1 = ica(kc,ks,ig)
      lc1 = ica(lc,ks,ig)
      kc2 = max (kc1, lc1)
      lc2 = min (kc1, lc1)
c     do 704 icc = ncc, 1, -1
      icc = ncc
701   if (icc .eq. 0) go to 700
      if (jc1 - icb(2,1,icc)) 704, 702, 700
702   if (kc2 - icb(3,1,icc)) 704, 703, 700
703   if (lc2 - icb(4,1,icc)) 704, 708, 700
704   icc = icc - 1
      go to 701
700   if (jc1.ne.jc.or.kc2.ne.kc.or.lc2.ne.lc) go to 710
c     do --- igw1 = igw, 1, -1
      igw1 = igw
705   if (igw1 .eq. 0) go to 706
      if (kc1.eq.icb(3,igw1,ncc+1).and.lc1.eq.icb(4,igw1,ncc+1)) go
     cto 710
      igw1 = igw1 - 1
      go to 705
706   igw = igw + 1
      icb (1,igw,ncc+1) = ic
      icb (2,igw,ncc+1) = jc1
      icb (3,igw,ncc+1) = kc1
      icb (4,igw,ncc+1) = lc1
710   continue
      ncc = ncc + 1
      nfct(ncc) = ic
      ngw(ncc) = igw
      go to 720
708   nfct(icc) = nfct(icc) + ic
720   continue
      go to 900
c        center combinations for is>js>ks>ls
36    do 820 jc = 1, jcu
      do 820 kc = 1, kcu
      do 820 lc = 1, lcu
      do 810 ig = 1, ng
      if (ica(ic,is,ig) .ne. ic) go to 810
      jc1 = ica(jc,js,ig)
      kc1 = ica(kc,ks,ig)
      lc1 = ica(lc,ls,ig)
c     do 804 icc = ncc, 1, -1
      icc = ncc
801   if (icc .eq. 0) go to 810
      if (jc1 - icb(2,1,icc)) 804, 802, 810
802   if (kc1 - icb(3,1,icc)) 804, 803, 810
803   if (lc1 - icb(4,1,icc)) 804, 808, 810
804   icc = icc - 1
      go to 801
810   continue
      ncc = ncc + 1
      icb (1, 1, ncc) = ic
      icb (2, 1, ncc) = jc
      icb (3, 1, ncc) = kc
      icb (4, 1, ncc) = lc
      nfct(ncc) = ic
      go to 820
808   nfct(icc) = nfct(icc) + ic
820   continue
  900 IF(NCC.GT.MCCU) THEN
        WRITE(IOUT,988) NCC
        STOP
  988   FORMAT(' MCCU TOO SMALL',I4)
      END IF
      return
      end
