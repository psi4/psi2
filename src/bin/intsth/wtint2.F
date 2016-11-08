      subroutine wtint2(a,il,iprst,iscm,ibfld,nsctld)
*mdc*on vax ibm cray
      implicit real*8 (a-h,o-z)
*mdc*off
      logical eprsf, esfc, esfcij, esfckl, esfi, esfk
      COMMON /OUT/    IOUT,ICHECK,IFLAG3,NSECT,NINTS
      COMMON /TAPE/ NTAPE
      common /parmr/ cutoff, tol
      common /parmi/ mblu, mccu, mconu, mcu, mcxu,
     1   mrcru, mstu, msu, nblu, ng, ns, nst
      common /lgcl/ eprsf(3), esfi(3), esfk(3)
      common /ntgr/itl,itu,jtl,jtu,ktl,ktu,ltl,ltu,jcu,kcu,lcu,inx,
     1jnx,knx,lnx,nwt,nwt1(3),ijsf(3),klsf(3),icxs(2,3),kcxs(2,3),
     2npri(2,3),nprk(2,3),mdum(6)
c
      COMMON /OUTBUF/ IBLK, IBUF, ALBL(2235), AI(2235)
      common /ikpr/ npriri(2,3,8), nprirk(2,3,8)
      common/misc/xij,xijm,yij,yijm,zij,zijm,xkl,xklm,ykl,yklm,zkl,zklm,
     1icons,ircru,jcons,jrcru,kcons,krcru,lcons,lrcru,esfc,esfcij,esfckl
      dimension a(*), il(*), iprst(*)
      DIMENSION LBL(4470)
      EQUIVALENCE (ALBL,LBL)
      DATA IDTR/2235/
c
      ibl = 0
      ipr = iprst(ijsf(iscm))
      if(.not.eprsf(iscm)) go to 10
      kprib=ipr
      go to 20
   10 kprib=iprst(klsf(iscm))
   20 do 170 ircr=1,ircru
      if(esfi(iscm)) jrcru=ircr
      do 170 jrcr=1,jrcru
      if((.not.esfi(iscm)).or.ircr.ne.jrcr) go to 50
      iesfb=1
      go to 60
   50 iesfb=2
   60 iprkb=ipr
      kpr=kprib
      if(eprsf(iscm)) krcru=ircr
      do 170 krcr=1,krcru
      if(.not.eprsf(iscm)) go to 94
      if(ircr.ne.krcr) go to 90
      lrcru=jrcr
      go to 100
   90 if(.not.esfk(iscm)) go to 92
      lrcru=krcr
      go to 100
   92 lrcru=jrcru
      go to 100
   94 if(esfk(iscm)) lrcru=krcr
  100 do 170 lrcr=1,lrcru
      if((.not.esfk(iscm)).or.krcr.ne.lrcr) go to 110
      kesfb=1
      go to 120
  110 kesfb=2
  120 ipr=iprkb
      do 160 ist=1,nst
      npairi = npriri(iesfb,iscm,ist)
      npairk = nprirk(kesfb,iscm,ist)
      if(npairi.ne.0) go to 30
      kpr=kpr+npairk
      go to 160
   30 if(npairk.ne.0) go to 40
      ipr=ipr + npairi
      go to 160
   40 kpri=kpr
      do 150 ipair=1,npairi
      ipr=ipr + 1
      if((eprsf(iscm).and.ircr.eq.krcr).and.jrcr.eq.lrcr) npairk=ipair
      kpr=kpri
      do 140 kpair=1,npairk
      kpr=kpr + 1
      ibl=ibl + 1
ctph  if(ibfld.ge.1)write(6,*)' ibl cutoff',ibl,cutoff
      if(abs(a(ibl)).le.cutoff) go to 140
      IF(IBUF.EQ.IDTR) THEN
        IBFLD=IBFLD+1
        NSCTLD=NSCTLD+I2SEC(8942)
        CALL SWRIT(NTAPE,IBLK,8942)
C       WRITE (NTAPE) IBLK, IBUF, ALBL, AI
        IBUF=0
      END IF
      ibuf=ibuf+1
      ijlb=il(ipr)
      kllb=il(kpr)
c
      if(kllb.ge.ijlb) then
        LBL(2*IBUF-1)=KLLB
        LBL(2*IBUF  )=ISHFT(IJLB,2)
      else
        LBL(2*IBUF-1)=IJLB
        LBL(2*IBUF  )=ISHFT(KLLB,2)
      endif
      AI(IBUF) = A(IBL)
c
  140 continue
  150 continue
  160 continue
  170 continue
      return
      end
