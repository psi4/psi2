      subroutine aoso2e(a,cx,iscm,val)
*mdc*on vax ibm cray
      implicit real*8 (a-h,o-z)
*mdc*off
      common /parmi/ mblu, mccu, mconu, mcu, mcxu,
     1   mrcru, mstu, msu, nblu, ng, ns, nst
      common /ntgr/itl,itu,jtl,jtu,ktl,ktu,ltl,ltu,jcu,kcu,lcu,inx,
     1jnx,knx,lnx,nwt,nwt1(3),ijsf(3),klsf(3),icxs(2,3),kcxs(2,3),
     2npri(2,3),nprk(2,3),iesfb,kesfb,ircru,jrcru,krcru,lrcru
      common /ikpr/ npriri(2,3,8), nprirk(2,3,8)
      common /sf3ao/ ibl
      dimension a(*), cx(*)
      icx = icxs(iesfb,iscm)
      jcx = kcxs(iesfb,iscm)
      if(icx.ne.jcx) go to 24
      assign 72 to iaoe
      go to 32
   24 assign 64 to iaoe
   32 do 92 ist=1,nst
      npairi = npriri(iesfb,iscm,ist)
      if(npairi.eq.0) goto 92
      kcxst = jcx
      lcxst = icx
      do 88 ijbl=1,npairi
      icx=icx+1
      jcx=jcx+1
      if(cx(icx)) 56,40,56
   40 if(cx(jcx)) 56,48,56
   48 ibl=ibl+ijbl
      go to 88
   56 kcx = kcxst
      lcx = lcxst
      do 80 klbl=1,ijbl
      kcx=kcx+1
      prd = cx(icx) * cx(kcx)
      go to iaoe, (64,72)
   64 lcx=lcx+1
      prd=prd+cx(jcx)*cx(lcx)
   72 ibl=ibl+1
   80 a(ibl)=a(ibl)+prd*val
   88 continue
   92 continue
      return
      end
