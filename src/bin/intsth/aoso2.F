      subroutine aoso2(a,cx,iscm,val)
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
      kcxst = kcxs(kesfb,iscm)
      do 56 ist=1,nst
      npairi = npriri(iesfb,iscm,ist)
      npairk = nprirk(kesfb,iscm,ist)
      if(npairi.eq.0) goto 52
      if(npairk.ne.0) goto 20
      icx=icx+npairi
      go to 56
   20 do 48 ijbl=1,npairi
      icx=icx+1
      if(cx(icx)) 32,48,32
   32 prd=cx(icx)*val
      do 40 klbl=1,npairk
   40 a(ibl+klbl)=a(ibl+klbl)+cx(kcxst+klbl)*prd
   48 ibl=ibl+npairk
   52 kcxst = kcxst + npairk
   56 continue
      return
      end
