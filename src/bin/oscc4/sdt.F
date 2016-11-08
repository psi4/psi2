      subroutine sdt(t1,t2,din,ein,fin,fock,
     .               orbsym,flov,fzo,fzv,toff,tadd,zlx,ipq,
     .               itr,itv,eoff,eadd,foff,fadd,doff,dadd,
     .               no,nv,nirred,xet,xet4,trstr,
     .               map,spin,ssov,omap)
      implicit integer(a-z)
      real*8 din(*),ein(*),fin(*),t1(no,nv),t2(*),fock(*)
      integer orbsym(*),flov(8,2,4),fzo(*),fzv(*),
     .        toff(no,no),tadd(nv,nv),zlx(nv,nv),
     .        itr(*),itv(*),doff(*),dadd(*),eoff(*),eadd(*),
     .        foff(*),fadd(*),ipq(*),
     .        map(*),spin(*),ssov(2,4),spen(2),omap(*)
      character*4 trstr
      real*8 val,valf,vale,val1,val2,delta,xet,xet4,fac,xval,aval,bval,
     .       third,valv,valw,ival,jval,kval,sign,xet1,xet2,xet3
c
      write(4,*)'  '
      write(4,*)'Entering ccsd(t) energy calculation  '
      write(6,*)'  '
      write(6,*)'Entering ccsd(t) energy calculation  '
c
      third=1.0d0/3.0d0
      xet=0.0d0
      xet1=0.0d0
      xet2=0.0d0
      xet3=0.0d0
      xet4=0.0d0
c
      if(trstr.eq.'YES ')then
         call ffile(9,'trstr',1)
         backspace 9
         backspace 9
         backspace 9
         read(9,'(i4)')lmo
         read(9,'(e20.15)')xet
         write(6,*)'restarting calculation from virtual mo',lmo
         write(6,'(a32,f20.12)')'restarting calculation with xe',xet
         write(6,*)
      else
         call ffile(9,'trstr',0)
         lmo=0
      endif
c
      do 950 a=1,nv
         if(a.le.lmo)go to 950
         if(fzv(a).eq.1)go to 950
         aspn=spin(a+no)
         xa=map(a+no)
         asym=orbsym(a+no)
         aa=ipq(a+no)+a+no
         do 850 b=1,a-1
            if(fzv(b).eq.1)go to 850
            bspn=spin(b+no)
            xb=map(b+no)
            bsym=orbsym(b+no)
            absym=ieor(asym,bsym)
            bb=ipq(b+no)+b+no
            do 750 c=1,b-1 
               if(fzv(c).eq.1)go to 750
               cspn=spin(c+no)
               xc=map(c+no)
               cc=ipq(c+no)+c+no
               csym=orbsym(c+no)
               abcsym=ieor(absym,csym)
c
               case=aspn+bspn+cspn-2
c
               do 650 i=1,no
                  if(fzo(i).eq.1)go to 650
                  ispn=spin(i)
                  if(ispn.ne.aspn)go to 650
                  xi=map(i)
                  isym=orbsym(i)
                  ia=ipq(max0(xi,xa))+min0(xi,xa)
                  ib=ipq(max0(xi,xb))+min0(xi,xb)
                  ic=ipq(max0(xi,xc))+min0(xi,xc)
                  iasym=ieor(isym,asym)
                  ibsym=ieor(isym,bsym)
                  icsym=ieor(isym,csym)
                  ii=ipq(i)+i
                  do 550 j=1,i-1
                     if(fzo(j).eq.1)go to 550
                     jspn=spin(j)
                     if(jspn.ne.bspn)go to 550
                     xj=map(j)
                     jsym=orbsym(j)
                     ijsym=ieor(isym,jsym)
                     jasym=ieor(jsym,asym)
                     jbsym=ieor(jsym,bsym)
                     jcsym=ieor(jsym,csym)
                     ja=ipq(max0(xj,xa))+min0(xj,xa)
                     jb=ipq(max0(xj,xb))+min0(xj,xb)
                     jc=ipq(max0(xj,xc))+min0(xj,xc)
                     jj=ipq(j)+j
                     do 450 k=1,j-1
                        if(fzo(k).eq.1)go to 450
                        kspn=spin(k)
                        if(kspn.ne.cspn)go to 450
                        xk=map(k)
                        ksym=orbsym(k)
                        ijksym=ieor(ijsym,ksym)
                        if(ijksym.ne.abcsym)go to 450
                        kasym=ieor(ksym,asym)
                        kbsym=ieor(ksym,bsym)
                        kcsym=ieor(ksym,csym)
                        ka=ipq(max0(xk,xa))+min0(xk,xa)
                        kb=ipq(max0(xk,xb))+min0(xk,xb)
                        kc=ipq(max0(xk,xc))+min0(xk,xc)
                        kk=ipq(k)+k
c
      valv=0.0d0
      valw=0.0d0
      valf=0.0d0
      vale=0.0d0
c
      if (case.eq.1.or.case.eq.4)then
      spn=ispn
c
      jbkc=doff(max0(jb,kc))+dadd(min0(jb,kc))
      jckb=doff(max0(jc,kb))+dadd(min0(jc,kb))
      ibkc=doff(max0(ib,kc))+dadd(min0(ib,kc))
      ickb=doff(max0(ic,kb))+dadd(min0(ic,kb))
      icjb=doff(max0(ic,jb))+dadd(min0(ic,jb))
      ibjc=doff(max0(ib,jc))+dadd(min0(ib,jc))
      jakc=doff(max0(ja,kc))+dadd(min0(ja,kc))
      jcka=doff(max0(jc,ka))+dadd(min0(jc,ka))
      iakc=doff(max0(ia,kc))+dadd(min0(ia,kc))
      icka=doff(max0(ic,ka))+dadd(min0(ic,ka))
      icja=doff(max0(ic,ja))+dadd(min0(ic,ja))
      iajc=doff(max0(ia,jc))+dadd(min0(ia,jc))
      jbka=doff(max0(jb,ka))+dadd(min0(jb,ka))
      jakb=doff(max0(ja,kb))+dadd(min0(ja,kb))
      ibka=doff(max0(ib,ka))+dadd(min0(ib,ka))
      iakb=doff(max0(ia,kb))+dadd(min0(ia,kb))
      iajb=doff(max0(ia,jb))+dadd(min0(ia,jb))
      ibja=doff(max0(ib,ja))+dadd(min0(ib,ja))

      valv=(din(jbkc)-din(jckb))*t1(i,a)
     .    -(din(ibkc)-din(ickb))*t1(j,a)
     .    -(din(icjb)-din(ibjc))*t1(k,a)
     .    -(din(jakc)-din(jcka))*t1(i,b)
     .    +(din(iakc)-din(icka))*t1(j,b)
     .    +(din(icja)-din(iajc))*t1(k,b)
     .    -(din(jbka)-din(jakb))*t1(i,c)
     .    +(din(ibka)-din(iakb))*t1(j,c)
     .    +(din(iajb)-din(ibja))*t1(k,c)
c               
c
c     lf2=lf
c     if(lf.ge.c)lf2=c-1
c     do 351 f=ff,lf2
c     ff2=ff
c     if(ff.le.c)ff2=c+1
c     do 361 f=ff2,lf
c
c
      fsym=ieor(bsym,iasym)+1
      ff=flov(fsym,spn,3)-no
      lf=flov(fsym,spn,4)-no
      do 351 f=ff,lf
         if(f.eq.c)go to 351
         xf=map(f+no)
         fb=ipq(max0(xf,xb))+min0(xf,xb)
         fa=ipq(max0(xf,xa))+min0(xf,xa)
         iafb=foff(max0(ia,fb))+fadd(min0(ia,fb))
         ibfa=foff(max0(ib,fa))+fadd(min0(ib,fa))
         ival=+fin(iafb)-fin(ibfa)
         jkfc=toff(j,k)+tadd(f,c)
         val=+ival*t2(jkfc)
         sign=+1.0d0
         if(c.gt.f)sign=-1.0d0
         valf=valf+val*sign
  351 continue
c
      fsym=ieor(bsym,jasym)+1
      ff=flov(fsym,spn,3)-no
      lf=flov(fsym,spn,4)-no
      do 352 f=ff,lf
         if(f.eq.c)go to 352
         xf=map(f+no)
         fb=ipq(max0(xf,xb))+min0(xf,xb)
         fa=ipq(max0(xf,xa))+min0(xf,xa)
         jafb=foff(max0(ja,fb))+fadd(min0(ja,fb))
         jbfa=foff(max0(jb,fa))+fadd(min0(jb,fa))
         jval=-fin(jafb)+fin(jbfa)
         ikfc=toff(i,k)+tadd(f,c)
         val=+jval*t2(ikfc)
         sign=+1.0d0
         if(c.gt.f)sign=-1.0d0
         valf=valf+val*sign
  352 continue
c
      fsym=ieor(bsym,kasym)+1
      ff=flov(fsym,spn,3)-no
      lf=flov(fsym,spn,4)-no
      do 353 f=ff,lf
         if(f.eq.c)go to 353
         xf=map(f+no)
         fb=ipq(max0(xf,xb))+min0(xf,xb)
         fa=ipq(max0(xf,xa))+min0(xf,xa)
         kafb=foff(max0(ka,fb))+fadd(min0(ka,fb))
         kbfa=foff(max0(kb,fa))+fadd(min0(kb,fa))
         kval=+fin(kafb)-fin(kbfa)
         ijfc=toff(i,j)+tadd(f,c)
         val=+kval*t2(ijfc)
         sign=+1.0d0
         if(c.gt.f)sign=-1.0d0
         valf=valf+val*sign
  353 continue
c
      fsym=ieor(bsym,icsym)+1
      ff=flov(fsym,spn,3)-no
      lf=flov(fsym,spn,4)-no
      do 354 f=ff,lf
         if(f.eq.a)go to 354
         xf=map(f+no)
         fb=ipq(max0(xf,xb))+min0(xf,xb)
         fc=ipq(max0(xf,xc))+min0(xf,xc)
         icfb=foff(max0(ic,fb))+fadd(min0(ic,fb))
         ibfc=foff(max0(ib,fc))+fadd(min0(ib,fc))
         ival=-fin(icfb)+fin(ibfc)
         jkfa=toff(j,k)+tadd(f,a)
         val=+ival*t2(jkfa)
         sign=+1.0d0
         if(a.gt.f)sign=-1.0d0
         valf=valf+val*sign
  354 continue
c
      fsym=ieor(bsym,jcsym)+1
      ff=flov(fsym,spn,3)-no
      lf=flov(fsym,spn,4)-no
      do 355 f=ff,lf
         if(f.eq.a)go to 355
         xf=map(f+no)
         fb=ipq(max0(xf,xb))+min0(xf,xb)
         fc=ipq(max0(xf,xc))+min0(xf,xc)
         jcfb=foff(max0(jc,fb))+fadd(min0(jc,fb))
         jbfc=foff(max0(jb,fc))+fadd(min0(jb,fc))
         jval=+fin(jcfb)-fin(jbfc)
         ikfa=toff(i,k)+tadd(f,a)
         val=+jval*t2(ikfa)
         sign=+1.0d0
         if(a.gt.f)sign=-1.0d0
         valf=valf+val*sign
  355 continue
c
      fsym=ieor(bsym,kcsym)+1
      ff=flov(fsym,spn,3)-no
      lf=flov(fsym,spn,4)-no
      do 356 f=ff,lf
         if(f.eq.a)go to 356
         xf=map(f+no)
         fb=ipq(max0(xf,xb))+min0(xf,xb)
         fc=ipq(max0(xf,xc))+min0(xf,xc)
         kcfb=foff(max0(kc,fb))+fadd(min0(kc,fb))
         kbfc=foff(max0(kb,fc))+fadd(min0(kb,fc))
         kval=-fin(kcfb)+fin(kbfc)
         ijfa=toff(i,j)+tadd(f,a)
         val=+kval*t2(ijfa)
         sign=+1.0d0
         if(a.gt.f)sign=-1.0d0
         valf=valf+val*sign
  356 continue
c
      fsym=ieor(csym,iasym)+1
      ff=flov(fsym,spn,3)-no
      lf=flov(fsym,spn,4)-no
      do 357 f=ff,lf
         if(f.eq.b)go to 357
         xf=map(f+no)
         fc=ipq(max0(xf,xc))+min0(xf,xc)
         fa=ipq(max0(xf,xa))+min0(xf,xa)
         iafc=foff(max0(ia,fc))+fadd(min0(ia,fc))
         icfa=foff(max0(ic,fa))+fadd(min0(ic,fa))
         ival=-fin(iafc)+fin(icfa)
         jkfb=toff(j,k)+tadd(f,b)
         val=+ival*t2(jkfb)
         sign=+1.0d0
         if(b.gt.f)sign=-1.0d0
         valf=valf+val*sign
  357 continue
c
      fsym=ieor(csym,jasym)+1
      ff=flov(fsym,spn,3)-no
      lf=flov(fsym,spn,4)-no
      do 358 f=ff,lf
         if(f.eq.b)go to 358
         xf=map(f+no)
         fc=ipq(max0(xf,xc))+min0(xf,xc)
         fa=ipq(max0(xf,xa))+min0(xf,xa)
         jafc=foff(max0(ja,fc))+fadd(min0(ja,fc))
         jcfa=foff(max0(jc,fa))+fadd(min0(jc,fa))
         jval=+fin(jafc)-fin(jcfa)
         ikfb=toff(i,k)+tadd(f,b)
         val=+jval*t2(ikfb)
         sign=+1.0d0
         if(b.gt.f)sign=-1.0d0
         valf=valf+val*sign
  358 continue
c
      fsym=ieor(csym,kasym)+1
      ff=flov(fsym,spn,3)-no
      lf=flov(fsym,spn,4)-no
      do 359 f=ff,lf
         if(f.eq.b)go to 359
         xf=map(f+no)
         fc=ipq(max0(xf,xc))+min0(xf,xc)
         fa=ipq(max0(xf,xa))+min0(xf,xa)
         kafc=foff(max0(ka,fc))+fadd(min0(ka,fc))
         kcfa=foff(max0(kc,fa))+fadd(min0(kc,fa))
         kval=-fin(kafc)+fin(kcfa)
         ijfb=toff(i,j)+tadd(f,b)
         val=+kval*t2(ijfb)
         sign=+1.0d0
         if(b.gt.f)sign=-1.0d0
         valf=valf+val*sign
  359 continue
c
c
      msym=ieor(jsym,iasym)+1
      fm=flov(msym,spn,1)
      lm=flov(msym,spn,2)
      do 251 m=fm,lm
         if(m.eq.k)go to 251
         xm=map(m)
         mi=ipq(max0(xm,xi))+min0(xm,xi)
         mj=ipq(max0(xm,xj))+min0(xm,xj)
         mjia=eoff(max0(mj,ia))+eadd(min0(mj,ia))
         mija=eoff(max0(mi,ja))+eadd(min0(mi,ja))
         val=-ein(mjia)+ein(mija)
         sign=1.0d0
         if(k.gt.m)sign=-1.0d0
         mkbc=toff(m,k)+tadd(b,c)
         vale=vale+val*t2(mkbc)*sign
  251 continue
c
      msym=ieor(jsym,kasym)+1
      fm=flov(msym,spn,1)
      lm=flov(msym,spn,2)
      do 252 m=fm,lm
         if(m.eq.i)go to 252
         xm=map(m)
         mk=ipq(max0(xm,xk))+min0(xm,xk)
         mj=ipq(max0(xm,xj))+min0(xm,xj)
         mjka=eoff(max0(mj,ka))+eadd(min0(mj,ka))
         mkja=eoff(max0(mk,ja))+eadd(min0(mk,ja))
         val=+ein(mjka)-ein(mkja)
         sign=1.0d0
         if(i.gt.m)sign=-1.0d0
         mibc=toff(m,i)+tadd(b,c)
         vale=vale+val*t2(mibc)*sign
  252 continue
c
      msym=ieor(ksym,iasym)+1
      fm=flov(msym,spn,1)
      lm=flov(msym,spn,2)
      do 253 m=fm,lm
         if(m.eq.j)go to 253
         xm=map(m)
         mi=ipq(max0(xm,xi))+min0(xm,xi)
         mk=ipq(max0(xm,xk))+min0(xm,xk)
         mkia=eoff(max0(mk,ia))+eadd(min0(mk,ia))
         mika=eoff(max0(mi,ka))+eadd(min0(mi,ka))
         val=+ein(mkia)-ein(mika)
         sign=1.0d0
         if(j.gt.m)sign=-1.0d0
         mjbc=toff(m,j)+tadd(b,c)
         vale=vale+val*t2(mjbc)*sign
  253 continue
c
      msym=ieor(jsym,ibsym)+1
      fm=flov(msym,spn,1)
      lm=flov(msym,spn,2)
      do 254 m=fm,lm
         if(m.eq.k)go to 254
         xm=map(m)
         mi=ipq(max0(xm,xi))+min0(xm,xi)
         mj=ipq(max0(xm,xj))+min0(xm,xj)
         mjib=eoff(max0(mj,ib))+eadd(min0(mj,ib))
         mijb=eoff(max0(mi,jb))+eadd(min0(mi,jb))
         val=+ein(mjib)-ein(mijb)
         sign=1.0d0
         if(k.gt.m)sign=-1.0d0
         mkac=toff(m,k)+tadd(a,c)
         vale=vale+val*t2(mkac)*sign
  254 continue
c
      msym=ieor(jsym,kbsym)+1
      fm=flov(msym,spn,1)
      lm=flov(msym,spn,2)
      do 255 m=fm,lm
         if(m.eq.i)go to 255
         xm=map(m)
         mk=ipq(max0(xm,xk))+min0(xm,xk)
         mj=ipq(max0(xm,xj))+min0(xm,xj)
         mjkb=eoff(max0(mj,kb))+eadd(min0(mj,kb))
         mkjb=eoff(max0(mk,jb))+eadd(min0(mk,jb))
         val=-ein(mjkb)+ein(mkjb)
         sign=1.0d0
         if(i.gt.m)sign=-1.0d0
         miac=toff(m,i)+tadd(a,c)
         vale=vale+val*t2(miac)*sign
  255 continue
c
      msym=ieor(ksym,ibsym)+1
      fm=flov(msym,spn,1)
      lm=flov(msym,spn,2)
      do 256 m=fm,lm
         if(m.eq.j)go to 256
         xm=map(m)
         mi=ipq(max0(xm,xi))+min0(xm,xi)
         mk=ipq(max0(xm,xk))+min0(xm,xk)
         mkib=eoff(max0(mk,ib))+eadd(min0(mk,ib))
         mikb=eoff(max0(mi,kb))+eadd(min0(mi,kb))
         val=-ein(mkib)+ein(mikb)
         sign=1.0d0
         if(j.gt.m)sign=-1.0d0
         mjac=toff(m,j)+tadd(a,c)
         vale=vale+val*t2(mjac)*sign
  256 continue
c
      msym=ieor(jsym,icsym)+1
      fm=flov(msym,spn,1)
      lm=flov(msym,spn,2)
      do 257 m=fm,lm
         if(m.eq.k)go to 257
         xm=map(m)
         mi=ipq(max0(xm,xi))+min0(xm,xi)
         mj=ipq(max0(xm,xj))+min0(xm,xj)
         mjic=eoff(max0(mj,ic))+eadd(min0(mj,ic))
         mijc=eoff(max0(mi,jc))+eadd(min0(mi,jc))
         val=-ein(mjic)+ein(mijc)
         sign=1.0d0
         if(k.gt.m)sign=-1.0d0
         mkab=toff(m,k)+tadd(a,b)
         vale=vale+val*t2(mkab)*sign
  257 continue
c
      msym=ieor(jsym,kcsym)+1
      fm=flov(msym,spn,1)
      lm=flov(msym,spn,2)
      do 258 m=fm,lm
         if(m.eq.i)go to 258
         xm=map(m)
         mk=ipq(max0(xm,xk))+min0(xm,xk)
         mj=ipq(max0(xm,xj))+min0(xm,xj)
         mjkc=eoff(max0(mj,kc))+eadd(min0(mj,kc))
         mkjc=eoff(max0(mk,jc))+eadd(min0(mk,jc))
         val=+ein(mjkc)-ein(mkjc)
         sign=1.0d0
         if(i.gt.m)sign=-1.0d0
         miab=toff(m,i)+tadd(a,b)
         vale=vale+val*t2(miab)*sign
  258 continue
c
      msym=ieor(ksym,icsym)+1
      fm=flov(msym,spn,1)
      lm=flov(msym,spn,2)
      do 259 m=fm,lm
         if(m.eq.j)go to 259
         xm=map(m)
         mi=ipq(max0(xm,xi))+min0(xm,xi)
         mk=ipq(max0(xm,xk))+min0(xm,xk)
         mkic=eoff(max0(mk,ic))+eadd(min0(mk,ic))
         mikc=eoff(max0(mi,kc))+eadd(min0(mi,kc))
         val=+ein(mkic)-ein(mikc)
         sign=1.0d0
         if(j.gt.m)sign=-1.0d0
         mjab=toff(m,j)+tadd(a,b)
         vale=vale+val*t2(mjab)*sign
  259 continue
c
  997 continue
      endif
c
c >>>>>>>>> case =2  <<<<<<<<<<<<<<<
c
      if (case.eq.2)then     
c
      jbkc=doff(max0(jb,kc))+dadd(min0(jb,kc))
      jckb=doff(max0(jc,kb))+dadd(min0(jc,kb))
      iakc=doff(max0(ia,kc))+dadd(min0(ia,kc))
      iajc=doff(max0(ia,jc))+dadd(min0(ia,jc))
      iakb=doff(max0(ia,kb))+dadd(min0(ia,kb))
      iajb=doff(max0(ia,jb))+dadd(min0(ia,jb))

      valv=(din(jbkc)-din(jckb))*t1(i,a)
     .    +din(iakc)*t1(j,b)
     .    -din(iajc)*t1(k,b)
     .    -din(iakb)*t1(j,c)
     .    +din(iajb)*t1(k,c)
c               
c
      fsym=ieor(bsym,iasym)+1
      ff=flov(fsym,1,3)-no
      lf=flov(fsym,1,4)-no
      do 551 f=ff,lf
         if(f.eq.c)go to 551
         xf=map(f+no)
         fb=ipq(max0(xf,xb))+min0(xf,xb)
         iafb=foff(max0(ia,fb))+fadd(min0(ia,fb))
         ival=+fin(iafb)
         jkfc=toff(j,k)+tadd(f,c)
         val=+ival*t2(jkfc)
         sign=+1.0d0
         if(c.gt.f)sign=-1.0d0
         valf=valf+val*sign
  551 continue
c
      fsym=ieor(bsym,jasym)+1
      ff=flov(fsym,2,3)-no
      lf=flov(fsym,2,4)-no
      do 552 f=ff,lf
         xf=map(f+no)
         fa=ipq(max0(xf,xa))+min0(xf,xa)
         jbfa=foff(max0(jb,fa))+fadd(min0(jb,fa))
         jval=+fin(jbfa)
         ikfc=toff(i,k)+tadd(f,c)
         val=+jval*t2(ikfc)
         valf=valf+val
  552 continue
c
      fsym=ieor(bsym,kasym)+1
      ff=flov(fsym,2,3)-no
      lf=flov(fsym,2,4)-no
      do 553 f=ff,lf
         xf=map(f+no)
         fa=ipq(max0(xf,xa))+min0(xf,xa)
         kbfa=foff(max0(kb,fa))+fadd(min0(kb,fa))
         kval=-fin(kbfa)
         ijfc=toff(i,j)+tadd(f,c)
         val=+kval*t2(ijfc)
         valf=valf+val
  553 continue
c
      fsym=ieor(bsym,jcsym)+1
      ff=flov(fsym,1,3)-no
      lf=flov(fsym,1,4)-no
      do 555 f=ff,lf
         xf=map(f+no)
         fb=ipq(max0(xf,xb))+min0(xf,xb)
         fc=ipq(max0(xf,xc))+min0(xf,xc)
         jcfb=foff(max0(jc,fb))+fadd(min0(jc,fb))
         jbfc=foff(max0(jb,fc))+fadd(min0(jb,fc))
         jval=-fin(jcfb)+fin(jbfc)
         ikaf=toff(i,k)+tadd(a,f)
         val=+jval*t2(ikaf)
         valf=valf+val
  555 continue
c
      fsym=ieor(bsym,kcsym)+1
      ff=flov(fsym,1,3)-no
      lf=flov(fsym,1,4)-no
      do 556 f=ff,lf
         xf=map(f+no)
         fb=ipq(max0(xf,xb))+min0(xf,xb)
         fc=ipq(max0(xf,xc))+min0(xf,xc)
         kcfb=foff(max0(kc,fb))+fadd(min0(kc,fb))
         kbfc=foff(max0(kb,fc))+fadd(min0(kb,fc))
         kval=+fin(kcfb)-fin(kbfc)
         ijaf=toff(i,j)+tadd(a,f)
         val=+kval*t2(ijaf)
         valf=valf+val
  556 continue
c
      fsym=ieor(csym,iasym)+1
      ff=flov(fsym,1,3)-no
      lf=flov(fsym,1,4)-no
      do 557 f=ff,lf
         if(f.eq.b)go to 557
         xf=map(f+no)
         fc=ipq(max0(xf,xc))+min0(xf,xc)
         iafc=foff(max0(ia,fc))+fadd(min0(ia,fc))
         ival=-fin(iafc)
         jkfb=toff(j,k)+tadd(f,b)
         val=+ival*t2(jkfb)
         sign=+1.0d0
         if(b.gt.f)sign=-1.0d0
         valf=valf+val*sign
  557 continue
c
      fsym=ieor(csym,jasym)+1
      ff=flov(fsym,2,3)-no
      lf=flov(fsym,2,4)-no
      do 558 f=ff,lf
         xf=map(f+no)
         fa=ipq(max0(xf,xa))+min0(xf,xa)
         jcfa=foff(max0(jc,fa))+fadd(min0(jc,fa))
         jval=-fin(jcfa)
         ikfb=toff(i,k)+tadd(f,b)
         val=+jval*t2(ikfb)
         valf=valf+val
  558 continue
c
      fsym=ieor(csym,kasym)+1
      ff=flov(fsym,2,3)-no
      lf=flov(fsym,2,4)-no
      do 559 f=ff,lf
         xf=map(f+no)
         fa=ipq(max0(xf,xa))+min0(xf,xa)
         kcfa=foff(max0(kc,fa))+fadd(min0(kc,fa))
         kval=+fin(kcfa)
         ijfb=toff(i,j)+tadd(f,b)
         val=+kval*t2(ijfb)
         valf=valf+val
  559 continue
c
c
      msym=ieor(jsym,iasym)+1
      fm=flov(msym,1,1)
      lm=flov(msym,1,2)
      do 451 m=fm,lm
         if(m.eq.k)go to 451
         xm=map(m)
         mj=ipq(max0(xm,xj))+min0(xm,xj)
         mjia=eoff(max0(mj,ia))+eadd(min0(mj,ia))
         val=-ein(mjia)
         sign=1.0d0
         if(k.gt.m)sign=-1.0d0
         mkbc=toff(m,k)+tadd(b,c)
         vale=vale+val*t2(mkbc)*sign
  451 continue
c
      msym=ieor(ksym,iasym)+1
      fm=flov(msym,1,1)
      lm=flov(msym,1,2)
      do 453 m=fm,lm
         if(m.eq.j)go to 453
         xm=map(m)
         mk=ipq(max0(xm,xk))+min0(xm,xk)
         mkia=eoff(max0(mk,ia))+eadd(min0(mk,ia))
         val=+ein(mkia)
         sign=1.0d0
         if(j.gt.m)sign=-1.0d0
         mjbc=toff(m,j)+tadd(b,c)
         vale=vale+val*t2(mjbc)*sign
  453 continue
c
      msym=ieor(jsym,ibsym)+1
      fm=flov(msym,2,1)
      lm=flov(msym,2,2)
      do 454 m=fm,lm
         xm=map(m)
         mi=ipq(max0(xm,xi))+min0(xm,xi)
         mijb=eoff(max0(mi,jb))+eadd(min0(mi,jb))
         val=-ein(mijb)
         mkac=toff(m,k)+tadd(a,c)
         vale=vale+val*t2(mkac)
  454 continue
c
      msym=ieor(jsym,kbsym)+1
      fm=flov(msym,1,1)
      lm=flov(msym,1,2)
      do 455 m=fm,lm
         xm=map(m)
         mk=ipq(max0(xm,xk))+min0(xm,xk)
         mj=ipq(max0(xm,xj))+min0(xm,xj)
         mjkb=eoff(max0(mj,kb))+eadd(min0(mj,kb))
         mkjb=eoff(max0(mk,jb))+eadd(min0(mk,jb))
         val=+ein(mjkb)-ein(mkjb)
         imac=toff(i,m)+tadd(a,c)
         vale=vale+val*t2(imac)
  455 continue
c
      msym=ieor(ksym,ibsym)+1
      fm=flov(msym,2,1)
      lm=flov(msym,2,2)
      do 456 m=fm,lm
         xm=map(m)
         mi=ipq(max0(xm,xi))+min0(xm,xi)
         mikb=eoff(max0(mi,kb))+eadd(min0(mi,kb))
         val=+ein(mikb)
         mjac=toff(m,j)+tadd(a,c)
         vale=vale+val*t2(mjac)
  456 continue
c
      msym=ieor(jsym,icsym)+1
      fm=flov(msym,2,1)
      lm=flov(msym,2,2)
      do 457 m=fm,lm
         xm=map(m)
         mi=ipq(max0(xm,xi))+min0(xm,xi)
         mijc=eoff(max0(mi,jc))+eadd(min0(mi,jc))
         val=+ein(mijc)
         mkab=toff(m,k)+tadd(a,b)
         vale=vale+val*t2(mkab)
  457 continue
c
      msym=ieor(jsym,kcsym)+1
      fm=flov(msym,1,1)
      lm=flov(msym,1,2)
      do 458 m=fm,lm
         xm=map(m)
         mk=ipq(max0(xm,xk))+min0(xm,xk)
         mj=ipq(max0(xm,xj))+min0(xm,xj)
         mjkc=eoff(max0(mj,kc))+eadd(min0(mj,kc))
         mkjc=eoff(max0(mk,jc))+eadd(min0(mk,jc))
         val=-ein(mjkc)+ein(mkjc)
         imab=toff(i,m)+tadd(a,b)
         vale=vale+val*t2(imab)
  458 continue
c
      msym=ieor(ksym,icsym)+1
      fm=flov(msym,2,1)
      lm=flov(msym,2,2)
      do 459 m=fm,lm
         xm=map(m)
         mi=ipq(max0(xm,xi))+min0(xm,xi)
         mikc=eoff(max0(mi,kc))+eadd(min0(mi,kc))
         val=-ein(mikc)
         mjab=toff(m,j)+tadd(a,b)
         vale=vale+val*t2(mjab)
  459 continue
c
  991 continue
      endif
c
c
c >>>>>>>>> case =3  <<<<<<<<<<<<<<<
c
      if (case.eq.3)then     
      jbkc=doff(max0(jb,kc))+dadd(min0(jb,kc))
      ibkc=doff(max0(ib,kc))+dadd(min0(ib,kc))
      jakc=doff(max0(ja,kc))+dadd(min0(ja,kc))
      iakc=doff(max0(ia,kc))+dadd(min0(ia,kc))
      iajb=doff(max0(ia,jb))+dadd(min0(ia,jb))
      ibja=doff(max0(ib,ja))+dadd(min0(ib,ja))
c
      valv= din(jbkc)           *t1(i,a)
     .    - din(ibkc)           *t1(j,a)
     .    - din(jakc)           *t1(i,b)
     .    + din(iakc)           *t1(j,b)
     .    +(din(iajb)-din(ibja))*t1(k,c)
 
c
      fsym=ieor(bsym,iasym)+1
      ff=flov(fsym,2,3)-no
      lf=flov(fsym,2,4)-no
      do 751 f=ff,lf
         xf=map(f+no)
         fb=ipq(max0(xf,xb))+min0(xf,xb)
         fa=ipq(max0(xf,xa))+min0(xf,xa)
         iafb=foff(max0(ia,fb))+fadd(min0(ia,fb))
         ibfa=foff(max0(ib,fa))+fadd(min0(ib,fa))
         ival=+fin(iafb)-fin(ibfa)
         jkfc=toff(j,k)+tadd(f,c)
         val=+ival*t2(jkfc)
         valf=valf+val
  751 continue
c
      fsym=ieor(bsym,jasym)+1
      ff=flov(fsym,2,3)-no
      lf=flov(fsym,2,4)-no
      do 752 f=ff,lf
         xf=map(f+no)
         fb=ipq(max0(xf,xb))+min0(xf,xb)
         fa=ipq(max0(xf,xa))+min0(xf,xa)
         jafb=foff(max0(ja,fb))+fadd(min0(ja,fb))
         jbfa=foff(max0(jb,fa))+fadd(min0(jb,fa))
         jval=-fin(jafb)+fin(jbfa)
         ikfc=toff(i,k)+tadd(f,c)
         val=+jval*t2(ikfc)
         valf=valf+val
  752 continue
c
      fsym=ieor(bsym,icsym)+1
      ff=flov(fsym,1,3)-no
      lf=flov(fsym,1,4)-no
      do 754 f=ff,lf
         xf=map(f+no)
         fc=ipq(max0(xf,xc))+min0(xf,xc)
         ibfc=foff(max0(ib,fc))+fadd(min0(ib,fc))
         ival=-fin(ibfc)
         jkaf=toff(j,k)+tadd(a,f)
         val=+ival*t2(jkaf)
         valf=valf+val
  754 continue
c
      fsym=ieor(bsym,jcsym)+1
      ff=flov(fsym,1,3)-no
      lf=flov(fsym,1,4)-no
      do 755 f=ff,lf
         xf=map(f+no)
         fc=ipq(max0(xf,xc))+min0(xf,xc)
         jbfc=foff(max0(jb,fc))+fadd(min0(jb,fc))
         jval=+fin(jbfc)
         ikaf=toff(i,k)+tadd(a,f)
         val=+jval*t2(ikaf)
         valf=valf+val
  755 continue
c
      fsym=ieor(bsym,kcsym)+1
      ff=flov(fsym,2,3)-no
      lf=flov(fsym,2,4)-no
      do 756 f=ff,lf
         if(f.eq.a)go to 756
         xf=map(f+no)
         fb=ipq(max0(xf,xb))+min0(xf,xb)
         kcfb=foff(max0(kc,fb))+fadd(min0(kc,fb))
         kval=-fin(kcfb)
         ijfa=toff(i,j)+tadd(f,a)
         val=+kval*t2(ijfa)
         sign=+1.0d0
         if(a.gt.f)sign=-1.0d0
         valf=valf+val*sign
  756 continue
c
      fsym=ieor(csym,iasym)+1
      ff=flov(fsym,1,3)-no
      lf=flov(fsym,1,4)-no
      do 757 f=ff,lf
         xf=map(f+no)
         fc=ipq(max0(xf,xc))+min0(xf,xc)
         iafc=foff(max0(ia,fc))+fadd(min0(ia,fc))
         ival=+fin(iafc)
         jkbf=toff(j,k)+tadd(b,f)
         val=+ival*t2(jkbf)
         valf=valf+val
  757 continue
c
      fsym=ieor(csym,jasym)+1
      ff=flov(fsym,1,3)-no
      lf=flov(fsym,1,4)-no
      do 758 f=ff,lf
         xf=map(f+no)
         fc=ipq(max0(xf,xc))+min0(xf,xc)
         jafc=foff(max0(ja,fc))+fadd(min0(ja,fc))
         jval=-fin(jafc)
         ikbf=toff(i,k)+tadd(b,f)
         val=+jval*t2(ikbf)
         valf=valf+val
  758 continue
c
      fsym=ieor(csym,kasym)+1
      ff=flov(fsym,2,3)-no
      lf=flov(fsym,2,4)-no
      do 759 f=ff,lf
         if(f.eq.b)go to 759
         xf=map(f+no)
         fa=ipq(max0(xf,xa))+min0(xf,xa)
         kcfa=foff(max0(kc,fa))+fadd(min0(kc,fa))
         kval=+fin(kcfa)
         ijfb=toff(i,j)+tadd(f,b)
         val=+kval*t2(ijfb)
         sign=+1.0d0
         if(b.gt.f)sign=-1.0d0
         valf=valf+val*sign
  759 continue
c
c
      msym=ieor(jsym,iasym)+1
      fm=flov(msym,2,1)
      lm=flov(msym,2,2)
      do 651 m=fm,lm
         xm=map(m)
         mi=ipq(max0(xm,xi))+min0(xm,xi)
         mj=ipq(max0(xm,xj))+min0(xm,xj)
         mjia=eoff(max0(mj,ia))+eadd(min0(mj,ia))
         mija=eoff(max0(mi,ja))+eadd(min0(mi,ja))
         val=-ein(mjia)+ein(mija)
         mkbc=toff(m,k)+tadd(b,c)
         vale=vale+val*t2(mkbc)
  651 continue
c
      msym=ieor(jsym,kasym)+1
      fm=flov(msym,1,1)
      lm=flov(msym,1,2)
      do 652 m=fm,lm
         xm=map(m)
         mk=ipq(max0(xm,xk))+min0(xm,xk)
         mkja=eoff(max0(mk,ja))+eadd(min0(mk,ja))
         val=+ein(mkja)
         imbc=toff(i,m)+tadd(b,c)
         vale=vale+val*t2(imbc)
  652 continue
c
      msym=ieor(ksym,iasym)+1
      fm=flov(msym,1,1)
      lm=flov(msym,1,2)
      do 653 m=fm,lm
         xm=map(m)
         mk=ipq(max0(xm,xk))+min0(xm,xk)
         mkia=eoff(max0(mk,ia))+eadd(min0(mk,ia))
         val=-ein(mkia)
         jmbc=toff(j,m)+tadd(b,c)
         vale=vale+val*t2(jmbc)
  653 continue
c
      msym=ieor(jsym,ibsym)+1
      fm=flov(msym,2,1)
      lm=flov(msym,2,2)
      do 654 m=fm,lm
         xm=map(m)
         mi=ipq(max0(xm,xi))+min0(xm,xi)
         mj=ipq(max0(xm,xj))+min0(xm,xj)
         mjib=eoff(max0(mj,ib))+eadd(min0(mj,ib))
         mijb=eoff(max0(mi,jb))+eadd(min0(mi,jb))
         val=+ein(mjib)-ein(mijb)
         mkac=toff(m,k)+tadd(a,c)
         vale=vale+val*t2(mkac)
  654 continue
c
      msym=ieor(jsym,kbsym)+1
      fm=flov(msym,1,1)
      lm=flov(msym,1,2)
      do 655 m=fm,lm
         xm=map(m)
         mk=ipq(max0(xm,xk))+min0(xm,xk)
         mkjb=eoff(max0(mk,jb))+eadd(min0(mk,jb))
         val=-ein(mkjb)
         imac=toff(i,m)+tadd(a,c)
         vale=vale+val*t2(imac)
  655 continue
c
      msym=ieor(ksym,ibsym)+1
      fm=flov(msym,1,1)
      lm=flov(msym,1,2)
      do 656 m=fm,lm
         xm=map(m)
         mk=ipq(max0(xm,xk))+min0(xm,xk)
         mkib=eoff(max0(mk,ib))+eadd(min0(mk,ib))
         val=+ein(mkib)
         jmac=toff(j,m)+tadd(a,c)
         vale=vale+val*t2(jmac)
  656 continue
c
      msym=ieor(jsym,kcsym)+1
      fm=flov(msym,2,1)
      lm=flov(msym,2,2)
      do 658 m=fm,lm
         if(m.eq.i)go to 658
         xm=map(m)
         mj=ipq(max0(xm,xj))+min0(xm,xj)
         mjkc=eoff(max0(mj,kc))+eadd(min0(mj,kc))
         val=+ein(mjkc)
         sign=1.0d0
         if(i.gt.m)sign=-1.0d0
         miab=toff(m,i)+tadd(a,b)
         vale=vale+val*t2(miab)*sign
  658 continue
c
      msym=ieor(ksym,icsym)+1
      fm=flov(msym,2,1)
      lm=flov(msym,2,2)
      do 659 m=fm,lm
         if(m.eq.j)go to 659
         xm=map(m)
         mi=ipq(max0(xm,xi))+min0(xm,xi)
         mikc=eoff(max0(mi,kc))+eadd(min0(mi,kc))
         val=-ein(mikc)
         sign=1.0d0
         if(j.gt.m)sign=-1.0d0
         mjab=toff(m,j)+tadd(a,b)
         vale=vale+val*t2(mjab)*sign
  659 continue
  992 continue
c
      endif
c
      valw=vale+valf
      delta=fock(ii)+fock(jj)+fock(kk)-fock(aa)-fock(bb)-fock(cc)
      xval=(valw+valv)*valw/delta
      xet =xet +xval
      if(case.eq.1)xet1=xet1+xval
      if(case.eq.2)xet2=xet2+xval
      if(case.eq.3)xet3=xet3+xval
      if(case.eq.4)xet4=xet4+xval
  450 continue
  550 continue
  650 continue
c
  750 continue
  850 continue
      write(6,678)omap(a+no),aspn,xet
      write(4,678)omap(a+no),aspn,xet
c     write(4,679)xet1,xet2,xet3,xet4
c     write(6,679)xet1,xet2,xet3,xet4
  678 format(' done with virtual mo =',i4,' spin',i4,' energy=',f20.12)
  679 format(4f18.10)
C      call timit(4,6)
      close(unit=4,err=51)
   51 continue
      call ffile(4,'term',1)
c
      write(9,'(i4)')a
      write(9,'(f20.15)')xet
      close(unit=9,err=59)
   59 continue
      call ffile(9,'trstr',1)
      
  950 continue
c
      write(4,*)
      write(6,*)
      write(4,679)xet1,xet2,xet3,xet4
      write(6,679)xet1,xet2,xet3,xet4
c
      write(6,*)
      write(6,*)'triples contribution to ccsd(t)',xet
c     write(6,*)'triples contribution to mp4    ',xet4
      write(6,*)
c     call flush(6)
c
      return
      end
