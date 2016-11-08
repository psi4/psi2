      subroutine sdtg(t1,t2,din,ein,fin,fock,w,v,w1,v1,
     .                orbsym,flov,fzo,fzv,toff,tadd,zlx,ipq,
     .                itr,itv,eoff,eadd,foff,fadd,doff,dadd,
     .                no,nv,nirred,xet,xet4,tra,tva,
     .                j1,j2,moo,nvv,dden,eden,fden,
     .                nshov,nslov,nslvo,ndimt2,j22,ffden,
     .                nthrd)
      implicit integer(a-z)
      real*8 din(*),ein(*),fin(*),t1(no,nv),t2(*),fock(*),
     .       moo(no,no,nthrd),nvv(nv,nv,nthrd)
      integer orbsym(*),flov(nirred,4),fzo(*),fzv(*),
     .        toff(no,no,2),tadd(nv,nv),zlx(nv,nv),
     .        itr(*),itv(*),doff(*),dadd(*),eoff(*),eadd(*),
     .        foff(*),fadd(*),ipq(*),tra(no,no),tva(nv,nv)
      real*8 val,valf,vale,val1,val2,delta,xet,xet4,fac,xval,aval,bval,
     .       third,qvalw,qvalv,avalw,avalv,bvalw,bvalv,
     .       qval,dval,delta3,dvale,deli,delj,delc,deld,dvalj,dvalf,
     .       pvalw,rvalv,pvalv,dval1,rvalw,half,term,
     .       cvalw,rvalw1,rvalw2,rvalw3,val3
      real*8 w(no,no,no,nthrd),v(no,no,no,nthrd)
      real*8 dden(nshov,nthrd),eden(nslov,nthrd),fden(no,nv,nv),
     .       j1(no,nv,nthrd),j2(ndimt2,nthrd),j22(nshov)
      real*8 ffden(nslvo,nthrd),w1(nv,nv,nv,nthrd),v1(nv,nv,nv,nthrd)
      integer prcntr
c
      half=1.0d0/2.0d0
      third=1.0d0/3.0d0
c
      call zero(moo,no*no*nthrd)
      call zero(nvv,nv*nv*nthrd)
      call zero(dden,nshov*nthrd)
      call zero(eden,nslov*nthrd)
      call zero(ffden,nslvo*nthrd)
      call zero(j2,ndimt2*nthrd) 
      call zero(j1,no*nv*nthrd)
c
      write(4,*)'  '
      write(4,*)'Entering  new ccsd(t) gradient  calculation  '
      write(6,*)'  '
      write(6,*)'Entering  new ccsd(t) gradient  calculation  '
#ifdef SGI_PARALLEL
C$    write(6,*)'parallel sgi version'
C$    write(6,*)'running with ', nthrd, ' threads'
#endif
c
      xet=0.0d0
c
#ifdef SGI_PARALLEL
C$DOACROSS local(
C$&        aa abcsym,acsym,af,ai,aibf,aijm,
C$&        aikm,asym,avalv,avalw,bb,bcsym,bf,
C$&        bi,bj,bjaf,bjcf,bjim,bjkm,bsym,
C$&        bvalw,cc,cf,ci,ck,ckaf,ckbf,ckim,
C$&        ckjm,csym,cvalw,delta,dval,dvalf,
C$&        eoffia,eoffjb,eoffkc,faddia,faddjb,
C$&        faddkc,ff,ff2,fl,fm,fm2,fsym,ia,
C$&        iabf,iacf,iajb,iajm,iakm,iasym,ii,
C$&        ijaf,ijksym,ijsym,ikaf,iksym,im,
C$&        imab,imac,isym,jb,jbaf,jbcf,jbim,
C$&        jbkc,jbkm,jbsym,jibf,jj,jkbf,jksym,
C$&        jm,jmbc,jsym,kc,kcaf,kcai,kcbf,
C$&        kcia,kcim,kcjm,kcsym,kicf,kjcf,kk,
C$&        km,ksym,lf,lf2,ll,lm,lm2,miba,mica,
C$&        mjab,mjcb,mkac,mkbc,msym,mythrd,
C$&        pvalv,pvalw,qvalw,rvalv,rvalw,
C$&        rvalw1,rvalw2,rvalw3,term,val,val1,
C$&        val2,val3,vale,valf,xval,zlab,
C$&        zlac,zlaf,zlba,zlbc,zlbf,zlca,zlcb,
C$&        zlcf,absym,aa,abcsym,aicf,l,
C$&        a,b,c,i,j,k,f,m
C$&        ),
C$&        reduction(xet),
C$&        mp_schedtype=dynamic,
C$&        chunk=1
#endif
      do 3290 a=1,nv
#ifdef SGI_PARALLEL
         mythrd = mp_my_threadnum() + 1
#else
         mythrd = 1
#endif
         if(fzv(a).eq.1)go to 3290
         asym=orbsym(a+no)
         aa=ipq(a+no)+a+no
         do 3280 b=1,a 
            if(fzv(b).eq.1)go to 3280
            bsym=orbsym(b+no)
            absym=ieor(asym,bsym)
            zlba=zlx(b,a)
            zlab=zlx(a,b)
            bb=ipq(b+no)+b+no
            do 3470 c=1,b 
               if(fzv(c).eq.1)go to 3470
               zlbc=zlx(b,c)
               zlcb=zlx(c,b)
               zlac=zlx(a,c)
               zlca=zlx(c,a)
               cc=ipq(c+no)+c+no
               csym=orbsym(c+no)
               abcsym=ieor(absym,csym)
               bcsym=ieor(bsym,csym)
               acsym=ieor(asym,csym)
               call zero(w(1,1,1,mythrd),no*no*no)
               call zero(v(1,1,1,mythrd),no*no*no)
c
               do 3360 i=1,no
                  if(fzo(i).eq.1)go to 3360
                  isym=orbsym(i)
                  ia=itr(a)+i
                  iasym=ieor(isym,asym)
                  ii=ipq(i)+i
                  faddia=fadd(ia)
                  eoffia=eoff(ia)
                  do 3350 j=1,no
                     if(fzo(j).eq.1)go to 3350
                     jsym=orbsym(j)
                     ijsym=ieor(isym,jsym)
                     jbsym=ieor(jsym,bsym)
                     jb=itr(b)+j
                     jj=ipq(j)+j
                     faddjb=fadd(jb)
                     eoffjb=eoff(jb)
                     do 3340 k=1,no
                        if(fzo(k).eq.1)go to 3340
                        ksym=orbsym(k)
                        ijksym=ieor(ijsym,ksym)
                        if(ijksym.ne.abcsym)go to 3340
                        kcsym=ieor(ksym,csym)
                        kc=itr(c)+k
                        faddkc=fadd(kc)
                        eoffkc=eoff(kc)
                        kk=ipq(k)+k
c
      delta=fock(ii)+fock(jj)+fock(kk)-fock(aa)-fock(bb)-fock(cc)
c
      iajb=doff(max0(ia,jb))+dadd(min0(ia,jb))
      jbkc=doff(max0(kc,jb))+dadd(min0(kc,jb))
      kcia=doff(max0(kc,ia))+dadd(min0(kc,ia))
      val=din(iajb)*t1(k,c)+din(jbkc)*t1(i,a)+din(kcia)*t1(j,b)
      v(i,j,k,mythrd)=v(i,j,k,mythrd)+val/delta
c               
      valf=0.0d0
      vale=0.0d0
c
      fsym=ieor(bsym,iasym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.gt.b)lf2=b
      do 3351 f=ff,lf2
         bf=itv(b)+f
         iabf=foff(bf)+faddia
         zlcf=zlx(c,f)
         kjcf=toff(k,j,zlcf)+tadd(c,f)
         valf=valf+fin(iabf)*t2(kjcf)
 3351 continue
      ff2=ff
      if(ff.le.b)ff2=b+1
      do 3361 f=ff2,lf
         bf=itv(f)+b
         iabf=foff(bf)+faddia
         kjcf=toff(k,j,2)+tadd(c,f)
         valf=valf+fin(iabf)*t2(kjcf)
 3361 continue
c
      fsym=ieor(csym,jbsym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.gt.c)lf2=c
      do 3352 f=ff,lf2
         cf=itv(c)+f
         jbcf=foff(cf)+faddjb
         ikaf=toff(i,k,1)+tadd(a,f)
         valf=valf+fin(jbcf)*t2(ikaf)
 3352 continue
      ff2=ff
      if(ff.le.c)ff2=c+1
      do 3362 f=ff2,lf
         cf=itv(f)+c
         jbcf=foff(cf)+faddjb
         zlaf=zlx(a,f)
         ikaf=toff(i,k,zlaf)+tadd(a,f)
         valf=valf+fin(jbcf)*t2(ikaf)
 3362 continue
c
      fsym=ieor(asym,kcsym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.gt.a)lf2=a
      do 3353 f=ff,lf2
         af=itv(a)+f
         kcaf=foff(af)+faddkc
         zlbf=zlx(b,f)
         jibf=toff(j,i,zlbf)+tadd(b,f)
         valf=valf+fin(kcaf)*t2(jibf)
 3353 continue
      ff2=ff
      if(ff.le.a)ff2=a+1
      do 3363 f=ff2,lf
         af=itv(f)+a
         kcaf=foff(af)+faddkc
         jibf=toff(j,i,2)+tadd(b,f)
         valf=valf+fin(kcaf)*t2(jibf)
 3363 continue
c
      fsym=ieor(csym,iasym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.gt.c)lf2=c
      do 3354 f=ff,lf2
         cf=itv(c)+f
         iacf=foff(cf)+faddia
         jkbf=toff(j,k,1)+tadd(b,f)
         valf=valf+fin(iacf)*t2(jkbf)
 3354 continue
      ff2=ff
      if(ff.le.c)ff2=c+1
      do 3364 f=ff2,lf
         cf=itv(f)+c
         iacf=foff(cf)+faddia
         zlbf=zlx(b,f)
         jkbf=toff(j,k,zlbf)+tadd(b,f)
         valf=valf+fin(iacf)*t2(jkbf)
 3364 continue
c
      fsym=ieor(asym,jbsym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.gt.a)lf2=a
      do 3355 f=ff,lf2
         af=itv(a)+f
         jbaf=foff(af)+faddjb
         zlcf=zlx(c,f)
         kicf=toff(k,i,zlcf)+tadd(c,f)
         valf=valf+fin(jbaf)*t2(kicf)
 3355 continue
      ff2=ff
      if(ff.le.a)ff2=a+1
      do 3365 f=ff2,lf
         af=itv(f)+a
         jbaf=foff(af)+faddjb
         kicf=toff(k,i,2)+tadd(c,f)
         valf=valf+fin(jbaf)*t2(kicf)
 3365 continue
c
      fsym=ieor(bsym,kcsym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.gt.b)lf2=b
      do 3356 f=ff,lf2
         bf=itv(b)+f
         kcbf=foff(bf)+faddkc
         ijaf=toff(i,j,1)+tadd(a,f)
         valf=valf+fin(kcbf)*t2(ijaf)
 3356 continue
      ff2=ff
      if(ff.le.b)ff2=b+1
      do 3366 f=ff2,lf
         bf=itv(f)+b
         kcbf=foff(bf)+faddkc
         zlaf=zlx(a,f)
         ijaf=toff(i,j,zlaf)+tadd(a,f)
         valf=valf+fin(kcbf)*t2(ijaf)
 3366 continue
c
      msym=ieor(jsym,iasym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      lm2=lm
      if(lm.gt.j)lm2=j
      do 3251 m=fm,lm2
         jm=itr(j)+m
         iajm=eoffia+eadd(jm)
         mkbc=toff(m,k,zlbc)+tadd(b,c)
         vale=vale+ein(iajm)*t2(mkbc)
 3251 continue
      fm2=fm
      if(fm.le.j)fm2=j+1
      do 3261 m=fm2,lm
         jm=itr(m)+j
         iajm=eoffia+eadd(jm)
         mkbc=toff(m,k,zlbc)+tadd(b,c)
         vale=vale+ein(iajm)*t2(mkbc)
 3261 continue
c
      msym=ieor(ksym,jbsym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      lm2=lm
      if(lm.gt.k)lm2=k
      do 3252 m=fm,lm2
         km=itr(k)+m
         jbkm=eoffjb+eadd(km)
         mica=toff(m,i,zlca)+tadd(c,a)
         vale=vale+ein(jbkm)*t2(mica)
 3252 continue
      fm2=fm
      if(fm.le.k)fm2=k+1
      do 3262 m=fm2,lm
         km=itr(m)+k
         jbkm=eoffjb+eadd(km)
         mica=toff(m,i,zlca)+tadd(c,a)
         vale=vale+ein(jbkm)*t2(mica)
 3262 continue
c
      msym=ieor(isym,kcsym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      lm2=lm
      if(lm.gt.i)lm2=i
      do 3253 m=fm,lm2
         im=itr(i)+m
c        im=tra(i,m)
         kcim=eoffkc+eadd(im)
         mjab=toff(m,j,zlab)+tadd(a,b)
         vale=vale+ein(kcim)*t2(mjab)
 3253 continue
      fm2=fm
      if(fm.le.i)fm2=i+1
      do 3263 m=fm2,lm
         im=itr(m)+i
c        im=tra(i,m)
         kcim=eoffkc+eadd(im)
         mjab=toff(m,j,zlab)+tadd(a,b)
         vale=vale+ein(kcim)*t2(mjab)
 3263 continue
c
      msym=ieor(ksym,iasym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      lm2=lm
      if(lm.gt.k)lm2=k
      do 3254 m=fm,lm2
         km=itr(k)+m
         iakm=eoffia+eadd(km)
         mjcb=toff(m,j,zlcb)+tadd(c,b)
         vale=vale+ein(iakm)*t2(mjcb)
 3254 continue
      fm2=fm
      if(fm.le.k)fm2=k+1
      do 3264 m=fm2,lm
         km=itr(m)+k
         iakm=eoffia+eadd(km)
         mjcb=toff(m,j,zlcb)+tadd(c,b)
         vale=vale+ein(iakm)*t2(mjcb)
 3264 continue
c
      msym=ieor(isym,jbsym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      lm2=lm
      if(lm.gt.i)lm2=i
      do 3255 m=fm,lm2
         im=itr(i)+m
         jbim=eoffjb+eadd(im)
         mkac=toff(m,k,zlac)+tadd(a,c)
         vale=vale+ein(jbim)*t2(mkac)
 3255 continue
      fm2=fm
      if(fm.le.i)fm2=i+1
      do 3265 m=fm2,lm
         im=itr(m)+i
         jbim=eoffjb+eadd(im)
         mkac=toff(m,k,zlac)+tadd(a,c)
         vale=vale+ein(jbim)*t2(mkac)
 3265 continue
c
      msym=ieor(jsym,kcsym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      lm2=lm
      if(lm.gt.j)lm2=j
      do 3256 m=fm,lm2
         jm=itr(j)+m
         kcjm=eoffkc+eadd(jm)
         miba=toff(m,i,zlba)+tadd(b,a)
         vale=vale+ein(kcjm)*t2(miba)
 3256 continue
      fm2=fm
      if(fm.le.j)fm2=j+1
      do 3266 m=fm2,lm
         jm=itr(m)+j
         kcjm=eoffkc+eadd(jm)
         miba=toff(m,i,zlba)+tadd(b,a)
         vale=vale+ein(kcjm)*t2(miba)
 3266 continue
c
      w(i,j,k,mythrd)=w(i,j,k,mythrd)+(valf-vale)/delta
 3340 continue
 3350 continue
 3360 continue
c3370 continue
c
c           do 3470 c=1,b 
c              if(fzv(c).eq.1)go to 3470
c              csym=orbsym(c+no)
c              bcsym=ieor(bsym,csym)
c              acsym=ieor(asym,csym)
c              abcsym=ieor(absym,csym)
c              cc=ipq(c+no)+c+no
c              zlbc=zlx(b,c)
               do 3460 i=1,no
                  if(fzo(i).eq.1)go to 3460
                  isym=orbsym(i)
                  ia=itr(a)+i
                  ai=ia
                  bi=itr(b)+i
                  ci=itr(c)+i
                  iasym=ieor(isym,asym)
                  ii=ipq(i)+i
                  do 3450 j=1,no
                     if(fzo(j).eq.1)go to 3450
                     jsym=orbsym(j)
                     ijsym=ieor(isym,jsym)
                     jbsym=ieor(jsym,bsym)
                     jb=itr(b)+j
                     bj=jb
                     jj=ipq(j)+j
                     do 3440 k=1,no
                        if(fzo(k).eq.1)go to 3440
                        ksym=orbsym(k)
                        ijksym=ieor(ijsym,ksym)
                        iksym=ieor(isym,ksym)
                        jksym=ieor(jsym,ksym)
                        if(ijksym.ne.abcsym)go to 3440
                        kcsym=ieor(ksym,csym)
                        kc=itr(c)+k
                        ck=kc
                        kk=ipq(k)+k
      delta=fock(ii)+fock(jj)+fock(kk)-fock(aa)-fock(bb)-fock(cc)
c
      avalw=4.d0*w(i,j,k,mythrd)+  w(k,i,j,mythrd)+  w(j,k,i,mythrd) 
      bvalw= w(i,j,k,mythrd)+4.0d0*w(k,i,j,mythrd)+  w(j,k,i,mythrd) 
      cvalw= w(i,j,k,mythrd)+  w(k,i,j,mythrd)+4.0d0*w(j,k,i,mythrd) 
      pvalw=2.d0*(w(k,j,i,mythrd)+w(i,k,j,mythrd)+w(j,i,k,mythrd))
c     qvalw=avalw-bvalw
      qvalw=avalw-pvalw
      rvalw=avalw-pvalw
      rvalw1=avalw-pvalw
      rvalw2=bvalw-pvalw
      rvalw3=cvalw-pvalw
c
      avalv=4.d0*v(i,j,k,mythrd)+v(k,i,j,mythrd)+v(j,k,i,mythrd) 
      pvalv=2.d0*(v(k,j,i,mythrd)+v(i,k,j,mythrd)+v(j,i,k,mythrd))
      rvalv=avalv-pvalv
c
      if(a.eq.b.or.b.eq.c)then
         qvalw=qvalw*half 
         rvalw1=rvalw1*half 
         rvalw2=rvalw2*half 
         rvalw3=rvalw3*half 
         rvalv=rvalv*half
      endif
c
      val1=w(i,j,k,mythrd)+v(i,j,k,mythrd)
      xval=val1*qvalw*delta
      xet =xet +xval
c
      dval=rvalw1+rvalw1+rvalv
      dval=dval+dval
c
      if(ksym.eq.csym)then
      iajb=doff(max0(ia,jb))+dadd(min0(ia,jb))
      dden(iajb,mythrd)=dden(iajb,mythrd)+t1(k,c)*rvalw1
      j1(k,c,mythrd)=j1(k,c,mythrd)+din(iajb)*rvalw1
      endif
c
      if(isym.eq.asym)then
      jbkc=doff(max0(jb,kc))+dadd(min0(jb,kc))
      dden(jbkc,mythrd)=dden(jbkc,mythrd)+t1(i,a)*rvalw1
      j1(i,a,mythrd)=j1(i,a,mythrd)+din(jbkc)*rvalw1
      endif
c
      if(jsym.eq.bsym)then
      kcai=doff(max0(kc,ai))+dadd(min0(kc,ai))
      dden(kcai,mythrd)=dden(kcai,mythrd)+t1(j,b)*rvalw1
      j1(j,b,mythrd)=j1(j,b,mythrd)+din(kcai)*rvalw1
      endif
c
      msym=ieor(bcsym,ksym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      do 3410 m=fm,lm
         mkbc=toff(m,k,1)+tadd(b,c)
         jm=itr(max0(j,m))+min0(j,m)
         aijm=eoff(ai)+eadd(jm)
         eden(aijm,mythrd)=eden(aijm,mythrd)-t2(mkbc)*dval
 3410 continue
      msym=ieor(bcsym,jsym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      do 3411 m=fm,lm
         jmbc=toff(j,m,1)+tadd(b,c)
         km=itr(max0(k,m))+min0(k,m)
         aikm=eoff(ai)+eadd(km)
         eden(aikm,mythrd)=eden(aikm,mythrd)-t2(jmbc)*dval
 3411 continue
      msym=ieor(acsym,ksym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      do 3412 m=fm,lm
         mkac=toff(m,k,1)+tadd(a,c)
         im=itr(max0(i,m))+min0(i,m)
         bjim=eoff(bj)+eadd(im)
         eden(bjim,mythrd)=eden(bjim,mythrd)-t2(mkac)*dval
 3412 continue
      msym=ieor(acsym,isym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      do 3413 m=fm,lm
         imac=toff(i,m,1)+tadd(a,c)
         km=itr(max0(k,m))+min0(k,m)
         bjkm=eoff(bj)+eadd(km)
         eden(bjkm,mythrd)=eden(bjkm,mythrd)-t2(imac)*dval
 3413 continue
      msym=ieor(absym,isym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      do 3414 m=fm,lm
         imab=toff(i,m,1)+tadd(a,b)
         jm=itr(max0(j,m))+min0(j,m)
         ckjm=eoff(ck)+eadd(jm)
         eden(ckjm,mythrd)=eden(ckjm,mythrd)-t2(imab)*dval
 3414 continue
      msym=ieor(absym,jsym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      do 3415 m=fm,lm
         mjab=toff(m,j,1)+tadd(a,b)
         im=itr(max0(i,m))+min0(i,m)
         ckim=eoff(ck)+eadd(im)
         eden(ckim,mythrd)=eden(ckim,mythrd)-t2(mjab)*dval
 3415 continue
c
      msym=ieor(bcsym,ksym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      if(b.ge.c)then
      do 3421 m=fm,lm
         mkbc=toff(m,k,1)+tadd(b,c)
         jm=itr(max0(j,m))+min0(j,m)
         aijm=eoff(ai)+eadd(jm)
         j2(mkbc,mythrd)=j2(mkbc,mythrd)-ein(aijm)*dval
 3421 continue
      endif
      if(b.le.c)then
      do 3422 m=fm,lm
         mkbc=toff(m,k,2)+tadd(b,c)
         jm=itr(max0(j,m))+min0(j,m)
         aijm=eoff(ai)+eadd(jm)
         j2(mkbc,mythrd)=j2(mkbc,mythrd)-ein(aijm)*dval
 3422 continue
      endif
c
      msym=ieor(acsym,isym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      if(a.ge.c)then
      do 3423 m=fm,lm
         imac=toff(i,m,1)+tadd(a,c)
         km=itr(max0(k,m))+min0(k,m)
         bjkm=eoff(bj)+eadd(km)
         j2(imac,mythrd)=j2(imac,mythrd)-ein(bjkm)*dval
 3423 continue
      endif
      if(a.le.c)then
      do 3424 m=fm,lm
         imac=toff(i,m,2)+tadd(a,c)
         km=itr(max0(k,m))+min0(k,m)
         bjkm=eoff(bj)+eadd(km)
         j2(imac,mythrd)=j2(imac,mythrd)-ein(bjkm)*dval
 3424 continue
      endif
c
      msym=ieor(absym,jsym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      if(a.ge.b)then
      do 3425 m=fm,lm
         mjab=toff(m,j,1)+tadd(a,b)
         im=itr(max0(i,m))+min0(i,m)
         ckim=eoff(ck)+eadd(im)
         j2(mjab,mythrd)=j2(mjab,mythrd)-ein(ckim)*dval
 3425 continue
      endif
      if(a.le.b)then
      do 3426 m=fm,lm
         mjab=toff(m,j,2)+tadd(a,b)
         im=itr(max0(i,m))+min0(i,m)
         ckim=eoff(ck)+eadd(im)
         j2(mjab,mythrd)=j2(mjab,mythrd)-ein(ckim)*dval
 3426 continue
      endif
c
      msym=ieor(bcsym,jsym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      if(b.ge.c)then
      do 3441 m=fm,lm
         jmbc=toff(j,m,1)+tadd(b,c)
         km=itr(max0(k,m))+min0(k,m)
         aikm=eoff(ai)+eadd(km)
         j2(jmbc,mythrd)=j2(jmbc,mythrd)-ein(aikm)*dval
 3441 continue
      endif
      if(b.le.c)then
      do 3442 m=fm,lm
         jmbc=toff(j,m,2)+tadd(b,c)
         km=itr(max0(k,m))+min0(k,m)
         aikm=eoff(ai)+eadd(km)
         j2(jmbc,mythrd)=j2(jmbc,mythrd)-ein(aikm)*dval
 3442 continue
      endif
c
      msym=ieor(acsym,ksym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      if(a.ge.c)then
      do 3443 m=fm,lm
         mkac=toff(m,k,1)+tadd(a,c)
         im=itr(max0(i,m))+min0(i,m)
         bjim=eoff(bj)+eadd(im)
         j2(mkac,mythrd)=j2(mkac,mythrd)-ein(bjim)*dval
 3443 continue
      endif
      if(a.le.c)then
      do 3444 m=fm,lm
         mkac=toff(m,k,2)+tadd(a,c)
         im=itr(max0(i,m))+min0(i,m)
         bjim=eoff(bj)+eadd(im)
         j2(mkac,mythrd)=j2(mkac,mythrd)-ein(bjim)*dval
 3444 continue
      endif
c
      msym=ieor(absym,isym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      if(a.ge.b)then
      do 3445 m=fm,lm
         imab=toff(i,m,1)+tadd(a,b)
         jm=itr(max0(j,m))+min0(j,m)
         ckjm=eoff(ck)+eadd(jm)
         j2(imab,mythrd)=j2(imab,mythrd)-ein(ckjm)*dval
 3445 continue
      endif
      if(a.le.b)then
      do 3446 m=fm,lm
         imab=toff(i,m,2)+tadd(a,b)
         jm=itr(max0(j,m))+min0(j,m)
         ckjm=eoff(ck)+eadd(jm)
         j2(imab,mythrd)=j2(imab,mythrd)-ein(ckjm)*dval
 3446 continue
      endif
c
c
      fsym=ieor(jksym,csym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.ge.c)lf2=c
      do 3721 f=ff,lf2
         kjcf=toff(k,j,1)+tadd(c,f)
         bf=itv(max0(b,f))+min0(b,f)
         aibf=foff(bf)+fadd(ai)
         j2(kjcf,mythrd)=j2(kjcf,mythrd)+fin(aibf)*dval
 3721 continue
      ff2=ff
      if(ff.le.c)ff2=c
      do 3722 f=ff2,lf
         kjcf=toff(k,j,2)+tadd(c,f)
         bf=itv(max0(b,f))+min0(b,f)
         aibf=foff(bf)+fadd(ai)
         j2(kjcf,mythrd)=j2(kjcf,mythrd)+fin(aibf)*dval
 3722 continue
c
      fsym=ieor(iksym,asym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.ge.a)lf2=a
      do 3723 f=ff,lf2
         ikaf=toff(i,k,1)+tadd(a,f)
         cf=itv(max0(c,f))+min0(c,f)
         bjcf=foff(cf)+fadd(bj)
         j2(ikaf,mythrd)=j2(ikaf,mythrd)+fin(bjcf)*dval
 3723 continue
      ff2=ff
      if(ff.le.a)ff2=a
      do 3724 f=ff2,lf
         ikaf=toff(i,k,2)+tadd(a,f)
         cf=itv(max0(c,f))+min0(c,f)
         bjcf=foff(cf)+fadd(bj)
         j2(ikaf,mythrd)=j2(ikaf,mythrd)+fin(bjcf)*dval
 3724 continue
c
      fsym=ieor(ijsym,bsym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.ge.b)lf2=b
      do 3725 f=ff,lf2
         jibf=toff(j,i,1)+tadd(b,f)
         af=itv(max0(a,f))+min0(a,f)
         ckaf=foff(af)+fadd(ck)
         j2(jibf,mythrd)=j2(jibf,mythrd)+fin(ckaf)*dval
 3725 continue
      ff2=ff
      if(ff.le.b)ff2=b
      do 3726 f=ff2,lf
         jibf=toff(j,i,2)+tadd(b,f)
         af=itv(max0(a,f))+min0(a,f)
         ckaf=foff(af)+fadd(ck)
         j2(jibf,mythrd)=j2(jibf,mythrd)+fin(ckaf)*dval
 3726 continue
c
c
      fsym=ieor(iksym,csym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.ge.c)lf2=c
      do 3741 f=ff,lf2
         kicf=toff(k,i,1)+tadd(c,f)
         af=itv(max0(a,f))+min0(a,f)
         bjaf=foff(af)+fadd(bj)
         j2(kicf,mythrd)=j2(kicf,mythrd)+fin(bjaf)*dval
 3741 continue
      ff2=ff
      if(ff.le.c)ff2=c
      do 3742 f=ff2,lf
         kicf=toff(k,i,2)+tadd(c,f)
         af=itv(max0(a,f))+min0(a,f)
         bjaf=foff(af)+fadd(bj)
         j2(kicf,mythrd)=j2(kicf,mythrd)+fin(bjaf)*dval
 3742 continue
c
      fsym=ieor(ijsym,asym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.ge.a)lf2=a
      do 3743 f=ff,lf2
         ijaf=toff(i,j,1)+tadd(a,f)
         bf=itv(max0(b,f))+min0(b,f)
         ckbf=foff(bf)+fadd(ck)
         j2(ijaf,mythrd)=j2(ijaf,mythrd)+fin(ckbf)*dval
 3743 continue
      ff2=ff
      if(ff.le.a)ff2=a
      do 3744 f=ff2,lf
         ijaf=toff(i,j,2)+tadd(a,f)
         bf=itv(max0(b,f))+min0(b,f)
         ckbf=foff(bf)+fadd(ck)
         j2(ijaf,mythrd)=j2(ijaf,mythrd)+fin(ckbf)*dval
 3744 continue
c
      fsym=ieor(jksym,bsym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.ge.b)lf2=b
      do 3745 f=ff,lf2
         jkbf=toff(j,k,1)+tadd(b,f)
         cf=itv(max0(c,f))+min0(c,f)
         aicf=foff(cf)+fadd(ai)
         j2(jkbf,mythrd)=j2(jkbf,mythrd)+fin(aicf)*dval
 3745 continue
      ff2=ff
      if(ff.le.b)ff2=b
      do 3746 f=ff2,lf
         jkbf=toff(j,k,2)+tadd(b,f)
         cf=itv(max0(c,f))+min0(c,f)
         aicf=foff(cf)+fadd(ai)
         j2(jkbf,mythrd)=j2(jkbf,mythrd)+fin(aicf)*dval
 3746 continue
c
      fsym=ieor(jksym,csym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      do 3431 f=ff,lf
         zlcf=zlx(c,f)
         kjcf=toff(k,j,zlcf)+tadd(c,f)
         dvalf=dval*t2(kjcf)
         bf=itv(max0(b,f))+min0(b,f)
         aibf=foff(bf)+fadd(ai)
         ffden(aibf,mythrd)=ffden(aibf,mythrd)+dvalf
 3431 continue
      fsym=ieor(jksym,bsym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      do 3432 f=ff,lf
         zlbf=zlx(b,f)
         jkbf=toff(j,k,zlbf)+tadd(b,f)
         dvalf=dval*t2(jkbf)
         cf=itv(max0(c,f))+min0(c,f)
         aicf=foff(cf)+fadd(ai)
         ffden(aicf,mythrd)=ffden(aicf,mythrd)+dvalf
 3432 continue
      fsym=ieor(iksym,asym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      do 3433 f=ff,lf
         zlaf=zlx(a,f)
         ikaf=toff(i,k,zlaf)+tadd(a,f)
         dvalf=dval*t2(ikaf)
         cf=itv(max0(c,f))+min0(c,f)
         jbcf=foff(cf)+fadd(jb)
         ffden(jbcf,mythrd)=ffden(jbcf,mythrd)+dvalf
 3433 continue
      fsym=ieor(iksym,csym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      do 3434 f=ff,lf
         zlcf=zlx(c,f)
         kicf=toff(k,i,zlcf)+tadd(c,f)
         dvalf=dval*t2(kicf)
         af=itv(max0(a,f))+min0(a,f)
         jbaf=foff(af)+fadd(jb)
         ffden(jbaf,mythrd)=ffden(jbaf,mythrd)+dvalf
 3434 continue
      fsym=ieor(ijsym,bsym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      do 3435 f=ff,lf
         zlbf=zlx(b,f)
         jibf=toff(j,i,zlbf)+tadd(b,f)
         dvalf=dval*t2(jibf)
         af=itv(max0(a,f))+min0(a,f)
         kcaf=foff(af)+fadd(kc)
         ffden(kcaf,mythrd)=ffden(kcaf,mythrd)+dvalf
 3435 continue
      fsym=ieor(ijsym,asym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      do 3436 f=ff,lf
         zlaf=zlx(a,f)
         ijaf=toff(i,j,zlaf)+tadd(a,f)
         dvalf=dval*t2(ijaf)
         bf=itv(max0(b,f))+min0(b,f)
         kcbf=foff(bf)+fadd(kc)
         ffden(kcbf,mythrd)=ffden(kcbf,mythrd)+dvalf
 3436 continue
c
c >>> moo
c
      fl=flov(ksym+1,1)
      ll=flov(ksym+1,2)
      do 3530 l=fl,ll
      val1=w(i,j,l,mythrd)+v(i,j,l,mythrd)
      val3=w(j,l,i,mythrd)+v(j,l,i,mythrd)
      val2=w(l,i,j,mythrd)+v(l,i,j,mythrd)
      term=+val1*rvalw1+val2*rvalw2+val3*rvalw3
      moo(k,l,mythrd)=moo(k,l,mythrd)+term
c     moo(l,k)=moo(l,k)+term
 3530 continue
c
 3440 continue
 3450 continue
 3460 continue
 3470 continue
c
 3280 continue
#ifndef SGI_PARALLEL
c     write(6,*)' done with virtual mo =',a
      write(4,*)' done with virtual mo =',a
C      if (prcntr('IS_ON BRIEF').eq.0) call timit(4,6)
      close(unit=4,err=351)
  351 continue
      call ffile(4,'term',1)
c     write(6,'(8f10.6)') fden
#endif
 3290 continue
c
c sum the contributions from each thread
c
      do 4590 ith=2,nthrd
         do 4580 j=1,nslvo
            ffden(j,1) = ffden(j,1) + ffden(j,ith)
 4580    continue
 4590 continue
      do 4390 ith=2,nthrd
         do 4380 j=1,nslov
            eden(j,1) = eden(j,1) + eden(j,ith)
 4380    continue
 4390 continue
      do 4690 ith=2,nthrd
         do 4680 i=1,nv
            do 4670 j=1,no
               j1(j,i,1) = j1(j,i,1) + j1(j,i,ith)
 4670       continue
 4680    continue
 4690 continue
      do 4790 ith=2,nthrd
         do 4780 j=1,ndimt2
            j2(j,1) = j2(j,1) + j2(j,ith)
 4780    continue
 4790 continue
c
c    symmetrize fden and write to file74 in appropiate fromat for ccqg
c
      write(6,*)'new fden'
      itap74=74
      call rfile(itap74)
      i74=1
      do 3890 a=1,nv
         if(fzv(a).eq.1)go to 3890
         asym=orbsym(a+no)
         call zero(fden,no*nv*nv)
      do 3880 i=1,no
         if(fzo(i).eq.1)go to 3880
         isym=orbsym(i)
         aisym=ieor(asym,isym)
         ai=itr(a)+i
      do 3870 b=1,nv
         if(fzv(b).eq.1)go to 3870
         bsym=orbsym(b+no)
      do 3860 c=1,nv
         if(fzv(c).eq.1)go to 3860
         csym=orbsym(c+no)
         bcsym=ieor(bsym,csym)
         if(aisym.ne.bcsym)go to 3860
         bc=itv(max0(b,c))+min0(b,c)
         aibc=foff(bc)+fadd(ai)
         val=ffden(aibc,1)
c        write(6,'(i6,f10.6)') aibc,val
         if(b.eq.c)val=val+val
         fden(i,b,c)=val
         fden(i,c,b)=val
 3860 continue
 3870 continue
 3880 continue
      call wwritw(itap74,fden,intowp(no*nv*nv),i74,i74)
 3890 continue
      call rclose(itap74,3)
c
      xet=xet*2.0d0
c
      write(6,*)
      write(6,*)'triples contribution to ccsd(t)',xet
      write(6,*)
      call fioflu(6)
c
c >>> avoid symmetryzing eden inside n**7 loop !
c
      do 3070 a=1,nv
         if(fzv(a).eq.1)go to 3070
         asym=orbsym(a+no)
         do 3050 i=1,no
            if(fzo(i).eq.1)go to 3050
            isym=orbsym(i)
            if(isym.ne.asym)go to 3050
            ai=itr(a)+i
            do 3030 j=1,no
               if(fzo(j).eq.1)go to 3030
               jj=itr(j)+j
               aijj=eoff(ai)+eadd(jj)
               eden(aijj,1)= (eden(aijj,1)+eden(aijj,1))
 3030       continue
 3050    continue
 3070 continue
c
c
c >>> write file70 containing j1 & j2 (read by initz.f in cczv)
c
      do 1997 i=1,no
      do 1996 a=1,nv
         j1(i,a,1)=j1(i,a,1)*2.0d0
 1996 continue
 1997 continue
c
c     write(6,*)' new j1 length is',no*nv
c     write(6,'(8f10.6)')  j1 
c     write(6,*)' new j2 length is',ndimt2
c     write(6,'(8f10.6)')  j2 
      call fioflu(6)
c
c     write(6,*)'now writing file70'
      itap70=70
      call rfile(itap70)
      call wwritw(itap70,j1,intowp(no*nv),1,i70)
      call wwritw(itap70,j2,intowp(ndimt2),i70,i70)
      call rclose(itap70,3)
c
c >>> form triples again !!!
c
      xet=0.0d0
c
#ifdef SGI_PARALLEL
C$DOACROSS local(
C$&        mythrd,isym,ii,jsym,jj,ijsym,ksym,kk,
C$&        ijksym,asym,ia,iasym,aa,faddia,eoffia,
C$&        bsym,absym,zlba,zlab,jbsym,jb,faddjb,
C$&        eoffjb,csym,zlbc,zlcb,zlac,zlca,
C$&        abcsym,af,ai,avalw,
C$&        bb,bcsym,bf,bvalw,cc,cf,cvalw,
C$&        delta,eoffkc,faddkc,fd,ff,ff2,fm,fm2,
C$&        fsym,iabf,iacf,iajb,iajm,iakm,
C$&        ijaf,ikaf,im,jbaf,jbcf,jbim,jbkc,
C$&        jbkm,jibf,jkbf,jm,kc,kcaf,kcbf,
C$&        kcia,kcim,kcjm,kcsym,kicf,kjcf,km,ld,lf,
C$&        lf2,lm,lm2,miba,mica,mjab,mjcb,mkac,
C$&        mkbc,msym,pvalw,qvalw,rvalw1,rvalw2,
C$&        rvalw3,term,val,val1,val2,val3,vale,
C$&        valf,xval,zlaf,zlbf,zlcf,d,
C$&        i,j,k,a,b,c,f,m
C$&        ),
C$&        reduction(xet),
C$&        mp_schedtype=dynamic,
C$&        chunk=1
#endif
      do 5290 i=1,no
#ifdef SGI_PARALLEL
         mythrd = mp_my_threadnum() + 1
#else
         mythrd = 1
#endif
         if(fzo(i).eq.1)go to 5290
         isym=orbsym(i)
         ii=ipq(i)+i
         do 5280 j=1,i 
            if(fzo(j).eq.1)go to 5280
            jsym=orbsym(j)
            jj=ipq(j)+j
            ijsym=ieor(isym,jsym)
            do 5270 k=1,j
               if(fzo(k).eq.1)go to 5270
               ksym=orbsym(k)
               kk=ipq(k)+k
               ijksym=ieor(ijsym,ksym)
               call zero(w1(1,1,1,mythrd),nv*nv*nv)
               call zero(v1(1,1,1,mythrd),nv*nv*nv)
c
               do 5360 a=1,nv
                  if(fzv(a).eq.1)go to 5360
                  asym=orbsym(a+no)
                  ia=itr(a)+i
                  iasym=ieor(isym,asym)
                  aa=ipq(a+no)+a+no
                  faddia=fadd(ia)
                  eoffia=eoff(ia)
                  do 5350 b=1,nv
                     if(fzv(b).eq.1)go to 5350
                     bsym=orbsym(b+no)
                     absym=ieor(asym,bsym)
                     zlba=zlx(b,a)
                     zlab=zlx(a,b)
                     bb=ipq(b+no)+b+no
                     jbsym=ieor(jsym,bsym)
                     jb=itr(b)+j
                     faddjb=fadd(jb)
                     eoffjb=eoff(jb)
                     do 5340 c=1,nv
                        if(fzv(c).eq.1)go to 5340
                        csym=orbsym(c+no)
                        zlbc=zlx(b,c)
                        zlcb=zlx(c,b)
                        zlac=zlx(a,c)
                        zlca=zlx(c,a)
                        cc=ipq(c+no)+c+no
                        abcsym=ieor(absym,csym)
                        if(ijksym.ne.abcsym)go to 5340
                        kcsym=ieor(ksym,csym)
                        kc=itr(c)+k
                        faddkc=fadd(kc)
                        eoffkc=eoff(kc)
c
      delta=fock(ii)+fock(jj)+fock(kk)-fock(aa)-fock(bb)-fock(cc)
c
      iajb=doff(max0(ia,jb))+dadd(min0(ia,jb))
      jbkc=doff(max0(kc,jb))+dadd(min0(kc,jb))
      kcia=doff(max0(kc,ia))+dadd(min0(kc,ia))
      val=din(iajb)*t1(k,c)+din(jbkc)*t1(i,a)+din(kcia)*t1(j,b)
      v1(a,b,c,mythrd)=v1(a,b,c,mythrd)+val/delta
c               
      valf=0.0d0
      vale=0.0d0
c
      fsym=ieor(bsym,iasym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.gt.b)lf2=b
      do 5351 f=ff,lf2
         bf=itv(b)+f
         iabf=foff(bf)+faddia
         zlcf=zlx(c,f)
         kjcf=toff(k,j,zlcf)+tadd(c,f)
         valf=valf+fin(iabf)*t2(kjcf)
 5351 continue
      ff2=ff
      if(ff.le.b)ff2=b+1
      do 5361 f=ff2,lf
         bf=itv(f)+b
         iabf=foff(bf)+faddia
         zlcf=zlx(c,f)
         kjcf=toff(k,j,zlcf)+tadd(c,f)
         valf=valf+fin(iabf)*t2(kjcf)
 5361 continue
c
      fsym=ieor(csym,jbsym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.gt.c)lf2=c
      do 5352 f=ff,lf2
         cf=itv(c)+f
         jbcf=foff(cf)+faddjb
         zlaf=zlx(a,f)
         ikaf=toff(i,k,zlaf)+tadd(a,f)
         valf=valf+fin(jbcf)*t2(ikaf)
 5352 continue
      ff2=ff
      if(ff.le.c)ff2=c+1
      do 5362 f=ff2,lf
         cf=itv(f)+c
         jbcf=foff(cf)+faddjb
         zlaf=zlx(a,f)
         ikaf=toff(i,k,zlaf)+tadd(a,f)
         valf=valf+fin(jbcf)*t2(ikaf)
 5362 continue
c
      fsym=ieor(asym,kcsym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.gt.a)lf2=a
      do 5353 f=ff,lf2
         af=itv(a)+f
         kcaf=foff(af)+faddkc
         zlbf=zlx(b,f)
         jibf=toff(j,i,zlbf)+tadd(b,f)
         valf=valf+fin(kcaf)*t2(jibf)
 5353 continue
      ff2=ff
      if(ff.le.a)ff2=a+1
      do 5363 f=ff2,lf
         af=itv(f)+a
         kcaf=foff(af)+faddkc
         zlbf=zlx(b,f)
         jibf=toff(j,i,zlbf)+tadd(b,f)
         valf=valf+fin(kcaf)*t2(jibf)
 5363 continue
c
      fsym=ieor(csym,iasym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.gt.c)lf2=c
      do 5354 f=ff,lf2
         cf=itv(c)+f
         iacf=foff(cf)+faddia
         zlbf=zlx(b,f)
         jkbf=toff(j,k,zlbf)+tadd(b,f)
         valf=valf+fin(iacf)*t2(jkbf)
 5354 continue
      ff2=ff
      if(ff.le.c)ff2=c+1
      do 5364 f=ff2,lf
         cf=itv(f)+c
         iacf=foff(cf)+faddia
         zlbf=zlx(b,f)
         jkbf=toff(j,k,zlbf)+tadd(b,f)
         valf=valf+fin(iacf)*t2(jkbf)
 5364 continue
c
      fsym=ieor(asym,jbsym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.gt.a)lf2=a
      do 5355 f=ff,lf2
         af=itv(a)+f
         jbaf=foff(af)+faddjb
         zlcf=zlx(c,f)
         kicf=toff(k,i,zlcf)+tadd(c,f)
         valf=valf+fin(jbaf)*t2(kicf)
 5355 continue
      ff2=ff
      if(ff.le.a)ff2=a+1
      do 5365 f=ff2,lf
         af=itv(f)+a
         jbaf=foff(af)+faddjb
         zlcf=zlx(c,f)
         kicf=toff(k,i,zlcf)+tadd(c,f)
         valf=valf+fin(jbaf)*t2(kicf)
 5365 continue
c
      fsym=ieor(bsym,kcsym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.gt.b)lf2=b
      do 5356 f=ff,lf2
         bf=itv(b)+f
         kcbf=foff(bf)+faddkc
         zlaf=zlx(a,f)
         ijaf=toff(i,j,zlaf)+tadd(a,f)
         valf=valf+fin(kcbf)*t2(ijaf)
 5356 continue
      ff2=ff
      if(ff.le.b)ff2=b+1
      do 5366 f=ff2,lf
         bf=itv(f)+b
         kcbf=foff(bf)+faddkc
         zlaf=zlx(a,f)
         ijaf=toff(i,j,zlaf)+tadd(a,f)
         valf=valf+fin(kcbf)*t2(ijaf)
 5366 continue
c
      msym=ieor(jsym,iasym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      lm2=lm
      if(lm.gt.j)lm2=j
      do 5251 m=fm,lm2
         jm=itr(j)+m
         iajm=eoffia+eadd(jm)
         mkbc=toff(m,k,zlbc)+tadd(b,c)
         vale=vale+ein(iajm)*t2(mkbc)
 5251 continue
      fm2=fm
      if(fm.le.j)fm2=j+1
      do 5261 m=fm2,lm
         jm=itr(m)+j
         iajm=eoffia+eadd(jm)
         mkbc=toff(m,k,zlbc)+tadd(b,c)
         vale=vale+ein(iajm)*t2(mkbc)
 5261 continue
c
      msym=ieor(ksym,jbsym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      lm2=lm
      if(lm.gt.k)lm2=k
      do 5252 m=fm,lm2
         km=itr(k)+m
         jbkm=eoffjb+eadd(km)
         mica=toff(m,i,zlca)+tadd(c,a)
         vale=vale+ein(jbkm)*t2(mica)
 5252 continue
      fm2=fm
      if(fm.le.k)fm2=k+1
      do 5262 m=fm2,lm
         km=itr(m)+k
         jbkm=eoffjb+eadd(km)
         mica=toff(m,i,zlca)+tadd(c,a)
         vale=vale+ein(jbkm)*t2(mica)
 5262 continue
c
      msym=ieor(isym,kcsym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      lm2=lm
      if(lm.gt.i)lm2=i
      do 5253 m=fm,lm2
         im=itr(i)+m
c        im=tra(i,m)
         kcim=eoffkc+eadd(im)
         mjab=toff(m,j,zlab)+tadd(a,b)
         vale=vale+ein(kcim)*t2(mjab)
 5253 continue
      fm2=fm
      if(fm.le.i)fm2=i+1
      do 5263 m=fm2,lm
         im=itr(m)+i
c        im=tra(i,m)
         kcim=eoffkc+eadd(im)
         mjab=toff(m,j,zlab)+tadd(a,b)
         vale=vale+ein(kcim)*t2(mjab)
 5263 continue
c
      msym=ieor(ksym,iasym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      lm2=lm
      if(lm.gt.k)lm2=k
      do 5254 m=fm,lm2
         km=itr(k)+m
         iakm=eoffia+eadd(km)
         mjcb=toff(m,j,zlcb)+tadd(c,b)
         vale=vale+ein(iakm)*t2(mjcb)
 5254 continue
      fm2=fm
      if(fm.le.k)fm2=k+1
      do 5264 m=fm2,lm
         km=itr(m)+k
         iakm=eoffia+eadd(km)
         mjcb=toff(m,j,zlcb)+tadd(c,b)
         vale=vale+ein(iakm)*t2(mjcb)
 5264 continue
c
      msym=ieor(isym,jbsym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      lm2=lm
      if(lm.gt.i)lm2=i
      do 5255 m=fm,lm2
         im=itr(i)+m
         jbim=eoffjb+eadd(im)
         mkac=toff(m,k,zlac)+tadd(a,c)
         vale=vale+ein(jbim)*t2(mkac)
 5255 continue
      fm2=fm
      if(fm.le.i)fm2=i+1
      do 5265 m=fm2,lm
         im=itr(m)+i
         jbim=eoffjb+eadd(im)
         mkac=toff(m,k,zlac)+tadd(a,c)
         vale=vale+ein(jbim)*t2(mkac)
 5265 continue
c
      msym=ieor(jsym,kcsym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      lm2=lm
      if(lm.gt.j)lm2=j
      do 5256 m=fm,lm2
         jm=itr(j)+m
         kcjm=eoffkc+eadd(jm)
         miba=toff(m,i,zlba)+tadd(b,a)
         vale=vale+ein(kcjm)*t2(miba)
 5256 continue
      fm2=fm
      if(fm.le.j)fm2=j+1
      do 5266 m=fm2,lm
         jm=itr(m)+j
         kcjm=eoffkc+eadd(jm)
         miba=toff(m,i,zlba)+tadd(b,a)
         vale=vale+ein(kcjm)*t2(miba)
 5266 continue
c
      w1(a,b,c,mythrd)=w1(a,b,c,mythrd)+(valf-vale)/delta
 5340 continue
 5350 continue
 5360 continue
c
               do 5460 a=1,nv
                  if(fzv(a).eq.1)go to 5460
                  asym=orbsym(a+no)
                  ia=itr(a)+i
                  ai=ia
                  iasym=ieor(isym,asym)
                  aa=ipq(a+no)+a+no
                  do 5450 b=1,nv
                     if(fzv(b).eq.1)go to 5450
                     bsym=orbsym(b+no)
                     absym=ieor(asym,bsym)
                     jbsym=ieor(jsym,bsym)
                     jb=itr(b)+j
                     bb=ipq(b+no)+b+no
                     do 5440 c=1,nv
                        if(fzv(c).eq.1)go to 5440
                        csym=orbsym(c+no)
                        abcsym=ieor(absym,csym)
                        bcsym=ieor(bsym,csym)
                        if(ijksym.ne.abcsym)go to 5440
                        kcsym=ieor(ksym,csym)
                        kc=itr(c)+k
                        cc=ipq(c+no)+c+no
c
      delta=fock(ii)+fock(jj)+fock(kk)-fock(aa)-fock(bb)-fock(cc)
c
      avalw=4.d0*w1(a,b,c,mythrd)+  w1(c,a,b,mythrd)+  w1(b,c,a,mythrd) 
      bvalw= w1(a,b,c,mythrd)+4.0d0*w1(c,a,b,mythrd)+  w1(b,c,a,mythrd) 
      cvalw= w1(a,b,c,mythrd)+  w1(c,a,b,mythrd)+4.0d0*w1(b,c,a,mythrd) 
      pvalw=2.d0*(w1(c,b,a,mythrd)+w1(a,c,b,mythrd)+w1(b,a,c,mythrd))
      qvalw=avalw-pvalw
c     rvalw=avalw-pvalw
      rvalw1=avalw-pvalw
      rvalw2=bvalw-pvalw
      rvalw3=cvalw-pvalw
c
      if(i.eq.j.or.j.eq.k)then
         qvalw=qvalw*half 
         rvalw1=rvalw1*half 
         rvalw2=rvalw2*half 
         rvalw3=rvalw3*half 
      endif
c
      val1=w1(a,b,c,mythrd)+v1(a,b,c,mythrd)
      xval=val1*qvalw*delta
      xet =xet +xval
c
c >>> nvv
c
      fd=flov(csym+1,3)-no
      ld=flov(csym+1,4)-no
      do 5660 d=fd,ld
         val1=w1(a,b,d,mythrd)+v1(a,b,d,mythrd)
         val3=w1(b,d,a,mythrd)+v1(b,d,a,mythrd)
         val2=w1(d,a,b,mythrd)+v1(d,a,b,mythrd)
         term=+val1*rvalw1+val2*rvalw2+val3*rvalw3
         nvv(c,d,mythrd)=nvv(c,d,mythrd)+term
c        nvv(d,c)=nvv(d,c)+term
 5660 continue
c
 5440 continue
 5450 continue
 5460 continue
c
 5270 continue
 5280 continue
#ifndef SGI_PARALLEL
c     write(6,*)' done with virtual mo =',a
      write(4,*)' done with occupied mo =',i
C      if (prcntr('IS_ON BRIEF').eq.0) call timit(4,6)
      close(unit=4,err=51)
   51 continue
      call ffile(4,'term',1)
#endif
 5290 continue
c
      xet=xet*2.0d0
c
      write(6,*)
      write(6,*)' 2nd triples contribution to ccsd(t)',xet
      write(6,*)
      call fioflu(6)
c
c sum the contributions from each thread
c
      do 4090 ith=2,nthrd
         do 4080 j=1,no
            do 4070 i=1,no
               moo(i,j,1)= moo(i,j,1) + moo(i,j,ith)
 4070       continue
 4080    continue
 4090 continue
      do 4190 ith=2,nthrd
         do 4180 j=1,nv
            do 4170 i=1,nv
               nvv(i,j,1)= nvv(i,j,1) + nvv(i,j,ith)
 4170       continue
 4180    continue
 4190 continue
      do 4290 ith=2,nthrd
         do 4280 j=1,nshov
            dden(j,1) = dden(j,1) + dden(j,ith)
 4280    continue
 4290 continue
c
c >>> write file76  (read by ccqg)
c
      do 1990 i=1,no
         do 1980 j=1,no
            moo(i,j,1)=moo(i,j,1)*2.0d0
 1980    continue
 1990 continue
c
      do 1970 i=1,nv
         do 1960 j=1,nv
            nvv(i,j,1)=nvv(i,j,1)*2.0d0 
 1960    continue
 1970 continue
c
      do 1995 i=1,nshov
         dden(i,1)=dden(i,1)*2.0d0
 1995 continue
c
c     write(6,*)' new moo length is',no*no
c     write(6,'(8f10.6)')  moo 
c     call matout(moo,no,no,no,no,6)
c     write(6,*)' new nvv length is',nv*nv
c     write(6,'(8f10.6)')  nvv 
c     call matout(nvv,nv,nv,nv,nv,6)
c     write(6,*)' new dden length is',nshov
c     write(6,'(8f10.6)')  dden 
c     write(6,*)' new eden length is',nslov
c     write(6,'(8f10.6)')  eden 
c
      itap76=76
      call rfile(itap76)
      call wwritw(itap76,moo,intowp(no*no),  1,i76)
      call wwritw(itap76,nvv,intowp(nv*nv),i76,i76)
      call wwritw(itap76,dden,intowp(nshov),i76,i76)
      call wwritw(itap76,eden,intowp(nslov),i76,i76)
      call rclose(itap76,3)
c
      return
      end
