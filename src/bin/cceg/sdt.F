      subroutine sdt(t1,t2,din,ein,fin,fock,v,w,
     .               orbsym,flov,fzo,fzv,toff,tadd,zlx,ipq,
     .               itr,itv,eoff,eadd,foff,fadd,doff,dadd,
     .               no,nv,nirred,xet,xet4,tra,tva,trstr,nthrd)
      implicit integer(a-z)
      real*8 din(*),ein(*),fin(*),t1(no,nv),t2(*),fock(*)
      integer orbsym(*),flov(nirred,4),fzo(*),fzv(*),
     .        toff(no,no,2),tadd(nv,nv),zlx(nv,nv),
     .        itr(*),itv(*),doff(*),dadd(*),eoff(*),eadd(*),
     .        foff(*),fadd(*),ipq(*),tra(no,no),tva(nv,nv)
      character*4 trstr
      real*8 val,valf,vale,val1,val2,delta,xet,xet4,fac,xval,aval,bval,
     .       third
      real*8 w(no,no,no,nthrd),v(no,no,no,nthrd)
      integer prcntr
c
      write(4,*)'  '
      write(4,*)'Entering ccsd(t) energy calculation  '
      write(6,*)'  '
      write(6,*)'Entering ccsd(t) energy calculation  '
#ifdef SGI_PARALLEL
C$    write(6,*)'parallel sgi version -- no restarts'
C$    write(6,*)'running with ', nthrd, ' threads'
#endif
c
      third=1.0d0/3.0d0
      xet=0.0d0
      xet4=0.0d0
c
      if(trstr.eq.'YES ')then
         call ffile(9,'trstr',1)
         backspace 9
         backspace 9
         backspace 9
         read(9,'(i4)')lmo
         read(9,'(e20.15)')xet
         nextmo = lmo + 1
         write(6,*)'restarting calculation from virtual mo',nextmo
         write(6,*)'restarting calculation with xe',xet
      else
         call ffile(9,'trstr',0)
         nextmo=1
      endif
c
#ifdef SGI_PARALLEL
C$DOACROSS local(a,asym,aa,b,bsym,absym,zlba,zlab,bb,c,
C$&              zlbc,zlcb,zlac,zlca,cc,csym,abcsym,
C$&              i,isym,ia,iasym,ii,faddia,eoffia,
C$&              j,jsym,ijsym,jbsym,jb,jj,faddjb,eoffjb,
C$&              k,ksym,ijksym,kcsym,kc,kk,faddkc,eoffkc,
C$&              iajb,jbkc,kcia,val,valf,vale,fsym,ff,lf,lf2,
C$&              f,bf,iabf,zlcf,kjcf,cf,jbcf,ikaf,
C$&              ff2,zlaf,af,kcaf,zlbf,jibf,iacf,jkbf,
C$&              jbaf,kicf,kcbf,ijaf,msym,
C$&              fm,lm,lm2,m,jm,iajm,mkbc,fm2,jbkm,mica,
C$&              im,kcim,mjab,km,iakm,mjcb,jbim,mkac,kcjm,
C$&              miba,
C$&              aval,bval,val2,delta,val1,xval,mythrd
C$&              ),
C$&        reduction(xet),
C$&        mp_schedtype=dynamic,
C$&        chunk=1
#endif
      do 950 a=nextmo,nv
#ifdef SGI_PARALLEL
         mythrd = mp_my_threadnum() + 1
#else
         mythrd = 1
#endif
         if(fzv(a).eq.1)go to 950
         asym=orbsym(a+no)
         aa=ipq(a+no)+a+no
         do 850 b=1,a 
            if(fzv(b).eq.1)go to 850
            bsym=orbsym(b+no)
            absym=ieor(asym,bsym)
            zlba=zlx(b,a)
            zlab=zlx(a,b)
            bb=ipq(b+no)+b+no
            do 750 c=1,b 
               if(fzv(c).eq.1)go to 750
c              if(c.eq.a)go to 750
               zlbc=zlx(b,c)
               zlcb=zlx(c,b)
               zlac=zlx(a,c)
               zlca=zlx(c,a)
               cc=ipq(c+no)+c+no
               csym=orbsym(c+no)
               abcsym=ieor(absym,csym)
c
               call zero(w(1,1,1,mythrd),no*no*no)
               call zero(v(1,1,1,mythrd),no*no*no)
c
               do 650 i=1,no
                  if(fzo(i).eq.1)go to 650
                  isym=orbsym(i)
                  ia=itr(a)+i
                  iasym=ieor(isym,asym)
                  ii=ipq(i)+i
                  faddia=fadd(ia)
                  eoffia=eoff(ia)
                  do 550 j=1,no
                     if(fzo(j).eq.1)go to 550
                     jsym=orbsym(j)
                     ijsym=ieor(isym,jsym)
                     jbsym=ieor(jsym,bsym)
                     jb=itr(b)+j
                     jj=ipq(j)+j
                     faddjb=fadd(jb)
                     eoffjb=eoff(jb)
                     do 450 k=1,no
                        if(fzo(k).eq.1)go to 450
                        ksym=orbsym(k)
c                       if(k.eq.i)go to 450
                        ijksym=ieor(ijsym,ksym)
                        if(ijksym.ne.abcsym)go to 450
                        kcsym=ieor(ksym,csym)
                        kc=itr(c)+k
                        kk=ipq(k)+k
                        faddkc=fadd(kc)
                        eoffkc=eoff(kc)
c
      iajb=doff(max0(ia,jb))+dadd(min0(ia,jb))
      jbkc=doff(max0(kc,jb))+dadd(min0(kc,jb))
      kcia=doff(max0(kc,ia))+dadd(min0(kc,ia))
      val=din(iajb)*t1(k,c)+din(jbkc)*t1(i,a)+din(kcia)*t1(j,b)
      v(i,j,k,mythrd)=v(i,j,k,mythrd)+val
c               
      valf=0.0d0
      vale=0.0d0
c
      fsym=ieor(bsym,iasym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.gt.b)lf2=b
      do 351 f=ff,lf2
         bf=itv(b)+f
         iabf=foff(bf)+faddia
         zlcf=zlx(c,f)
         kjcf=toff(k,j,zlcf)+tadd(c,f)
         valf=valf+fin(iabf)*t2(kjcf)
  351 continue
      ff2=ff
      if(ff.le.b)ff2=b+1
      do 361 f=ff2,lf
         bf=itv(f)+b
         iabf=foff(bf)+faddia
         kjcf=toff(k,j,2)+tadd(c,f)
         valf=valf+fin(iabf)*t2(kjcf)
  361 continue
c
      fsym=ieor(csym,jbsym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.gt.c)lf2=c
      do 352 f=ff,lf2
         cf=itv(c)+f
         jbcf=foff(cf)+faddjb
         ikaf=toff(i,k,1)+tadd(a,f)
         valf=valf+fin(jbcf)*t2(ikaf)
  352 continue
      ff2=ff
      if(ff.le.c)ff2=c+1
      do 362 f=ff2,lf
         cf=itv(f)+c
         jbcf=foff(cf)+faddjb
         zlaf=zlx(a,f)
         ikaf=toff(i,k,zlaf)+tadd(a,f)
         valf=valf+fin(jbcf)*t2(ikaf)
  362 continue
c
      fsym=ieor(asym,kcsym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.gt.a)lf2=a
      do 353 f=ff,lf2
         af=itv(a)+f
         kcaf=foff(af)+faddkc
         zlbf=zlx(b,f)
         jibf=toff(j,i,zlbf)+tadd(b,f)
         valf=valf+fin(kcaf)*t2(jibf)
  353 continue
      ff2=ff
      if(ff.le.a)ff2=a+1
      do 363 f=ff2,lf
         af=itv(f)+a
         kcaf=foff(af)+faddkc
         jibf=toff(j,i,2)+tadd(b,f)
         valf=valf+fin(kcaf)*t2(jibf)
  363 continue
c
      fsym=ieor(csym,iasym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.gt.c)lf2=c
      do 354 f=ff,lf2
         cf=itv(c)+f
         iacf=foff(cf)+faddia
         jkbf=toff(j,k,1)+tadd(b,f)
         valf=valf+fin(iacf)*t2(jkbf)
  354 continue
      ff2=ff
      if(ff.le.c)ff2=c+1
      do 364 f=ff2,lf
         cf=itv(f)+c
         iacf=foff(cf)+faddia
         zlbf=zlx(b,f)
         jkbf=toff(j,k,zlbf)+tadd(b,f)
         valf=valf+fin(iacf)*t2(jkbf)
  364 continue
c
      fsym=ieor(asym,jbsym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.gt.a)lf2=a
      do 355 f=ff,lf2
         af=itv(a)+f
         jbaf=foff(af)+faddjb
         zlcf=zlx(c,f)
         kicf=toff(k,i,zlcf)+tadd(c,f)
         valf=valf+fin(jbaf)*t2(kicf)
  355 continue
      ff2=ff
      if(ff.le.a)ff2=a+1
      do 365 f=ff2,lf
         af=itv(f)+a
         jbaf=foff(af)+faddjb
         kicf=toff(k,i,2)+tadd(c,f)
         valf=valf+fin(jbaf)*t2(kicf)
  365 continue
c
      fsym=ieor(bsym,kcsym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      lf2=lf
      if(lf.gt.b)lf2=b
      do 356 f=ff,lf2
         bf=itv(b)+f
         kcbf=foff(bf)+faddkc
         ijaf=toff(i,j,1)+tadd(a,f)
         valf=valf+fin(kcbf)*t2(ijaf)
  356 continue
      ff2=ff
      if(ff.le.b)ff2=b+1
      do 366 f=ff2,lf
         bf=itv(f)+b
         kcbf=foff(bf)+faddkc
         zlaf=zlx(a,f)
         ijaf=toff(i,j,zlaf)+tadd(a,f)
         valf=valf+fin(kcbf)*t2(ijaf)
  366 continue
c
      msym=ieor(jsym,iasym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      lm2=lm
      if(lm.gt.j)lm2=j
      do 251 m=fm,lm2
         jm=itr(j)+m
         iajm=eoffia+eadd(jm)
         mkbc=toff(m,k,zlbc)+tadd(b,c)
         vale=vale+ein(iajm)*t2(mkbc)
  251 continue
      fm2=fm
      if(fm.le.j)fm2=j+1
      do 261 m=fm2,lm
         jm=itr(m)+j
         iajm=eoffia+eadd(jm)
         mkbc=toff(m,k,zlbc)+tadd(b,c)
         vale=vale+ein(iajm)*t2(mkbc)
  261 continue
c
      msym=ieor(ksym,jbsym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      lm2=lm
      if(lm.gt.k)lm2=k
      do 252 m=fm,lm2
         km=itr(k)+m
         jbkm=eoffjb+eadd(km)
         mica=toff(m,i,zlca)+tadd(c,a)
         vale=vale+ein(jbkm)*t2(mica)
  252 continue
      fm2=fm
      if(fm.le.k)fm2=k+1
      do 262 m=fm2,lm
         km=itr(m)+k
         jbkm=eoffjb+eadd(km)
         mica=toff(m,i,zlca)+tadd(c,a)
         vale=vale+ein(jbkm)*t2(mica)
  262 continue
c
      msym=ieor(isym,kcsym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      lm2=lm
      if(lm.gt.i)lm2=i
      do 253 m=fm,lm2
         im=itr(i)+m
c        im=tra(i,m)
         kcim=eoffkc+eadd(im)
         mjab=toff(m,j,zlab)+tadd(a,b)
         vale=vale+ein(kcim)*t2(mjab)
  253 continue
      fm2=fm
      if(fm.le.i)fm2=i+1
      do 263 m=fm2,lm
         im=itr(m)+i
c        im=tra(i,m)
         kcim=eoffkc+eadd(im)
         mjab=toff(m,j,zlab)+tadd(a,b)
         vale=vale+ein(kcim)*t2(mjab)
  263 continue
c
      msym=ieor(ksym,iasym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      lm2=lm
      if(lm.gt.k)lm2=k
      do 254 m=fm,lm2
         km=itr(k)+m
         iakm=eoffia+eadd(km)
         mjcb=toff(m,j,zlcb)+tadd(c,b)
         vale=vale+ein(iakm)*t2(mjcb)
  254 continue
      fm2=fm
      if(fm.le.k)fm2=k+1
      do 264 m=fm2,lm
         km=itr(m)+k
         iakm=eoffia+eadd(km)
         mjcb=toff(m,j,zlcb)+tadd(c,b)
         vale=vale+ein(iakm)*t2(mjcb)
  264 continue
c
      msym=ieor(isym,jbsym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      lm2=lm
      if(lm.gt.i)lm2=i
      do 255 m=fm,lm2
         im=itr(i)+m
         jbim=eoffjb+eadd(im)
         mkac=toff(m,k,zlac)+tadd(a,c)
         vale=vale+ein(jbim)*t2(mkac)
  255 continue
      fm2=fm
      if(fm.le.i)fm2=i+1
      do 265 m=fm2,lm
         im=itr(m)+i
         jbim=eoffjb+eadd(im)
         mkac=toff(m,k,zlac)+tadd(a,c)
         vale=vale+ein(jbim)*t2(mkac)
  265 continue
c
      msym=ieor(jsym,kcsym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      lm2=lm
      if(lm.gt.j)lm2=j
      do 256 m=fm,lm2
         jm=itr(j)+m
         kcjm=eoffkc+eadd(jm)
         miba=toff(m,i,zlba)+tadd(b,a)
         vale=vale+ein(kcjm)*t2(miba)
  256 continue
      fm2=fm
      if(fm.le.j)fm2=j+1
      do 266 m=fm2,lm
         jm=itr(m)+j
         kcjm=eoffkc+eadd(jm)
         miba=toff(m,i,zlba)+tadd(b,a)
         vale=vale+ein(kcjm)*t2(miba)
  266 continue
c
      w(i,j,k,mythrd)=w(i,j,k,mythrd)+valf-vale
  450 continue
  550 continue
  650 continue
               do 655 i=1,no
                  if(fzo(i).eq.1)go to 655
                  isym=orbsym(i)
                  ia=itr(a)+i
                  iasym=ieor(isym,asym)
                  ii=ipq(i)+i
                  do 555 j=1,no
                     if(fzo(j).eq.1)go to 555
                     jsym=orbsym(j)
                     ijsym=ieor(isym,jsym)
                     jbsym=ieor(jsym,bsym)
                     jb=itr(b)+j
                     jj=ipq(j)+j
                     do 455 k=1,no
                        if(fzo(k).eq.1)go to 455
                        ksym=orbsym(k)
                        ijksym=ieor(ijsym,ksym)
                        if(ijksym.ne.abcsym)go to 455
                        kcsym=ieor(ksym,csym)
                        kc=itr(c)+k
                        kk=ipq(k)+k
c
      aval=4.d0*w(i,j,k,mythrd)+w(k,i,j,mythrd)+w(j,k,i,mythrd) 
      bval=w(k,j,i,mythrd)+w(i,k,j,mythrd)+w(j,i,k,mythrd)
      bval=bval+bval
      val2=aval-bval
      if(a.eq.b.or.b.eq.c)then
         aval=4.0d0*w(i,j,k,mythrd)+w(k,i,j,mythrd)+w(j,k,i,mythrd)
          bval=4.0d0*w(k,j,i,mythrd)+w(i,k,j,mythrd)+w(j,i,k,mythrd)
         val2=(aval-bval)*third
      endif
      delta=fock(ii)+fock(jj)+fock(kk)-fock(aa)-fock(bb)-fock(cc)
      val1=w(i,j,k,mythrd)+v(i,j,k,mythrd)
      xval=val1*val2/delta
      xet =xet +xval
c     write(6,'(3i3,5f12.8)')ai,bj,ck,val1,aval,bval,val2,xval
cges triples contribution to mp4 may be obtained with formula below
c    iff w is contsructed with mp2 t2's
c     xet4=xet4+(w(i,j,k,mythrd)               )*val2/delta
  455 continue
  555 continue
  655 continue
c
  750 continue
  850 continue

#ifndef SGI_PARALLEL
      write(6,*)' done with virtual mo =',a,' energy=',xet
      write(4,*)' done with virtual mo =',a,' energy=',xet
C      if (prcntr('IS_ON BRIEF').eq.0) call timit(4,6)
      call fioflu(6)
#ifndef SGI_PARALLEL
      close(unit=4,err=51)
   51 continue
      call ffile(4,'term',1)
#else
      call fioflu(4)
#endif
c
         write(9,'(i4)')a
         write(9,'(f20.15)')xet
#endif
  950 continue
c
      xet=xet*2.0d0
c
c     xet4=xet4*2.0d0
c
      write(6,*)
      write(6,*)'triples contribution to ccsd(t)',xet
c     write(6,*)'triples contribution to mp4    ',xet4
      write(6,*)
      call fioflu(6)
c
      return
      end
