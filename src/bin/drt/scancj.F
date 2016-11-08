      subroutine scancj(ntype,bfnum,bfsym,bfcode)
C
C
      IMPLICIT INTEGER (A-Z)
cets050191 yeah, i know, it's not very elegant, but it works for
c          most things.
c
#include <error.h>
      integer errcod, frdia1, frdc, fcnt, fcnt1, frdca2
      character*6 dertyp
      character*8 opentp, cityp
      logical res,mr,mref
      integer numint
      character*1 multrf,valenc,clef,junk1,bfkey,getkey
      character*3 codes,words
      character*4 label
      integer docc(8),socc(8),fdocc(8),fuocc(8)
      integer valdoc(8),valuoc(8),valsoc(8)
      integer uocc(8),tocc(8),valorb(8)
      integer mrorbs(100)
      character*3 refcod(500)
      integer ITEMP(200),OCCTEMP(8)
      integer NSYMHF,point
      character*4 LAB(8),HFLAB(8)

C
      COMMON /CBFKEY/ BFKEY(1000)
      COMMON /DRTINF/ NA,NB,NS,NESPEC,MAXB,LEVFRM,LEVVAL,LEVOPN,LEVMUL
     #,               LEVOCC,SPEC,SSPESH,VAL
      COMMON /DIMENS/ NBF,NSYM,NORBS,NROWSP,NROW4P,NROWS,NROWS4
     #,               NLEVS,NREFS,NROWOC,NROW4O,NWKS,NWKSOC,NLEVOC
     #,               ORBFRM,SYMORB,NUMIJ,NGROUP,NUMINT,NMAX,NSPC,NVREF
     #,               NIJVIR
      COMMON /DRTCOD/ NCODES,DELA(9),DELB(9),DELELE(9)
     #,               NTYPES,VIRTUL,OCCUPD,VALOCC,RESCOR,RESVIR,FROZEN
     #,               VALVIR,OPENSH,MULTI,SPESHL,MULTRF,VALENC
     #,               WORDS(6,9),CODES(9)
      COMMON /CODE/   FZC, FZV, COR, VIR, DOC, UOC, ALP, BET, SPE
      COMMON /TAPES/  OUT,ERROUT,INPUT,DRTTAP
      COMMON /OPTION/ OPTION(16),IDENT(26),LABEL(26)
C
      DIMENSION NTYPE(NTYPES),BFNUM(NTYPES,NBF)
      DIMENSION BFSYM(NBF),BFCODE(NREFS,NBF)

      IEND=1
      MPOINT=200
      MCONST=200
      MCALCS=100
      ITAP30=30

C
      WRITE (OUT,11)
   11 FORMAT (//,T28,'**** ORBITAL INFORMATION ****',//
     #,      T21,'BASIS FCTS',T33,'#',T35,'KEY',T42,'TYPE OF ORBITAL'
     #,      T59,'SYM',/,T21,'----------',T33,'-',T35,'---',T42
     #,      '---------------',T59,'---')
c
      BF=0

      dertyp='NONE  '
      errcod = frdc('DERTYPE',dertyp)

CeebCAS
      cityp = 'CI      '
      errcod = frdc('WFN',cityp)
CeebCAS

      opentp='NONE    '
      errcod = frdc('OPENTYPE',opentp)

      nexcit=2
      errcod = frdi('EX_LVL',nexcit)

      res=.false.
      if(dertyp.ne.'NONE  ') res=.true.

CeebCAS
      if (cityp.eq.'CASSCF  ') res=.true.
CeebCAS

      do 110 i=1,100
 110  mrorbs(i)=0

      if(nrefs.gt.1) then
        errcod = fcnt('MULTI_ORBS',nmro)
        if(errcod.ne.EOK) then
          write(6,*)' nrefs =',nrefs
          write(6,*)' you must supply the numbers of the mos which'
          write(6,*)' differ in each reference with the'
          write(6,*)' MULTI_ORBS vector'
cets one for the Welsh, and one for the Irish
          call pabort
          call qabort
        endif

        errcod = frdia1('MULTI_ORBS',mrorbs,nmro)

        errcod = fcnt('REFERENCES',nr)
        if(errcod.ne.EOK.or.nr.ne.nrefs) then
          write(6,*)' number of REFERENCES ne NREFS'
          call pabort
          call qabort
        endif

        errcod = frdca2('REFERENCES',refcod,nmro,nrefs)
        if(errcod.ne.EOK) then
          write(6,*)' problem reading REFERENCES'
          call pabort
          call qabort
        endif
      endif

      do 100 i=1,8
        docc(i)=0
        socc(i)=0
        fdocc(i)=0
        fuocc(i)=0
        valdoc(i)=0
        valsoc(i)=0
        valuoc(i)=0
        uocc(i)=0
 100  continue

      DO 7000 I=1,MPOINT
         ITEMP(I)=0
 7000 continue
      call rfile(itap30)
      call wreadw (itap30,itemp,mconst,101,junk)
      
      NIRRED = ITEMP(28)
      NSYMHF = ITEMP(41)
      MXCOEF = ITEMP(42)
      MSFRU = ITEMP(18)
      IOPEN = ITEMP(43)
      call wreadw (itap30,itemp,mpoint,301,junk)
      call wreadw (itap30,lab,nirred,itemp(16),junk)
      call wreadw (itap30,itemp,mcalcs,501,junk)
      call wreadw (itap30,point,1,itemp(1)+60,junk)
      point = point+(2*mxcoef)+(2*msfru)
      call wreadw (itap30,hflab,nsymhf,point,point)
      call wreadw (itap30,docc,nsymhf,point+nsymhf,point)
CTDC.  This was "if(iopen.eq.1) ..."; I don't know why this worked before.
      if(iopen.ne.0) call wreadw (itap30,socc,nsymhf,point,point)
      M=0
      if(nirred.ne.nsymhf) THEN
         do 7100 I=1,nirred
            if(lab(i+m).ne.hflab(i)) THEN
              do 7110 k=nirred,i+1,-1
                 docc(k)=docc(k-1)
                 socc(k)=socc(k-1)
 7110          continue
               docc(i)=0
               socc(i)=0
               m=m+1
            endif
 7100    continue
      endif
      call rclose(ITAP30,3)

Crawdad.9.19.93 Check length of frozen_docc and frozen_uocc vectors.
      nfrz = 0
      errcod = fcnt('FROZEN_DOCC',nfrz)
      if(nfrz.ne.nsym.and.nfrz.ne.0) then
		write(6,*) 'Hmmm...your FROZEN_DOCC array is too small.'
		write(6,*) 'It is',nfrz,'when it should be',nsym
		call pabort
		call qabort
      endif
Crawdad.9.19.93 Check length of frozen_docc and frozen_uocc vectors.
      errcod = frdia1('FROZEN_DOCC',fdocc,nsym)
Crawdad.9.19.93 Check length of frozen_docc and frozen_uocc vectors.
      nfrz = 0
      errcod = fcnt('FROZEN_UOCC',nfrz)
      if(nfrz.ne.nsym.and.nfrz.ne.0) then
		write(6,*) 'Hmmm...your FROZEN_UOCC array is too small.'
		write(6,*) 'It is',nfrz,'when it should be',nsym
		call pabort
		call qabort
      endif
Crawdad.9.19.93 Check length of frozen_docc and frozen_uocc vectors.
      errcod = frdia1('FROZEN_UOCC',fuocc,nsym)
      if(option(3).gt.0) then
        errcod = frdia1('VAL_ORB',valorb,nsym)

        do 300 i=1,nsym
          if(socc(i).le.valorb(i)) then
            valsoc(i)=socc(i)
          else
            valsoc(i)=valorb(i)
          endif
          valorb(i)=valorb(i)-valsoc(i)
  300   continue

        do 301 i=1,nsym
          if((docc(i)-fdocc(i)).le.valorb(i)) then
            valdoc(i)=docc(i)-fdocc(i)
          else
            valdoc(i)=valorb(i)
          endif
          valuoc(i)=valorb(i)-valdoc(i)
  301   continue
      endif

      errcod = frdia1('UOCC',uocc,nsym)
      if(errcod.ne.EOK) then
        call rfile(34)
        pos=26+2*nsym
        call wreadw(34,tocc,nsym,pos,junk)
        call rclose(34,3)
        
        do 102 i=1,nsym
          uocc(i)=tocc(i)-docc(i)-socc(i)-fuocc(i)-valuoc(i)
  102   continue
      endif

      ndocc=0
      nsocc=0
      nfdocc=0
      nfuocc=0
      nvdocc=0
      nvuocc=0
      nuocc=0
      
      do 101 i=1,nsym
        docc(i) = docc(i)-fdocc(i)-valdoc(i)
        ndocc = ndocc+docc(i)
        socc(i) = socc(i)-valsoc(i)
        nsocc = nsocc+socc(i)
        nfdocc = nfdocc+fdocc(i)
        nfuocc = nfuocc+fuocc(i)
        nvdocc = nvdocc+valdoc(i)
        nvsocc = nvsocc+valsoc(i)
        nvuocc = nvuocc+valuoc(i)
        nuocc = nuocc+uocc(i)
 101  continue

      if(opentp.eq.'TWOCON  ') then
        if(nrefs.gt.1) then
          write(6,*)'nrefs=',nrefs
          write(6,*)'TWOCON is a special case, so you should not'
          write(6,*)' use multiple references'
          call qabort
        else if((nvdocc+nvsocc+nvuocc).gt.0) then
          write(6,*) 'use old style input for GVB multireference ci'
          call qabort
        endif
      endif

      nvalor=0
      if((nvdocc+nvsocc+nvuocc).gt.0) nvalor=1

      ntype(occupd) = ndocc
      ntype(virtul) = nuocc
      ntype(valocc) = nvdocc+nvsocc
      ntype(valvir) = nvuocc
      if(res) then
        ntype(rescor) = nfdocc
        ntype(resvir) = nfuocc
      else
        ntype(frozen) = nfdocc+nfuocc
      endif

      if(opentp.ne.'NONE    '.and.opentp.ne.'TWOCON  '.and.
     +                           (nvalor.gt.0.or.nexcit.gt.2)) then
        ntype(opensh) = nsocc
      else if(opentp.eq.'HIGHSPIN') then
        ntype(opensh) = nsocc
      else if(opentp.eq.'SINGLET '.or.opentp.eq.'TWOCON  ') then
        ntype(speshl)=nsocc
      endif

      norbs = ntype(occupd)+ntype(virtul)+ntype(valocc)+ntype(valvir)
     1      + ntype(rescor)+ntype(resvir)+ntype(opensh)+ntype(speshl)
C
      nfz=0
      nfc=0
      nd=0
      nvd=0
      ns=0
      nvu=0
      nu=0
      nfu=0
      nmul=0
      noporb=0

      do 200 i=1,nsym
        if(res) then
          if(fdocc(i).gt.0) write(out,1)bf+1,bf+fdocc(i),fdocc(i),' ',
     +                                   (words(j,cor),j=1,6),i
          do 201 j=1,fdocc(i)
            bf=bf+1
            nfc=nfc+1
            bfnum(rescor,nfc)=bf
            bfsym(bf)=i-1
            bfkey(bf)=' '
            bfcode(1,bf)=cor
  201     continue
        else
          if(fdocc(i).gt.0) write(out,1)bf+1,bf+fdocc(i),fdocc(i),' ',
     +                                   (words(j,fzc),j=1,6),i
          do 202 j=1,fdocc(i)
            bf=bf+1
            nfz=nfz+1
            bfnum(frozen,nfz)=bf
            bfsym(bf)=i-1
            bfkey(bf)=' '
            bfcode(1,bf)=fzc
  202     continue
        endif
        if(docc(i).gt.0) write(out,1)bf+1,bf+docc(i),docc(i),' ',
     +                                   (words(j,doc),j=1,6),i
        do 203 j=1,docc(i)
          bf=bf+1
          mr=mref(bf,mrorbs,nmro,nrefs)
          if(mr) then
            nmul=nmul+1
            bfnum(multi,nmul)=bf
            bfkey(bf)='/'
          else
            nd=nd+1
            bfnum(occupd,nd)=bf
            bfkey(bf)=' '
          endif
          bfsym(bf)=i-1
          bfcode(1,bf)=doc
  203   continue
        if(valdoc(i).gt.0) write(out,1)bf+1,bf+valdoc(i),valdoc(i),'%',
     +                                   (words(j,doc),j=1,6),i
        do 204 j=1,valdoc(i)
          bf=bf+1
          nvd=nvd+1
          bfnum(valocc,nvd)=bf
          bfsym(bf)=i-1
          bfkey(bf)='%'
          bfcode(1,bf)=doc
  204   continue
        if(opentp.eq.'TWOCON  '.or.
     +    (opentp.eq.'SINGLET '.and.nrefs.eq.1.and.
     +                    nexcit.le.2.and.nvalor.eq.0)) then
          if(socc(i).gt.0) write(out,1)bf+1,bf+socc(i),socc(i),' ',
     +                                   (words(j,spe),j=1,6),i
          do 225 j=1,socc(i)
            sspesh=ieor(i-1,sspesh)
            bf=bf+1
            ns=ns+1
            bfnum(speshl,ns)=bf
            bfkey(bf)=' '
            bfcode(1,bf)=spe
            bfsym(bf)=i-1
  225     continue
        else if(opentp.eq.'SINGLET '.and.noporb.eq.0) then
          if(socc(i).gt.0) write(out,1)bf+1,bf+socc(i),socc(i),' ',
     +                                   (words(j,bet),j=1,6),i
          option(4)=1
          do 215 j=1,socc(i)
            noporb=1
            bf=bf+1
            mr=mref(bf,mrorbs,nmro,nrefs)
            if(mr) then
              nmul=nmul+1
              bfnum(multi,nmul)=bf
              bfkey(bf)='/'
            else
              ns=ns+1
              bfnum(opensh,ns)=bf
              bfkey(bf)=' '
            endif
            bfcode(1,bf)=bet
            bfsym(bf)=i-1
  215     continue
        else
          if(socc(i).gt.0) write(out,1)bf+1,bf+socc(i),socc(i),' ',
     +                                   (words(j,alp),j=1,6),i
          do 206 j=1,socc(i)
            bf=bf+1
            mr=mref(bf,mrorbs,nmro,nrefs)
            if(mr) then
              nmul=nmul+1
              bfnum(multi,nmul)=bf
              bfkey(bf)='/'
            else
              ns=ns+1
              bfnum(opensh,ns)=bf
              bfkey(bf)=' '
            endif
            bfsym(bf)=i-1
            bfcode(1,bf)=alp
  206     continue
        endif
        if(opentp.eq.'TWOCON  '.or.opentp.eq.'SINGLET '.and.
     +                                             noporb.eq.0) then
          if(valsoc(i).gt.0)write(out,1)bf+1,bf+valsoc(i),valsoc(i),'%',
     +                                   (words(j,bet),j=1,6),i
          do 217 j=1,valsoc(i)
            noporb=1
            bf=bf+1
            nvd=nvd+1
            bfnum(valocc,nvd)=bf
            bfsym(bf)=i-1
            bfkey(bf)='%'
            bfcode(1,bf)=bet
  217     continue
        else
          if(valsoc(i).gt.0)write(out,1)bf+1,bf+valsoc(i),valsoc(i),'%',
     +                                   (words(j,alp),j=1,6),i
          do 207 j=1,valsoc(i)
            bf=bf+1
            nvd=nvd+1
            bfnum(valocc,nvd)=bf
            bfsym(bf)=i-1
            bfkey(bf)='%'
            bfcode(1,bf)=alp
  207     continue
        endif
        if(valuoc(i).gt.0)write(out,1)bf+1,bf+valuoc(i),valuoc(i),'%',
     +                                   (words(j,uoc),j=1,6),i
        do 208 j=1,valuoc(i)
          bf=bf+1
          nvu=nvu+1
          bfnum(valvir,nvu)=bf
          bfsym(bf)=i-1
          bfkey(bf)='%'
          bfcode(1,bf)=uoc
  208   continue
        if(uocc(i).gt.0) write(out,1)bf+1,bf+uocc(i),uocc(i),' ',
     +                                   (words(j,uoc),j=1,6),i
        do 209 j=1,uocc(i)
          bf=bf+1
          mr=mref(bf,mrorbs,nmro,nrefs)
          if(mr) then
            nmul=nmul+1
            bfnum(multi,nmul)=bf
            bfkey(bf)='/'
          else
            nu=nu+1
            bfnum(virtul,nu)=bf
            bfkey(bf)=' '
          endif
          bfsym(bf)=i-1
          bfcode(1,bf)=uoc
  209   continue
        if(res) then
          if(fuocc(i).gt.0) write(out,1)bf+1,bf+fuocc(i),fuocc(i),' ',
     +                                   (words(j,vir),j=1,6),i
          do 211 j=1,fuocc(i)
            bf=bf+1
            nfu=nfu+1
            bfnum(resvir,nfu)=bf
            bfsym(bf)=i-1
            bfkey(bf)=' '
            bfcode(1,bf)=vir
  211     continue
        else
          if(fuocc(i).gt.0) write(out,1)bf+1,bf+fuocc(i),fuocc(i),' ',
     +                                   (words(j,fzv),j=1,6),i
          do 212 j=1,fuocc(i)
            bf=bf+1
            nfz=nfz+1
            bfnum(frozen,nfz)=bf
            bfsym(bf)=i-1
            bfkey(bf)=' '
            bfcode(1,bf)=fzv
  212     continue
        endif
  200 continue

    1 FORMAT (T22,I3,'-',I3,I4,2X,A1,T43,6A3,T61,I1)

      spec = ntype(speshl)
      if(nrefs.eq.1) return

      ntype(multi)=nmro

      write(out,12) 1
      do 400 i=1,nmro
        code=doc
        type=occupd
        bf=mrorbs(i)
        if(refcod(i).eq.'UOC') then
          type=virtul
          code=uoc
        else if(refcod(i).eq.'ALP') then
          type=opensh
          if(spec.eq.2) type=speshl
          code=alp
        else if(refcod(i).eq.'BET') then
          type=opensh
          if(spec.eq.2) type=speshl
          code=bet
        endif
        ntype(type)=ntype(type)-1
        bfcode(1,bf)=code
        write (out,13) bf,(words(j,bfcode(1,bf)),j=1,6),bfsym(bf)+1
 400  continue

      spec = 0
      mri=nmro
 
      do 401 i=2,nrefs
        write(out,12) i
        do 402 j=1,nbf
          bfcode(i,j)=bfcode(1,j)
 402    continue
        do 403 j=1,nmro
          mri=mri+1
          bf=mrorbs(j)
          code=doc
          if(refcod(mri).eq.'UOC') then
            code=uoc
          else if(refcod(mri).eq.'ALP') then
            code=alp
          else if(refcod(mri).eq.'BET') then
            code=bet
          endif
          bfcode(i,bf)=code
        write (out,13) bf,(words(k,bfcode(i,bf)),k=1,6),bfsym(bf)+1
 403    continue
 401  continue

   12 FORMAT (//,T30,'REFERENCE #',I3,//,T25,'BASIS FCT',T35,'TYPE OF
     * ',       'ORBITAL',T54,'SYM',/,T25,'---------',T35,'--------'
     #   ,       '-------',T54,'---')
   13 FORMAT (T25,I5,5X,6A3,2X,I1)

      return
      
      end

      logical function mref(bf,mrorbs,nmro,nrefs)
      implicit integer(a-z)
      dimension mrorbs(nmro)

      if(nrefs.eq.1) then
        mref=.false.
        return
      else
        do 10 i=1,nmro
          if(mrorbs(i).eq.bf) then
            mref=.true.
            return
          endif
  10    continue
      endif
      mref=.false.
      return
      end

















