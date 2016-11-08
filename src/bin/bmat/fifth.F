C     //////////////////////////////////////////////////////////
      subroutine fifth(inp,iout,inp2,ipun,igmup,nsim,nsym,nstre,nbnd,
     1  ntor,nlin1,nlin2,ninv,nout,ndisp,ndum,napdm,nfix,
     2  stype,satom,symms,symmc,ccin,ccc,symb,qqq,fii,
     3  b,ww,l,m,ia,ifix,ibohr,vmode,c,cc,f,fi,xa,ya,za,xm,xy,qq,qq1,
     4  nv,he,hn,n2v,dv,n3v,grc,gec,qd,fd,a,indi,ifstre,geom,geobf,
     5  core,lbls)
      IMPLICIT REAL*8 (A-H,O-Z)
#include <error.h>
      integer errcod, frdi, frdboo, frdfa2, frdfa1, frdia1, frdc, fcnt
      integer frdca1
      integer coallc

      integer core(1)
      real*8 symmc(nsym,100),ccin(4,ndisp),ccc(nsym)
      real*8 qqq(nsym),fii(nsym),b(3*napdm,nsym),ww(2),vmode(nsym)
      real*8 c(nsym*nsym),cc(nsym),f(3*napdm),fi(nsym)
      real*8 xa(napdm),ya(napdm),za(napdm),xm(3*napdm),xy(3*napdm)
      real*8 qq(nsym),qq1(nsym),nv(nsym),he(nsym,nsym)
      real*8 hn(nsym,nsym),n2v(nsym),dv(nsym),n3v(nsym),grc(nsym)
      real*8 gec(nsym),qd(nsym,nsym),fd(nsym,nsym),a(nsym*nsym)
      real*8 geom(3,na),gtmp(1500)

      integer satom(nsim,4),symms(nsym,0:100)
      integer l(nsym),m(nsym),ia(na),ifix(nfix),ibohr(nstre)
      integer indi(4,ndisp)

      integer inp,iout,inp2,ipun,igmup
      integer nsim,nsym,nstre,nbnd,ntor,nlin1,nlin2,ninv,nout,ndisp
      integer ndum,napdm,nfix

      character*5 symb(na)
      character*7 mftyp
      character*4 stype(nsim),updt
      character*80 geobf(2000),lbls(nsym)
      logical ifstre(nsym)

      INTEGER P,Q,NR
      CHARACTER WO*4,WK(27)*4,W1*2
      INTEGER INDT(4)
      LOGICAL WRIT,SFIRST,DISP
      DIMENSION CIN(4),ETA(4),ENERGY(4)
      DIMENSION CI(9),AE(3)
      COMMON /ENER/ TITLE(16)
      COMMON /GANZ/ NA,LOPT(27)
      COMMON /IO/ IIN,IIOUT,IINP2,IIPUN,IIGMUP
      COMMON /IOP/ IOP(50)
      COMMON/DISCON/SSS(200),SFIRST
      COMMON /BMAT/ BL(1)
      PARAMETER(ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0)
      PARAMETER(THREE=3.0D0,FOUR=4.0D0,FIVE=5.0D0,TEN=10.0D0)
      PARAMETER(ANG=1.889726664D0,GDYN=8.2388575D0)
      PARAMETER(FSUMT=1.0D-7,CRTTOL=1.0D-12,MXITER=40)
      PARAMETER(RADC=57.29577951308D0)
      DATA (WK(I),I=1,6)/'CARD','FMAT','DISP','OLDF','PUNC','PRIN'/
      DATA (WK(I),I=7,12)/'ANGS','FIXC','FINT','FLT1','EXFI','BOHR'/
      DATA (WK(I),I=13,18)/'FLT2','GDYN','OLDH','DEGR','DFLP','BFGS'/
      DATA (WK(I),I=19,23)/'POWL','EIGF','DIIS','MODE','DUMB'/
      DATA (WK(I),I=24,27)/'NOUP','NOEX','MURT','NTRY'/
      DATA XNIST/2HN=/
      REWIND INP2
      BILION = 1000000000.0D0
      SFIRST=.TRUE.
      DISP=.TRUE.
      IIN=INP
      IIOUT=IOUT
      IINP2=INP2
      IIPUN=IPUN
      IIGMUP=IGMUP
      IOP(5)=1
      ISCALE=0
      NMAX=napdm
      NCMAX=3*NMAX
      NQMAX=nsym
C
C     DIMENSION OF THE B MATRIX MUST CONFORM TO B(NCMAX,NQMAX)
C.................................................................
C     READ IN INPUT OPTIONS                                      .
C.................................................................
      lopt(1)=1
      lopt(3)=ndisp
      lopt(8)=nfix
      lopt(23)=ndum
Crawdad Adding input option for reading geometry from geom.dat
      iop(50) = 0
      errcod = frdboo('READ_GEOM',iop(50))
Crawdad is iop(50) available?  I think so.

      lopt(2)=1
      errcod = frdboo('FMAT',lopt(2))

      lopt(4)=1
      errcod = frdboo('WRT_BMAT',lopt(5))

      lopt(6)=1
      errcod = frdboo('ANGSTROMS',lopt(7))
      if(lopt(7).gt.0) 
     +  write(iout,*) ' ANGSTROM units will be used'

      errcod = fcnt('INT_FORCES',lopt(9))
      if(lopt(9).gt.0) 
     +  write(iout,*) ' Internal forces will be read from input'

      if(ndisp.gt.0) then
        errcod = frdboo('BOHR',lopt(12))
        if(lopt(12).gt.0) 
     +    write(iout,*) ' BOHR units will be used in displacements'
        errcod = frdboo('DEGR',lopt(16))
        if(lopt(16).gt.0) 
     +    write(iout,*) ' Degree units will be used in displacements'
      endif

      errcod = frdboo('GDYN',lopt(14))
      if(lopt(14).gt.0) 
     +  write(iout,*) ' Forces are assumed to be in mdyn units'

      updt=' '
      errcod = frdc('UPDATE',updt)
      if(updt.eq.'DFLP') then
        lopt(17)=1
        write(iout,*) ' DFLP update will be used'
      else if(updt.eq.'BFGS') then
        lopt(18)=1
        write(iout,*) ' BFGS update will be used'
      else if(updt.eq.'POWL') then
        lopt(19)=1
        write(iout,*) ' POWL update will be used'
      else if(updt.eq.'MURT') then
        lopt(26)=1
        write(iout,*) ' POWL update will be used'
      endif

c     errcod = frdboo('DIIS',lopt(21))

      errcod = frdboo('NOUP',lopt(24))
      if(lopt(24).gt.0) 
     +  write(iout,*) ' Hessian will not be updated'

c     errcod = frdboo('NOEX',lopt(25))

      errcod = frdi('NTRY',lopt(27))

      if(nfix.gt.0) then
        errcod = frdia1('FIX_COORD',ifix,nfix)
        write(iout,*) ' FIX_COORD array',ifix
      endif
        
      nflt1=0
      errcod = fcnt('FLT1',nflt1)
      if(nflt1.ne.0) then
        errcod = frdfa1('FLT1',eta,nflt1)
        lopt(10)=1
        write(iout,*) ' FLT1 array',eta
      endif
       
      nflt2=0
      errcod = fcnt('FLT2',nflt2)
      if(nflt2.ne.0) then
        errcod = frdfa1('FLT2',eta,nflt)
        errcod = frdfa1('FLT2E',energy,nflt)
        lopt(10)=1
        write(iout,*) ' FLT2 array',eta
        write(iout,*) ' FLT2E array',energy
      endif
       
      errcod = frdboo('EIGF',lopt(20))
      if(lopt(20).ne.0) then
        write(iout,*) ' Eigenvector following will be used'
        errcod = frdboo('MINIMUM',iop(5))
        if(iop(5).eq.0) 
     +   write(6,*) ' Transition state search'
        errcod = frdi('MAX_STEP',iop(8))
        errcod = frdi('UPDT',iop(13))
        errcod = frdi('EIG_MAX',iop(16))
        errcod = frdi('EIG_MIN',iop(17))
        errcod = frdboo('NWT_RPH',iop(19))
        errcod = frdboo('NO_PRINT',iop(33))
        errcod = frdboo('DEBUG',iop(34))

        lopt(22)=0
        errcod = frdboo('MODEF',lopt(22))
        if(lopt(22).ne.0) then
          mode1=0
          mftyp='HESSIAN'
          errcod = frdc('MF_METHOD',mftyp)
          if(mftyp.eq.'HESSIAN') mode1=0
          if(mftyp.eq.'INTCO  ') mode1=1

          mode2=1
          errcod = frdi('MF_MODE',mode2)
        endif
      endif

      errcod = frdca1('SYMBOLS',symb,na)
      errcod = frdia1('AT_NUMBER',ia,na)
      write(6,*) ' '
        
C..................................................................
C     CHECK THE NUMBER OF ITERATIONS                              .
C..................................................................
   40 CONTINUE
      read(inp2,2070,end=979) ngeo,numit
      goto 978
  979 write(inp2,2070) 0,0
      rewind inp2
      goto 40
  978 continue
C..................................................................
C     SET THE DEFAULTS HERE                                       .
C..................................................................
      IF(LOPT(3).GT.0.OR.LOPT(10).GT.0.OR.LOPT(13).GT.0) DISP=.FALSE.
      IF(LOPT(10).GT.0) THEN
         LOPT(15)=1
         IF(LOPT(24).NE.0) LOPT(15)=0
      ELSE
         LOPT(4)=1
         IF(LOPT(24).NE.0) LOPT(4)=0
      ENDIF
      IF(LOPT(25).GT.0) THEN
         LOPT(11)=0
      ELSE
         LOPT(11)=1
      ENDIF
      IF(LOPT(17).EQ.0.AND.LOPT(18).EQ.0.AND.LOPT(19).EQ.0.AND.
     &   LOPT(26).EQ.0) THEN
         IF(IOP(5).EQ.1) LOPT(18)=1
         IF(IOP(5).EQ.0) LOPT(19)=1
      ENDIF
C..................................................................
C     TURN OFF HESSIAN UPDATING IN OPTEFC UNLESS BFGS WITH        .
C     POSITIVE DEFINITENESS IS GUARANTEED.                        .
C..................................................................
      IF(LOPT(18).EQ.1 .AND. IOP(13).EQ.2) THEN
         LOPT(24) = 1
      ELSE
         IOP(13) = -1
      ENDIF
      IF(LOPT(11).EQ.0.OR.NUMIT.LT.1) GO TO 60
C..................................................................
C     CHECK THE NUMBER OF ATOMS IN INP2 IS EQUAL NA (EXFI OPTION) .
C..................................................................
      DO 50 I=1,NA
      READ(INP2,2080) W1
      IF(W1.NE.'N=') THEN
         WRITE(IOUT,2090) W1
         RETURN
      END IF
   50 CONTINUE
   60 CONTINUE
C..................................................................
C     read geometry in from file11.dat
c     ffile should position us at the bottom of file11
C..................................................................
Crawdad Altering read from file11 or geom.dat
      if(iop(50).eq.0) then
        do 101 i=1,2*na+1
  101   backspace 11
      endif
      WRITE( IOUT,2110)
      errcod = frdfa1('MASSES',xy,na)
      if(iop(50).gt.0) then
Crawdad Do I need this?
          errcod = fwkadd(8,'DUMMY')
          errcod = frdfa2(':GEOMETRY',gtmp,3,na)  
          if (errcod.ne.EOK) then
            write(6,*) 'ERROR: problem reading geometry from geom.dat'
            call pabort
            call qabort
          endif
          ij=1
          do 70 j=1,na
              ij=((j-1)*3)+1
              xa(j)=gtmp(ij)
              ya(j)=gtmp(ij+1)
              za(j)=gtmp(ij+2)
 70        continue
      endif
      DO 90 I=1,NA
        if(iop(50).eq.0) then
          read(11,2141) xa(i),ya(i),za(i)
        endif
        i3=i*3
        xm(i3)=xy(i)
C..................................................................
C     DEFAULT MASS IS ONE                                         .
C..................................................................
        if(xm(i3).eq.zero) xm(i3)=one
        if (lopt(7).eq.1) go to 80
C..................................................................
C     CONVERT BOHR TO ANGSTROMS                                   .
C..................................................................
        XA(I)=XA(I)/ANG
        YA(I)=YA(I)/ANG
        ZA(I)=ZA(I)/ANG
   80   CONTINUE
C..................................................................
C     CONVERT ANGSTROMS TO BOHR FOR PRINTING                      .
C..................................................................
        XXWA=XA(I)*ANG
        YYWA=YA(I)*ANG
        ZZWA=ZA(I)*ANG
        write(iout,2130) i,ia(i),xxwa,yywa,zzwa,xm(i3)
C..................................................................
C     PUT INVERSE MASS IN XM                                      .
C..................................................................
        xm(i3)=one/xm(i3)
        xm(i3-1)=xm(i3)
        xm(i3-2)=xm(i3)
   90 CONTINUE
C..................................................................
      n3n=3*na
      nad3=3*napdm
      IF (LOPT(3).GT.0) GO TO 130
      IF (LOPT(9).GT.0) GO TO 190
C..................................................................
C     read forces in from file11.dat                              .
C..................................................................
      do 102 i=1,nad3
  102 f(i)=0.0d0
      read (11,2141) (f(i),i=1,n3n)
      if(lopt(14).eq.1) go to 110
 2141 format (20x,3f20.10)
C..................................................................
C     CONVERT GRADIENTS IN ATOMIC UNITS TO FORCES IN MDYN         .
C..................................................................
      do 100 i=1,n3n
         f(i)=-f(i)*gdyn
  100 continue
  110 CONTINUE
      WRITE(IOUT,2150)
      WRITE(IOUT,2160) (F(I),I=1,n3n)
C..................................................................
C     CHECK THAT FORCES ADD TO ZERO                               .
C..................................................................
      XFQX=F(1)
      DO 120 I=2,n3n
  120 XFQX=XFQX+F(I)
      IF(DABS(XFQX).GT.FSUMT) WRITE (IOUT,2170) XFQX
      IF(LOPT(3).LE.0) GO TO 190
C..................................................................
C     READ IN DISPLACEMENTS                                       .
C..................................................................
  130 CONTINUE
      call rddisp(indi,ccin,ndisp)
      IF(LOPT(12).NE.1) GO TO 190
C..................................................................
C     CHANGE DISPLACEMENTS TO BOHR INSTEAD OF ANGSTROMS           .
C..................................................................
      nbohr=0
      do 150 i=1,nsym
        k=symms(i,1)
        if(stype(k).eq.'STRE') then
          nbohr=nbohr+1
          ibohr(nbohr)=i
        endif
  150 continue
      DO 180  K=1,nbohr
        DO 180  I=1,ndisp
          DO 180  J=1,4
            IF(INDI(J,I).EQ.IBOHR(K))THEN
              CCIN(J,I)=CCIN(J,I)/ANG
            ELSE IF (LOPT(16).EQ.1) THEN
              CCIN(J,I)=CCIN(J,I)/RADC
            END IF
  180 CONTINUE
C....................................................................
C     PUT DUMMY ATOM COORDINATES IN XA,YA,ZA (NOT COMPLETELY TESTED).
C....................................................................
  190 CONTINUE
      if(ndum.gt.0) then
        errcod = frdfa2('DUMMY_ATOMS',geom,3,ndum)
        if(errcod.ne.EOK) then
          write(6,*)'trouble reading in DUMMY_ATOMS'
          call pabort
          call qabort
        endif
      endif
      do 200 i=1,ndum
        i3=3*(na+i)
        dunorm = dsqrt(geom(1,i)**2 + geom(2,i)**2 + geom(3,i)**2)
        xa(na+i)=geom(1,i)*bilion/dunorm
        ya(na+i)=geom(2,i)*bilion/dunorm
        za(na+i)=geom(3,i)*bilion/dunorm
        xm(i3-2)=1.0d0
        xm(i3-1)=1.0d0
        xm(i3)=1.0d0
  200 continue
CTPH  ENDIF
      WRIT=.FALSE.
      IF (LOPT(6).EQ.1) WRIT=.TRUE.
C....................................................................
C     CHECK THE TYPE OF INTERNAL COORDINATE FOR OPTEFC CONVERSION   .
C     OF ANGSTROMS TO BOHRS                                         .
C....................................................................
      do 205 i=1,nsym
        k=symms(i,1)
        ifstre(i)=stype(k).eq.'STRE'
  205 continue
C....................................................................
C     MAKE B MATRIX                                                 .
C....................................................................
      call machb(b,stype,satom,symms,symmc,nsim,nsym,nstre,nbnd,
     1  ntor,nout,nlin1,nlin2,ninv,nad3,xa,ya,za,qq,napdm,iout,writ,
     2  .false.)
      IF(LOPT(5).EQ.0) GO TO 220
C....................................................................
C     PUNCH B MATRIX                                                .
C....................................................................
      WRITE(IPUN,2200)
      WRITE (IPUN,2210)
      DO 210 I=1,nsym
          WRITE (IPUN,2220) I
          WRITE (IPUN,2230) (B(J,I),J=1,nad3)
  210 CONTINUE
  220 CONTINUE
C....................................................................
C     INVERT B*M*B+  (SEE PULAY'S ARTICLE IN VOL 4, P 166)          .
C....................................................................
      I1=0
      DO 240 I=1,nsym
      DO 240 J=1,nsym
      I1=I1+1
      S=ZERO
      DO 230 K=1,nad3
  230 S=S+B(K,I)*B(K,J)*XM(K)
  240 C(I1)=S
      TOL=1.0D-8
      CALL OSINV1 (C,nsym,D,TOL,L,M)
C....................................................................
      WRITE (IOUT,2240) D
      IF (LOPT(3).GT.0) GO TO 1150
C....................................................................
C     FORM B*M*F                                                    .
C....................................................................
      WRITE( IOUT,2250)
      DO 260 I=1,nsym
         T=ZERO
         DO 250 J=1,nad3
  250      T=T+B(J,I)*F(J)*XM(J)
  260 CC(I)=T
C....................................................................
C     FORM (B*M*B+)-1 * B*M*F = INTERNAL FORCES                     .
C....................................................................
      if(lopt(9).gt.0) errcod = frdfa1('INT_FORCES',fi,lopt(9))
      ij=0
      do 300 i=1,nsym
        t=zero
        if(lopt(9).le.0) then
          do 270 j=1,nsym
            IJ=IJ+1
            T=T+C(IJ)*CC(J)
  270     continue
          fi(i)=t
        endif
C....................................................................
C     WRITE OUT COORDINATE NUMBER, INTERNAL COORDINATE AND FORCE    .
C....................................................................
        if(qq(i).gt.900000000.0d0) then
          write (iout,2265) i,qq(i),fi(i)
        else
          write (iout,2260) i,qq(i),fi(i)
        endif
        qqq(i)=qq(i)
        fii(i)=fi(i)
  300 continue
C....................................................................
C     PUNCH INTERNAL FORCES                                         .
C....................................................................
      WRITE(IPUN,2270) (FI(I),I=1,nsym)
      IF(LOPT(10).GE.1.AND.LOPT(13).LT.1) WRITE(IOUT,2280)
      IF(LOPT(13).GE.1) WRITE(IOUT,2290)
      IF (LOPT(2).LT.1) GO TO 1530
C....................................................................
C     read internal force constant matrix from fconst.dat or inp2   .
C....................................................................
      IF(LOPT(10).GE.1) THEN
           WRITE(IOUT,2300)
      ELSE
           WRITE(IOUT,2310)
      END IF
      IF(LOPT(11).GT.0.AND.NUMIT.GT.0) then
        IX=INP2
      else
        ix=9
        call ffile(ix,'fconst',0)
        rewind(ix)
      endif
      DO 320 I=1,nsym
        read(ix,2270)(f(j),j=1,i)
        write (iout,2320) (f(j),j=1,i)
C....................................................................
C     PUT FORCE CONSTANT IN SQUARE MATRIX C, TAKE THE INVERSE       .
C....................................................................
        DO 320 J=1,I
          IJ=(I-1)*nsym+J
          JI=(J-1)*nsym+I
          C(JI)=F(J)
          C(IJ)=C(JI)
  320 CONTINUE
      CALL OSINV1 (C,nsym,D,1.0D-8,L,M)
      IF(NUMIT.LT.1) GO TO 870
C....................................................................
C     READ IN OLD INTERNAL COORDINATES, FORCES, AND DISPLACEMENTS   .
C     FROM INP2 OR INP                                              .
CTPH  INP2 SHOULD HAVE THIS WHETHER EXFI IS USED OR NOT             .
C....................................................................
      ix=inp2
      DO 340 I=1,nsym
         READ(IX,2340) CC(I),F(I),GEC(I)
  340 CONTINUE
cets071891      MODE= 0
      IF(NUMIT.GT.0.AND.LOPT(11).GT.0) READ(INP2,2620) MODE
      IX=INP
      IF (LOPT(4).LT.1.AND.LOPT(15).LT.1) GO TO 870
C....................................................................
C     UPDATE THE HESSIAN UNLESS IT IS THE FIRST OPTIMIZATION CYCLE  .
C     ALL UPDATES ARE THE POWELL UPDATE THE INVERSE HESSIAN         .
C....................................................................
      WRITE(IOUT,2330)
CTPH  IF(LOPT(11).GT.0) IX=INP2
      DO 345 I=1,nsym
         IF(QQ(I).GT.900000000.0D0) THEN
           WRITE (IOUT,2265) I,CC(I),F(I)
         ELSE
           WRITE (IOUT,2260) I,CC(I),F(I)
         ENDIF
  345 CONTINUE
CTPH  IX=INP
C..................................................................
      IF(LOPT(15).GE.1) GO TO 820
      IF(LOPT(17) .EQ. 1) GOTO 420
      IF(LOPT(18) .EQ. 1) GOTO 520
      IF(LOPT(19) .EQ. 1) GOTO 620
C....................................................................
C     MURTAGH-SARGENT UPDATE                                        .
C....................................................................
      WRITE(IOUT,2350)
      DO 350 I=1,nsym
      GRC(I)=-FI(I)+F(I)
  350 CONTINUE
C 350 GEC(I)=QQ(I)-CC(I)
      IJ=0
      DO 360 I=1,nsym
      DO 360 J=1,nsym
      IJ=IJ+1
  360 HE(I,J)=C(IJ)
      DO 370 I=1,nsym
      NV(I)=0.D0
      DO 370 J=1,nsym
  370 NV(I)=NV(I)+HE(I,J)*GRC(J)
      DO 380 I=1,nsym
  380 NV(I)=GEC(I)-NV(I)
      DN=0.D0
      DO 390 I=1,nsym
  390 DN=DN+NV(I)*GRC(I)
      DO 400 I=1,nsym
      DO 400 J=1,nsym
  400 HE(I,J)=HE(I,J)+NV(I)*NV(J)/DN
      IJ=0
      DO 410 I=1,nsym
      DO 410 J=1,nsym
      IJ=IJ+1
  410 C(IJ)=HE(J,I)
      GOTO 820
  420 CONTINUE
C....................................................................
C     DAVIDON-FLETCHER-POWELL (DFP) UPDATE                          .
C....................................................................
      WRITE(IOUT,2360)
      DO 430 I=1,nsym
      GRC(I)=-FI(I)+F(I)
  430 CONTINUE
C 430 GEC(I)=QQ(I)-CC(I)
      IJ=0
      DO 440 I=1,nsym
      DO 440 J=1,nsym
      IJ=IJ+1
  440 HE(I,J)=C(IJ)
      DN=0.D0
      DO 450 I=1,nsym
  450 DN=DN+GEC(I)*GRC(I)
      DO 460 I=1,nsym
      DO 460 J=1,nsym
  460 HN(I,J)=GEC(I)*GEC(J)/DN
      DO 470 I=1,nsym
      DV(I)=0.D0
      DO 470 J=1,nsym
  470 DV(I)=DV(I)+HE(I,J)*GRC(J)
      DN=0.D0
      DO 480 I=1,nsym
  480 DN=DN+GRC(I)*DV(I)
      DO 490 I=1,nsym
      N2V(I)=0.D0
      DO 490 J=1,nsym
  490 N2V(I)=N2V(I)+GRC(I)*HE(I,J)
      DO 500 I=1,nsym
      DO 500 J=1,nsym
  500 HE(I,J)=HE(I,J)+HN(I,J)-DV(I)*N2V(J)/DN
      IJ=0
      DO 510 I=1,nsym
      DO 510 J=1,nsym
      IJ=IJ+1
  510 C(IJ)=HE(I,J)
      GOTO 820
  520 CONTINUE
C....................................................................
C     BFGS UPDATE                                                   .
C....................................................................
      WRITE(IOUT,2370)
      DO 530 I=1,nsym
      GRC(I)=-FI(I)+F(I)
  530 CONTINUE
C 530 GEC(I)=QQ(I)-CC(I)
      IJ=0
      DO 540 I=1,nsym
      DO 540 J=1,nsym
      IJ=IJ+1
  540 HE(I,J)=C(IJ)
      DN=0.D0
      DO 550 I=1,nsym
  550 DN=DN+GEC(I)*GRC(I)
      DO 560 I=1,nsym
      NV(I)=0.D0
      DO 560 J=1,nsym
  560 NV(I)=NV(I)+GRC(J)*HE(I,J)
      SU=0.D0
      DO 570 I=1,nsym
  570 SU=SU+NV(I)*GRC(I)
      SU=1.D0+SU/DN
      DO 580 I=1,nsym
      DO 580 J=1,nsym
  580 HN(I,J)=HN(I,J)+SU/DN*GEC(I)*GEC(J)
      DO 590 I=1,nsym
      N2V(I)=0.D0
      DO 590 J=1,nsym
  590 N2V(I)=N2V(I)+GRC(J)*HE(I,J)
      DO 600 I=1,nsym
      DO 600 J=1,nsym
  600 HE(I,J)=HE(I,J)+HN(I,J)-GEC(I)*N2V(J)/DN-NV(I)*GEC(J)/DN
      IJ=0
      DO 610 I=1,nsym
      DO 610 J=1,nsym
      IJ=IJ+1
  610 C(IJ)=HE(I,J)
      GO TO 820
  620 CONTINUE
C......................................................................
C     POWELL UPDATE (SYMMETRIC BROYDEN)                               .
C     INVERT C TO REGAIN THE HESSIAN, UPDATE, THEN RESTORE THE INVERSE.
C......................................................................
      CALL OSINV1 (C,nsym,D,1.0D-8,L,M)
      WRITE(IOUT,2380)
      IJ=0
      DO 630 I=1,nsym
      DO 630 J=1,nsym
      IJ=IJ+1
  630 HE(I,J)=C(IJ)
CTPH  DO 640 I=1,nsym
CTPH    GEC(I)=QQ(I)-CC(I)
C 640 CONTINUE
      DO 660 I=1,nsym
      N2V(I)=0.0D0
         DO 650 J=1,nsym
            N2V(I)=N2V(I)+HE(I,J)*GEC(J)
  650    CONTINUE
      GRC(I)=-FI(I)+F(I)-N2V(I)
  660    CONTINUE
      DN=0.D0
      DO 670 I=1,nsym
  670 DN=DN+GEC(I)*GEC(I)
      SU=0.D0
      DO 680 I=1,nsym
  680 SU=SU+GEC(I)*GRC(I)
      SU=SU/DN
      DABSSU=DABS(SU)
C....................................................................
C     YOU MAY NOT WANT TO DO THE UPDATE IF THE DISPLACEMENT IS SMALL.
C....................................................................
CTPH  IF(DN.LT.0.000001D0) GO TO 700
      DO 690 I=1,nsym
      DO 690 J=1,I
      TEMP=GEC(I)*GRC(J)+GEC(J)*GRC(I)-GEC(I)*SU*GEC(J)
      HE(I,J)=HE(I,J)+TEMP/DN
      HE(J,I)=HE(I,J)
  690 CONTINUE
  700 CONTINUE
      IJ=0
      DO 710 I=1,nsym
      DO 710 J=1,nsym
      IJ=IJ+1
  710 C(IJ)=HE(I,J)
      CALL OSINV1 (C,nsym,D,1.0E-8,L,M)
C....................................................................
C     THE OPTIMALLY CONDITIONED UPDATE SHOULD GO HERE               .
C....................................................................
  820 IF(LOPT(4).GE.1) GO TO 870
C....................................................................
C     UPDATE THE FLETCHER POWELL FORCE CONSTANTS                    .
C....................................................................
      IJ=0
      DO 840  I=1,nsym
      S=ZERO
      DO 830  J=1,nsym
      IJ=IJ+1
  830 S=S-C(IJ)*F(J)
  840 CCC(I)=S
      S=ZERO
      S1=ZERO
      DO 850  I=1,nsym
      S=S-CC(I)*F(I)
  850 S1=S1+CCC(I)*F(I)
      IJ=0
      DO 860  I=1,nsym
      DO 860  J=1,nsym
      IJ=IJ+1
  860 C(IJ)=C(IJ)+CC(I)*CC(J)/S+CCC(I)*CCC(J)/S1
  870 CONTINUE
C....................................................................
C     CALCULATE THE NEW COORDINATES IN ONE OF THREE WAYS:           .
C    (A) EIGENVECTOR FOLLOWING  (B) DIIS  (C) NEWTON-RAPHSON        .
C....................................................................
      IF(LOPT(20).GE.1) THEN
C....................................................................
C    (A)  EIGENVECTOR FOLLOWING, EFC IS AN INTERFACE BETWEEN        .
C         BMATIN6 AND OPTEFC WRITTEN BY BRIAN YATES                 .
C         THE HESSIAN IS NEEDED INSTEAD OF ITS INVERSE              .
C....................................................................
CTPH  MODE = 0
         IF(LOPT(22).EQ.1 .AND. NUMIT.LT.1) THEN
            IF(MODE1.EQ.1) MODE = MODE2
            IF(MODE1.EQ.0) MODE = MODE2 + nsym
         ENDIF
         write(6,*) ' mode = ',mode
         CALL OSINV1 (C,nsym,D,1.0D-8,L,M)
         CALL EFC(NUMIT,nsym,F,QQ,FI,C,MODE,QQ1,IFSTRE,GEC,vmode)
         CALL OSINV1 (C,nsym,D,1.0D-8,L,M)
      ELSE
         IF(LOPT(21).GT.0) NGEO=NGEO+1
         IF(NGEO.GT.1.AND.LOPT(21).GT.0) THEN
C....................................................................
C (B) DIIS WITH THE POSSIBILITY OF UPDATING THE HESSIAN; THEREFORE  .
C     THE ERROR VECTORS MUST BE RECALCULATED FROM THE PREVIOUS AND  .
C     PRESENT GRADIENTS, USING THE NEW HESSIAN  E = (H)-1*G         .
C....................................................................
C....................................................................
C     READ IN PAST COORDINATES AND GRADIENTS INTO QD AND FD         .
C....................................................................
            IF(NGEO.GE.nsym) NGEO = nsym-1
            DO 880 J=1,NGEO-1
               DO 880 I=1,nsym
                  READ(INP2,2390) QD(I,J),FD(I,J)
  880       CONTINUE
            M1=NGEO+1
            MM1=nsym
            MNQ=NMAX
            DO 890 I=1,nsym
               FD(I,NGEO)=FI(I)
               QD(I,NGEO)=QQ(I)
  890       CONTINUE
      CALL GDIIS(MM1,NQMAX,M1,nsym,QD,FD,C,A,NV,N2V,QQ,FI,DV,XLAM,IER)
CTPH.................................................................
C     EXTRAPOLATED COORDINATES AND FORCES IN QQ AND FI, THEN RELAX  .
C     PRINTING OF DIIS MATRIX FOR DEBUGGING, TO BE CHANGED LATER    .
CTPH.................................................................
            WRITE(IOUT,2401)
 2401       FORMAT('INVERSE HESSIAN')
            DO 901 I=1,nsym
               J1=(I-1)*nsym+1
               JN=(I-1)*nsym+I
               WRITE(IOUT,2270) (C(J), J=J1,JN)
  901       CONTINUE
            WRITE(IOUT,2400)
            DO 900 I=1,M1
               J1=(I-1)*M1+1
               JN=(I-1)*M1+I
               WRITE(IOUT,2270) (A(J), J=J1,JN)
  900       CONTINUE
            IF (IER.EQ.0) GO TO 910
CTPH.................................................................
C     ERROR SECTION, TO BE CHANGED LATER TO STOP PROGRAM            .
CTPH.................................................................
            WRITE(IOUT,2405)
            RETURN
  910       CONTINUE
C....................................................................
C     PUT INTERNAL DISPLACEMENTS AND OLD INT. COORDINATES IN QQ1, QQ.
C....................................................................
CC    DO 930 I=1,NQ
CC       QQ1(I)=QQ(I)-QD(I,NGEO)
CC       QQ(I)=QD(I,NGEO)
CC930 CONTINUE
            WRITE(IOUT,2320) (QQ(I),I=1,nsym)
            WRITE(IOUT,2320) (FI(I),I=1,nsym)
            WRITE(IOUT,2410) (DV(I),I=1,NGEO)
            WRITE(IOUT,2420) XLAM
         ENDIF
C....................................................................
C    (C)  NEWTON RAPHSON STEP WHICH IS THE DEFAULT                  .
C....................................................................
         IJ=0
         DO 950 I=1,nsym
            S=ZERO
            DO 940 J=1,nsym
               IJ=IJ+1
               S=S+C(IJ)*FI(J)
  940       CONTINUE
            QQ1(I)=S
  950    CONTINUE
      ENDIF
C....................................................................
C     CALCULATE FLAMBDA FOR FLETCHER PROCEDURE, PRINT NEW INTERNAL  .
C     COORDINATES AND DISPLACEMENTS, PUT DISPLACEMENTS IN GEC FOR   .
C     WRITING TO RESUL2 LATER                                       .
C....................................................................
      IF(LOPT(10).GE.1) THEN
      WRITE(IOUT,2430)
      ELSE
      WRITE(IOUT,2440)
      END IF
      FLAMDA=0.0D0
      DO 960 I=1,nsym
         FLAMDA=FLAMDA+QQ1(I)*FII(I)
         IF(NGEO.GT.1.AND.LOPT(21).GT.0) QQ1(I)=QQ(I)+QQ1(I)-QD(I,NGEO)
         GEC(I)=QQ1(I)
         CC(I)=GEC(I)
  960 CONTINUE
C....................................................................
C     SET CERTAIN DISPLACEMENTS TO ZERO (FIXC OPTION)               .
C....................................................................
      IF(NFIX.GT.0) THEN
        DO 970  I=1,NFIX
          CC(IFIX(I))=ZERO
          QQ1(IFIX(I))=ZERO
          GEC(IFIX(I))=ZERO
  970   CONTINUE
      ENDIF
      DO 975 I=1,nsym
         IF(QQ(I).GT.900000000.0D0) THEN
           WRITE (IOUT,2266) I,QQ1(I),QQ(I)+QQ1(I)
         ELSE
           WRITE (IOUT,2260) I,QQ1(I),QQ(I)+QQ1(I)
         ENDIF
  975 CONTINUE
C....................................................................
C     REGAIN THE HESSIAN FROM THE INVERSE                           .
C....................................................................
  980 CALL OSINV1 (C,nsym,D,1.0E-8,L,M)
      IF(LOPT(4).GT.0.OR.LOPT(15).GT.0) GO TO 1010
C....................................................................
C     CHECK THE FORCE CONSTANTS TO SEE WHICH FORMAT TO USE, PRINT   .
C     IN FILE INP2                                                  .
C....................................................................
      REWIND INP2
      JTRI=nsym*(nsym+1)/2
      IBIG=0
      DO 990 J=1,JTRI
         IF(C(J).LE.-10.0D0.OR.C(J).GE.100.0D0) IBIG=1
  990 CONTINUE
      DO 995 J=1,JTRI
         IF(C(J).LE.-100.0D0.OR.C(J).GE.1000.0D0) IBIG=2
  995 CONTINUE
      IJ=1
      DO 1000 I=1,nsym
         IJ1=IJ+I-1
         IF(IBIG.EQ.2) THEN
            WRITE(INP2,2325)(C(J),J=IJ,IJ1)
         ELSE IF(IBIG.EQ.1) THEN
            WRITE(INP2,2320)(C(J),J=IJ,IJ1)
         ELSE
            WRITE(INP2,2450)(C(J),J=IJ,IJ1)
         ENDIF
         IJ=IJ+nsym
 1000 CONTINUE
C....................................................................
      IF(LOPT(10).GT.0) GO TO 1040
      IF(LOPT(4).LT.1) GO TO 1210
C....................................................................
C     CHECK THE FORCE CONSTANTS TO SEE WHICH FORMAT TO USE, PRINT   .
C     IN FILE IPUN,IOUT,INP2                                        .
C....................................................................
 1010 IF(NUMIT.GT.0) WRITE(IOUT,2460)
      REWIND INP2
      JTRI=nsym*(nsym+1)/2
      IBIG=0
      DO 1020 J=1,JTRI
         IF(C(J).LE.-10.0D0.OR.C(J).GE.100.0D0) IBIG=1
 1020 CONTINUE
      DO 1025 J=1,JTRI
         IF(C(J).LE.-100.0D0.OR.C(J).GE.1000.0D0) IBIG=2
 1025 CONTINUE
      IJ=1
      IF(LOPT(5).NE.0) WRITE(IPUN,2310)
      DO 1030 I=1,nsym
         IJ1=IJ+I-1
         IF(IBIG.EQ.2) THEN
       IF(NUMIT.GE.1.AND.LOPT(24).EQ.0) WRITE(IOUT,2325) (C(J),J=IJ,IJ1)
            WRITE(INP2,2325)(C(J),J=IJ,IJ1)
            IF (LOPT(5).GT.0) WRITE(IPUN,2325) (C(J),J=IJ,IJ1)
         ELSE IF(IBIG.EQ.1) THEN
       IF(NUMIT.GE.1.AND.LOPT(24).EQ.0) WRITE(IOUT,2320) (C(J),J=IJ,IJ1)
            WRITE(INP2,2320)(C(J),J=IJ,IJ1)
            IF (LOPT(5).GT.0) WRITE(IPUN,2320) (C(J),J=IJ,IJ1)
         ELSE
       IF(NUMIT.GE.1.AND.LOPT(24).EQ.0) WRITE(IOUT,2450) (C(J),J=IJ,IJ1)
            WRITE(INP2,2450)(C(J),J=IJ,IJ1)
            IF (LOPT(5).GT.0) WRITE(IPUN,2450) (C(J),J=IJ,IJ1)
         ENDIF
         IJ=IJ+nsym
 1030 CONTINUE
      IF(LOPT(4).GE.1) GO TO 1210
C....................................................................
C     FIRST PART OF THE FLETCHER POWELL PROCEDURE (FLT1)            .
C....................................................................
 1040 DO 1050  I=1,nsym
 1050 CCC(I)=QQ1(I)
      IF(LOPT(13).EQ.1) GO TO 1080
      IFLT=1
C....................................................................
C     CHECK TO SEE IF LAST STEP ALONG FLETCHER-POWELL VECTOR IS MADE.
C....................................................................
 1060 IF(ETA(IFLT).EQ.ZERO) GO TO 1450
      DO 1070  I=1,nsym
      QQ(I)=QQQ(I)
      QQ1(I)=ETA(IFLT)*CCC(I)
 1070 CC(I)=QQ1(I)
      WRITE(IOUT,2470) ETA(IFLT)
      GO TO 1210
C....................................................................
C     SECOND PART OF THE FLETCHER POWELL PROCEDURE (FLT2)           .
C....................................................................
 1080 IJ=1
      IF(ENERGY(3).GE.0.0D0) GO TO 1120
      DO 1090 I=1,3
      CI(IJ)=1
      CI(IJ+1)=ETA(I)
      CI(IJ+2)=ETA(I)**2
 1090 IJ=IJ+3
      CALL OSINV1(CI,3,D,1.0E-8,L,M)
      IJ=0
      DO 1110 I=1,3
      S=ZERO
      DO 1100 J=1,3
      IJ=IJ+1
 1100 S=S+CI(IJ)*ENERGY(J)
 1110 AE(I)=S
      ETA(4)=-AE(2)/(TWO*AE(3))
      ENERGY(4)=ETA(4)*(ETA(4)*AE(3)+AE(2)) +AE(1)
      GO TO 1130
 1120 WRITE(IOUT,*) 'FLAMDA=',FLAMDA
      FLAMDA=FLAMDA*ANG/GDYN
      ELAMDA=(ENERGY(1)-ENERGY(2))/ETA(1)+FLAMDA
      ETA(4)=FLAMDA*ETA(1)/(2.0D0*ELAMDA)
      ENERGY(4)=ENERGY(2)-FLAMDA*ETA(4)/2.0D0
 1130 IFLT=4
      DO 1140 I=1,nsym
      QQ(I)=QQQ(I)
      QQ1(I)=ETA(IFLT)*CCC(I)
 1140 CC(I)=QQ1(I)
      WRITE(IOUT,2480) ENERGY(4)
      WRITE(IOUT,2490) ETA(4)
      GO TO 1210
 1150 CONTINUE
C....................................................................
C     WRITE OUT ORIGINAL INTERNAL COORDINATES, DISPLACEMENTS        .
C     FOR ONE DISTORTED GEOMETRY AT A TIME (DISP OPTION)            .
C     THE PROGRAM LOOPS BACK TO 1160 UNTIL II = # DISPLACEMENTS     .
C....................................................................
      WRITE(IOUT,2500)
      WRITE(IOUT,2510) (QQ(J),J=1,nsym)
      II=0
 1160 II=II+1
         DO 1170 J=1,4
            INDT(J)=INDI(J,II)
 1170       CIN(J)=CCIN(J,II)
      WRITE(IOUT,2200)
      WRITE(IOUT,2520) (INDT(J),CIN(J),J=1,4)
      DO 1180 J=1,nsym
         CC(J)=ZERO
 1180 CONTINUE
      DO 1190 J=1,4
      J1=INDT(J)
        IF (J1.LE.0) GO TO 1190
      CC(J1) = CIN(J)
 1190 CONTINUE
      DO 1200 J=1,nsym
 1200   QQ1(J)=CC(J)
      GO TO 1250
C....................................................................
C     MAKE B MATRIX, THEN (B*M*B+)-1  IN C                          .
C....................................................................
 1210 CONTINUE
 1220 CONTINUE
      call machb(b,stype,satom,symms,symmc,nsim,nsym,nstre,nbnd,
     1  ntor,nout,nlin1,nlin2,ninv,nad3,xa,ya,za,qq,napdm,iout,.false.,
     2  .false.)
      I1=0
      DO 1240 I=1,nsym
      DO 1240 J=1,nsym
      I1=I1+1
      S=ZERO
      DO 1230 K=1,nad3
 1230 S=S+B(K,I)*B(K,J)*XM(K)
 1240 C(I1)=S
      CALL OSINV1(C,nsym,D,TOL,L,M)
 1250 CONTINUE
C....................................................................
C     NOW TRANSFORM THE INTERNAL DISPLACEMENTS TO CARTESIAN         .
C     DISPLACEMENTS BY DELTA X = (B*M*B+)-1 * B*M * (DELTA Q)       .
C     THIS IS THE SAME TRANSFORMATION USED TO CONVERT CARTESIAN     .
C     FORCES TO INTERNAL FORCES, AND CAN BE USED BECAUSE THE        .
C     DERIVATIVE WITH RESPECT TO A COORDINATE AND A COORDINATE ARE  .
C     CONTRAVARIANT                                                 .
C....................................................................
         ISCALE=ISCALE+1
         IREP=0
         I3=0
         DO 1260 I=1,napdm
            XY(I)=XA(I)
         J=I+napdm
         XY(J)=YA(I)
         J=J+napdm
         XY(J)=ZA(I)
 1260    CONTINUE
C...........................................................
 1270    CONTINUE
         IREP=IREP+1
         IJ=0
         DO 1290 I=1,nsym
            S=ZERO
         DO 1280 J=1,nsym
               IJ=IJ+1
               S=S+C(IJ)*QQ1(J)
 1280 CONTINUE
         FI(I)=S
 1290 CONTINUE
C...........................................................
         DO 1310 I=1,nad3
            S=ZERO
            DO 1300 J=1,nsym
               S=S+B(I,J)*FI(J)
 1300       CONTINUE
         QQ1(I)=S*XM(I)
 1310    CONTINUE
C....................................................................
C     UPDATE THE CARTESIAN COORDINATES                              .
C....................................................................
         I3=0
         DO 1320 I=1,napdm
            I3=I3+3
            XY(I)=XY(I)+QQ1(I3-2)
            J=I+napdm
            XY(J)=XY(J)+QQ1(I3-1)
            J=J+napdm
            XY(J)=XY(J)+QQ1(I3)
 1320   CONTINUE
C....................................................................
C     RECALCULATE THE INTERNALS AT THE NEW CARTESIAN COORDINATES    .
C     USING THE SAME OLD B MATRIX                                   .
C....................................................................
         I1=napdm+1
         I2=2*napdm+1
 1330 CONTINUE
         call machb(b,stype,satom,symms,symmc,nsim,nsym,nstre,nbnd,
     1    ntor,nout,nlin1,nlin2,ninv,nad3,xy,xy(i1),xy(i2),qq1,napdm,
     2    iout,.false.,.true.)
C....................................................................
C     THE INTERNAL COORDINATES CORRESPONDING TO THE PRESENT         .
C     CARTESIAN COORDINATES ARE NOW IN QQ1                          .
C....................................................................
          WRITE (IOUT,2530) IREP
         WRITE  (IOUT,2510)  (QQ1(I),I=1,nsym)
         QMAX=ZERO
         DO 1340 I=1,nsym
            QQ1(I)=QQ(I)+CC(I)-QQ1(I)
            IF  (QMAX.LT.DABS(QQ1(I)))  QMAX=ABS(QQ1(I))
 1340    CONTINUE
         IF (IREP.LT.MXITER.AND.QMAX.GT.CRTTOL) GO TO 1270
C....................................................................
C     IF QQ1 MATCHES THE DESIRED FINAL SET OF INTERNAL COORDINATES  .
C     GO ON, ELSE DO ANOTHER ITERATIVE CYCLE                        .
C....................................................................
         IF (QMAX.GT.CRTTOL) THEN
            WRITE(IOUT,2540) MXITER
C....................................................................
C     IF DISPLACEMENT FAILS, SCALE AND TRY ONCE MORE                .
C....................................................................
            IF(DISP.AND.ISCALE.LT.LOPT(27)) THEN
               WRITE(IOUT,2545)
               DO 1345 I=1,nsym
                  QQ1(I)=0.5D0*CC(I)
                  CC(I) = QQ1(I)
 1345          CONTINUE
               GO TO 1250
            ENDIF
         ENDIF
         WRITE(IOUT,2550)
 1350    CONTINUE
C....................................................................
C     PRINT THE DISPLACEMENTS IN CARTESIAN COORDINATES              .
C....................................................................
         I3=0
         DO 1360 I=1,NA
            I3=I3+3
            YM=ONE/XM(I3)
            J=I+napdm
            J1=J+napdm
            X=XY(I)-XA(I)
            Y=XY(J)-YA(I)
            Z=XY(J1)-ZA(I)
            WRITE (IOUT,2560) I,IA(I),X,Y,Z,YM
 1360    CONTINUE
C....................................................................
C     PRINT THE NEW CARTESIAN COORDINATES                           .
C....................................................................
         WRITE(IOUT,2570)
      I3=0
      DO 1370 I=1,NA
           J=I+napdm
           J1=J+napdm
           I3=I3+3
           YM=ONE/XM(I3)
           WRITE(IOUT,2560) I,IA(I),XY(I),XY(J),XY(J1),YM
 1370 CONTINUE
C.................................................................
      IF(LOPT(2).EQ.1) GO TO 1390
 1380 CONTINUE
      IF(II.LT.ndisp) GO TO 1420
 1390 REWIND INP2
      IF(LOPT(3).GT.0) GO TO 1410
C....................................................................
C     READ FORCE CONSTANT MATRIX BACK IN FROM INP2                  .
C....................................................................
      JJ=1
      DO 1400 I=1,nsym
      IJ1=JJ + I -1
      READ(INP2,2450)(C(J),J=JJ,IJ1)
      JJ=JJ+nsym
 1400 CONTINUE
C.................................................................
 1410 REWIND INP2
 1420 WRITE(IOUT,2580)
C....................................................................
C     WRITE NEW GEOMETRY IN ATOMIC UNITS TO IGMUP                   .
C....................................................................
      call wrtgeo(xy,geobf,na,nad3,numit,igmup,iout,ang,ii,ndisp,
     +  lbls,nsym,ccin,indi)
 1450 continue
C....................................................................
C     IF NOT THE LAST DISPLACEMENT USING DISP OPTION, DO ANOTHER    .
C....................................................................
      IF(II.LT.ndisp) GO TO 1160
C....................................................................
C     IF NOT THE LAST STEP ALONG FLETCHER-POWELL VECTOR (FLT1) DO   .
C     THE NEXT                                                      .
C....................................................................
      IF(LOPT(10).GE.1.AND.IFLT.LT.4) THEN
           IFLT=IFLT+1
           GO TO 1060
      END IF
      IF(LOPT(3).GT.0) GO TO 1520
C....................................................................
C  WRITE OUT NGEO AND NUMIT                                         .
C....................................................................
      REWIND INP2
      WRITE(INP2,2070) NGEO,NUMIT+1
C....................................................................
C     WRITE CARTESIAN COORDINATES IN ATOMIC UNITS IN INP2           .
C....................................................................
      DO 1460 I=1,NA
      J=I+napdm
      J1=J+napdm
      X=XY(I)*ANG
      Y=XY(J)*ANG
      Z=XY(J1)*ANG
      WRITE(INP2,2610)SYMB(I),IA(I),X,Y,Z
 1460 CONTINUE
C....................................................................
C     WRITE FORCE CONSTANTS IN INTERNAL COORDINATES IN INP2, AFTER  .
C     CHECKING TO SEE WHICH FORMAT TO USE                           .
C....................................................................
      JTRI=nsym*(nsym+1)/2
      IBIG=0
      DO 1470 J=1,JTRI
         IF(C(J).LE.-10.0D0.OR.C(J).GE.100.0D0) IBIG=1
 1470 CONTINUE
      DO 1475 J=1,JTRI
         IF(C(J).LE.-100.0D0.OR.C(J).GE.1000.0D0) IBIG=2
 1475 CONTINUE
      JJ=1
      DO 1480 I=1,nsym
         IJ1=JJ+I-1
         IF(IBIG.EQ.2) THEN
            WRITE(INP2,2325)(C(J),J=JJ,IJ1)
         ELSE IF(IBIG.EQ.1) THEN
            WRITE(INP2,2320)(C(J),J=JJ,IJ1)
         ELSE
            WRITE(INP2,2450)(C(J),J=JJ,IJ1)
         ENDIF
         JJ=JJ+nsym
 1480 CONTINUE
C....................................................................
C     WRITE OUT INTERNAL FORCES, COORDINATES, AND DISPLACEMENTS;    .
C     THEN MODE1
C....................................................................
      DO 1490 I=1,nsym
        IF(QQQ(I).LT.999.0D0) THEN
          WRITE(INP2,2340) QQQ(I),FII(I),CC(I)
        ELSE
          WRITE(INP2,2345) QQQ(I),FII(I),CC(I)
        ENDIF
 1490 CONTINUE
      WRITE(INP2,2620) MODE
C....................................................................
C     WRITE OUT INTERNAL COORDINATES AND FORCES FOR ALL PREVIOUS    .
C     GEOMETRIES, STARTING WITH THE OLDEST FOR GDIIS                .
C....................................................................
      IF(LOPT(21).GT.0) THEN
         IF(NUMIT.GE.1) THEN
            JS=1
            IF(NGEO.GE.nsym-1) JS=2
            DO 1500 J=JS,NGEO
               DO 1500 I=1,nsym
                 IF(QD(I,J).LT.999.0D0) THEN
                   WRITE(INP2,2390) QD(I,J),FD(I,J)
                 ELSE
                   WRITE(INP2,2395) QD(I,J),FD(I,J)
                 ENDIF
 1500       CONTINUE
         ELSE
            DO 1505 I=1,nsym
              IF(QQQ(I).LT.999.0D0) THEN
                WRITE(INP2,2390) QQQ(I),FII(I)
              ELSE
                WRITE(INP2,2395) QQQ(I),FII(I)
              ENDIF
 1505       CONTINUE
         ENDIF
      ENDIF
C....................................................................
C     WRITE OUT THE MODE NUMBER OF THE MODE MOST LIKE THE PREVIOUS  .
C     ONE, THEN THE EIGENMODE OF THE HESSIAN ITSELF                 .
C....................................................................
      IF(LOPT(20).GT.0) THEN
         DO 1510 I=1,nsym
            WRITE(INP2,2630) VMODE(I)
 1510 CONTINUE
      ENDIF
 1520 CONTINUE
      REWIND INP2
 1530 RETURN
C
 2000 FORMAT (A4,1X,I5)
 2010 FORMAT(1X,A4,2X,13HOPTION IS ON ,I5)
 2020 FORMAT(A4,6X,4F6.0)
 2025 FORMAT(1X,A4,2X,13HOPTION IS ON ,4F6.0)
 2030 FORMAT(A4,8X,3(F6.0,F16.10))
 2035 FORMAT(1X,A4,8X,3(F6.3,F16.10))
 2040 FORMAT(A4,1X,9I5)
 2045 FORMAT(1X,A4,2X,13HOPTION IS ON ,9I5)
 2050 FORMAT(A4,1X,2I5)
 2055 FORMAT(1X,A4,2X,13HOPTION IS ON ,2I5)
 2060 FORMAT(A4,12X,3F16.10)
 2065 FORMAT(1X,A4,1X,4H IS ,3(F16.2,2X))
 2070 FORMAT(2I5)
 2080 FORMAT(1X,A2)
 2090 FORMAT(1X,'EXFI READ ERROR--',A2)
 2100 FORMAT (1X,'TOO MANY NUCLEI, EXPECTING ',I5)
 2105 FORMAT (1X,'CARD OPTION MISSING, PLEASE TELL ME HOW MANY ATOMS
     &        THERE ARE!')
 2110 FORMAT (20X,28HNUCLEAR COORDINATES IN BOHRS    ,/)
 2120 FORMAT(1X,A5,2X,I2,6X,3F16.12,F10.5)
 2130 FORMAT (3X,I2,3X,I2,3(F16.12,2X),F12.6)
 2140 FORMAT (3F16.12)
 2150 FORMAT (//,1X,6HFORCES,/)
 2160 FORMAT (1X,3F16.12)
 2170 FORMAT (/,1X,26HFORCES DO NOT VANISH, SUM=,F15.7,/)
 2180 FORMAT (4(I2,2X,F12.8,2X))
 2190 FORMAT (A4,22X,A2)
 2200 FORMAT(///)
 2210 FORMAT (1X,9HB MATRIX ,//)
 2220 FORMAT ('I=',I2)
 2230 FORMAT (5F16.12)
 2240 FORMAT(/,5X,16H DETERMINANT  = , E16.5)
 2250 FORMAT (/,1X,32HINTERNAL COORDINATES AND FORCES ,/)
 2260 FORMAT (2X,I2,4X,2F16.12)
 2265 FORMAT (2X,I2,4X,F16.4,F16.12)
 2266 FORMAT (2X,I2,4X,F16.12,F16.4)
 2270 FORMAT (8F10.7)
 2280 FORMAT(//,1X,'FLETCHER-POWELL OPTIMIZATION  PART 1',/,120('-'))
 2290 FORMAT(//,1X,'FLETCHER-POWELL OPTIMIZATION  PART 2',/,120('-'))
 2300 FORMAT (1X,'FLETCHER POWELL FORCE CONSTANTS',/)
 2310 FORMAT (1X,16HFORCE CONSTANTS ,/)
 2320 FORMAT (8F10.6)
 2325 FORMAT (8F10.5)
 2330 FORMAT (/1X,24HOLD GEOMETRY AND FORCES ,/)
 2340 FORMAT(3F16.12)
 2345 FORMAT(F16.4,2F16.12)
 2350 FORMAT (1X,'MURTAGH-SARGENT HESSIAN UPDATE',/)
 2360 FORMAT (1X,'DAVIDON-FLETCHER-POWELL HESSIAN UPDATE',/)
 2370 FORMAT (1X,'BFGS HESSIAN UPDATE',/)
 2380 FORMAT (1X,'POWELL HESSIAN UPDATE',/)
 2390 FORMAT(2F16.12)
 2395 FORMAT(F16.4,F16.12)
 2400 FORMAT(1X,'DIIS MATRIX')
 2405 FORMAT(/1X,'ERROR IN DIIS',/1X,'PROBABLY GRADIENTS ARE PARALLEL',
     1 /)
 2410 FORMAT(/1X,'DIIS COEFFICIENTS',/1X,1(10F10.4,/))
 2420 FORMAT(/1X,'LAMBDA =',E12.4)
 2430 FORMAT(/1X,'FLETCHER POWELL DIRECTION VECTOR AND NEW INTERNAL COOR
     1DINATES FOR ETA=1',/)
 2440 FORMAT(/1X,44HGEOMETRY CHANGE AND NEW INTERNAL COORDINATES,/)
 2450 FORMAT (8F10.7)
 2460 FORMAT (/,1X,31HIMPROVED FORCE CONSTANT MATRIX ,/)
 2470 FORMAT(80('-'),/,1X,'ETA= ',F14.10,/)
 2480 FORMAT(80('-'),/,1X,'PREDICTED ENERGY = ',F16.10)
 2490 FORMAT(1X,'ETA= ',F14.10,/)
 2500 FORMAT(/1X,30HORIGINAL INTERNAL COORDINATES ,/)
 2510 FORMAT (2X,3F24.12)
 2520 FORMAT(120('-'),/,' INTERNAL DISPLACEMENT',/,4(I2,2X,F12.8,4X),/)
 2530 FORMAT (1X,5HSTEP=,I4)
 2540 FORMAT(/,1X,'*****CAUTION,WITHIN',I3,'STEPS NO CONVERGENCE*****'/)
 2545 FORMAT(/,1X,'***** TRYING TO SCALE BACK DISPLACEMENT BY HALF **'/)
 2550 FORMAT (/,1X,'  CARTESIAN DISPLACEMENTS AND NEW CARTESIAN COORDINA
     1TES IN ANGSTROMS',/)
 2560 FORMAT(1X,I3,4H   N,I4,3F16.12,F12.6)
 2570 FORMAT(//)
 2580 FORMAT(/' CARTESIAN COORDINATES IN A.U.'/)
 2590 FORMAT(3F20.10)
 2600 FORMAT(1X,I3,4H   N,I4,3F16.12)
 2610 FORMAT(1X,'N=',A4,4X,I2,3F16.12)
 2620 FORMAT(I5)
 2630 FORMAT(8F10.4)
      END
