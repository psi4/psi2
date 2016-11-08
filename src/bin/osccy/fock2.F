      SUBROUTINE FOCK2 (NT,HH,FOCK,EREF,ENUC,EIG,
     .                  IPQ,ITR,ITV,NO,NV,xIN,
     .                  iOFF,iADD,
     .                  ORBSYM,FLOV,NIRRED,I79,FZO,FZV,
     .                  ntmo,spin,map,cctscf,nfzo,nfzv)
      IMPLICIT INTEGER (A-Z)
      REAL*8 HH(*),FOCK(*),xIN(*),EIG(*)
      REAL*8 EREF,ENUC
      INTEGER IPQ(*),ITR(*),ITV(*),FZO(no),FZV(nv),
     .        iOFF(*),iADD(*),nfzo(*),nfzv(*),
     .        ORBSYM(nt),FLOV(8,2,4),
     .        spin(nt),map(nt),cctscf(*)
C
      NT2=NT*(NT+1)/2
      ntmo2=ntmo*(ntmo+1)/2
c     write(6,*)' ENTERING FOCK2'
c     write(6,*)' 100 integrals'
c     write(6,*)(xin(jkh),jkh=1,100)
c     write(6,*)'orbsym'
c     write(6,*) orbsym 
c     write(6,*)' spin' 
c     write(6,*)  spin  
c     write(6,*)' map' 
c     write(6,*)  map  
c     write(6,*)' fzo' 
c     write(6,*)  fzo  
c     write(6,*)' fzv' 
c     write(6,*)  fzv  
C
C >>> READ ONE-E MO INTS FROM FILE79
C
      ITAP79=79
COK   CALL RFILE(ITAP79)
      CALL WREADW(ITAP79,HH,INTOWP(NTmo2),1,PON79)
C
c     WRITE(6,*)'ONE E HAMILTONIAN MO BASIS'
c     CALL PRINT(HH,NTmo2,NTmo,6)
C
C >>> ADD FROZEN ORBITALS CONTRIBUTION OBTAINED IN RDINTS
C
C
C     WRITE(6,*)'FROZEN MOS CONTRIBUTION '
C     WRITE(6,*)'  '
C     CALL PRINT(FOCK,NT2,NT,6)
C
c     DO 220 IJ=1,NT2
c        HH(IJ)=HH(IJ)+FOCK(IJ)
c 220 CONTINUE
C
C >>> 2 ELECTRON CONTRIBUTION TO FOCK MATRIX
C
      CALL ZERO (FOCK,NT2)
c
      DO 248  I=1,no
c        IF(FZO(I).EQ.1)GO TO 248
         ISYM=ORBSYM(I)
         ispn=spin(i)
         xi=map(i)
         FM=FLOV(ISYM+1,ispn,1)-nfzo(isym+1)
         LM=FLOV(ISYM+1,ispn,2)
         IF(LM.GT.I)LM=I
         DO 246  M=FM,LM
            xm=map(m)
            xIM=ipq(MAX0(xI,xm))+MIN0(xI,xm)
            im=ipq(i)+m
            DO 240 J=1,no
c              IF(FZO(J).EQ.1)GO TO 240
               jspn=spin(j)
               jsym=orbsym(j)
               xj=map(j)
               xJJ=ipq(xJ)+xJ
               IMJJ=ioff(MAX0(xIM,xJJ))+iADD(MIN0(xIM,xJJ))
               FOCK(IM)=FOCK(IM)+xIN(IMJJ)
               if(jspn.eq.ispn)then
                  xIJ=ipq(MAX0(xI,xJ))+MIN0(xI,xJ)
                  xMJ=ipq(MAX0(xM,xJ))+MIN0(xM,xJ)
                  IJMJ=iOFF(MAX0(xIJ,xMJ))+iADD(MIN0(xIJ,xMJ))
                  FOCK(IM)=FOCK(IM)-xIN(IJMJ)
               endif
  240       CONTINUE
  246    CONTINUE
  248 CONTINUE
c
      DO 348 I=1,nv
c        IF(FZv(I).EQ.1)GO TO 348
         ISYM=ORBSYM(I+no)
         ispn=spin(i+no)
         xi=map(i+no)
         FM=FLOV(ISYM+1,ispn,1)-nfzo(isym+1)
         LM=FLOV(ISYM+1,ispn,2)
         DO 346  M=FM,LM
            xm=map(m)
            xIM=ipq(MAX0(xI,xm))+MIN0(xI,xm)
            im=ipq(i+no)+m
            DO 340 J=1,no
c              IF(FZO(J).EQ.1)GO TO 340
               jsym=orbsym(j)
               jspn=spin(j)
               xj=map(j)
               xJJ=ipq(xJ)+xJ
               IMJJ=ioff(MAX0(xIM,xJJ))+iADD(MIN0(xIM,xJJ))
               FOCK(IM)=FOCK(IM)+xIN(IMJJ)
               if(jspn.eq.ispn)then
                  xIJ=ipq(MAX0(xI,xJ))+MIN0(xI,xJ)
                  xMJ=ipq(MAX0(xM,xJ))+MIN0(xM,xJ)
                  IJMJ=iOFF(MAX0(xIJ,xMJ))+iADD(MIN0(xIJ,xMJ))
                  FOCK(IM)=FOCK(IM)-xIN(IJMJ)
               endif
  340       CONTINUE
  346    CONTINUE
  348 CONTINUE
c
      DO 448  I=1,nv
c        IF(FZv(I).EQ.1)GO TO 448
         ISYM=ORBSYM(I+no)
         ispn=spin(i+no)
         xi=map(i+no)
         FM=FLOV(ISYM+1,ispn,3)-no
         LM=FLOV(ISYM+1,ispn,4)-no+nfzv(isym+1)
         IF(LM.GT.I)LM=I
         DO 446  M=FM,LM
            xm=map(m+no)
            xIm=ipq(MAX0(xI,xm))+MIN0(xI,xm)
            im=ipq(i+no)+m+no
            DO 440 J=1,no
c              IF(FZO(J).EQ.1)GO TO 440
               jspn=spin(j)
               jsym=orbsym(j)
               xj=map(j)
               xJJ=ipq(xJ)+xJ
               IMJJ=ioff(MAX0(xIM,xJJ))+iADD(MIN0(xIM,xJJ))
               FOCK(IM)=FOCK(IM)+xIN(IMJJ)
               if(jspn.eq.ispn)then
                  xIJ=ipq(MAX0(xI,xJ))+MIN0(xI,xJ)
                  xMJ=ipq(MAX0(xM,xJ))+MIN0(xM,xJ)
                  IJMJ=iOFF(MAX0(xIJ,xMJ))+iADD(MIN0(xIJ,xMJ))
                  FOCK(IM)=FOCK(IM)-xIN(IJMJ)
               endif
  440       CONTINUE
  446    CONTINUE
  448 CONTINUE
c
c     WRITE(6,*)'2e contr. to FOCK MATRIX IN so BASIS'
c     WRITE(6,*)'  '
c     CALL PRINT(FOCK,NT2,NT,6)
c
      DO 221 I=1,NT
         ispn=spin(i)
         xi=map(i)
         do 220 j=1,i
            jspn=spin(j)
            if(ispn.ne.jspn)go to 220
            xj=map(j)
            xij=ipq(max0(xi,xj))+min0(xi,xj)
            ij=ipq(i)+j
            fock(ij)=fock(ij)+hh(xij)
  220    CONTINUE
  221 CONTINUE
c
c     WRITE(6,*)'FOCK MATRIX IN so BASIS'
c     WRITE(6,*)'  '
c     CALL PRINT(FOCK,NT2,NT,6)
C
      EREF=0.0D0
      DO 700 I=1,no
         xi=map(i)
         xii=ipq(xi)+xi
         FII=IPQ(I)+I
         EREF=EREF+FOCK(FII)+hh(xii)
  700 CONTINUE
      EREF=EREF*0.5d0+ENUC
      WRITE(6,*)' '
      WRITE(6,638)EREF
  638 FORMAT(/,' REFERENCE ENERGY  =  ',F20.12)
C
C >>> WRITE FOCK MATRIX IN FILE79 AFTER THE ONE-E MO INTS
C
c     PON79=1+INTOWP(NT2)
c     CALL WWRITW(ITAP79,FOCK,INTOWP(NT2),PON79,PON79)
c     I79=PON79
COK   CALL RCLOSE(ITAP79,3)
C
      RETURN
      END
