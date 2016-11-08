      subroutine readmo(cmo,eigocc,msfru,
     &                 su,usotao,eig,mxnao,mxntri)
C     This subroutine reads MOs from file30, and change the occupied
C  MOs into the form in terms of AOs from SOs.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 SU(mxntri),USOTAO(mxnao,mxnao),EIG(mxnao)
      INTEGER ITEMP(200),NLAMDA(8),NC(8),NO(8)
      logical change
      common/funcs/natom,n3n,nbfao,nbatri,nocc,mb(20)
      real*8 cmo(mxnao,mxnao)
      real*8 eigocc(mxnao)
      integer msfru
C
      ITAP30=30
      CALL SREW(ITAP30)
      IEND=1
      MPOINT=200
      MCONST=200
      MCALCS=100
      DO 1 I=1,MPOINT
         ITEMP(I)=0
    1 CONTINUE
C
      CALL WREADW (ITAP30,ITEMP,MCONST,101,JUNK)
      MSFRU = ITEMP(18)
      NAT = ITEMP(19)
      NUM = ITEMP(22)
      NSYMHF = ITEMP(41)
      MXCOEF = ITEMP(42)
      IOPEN = ITEMP(43)
C
C     ----- array for pointers -----      (always from pointer 301)
      CALL WREADW (ITAP30,ITEMP,MPOINT,301,IEND)
      ITMP29 = ITEMP(29)
C    itmp29 is the pointer to AO-SO transformation matrix.
C
C     ----- pointers to calculations -----  (always from pointer 501)
      CALL WREADW (ITAP30,ITEMP,MCALCS,501,IEND)
      ITMP1=ITEMP(1)
C     itmp1 is the pointer to directory for first calculation.
C
      IF (ITMP29.EQ.0) GO TO 99
      CALL WREADW(ITAP30,USOTAO,INTOWP(NUM*MSFRU),ITMP29,JUNK)
C     skipped 60 points, which are for title.
      ITMP61=ITMP1+60
C     ----- zero pointers ot vectors etc. -----
      CALL WREADW (ITAP30,ITEMP,20,ITMP61,IEND)
C   ITEMP(1) is the pointer to SCF vectors, eigenvalues, etc.
      IF (ITEMP(1).EQ.0) GO TO 99
C
C  SU: symmetry blocked SCF vector;  EIG: orbital energies.
C  NLAMDA: orbital # per irrep;  NC: closed shell orbital # per irrep.
      IEND=IEND+INTOWP(3*NAT)
      CALL WREADW(ITAP30,SU,INTOWP(MXCOEF),ITEMP(1),IEND)
      CALL WREADW(ITAP30,EIG,INTOWP(MSFRU),IEND,IEND)
      iend=iend+nsymhf
      CALL WREADW(ITAP30,NLAMDA,NSYMHF,IEND,IEND)
      CALL WREADW(ITAP30,NC,NSYMHF,IEND,IEND)
      IF(IOPEN.NE.0) CALL WREADW(ITAP30,NO,NSYMHF,IEND,IEND)
C
      nocc=0
      nscc=0
      nmo=0
      do 110 i=1,nsymhf
  110 nocc=nocc+nc(i)
      do 120 i=1,nsymhf
  120 nscc=nscc+no(i)
      do 130 i=1,nsymhf
  130 nmo=nmo+nlamda(i)
      call cmoao(nsymhf,num,msfru,nc,no,nlamda,usotao,su,cmo,mxnao,
     &  iopen)
c
c    pick up the eigenvalues of the OCCUPIED MOs.
      do 150 nblc=1,nsymhf
      i1=1
      if (nblc.ne.1) i1=i2+1
      i2=i1-1+nc(nblc)
      k1=1
      if (nblc.ne.1) k1=k2+1
      k2=k1-1+nlamda(nblc)
      k=k1-1
      do 140 i=i1,i2
      k=k+1
      eigocc(i)=eig(k)
  140 continue
  150 continue
c
      if(iopen.eq.0) goto 159
c    pick up the eigenvalues of the SINGLE OCCUPIED MOs.
      k2=0
      do 151 nblc=1,nsymhf
      i1=i2+1
      i2=i1-1+no(nblc)
      k=k2+nc(nblc)
      k2=k2+nlamda(nblc)
      WRITE(6,*) 'SINGLE',I1,I2,K,K2
      if (i2.lt.i1) goto 151
      do 141 i=i1,i2
      k=k+1
      eigocc(i)=eig(k)
  141 continue
  151 continue
c
  159 continue
c
c    pick up the eigenvalues of the VIRTUAL OCCUPIED MOs.
      k2=0
      do 152 nblc=1,nsymhf
      i1=i2+1
      i2=i1-1+nlamda(nblc)-nc(nblc)-no(nblc)
      k=k2+nc(nblc)+no(nblc)
      WRITE(6,*) 'VIRTUAL',I1,I2,K,K2
      k2=k2+nlamda(nblc)
      do 142 i=i1,i2
      k=k+1
      eigocc(i)=eig(k)
  142 continue
  152 continue
c
c    reorder the OCCUPIED MOs.
      do 180 i=1,nocc-1
      change=.false.
      eigen=eigocc(i)
      do 160 j=i+1,nocc
      if (eigocc(j).lt.eigocc(i)) then
      eigocc(i)=eigocc(j)
      jj=j
      change=.true.
      endif
  160 continue
      if (change) then
      eigocc(jj)=eigen
      do 170 k=1,num
      a=cmo(i,k)
      cmo(i,k)=cmo(jj,k)
      cmo(jj,k)=a
  170 continue
      endif
  180 continue
c
c    reorder the VIRTUAL MOs.
      change=.false.
      nto=nscc+nocc+1
      do 181 i=nto,nmo-1
      eigen=eigocc(i)
      do 161 j=i+1,nmo
      if (eigocc(j).lt.eigocc(i)) then
      eigocc(i)=eigocc(j)
      jj=j
      change=.true.
      endif
  161 continue
      if (change) then
      eigocc(jj)=eigen
      do 171 k=1,num
      a=cmo(i,k)
      cmo(i,k)=cmo(jj,k)
      cmo(jj,k)=a
  171 continue
      endif
  181 continue
C
      GO TO 100
C
   99 WRITE (6,*) '          '
      WRITE (6,*) 'FAILED!  IT IS SUPPOSED TO RUN SCF BEFORE THIS JOB.'
  100 RETURN
c
      END
