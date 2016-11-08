C====================================================================
C  Construct natural orbitals for correlated wave functions
C     Diagonalize the one particle density matrix
c    Last updated on January 22, 2003 by Dr. Yukio Ymaguchi 
C====================================================================
      SUBROUTINE NORBIT(ESO,EAO,ENO,PMO,PNO,OCCNO,
     *                  PROOT,QROOT,RROOT,FV1,FV2,
     *                  PTEMP,QTEMP,RTEMP,
     *                  PEIV,QEIV,REIV,SEIV,ENEIV)
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension eso(nbfso,nbasis),eao(nbfao,nbasis),eno(naoso,nbasis)
      dimension pmo(ntri),pno(ntri)
      dimension occno(nbasis)
      dimension proot(nbasis),qroot(nbasis),rroot(nbasis)
      dimension fv1(nbasis),fv2(nbasis)
      dimension ptemp(ntri),qtemp(ntri),rtemp(ntri)
      dimension peiv(naoso,nbasis),qeiv(naoso,nbasis),reiv(naoso,nbasis)
      dimension seiv(naoso,nbasis),eneiv(naoso*nbasis)
      COMMON/BASIS/NBASIS,NBFAO,NBFSO,NBATRI,NTRI
      COMMON/FUNCS/NATOM,N3N,NSYM
      common/mobas/iaoso,naoso
      COMMON/SIGNS/IOFF(500),IPRNT
      common/symmt/mosym(256),nlamda(10)
      common/orbit/ntype(256),ncore(256),nvirt(256),nactv(256)
      DATA a00 / 0.0D+00 /
      data dlimit / 1.0D-6 /
    1 FORMAT(//,2X,' THE ONE PARTICLE DENSITY MATRIX (MO BASIS)',/)
    2 FORMAT(//,2X,' THE ESO MATRIX',/)
    3 FORMAT(//,2X,' THE EAO MATRIX',/)
    4 FORMAT(//,2X,' EXTRACTED PMO MATRIX FOR ISYM = ',I5,/)
    5 FORMAT(//,2X,' EXTRACTED ESO/EAO MATRIX FOR ISYM = ',I5,/)
    6 FORMAT(//,2X,' EXTRACTED PMO MATRIX FOR ACTIVE MOS',/)
    7 FORMAT(//,2X,' EXTRACTED ESO/EAO MATRIX FOR ACTIVE MOS',/)
    8 FORMAT(//,2X,' EIGENVECTORS FOR ACTIVE MOS',/)
    9 FORMAT(//,2X,' NATURAL ORBITAL MATRIX FOR ISYM = ',I5,/)
   10 FORMAT(//,2X,' REORDERED NATURAL ORBITAL MATRIX FOR ISYM = ',I5,/)
   11 FORMAT(//,2X,' TRANSFERED NATURAL ORBITAL AFTER ISYM = ',I5,/)
   12 FORMAT(//,2X,' THE FINAL NATURAL ORBITAL MATRIX',/)
   13 FORMAT(//,2X,' THE ONE PARTICLE DENSITY MATRIX (NO BASIS)',/)
C                                                                       CPC05250
      if(iprnt.ge.2) then
        WRITE(6,1)
        CALL PRINT(PMO,NTRI,NBASIS,6)
      end if
      if(iaoso.eq.1) then
        if(iprnt.ge.2) then
          write(6,2)
          call matout(eso,nbfso,nbasis,nbfso,nbasis,6)
        end if
      end if
      if(iaoso.eq.2) then
        if(iprnt.ge.2) then
          write(6,3)
          call matout(eao,nbfao,nbasis,nbfao,nbasis,6)
        end if
      end if
      if(iprnt.ge.2) then
        write(6,*) '   '
        write(6,*) ' nlamda = ',(nlamda(i),i=1,nsym)
        write(6,*) '   '
      end if
C
C  Loop over each symmetry
      call zero(eno,naoso*nbasis)
      call zero(pno,ntri)
      call zero(occno,nbasis)
C
      nstot=0
      do 120 ISYM=1,NSYM
      ns=nlamda(ISYM)
      nstri=ioff(ns+1)
      if(iprnt.ge.2) then
        write(6,*) '  ISYM   = ',isym
        write(6,*) '  nbasis = ',nbasis
        write(6,*) '  nbfso  = ',nbfso
        write(6,*) '  nbfao  = ',nbfao
        write(6,*) '  iaoso  = ',iaoso
        write(6,*) '  naoso  = ',naoso
        write(6,*) '  nstot  = ',nstot
        write(6,*) '  ns     = ',ns
        write(6,*) '  nstri = ',nstri
      end if
      if(ns.eq.0) go to 120
c
C  Extract OPDM elements for ISYM symmetry
      call zero(ptemp,nstri)
      call zero(proot,ns)
      do 101 i=1,ns
      ii=nstot+i
      do 101 j=1,i
      jj=nstot+j
      ij=ioff(i)+j
      iijj=ioff(ii)+jj
      ptemp(ij)=pmo(iijj)
  101 continue
      if(iprnt.ge.2) then
        write(6,4) isym
        call print(ptemp,nstri,ns,6)
      end if
      do 102 i=1,ns
      ii=ioff(i+1)
      proot(i)=ptemp(ii)
  102 continue
c
c  Extract ESO/EAO elements for ISYM symmetry
      call zero(peiv,naoso*nbasis)
      if(iaoso.eq.1) then
        do 103 i=1,naoso
        do 103 j=1,ns
        jj=nstot+j
        peiv(i,j)=eso(i,jj)
  103   continue
      end if
      if(iaoso.eq.2) then
        do 104 i=1,naoso
        do 104 j=1,ns
        jj=nstot+j
        peiv(i,j)=eao(i,jj)
  104   continue
      end if
      if(iprnt.ge.2) then
        write(6,5) isym
        call matout(peiv,naoso,nbasis,naoso,ns,6)
      end if
c
c  Classify frozen core, deleted virtual, and active orbitals
      call izero(ncore,ns)
      call izero(nvirt,ns)
      call izero(nactv,ns)
      call izero(ntype,ns)
      nact=0
      ncor=0
      nvir=0
      do 105 i=1,ns
      valu=proot(i)
      if(dabs(valu-2.0d+00).lt.dlimit) then
        ncor=ncor+1
        ncore(ncor)=i
        ntype(i)=0
        go to 105
      end if
      if(dabs(valu).lt.dlimit) then
        nvir=nvir+1
        nvirt(nvir)=i
        ntype(i)=0
        go to 105
      end if
      nact=nact+1
      nactv(nact)=i
      ntype(i)=1
  105 continue
      nactri=ioff(nact+1)
      if(iprnt.ge.2) then
        write(6,*) '  ncor   = ',ncor
        write(6,*) '  nvir   = ',nvir
        write(6,*) '  nact   = ',nact
        write(6,*) '  nactri = ',nactri
      end if
      if(nact.eq.0) go to 200
      if(nact.eq.1) then
        ii=nactv(1)
        ntype(ii)=0
        go to 200
      end if
c
c  Further extract active orbitals
C    for nact >= 2
      call zero(qtemp,nactri)
      do 106 i=1,nact
      ii=nactv(i)
      do 106 j=1,i
      jj=nactv(j)
      ij=ioff(i)+j
      iijj=ioff(ii)+jj
      qtemp(ij)=ptemp(iijj)
  106 continue
      if(iprnt.ge.2) then
        write(6,6)
        call print(qtemp,nactri,nact,6)
      end if
C
c  Further extract active orbitals
      call zero(qeiv,naoso*nbasis)
      do 107 i=1,naoso
      do 107 j=1,nact
      jj=nactv(j)
      qeiv(i,j)=peiv(i,jj)
  107 continue
      if(iprnt.ge.2) then
        write(6,7)
        call matout(qeiv,naoso,nbasis,naoso,nact,6)
      end if
c
c  Diagonalize the extracted OPDM matrix
      call rsp(nact,nact,nactri,qtemp,qroot,1,eneiv,fv1,fv2)
      if(iprnt.ge.2) then
        write(6,8)
        call vecout(eneiv,qroot,nact,nact,nact,nact,6)
      end if
C
C  Construct the natural orbitals for ISYM symmetry
      call zero(reiv,naoso*nbasis)
      do 109 i=1,naoso
      do 109 j=1,nact
      jj=nact*(j-1)
      valu=a00
      do 108 k=1,nact
      kj=jj+k
      valu=valu+qeiv(i,k)*eneiv(kj)
  108 continue
      reiv(i,j)=valu
  109 continue
      if(iprnt.ge.2) then
        write(6,9) isym
        call vecout(reiv,qroot,naoso,nbasis,naoso,nact,6)
      end if
C
C  Reorder the natural orbitals
C     in order of occupation numbers
      call zero(rroot,nbasis)
      call zero(seiv,naoso*nbasis)
      do 110 i=1,naoso
      do 110 j=1,nact
      jj=nact-j+1
      seiv(i,j)=reiv(i,jj)
  110 continue
      do 111 i=1,nact
      ii=nact-i+1
      rroot(i)=qroot(ii)
  111 continue
      if(iprnt.ge.2) then
        write(6,10) isym
        call vecout(seiv,rroot,naoso,nbasis,naoso,nact,6)
      end if
c
c  Transfer the natural orbitals for ISYM symmetry to a whole matrix
  200 continue
      mact=0
      do 115 imo=1,ns
      mm=ntype(imo)
      ii=nstot+imo
      if(iprnt.ge.2) then
        write(6,*) ' imo = ',imo
        write(6,*) ' mm  = ',mm
        write(6,*) ' na  = ',na
        write(6,*) ' ii  = ',ii
        write(6,*) ' ns  = ',ns
      end if
      if(mm.eq.0) then
        do 112 i=1,naoso
        eno(i,ii)=peiv(i,imo)
  112   continue
        occno(ii)=proot(imo)
      end if
C  Only for nact >= 2
      if(mm.ne.0) then
        mact=mact+1
        na=nactv(mact)
        iii=nstot+na
        do 113 i=1,naoso
        eno(i,iii)=seiv(i,mact)
  113   continue
        occno(iii)=rroot(mact)
      end if
  115 continue
      if(iprnt.ge.2) then
        write(6,11) isym
        call vecout(eno,occno,naoso,nbasis,naoso,nbasis,6)
      end if
C
      nstot=nstot+ns
  120 CONTINUE
c
c  Final natural orbitals
      if(iprnt.ge.2) then
        write(6,12)
        call vecout(eno,occno,naoso,nbasis,naoso,nbasis,6)
      end if
C
C  One particle density matrix in the natural orbital basis
      do 125 i=1,nbasis
      ii=ioff(i+1)
      pno(ii)=occno(i)
  125 continue
      if(iprnt.ge.2) then
        write(6,13)
        call print(pno,ntri,nbasis,6)
      end if
C
      RETURN
      END
