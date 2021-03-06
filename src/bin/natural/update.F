      subroutine update(eno,eig,ee,ess)
      implicit real*8 (a-h,o-z)
      dimension i30(200),nlamda(10)
      common/basis/nbasis,nbfao,nbfso,nbatri,ntri
      common/funcs/natom,n3n,nsym
      common/signs/ioff(500),iprnt
      dimension eno(nbasis,nbasis),eig(nbasis)
      dimension ee(nbasis*nbasis),ess(nbasis,nbasis)
      data a00 / 0.0d+00 /
    1 format(//,2x,' The old Eso Matrix in file30'/)
    2 format(//,2x,' The new Eno Matrix'/)
    3 format(//,2x,' The Eno Matrix has been replaced in file30'/)
C
C   GET CONSTANTS FROM TAPE30
      ITAP30=30
      CALL RFILE(ITAP30)
      CALL WREADW(ITAP30,I30,200,101,JUNK)
C
      IEND    = I30(1)
      MPOINT  = I30(2)
      MCONST  = I30(3)
      MCALCS  = I30(4)
      NCALCS  = I30(5)
      MXCOEF  = I30(42)
      write(3,*) '  iend   = ',iend
      write(3,*) '  mpoint = ',mpoint
      write(3,*) '  mconst = ',mconst
      write(3,*) '  mcalcs = ',mcalcs
      write(3,*) '  ncalcs = ',ncalcs
      write(3,*) '  mxcoef = ',mxcoef
c
C   Find location of eigenvectors
      IPOS=101+MCONST+MPOINT+NCALCS-1
      write(3,*) '  ipos   = ',ipos
      CALL WREADW(ITAP30,LOCCAL,1,IPOS,JUNK)
      IPOS=LOCCAL+60
      write(3,*) '  loccal = ',loccal
      CALL WREADW(ITAP30,I30,20,IPOS,JUNK)
      LOCVEC=I30(1)
      locsav=i30(1)
      write(3,*) '  locvec = ',locvec
C
C   READ IN SO-MO EIGENVECTORS
      CALL WREADW(ITAP30,EE,MXCOEF*2,LOCVEC,LOCVEC)
      CALL WREADW(ITAP30,EIG,NBASIS*2,LOCVEC,LOCVEC)
      LOCVEC=LOCVEC+NSYM
      CALL WREADW(ITAP30,NLAMDA,NSYM,LOCVEC,LOCVEC)
C
C   FORM SO-MO EIGENVECTOR MATRIX
      call zero(ess,nbasis*nbasis)
      I=0
      J=0
      DO 102 ISYM=1,NSYM
      N=NLAMDA(ISYM)
      DO 101 K=J+1,J+N
      DO 101 L=J+1,J+N
      I=I+1
      ESS(L,K)=EE(I)
  101 CONTINUE
      J=J+N
  102 CONTINUE
C
      if(iprnt.ge.1) then
        write(6,1)
        call eivout(ess,eig,nbasis,nbasis,nbasis,nbasis,6)
      end if
C
      if(iprnt.ge.1) then
        write(6,2)
        call eivout(eno,eig,nbasis,nbasis,nbasis,nbasis,6)
      end if
c
C   Store SO-MO (ENO) EIGENVECTOR in file30 format
      call zero(ee,nbasis*nbasis)
      I=0
      J=0
      DO 105 ISYM=1,NSYM
      N=NLAMDA(ISYM)
      DO 104 K=J+1,J+N
      DO 104 L=J+1,J+N
      I=I+1
      EE(I)=ENO(L,K)
  104 CONTINUE
      J=J+N
  105 CONTINUE
      CALL Wwritw(ITAP30,EE,MXCOEF*2,LOCsav,junk)
      write(6,3)
C
      call rclose(itap30,3)
      return
      end
