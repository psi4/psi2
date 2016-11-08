      PROGRAM LOCAPM
C  This program localizes M.O. by Pipek-Mezey Method.
C  It reads overlap basis function and CMO coefficient matrix C from
C     FILE30. which is produced by program SCF.
C  Overlap integrals in terms of AO's (matrix S) are calculated in this
C     program.
C
C  No input required.
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      integer mxntri,mxnao
      parameter (mxnao=255,mxntri=32640)
      real*8 E(2,mxnao), temp(40,mxnao)
c mxnao is used to allocate storage for both ao and so arrays
      dimension s(mxntri),sa(mxnao,mxnao)
      dimension su(mxntri),usotao(mxnao,mxnao),eig(mxnao)
c s is the ao overlap matrix (trigonal)
c s is later dimensioned as nbatri so idims > = nbatri
c
c  common to basset, derao and pm;
      common/funcs/natom,n3n,nbfao,nbatri,nocc,mb(20)
c
c  common to basset and derao;
      common/coord/chgs(112),c(3,112),ict(112),ns(112),ks(112)
      common/gauss/zs(200),cs(200)
      common/nshel/knuc(200),ktype(200),kprim(200),kloc(200),
     1 kmin(200),kmax(200),kng(200),nshell,nprim
      common/signs/ipq(500),iprnt
c
c  common to derao and ecal routines ;
      common/exyzs/ex(7,7,13),ey(7,7,13),ez(7,7,13)
c
c  common to readmo;
cccccclj      common/paramo/cmo(mxnao,mxnao),eigocc(mxnao),msfru
      real*8 cmo(mxnao,mxnao)
      real*8 eigocc(mxnao)
      integer msfru
c
      integer errcod,prdup,prcntr
c
 500  format(//1x,'not enough memory in routine overlap for s matrix'/)
c
      call psinit('APPEND')
      call tstart(6)

      errcod = prcntr('ON BRIEF')
      errcod = prdup()

      iprnt = 0
c
c   read in basis set information
      call basset
c
c   calculate overlap integrals
      if(nbatri.lt.mxntri)then
        call derao(s)
      else
        write(6,500)
      endif
c
c     change overlap integrals from triangle form to a matrix
      jj=0
      do 20 i=1,nbfao
      do 20 j=1,i
      jj=jj+1
      sa(j,i)=s(jj)
      sa(i,j)=s(jj)
   20 continue
c
c     read MO's from file30
      CALL READMO(cmo,eigocc,msfru,su,usotao,eig,mxnao,mxntri)
c
C     do 30 i=1,nbfao
C  30 write (6,10) (sa(i,j),j=1,nbfao)
C     do 40 j=1,nocc
C  40 WRITE (6,10) (CMO(J,K),K=1,msfru)
C     WRITE (6,10) (eigocc(J),j=1,nocc)
C  10 format (7f10.6)
C
C     localization
      CALL PM(SA,CMO,EIGOCC,e,temp,mxnao,msfru)
C
      CALL TSTOP(6)
      call psdone()
      STOP
C
      END
