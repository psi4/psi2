      SUBROUTINE Orbs(q,ncore,nvirt,nolst)
      IMPLICIT REAL*8 (A-H,O-Z)                                         SMM13270
      DIMENSION Q(NTRI)
      COMMON/BASIS/NBASIS,NTRI,NTOT,NSYM
      common/posit/ioemc,iopmc,ilgmc
      COMMON/SIGNS/IOFF(32897),IPRNT    
      data a00,two / 0.0d+00 , 2.0d+00 /
      data xlim / 1.0d-08 /
    1 FORMAT(//,2X,' One particle density matrix in Orbs'/)
C                                                          
C   Read in the one-PDM
      call mread(q,iopmc)
      if(iprnt.gt.4) then
        write(6,1)
        call print(q,ntri,nbasis,6)
      end if
C  
C   Find numbers of the frozen core and deleted virtual orbitals
      ncore=0
      nvirt=0
      nolst=0
      do 101 i=1,nbasis
      ii=ioff(i+1)
      valu=q(ii)
      if(dabs(valu-two).lt.xlim) then
        ncore=ncore+1
      end if
      if(dabs(valu-a00).lt.xlim) then
        nvirt=nvirt+1
      end if
      if(dabs(valu).gt.xlim) then
        nolst=i
      end if
  101 continue
      write(3,*) ' ncore in orbs = ',ncore
      write(3,*) ' nvirt in orbs = ',nvirt
      write(3,*) ' nolst in orbs = ',nolst
      return
      end
