      SUBROUTINE RD78 (MOINTS)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 MOINTS(*),BUF(51200)
      integer IBUF(2*51200)
      equivalence(buf,ibuf)
c
      ind(i,j)=i*(i-1)/2+j
c     write(6,*)'intowp(1)',intowp(1)
      idisk=78
      jout=6
      CALL rfile(IDISK)
      call sread(idisk,lenbuf,1)
c     WRITE (JOUT,*) ' LENBUF=',LENBUF
c     call flush(6)
      ITEMP=255
      INTLEN=(INTOWP(LENBUF)-2)/INTOWP(1)
      MAXVAL=INTOWP(INTLEN)/(1+INTOWP(1))
      IBOFF=(MAXVAL+2)/INTOWP(1)
  111 CALL SREAD(IDISK,iBUF,intowp(LENBUF))
      IFLG=IBUF(1)
      MBUF=IBUF(2)
c     WRITE (JOUT,*) ' IFLG , mbuf,',IFLG,mbuf   
      DO 12 II=1,MBUF
      IJKL=IBUF(II+2)
         I=ISHFT(IJKL,-24)
         J=IAND(ITEMP,ISHFT(IJKL,-16))
         IJ=ind(i,j)
         K=IAND(ITEMP,ISHFT(IJKL,- 8))
         L=IAND(ITEMP,IJKL)
         KL=ind(k,l)
         IJKL=ind(ij,kl)
         MOINTS(IJKL)=BUF(II+IBOFF)
c        write(6,678)i,j,k,l,ijkl,moints(ijkl)
  678 format(5i6,f20.12)
   12 CONTINUE
      IF (IFLG.EQ.0) GOTO 111
c     call flush(6)
      call rclose(idisk,3)
c
      RETURN
      END
