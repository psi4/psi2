      function i21d(ia,jb)
      common/signs/ioff(32897),iprnt
c
      i21d=ioff(max0(ia,jb))+min0(ia,jb)
c
      return
      end

      subroutine qflin(a,idim,in,im,det,ip100,im100)
      implicit real*8 (a-h,o-z)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c   Note : The following three variables are defined         :
c     in quadruple precision                                 :
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c     Last update on October 19, 1993 by Dr. Yukio Yamaguchi :
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      REAL*16 DET,D,SIGN,Q100,QM100,QP100
C
C     LINEAR SIMULTANEOUS EQUATION
C
C     A(IN*IN) * X(IN*IM) = B(IN*IM)
C
C     A & B SHOULD BE STORED ON A(IN*(IN+IM))
C     SOLUTION X WILL BE STORED ON B PART IN DIMENSION A.
C
      DIMENSION A(IDIM,1)
      data a00,ONE / 0.0D+00 , 1.0D+00 /
      data qm100,qp100 / 1.0q-100 , 1.0q+100 /
      data q100 / 100.0q+00 /
C
      ip100=0
      im100=0
      N=IN
      NR=IM
      JMAX=N+NR
cyy   sign=1.0d+00
      SIGN=1.0Q+00
C M IS THE STAGE OF ELIMINATION
      DO 49 M=1,N
      TEMP=a00
      DO 41 I=M,N
      IF(M.GT.1)A(I,M)=A(I,M)-DOTT(A(I,1),IDIM,A(1,M),1,M-1)
      AVAL=A(I,M)
      IF(DABS(AVAL).LE.TEMP)GOTO 41
      TEMP=DABS(AVAL)
      IMAX=I
 41   CONTINUE
      IF(TEMP.LE.a00)GOTO 999
      IF(IMAX.EQ.M)GOTO 45
      SIGN=-SIGN
      DO 44 J=1,JMAX
      STOR=A(M,J)
      A(M,J)=A(IMAX,J)
      A(IMAX,J)=STOR
 44   CONTINUE
 45   CONTINUE
      JJ=M+1
      IF(JJ.GT.JMAX)GOTO 49
      IF(M.GT.1)GOTO 47
      DO 46 J=JJ,JMAX
      A(1,J)=A(1,J)/A(1,1)
 46   CONTINUE
      D=A(1,1)
      GOTO 49
 47   CONTINUE
      DO 48 J=JJ,JMAX
      A(M,J)=(A(M,J)-DOTT(A(M,1),IDIM,A(1,J),1,M-1))/A(M,M)
 48   CONTINUE
      D=D*A(M,M)
c
c======================================================
c   Adjust the magnitude of determinant value
c    to avoid under- or over-flows
c======================================================
      if(qabs(d).gt.qp100) then
        d=d/q100
        ip100=ip100+1
      end if
      if(qabs(d).lt.qm100) then
        d=d*q100
        im100=im100+1
      end if
c======================================================
c
 49   CONTINUE
      IF(NR.EQ.0) RETURN
      DO 59 I=1,NR
      NPI=N+I
      DO 58 K=2,N
      J=N+1-K
      A(J,NPI)=A(J,NPI)-DOTT(A(J,J+1),IDIM,A(J+1,NPI),1,K-1)
 58   CONTINUE
 59   CONTINUE
      DET=D*SIGN
      RETURN
c
C ON ZERO PIVOT, SET DET=0.AND RETURN TO CALLING PROGRAM NOV 1972
cyy   999 det=0.0d+00
 999  DET=0.0Q+00
      RETURN
      END
