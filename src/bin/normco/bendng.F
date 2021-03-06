      SUBROUTINE BENDNG(B,EX,EY,EZ,NXYZ,NNE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION B(NXYZ,NXYZ),EX(NNE,NNE),EY(NNE,NNE),EZ(NNE,NNE)
      COMMON/VIB103/NST,NBND,NLIB,NDEF,NTORS
      COMMON/VIB202/R(1275)
      COMMON/VIB205/IOFF(150)
      COMMON/VIB302/KA(150),LA(150),MA(150)
      COMMON/VIB306/AGLB(150),AGLT(150),AGLO(150),AGLR(150)
C
C   THE CALCULATION OF S VECTORS FOR ANGLE BENDINGS
      DO 101 I=1,NBND
      II=I+NST
      K=KA(I)
      L=LA(I)
      M=MA(I)
      KX=3*K-2
      KY=3*K-1
      KZ=3*K
      LX=3*L-2
      LY=3*L-1
      LZ=3*L
      MX=3*M-2
      MY=3*M-1
      MZ=3*M
      AG=AGLB(I)
      LK=IOFF(MAX0(L,K))+MIN0(L,K)
      LM=IOFF(MAX0(L,M))+MIN0(L,M)
      B(II,KX)=(DCOS(AG)*EX(L,K)-EX(L,M))/(R(LK)*DSIN(AG))
      B(II,KY)=(DCOS(AG)*EY(L,K)-EY(L,M))/(R(LK)*DSIN(AG))
      B(II,KZ)=(DCOS(AG)*EZ(L,K)-EZ(L,M))/(R(LK)*DSIN(AG))
      B(II,MX)=(DCOS(AG)*EX(L,M)-EX(L,K))/(R(LM)*DSIN(AG))
      B(II,MY)=(DCOS(AG)*EY(L,M)-EY(L,K))/(R(LM)*DSIN(AG))
      B(II,MZ)=(DCOS(AG)*EZ(L,M)-EZ(L,K))/(R(LM)*DSIN(AG))
      B(II,LX)=-B(II,KX)-B(II,MX)
      B(II,LY)=-B(II,KY)-B(II,MY)
  101 B(II,LZ)=-B(II,KZ)-B(II,MZ)
      RETURN
      END
