      SUBROUTINE LINBND(B,EX,EY,EZ,NXYZ,NNE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION B(NXYZ,NXYZ),EX(NNE,NNE),EY(NNE,NNE),EZ(NNE,NNE)
      COMMON/VIB103/NST,NBND,NLIB,NDEF,NTORS
      COMMON/VIB201/CHG(50),X(50),Y(50),Z(50),W(50)
      COMMON/VIB202/R(1275)
      COMMON/VIB205/IOFF(150)
      COMMON/VIB303/KB(150),LB(150),MB(150)
      DATA PI / 3.1415926536D+00 /
      DATA ONE / 1.0D+00 /
C
      NLHF=NLIB/2
C   N IS A DUMMY ATOM
      DO 101 I=1,NLHF
      I1=2*I-1+NST+NBND
      I2=I1+1
      K=KB(I)
      L=LB(I)
      M=MB(I)
      KX=3*K-2
      KY=3*K-1
      KZ=3*K
      LX=3*L-2
      LY=3*L-1
      LZ=3*L
      MX=3*M-2
      MY=3*M-1
      MZ=3*M
      XN=X(M)+ONE
      YN=Y(M)+ONE
      ZN=Z(M)+ONE
      RMN=DSQRT(3.0D+00)
      ML=IOFF(MAX0(M,L))+MIN0(M,L)
      RML=R(ML)
      KL=IOFF(MAX0(K,L))+MIN0(K,L)
      RNL=DSQRT((XN-X(L))**2+(YN-Y(L))**2+(ZN-Z(L))**2)
      SQCOS=(RNL**2-RMN**2-RML**2)/(2.0D+00*RML*RMN)
      AGM=PI-DACOS(SQCOS)
      EXMN=ONE/RMN
      EYMN=ONE/RMN
      EZMN=ONE/RMN
      B(I1,KX)=(DCOS(AGM)*EX(M,L)-EXMN)/(R(KL)*DSIN(AGM))
      B(I1,KY)=(DCOS(AGM)*EY(M,L)-EYMN)/(R(KL)*DSIN(AGM))
      B(I1,KZ)=(DCOS(AGM)*EZ(M,L)-EZMN)/(R(KL)*DSIN(AGM))
      B(I1,MX)=(DCOS(AGM)*EX(M,L)-EXMN)/(R(ML)*DSIN(AGM))
      B(I1,MY)=(DCOS(AGM)*EY(M,L)-EYMN)/(R(ML)*DSIN(AGM))
      B(I1,MZ)=(DCOS(AGM)*EZ(M,L)-EZMN)/(R(ML)*DSIN(AGM))
      B(I1,LX)=-B(I1,KX)-B(I1,MX)
      B(I1,LY)=-B(I1,KY)-B(I1,MY)
      B(I1,LZ)=-B(I1,KZ)-B(I1,MZ)
      B(I2,KX)=-(EYMN*EZ(M,L)-EZMN*EY(M,L))/(R(KL)*DSIN(AGM))
      B(I2,KY)=-(EZMN*EX(M,L)-EXMN*EZ(M,L))/(R(KL)*DSIN(AGM))
      B(I2,KZ)=-(EXMN*EY(M,L)-EYMN*EX(M,L))/(R(KL)*DSIN(AGM))
      B(I2,MX)=-(EYMN*EZ(M,L)-EZMN*EY(M,L))/(R(ML)*DSIN(AGM))
      B(I2,MY)=-(EZMN*EX(M,L)-EXMN*EZ(M,L))/(R(ML)*DSIN(AGM))
      B(I2,MZ)=-(EXMN*EY(M,L)-EYMN*EX(M,L))/(R(ML)*DSIN(AGM))
      B(I2,LX)=-B(I2,KX)-B(I2,MX)
      B(I2,LY)=-B(I2,KY)-B(I2,MY)
      B(I2,LZ)=-B(I2,KZ)-B(I2,MZ)
  101 CONTINUE
      RETURN
      END
