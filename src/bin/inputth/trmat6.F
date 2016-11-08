      SUBROUTINE TRMAT6(PTR,DTR,FTR,GTR,T)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV
      COMMON /INOUT/  IOUT,ICHECK,ITAP30
      COMMON /MAXDIM/ MXAT,MXTRN,MXPR,MXSH,MXSYM,MXCF,MXT,MXISO,MXPSYM
      COMMON /OUTPUT/ NPRINT,NGEOM,NORMF,NORMP
      COMMON /SYMTRY/ NT
      COMMON /TRANSF/ XSMAL,YSMAL,ZSMAL,XNEW,YNEW,ZNEW,XP,YP,ZP
      COMMON /FRAME/  U1,U2,U3,V1,V2,V3,W1,W2,W3,X0,Y0,Z0
C
      REAL*8 PTR(3,1),DTR(6,1),FTR(10,1),GTR(15,1),T(MXT)
C
      DIMENSION PNAME(3),DNAME(6)
C
      PARAMETER (ONE=1.0D+00)
C     PARAMETER (DSQRT3=SQRT(3.0D+00))
      DATA PNAME /4HX   ,4HY   ,4HZ   /
      DATA DNAME /4HXX  ,4HYY  ,4HZZ  ,4HXY  ,4HXZ  ,4HYZ  /
      DSQRT3=SQRT(3.0D+00)
C
 8869 FORMAT(//)
 8863 FORMAT(/,' TRANSFORMATION OF THE BASIS FUNCTIONS',/)
 8862 FORMAT(8X,10(3X,A4,3X))
 8861 FORMAT(2X,A4,2X,10F10.6)
 8860 FORMAT(1H1)
 8858 FORMAT(/,21X,'TRANSFORMATION NUMBER',I4,/)
 8855 FORMAT(/)
C
C  CALCULATE TRANSFORMS OF P AND D FUNCTIONS FOR ALL SYMMETRY OPERATIONS
C
      X=X0+ONE
      Y=Y0
      Z=Z0
      CALL LOCAL(X,Y,Z,XS,YS,ZS)
      XSMAL=XS
      YSMAL=YS
      ZSMAL=ZS
      DO 1760 IT=1,NT
        NN=9*(IT-1)
        CALL TRANS(T,NN)
        CALL ROT
        N=3*(IT-1)
        PTR(1,N+1)=XP-X0
        PTR(2,N+1)=YP-Y0
        PTR(3,N+1)=ZP-Z0
 1760 CONTINUE
      X=X0
      Y=Y0+ONE
      Z=Z0
      CALL LOCAL(X,Y,Z,XS,YS,ZS)
      XSMAL=XS
      YSMAL=YS
      ZSMAL=ZS
      DO 1770 IT=1,NT
        NN=9*(IT-1)
        CALL TRANS(T,NN)
        CALL ROT
        N=3*(IT-1)
        PTR(1,N+2)=XP-X0
        PTR(2,N+2)=YP-Y0
        PTR(3,N+2)=ZP-Z0
 1770 CONTINUE
      X=X0
      Y=Y0
      Z=Z0+ONE
      CALL LOCAL(X,Y,Z,XS,YS,ZS)
      XSMAL=XS
      YSMAL=YS
      ZSMAL=ZS
      DO 1780 IT=1,NT
        NN=9*(IT-1)
        CALL TRANS(T,NN)
        CALL ROT
        N=3*(IT-1)
        PTR(1,N+3)=XP-X0
        PTR(2,N+3)=YP-Y0
        PTR(3,N+3)=ZP-Z0
 1780 CONTINUE
      DO 1820 IT=1,NT
        NP=3*(IT-1)
        ND=6*(IT-1)
        NF=10*(IT-1)
        NG=15*(IT-1)
        DO 1788 I=1,6
          GO TO (1781,1782,1783,1784,1785,1786) ,I
 1781     J=1
          K=1
          GO TO 1787
 1782     J=2
          K=2
          GO TO 1787
 1783     J=3
          K=3
          GO TO 1787
 1784     J=1
          K=2
          GO TO 1787
 1785     J=1
          K=3
          GO TO 1787
 1786     J=2
          K=3
 1787     DTR(1,ND+I)=PTR(1,NP+J)*PTR(1,NP+K)
          DTR(2,ND+I)=PTR(2,NP+J)*PTR(2,NP+K)
          DTR(3,ND+I)=PTR(3,NP+J)*PTR(3,NP+K)
          DTR(4,ND+I)=PTR(1,NP+J)*PTR(2,NP+K)+PTR(2,NP+J)*PTR(1,NP+K)
          DTR(5,ND+I)=PTR(1,NP+J)*PTR(3,NP+K)+PTR(3,NP+J)*PTR(1,NP+K)
          DTR(6,ND+I)=PTR(2,NP+J)*PTR(3,NP+K)+PTR(3,NP+J)*PTR(2,NP+K)
 1788   CONTINUE
        DO 2988 I=1,10
          GO TO (2781,2782,2783,2784,2785,2786,2787,2788,2789,2790) ,I
 2781     J=1
          K=1
          L=1
          GO TO 2791
 2782     J=2
          K=2
          L=2
          GO TO 2791
 2783     J=3
          K=3
          L=3
          GO TO 2791
 2784     J=1
          K=1
          L=2
          GO TO 2791
 2785     J=1
          K=1
          L=3
          GO TO 2791
 2786     J=1
          K=2
          L=2
          GO TO 2791
 2787     J=2
          K=2
          L=3
          GO TO 2791
 2788     J=1
          K=3
          L=3
          GO TO 2791
 2789     J=2
          K=3
          L=3
          GO TO 2791
 2790     J=1
          K=2
          L=3
          GO TO 2791
 2791     FTR(1,NF+I)=PTR(1,NP+J)*PTR(1,NP+K)*PTR(1,NP+L)
          FTR(2,NF+I)=PTR(2,NP+J)*PTR(2,NP+K)*PTR(2,NP+L)
          FTR(3,NF+I)=PTR(3,NP+J)*PTR(3,NP+K)*PTR(3,NP+L)
          FTR(4,NF+I)=PTR(1,NP+J)*PTR(1,NP+K)*PTR(2,NP+L)
     1               +PTR(1,NP+J)*PTR(2,NP+K)*PTR(1,NP+L)
     2               +PTR(2,NP+J)*PTR(1,NP+K)*PTR(1,NP+L)
C
          FTR(5,NF+I)=PTR(1,NP+J)*PTR(1,NP+K)*PTR(3,NP+L)
     1               +PTR(3,NP+J)*PTR(1,NP+K)*PTR(1,NP+L)
     2               +PTR(1,NP+J)*PTR(3,NP+K)*PTR(1,NP+L)
C
          FTR(6,NF+I)=PTR(1,NP+J)*PTR(2,NP+K)*PTR(2,NP+L)
     1               +PTR(2,NP+J)*PTR(1,NP+K)*PTR(2,NP+L)
     2               +PTR(2,NP+J)*PTR(2,NP+K)*PTR(1,NP+L)
C
          FTR(7,NF+I)=PTR(3,NP+J)*PTR(2,NP+K)*PTR(2,NP+L)
     1               +PTR(2,NP+J)*PTR(3,NP+K)*PTR(2,NP+L)
     2               +PTR(2,NP+J)*PTR(2,NP+K)*PTR(3,NP+L)
C
          FTR(8,NF+I)=PTR(1,NP+J)*PTR(3,NP+K)*PTR(3,NP+L)
     1               +PTR(3,NP+J)*PTR(1,NP+K)*PTR(3,NP+L)
     2               +PTR(3,NP+J)*PTR(3,NP+K)*PTR(1,NP+L)
C
          FTR(9,NF+I)=PTR(2,NP+J)*PTR(3,NP+K)*PTR(3,NP+L)
     1               +PTR(3,NP+J)*PTR(2,NP+K)*PTR(3,NP+L)
     2               +PTR(3,NP+J)*PTR(3,NP+K)*PTR(2,NP+L)
C
          FTR(10,NF+I)=PTR(1,NP+J)*PTR(2,NP+K)*PTR(3,NP+L)
     1                +PTR(3,NP+J)*PTR(1,NP+K)*PTR(2,NP+L)
     2                +PTR(3,NP+J)*PTR(2,NP+K)*PTR(1,NP+L)
     3                +PTR(2,NP+J)*PTR(1,NP+K)*PTR(3,NP+L)
     4                +PTR(2,NP+J)*PTR(3,NP+K)*PTR(1,NP+L)
     5                +PTR(1,NP+J)*PTR(3,NP+K)*PTR(2,NP+L)
C
 2988   CONTINUE
        DO 5988 I=1,15
          GO TO (5781,5782,5783,5784,5785,5786,5787,5788,5789,5790,
     1           6781,6782,6783,6784,6785), I
 5781     J=1
          K=1
          L=1
          M=1
          GO TO 5791
 5782     J=2
          K=2
          L=2
          M=2
          GO TO 5791
 5783     J=3
          K=3
          L=3
          M=3
          GO TO 5791
 5784     J=1
          K=1
          L=1
          M=2
          GO TO 5791
 5785     J=1
          K=1
          L=1
          M=3
          GO TO 5791
 5786     J=1
          K=2
          L=2
          M=2
          GO TO 5791
 5787     J=2
          K=2
          L=2
          M=3
          GO TO 5791
 5788     J=1
          K=3
          L=3
          M=3
          GO TO 5791
 5789     J=2
          K=3
          L=3
          M=3
          GO TO 5791
 5790     J=1
          K=1
          L=2
          M=2
          GO TO 5791
 6781     J=1
          K=1
          L=3
          M=3
          GO TO 5791
 6782     J=2
          K=2
          L=3
          M=3
          GO TO 5791
 6783     J=1
          K=1
          L=2
          M=3
          GO TO 5791
 6784     J=1
          K=2
          L=2
          M=3
          GO TO 5791
 6785     J=1
          K=2
          L=3
          M=3
          GO TO 5791
 5791     GTR(1,NG+I)=PTR(1,NP+J)*PTR(1,NP+K)*PTR(1,NP+L)*PTR(1,NP+M)
          GTR(2,NG+I)=PTR(2,NP+J)*PTR(2,NP+K)*PTR(2,NP+L)*PTR(2,NP+M)
          GTR(3,NG+I)=PTR(3,NP+J)*PTR(3,NP+K)*PTR(3,NP+L)*PTR(3,NP+M)
          GTR(4,NG+I)=PTR(1,NP+J)*PTR(1,NP+K)*PTR(1,NP+L)*PTR(2,NP+M)
     1               +PTR(1,NP+J)*PTR(2,NP+K)*PTR(1,NP+L)*PTR(1,NP+M)
     2               +PTR(2,NP+J)*PTR(1,NP+K)*PTR(1,NP+L)*PTR(1,NP+M)
     3               +PTR(1,NP+J)*PTR(1,NP+K)*PTR(2,NP+L)*PTR(1,NP+M)
C
          GTR(5,NG+I)=PTR(1,NP+J)*PTR(1,NP+K)*PTR(1,NP+L)*PTR(3,NP+M)
     1               +PTR(1,NP+J)*PTR(3,NP+K)*PTR(1,NP+L)*PTR(1,NP+M)
     2               +PTR(3,NP+J)*PTR(1,NP+K)*PTR(1,NP+L)*PTR(1,NP+M)
     3               +PTR(1,NP+J)*PTR(1,NP+K)*PTR(3,NP+L)*PTR(1,NP+M)
C
          GTR(6,NG+I)=PTR(2,NP+J)*PTR(2,NP+K)*PTR(2,NP+L)*PTR(1,NP+M)
     1               +PTR(2,NP+J)*PTR(1,NP+K)*PTR(2,NP+L)*PTR(2,NP+M)
     2               +PTR(1,NP+J)*PTR(2,NP+K)*PTR(2,NP+L)*PTR(2,NP+M)
     3               +PTR(2,NP+J)*PTR(2,NP+K)*PTR(1,NP+L)*PTR(2,NP+M)
C
          GTR(7,NG+I)=PTR(2,NP+J)*PTR(2,NP+K)*PTR(2,NP+L)*PTR(3,NP+M)
     1               +PTR(2,NP+J)*PTR(3,NP+K)*PTR(2,NP+L)*PTR(2,NP+M)
     2               +PTR(3,NP+J)*PTR(2,NP+K)*PTR(2,NP+L)*PTR(2,NP+M)
     3               +PTR(2,NP+J)*PTR(2,NP+K)*PTR(3,NP+L)*PTR(2,NP+M)
C
          GTR(8,NG+I)=PTR(3,NP+J)*PTR(3,NP+K)*PTR(3,NP+L)*PTR(1,NP+M)
     1               +PTR(3,NP+J)*PTR(1,NP+K)*PTR(3,NP+L)*PTR(3,NP+M)
     2               +PTR(1,NP+J)*PTR(3,NP+K)*PTR(3,NP+L)*PTR(3,NP+M)
     3               +PTR(3,NP+J)*PTR(3,NP+K)*PTR(1,NP+L)*PTR(3,NP+M)
C
          GTR(9,NG+I)=PTR(3,NP+J)*PTR(3,NP+K)*PTR(3,NP+L)*PTR(2,NP+M)
     1               +PTR(3,NP+J)*PTR(2,NP+K)*PTR(3,NP+L)*PTR(3,NP+M)
     2               +PTR(2,NP+J)*PTR(3,NP+K)*PTR(3,NP+L)*PTR(3,NP+M)
     3               +PTR(3,NP+J)*PTR(3,NP+K)*PTR(2,NP+L)*PTR(3,NP+M)
C
          GTR(10,NG+I)=PTR(1,NP+J)*PTR(1,NP+K)*PTR(2,NP+L)*PTR(2,NP+M)
     1                +PTR(1,NP+J)*PTR(2,NP+K)*PTR(1,NP+L)*PTR(2,NP+M)
     2                +PTR(1,NP+J)*PTR(2,NP+K)*PTR(2,NP+L)*PTR(1,NP+M)
     3                +PTR(2,NP+J)*PTR(2,NP+K)*PTR(1,NP+L)*PTR(1,NP+M)
     4                +PTR(2,NP+J)*PTR(1,NP+K)*PTR(2,NP+L)*PTR(1,NP+M)
     5                +PTR(2,NP+J)*PTR(1,NP+K)*PTR(1,NP+L)*PTR(2,NP+M)
C
          GTR(11,NG+I)=PTR(1,NP+J)*PTR(1,NP+K)*PTR(3,NP+L)*PTR(3,NP+M)
     1                +PTR(1,NP+J)*PTR(3,NP+K)*PTR(1,NP+L)*PTR(3,NP+M)
     2                +PTR(1,NP+J)*PTR(3,NP+K)*PTR(3,NP+L)*PTR(1,NP+M)
     3                +PTR(3,NP+J)*PTR(3,NP+K)*PTR(1,NP+L)*PTR(1,NP+M)
     4                +PTR(3,NP+J)*PTR(1,NP+K)*PTR(3,NP+L)*PTR(1,NP+M)
     5                +PTR(3,NP+J)*PTR(1,NP+K)*PTR(1,NP+L)*PTR(3,NP+M)
C
          GTR(12,NG+I)=PTR(2,NP+J)*PTR(2,NP+K)*PTR(3,NP+L)*PTR(3,NP+M)
     1                +PTR(2,NP+J)*PTR(3,NP+K)*PTR(2,NP+L)*PTR(3,NP+M)
     2                +PTR(2,NP+J)*PTR(3,NP+K)*PTR(3,NP+L)*PTR(2,NP+M)
     3                +PTR(3,NP+J)*PTR(3,NP+K)*PTR(2,NP+L)*PTR(2,NP+M)
     4                +PTR(3,NP+J)*PTR(2,NP+K)*PTR(3,NP+L)*PTR(2,NP+M)
     5                +PTR(3,NP+J)*PTR(2,NP+K)*PTR(2,NP+L)*PTR(3,NP+M)
C
          GTR(13,NG+I)=PTR(1,NP+J)*PTR(1,NP+K)*PTR(2,NP+L)*PTR(3,NP+M)
     1                +PTR(1,NP+J)*PTR(1,NP+K)*PTR(3,NP+L)*PTR(2,NP+M)
     2                +PTR(1,NP+J)*PTR(2,NP+K)*PTR(1,NP+L)*PTR(3,NP+M)
     3                +PTR(1,NP+J)*PTR(3,NP+K)*PTR(1,NP+L)*PTR(2,NP+M)
     4                +PTR(1,NP+J)*PTR(3,NP+K)*PTR(2,NP+L)*PTR(1,NP+M)
     5                +PTR(1,NP+J)*PTR(2,NP+K)*PTR(3,NP+L)*PTR(1,NP+M)
     6                +PTR(2,NP+J)*PTR(3,NP+K)*PTR(1,NP+L)*PTR(1,NP+M)
     7                +PTR(3,NP+J)*PTR(2,NP+K)*PTR(1,NP+L)*PTR(1,NP+M)
     8                +PTR(3,NP+J)*PTR(1,NP+K)*PTR(2,NP+L)*PTR(1,NP+M)
     9                +PTR(2,NP+J)*PTR(1,NP+K)*PTR(3,NP+L)*PTR(1,NP+M)
     A                +PTR(3,NP+J)*PTR(1,NP+K)*PTR(1,NP+L)*PTR(2,NP+M)
     B                +PTR(2,NP+J)*PTR(1,NP+K)*PTR(1,NP+L)*PTR(3,NP+M)
C
          GTR(14,NG+I)=PTR(2,NP+J)*PTR(2,NP+K)*PTR(1,NP+L)*PTR(3,NP+M)
     1                +PTR(2,NP+J)*PTR(2,NP+K)*PTR(3,NP+L)*PTR(1,NP+M)
     2                +PTR(2,NP+J)*PTR(1,NP+K)*PTR(2,NP+L)*PTR(3,NP+M)
     3                +PTR(2,NP+J)*PTR(3,NP+K)*PTR(2,NP+L)*PTR(1,NP+M)
     4                +PTR(2,NP+J)*PTR(3,NP+K)*PTR(1,NP+L)*PTR(2,NP+M)
     5                +PTR(2,NP+J)*PTR(1,NP+K)*PTR(3,NP+L)*PTR(2,NP+M)
     6                +PTR(1,NP+J)*PTR(3,NP+K)*PTR(2,NP+L)*PTR(2,NP+M)
     7                +PTR(3,NP+J)*PTR(1,NP+K)*PTR(2,NP+L)*PTR(2,NP+M)
     8                +PTR(3,NP+J)*PTR(2,NP+K)*PTR(1,NP+L)*PTR(2,NP+M)
     9                +PTR(1,NP+J)*PTR(2,NP+K)*PTR(3,NP+L)*PTR(2,NP+M)
     A                +PTR(3,NP+J)*PTR(2,NP+K)*PTR(2,NP+L)*PTR(1,NP+M)
     B                +PTR(1,NP+J)*PTR(2,NP+K)*PTR(2,NP+L)*PTR(3,NP+M)
C
          GTR(15,NG+I)=PTR(3,NP+J)*PTR(3,NP+K)*PTR(1,NP+L)*PTR(2,NP+M)
     1                +PTR(3,NP+J)*PTR(3,NP+K)*PTR(2,NP+L)*PTR(1,NP+M)
     2                +PTR(3,NP+J)*PTR(1,NP+K)*PTR(3,NP+L)*PTR(2,NP+M)
     3                +PTR(3,NP+J)*PTR(2,NP+K)*PTR(3,NP+L)*PTR(1,NP+M)
     4                +PTR(3,NP+J)*PTR(2,NP+K)*PTR(1,NP+L)*PTR(3,NP+M)
     5                +PTR(3,NP+J)*PTR(1,NP+K)*PTR(2,NP+L)*PTR(3,NP+M)
     6                +PTR(1,NP+J)*PTR(2,NP+K)*PTR(3,NP+L)*PTR(3,NP+M)
     7                +PTR(2,NP+J)*PTR(1,NP+K)*PTR(3,NP+L)*PTR(3,NP+M)
     8                +PTR(2,NP+J)*PTR(3,NP+K)*PTR(1,NP+L)*PTR(3,NP+M)
     9                +PTR(1,NP+J)*PTR(3,NP+K)*PTR(2,NP+L)*PTR(3,NP+M)
     A                +PTR(2,NP+J)*PTR(3,NP+K)*PTR(3,NP+L)*PTR(1,NP+M)
     B                +PTR(1,NP+J)*PTR(3,NP+K)*PTR(3,NP+L)*PTR(2,NP+M)
C
 5988   CONTINUE
 1820 CONTINUE
      IF(NORMF.EQ.1.AND.NORMP.EQ.1) GO TO 1831
      DO 1830 IT=1,NT
        ND=6*(IT-1)
        NF=10*(IT-1)
        NG=15*(IT-1)
        DO 1822 I=1,6
          IF(I.GT.3) GO TO 1821
          DTR(4,ND+I)=DTR(4,ND+I)/DSQRT3
          DTR(5,ND+I)=DTR(5,ND+I)/DSQRT3
          DTR(6,ND+I)=DTR(6,ND+I)/DSQRT3
          GO TO 1822
 1821     DTR(1,ND+I)=DTR(1,ND+I)*DSQRT3
          DTR(2,ND+I)=DTR(2,ND+I)*DSQRT3
          DTR(3,ND+I)=DTR(3,ND+I)*DSQRT3
 1822   CONTINUE
CJ
        DO 2822 I=1,10
          IF(I.LE.3) THEN
            FTR(4,NF+I)=FTR(4,NF+I)*DSQRT(3.0D0/15.0D0)
            FTR(5,NF+I)=FTR(5,NF+I)*DSQRT(3.0D0/15.0D0)
            FTR(6,NF+I)=FTR(6,NF+I)*DSQRT(3.0D0/15.0D0)
            FTR(7,NF+I)=FTR(7,NF+I)*DSQRT(3.0D0/15.0D0)
            FTR(8,NF+I)=FTR(8,NF+I)*DSQRT(3.0D0/15.0D0)
            FTR(9,NF+I)=FTR(9,NF+I)*DSQRT(3.0D0/15.0D0)
            FTR(10,NF+I)=FTR(10,NF+I)*DSQRT(1.0D0/15.0D0)
          ELSE IF (I.LE.9) THEN
            FTR(1,NF+I)=FTR(1,NF+I)*DSQRT(15.0D0/3.0D0)
            FTR(2,NF+I)=FTR(2,NF+I)*DSQRT(15.0D0/3.0D0)
            FTR(3,NF+I)=FTR(3,NF+I)*DSQRT(15.0D0/3.0D0)
            FTR(10,NF+I)=FTR(10,NF+I)*DSQRT(1.0D0/3.0D0)
          ELSE
            FTR(1,NF+I)=FTR(1,NF+I)*DSQRT(15.0D0)
            FTR(2,NF+I)=FTR(2,NF+I)*DSQRT(15.0D0)
            FTR(3,NF+I)=FTR(3,NF+I)*DSQRT(15.0D0)
            FTR(4,NF+I)=FTR(4,NF+I)*DSQRT(3.0D0)
            FTR(5,NF+I)=FTR(5,NF+I)*DSQRT(3.0D0)
            FTR(6,NF+I)=FTR(6,NF+I)*DSQRT(3.0D0)
            FTR(7,NF+I)=FTR(7,NF+I)*DSQRT(3.0D0)
            FTR(8,NF+I)=FTR(8,NF+I)*DSQRT(3.0D0)
            FTR(9,NF+I)=FTR(9,NF+I)*DSQRT(3.0D0)
          END IF
 2822   CONTINUE
CJ
        DO 8822 I=1,15
          IF(I.LE.3) THEN
            GTR(4,NG+I)=GTR(4,NG+I)*DSQRT(1.0D0/7.0D0)
            GTR(5,NG+I)=GTR(5,NG+I)*DSQRT(1.0D0/7.0D0)
            GTR(6,NG+I)=GTR(6,NG+I)*DSQRT(1.0D0/7.0D0)
            GTR(7,NG+I)=GTR(7,NG+I)*DSQRT(1.0D0/7.0D0)
            GTR(8,NG+I)=GTR(8,NG+I)*DSQRT(1.0D0/7.0D0)
            GTR(9,NG+I)=GTR(9,NG+I)*DSQRT(1.0D0/7.0D0)
            GTR(10,NG+I)=GTR(10,NG+I)*DSQRT(1.0D0/35.0D0)
            GTR(11,NG+I)=GTR(11,NG+I)*DSQRT(1.0D0/35.0D0)
            GTR(12,NG+I)=GTR(12,NG+I)*DSQRT(1.0D0/35.0D0)
            GTR(13,NG+I)=GTR(13,NG+I)*DSQRT(1.0D0/35.0D0)
            GTR(14,NG+I)=GTR(14,NG+I)*DSQRT(1.0D0/35.0D0)
            GTR(15,NG+I)=GTR(15,NG+I)*DSQRT(1.0D0/35.0D0)
          ELSE IF (I.LE.9) THEN
            GTR(1,NG+I)=GTR(1,NG+I)*DSQRT(7.0D0)
            GTR(2,NG+I)=GTR(2,NG+I)*DSQRT(7.0D0)
            GTR(3,NG+I)=GTR(3,NG+I)*DSQRT(7.0D0)
            GTR(10,NG+I)=GTR(10,NG+I)*DSQRT(1.0D0/5.0D0)
            GTR(11,NG+I)=GTR(11,NG+I)*DSQRT(1.0D0/5.0D0)
            GTR(12,NG+I)=GTR(12,NG+I)*DSQRT(1.0D0/5.0D0)
            GTR(13,NG+I)=GTR(13,NG+I)*DSQRT(1.0D0/5.0D0)
            GTR(14,NG+I)=GTR(14,NG+I)*DSQRT(1.0D0/5.0D0)
            GTR(15,NG+I)=GTR(15,NG+I)*DSQRT(1.0D0/5.0D0)
          ELSE
            GTR(1,NG+I)=GTR(1,NG+I)*DSQRT(35.0D0)
            GTR(2,NG+I)=GTR(2,NG+I)*DSQRT(35.0D0)
            GTR(3,NG+I)=GTR(3,NG+I)*DSQRT(35.0D0)
            GTR(4,NG+I)=GTR(4,NG+I)*DSQRT(5.0D0)
            GTR(5,NG+I)=GTR(5,NG+I)*DSQRT(5.0D0)
            GTR(6,NG+I)=GTR(6,NG+I)*DSQRT(5.0D0)
            GTR(7,NG+I)=GTR(7,NG+I)*DSQRT(5.0D0)
            GTR(8,NG+I)=GTR(8,NG+I)*DSQRT(5.0D0)
            GTR(9,NG+I)=GTR(9,NG+I)*DSQRT(5.0D0)
          END IF
 8822   CONTINUE
CJ
 1830 CONTINUE
 1831 CONTINUE
C
C     ----- PRINT MATRICES IF NPRINT.EQ.5 -----
C
      IF(NPRINT.NE.5) GO TO 2500
      WRITE(IW,8863)
      DO 2155 IT=1,NT
        WRITE(IW,8860)
        WRITE(IW,8858) IT
        NP=3*(IT-1)
        WRITE(IW,8862) (PNAME(J),J=1,3)
        WRITE(IW,8855)
        DO 2150 I=1,3
 2150   WRITE(IW,8861) PNAME(I),(PTR(I,NP+J),J=1,3)
        WRITE(IW,8869)
        ND=6*(IT-1)
        WRITE(IW,8862) (DNAME(J),J=1,6)
        WRITE(IW,8855)
        DO 2151 I=1,6
 2151   WRITE(IW,8861) DNAME(I),(DTR(I,ND+J),J=1,6)
        WRITE(IW,8869)
 2155 CONTINUE
 2500 CONTINUE
      RETURN
      END
