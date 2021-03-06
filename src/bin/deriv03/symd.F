      SUBROUTINE SYMD(D1T)
      IMPLICIT REAL*8 (A-H,O-Z)
cets030691
#include <error.h>
      integer ip, prcntr
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOM,N3N,NATRI,NSYM
CTPH  COMMON/GAUSC/T(36),PTR(3,3,8),DTR(6,6,8)
      COMMON/GAUSC/T(225),PTR(3,3,8),DTR(6,6,8),FTR(10,10,8),
     1             GTR(15,15,8)
      COMMON/GAUSD/INVT(8),ISO(200,8),LBLAT(2,50),ICT(50,8)
      COMMON/PARA1/NSHELL,NT,MXT,NISO
      COMMON/TOLER/DLIM
      DIMENSION D1T(3,NATOM)
      DATA ZERO,ONE / 0.0D+00 , 1.0D+00 /
    1 FORMAT(//,' GRADIENT BEFORE SYMDE')
    2 FORMAT (/,' ICT:',(15I4))
    3 FORMAT (/,' INVT:',(15I4))
    4 FORMAT (//,' TRANSFORMATION MATRICES:',/,(3(/,1X,3F15.6),//))
C
      WRITE(3,1)
      CALL MATOUT(D1T,3,NATOM,3,NATOM,3)
      ip=prcntr('IS_ON DEBUG')
      if(ip.eq.1) then
        WRITE(6,1)
        CALL MATOUT(D1T,3,NATOM,3,NATOM,6)
        WRITE(6,2) ICT
        WRITE(6,3) INVT
        WRITE(6,4) (((PTR(IQ,JQ,KQ),IQ=1,3),JQ=1,3),KQ=1,NT)
      endif
      WRITE(3,2) ICT
      WRITE(3,3) INVT
      WRITE(3,4) (((PTR(IQ,JQ,KQ),IQ=1,3),JQ=1,3),KQ=1,NT)
C
      IF(NT.EQ.1) RETURN
C
C   SYMMETRIZE GRADIENT VECTOR
C
      DO 105 IC=1,NATOM
      DO 101 IT=1,NT
      IF(ICT(IC,IT).GT.IC) GO TO 105
  101 CONTINUE
      DEDX=ZERO
      DEDY=ZERO
      DEDZ=ZERO
      DO 102 IT=1,NT
      ICNU=ICT(IC,IT)
      DEDXP=D1T(1,ICNU)
      DEDYP=D1T(2,ICNU)
      DEDZP=D1T(3,ICNU)
      DEDX=DEDX+DEDXP*PTR(1,1,IT)+DEDYP*PTR(2,1,IT)+DEDZP*PTR(3,1,IT)
      DEDY=DEDY+DEDXP*PTR(1,2,IT)+DEDYP*PTR(2,2,IT)+DEDZP*PTR(3,2,IT)
      DEDZ=DEDZ+DEDXP*PTR(1,3,IT)+DEDYP*PTR(2,3,IT)+DEDZP*PTR(3,3,IT)
  102 CONTINUE
      D1T(1,IC)=DEDX
      D1T(2,IC)=DEDY
      D1T(3,IC)=DEDZ
      DO 104 IT=1,NT
      ICNU=ICT(IC,IT)
      IF(ICNU.EQ.IC) GO TO 104
      IF(IT.EQ.NT) GO TO 201
      IT1=IT+1
      DO 103 JT=IT1,NT
      IF(ICT(IC,JT).EQ.ICNU) GO TO 104
  103 CONTINUE
  201 CONTINUE
      JT=INVT(IT)
      D1T(1,ICNU)=D1T(1,IC)*PTR(1,1,JT)+D1T(2,IC)*PTR(2,1,JT)
     1                                 +D1T(3,IC)*PTR(3,1,JT)
      D1T(2,ICNU)=D1T(1,IC)*PTR(1,2,JT)+D1T(2,IC)*PTR(2,2,JT)
     1                                 +D1T(3,IC)*PTR(3,2,JT)
      D1T(3,ICNU)=D1T(1,IC)*PTR(1,3,JT)+D1T(2,IC)*PTR(2,3,JT)
     1                                 +D1T(3,IC)*PTR(3,3,JT)
  104 CONTINUE
  105 CONTINUE
C
      FACT=ONE/DFLOAT(NT)
      DO 106 J=1,NATOM
      DO 106 I=1,3
      D1T(I,J)=D1T(I,J)*FACT
  106 CONTINUE
C
      RETURN
      END
