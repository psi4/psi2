C     ////////////////////////////////////////////////////////////////
      SUBROUTINE TRN1E(SU,SFQ,FQ,SCRTCH,N,N0,NIJ,NMAT)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DOUBLE PRECISION E,SU(N,N0),SFQ(NMAT),FQ(NIJ),SCRTCH(N)
      INTEGER P,Q,PQ
C
      DO 10 K=1,NIJ
  10  FQ(K)=0.0D+00
      PQ=0
      DO 20 P=1,N
        DO 30 II=1,N
          SCRTCH(II)=0.0D+00
 30     CONTINUE
        E=1.0D+00
        DO 40 Q=1,P
          PQ=PQ+1
          IF (P.EQ.Q) E=0.5D+00
          F=SFQ(PQ)*E
          DO 50 J=1,N0
            SCRTCH(J)=SCRTCH(J)+SU(Q,J)*F
 50       CONTINUE
 40     CONTINUE
        IJ=0
        DO 60 I=1,N0
          SUPI=SU(P,I)
          SCRI=SCRTCH(I)
          DO 70 J=1,I
            IJ=IJ+1
            FQ(IJ)=FQ(IJ)+SCRTCH(J)*SUPI+SCRI*SU(P,J)
 70       CONTINUE
 60     CONTINUE
 20   CONTINUE
      RETURN
      END
