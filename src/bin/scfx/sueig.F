C     /////////////////////////////////////////////////////////////////
      SUBROUTINE SUEIG(SU,EIG,SCRTCH,N,N0,A,CM,IRORDR,NIJ)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DOUBLE PRECISION SU(N,N0),EIG(N0,N0),SCRTCH(N,N0)
      DOUBLE PRECISION CM(2,N0),A(NIJ)
C
      CALL EBC(SCRTCH,SU,EIG,N,N0,N0)
      DO 10 I=1,N
        DO 10 J=1,N0
  10  SU(I,J)=SCRTCH(I,J)
C
C  SECTION TO ORDER VECTORS
C
      IF (IRORDR .EQ. 1) GOTO 160
C     PRINT,'  ORDERING VECTORS IN SUEIG'
      DO 150 I=1,N0
        II=I*(I+1)/2
        DO 150 J=I,N0
          JJ=J*(J+1)/2
C               IF (ABS(CM(2,J)).GT.1.0D-05 .AND. IRORDR.LT.2) GOTO 149
C               IF (ABS(CM(2,I)).GT.1.0D-05 .AND. IRORDR.LT.2) GOTO 149
          IF (A(II).GT.A(JJ)) THEN
            TEMP=A(II)
            A(II)=A(JJ)
            A(JJ)=TEMP
            DO 140 K=1,N
              TEMP=SU(K,I)
              SU(K,I)=SU(K,J)
  140       SU(K,J)=TEMP
          END IF
  149     CONTINUE
  150 CONTINUE
  160 CONTINUE
      RETURN
      END
