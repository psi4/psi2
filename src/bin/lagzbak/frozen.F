C
C-----------------------------------------------------------------------
C
      SUBROUTINE FROZEN (FLOV,NLAMDA,NC,NFZO,NFZV,NDOC,NUOC,FZT,
     .                   ITYP,NIRRED,NO,NV,IW,NOACT,NVACT,NSYMHF)
      IMPLICIT INTEGER(A-Z)
      INTEGER FLOV(NIRRED,2),NFZO(NIRRED),NLAMDA(*),NC(*),
     .        NFZV(NIRRED),NDOC(NIRRED),NUOC(NIRRED),FZT(*)
      CHARACTER*4 ITYP(NIRRED)
C
      NOACT=0
      NVACT=0
      DO 127 I=1,NSYMHF
         NDOC(I)=NC(I)-NFZO(I)
         NUOC(I)=NLAMDA(I)-NC(I)-NFZV(I)
         NOACT=NOACT+NDOC(I)
         NVACT=NVACT+NUOC(I)
  127 CONTINUE
C
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)' >>>  FROZEN & ACTIVE ORBITALS   <<<'
      WRITE(6,*)
C
      DO 129 I=1,NSYMHF
         WRITE(IW,629)ITYP(I),NFZO(I),NDOC(I),NUOC(I),NFZV(I)
  129 CONTINUE
  629 FORMAT(2X,A4,3X,I3,'COR',I3,'DOC',I3,'UOC',I3,'VIR')
C
      WRITE(6,*)
      WRITE(6,*)' NOACT=',NOACT,'    NVACT=',NVACT
      WRITE(6,*)
      WRITE(6,*)'  NO  =',NO   ,'     NV  =',NV
      WRITE(6,*)
C
      JJ=0
      DO 140 IS=1,NSYMHF
         DO 131 I=1,NFZO(IS)
            JJ=JJ+1
            FZT(JJ)=1
  131    CONTINUE
         DO 132 I=1,NDOC(IS)
            JJ=JJ+1
            FZT(JJ)=0
  132    CONTINUE
         DO 134 I=1,NUOC(IS)
            JJ=JJ+1
            FZT(JJ)=0
  134    CONTINUE
         DO 135 I=1,NFZV(IS)
            JJ=JJ+1
            FZT(JJ)=1
  135    CONTINUE
  140 CONTINUE
C
      DO 170 I=1,NSYMHF
         FLOV(I,1)=FLOV(I,1)+NFZO(I)
         FLOV(I,2)=FLOV(I,2)-NFZV(I)
  170 CONTINUE
C
C     WRITE(6,*)'ISYM,FT,LT'
C     DO 225 ISYM = 1,NSYMHF
C        FT=FLOV(ISYM,1)
C        LT=FLOV(ISYM,2)
C     WRITE(6,*)ISYM,FT,LT
C 225 CONTINUE
C
      RETURN
      END
