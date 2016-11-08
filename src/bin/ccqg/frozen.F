C
C-----------------------------------------------------------------------
C
      SUBROUTINE FROZEN (FLOV,NOSYM,NVSYM,NFZO,NFZV,NDOC,NUOC,FZO,FZV,
     .                   ITYP,NIRRED,NO,NV,IW,NOACT,NVACT,NSYMHF)
      IMPLICIT INTEGER(A-Z)
      INTEGER FLOV(NIRRED,4),NOSYM(NIRRED),NVSYM(NIRRED),NFZO(NIRRED),
     .        NFZV(NIRRED),NDOC(NIRRED),NUOC(NIRRED),FZO(NO),FZV(NV),
     .        ITYP(NIRRED)
C
      DO 125 I=1,NIRRED
         NOSYM(I)=FLOV(I,2)-FLOV(I,1)+1
         IF(FLOV(I,2).LT.FLOV(I,1))NOSYM(I)=0
         NVSYM(I)=FLOV(I,4)-FLOV(I,3)+1
         IF(FLOV(I,4).LT.FLOV(I,3))NVSYM(I)=0
  125 CONTINUE
C
      NOACT=0
      NVACT=0
      DO 127 I=1,NIRRED
         NDOC(I)=NOSYM(I)-NFZO(I)
         NUOC(I)=NVSYM(I)-NFZV(I)
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
         WRITE(6,629)ITYP(I),NFZO(I),NDOC(I),NUOC(I),NFZV(I)
  129 CONTINUE
  629 FORMAT(2X,A4,3X,I3,'COR',I3,'DOC',I3,'UOC',I3,'VIR')
C
      WRITE(6,*)
      WRITE(6,*)' NOACT=',NOACT,'    NVACT=',NVACT
      WRITE(6,*)
      WRITE(6,*)'  NO  =',NO   ,'     NV  =',NV
      WRITE(6,*)
C
      JO=0
      DO 133 IS=1,NIRRED
         DO 131 I=1,NFZO(IS)
            JO=JO+1
            FZO(JO)=1
  131    CONTINUE
         DO 132 I=1,NDOC(IS)
            JO=JO+1
            FZO(JO)=0
  132    CONTINUE
  133 CONTINUE
C
      WRITE(6,*)' FZO'
      WRITE(6,*)  FZO
C
      JV=0
      DO 136 IS=1,NIRRED
         DO 134 I=1,NUOC(IS)
            JV=JV+1
            FZV(JV)=0
  134    CONTINUE
         DO 135 I=1,NFZV(IS)
            JV=JV+1
            FZV(JV)=1
  135    CONTINUE
  136 CONTINUE
C
      WRITE(6,*)' FZV'
      WRITE(6,*)  FZV
C
      DO 170 I=1,NIRRED
         FLOV(I,1)=FLOV(I,1)+NFZO(I)
         FLOV(I,4)=FLOV(I,4)-NFZV(I)
  170 CONTINUE
C
      RETURN
      END
