C
C---------------------------------------------------------------------
C
      SUBROUTINE FROZEN (FLOV,NOSYM,NVSYM,NFZO,NFZV,NDOC,NUOC,FZO,FZV,
     .                   ITYP,NIRRED,NO,NV,IW,NOACT,NVACT,NSYMHF)
      IMPLICIT INTEGER(A-Z)
      INTEGER FLOV(NIRRED,4),NOSYM(NIRRED),NVSYM(NIRRED),NFZO(NIRRED),
     .        NFZV(NIRRED),NDOC(NIRRED),NUOC(NIRRED),FZO(NO),FZV(NV)
      CHARACTER*4 ITYP(NIRRED)
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
      WRITE(6,*)' >>>  FROZEN & ACTIVE ORBITALS   <<<'
C
      DO 129 I=1,NSYMHF
         WRITE(IW,629)ITYP(I),NFZO(I),NDOC(I),NUOC(I),NFZV(I)
  129 CONTINUE
  629 FORMAT(2X,A4,3X,I3,'COR',I3,'DOC',I3,'UOC',I3,'VIR')
C
      WRITE(6,*)
      write(6,1010) 'Active',noact,nvact
      write(6,1010) 'All   ',no,nv
 1010 format(1X,A,' orbitals: Occupied: ',I4,', Unoccupied: ',I4)
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
      DO 170 I=1,NIRRED
         FLOV(I,1)=FLOV(I,1)+NFZO(I)
         FLOV(I,4)=FLOV(I,4)-NFZV(I)
  170 CONTINUE
C
C     WRITE(6,*)'ISYM,FO,LO,FV,LV'
C     DO 225 ISYM = 1,NIRRED
C        FO=FLOV(ISYM,1)
C        LO=FLOV(ISYM,2)
C        FV=FLOV(ISYM,3)
C        LV=FLOV(ISYM,4)
C     WRITE(6,*)ISYM,FO,LO,FV,LV
C 225 CONTINUE
C
      RETURN
      END
