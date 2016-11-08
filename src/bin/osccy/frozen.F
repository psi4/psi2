      SUBROUTINE FROZEN (FLOV,NOSYM,NVSYM,NFZO,NFZV,NDOC,NUOC,FZO,FZV,
     .                   ITYP,NIRRED,NO,NV,IW,NOACT,NVACT,NSYMHF,
     .                   nasym,mofz)
      IMPLICIT INTEGER(A-Z)
      INTEGER FLOV(8,2,4),NOSYM(NIRRED),NVSYM(NIRRED),NFZO(NIRRED),
     .        NFZV(NIRRED),NDOC(NIRRED),NUOC(NIRRED),FZO(NO),FZV(NV),
     .        nasym(nirred),mofz(*)
      CHARACTER*4 ITYP(NIRRED)
C
      DO 125 I=1,NIRRED
         NOSYM(I)=FLOV(I,2,2)-FLOV(I,2,1)+1
         IF(FLOV(I,2,2).LT.FLOV(I,2,1))NOSYM(I)=0
         NVSYM(I)=FLOV(I,1,4)-FLOV(I,1,3)+1
         IF(FLOV(I,1,4).LT.FLOV(I,1,3))NVSYM(I)=0
         NaSYM(I)=FLOV(I,1,2)-FLOV(I,1,1)-flov(i,2,2)+flov(i,2,1)
  125 CONTINUE
C
      NOACT=0
      NVACT=0
      DO 127 I=1,NIRRED
         NDOC(I)=NOSYM(I)-NFZO(I)
         NUOC(I)=NVSYM(I)-NFZV(I)
         NOACT=NOACT+NDOC(I)*2+nasym(i)
         NVACT=NVACT+NUOC(I)*2+nasym(i)
  127 CONTINUE
C
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)' >>>  FROZEN & ACTIVE ORBITALS   <<<'
      WRITE(6,*)
C
      DO 129 I=1,NSYMHF
         WRITE(IW,629)ITYP(I),NFZO(I),NDOC(I),nasym(i),NUOC(I),NFZV(I)
  129 CONTINUE
  629 FORMAT(2X,A4,3X,I3,'COR',I3,'DOC',I3,'ALP',i3,'UOC',I3,'VIR')
C
      WRITE(6,*)
      WRITE(6,*)' NOACT=',NOACT,'    NVACT=',NVACT
      WRITE(6,*)
      WRITE(6,*)' NO   =',NO   ,'    NV   =',NV 
      WRITE(6,*)
c     write(6,*)' NOTE: for allocation of integrals purpose open-shells'
c     write(6,*)'       are treated as virtuals ...                    '
C
      JO=0
      jmo=0
      DO 143 IS=1,NIRRED
         DO 131 I=1,NFZO(IS)
            JO=JO+1
            FZO(JO)=1
            jmo=jmo+1
            mofz(jmo)=1
  131    CONTINUE
         DO 132 I=1,NDOC(IS)
            JO=JO+1
            FZO(JO)=0
            jmo=jmo+1
            mofz(jmo)=0
  132    CONTINUE
         DO 133 I=1,Nasym(is)
            JO=JO+1
            FZO(JO)=0
  133    CONTINUE
  143 CONTINUE
      DO 153 IS=1,NIRRED
         DO 151 I=1,NFZO(IS)
            JO=JO+1
            FZO(JO)=1
  151    CONTINUE
         DO 152 I=1,NDOC(IS)
            JO=JO+1
            FZO(JO)=0
  152    CONTINUE
  153 CONTINUE
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
      DO 166 IS=1,NIRRED
         DO 164 I=1,NUOC(IS)+nasym(is)
            JV=JV+1
            FZV(JV)=0
            jmo=jmo+1
            mofz(jmo)=0
  164    CONTINUE
         DO 165 I=1,NFZV(IS)
            JV=JV+1
            FZV(JV)=1
            jmo=jmo+1
            mofz(jmo)=1
  165    CONTINUE
  166 CONTINUE
C
      DO 170 I=1,NIRRED
         FLOV(I,1,1)=FLOV(I,1,1)+NFZO(I)
         FLOV(I,1,4)=FLOV(I,1,4)-NFZV(I)
         FLOV(I,2,1)=FLOV(I,2,1)+NFZO(I)
         FLOV(I,2,4)=FLOV(I,2,4)-NFZV(I)
  170 CONTINUE
C
c     WRITE(6,*)'ISYM,FO,LO,FV,LV,alpha,beta'
      do 225,ispn=1,2
      DO 225 ISYM = 1,NIRRED
         FO=FLOV(ISYM,ispn,1)
         LO=FLOV(ISYM,ispn,2)
         FV=FLOV(ISYM,ispn,3)
         LV=FLOV(ISYM,ispn,4)
c     WRITE(6,639)ISYM,FO,LO,FV,LV
  639 format(i4,5x,4(i3,1x))
  225 CONTINUE
C
      RETURN
      END
