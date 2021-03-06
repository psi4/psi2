      SUBROUTINE ALLOC(IOFF,NSOFF,OFFA1,NSO,NSYMHF,BFDRT,IOUT)
C
C
C
      IMPLICIT INTEGER (A-Z)
C
C
      COMMON /OUT/    IOUT6,ICHECK,IDBG,ICNT
      COMMON /DIMS/   SYMORB,NUMIJ,NIJ,NBF,NORBS,NMAX,NBFORB,NSYM0
C
      DIMENSION IOFF(NIJ),NSOFF(NSYMHF),OFFA1(NSYMHF),NSO(NSYMHF)
      DIMENSION IOUT(NBF),BFDRT(NORBS)
C
C     PRINT,'  IN ALLOC:'
C     PRINT,'  NSYMHF=',NSYMHF,' NBF=',NBF,'  NORBS=',NORBS
C     IDBG2=1
CDB   IF ((IDBG.AND.2).NE.0) IDBG2=1
C  LOWER TRIANGLE OFFSET
      IOFF(1)=0
C     IF (IDBG2.NE.0)WRITE(IOUT6,'('        IOFF')')
      DO 1 I=1,NIJ-1
           IOFF(I+1)=IOFF(I)+I
    1 CONTINUE
C     IF(IDBG2.NE.0)PRINT,' IOFF= ',IOFF
C  SYMMETRY ORBITAL OFFSET AND ONE ELECTRON SYM BLOCK OFFSET
      OFFA1(1)=0
      NSOFF(1)=0
C     IF(IDBG2.NE.0)WRITE(IOUT6,'('     I  NSO NSOFF OFFA1')')
      DO 2 I=2,NSYMHF
           IF (NSO(I-1) .EQ. 0) THEN
                NSOFF(I)=NSOFF(I-1)
                OFFA1(I)=OFFA1(I-1)
           ELSE
                NSOFF(I)=NSO(I-1)+NSOFF(I-1)
                OFFA1(I)=OFFA1(I-1)+(IOFF(NSO(I-1))+NSO(I-1))
C     IF(IDBG2.NE.0)WRITE(IOUT6,'(4I5)')I,NSO(I),NSOFF(I),OFFA1(I)
           END IF
    2 CONTINUE
C BASIS FUNCTION TO DRT ORDERING ARRAY
      JJ=0
      DO 3 II=1,NBF
           IF (IOUT(II) .GT. 0) THEN
                JJ=JJ+1
                BFDRT(JJ)=IOUT(II)
           END IF
    3 CONTINUE
C     PRINT,'  BFDRT=',BFDRT
      IF (JJ .NE. NORBS) THEN
C     WRITE(IOUT6,'('  IN ALLOC, JJ=',I4,'  NORBS=',I4)')JJ,NORBS
           CALL MABORT
      END IF
      RETURN
      END
