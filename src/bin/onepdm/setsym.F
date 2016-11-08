C////////////////////////////////////////////////////////////////////
      SUBROUTINE SETSYM(PGROUP)
      CHARACTER PGROUP*3,SYMLBL*3
      INTEGER PRINT,MAX,NOFLG,PRPFLG,ROOTI,ROOTF,IRSYM,SYMOFF
      INTEGER NTOTS,SYMDIF
      COMMON /OPS/  PRINT,MAX,NOFLG,PRPFLG,ROOTI,ROOTF,IRSYM
      COMMON /SYMM/ SYMOFF(8),NTOTS(8,2),SYMDIF,SYMLBL(15)
      SYMLBL(1)='A'
      SYMLBL(2)='A'''
      SYMLBL(3)='A"'
      SYMLBL(4)='A1'
      SYMLBL(5)='A2'
      SYMLBL(6)='B1'
      SYMLBL(7)='B2'
      SYMLBL(8)=' AG'
      SYMLBL(9)='B1G'
      SYMLBL(10)='B2G'
      SYMLBL(11)='B3G'
      SYMLBL(12)=' AU'
      SYMLBL(13)='B1U'
      SYMLBL(14)='B2U'
      SYMLBL(15)='B3U'
      IF(PGROUP.EQ.' C1'.OR.PGROUP.EQ.' C2') THEN
      IRSYM=1
      SYMLBL(2)='A'
      SYMLBL(3)='B'
      ELSE IF (PGROUP.EQ.' CS') THEN
      IRSYM=2
      ELSE IF (PGROUP.EQ.' CI') THEN
      SYMLBL(2)='AG'
      SYMLBL(3)='AU'
      IRSYM=2
      ELSE IF (PGROUP.EQ.'C2V') THEN
      IRSYM=2
      ELSE IF (PGROUP.EQ.' D2') THEN
      SYMLBL(4)=' A'
      SYMLBL(5)='B1'
      SYMLBL(6)='B2'
      SYMLBL(7)='B3'
      IRSYM=2
      ELSE IF (PGROUP.EQ.'C2H') THEN
      SYMLBL(4)='AG'
      SYMLBL(5)='BG'
      SYMLBL(6)='AU'
      SYMLBL(7)='BU'
      IRSYM=2
      ELSE IF (PGROUP.EQ.'D2H') THEN
      IRSYM=3
      END IF
      RETURN
      END
