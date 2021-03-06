      SUBROUTINE GETDRT(KADD,LADD,IJADD,IJGRP,BFSYM,ORBTBF,IOUT,ORBSYM
     #,                 END,NGROUP,NREFS)
C
C  READ IN DRT INFORMATION
C
      IMPLICIT INTEGER (A-Z)
C
C
      COMMON /TAPES/ IT34,ITAPE8,IT47,IT52,IT91,IT93,IT95
      COMMON /DIMS/ SYMORB,NUMIJ,NIJ,NBF,NORBS,NMAX,NBFORB,NSYM0
      COMMON /OUT/    IOUT6,ICHECK,IDBG,ICNT
      DIMENSION KADD(SYMORB),LADD(SYMORB),IJADD(NUMIJ),IJGRP(NUMIJ)
      DIMENSION BFSYM(NBF),IOUT(NBF),ORBTBF(NORBS),ORBSYM(14   )
C
      CALL WREADW(ITAPE8,KADD,SYMORB,END,END)
      CALL WREADW(ITAPE8,LADD,SYMORB,END,END)
      CALL WREADW(ITAPE8,IJADD,NUMIJ,END,END)
      CALL WREADW(ITAPE8,IJGRP,NUMIJ,END,END)
      END=END+4*NORBS+NGROUP
      CALL WREADW(ITAPE8,IOUT,NBF,END,END)
      CALL WREADW(ITAPE8,ORBSYM,NORBS,END,END)
      CALL WREADW(ITAPE8,BFSYM,NBF,END,END)
      END=END+NREFS*NBF
      CALL WREADW(ITAPE8,ORBTBF,NORBS,END,END)
CDB   PRINT *,' ORBSYM= ',ORBSYM
CDB   PRINT *,' LADD= ',LADD
CDB   PRINT *,' KADD= ',KADD
CDB   PRINT *,' IJADD= ',IJADD
CDB   PRINT *,' ORBTBF= ',ORBTBF
CDB   PRINT *,' BFSYM= ',BFSYM
CDB   PRINT *,' IOUT=',IOUT
CDB   PRINT *,' IJGRP',IJGRP
      RETURN
      END
