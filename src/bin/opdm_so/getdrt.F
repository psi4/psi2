      SUBROUTINE GETDRT(KADD,LADD,IJADD,IJGRP,IEND,NGROUP)              DRT01860
      IMPLICIT REAL*8 (A-H,O-Z)                                         DRT01870
      DIMENSION KADD(ISYMRB),LADD(ISYMRB),IJADD(NUMIJ),IJGRP(NUMIJ)     DRT01880
      COMMON/DIMES/ISYMRB,NUMIJ,NBF,NORBS,NMAX                          DRT01890
      COMMON/ORDCI/IDRSYM(256),IORDRT(256)                              DRT01900
      COMMON/SIGNS/IOFF(256),IPRNT                                      DRT01910
      COMMON/TAPES/ITAP53,ITAP58 
    1 FORMAT(//,2X,' INFORMATION IN GETDRT'/)                           DRT01930
    2 FORMAT(2X,20I5)                                                   DRT01940
C                                                                       DRT01950
      CALL WREADW(ITAP58,KADD,ISYMRB,IEND,IEND)                         DRT01960
      CALL WREADW(ITAP58,LADD,ISYMRB,IEND,IEND)                         DRT01970
      CALL WREADW(ITAP58,IJADD,NUMIJ,IEND,IEND)                         DRT01980
      CALL WREADW(ITAP58,IJGRP,NUMIJ,IEND,IEND)                         DRT01990
      IEND=IEND+4*NORBS+NGROUP+NBF                                      DRT02000
      CALL WREADW(ITAP58,IDRSYM,NORBS,IEND,IEND)                        DRT02010
      IEND=IEND+NBF*2                                                   DRT02020
      CALL WREADW(ITAP58,IORDRT,NORBS,IEND,IEND)                        DRT02030
C                                                                       DRT02040
      IF(IPRNT.LE.4) GO TO 201                                          DRT02050
      WRITE(6,1)                                                        DRT02060
      WRITE(6,*) ' KADD ',KADD                                          DRT02070
      WRITE(6,*) ' LADD ',LADD                                          DRT02080
      WRITE(6,*) ' IJADD ',IJADD                                        DRT02090
      WRITE(6,*) ' IJGRP ',IJGRP                                        DRT02100
      WRITE(6,*) ' IDRSYM '                                             DRT02110
      WRITE(6,2) (IDRSYM(I),I=1,NBF)                                    DRT02120
      WRITE(6,*) ' IORDRT '                                             DRT02130
      WRITE(6,2) (IORDRT(I),I=1,NBF)                                    DRT02140
C                                                                       DRT02150
  201 CONTINUE                                                          DRT02160
      RETURN                                                            DRT02170
      END                                                               DRT02180
