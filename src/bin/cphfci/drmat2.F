ctph change several numij to nbatri, and nbf to nbfao
      SUBROUTINE DRMAT2(VEC,SCR1,SCR2,SS,BUF1,BUF2,IOFF,JT1C)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer norbs,numij,ndf,nbset
ctph Change dimensions of the following line
      DIMENSION VEC(NBFAO,NORBS),SCR1(NORBS),SCR2(NBATRI),
     1BUF1(NBATRI,NDF),BUF2(NBATRI,NDF,NBSET),SS(NBATRI)
      INTEGER IOFF(NUMIJ),JT1C(1)
      COMMON /TAPES/ NOUT,ICHK,IT30,IT42,IT47,IT52,IT58,IT94,IT96,IT98
      COMMON /DIMSG/ ISYMRB,NUMIJ,NBF,NMAX,NBFORB,NSYM,NDF,NORBS
      COMMON /FUNCS/ IOCC,JOCC,KOCC,NTYPES,NBSET,NAT
ctph add common statements
      COMMON /GRAD/ E1A(150)
      COMMON /DIMAO/ NBFAO,NBATRI
      COMMON /LOC42/ IA42(100),IBF2AT(284)
      COMMON /IEXT/ IPRINT,IGRP,ICIDIP
C
CTJL  WRITE(ICHK,2222)
C2222 FORMAT(/'  IN DERMAT')
      CALL RFILE(IT94)
CIBM  CALL SREW(IT94)
cets030691      CALL WREADW(IT42,IBF2AT,NORBS,IA42(11),JUNK)
      call srew(it42)
      call sread(it42,ibf2at,intowp(nbf))
c position file42 to beginning of s integrals
ctph read in forces into e1a
      call sread(it42,E1A,intowp(ndf))
      call sread(it42,buf1,intowp(ndf*ndf))
      call zero(buf1,nbatri*ndf)
C
C READ IN IT42 DATA; TRANSFORM TO MO BASIS; STORE ON IT94
C
C READ S DERIVATIVE INTEGRALS
cets030691      CALL WREADW(IT42,BUF1,INTOWP(nbatri*3),IA42(12),JUNK)
      do 100 i=1,3
        call sread(it42,buf1(1,i),intowp(nbatri))
c       write(6,*) 'ss i=',i
c       call print(buf1(1,i),nbatri,nbfao,6)
 100  continue
      CALL SAMAT(BUF1,nbatri,SS,NUMIJ,VEC,IOFF,SCR1,SCR2,JT1C)
C
C READ H DERIVATIVE INTEGRALS
cets030691      CALL WREADW(IT42,BUF1,INTOWP(nbatri*NDF),IA42(14),JUNK)
      call zero(buf1,nbatri*ndf)
      do 102 i=1,ndf
        call sread(it42,buf1(1,i),intowp(nbatri))
c       write(6,*) 'h i=',i
c       call print(buf1(1,i),nbatri,nbfao,6)
  102 continue
      IS=JT1C(NDF)+INTOWP(NUMIJ)
      DO 107 IABC=1,NDF
        IF (IAND(IPRINT,4) .NE. 0) THEN
          WRITE(NOUT,1001) IABC
 1001     FORMAT(/'  THE HA MATRIX FOR DEGREE OF FREEDOM ',I5)
          CALL PRINT(BUF1(1,IABC),nbatri,nbfao,NOUT)
        END IF
        CALL MOCONV(BUF1(1,IABC),nbatri,SS,NUMIJ,VEC,SCR1,SCR2)
        JT1C(IABC+NDF)=IS
        CALL WWRITW(IT94,SS,INTOWP(NUMIJ),IS,IS)
        IF (IAND(IPRINT,4) .NE. 0) THEN
          WRITE(NOUT,1002) IABC
 1002     FORMAT(/'  THE HM MATRIX FOR DEGREE OF FREEDOM ',I5)
          CALL PRINT(SS,NUMIJ,NORBS,NOUT)
          WRITE(NOUT,1003) IABC,JT1C(IABC+NDF),IS
 1003     FORMAT(/'  HM FOR DEGREE ',I5,'  WRITTEN TO WORD ',I8,
     1    '  NEXT WORD IS ',I8)
        END IF
 107  CONTINUE
C
C  READ 2E DERIVATIVE INTEGRALS (T MATRIX)
CTJL  WRITE(*,*) ' IA42(17),NBSET,NDF,nbatri = ',
CTJL 1    IA42(17),NBSET,NDF,nbatri
cets030691     CALL WREADW(IT42,BUF2,INTOWP(NBSET*NDF*nbatri),IA42(17),JUNK)
      do 110 i=1,nbset
        do 111 j=1,ndf
          call sread(it42,buf2(1,j,i),intowp(nbatri))
c         write(6,*) 'bmat i j',i,j
c         call print(buf2(1,j,i),nbatri,nbfao,6)
 111    continue
 110  continue
c
      IJT=2*NDF
      DO 108 IABC = 1,NDF
        DO 108 ISET = 1,NBSET
          IF (IAND(IPRINT,8) .NE. 0) THEN
            WRITE(NOUT,1004) ISET,IABC
 1004       FORMAT(/'  THE TA MATRIX FOR SET ',I5,'  DEGREE ',
     *      I5)
            CALL PRINT(BUF2(1,IABC,ISET),nbatri,nbfao,NOUT)
          END IF
          CALL MOCONV(BUF2(1,IABC,ISET),nbatri,SS,NUMIJ,VEC,SCR1,
     *    SCR2)
          IJT=IJT+1
          JT1C(IJT)=IS
          CALL WWRITW(IT94,SS,INTOWP(NUMIJ),IS,IS)
          IF (IAND(IPRINT,8) .NE. 0) THEN
            WRITE(NOUT,1005) ISET,IABC
 1005       FORMAT(/'  THE TM MATRIX FOR SET ',I5,'  DEGREE ',
     *      I5)
            CALL PRINT(SS,NUMIJ,NORBS,NOUT)
            WRITE(NOUT,1006) IABC,IJT,JT1C(IJT),IS
 1006       FORMAT(/'  TM FOR DEGREE ',I5,'  IN JT1C ',I5/
     1      '  WRITTEN TO WORD',I5,'  NEXT WORD IS ',I8)
          END IF
 108  CONTINUE
      RETURN
      END
