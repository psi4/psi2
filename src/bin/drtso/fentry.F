      subroutine fentry(cc,ia,maxcor)
C**********************************************************************
C     PROGRAM DRTSO                                                   *
C**********************************************************************
C   THIS PROGRAM REORDERS MO INTEGRALS AND TWO-PDMS IN DRT ORDER      *
C   TO THOSE IN SORTED-SCF ORDER FOR CI AND MCSCF WAVEFUNCTIONS       *
C**********************************************************************
C*   Last updated on Febuary 2, 1995 by Evan Bolton                   *
C   Reason:  To put in PSI 2.0 input parsing.  See comments Ceeb      *
C**********************************************************************
C*  Last updated on November 05, 1993 by Dr. Yukio Yamaguchi          *
c*  Dafault inputs were added                                         * 
C**********************************************************************
C*  Last updated on July 12, 1993 by Dr. Yukio Yamaguchi              *
c*  Statements ``call tstart'' and ``call initmf''                    *
c*  have been move to after ``call psinit'' line                      *
C**********************************************************************
C*  Last updated on August 29, 1992 by Dr. Yukio Yamaguchi            *
c*  Modification for the Unix Version                                 *
C**********************************************************************
C*  LAST UPDATED ON MAY 02, 1988 BY YUKIO YAMAGUCHI                   *
C**********************************************************************      
C                                                                  
      IMPLICIT REAL*8 (A-H,O-Z)                                
      dimension cc(maxcor),ia(maxcor*2)
      character*8 adum,bdum

Ceeb For New Style Input
      logical oldinp
      integer errcod
      integer fwkclr, fwkadd, frdi

      DIMENSION LBLDRT(26),NDUM(20)                                     DRT00090
      COMMON/CALIF/LPARA(1024),APARA(1024)                              DRT00120
      COMMON/DIMES/ISYMRB,NUMIJ,NBF,NORBS,NMAX                          DRT00130
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          DRT00150
      COMMON/MFSEC/MASTER,NSECT                                         DRT00160
      COMMON/ORDCI/IDRSYM(256),IORDRT(256)                              DRT00170
      COMMON/ORDSO/IORSO(256),IDRTSO(256)                               DRT00180
      common/posit/ioeci,ioemc,iopci,iopmc
      COMMON/SIGNS/IOFF(256),IPRNT                                      DRT00190
      COMMON/TAPES/ITAP35,ITAP52,ITAP53,ITAP56,ITAP58                   DRT00200
    1 FORMAT(//,2X,' THE DRT-SO PROGRAM  (Version of 02/02/95)',/,11x,
     1 'Written by Dr. Yukio Yamaguchi',//,11X, 'Modifications for ',
     2 'PSI 2.0 use by Evan Bolton',///)
    2 FORMAT(A8,2X,A8,2X)
    3 FORMAT(5I5)                                                       DRT00240
    4 FORMAT(2X,' DRT LABEL = ',26A3)                                   DRT00250
    5 FORMAT(/,2X,' PARAMETERS'/           
     1 2X,' NBF     = ',I8/                                             DRT00270
     2 2X,' NSYM    = ',I8/                                             DRT00280
     3 2X,' NORBS   = ',I8/                                             DRT00290
     4 2X,' NWKS    = ',I8/                                             DRT00300
     5 2X,' NMAX    = ',I8/                                             DRT00310
     6 2X,' NUMIJ   = ',I8/                                             DRT00320
     7 2X,' NGROUP  = ',I8/                                             DRT00330
     8 2X,' ISYMRB  = ',I8/                                             DRT00340
     9 2X,' IPRNT   = ',I8/)                                            DRT00370
    6 FORMAT(//,2X,' NUMBER OF ORBITALS DOES NOT MATCH',                DRT00380
     1             ' NBF = ',I5,' ISUM = ',I5/)                         DRT00390
    7 FORMAT(//,2X,' I , IORDRT , IORSO , IDRTSO'/)                     DRT00400
    8 FORMAT(2X,5I5)                                                    DRT00410
    9 FORMAT(//,2X,' WRONG INPUT',A8/)                                  DRT00420
C                                                                       DRT00430
C     ITAP35 = CONTAINS THE TWO ELECTRON INTEGRALS (SORTED SO)          DRT00480
C     ITAP52 = CONTAINS THE MO ONE AND TWO ELECTRON INTEGRALS (DRT)     DRT00490
C     ITAP53 = CONTAINS THE ONE AND TWO PDM'S IN THE MO BASIS (DRT)     DRT00500
C     ITAP56 = CONTAINS TWO PDM'S IN THE MO BASIS (SORTED SO)           DRT00510
C     ITAP58 = CONTAINS THE DRT INFORMATION                             DRT00520
C                                                                       DRT00530
      ITAPE3=3                                                          DRT00540
      INPUT=5                                                           DRT00550
      ITAP35=35                                                         DRT00560
      ITAP52=52                                                         DRT00570
      ITAP53=53                                                         DRT00580
      ITAP56=56                                                         DRT00590
      ITAP58=58                                                         DRT00600
c
c=====================================================================
C    Followings are the positions of arrays in the master file
c=====================================================================
      ioeci=16
      ioemc=23
      iopci=35
      iopmc=39
c=====================================================================
c
      call psinit('APPEND')
      CALL TSTART(6)                 
      CALL INITMF(1)              
C                                                                       DRT00620
      WRITE(6,1)                                                        DRT00650
      write(6,*) '  maxcor  = ',maxcor
c
      CALL LOCATE(INPUT,'# DRTSO ##',IERR)                              DRT00630

Ceeb  For New Style Input
      oldinp = .false.
      iprnt=0

      IF (IERR.eq.0) THEN
        oldinp = .true.
        write(6,*) 
        write(6,*) 'WARNING: using old style ''# DRTSO ##'''
      END IF

      if(oldinp) then
        READ(5,2) adum,bdum
        READ(5,3) IPRNT           
      else
        errcod = fwkclr()
        errcod = fwkadd(5,'DRTSO')
        errcod = fwkadd(5,'DEFAULT')

        errcod = frdi('PRINT',iprnt)
      end if

Ceeb End New Style Input

      IOFF(1)=0                                                         DRT00870
      DO 101 I=1,255                                                    DRT00880
      IOFF(I+1)=IOFF(I)+I                                               DRT00890
  101 continue
C                                                                       DRT00900
C   READ DIMENSIONS FROM DRT TAPE                                       DRT00910
C                                                                       DRT00920
      CALL RFILE(ITAP58)                                                DRT00930
      CALL WREADW(ITAP58,DRTVER,2,1,IEND)                               DRT00940
      CALL WREADW(ITAP58,LBLDRT,26,IEND,IEND)                           DRT00950
      CALL WREADW(ITAP58,NBF,1,IEND,IEND)                               DRT00960
      CALL WREADW(ITAP58,NSYM,1,IEND,IEND)                              DRT00970
      CALL WREADW(ITAP58,NORBS,1,IEND,IEND)                             DRT00980
      CALL WREADW(ITAP58,NDUM,6,IEND,IEND)                              DRT00990
      CALL WREADW(ITAP58,NWKS,1,IEND,IEND)                              DRT01000
      CALL WREADW(ITAP58,NDUM,3,IEND,IEND)                              DRT01010
      CALL WREADW(ITAP58,ISYMRB,1,IEND,IEND)                            DRT01020
      CALL WREADW(ITAP58,NUMIJ,1,IEND,IEND)                             DRT01030
      CALL WREADW(ITAP58,NGROUP,1,IEND,IEND)                            DRT01040
      CALL WREADW(ITAP58,NUMINT,1,IEND,IEND)                            DRT01050
      CALL WREADW(ITAP58,NMAX,1,IEND,IEND)                              DRT01060
      IEND=IEND+3                                                       DRT01070
C                                                                       DRT01080
c*    WRITE(6,4) LBLDRT                                                 DRT01090
      WRITE(6,5) NBF,NSYM,NORBS,NWKS,NMAX,NUMIJ,NGROUP,ISYMRB,iprnt
C                                                                       DRT01120
      IC1=1                                                             DRT01130
      IA1=1                                                             DRT01140
      IC2=IC1+ISYMRB                                                    DRT01150
      IA2=IC2+IC2-1                                                     DRT01160
      IC3=IC2+ISYMRB                                                    DRT01170
      IA3=IC3+IC3-1                                                     DRT01180
      IC4=IC3+NUMIJ                                                     DRT01190
      IA4=IC4+IC4-1                                                     DRT01200
      ICSAVE=IC4+NUMIJ                                                  DRT01210
C                                                                       DRT01220
C.................KADD    LADD    IJADD   IJGRP...............          DRT01230
      CALL GETDRT(IA(IA1),IA(IA2),IA(IA3),IA(IA4),IEND,NGROUP)          DRT01240
C                                                                       DRT01250
C   FORM AN ARRAY TO RELATE DRT AND SO ORDERINGS                        DRT01260
      CALL Mwrit(IORDRT,5)        
      CALL Mread(IORSO,4)    
      ISUM=0                                                            DRT01280
      DO 103 I=1,NBF                                                    DRT01290
      II=IORDRT(I)                                                      DRT01300
      DO 102 J=1,NBF                                                    DRT01310
      JJ=IORSO(J)                                                       DRT01320
      IF(II.NE.JJ) GO TO 102                                            DRT01330
      ISUM=ISUM+1                                                       DRT01340
      IDRTSO(I)=J                                                       DRT01350
  102 CONTINUE                                                          DRT01360
  103 CONTINUE                                                          DRT01370
      IF(ISUM.EQ.NBF) GO TO 301                                         DRT01380
      WRITE(6,6) NBF,ISUM                                               DRT01390
      STOP                                                              DRT01400
  301 CONTINUE                                                          DRT01410
      IF(IPRNT.LE.4) GO TO 302                                          DRT01420
      WRITE(6,7)                                                        DRT01430
      DO 104 I=1,NBF                                                    DRT01440
  104 WRITE(6,8) I,IORDRT(I),IORSO(I),IDRTSO(I)                         DRT01450
C                                                                       DRT01460
C   READ IN MO INTEGRALS IN DRT ORDER AND REORDER THEM IN SORTED-SCF    DRT01470
C   ORDER                                                               DRT01480
C                                                                       DRT01490
  302 CONTINUE                                                          DRT01500
      IC5=ICSAVE                                                        DRT01530
      IC6=IC5+NUMIJ                                                     DRT01540
      IC7=IC6+NMAX                                                      DRT01550
      IA7=IC7+IC7-1                                                     DRT01560
      ICMAX=IC7+MAXBF2                                                  DRT01570
C................KADD    LADD    IJADD   IJGRP   HQ      BUFI....       DRT01580
      CALL GETIT(IA(IA1),IA(IA2),IA(IA3),IA(IA4),CC(IC5),CC(IC6),       DRT01590
C.............. .LBLIO   BUFO..................                         DRT01600
     1           IA(IA7),CC(IC7),ITAP52,ITAP35)                         DRT01610
C                                                                       DRT01620
C   READ IN PDM MATRICES IN DRT ORDER AND REORDER THEM IN SORTED-SCF    DRT01630
C   ORDER                                                               DRT01640
C                                                                       DRT01650
      IC5=ICSAVE                                                        DRT01700
      IC6=IC5+NUMIJ                                                     DRT01710
      IC7=IC6+NMAX                                                      DRT01720
      IA7=IC7+IC7-1                                                     DRT01730
      ICMAX=IC7+MAXBF2                                                  DRT01740
C................KADD    LADD    IJADD   IJGRP   HQ      BUFI....       DRT01750
      CALL GETIT(IA(IA1),IA(IA2),IA(IA3),IA(IA4),CC(IC5),CC(IC6),       DRT01760
C.............. .LBLIO   BUFO..................                         DRT01770
     1           IA(IA7),CC(IC7),ITAP53,ITAP56)                         DRT01780
C                                                                       DRT01790
  400 CONTINUE                                                          DRT01800
      CALL RCLOSE(MASTER,3)                                             DRT01810
      CALL RCLOSE(ITAP58,3)                                             DRT01820
      CALL TSTOP(6)                                                     DRT01830
      STOP                                                              DRT01840
      END                                                               DRT01850
