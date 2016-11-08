      subroutine fentry(cc,ia,maxcor)
C**********************************************************************
C     PROGRAM OPDM_SO                                                 *
C**********************************************************************
C   THIS PROGRAM REORDERS ONE-PDM IN DRT ORDER                        *
C   TO THOSE IN PITZER-SCF ORDER FOR CI AND MCSCF WAVEFUNCTIONS       *
C======================================================================
C  Modification for LINUX symstem
C  Subroutine Mwrit has been commetned out
C  on September 07, 2005 by Dr. Yukio Yamaguchi
C======================================================================
C**********************************************************************
C*  Last updated on January 09, 2003 by Dr. Yukio Yamaguchi           *
C**********************************************************************
C                                                                  
      IMPLICIT REAL*8 (A-H,O-Z)                                
      dimension cc(maxcor),ia(maxcor*2)
Ceeb For New Style Input
      integer errcod
      integer fwkclr, fwkadd, frdi
      DIMENSION LBLDRT(26),NDUM(20)                                     DRT00090
      COMMON/CALIF/LPARA(1024),APARA(1024)                              DRT00120
      COMMON/DIMES/ISYMRB,NUMIJ,NBF,NORBS,NMAX                          DRT00130
      COMMON/MFSEC/MASTER,NSECT                                         DRT00160
      COMMON/ORDCI/IDRSYM(256),IORDRT(256)                              DRT00170
      COMMON/ORDSO/IORSO(256),IDRTSO(256)                               DRT00180
      common/posit/iopci,iopmc
      COMMON/SIGNS/IOFF(256),IPRNT                                      DRT00190
      COMMON/TAPES/ITAP53,ITAP58   
    1 FORMAT(//,2X,' THE OPDM-SO PROGRAM  (Version of 09/07/05)',/,11x,
     1 'Written by Dr. Yukio Yamaguchi',///)
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
C     ITAP53 = CONTAINS THE ONE AND TWO PDM'S IN THE MO BASIS (DRT)     DRT00500
C     ITAP58 = CONTAINS THE DRT INFORMATION                             DRT00520
C                                                                       DRT00530
      ITAPE3=3                                                          DRT00540
      INPUT=5                                                           DRT00550
      ITAP53=53                                                         DRT00580
      ITAP58=58                                                         DRT00600
c
c=====================================================================
C    Followings are the positions of arrays in the master file
c=====================================================================
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
      write(3,*) ' This is a Chek file'
      write(3,*) '  maxcor  = ',maxcor
      write(3,*) '  master  = ',master
      write(3,*) '  nsect   = ',nsect
C
      errcod = fwkclr()
      errcod = fwkadd(5,'OPDM_SO')
      errcod = fwkadd(5,'DEFAULT')
      iprnt=0
      errcod = frdi('PRINT',iprnt)
C
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
      write(3,*) ' Now you are before getdrt'
C                                                                       DRT01220
C.................KADD    LADD    IJADD   IJGRP...............          DRT01230
      CALL GETDRT(IA(IA1),IA(IA2),IA(IA3),IA(IA4),IEND,NGROUP)          DRT01240
      write(3,*) ' Now you are after getdrt'
C                                                                       DRT01250
C   FORM AN ARRAY TO RELATE DRT AND SO ORDERINGS                        DRT01260
      write(3,*) ' Now you are before Mwrit'
cyy9-07-05 CALL Mwrit(IORDRT,5)        
      write(3,*) ' Now you are after Mwrit'
      write(3,*) ' Now you are before Mread'
      CALL Mread(IORSO,4)    
      write(3,*) ' Now you are after Mread'
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
C   READ IN PDM MATRICES IN DRT ORDER AND REORDER THEM IN SORTED-SCF    DRT01630
C   ORDER                                                               DRT01640
C                                                                       DRT01650
  302 CONTINUE                                                          DRT01500
      IC5=ICSAVE                                                        DRT01700
      IC6=IC5+NUMIJ                                                     DRT01710
      ICMAX=IC6+NMAX   
      write(3,*) ' Now you are before getit'
C................KADD    LADD    IJADD   IJGRP   HQ      BUFI....
      CALL GETIT(IA(IA1),IA(IA2),IA(IA3),IA(IA4),CC(IC5),CC(IC6))
      write(3,*) ' Now you are after getit'
C
C                                                                       DRT01790
  400 CONTINUE                                                          DRT01800
      CALL RCLOSE(MASTER,3)                                             DRT01810
      CALL RCLOSE(ITAP58,3)                                             DRT01810
      CALL TSTOP(6)                                                     DRT01830
      STOP                                                              DRT01840
      END                                                               DRT01850
