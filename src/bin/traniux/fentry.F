      subroutine fentry(cc,ia,maxcor)
C**********************************************************************
C* 09/14/92      PROGRAM TRANIUX                                      *
C**********************************************************************
C   THE TRANSFORMATION PROGRAM                                        *
C   THIS PROGRAM TRANSFORMS INTEGRALS FROM PITZER'S SO TO MO BASIS    *
C   FOR ONLY ``OCCUPIED'' AND ``INDEPENDENT PAIRS''                   *
C**********************************************************************
C   Last updated on January 14, 1995 by Evan Bolton                   *
C   Reason:  Forgot to claim lib parse functions as integer           *
C**********************************************************************
C   Last updated on January 10, 1995 by Evan Bolton                   *
C   Reason:  To put in PSI 2.0 input parsing.  See comments Ceeb      *
C**********************************************************************
C   Last updated on September 10, 1993 by Dr. Yukio Yamaguchi         *
C   ``Maxbas'' has been reset to 256                                  *
c   ``Nbfao'' and ``Nbatri'' have been eliminated to avoid confusion  *
C**********************************************************************
C   Last updated on August 03, 1993 by Dr. Yukio Yamaguchi            *
C   ``call tstart'' , ''call initmf''  and ``maxbas'' were moved to   *
c   after ``call psinit''                                             *
C**********************************************************************
C   Last updated on September 14, 1992 by Dr. Yukio Yamaguchi         *
C   Modification for the Unix Version - Compatible to Psi-2           *
C**********************************************************************
C*  LAST UPDATED ON OCTOBER 29, 1989 BY YUKIO YAMAGUCHI               *
C**********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension cc(maxcor),ia(maxcor*2)
      DIMENSION NDEG(10),LORD(256),IPZSYM(256)
Ceeb Needed For Psi2.0 Input
      integer errcod, inold
      integer fwkclr, fwkadd, frdi
Ceeb
      COMMON/BASIS/NBASIS,NTRI
      COMMON/CALIF/LPARA(1024),APARA(1024)
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4
      COMMON/MFSEC/MASTER,NSECT
      COMMON/SDATA/ISMAT(8,8)
      COMMON/SIGNS/IOFF(257),IPRNT
      COMMON/SYMMS/MOSYM(256),MOVIR(256)
      COMMON/TAPES/ITAP35,ITAP36
    1 FORMAT(//,2X,' TRANIUX:  A transformation program for ',
     1 '''Occupied'' and ''Independent Pairs''  (Version of 07/30/03)',
     2 /,13X,'Alterations of TRANS for Stability Analysis by Dr. ',
     3 'Yukio Yamaguchi',//,13X, 'Modifications for PSI 2.0 use by ',
     4 'Evan Bolton'//)
    2 FORMAT(5I5)
    3 FORMAT(//,2X,' PARAMETERS'/
     1 2X,' NBASIS = ',I8/
     2 2X,' NTRI   = ',I8/
     3 2X,' NSYM   = ',I8/
     4 2X,' NST    = ',I8/
     5 2X,' IOCC   = ',I8/
     6 2X,' JOCC   = ',I8/
     7 2X,' KOCC   = ',I8/
     8 2X,' IPRNT  = ',I8/)
    4 FORMAT(//,2X,' SO-MO EIGENVECTOR MATRIX'/)
    5 FORMAT(//,2X,' MO ORDERING FOR SORTED SCF'//2X,20I5/)
    6 FORMAT(//,2X,' NBASIS EXCEEDS MAXBAS, NBASIS = ',I8,
     1             ' MAXBAS = ',I8/)
    7 FORMAT(//,2X,' REQUIRED MEMORY, ICMAX = ',I8,' MAXCOR = ',I8/)
    8 FORMAT(//,2X,' REQUIRED MEMORY EXCEEDS MAXCOR',' ICMAX = ',I8,
     1             ' MAXCOR = ',I8/)
C
C   ITAP36 = TWO ELECTRON SO INTEGRALS IN PITZER'S SCF
C   ITAP35 = TWO ELECTRON MO INTEGRALS IN PITZER'S SCF
C
      ITAPE3=3
      ITAPE6=6
      INPUT=5
      ITAP35=35
      ITAP36=36
C

      call psinit('APPEND')
      CALL TSTART(6)
      CALL INITMF(1)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C:::  DEFINE MACHINE DEPENDENT PARAMETERS
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      MAXBAS=256

      iprnt = 0

      CALL LOCATE(INPUT,'# TRANS ###',IERR)

Ceeb  Check to see if using old style input
      if (ierr.ne.0) inold = 0
      if (ierr.eq.0) inold = 1

      if (inold.eq.1) then
        write(6,*) 'WARNING: using old style ''# TRANS ##'''
      endif
Ceeb

C
      WRITE(6,1)
      WRITE(3,1)

Ceeb  Read in Old style input
      if (inold.eq.1) then
        READ(5,2) IPRNT
      endif
Ceeb  Finished reading in Old Style input

Ceeb  Read in New Style input
      if (inold.eq.0) then
        errcod = fwkclr()
        errcod = fwkadd(5,'TRANS')
        errcod = fwkadd(5,'DEFAULT')
        errcod = frdi('PRINT',iprnt)
      endif
Ceeb  Finished reading in New Style input

C
      IOFF(1)=0
      DO 101 I=1,256
      IOFF(I+1)=IOFF(I)+I
  101 CONTINUE
C
C   READ IN PARAMETERS
      NBASIS=LPARA(3)
      NTRI=LPARA(4)
      IOCC=LPARA(7)
      KOCC=LPARA(8)
      JOCC=LPARA(9)
      NSYM=LPARA(13)
      NST=LPARA(30)
      WRITE(6,3) NBASIS,NTRI,NSYM,NST,IOCC,JOCC,KOCC,IPRNT
      if(nbasis.gt.maxbas) then
        write(6,*) '  Nbasis cxceeds Maxbas'
        write(6,*) '  Nbasis = ',nbasis,' Maxbas = ',maxbas
      end if
C
C   READ IN SO-MO EIGENVECTOR MATRIX IN PITZER'S SCF
      IC1=1
      CALL MREAD(CC(IC1),19)
      IF(IPRNT.LE.2) GO TO 301
      WRITE(6,4)
      call matout(cc(ic1),nbasis,nbasis,nbasis,nbasis,6)
C
C   ASSIGN SYMMETRY TYPE FOR EACH MOLECULAR ORBITAL
  301 CONTINUE
      DO 102 I=1,NST
      NDEG(I)=LPARA(I+30)
  102 CONTINUE
C*    WRITE(6,*) (NDEG(I),I=1,NST)
      WRITE(3,*) (NDEG(I),I=1,NST)
      CALL MREAD(LORD,4)
C*    WRITE(6,5) (LORD(I),I=1,NBASIS)
      WRITE(3,5) (LORD(I),I=1,NBASIS)
      CALL IZERO(IPZSYM,nbasis)
      MM=0
      DO 104 I=1,NST
      NN=NDEG(I)
      IF(NN.EQ.0) GO TO 104
      DO 103 J=1,NN
      MM=MM+1
      IPZSYM(MM)=I
  103 CONTINUE
  104 CONTINUE
C*    WRITE(6,*) (IPZSYM(I),I=1,NBASIS)
      WRITE(3,*) (IPZSYM(I),I=1,NBASIS)
      CALL IZERO(MOSYM,nbasis)
      DO 105 I=1,NBASIS
      MM=LORD(I)
      NN=IPZSYM(MM)
      MOSYM(I)=NN
  105 CONTINUE
C*    WRITE(6,*) (MOSYM(I),I=1,NBASIS)
      WRITE(3,*) (MOSYM(I),I=1,NBASIS)
C
C   DISTINGUISH BETWEEN OCCUPIED AND VIRTUAL ORBITALS
      CALL IZERO(MOVIR,nbasis)
      DO 106 I=1,NBASIS
      IF(I.LE.JOCC) GO TO 106
      MOVIR(I)=1
  106 CONTINUE
C*    WRITE(6,*) (MOVIR(I),I=1,NBASIS)
      WRITE(3,*) (MOVIR(I),I=1,NBASIS)
C
C   TRANSFORM SO INTEGRALS INTO MO INTEGRALS IN SORTED SCF
C::::::::::::::::::::::::::::
C:::ONE ELECTRON INTEGRALS:::
C::::::::::::::::::::::::::::
      WRITE(3,21)
   21 FORMAT(//,2X,' NOW YOU ARE IN IJMO'/)
      IC1=1
      IC2=IC1+NBASIS*NBASIS
      IC3=IC2+NTRI
      IC4=IC3+NTRI
      IC5=IC4+NBASIS*NBASIS
      ICMAX=IC5+NBASIS*NBASIS
C...............ESO     HSO     HMO     U       T.......
      CALL IJMO(CC(IC1),CC(IC2),CC(IC3),CC(IC4),CC(IC5))
C
C::::::::::::::::::::::::::::::::::::::::::
C:::TWO ELECTRON INTEGRALS WITH SYMMETRY:::
C::::::::::::::::::::::::::::::::::::::::::
      WRITE(3,22)
   22 FORMAT(//,2X,' NOW YOU ARE IN SMIJKL'/)
      IC1=1
      IC2=IC1+NBASIS*NBASIS
      IC3=IC2+NBASIS
      IC4=IC3+NTRI
      IC5=IC4+NTRI*NBASIS
      IA5=IC5+IC5-1
      IC6=IC5+MAXBF2
      IA6=IC6+IC6-1
      ICMAX=IC6+MAXBF2
      IF(ICMAX.GT.MAXCOR) GO TO 399
      WRITE(3,7) ICMAX,MAXCOR
C.................EV      IJKU    IJTU    ISTU    LBLI    BUFI....
      CALL SMIJKL(CC(IC1),CC(IC2),CC(IC3),CC(IC4),IA(IA5),CC(IC5),
C.................LBLO    BUFO....
     1            IA(IA6),CC(IC6))
      GO TO 400
C
  399 CONTINUE
      WRITE(6,8) ICMAX,MAXCOR
      STOP
  400 CONTINUE
C
      CALL RCLOSE(MASTER,3)
      CALL RCLOSE(ITAP35,3)
C
      CALL TSTOP(3)
      CALL TSTOP(6)
C
      STOP
      END
