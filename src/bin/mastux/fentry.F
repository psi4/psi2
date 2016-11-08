      subroutine fentry(cc,ia,maxcor)
C**********************************************************************
C 09/14/92       PROGRAM MASTUX
C**********************************************************************
C   THIS PROGRAM FORMS THE MASTER FILE
C   USING THE PITZER'S SCF INFORMATION
C**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
C**********************************************************************
C   Last updated on Febuary 03, 1995 by Evan Bolton                   *
C   Reason:  Addition of CASSCF WFN type required the hard wiring of  *
C              of the SCFTYP, DERTYP, and idel34.  Those variables    *
C              must be GRSCF, SECOND, and 0, respectively.            *
C              See comments labeled CeebCAS.                          *
C**********************************************************************
C   Last updated on January 14, 1995 by Evan Bolton                   *
C   Reason:  Alter routine which calculates the number of open shells.*
C		Major problem if degeneracy.  Correction as per YY's  *
C		suggestions.  See comments Ceeb			      *
C**********************************************************************
C   Last updated on January 10, 1995 by Evan Bolton                   *
C   Reason:  To put in PSI 2.0 input parsing.  See comments Ceeb      *
C**********************************************************************
C   Last updated on October 18, 1993 by Dr. Yukio Yamaguchi           *
C*  Modification for subroutine ``sort'' for TCSCF wavefunctions      *
C**********************************************************************
C   Last updated on July 12, 1993 by Dr. Yukio Yamaguchi              *
C*  Statements ``call tstart'' and ``call initmf''                    *
c*  have been moved to below ``call psinit'' line                     *
C**********************************************************************
C   Last updated on March 17, 1993 by Dr. Yukio Yamaguchi             *
C   Modification for sorting open-shell energies                      *
c   in Subroutine Sort                                                *
C**********************************************************************
C   Last updated on September 17, 1992 by Dr. Yukio Yamaguchi         *
C   Modification for the Unix Version - Compatible to Psi-2           *
C**********************************************************************
C*  LAST UPDATED ON OCTOBER 26, 1989 BY YUKIO YAMAGUCHI               *
C**********************************************************************
C*  BY:  TRACY HAMILTON                                               *
C*  DATE:  JANUARY 16, 1989                                           *
C*  REASON: CHANGE EAO DIMENSION TO AO,SO, NEED TO CHECK ALL DIMENSION*
C*          ADD ALP AND BET TO MASTER FILE, DRT ORDERED AO-MO VECTOR  *
C**********************************************************************
C*  BY:  RICHARD REMINGTON                              cMachineSize  *
C*  DATE:  MARCH  12,  1988                                           *
C*  REASON: DECREASE CORE TO RUN IN 7MB ON 9370                       *
C**********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension cc(maxcor),ia(maxcor*2)
      INTEGER BLABEL(20)
      INTEGER ND(10),MNL(142),MS(142),MTYPE(2,20),ITYP(20)
      DIMENSION A30(200)
      CHARACTER*8 SCFTYP,CITYP,DERTYP
      CHARACTER*8 TAPE,OPTYPE
      LOGICAL TAPE36
      COMMON/BASIS/NBASIS,NTRI,NST,NSYM,NBFAO,NBFSO,NBATRI
      COMMON/CALIF/LPARA(1024),APARA(1024)
      COMMON/COUPL/ALPC(15),BETC(15),ALPT(15),BETT(15)
      COMMON/ENRGY/ESCF,ENUC
      COMMON/FLAGS/MCONST,MPOINT,NCALCS,MXCOEF
      COMMON/FUNCS/NTYPES,NATOM,N3N
      COMMON/GRSCF/FOCC(10),NSORB(10)
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4
      COMMON/LOCAT/LOCS(1024)
      COMMON/MFSEC/MFILE,NSECT
      COMMON/OPEN1/NUNIQ,LOPN(5),MOPN(5,5)
      COMMON/OPEN2/GOCC(10)
      COMMON/SIGNS/IOFF(256),IPRNT
      COMMON/SYMMS/IDEG(10),NDEG(10),LSYM(10)
      COMMON/TAPES/ITAP34,ITAP36
      COMMON/TEMPS/I30(400)
      COMMON/TYPES/SCFTYP,CITYP,DERTYP
      COMMON/CI101/IOPEN,IOCC,JOCC,KOCC
      EQUIVALENCE (I30,A30)
ceeb  For Psi 2.0 input
#include <error.h>
      integer frdc, frdi, frdf, fwkclr, fwkadd, frdboo, frdpos
      integer fcnt, frdia1
      integer errcod, isocc(8)
      logical oldinp
Ceeb
    1 FORMAT(//,2X,' MASTUX: A MASTER FILE program ',
     1 ' (Version of 01/14/95)',/,11X,
     2 'Alterations of MASTER for Stability Analysis by Dr. Yukio ',
     3 'Yamaguchi',//,11X, 'Modifications for PSI 2.0 use by ',
     4 'Evan Bolton'///)
    2 FORMAT(3(A8,2X))
    3 FORMAT(5I5)
    4 FORMAT(2X,' BLABEL = ',20A4/)
    5 FORMAT(2X,' ENUC   = ',F15.8/)
    6 FORMAT(2X,' NST    = ',I5/2X,' NDEG   =',8I5/)
    7 FORMAT(A8,2X,I5)
    8 FORMAT(F10.5,10I5)
    9 FORMAT(A8,2X)
   10 FORMAT(2X,' NBASIS = ',I8/
     1       2X,' NTRI   = ',I8/
     2       2X,' NST    = ',I8/
     3       2X,' NSYM   = ',I8/
     4       2X,' NBFAO  = ',I8/
     5       2X,' NBFSO  = ',I8/
     6       2X,' NBATRI = ',I8/
     7       2X,' NATOM  = ',I8/
     8       2X,' N3N    = ',I8/
     9       2X,' IOPEN  = ',I8/
     A       2X,' IPRNT  = ',I8/
     B       2X,' SCFTYP = ',A8/
     C       2X,' CITYP  = ',A8/
     D       2X,' DERTYP = ',A8/
     E       2X,' TAPE   = ',A8/
     F       2X,' idel34 = ',i8)
   11 FORMAT(//,2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)
   12 FORMAT(//,2X,' REQUIRED MEMORY EXCEEDS MAXCOR '/
     1          2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)
C
C   ITAP34 = SO INTEGRALS
C   ITAP36 = TWO ELECTRON SO INTEGRALS
C
      call psinit('APPEND')
      CALL TSTART(6)
      CALL INITMF(0)

      ITAPE3=3
      ITAPE6=6
      ITAP30=30
      ITAP34=34
      ITAP36=36
      TAPE36=.TRUE.
      INPUT=5
C
      WRITE(6,1)
      WRITE(3,1)
      write(6,*) '  maxcor = ',maxcor

      CALL LOCATE(INPUT,'# MASTER #',IERR)

      idel34 = 1

Ceeb
      if (ierr.eq.0) then
        oldinp=.true.
        READ(5,2) SCFTYP,CITYP,DERTYP
        READ(5,3) IPRNT,idel34
      else
        errcod = fwkclr()
        errcod = fwkadd(5,'MASTER')
        errcod = fwkadd(5,'DEFAULT')

        oldinp=.false.

        cityp  = 'SCF     '
        errcod = frdc('WFN',cityp)

        dertyp = 'SECOND  '
        if(cityp.eq.'CI      ') dertyp = 'FIRST   '
        errcod = frdc('DERTYPE',dertyp)

        optype = 'NONE    '
        errcod = frdc('OPENTYPE',optype)

        iprnt = 0
        errcod = frdi('PRINT',iprnt)

        idel34 = 1
        if(cityp.eq.'CI      ') idel34 = 0
        errcod = frdboo('DELETE34',idel34)

        if (cityp.eq.'GVBCI   '.or.optype.eq.'TWOCON  ') then
          scftyp = 'TCSCF   '
        else if(cityp.eq.'CI      '.or.optype.eq.'HIGHSPIN'.or.
     +                                 optype.eq.'SINGLET ') then
          scftyp = 'GRSCF   '
        else
          scftyp = 'GRSCF   '
        endif

CeebCAS Begin CASSCF Fix
        if (cityp.eq.'CASSCF  ') then
          scftyp = 'GRSCF   '
          dertyp = 'SECOND  '
          idel34 = 0
        endif
CeebCAS END CASSCF Fix
      endif

Ceeb

C
      IOFF(1)=0
      DO 101 I=1,255
  101 IOFF(I+1)=IOFF(I)+I
C
C   READ IN INTEGRAL INFORMATION FROM 34 FILE
      CALL RFILE(ITAP34)
      CALL WREADW(ITAP34,NKIND,1,1,IEND)
c*    write(6,*) ' nkind = ',nkind
      CALL WREADW(ITAP34,ICLOS,1,IEND,IEND)
c*    write(6,*) ' iclos = ',iclos
      CALL WREADW(ITAP34,BLABEL,20,IEND,IEND)
c*    write(6,*) ' blabel = ',blabel
      CALL WREADW(ITAP34,ENUC,2,IEND,IEND)
c*    write(6,*) ' enuc = ',enuc
      CALL WREADW(ITAP34,NST,1,IEND,IEND)
c*    write(6,*) ' nst = ',nst
      CALL WREADW(ITAP34,ND,NST,IEND,IEND)
c*    write(6,*) ' nd = ',nd
      CALL WREADW(ITAP34,ITYP,NST,IEND,IEND)
c*    write(6,*) ' ityp = ',ityp
      CALL WREADW(ITAP34,NDEG,NST,IEND,IEND)
c*    write(6,*) ' ndeg = ',ndeg
      CALL WREADW(ITAP34,NS,1,IEND,IEND)
c*    write(6,*) ' ns = ',ns
      CALL WREADW(ITAP34,MTYPE,2*NS,IEND,IEND)
      CALL WREADW(ITAP34,LAPU,1,IEND,IEND)
      CALL WREADW(ITAP34,MS,LAPU,IEND,IEND)
      CALL WREADW(ITAP34,MNL,LAPU,IEND,IPARA)
C     WRITE(6,30) IPARA
C  30 FORMAT(2X,' IPARA = ',I10/)
      ISADR=I2SEC(IPARA)+1
      CALL RSETSA(ITAP34,ISADR)
C
      WRITE(6,4) BLABEL
C     WRITE(3,4) BLABEL
C     WRITE(6,5) ENUC
      WRITE(3,5) ENUC
C     WRITE(6,6) NST,(NDEG(I),I=1,NST)
      WRITE(3,6) NST,(NDEG(I),I=1,NST)
C
      NSYM=0
      NBASIS=0
      DO 102 I=1,NST
      IF(NDEG(I).GT.0) NSYM=NSYM+1
      NBASIS=NBASIS+NDEG(I)
  102 CONTINUE
      IDEG(1)=0
      LSYM(1)=0
      IF(NSYM.EQ.1) GO TO 203
      II=0
      DO 103 I=2,NST
      IF(NDEG(I).LE.0) GO TO 202
  201 II=II+1
      NN=NDEG(II)
      IF(NN.EQ.0) GO TO 201
      IDEG(I)=IDEG(I-1)+NN
      LSYM(I)=LSYM(I-1)+IOFF(NN+1)
      GO TO 103
  202 LSYM(I)=LSYM(I-1)
      IDEG(I)=IDEG(I-1)
  103 CONTINUE
  203 NTRI=IOFF(NBASIS+1)
C
C   READ IN PARAMETERS FROM TAPE30
      CALL RFILE(ITAP30)
      CALL WREADW(ITAP30,I30,200,101,JUNK)
      IEND=I30(1)
      MPOINT=I30(2)
      MCONST=I30(3)
      MCALCS=I30(4)
      NCALCS=I30(5)
      NBFSO=I30(18)
      NATOM=I30(19)
      NBFAO=I30(22)
      NBATRI=IOFF(NBFAO+1)
      MXCOEF=I30(42)
      IOPEN=I30(43)
      N3N=NATOM*3
C
C   INPUT SHELL INFORMATION ONLY FOR A GRSCF OPEN-SHELL CASE
C   UNNECESSARY FOR A TCSCF CASE

Ceeb
      if(oldinp) then
        CALL LOCATE(INPUT,'# MASTER #',IERR)
        READ(5,2) SCFTYP,CITYP,DERTYP
        READ(5,3) IPRNT,idel34
        IF(SCFTYP.NE.'GRSCF   ') GO TO 205
        IF(IOPEN.EQ.0) GO TO 205
        READ(5,7) OPTYPE,NUNIQ
        DO 105 I=1,NUNIQ
        READ(5,8) GOCC(I),LL,(MOPN(I,J),J=1,LL)
        LOPN(I)=LL
  105   CONTINUE
  205   CONTINUE
        READ(5,9) TAPE
        IF(TAPE.EQ.'NOTAPE36') TAPE36=.FALSE.
        IF(TAPE.EQ.'NOFILE36') TAPE36=.FALSE.
        IF(TAPE36) TAPE='  TAPE36'
      else
        if(optype.eq.'SINGLET ') then
          nuniq=2
          gocc(1)=1.0d0
          gocc(2)=1.0d0
          lopn(1)=1
          lopn(2)=1
          mopn(1,1) = 1
          mopn(2,1) = 2
        else if(optype.eq.'HIGHSPIN') then
          nuniq=1
          gocc(1)=1.0d0

Ceeb	Calculates the number of openshells properly.

          errcod = fwkclr()
          errcod = fwkadd(5,'MASTER')
          errcod = fwkadd(5,'SCF')
          errcod = fwkadd(5,'DEFAULT')

          errcod = fcnt('SOCC',nsocc)

          if (errcod.ne.EOK) then
            write(6,*) 'ERROR: Problems reading the ''SOCC'' array'
            write(6,*) ' '
            call mabort
          endif

          if (nsocc.gt.8) then
            write(6,*) 'ERROR: ''SOCC'' array is too large!'
            write(6,*) ' '
            call mabort
          endif

          errcod = frdia1('SOCC',isocc,nsocc)

          if (errcod.ne.EOK) then
            write(6,*) 'ERROR: Problems reading the ''SOCC'' array'
            write(6,*) ' '
            call mabort
          endif

          do 5001 i=1,nsocc
 5001       insocc = insocc + isocc(i)

Ceeb          fmm = dsqrt(2.0d0*dfloat(iopen)+0.25d0) - 0.5d0
Ceeb          lopn(1) = dint(fmm)

          lopn(1) = insocc
          do 5000 i=1,lopn(1)
 5000       mopn(1,i) = i

        else if(optype.eq.'SPECIAL ') then
          write(6,*) 'cannot do special cases'
          call qabort
        endif
      endif

Ceeb

      WRITE(6,10) NBASIS,NTRI,NST,NSYM,NBFAO,NBFSO,NBATRI,NATOM,N3N,
     1            IOPEN,IPRNT,SCFTYP,CITYP,DERTYP,TAPE,idel34

CeebCAS  Now that we have shown that the CASSCF is the type lets 
CeebCAS   switch it to SCF so it will properly set things up as
CeebCAS   per a SCF SECOND type.
      if (cityp.eq.'CASSCF  ') then
        cityp = 'SCF     '
      endif
CeebCAS   End CASSCF Fix

C
C   FORM MASTER FILE
      WRITE(3,20)
   20 FORMAT(//,2X,' NOW YOU ARE IN FORMMS'/)
      CALL FORMMS(CC,IA,MAXCOR)
C
C   CREATE SO INTEGRAL FILES
      WRITE(3,21)
   21 FORMAT(//,2X,' NOW YOU ARE IN RDONE'/)
      IBUFHF=4471
      IC1=1
      IC2=IC1+NTRI
      IC3=IC2+NTRI
      IC4=IC3+NTRI
      IC5=IC4+NTRI
      IA5=IC5+IC5-1
      ICMAX=IC5+IBUFHF
      WRITE(3,11) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 399
C................S       H       T       V       LBLI    STVI....
      CALL RDONE(CC(IC1),CC(IC2),CC(IC3),CC(IC4),IA(IA5),CC(IC5),
C.......................
     1           IBUFHF)
      CALL MWRIT(CC(IC1),14)
      CALL MWRIT(CC(IC2),15)
      CALL MWRIT(CC(IC1),21)
      CALL MWRIT(CC(IC2),22)
C
      IF(.NOT.TAPE36) GO TO 300
      WRITE(3,22)
   22 FORMAT(//,2X,' NOW YOU ARE IN RDTWO'/)
      IBUFHF=4471
      IC1=1
      IA1=1
      IC2=IC1+IBUFHF
      IA2=IC2+IC2-1
C................LBLI    PKI     LBLO    BUFO...........
      CALL RDTWO(IA(IA1),CC(IC1),IA(IA2),CC(IC2),IBUFHF)
C
C   READ IN SCF INFORMATION FROM TAPE30
  300 CONTINUE
      WRITE(3,23)
   23 FORMAT(//,2X,' NOW YOU ARE IN SCFINP'/)
      IC1=1
      IC2=IC1+NBASIS
      IC3=IC2+NBASIS
      IC4=IC3+NBFSO*NBASIS
      IC5=IC4+NBFAO*NBFSO
      IC6=IC5+NBFAO*NBASIS
      IC7=IC6+NBFAO*NBFAO
      IA7=IC7+IC7-1
      IC8=IC7+NATOM
      IA8=IC8+IC8-1
      ICMAX=IC8+NATOM*3
      WRITE(3,11) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 399
C.................EIG     OCC     ESO     UAS     EAO     EE......
      CALL SCFINP(CC(IC1),CC(IC2),CC(IC3),CC(IC4),CC(IC5),CC(IC6),
C.................CHGS    XX......
     1            CC(IC7),CC(IC8))
      CALL MWRIT(CC(IC1),10)
      CALL MWRIT(CC(IC2),11)
      CALL MWRIT(CC(IC3),12)
      CALL MWRIT(CC(IC4),9)
      CALL MWRIT(CC(IC5),13)
      CALL MWRIT(CC(IC7),8)
C
C    SORT EIGENVECTORS IF NECESSARY
      WRITE(3,24)
   24 FORMAT(//,2X,' NOW YOU ARE IN SORT'/)
      CALL MREAD(LPARA,2)
      IC1=1
      IC2=IC1+NBASIS
      IC3=IC2+NBASIS
      IC4=IC3+NBFSO*NBASIS
      IC5=IC4+NBFAO*NBFSO
      IC6=IC5+NBFAO*NBASIS
      IC7=IC6+NBFAO*NBFAO
      IC8=IC7+NBASIS
      IC9=IC8+NBASIS
      IA9=IC9+IC9-1
      IC10=IC9+NBASIS
      IA10=IC10+IC10-1
      IC11=IC10+NBASIS
      IA11=IC11+IC11-1
      IC12=IC11+NBASIS
      IA12=IC12+IC12-1
      ICMAX=IC12+NBASIS
      WRITE(3,11) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 399
C...............EIG     OCC     ESO     UAS     EAO     EE......
      CALL SORT(CC(IC1),CC(IC2),CC(IC3),CC(IC4),CC(IC5),CC(IC6),
C...............EGN     OCN.....
     1          CC(IC7),CC(IC8),
C...............KORD    LORD     MORD     NORD.....
     2          IA(IA9),IA(IA10),IA(IA11),IA(IA12))
      CALL MWRIT(CC(IC1),17)
      CALL MWRIT(CC(IC2),18)
      CALL MWRIT(CC(IC3),19)
      CALL MWRIT(CC(IC5),20)
      CALL MWRIT(LPARA,2)
      CALL MWRIT(I30,3)
      CALL MWRIT(IA(IA10),4)
C
C    FORM LAGRANGIAN AND DENSITY MATRICES FOR CLSCF
      IF(SCFTYP.NE.'CLSCF   ') GO TO 206
      WRITE(3,25)
   25 FORMAT(//,2X,' NOW YOU ARE IN SCFCL'/)
      IC1=1
      IC2=IC1+NBASIS
      IC3=IC2+NBASIS
      IC4=IC3+NBFSO*NBASIS
      IC5=IC4+NBFAO*NBFSO
      IC6=IC5+NBFAO*NBASIS
      IC7=IC6+NBATRI
      IC8=IC7+NBATRI
      ICMAX=IC8+NBATRI
      WRITE(3,11) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 399
C................EIG     OCC     EAO     WM      WA      DC......
      CALL SCFCL(CC(IC1),CC(IC2),CC(IC5),CC(IC6),CC(IC7),CC(IC8))
      GO TO 210
C
C   FORM HIGH-SPIN OPEN-SHELL MO
  206 IF(SCFTYP.NE.'HSOP    ') GO TO 207
      WRITE(3,26)
   26 FORMAT(//,2X,' NOW YOU ARE IN SCFHS'/)
      IC1=1
      IC2=IC1+NBASIS
      IC3=IC2+NBASIS
      IC4=IC3+NBFSO*NBASIS
      IC5=IC4+NBFAO*NBFSO
      IC6=IC5+NBFAO*NBASIS
      IC7=IC6+NBFAO*NBFAO
      IC8=IC7+NBFAO*NBFAO
      IC9=IC8+NBATRI
      IC10=IC9+NBATRI
      IC11=IC10+NBATRI
      IC12=IC11+NBATRI
      IC13=IC12+NBATRI
      IC14=IC13+NBATRI
      IC15=IC14+NBATRI
      IC16=IC15+NBATRI
      IC17=IC16+NBATRI
      IC18=IC17+NBFAO*NBFAO
      IC19=IC18+NBFAO*NBFAO
      IA19=IC19+IC19-1
      IC20=IC19+MAXBF2
      ICMAX=IC20+NBFAO*4
      WRITE(3,11) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 399
C................EIG     OCC     ESO     UAS     EAO     U.......
      CALL SCFHS(CC(IC1),CC(IC2),CC(IC3),CC(IC4),CC(IC5),CC(IC6),
C................T       H       DC      DO       FC       FF.......
     1           CC(IC7),CC(IC8),CC(IC9),CC(IC10),CC(IC11),CC(IC12),
C................FK       FT       WM       WA       EE       EZ.......
     2           CC(IC13),CC(IC14),CC(IC15),CC(IC16),CC(IC17),CC(IC18),
C................LBLI     BUFI     CC.......
     3           IA(IA19),CC(IC19),CC(IC20))
      CALL MWRIT(CC(IC1),17)
      CALL MWRIT(CC(IC2),18)
      CALL MWRIT(CC(IC3),19)
      CALL MWRIT(CC(IC5),20)
      GO TO 210
C
C   FORM GENERALIZED OPEN-SHELL MO
  207 IF(SCFTYP.NE.'GRSCF   '.AND.SCFTYP.NE.'TCSCF   ') GO TO 210
      IF(CITYP.EQ.'SCF     '.AND.DERTYP.EQ.'FIRST   ') GO TO 208
      IF(CITYP.EQ.'SCF     '.AND.DERTYP.EQ.'SECOND  ') GO TO 208
      WRITE(3,27)
   27 FORMAT(//,2X,' NOW YOU ARE IN SCFGR'/)
      IC1=1
      IC2=IC1+NBASIS
      IC3=IC2+NBASIS
      IC4=IC3+NBFSO*NBASIS
      IC5=IC4+NBFAO*NBFSO
      IC6=IC5+NBFAO*NBASIS
      IC7=IC6+NBFAO*NBFAO
      IC8=IC7+NBFAO*NBFAO
      IC9=IC8+NBATRI
      IC10=IC9+NBATRI
      IC11=IC10+NBATRI
      IC12=IC11+NBATRI
      IC13=IC12+NBATRI
      IC14=IC13+NBFAO*NBFAO
      IC15=IC14+NBFAO*NBFAO
      IA15=IC15+IC15-1
      IC16=IC15+MAXBF2
      ICMAX=IC16+NBFAO*4
      WRITE(3,11) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 399
C................EIG     OCC     ESO     UAS     EAO     U.......
      CALL SCFGR(CC(IC1),CC(IC2),CC(IC3),CC(IC4),CC(IC5),CC(IC6),
C................T       H       DC      FC       FM       FF.......
     1           CC(IC7),CC(IC8),CC(IC9),CC(IC10),CC(IC11),CC(IC12),
C................EE       EZ       LBLI     BUFI     CC.......
     2           CC(IC13),CC(IC14),IA(IA15),CC(IC15),CC(IC16))
      CALL MWRIT(CC(IC1),17)
      CALL MWRIT(CC(IC2),18)
      CALL MWRIT(CC(IC3),19)
      CALL MWRIT(CC(IC5),20)
C
C   FORM ALPA AND BETA MATRICES
  208 CONTINUE
      WRITE(3,28)
   28 FORMAT(//,2X,' NOW YOU ARE IN ABMAT'/)
      IC1=1
      IC2=IC1+NBASIS
      IC3=IC2+NBASIS
      IC4=IC3+NBFSO*NBASIS
      IC5=IC4+NBFAO*NBFSO
      IC6=IC5+NBFAO*NBASIS
      IC7=IC6+NBFAO*NBFAO
      IC8=IC7+NBFAO*NBFAO
      IC9=IC8+NBATRI
      IC10=IC9+NBATRI
      IC11=IC10+NBATRI
      IC12=IC11+NBATRI
      IC13=IC12+NBATRI
      IA13=IC13+IC13-1
      ICMAX=IC13+NBASIS
      WRITE(3,11) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 399
C................ALPA    BETA     ALP      BET      LORD.....
      CALL ABMAT(CC(IC9),CC(IC10),CC(IC11),CC(IC12),IA(IA13))
C
C   FORM LAGRANGIAN MATRICES
      WRITE(3,29)
   29 FORMAT(//,2X,' NOW YOU ARE IN ZETMAT'/)
      IC1=1
      IC2=IC1+NBASIS
      IC3=IC2+NBASIS
      IC4=IC3+NBFSO*NBASIS
      IC5=IC4+NBFAO*NBFSO
      IC6=IC5+NBFAO*NBASIS
      IC7=IC6+NBFAO*NBFAO
      IC8=IC7+NBFAO*NBFAO
      IC9=IC8+NBATRI
      IC10=IC9+NBATRI
      IC11=IC10+NBATRI
      IC12=IC11+NBATRI*NTYPES
      IC13=IC12+NBATRI*NTYPES
      IC14=IC13+NBATRI*NTYPES
      IC15=IC14+NBATRI*NTYPES
      IC16=IC15+NBATRI
      IC17=IC16+NBATRI
      IC18=IC17+NBATRI
      IC19=IC18+NBATRI
      IA19=IC19+IC19-1
      ICMAX=IC19+MAXBF2
      WRITE(3,11) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 399
C.................EIG     OCC     ESO     EAO     U       T.......
      CALL ZETMAT(CC(IC1),CC(IC2),CC(IC3),CC(IC5),CC(IC6),CC(IC7),
C.................H       ALPA    BETA     DA       DB       ZT.......
     1            CC(IC8),CC(IC9),CC(IC10),CC(IC11),CC(IC12),CC(IC13),
C.................ZETA     EPS      EPA      DT       XA       LBLI.....
     2            CC(IC14),CC(IC15),CC(IC16),CC(IC17),CC(IC18),IA(IA19),
C.................BUFI.....
     3            CC(IC19))
C
  210 CONTINUE
      GO TO 400
  399 WRITE(3,12) ICMAX,MAXCOR
  400 CONTINUE
C
Ceeb
      CALL RCLOSE(ITAP36,3)
      CALL RCLOSE(MFILE,3)
      if (idel34.eq.1) then
         call rclose(itap34,4)
      else
         call rclose(itap34,3)
      endif
Ceeb

Ceeb  write(6,*) (nsorb(i),i=1,10)

      CALL TSTOP(6)
C
      STOP
      END
