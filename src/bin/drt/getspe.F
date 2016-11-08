      SUBROUTINE GETSPE(SPC,oldinp)
C
C***********************************************************************
C     READ IN A CODE TELLING WHAT TYPE OF SPECIAL ORBITALS TO HANDLE,  *
C     OR READ IN THE MATRIX ITSELF UNDER THE 'MAT' OPTION.             *
C***********************************************************************
C
      IMPLICIT INTEGER (A-Z)
#include <error.h>
      INTEGER NUMINT
      integer errcod, frdc
      logical oldinp
      character*8 optyp
      CHARACTER*1 MULTRF,VALENC
      CHARACTER*3 GVBW,OSSW,SPCODE,CODES,WORDS
C
      COMMON /TAPES/  OUT,ERROUT,INPUT,DRTTAP
      COMMON /DIMENS/ NBF,NSYM,NORBS,NROWSP,NROW4P,NROWS,NROWS4
     #,               NLEVS,NREFS,NROWOC,NROW4O,NWKS,NWKSOC,NLEVOC
     #,               ORBFRM,SYMORB,NUMIJ,NGROUP,NUMINT,NMAX,NSPC,NVREF
     #,               NIJVIR
      COMMON /DRTINF/ NA,NB,NS,NESPEC,MAXB,LEVFRM,LEVVAL,LEVOPN,LEVMUL
     #,               LEVOCC,SPEC,SSPESH,VAL
      COMMON /DRTCOD/ NCODES,DELA(9),DELB(9),DELELE(9)
     #,               NTYPES,VIRTUL,OCCUPD,VALOCC,RESCOR,RESVIR,FROZEN
     #,               VALVIR,OPENSH,MULTI,SPESHL,MULTRF,VALENC
     #,               WORDS(6,9),CODES(9)
C
      DIMENSION SPC(NSPC),GVB(16),GVBW(6),OSS(16),OSSW(6)
C
      DATA GVB/ 0,0,0,0, 0,1,1,1, 0,1,1,1, 0,1,1,2 /
      DATA OSS/ 0,0,0,1, 0,1,0,1, 0,0,1,1, 1,1,1,2 /
      DATA GVBW/'GVB',' OR','BIT','AL ','   ','   '/
      DATA OSSW/'OPE','N S','HEL','L S','ING','LET'/
C
      if(oldinp) then
        READ (INPUT,1,END=99) SPCODE
      else
        optyp = 'NONE    '
        errcod = frdc('OPENTYPE',optyp)

        spcode=' '
        if(optyp.eq.'SINGLET ') spcode='OSS'
        if(optyp.eq.'TWOCON  ') spcode='GVB'
      endif   

    1 FORMAT (A3)
      IF (SPCODE.NE.'GVB') GO TO 6
      IF (SPEC.EQ.2) GO TO 3
      WRITE (OUT,2) SPEC
      WRITE (ERROUT,2) SPEC
    2 FORMAT (//,' CALLING FOR GVB REFERENCE WITH',I3,' SPECIAL '
     #,          'ORBITALS: IMPOSSIBLE',//)
      CALL mabort
    3 CONTINUE
      DO 4 I=1,16
         SPC(I)=GVB(I)
    4 CONTINUE
      DO 5 I=1,6
         WORDS(I,9)=GVBW(I)
    5 CONTINUE
      NESPEC=2
      NA=NA+1
      WRITE (OUT,16)
   16 FORMAT (//,T30,'** GVB TWO-REFERENCE CI **')
      RETURN
C
    6 CONTINUE
      IF (SPCODE.NE.'OSS') GO TO 11
      IF (SPEC.EQ.2) GO TO 8
      WRITE (OUT,7) SPEC
      WRITE (ERROUT,7) SPEC
    7 FORMAT (//,' CALLING FOR AN OPEN-SHELL SINGLET WITH',I3,' SPECIAL'
     #,          ' ORBITALS: IMPOSSIBLE',//)
      CALL mabort
    8 CONTINUE
      DO 9 I=1,16
         SPC(I)=OSS(I)
    9 CONTINUE
      DO 10 I=1,6
         WORDS(I,9)=OSSW(I)
   10 CONTINUE
      NESPEC=2
      NA=NA+1
      NS=IEOR(NS,SSPESH)
      WRITE (OUT,17)
   17 FORMAT (//,T30,'** OPEN SHELL SINGLET REFERENCE **')
      RETURN
C
   11 CONTINUE
      IF (SPCODE.NE.'MAT') GO TO 90
      READ (INPUT,12,END=99) NESPEC,B,S
   12 FORMAT (3I5)
      READ (INPUT,13,END=99) (SPC(I),I=1,NSPC)
   13 FORMAT (4(4I1,1X))
      WRITE (OUT,14) NESPEC,B,S,(SPC(I),I=1,NSPC)
   14 FORMAT (//,T25,'** USING INPUT MATRIX FOR SPECIAL ORBITALS **'
     #,       //,T30,'NUMBER OF ELECTRONS:',I2
     #,        /,T30,'     SPIN TIMES TWO:',I2
     #,        /,T30,'           SYMMETRY:',I2
     #,       //,(T30,4(4I1,1X)))
      NA=NA+(NESPEC-B)/2
      NB=NB+B
      NS=IEOR(NS,S-1)
      RETURN
C
   90 CONTINUE
      WRITE (ERROUT,91) SPCODE
   91 FORMAT (//,' NO CODE AVAILABLE YET FOR SPECIAL CODE OF ',A3,//)
      CALL mabort
   99 CONTINUE
      WRITE (OUT,100)
      WRITE (ERROUT,100)
  100 FORMAT (//,' COULD NOT READ IN ALL INFORMATION FOR',I3,' SPECIAL '
     #,          'ORBITALS',/,' PLEASE CHECK YOUR INPUT DECK.',//)
      CALL mabort
C
      END
