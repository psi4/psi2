      SUBROUTINE START(A,IA,MAXCOR)
      IMPLICIT INTEGER(A-Z)
      INTEGER IA(MAXCOR)
      REAL*8   A(1)
      integer dist
      integer prcntr
      COMMON /MAXDIM/ MXAT,MXTRN,MXPR,MXSH,MXSYM,MXCF,MXT,MXISO,MXPSYM
      COMMON /PSIZE/  MSU,KAORDS,MCONSU,MGCSU,MSTU,MRU,MCSU,MCTU,MCONU,
     1                MCU,MSFU,MGU,MSFRU,MNSFU,MPRU,MCXU,MCCU,MBLU,MRCRU
      COMMON /INFOA/  NAT,NUM,NQW,NX
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV
      COMMON /INOUT/  IOUT,ICHECK,ITAP30
      COMMON /NSHEL/  NSHELL,MSHELL
C***********************************************************************
C     ALLOCATE CORE FOR MOLE.                                          *
C***********************************************************************
      ZAN=1
      C=ZAN+MXAT
      EX=C+3*MXAT
      C1=EX+MXPR
      CSPD=C1+MRCRU
      NGNCN=WPADTI(CSPD+MRCRU*MXPR*5)
      KSTART=NGNCN+MXSH
      KATOM=KSTART+MXSH
      KTYPE=KATOM+MXSH
      KNG=KTYPE+MXSH
      KLOC=KNG+MXSH
      KMIN=KLOC+MXSH
      KMAX=KMIN+MXSH
      NPERSH=KMAX+MXSH
      NEWSH=NPERSH+5
      TITLE=NEWSH+MXSH*MXTRN
      T=IADTWP(TITLE+20)
      INVT=WPADTI(T+MXT)
      ISO=INVT+MXTRN
      ICT=ISO+MXSH*MXISO
      NS=ICT+MXAT*MXTRN
      KS=NS+MXAT
      LBLAT=KS+MXAT
      CSINP=IADTWP(LBLAT+2*MXAT)
      NEED=WPADTI(CSINP+MRCRU*MXPR*5)
C***********************************************************************
C     CHECK IF THERE IS ENOUGH CORE FOR PROGRAM MOLE.                  *
C     ZERO OUT THE AMOUNT OF CORE THAT IS NEEDED.                      *
C***********************************************************************
ctph  WRITE(IOUT,546) NEED,MAXCOR
  546 FORMAT(//' NEED=',I12,'   MAXCOR=',I12//)
      IF (NEED.GT.MAXCOR) then
        WRITE(IOUT,546) NEED,MAXCOR
        write(6,*) 'ERROR: NOT ENOUGH MEMORY FOR MOLE'
        call qabort
        endif
      NTOZER=NEED-WPADTI(ZAN)+1
      CALL IZERO(A(WPADTI(ZAN)),NTOZER)
C***********************************************************************
C     MOLE READS IN THE BASIS SET AND OPTIONS, AND GENERATES ALL OF    *
C     THE ATOMS FROM A SYMMETRY UNIQUE LIST.                           *
C***********************************************************************
      CALL MOLE( A(ZAN),A(C),A(EX),A(C1),A(CSPD),IA(NGNCN),IA(NPERSH),
     1  IA(KSTART),IA(KATOM),IA(KTYPE),IA(KNG),IA(KLOC),IA(KMIN),
     2  IA(KMAX),IA(NEWSH),IA(TITLE), A(T),IA(INVT),IA(ISO),IA(ICT),
     3  IA(NS),IA(KS),IA(LBLAT), A(CSINP), ia(need), maxcor-need)
C***********************************************************************
C     CALCULATE THE INTERNUCLEAR DISTANCES.                            *
C***********************************************************************
c Print out various distances if the user has requested it.
      if (prcntr('IS_ON BOHR').ne.0) then
        write(6,*)
        write(6,*) 'The interatomic distances in bohr:'
        call coinit(ia(need),maxcor-need)
        dist = coallc(ia(need),(nat*(nat+1))/2)
        CALL atomd(A(C),mxat,nat,ia(need-1+dist),1.D+00)
        endif
      if (prcntr('IS_ON ANGSTROM').ne.0) then
        write(6,*)
        write(6,*) 'The interatomic distances in angstroms:'
        call coinit(ia(need),maxcor-need)
        dist = coallc(ia(need),(nat*(nat+1))/2)
        CALL atomd(A(C),mxat,nat,ia(need-1+dist),0.529177249D+00)
        endif
      if ((prcntr('IS_ON RADIAN').ne.0).and.(nat.gt.2)) then
        write(6,*)
        write(6,*) 'The bond angles in radians:'
        CALL bonda(A(C),mxat,nat,1.D+00)
        endif
      if ((prcntr('IS_ON DEGREE').ne.0).and.(nat.gt.2)) then
        write(6,*)
        write(6,*) 'The bond angles in degrees:'
        CALL bonda(A(C),mxat,nat,57.29577951D+00)
        endif
C***********************************************************************
C     ALLOCATE CORE FOR SYMORB.                                        *
C***********************************************************************
      PTR=IADTWP(NEED)
      DTR=PTR+3*3*MXTRN
      FTR=WPADTI(DTR+6*6*MXTRN)
      GTR=WPADTI(FTR+10*10*MXTRN)
      ISC=WPADTI(GTR+15*15*MXTRN)
      IPC=ISC+MXSH*MXPSYM
      SB=IPC+MXSH*MXPSYM
      CC=IADTWP(SB+MXSYM)
      LAB=WPADTI(CC+MXCF)
      LOC1=LAB+MXCF
      LOC2=LOC1+MXAT*5
      KLASS=LOC2+MXAT*5
      LAMBDA=KLASS+MXAT
      ISGRP=LAMBDA+MXSYM
      IADDR=ISGRP+MXAT
      IPA=IADDR+MXSYM
      GEN1=IADTWP(IPA+MXSYM)
      GEN2=GEN1+MXSYM*3*3
      GEN3=GEN2+MXSYM*3*3
      GAM=GEN3+MXSYM*3*3
      TRACEP=GAM+MXTRN**2
      TRACED=TRACEP+MXTRN
      TRACEF=TRACED+MXTRN
      TRACEG=TRACEF+MXTRN
      CHARR=TRACEG+MXTRN
      NEED=WPADTI(CHARR+MXTRN*MXSYM)
C***********************************************************************
C     CHECK IF THERE IS ENOUGH CORE FOR PROGRAM SYMORB.                *
C     ZERO OUT THE AMOUNT OF CORE THAT IS NEEDED.                      *
C***********************************************************************
CTPH  WRITE(IOUT,546) NEED,MAXCOR
      IF (NEED.GT.MAXCOR) then
        WRITE(IOUT,546) NEED,MAXCOR
        write(6,*) 'ERROR: NOT ENOUGH MEMORY FOR SYMORB'
        call qabort
        endif
      NTOZER=NEED-WPADTI(PTR)+1
      CALL IZERO(IA(WPADTI(PTR)),NTOZER)
C***********************************************************************
C     SYMORB GENERATES THE SYMMETRY ADAPTED BASIS FUNCTIONS.           *
C***********************************************************************
      CALL SYMORB(A(T),IA(INVT),IA(ISO),A(PTR),A(DTR),IA(KATOM),
     1  IA(KTYPE),IA(ISC),IA(IPC),IA(SB),A(CC),IA(LAB),IA(LOC1),
     2  IA(LOC2),IA(KLASS),IA(LAMBDA),IA(IADDR),IA(IPA), A(GEN1),
     3  A(GEN2),A(GEN3),IA(KMIN),IA(ICT),A(GAM),A(TRACEP),A(TRACED),
     4  A(CHARR),IA(ISGRP),A(FTR),A(GTR),A(TRACEF),A(TRACEG))
C

      if (prcntr('IS_ON BRIEF').eq.0) then
        WRITE (IW,85) MSU,KAORDS,MCONSU,MGCSU,MSTU,MRU,MCSU,MCTU,
     1    MCONU,MCU,MSFU,MGU,MSFRU,MRCRU
   85   FORMAT(//'    MSU:',I5/' KAORDS:',I5/' MCONSU:',I5/
     1           '  MGCSU:',I5/'   MSTU:',I5/'    MRU:',I5/
     2           '   MCSU:',I5/'   MCTU:',I5/'  MCONU:',I5/
     3           '    MCU:',I5/'   MSFU:',I5/'    MGU:',I5/
     4           '  MSFRU:',I5/'  MRCRU:',I5///)
        endif
      WRITE(IW,8877) MSHELL,MSFRU,NAT
C***********************************************************************
C           WRITE OUT HONDO ARRAYS TO FILE 10 NOW 30                   *
C           SET UP ADDRESSES FOR REAL VARIABLES IN INTEGER ARRAYS      *
C           THE SEL-87 REQUIRES USE OF INTEGER ARRAYS FOR WORD-        *
C           ADDRESSABLE RANDOM ACCESS I/O                              *
C***********************************************************************
      RZAN=WPADTI(ZAN)
      RC=WPADTI(C)
      REX=WPADTI(EX)
      RCSPD=WPADTI(CSPD)
      RCC=WPADTI(CC)
      RT=WPADTI(T)
      RPTR=WPADTI(PTR)
      RDTR=WPADTI(DTR)
      RFTR=WPADTI(FTR)
      RGTR=WPADTI(GTR)
      MXPRD=INTOWP(MXPR)
C***********************************************************************
C    SET MXPR TO THE TRUE NUMBER OF PRIMITIVES, MXPRO IS ORIGINAL MXPR *
C***********************************************************************
      MXPRO=MXPR
      MXPR=IA(KSTART+NSHELL-1)+IA(KNG+NSHELL-1)-1
C***********************************************************************
C     FIND THE LARGEST NUMBER OF GENERAL CONTRACTIONS.                 *
C***********************************************************************
      IRCRU=0
      DO 110 I=1,NSHELL
         IRCRU=MAX(IRCRU,IA(NGNCN+I-1))
  110 CONTINUE
C***********************************************************************
C  CSPD(K,J,I) WITH DIMENSION CSPD(IRCRU,MXPR,5)=                      *
C     CSPD(K,J,I) WITH DIMENSION CSPD(MRCRU,MXPRO,5)                   *
C     WHERE MRCRU AND MXPRO ARE THE MAXIMUM POSSIBLE VALUES AND        *
C     IRCRU AND MXPR ARE THE ACTUAL MAXIMUM VALUES FOR THE BASIS SET.  *
C***********************************************************************
      DO 120 I=1,5
      DO 120 J=1,MXPR
      DO 120 K=1,IRCRU
        A(CSPD+(IRCRU*MXPR)*(I-1)+IRCRU*(J-1)+K-1)=
     1  A(CSPD+(MRCRU*MXPRO)*(I-1)+MRCRU*(J-1)+K-1)
  120 CONTINUE
      MRCRU=IRCRU
C***********************************************************************
C     WRITE OUT ARRAYS TO FILE 30.                                     *
C***********************************************************************
      CALL OUT(IA(RZAN),IA(RC),IA(ICT),IA(NS),IA(KS),IA(REX),IA(RCSPD),
     1  IA(NGNCN),IA(KSTART),IA(KATOM),IA(KTYPE),IA(KNG),IA(KLOC),
     2  IA(KMIN),IA(KMAX),IA(NPERSH),IA(ISC),IA(IPC),IA(RCC),IA(LAB),
     3  IA(LOC1),IA(LOC2),IA(KLASS),IA(LAMBDA),IA(IADDR),IA(IPA),IA(RT),
     4  IA(INVT),IA(ISO),IA(TITLE),IA(LBLAT),IA(SB),IA(RPTR),
     5  IA(RDTR),MXPRD,IA(RFTR),IA(RGTR))
C
      RETURN
 8877 FORMAT(/,' TOTAL NUMBER OF SHELLS ',13X,'=',I5,/,
     1  ' TOTAL NUMBER OF BASIS FUNCTIONS     =',I5,/,
     1  ' TOTAL NUMBER OF ATOMS',15X,'=',I5)
      END
