      SUBROUTINE PITZER(RA,IA,MAXCOR)
C
C
C  SCF PROGRAM WRITTEN BY R.M.PITZER ET AL, MODIFIED BY PAUL SAXE,
C  GEORGE FITZGERALD, MARK VINCENT AND YUKIO YAMAGUCHI AT BERKELEY.
C
C
C  1 JULY 1983  MODIFIED INPUT SECTIONS TO SIMPLIFIY INPUT AND USE
C               INFORMATION ON HEADER RECORD OF INTEGRALS.
C                                                           PWS
C
C
C--------------------------------------------------------------------
C
      IMPLICIT INTEGER (A-Z)
C
C
cets022791
#include <error.h>
      integer frdc,frdi,errcod
c
      REAL*8 RA(1),SECS
      CHARACTER*6 CALTYP,DERTYP,CITYP
      INTEGER IA(MAXCOR)
      CHARACTER*1 P,PP
      INTEGER TIME
      LOGICAL IPOPLE
      REAL*8 REPNUC,RTIME,RTIMD
      REAL*8 FMT,TEIMEN
      REAL*8 Z
      INTEGER BLABEL,ALABEL,PRINT
      INTEGER NUCREP(1),FCB44(16),FCB10(16)
C
      COMMON /T10/ I10(200)
      COMMON /DERT/ CALTYP,DERTYP,CITYP
      COMMON /PRINTQ/ PRINT
      COMMON /KIND/ NKIND,ICLOS
      COMMON /STIMEX/ KTIMES,ITIME5,ITIME6,ITIME7,ITIME8,ICONT5,ICONT6
     1,              ICONT7,ICONT8,ITIMED,ICONTD,ITEMPT
      COMMON /SIGNS/IOFF(200),MASK,IPRNT
      COMMON /SCTRS/ISECT,MXISEC
      COMMON /SYMMS/ND(10),LLSYM(10)
      COMMON /POP/ IPOPLE
      COMMON /RROB/ BTAG,POINTR,TAGMAX,RDCNT
      COMMON /INPT/  NSYM,IPMAX,NMAT,ITMAX,NSFMAX,NTOTAL,MXCOEF,MSFRU
      COMMON /SCFIN/ REPNUC,IOPEN,ISCF,IPRCT,IVECT,NBFLD,ICUT
      COMMON /TAPES/ NIN,IR,IOUT,ICHECK,IFLAG3,ITAPE,I30
      COMMON /T102/ NOUT,HALF,THIRD,MAXNO,OUTBUF,KINADR
      COMMON /HOLL/ BLABEL(20),ALABEL(26),FMT(10),JL
      COMMON /TIMX/ ITIME,ICONT0,JCOUNT,ITIME2,ICONT2,ITIME3,ICONT3
     1,             ITIME4,ICONT4,ITIME1,ICONT1
      COMMON /PASS1/IDUMC(7)
      COMMON /PASS2/IDUMD(11)
      COMMON /LIMIT/TEIMEN,MAXBUF
cets022791
      logical twocon
      COMMON /GVB/ ALPHA1,ALPHA2,ALPHA3,IJUP,TWOCON
C     COMMON Z(1)
C
C
C
      EQUIVALENCE (NUCREP,REPNUC)
C     EQUIVALENCE (Z,A)
C
C
      DIMENSION ITYP(20),NSO(20),MTYPE(2,20),MS(142),MNL(142)
CTJL  MASK=X'7FFFFF'
      ISECT=0
      TEIMEN=1.0D-10
      MAXBUF=3000
C
      IPRNT=0
      ITAPE=34
      IPOPLE=.FALSE.
      SIZE=8942
C1004 PSIZE=20736   PSIZE = N(DSHELLS)*6**4   THIS FIX IS FOR MORE D'S
cets050391      PSIZE=41472
      psize=82944
      PSIZE2=PSIZE/2
      MAXBUF=PSIZE
C     MAXCOR=1
      INBUF=1532
      OUTBUF=22394
      OUTBUF=(OUTBUF-2)/12*12+2
      NIN=14
      I30 = 30
      IR=5
      ICHECK = 3
C     CALL LFN(3,IFLAG3)
C     IF((IFLAG3.AND.'2000000).EQ.0) THEN
      IFLAG3=3
C     ELSE
C     IFLAG3=-1
C     END IF
C     CALL LFN(4,IVAL)
C     IF(IVAL.EQ.0) THEN
      IOUT=6
C     ELSE
C     IOUT=4
C     END IF
      NOUT=92
C
C     ----- OPEN OUTPUT FILES FOR APPEND ACCESS -----
C
CTJL  OPEN (UNIT=IOUT,OPENMODE='A')
CTJL  OPEN (UNIT=IFLAG3,OPENMODE='A')
C
C
C
      IOFF(1)=0
      DO 1 I=1,199
        IOFF(I+1)=IOFF(I)+I
    1 CONTINUE
    2 CONTINUE
C
C     ----- READ HEADER INFORMATION FROM INTEGRAL TAPE -----
C
      CALL RFILE(ITAPE)
CTJL  CALL SREW(ITAPE)     THIS IS DONE IN RFILE
CSEL  READ(ITAPE)NKIND,ICLOS,BLABEL,REPNUC,NST,(ND(I),I=1,NST),(ITYP(I),
CSEL 1I=1,NST), (NSO(I), I=1,NST), NS, (MTYPE(IS), IS=1,NS),LAPU,(MS(I)
CSEL 2, I=1,LAPU), (MNL(I), I=1,LAPU)
CSEL  CALL RSETSIA(ITAPE,10)
C
C
      CALL WREADW(ITAPE,NKIND,1,1,IEND)
      CALL WREADW(ITAPE,ICLOS,1,IEND,IEND)
      CALL WREADW(ITAPE,BLABEL,20,IEND,IEND)
      CALL WREADW(ITAPE,NUCREP(1),INTOWP(1),IEND,IEND)
      CALL WREADW(ITAPE,NST,1,IEND,IEND)
      CALL WREADW(ITAPE,ND,NST,IEND,IEND)
      CALL WREADW(ITAPE,ITYP,NST,IEND,IEND)
      CALL WREADW(ITAPE,NSO,NST,IEND,IEND)
      CALL WREADW(ITAPE,NS,1,IEND,IEND)
      CALL WREADW(ITAPE,MTYPE,2*NS,IEND,IEND)
      CALL WREADW(ITAPE,LAPU,1,IEND,IEND)
      CALL WREADW(ITAPE,MS,LAPU,IEND,IEND)
      CALL WREADW(ITAPE,MNL,LAPU,IEND,IEND)
C   THIS IS TO PLACE THE CURRENT FILE POINTER ON ITAPE TO A SECTOR
C   BOUNDARY.
      ITOT = 25 + INTOWP(1) + 2*NS + 3*NST + 2*LAPU
      ISADR = I2SEC(ITOT) + 1
      CALL RSETSA(ITAPE,ISADR)
C
C
C
CDBG  WRITE (3,523) NUCREP(1),NUCREP(2),NUCREP(3),REPNUC
C 523 FORMAT (' NUCREP:',3Z12,Z20)
CDB   WRITE (3,546)
C 546 FORMAT (//,' BACK FROM READING LABEL OF INTEGRALS ',//)
CDB   WRITE (3,548) NKIND,NCLOS,REPNUC,NST,NS,LAPU,(BLABEL(IQ),IQ=1,20)
C 548 FORMAT (/,' NKIND=',I5,' NCLOS=',I5,' REPNUC=',F10.4,' NST='
CDB  #,       I5,/,' NS=',I5,' LAPU=',I5,/,(1X,20A4))
CDB   WRITE (3,549) (ND(IQ),IQ=1,NST)
C 549 FORMAT (' ND:',(T10,10I5))
CDV   WRITE (3,550) (ITYP(IQ),IQ=1,NST)
C 550 FORMAT (' ITYP:',(T10,10A5))
CDV   WRITE (3,551) (NSO(IQ),IQ=1,NST)
C 551 FORMAT (' NSO:',(T10,10I5))
CD    WRITE (3,552) (MTYPE(1,IQ),MTYPE(2,IQ),IQ=1,NS)
C 552 FORMAT (' MTYPE:',(T10,5(2A4,2X)))
CD    WRITE (3,553) (MS(IQ),IQ=1,LAPU)
C 553 FORMAT (' MS:',(T10,10I5))
CD    WRITE (3,554) (MNL(IQ),IQ=1,LAPU)
C 554 FORMAT (' MNL:',(T10,10I5))
      IF (NKIND.NE.1 .AND. NKIND.NE.2) THEN
        WRITE (IOUT,900)
  900   FORMAT (14X,' ***** PROBLEMS WITH INTEGRAL FILE *****',/
     #  ,            14X,' ***** USE INTD2H OR INTSYM        *****')
        IF (IFLAG3.GT.0) WRITE (IFLAG3,900)
        STOP
      END IF
      ISECT=ISECT+1
      NSYM=0
      MXCOEF=0
      NSFMAX=0
      NMAT=0
      MSFRU=0
      DO 3 I=1,NST
        N=NSO(I)
        IF (N.LE.0) GO TO 3
        NSYM=NSYM+1
        IF (N.GT.NSFMAX) NSFMAX=N
        NMAT=NMAT+N*(N+1)/2
        MXCOEF=MXCOEF+N**2
        MSFRU=MSFRU+N
    3 CONTINUE
      NSYM2=NSYM
      MXCOF2=MXCOEF
      LLSYM(1)=0
      IF (NSYM .NE. 1)THEN
        II=0
        DO 4 I=2,NST
          IF (NSO(I).LE.0) THEN
            LLSYM(I)=LLSYM(I-1)
          ELSE
  201       II=II+1
            NN=NSO(II)
            IF (NN.EQ.0) GOTO 201
            LLSYM(I)=LLSYM(I-1)+IOFF(NN+1)
          END IF
    4   CONTINUE
      END IF
C
C
      IF (NSYM.GT.8 .AND. NKIND.EQ.1) THEN
        WRITE (IOUT,901) NSYM
  901   FORMAT (' ***** SCF ONLY RUNS D2H AND SUBGROUPS ****',I10)
        IF (IFLAG3.GT.0) WRITE (IFLAG3,901) NSYM
        STOP
      END IF
C
C     ----- ALLOCATE SOME CORE -----
C
      MAX=MAXCOR/INTOWP(1)
      CALL ZERO(RA,MAX)
C
      NORBS=1
      NC=NORBS+NSYM
      NO=NC+NST
      NLAMDA=NO+NST
      IJL=NLAMDA+NSYM
      IJU=IJL+NSYM
      ALPHA=IADTWP(IJU+NSYM)
      BETA=ALPHA+NSYM*(NSYM+1)/2
      SU=BETA+NSYM*(NSYM+1)/2
      SU2=SU+MXCOEF
      CM=SU2+MXCOEF
C
C     -----INTEGER EQUIVALENTS -----
C
      RALPHA=WPADTI(ALPHA)
      RBETA=WPADTI(BETA)
      RSU=WPADTI(SU)
      RCM=WPADTI(CM)
C
C
C
C
      CALL INPUT(IA(NC),IA(NO),RA(ALPHA),RA(BETA),IA(RSU),ITYP,NST,
     1RA(SU),RA(SU2),NSO)
C
C
C
      CALL ADDTO(ITYP,IA(NORBS),IA(NC),IA(NO),IA(NLAMDA),IA(IJL),IA(IJU)
     1,          RA(ALPHA),RA(BETA),RA(SU),RA(CM),NSO,NST)
C
C
C
C
      IF (IAND(PRINT,4096).NE.0) IPRNT=1
C
      D1OLD=CM+MSFRU*2
      DT=D1OLD+NMAT
      DO=DT+NMAT
      S=DO+NMAT
      T=S+NMAT
      V=T+NMAT
      C=V+NMAT
      H=C+NMAT
      UNIT=H+NMAT
      C1=UNIT+NMAT
      C2=C1+NSYM
      SD=C2+NSYM
      NV=WPADTI(SD+NSFMAX*(NSFMAX+1)/2)
      NX=NV+NSYM*2
      IDMAT=IADTWP(NX+NSFMAX)
C
C  PARAMETERS FOR SOCONV ROUTINES
C
      PA=IDMAT
      IF(IPOPLE) PA=IDMAT+8*NMAT+2
      VAL1=PA
      PB=PA+PSIZE2
      VAL2=PB
      PC=PB+PSIZE2
      LBLINT=PC+PSIZE2
      VAL0=LBLINT
      INEXT=WPADTI(LBLINT+PSIZE2)
      LBLI=IADTWP(INEXT+PSIZE2+128)
      OBUF=LBLI+SIZE/INTOWP(1)
C
C     ----- INTEGER EQUIVALENTS -----
C
      ROBUF=WPADTI(OBUF)
      RLBLIN=WPADTI(LBLINT)
      RLBLI =WPADTI(LBLI)
C
      NOPN=5
      IF (IOPEN.EQ.0)NOPN=3
      IF (WPADTI(OBUF)+NOPN*OUTBUF.GT.MAXCOR) THEN
        WRITE(IOUT,100) OBUF+OUTBUF*NOPN
        IF(IFLAG3.GT.0) WRITE(IFLAG3,100) OBUF+OUTBUF*NOPN
  100   FORMAT(' ******BLANK COMMON NOT LARGE ENOUGH, INCREASE TO',I7
     *  )
        STOP
      END IF
C
C
C
      CALL INONE(RA(S),RA(T),RA(V),RA(C),RA(H),NMAT,INBUF,IA(RLBLI)
     #,          RA(LBLI))
C
C
C
CTEST CALL SCF1(RA(CM),IA(NLAMDA),IA(NV),RA(SU),ITYP,RA(C1),RA(C2)
CTEST1,  IA(NX),IA(NORBS),RA(ALPHA),RA(BETA),IA(NC),IA(NO),RA(S)
CTEST2,      RA(H),RA(T),RA(V),RA(DT),RA(DO),RA(SD),RA(IDMAT),RA(D1OLD))
C
C
C
      NTOTAL=MSFRU
      NDIIS=8
      NFOCK=1
      IF (IOPEN.NE.0) NFOCK=2
C
      call wreadw(30,i10,200,101,junk)
      nxx=i10(24)
c
      SF1=PA
      SF11=SF1+NMAT
      SF12=SF11+NMAT
      SF21=SF12+NMAT
      SF22=SF21+NMAT
      FC=SF22+NMAT
      FO=FC+NSFMAX*(NSFMAX+1)/2
      EIG=FO+NSFMAX*(NSFMAX+1)/2
      SF2=EIG+NTOTAL
      FSAVE=SF2+NMAT
      ERROR=FSAVE+NMAT*NFOCK*NDIIS
      SMHALF=ERROR+NMAT*NDIIS
      U=SMHALF+NMAT
      EIGVAL=U+MXCOEF
      T1=EIGVAL+NTOTAL
      T2=T1+NSFMAX**2
      TRIANG=T2+NSFMAX**2
      BMAT=TRIANG+NMAT
      wmsiz=nsfmax*(nsfmax+1)/2
      wmot=bmat+(ndiis+1)*(ndiis+2)
      sd2=wmot+nxx
      buffer=sd2+msfru
      RBUFFR=WPADTI(BUFFER)
      REIG=WPADTI(EIG)
      TAGMAX=(MAXCOR-RBUFFR-1)/OUTBUF
      IF (TAGMAX .LT. 2) THEN
        WRITE (IOUT,902)
  902   FORMAT (' ***** CANNOT MAKE TWO BUFFERS FOR CONTR2 *****')
        IF (IFLAG3.GT.0) WRITE (IFLAG3,902)
        STOP
      END IF
C
C
C
      CALL INTWO (RA(ALPHA),RA(BETA),IA(IJL),IA(IJU),RA(VAL1),RA(VAL2)
     1,          IA(ROBUF),RA(OBUF),INBUF2,RA(PA),RA(PB),RA(PC)
     2,           IA(RLBLIN),RA(LBLINT))
C
C
      IF (MAXNO.GT.0) GO TO 298
C
C
      CALL INFIND (RA(PA),RA(PB),RA(PC),IA(RLBLIN),IA(INEXT),PSIZE2
     #,            PSIZE)
C
C
C
      CALL TWOI(RA(PA),RA(PB),RA(PC),IA(RLBLIN),IA(RLBLI),RA(LBLI)
     1,         RA(VAL1),RA(VAL2),ra(alpha),ra(beta),ia(ijl),ia(iju)
     #,         ia(robuf),ra(obuf),inbuf2,ia(inext),psize2,psize)
C
C
C
      CALL CLNUP(ra(alpha),ra(beta),ia(ijl),ia(iju),ra(val1),ra(val2),
     +           ia(robuf),ra(obuf),inbuf2,ra(pa),ra(pb),ra(pc),
     +           ia(rlblin),ra(lblint))
C
C
C
      CALL TAKON (RA(T),RA(V),RA(C),NMAT)
C
C
C
  298 CONTINUE
      WRITE(ICHECK,299)TAGMAX,OUTBUF,BUFFER,MAXNO
CRBR  IF(IFLAG3.GT.0)WRITE(ICHECK,299)TAGMAX,OUTBUF,BUFFER,MAXNO
299   FORMAT(/'  TOTAL OF ',I4,' BUFFERS AVAILABLE,'/'  EACH OF SIZE ',
     1I7,'  STARTING AT ADDRESS ',I7/,'  TOTAL OF ',I4,' BUFFERS ARE
     2 REQUIRED '/)
C
C      ----- REZERO CORE USED TO READ IN INTEGRALS -----
C
      MAX=MAXCOR/INTOWP(1)-SF1
      CALL ZERO(RA(SF1),MAX)
C
      CALL PRIME(IA(RBUFFR))
C
C
C
      CALL SCF1(RA(CM),IA(NLAMDA),IA(NV),RA(SU),ITYP,RA(C1),RA(C2)
     1,  IA(NX),IA(NORBS),RA(ALPHA),RA(BETA),IA(NC),IA(NO),RA(S)
     2,      RA(H),RA(T),RA(V),RA(DT),RA(DO),RA(SD),RA(IDMAT),RA(D1OLD)
     3,         RA(EIG),RA(SF1),RA(SF2),RA(SF11),RA(SF12),RA(SF21)
     1,         RA(SF22),RA(FC),RA(FO),IA(RBUFFR),RA(BUFFER)
     2,         IA(RCM),IA(RSU),IA(RALPHA),IA(RBETA),IA(REIG)
     3,         RA(UNIT)
     5,  RA(FSAVE),RA(ERROR),NFOCK,NDIIS,RA(SMHALF),RA(U),RA(EIGVAL)
     6,  RA(T1),RA(T2),RA(TRIANG),RA(BMAT),ra(wmot),ra(sd2),wmsiz)
C
CTJL  CALL CPUTIM(SECS)
CTJL  MINS=SECS/60
CTJL  HRS=MINS/60
CTJL  SECS=SECS-MINS*60
CTJL  MINS=MINS-HRS*60
CTJL  WRITE (IOUT,312) HRS,MINS,SECS
C 312 FORMAT (/,T15,'SCF COMPLETED IN',I4,':',I2,':',F6.3,///)
C
C  TFOCK ROUTINES FOR SCF FIRST DERIVATIVES
C
      CALL LOCATE(5,'# TFOCK ##',IERROR)
cets022791
      IF (IERROR .NE. 0) THEN
        dertyp = 'FIRST '
        cityp = 'SCF   '

        errcod = frdi('IPRINT2',iprint)
        if(errcod.ne.EOK) iprint = 0

        errcod = frdc('DERTYPE',dertyp)

        errcod = frdc('WFN',cityp)

        if(twocon) then
          caltyp = 'TCSCF '
        else if(cityp.ne.'SCF   '.or.iopen.eq.1) then
          caltyp = 'GRSCF '
        else
          caltyp = 'CLSCF '
        endif
      else
        read(5,922) caltyp,cityp,dertyp
 922    FORMAT(3(A6,4X))
        READ(5,923) IPRINT
 923    FORMAT(I5)
      END IF
c
      IF (CITYP.EQ.'      ') CITYP='SCF   '
      IF (DERTYP.EQ.'      ') DERTYP='FIRST '
      WRITE(IOUT,924) CALTYP,CITYP,DERTYP
      WRITE(ICHECK,924) CALTYP,CITYP,DERTYP
 924  FORMAT(3X,3(A8,2X))
C     IF (CITYP .NE. 'SCF   ') THEN
C       WRITE(IOUT,921) CITYP
C921    FORMAT(/'  ONLY ABLE TO PERFORM SCF DERIVATIVES.'/
C    1  /'  YOUR CITYP IS ',A6)
C       STOP
C     END IF
cets030291 no need for tfock if not a ci gradient
      if(cityp.ne.'CI    ') return
c
      ITAP30=30
CTJL  CALL RFILE(ITAP30)   THIS FILE ALREADY OPENED
      CALL WREADW(ITAP30,I10,200,101,JUNK)
      NBFSO=I10(18)
      NBFAO=I10(22)
cets030291 for now, can't do 5d gradients
      if(nbfso.ne.nbfao) then
        write(6,*) ' '
        write(6,*) '************ WARNING ******************'
        write(6,*) '  cannot do 5d ci gradients just yet'
        write(6,*) ' '
        return
      endif
c
      MPOINT=I10(2)
      MCONST=I10(3)
      NATOM=I10(19)
      JUNK=101+MCONST
      CALL WREADW(ITAP30,I10,MPOINT,JUNK,JUNK)
      LOCU=I10(29)
      EIG=D1OLD
      SOMOV=EIG+NTOTAL
      SOAOV=SOMOV+NBFSO*NBFSO
      AOMOV=SOAOV+NBFSO*NBFAO
      BUF=AOMOV+NBFSO*NBFSO
      WMO=BUF+NBFSO*NBFSO
      WMO2=WMO+NSFMAX*(NSFMAX+1)/2
      WMOT=WMO2+NSFMAX*(NSFMAX+1)/2
      NBATRI=NBFAO*(NBFAO+1)/2
      WAO=WMOT+NBATRI
      DC=WAO+NBATRI
      DO=DC+NBATRI
      DO2=DO+NBATRI
      SCRTCH=DO2+NBFSO*NBFSO
      ALPHA=SCRTCH+NBFAO*NBFAO
      BETA=ALPHA+NBATRI
      ITOP=BETA+NBATRI
      IF (ITOP.GT.MAXCOR) THEN
        WRITE(IOUT,905)
 905    FORMAT(/'  ITOP=',I6,'  MAXCOR=',I6)
        STOP
      END IF
      ISOAOV=WPADTI(SOAOV)
C
      DO 906 III=SOMOV,ITOP
 906  RA(III)=0.0D+00
C
      WMSIZ = NSFMAX*(NSFMAX + 1)/2
CRBRBRBRBRBRBRBR
CC    WRITE(ICHECK,*) ' WMSIZ = ',WMSIZ
CC    WRITE(ICHECK,*) 'CM,SU,SOMOV,SOAOV,AOMOV,
CC   1  WMO,WMO2,WMOT,WAO,DC,DO,DO2,
CC   2  SCRTCH,NBFSO,NBFAO,LOCU,NBATRI,ISOAOV,NLAMDA,
CC   3  ITAP30,NSYM2,MXCOF2,IPRINT,EIG,NTOTAL,NATOM,ALPHA,
CC   4  BETA,BUF,WMSIZ ',
CC   3                 CM,SU,SOMOV,SOAOV,AOMOV,
CC   1  WMO,WMO2,WMOT,WAO,DC,DO,DO2,
CC   2  SCRTCH,NBFSO,NBFAO,LOCU,NBATRI,ISOAOV,NLAMDA,
CC   3  ITAP30,NSYM2,MXCOF2,IPRINT,EIG,NTOTAL,NATOM,ALPHA,
CC   4  BETA,BUF,WMSIZ
C
      CALL TFOCK(RA(CM),RA(SU),RA(SOMOV),RA(SOAOV),RA(AOMOV),
     1RA(WMO),RA(WMO2),RA(WMOT),RA(WAO),RA(DC),RA(DO),RA(DO2),
     2RA(SCRTCH),NBFSO,NBFAO,LOCU,NBATRI,IA(ISOAOV),IA(NLAMDA),
     3ITAP30,NSYM2,MXCOF2,IPRINT,RA(EIG),NTOTAL,NATOM,RA(ALPHA),
     4RA(BETA),RA(BUF),WMSIZ)
      RETURN
      END
