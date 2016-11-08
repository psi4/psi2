      subroutine fentry(ic,cc,maxcor)
C
C >>>  LAST UPDATE september 19,1989 <<<
C
C***********************************************************************
C  NOTICE OF PROGRAM MODIFICATIONS                                     *
C      BY:  GUSTAVO E. SCUSERIA                                        *
C                                                                      *
c    date:  september 19, 1989
c  reason:  avoid transforming frozen mo's
c
C    DATE:  SEPTEMBER 16, 1988                                         *
C  REASON:  FORCE NBFAO=NBF FOR NEW GCINTS CODE.                       *
C                                                                      *
C    DATE:  APRIL 13, 1989                                             *
C  REASON:  CHANGE READ30 ROUTINE TO NEW VERSION (NEW SYMMETRY ORDER)  *
C                                                                      *
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 CC(1)
      INTEGER IC(maxcor*2),I30(200)
      INTEGER WPADTI,SEC2I,SPACE,SPACE2,nfzo(8),nfzv(8),ndoc(8),nuoc(8)
      character*80 label
      COMMON/BASIS/NBASIS,MTRI,NSYM
      call drum
      CALL TSTART(6)
C
      WRITE(6,6000)
 6000 FORMAT(//,63('-'),/,63('-'),/,
     .' CCTRANSL2 : integral transformation for large CC calculations'/,
     .'   written by Andy Scheiner & Gustavo Scuseria ',/,63('-'),/,
     .   63('-'),//)
C
    1 FORMAT ('                                                ')
    2 FORMAT ('************************************************')
  121 FORMAT (' NOT ENOUGH CC TO DO FINAL 2-ELECTRON HALF TRANSFORM')
C
      INPUT=5
      JOUT=6
      KOUT=8
      ITAP30=30
      ITAP71=71
      ITAP72=72
      ITAP77=77
      ITAP78=78
      ITAP34=34
      IBUFHF=4471
      IFLAG=1
      ISORT=1
      JBUF=0
      INTBUF=SEC2I(100)/INTOWP(1)
c     MAXBKT=120
c
c    change maximum number of buckets to 200
c
      MAXBKT=200
C
      CALL RFILE(ITAP71)
      CALL RFILE(ITAP72)
      CALL RFILE(ITAP77)
      CALL RFILE(ITAP78)
      CALL RFILE(ITAP30)
C
C >>> READ POINTERS IN FILE30
C
      CALL WREADW(ITAP30,I30,200,101,LWORD)
      MPOINT = I30( 2)
      MCONST = I30( 3)
      MCALCS = I30( 4)
      NCALCS = I30( 5)
      NBF    = I30(18)
      NAT    = I30(19)
      NBFAO  = I30(22)
      NSHELL = I30(27)
      NIRRED = I30(28)
      NPRIM  = I30(32)
      NSYMHF = I30(41)
      MXCOEF = I30(42)
      IOPEN  = I30(43)
C
C >>> FORCE NBFAO = NBF FOR NEW INTEGRALS CODE
C
      NBFAO=NBF
      NTRI=(NBF*(NBF+1))/2
      NTRIAO=(NBFAO*(NBFAO+1))/2
      MTRI=NTRIAO
C
C >>> ZERO OUT THE CC
C
      CALL ZERO(CC,MAXCOR)
C
C >>> DYNAMICAL ALLOCATION
C
      K1     = 1
      K2     = K1 + NBF*NBFAO
      K3     = WPADTI(K2 + NBF)
      KTPO   = K3 + NSYMHF
      KTYP   = KTPO + NBF
      KFLOV  = KTYP + NSYMHF
      KORBSM = KFLOV + 2*NIRRED
      KNMSM  = KORBSM + NBF
      KIWRK  = KNMSM + NIRRED
      KKLPT  = KIWRK + NIRRED
      KIJPT  = KKLPT + NTRI
      KOFF1  = KIJPT + NTRI
      KOFF2  = KOFF1 + NTRI
      K4     = KOFF2 + NTRI
      K5     = K4 + MCALCS
      fzt    = k5 + nsymhf
      itop   = fzt + nbf
      K6     = IADTWP( itop )
      K7     = K6 + MXCOEF
      KTOP   = K7 + NBF*NBF 
C
C >>> READ INPUT FILE
C
      call rdinpt(label,nfzo,nfzv,nirred)

C
      CALL READ30 (ITAP30,MPOINT,MCONST,MCALCS,NCALCS,NBF,NBFAO,NSYMHF,
     .             MXCOEF,CC(K1),CC(K2),IC(K3),IC(K4),CC(K6)
     .             ,JOUT,IC(K5),NO,IC(KTPO),CC(K7),IC(KTYP),
     .             NIRRED,IC(KFLOV),IC(KORBSM),NTRI,IC(KNMSM),
     .             IC(KIWRK),IC(KKLPT),IC(KIJPT),IC(KOFF1),
     .             IC(KOFF2))
      NV=NBF-NO
C
C >>> SET FROZEN ARRAYS AND ACTIVE ORBITALS PER SYMMETRY
C
      nlamda=k3
      nc=k5
      CALL FROZEN (ic(kFLOV),ic(k3),ic(k5),NFZO,NFZV,NDOC,NUOC,IC(FZT),
     .             ic(kTYP),NIRRED,NO,NV,JOUT,NOACT,NVACT,NSYMHF)
      LFZO=NO-NOACT
      LFZV=NV-NVACT
C
C     WRITE(*,*)  '  NV = ',NV
C     WRITE(JOUT,*)
C     WRITE(JOUT,*) ' THE SCF EIGENVECTORS '
C     WRITE(JOUT,*)
C     CALL MATOUT (CC(K1),NBFAO,NBF,NBFAO,NBF,JOUT)
C     WRITE (JOUT,1)
C     WRITE (JOUT,2)
C     WRITE (JOUT,*) ' ************  EIGENVECTOR  *************'
C     WRITE (JOUT,2)
C     WRITE (JOUT,1)
C     DO 4 I=1,196
C        WRITE (JOUT,*) I,CC(I)
C   4 CONTINUE
C     WRITE (JOUT,1)
C     WRITE (JOUT,2)
C     WRITE (JOUT,*) ' ************  EIGENVALUES  *************'
C     WRITE (JOUT,2)
C     WRITE (JOUT,1)
C     DO 5 I=197,210
C        WRITE (JOUT,*) I,CC(I)
C   5 CONTINUE
C     WRITE (JOUT,1)
C     WRITE (JOUT,2)
C     WRITE (JOUT,*) ' ************  NLAMDA  *************'
C     WRITE (JOUT,2)
C     WRITE (JOUT,1)
C     WRITE (JOUT,*) IC (421)
C     WRITE (JOUT,1)
C     WRITE (JOUT,2)
C     WRITE (JOUT,*) ' ************  ITEMP  *************'
C     WRITE (JOUT,2)
C     WRITE (JOUT,1)
C     DO 6 I=422,521
C        WRITE (JOUT,*) I,IC(I)
C   6 CONTINUE
C     WRITE (JOUT,1)
C     WRITE (JOUT,2)
C     WRITE (JOUT,*) ' ************  XVEC  *************'
C     WRITE (JOUT,2)
C     WRITE (JOUT,1)
C     DO 7 I=262,457
C        WRITE (JOUT,*) I,CC(I)
C   7 CONTINUE
C
      J1   = IADTWP(K4)
      J2   = J1 + NTRIAO
      J3   = J2 + NTRIAO
      J4   = J3 + NTRIAO
      J5   = J4 + NTRIAO
      J6   = WPADTI(J5)
      JTOP = J5 + IBUFHF
      WRITE (JOUT,*) ' before reading file34 JTOP=',JTOP
      IF (JTOP.GT.MAXCOR) THEN
         WRITE (JOUT,*) ' MAXCOR=',MAXCOR
         WRITE (JOUT,*) ' INCREASE SIZE OF MAXCOR'
         stop
      END IF
C
      CALL READ34 (CC(J1),CC(J2),CC(J3),CC(J4),CC(J5),IC(J6),
     .             IBUFHF,ITAP34,JOUT)
C
      LENINT=MAXCOR-IBUFHF-J2-3*NBFAO*NBFAO-(MAXBKT+1)/INTOWP(1)
C     LENINT=MAXCOR-INTBUF-J2-3*NBFAO*NBFAO-(MAXBKT+1)/INTOWP(1)
      ITRI=LENINT/NTRIAO
      WRITE (JOUT,*) ' LENINT=',LENINT
      WRITE (JOUT,*) ' NTRIAO=',NTRIAO
      WRITE (JOUT,*) ' ITRI=lenint/ntriao',ITRI
CGES
      IF(LENINT.LE.0)THEN
         WRITE(6,*)' LENINT.LE.0'
         WRITE(6,*)' INCREASE MAXCOR'
         stop
      ENDIF
CGES
         ISORT=2
         IF (ITRI.GT.NTRIAO) ITRI=NTRIAO
         NBUCK=NTRIAO/ITRI
         MODTRI=MOD(NTRIAO,ITRI)
         IF (ITRI*NBUCK.NE.NTRIAO) NBUCK=NBUCK+1
         IBKTSP=MAXCOR-IBUFHF-J2-(MAXBKT+1)/INTOWP(1)
         LENGTH=IBKTSP/NBUCK
         NIOBF=INTOWP(LENGTH)/SEC2I(1)
cges
         IF (NIOBF.GT.100) NIOBF=100
cges
         LENGTH=NIOBF*(SEC2I(1)/INTOWP(1))
         SPACE=MAXCOR-LENGTH-3*NBFAO*NBFAO-J2-(MAXBKT+1)/INTOWP(1)
         IF (NTRIAO*ITRI.GT.SPACE) THEN
            LENGTH=LENGTH-(NTRIAO*ITRI-SPACE)
            NIOBF=INTOWP(LENGTH)/SEC2I(1)
            LENGTH=NIOBF*(SEC2I(1)/INTOWP(1))
            SPACE=MAXCOR-LENGTH-3*NBFAO*NBFAO-J2-(MAXBKT+1)/INTOWP(1)
         END IF
         IF (SPACE.GT.NTRIAO*NTRIAO) SPACE=NTRIAO*NTRIAO
CGES
      IF(SPACE.LE.0)THEN
         WRITE(6,*)' SPACE.LE.0'
         WRITE(6,*)' INCREASE MAXCOR OR DECREASE NIOBF'
         stop
      ENDIF
CGES
C
         N1=J2
         IN1=WPADTI(N1)
         N2=IADTWP(IN1+MAXBKT)
         IN2=WPADTI(N2)
         N3=N2+NBUCK*LENGTH
         IN3=WPADTI(N3)
         NTOP=N3+IBUFHF-1
         WRITE (JOUT,*) ' before first sort nTOP=',nTOP
         WRITE (JOUT,*) '       NBUCK =',NBUCK
         WRITE (JOUT,*) '      LENGTH =',LENGTH
         WRITE (JOUT,*) '      IBKTSP =',IBKTSP
         WRITE (JOUT,*) '       SPACE =',SPACE
         WRITE (JOUT,*) '       NIOBF =',NIOBF
         WRITE (JOUT,*) '        ITRI =', ITRI
         WRITE (JOUT,*) '      MAXBKT =', MAXBKT
         IF (nTOP.GT.MAXCOR) THEN
             WRITE (JOUT,*) ' MAXCOR=',MAXCOR
             WRITE (JOUT,*) ' INCREASE SIZE OF MAXCOR'
             stop
         END IF
         CALL ZERO(CC(N1),MAXCOR-N1)
c
        CALL SORT (CC(N2),IC(IN2),IC(IN3),CC(N3),NBUCK*LENGTH,
     .             IBUFHF,NBUCK,JOUT,ITAP34,ITRI,ITAP77,IC(IN1),
     .             MAXBKT,NIOBF,LENGTH,MAXVAL,IWRIT)
c
         WRITE (6,*) ' exiting sort 1 : file34 --> file77'                
cges     CALL RCLOSE(ITAP34,3)
         CALL RCLOSE(ITAP34,4)
C
         N6=N2
         N7=N6+SPACE
         IN7=WPADTI(N7)
         N8=N7+LENGTH
         N9=N8+NBFAO*NBFAO
         N10=N9+NBFAO*NBFAO
         NTOP=N10+NBFAO*NBFAO-1  
         WRITE (JOUT,*) '  before 1st htrans  nTOP=',nTOP
         IF (nTOP.GT.MAXCOR) THEN
             WRITE (JOUT,*) ' MAXCOR=',MAXCOR
             WRITE (JOUT,*) ' INCREASE SIZE OF MAXCOR'
             stop
         END IF
         CALL ZERO(CC(N6),MAXCOR-N6)
         JBUF=0
         DO 101 I=1,NBUCK
C           WRITE (JOUT,*) ' ENTERING BCKCHN 1'
C           WRITE (JOUT,*) ' NBUCK=',I
C           WRITE (JOUT,*) ' LENGTH=',LENGTH
C           WRITE (JOUT,*) ' SPACE=',SPACE
C           WRITE (JOUT,*) ' MAXBKT=',MAXBKT
            CALL BCKCHN(IC(IN1),I,ITAP77,CC(N7),IC(IN7),LENGTH,
     .                  MAXBKT,MAXVAL,CC(N6),SPACE,NTRIAO,JOUT,ITRI)
            NTRIMX=ITRI
            IF ((I.EQ.NBUCK).AND.(MODTRI.NE.0)) NTRIMX=MODTRI
            IFLAG=2
            CALL HTRANS (CC(N8),CC(N9),CC(N10),CC(N6),
     .                   CC(N7),IC(IN7),SPACE,
     .                   LENGTH,NBFAO,CC(K1),NBFAO,NBF,JOUT,NTRIAO,
     .                   NTRI,IFLAG,ISORT,ITAP78,NTRIMX,I,ITRI,JBUF,
     .                   NBUCK,IC(KORBSM),IC(KFLOV),NIRRED,IC(K3),
     .                   NSYMHF,IC(KKLPT),IC(KIJPT),NTRI,IC(KTPO),
     .                   IC(KOFF1),IC(KOFF2),NO,nfzo,nfzv)
            CALL ZERO (CC(N6),MAXCOR-N6)
  101    CONTINUE
         CALL RCLOSE(ITAP77,4)
         WRITE (6,*) ' exiting 1st half trans : file77 --> file78'
c
         CALL RFILE(ITAP77)
         CALL ZERO(CC(J2),MAXCOR-J2)
         NBUCK=NTRI/ITRI
         MODTRI=MOD(NTRI,ITRI)
         IF (ITRI*NBUCK.NE.NTRI) NBUCK=NBUCK+1
         IBKTSP=MAXCOR-LENGTH-J2-(MAXBKT+1)/INTOWP(1)
         LEN2=IBKTSP/NBUCK
         NIOBF=INTOWP(LEN2)/SEC2I(1)
c
         IF (NIOBF.GT.100) NIOBF=100
c
         LEN2=NIOBF*(SEC2I(1)/INTOWP(1))
         SPACE2=MAXCOR-LEN2-3*NBF*NBFAO-J2-(MAXBKT+1)/INTOWP(1)
         IF (NTRIAO*ITRI.GT.SPACE2) THEN
            LEN2=LEN2-(NTRIAO*ITRI-SPACE2)
            NIOBF=INTOWP(LEN2)/SEC2I(1)
            LEN2=NIOBF*(SEC2I(1)/INTOWP(1))
         SPACE2=MAXCOR-LEN2-3*NBF*NBFAO-J2-(MAXBKT+1)/INTOWP(1)
         END IF
         IF (SPACE2.GT.NTRI*NTRIAO) SPACE2=NTRI*NTRIAO
CGES
      IF(SPACE2.LE.0)THEN
         WRITE(6,*)' SPACE2.LE.0'
         WRITE(6,*)' INCREASE MAXCOR OR DECREASE NIOBF'
         stop
      ENDIF
CGES
C
         N1=J2
         IN1=WPADTI(N1)
         N2=IADTWP(IN1+MAXBKT)
         IN2=WPADTI(N2)
         N3=N2+LEN2*NBUCK
         IN3=WPADTI(N3)
         NTOP=N3+LENGTH-1
         WRITE (JOUT,*) ' before 2nd   sort nTOP=',nTOP
         WRITE (JOUT,*) '       NBUCK =',NBUCK
         WRITE (JOUT,*) '       LEN2  =',LEN2  
         WRITE (JOUT,*) '      IBKTSP =',IBKTSP
         WRITE (JOUT,*) '      SPACE2 =',SPACE2
         WRITE (JOUT,*) '       NIOBF =',NIOBF
         WRITE (JOUT,*) '        ITRI =', ITRI
         WRITE (JOUT,*) '      MAXBKT =', MAXBKT
         IF (nTOP.GT.MAXCOR) THEN
             WRITE (JOUT,*) ' MAXCOR=',MAXCOR
             WRITE (JOUT,*) ' INCREASE SIZE OF MAXCOR'
             stop
         END IF
c
         CALL ZERO(CC(N1),MAXCOR-N1)
         CALL SORT2 (ITAP77,ITAP78,CC(N3),IC(IN3),LENGTH,CC(N2),
     .               IC(IN2),NBUCK*LEN2,ITRI,IC(IN1),MAXBKT,
     .               JOUT,NBUCK,NIOBF,LEN2,JBUF,MAXVAL)
C
         CALL RCLOSE(ITAP78,4)
         WRITE (6,*) ' exiting sort 2 : file78 --> file77'                
c 
         CALL RFILE(ITAP78)
         CALL SWRIT(ITAP78,LEN2,1)
         CALL SREW(ITAP78)
         CALL SREAD(ITAP78,LEN2,1)
         N6=N2
         N7=N6+SPACE2
         IN7=WPADTI(N7)
         N8=N7+LEN2
         N9=N8+NBFAO*NBFAO
         N10=N9+NBFAO*NBFAO
         NTOP=N10+NBFAO*NBFAO-1
         WRITE (JOUT,*) '  before 2nd htrans  nTOP=',nTOP
         IF (nTOP.GT.MAXCOR) THEN
             WRITE (JOUT,*) ' MAXCOR=',MAXCOR
             WRITE (JOUT,*) ' INCREASE SIZE OF MAXCOR'
             stop
         END IF
c
         CALL ZERO(CC(N6),MAXCOR-N6)
         JBUF=0
         DO 111 I=1,NBUCK
C           WRITE (JOUT,*) ' LEN2=',LENGTH
C           WRITE (JOUT,*) ' SPACE2=',SPACE2
C           WRITE (JOUT,*) ' MAXBKT=',MAXBKT
            CALL BCKCHN(IC(IN1),I,ITAP77,CC(N7),IC(IN7),LEN2,
     .                  MAXBKT,MAXVAL,CC(N6),SPACE2,NTRIAO,JOUT,ITRI)
            NTRIMX=ITRI
            IF ((I.EQ.NBUCK).AND.(MODTRI.NE.0)) NTRIMX=MODTRI
            IFLAG=1
            CALL SETMBF(N8,N7+LEN2)
            CALL SETMBF(N9,N8+NBFAO*NBFAO)
            CALL SETMBF(N10,N9+NBFAO*NBFAO)
            CALL SETMBF(NTOP,N10+NBFAO*NBFAO)
C           WRITE(*,*) '  THIRD NTOP = ',NTOP
            CALL HTRANS (CC(N8),CC(N9),CC(N10),CC(N6),
     .                   CC(N7),IC(IN7),
     .                   SPACE2,LEN2,NBFAO,CC(K1),NBFAO,NBF,JOUT,NTRI,
     .                   NTRI,IFLAG,ISORT,ITAP78,NTRIMX,I,ITRI,JBUF,
     .                   NBUCK,IC(KORBSM),IC(KFLOV),NIRRED,IC(K3),
     .                   NSYMHF,IC(KKLPT),IC(KIJPT),NTRI,IC(KTPO),
     .                   IC(KOFF1),IC(KOFF2),NO,nfzo,nfzv)
            CALL ZERO(CC(N6),MAXCOR-N6)
  111    CONTINUE
         CALL RCLOSE(ITAP77,4)
         WRITE (6,*) ' exiting 2nd half trans file77 --> file78'                
c
      I1=J2
      II1=WPADTI(I1)
      I2=I1+LEN2
      I3=I2+NBFAO*NBFAO
      I4=I3+NBFAO*NBFAO
      ITOP=I4+NBFAO*NBFAO  
C
      CALL ZERO(CC(I1),MAXCOR-I1)
      ISORT=1
      IFLAG=1
      CALL HTRANS (CC(I2),CC(I3),CC(I4),CC(J1),
     .             CC(I1),IC(II1),NTRIAO*1,
     .             LEN2,NBFAO,CC(K1),NBFAO,NBF,JOUT,1,1,IFLAG,
     .             ISORT,ITAP78,1,1,1,JBUF,1,IC(KORBSM),IC(KFLOV),
     .             NIRRED,IC(K3),NSYMHF,IC(KKLPT),IC(KIJPT),NTRI,
     .             IC(KTPO),IC(KOFF1),IC(KOFF2),NO,nfzo,nfzv)
C
      write(6,*)' exit transformation of one-e ints'
      CALL RCLOSE(ITAP30,3)
      CALL RCLOSE(ITAP78,3)
      CALL TSTOP(6)
C
      END
