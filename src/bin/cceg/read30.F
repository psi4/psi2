C
C----------------------------------------------------------------------
C
      SUBROUTINE READ30(ITAP30,MPOINT,MCONST,MCALCS,NCALCS,NT,
     .                  NSYMHF,MXCOEF,EIG,NLAMDA,JOUT,NCL,NO,
     .                  WTEMP,SCFTCC,ITYP,ORBSYM,FLOV,NIRRED,
     .                  ENUC,ESCF,NAT,NSHELL,NPRIM,
     .                  CHAR,NORD,ZAN,COORD,
     .                  IOPEN,NOP,dist)
      IMPLICIT INTEGER(A-Z)
      CHARACTER*4 ITYP(NSYMHF),ICSYM
      CHARACTER*4 CHAR,LABEL(20)
      REAL*8 ENUC,ESCF,EE
      real*8 distx,disty,distz,dist(*)
      REAL*8 EIG(NT),WTEMP(NT),ZAN(NAT),COORD(3,NAT),A30(100)
      INTEGER NLAMDA(NSYMHF),NCL(NSYMHF),SCFTCC(NT),
     .        ORBSYM(NT),FLOV(NIRRED,4),SYMORB,NOP(8),
     .        ITEMP(100),IA30(200),IALBEL(20)
      integer errcod,frdc
      EQUIVALENCE (IA30,A30)
      nat2=nat*(nat-1)/2
C
C     MCALCS=100, MPOINT=200, MCONST=200
C
      JPOINT = 101 + MCONST + MPOINT
      CALL WREADW (ITAP30,ITEMP,MCALCS,JPOINT,JPOINT)
      LOCCAL = ITEMP(NCALCS)
      JPOINT = LOCCAL + 60
      CALL WREADW (ITAP30,LOCVEC,1,JPOINT,JPOINT)
      LOCVEC = LOCVEC + INTOWP(MXCOEF)
      CALL WREADW (ITAP30,EIG,INTOWP(NT),LOCVEC,LOCVEC)
      CALL WREADW (ITAP30,ITYP,NSYMHF,LOCVEC,LOCVEC)
      CALL WREADW (ITAP30,NLAMDA,NSYMHF,LOCVEC,LOCVEC)
      CALL WREADW (ITAP30,NCL,NSYMHF,LOCVEC,LOCVEC)
      IF(IOPEN.NE.0) THEN
         CALL WREADW (ITAP30,NOP,NSYMHF,LOCVEC,LOCVEC)
      ENDIF
C
C >>> CALCULATE THE NUMBER OF OCCUPIED ALPHA AND BETA ORBITALS
C
      NOB=0
      NOA=0
      DO 125 I=1,NSYMHF
         NOA=NOA+NCL(I)+NOP(I)
         NOB=NOB+NCL(I)
C        WRITE (6,*) ' I=',I,' NCL(I)=',NCL(I),' NOP(I)=',NOP(I)
  125 CONTINUE
      NO=NOB
C
C >>> FORM SCF TO CC ORDERING ARRAY
C
      ICNT=0
      IOF=0
      DO 258 I=1,NSYMHF
         IF(I.NE.1) IOF=IOF+NLAMDA(I-1)
         NOI=NCL(I)
         DO 257 J=1,NOI
            ICNT=ICNT+1
            IPT=IOF+J
            SCFTCC(IPT)=ICNT
  257    CONTINUE
  258 CONTINUE
C
      IOF = 0
      DO 358 I=1,NSYMHF
         IF(I.NE.1) IOF=IOF+NLAMDA(I-1)
         NVI=NLAMDA(I)-NCL(I)
         DO 357 J=1,NVI
            ICNT=ICNT+1
            IPT=IOF+NCL(I) + J
            SCFTCC(IPT)=ICNT
  357    CONTINUE
  358 CONTINUE
C
C >>> SORT THE MO ENERGIES TO CC ORDERING
C
      DO 459 I=1,NT
         IPT=SCFTCC(I)
         WTEMP(IPT)=EIG(I)
  459 CONTINUE
C
      DO 469 I=1,NT
         EIG(I)=WTEMP(I)
C        WRITE(6,*)'I=',I,'E(I)=',EIG(I)
  469 CONTINUE
C
C >>> GET SYMMETRY LABEL FROM INPUT
C
      CALL LOCATE (5,'# INPUT ###',IERR)
      IF(IERR.NE.0) THEN
         errcod = frdc('SYMMETRY',char)
         if (errcod.ne.0) then
           write(6,*) 'ERROR: problem reading ''SYMMETRY'' from input'
           call pabort
           call qabort
           endif
         call consym(char,nord)
      ELSE
         write(6,*) 'WARNING: using ''# INPUT ###'' for symmetry'
         READ(5,6004)LABEL
         READ(5,6004)LABEL
         READ(5,6005)CHAR,NORD
      END IF
C
C >>> FORM FIRST & LAST MO ARRAY FOR EACH IRREP
C
      DO 805 ISYM = 1,NIRRED
         FLOV(ISYM,1) = 0
         FLOV(ISYM,2) = -1
         FLOV(ISYM,3) = 0
         FLOV(ISYM,4) = -1
  805 CONTINUE
C
C >>>  AND FORM ORBITAL SYMMETRY ARRAY
C
      ICNT=0
      DO 705 ISYM=1,NSYMHF
         ICNTH=ICNT+1
         NOI=NCL(ISYM)
         DO 702 I=1,NOI
            ICNT=ICNT+1
            ORBSYM(ICNT)=ISYM-1
  702    CONTINUE
         FLOV(ISYM,1)=ICNTH
         FLOV(ISYM,2)=ICNT
  705 CONTINUE
      DO 715 ISYM=1,NSYMHF
         ICNTH=ICNT+1
         NVI=NLAMDA(ISYM)-NCL(ISYM)
         DO 712 I=1,NVI
            ICNT=ICNT+1
            ORBSYM(ICNT)=ISYM-1
  712    CONTINUE
         FLOV(ISYM,3)=ICNTH
         FLOV(ISYM,4)=ICNT
  715 CONTINUE
C
C >>> GET ENERGIES, COORDINATES AND ATOMIC NUMBERS FROM FILE30
C
      JUNK = 101 + MCONST
      CALL WREADW(ITAP30,IA30,MPOINT,JUNK,JUNK)
      CALL WREADW(ITAP30,ZAN, INTOWP(NAT),IA30(1),JUNK)
      JUNK = 101 + MCONST + MPOINT + NCALCS - 1
      CALL WREADW(ITAP30,LOCCAL,1,JUNK,JUNK)
      CALL WREADW(ITAP30,IA30,60,LOCCAL,LOCCAL)
      DO 500  ILBL=1,20
         IALBEL(ILBL)=IA30(ILBL)
  500 CONTINUE
      WRITE(6,*)
      WRITE(6,*)' >>>  LABEL FROM FILE30  <<<'
      WRITE(6,905)IALBEL
  905 FORMAT(2X,20A4)
      LOCCAL=LOCCAL+20
      CALL WREADW(ITAP30,COORD,INTOWP(3*NAT),LOCCAL,LOCCAL)
      WRITE(6,*)
      WRITE(6,*)' >>>  GEOMETRY FROM FILE30  <<<'
      DO 523 I=1,NAT
      WRITE(6,909)I,ZAN(I),COORD(1,I),COORD(2,I),COORD(3,I)
      WRITE(10,909)I,ZAN(I),COORD(1,I),COORD(2,I),COORD(3,I)
  909 FORMAT(1X,I4,3X,F5.2,5X,3(4X,F13.10))
  523 CONTINUE
c
c     ij=0
c     do 525 i=1,nat
c        do 524 j=1,i-1
c           distx=(coord(1,i)-coord(1,j))**2.0d0
c           disty=(coord(2,i)-coord(2,j))**2.0d0
c           distz=(coord(3,i)-coord(3,j))**2.0d0
c           ij=ij+1
c           dist(ij) =dsqrt(distx+disty+distz)*0.529177d0
c 524    continue
c 525 continue
c           write(10,6078)(dist(ij),ij=1,nat2)
c6078 format(6f12.8)
C
C >>> GET ENERGIES AND WRITE OUT
C
      CALL WREADW(ITAP30,IA30,INTOWP(10),LOCCAL,LOCCAL)
      ESCF=A30(2)
      ENUC=A30(1)
      EE=ESCF-ENUC
      WRITE(6,*)
      WRITE(6,*)' >>>  ENERGIES FROM FILE30  <<<'
      WRITE(6,901) ENUC
  901 FORMAT(1X,' NUCLEAR REPULSION ENERGY =',F20.12)
      WRITE(6,902) ESCF
  902 FORMAT(1X,' SCF ENERGY               =',F20.12)
      WRITE(6,903) EE
  903 FORMAT(1X,' TOTAL ELECTRONIC  ENERGY =',F20.12)
C
 6004 FORMAT(20A4)
 6005 FORMAT(A4,5X,I1)
C
C >>> WRITE OUT MO ENERGIES FROM FILE30
C
      WRITE(6,*)
      WRITE(6,*)' >>>  SCF EIGENVALUES FROM FILE30  <<<'
 6021 FORMAT(4(I3,1X,A3,1X,F10.6,1X))
      WRITE(6,6021) (I,ITYP(ORBSYM(I)+1),EIG(I),I=1,NT)
C
      CALL RCLOSE(ITAP30,3)
C
      RETURN
      END
