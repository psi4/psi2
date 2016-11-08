      SUBROUTINE READ30(ITAP30,MPOINT,MCONST,MCALCS,NCALCS,NT,
     .                  NSYMHF,MXCOEF,EIG,NLAMDA,JOUT,NCL,NO,
     .                  WTEMP,SCFTCC,ITYP,ORBSYM,FLOV,NIRRED,
     .                  ENUC,ESCF,NAT,NSHELL,NPRIM,
     .                  CHAR,NORD,ZAN,COORD,cctscf,
     .                  IOPEN,NOP,map,spin,ntmo,mosym,dist,
     .                  nocs,nvcs,f1,g1,omap,nalp)
      IMPLICIT INTEGER(A-Z)
      CHARACTER*4 ITYP(NSYMHF),ICSYM
      CHARACTER*4 CHAR,LABEL(20)
      REAL*8 ENUC,ESCF,EE,distx,disty,distz,dist(*)
      REAL*8 EIG(nt),WTEMP(nt),ZAN(NAT),COORD(3,NAT),A30(100)
      INTEGER NLAMDA(*),NCL(*),SCFTCC(nt),cctscf(nt),
     .        ORBSYM(NT),FLOV(8,2,4),SYMORB,NOP(*),
     .        ITEMP(100),IA30(200),IALBEL(20),
     .        map(*),spin(*),mosym(*),f1(nt),g1(nt),omap(nt)
      integer errcod,frdc
      EQUIVALENCE (IA30,A30)
c
c     most arrays are nt long but only ntmo are needed in many cases
c
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
      CALL WREADW (ITAP30,EIG,INTOWP(NTmo),LOCVEC,LOCVEC)
      CALL WREADW (ITAP30,ITYP,NSYMHF,LOCVEC,LOCVEC)
      CALL WREADW (ITAP30,NLAMDA,NSYMHF,LOCVEC,LOCVEC)
      CALL WREADW (ITAP30,NCL,NSYMHF,LOCVEC,LOCVEC)
      IF(IOPEN.NE.0) THEN
         CALL WREADW (ITAP30,NOP,NSYMHF,LOCVEC,LOCVEC)
      ENDIF
C
C >>> CALCULATE THE NUMBER OF OCCUPIED ALPHA AND BETA ORBITALS
C
      write(6,*)
      write(6,*)'occupations from file30'
      NOB=0
      NOA=0
      nalp=0
      nt30=0
      DO 125 I=1,NSYMHF
         nt30=nt30+nlamda(i)
         NOA=NOA+NCL(I)+NOP(I)
         NOB=NOB+NCL(I)
         nalp=nalp+nop(i)
         WRITE (6,*) ' I=',I,' NCL(I)=',NCL(I),' NOP(I)=',NOP(I)
  125 CONTINUE
      NO=noa+NOB
      if(nt30.ne.ntmo)then
         write(6,*)' nt30.ne.ntmo  error... '
         stop
      endif
      nv=nt-no
c nalp is the number of open-shell orbitals
c number of closed-shell beta orbitals is equal to # doc beta orbitals
      NOcs=NOB
      nvcs=ntmo-nocs-nalp
      write(6,*)
      write(6,*)'NOTE ! nocs+nalp+nvcs=ntmo '
      write(6,*)'nocs=',nocs
      write(6,*)'nalp=',nalp
      write(6,*)'nvcs=',nvcs
      write(6,*)'ntmo=',ntmo
c
      write(6,*)' '
      write(6,*)' the total number occ sos  =',no
      write(6,*)' the total number vir sos  =',nv
      write(6,*)' the total number  of sos  =',nt
      write(6,*)' '
C
C >>> FORM SCF TO CC ORDERING ARRAY
c     this will implicitily map open-shells as virtual mo's
C
      ICNT=0
      IOF=0
      DO 258 I=1,NSYMHF
         IF(I.NE.1) IOF=IOF+NLAMDA(I-1)
         NOI=NCL(I)
         DO 255 J=1,NOI
            ICNT=ICNT+1
            IPT=IOF+J
            SCFTCC(IPT)=ICNT
  255    CONTINUE
  258 CONTINUE
c
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
c
c >>> SORT THE MO ENERGIES TO CC ORDERING
C
      DO 459 I=1,ntmo
         IPT=SCFTCC(I)
         cctscf(ipt)=i
         WTEMP(IPT)=EIG(I)
  459 CONTINUE
C
      DO 469 I=1,ntmo
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
C >>> FORM FIRST & LAST MO ARRAY FOR EACH IRREP and spin
C
      do 804 ispin=1,2
         DO 805 ISYM = 1,NIRRED
            FLOV(ISYM,ispin,1) = 0
            FLOV(ISYM,ispin,2) = -1
            FLOV(ISYM,ispin,3) = 0
            FLOV(ISYM,ispin,4) = -1
  805    CONTINUE
  804 CONTINUE
C
C >>>  AND FORM ORBITAL SYMMETRY ARRAY
C
      ICNT=0
      mocnt=0
c     write(6,*)'occ alpha'
      DO 705 ISYM=1,NSYMHF
         ICNTH=ICNT+1
         NOI=NCL(ISYM)
         nai=nop(isym)
         DO 702 I=1,NOI+nai
            ICNT=ICNT+1
            ORBSYM(ICNT)=ISYM-1
  702    CONTINUE
         FLOV(ISYM,1,1)=ICNTH
         FLOV(ISYM,1,2)=ICNT
c      write(6,*)'isym',isym,icnth,icnt
  705 CONTINUE
c
c     write(6,*)'occ beta '
      DO 725 ISYM=1,NSYMHF
         ICNTH=ICNT+1
         NOI=NCL(ISYM)
         DO 703 I=1,NOI
            ICNT=ICNT+1
            mocnt=mocnt+1
            ORBSYM(ICNT)=ISYM-1
            mosym(mocnt)=isym-1
  703    CONTINUE
         FLOV(ISYM,2,1)=ICNTH
         FLOV(ISYM,2,2)=ICNT
c      write(6,*)'isym',isym,icnth,icnt
  725 CONTINUE
c
c     write(6,*)'vir alpha'
      DO 715 ISYM=1,NSYMHF
         ICNTH=ICNT+1
         NVI=NLAMDA(ISYM)-NCL(ISYM)-nop(isym)
         DO 712 I=1,NVI
            ICNT=ICNT+1
            ORBSYM(ICNT)=ISYM-1
  712    CONTINUE
         FLOV(ISYM,1,3)=ICNTH
         FLOV(ISYM,1,4)=ICNT
c      write(6,*)'isym',isym,icnth,icnt
  715 CONTINUE
c
c     write(6,*)'vir beta '
      DO 735 ISYM=1,NSYMHF
         ICNTH=ICNT+1
         NVI=NLAMDA(ISYM)-NCL(ISYM)
         DO 713 I=1,NVI
            ICNT=ICNT+1
            mocnt=mocnt+1
            ORBSYM(ICNT)=ISYM-1
            mosym(mocnt)=isym-1
  713    CONTINUE
         FLOV(ISYM,2,3)=ICNTH
         FLOV(ISYM,2,4)=ICNT
c      write(6,*)'isym',isym,icnth,icnt
  735 CONTINUE
c     
c     do 998 isym=1,nsymhf
c 998 write(6,*)'isym=',isym,'flocc,alpha',flov(isym,1,1),flov(isym,1,2)
c     do 999 isym=1,nsymhf
c 999 write(6,*)'isym=',isym,'flocc, beta',flov(isym,2,1),flov(isym,2,2)
c     do 996 isym=1,nsymhf
c 996 write(6,*)'isym=',isym,'flvir,alpha',flov(isym,1,3),flov(isym,1,4)
c     do 997 isym=1,nsymhf
c 997 write(6,*)'isym=',isym,'flvir, beta',flov(isym,2,3),flov(isym,2,4)
C
C >>>  form map and spin
C
      ICNT=0
      shift=0
      opcnt=0
      rcnt=0
      DO 755 ISYM=1,NSYMHF
         NOI=NCL(ISYM)
         DO 752 I=1,NOI
            ICNT=ICNT+1
            rcnt=rcnt+1
            omap(ICNT)=rcnt
            map(icnt)=rcnt+nvcs
            f1(map(icnt))=omap(icnt)
            g1(omap(icnt))=map(icnt)
            spin(icnt)=1
  752    CONTINUE
         ocnt=0
         nai=nop(isym)
         DO 753 I=1,NaI
            ICNT=ICNT+1
            ocnt=ocnt+1
            opcnt=opcnt+1
            omap(ICNT)=nob+shift+ocnt
            map(icnt)=nvcs+nocs+opcnt
            f1(map(icnt))=omap(icnt)
            g1(omap(icnt))=map(icnt)
            spin(icnt)=1
  753    CONTINUE
         shift=shift+nlamda(isym)-ncl(isym)
  755 CONTINUE
      rcnt=0
      DO 756 ISYM=1,NSYMHF
         NOI=NCL(ISYM)
         DO 754 I=1,NOI
            ICNT=ICNT+1
            rcnt=rcnt+1
            omap(ICNT)=rcnt
            map(icnt)=rcnt+nvcs
            spin(icnt)=2
  754    CONTINUE
  756 CONTINUE
c
      shift=0
      rcnt=0
      DO 765 ISYM=1,NSYMHF
         NvI=nlamda(isym)-NCL(ISYM)-nop(isym)
         shift=shift+nop(isym)
         DO 762 I=1,NvI
            ICNT=ICNT+1
            rcnt=rcnt+1
            omap(ICNT)=rcnt+shift+nob
            map(icnt)=rcnt
            f1(map(icnt))=omap(icnt)
            g1(omap(icnt))=map(icnt)
            spin(icnt)=1
  762    CONTINUE
  765 CONTINUE
      rcnt=0
      DO 766 ISYM=1,NSYMHF
         NvI=nlamda(isym)-NCL(ISYM)
         DO 764 I=1,NvI
            ICNT=ICNT+1
            rcnt=rcnt+1
            omap(ICNT)=rcnt+nob
            map(ICNT)=g1(omap(icnt))
            spin(icnt)=2
  764    CONTINUE
  766 CONTINUE
c
c     do 769 i=1,nt
c        write(6,675)i,spin(i),omap(i),map(i)
c 769 continue
c 675 format('i=',i4,'   spin=',i4,'    omap=',i4,'     map=',i4)
c
c     write(6,*)'i,f1(i),g1(i)'
c     do 1390 i=1,ntmo
c        write(6,'(3i4)')i,f1(i),g1(i)
c1390 continue
c
c >>> re-map map using g1
c
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
      WRITE(6,*)
      WRITE(6,*)' >>>  LABEL FROM FILE30  <<<'
      WRITE(6,*)
      WRITE(6,905)IALBEL
  905 FORMAT(2X,20A4)
      LOCCAL=LOCCAL+20
      CALL WREADW(ITAP30,COORD,INTOWP(3*NAT),LOCCAL,LOCCAL)
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)' >>>  GEOMETRY FROM FILE30  <<<'
      WRITE(6,*)
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
      WRITE(6,*)
      WRITE(6,*)' >>>  ENERGIES FROM FILE30  <<<'
      WRITE(6,*)
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
      WRITE(6,*)
      WRITE(6,*)' >>>  SCF EIGENVALUES FROM FILE30  <<<'
      WRITE(6,*)
 6021 FORMAT(4(I3,1X,A3,1X,F10.6,1X))
      WRITE(6,6021) (I,ITYP(mosym(I)+1),EIG(I),I=1,NTmo)
C
      CALL RCLOSE(ITAP30,3)
C
      RETURN
      END
