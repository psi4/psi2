c************************************************************
c
      SUBROUTINE BASSET
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION I30(200),A30(100),CONTR(1000)
      COMMON/COORD/CHGS(112),XX(3,112),ict(112),ns(112),ks(112)
      common/funcs/natom,n3n,nbfao,nbatri,nocc,mb(20)
      COMMON/GAUSS/ZS(200),CS(200)
      COMMON/NSHEL/KNUC(200),KTYPE(200),KPRIM(200),KLOC(200),
     1 KMIN(200),KMAX(200),KNG(200),NSHELL,NPRIM
      COMMON/SIGNS/IPQ(500),IPRNT
      EQUIVALENCE (I30,A30)
      DATA ZERO0 / 0.0D+00 /
    1 FORMAT(//,2X,' NPRIM  IS ',I5)
    2 FORMAT(2X,I5,2F20.10)
    3 FORMAT(//,2X,' NSHELL IS ',I5)
    4 FORMAT(2X,8I5)
    5 FORMAT(//,2X,' NATOM IS ',I5)
    6 FORMAT(2X,I5,3F16.10,F8.4,5X,2I5)
    8 FORMAT(/,2X,' PARAMETERS'/
     * 2X,' NBFAO  = ',I8/
     * 2X,' NBATRI = ',I8/
     * 2X,' NATOM  = ',I8/
     * 2X,' N3N    = ',I8/
     * 2X,' NSHELL = ',I8/
     * 2X,' NPRIM  = ',I8/
     * 2X,' IPRNT  = ',I8/)
    9 FORMAT(//,3X,20A4)
   10 FORMAT(//,2X,' ENUC  = ',F20.10/
     1          2X,' ESCF  = ',F20.10/)
C
C
      ITAP30=30
c
      ipq(1)=0
      DO 100 i=1,499
  100 ipq(i+1)=ipq(i)+i
c
      CALL RFILE(ITAP30)
      CALL WREADW(ITAP30,I30,200,101,JUNK)
      IEND=I30(1)
      MPOINT=I30(2)
      MCONST=I30(3)
      MCALCS=I30(4)
      NCALCS=I30(5)
      NATOM=I30(19)
      NBFAO=I30(22)
      NSHELL=I30(27)
      NT=I30(29)
      NPRIM=I30(32)
      N3N=NATOM*3
      nbatri=((nbfao+1)*nbfao)/2
c
      if(iprnt.gt.0)write(6,8) nbfao,nbatri,natom,n3n,nshell,
     1    nprim,iprnt
C
C   READ IN LABELS
      JUNK=101+MCONST+MPOINT+NCALCS-1
      CALL WREADW(ITAP30,LOCCAL,1,JUNK,JUNK)
      CALL WREADW(ITAP30,I30,60,LOCCAL,LOCCAL)
      WRITE(6,9) (I30(I),I=1,20)
C
C   READ IN GEOMETRY
      LOCCAL=LOCCAL+20
      CALL WREADW(ITAP30,XX,intowp(NATOM*3),LOCCAL,LOCCAL)
C
C   READ IN ENUC AND ESCF
      CALL WREADW(ITAP30,I30,20,LOCCAL,LOCCAL)
      ENUC=A30(1)
      ESCF=A30(2)
      WRITE(6,10) ENUC,ESCF
c
      JUNK=101+MCONST
      CALL WREADW(ITAP30,I30,MPOINT,JUNK,JUNK)
C
      CALL ZERO(CONTR,1000)
      ITYFC =  I30(7)-I30(6)
C
      WRITE(3,*) ' I30(6)   = ',I30(6)
      WRITE(3,*) ' I30(7)   = ',I30(7)
      WRITE(3,*) ' ITYFC    = ',ITYFC
      WRITE(3,*) ' NPRIM    = ',NPRIM
      WRITE(3,*) ' NPRIM*6  = ',NPRIM*6
C
C   READ IN BASIS SET
      CALL WREADW(ITAP30,CHGS,intowp(NATOM),I30(1),JUNK)
      CALL WREADW(ITAP30,ICT,NT*NATOM,I30(2),JUNK)
      CALL WREADW(ITAP30,NS,NATOM,I30(3),JUNK)
      CALL WREADW(ITAP30,KS,NATOM,I30(4),JUNK)
      CALL WREADW(ITAP30,ZS,intowp(NPRIM),I30(5),JUNK)
      CALL WREADW(ITAP30,CONTR,ITYFC,I30(6),JUNK)
      CALL WREADW(ITAP30,KPRIM,NSHELL,I30(7),JUNK)
      CALL WREADW(ITAP30,KNUC,NSHELL,I30(8),JUNK)
      CALL WREADW(ITAP30,KTYPE,NSHELL,I30(9),JUNK)
      CALL WREADW(ITAP30,KNG,NSHELL,I30(10),JUNK)
      CALL WREADW(ITAP30,KLOC,NSHELL,I30(11),JUNK)
      CALL WREADW(ITAP30,KMIN,NSHELL,I30(12),JUNK)
      CALL WREADW(ITAP30,KMAX,NSHELL,I30(13),JUNK)

      write(6,*) 'The geometry and charges:'
      DO 105 I=1,NATOM
  105 WRITE(6,6) I,XX(1,I),XX(2,I),XX(3,I),CHGS(I)
C
C   MB(I) is the # of basis functions on the Ith atom.
      do 80 i=1,natom
      j1=ks(i)
      j2=ks(i)+ns(i)-1
      mb(i)=0
      do 70 j=j1,j2
      go to (61,62,63,64,65) ktype(j)
   61 ktn=1
      go to 66
   62 ktn=3
      go to 66
   63 ktn=6
      go to 66
   64 ktn=10
      go to 66
   65 ktn=15
   66 mb(i)=mb(i)+ktn
   70 continue
   80 continue
c
      IJ=0
      DO 102 I=1,5
      DO 101 J=1,NPRIM
      IJ=IJ+1
      VALU=CONTR(IJ)
      IF(VALU.EQ.ZERO0) GO TO 101
      CS(J)=VALU
  101 CONTINUE
  102 CONTINUE
C
      if(iprnt.gt.0)then
      WRITE(6,1) NPRIM
      DO 103 I=1,NPRIM
  103 WRITE(6,2) I,ZS(I),CS(I)
      WRITE(6,3) NSHELL
      DO 104 I=1,NSHELL
  104 WRITE(6,4) I,KNUC(I),KTYPE(I),KPRIM(I),KLOC(I),KMIN(I),KMAX(I),
     1 KNG(I)
      WRITE(6,5) NATOM
      endif
C
      END
