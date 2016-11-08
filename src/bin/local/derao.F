c
c*********************************************************************
c
      subroutine derao(s)
c
      implicit real*8 (a-h,o-z)
      logical iandj
      parameter(jnmax=36)
      dimension s(nbatri)
      dimension sx(jnmax),sy(jnmax),sz(jnmax)
      dimension ijx(225),ijy(225),ijz(225),ijpos(225)
      dimension ix(84),iy(84),iz(84),jx(84),jy(84),jz(84)
c
      COMMON/COORD/CHGS(112),C(3,112),ict(112),ns(112),ks(112)
      COMMON/EXYZS/EX(7,7,13),EY(7,7,13),EZ(7,7,13)
      COMMON/FUNCS/NATOM,N3N,nbfao,nbatri
      COMMON/GAUSS/ZS(200),CS(200)
      COMMON/NSHEL/KNUC(200),KTYPE(200),KSTART(200),KLOC(200),
     1 KMIN(200),KMAX(200),KNG(200),NSHELL,NPRIM
      COMMON/SIGNS/IPQ(500),IPRNT
c
      DATA ITOL / 10 /
      DATA PI32 / 5.56832799683170D+00 /
      DATA DEBYE / 2.541765480D+00 /
      DATA ZERO0,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /
      DATA NQLIM / 7 /
C INDEX DATA
C FIRST ROW CONTAINS S,P, AND D INDICES
C SECOND ROW CONTAINS F
C THIRD ROW CONTAINS G
C FOURTH AND FIFTH ROWS CONTAIN  H
C SIXTH AND SEVENTH ROWS CONTAIN  I
C
      DATA IX / 1, 8, 1, 1,15, 1, 1, 8, 8, 1,
     *         22, 1, 1,15,15, 8, 1, 8, 1, 8,
     *         29, 1, 1,22,22, 8, 1, 8, 1,15,15, 1,15, 8, 8,
     *         36, 1, 1,29,29, 8, 1, 8, 1,22,22,15, 1,15, 1,22, 8, 8,15,
     *          8,15,
     *         43, 1, 1,36,36, 8, 1, 8, 1,29,29,15, 1,15, 1,22,22, 1,29,
     *          8, 8,22,22,15, 8,15, 8,15/
      DATA JX / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0,
     *          3, 0, 0, 2, 2, 1, 0, 1, 0, 1,
     *          4, 0, 0, 3, 3, 1, 0, 1, 0, 2, 2, 0, 2, 1, 1,
     *          5, 0, 0, 4, 4, 1, 0, 1, 0, 3, 3, 2, 0, 2, 0, 3, 1, 1, 2,
     *          1, 2,
     *          6, 0, 0, 5, 5, 1, 0, 1, 0, 4, 4, 2, 0, 2, 0, 3, 3, 0, 5,
     *          1, 1, 3, 3, 2, 1, 2, 1, 2/
      DATA IY / 1, 1, 8, 1, 1,15, 1, 8, 1, 8,
     *          1,22, 1, 8, 1,15,15, 1, 8, 8,
     *          1,29, 1, 8, 1,22,22, 1, 8,15, 1,15, 8,15, 8,
     *          1,36, 1, 8, 1,29,29, 1, 8,15, 1,22,22, 1,15, 8,22, 8,15,
     *         15, 8,
     *          1,43, 1, 8, 1,36,36, 1, 8,15, 1,29,29, 1,15,22, 1,22, 8,
     *         29, 8,15, 8,22,22, 8,15,15/
      DATA JY / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1,
     *          0, 3, 0, 1, 0, 2, 2, 0, 1, 1,
     *          0, 4, 0, 1, 0, 3, 3, 0, 1, 2, 0, 2, 1, 2, 1,
     *          0, 5, 0, 1, 0, 4, 4, 0, 1, 2, 0, 3, 3, 0, 2, 1, 3, 1, 2,
     *          2, 1,
     *          0, 6, 0, 1, 0, 5, 5, 0, 1, 2, 0, 4, 4, 0, 2, 3, 0, 3, 1,
     *          5, 1, 2, 1, 3, 3, 1, 2, 2/
      DATA IZ / 1, 1, 1, 8, 1, 1,15, 1, 8, 8,
     *          1, 1,22, 1, 8, 1, 8,15,15, 8,
     *          1, 1,29, 1, 8, 1, 8,22,22, 1,15,15, 8, 8,15,
     *          1, 1,36, 1, 8, 1, 8,29,29, 1,15, 1,15,22,22, 8, 8,22, 8,
     *         15,15,
     *          1, 1,43, 1, 8, 1, 8,36,36, 1,15, 1,15,29,29, 1,22,22, 8,
     *          8,29, 8,15, 8,15,22,22,15/
      DATA JZ / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1,
     *          0, 0, 3, 0, 1, 0, 1, 2, 2, 1,
     *          0, 0, 4, 0, 1, 0, 1, 3, 3, 0, 2, 2, 1, 1, 2,
     *          0, 0, 5, 0, 1, 0, 1, 4, 4, 0, 2, 0, 2, 3, 3, 1, 1, 3, 1,
     *          2, 2,
     *          0, 0, 6, 0, 1, 0, 1, 5, 5, 0, 2, 0, 2, 4, 4, 0, 3, 3, 1,
     *          1, 4, 1, 2, 1, 2, 3, 3, 2/
c
      TOL=ITOL*2.30258D+00
      call zero(s,nbatri)
C
C:::::::::::::::::::::::::::::
C:::::----- I SHELL -----:::::
C:::::::::::::::::::::::::::::
C
      DO 6000 ISH=1,NSHELL
C
      IATOM=KNUC(ISH)
      AX=C(1,IATOM)
      AY=C(2,IATOM)
      AZ=C(3,IATOM)
      IGMIN=KSTART(ISH)
      IGMAX=IGMIN+KNG(ISH)-1
      LIT=KTYPE(ISH)
      MINI=KMIN(ISH)
      MAXI=KMAX(ISH)
      LOCI=KLOC(ISH)
C
C:::::::::::::::::::::::::::::
C:::::----- J SHELL -----:::::
C:::::::::::::::::::::::::::::
C
      DO 5000 JSH=1,ISH
C
      JATOM=KNUC(JSH)
      BX=C(1,JATOM)
      BY=C(2,JATOM)
      BZ=C(3,JATOM)
      JGMIN=KSTART(JSH)
      JGMAX=JGMIN+KNG(JSH)-1
      LJT=KTYPE(JSH)
      MINJ=KMIN(JSH)
      MAXJ=KMAX(JSH)
      LOCJ=KLOC(JSH)
      IANDJ=ISH.EQ.JSH
C
C     ----- WORK OUT INDICES FOR COMBINING 2D INTEGRALS -----
C
      MIJ=0
      JMAX=MAXJ
      II=-1
      DO 104 I=MINI,MAXI
      II=II+1
      JJ=-1
      NX=IX(I)
      NY=IY(I)
      NZ=IZ(I)
      IF(IANDJ) JMAX=I
      DO 103 J=MINJ,JMAX
      MIJ=MIJ+1
      JJ=JJ+1
      IJX(MIJ)=NX+JX(J)
      IJY(MIJ)=NY+JY(J)
      IJZ(MIJ)=NZ+JZ(J)
      IDXI=LOCI+II
      IDXJ=LOCJ+JJ
      IDXS=IPQ(MAX0(IDXI,IDXJ))+MIN0(IDXI,IDXJ)
      IJPOS(MIJ)=IDXS
  103 CONTINUE
  104 CONTINUE
C
      XAB=AX-BX
      YAB=AY-BY
      ZAB=AZ-BZ
      RAB=XAB*XAB+YAB*YAB+ZAB*ZAB
      MAXIJ=MAX0(LIT+2,LJT+2)
C
      DO 4000 IG=IGMIN,IGMAX
      CSI=CS(IG)*PI32
      AI=ZS(IG)
      AXI=AI*AX
      AYI=AI*AY
      AZI=AI*AZ
C
      IF(IANDJ) JGMAX=IG
      DO 3000 JG=JGMIN,JGMAX
      CSJ=CS(JG)
      BJ=ZS(JG)
      BXJ=BJ*BX
      BYJ=BJ*BY
      BZJ=BJ*BZ
C
      PP=AI+BJ
      TPP=ONE/PP
      TEMP=(AI*BJ*RAB)*TPP
C
C     ----- IGNORE INTEGRAL IF PRE-EXPONENTIAL TOO SMALL -----
C
      IF(TEMP.GT.TOL) GO TO 3000
C
C     ----- DENSITY FACTOR -----
C
      TP=DSQRT(TPP)
      DIJ=CSI*CSJ*DEXP(-TEMP)*TPP
      IF(IANDJ.AND.IG.NE.JG) DIJ=DIJ+DIJ
      PFAC=DIJ*TP
      PTWO=TPP*HALF
C
C     FORM SOME PRODUCTS NEEDED IN INTEGRAL EVALUATION
C
      PX=(AXI+BXJ)*TPP
      PY=(AYI+BYJ)*TPP
      PZ=(AZI+BZJ)*TPP
      PAX=PX-AX
      PAY=PY-AY
      PAZ=PZ-AZ
      PBX=PX-BX
      PBY=PY-BY
      PBZ=PZ-BZ
C
      CALL ECAL(PAX,PAY,PAZ,PBX,PBY,PBZ,PTWO,MAXIJ)
C
      A12=AI+AI
      B12=BJ+BJ
C
      IN=-NQLIM
c
      DO 106 NI=1,LIT
      IN=IN+NQLIM
      A11=DBLE(NI-1)
      DO 105 NJ=1,LJT
      JN=IN+NJ
      B11=DBLE(NJ-1)
c
c   overlap integrals
      sx(jn)=ex(ni,nj,1)
      sy(jn)=ey(ni,nj,1)
      sz(jn)=ez(ni,nj,1)
c
 105  continue
 106  continue
c
      do 107 i=1,mij
      ij=ijpos(i)
      nx=ijx(i)
      ny=ijy(i)
      nz=ijz(i)
c
c   form s(i,j) ao overlap integrals
c
      s(ij)=s(ij)+ sx(nx)*sy(ny)*sz(nz)*pfac
c
 107  continue
c
c    finish loops over shells
 3000 continue
 4000 continue
 5000 continue
 6000 continue
c
c   print out overlap integrals
c
      if(iprnt.gt.0)then
      write(6,500)
      call print(s,nbatri,nbfao,6)
      endif
c
 500  format(/1x,'overlap matrix')
c
      end
