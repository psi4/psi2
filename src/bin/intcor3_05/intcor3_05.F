      subroutine fentry(cc,ia,maxcor)
c**********************************************************************
C     PROGRAM INTCOR3_05                                              *
c     This program converts gradientrs and Hessian                    *
c       from Cartesian to internal coordinate system                  *
c**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
c**********************************************************************
C   Last updated on August 02, 2005 by Dr. Yukio Yamaguchi            *
c**********************************************************************
C   Last updated on April 30, 2003 by Dr. Yukio Yamaguchi             *
c**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension cc(maxcor),ia(maxcor*2)
      character*1 maxis, c2axis
      character*2 mplane
      character*3 grad, hess
      character*6 prog, dipole
      character*8 moltyp, unit
      character*12 wfn
      character*80 title
      common/com101/natoms
      common/com102/iprnt
      common/com110/escf
      common/com201/chg(10),w(10)
      common/com202/x(10),y(10),z(10)
      common/com205/r(10,10)
      common/com206/kl(2),lm(2),klm(3)
      common/com208/maxis,c2axis,mplane
      common/com209/ilin,isym
      common/com210/grad,hess
      common/vib101/n3n,nints
      common/vib102/nstr,nbnd,nlib
      common/vib103/engx(10),engy(10),engz(10)
      common/vib104/rint(10),sint(10),engr(10),engs(10)
      common/vib201/kr(10),lr(10)
      common/vib202/ka(10),la(10),ma(10),agla(10)
      common/vib203/kb(10),lb(10),mb(10),aglb(10)
c     data bohr / 0.529177249d+00 /
      data bohr / 0.52917706d+00 /
c==================================================================
c  For psi 2.0 input
#include <error.h>
      integer frdc, frdi, frdf, fwkclr, fwkadd, frdboo, frdpos
      integer fcnt, frdia1
      integer errcod
c==================================================================
    1 FORMAT(//,2X,' THE INTERNAL ENERGY GRADIENT AND HESSIAN PROGRAM'/
     1          2x,'              (Version of 08/02/2005)'/)
    2 FORMAT(A80)                                                       FDG00200
    3 FORMAT(3x,a80)
    4 FORMAT(I5,f20.10)
    5 FORMAT(//,2X,' PARAMETERS'/
     1 2X,' PROGRAM   = ',A6/
     2 2X,' WFN       = ',A12/
     3 2X,' UNIT      = ',A8/
     4 2X,' MOL_TYPE  = ',A8/
     5 2X,' MOL_AXIS  = ',A1/
     6 2X,' C2_AXIS   = ',A1/
     7 2X,' MOL_PLANE = ',A2/
     8 2X,' GRADIENT  = ',A3/
     9 2X,' HESSIAN   = ',A3/
     a 2X,' IPRNT     = ',I5/
     b 2X,' NATOMS    = ',I5/
     c 2X,' N3N       = ',I5/
     d 2X,' NSTR      = ',I5/
     e 2X,' NBND      = ',I5/
     f 2X,' NLIB      = ',I5/
     g 2X,' NINTS     = ',I5/
     h 2X,' ILIN      = ',I5/
     i 2X,' ISYM      = ',I5/)
    6 FORMAT(//,2X,' THE REFERENCE GEOMETRY (IN BOHR)'/
     1 5X,4H NO.,11x,7H CHARGE,13X,2H X,18X,2H Y,18X,2H Z/)               FDG00380
    7 FORMAT(4F20.10)  
    8 FORMAT(2X,I7,5X,4F20.10)                                          FDG00390
    9 FORMAT(//,2X,' CARTESIAN ENERGY GRADIENTS (IN HARTREE/BOHR)'/
     1 5X,4H NO.,11x,7H CHARGE,13X,2H X,18X,2H Y,18X,2H Z/)               FDG00380
C                                                                       FDG00570
      ITAP11=11                                                         FDG00600
      ITAP12=12                                                         FDG00600
      ITAP15=15                                                         FDG00600
      ITAP16=16                                                         FDG00600
C                                                                       FDG00620
      call psinit('APPEND')
      CALL TSTART(6)                                                    FDG00550
      call ffile(itap11,'file11',0)
      call ffile(itap12,'file12',0)
      call ffile(itap15,'file15',0)
      call ffile(itap16,'hess',0)
      WRITE(6,1)                                                        FDG00630
      write(6,*) '  maxcor    = ',maxcor
C
C=============================================================
C   MOL_TYPE
C=============================================================
C    Linear asymmetric molecules
C     A(1) atom is located in the 'origin'
C     B(2) atom is located along the 'specified axis (maxis)'
C     C(3) atom is located along the 'specified axis (maxis)'
C
C    Linear symmetric molecules
C     A(1) atom is located in the origin bisecting 'C2' axis
C     B(2) atom is located along the 'specified axis (maxis)'
C     B(3) atom is located along the 'specified axis (maxis)'
C
C    Bent asymmetric molecules
C     A(1) atom is located in the 'origin'
C     B(2) atom is located along the `specified axis (maxis)'
C     C(3) atom is located in the `specified plane (mplane)'
C
C    Bent symmetric molecules
C     A(1) atom is located in the 'origin' bisecting 'C2' axis
C     B(2) atom is located in the 'specified plane (mplane)'
C     B(3) atom is located in the 'specified plane (mplane)'

C=============================================================
C
C     A(1)-B(2)-C(3)
C    'abc_lin ' : R1 = r(ab) = r(1,2), R2 = r(bc) = r(2,3)
C    'abc_bent' : R1 = r(ab) = r(1,2), R2 = r(bc) = r(2,3)
C               :    R3 = <abc = <(1,2,3)
C
C     C(3)-A(1)-B(2)
C    'cab_lin ' : R1 = r(ab) = r(1,2), R2 = r(ac) = r(3,1)
C    'cab_bent' : R1 = r(ab) = r(1,2), R2 = r(ac) = r(3,1)
C               :    R3 = <cab = <(3,1,2)
C
C     B(2)-A(1)-C(3)
C    'bac_lin ' : R1 = r(ab) = r(2,1), R2 = r(ac) = r(1,3)
C    'bac_bent' : R1 = r(ab) = r(2,1), R2 = r(ac) = r(1,3)
C               :    R3 = <bac = <(2,1,3)
C
C     B(2)-A(1)-B(3)
C    'bab_lin ' : R1 = r(ba) = r(2,1), R2 = r(ac) = r(1,3)
C    'bab_bent' : R1 = r(ab) = r(2,1), R2 = r(ac) = r(1,3)
C               :    R3 = <bab = <(2,1,3)
C
C=============================================================
C
      errcod    = fwkclr()
      errcod    = fwkadd(5,'INTCOR3')
      errcod    = fwkadd(5,'DEFAULT')
      errcod    = frdi('PRINT',iprnt)
      unit      = 'ATOMIC  '
      errcod    = frdc('UNIT',unit)
      prog      = 'PSI  '
      errcod    = frdc('PROGRAM',prog)
      wfn       = 'SCF  '
      errcod    = frdc('WFN',wfn)
      grad      = 'YES'
      errcod    = frdc('GRADIENT',grad)
      hess      = 'YES'
      errcod    = frdc('HESSIAN',hess)
C
      errcod = fwkadd(5,'MOLECULE')
      errcod = frdc('MOL_TYPE',moltyp)
      errcod = frdc('MOL_AXIS',maxis)
      c2axis = 'N'
      errcod = frdc('C2_AXIS',c2axis)
      errcod = frdc('MOL_PLANE',mplane)
      call izero(kl,2)
      call izero(lm,2)
      call izero(klm,3)
      errcod = frdia1('RKL',kl,2)
      errcod = frdia1('RLM',lm,2)
      errcod = frdia1('AKLM',klm,3)
      errcod = frdia1('BOND_1',kl,2)
      errcod = frdia1('BOND_2',lm,2)
      errcod = frdia1('ANGLE',klm,3)
C
      REWIND ITAP11  
      READ(ITAP11,2) title
      write(6,3) title
      READ(ITAP11,4) NATOMS,escf
      n3n = natoms*3
C
      kr(1) = kl(1)
      lr(1) = kl(2)
      kr(2) = lm(1)
      lr(2) = lm(2)
      ka(1) = klm(1)
      la(1) = klm(2)
      ma(1) = klm(3)
      kb(1) = klm(1)
      lb(1) = klm(2)
      mb(1) = klm(3)
C
      nstr  = 0
      nbnd  = 0
      nlib  = 0
      ilin  = 0
      isym  = 0
      if(natoms.eq.2) then
        nstr = 1
      end if
      if(natoms.eq.3) then
        nstr = 2
        nbnd = 1
        nlib = 0
        ilin = 0
      end if
      linear = index(moltyp(1:),'LIN')
      if(linear.ne.0) then
        nstr = 2
        nbnd = 0
        nlib = 2
        ilin = 1
      end if
      isymty = index(moltyp(1:),'BAB')
      if(isymty.ne.0) then
        isym = 1
      end if
      nints = nstr + nbnd + nlib
C
      WRITE(6,5) PROG,WFN,UNIT,MOLTYP,MAXIS,C2AXIS,MPLANE,GRAD,HESS,
     *           iprnt,natoms,n3n,nstr,nbnd,nlib,nints,ilin,isym
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::
C   Read in geometries from `file11.dat' file
C::::::::::::::::::::::::::::::::::::::::::::::::::::::
      REWIND ITAP11                                                     FDG00800
      READ(ITAP11,2) TITLE                                              FDG00810
      READ(ITAP11,4) NATOMS,escf
      WRITE(6,6)    
      DO 101 I=1,NATOMS                                                 FDG00830
      READ(ITAP11,7) CHG(I),X(I),Y(I),Z(I)
      WRITE(6,8) I,CHG(I),X(I),Y(I),Z(I)
  101 CONTINUE                                                          FDG00850
      call dist
      WRITE(6,9)
      DO 102 I=1,NATOMS                                                 FDG00830
      READ(ITAP11,7) dummy,engx(i),engy(i),engz(i)
      WRITE(6,8) i,chg(i),engx(i),engy(i),engz(i)
  102 CONTINUE                                                          FDG00850
C                                                                       FDG01130
C===========================================
C  Dynamic allocation
C===========================================
      n3nsq = n3n*n3n
      natsq = natoms*natoms
      ic1   = 1
      ic2   = ic1 + n3nsq
      ic3   = ic2 + n3nsq
      ic4   = ic3 + n3nsq
      ic5   = ic4 + n3nsq
      ic6   = ic5 + n3nsq
      ic7   = ic6 + n3n
      ia7   = ic7 + ic7 - 1
      ic8   = ic7 + n3n
      ic9   = ic8 + natsq
      ic10  = ic9 + natsq
      icmax = ic10 + natsq
      if(icmax.gt.maxcor) go to 400
C===========================================
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C   Form the B matrix
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C...............B       BIN     AA      BB      FX      FF......
      call bmat(cc(ic1),cc(ic2),CC(ic3),cc(ic4),cc(ic5),cc(ic5),
C...............PIVOT   INDEXX  EX      EY      EZ.......
     *          cc(ic6),ia(ia7),cc(ic8),cc(ic9),cc(ic10))
C
      go to 500
C
  400 continue
      write(6,*) '  '
      write(6,*) '  ***************************'
      write(6,*) '  ICMAX exceeds MAXCOR !!!   '
      write(6,*) '  ICMAX   = ',icmax
      write(6,*) '  MAXCOR  = ',maxcor
      write(6,*) '  ---> STOP !!!              '
      write(6,*) '  ***************************'
      write(6,*) '  '
C
  500 continue
C
      CALL TSTOP(6)                                                     FDG01730
      STOP                                                              FDG01750
      END                                                               FDG01760
      subroutine dist
      IMPLICIT REAL*8 (A-H,O-Z)                                         NOR05390
      common/com101/natoms
      common/com202/x(10),y(10),z(10)
      common/com205/r(10,10)
    1 FORMAT(//,2X,' INTERATOMIC DISTANCE MATRIX'/)                     NOR05450
C                                                                       NOR05460
      call zero(r,100)
      DO 101 I=1,NATOMS
      DO 101 J=1,NATOMS
      IF(I.EQ.J) GO TO 101                                              NOR05520
      XD=X(I)-X(J)                                                      NOR05530
      YD=Y(I)-Y(J)                                                      NOR05540
      ZD=Z(I)-Z(J)                                                      NOR05550
      valu=DSQRT(XD*XD+YD*YD+ZD*ZD)   
      r(i,j)=valu
  101 CONTINUE                                                          NOR05570
      WRITE(6,1)                                                        NOR05590
      CALL MATOUT(R,10,10,NATOMS,NATOMS,6)
c
      return
      end
      SUBROUTINE BMAT(B,BIN,AA,BB,FX,FF,PIVOT,INDEXX,EX,EY,EZ) 
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT02220
      DIMENSION B(N3N,N3N),BIN(N3N,N3N),AA(N3N,N3N),BB(N3N,N3N)
      DIMENSION FX(N3N,N3N),FF(N3N,N3N)
      DIMENSION PIVOT(N3N),INDEXX(N3N*2)
      DIMENSION EX(NATOMS,NATOMS),EY(NATOMS,NATOMS),EZ(NATOMS,NATOMS)
      DIMENSION FS(4,4)
      character*3 grad, hess
      common/com101/natoms
      common/com102/iprnt
      common/com110/escf
      common/com201/chg(10),w(10)
      common/com202/x(10),y(10),z(10)
      common/com205/r(10,10)
      common/com209/ilin,isym
      common/com210/grad,hess
      COMMON/VIB101/N3N,NINTS
      COMMON/VIB102/NSTR,NBND,NLIB
      COMMON/VIB103/ENGX(10),ENGY(10),ENGZ(10)                          INT02330
      COMMON/VIB104/RINT(10),SINT(10),ENGR(10),ENGS(10)
      COMMON/VIB201/KR(10),LR(10)
      COMMON/VIB202/KA(10),LA(10),MA(10),AGLA(10)
      COMMON/VIB203/KB(10),LB(10),MB(10),AGLB(10) 
      DATA A00,ONE,TWO / 0.0D+00 , 1.0D+00 , 2.0D+00 /
      DATA PI / 3.1415926536D+00 /                                      INT02450
    1 FORMAT(//)                                                        INT02470
    2 FORMAT(2X,' PARAMETERS FOR THE CALCULATION'//                     INT02490
     * 2X,' IPRNT  = ',I5/      
     * 2X,' NSTR   = ',I5/                                              INT02500
     * 2X,' NBND   = ',I5/                                              INT02510
     * 2X,' NLIB   = ',I5/                                              INT02520
     * 2X,' NINTS  = ',I5/                                              INT02600
     * 2X,' NATOMS = ',I5/                                              INT02610
     * 2X,' N3N    = ',I5)                                              INT02620
    3 FORMAT(//,2X,' EX MATRIX'/)                                       INT02830
    4 FORMAT(//,2X,' EY MATRIX'/)                                       INT02840
    5 FORMAT(//,2X,' EZ MATRIX'/)                                       INT02850
    6 FORMAT(//,2X,21H INTERNAL COORDINATES,50X,11H COORDINATE/         INT02640
     1 73X,4H NO.)                                                      INT02650
    7 FORMAT(2X,' BOND LENGTH (IN BOHR)'/
     1 2X,20H   R   K-L    K    L/)                                     INT02670
    8 FORMAT(2X,F10.5,2I5,51X,2H (,I2,1H))                              INT02680
    9 FORMAT(/,2X,' BOND ANGLE (IN DEGREES)'/
     1 2X,25H ANG-KLM      K    L    M/)                                INT02710
   10 FORMAT(2X,F10.5,3I5,46X,2H (,I2,1H))                              INT02720
   11 FORMAT(/,2X,' LINEAR BOND ANGLE (IN DEGREES)'/ 
     1 2X,8H ANG-KLM,5X,2H K,3X,2H L,3X,2H M/)                          INT02740
   12 FORMAT(//,2X,' B MATRIX'/)                                        INT02860
   13 FORMAT(//,2X,' B INVERSE MATRIX'/)                                INT02870
   14 FORMAT(//,2X,' DETERM OF B MATRIX = ',D15.7)                      INT02880
   15 FORMAT(//,2X,' B * B(-1) = E'/)                                   INT02890
   16 FORMAT(//,2X,' BB MATRIX'/)
   17 FORMAT(//,2X,' INTERNAL ENERGY GRADIENTS '/                       INT02900
     *2x,' (BOND LENGTH IN HARTREE/BOHR AND BOND ANGLE IN HARTREE/RADIAN
     *)'/5X,4H NO.,10X,7H COORD.,13X,6H GRAD./)
   18 FORMAT(2X,I7,5X,2F20.10)                                          INT02920
   19 FORMAT(//,2X,' SYMMETRY INTERNAL ENERGY GRADIENTS '/    
     *2x,' (BOND LENGTH IN HARTREE/BOHR AND BOND ANGLE IN HARTREE/RADIAN
     *)'/5X,4H NO.,10X,7H COORD.,13X,6H GRAD./)
   20 FORMAT(2I5)                                                       INT02930
   21 FORMAT(3F20.10)                                                   INT02940
   22 FORMAT(//,2X,' FX MATRIX (IN HARTREE/BOHR+2)'/)
   23 FORMAT(//,2X,' BIN * FX  MATRIX'/)                                INT02960
   24 FORMAT(//,2X,' F MATRIX'/                                         INT02970
     * 2x,' [STR-STR IN HARTREE/(BOHR+2), STR-BEND IN HARTREE/(BOHR.RADI
     *AN), and BEND-BEND IN HARTREE/(RADIAN+2)]'/)
   25 FORMAT(//,2X,' FS MATRIX'/   
     * 2x,' [STR-STR IN HARTREE/(BOHR+2), STR-BEND IN HARTREE/(BOHR.RADI
     *AN), and BEND-BEND IN HARTREE/(RADIAN+2)]'/)
   26 FORMAT(2I5,F20.10) 
   27 FORMAT(4F20.10)                                                   INT02990
   28 FORMAT(2F20.10)                                                   INT03000
   29 FORMAT(3F20.10)                                                   INT03010
C                       
      AR = PI/(180.0D+00)
      rsq2 = one / (dsqrt(two))
C
      if(iprnt.ge.3) then
        WRITE(6,1) 
        WRITE(6,2) IPRNT,NSTR,NBND,NLIB,NINTS,NATOMS,N3N
      end if
      n3nsq = n3n*n3n
      natsq = natoms*natoms
C                                                                       INT03110
C   THE CALCULATION OF B VECTORS                                        INT03120
C   THE CALCULATION OF UNIT VECTORS                                     INT03130
      call zero(ex,natsq)
      call zero(ey,natsq)
      call zero(ez,natsq)
      call zero(b,n3nsq)
c
      DO 101 I=1,natoms
      DO 101 J=1,natoms
      if(i.eq.j) go to 101
      EX(I,J)=-(X(I)-X(J))/R(I,J) 
      EY(I,J)=-(Y(I)-Y(J))/R(I,J)
      EZ(I,J)=-(Z(I)-Z(J))/R(I,J) 
  101 CONTINUE                                                          INT03320
C
C   PRINT OUT E VECTORS                                                 INT04650
      if(iprnt.ge.3) then
        WRITE(6,3) 
        CALL MATOUT(EX,NATOMS,NATOMS,NATOMS,NATOMS,6)    
        WRITE(6,4)                                  
        CALL MATOUT(EY,NATOMS,NATOMS,NATOMS,NATOMS,6)  
        WRITE(6,5)                                
        CALL MATOUT(EZ,NATOMS,NATOMS,NATOMS,NATOMS,6) 
      end if
C                                                                       INT04730
C  FOR STRETCHINGS
      WRITE(6,6)                                                        INT03350
      WRITE(6,7)                                                        INT03360
      DO 102 I=1,NSTR
      K=KR(I)                                                           INT03380
      L=LR(I)                                                           INT03390
      II=I                                                              INT03400
      RINT(II)=R(K,L)
      WRITE(6,8) R(K,L),KR(I),LR(I),II  
  102 CONTINUE                                                          INT03440
      CALL STRNG(B,EX,EY,EZ)   
      IF(NBND.EQ.0) GO TO 201                                           INT03460
C
C  FOR BENDINGS
      WRITE(6,9)                                                        INT03490
      DO 103 I=1,NBND                                                   INT03500
      K=KA(I)                                                           INT03510
      L=LA(I)                                                           INT03520
      M=MA(I)                                                           INT03530
      II=I+NSTR                                                         INT03540
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::INT03550
C   BOND ANGLE CALCULATION (WILSON-DECIOUS-CROSS  P.57, EQ. 8)          INT03560
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::INT03570
      AKLM=EX(L,K)*EX(L,M)+EY(L,K)*EY(L,M)+EZ(L,K)*EZ(L,M)              INT03580
      AGLA(I)=DACOS(AKLM)                                               INT03590
      ANG=AGLA(I)/AR                                                    INT03600
      RINT(II)=AGLA(I)                                                  INT03610
      WRITE(6,10) ANG,KA(I),LA(I),MA(I),II                              INT03620
  103 CONTINUE                                                          INT03630
      CALL BENDNG(B,EX,EY,EZ) 
  201 continue
      IF(NLIB.EQ.0) GO TO 202                                           INT03650
C
C  FOR LINEAR BENDINGS
      NLHF=NLIB/2                                                       INT03670
      WRITE(6,11)                                                       INT03690
      DO 104 I=1,NLHF                                                   INT03700
      I1=2*I-1+NSTR+NBND
      I2=I1+1                                                           INT03720
      K=KB(I)                                                           INT03730
      L=LB(I)                                                           INT03740
      M=MB(I)                                                           INT03750
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::INT03760
C   BOND ANGLE CALCULATION (WILSON-DECIOUS-CROSS  P.57, EQ. 8)          INT03770
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::INT03780
      AKLM=EX(L,K)*EX(L,M)+EY(L,K)*EY(L,M)+EZ(L,K)*EZ(L,M)              INT03790
      IF(DABS(AKLM).GT.ONE) AKLM=DSIGN(ONE,AKLM)                        INT03800
      AGLB(I)=DACOS(AKLM)                                               INT03810
      ANG=AGLB(I)/AR                                                    INT03820
      RINT(I1)=AGLB(I)                                                  INT03830
      RINT(I2)=AGLB(I)                                                  INT03840
      WRITE(6,10) ANG,KB(I),LB(I),MB(I),I1                              INT03850
      WRITE(6,10) ANG,KB(I),LB(I),MB(I),I2                              INT03860
  104 CONTINUE                                                          INT03870
      CALL LINBND(B,EX,EY,EZ)   
  202 continue
C                                                                       INT04640
      if(iprnt.ge.3) then
        WRITE(6,12)    
        CALL MATOUT(B,N3N,N3N,NINTS,N3N,6)   
      end if
C
C   THE CALCULATION OF INVERSE B MATRIX                                 INT04770
      call zero(bb,n3nsq)
      DO 107 I=1,NINTS
      DO 107 J=1,NINTS
      valu=a00
      DO 106 K=1,N3N                                                    INT04810
      valu=valu+B(I,K)*B(J,K)   
  106 CONTINUE                                                          INT04830
      BB(i,j)=valu
  107 CONTINUE                                                          INT04830
      call zero(bin,n3nsq)
      DO 108 I=1,NINTS
      DO 108 J=1,NINTS
      BIN(I,J)=BB(I,J)                                                  INT04860
  108 continue
C
      CALL MATINV(BIN,AA,PIVOT,INDEXX,NINTS,NINTS,DETERM,N3N)
      if(iprnt.ge.3) then
        WRITE(6,13) 
        CALL MATOUT(BIN,N3N,N3N,NINTS,NINTS,6)   
      end if
      WRITE(6,14) DETERM                                                INT04910
C
      call zero(aa,n3nsq)
      DO 110 I=1,NINTS
      DO 110 J=1,NINTS
      valu=A00 
      DO 109 K=1,NINTS
      valu=valu+BIN(I,K)*BB(K,J)   
  109 CONTINUE                                                          INT04970
      AA(i,j)=valu
  110 CONTINUE                                                          INT04970
      if(iprnt.ge.3) then
        WRITE(6,15)    
        CALL MATOUT(AA,N3N,N3N,NINTS,NINTS,6)   
      end if
C
      call zero(bb,n3nsq)
      DO 112 I=1,NINTS
      DO 112 J=1,N3N                                                    INT05020
      valu=A00
      DO 111 K=1,NINTS
      valu=valu+BIN(I,K)*B(K,J)    
  111 CONTINUE                                                          INT05060
      BB(I,J)=valu
  112 CONTINUE                                                          INT05060
      if(iprnt.gt.3) then
        write(6,16)
        call matout(bb,n3n,n3n,nints,n3n,6)
      end if
C                                                                       INT05070
C***********************************************************            INT05080
C***TRANSFORM GRADIENTS FROM CARTESIAN TO INTERNAL SYSTEM***            INT05090
C***********************************************************            INT05100
      if(grad.ne.'YES') go to 205
      DO 116 I=1,NINTS                                                  INT05110
      valu=A00
      DO 115 J=1,NATOMS
      JX=3*J-2                                                          INT05140
      JY=3*J-1                                                          INT05150
      JZ=3*J                                                            INT05160
      valu=valu+ENGX(J)*BB(I,JX)+ENGY(J)*BB(I,JY)+ENGZ(J)*BB(I,JZ)
  115 CONTINUE                                                          INT05190
      ENGR(I)=valu
  116 CONTINUE                                                          INT05200
      WRITE(6,17)                                                       INT05210
      DO 117 I=1,NINTS                                                  INT05220
      WRITE(6,18) I,RINT(I),ENGR(I)                                     INT05230
  117 CONTINUE                                                          INT05240
C
C   Adjustment for symmetry coordinates
      if(isym.ne.0) then
        SINT(1) = rsq2 * ( RINT(1) + RINT(2) )
        SINT(2) = RINT(3)
        SINT(3) = rsq2 * ( RINT(1) - RINT(2) )
        ENGS(1) = rsq2 * ( ENGR(1) + ENGR(2) )
        ENGS(2) = ENGR(3)
        ENGS(3) = rsq2 * ( ENGR(1) - ENGR(2) )
        ENGS(4) = ENGS(2)
        WRITE(6,19)        
        DO 118 I=1,NINTS  
        WRITE(6,18) I,SINT(I),ENGS(I)  
  118   CONTINUE
      end if
C
C***********************************************************************INT05260
C***TRANSFORM FORCE CONSTANT MATRIX FROM CARTESIAN TO INTERNAL SYSTEM***INT05270
C***********************************************************************INT05280
  205 continue
      if(hess.ne.'YES') go to 210
C   READ IN FORCE CONSTANT MATRIX                                       INT05300
      call zero(fx,n3nsq)
      itap15 = 15
      REWIND ITAP15                                                     INT05310
      READ(ITAP15,20) NATOMX,NNC 
      READ(ITAP15,21) ((FX(I,J),J=1,N3N),I=1,N3N)                       INT05330
      REWIND ITAP15                                                     INT05310
      if(iprnt.gt.3) then
       WRITE(6,22)                        
       CALL MATOUT(FX,N3N,N3N,N3N,N3N,6)     
      end if
c
c   TRANSFORM FX INTO F                                                 INT05370
      call zero(aa,n3nsq)
      DO 121 I=1,NINTS 
      DO 121 J=1,N3N                                                    INT05390
      valu=A00
      DO 120 K=1,N3N                                                    INT05410
      valu=valu+BB(I,K)*FX(K,J) 
  120 CONTINUE                                                          INT05430
      AA(I,J)=valu
  121 CONTINUE                                                          INT05440
      if(iprnt.gt.3) then
        WRITE(6,23)     
        CALL MATOUT(AA,N3N,N3N,NINTS,N3N,6)    
      end if
C
      call zero(ff,n3nsq)
      DO 123 I=1,NINTS
      DO 123 J=1,NINTS
      valu=A00
      DO 122 K=1,N3N                                                    INT05510
      valu=valu+AA(I,K)*BB(J,K)  
  122 CONTINUE                                                          INT05530
      FF(I,j)=valu
  123 CONTINUE                                                          INT05540
      WRITE(6,24)    
      CALL MATOUT(FF,N3N,N3N,NINTS,NINTS,6)    
      write(6,*) ' '
      if(isym.ne.0) then
        CALL ZERO(FS,16)
        FS(1,1) = FF(1,1) + FF(1,2)
        FS(1,2) = FF(1,3) * DSQRT(TWO)
        FS(2,1) = FS(1,2)
        FS(2,2) = FF(3,3)
        FS(3,3) = FF(1,1) - FF(1,2)
        FS(4,4) = FS(2,2)
        WRITE(6,25)    
        CALL MATOUT(FS,4,4,NINTS,NINTS,6)    
      end if
C
C   STORE GRADIENT INFORMATION IN TAPE12
  210 CONTINUE 
      if(grad.ne.'YES') go to 211
      itap12 = 12
      REWIND ITAP12                                                     INT05620
      WRITE(ITAP12,26) NATOMS,NINTS,ESCF
      DO 125 I=1,NATOMS
      WRITE(ITAP12,27) CHG(I),X(I),Y(I),Z(I)                            INT05650
  125 CONTINUE                                                          INT05660
      DO 126 I=1,NATOMS                                                 INT05670
      WRITE(ITAP12,27) CHG(I),ENGX(I),ENGY(I),ENGZ(I)                   INT05680
  126 CONTINUE                                                          INT05690
      if(isym.eq.0) then
        DO 127 I=1,NINTS     
        WRITE(ITAP12,28) RINT(I),ENGR(I) 
  127   CONTINUE                        
      end if
      if(isym.ne.0) then
        DO 128 I=1,NINTS
        WRITE(ITAP12,28) SINT(I),ENGS(I)
  128   CONTINUE
      end if
      REWIND ITAP12                                                     INT05760
C
C   STORE HESSIAN INFORMATION IN TAPE16 (HESS.DAT)
  211 continue
      if(hess.eq.'YES') then
        itap16 = 16
        rewind itap16
        if(isym.eq.0) then
          WRITE(ITAP16,29) ((FF(I,J),J=1,NINTS),I=1,NINTS)   
        end if
        if(isym.ne.0) then
          WRITE(ITAP16,29) ((FS(I,J),J=1,NINTS),I=1,NINTS)   
        end if
        rewind itap16
      end if
C                                                                       INT05770
      RETURN                                                            INT05780
      END                                                               INT05790
      SUBROUTINE STRNG(B,EX,EY,EZ)     
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT05810
      DIMENSION B(N3N,N3N)
      DIMENSION EX(NATOMS,NATOMS),EY(NATOMS,NATOMS),EZ(NATOMS,NATOMS)
      common/com101/natoms
      COMMON/VIB101/N3N,NINTS
      COMMON/VIB102/NSTR,NBND,NLIB
      COMMON/VIB201/KR(10),LR(10)
C                                                                       INT05870
C   THE CALCULATION OF S VECTORS FOR STRETCHINGS                        INT05880
      DO 101 I=1,NSTR                                                   INT05890
      K=KR(I)                                                           INT05900
      L=LR(I)                                                           INT05910
      KX=3*K-2                                                          INT05920
      KY=3*K-1                                                          INT05930
      KZ=3*K                                                            INT05940
      LX=3*L-2                                                          INT05950
      LY=3*L-1                                                          INT05960
      LZ=3*L                                                            INT05970
      B(I,KX)=EX(L,K)                                                   INT05980
      B(I,KY)=EY(L,K)                                                   INT05990
      B(I,KZ)=EZ(L,K)                                                   INT06000
      B(I,LX)=EX(K,L)                                                   INT06010
      B(I,LY)=EY(K,L)                                                   INT06020
      B(I,LZ)=EZ(K,L)                                                   INT06030
  101 continue
c
      RETURN                                                            INT06040
      END                                                               INT06050
      SUBROUTINE BENDNG(B,EX,EY,EZ)
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT06070
      DIMENSION B(N3N,N3N)
      DIMENSION EX(NATOMS,NATOMS),EY(NATOMS,NATOMS),EZ(NATOMS,NATOMS) 
      common/com101/natoms
      common/com205/r(10,10)
      COMMON/VIB101/N3N,NINTS
      COMMON/VIB102/NSTR,NBND,NLIB
      COMMON/VIB202/KA(10),LA(10),MA(10),AGLA(10)
C                                                                       INT06150
C   THE CALCULATION OF S VECTORS FOR ANGLE BENDINGS                     INT06160
      DO 101 I=1,NBND                                                   INT06170
      II=I+NSTR
      K=KA(I)                                                           INT06190
      L=LA(I)                                                           INT06200
      M=MA(I)                                                           INT06210
      KX=3*K-2                                                          INT06220
      KY=3*K-1                                                          INT06230
      KZ=3*K                                                            INT06240
      LX=3*L-2                                                          INT06250
      LY=3*L-1                                                          INT06260
      LZ=3*L                                                            INT06270
      MX=3*M-2                                                          INT06280
      MY=3*M-1                                                          INT06290
      MZ=3*M                                                            INT06300
      AG=AGLA(I)                                                        INT06310
      B(II,KX)=(DCOS(AG)*EX(L,K)-EX(L,M))/(R(L,K)*DSIN(AG))    
      B(II,KY)=(DCOS(AG)*EY(L,K)-EY(L,M))/(R(L,K)*DSIN(AG))   
      B(II,KZ)=(DCOS(AG)*EZ(L,K)-EZ(L,M))/(R(L,K)*DSIN(AG))  
      B(II,MX)=(DCOS(AG)*EX(L,M)-EX(L,K))/(R(L,M)*DSIN(AG)) 
      B(II,MY)=(DCOS(AG)*EY(L,M)-EY(L,K))/(R(L,M)*DSIN(AG))
      B(II,MZ)=(DCOS(AG)*EZ(L,M)-EZ(L,K))/(R(L,M)*DSIN(AG))
      B(II,LX)=-B(II,KX)-B(II,MX)                                       INT06400
      B(II,LY)=-B(II,KY)-B(II,MY)                                       INT06410
      B(II,LZ)=-B(II,KZ)-B(II,MZ)                                       INT06420
  101 continue
c
      RETURN                                                            INT06430
      END                                                               INT06440
      SUBROUTINE LINBND(B,EX,EY,EZ) 
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT06460
      DIMENSION B(N3N,N3N)
      DIMENSION EX(NATOMS,NATOMS),EY(NATOMS,NATOMS),EZ(NATOMS,NATOMS)
      character*1 maxis, c2axis
      character*2 mplane
      common/com101/natoms
      common/com205/r(10,10)
      common/com208/maxis,c2axis,mplane
      COMMON/VIB101/N3N,NINTS
      COMMON/VIB102/NSTR,NBND,NLIB
      COMMON/VIB203/KB(10),LB(10),MB(10),AGLB(10) 
      DATA ONE / 1.0D+00 /
      DATA PI / 3.1415926536D+00 /                                      INT06560
C                                                                       INT06570
      NLHF=NLIB/2                                                       INT06590
C   N IS A DUMMY ATOM                                                   INT06600
      DO 101 I=1,NLHF                                                   INT06610
      I1=2*I-1+NSTR+NBND 
      I2=I1+1                                                           INT06630
      K=KB(I)                                                           INT06640
      L=LB(I)                                                           INT06650
      M=MB(I)                                                           INT06660
      KX=3*K-2                                                          INT06670
      KY=3*K-1                                                          INT06680
      KZ=3*K                                                            INT06690
      LX=3*L-2                                                          INT06700
      LY=3*L-1                                                          INT06710
      LZ=3*L                                                            INT06720
      MX=3*M-2                                                          INT06730
      MY=3*M-1                                                          INT06740
      MZ=3*M                                                            INT06750
      rkl=one/R(K,L)
      rlm=one/R(L,M)
C
C======================================
      if(maxis.eq.'X') then
C  R(Y)
        B(I1,KZ) = rkl
        B(I1,MZ) = rlm
        B(I1,LZ) = -rkl-rlm
C  R(Z)
        B(I2,KY) = rkl
        B(I2,MY) = rlm
        B(I2,LY) = -rkl-rlm
      end if
C======================================
C
C======================================
      if(maxis.eq.'Y') then
C  R(Z)
        B(I1,KX) = rkl
        B(I1,MX) = rlm
        B(I1,LX) = -rkl-rlm
C  R(X)
        B(I2,KZ) = rkl
        B(I2,MZ) = rlm
        B(I2,LZ) = -rkl-rlm
      end if
C======================================
C
C======================================
      if(maxis.eq.'Z') then
C  R(Y)
        B(I1,KX) = rkl
        B(I1,MX) = rlm
        B(I1,LX) = -rkl-rlm
C  R(X)
        B(I2,KY) = rkl
        B(I2,MY) = rlm
        B(I2,LY) = -rkl-rlm
      end if
C======================================
C
  101 CONTINUE                                                          INT07070
c
      RETURN                                                            INT07080
      END                                                               INT07090
