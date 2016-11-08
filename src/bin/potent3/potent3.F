      subroutine fentry(cc,ia,maxcor)
c**********************************************************************
C     PROGRAM POTENT3                                                 *
c     This program will construct an analytic potential               *
c       for a triatomic molecule                                      *
c**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
c**********************************************************************
C   Last updated on April 18, 2003 by Dr. Yukio Yamaguchi             *
c**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (max_cal=5000)
      parameter (max_ijk=1000)
      dimension cc(maxcor),ia(maxcor*2)
      dimension RI0(max_cal),RJ0(max_cal),PK0(max_cal)
      dimension RI(max_cal),RJ(max_cal),PK(max_cal)
      character*2 wover
      character*2 kadjst
      character*3 legdre, spf
      character*6 padjst
      character*8 moltyp, unit
      character*12 wfn
      character*80 title, label
      common/com101/natoms
      common/com102/iprnt
      common/com103/itap10
      common/com105/apir,aglb,ang
      common/com201/chg(10),w(10)
      common/com202/xold(10),yold(10),zold(10)
      common/com203/xnew(10),ynew(10),znew(10)
      common/com204/xa(10),ya(10),za(10)
      common/com205/r(10,10)
      common/com206/kl(2),lm(2),klm(3)
      common/pot101/no_cals,ndim
      common/pot102/ipower,jpower,kpower,max_power
      common/pot103/del_max
      common/pot105/eref
      common/pot106/k_pref
      data a00,one,two / 0.0d+00 , 1.0d+00 , 2.0d+00 /
c     data bohr / 0.529177249d+00 /
      data bohr / 0.52917706d+00 /
      data pi / 3.1415926536d+00 /
      data wave / 219474.63d+00 /
c==================================================================
c  For psi 2.0 input
#include <error.h>
      integer frdc, frdi, frdf, fwkclr, fwkadd, frdboo, frdpos
      integer fcnt, frdia1
      integer errcod
c==================================================================
    1 FORMAT(//,2X,' THE ANALYTIC POTENTIAL PROGRAM FOR A TRIATOMIC MOLE
     1CULE'/
     2          2x,'            <To test an analytic potential>'/
     3          2x,'              (Version of 04/18/2003)'/)
    2 FORMAT(A80)                                                       FDG00200
    3 FORMAT(I5,F20.10,F15.10,2F20.10)                                  FDG00210
    4 FORMAT(I5)
    5 FORMAT(//,2X,' PARAMETERS'/
     1 2X,' WFN       = ',A12/
     2 2X,' UNIT      = ',A8/
     3 2X,' MOL_TYPE  = ',A8/
     4 2X,' IPRNT     = ',I6/
     5 2X,' NATOMS    = ',I6/
     6 2X,' N3N       = ',I6/
     7 2X,' LEGENDRE  = ',A3/
     8 2X,' SPF       = ',A3/
     9 2X,' NO_CALS   = ',I6/
     a 2X,' I_POWER   = ',I6/ 
     b 2X,' J_POWER   = ',I6/ 
     c 2X,' K_POWER   = ',I6/
     d 2X,' MAX_POWER = ',I6/ 
     e 2X,' NDIM      = ',I6/ 
     f 2X,' MAX_DELTA = ',I6/ 
     g 2X,' K_ADJUST  = ',A2/ 
     h 2X,' P_ADJUST  = ',A6/)
    6 FORMAT(//,2X,' THE REFERENCE GEOMETRY IN TAPE11'/
     1 39X,' Cartesian X in a.u.'/
     2 4X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z/)               FDG00380
    7 FORMAT(4F20.10)  
    8 FORMAT(3X,I5,1X,4F16.10)                                          FDG00390
    9 format(//,2x,' Equilibrium bond distances'/
     1          2x,' R_i(e)        = ',f20.10/
     2          2x,' R_j(e)        = ',f20.10/)
   10 FORMAT(//,2X,' THE REFERENCE GEOMETRY IN TAPE11'/
     1 39X,' Cartesian X in angstrom'/
     2 4X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z/)               FDG00380
   11 format(//,2x,' PERTURBED GEOMETRIES'/
     1 4X,4H NO.,6X,4H R_i,14X,4H R_j,14X,4H P_k/
     2 14X,7H (Bohr),11X,7H (Bohr),11X,9H (radian)/)
   12 format(i5,5x,3f15.10)
   13 format(2x,i5,3x,f15.10,3x,f15.10,3x,f15.10)
   14 format(//,2x,' PERTURBED GEOMETRIES'/
     1 4X,4H NO.,6X,7H R_i-Re,11X,7H R_j-Re,11X,
     2 7H P_k-Pe/
     3 14X,7H (Bohr),11X,7H (Bohr),11X,9H (radian)/)
   15 format(//,2x,' PERTURBED GEOMETRIES'/
     1 4X,4H NO.,6X,7H R_i-Re,11X,7H R_j-Re,11X,
     2 7H Pi-P_k/
     3 14X,7H (Bohr),11X,7H (Bohr),11X,9H (radian)/)
   16 format(//,2x,' PERTURBED GEOMETRIES'/
     1 4X,4H NO.,6X,15H (R_i-Re)/(R_i),3X,15H (R_j-Re)/(R_j),3X,
     2 7H P_k-Pe/
     3 14X,6H ( - ),12X,6H ( - ),12X,9H (radian)/)
   17 format(//,2x,' PERTURBED GEOMETRIES AND ENERGIES'/
     1 4X,4H NO.,6X,15H (R_i-Re)/(R_i),3X,15H (R_j-Re)/(R_j),3X,
     2 7H Pi-P_k/
     3 14X,6H ( - ),12X,6H ( - ),12X,9H (radian)/)
C                                                                       FDG00570
      INPUT=5                                                           FDG00580
      ITAPE6=6                                                          FDG00590
      ITAP63=63
      ITAP71=71
      ITAP11=11                                                         FDG00600
      ITAP30=30                                                         FDG00610
C                                                                       FDG00620
      call psinit('APPEND')
      CALL TSTART(6)                                                    FDG00550
      call ffile(itap11,'file11',0)
      call ffile(itap63,'potent',0)
      call ffile(itap71,'file71',0)
      WRITE(6,1)                                                        FDG00630
      write(6,*) '  maxcor    = ',maxcor
      write(6,*) '  max_cal   = ',max_cal
      write(6,*) '  max_ijk   = ',max_ijk
C
C=============================================================
C   MOL_TYPE
C=============================================================
C    Asymmetric molecules
C     A(1) atom is located in the `origin'
C     B(2) atom is located along the `specified axis'
C     C(3) atom is located along the `specified axis'
C               or in the `specified plane'
C
C    Symmetric molecules
C     A(1) atom is located in the `origin'
C     B(2) atom is located along the `specified axis'
C                  or in the `specified plane'
C     B(3) atom is located along the `specified axis'
C                  or in the `specified plane'
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
C     f_{ijk} = R_{i} * R_{j} * P_{k}
C=============================================================
C   Power for each function starts from a value `0'
C     I_POWER   :  0, 1, 2, 3, ... ipower
C     J_POWER   :  0, 1, 2, 3, ... jpower
C     K_POWER   :  0, 1, 2, 3, ... kpower
C
C   Thus, the possible dimension of f_{ijk} is 
C         ndim = (ipower+1) * (jpower+1) * (kpower+1)
C   However, the folowing restriction is inforced
C         max_power <= i + j + k
C=============================================================
C
      errcod    = fwkclr()
      errcod    = fwkadd(5,'POTENT3')
      errcod    = fwkadd(5,'DEFAULT')
      errcod    = frdi('PRINT',iprnt)
      unit      = 'ATOMIC  '
      errcod    = frdc('UNIT',unit)
      wfn       = 'SCF  '
      errcod    = frdc('WFN',wfn)
      ipower    = 3
      errcod    = frdi('I_POWER',ipower)
      jpower    = 3
      errcod    = frdi('J_POWER',jpower)
      kpower    = 3
      errcod    = frdi('K_POWER',kpower)
      max_power = max0(ipower,jpower,kpower)
      errcod    = frdi('MAX_POWER',max_power)
      max_delta = 5000
      errcod    = frdi('MAX_DELTA',max_delta)
      del_max   = dfloat(max_delta)
      kadjst    = '2K'
      errcod    = frdc('K_ADJUST',kadjst)
      if(kadjst.eq.'2K') then
        k_pref  = 2
      end if
      if(kadjst.eq.'1K') then
        k_pref  = 1
      end if
      padjst    = 'EQUIL '
      errcod    = frdc('P_ADJUST',padjst)
      legdre    = 'NO '
      errcod    = frdc('LEGENDRE',legdre)
      spf       = 'NO '
      errcod    = frdc('SPF',spf)
C
      errcod = fwkadd(5,'MOLECULE')
      errcod = frdc('MOL_TYPE',moltyp)
      errcod = frdia1('RKL',kl,2)
      errcod = frdia1('RLM',lm,2)
      errcod = frdia1('AKLM',klm,3)
      errcod = frdia1('BOND_1',kl,2)
      errcod = frdia1('BOND_2',lm,2)
      errcod = frdia1('ANGLE',klm,3)
C
    
      ijk = 0
      do 101 i=1,ipower+1
      do 101 j=1,jpower+1
      do 101 k=1,kpower+1
      mpower = (i-1) + (j-1) + (k-1)
      if(mpower.gt.max_power) go to 101
      ijk = ijk + 1
  101 continue
      ndim = ijk
C
      REWIND ITAP11                                                     FDG00800
      READ(ITAP11,2) TITLE                                              FDG00810
      READ(ITAP11,3) NATOMS,ENERG0,DIPX,DIPY,DIPZ
C
      REWIND ITAP71                                                     FDG00800
      READ(ITAP71,4) NO_CALS
C
      n3n = natoms*3
      WRITE(6,5) WFN,UNIT,MOLTYP,IPRNT,NATOMS,N3N,
     *           legdre,spf,no_cals,ipower,jpower,kpower,max_power,
     *           ndim,max_delta,kadjst,padjst
C
      max_dim =  (ipower+1) * (jpower+1) * (kpower+1)
      if(ndim.gt.max_dim) then
        write(6,*) '  '
        write(6,*) '  ***************************'
        write(6,*) '  NDIM exceeds MAX_DIM !!!'
        write(6,*) '  NDIM    = ',ndim
        write(6,*) '  MAX_DIM = ',max_dim
        write(6,*) '  ---> STOP !!!'
        write(6,*) '  ***************************'
        write(6,*) '  '
        stop
      end if
      if(max_power.gt.max_dim) then
        write(6,*) '  '
        write(6,*) '  *****************************'
        write(6,*) '  MAX_POWER exceeds MAX_DIM !!!'
        write(6,*) '  MAX_POWER = ',max_power
        write(6,*) '  MAX_DIM   = ',max_dim
        write(6,*) '  ---> STOP !!!'
        write(6,*) '  *****************************'
        write(6,*) '  '
        stop
      end if
      if(no_cals.gt.max_cal) then
        write(6,*) '  '
        write(6,*) '  ***************************'
        write(6,*) '  NO_CALS exceeds MAX_CAL !!!'
        write(6,*) '  NO_CALS = ',no_cals
        write(6,*) '  MAX_CAL = ',max_cal
        write(6,*) '  ---> STOP !!!'
        write(6,*) '  ***************************'
        write(6,*) '  '
        stop
      end if
      if(ndim.gt.max_ijk) then
        write(6,*) '  '
        write(6,*) '  ***************************'
        write(6,*) '  NDIM exceeds MAX_IJK !!!'
        write(6,*) '  NDIM    = ',ndim
        write(6,*) '  MAX_IJK = ',max_ijk
        write(6,*) '  ---> STOP !!!'
        write(6,*) '  ***************************'
        write(6,*) '  '
        stop
      end if
C                                                                       FDG01130
C::::::::::::::::::::::::::::::::::::::::::::::::::::::
C   Read in geometries from `file11.dat' file
C::::::::::::::::::::::::::::::::::::::::::::::::::::::
      REWIND ITAP11                                                     FDG00800
      READ(ITAP11,2) TITLE                                              FDG00810
      READ(ITAP11,3) NATOMS,ENERG0,DIPX,DIPY,DIPZ
      WRITE(6,6)    
      DO 102 I=1,NATOMS                                                 FDG00830
      READ(ITAP11,7) CHG(I),XOLD(I),YOLD(I),ZOLD(I)
      WRITE(6,8) I,CHG(I),XOLD(I),YOLD(I),ZOLD(I)
  102 CONTINUE                                                          FDG00850
      call dist
      rie = r(kl(1),kl(2))
      rje = r(lm(1),lm(2))
      write(6,9) rie,rje
      call angl
      pke = aglb
C
      WRITE(6,10)
      DO 103 I=1,NATOMS                                                 FDG01160
      xa(i)=xold(i) * bohr
      ya(i)=yold(i) * bohr
      za(i)=zold(i) * bohr
      WRITE(6,8) I,CHG(I),XA(I),YA(I),ZA(I)
  103 CONTINUE                                                          FDG01180
      call length(xa,ya,za)
C                                                                       FDG01130
C::::::::::::::::::::::::::::::::::::::::::::::::::::::
C   Read in geometrical parameters
C::::::::::::::::::::::::::::::::::::::::::::::::::::::
      REWIND ITAP71                                                     FDG00800
      READ(ITAP71,4) NO_CALS
      write(6,11)
      do 105 i=1,no_cals
      read(itap71,12) nn,rri,rrj,ppk
      ri0(i)=rri
      rj0(i)=rrj
      pk0(i)=ppk
      write(6,13) nn,rri,rrj,ppk
  105 continue
      rewind itap71
C
      if(padjst.eq.'EQUIL ' .and. spf.eq.'NO ') then
        write(6,14)
        do 106 i=1,no_cals
        ri(i)=ri0(i)-rie
        rj(i)=rj0(i)-rje
        pk(i)=pk0(i)-pke
        write(6,13) i,ri(i),rj(i),pk(i)
  106   continue
      end if
C
      if(padjst.eq.'LINEAR' .and. spf.eq.'NO ') then
        write(6,15)
        do 107 i=1,no_cals
        ri(i)=ri0(i)-rie
        rj(i)=rj0(i)-rje
        pk(i)=pi-pk0(i)
        write(6,13) i,ri(i),rj(i),pk(i)
  107   continue
      end if
C
      if(padjst.eq.'EQUIL ' .and. spf.eq.'YES') then
        write(6,16)
        do 108 i=1,no_cals
        ri(i) = ( ri0(i)-rie ) / ( ri0(i) )
        rj(i) = ( rj0(i)-rje ) / ( rj0(i) )
        pk(i) = pk0(i)-pke
        write(6,13) i,ri(i),rj(i),pk(i)
  108   continue
      end if
C
      if(padjst.eq.'LINEAR' .and. spf.eq.'YES') then
        write(6,17)
        do 109 i=1,no_cals
        ri(i) = ( ri0(i)-rie ) / ( ri0(i) )
        rj(i) = ( rj0(i)-rje ) / ( rj0(i) )
        pk(i) = pi-pk0(i)
        write(6,13) i,ri(i),rj(i),pk(i)
  109   continue
      end if
C
C===========================================
C  Dynamic allocation
C===========================================
C  (1)  Fijk (ndim*no_cals)
C  (2)  Cijk (ndim)
C  (3)  Diff (no_cals)
C===========================================
      ic1   = 1
      ic2   = ic1 + ndim * no_cals
      ic3   = ic2 + ndim 
      icmax = ic3 + no_cals
      if(icmax.gt.maxcor) go to 400
C===========================================
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::
C   Form the Fijk matrix
C::::::::::::::::::::::::::::::::::::::::::::::::::::::
      if(legdre.eq.'NO ') then
C....................Fijk.............
        call mk_fijk(cc(ic1),ri,rj,pk)
      end if
C
      if(legdre.eq.'YES') then
C.....................Fijk.............
        call mk_fijkL(cc(ic1),ri,rj,pk)
      end if
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::
C   Read in the Cijk vector
C::::::::::::::::::::::::::::::::::::::::::::::::::::::
C..................Cijk....
      call rd_cijk(cc(ic2))
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::
C   Test the analytic potential
C::::::::::::::::::::::::::::::::::::::::::::::::::::::
C................Fijk    Cijk    Diff....
      call ptest(cc(ic1),cc(ic2),cc(ic3))
      go to 500
C
  400 continue
      write(6,*) '  '
      write(6,*) '  ***************************'
      write(6,*) '  ICMAX exceeds MAXCOR !!!'
      write(6,*) '  ICMAX   = ',icmax
      write(6,*) '  MAXCOR  = ',maxcor
      write(6,*) '  ---> STOP !!!'
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
      subroutine length(x,y,z)
      IMPLICIT REAL*8 (A-H,O-Z)                              
      dimension x(10),y(10),z(10)
      dimension r(10,10)
      common/com101/natoms
    1 FORMAT(//,2X,' INTERATOMIC DISTANCE MATRIX'/)           
C                                                            
      call zero(r,100)
      DO 101 I=1,NATOMS
      DO 101 J=1,NATOMS
      IF(I.EQ.J) GO TO 101                                  
      XD=X(I)-X(J)                                         
      YD=Y(I)-Y(J)                                        
      ZD=Z(I)-Z(J)                                       
      valu=DSQRT(XD*XD+YD*YD+ZD*ZD)
      r(i,j)=valu
  101 CONTINUE                                          
      WRITE(6,1)                                       
      CALL MATOUT(R,10,10,NATOMS,NATOMS,6)
c
      return
      end
      subroutine angl
      implicit real*8 (a-h,o-z)
      dimension ex(10,10),ey(10,10),ez(10,10)
      common/com101/natoms
      common/com102/iprnt
      common/com105/apir,aglb,ang
      common/com202/x(10),y(10),z(10)
      common/com205/r(10,10)
      common/com206/kl(2),lm(2),klm(3)
      data pi / 3.1415926536D+00 /
    1 format(//,2x,' Ex matrix'/)
    2 format(//,2x,' Ey matrix'/)
    3 format(//,2x,' Ez matrix'/)
c
      apir=pi/180.0D+00
      call zero(ex,100)
      call zero(ey,100)
      call zero(ez,100)
      DO 101 I=1,NATOMS
      DO 101 J=1,NATOMS
      if(i.eq.j) go to 101
      EX(I,J)=-(X(I)-X(J))/R(I,J)
      EY(I,J)=-(Y(I)-Y(J))/R(I,J)
      EZ(I,J)=-(Z(I)-Z(J))/R(I,J)
  101 CONTINUE                              
      if(iprnt.ge.5) then
        write(6,1)
        call matout(ex,10,10,natoms,natoms,6)
        write(6,2)
        call matout(ey,10,10,natoms,natoms,6)
        write(6,3)
        call matout(ez,10,10,natoms,natoms,6)
      end if
c
C  Bond angle is defined as K-L-M
      k=klm(1)
      l=klm(2)
      m=klm(3)
      AKLM=EX(L,K)*EX(L,M)+EY(L,K)*EY(L,M)+EZ(L,K)*EZ(L,M)
      AGLB=DACOS(AKLM)
      ANG=(AGLB)/(apir)
      write(6,*) ' '
      WRITE(6,*) '  Bond angle (k,l,m)    = ',k,l,m
      WRITE(6,*) '  Bond angle in radians = ',aglb
      WRITE(6,*) '  Bond angle in degrees = ',ang
c
      return
      end
      subroutine mk_fijk(fijk,ri,rj,pk)
      implicit real*8 (a-h,o-z)
      dimension ri(no_cals),rj(no_cals),pk(no_cals)
      dimension fijk(ndim,no_cals)
      common/com102/iprnt
      common/pot101/no_cals,ndim
      common/pot102/ipower,jpower,kpower,max_power
      common/pot106/k_pref
    1 format(//,2x,' The Fijk matrix'/)
C
C  Form the Fijk matrix
      if(iprnt.ge.5) then
        write(6,*) '  '
        write(6,*) ' iprnt in mk_fijk     = ',iprnt
        write(6,*) ' no_cals in mk_fijk   = ',no_cals
        write(6,*) ' ndim in mk_fijk      = ',ndim
        write(6,*) ' ipower in mk_fijk    = ',ipower
        write(6,*) ' jpower in mk_fijk    = ',jpower
        write(6,*) ' kpower in mk_fijk    = ',kpower
        write(6,*) ' max_power in mk_fijk = ',max_power
        write(6,*) ' k_pref in mk_fijk    = ',k_pref
      end if
C
      if(iprnt.ge.15) then
        write(6,*) '  '
        write(6,*) ' ipts, ijk, i, j, k'
        write(6,*) ' rri, rrj, ppk'
        write(6,*) ' rivalu, rjvalu, pkvalu, fvalu'
        write(6,*) '  '
      end if
C
      call zero(fijk,ndim*no_cals)
      do 105 ipts = 1, no_cals
      rri = ri(ipts)
      rrj = rj(ipts)
      ppk = pk(ipts)
C
        ijk = 0
        do 103 i = 1, ipower+1
        rivalu = rri ** (i-1)
        do 102 j = 1, jpower+1
        rjvalu = rrj ** (j-1)
        rijval = rivalu * rjvalu
        do 101 k = 1, kpower+1
C
        mpower = (i-1) + (j-1) + (k-1)
        if(mpower.gt.max_power) go to 101
        kvalu = k_pref*(k-1)
        pkvalu = ppk ** (kvalu)
        ijk = ijk + 1
        fvalu = rijval * pkvalu
        if(iprnt.ge.15) then
          write(6,*) ipts, ijk, i, j, k
          write(6,*) rri,rrj,ppk
          write(6,*) rivalu,rjvalu,pkvalu,fvalu
        end if
        fijk(ijk,ipts) = fvalu
  101   continue
  102   continue
  103   continue
  105 continue
C
      if(iprnt.ge.3) then
        write(6,*) ' '
        write(6,*) '  ijk in mk_fijk = ',ijk
        write(6,*) '  '
      end if
      if(iprnt.ge.5) then
        write(6,1)
        call matout(fijk,ndim,no_cals,ndim,no_cals,6)
      end if
C
      return
      end
      subroutine mk_fijkL(fijk,ri,rj,pk)
      implicit real*8 (a-h,o-z)
      dimension ri(no_cals),rj(no_cals),pk(no_cals)
      dimension fijk(ndim,no_cals)
      common/com102/iprnt
      common/pot101/no_cals,ndim
      common/pot102/ipower,jpower,kpower,max_power
      common/pot106/k_pref
    1 format(//,2x,' The Fijk matrix'/)
C
C  Form the Fijk matrix
      if(iprnt.ge.5) then
        write(6,*) '  '
        write(6,*) ' iprnt in mk_fijkL     = ',iprnt
        write(6,*) ' no_cals in mk_fijkL   = ',no_cals
        write(6,*) ' ndim in mk_fijkL      = ',ndim
        write(6,*) ' ipower in mk_fijkL    = ',ipower
        write(6,*) ' jpower in mk_fijkL    = ',jpower
        write(6,*) ' kpower in mk_fijkL    = ',kpower
        write(6,*) ' max_power in mk_fijkL = ',max_power
        write(6,*) ' k_pref in mk_fijkL    = ',k_pref
      end if
C
      if(iprnt.ge.15) then
        write(6,*) '  '
        write(6,*) ' ipts, ijk, i, j, k'
        write(6,*) ' rri, rrj, ppk'
        write(6,*) '  '
      end if
C
      call zero(fijk,ndim*no_cals)
      do 105 ipts = 1, no_cals
      rri = ri(ipts)
      rrj = rj(ipts)
      ppk = pk(ipts)
C
        ijk = 0
        do 103 i = 1, ipower+1
        rivalu = rri ** (i-1)
        do 102 j = 1, jpower+1
        rjvalu = rrj ** (j-1)
        rijval = rivalu * rjvalu
        do 101 k = 1, kpower+1
C
        mpower = (i-1) + (j-1) + (k-1)
        if(mpower.gt.max_power) go to 101
C*      kvalu  = k_pref*(k-1)
C*      pkvalu = pcal(ppk,kvalu)
        kvalu  = k-1
        pktemp = pcal(ppk,kvalu)
        pkvalu = pktemp ** (k_pref)
        ijk = ijk + 1
        fvalu = rijval * pkvalu
        if(iprnt.ge.15) then
          write(6,*) ipts, ijk, i, j, k
          write(6,*) rri,rrj,ppk
          write(6,*) rivalu,rjvalu,pkvalu,fvalu
        end if
        fijk(ijk,ipts) = fvalu
  101   continue
  102   continue
  103   continue
  105 continue
C
      if(iprnt.ge.3) then
        write(6,*) ' '
        write(6,*) '  ijk in mk_fijk = ',ijk
        write(6,*) '  '
      end if
      if(iprnt.ge.5) then
        write(6,1)
        call matout(fijk,ndim,no_cals,ndim,no_cals,6)
      end if
C
      return
      end
      function pcal(pvalu,lpower)
      implicit real*8 (a-h,o-z)
      dimension q(100)
      data one, two / 1.0d+00 , 2.0d+00 /
    1 format(2x,i5,2f20.10)
C
C  Form the Legendre polynmials
      q(1) = one
      if(lpower.eq.0) go to 200
C
      zvalu = dcos(pvalu)
      q(2)  = zvalu
      if(lpower.eq.1) go to 200
C
      llim = lpower-1
      do 101 l = 1, llim
      rl     = dfloat(l)
      fact   = two*rl + one
      denom  = rl + one
      valu   = fact*zvalu*q(l+1) - rl*q(l)
      q(l+2) = valu / ( denom )
  101 continue
C
  200 continue
      pcal = q(lpower+1)
C
      return
      end
      subroutine rd_cijk(cijk)
      implicit real*8 (a-h,o-z)
      dimension cijk(ndim)
      common/com102/iprnt
      common/pot101/no_cals,ndim
      common/pot105/eref
    1 format(i5,15x,f20.10)
    2 format(//,2x,' Cijk vector'/)
    3 format(4i5,f20.10)
    4 format(2x,i5,2x,'(',i3,i3,i3,')',5x,f20.10)
C
C  Read in the Cijk vector from file63 (potent.dat) file
      itap63 = 63
      if(iprnt.ge.5) then
        write(6,*) '  '
        write(6,*) ' iprnt in rd_cijk     = ',iprnt
        write(6,*) ' no_cals in rd_cijk   = ',no_cals
        write(6,*) ' ndim in rd_cijk      = ',ndim
        write(6,*) ' itap63 in rd_cijk    = ',itap63
      end if
C
      rewind itap63
      read(itap63,1) nndim,eref
      if(nndim.ne.ndim) then
        write(6,*) '  '
        write(6,*) ' nndim is not equal to ndim ---> STOP!'
        write(6,*) ' nndim = ',nndim
        write(6,*) ' ndim  = ',ndim
        write(6,*) '  '
      end if
C
      call zero(cijk,ndim)
      write(6,2)
      DO 101 I=1,ndim
      read(itap63,3) ijk,ii,jj,kk,valu
      cijk(I)=valu
      write(6,4) ijk,ii,jj,kk,valu
  101 CONTINUE   
      rewind itap63
C
      return
      end
      subroutine ptest(fijk,cijk,diff)
      implicit real*8 (a-h,o-z)
      dimension fijk(ndim,no_cals),cijk(ndim)
      dimension diff(no_cals)
      character*2 wover
      common/com102/iprnt
      common/pot101/no_cals,ndim
      common/pot103/del_max
      common/pot105/eref
      data a00 / 0.0d+00 /
      data wave / 219474.63d+00 /
    1 format(//,2x,' The Cijk vector'/)
    2 format(2x,i5,5x,f20.10)
    3 format(//,2x,' The Fijk matrix'/)
    4 format(//,2x,' Reference energy = ',f20.10)
    5 format(//,2x,' IPTS          ANALYTICAL',
     1 10X,'DELTA',19x,'DELTA'/
     2 16X,10H (hartree),10x,10H (hartree),
     3 14x,7H (CM-1),/)
    6 format(2x,i5,2f20.10,4x,f20.10,a2)
C
C  Evaluate an omega value
      if(iprnt.ge.5) then
        write(6,*) '  '
        write(6,*) ' iprnt in ptest   = ',iprnt
        write(6,*) ' no_cals in ptest = ',no_cals
        write(6,*) ' ndim in ptest    = ',ndim
      end if
C
      if(iprnt.ge.10) then
        write(6,1)
        do 101 i=1,ndim
        valu = cijk(i)
        write(6,2) i,valu
  101   continue
        write(6,3)
        call matout(fijk,ndim,no_cals,ndim,no_cals,6)
      end if
C
      call zero(diff,no_cals)
      write(6,4) eref
      write(6,5)
      do 103 ipts = 1, no_cals
C
      valu = a00
      do 102 ijk = 1, ndim
      valu = valu + cijk(ijk) * fijk(ijk,ipts)
  102 continue
      delta = valu - eref
      diff(ipts) = delta
      delw = delta * wave
      wover = '  '
      if(dabs(delw).gt.del_max) then
        wover = ' *'
      end if
      write(6,6) ipts,valu,delta,delw,wover
  103 continue
      write(6,*) ' '
C
      return
      end
