      subroutine fentry(cc,ia,maxcor)
c**********************************************************************
c     PROGRAM CS3_MP06                                                *
c**********************************************************************
c   Last updated on January 05, 2007 by Dr. Yukio Yamaguchi           *
c   Modified for a bent C_{s} type molecule                           *
c**********************************************************************
c   Last updated on January 11, 2006 by Dr. Yukio Yamaguchi           *
c   Adjustment to the MOLPRO program package                          *
c**********************************************************************
c   This program determines force constants and frequencies           *
c   for a bent C_{s} type molecule                                    *
c**********************************************************************
      implicit real*8 (a-h,o-z)
      dimension cc(maxcor),ia(maxcor*2)
      parameter (max_cal=25)
      dimension ri(max_cal),rj(max_cal),pk(max_cal),energ(max_cal)
      dimension ti(max_cal),tj(max_cal),tk(max_cal),tener(max_cal)
      dimension ATM(36)
      dimension WISO(20,3)
      dimension GXX(6),FXX(3,3)
      character*3  p_side
      character*6  ref_en
      character*8  moltyp
      character*8  unit_r, unit_a
      character*12 vibtyp
      character*12 wfn
      character*80 atitle
      common/com101/natoms
      common/com102/iprnt
      common/com201/chg(10),w(10)
      common/com202/x(10),y(10),z(10)
      common/com205/r(10,10)
      common/com206/kl(2),lm(2),klm(3)
      common/vib101/wave
      common/vib102/pertr,perta
      common/vib104/cosval
      data xpertr, xperta / 0.005d+00, 0.01d+00 /
      data one,two / 1.0d+00 , 2.0d+00  /
      DATA ATM / 1.007825D+00 ,  4.00260D+00 ,  7.01600D+00 ,
     *           9.01218D+00  , 11.00931D+00 , 12.00000D+00 ,           NOR00620
     *           14.00307D+00 , 15.99491D+00 , 18.99840D+00 ,           NOR00630
     *           19.99244D+00 , 22.98977D+00 , 23.98504D+00 ,           NOR00640
     *           26.98154D+00 , 27.97693D+00 , 30.97376D+00 ,           NOR00650
     *           31.97207D+00 , 34.96885D+00 , 39.96238D+00 ,           NOR00660
     *           38.96371D+00 , 39.96259D+00 , 44.95591D+00 ,           NOR00670
     *           47.94795D+00 , 50.94396D+00 , 51.94051D+00 ,           NOR00680
     *           54.93805D+00 , 55.93494D+00 , 58.93320D+00 ,           NOR00690
     *           57.93535D+00 , 62.92960D+00 , 63.92915D+00 ,           NOR00700
     *           68.92558D+00 , 73.92118D+00 , 74.92160D+00 ,           NOR00710
     *           79.91652D+00 , 78.91834D+00 , 83.91151D+00   /         NOR00720
      DATA PH,BK,AVN / 6.626176D+00 , 1.380662D+00 , 6.022045D+00 /     NOR00730
      DATA A0,HE,EV / 0.529177249D+00 , 4.359813653D+00 , 1.6021892D+00/NOR00740
      DATA CL,TEMP,DANG3 / 2.99792458D+00 , 298.15D+00 , 1.4409726D+01 /NOR00750
      DATA PI,ALN / 3.1415926536D+00 , 2.3025850930D+00 /               NOR00760
c==================================================================
c  For psi 2.0 input
#include <error.h>
      integer frdc, frdi, frdf, fwkclr, fwkadd, frdboo, frdpos
      integer fcnt, frdia1
      integer errcod
c==================================================================
    1 FORMAT(/,2X,' ============================================='/
     1        ,2X,' THE NUMERICAL DIFFERENTIATION PROGRAM        '/
     2         2x,' to determine force constants and frequencies '/
     3         2x,' for a linear C_{s} type triatomic molecule   '/
     4         2x,' adjusted to the MOLPRO program package       '/
     5         2x,'         (Version of 01/05/2007)              '/ 
     6         2x,'              (16:30 PM)                      '/ 
     7        ,2X,' ============================================='/)
    2 FORMAT(A80)                                                       MK100320
    3 FORMAT(3X,A80/3X,'(title in file11)'/) 
    4 FORMAT(3X,A80/3X,'(title in file61)') 
    5 FORMAT(//,2X,' PARAMETERS'/
     1 2X,' WFN        = ',A12/
     2 2X,' UNIT_R     = ',A8/
     3 2X,' UNIT_A     = ',A8/
     4 2X,' REF_ENERGY = ',A6/
     5 2X,' MOL_TYPE   = ',A8/
     6 2X,' VIB_TYPE   = ',A12/
     7 2X,' PERT_SIDES = ',A3/
     8 2X,' PERTR      = ',f10.5/
     9 2X,' PERTA      = ',f10.5/
     a 2X,' IPRNT      = ',I5/
     b 2X,' NATOMS     = ',I5/
     c 2X,' N3N        = ',I5/
     d 2X,' IVTYPE     = ',I5/
     e 2X,' ISOTOP     = ',I5/
     f 2X,' NPOINTS    = ',I5/
     g 2X,' FORMULA    = ',I5/)
    6 FORMAT(//,2X,' REFERENCE GEOMETRY IN ATOMIC UNIT (from file11)'/
     1 5X,4H NO.,11X,7H CHARGE,13X,2H X,18X,2H Y,18X,2H Z/)             MK100410
    7 FORMAT(10X,F8.2,3F15.9) 
    8 FORMAT(2X,I7,5X,4F20.10)                                          MK100430
    9 format(//,2x,' Equilibrium bond distances in atomic unit (from fil
     1e11)'/)
   10 format(4F20.10)
   11 format(//,2x,' PERTURBED GEOMETRIES AND ENERGIES'/
     1 4X,4H NO.,6X,4H R_i,14X,4H R_j,14X,4H P_k,19X,7H Energy/
     2 14X,7H (Bohr),11X,7H (Bohr),11X,9H (radian),14X,10H (hartree)/)
   12 format(//,2x,' PERTURBED GEOMETRIES AND ENERGIES'/
     1 4X,4H NO.,6X,4H R_i,14X,4H R_j,14X,4H P_k,19X,7H Energy/
     2 14X,7H (Bohr),11X,7H (Bohr),11X,9H (degree),14X,10H (hartree)/)
   13 format(//,2x,' PERTURBED GEOMETRIES AND ENERGIES'/
     1 4X,4H NO.,6X,4H R_i,14X,4H R_j,14X,4H P_k,19X,7H Energy/
     2 14X,7H (A   ),11X,7H (A   ),11X,9H (radian),14X,10H (hartree)/)
   14 format(//,2x,' PERTURBED GEOMETRIES AND ENERGIES'/
     1 4X,4H NO.,6X,4H R_i,14X,4H R_j,14X,4H P_k,19X,7H Energy/
     2 14X,7H (A   ),11X,7H (A   ),11X,9H (degree),14X,10H (hartree)/)
   15 format(//,2x,' PERTURBED GEOMETRIES AND ENERGIES'/
     1 4X,4H NO.,6X,7H R_i-Re,11X,7H R_j-Re,11X,
     2 7H P_k-Pe,16X,7H Energy,/
     3 14X,7H (Bohr),11X,7H (Bohr),11X,9H (radian),14X,10H (hartree)/)
   16 format(2x,i5,3x,f15.10,3x,f15.10,3x,f15.10,5x,f20.12,2x,f15.6)
   17 format(f20.10)
   18 FORMAT(///,2X,' REFERENCE GEOMETRY IN ANGSTROM'/
     1 5X,4H NO.,11X,12H ATOMIC MASS,8X,2H X,18X,2H Y,18X,2H Z/)   
   19 format(//,2x,' Equilibrium bond distances in angstrom (from file11
     1)'/
     2          2x,' Bond_1(e)       = ',f20.10/
     3          2x,' Bond_2(e)       = ',f20.10/)
   20 format(//,2x,' Equilibrium bond distances in atomic unit (from fil
     1e61)'/
     2          2x,' Bond_1(e)       = ',f20.10/
     3          2x,' Bond_2(e)       = ',f20.10/)
   21 format(   2x,' Pertr in bohr   = ',f20.10/
     1          2x,' Perta in radian = ',f20.10)
C                                                
      call psinit('APPEND')
      CALL TSTART(6)                                                    NOR01250
      input=5
      itap11=11
      itap61=61
      call ffile(itap11,'file11',0)
      call ffile(itap61,'vibener',0)
c
      write(6,1)
      write(6,*) '  '
      write(6,*) '  maxcor     = ',maxcor
      write(6,*) '  '
C
C   PHYSICAL CONSTANT
      WAVE=1.0D+04*DSQRT(AVN)/(2.0D+00*PI*CL)                           NOR01300
C
C=============================================================
C   MOL_TYPE
C=============================================================
C
C     A(1)-B(2)-C(3)   C_{s} point group symmetry
C    'abc_bent'
C     R1 = R(1,2)
C     R2 = R(2,3)
C     R3 = ANG(1,2,3)
C
C     B(2)-A(1)-C(3)   C_{s} point group symmetry
C    'bac_bent'
C     R1 = R(1,2)
C     R2 = R(1,3)
C     R3 = ANG(2,1,3)
C
C     C(3)-A(1)-B(2)   C_{s} point group symmetry
C    'cab_bent'
C     R1 = R(1,2)
C     R2 = R(1,3)
C     R3 = ANG(2,1,3)
C
C=============================================================
C
      errcod = fwkclr()
      errcod = fwkadd(5,'CS3')
      errcod = fwkadd(5,'DEFAULT')
      iprnt  = 0
      errcod = frdi('PRINT',iprnt)
      isotop = 1
      errcod = frdi('ISOTOP',isotop)
      wfn    = 'SCF         '
      errcod = frdc('WFN',wfn)
C
      errcod = fwkadd(5,'PERTURB')
      unit_r = 'ATOMIC  '
      errcod = frdc('UNIT_R',unit_r)
      unit_a = 'RADIAN  '
      errcod = frdc('UNIT_A',unit_a)
      pertr  = xpertr
      errcod = frdf('PERT_R',pertr)
      perta  = xperta
      errcod = frdf('PERT_A',perta)
      p_side = 'TWO'
      errcod = frdc('PERT_SIDES',p_side)
      ref_en = 'TOP   '
      errcod = frdc('REF_ENERGY',ref_en)
      vibtyp = 'ALL         '
      errcod = frdc('VIB_TYPE',vibtyp)
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
C   READ IN THE REFERENCE GEOMETRY FROM FILE11
      REWIND ITAP11                                                     MK100850
      READ(ITAP11,2) ATITLE                                             MK100860
      WRITE(6,3) ATITLE                                                 MK100870
      natoms = 0
  200 continue
      READ(ITAP11,2,end=201) ATITLE
      natoms = natoms + 1
      go to 200
  201 continue
C
C   READ IN ENERGY POINTS FROM FILE61
      REWIND ITAP61                                                     MK100850
      READ(ITAP61,2) ATITLE                                             MK100860
      WRITE(6,4) ATITLE                                                 MK100870
      READ(ITAP61,2) ATITLE                                             MK100860
      npoints = 0
  202 continue
      READ(ITAP61,2,end=203) ATITLE
      npoints = npoints + 1
      go to 202
  203 continue
C
      nabc = index(moltyp,'ABC_BENT')
      nbac = index(moltyp,'BAC_BENT')
      ncab = index(moltyp,'CAB_BENT')
      nmax = max(nabc, nbac, ncab)
      if(nmax .eq. 0) then
        write(6,*) '  '
        write(6,*) '  ******************************'
        write(6,*) '  You cannot have this type of  '
        write(6,*) '  MOL_TYPE                      '
        write(6,*) '  MOL_TYPE = ',moltyp
        write(6,*) '  NABC     = ',nabc
        write(6,*) '  NBAC     = ',nbac
        write(6,*) '  NCAB     = ',ncab
        write(6,*) '  NMAX     = ',nmax
        write(6,*) '  ---> STOP !!!                 '
        write(6,*) '  ******************************'
        write(6,*) '  '
        stop
      end if
C
      if(natoms .ne. 3) then
        write(6,*) '  '
        write(6,*) '  ******************************'
        write(6,*) '  You cannot have this number of'
        write(6,*) '  NATOMS                        '
        write(6,*) '  NATOMS   = ',natoms
        write(6,*) '  ---> STOP !!!                 '
        write(6,*) '  ******************************'
        write(6,*) '  '
        stop
      end if
C
      if(npoints .le. 1 .or. npoints .gt. 25) then
        write(6,*) '  '
        write(6,*) '  ******************************'
        write(6,*) '  You cannot have a number      '
        write(6,*) '  for NPOINTS                   '
        write(6,*) '  NPOINTS = ',npoints
        write(6,*) '  ---> STOP !!!                 '
        write(6,*) '  ******************************'
        write(6,*) '  '
        stop
      end if
C
      n3n=natoms*3
      if(isotop .eq. 0) then
        isotop = 1
      end if
      nform = 5
C
      WRITE(6,5) WFN,UNIT_R,UNIT_A,REF_EN,MOLTYP,VIBTYP,
     * P_SIDE,PERTR,PERTA,
     * IPRNT,NATOMS,N3N,IVTYPE,ISOTOP,NPOINTS,NFORM
C
      REWIND ITAP11                                                     MK100850
      READ(ITAP11,2) ATITLE                                             MK100860
      WRITE(6,6)                                                        MK100910
      DO 101 I=1,NATOMS                                       
      READ(ITAP11,7) CHG(I),X(I),Y(I),Z(I)                              MK100930
      WRITE(6,8) I,CHG(I),X(I),Y(I),Z(I)  
  101 CONTINUE                                                          MK100950
      call dist
      rie = r(kl(1),kl(2))
      rje = r(lm(1),lm(2))
      write(6,9) 
      WRITE(6,*) '  Bond_1 (k,l)    = ',kl(1),kl(2),rie
      WRITE(6,*) '  Bond_2 (l,m)    = ',lm(1),lm(2),rje
C
      call angl
      apir=pi/180.0D+00
C
C================================================
c   Obtain perturbed energies from file61
C================================================
C
      call zero(ti,max_cal)
      call zero(tj,max_cal)
      call zero(tk,max_cal)
      call zero(tener,max_cal)
      call zero(ri,max_cal)
      call zero(rj,max_cal)
      call zero(pk,max_cal)
      call zero(energ,max_cal)
C
      rewind itap61  
      read(itap61,2) atitle
      read(itap61,2) atitle
      do 102 i=1,npoints
      read(itap61,10) ti(i),tj(i),tk(i),tener(i)
  102 continue
      rewind itap61
C
      if(ref_en .eq. 'TOP   ') then
        do 103 i=1,npoints
        ri(i)=ti(i)
        rj(i)=tj(i)
        pk(i)=tk(i)
        energ(i)=tener(i)
  103   continue
      end if
      if(ref_en .eq. 'BOTTOM') then
        ri(1)=ti(npoints)
        rj(1)=tj(npoints)
        pk(1)=tk(npoints)
        energ(1)=tener(npoints)
        do 104 i=1,npoints-1
        ri(i+1)=ti(i)
        rj(i+1)=tj(i)
        pk(i+1)=tk(i)
        energ(i+1)=tener(i)
  104   continue
      end if
C
      write(6,*) ' '
      if(unit_r.eq.'ATOMIC  ' .and. unit_a.eq.'RADIAN  ') then
        write(6,11)
        do 105 i=1,npoints
        write(6,16) i,ri(i),rj(i),pk(i),energ(i)
  105   continue
      end if
      if(unit_r.eq.'ATOMIC  ' .and. unit_a.eq.'DEGREE  ') then
        perta = perta * apir
        write(6,12)
        do 106 i=1,npoints
        write(6,16) i,ri(i),rj(i),pk(i),energ(i)
  106   continue
        write(6,11)
        do 107 i=1,npoints
        pk(i) = pk(i) * apir
        write(6,16) i,ri(i),rj(i),pk(i),energ(i)
  107   continue
      end if
      if(unit_r.eq.'ANGSTROM' .and. unit_a.eq.'RADIAN  ') then
        pertr = pertr / ( a0 )
        write(6,13)
        do 108 i=1,npoints
        write(6,16) i,ri(i),rj(i),pk(i),energ(i)
  108   continue
        write(6,11)
        do 109 i=1,npoints
        ri(i) = ri(i) / (a0)
        rj(i) = rj(i) / (a0)
        write(6,16) i,ri(i),rj(i),pk(i),energ(i)
  109   continue
      end if
      if(unit_r.eq.'ANGSTROM' .and. unit_a.eq.'DEGREE  ') then
        pertr = pertr / ( a0 )
        perta = perta * apir
        write(6,14)
        do 110 i=1,npoints
        write(6,16) i,ri(i),rj(i),pk(i),energ(i)
  110   continue
        write(6,11)
        do 111 i=1,npoints
        ri(i) = ri(i) / (a0)
        rj(i) = rj(i) / (a0)
        pk(i) = pk(i) * apir
        write(6,16) i,ri(i),rj(i),pk(i),energ(i)
  111   continue
      end if
C
      ri0 = ri(1)
      rj0 = rj(1)
      pk0 = pk(1)
      write(6,15)
      do 112 i=1,npoints
      ri(i)=ri(i)-ri0
      rj(i)=rj(i)-rj0
      pk(i)=pk(i)-pk0
      write(6,16) i,ri(i),rj(i),pk(i),energ(i)
  112 continue
C
C==============================================================
C   READ IN ISOTOPIC ATOMIC MASS
C   Isotopic atomic weights are read in from "input.dat" file
C==============================================================
      if(isotop.gt.1) then
        input = 5
        CALL LOCATE(INPUT,'# ISOTOPE #',IERR)
        DO 115 ISO=2,ISOTOP
        READ(5,17) (WISO(ISO,I),I=1,NATOMS)
  115   CONTINUE
      end if
C
C   VIBRATIONAL ANALYSIS FOR ISOTOPIC MOLECULES                         NOR03670
      NISO=0                                                            NOR03700
  400 continue
      NISO=NISO+1                                                       NOR03710
      IF(NISO.EQ.1) GO TO 401                                           NOR03720
      IF(NISO.GT.ISOTOP) GO TO 405                                      NOR03730
      write(6,*) ' '
      write(6,*) ' '
      write(6,*) ' ================================='
      write(6,*) ' Isotopic molecule number ',niso
      write(6,*) ' ================================='
      DO 116 I=1,NATOMS
      W(I)=WISO(NISO,I)                                                 NOR03800
  116 CONTINUE                                                          NOR03820
C
  401 continue
      write(6,18)
      DO 117 I=1,NATOMS
      if(niso.eq.1) then
        ii=chg(i)
        w(i)=atm(ii)
        x(i) = x(i) * A0
        y(i) = y(i) * A0
        z(i) = z(i) * A0
      end if
      WRITE(6,8) I,W(I),X(I),Y(I),Z(I)
  117 CONTINUE                         
      call dist
      rie = r(kl(1),kl(2))
      rje = r(lm(1),lm(2))
      write(6,19) rie,rje
      call angl
      write(6,20) ri0,rj0
      write(6,21) pertr,perta
C
C========================================
c   Determine the G matrix elements
c   for the linear C_{infv} type molecule
C========================================
c
      k     = klm(1)
      l     = klm(2)
      m     = klm(3)
      rwk   = one/w(k)
      rwl   = one/w(l)
      rwm   = one/w(m)
      rkl   = one/(rie)
      rml   = one/(rje)
      rklsq = one/(rie * rie)
      rmlsq = one/(rje * rje)
      rklml = one/(rie * rje)
      sinval = dsqrt( one - cosval * cosval )
c   For the R1 mode
      g11 = rwk + rwl
c   For the R2 mode
      g22 = rwk + rwl
c   For the R3 mode
      g33 = rklsq * rwk + rmlsq * rwm
     1    + ( rklsq + rmlsq - two * cosval * rklml ) * rwl
c   For the R1 X R2 interaction
      g12 = rwl * cosval
c   For the R1 X R3 interaction
      g13 = - rml * sinval * rwl
c   For the R2 X R3 interaction
      g23 = - rkl * sinval * rwl
C
      write(6,*) '  '
      write(6,*) '  '
      write(6,*) '  '
      write(6,*) ' :::::::::::::::::::::::::::::::::::'
      write(6,*) ' :  The G matrix element           :'
      write(6,*) ' :::::::::::::::::::::::::::::::::::'
      write(6,*) '  '
      write(6,*) ' The G matrix elements for the linear molecule'
      write(6,*) ' [in spectroscopic units]'
      write(6,*) '  '
      write(6,*) ' g11       = ',g11
      write(6,*) ' g22       = ',g22
      write(6,*) ' g33       = ',g33
      write(6,*) ' g12       = ',g12
      write(6,*) ' g13       = ',g13
      write(6,*) ' g23       = ',g23
      write(6,*) '  '
C
C  For all modes
      GXX(1) = g11
      GXX(2) = g12
      GXX(3) = g22
      GXX(4) = g13
      GXX(5) = g23
      GXX(6) = g33
      if(niso.eq.1 ) then
        call mk_force(energ,fxx)
      end if
      call normco(gxx,fxx)
C
C=====================================
C  END OF AN ISOTOPIC MOLECULE
C=====================================
      go to 400
C
  405 continue
      call tstop(6)
      stop
      end
      subroutine dist
      implicit real*8 (a-h,o-z)
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
      write(6,1)
      call matout(r,10,10,natoms,natoms,6)
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
      common/vib104/cosval
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
      if(iprnt.ge.3) then
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
      cosval=aklm
      ANG=(AGLB)/(apir)
      write(6,*) ' '
      WRITE(6,*) '  Bond angle (k,l,m)    = ',k,l,m
      WRITE(6,*) '  Bond angle in radians = ',aglb
      WRITE(6,*) '  Bond angle in degrees = ',ang
      write(6,*) ' '
c
      return
      end
