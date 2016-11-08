      subroutine fentry(cc,ia,maxcor)
c**********************************************************************
C     PROGRAM EXTR3D_MP0703                                           *
c     This program will extract energy information                    *
c       for an analytic potential of a triatomic molecule             *
C     Adjusted to the MOLPRO programs                                 *
C=====================================================================*
C  Structure of the program                                           *
C=====================================================================*
C  A. Read in the reference geometry from "file11.dat" file           *
C  B. Read in perturbed geometrical parameters, dipole moments,       *
C       and total energies from "file81.dat" file                     *
C  C. Extract the energy points less than MAX_DELTA (in cm-1)         *
C       and within criteria                                           *
C  D. Store the extracted energy points                               *
C         in "extr81.dat (file87)" file                               *
C   (Extracted energy points are less than MAX_DELTA                  *
C      and R1, R2, and R3 values are within criteria)                 *
C=====================================================================*
C*   NOTICE OF PROGRAM MODIFICATION                                   *
c**********************************************************************
C   Last updated on 03-29-2007 by Dr. Yukio Yamaguchi                 *
c**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (max_cal=10000)
      dimension cc(maxcor),ia(maxcor*2)
      dimension RI0(max_cal),RJ0(max_cal),PK0(max_cal)
      dimension DIPZ0(max_cal),DIPX0(max_cal),ENER0(max_cal)
      dimension RIT(max_cal),RJT(max_cal),PKT(max_cal)
      dimension DIPZT(max_cal),DIPXT(max_cal),ENERT(max_cal)
      character*2 wover, wr1, wr2, wr3
      character*3 str_only
      character*6 prog
      character*6 padjst
      character*8 unit_r, unit_a, unit_d
      character*8 moltyp
      character*12 wfn
      character*80 atitle
      character*80 title1, title2
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
      data a00,one,two / 0.0d+00 , 1.0d+00 , 2.0d+00 /
      data bohr / 0.529177249d+00 /
      data pi / 3.1415926536d+00 /
      data wave / 219474.63d+00 /
      data rrmax, rrmin / 5.0d+00 , 0.50d+00 /
      data ppmax, ppmin / 4.0d+00 , 0.10d+00 /
c==================================================================
c  For psi 2.0 input
#include <error.h>
      integer frdc, frdi, frdf, fwkclr, fwkadd, frdboo, frdpos
      integer fcnt, frdia1
      integer errcod
c==================================================================
    1 FORMAT(//,2X,' THE ANALYTIC POTENTIAL PROGRAM FOR A TRIATOMIC MOLE
     1CULE'/
     2          2x,'          EXTRACTING NECESSARY ENERGY POINTS'/
     3          2x,'              (Version of 03/29/2007)'/)
    2 FORMAT(A80)                                                       FDG00200
    3 FORMAT(3x,a80/3x,'(title in file11)')
    4 FORMAT(3x,a80/3x,'(title in file81)')
    5 FORMAT(//,2X,' PARAMETERS'/
     1 2X,' PROG      = ',A6/
     2 2X,' WFN       = ',A12/
     3 2X,' UNIT_R    = ',A8/
     4 2X,' UNIT_A    = ',A8/
     5 2X,' UNIT_D    = ',A8/
     6 2X,' MOL_TYPE  = ',A8/
     7 2X,' STR_ONLY  = ',A3/
     8 2X,' IPRNT     = ',I6/
     9 2X,' NATOMS    = ',I6/
     1 2X,' N3N       = ',I6/
     2 2X,' NUM_PTS   = ',I6/
     3 2X,' MAX_DELTA = ',I6/ 
     4 2X,' P_ADJUST  = ',A6/
     5 2X,' R1_MAX    = ',F10.6/
     6 2X,' R1_MIN    = ',F10.6/
     7 2X,' R2_MAX    = ',F10.6/
     8 2X,' R2_MIN    = ',F10.6/
     9 2X,' R3_MAX    = ',F10.6/
     1 2X,' R3_MIN    = ',F10.6/)
    6 FORMAT(//,2X,' THE REFERENCE GEOMETRY IN TAPE11'/
     1 39X,' Cartesian X in a.u.'/
     2 4X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z/)               FDG00380
    7 FORMAT(4F20.10)  
    8 FORMAT(3X,I5,1X,4F16.10)                                          FDG00390
    9 FORMAT(10X,F8.2,3F15.9)
   10 format(//,2x,' Equilibrium bond distances'/
     1          2x,' R_i(e)        = ',f20.10/
     2          2x,' R_j(e)        = ',f20.10/)
   11 FORMAT(//,2X,' THE REFERENCE GEOMETRY IN TAPE11'/
     1 39X,' Cartesian X in angstrom'/
     2 4X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z/)               FDG00380
   12 FORMAT(I5)
   13 format(//,2x,' PERTURBED GEOMETRIES, DIPOLE MOMENTS, AND ENERGIES'
     1 /4X,4H NO.,3X,4H R_i,11X,4H R_j,11X,4H P_k,14X,5H DIPZ,
     2 15X,5H DIPX,15X,7H Energy,14X,7H RELAT./
     3 11X,7H (Bohr),8X,7H (Bohr),8X,9H (radian),9X,8H (Debye),
     4 12X,8H (Debye),12X,10H (hartree),11X,7H (CM-1),/)
   14 format(i5,3f15.10,3f20.12)
   15 format(f7.1,3f15.10,3f20.12)
   16 format(2x,i5,3f15.10,3f20.12,f15.6,4a2)
   17 format(//,2x,' PERTURBED GEOMETRIES, DIPOLE MOMENTS, AND ENERGIES'
     1 /4X,4H NO.,3X,7H R_i-Re,8X,7H R_j-Re,8X,7H P_k-Pe,
     2 11X,5H DIPZ,15X,5H DIPX,15X,7H Energy,14X,7H RELAT./
     3 11X,7H (Bohr),8X,7H (Bohr),8X,9H (radian),9X,8H (Debye),
     4 12X,8H (Debye),12X,10H (hartree),11X,7H (CM-1),/)
   18 format(//,2x,' PERTURBED GEOMETRIES, DIPOLE MOMENTS, AND ENERGIES'
     1 /4X,4H NO.,3X,7H R_i-Re,8X,7H R_j-Re,8X,7H Pi-P_k,
     2 11X,5H DIPZ,15X,5H DIPX,15X,7H Energy,14X,7H RELAT./
     3 11X,7H (Bohr),8X,7H (Bohr),8X,9H (radian),9X,8H (Debye),
     4 12X,8H (Debye),12X,10H (hartree),11X,7H (CM-1),/)
   19 format(//,2x,' PERTURBED GEOMETRIES, DIPOLE MOMENTS, AND ENERGIES'
     1 /4X,4H NO.,3X,7H R_i-Re,8X,7H R_j-Re,8X,7H P_k   ,
     2 11X,5H DIPZ,15X,5H DIPX,15X,7H Energy,14X,7H RELAT./
     3 11X,7H (Bohr),8X,7H (Bohr),8X,9H (radian),9X,8H (Debye),
     4 12X,8H (Debye),12X,10H (hartree),11X,7H (CM-1),/)
   20 format(//,2x,' PERTURBED GEOMETRIES, DIPOLE MOMENTS, AND ENERGIES'
     1 ,' (EXTRACTED)', 
     2 /4X,4H NO.,3X,4H R_i,11X,4H R_j,11X,4H P_k,14X,5H DIPZ,
     3 15X,5H DIPX,15X,7H Energy,14X,7H RELAT./
     4 11X,7H (Bohr),8X,7H (Bohr),8X,9H (radian),9X,8H (Debye),
     5 12X,8H (Debye),12X,10H (hartree),11X,7H (CM-1),/)
   21 format(//,2x,' PERTURBED GEOMETRIES, DIPOLE MOMENTS, AND ENERGIES'
     1 ,' (EXTRACTED)',
     2 /4X,4H NO.,3X,7H R_i-Re,8X,7H R_j-Re,8X,7H P_k-Pe,11X,5H DIPZ,
     3 15X,5H DIPX,15X,7H Energy,14X,7H RELAT./
     4 11X,7H (Bohr),8X,7H (Bohr),8X,9H (radian),9X,8H (Debye),
     5 12X,8H (Debye),12X,10H (hartree),11X,7H (CM-1),/)
   22 format(//,2x,' PERTURBED GEOMETRIES, DIPOLE MOMENTS, AND ENERGIES'
     1 ,' (EXTRACTED)',
     2 /4X,4H NO.,3X,7H R_i-Re,8X,7H R_j-Re,8X,7H Pi-P_k,11X,5H DIPZ,
     3 15X,5H DIPX,15X,7H Energy,14X,7H RELAT./
     4 11X,7H (Bohr),8X,7H (Bohr),8X,9H (radian),9X,8H (Debye),
     5 12X,8H (Debye),12X,10H (hartree),11X,7H (CM-1),/)
   23 format(//,2x,' PERTURBED GEOMETRIES, DIPOLE MOMENTS, AND ENERGIES'
     1 ,' (EXTRACTEDD)',
     2 /4X,4H NO.,3X,7H R_i-Re,8X,7H R_j-Re,8X,7H P_k   ,11X,5H DIPZ,
     3 15X,5H DIPX,15X,7H Energy,14X,7H RELAT./
     4 11X,7H (Bohr),8X,7H (Bohr),8X,9H (radian),9X,8H (Debye),
     5 12X,8H (Debye),12X,10H (hartree),11X,7H (CM-1),/)
   24 format(2x,i5,3x,f15.10,3x,f15.10,3x,f15.10,5x,f20.12)
   26 format(/,2x,' Number of original sets  = ',i5)
   27 format(2x,' Number of deleted sets   = ',i5)
   28 format(2x,' Number of extracted sets = ',i5/)
C                                                                       FDG00570
      INPUT=5                                                           FDG00580
      ITAPE6=6                                                          FDG00590
      ITAP11=11                                                         FDG00600
      ITAP81=81
      ITAP87=87
C                                                                       FDG00620
      call psinit('APPEND')
      CALL TSTART(6)                                                    FDG00550
      call ffile(itap11,'file11',0)
      call ffile(itap81,'file81',0)
      call ffile(itap87,'extr81',0)
      WRITE(6,1)                                                        FDG00630
      write(6,*) '  '
      write(6,*) '  MAXCOR    = ',maxcor
      write(6,*) '  MAX_CAL   = ',max_cal
C
C=============================================================
C   MOL_TYPE
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
      errcod    = fwkadd(5,'EXTR3D_MP0703')
      errcod    = fwkadd(5,'DEFAULT')
      errcod    = frdi('PRINT',iprnt)
      r1max     = rrmax
      errcod    = frdf('R1_MAX',r1max)
      r1min     = rrmin
      errcod    = frdf('R1_MIN',r1min)
      r2max     = rrmax
      errcod    = frdf('R2_MAX',r2max)
      r2min     = rrmin
      errcod    = frdf('R2_MIN',r2min)
      r3max     = ppmax
      errcod    = frdf('R3_MAX',r3max)
      r3min     = ppmin
      errcod    = frdf('R3_MIN',r3min)
C
      errcod    = fwkadd(5,'POTENT3')
      unit_r    = 'ATOMIC  '
      errcod    = frdc('UNIT_R',unit_r)
      unit_a    = 'RADIAN  '
      errcod    = frdc('UNIT_A',unit_a)
      unit_d    = 'ATOMIC  '
      errcod    = frdc('UNIT_D',unit_d)
      prog      = 'PSI   '
      errcod    = frdc('PROGRAM',prog)
      wfn       = 'SCF  '
      errcod    = frdc('WFN',wfn)
      str_only  = 'NO '
      errcod    = frdc('STR_ONLY',str_only)
      max_delta = 5000
      errcod    = frdi('MAX_DELTA',max_delta)
      del_max   = dfloat(max_delta)
      padjst    = 'EQUIL '
      errcod    = frdc('P_ADJUST',padjst)
C
      errcod    = fwkadd(5,'MOLECULE')
      errcod    = frdc('MOL_TYPE',moltyp)
      errcod    = frdia1('RKL',kl,2)
      errcod    = frdia1('RLM',lm,2)
      errcod    = frdia1('AKLM',klm,3)
      errcod    = frdia1('BOND_1',kl,2)
      errcod    = frdia1('BOND_2',lm,2)
      errcod    = frdia1('ANGLE',klm,3)
C
C  Read in information from file11
      REWIND ITAP11                                                     FDG00800
      READ(ITAP11,2) ATITLE
      write(6,*) ' '
      write(6,3) atitle
      natoms = 0
  200 continue
      read(itap11,2,end=201) atitle
      natoms = natoms + 1
      go to 200
  201 continue
C
C  Read in information from file81
      REWIND ITAP81                                                     FDG00800
      READ(ITAP81,2) atitle
      write(6,*) ' '
      write(6,4) atitle
      READ(ITAP81,2) atitle
      READ(ITAP81,2) atitle
      npoints = 0
  202 continue
      read(itap81,2,end=203) atitle
      npoints = npoints + 1
      go to 202
  203 continue
C
      num_pts = npoints
      n3n = natoms*3
      WRITE(6,5) PROG,WFN,UNIT_R,UNIT_A,UNIT_D,MOLTYP,STR_ONLY,
     *           IPRNT,NATOMS,N3N,
     *           num_pts,max_delta,padjst,
     *           r1max,r1min,r2max,r2min,r3max,r3min
C
      if(num_pts.gt.max_cal) then
        write(6,*) '  '
        write(6,*) '  ***************************'
        write(6,*) '  NUM_PTS exceeds MAX_CAL !!!'
        write(6,*) '  NUM_PTS = ',num_pts
        write(6,*) '  MAX_CAL = ',max_cal
        write(6,*) '  ---> STOP !!!'
        write(6,*) '  ***************************'
        write(6,*) '  '
        stop
      end if
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C   Read in the reference geometry from "file11.dat" file
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      REWIND ITAP11                                                     FDG00800
      READ(ITAP11,2) ATITLE  
      if(prog.eq.'PSI   ') then
        WRITE(6,6)    
        DO 101 I=1,NATOMS                           
        READ(ITAP11,7) CHG(I),XOLD(I),YOLD(I),ZOLD(I)
        WRITE(6,8) I,CHG(I),XOLD(I),YOLD(I),ZOLD(I)
  101   CONTINUE                                   
      end if
C
      if(prog.eq.'MOLPRO') then
        WRITE(6,6)
        DO 102 I=1,NATOMS
        READ(ITAP11,9) CHG(I),XOLD(I),YOLD(I),ZOLD(I)
        WRITE(6,8) I,CHG(I),XOLD(I),YOLD(I),ZOLD(I)
  102   CONTINUE
      end if
C
      call dist
      rie = r(kl(1),kl(2))
      rje = r(lm(1),lm(2))
      write(6,10) rie,rje
      call angl
      pke = aglb
C
C  Cartesian coordinates in Angstrom
      WRITE(6,11)
      DO 103 I=1,NATOMS
      xa(i)=xold(i) * bohr
      ya(i)=yold(i) * bohr
      za(i)=zold(i) * bohr
      WRITE(6,8) I,CHG(I),XA(I),YA(I),ZA(I)
  103 CONTINUE
      call length(xa,ya,za)
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C   Read in perturbed geometrical parameters, dipole moments,
C      and total energies from "file81.dat" file
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      REWIND ITAP81                                                     FDG00800
      READ(ITAP81,2) title1
      READ(ITAP81,2) title2
      READ(ITAP81,12) num_temp
C
      if(num_temp.ne.num_pts) then
        write(6,*) '  '
        write(6,*) '  *********************************'
        write(6,*) '  NUM_TEMP not equal to NUM_PTS !!!'
        write(6,*) '  NUM_TEPM = ',num_temp
        write(6,*) '  NUM_PTS  = ',num_pts
        write(6,*) '  ---> STOP !!!'
        write(6,*) '  *********************************'
        write(6,*) '  '
        stop
      end if
C
C===========================================================
C  Extract energy points less than MAX_DELTA (in cm-1)
C===========================================================
C  For standard units
C  (bond lengths in atomic unit and bond angles in radians)
C===========================================================
      write(6,13)
      icount = 0
      ncount = 0
      do 104 i=1,num_pts
        if(prog.eq.'PSI   ') then
          read(itap81,14) nn,rri,rrj,ppk,dipz,dipx,evalu
        end if
        if(prog.eq.'MOLPRO') then
          read(itap81,15) xnn,rri,rrj,ppk,dipz,dipx,evalu
          nn = int(xnn)
        end if

      iflag = 0
      if(i.eq.1) then
        eref = evalu
      end if
      dele = evalu - eref
      delw = dele * wave
      wover = '  '
      if(delw.gt.del_max) then
        wover = ' *'
        iflag = iflag + 1
        ncount = ncount + 1
      end if

      wr1     = '  '
      if(rri .gt. r1max) then
        iflag = iflag + 1
        ncount = ncount + 1
        wr1 = ' +'
      end if
      if(rri .lt. r1min) then
        iflag = iflag + 1
        ncount = ncount + 1
        wr1 = ' -'
      end if

      wr2     = '  '
      dr2_max = r2max - rrj
      dr2_min = r2min - rrj
      if(rrj .gt. r2max) then
        iflag = iflag + 1
        ncount = ncount + 1
        wr2 = ' +'
      end if
      if(rrj .lt. r2min) then
        iflag = iflag + 1
        ncount = ncount + 1
        wr2 = ' -'
      end if

      wr3     = '  '
      if(ppk .gt. r3max) then
        iflag = iflag + 1
        ncount = ncount + 1
        wr3 = ' +'
      end if
      if(ppk .lt. r3min) then
        iflag = iflag + 1
        ncount = ncount + 1
        wr3 = ' -'
      end if

      write(6,16) nn,rri,rrj,ppk,dipz,dipx,evalu,delw,wover,
     1            wr1,wr2,wr3
      ri0(i)   = rri
      rj0(i)   = rrj
      pk0(i)   = ppk
      dipz0(i) = dipz
      dipx0(i) = dipx
      ener0(i) = evalu

      if(iflag.ne.0) go to 104

      icount = icount + 1
      rit(icount)   = rri
      rjt(icount)   = rrj
      pkt(icount)   = ppk
      dipzt(icount) = dipz
      dipxt(icount) = dipx
      enert(icount) = evalu
  104 continue
      if(iprnt.ge.3) then
        write(6,*) '  '
        write(6,*) '  Number of deleted points = ',ncount
        write(6,*) '  '
      end if
C
C=================================================
C Case (1) : equilibrium angle reference
C    a) r_{i} = r_{i} - r_{i}^{e}
C    b) r_{j} = r_{j} - r_{j}^{e}
C    c) p_{k} = p_{k} - p_{k}^{e}
C=================================================
      if(padjst.eq.'EQUIL ') then
        write(6,17)
        do 105 i=1,num_pts
        ridif = ri0(i)-rie
        rjdif = rj0(i)-rje
        pkdif = pk0(i)-pke
        dipz  = dipz0(i)
        dipx  = dipx0(i)
        evalu = ener0(i)
        del   = evalu - eref
        delw  = del * wave
        wover = '  '
        wr1   = '  '
        wr2   = '  '
        wr3   = '  '
        if(delw.gt.del_max) then
          wover = ' *'
        end if
        write(6,16) i,ridif,rjdif,pkdif,dipz,dipx,evalu,delw,wover,
     1              wr1,wr2,wr3
  105   continue
      end if
C
C=================================================
C Case (2) : linear angle reference
C    a) r_{i} = r_{i} - r_{i}^{e}
C    b) r_{j} = r_{j} - r_{j}^{e}
C    c) p_{k} =   pi  - p_{k}
C=================================================
      if(padjst.eq.'LINEAR') then
        write(6,18)
        do 106 i=1,num_pts
        ridif = ri0(i)-rie
        rjdif = rj0(i)-rje
        pkdif = pi-pk0(i)
        dipz  = dipz0(i)
        dipx  = dipx0(i)
        evalu = ener0(i)
        del   = evalu - eref
        delw  = del * wave
        wover = '  '
        wr1   = '  '
        wr2   = '  '
        wr3   = '  '
        if(delw.gt.del_max) then
          wover = ' *'
        end if
        write(6,16) i,ridif,rjdif,pkdif,dipz,dipx,evalu,delw,wover
  106   continue
      end if
C
C=================================================
C Case (3) : no angle reference
C    a) r_{i} = r_{i} - r_{i}^{e}
C    b) r_{j} = r_{j} - r_{j}^{e}
C    c) p_{k} = p_{k}
C=================================================
      if(padjst.eq.'NONE  ') then
        write(6,19)
        do 107 i=1,num_pts
        ridif = ri0(i)-rie
        rjdif = rj0(i)-rje
        pkdif = pk0(i)
        dipz  = dipz0(i)
        dipx  = dipx0(i)
        evalu = ener0(i)
        del   = evalu - eref
        delw  = del * wave
        wover = '  '
        wr1   = '  '
        wr2   = '  '
        wr3   = '  '
        if(delw.gt.del_max) then
          wover = ' *'
        end if
        write(6,16) i,ridif,rjdif,pkdif,dipz,dipx,evalu,delw,wover,
     1              wr1,wr2,wr3
  107   continue
      end if
C
C=================================================
C  For extracted energy points less than MAX_DELTA
C=================================================
      write(6,20)
      do 110 i=1,icount
      rri   = rit(i)
      rrj   = rjt(i)
      ppk   = pkt(i)
      dipz  = dipzt(i)
      dipx  = dipxt(i)
      evalu = enert(i)
      del   = evalu - eref
      delw  = del * wave
      wover = '  '
      wr1   = '  '
      wr2   = '  '
      wr3   = '  '
      if(delw.gt.del_max) then
        wover = ' *'
      end if
      write(6,16) i,rri,rrj,ppk,dipz,dipx,evalu,delw,wover,
     1              wr1,wr2,wr3
  110 continue
C
C=================================================
C Case (1) : equilibrium angle reference
C    a) r_{i} = r_{i} - r_{i}^{e}
C    b) r_{j} = r_{j} - r_{j}^{e}
C    c) p_{k} = p_{k} - p_{k}^{e}
C=================================================
      if(padjst.eq.'EQUIL ') then
        write(6,21)
        do 111 i=1,icount
        ridif = rit(i)-rie
        rjdif = rjt(i)-rje
        pkdif = pkt(i)-pke
        dipz  = dipzt(i)
        dipx  = dipxt(i)
        evalu = enert(i)
        del   = evalu - eref
        delw  = del * wave
        wover = '  '
        wr1   = '  '
        wr2   = '  '
        wr3   = '  '
        if(delw.gt.del_max) then
          wover = ' *'
        end if
        write(6,16) i,ridif,rjdif,pkdif,dipz,dipx,evalu,delw,wover,
     1              wr1,wr2,wr3
  111   continue
        write(6,*) '  '
      end if
C
C=================================================
C Case (2) : linear angle reference
C    a) r_{i} = r_{i} - r_{i}^{e}
C    b) r_{j} = r_{j} - r_{j}^{e}
C    c) p_{k} =   pi  - p_{k}
C=================================================
      if(padjst.eq.'LINEAR') then
        write(6,22)
        do 112 i=1,icount
        ridif = rit(i)-rie
        rjdif = rjt(i)-rje
        pkdif = pi-pkt(i)
        dipz  = dipzt(i)
        dipx  = dipxt(i)
        evalu = enert(i)
        del   = evalu - eref
        delw  = del * wave
        wover = '  '
        wr1   = '  '
        wr2   = '  '
        wr3   = '  '
        if(delw.gt.del_max) then
          wover = ' *'
        end if
        write(6,16) i,ridif,rjdif,pkdif,dipz,dipx,evalu,delw,wover,
     1              wr1,wr2,wr3
  112   continue
        write(6,*) '  '
      end if
C
C=================================================
C Case (3) : no angle reference
C    a) r_{i} = r_{i} - r_{i}^{e}
C    b) r_{j} = r_{j} - r_{j}^{e}
C    c) p_{k} = p_{k}
C=================================================
      if(padjst.eq.'NONE  ') then
        write(6,23)
        do 113 i=1,icount
        ridif = rit(i)-rie
        rjdif = rjt(i)-rje
        pkdif = pkt(i)
        dipz  = dipzt(i)
        dipx  = dipxt(i)
        evalu = enert(i)
        del   = evalu - eref
        delw  = del * wave
        wover = '  '
        wr1   = '  '
        wr2   = '  '
        wr3   = '  '
        if(delw.gt.del_max) then
          wover = ' *'
        end if
        write(6,16) i,ridif,rjdif,pkdif,dipz,dipx,evalu,delw,wover,
     1              wr1,wr2,wr3
  113   continue
        write(6,*) '  '
      end if
C
C##################################################################
C  Store the extracted information in "extract81.dat (file87)" file
C##################################################################
      rewind itap87
      write(itap87,2) title1
      write(itap87,2) title2
      write(itap87,12) icount
      do 117 i=1,icount
      rri   = rit(i)
      rrj   = rjt(i)
      ppk   = pkt(i)
      dipz  = dipzt(i)
      dipx  = dipxt(i)
      evalu = enert(i)
      if(prog.eq.'PSI   ') then
        write(itap87,14) i,rri,rrj,ppk,dipz,dipx,evalu
      end if
      if(prog.eq.'MOLPRO') then
        xii = i
        write(itap87,15) xii,rri,rrj,ppk,dipz,dipx,evalu
      end if
  117 continue
      rewind itap87
      write(6,26) num_pts
      write(6,27) ncount
      write(6,28) icount
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
