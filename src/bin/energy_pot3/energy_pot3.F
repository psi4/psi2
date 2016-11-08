      subroutine fentry(cc,ia,maxcor)
c*    program energy_pot3
C**********************************************************************
c   This program extracts an energy of the given wave function        *
c        and add to the file61                                        *
C**********************************************************************
C*  Last updated on March 05, 2003 by Dr. Yukio Yamaguchi             *
C**********************************************************************
C                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                    
      dimension cc(maxcor),ia(maxcor*2)                                
      parameter (max_cal=5000)
      character*8 moltyp, unit
      character*8 wfn, opentyp
      character*80 title1, title2
      dimension ipts(max_cal),energ(max_cal)
      dimension ri(max_cal),rj(max_cal),pk(max_cal)
      common/enr101/energy
      common/enr102/wfn
      common/com101/natoms
      common/com102/iprnt
      common/com103/itap10
      common/com105/apir,aglb,ang
      common/com201/chg(10),w(10)
      common/com202/x(10),y(10),z(10)
      common/com205/r(10,10)
      common/com206/kl(2),lm(2),klm(3)
c==================================================================
c  For psi 2.0 input
#include <error.h>
      integer frdc, frdi, frdf, fwkclr, fwkadd, frdboo, frdpos
      integer fcnt, frdia1
      integer errcod
c==================================================================
    1 FORMAT(/,2X,' A program to extract the total energy '/
     1          2X,'         (VERSION OF 03/05/03)          '/)        
    2 format(A80)
    3 FORMAT(I5)
    4 FORMAT(//,2X,' PARAMETERS'/
     1 2X,' WFN       = ',A8/
     2 2X,' OPEN_TYPE = ',A8/
     3 2X,' MOL_TYPE  = ',A8/
     4 2X,' UNIT      = ',A8/
     5 2X,' IPRNT     = ',I6/
     6 2X,' NATOMS    = ',I6/
     7 2X,' N3N       = ',I6/
     8 2x,' NPERTS    = ',I6/
     9 2x,' NPERTP    = ',I6/)
    5 format(//,2x,' Current Bond distances'/
     1          2x,' R_i(e)        = ',f20.10/
     2          2x,' R_j(e)        = ',f20.10)
    6 format(//,2x,' PERTURBED GEOMETRIES AND ENERGIES'/
     1 4X,4H NO.,6X,4H R_i,14X,4H R_j,14X,4H P_k,19X,7H Energy,/
     2 14X,7H (Bohr),11X,7H (Bohr),11X,9H (radian),14X,10H (hartree),/)
    7 format(i5,5x,3f15.10,5x,f20.12)
    8 format(2x,i5,3x,f15.10,3x,f15.10,3x,f15.10,4x,f25.16,2x,f15.6)
c
      call psinit('APPEND')
      call tstart(6)
      input=5
      itap10=10
      itap61=61
      call ffile(itap10,'ccteng',1)
      call ffile(itap61,'file61',1)
c
      write(6,1)
      write(6,*) '  maxcor    = ',maxcor
      write(6,*) '  max_cal   = ',max_cal
c
      errcod  = fwkclr()
      errcod  = fwkadd(5,'ENERGY_POT3')
      errcod  = fwkadd(5,'DEFAULT')
      errcod  = frdi('PRINT',iprnt)
      unit    = 'ATOMIC  '
      errcod  = frdc('UNIT',unit)
      wfn     = 'SCF     '
      errcod  = frdc('WFN',wfn)
      opentyp = 'NONE    '
      errcod  = frdc('OPENTYPE',opentyp)
c
      errcod  = fwkadd(5,'MOLECULE')
      errcod  = frdc('MOL_TYPE',moltyp)
      errcod  = frdia1('RKL',kl,2)
      errcod  = frdia1('RLM',lm,2)
      errcod  = frdia1('AKLM',klm,3)
      errcod  = frdia1('BOND_1',kl,2)
      errcod  = frdia1('BOND_2',lm,2)
      errcod  = frdia1('ANGLE',klm,3)
 
C   Find the number of energy points stored in file61
      call izero(ipts,max_cal)
      call zero(ri,max_cal)
      call zero(rj,max_cal)
      call zero(pk,max_cal)
      call zero(energ,max_cal)
      rewind itap61
      read(itap61,2,err=200) title1
      read(itap61,2,err=200) title2
      read(itap61,3,err=200) nperts
      if(nperts.ne.0) go to 201
  200 continue
      nperts = 0
      title1 = 'This is a provided label (1).'
      title2 = 'This is a provided label (2).'
  201 continue
      npertp=nperts+1
C
C   Read in geometry in file30
      call read30
      n3n = natoms*3
      write(6,4) wfn,opentyp,moltyp,unit,
     *           iprnt,natoms,n3n,nperts,npertp
c
      call dist
      rie = r(kl(1),kl(2))
      rje = r(lm(1),lm(2))
      write(6,5) rie,rje
      call angl
C
      if(nperts.eq.0) go to 205
      write(6,6)
      do 101 i=1,nperts
      read(itap61,7) ii,rri,rrj,ppk,valu
      write(6,8) ii,rri,rrj,ppk,valu
      ipts(i)=ii
      ri(i)=rri
      rj(i)=rrj
      pk(i)=ppk
      energ(i)=valu
  101 continue
  205 continue
c
c====================================
c=  SCF wave functions
c====================================
      if(wfn.eq.'SCF     ') then
        call scfwfn
      end if
      if(wfn.eq.'TCSCF   ') then
        call scfwfn
      end if
      if(wfn.eq.'TCSCFX  ') then
        call scfwfn
      end if
c
c====================================
c=  CI wave functions
c====================================
      if(wfn.eq.'CI      ') then
        call ciwfn
      end if
      if(wfn.eq.'SOCI    ') then
        call ciwfn
      end if
      if(wfn.eq.'CASSCF  ') then
        call ciwfn
      end if
      if(wfn.eq.'SACASSCF') then
        call ciwfn
      end if
c
c====================================
c=  CCSD and CCSD(T) wave functions
c====================================
      if(wfn.eq.'CCSD    ') then
        if(opentyp.eq.'NONE    ') then 
          call clccwf
        else
          call osccwf
        end if
      end if
      if(wfn.eq.'CCSD_T  ') then
        if(opentyp.eq.'NONE    ') then 
          call clccwf
        else
          call osccwf
        end if
      end if
c
c  Add a new geometry and energy in file61
      ipts(npertp)=npertp
      ri(npertp)=rie
      rj(npertp)=rje
      pk(npertp)=aglb
      energ(npertp)=energy
      write(6,*) ' '
      write(6,8) ipts(npertp),ri(npertp),rj(npertp),pk(npertp),
     *           energ(npertp)
      write(6,*) ' '
C
      rewind itap61
      write(itap61,2) title1
      write(itap61,2) title2
      write(itap61,3) npertp
      do 102 i=1,npertp
      write(itap61,7) ipts(i),ri(i),rj(i),pk(i),energ(i)
  102 continue
      rewind itap61
c
      CALL TSTOP(6)                                                 
      STOP                                                         
      END                                                         
      subroutine read30
      implicit real*8 (a-h,o-z)
      dimension i30(200),coord(3,10)
      common/com101/natoms
      common/com201/chg(10),w(10)
      common/com202/x(10),y(10),z(10)
    1 FORMAT(//,2X,' PARAMETERS IN FILE30'/
     1 2X,' MPOINT    = ',I5/  
     2 2X,' MCONST    = ',I5/ 
     3 2X,' NCALCS    = ',I5/
     4 2X,' NATOMS    = ',I5/ 
     5 2X,' N3N       = ',I5/
     6 2X,' IGEOP     = ',I5/)
    2 FORMAT(//,2X,' GEOMETRY IN TAPE30'/       
     1 39X,' CURRENT X in a.u.'/   
     2 4X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z/)               FDG00460
    3 FORMAT(3X,I5,1X,4F16.10)                                          FDG00470
C
C   READ IN GEOMETRY FROM TAPE30                                        FDG00690
      itap30 = 30
      CALL RFILE(ITAP30)                                                FDG00700
      CALL WREADW(ITAP30,I30,200,101,JUNK)                              FDG00710
      MPOINT=I30(2)                                                     FDG00720
      MCONST=I30(3)                                                     FDG00730
      NCALCS=I30(5)                                                     FDG00740
      NATOMS=I30(19)                                                     FDG00750
      N3N=NATOMS*3
      IGEOP=100+MCONST+MPOINT+NCALCS                                    FDG00770
C                                                                       FDG00780
      WRITE(6,1) MPOINT,MCONST,NCALCS,NATOMS,N3N,IGEOP
C                                                                       FDG01130
C   READ IN NUCLEAR CHARGE                                              MAS11050
      JUNK=101+MCONST                                                   MAS11020
      CALL WREADW(ITAP30,I30,MPOINT,JUNK,JUNK)                          MAS11030
      CALL WREADW(ITAP30,CHG,NATOMS*2,I30(1),JUNK)
C                                                                       FDG01130
C   READ IN THE PRESENT GEOMETRY                                        FDG01260
      CALL WREADW(ITAP30,LOCCAL,1,IGEOP,JUNK)                           FDG01270
      IGEOP=LOCCAL+80                                                   FDG01280
      CALL WREADW(ITAP30,COORD,N3N*2,IGEOP,JUNK)                        FDG01290
      WRITE(6,2)
      do 101 i=1,natoms
      x(i)=coord(1,i)
      y(i)=coord(2,i)
      z(i)=coord(3,i)
      WRITE(6,3) I,CHG(I),x(i),y(i),z(i)
  101 continue
      CALL RCLOSE(ITAP30,3)
C
      return
      end
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
      subroutine scfwfn
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension i30(200),a30(100),xx(300)
      equivalence (i30,a30)
      common/enr101/energy
      common/com102/iprnt
c
c=======================================================================
c   Obtain SCF information from master file 30
c=======================================================================
      itap30=30
      call rfile(itap30)
      CALL WREADW(ITAP30,I30,200,101,JUNK)                              MAS01580
      IEND=I30(1)                                                       MAS01590
      MPOINT=I30(2)                                                     MAS01600
      MCONST=I30(3)                                                     MAS01610
      MCALCS=I30(4)                                                     MAS01620
      NCALCS=I30(5)                                                     MAS01630
      NBFSO=I30(18)                                                     MAS01640
      NATOM=I30(19)                                                     MAS01650
      NBFAO=I30(22)                                                     MAS01660
      MXCOEF=I30(42)                                                    MAS01680
      IOPEN=I30(43)                                                     MAS01690
c
C   READ IN SCF ENERGY                                                  MAS11240
      JUNK=101+MCONST+MPOINT+NCALCS-1                                   MAS11160
      CALL WREADW(ITAP30,LOCCAL,1,JUNK,JUNK)                            MAS11170
      if(iprnt.gt.5) then
        write(6,*) '  loccal 1 = ',loccal
      end if
      LOCCAL=LOCCAL+80+natom*6
      if(iprnt.gt.5) then
        write(6,*) '  loccal 2 = ',loccal
      end if
      CALL WREADW(ITAP30,I30,20,LOCCAL,LOCCAL)                          MAS11250
      ENUC=A30(1)                                                       MAS11260
      ESCF=A30(2)                                                       MAS11270
C
      if(iprnt.ge.3) then
        write(6,*) '   '
        write(6,*) '   '
        write(6,*) '  enuc   = ',enuc
        write(6,*) '  escf   = ',escf
        write(6,*) '   '
      end if
      energy=escf
c
      call rclose(itap30,3)
      return
      end
      subroutine ciwfn
      IMPLICIT REAL*8 (A-H,O-Z)
      integer intlbl(26),scflbl(26)
      common/enr101/energy
    1 format(26a3)
c
c=======================================================================
c   Obtain CI information from CI vector file (54)
c=======================================================================
      itap54=54
      call rfile(itap54)
      call dcdlbl(itap54,intlbl,scflbl,nrfwlk,ix,nwks,eguess,eci,czero)
C                  
c     write(6,*) '  Label from Integral'
c     write(6,1) (intlbl(i),i=1,26)
c     write(6,*) '  Label from Scf'
c     write(6,1) (scflbl(i),i=1,26)
      write(6,*) '  nrfwlk = ',nrfwlk
      write(6,*) '  ix     = ',ix
      write(6,*) '  nwks   = ',nwks
      write(6,*) '  eguess = ',eguess
      write(6,*) '  eci    = ',eci
      write(6,*) '  czero  = ',czero
      write(6,*) '  '
      energy=eci
c
      call rclose(itap54,3)
      return
      end
      subroutine dcdlbl(itape,na,nb,nc,nd,ne,rf,rg,rh)
      implicit real*8 (a-h,o-z)
      dimension na(26),nb(26),n(112),rn(10)
      equivalence (n(57),rn(1))
c
      call sread(itape,n,112)
      do 101 i=1,26
      na(i)=n(i)
      nb(i)=n(i+26)
  101 continue
      nc=n(53)
      nd=n(54)
      ne=n(55)
      rf=rn(1)
      rg=rn(2)
      rh=rn(3)
c
      return
      end
      subroutine clccwf
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension i30(200),a30(100)
      equivalence (i30,a30)
      common/enr101/energy
c
c=======================================================================
c   Obtain SCF information from master file 30
c=======================================================================
      itap30=30
      call rfile(itap30)
      CALL WREADW(ITAP30,I30,200,101,JUNK)                              MAS01580
      IEND=I30(1)                                                       MAS01590
      MPOINT=I30(2)                                                     MAS01600
      MCONST=I30(3)                                                     MAS01610
      MCALCS=I30(4)                                                     MAS01620
      NCALCS=I30(5)                                                     MAS01630
      NBFSO=I30(18)                                                     MAS01640
      NATOM=I30(19)                                                     MAS01650
      NBFAO=I30(22)                                                     MAS01660
      MXCOEF=I30(42)                                                    MAS01680
      IOPEN=I30(43)                                                     MAS01690
      call rclose(itap30,3)
c
      itap79=79
      call rfile(itap79)
      nt=i30(18)
      nt2=nt*(nt+1)/2
      iadr=1+intowp(nt2)+intowp(nt2)
      write(6,*) '  nt     = ',nt
      write(6,*) '  nt2    = ',nt2
      write(6,*) '  iadr   = ',iadr
      call wreadw(itap79,etot69,intowp(1),iadr,junk)
      write(6,*) '  etot69 = ',etot69
      write(6,*) '  '
      energy=etot69
c
      call rclose(itap79,3)
      return
      end
      subroutine osccwf
      IMPLICIT REAL*8 (A-H,O-Z)
      character*8 wfn
      common/enr101/energy
      common/enr102/wfn
    1 format(4x,f25.16)
c
      itap10=10
      rewind itap10
c
      if(wfn.eq.'ccsd    '.or.wfn.eq.'CCSD    ') then
        read(itap10,1) escf
        read(itap10,1) eccsd
        write(6,*) '  escf   = ',escf
        write(6,*) '  eccsd  = ',eccsd
        write(6,*) '  '
        energy=eccsd
      end if
c
      if(wfn.eq.'ccsd_t  '.or.wfn.eq.'CCSD_T  ') then
        read(itap10,1) escf
        read(itap10,1) eccsd
        read(itap10,1) eccsdt
        write(6,*) '  escf   = ',escf
        write(6,*) '  eccsd  = ',eccsd
        write(6,*) '  eccsdt = ',eccsdt
        write(6,*) '  '
        energy=eccsdt
      end if
      rewind itap10
c
      return
      end
