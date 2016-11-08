      program grdtol
c
c	written by george -croag- vacek
c	Thu Sep  1 13:52:22 EDT 1994
c
c	some ideas and routines may have been taken from
c	bmat and-or gradconv
c
c	i've used mostly double precision throughout, but
c	single precision should be enough.

c
c	notes to croag
c
c	rip off bmat intcoord subroutines
c

      implicit integer (a-z)
#include <error.h>
      integer prcntr,prdup
      integer i11,iout,length,natom,maxatm
      integer frstgr,minima,ip,noint
      integer errcod,frdboo,frdi,int1,int2,i,j
      real*8 cgrmsi,cgmaxi,cgrmsc,cgmaxc,cdmaxc,cdrmsc,enrg0,enrg1
      character*1 garbage
c	current natom limit is 255
      parameter (maxatm=255)
      real*8 x0(maxatm,3),g0(maxatm,3),grms0,gmax0,cmxnrg
      real*8 x1(maxatm,3),g1(maxatm,3),grms1,gmax1
      real*8 dtemp,drms,dmax

c	Cartesian gradients are in file11
c	internal gradients will be determined
      iout = 6
      i11 = 11

c	note that if you call unit11 by the name file11.dat, the
c	program will open a file11.dat.d which is not right.
c	the 0 flag will open to the top of file11 rather than 
c	the bottom as a 1 flag would

      call psinit('APPEND')
      call ffile(i11,' ',0)

      call tstart(iout)

c
c  print control options
c
c have three modes for now, DEBUG, VERBOSE, and BRIEF
c       BRIEF  prints almost nothing
c       VERBOSE  prints some info on the various values 
c       DEBUG  prints everything I can find
c
      errcod = prcntr('ON BRIEF')
      errcod = prcntr('OFF VERBOSE DEBUG')
      errcod = prdup()
      if(prcntr('IS_ON DEBUG').eq.1) errcod = prcntr('ON VERBOSE')
      if(prcntr('IS_ON VERBOSE').eq.1) errcod = prcntr('OFF BRIEF')

      write(iout,*)
      write(iout,10) 'Convergence Criteria Test Program'
      write(iout,10) ' Written by George -Croag- Vacek '
      write(iout,10) '           PSI2.0.8'
      write(iout,*)
 10   format(2x,t20,a)

c	set the defaults
c	read any input information and override defaults if necessary
c	crit.grad-or-disp.rms-or-max.internal-or-Cartesian
c	these could be adjusted some for wfn=, etc.

c	default assumes a minimum is searched for, energy should go down
c	default is to do only a Cartesian convergence check

      minima=1
      errcod = frdboo('MINIMUM',minima)
      noint=1
      errcod = frdboo('NO_INT',noint)
      
      if(prcntr('IS_ON DEBUG').eq.1) then
        write(iout,*)
        write(iout,*) ' MINIMUM = ',minima
        write(iout,*) ' NO_INT = ',noint
      endif

c	defaults:
c	  Energy = 10-10 
c	  Internal RMS GRADIENT = 10-6
c	  Cartesian RMS GRADIENT = 10-6
C	  MAXIMUM GRADIENT = 1.5 * RMS GRADIENT
C         RMS DISPLACEMENT = 4.0 * Cartesian RMS GRADIENT
C         MAX DISPLACEMENT = 1.5 * RMS DISPLACEMENT

      int1 = 10
      errcod = frdi('MAX_NRG',int1)
      cmxnrg = 1.0/10**int1

      int1 = 6
      int2 = 6
      errcod = frdi('RMS_GRAD_INT',int1)
      errcod = frdi('RMS_GRAD_CART',int2)
      cgrmsi = 1.0/10**int1
      cgrmsc = 1.0/10**int2

      int1 = 0
      int2 = 0
      errcod = frdi('MAX_GRAD_INT',int1)
      errcod = frdi('MAX_GRAD_CART',int2)

      if (int1.ne.0) then
	cgmaxi = 1.0/10**int1
      else
	cgmaxi = 1.5* cgrmsi
      endif

      if (int2.ne.0) then
	cgmaxc = 1.0/10**int2
      else
	cgmaxc = 1.5* cgrmsc
      endif

      int1 = 0
      int2 = 0
      errcod = frdi('RMS_DISP_CART',int1)
      errcod = frdi('MAX_DISP_CART',int2)

      if (int1.ne.0) then
        cdrmsc = 1.0/10**int1
      else
        cdrmsc = 4.0* cgrmsc
      endif

      if (int2.ne.0) then
        cdmaxc = 1.0/10**int2
      else
        cdmaxc = 1.5* cdrmsc
      endif

      write(iout,*)
      write(iout,10) 'GradTol will use the following criteria:'
      write(iout,*)
      write(iout,20) '  rms Cartesian grad = ',cgrmsc
      write(iout,20) '  max Cartesian grad = ',cgmaxc
      write(iout,20) '  rms Cartesian disp = ',cdrmsc
      write(iout,20) '  max Cartesian disp = ',cdmaxc
      write(iout,20) '  rms internal grad  = ',cgrmsi
      write(iout,20) '  max internal grad  = ',cgmaxi
      write(iout,*)
 20   format(2x,t20,a,e10.3)

c	read geometry and cartesian gradients from file11
c	we need to get to the bottom (reading a character at a time) and backup
c	length is number of lines in file11

      length = 0
100   read(i11,'(1A1)',end=110) garbage
	length = length +1
	goto 100
110   continue

c	since there is a coord and grad for each atom, we can backup 2
c	lines at a time.  backup until we find data.
c	one could count the number of backup iterations to get natom as well.

      backspace i11
      backspace i11
      read(i11,'(i5)')natom
      if(natom.eq.0) goto 110
      write(iout,*) '  Number of atoms is ', natom

c	check that this is not too many atoms
      if(natom.gt.maxatm) then
	write(iout,*)
	write(iout,*) ' GradTol cannot currently handle that many atoms.'
	write(iout,*)' Please adjust the program (1 line) and recompile.'
#if defined(NO_EXIT_ROUTINE)
          stop 111
#else
          call exit(111)
#endif
      endif

c	and we need to backup once more, to get to the start of that
c	same line of data

      backspace i11

c	lets check if more than one gradient point is in file11.
c	keep option to override

      if(length.eq.2*natom+2) then
	frstgr = 1
      else
	frstgr = 0
      endif
      errcod = frdboo('FIRST_RUN',frstgr)

c	check printing demands
      ip = 0
      ip = prcntr('IS_ON DEBUG')
      if(ip.eq.1) then
        write(iout,*)
        write(iout,*) '  Beginning to read file11 data. '
      endif

c   	if it exists, backup more, and read the previous gradient point data

      if(frstgr.eq.0) then
	do 120 i=1,2*natom+2
	  backspace i11
120	continue

c	read the energy of the previous grad point

	enrg0 = 0.0d+00
        read(i11,30) enrg0
30	format(5x,f20.10)
      	if(ip.eq.1) then
          write(iout,*)
          write(iout,35) ' Energy of the previous point is ',enrg0
      	endif
35      format(2x,a,f20.10)

c	read the Cart coord of the previous grad point
	do 130 i=1,natom
	  x0(i,1) =0.0d+00
	  x0(i,2) =0.0d+00
	  x0(i,3) =0.0d+00
	  read(i11,40) (x0(i,j), j=1,3)
130	continue
40	format(20x,3f20.10)
      	if(ip.eq.1) then
          write(iout,*)
          write(iout,*) '  Geometry of the previous point is:'
	  do 135 i=1,natom
	    write(iout,40) (x0(i,j), j=1,3)
135	  continue
      	endif

c	read the Cart grad of the previous grad point
	do 140 i=1,natom
	  g0(i,1) =0.0d+00
	  g0(i,2) =0.0d+00
	  g0(i,3) =0.0d+00
	  read(i11,40) (g0(i,j), j=1,3)
140	continue
      	if(ip.eq.1) then
          write(iout,*)
          write(iout,*) '  Gradient of the previous point is:'
	  do 145 i=1,natom
	    write(iout,40) (g0(i,j), j=1,3)
145	  continue
	endif

c	skip a line for the title
	read(i11,*)
      endif
      
c	read the data for the current gradient point

      enrg1 = 0.0d+00
      read(i11,30) enrg1
      if(ip.eq.1) then
        write(iout,*)
        write(iout,35) ' Energy of the current point is ',enrg1
      endif

c       read the Cart coord of the current grad point
      do 150 i=1,natom
	x1(i,1) =0.0d+00
	x1(i,2) =0.0d+00
	x1(i,3) =0.0d+00
        read(i11,40) (x1(i,j), j=1,3)
150   continue
      if(ip.eq.1) then
        write(iout,*)
        write(iout,*) '  Geometry of the current point is:'
        do 155 i=1,natom
          write(iout,40) (x1(i,j), j=1,3)
155     continue
      endif

c       read the Cart grad of the current grad point
      do 160 i=1,natom
	g1(i,1) =0.0d+00
	g1(i,2) =0.0d+00
	g1(i,3) =0.0d+00
        read(i11,40) (g1(i,j), j=1,3)
160   continue
      if(ip.eq.1) then
        write(iout,*)
        write(iout,*) '  Gradient of the current point is:'
        do 165 i=1,natom
          write(iout,40) (g1(i,j), j=1,3)
165     continue
      endif

c	ok, that is all the data to be read, so now there is
c	only the analysis left

c	if it is the first gradient, it should still be somewhat
c	reasonable - test it.

      if(frstgr.eq.1) then
	enrg0 = 0.0d+00
	grms0 = 0.75d+00
	gmax0 = 1.0d+00
      else
	grms0 = 0.0d+00
	gmax0 = 0.0d+00
	do 170 i = 1, natom
	  do 180 j=1,3
	    if(dabs(g0(i,j)).gt.gmax0) gmax0 = dabs(g0(i,j))
	    grms0 = grms0 + g0(i,j)*g0(i,j)
180	  continue
170	continue
	grms0 = dsqrt(grms0/(natom*3))
      endif

      grms1 = 0.0d+00
      gmax1 = 0.0d+00
      do 190 i=1,natom
	do 200 j=1,3
	  if(dabs(g1(i,j)).gt.gmax1) gmax1 = dabs(g1(i,j))
	  grms1 = grms1 + g1(i,j)*g1(i,j)
200	continue
190   continue
      grms1 = dsqrt(grms1/(natom*3))

      if(frstgr.eq.0) then
	write(iout,*)
        write(iout,10) 'For the previous geometry point: '
        write(iout,*)
        write(iout,20) '  rms Cartesian grad = ',grms0
        write(iout,20) '  max Cartesian grad = ',gmax0
      endif

      write(iout,*)
      write(iout,10) 'For the current geometry point: '
      write(iout,*)
      write(iout,20) '  rms Cartesian grad = ',grms1
      write(iout,20) '  max Cartesian grad = ',gmax1

c	check that the results have improved

      if(minima.eq.1) then
	if((enrg0 - enrg1).lt.-cmxnrg) then
	  write(iout,*)
	  write(iout,*) '  Warning:  You may have an error.'
	  write(iout,*) '  You are looking for a minimum, but your'
	  write(iout,*) '  total energy is not decreasing.'
	else if(ip.eq.1)   then
          write(iout,*)
          write(iout,*) '  Energy looks ok.'
	endif
	if((grms1.gt.grms0+cgrmsc).or.(gmax1.gt.gmax0+cgmaxc)) then
	  write(iout,*)
          write(iout,*) '  Calculation aborted.'
	  write(iout,*) '  You are looking for a minimum, but your'
	  write(iout,*) '  gradients show an unreasonable increase.'
	  call tstop(iout)
	  call psdone()
#if defined(NO_EXIT_ROUTINE)
          stop 111
#else
          call exit(111)
#endif
	else if(ip.eq.1) then
          write(iout,*)
          write(iout,*) '  Cartesian gradients look ok.'
	endif
      endif

c	do a check -kind of- to see if it is the same geometry point
c	twice in a row

      if(((grms1.eq.grms0).and.(gmax1.eq.gmax0)).and.(enrg1.eq.enrg0))
     &  then
	write(iout,*)
	write(iout,*) '  Calculation aborted.'
	write(iout,*) '  Why would someone want to do the same '
	write(iout,*) '  geometry point twice? '
	call tstop(iout)
	call psdone()
#if defined(NO_EXIT_ROUTINE)
        stop 111
#else
        call exit(111)
#endif
      endif

c	check to see if it is converged.  at any point which convergence
c	fails, stop checking and allow psi to continue uninterupted.

      if((frstgr.eq.0).and.(dabs(enrg0-enrg1).gt.cmxnrg)) then
	write(iout,*)
	write(iout,*) '  Energy is still changing.'
	write(iout,*) '  Calculation continuing.'
	goto 999
      else if(gmax1.gt.cgmaxc) then
        write(iout,*)
        write(iout,*) '  Largest Cartesian gradient is too big.'
        write(iout,*) '  Calculation continuing.'
        goto 999
      else if(grms1.gt.cgrmsc) then
        write(iout,*)
        write(iout,*) '  RMS Cartesian gradient is too big.'
        write(iout,*) '  Calculation continuing.'
        goto 999
      endif

c	time to check the Cartesian geometry change

      if(frstgr.eq.0) then
	drms = 0.0d+00
	dmax = 0.0d+00
	dtemp = 0.0d+00
	do 210 i=1,natom
	  do 220 j=1,3
	    dtemp = dabs(x1(i,j)-x0(i,j)) 
	    if(dtemp.gt.dmax) dmax = dtemp
	    drms = dtemp*dtemp + drms
220	  continue
210	continue
	drms = dsqrt(drms/(natom*3))
	if(dmax.gt.cdmaxc) then
          write(iout,*)
          write(iout,*) '  Largest Cartesian geometry change is '
	  write(iout,*) '  still quite large.'
          write(iout,*) '  Calculation continuing.'
          goto 999
	elseif(drms.gt.cdrmsc) then
          write(iout,*)
          write(iout,*) '  RMS Cartesian geometry change is '
          write(iout,*) '  still quite large.'
          write(iout,*) '  Calculation continuing.'
          goto 999
	endif
      endif

c	currently there is no capacity for internal coordinate analysis

      if(noint.eq.0) then
	write(iout,*)
        write(iout,*) '  Sorry, this program can not do internal '
	write(iout,*) '  coordinate analysis right now.'
      endif

c	so it passed all tests.  it is converged.

      write(iout,*)
      write(iout,*) '  Calculation is converged.'
      write(iout,*)
      call tstop(iout)
      call psdone()
#if defined(NO_EXIT_ROUTINE)
          stop 11
#else
          call exit(11)
#endif

c	it can only get to the following lines if things are not 
c	converged and it hits a goto 999

999   continue
      call tstop(iout)
      call psdone()
      stop
      end
