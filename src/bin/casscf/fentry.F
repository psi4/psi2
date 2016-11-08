      subroutine fentry(cc,ia,maxcor)
C*    PROGRAM CASSCF
C**********************************************************************
C   A Small Complete Active SCF (CASSCF) Program                      *
C**********************************************************************
C*   Last updated on Febuary 7, 1995 by Evan Bolton                   *
C*  Reason:  Better error checking, etc.       See comments Ceeb      *
C**********************************************************************
C*   Last updated on Febuary 2, 1995 by Evan Bolton                   *
C*  Reason:  To put in PSI 2.0 input parsing.  See comments Ceeb      *
C**********************************************************************
C*  Last updated on November 10, 1993 by Dr. Yukio Yamaguchi          *
c*  An option for a smallest determinant value is added               *
C**********************************************************************
C*  This program loads two big (MO integrals and two-PDM) files       *
c*  on the core memory                                                *
C**********************************************************************
C*  Last updated on August 30, 1993 by Dr. Yukio Yamaguchi            *
C*  Function `YMAT' was replaced by Function `Y4MAT'                  *
C**********************************************************************
C*  Last updated on August 21, 1993 by Dr. Yukio Yamaguchi            *
C*  The dimension of ioff was increased to 256*(256+1)/2 +1 = 32897   *
c**********************************************************************
C*  Last updated on July 14, 1993 by Dr. Yukio Yamaguchi              *
C**********************************************************************
C                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                    
      dimension cc(maxcor),ia(maxcor*2)                                
      dimension indmc(300),gmcs(300),emcs(300)

Ceeb Needed For Psi2.0 Input
      integer errcod, inold
      integer isizof, fwkclr, fwkadd, frdi, frdf
      logical kill
Ceeb

      COMMON/BASIS/NBASIS,NTRI,NTOT,NSYM
      COMMON/CALIF/LPARA(1024),APARA(1024)                              
      COMMON/CIENR/CIELEC,CITOT                                         
      COMMON/CLIMT/cmlim,dlamb
      COMMON/ENRGY/ESCF,ENUC                                            
      COMMON/GRADS/GLIM,GMNORM
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          
      COMMON/MFSEC/MASTER,NSECT                                         
      common/posit/ioemc,iopmc,ilgmc
      COMMON/SIGNS/IOFF(32897),IPRNT 
      COMMON/TAPES/ITAP35,ITAP56
      COMMON/WORKS/jtape1,jtape2
      COMMON/CI101/NOCMO,NOCTRI,NOCTOT
      COMMON/CI102/NIND,NINTRI
      COMMON/CI103/MXMIC,IDLIM
      DATA A00,TEN / 0.0D+00 , 10.0D+00 /                               
    1 FORMAT(//,2X,' A SMALL CASSCF PROGRAM  (VERSION OF 02/03/95)',
     1 /,11x,
     2 'Written by Dr. Yukio Yamaguchi',//,11X, 'Modifications for ',
     3 'PSI 2.0 use by Evan Bolton',///)
    2 FORMAT(14I5)                                                     
    3 FORMAT(//,2X,' PARAMETERS'/                                      
     1 2X,' ECONVERGENCE = ',I8/  
     2 2X,' DCONVERGENCE = ',I8/                                         
     3 2X,' GCONVERGENCE = ',I8/                                       
     4 2X,' MAXITERATION = ',I8/
     5 2X,' MOINCREMENT  = ',F8.3/                                        
     6 2X,' HESSHIFT     = ',F8.3/                                      
     7 2X,' MICROITER    = ',I8/                                   
     8 2X,' DLIMIT       = ',I8/            
     9 2X,' PRINT        = ',I8/)
    4 FORMAT(2X,' NBASIS = ',I8/
     1 2X,' NTRI   = ',I8/                               
     2 2X,' NTOT   = ',I8/                                              SMM00450
     3 2X,' NSYM   = ',I8/                                              SMM00470
     4 2X,' NST    = ',I8/                                              SMM00480
     5 2X,' IOCC   = ',I8/                                              SMM00490
     6 2X,' JOCC   = ',I8/                                              SMM00500
     7 2X,' KOCC   = ',I8/                                              SMM00510
     8 2X,' NCORE  = ',I8/                                              SMM00670
     9 2X,' NVIRT  = ',I8/                                              SMM00670
     a 2X,' NOLST  = ',I8/                                              SMM00670
     b 2X,' NACMO  = ',I8/                                              SMM00670
     c 2X,' NOCMO  = ',I8/                                              SMM00630
     d 2X,' NOCTRI = ',I8/                                              SMM00630
     e 2X,' NOCTOT = ',I8)                                              SMM00660
    7 FORMAT(//,2X,' PREVIOUS MCSCF MOLECULAR ORBITALS'/) 
    8 FORMAT(//,2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)                SMM00710
    9 FORMAT(//,2X,' REQUIRED MEMORY EXCEEDS MAXCOR'/                   SMM00720
     1          2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)                SMM00730
   10 FORMAT(//,2X,' NBASIS IS TOO LARGE'/                              SMM00740
     1          2X,' NABSIS = ',I5,' MAXBAS = ',I5/)                    SMM00750
   11 FORMAT(//,2X,' MCITER',4X,' I.P.',11X,' DELE',18X,' EMCSCF'/)     SMM00760
   12 FORMAT(//,2X,' MCITER',4X,' I.P.',11X,' GMNORM ',15X, ' DELE',
     1         18X,' EMCSCF'/)
   13 FORMAT(3X,I3,6X,I5,3X,F25.15,3X,F25.15,3x,f25.15)
   14 FORMAT(//,2x,' CONVERGENCE CRITERIA NOT MET'// 
     1          2x,' MCSCF PROCEDURE TO BE CONTINUED'/)
   15 FORMAT(//,2X,' NO MO INDEPENDENT PAIRS'/)                         SMM00810
   16 FORMAT(//,2X,' SPECIFIED MCSCF ENERGY CONVERGENCE ACHIEVED'/) 
   17 FORMAT(//,2X,' SPECIFIED MCSCF DENSITY CONVERGENCE ACHIEVED'/)
   19 FORMAT(//,2X,' WARNING:  You have exceeded 20 MCSCF iterations'
     1 //,5X,'Please check your MCSCF expansion as something may be ',
     2 'wrong.'/)
   20 FORMAT(//,2X,' LATEST MCSCF MOLECULAR ORBITALS'/)
   21 format(2i5)
   22 format(2i5,2f25.15)
   23 format(2x,'-----',5X,'-----',3X,'-------------------------',3X,
     1          '-------------------------',3X,
     2          '-------------------------',3X)
   24 format(/2x,'#ITER',6X,'#IP',6X,'RMS of MCSCF Gradient',6x,
     1          'Change in CASSCF Energy',10x,'CASSCF Energy')
C                                                                       SMM00860
      call psinit('APPEND')
      CALL TSTART(6)                                                    SMM00940
      CALL INITMF(1)                                                    SMM00960

C:::::::::::::::::::::::::::::::::::::::::
C:::Def:ine machine dependent variables:::
C:::::::::::::::::::::::::::::::::::::::::
      maxbas=256
      maxsiz=maxbas*(maxbas+1)/2
      maxscf=300

Ceeb
      kill = .false.

c
C   ITAP35 = TWO ELECTRON MO INTEGRALS                                  SMM00870
C   ITAP56 = TWO PARTICLE DENSITY MATRIX                                SMM00890
C   jtape91 = WORKING FILE
C   jtape92 = WORKING FILE 
C                                                                       SMM00930
      ITAPE3=3                                                          SMM00980
      ITAPE4=4                                                          SMM00980
      INPUT=5                                                           SMM00990
      ITAPE6=6                                                          SMM01000
      ITAP14=14                                                         SMM01000
      ITAP35=35                                                         SMM01020
      ITAP56=56                                                         SMM01040
      jtape1=91                                                         SMM01050
      jtape2=92                                                         SMM01060
c====================================================================
c   Followings are the position of arrays in the master file
c====================================================================
      ioemc=23
      iopmc=39
      ilgmc=37
c====================================================================
c

      call ffile(itape4,'file4',0)
      call ffile(itap14,'file14',0)
      WRITE(6,1)                                                        SMM01170
      WRITE(3,1)                                                        SMM01180
      WRITE(4,1)                                                        SMM01190
      write(6,*) '  maxcor = ',maxcor

      CALL LOCATE(INPUT,'# CASSCF #',IERR)

Ceeb  Check to see if using old style input
      if (ierr.ne.0) inold = 0
      if (ierr.eq.0) inold = 1

      if (inold.eq.1) then
        write(6,*)
        write(6,*) 'WARNING: using old style ''# CASSCF #'''
      endif
Ceeb

c
c   Parameters from the master file
      nbasis=lpara(3)
      ntri=lpara(4)
      if(nbasis.gt.maxbas) then
        write(6,*) '  Nbasis exceeds Maxbas -----> Stop ! '
        write(6,*) '  Nbasis = ',nbasis,' Maxbas = ',maxbas
        write(6,*) '  '
        write(6,*) '  '
        call mabort
        stop
      end if
      if(ntri.gt.maxsiz) then
        write(6,*) '  Ntri exceeds Maxsiz -----> Stop ! '
        write(6,*) '  Ntri = ',ntri,' Maxsiz = ',maxsiz
        write(6,*) '  '
        write(6,*) '  '
        call mabort
        stop
      end if
c
c   Form the IOFF array
      IOFF(1)=0                                          
      DO 101 I=1,ntri
      IOFF(I+1)=IOFF(I)+I                                          
  101 CONTINUE                                                

c   
c======================================================================
c   Ieconv  =  tolerance for the MCSCF energy
c   Idconv  =  tolerance for the density matrix
c   Igconv  =  tolerance for the Lagrangian matrix
c   Maxit   =  maximum number of MCSCF iterations
c   Nmlim   =  Increment for MO changes
c   Nlamb   =  Increment for level shift
c   Mxmic   =  maximum number of micro iterations
c   Idlim   =  the smallest positive determinant value allowed
c   Iprint  =  Print option
c======================================================================
C
Ceeb  Read in Old style input
      if (inold.eq.1) then
C
C   READIN PARAMETERS
        read(5,2) ieconv,idconv,igconv,maxit,nmlim,nlamb,mxmic,idlim,iprnt

c======================================================================
c   Note maximum increment in MO changes and level shift       
c   Maxmum changes in MO coefficients are [ tlimit / nmlim ]           
c   Maxmum change in level shift is [ tlimit / nlamb ]                 
c      where tlimit = 2.00  in this version
c======================================================================
c   The smallest positive determinant value allowed before level shift
c     is [ dlim = 1.00 / 10.0**(idlim) ]
c======================================================================
        tlimit=2.00d+00
        if(nmlim.eq.0) nmlim=20
        if(nlamb.eq.0) nlamb=20
        cmlim=tlimit/dfloat(nmlim)
        dlamb=tlimit/dfloat(nlamb)
      endif
Ceeb  Finished reading in Old Style input

Ceeb  Read in New Style input
      if (inold.eq.0) then
        errcod = fwkclr()
        errcod = fwkadd(5,'CASSCF')
        errcod = fwkadd(5,'DEFAULT')

        errcod = frdi('ECONVERGENCE',ieconv)
        errcod = frdi('DCONVERGENCE',idconv)
        errcod = frdi('GCONVERGENCE',igconv)

Ceeb        errcod = frdi('CASITER',maxit)

        cmlim = 0.010D+00
        errcod = frdf('MOINCR',cmlim)

        dlamb = 0.010D+00
        errcod = frdf('HESSHIFT',dlamb)

        errcod = frdi('MICROITER',mxmic)
        errcod = frdi('DLIMIT',idlim)
        errcod = frdi('PRINT',iprnt)
      endif
Ceeb  Finished reading in New Style input

Ceeb  This sets up the Defaults properly after reading the input
      if(ieconv.eq.0) ieconv=8
      if(idconv.eq.0) idconv=8
      if(igconv.eq.0) igconv=8
      ECONV=TEN**(-IECONV)                                             
      DCONV=TEN**(-IDCONV)                                         
      GLIM=TEN**(-IGCONV)
      if(maxit.eq.0) maxit=20

      if(mxmic.eq.0) mxmic=30
      if(idlim.eq.0) idlim=5

c======================================================================
C                                                                       SMM01350
C   READIN SCF INFORMATION FROM THE MASTER FILE                         SMM01450
      NTOT=IOFF(NTRI+1)                                                 SMM01480
      IOCC=LPARA(7)                                                     SMM01490
      KOCC=LPARA(8)                                                     SMM01500
      JOCC=LPARA(9)                                                     SMM01510
      NSYM=LPARA(13)                                                    SMM01540
      NST=LPARA(30)                                                     SMM01570
      ENUC=APARA(1)                                                     SMM01580
      ESCF=APARA(2)                                                     SMM01590
C                                                                       SMM01680
C   Find numbers of frozen core and deleted vitual oribtals
      ic1=1
c...............q.........................
      call orbs(cc(ic1),ncore,nvirt,nolst)
      nacmo=nbasis-ncore-nvirt
c======================================================================
c   Note
c    A parameter ``NOCMO''  
c      is a number for the highest occupied molecular orbital
c      However, unoccupied orbitals may exist below NOCMO
c       in present file structure)
c======================================================================
      NOCMO=NCORE+NACMO                                                
      nocmo=max0(nocmo,nolst)
      NOCTRI=IOFF(NOCMO+1)                                         
      NOCTOT=IOFF(NOCTRI+1)
C
      write(6,3) ieconv,idconv,igconv,maxit,cmlim,dlamb,mxmic,idlim,
     1           iprnt
      WRITE(6,4) NBASIS,NTRI,NTOT,NSYM,NST,IOCC,JOCC,KOCC,
     1           ncore,nvirt,nolst,nacmo,nocmo,noctri,noctot

Ceeb Flush out parameters before things get started...
      call fioflu(6)
C                                                                      
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C  DYNAMIC ALLOCATION                                                
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C     IC1  = EIG(NBASIS)                                          
C     IC2  = OCC(NBASIS)                                      
C     IC3  = ESO(NBASIS,NBASIS)                           
C     IC4  = H(NTRI)                                               
C     IC5  = Q(NTRI)             
C     IC6  = X(NBASIS,NBASIS)                               
C     IC7  = NIJ(NTRI,2)
C     IC8  = GMO(NTRI)
C     IC9  = RMO(NTRI)
C     IC10 = G(NOCTOT)
C     IC11 = CC(NTOT)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      IC1=1                                                          
      IC2=IC1+NBASIS                                          
      IC3=IC2+NBASIS                                    
      IC4=IC3+NBASIS*NBASIS                          
      IC5=IC4+NTRI                                           
      IC6=IC5+NTRI          
      IC7=IC6+NBASIS*NBASIS        
      IA7=IC7+IC7-1
      IC8=IC7+NTRI
      IC9=IC8+NTRI
      ickeep=ic9+ntri
      IC10=ickeep
      IC11=IC10+NOCTOT
      icsave=ic11+ntot
      write(3,*) ' ic1    = ',ic1
      write(3,*) ' ic2    = ',ic2
      write(3,*) ' ic3    = ',ic3
      write(3,*) ' ic4    = ',ic4
      write(3,*) ' ic5    = ',ic5
      write(3,*) ' ic6    = ',ic6
      write(3,*) ' ic7    = ',ic7
      write(3,*) ' ia7    = ',ia7
      write(3,*) ' ic8    = ',ic8
      write(3,*) ' ic9    = ',ic9
      write(3,*) ' ickeep = ',ickeep
      write(3,*) ' ic10   = ',ic10
      write(3,*) ' ic11   = ',ic11
      write(3,*) ' icsave = ',icsave
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::SMM02480
C                                                                       SMM02490
C   READ IN SO-MO EIGENVECTORS                                          SMM02500
      CALL MREAD(CC(IC1),17)                                            SMM02510
      CALL MREAD(CC(IC2),18)                                            SMM02520
      CALL MREAD(CC(IC3),19)                                            SMM02530
      IF(IPRNT.LE.2) GO TO 301                                          SMM02540
      WRITE(6,7)                                                        SMM02550
      CALL EIGOUT(CC(IC3),CC(IC1),CC(IC2),NBASIS,NBASIS,NBASIS,NBASIS,6)SMM02560
C                                                                       SMM02570
C=======================================================================SMM02580
C==  START THE MCSCF PROCEDURE                                        ==SMM02590
C=======================================================================SMM02600
  301 CONTINUE                                                          SMM02610
C                                                                       SMM02620
C   Read in information for previous calculations
      rewind itap14
Ceeb  read(itap14,21) mciter
      read(itap14,21,end=302) mciter
      if(mciter.eq.0) go to 302
      do 102 i=1,mciter
      read(itap14,22) ijunk,indmc(i),gmcs(i),emcs(i)
      if (ijunk.ne.i) then
         write(6,*) 'Bad CASSCF restart file'
         call mabort
      endif
  102 continue
      rewind itap14
c
  302 continue
      MCITER=MCITER+1                                                   SMM02650
      WRITE(3,25) MCITER                                                SMM02660
      WRITE(4,25) MCITER                                                SMM02670
      WRITE(6,25) MCITER                                                SMM02670
   25 FORMAT(//,2X,' MCITER = ',I5/)                                    SMM02680
      if(mciter.gt.maxscf) then
        write(6,*) '  mciter exceeds maxscf -----> Stop ! '
        write(6,*) '  Mciter = ',mciter,' Maxscf = ',maxscf
        write(6,*) '  '
        write(6,*) '  '
        call mabort
        stop
      end if
C                                                                      
C   CALCULATE THE CI (MCSCF) ENERGY AS A TEST
      WRITE(3,30)                                                      
   30 FORMAT(//,2X,' NOW YOU ARE IN ENERGY'/)                        
      ic12=icsave
      ia12=ic12+ic12-1
      ICMAX=IC12+MAXBF2   
      write(3,*) ' ic12 = ',ic12
      write(3,*) ' ia12 = ',ia12
      WRITE(3,8) ICMAX,MAXCOR                                          
      maxuse=icmax
      IF(ICMAX.GT.MAXCOR) GO TO 399                                   
C.................H       Q       G        LBLI     BUFI.....
      CALL ENERGY(CC(IC4),CC(IC5),CC(IC10),IA(IA12),CC(IC12))
C                                                                      
C   FORM THE LAGRANGIAN MATRIX                                         
      WRITE(3,31)                                                   
   31 FORMAT(//,2X,' NOW YOU ARE IN LAGRAN'/)                   
      ic12=icsave
      IA12=IC12+IC12-1             
      ICMAX=IC12+MAXBF2       
      write(3,*) ' ic12 = ',ic12
      write(3,*) ' ia12 = ',ia12
      maxuse=max0(icmax,maxuse)
      WRITE(3,8) ICMAX,MAXCOR                                           SMM03990
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     SMM04000
C.................H       Q       X       G        CC.......   
      CALL LAGRAN(CC(IC4),CC(IC5),CC(IC6),CC(IC10),CC(IC11), 
C.................LBLI     BUFI.....          
     1            IA(IA12),CC(IC12))        
C                                                                       SMM04050
C   FORM THE MO GRADIENT VECTOR                                         SMM04060
      WRITE(3,32)                                                       SMM04070
   32 FORMAT(//,2X,' NOW YOU ARE IN MOGRAD'/)                           SMM04080
      icmax=icsave
      maxuse=max0(icmax,maxuse)
      WRITE(3,8) ICMAX,MAXCOR                                           SMM04130
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     SMM04140
C.................X       NIJ     GMO.....
      CALL MOGRAD(CC(IC6),IA(IA7),CC(IC8))
      if(nind.gt.maxsiz) then
        write(6,*) '  Nind exceeds maxsiz -----> Stop ! '
        write(6,*) '  Nind = ',nind,' Maxsiz = ',maxsiz
        write(6,*) '  '
        write(6,*) '  '
        call mabort
        stop
      end if
c
      IF(NIND.EQ.0) GO TO 510                                           SMM04170
      nintri=ioff(nind+1)
C                                                                       SMM04570
C   FORM THE A11 MATRIX                                                 SMM04580
      WRITE(3,33)                                                       SMM04590
   33 FORMAT(//,2X,' NOW YOU ARE IN A11MAT'/)                           SMM04600
      ic12=icsave
      IA12=IC12+IC12-1
      ICMAX=IC12+MAXBF2
      write(3,*) ' ic12 = ',ic12
      write(3,*) ' ia12 = ',ia12
      maxuse=max0(icmax,maxuse)
      WRITE(3,8) ICMAX,MAXCOR                                           SMM04660
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     SMM04670
C.................H       Q       X       NIJ     G        CC.......
      CALL A11MAT(CC(IC4),CC(IC5),CC(IC6),IA(IA7),CC(IC10),CC(IC11),
C.................LBLI     BUFI.....     
     1            IA(IA12),CC(IC12))
C                                                                       SMM04720
C   SOLVE THE SIMULTANEOUS EQUATIONS
      WRITE(3,34)                                                       SMM05230
   34 FORMAT(//,2X,' NOW YOU ARE IN SIMULT'/)
      ic10=ickeep
      ic11=ic10+nind*(nind+1)
      ia11=ic11+ic11-1
      icmax=ic11+maxbf2
      maxuse=max0(icmax,maxuse)
      WRITE(3,8) ICMAX,MAXCOR                                           SMM05280
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     SMM05290
C.................GMO     RMO     CC       LBLI     BUFI.....
      CALL SIMULT(CC(IC8),CC(IC9),CC(IC10),IA(IA11),CC(IC11))
C                                                                       SMM05340
C   FIND AN IMPROVED MO EIGENVECTOR MATRIX                              SMM05700
      WRITE(3,35)                                                       SMM05710
   35 FORMAT(//,2X,' NOW YOU ARE IN UNITRY'/)                           SMM05720
      ic10=ickeep
      ic11=ic10+nbasis*nbasis
      ic12=ic11+nbasis*nbasis
      ic13=ic12+nbasis*nbasis
      ic14=ic13+nbasis
      ia14=ic14+ic14-1
      icmax=ic14+nbasis
      maxuse=max0(icmax,maxuse)
      WRITE(3,8) ICMAX,MAXCOR                                           SMM05800
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     SMM05810
C.................ESO     NIJ     RMO     R        U        W........
      CALL UNITRY(CC(IC3),IA(IA7),CC(IC9),CC(IC10),CC(IC11),CC(IC12),
C.................VEC      NORD.....           
     1            CC(IC13),IA(IA14))       
c
c   Update the eigenvectors in file30
      write(3,36)
   36 format(//,2x,' now you are in update'/)
      ic10=ickeep
      ic11=ic10+nbasis
      ic12=ic11+nbasis*nbasis
      icmax=ic12+nbasis*nbasis
      maxuse=max0(icmax,maxuse)
      WRITE(3,8) ICMAX,MAXCOR                                           SMM05800
      IF(ICMAX.GT.MAXCOR) GO TO 399                                     SMM05810
c.................eso     eig      ee       ess......
      call update(cc(ic3),cc(ic10),cc(ic11),cc(ic12))
c
C=======================================================================SMM05980
C==  END OF THE MCSCF PROCEDURE                                        =SMM05990
C=======================================================================SMM06000
      indmc(mciter)=nind
      gmcs(mciter)=gmnorm
      emcs(mciter)=citot
      apara(8)=citot
      call mwrit(apara,3)

      write(6,24)
      write(6,23)
      write(4,24)
      write(4,23)
      do 103 iii=1,mciter
      if(iii.eq.1) then
        dele=a00
      else
        dele=emcs(iii)-emcs(iii-1)
      end if
      WRITE(4,13) iii,indmc(iii),gmcs(iii),dele,emcs(iii)
      WRITE(6,13) iii,indmc(iii),gmcs(iii),dele,emcs(iii)
  103 continue
c
c   Store current information on file14
      rewind itap14
      write(itap14,21) mciter
      do 104 i=1,mciter
      write(itap14,22) i,indmc(i),gmcs(i),emcs(i)
  104 continue
      rewind itap14
c
      IF(MCITER.GE.2) THEN                                              SMM06110
        IF(DABS(DELE).LE.ECONV) GO TO 511                               SMM06120
        IF(GMNORM.LE.DCONV) GO TO 512                                   SMM06130
        IF(MCITER.GE.MAXIT) GO TO 514
      END IF   
      write(6,14)
Ceeb  GO TO 400                                                         SMM06170
      GO TO 401                                                         SMM06170
C
C===============================================
C   The MCSCF procedures has been terminated
C===============================================
  510 CONTINUE                                                          SMM06180
      WRITE(6,15)                                                       SMM06190
      GO TO 400                                                         SMM06200
  511 CONTINUE                                                          SMM06210
      WRITE(6,16)                                                       SMM06220
      GO TO 400                                                         SMM06230
  512 CONTINUE                                                          SMM06240
      WRITE(6,17)                                                       SMM06250
      GO TO 400                                                         SMM06260
  514 CONTINUE                                                          SMM06300
      WRITE(6,19)                                                       SMM06310
Ceeb  GO TO 400                                                         SMM06320
      GO TO 401                                                         SMM06320
C                                                                       SMM06330
  399 CONTINUE                                                          SMM06340
      WRITE(6,9) ICMAX,MAXCOR                                           SMM06350
      WRITE(3,9) ICMAX,MAXCOR                                           SMM06360

Ceeb  Notification of amount of memory required
      Ieebtmp = icmax * isizof('REAL_IN_BYTES') / 1000 + 1

      write(6,*) '  '
      write(6,*) '  Need approximately ',Ieebtmp,' KB of memory.'
      write(6,*) '  '

Ceeb "Proper" job termination for lack of memory
      kill = .true.
Ceeb

      go to 600 
C
C================================================
C   Print out the MOs for a meaningful result
C================================================
  400 CONTINUE                                                          SMM06370

Ceeb "Proper" job termination for Convergence criteria and No MO pairs
      kill = .true.
c
  401 CONTINUE                                                          SMM06370
Ceeb
      IF(IPRNT.LE.2) GO TO 600
      WRITE(6,20)
      CALL MATOUT(CC(IC3),NBASIS,NBASIS,NBASIS,NBASIS,6)
C
  600 CONTINUE
      write(6,*) '  '
      write(6,*) '  largest core used = ',maxuse
      write(6,*) '  '
      write(6,*) '  '
      CALL RCLOSE(MASTER,3)                                             SMM06380

      if (kill) then
        write(6,*) ' '
        write(6,*) ' '
        write(6,*) ' MCSCF Procedure Terminated due to CONVERGENCE'
        write(6,*) '   or PROBLEM.  Please closely read the output'
        write(6,*) '   and take preventative measures before '
        write(6,*) '   continuing this expensive procedure.'
        write(6,*) ' '
        write(6,*) '  '
        write(6,*) '  '
        call mabort 
      endif
C                                                                       SMM06400
      CALL TSTOP(6)                                                     SMM06410
      STOP                                                              SMM06420
      END                                                               SMM06430
