      SUBROUTINE MAININ(IZ,Z,MAXCOR)                                      
      IMPLICIT INTEGER (A-Z)                                            
cets030291
#include <error.h>
      integer errcod, frdi
      integer ip, prcntr, prdup
      real*8 Z(MAXCOR),FOCC,DRTVER,APARA                                
      real*8 ci1,ci2,citr,h11,h12,h22,one,derc1,derc2
      INTEGER IZ(1),LABEL(26),ia47(192),cicoef(2)                       
      LOGICAL CCSD                                                      
      COMMON /TAPES/ NOUT,ICHK,IT30,IT42,IT47,IT52,IT58,IT94,IT96,IT98  
      COMMON /DIMAO/ NBFAO,NBATRI                                       
      COMMON /DIMSG/ SYMORB,NUMIJ,NBF,NMAX,NBFORB,NSYM,NDF,NORBS        
      COMMON /POINT/ FOCC(5),IPOINT(5),JPOINT(5),NSORB(5)               
      common /ptrs/ motyp(500),mopos(5)
      common /inpspe/ kspe,lspe(5)
      common /tcscf/ derc1,derc2,ci1,ci2,h11,h12,h22,one,itcscf
      COMMON /FUNCS/ IOCC,JOCC,KOCC,NTYPES,NBSET,NAT                    
      COMMON /OCCS/ NSPE,NSA,NSB,NALP,NOP,NDOC,LFERM,NUOC,IJIND,IJIND3  
      COMMON/OCCS2/NUMCOR,NCOR,FDOC,CORDOC,NUMVIR,FVIR,LUOC,NUMDEP,VUOC 
      COMMON /LOC42/ IA42(100),IBF2AT(284)                              
      COMMON /IEXT/ IPRINT,IGRP,ICIDIP                                  
      COMMON /CALCJ/ ICALT,IDERT,ICIT                                   
      COMMON/DEP/IJDC,IFZC,IFZV,IZDIM                                   
      COMMON /CALIF/ LPARA(1024),APARA(1024)                            
      equivalence (cicoef,citr)
C                                                                       
      WRITE(NOUT,1000)                                                  
 1000 FORMAT(//10X,'THE GENERAL COUPLED PERTURBED HARTREE-FOCK PROGRAM',
     +       /,10X,'          BASED ON LUCK AND A GOOD GUESS         ', 
     +      //,10X,'            GRCPHF PROGRAM WRITTEN BY            ', 
     +       /,10X,'          PATIENT BUT FRUSTRATED PEOPLE          ', 
     +      //,10X,'  hacked up to use less memory and handle tcscf  ',
     +      //,10X,'  version for use with derivth and new lag and  ',
     +       /,10X,'         backtr programs written by TPH ',/
     +       /,10X,'              There can be only one!         ',/
     +        ,10X,'              version: July 13, 1991         ',//)

C************************                                               
C  READ INPUT FILE.     *                                               
C************************                                               
      it11 = 11
      it15 = 15
      call ffile(it11,' ',1)
      call ffile(it15,' ',0)
      INPUT=5                                                           
cets030691
      errcod = prcntr('OFF BRIEF DEBUG')
      errcod = prdup()
      CALL LOCATE(INPUT,'# GRCPHF #',IERR)                              
      IF (IERR .EQ. 0) THEN                                             
        READ(INPUT,'(3I5)') IPRINT,ICIDIP,ICCSD                         
      ELSE                                                              
cets030291        WRITE(NOUT,990)
        ICCSD = 0                                                       
        IPRINT = 0                                                      
        ICIDIP = 1                                                      
        errcod = frdi('IPRINT',iprint)
 990    FORMAT(/' # GRCPHF # NOT FOUND IN INPUT FILE; USING DEFAULTS')  
      END IF                                                            
      CCSD=.FALSE.                                                      
      IF(ICCSD.NE.0) CCSD=.TRUE.                                        
C*********************                                                  
C     OPEN FILES.    *                                                  
C*********************                                                  
      CALL INITMF(1)                                                    
      CALL RFILE(IT30)                                                  
      CALL RFILE(IT42)                                                  
      CALL RFILE(IT47)                                                  
      call wreadw(it47,ia47,192,1,junk)
      cicoef(1) = ia47(14)
      cicoef(2) = ia47(15)
      ci1 = citr
      cicoef(1) = ia47(16)
      cicoef(2) = ia47(17)
      ci2 = citr
cets030691
      ip = prcntr('IS_ON BRIEF')
      if(ip.eq.0) then
       write(6,*) ' '
       write(6,*) ' ci coefficients '
       write(6,*) ' ci1 = ',ci1,' ci2 = ',ci2
       write(6,*) ' '
      endif
      CALL RFILE(IT58)                                                  
C*********************************************                          
C     GET CONSTANTS FROM MASTER FILE 40.     *                          
C*********************************************                          
      IOCC=LPARA(7)                                                     
      JOCC=LPARA(8)                                                     
      KOCC=LPARA(9)                                                     
      NBF=LPARA(12)                                                     
      NBFAO=LPARA(11)                                                   
      NAT=LPARA(10)                                                     
      NTYPES=LPARA(18)                                                  
      ICALT=LPARA(19)                                                   
      ICIT=LPARA(20)                                                    
      IDERT=LPARA(21)                                                   
      NBATRI=NBFAO*(NBFAO+1)/2                                          
      NDF=NAT*3                                                         
      NIJ=NBF*(NBF+1)/2                                                 
      NBSET=NTYPES+1                                                    
C*********************                                                  
C  CHECK FOR ERRORS  *                                                  
C*********************                                                  
      IF (ICALT .LE. 3) THEN                                            
        ITCSCF=0                                                        
      ELSE IF (ICALT .EQ. 4) THEN                                       
        NBSET=7                                                         
        ITCSCF=1                                                        
        write(6,*) ' '
        write(6,*) '  !!! For TCSCF you better check your
     1 gradients carefully !!!'
        write(6,*) ' '
      ELSE                                                              
        WRITE(NOUT,1001) ICALT                                          
        CALL ABORTS                                                     
 1001   FORMAT(/'  ERROR: THIS PROGRAM DOES GENERAL SCF OR TCSCF'       
     1  /'  CAN NOT HANDLE YOUR CALTYP OF :',I5)                        
      END IF                                                            
      IF (ICIT.EQ.2) THEN                                               
        IF (IDERT.NE.1) THEN                                            
          WRITE(NOUT,1004)                                              
1004      FORMAT(/' ERROR: CAN ONLY DO FIRST DERIVATIVES OF CI OR CC')  
          CALL ABORTS                                                   
        END IF                                                          
      END IF                                                            
      IF(ICIT.GT.2) THEN                                                
         WRITE(NOUT,1005)                                               
1005     FORMAT(/' ERROR: ONLY CI OR CC DERIVATIVES CAN BE DONE.')      
         CALL ABORTS                                                    
      END IF                                                            
C*********************************************                          
C     GET CONSTANTS FROM MASTER FILE 58.     *                          
C*********************************************                          
      CALL WREADW(IT58,DRTVER,INTOWP(1),1,END)                          
      CALL WREADW(IT58,LABEL,26,END,END)                                
      CALL WREADW(IT58,NBF,1,END,END)                                   
      CALL WREADW(IT58,NSYM,1,END,END)                                  
      CALL WREADW(IT58,NORBS,1,END,END)                                 
      END=END+3                                                         
      CALL WREADW(IT58,NREFS,1,END,END)                                 
      END=END+5                                                         
      CALL WREADW(IT58,ORBSYM,1,END,END)                                
      CALL WREADW(IT58,SYMORB,1,END,END)                                
      CALL WREADW(IT58,NUMIJ,1,END,END)                                 
      CALL WREADW(IT58,NGROUP,1,END,END)                                
      END=END+1                                                         
      CALL WREADW(IT58,NMAX,1,END,END)                                  
cets030691
      ip = prcntr('IS_ON DEBUG')
      if(ip.eq.1) then
      WRITE(NOUT,1008) DRTVER,LABEL,NBF,NSYM,NORBS,NREFS,ORBSYM,NUMIJ,  
     1NGROUP,NMAX                                                       
 1008 FORMAT(/'  READING DATA FROM DRT TAPE, VERSION=',F10.6            
     1/'  LABEL FROM DRT: ',26A3/                                       
     2/'  DRT INFORMATION:'/'  NBF=',I5,'  NSYM=',I5,'  NORBS=',I5,     
     3/'  NREFS=',I5,'  ORBSYM=',I5,' NUMIJ=',I5,'  NGROUP=',I5,        
     4'  NMAX=',I5)                                                     
      endif
      IF (NORBS .NE. NBF) THEN                                          
        WRITE(NOUT,1009) NORBS,NBF                                      
        CALL ABORTS                                                     
 1009   FORMAT(/'  THIS DRT CONTAINS FROZEN ORBITALS.'                  
     1  /'  NORBS=',I6,'   NBF=',I6)                                    
      END IF                                                            
C                                                                       
      IF (NREFS .NE. 1) THEN                                            
        WRITE(NOUT,1011) NREFS                                          
        CALL ABORTS                                                     
 1011   FORMAT(/'  THE NUMBER OF REFERENCES IS: ',I5,                   
     1  /'  THIS PROGRAM HANDLES ONLY SINGLE REFERENCES!')              
      END IF                                                            
C*********************                                                  
C    ZERO OUT CORE   *                                                  
C*********************                                                  
      core = intowp(maxcor)
      do 10 i=1,core
       iz(i) = 0
   10 continue
      CALL COINIT(IZ,core)                                    
C************************************                                   
C    ALLOCATE CORE FOR DRT ARRAYS   *                                   
C************************************                                   
      IOFF = COALLC(IZ,isizof('INT')*NUMIJ)                             
      IJORD = COALLC(IZ,isizof('INT')*NUMIJ)                            
      IJDEP = COALLC(IZ,isizof('INT')*NUMIJ)                            
      KADD = COALLC(IZ,isizof('INT')*SYMORB)                            
      LADD = COALLC(IZ,isizof('INT')*SYMORB)                            
      IJADD = COALLC(IZ,isizof('INT')*NUMIJ)                            
      IJGRP = COALLC(IZ,isizof('INT')*NUMIJ)                            
      ORBSYM = COALLC(IZ,isizof('INT')*NORBS)                           
      DOFF = COALLC(IZ,isizof('INT')*NUMIJ)                             
      DADD = COALLC(IZ,isizof('INT')*NUMIJ)                             
      JT1C = COALLC(IZ,isizof('INT')*(4+NBSET)*NDF)                     
      OCC = COALLC(IZ,isizof('REAL')*NORBS)                             
      ALPHA = COALLC(IZ,isizof('REAL')*NUMIJ)                           
      BETA = COALLC(IZ,isizof('REAL')*NUMIJ)                            
      EPS = COALLC(IZ,isizof('REAL')*NUMIJ)                             
      VEC = COALLC(IZ,isizof('REAL')*NBFAO*NORBS)                       
      EAVF = COALLC(IZ,isizof('REAL')*NUMIJ)                            
      INTS = COALLC(IZ,isizof('REAL')*NMAX)                             
      EPA = COALLC(IZ,isizof('REAL')*NORBS*NORBS)                       
      if(itcscf.eq.1) then
       ha11 = COALLC(IZ,isizof('REAL')*ndf)                       
       ha12 = COALLC(IZ,isizof('REAL')*ndf)                       
       ha22 = COALLC(IZ,isizof('REAL')*ndf)                       
      else
       ha11 = 1
       ha12 = 1
       ha22 = 1
      endif
C                                                                       
      ORBTBF = COALLC(IZ,isizof('INT')*NORBS)                           
      DROL = COALLC(IZ,isizof('INT')*NORBS)                             
      IOUT = COALLC(IZ,isizof('INT')*NORBS)                             
      JCODE = COALLC(IZ,isizof('INT')*NORBS)                            
C                                                                       
      CALL GETDRT(IZ(KADD),IZ(LADD),IZ(IJADD),IZ(IJGRP),IZ(ORBTBF),     
     1 IZ(IOUT),IZ(ORBSYM),IZ(JCODE),NGROUP)                            
      CALL ALLOCA(IZ(ORBTBF),IZ(JCODE),IZ(IJORD),IZ(IOFF),IZ(IJDEP),    
     .           IZ(ORBSYM),IZ(DOFF),IZ(DADD),NMA1,itcscf)      
C                                                                       
      IF (ITCSCF .EQ. 1) IJIND=IJIND+2                                  
      IF (ITCSCF .EQ. 1) NMA1 = NMA1 + 2                                
      IJIND2=IJIND                                                      
      IJIND3=IJIND2+NDF                                                 
C*************************************                                  
C    READ IN DATA FOR SHELL TYPES    *                                  
C*************************************                                  
      DO 200 I=1,5                                                      
        FOCC(I)=APARA(40+I)                                             
        NSORB(I)=LPARA(40+I)                                            
  200 CONTINUE                                                          
cets030691
      ip = prcntr('IS_ON DEBUG')
      if(ip.eq.1) then
       WRITE(ICHK,1050) (FOCC(III), III=1,5), (NSORB(JJJ),JJJ=1,5)       
 1050  FORMAT(/'  FOCC:',5F5.3/'  NSORB:',5I5)                           
      endif
      CALL REGIST                                                       
C********************************************                           
C    ALLOCATE CORE FOR PRELIMINARY ARRAYS   *                           
C********************************************                           
      EIG = COALLC(IZ,isizof('REAL')*NORBS)                             
      SCR = COALLC(IZ,isizof('REAL')*NUMIJ)                             
C----------------------------------------------------------             
C    GET SOME ARRAYS FROM IT40 AND ROTATE TO DRT ORDER    |             
C----------------------------------------------------------             
      CALL READ40(IZ(IOFF),IZ(DROL),IZ(ORBTBF),IZ(EIG),IZ(OCC),         
     1            IZ(ALPHA),IZ(BETA),IZ(EPS),IZ(VEC),IZ(SCR),CCSD,      
     2            itcscf,nsa,nsb,ci1,ci2)
      WRITE(*,*)'EIGENVECTOR FROM TAPE 40'                              
      CALL COFREE(IZ,SCR)                                               
      CALL COFREE(IZ,EIG)                                               
C-------------------------------------------------------------------    
C    GET DERIVATIVE INTEGRALS FROM TAPE42; TRANSFORM TO MO BASIS   |    
C-------------------------------------------------------------------    
      SS = COALLC(IZ,isizof('REAL')*NBATRI)                              
      SCR1 = COALLC(IZ,isizof('REAL')*NORBS)                            
      SCR2 = COALLC(IZ,isizof('REAL')*NBATRI)                           
      BUF = COALLC(IZ,isizof('REAL')*NBSET*NDF*NBATRI)                  
      if(itcscf.eq.1) then
       CALL DRMAT2(IZ(VEC),IZ(SCR1),IZ(SCR2),IZ(SS),IZ(BUF),IZ(BUF),     
     1             IZ(IOFF),IZ(JT1C))                                    
      else
       CALL DERMAT(IZ(VEC),IZ(SCR1),IZ(SCR2),IZ(SS),IZ(BUF),IZ(BUF),     
     1             IZ(IOFF),IZ(JT1C))                                    
      endif
      CALL COFREE(IZ,BUF)                                               
      CALL COFREE(IZ,SCR2)                                              
      CALL COFREE(IZ,SCR1)                                              
      CALL COFREE(IZ,SS)                                                
C-------------------------------------------------------------          
C     DROL IS SCRATCHED TO C1 SYMMETRY FOR COUPLED CLUSTER   |          
C-------------------------------------------------------------          
      IF(CCSD)THEN                                                      
       DO 9156 J=1,NORBS                                                
 9156  IZ(DROL+J-1)=J                                                   
      ENDIF                                                             
C------------------------------------                                   
C    CONSTRUCT DERIV FOCK MATRIX    |                                   
C    LSTEPA IS SAME AS EPA          |                                   
C------------------------------------                                   
      TA = COALLC(IZ,isizof('REAL')*NUMIJ)                              
      HA = COALLC(IZ,isizof('REAL')*NUMIJ)                              
      BUF = COALLC(IZ,isizof('REAL')*NBSET*NUMIJ)                       
      if(itcscf.eq.1) then
       CALL FAMAT2(IZ(EPA),IZ(TA),IZ(HA),IZ(BUF),IZ(IOFF),IZ(JT1C),       
     1             IZ(IOUT),IZ(DROL),IZ(EPA),iz(alpha),iz(beta))          
       CALL HAMAT(IZ(HA),IZ(BUF),IZ(IOFF),IZ(JT1C),iz(iout),iz(drol),
     1            iz(ha11),iz(ha22),iz(ha12))
      else
       CALL FAMAT(IZ(EPA),IZ(TA),IZ(HA),IZ(BUF),IZ(IOFF),IZ(JT1C),       
     1            IZ(IOUT),IZ(DROL),IZ(EPA))                             
      endif
      CALL COFREE(IZ,BUF)                                               
      CALL COFREE(IZ,HA)                                                
      CALL COFREE(IZ,TA)                                                
      CALL COFREE(IZ,JCODE)                                             
      CALL COFREE(IZ,IOUT)                                              
      CALL COFREE(IZ,DROL)                                              
c      CALL COFREE(IZ,ORBTBF)                                            
C----------------------------                                           
C    CONSTRUCT ETA,ZETA     |                                           
C----------------------------                                           
      ZETA = COALLC(IZ,isizof('REAL')*NUMIJ*(NTYPES+1))                      
      CALL GETINT(IZ(EPS),IZ(ZETA),IZ(OCC),IZ(ALPHA),IZ(BETA),          
     1            IZ(KADD),IZ(LADD),IZ(IJADD),IZ(IJGRP),IZ(ORBSYM),     
     2            IZ(INTS),IZ(IOFF),IZ(EAVF),CCSD)                      
C--------------------------------------------------------------         
C    CONSTRUCT A AND B0 MATRICES IN 'CC' AND SOLVE A*U=B0    |         
C--------------------------------------------------------------         
      IDIMCC=NMA1*(NMA1+1)                                              
      WRITE(*,*)'       '                                               
      WRITE(*,*)'DIMENSION OF ZZ IS   ',IZDIM                           
      WRITE(*,*)'DIMENSION OF B LAMDA ',IJIND*NDF                       
      WRITE(*,*)'DIMENSION OF A MATRIX',NMA1*NMA1                       
      WRITE(*,*)'DIMENSION OF Z VECTOR',NMA1                            
      ip = prcntr('IS_ON DEBUG')
      if(ip.eq.1) then
       write(6,*) 'maxcor', maxcor
       write(6,*) 'ijind', ijind  
       write(6,*) 'numdep', numdep
       write(6,*) 'ndf',ndf       
       write(6,*) 'idimcc', idimcc
      endif
C                                                                       
      C = COALLC(IZ,isizof('REAL')*IJIND*NDF)                           
      EX = COALLC(IZ,isizof('REAL')*IJIND)                              
      EX2 = COALLC(IZ,isizof('REAL')*IJIND)                             
      RS = COALLC(IZ,isizof('REAL')*NUMDEP*NDF)                         
      LG = COALLC(IZ,isizof('REAL')*NORBS*NORBS)                        
      if(itcscf.eq.1) then
       e11 = coallc(iz,isizof('REAL')*norbs*norbs)
       e12 = coallc(iz,isizof('REAL')*norbs*norbs)
       e22 = coallc(iz,isizof('REAL')*norbs*norbs)
      else
       e11 = 1
       e12 = 1
       e22 = 1
      endif
      BOX = COALLC(IZ,isizof('REAL')*NDF)                               
c     IF(BOX+INTOWP((NUMIJ+1)*NDF+IDIMCC)+40.LT.CORE) THEN   
       ICASE = 1                                                        
       CC = COALLC(IZ,isizof('REAL')*IDIMCC)                            
       DB = CC                                                          
       ITOP = CC + 40 + INTOWP(IDIMCC)                                  
c     ELSE IF(BOX+INTOWP((NUMIJ+1)*NDF+20*IZDIM)                         
c    1        +40.LT.CORE) THEN                              
c      ICASE = 2                                                        
c      DB = COALLC(IZ,isizof('REAL')*20*IZDIM)                          
c      CC = DB                                                          
c      ITOP = DB + 40 + INTOWP(20*IZDIM)                                
c     ELSE                                                              
      IF(ITOP.GT.CORE) THEN 
       WRITE(6,*) 'NOT ENOUGH CORE - ERROR'                             
       WRITE(NOUT,1030) ITOP/2,CORE/2                                      
       CALL ABORTS                                                      
      ENDIF                                                             
      WRITE(NOUT,1030) ITOP/2,CORE/2                                      
 1030 FORMAT(/'  PROGRAM REQUESTED ',I8,'  REAL WORDS BEFORE B0MAT.'         
     1  /'  THERE ARE ',I8,'  REAL WORDS AVAILABLE.')                        
      SS = COALLC(IZ,isizof('REAL')*NUMIJ*NDF)                          
      CALL B0MAT(IZ(ALPHA),IZ(BETA),IZ(EPS),IZ(EPA),IZ(SS),             
     1           IZ(ZETA),IZ(CC),IZ(OCC),IZ(INTS),                      
     2           IZ(IJORD),IZ(IOFF),IZ(KADD),IZ(LADD),IZ(IJGRP),        
     3           IZ(IJADD),IZ(ORBSYM),IZ(JT1C),                         
     4           IZ(LG),IZ(BOX),IZ(C),                                  
     5           IZ(EX),IZ(RS),IZ(IJDEP),                               
     6           IZ(EAVF),IZ(EPA),IZ(EX2),CCSD,                         
     7           IZ(DOFF),IZ(DADD),NMA1,IDIMCC,IZ(DB),ICASE,            
     8           iz(ha11),iz(ha22),iz(ha12),iz(e11),iz(e22),iz(e12))
      CALL COFREE(IZ,SS)                                                
C--------------------------------                                       
C    CALCULATE DIPOLE MOMENT    |                                       
C--------------------------------                                       
c      ZZ = COALLC(IZ,isizof('REAL')*(ijind+IZDIM))                              
cedtool020998
      zz = COALLC(iz,isizof('REAL')*(numij))
      CIDT = COALLC(IZ,isizof('REAL')*3)                                
      SC1 = COALLC(IZ,isizof('REAL')*NBFAO*NBFAO)                       
      SC2 = COALLC(IZ,isizof('REAL')*NBFAO*NBFAO)                       
      DIP = COALLC(IZ,isizof('REAL')*NBATRI)                            
c      IF(ICIDIP.EQ.1) CALL CIDIP(IZ(DIP),IZ(VEC),IZ(SC1),IZ(SC2),IZ(CC),
c     1                           IZ(IOFF),IZ(IJORD),IZ(IJDEP),IZ(EAVF), 
c     2                           IZ(LG),IZ(OCC),IZ(ZZ),IZDIM,IZ(CIDT))  
cedtool020998
      IF(ICIDIP.EQ.1) CALL CIDIP(IZ(DIP),IZ(VEC),IZ(SC1),IZ(SC2),IZ(CC),
     1                           IZ(IOFF),IZ(IJORD),IZ(IJDEP),IZ(EAVF),
     2                           IZ(LG),IZ(OCC),IZ(ZZ),IZ(CIDT),
     3                           IZ(ORBTBF))
C--------------------------------                                       
C    PRINT FORCES IN FILE11     |                                       
C--------------------------------                                       
      if(iprint.ge.1) then
       call codump(iz)                                                   
      endif
      CALL PRN11(IZ(BOX),NDF,CCSD,IZ(CIDT))                             
C                                                                       
      write(6,*) '  '
      RETURN                                                            
      END                                                               
