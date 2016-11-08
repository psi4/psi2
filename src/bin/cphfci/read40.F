      SUBROUTINE READ40(IOFF,DROL,ORBTBF,EIG,OCC,ALPHA,BETA,            
     1EPS,VEC,SCR,CCSD,itcscf,nsa,nsb,ci1,ci2)                          
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      INTEGER IOFF(NUMIJ),ORBTBF(NORBS),DROL(NORBS),IREP(2)             
      DIMENSION EIG(NORBS),OCC(NORBS),ALPHA(NUMIJ),BETA(NUMIJ)          
      DIMENSION EPS(NUMIJ),VEC(NBFAO,NORBS),SCR(NUMIJ)                  
      LOGICAL CCSD                                                      
      EQUIVALENCE (IREP(1),REPN)                                        
      COMMON /CALIF/ LPARA(1024),APARA(1024)                            
      COMMON /LIMIT/ MAXBAS,MAXBUF,MAXBF2,MAXBF4                        
      COMMON /LOCAT/ LOCS(1024)                                         
      COMMON /MFSEC/ MFILE,NSECT                                        
      COMMON /REN/ REPNUC                                               
      COMMON /TAPES/ NOUT,ICHK,IT30,IT42,IT47,IT52,IT58,IT94,IT96,IT98  
      COMMON /DIMAO/ NBFAO,NBATRI                                       
      COMMON /DIMSG/ ISYMRB,NUMIJ,NBF,NMAX,NBFORB,NSYM,NDF,NORBS        
      COMMON /IEXT/ IPRINT,IGRP,ICIDIP                                  
C                                                                       
      REPNUC=APARA(1)                                                   
      WRITE(ICHK,3333)REPNUC                                            
      WRITE(NOUT,3333) REPNUC                                           
 3333 FORMAT(/'  THE NUCLEAR REPULSION ENERGY FROM TAPE40 =',F20.10)    
C                                                                       
CTPH  READ IN EIGENVALUES                                               
C                                                                       
      CALL MREAD(EIG,10)                                                
      IF (IAND(IPRINT,2) .NE. 0) THEN                                   
        WRITE(NOUT,101) 'EIG  '                                         
        WRITE(NOUT,102) (EIG(III),III=1,NORBS)                          
      END IF                                                            
C                                                                       
CTPH  READ IN EPS (SCF MO LAGRANGIAN), TEMPORARILY INTO ALPHA           
C                                                                       
      CALL MREAD(ALPHA,25)                                              
      IF (IAND(IPRINT,2) .NE. 0) THEN                                   
        WRITE(NOUT,101) 'EPS  '                                         
        CALL PRINT(ALPHA,NUMIJ,NORBS,NOUT)                              
      END IF                                                            
C                                                                       
CTPH  CHANGE EPS TO SORTED SCF ORDER USING DROL                         
C                                                                       
      CALL MREAD(DROL,4)                                                
      IF (IAND(IPRINT,2) .NE. 0) THEN                                   
        WRITE(NOUT,101) 'DROL '                                         
        WRITE(NOUT,*) (DROL(III),III=1,NORBS)                           
      END IF                                                            
      DO 20 I=1,NORBS                                                   
        M=DROL(I)                                                       
        DO 10 J=1,I                                                     
          N=DROL(J)                                                     
          IJ=IOFF(I)+J                                                  
          MN=IOFF(MAX(M,N))+MIN(M,N)                                    
          EPS(MN)=ALPHA(IJ)                                             
   10   CONTINUE                                                        
   20 CONTINUE                                                          
      IF (IAND(IPRINT,2) .NE. 0) THEN                                   
        WRITE(NOUT,101) 'EPS  '                                         
        CALL PRINT(EPS,NUMIJ,NORBS,NOUT)                                
      END IF                                                            
C                                                                       
CTPH  READ IN ALPHA                                                     
C                                                                       
      CALL MREAD(ALPHA,40)                                              
      IF (IAND(IPRINT,2) .NE. 0) THEN                                   
        WRITE(NOUT,101) 'ALPHA'                                         
        CALL PRINT(ALPHA,NUMIJ,NORBS,NOUT)                              
      END IF                                                            
C                                                                       
CTPH  READ IN AO-MO EIGENVECTOR (DRT ORDERED FROM MASTER FILE)          
C                                                                       
      CALL MREAD(VEC,42)                                                
      IF (IAND(IPRINT,2) .NE. 0) THEN                                   
        WRITE(NOUT,101) 'VEC  '                                         
        CALL MATOUT(VEC,NBFAO,NORBS,NBFAO,NORBS,NOUT)                   
      END IF                                                            
C                                                                       
      CALL MREAD(OCC,11)                                                
      IF (IAND(IPRINT,2) .NE. 0) THEN                                   
        WRITE(NOUT,101) 'OCC  '                                         
        WRITE(NOUT,102) (OCC(III),III=1,NORBS)                          
      END IF                                                            
C                                                                       
C RE-ORDER TO DRT ORDERING                                              
      CALL ORDER(EIG,SCR,ORBTBF,DROL,CCSD)                              
      IF (IAND(IPRINT,2) .NE. 0) THEN                                   
        WRITE(NOUT,103) 'EIG  '                                         
        WRITE(NOUT,102) (EIG(III),III=1,NORBS)                          
      END IF                                                            
C                                                                       
      CALL ORDER(OCC,SCR,ORBTBF,DROL,CCSD)                              
ctph  if(itcscf.eq.1) then
ctph   ci1 =  dsqrt(occ(nsb)*0.5d0)
ctph   ci2 = -dsqrt(occ(nsa)*0.5d0)
ctph   write(6,*) 'CI1=',ci1,'CI2=',ci2
ctph  endif
      IF (IAND(IPRINT,2) .NE. 0) THEN                                   
        WRITE(NOUT,103) 'OCC  '                                         
        WRITE(NOUT,102) (OCC(III),III=1,NORBS)                          
      END IF                                                            
C                                                                       
      CALL ORDER1(ALPHA,SCR,IOFF,ORBTBF,DROL,CCSD)                      
      IF (IAND(IPRINT,2) .NE. 0) THEN                                   
        WRITE(NOUT,103) 'ALPHA'                                         
        CALL PRINT(ALPHA,NUMIJ,NORBS,NOUT)                              
      END IF                                                            
C                                                                       
CTPH  READ IN BETA                                                      
C                                                                       
      CALL MREAD(BETA,41)                                               
      IF (IAND(IPRINT,2) .NE. 0) THEN                                   
        WRITE(NOUT,101) 'BETA '                                         
        CALL PRINT(BETA,NUMIJ,NORBS,NOUT)                               
      END IF                                                            
      CALL ORDER1(BETA,SCR,IOFF,ORBTBF,DROL,CCSD)                       
      if(itcscf.eq.1) then
       beta(numij-1) = ci1*ci2
       beta(numij) = 0.0d0
      endif
      IF (IAND(IPRINT,2) .NE. 0) THEN                                   
        WRITE(NOUT,103) 'BETA '                                         
        CALL PRINT(BETA,NUMIJ,NORBS,NOUT)                               
      END IF                                                            
C                                                                       
 101  FORMAT(3X,A5,'  ARRAY AS READ FROM FILE40:')                      
 102  FORMAT(5F15.6)                                                    
 103  FORMAT(3X,A5,'  ARRAY AFTER RE-ORDERING TO DRT:')                 
 104  FORMAT(10I5)                                                      
      RETURN                                                            
      END                                                               
