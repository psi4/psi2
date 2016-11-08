      SUBROUTINE REGIST                                                 
      IMPLICIT INTEGER (A-Z)                                            
      REAL*8 FOCC                                                       
      COMMON /TAPES/ NOUT,ICHK,IT30,IT42,IT47,IT52,IT58,IT94,IT96,IT98  
      COMMON /DIMSG/ SYMORB,NUMIJ,NBF,NMAX,NBFORB,NSYM,NDF,NORBS        
      COMMON /POINT/ FOCC(5),IPOINT(5),JPOINT(5),NSORB(5)               
      COMMON /FUNCS/ IOCC,JOCC,KOCC,NTYPES,NBSET,NAT                    
      COMMON /OCCS/ NSPE,NSA,NSB,NALP,NOP,NDOC,LFERM,NUOC,IJIND,IJIND3  
C                                                                       
CTJL  WRITE(ICHK,2222)                                                  
C2222 FORMAT(/'  IN REGIST')                                            
      IPOINT(1)=IOCC+1                                                  
      JPOINT(1)=1                                                       
      DO 10 I=2,NTYPES+1                                                
        JPOINT(I)=JPOINT(I-1)+NSORB(I-1)                                
 10   IPOINT(I)=IPOINT(I-1)+NSORB(I)                                    
      IT1=NTYPES+1                                                      
      WRITE(ICHK,101)                                                   
      WRITE(ICHK,99) (IPOINT(I), I=1,IT1)                               
      WRITE(ICHK,102)                                                   
      WRITE(ICHK,99) (JPOINT(I), I=1,IT1)                               
  99  FORMAT(2X,5I5)                                                    
 101  FORMAT('  THE IPOINT ARRAY:')                                     
 102  FORMAT('  THE JPOINT ARRAY:')                                     
      RETURN                                                            
      END                                                               
