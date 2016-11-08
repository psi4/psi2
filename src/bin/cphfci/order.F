      SUBROUTINE ORDER(A,SCR,ORBTBF,DROL,CCSD)                          
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      INTEGER ORBTBF(NORBS),DROL(NORBS)                                 
      DIMENSION A(NORBS),SCR(NORBS)                                     
      LOGICAL CCSD                                                      
      COMMON /DIMSG/ ISYMRB,NUMIJ,NBF,NMAX,NBFORB,NSYM,NDF,NORBS        
C                                                                       
      DO 101 I=1,NORBS                                                  
        SCR(I)=A(I)                                                     
 101  A(I)=0.0D+00                                                      
C                                                                       
      IF (CCSD)THEN                                                     
       DO 102 IOR=1,NORBS                                               
        IBF=DROL(ORBTBF(IOR))                                           
 102   A(IOR)=SCR(IBF)                                                  
      ELSE                                                              
       DO 103 IOR=1,NORBS                                               
        IBF=ORBTBF(IOR)                                                 
 103   A(IOR)=SCR(IBF)                                                  
      ENDIF                                                             
      RETURN                                                            
      END                                                               
