      SUBROUTINE ORDER2(A,SCR,ORBTBF)                                   
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION A(NBFAO,NORBS),SCR(NBFAO,NORBS)                         
      INTEGER ORBTBF(NORBS)                                             
      COMMON /DIMAO/ NBFAO,NBATRI                                       
      COMMON /DIMSG/ ISYMRB,NUMIJ,NBF,NMAX,NBFORB,NSYM,NDF,NORBS        
C                                                                       
      DO 101 I=1,NBFAO                                                  
        DO 101 J=1,NORBS                                                
          SCR(I,J)=A(I,J)                                               
 101  A(I,J)=0.0D+00                                                    
C                                                                       
      DO 102 IOR=1,NORBS                                                
        IBF=ORBTBF(IOR)                                                 
        DO 102 J=1,NBFAO                                                
 102  A(J,IOR)=SCR(J,IBF)                                               
C                                                                       
      RETURN                                                            
      END                                                               
