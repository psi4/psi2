      SUBROUTINE ORDER1(A,SCR,IOFF,ORBTBF,DROL,CCSD)                    
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION A(NUMIJ),SCR(NUMIJ)                                     
      INTEGER ORBTBF(NORBS),IOFF(NUMIJ),DROL(NORBS)                     
      LOGICAL CCSD                                                      
      COMMON /DIMSG/ ISYMRB,NUMIJ,NBF,NMAX,NBFORB,NSYM,NDF,NORBS        
C                                                                       
      DO 101 IJ=1,NUMIJ                                                 
        SCR(IJ)=A(IJ)                                                   
 101  A(IJ)=0.0D+00                                                     
C                                                                       
      IF(CCSD)THEN                                                      
       IJ=0                                                             
       DO 102 IOR=1,NORBS                                               
        IBF=DROL(ORBTBF(IOR))                                           
        DO 102 JOR=1,IOR                                                
         IJ=IJ+1                                                        
         JBF=DROL(ORBTBF(JOR))                                          
         IBM=MAX0(IBF,JBF)                                              
         JBM=MIN0(IBF,JBF)                                              
         IJBF=IOFF(IBM)+JBM                                             
 102   A(IJ)=SCR(IJBF)                                                  
      ELSE                                                              
       IJ=0                                                             
       DO 103 IOR=1,NORBS                                               
        IBF=ORBTBF(IOR)                                                 
        DO 103 JOR=1,IOR                                                
         IJ=IJ+1                                                        
         JBF=ORBTBF(JOR)                                                
         IBM=MAX0(IBF,JBF)                                              
         JBM=MIN0(IBF,JBF)                                              
         IJBF=IOFF(IBM)+JBM                                             
 103   A(IJ)=SCR(IJBF)                                                  
      ENDIF                                                             
      RETURN                                                            
      END                                                               
