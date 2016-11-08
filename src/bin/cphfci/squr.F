      SUBROUTINE SQUR(T,SQ,N)                                           
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION SQ(2),T(2)                                              
      IJ=0                                                              
      II=0                                                              
      DO 20 I=1,N                                                       
       JJ=0                                                             
       DO 10 J=1,I                                                      
        IJ=IJ+1                                                         
        SQ(II+J)=T(IJ)                                                  
        SQ(JJ+I)=T(IJ)                                                  
10     JJ=JJ+N                                                          
20    II=II+N                                                           
      RETURN                                                            
      END                                                               
