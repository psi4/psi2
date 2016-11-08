      SUBROUTINE MATOUT(A,NAD,NBD,M,N,IOUT)                             
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      CHARACTER*6 LINE(13)                                              
      DIMENSION A(NAD,NBD)                                              
    1 FORMAT(3X,6(7X,I5))                                              
    2 FORMAT(3X,13A6)                                                   
    3 FORMAT(2X,I3,2X,6F12.7)                                          
    4 FORMAT(/)                                                         
C                                                                       
      DATA LINE / 13*'------' /                                         
      II=0                                                              
      JJ=0                                                              
  200 II=II+1                                                           
      JJ=JJ+1                                                           
      KK=6*JJ                                                          
      NN=N                                                              
      IF(N.GT.KK) NN=KK                                                 
      LL=2*(NN-II+1)+1                                                  
      WRITE (IOUT,1) (I,I=II,NN)                                        
      WRITE (IOUT,2) (LINE(I),I=1,LL)                                   
      DO 101 I=1,M                                                      
      WRITE (IOUT,3) I,(A(I,J),J=II,NN)                                 
  101 CONTINUE                                                          
      IF(N.LE.KK) GO TO 201                                             
      WRITE (IOUT,4)                                                    
      II=KK                                                             
      GO TO 200                                                         
  201 RETURN                                                            
      END                                                               
