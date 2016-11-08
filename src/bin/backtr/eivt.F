      SUBROUTINE EIVT(EIV,WORK,EIVTR,NBF,NBFAO)                         
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION EIV(NBFAO,NBF),WORK(NBF,NBFAO),EIVTR(NBF,NBFAO)         
C*************************************                                  
C     TRANSPOSE AO-MO EIGENVECTOR    *                                  
C*************************************                                  
      DO 20 I=1,NBF                                                     
        DO 10 J=1,NBFAO                                                 
          WORK(I,J) = EIV(J,I)                                          
10      CONTINUE                                                        
20    CONTINUE                                                          
      DO 40 J=1,NBFAO                                                   
        DO 30 I=1,NBF                                                   
          EIVTR(I,J) = WORK(I,J)                                        
30      CONTINUE                                                        
40    CONTINUE                                                          
      RETURN                                                            
      END                                                               
