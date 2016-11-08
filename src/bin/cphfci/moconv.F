      SUBROUTINE MOCONV(SA,NNA,SM,NNM,EAO,U,T)                          
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      INTEGER P,Q,PQ                                                    
      DIMENSION SA(NNA),SM(NNM),EAO(NBFAO,NORBS),U(NORBS),T(NNA)        
      COMMON /DIMAO/ NBFAO,NBATRI                                       
      COMMON /DIMSG/ ISYMRB,NUMIJ,NBF,NMAX,NBFORB,NSYM,NDF,NORBS        
c                                                                       
      A0=0.0D+00                                                        
      A1=1.0D+00                                                        
      A1S2=0.50D+00                                                     
      DO 101 IJ=1,NNA                                                   
 101  T(IJ)=SA(IJ)                                                      
      DO 102 IJ=1,NNM                                                   
 102  SM(IJ)=A0                                                         
C                                                                       
      PQ=0                                                              
      DO 40 P=1,NBFAO                                                   
        DO 5 II=1,NORBS                                                 
 5      U(II)=A0                                                        
        E=A1                                                            
        DO 20 Q=1,P                                                     
          PQ=PQ+1                                                       
          IF (P .EQ. Q) E=A1S2                                          
          F=T(PQ)*E                                                     
          DO 10 J=1,NORBS                                               
 10       U(J)=U(J)+EAO(Q,J)*F                                          
 20     CONTINUE                                                        
        IJ=0                                                            
        DO 30 I=1,NORBS                                                 
          U1I=U(I)                                                      
          DO 25 J=1,I                                                   
            IJ=IJ+1                                                     
            SM(IJ)=SM(IJ)+U1I*EAO(P,J)+U(J)*EAO(P,I)                    
 25       CONTINUE                                                      
 30     CONTINUE                                                        
 40   CONTINUE                                                          
      RETURN                                                            
      END                                                               
