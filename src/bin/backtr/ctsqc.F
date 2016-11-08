      SUBROUTINE CTSQC(EIV,SQ,WORK,TRI,NBF,NBFAO,NBFTMP,NTRI)           
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION EIV(NBF,NBFAO),SQ(NBF,NBF),WORK(NBF,NBFTMP),TRI(NTRI)   
C************************************************************           
C     TRANSFORM A SQUARE ARRAY WITH C TRANSPOSE X SQ X C    *           
C     AND PUT INTO A PACKED TRIANGULAR ARRAY TRI.           *           
C************************************************************           
cets072091      DO 30 J=1,NBFTMP                                                
c        DO 20 I=1,NBF                                                     
c          sum = 0.0D0                                                   
c          DO 10 K=1,NBF                                                 
c            sum = sum + EIV(K,J)*SQ(K,I)                                
c10        CONTINUE                                                      
c          work(i,j) = sum                                               
c20      CONTINUE                                                          
c30    CONTINUE                                                        
      call mmult(eiv,nbf,1,sq,nbf,0,work,nbf,1,nbftmp,nbf,nbf,0)
      IJ = 0                                                            
      DO 60 I=1,NBFTMP                                                  
        DO 50 J=1,I                                                     
          IJ = IJ + 1                                                   
          sum = 0.0D0                                                   
          DO 40 K=1,NBF                                                 
            sum = sum + EIV(K,I)*WORK(K,J)                              
40        CONTINUE                                                      
          tri(ij) = sum                                                 
50      CONTINUE                                                        
60    CONTINUE                                                          
      RETURN                                                            
      END                                                               
