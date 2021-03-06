      SUBROUTINE FLINX(A,IDIM,IN,IM,DET)                                
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C     EXTENDED DUMMY A                                                  
C                                                                       
C     LINEAR SIMULTANEOUS EQUATION                                      
C                                                                       
C     A(IN*IN) * X(IN*IM) = B(IN*IM)                                    
C                                                                       
C     A & B SHOULD BE STORED ON A(IN*(IN+IM))                           
C     SOLUTION X WILL BE STORED ON B PART IN DIMENSION A.               
C                                                                       
      REAL*8 DET,D                                                     
      DIMENSION A(IDIM,1)                                               
      DATA ZERO,ONE / 0.0D+00 , 1.0D+00 /                               
C                                                                       
      N=IN                                                              
      NR=IM                                                             
      JMAX=N+NR                                                         
      SIGN=ONE                                                          
C M IS THE STAGE OF ELIMINATION                                         
      DO 49 M=1,N                                                       
        TEMP=ZERO                                                       
        DO 41 I=M,N                                                     
          IF(M.GT.1)A(I,M)=A(I,M)-DOT(A(I,1),IDIM,A(1,M),1,M-1)         
          AVAL=A(I,M)                                                   
          IF(DABS(AVAL).LE.TEMP)GOTO 41                                 
          TEMP=DABS(AVAL)                                               
          IMAX=I                                                        
 41     CONTINUE                                                        
        IF(TEMP.LE.ZERO)GOTO 999                                        
        IF(IMAX.EQ.M)GOTO 45                                            
        SIGN=-SIGN                                                      
        DO 44 J=1,JMAX                                                  
          STOR=A(M,J)                                                   
          A(M,J)=A(IMAX,J)                                              
          A(IMAX,J)=STOR                                                
 44     CONTINUE                                                        
 45     CONTINUE                                                        
        JJ=M+1                                                          
        IF(JJ.GT.JMAX)GOTO 49                                           
        IF(M.GT.1)GOTO 47                                               
        DO 46 J=JJ,JMAX                                                 
          A(1,J)=A(1,J)/A(1,1)                                          
 46     CONTINUE                                                        
        D=A(1,1)                                                        
        GOTO 49                                                         
 47     CONTINUE                                                        
        DO 48 J=JJ,JMAX                                                 
          A(M,J)=(A(M,J)-DOT(A(M,1),IDIM,A(1,J),1,M-1))/A(M,M)          
 48     CONTINUE                                                        
CTJL    D=D*A(M,M)                                                      
 49   CONTINUE                                                          
      IF(NR.EQ.0) RETURN                                                
      DO 59 I=1,NR                                                      
        NPI=N+I                                                         
        DO 58 K=2,N                                                     
          J=N+1-K                                                       
          A(J,NPI)=A(J,NPI)-DOT(A(J,J+1),IDIM,A(J+1,NPI),1,K-1)         
 58     CONTINUE                                                        
 59   CONTINUE                                                          
C***  IF(DABS(D).GE.1.0D+36) D=1.0D+36                                  
C***  IF(DABS(D).LE.1.0D-36) D=1.0D-36                                  
CTJL  DET=D*SIGN                                                        
      RETURN                                                            
C ON ZERO PIVOT, SET DET=0.AND RETURN TO CALLING PROGRAM NOV 1972       
 999  DET=ZERO                                                          
      RETURN                                                            
      END                                                               
