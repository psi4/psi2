      SUBROUTINE INITG(IOFF,IJOFF,IJBKT)                                
      IMPLICIT INTEGER (A-Z)                                            
      REAL*8 VAL,VAL2                                                   
C                                                                       
      DIMENSION IOFF(NIJ),IJOFF(NIJ),IJBKT(NIJ)                         
      COMMON /ORD/ IFAC1,IFAC2                                          
      COMMON/TAPES/ IOUT,ICHECK,IT30,IT52,IT53,IT58,IT91,IT93           
      COMMON /DIMS/ ISYMRB,NUMIJ,NIJ,NBF,NIJAO,NBFAO,NORBS,NMAX,NSYMHF  
      COMMON /SRT/ VAL,IADR,IBKT,NBKT,NGRPIJ,LENR                       
      COMMON /PRNT/ IPRINT                                              
C************************************************                       
C     MAKE INTEGER OFFSET ARRAY - IOFF(NIJ)     *                       
C************************************************                       
      DO 5 II=1,NIJ                                                     
        IOFF(II) = II*(II-1)/2                                          
    5 CONTINUE                                                          
C********************************************************************   
C     MAKE IJBKT ARRAY - WHICH BUCKET BASIS FUNCTION PAIR I,J IS IN *   
C     MAKE IJOFF ARRAY - OFFSET POSITION FOR PAIR I,J               *   
C********************************************************************   
      WRITE(ICHECK,10)                                                  
   10 FORMAT(' IN INITG ')                                              
      DO 20 IJ=1,NIJ                                                    
        IJBKT(IJ)=IJ/NGRPIJ                                             
        IF (NGRPIJ*IJBKT(IJ) .NE. IJ) IJBKT(IJ)=IJBKT(IJ)+1             
   20 CONTINUE                                                          
      DO 30 IJ=1,NIJ                                                    
        IJOFF(IJ)=(IJ-(IJBKT(IJ)-1)*NGRPIJ-1)*NIJ                       
   30 CONTINUE                                                          
      IF(IPRINT.GE.2) THEN                                              
        WRITE(6,*) 'IJBKT'                                              
        WRITE(6,35) (IJBKT(N),N=1,NIJ)                                  
        WRITE(6,*) 'IJOFF'                                              
        WRITE(6,35) (IJOFF(N),N=1,NIJ)                                  
   35 FORMAT(7I10)                                                      
      ENDIF                                                             
C                                                                       
      RETURN                                                            
      END
