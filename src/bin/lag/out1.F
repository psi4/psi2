      SUBROUTINE OUT1(F,X,ia47)                                         
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      integer ia47(192)
      REAL*8 F(NBF,NBF),X(NIJ)                                          
C                                                                       
      COMMON /DIMS/ ISYMRB,NUMIJ,NIJ,NBF,NIJAO,NBFAO,NORBS,NMAX,NSYMHF  
      COMMON/TAPES/ IOUT,ICHECK,IT30,IT52,IT53,IT58,IT91,IT93           
      COMMON /PRNT/ IPRINT                                              
C                                                                       
      it47 = 47
C*******************************************************************    
C     WRITE ONEPDM IN THE MO BASIS TO MASTER FILE 40 (TRIANG).     *    
C*******************************************************************    
      IJ = 0                                                            
      DO 20 IOR=1,NORBS                                                 
        DO 10 JOR=1,IOR                                                 
          IJ = IJ + 1                                                   
          X(IJ)=F(IOR,JOR)                                              
   10   CONTINUE                                                        
   20 CONTINUE                                                          
      CALL MWRIT(X,35)                                                  
ctph
ctph  mopdpt = ia47(119) + intowp(nbf*nbf) + 2*intowp(nijao)
ctph  mopdln = intowp(nij)
ctph  call wwritw(it47,mopdpt,1,123,junk)
ctph  call wwritw(it47,mopdln,1,173,junk)
      mopdpt = ia47(123)
      mopdln = ia47(173)
      call wwritw(it47,x,mopdln,mopdpt,junk)
ctph  
      IF(IPRINT.GE.2) THEN                                              
        WRITE(6,30)                                                     
   30   FORMAT(/1X,'ONE PDM IN MO BASIS AS WRITTEN TO 40'/)             
        CALL PRINT(X,NIJ,NBF,IOUT)                                      
      END IF                                                            
C**************                                                         
C     DONE    *                                                         
C**************                                                         
      RETURN                                                            
      END                                                               
