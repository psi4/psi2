      SUBROUTINE BACK1(TOLD,TNEW,XTP,EIV,SQ,MAPIJ,ILAST)                
C                                                                       
      IMPLICIT INTEGER (A-Z)                                            
      INTEGER MAPIJ(NIJAO),ILAST(NBKT)                                  
      REAL*8 TNEW(2),TOLD(2),EIV(NBF,NBFAO),                            
     .       SQ(NBFAO*NBFAO),XTP(NBFAO*NBFAO)                           
C                                                                       
      COMMON /DIMS/ NIJ,NBF,NIJAO,NBFAO,LENR,NST,NSHELL,AOSTRT,BKTSZ    
      COMMON/TAPES/ IOUT,ICHECK,IT30,IT55,IT91,IT95,IT96                
      COMMON /SRT/  NBKT,NGRPIJ,LGRPIJ,NBKT2,LGRP2,IBPT,IAPT,IGRP       
      IAPT=1                                                            
      IBPT=1                                                            
      GRPIJ = NGRPIJ                                                    
C*********************************                                      
C   START MAIN LOOP BUCKETS - 1  *                                      
C*********************************                                      
      DO 10 IBKT = 1,NBKT-1                                             
        IGRP = IBKT                                                     
        CALL TRANDM(TOLD,TNEW,EIV,GRPIJ,XTP,SQ,MAPIJ,TOLD,TNEW,         
     .              ILAST)                                              
   10 CONTINUE                                                          
C*******************************                                        
C     TRANSFORM LAST GROUP     *                                        
C*******************************                                        
      GRPIJ = LGRPIJ                                                    
      IGRP = NBKT                                                       
      CALL TRANDM(TOLD,TNEW,EIV,GRPIJ,XTP,SQ,MAPIJ,TOLD,TNEW,ILAST)     
      RETURN                                                            
      END                                                               
