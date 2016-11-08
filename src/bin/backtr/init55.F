      SUBROUTINE INIT55(AOPS,LNSR,IJOFF,KLOFF)                          
      IMPLICIT INTEGER (A-Z)                                            
      INTEGER AOPS(NSHELL),LNSR(NST),IJOFF(NST),KLOFF(NST)              
      COMMON/TAPES/ IOUT,ICHECK,IT30,IT55,IT91,IT95,IT96                
      COMMON /DIMS/ NIJ,NBF,NIJAO,NBFAO,LENR,NST,NSHELL,AOSTRT,BKTSZ    
C**********************************************                         
C     MAKE OFFSET ARRAY FOR SHELL PAIR K,L    *                         
C**********************************************                         
      NNST = 0                                                          
      KLOFF(1) = 0                                                      
      DO 20 S=1,NSHELL                                                  
        APS = AOPS(S)                                                   
        DO 10 T=1,S                                                     
          NNST = NNST + 1                                               
          LNSR(NNST) = APS*AOPS(T)                                      
          IF(NNST.NE.NST) KLOFF(NNST+1) = KLOFF(NNST) + LNSR(NNST)      
   10   CONTINUE                                                        
   20 CONTINUE                                                          
C**********************************************                         
C     MAKE OFFSET ARRAY FOR SHELL PAIR I,J    *                         
C**********************************************                         
      SRMN = 0                                                          
      IJOFF(1) = 0                                                      
      DO 40 ST = 1,NST                                                  
        LNST = LNSR(ST)                                                 
        IJOF = 0                                                        
        DO 30 MN = 1,ST                                                 
          SRMN = SRMN + 1                                               
          IJOF = IJOF + LNST*LNSR(MN)                                   
   30   CONTINUE                                                        
      IF(ST.NE.NST) IJOFF(ST+1) = IJOFF(ST) + IJOF                      
   40 CONTINUE                                                          
C*************************************************                      
C     WRITE THE ARRAYS OF POINTERS TO FILE 55    *                      
C*************************************************                      
      CALL WWRITW(IT55,AOPS,NSHELL,1,JUNK)                              
      CALL WWRITW(IT55,IJOFF,NST,JUNK,JUNK)                             
      CALL WWRITW(IT55,KLOFF,NST,JUNK,JUNK)                             
      CALL WWRITW(IT55,LNSR,NST,JUNK,JUNK)                              
c     WRITE(6,*) 'AOPS'                                                 
c     WRITE(6,*) AOPS                                                   
c     WRITE(6,*) 'IJOFF'                                                
c     WRITE(6,*) IJOFF                                                  
c     WRITE(6,*) 'KLOFF'                                                
c     WRITE(6,*) KLOFF                                                  
c     WRITE(6,*) 'LNSR'                                                 
c     WRITE(6,*) LNSR                                                   
C********************************************************               
C     WRITE THE NUMBER OF TPDM ELEMENTS.   AOSTRT IS    *               
C     POINTER FOR THE BEGINNING OF TPDM MATRIX.         *               
C********************************************************               
      TPDMEL = IJOFF(NST) + IJOF                                        
      WRITE(IOUT,*) 'THE NUMBER OF TPDM ELEMENTS IS ',TPDMEL            
      AOSTRT = JUNK                                                     
      WRITE(6,*) 'AOSTRT',AOSTRT                                        
C                                                                       
      RETURN                                                            
      END                                                               
