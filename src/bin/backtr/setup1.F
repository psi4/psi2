      SUBROUTINE SETUP1(AOPS,KTYPE,MAPIJ,IOFF,KLOC)                     
      IMPLICIT INTEGER (A-Z)                                            
      INTEGER AOPS(NSHELL),KTYPE(NSHELL),MAPIJ(NIJAO),KLOC(NSHELL)      
      INTEGER IOFF(NIJAO)                                               
      COMMON /DIMS/ NIJ,NBF,NIJAO,NBFAO,LENR,NST,NSHELL,AOSTRT,BKTSZ    
      COMMON/TAPES/ IOUT,ICHECK,IT30,IT55,IT91,IT95,IT96                
      COMMON /SRT/  NBKT,NGRPIJ,LGRPIJ,NBKT2,LGRP2,IBPT,IAPT,IGRP       
C**************************************************                     
C     MAKE AO PER SHELL ARRAY - AOPS(NSHELL)      *                     
C**************************************************                     
      DO 10 XSHL=1,NSHELL                                               
        MUTYP = KTYPE(XSHL)*(KTYPE(XSHL)+1)/2                           
        AOPS(XSHL) = MUTYP                                              
   10 CONTINUE                                                          
C**************************************************                     
C     MAKE INTEGER OFFSET ARRAY - IOFF(NIJAO)     *                     
C**************************************************                     
      DO 20 II=1,NIJAO                                                  
        IOFF(II)=II*(II-1)/2                                            
   20 CONTINUE                                                          
C****************************************************************       
C     MAKE MAP OF PAIR OF BASIS FUNCTIONS MU,NU TO POSITION     *       
C     IN FILE ORDERED BY PAIRS OF SHELLS                        *       
C****************************************************************       
      NWMUNU = 0                                                        
      DO 160 R = 1,NSHELL                                               
        FMU = KLOC(R)                                                   
        LMU = FMU + AOPS(R) - 1                                         
        DO 130 S = 1,R-1                                                
          FNU = KLOC(S)                                                 
          LNU = FNU + AOPS(S) - 1                                       
          DO 120 MU = FMU,LMU                                           
            DO 110 NU = FNU,LNU                                         
              MUNU = IOFF(MU) + NU                                      
              NWMUNU = NWMUNU + 1                                       
              MAPIJ(MUNU) = NWMUNU                                      
  110       CONTINUE                                                    
  120     CONTINUE                                                      
  130   CONTINUE                                                        
        DO 150 MU = FMU,LMU                                             
          DO 140 NU = FMU,MU                                            
            MUNU = IOFF(MU) + NU                                        
            NWMUNU = NWMUNU + 1                                         
            MAPIJ(MUNU) = NWMUNU                                        
  140     CONTINUE                                                      
  150   CONTINUE                                                        
  160 CONTINUE                                                          
C     WRITE(6,*) 'IOFF'                                                 
C     WRITE(6,*) IOFF                                                   
C                                                                       
      RETURN                                                            
      END                                                               
