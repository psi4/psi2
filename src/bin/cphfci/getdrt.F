      SUBROUTINE GETDRT(KADD,LADD,IJADD,IJGRP,ORBTBF,IOUT,ORBSYM,JCODE, 
     1NGROUP)                                                           
      IMPLICIT INTEGER (A-Z)                                            
      COMMON /TAPES/ NOUT,ICHK,IT30,IT42,IT47,IT52,IT58,IT94,IT96,IT98  
      COMMON /DIMSG/ ISYMRB,NUMIJ,NBF,NMAX,NBFORB,NSYM,NDF,NORBS        
      COMMON /FUNCS/ IOCC,JOCC,KOCC,NTYPES,NBSET,NAT                    
      COMMON /OCCS/ NSPE,NSA,NSB,NALP,NOP,NDOC,LFERM,NUOC,IJIND,IJIND3  
      COMMON /IEXT/ IPRINT,IGRP,ICIDIP                                  
      DIMENSION KADD(ISYMRB),LADD(ISYMRB),IJADD(NUMIJ),IJGRP(NUMIJ)     
      DIMENSION JCODE(NORBS),IOUT(NBF),ORBTBF(NORBS),ORBSYM(NORBS)      
C                                                                       
      WRITE(ICHK,2222)                                                  
 2222 FORMAT(/'  IN GETDRT')                                            
      END=48+INTOWP(1)                                                  
      CALL WREADW(IT58,KADD,ISYMRB,END,END)                             
      CALL WREADW(IT58,LADD,ISYMRB,END,END)                             
      CALL WREADW(IT58,IJADD,NUMIJ,END,END)                             
      CALL WREADW(IT58,IJGRP,NUMIJ,END,END)                             
      END=END+4*NORBS+NGROUP                                            
      CALL WREADW(IT58,IOUT,NBF,END,END)                                
      CALL WREADW(IT58,ORBSYM,NORBS,END,END)                            
      END=END+NBF                                                       
      CALL WREADW(IT58,JCODE,NBF,END,END)                               
      CALL WREADW(IT58,ORBTBF,NORBS,END,END)                            
COUT  WRITE(*,*)'ORBTBF FROM DRT'                                       
COUT  WRITE(*,*) ORBTBF                                                 
      IF (IAND(IPRINT,1) .NE. 0) THEN                                   
        WRITE(NOUT,*) '  KADD:'                                         
        WRITE(NOUT,710)  KADD                                           
        WRITE(NOUT,*) '  LADD:'                                         
        WRITE(NOUT,710)  LADD                                           
        WRITE(NOUT,*) '  IJADD:'                                        
        WRITE(NOUT,710)  IJADD                                          
        WRITE(NOUT,*) '  IJGRP:'                                        
        WRITE(NOUT,710)  IJGRP                                          
        WRITE(NOUT,*) '  IOUT:'                                         
        WRITE(NOUT,710)  IOUT                                           
        WRITE(NOUT,*) '  ORBSYM:'                                       
        WRITE(NOUT,710)  ORBSYM                                         
        WRITE(NOUT,*) '  JCODE:'                                        
        WRITE(NOUT,710)  JCODE                                          
        WRITE(NOUT,*) '  ORBTBF:'                                       
        WRITE(NOUT,710)  ORBTBF                                         
      END IF                                                            
 710  FORMAT(16I5)                                                      
C                                                                       
      RETURN                                                            
      END                                                               
