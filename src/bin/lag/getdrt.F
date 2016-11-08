      SUBROUTINE GETDRT(KADD,LADD,IJADD,IJGRP,ORBTBF,JOUT,ORBSYM,       
     .                  END,NGROUP)                                     
      IMPLICIT INTEGER (A-Z)                                            
      COMMON/TAPES/ IOUT,ICHECK,IT30,IT52,IT53,IT58,IT91,IT93           
      COMMON /DIMS/ ISYMRB,NUMIJ,NIJ,NBF,NIJAO,NBFAO,NORBS,NMAX,NSYMHF  
      COMMON /PRNT/ IPRINT                                              
      DIMENSION KADD(ISYMRB),LADD(ISYMRB),IJADD(NUMIJ),IJGRP(NUMIJ)     
      DIMENSION JOUT(NBF),ORBTBF(NORBS),ORBSYM(NORBS)                   
C                                                                       
      CALL WREADW(IT58,KADD,ISYMRB,END,END)                             
      CALL WREADW(IT58,LADD,ISYMRB,END,END)                             
      CALL WREADW(IT58,IJADD,NUMIJ,END,END)                             
      CALL WREADW(IT58,IJGRP,NUMIJ,END,END)                             
      END = END + 4*NORBS + NGROUP                                      
      CALL WREADW(IT58,JOUT,NBF,END,END)                                
      CALL WREADW(IT58,ORBSYM,NORBS,END,END)                            
      END = END + 2*NBF                                                 
      CALL WREADW(IT58,ORBTBF,NORBS,END,END)                            
      IF(IPRINT.GE.2) THEN                                              
        WRITE(6,*) 'KADD'                                               
        WRITE(6,10) (KADD(I),I=1,ISYMRB)                                
        WRITE(6,*) 'LADD'                                               
        WRITE(6,10) (LADD(I),I=1,ISYMRB)                                
        WRITE(6,*) 'IJADD'                                              
        WRITE(6,10) (IJADD(I),I=1,NUMIJ)                                
        WRITE(6,*) 'IJGRP'                                              
        WRITE(6,10) (IJGRP(I),I=1,NUMIJ)                                
        WRITE(6,*) 'JOUT'                                               
        WRITE(6,10) (JOUT(I),I=1,NBF)                                   
        WRITE(6,*) 'ORBSYM'                                             
        WRITE(6,10) (ORBSYM(I),I=1,NORBS)                               
        WRITE(6,*) 'ORBTBF'                                             
        WRITE(6,10) (ORBTBF(I),I=1,NORBS)                               
   10   FORMAT(7I10)                                                    
      ENDIF                                                             
C                                                                       
      RETURN                                                            
      END                                                               
