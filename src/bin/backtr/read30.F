      SUBROUTINE READ30(KTYPE,KLOC,I30)                                 
      IMPLICIT INTEGER (A-Z)                                            
      INTEGER KTYPE(NSHELL),KLOC(NSHELL),I30(200)                       
      COMMON /DIMS/ NIJ,NBF,NIJAO,NBFAO,LENR,NST,NSHELL,AOSTRT,BKTSZ    
      COMMON/TAPES/ IOUT,ICHECK,IT30,IT55,IT91,IT95,IT96                
C*********************************************************************  
C     THE POINTERS FOR THE BASIS SET INFO ON TAPE 30 START AT 301    *  
C     AND EXTEND FOR 200 INTEGER WORDS.                              *  
C*********************************************************************  
C     WRITE(6,*) 'IT30,301',IT30,301                                    
      CALL WREADW(IT30,I30,200,301,JUNK)                                
C                                                                       
      KTPPT  = I30( 9)                                                  
      KLOCPT = I30(11)                                                  
C     WRITE(6,*) 'IT30,KTPPT',IT30,KTPPT                                
C     WRITE(6,*) 'IT30,KLOCPT',IT30,KLOCPT                              
      CALL WREADW(IT30,KTYPE,NSHELL,KTPPT,JUNK)                         
      CALL WREADW(IT30,KLOC,NSHELL,KLOCPT,JUNK)                         
C                                                                       
      RETURN                                                            
      END                                                               
