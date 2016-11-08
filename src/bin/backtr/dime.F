      SUBROUTINE Dime(AOPS,mxdim1,mxdim2)                               
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      INTEGER AOPS(NSHELL),aostrt,bktsz                                 
      COMMON /DIMS/ NIJ,NBF,NIJAO,NBFAO,LENR,NST,NSHELL,AOSTRT,BKTSZ    
      COMMON/TAPES/ IOUT,ICHECK,IT30,IT55,IT91,IT95,IT96                
C                                                                       
      WRITE(ICHECK,10)                                                  
   10 FORMAT(/5X,'you are in space')                                    
C*****************************************************************      
C     subroutine to see how much space is needed for sort2.      *      
C*****************************************************************      
      mxdim1 = 0                                                        
      mxdim2 = 0                                                        
      IMU = 1                                                           
      DO 80 IR = 1,NSHELL                                               
        NMU = AOPS(IR)                                                  
        if(nmu.gt.mxdim2) mxdim2 = nmu                                  
        NTRI = (IMU+NMU-1)*(IMU+NMU)/2                                  
        INU = 1                                                         
        DO 70 IS = 1,IR                                                 
          IF(IR.EQ.IS) THEN                                             
            ISIZ = AOPS(IR)*(AOPS(IS) + 1) / 2                          
          ELSE                                                          
            ISIZ = AOPS(IR)*AOPS(IS)                                    
          ENDIF                                                         
          NNU = AOPS(IS)                                                
          LENGTH = ISIZ*NTRI                                            
          if(length.gt.mxdim1) mxdim1 = length                          
          INU = INU + NNU                                               
   70   CONTINUE                                                        
        IMU = IMU + NMU                                                 
   80 CONTINUE                                                          
      mxdim2 = mxdim2**4                                                
      write(6,*) 'MXDIM1 =',mxdim1                                      
      write(6,*) 'MXDIM2 =',mxdim2                                      
      RETURN                                                            
      END                                                               
