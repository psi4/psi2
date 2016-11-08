      SUBROUTINE SORT2(TPMOLD,TPMNEW,AOPS,ISIZE,IREAD,IOFF,MAPIJ,mxdim1,
     .                 mxdim2)                                          
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      INTEGER AOSTRT,AOPS(NSHELL),BKTSZ                                 
      COMMON/TAPES/ IOUT,ICHECK,IT30,IT55,IT91,IT95,IT96                
      COMMON /SRT/  NBKT,NGRPIJ,LGRPIJ,NBKT2,LGRP2,IBPT,IAPT,IGRP       
      COMMON /DIMS/ NIJ,NBF,NIJAO,NBFAO,LENR,NST,NSHELL,AOSTRT,BKTSZ    
      DIMENSION TPMOLD(mxdim1),TPMNEW(mxdim2),MAPIJ(NIJAO),             
     .          ISIZE(NST),IREAD(NST+1),IOFF(NIJAO)                     
C                                                                       
      WRITE(ICHECK,10)                                                  
   10 FORMAT(/5X,'IN SORT2')                                            
C*****************************************************************      
C     MAKE ARRAYS TO HELP DECIDE HOW MUCH OF FILE 96 TO READ     *      
C     IRSMX IS MAXIMUM ARRAYS OF DIM NIJAO THAT CAN BE IN CORE   *      
C*****************************************************************      
croag
croag   improved an error message a bit.
croag   Tue Dec  1 12:14:26 EST 1992
croag
      IRS = 0                                                           
      IPOS = 1                                                          
      IPOST = 1                                                         
      IREAD(1) = 1                                                      
      NREAD = 1                                                         
      IRSMX = LENR/NIJAO                                                
C     WRITE(6,*) 'IRSMX',IRSMX                                          
C     WRITE(6,*) 'LENR',LENR                                            
C     WRITE(6,*) 'IOFF',IOFF                                            
C     IOFSET = 0                                                        
      DO 30 IR = 1,NSHELL                                               
        DO 20 IS = 1,IR                                                 
          IRS = IRS + 1                                                 
          IF(IR.EQ.IS) THEN                                             
            ISIZ = AOPS(IR)*(AOPS(IS) + 1) / 2                          
          ELSE                                                          
            ISIZ = AOPS(IR)*AOPS(IS)                                    
          ENDIF                                                         
          ISIZE(IRS) = ISIZ                                             
          IPOS = IPOS + ISIZ                                            
C         WRITE(6,*) 'IRS,ISIZ,IPOS',IRS,ISIZ,IPOS                      
C         WRITE(6,*) 'IPOST',IPOST                                      
          IF(ISIZ.GT.IRSMX) THEN                                        
            write(6,*) 'isiz = ',isiz,' irsmx = ',irsmx
croag            WRITE(*,*) 'ERROR: CANNOT RUN SORT2 FOR THIS BASIS SIZE YET'
croag  -- new output follows
	    write(6,*) 'IRSMX IS MAXIMUM ARRAYS OF DIM NIJAO',
     &                  ' THAT CAN BE IN CORE'
	    write(6,*) 'ISIZ IS THE NECESSARY ARRAYS OF DIM',
     &			' NIJAO'  
	    write(6,*) 'TO CORRECT THIS ERROR ALLOCATE MORE MEMORY'
	    write(6,*) 'REQUIRED MEM= CURRENT MEM * ISIZ / IRSMX '
            STOP 1                                                     
          ENDIF                                                         
          IF(IPOS-IPOST.GT.IRSMX) THEN                                  
            IPOST = IPOS - ISIZ                                         
            NREAD = NREAD + 1                                           
            IREAD(NREAD) = IPOST                                        
C         WRITE(6,*) 'IPOST',IPOST                                      
          ENDIF                                                         
   20   CONTINUE                                                        
   30 CONTINUE                                                          
      IREAD(NREAD+1) = IPOS                                             
C*************************************                                  
C     LOOP OVER NUMBER OF SHELLS     *                                  
C*************************************                                  
C     WRITE(6,*) 'IREAD',IREAD                                          
      IBPT = AOSTRT                                                     
      IPT96 = 1                                                         
      NREAD = 1                                                         
      IRS = 0                                                           
      IRSIN = 0                                                         
      IMU = 1                                                           
      ICNT = 0                                                          
      DO 80 IR = 1,NSHELL                                               
        NMU = AOPS(IR)                                                  
        NTRI = (IMU+NMU-1)*(IMU+NMU)/2                                  
        INU = 1                                                         
        DO 70 IS = 1,IR                                                 
          NNU = AOPS(IS)                                                
          IRS = IRS + 1                                                 
          ISIZ1 = ISIZE(IRS)                                            
C********************************                                       
C     READ IN PART OF FILE96    *                                       
C********************************                                       
            IOFSET = NTRI * IRSIN                                       
            LENGTH = ISIZ1*NTRI                                         
C     WRITE(6,*) 'IT96,IPT96,LENGTH',IT96,IPT96,LENGTH                  
            CALL WREADW(IT96,TPMOLD,INTOWP(LENGTH),IPT96,IPTL)          
C     WRITE(6,*) 'TPMOLD IN IT96'                                       
C     WRITE(6,*) (TPMOLD(I),I=1,LENGTH)                                 
            IPT96 = IPTL                                                
C     WRITE(6,*) 'IOFSET=',IOFSET                                       
C*******************************************************                
C     SORT FILE 96 INTO ARRAYS OF PROPER DIMENSION     *                
C*******************************************************                
          IRSIN  = IRSIN  + ISIZ1                                       
          ISTR = 1                                                      
          IRHO = 1                                                      
          DO 60 IT = 1, IR                                              
            ITT = IT                                                    
            NRHO = AOPS(IT)                                             
            IF(IT.EQ.IR) ITT = IS                                       
            ISIG = 1                                                    
            DO 50 IV = 1, ITT                                           
              NSIG = AOPS(IV)                                           
              DO 40 MU = IMU,IMU + NMU - 1                              
                DO 35 NU = INU,INU + NNU - 1                            
                  MUNU = IOFF(MAX(MU,NU)) + MIN(MU,NU)                  
                  NWMUNU = NTRI*( MAPIJ(MUNU) - 1 ) - IOFSET            
                  DO 34 JRHO = IRHO,IRHO + NRHO -1                      
                    DO 33 JSIG = ISIG,ISIG + NSIG - 1                   
                      KLOFF = IOFF(MAX(JRHO,JSIG)) + MIN(JRHO,JSIG)     
C     WRITE(6,32) MU,NU,JRHO,JSIG,NWMUNU,KLOFF                          
   32 FORMAT('MU NU RHO SIG IJOFF KLOFF',5X,6I5)                        
                      ICNT = ICNT + 1                                   
                      TPMNEW(ICNT) = TPMOLD( NWMUNU + MAPIJ(KLOFF) )    
   33               CONTINUE                                            
   34             CONTINUE                                              
   35           CONTINUE                                                
   40         CONTINUE                                                  
              CALL WWRITW(IT55,TPMNEW,INTOWP(ICNT),IBPT,ILPT)           
              IBPT = ILPT                                               
C      IF(IR.LE.3) THEN                                                 
C               WRITE(6,*) ' OUT TO IT55'                               
C               WRITE(6,*) (TPMNEW(I),I=1,ICNT)                         
C      ENDIF                                                            
              ICNT = 0                                                  
              ISIG = ISIG + NSIG                                        
   50       CONTINUE                                                    
            IRHO = IRHO + NRHO                                          
   60     CONTINUE                                                      
C**************************                                             
C     WRITE TO FILE55     *                                             
C**************************                                             
          INU = INU + NNU                                               
   70   CONTINUE                                                        
        IMU = IMU + NMU                                                 
   80 CONTINUE                                                          
      WRITE(6,*) 'THERE ARE',(ILPT-AOSTRT)/2,' TPDM ELEMENTS WRITTEN '  
      RETURN                                                            
      END                                                               
