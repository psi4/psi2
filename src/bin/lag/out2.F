      SUBROUTINE OUT2(IOFF,F,X,EIV,ORBTBF,BUF,XAO,ia47)                 
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      INTEGER IOFF(NIJ),ORBTBF(NBF),ia47(192)                           
      REAL*8 X(NBF,NBF),F(NBF,NBF),XAO(NIJAO),                          
     1       EIV(NBFAO,NBF),BUF(NBF,NBFAO)                              
C                                                                       
      COMMON /DIMS/ ISYMRB,NUMIJ,NIJ,NBF,NIJAO,NBFAO,NORBS,NMAX,NSYMHF  
      COMMON/TAPES/ IOUT,ICHECK,IT30,IT52,IT53,IT58,IT91,IT93           
      COMMON /PRNT/ IPRINT                                              
      it47 = 47
C*******************************************************************    
C     READ ONEPDM FROM THE MO BASIS TO MASTER FILE 40 (TRIANG).    *    
C*******************************************************************    
      CALL MREAD(XAO,35)                                                
      IJ = 0                                                            
      DO 20 IOR=1,NORBS                                                 
        DO 10 JOR=1,IOR                                                 
          IJ = IJ + 1                                                   
          F(IOR,JOR) = XAO(IJ)                                          
          F(JOR,IOR) = XAO(IJ)                                          
   10   CONTINUE                                                        
   20 CONTINUE                                                          
      IF(IPRINT.GE.2) THEN                                              
        WRITE(*,30)                                                     
   30   FORMAT(/1X,'ONE PDM IN MO BASIS AS READ FROM 40'/)              
        CALL PRINT(XAO,NIJ,NBF,IOUT)                                    
      END IF                                                            
C*************************************************************          
C     WRITE THE LAGRANGIAN IN DRT ORDERING TO TAPE40.        *          
C     AND ALSO THE PITZER TO DRT ORDERING ARRAY, ORBTBF.     *          
C*************************************************************          
      CALL MWRIT(ORBTBF,5)                                              
      CALL MWRIT(X,33)                                                  
ctph
ctph  lagpt = ia47(119) + ia47(169)
ctph  laglen = intowp(nbf*nbf)
      lagpt = ia47(120)              
      laglen = ia47(170)              
ctph  call wwritw(it47,lagpt,1,120,junk)
ctph  call wwritw(it47,laglen,1,170,junk)
      call wwritw(it47,x,laglen,lagpt,junk)
ctph
      IF(IPRINT.GE.2) THEN                                              
        WRITE(IOUT,40)                                                  
   40   FORMAT(/5X,'REORDERED LAGRANGIAN AS WRITTEN TO 40'/)            
        CALL MATOUT(X,NBF,NBF,NBF,NBF,IOUT)                             
      END IF                                                            
C**************************************************************         
C     NOW PUT LAGRANGIAN INTO AO BASIS AND WRITE TO MASTER    *         
C**************************************************************         
      DO 70 II=1,NBFAO                                                  
        DO 60 I=1,NBF                                                   
          VAL=0.0D0                                                     
          DO 50 J=1,NBF                                                 
            VAL=VAL+EIV(II,J)*X(I,J)                                    
   50     CONTINUE                                                      
          BUF(I,II)=VAL                                                 
   60   CONTINUE                                                        
   70 CONTINUE                                                          
      IIJJ=0                                                            
      DO 100 II=1,NBFAO                                                 
        DO 90 JJ=1,II                                                   
          IIJJ=IIJJ+1                                                   
          VAL=0.0D0                                                     
          DO 80 I=1,NBF                                                 
            VAL=VAL+EIV(II,I)*BUF(I,JJ)                                 
   80     CONTINUE                                                      
          XAO(IIJJ)=VAL                                                 
   90   CONTINUE                                                        
  100 CONTINUE                                                          
      IF (IPRINT.GE.2) THEN                                             
        WRITE(IOUT,110)                                                 
  110   FORMAT(/5X,'AO LAGRANGIAN'/)                                    
        CALL PRINT(XAO,NIJAO,NBFAO,IOUT)                                
      END IF                                                            
      CALL MWRIT(XAO,32)                                                
ctph
ctph  iptao = lagpt + intowp(nbf*nbf)
ctph  lenao = intowp(nijao)
      iptao = ia47(121)              
      lenao = ia47(171)              
ctph  call wwritw(it47,iptao,1,121,junk)
ctph  call wwritw(it47,lenao,1,171,junk)
      call wwritw(it47,xao,lenao,iptao,junk)
ctph
C*********************************************************              
C     NOW PUT 1PDM INTO AO BASIS AND WRITE TO MASTER     *              
C*********************************************************              
      DO 140 II=1,NBFAO                                                 
        DO 130 I=1,NBF                                                  
          VAL=0.0D0                                                     
          DO 120 J=1,NBF                                                
            VAL=VAL+EIV(II,J)*F(I,J)                                    
  120     CONTINUE                                                      
          BUF(I,II)=VAL                                                 
  130   CONTINUE                                                        
  140 CONTINUE                                                          
      IIJJ=0                                                            
      DO 170 II=1,NBFAO                                                 
        DO 160 JJ=1,II                                                  
          IIJJ=IIJJ+1                                                   
          VAL=0.0D0                                                     
          DO 150 I=1,NBF                                                
            VAL=VAL+EIV(II,I)*BUF(I,JJ)                                 
  150     CONTINUE                                                      
          XAO(IIJJ)=VAL                                                 
  160   CONTINUE                                                        
  170 CONTINUE                                                          
      IF (IPRINT.GE.2) THEN                                             
        WRITE(IOUT,180)                                                 
  180   FORMAT(/1X,'AO ONE PDM'/)                                       
        CALL PRINT(XAO,NIJAO,NBFAO,IOUT)                                
      END IF                                                            
      CALL MWRIT(XAO,34)                                                
ctph
ctph  ipqao = lagpt + intowp(nbf*nbf)
ctph  lenao = intowp(nijao)
      ipqao = ia47(122)              
      lenao = ia47(172)              
ctph  call wwritw(it47,ipqao,1,122,junk)
ctph  call wwritw(it47,lenao,1,172,junk)
      call wwritw(it47,xao,lenao,ipqao,junk)
ctph
C**************                                                         
C     DONE    *                                                         
C**************                                                         
      RETURN                                                            
      END                                                               
