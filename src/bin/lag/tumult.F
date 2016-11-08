      SUBROUTINE tumult(ISYM,F,tpdm,ints,G,W,LENR)                                
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      INTEGER ISYM(NORBS)                                               
      REAL*8 F(NBF,NBF),G(LENR),W(LENR),tpdm(lenr),ints(lenr)           
C                                                                       
      COMMON/TAPES/ IOUT,ICHECK,IT30,IT52,IT53,IT58,IT91,IT93           
      COMMON /TUM/ IJSYM,IJSRT,ISUP,JSUP                                
      COMMON /DIMS/ ISYMRB,NUMIJ,NIJ,NBF,NIJAO,NBFAO,NORBS,NMAX,NSYMHF  
      COMMON /PRNT/ IPRINT                                              
C****************************************************                   
C     ARRAY G CONTAINS DENSITY MATRIX ELEMENTS      *                   
C     ARRAY W CONTAINS INTEGRALS                    *                   
C****************************************************                   
       IWPT=IJSRT                                                       
       IGPT=IJSRT                                                       
       IPT=IGPT                                                         
       IJ=0                                                             
       DO 30 I=1,NBF                                                    
         JPT=IWPT                                                       
         DO 20 J=1,NBF                                                  
           IJ=IJ+1                                                      
           IF(ISYM(I).NE.ISYM(J)) GO TO 20                              
           IK=IPT                                                       
           JK=JPT                                                       
           KMX=J-1                                                      
           IF(I.LT.J) KMX = I - 1                                       
           DO 10 K=1,KMX                                                
             IK=IK+1                                                    
             JK=JK+1                                                    
             F(I,J)=F(I,J)+G(IK)*W(JK)                                  
  10       CONTINUE                                                     
  20     JPT=JPT+J                                                      
         IPT=IPT+I                                                      
  30   CONTINUE                                                         
       IK=IGPT                                                          
       KPT=IWPT                                                         
       DO 60 K=1,NBF                                                    
         IPT=0                                                          
         LSM=IEOR((ISYM(K)-1),IJSYM)                                    
         DO 50 I=1,K                                                    
           IK=IK+1                                                      
           IF((ISYM(I)-1).NE.LSM) GO TO 50                              
           JK=KPT                                                       
           IJ=IPT                                                       
           DO 40 J=1,K                                                  
             JK=JK+1                                                    
             IJ=IJ+1                                                    
             F(I,J)=F(I,J)+G(IK)*W(JK)                                  
  40       CONTINUE                                                     
  50     IPT=IPT+NBF                                                    
         KPT=KPT+K                                                      
  60   CONTINUE                                                         
       IPT=0                                                            
       IK=IGPT                                                          
       DO 90 I=1,NBF                                                    
         LSM=IEOR((ISYM(I)-1),IJSYM)                                    
         KMX=I-1                                                        
         KPT=IWPT                                                       
         DO 80 K=1,KMX                                                  
           IK=IK+1                                                      
           IF((ISYM(K)-1).NE.LSM) GO TO 80                              
           JK=KPT                                                       
           IJ=IPT                                                       
           DO 70 J=1,K                                                  
             JK=JK+1                                                    
             IJ=IJ+1                                                    
             F(I,J)=F(I,J)+G(IK)*W(JK)                                  
  70       CONTINUE                                                     
  80     KPT=KPT+K                                                      
         IPT=IPT+NBF                                                    
         IK=IK+1                                                        
  90   CONTINUE                                                         
       JK=IWPT                                                          
       DO 120 J=1,NBF                                                   
         LSM=IEOR((ISYM(J)-1),IJSYM)                                    
         KMX=J-1                                                        
         KPT=IGPT                                                       
         DO 110 K=1,KMX                                                 
           JK=JK+1                                                      
           IF((ISYM(K)-1).NE.LSM) GO TO 110                             
           IK=KPT                                                       
           IJ=J                                                         
           DO 100 I=1,K                                                 
             IK=IK+1                                                    
             F(I,J)=F(I,J)+G(IK)*W(JK)                                  
             IJ=IJ+NBF                                                  
 100       CONTINUE                                                     
 110     KPT=KPT+K                                                      
 120   JK=JK+1                                                          
C                                                                       
       RETURN                                                           
       END                                                              
