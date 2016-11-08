      SUBROUTINE BACK2(TPMOLD,TPMNEW,MAPIJ,SQ,XTP,EIV,IPMOLD,IPMNEW,    
     .                 AOPS)                                            
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      INTEGER AOSTRT,IPMOLD(1),IPMNEW(1),BKTSZ,OFSET,AOPS(NSHELL)       
      COMMON/TAPES/ IOUT,ICHECK,IT30,IT55,IT91,IT95,IT96                
      COMMON /SRT/  NBKT,NGRPIJ,LGRPIJ,NBKT2,LGRP2,IBPT,IAPT,IGRP       
      COMMON /DIMS/ NIJ,NBF,NIJAO,NBFAO,LENR,NST,NSHELL,AOSTRT,BKTSZ    
      DIMENSION TPMOLD(LENR),TPMNEW(LENR),MAPIJ(NIJAO),                 
     .          SQ(NBFAO*NBFAO),XTP(NBFAO*NBFAO),EIV(NBF,NBFAO)         
C                                                                       
C************************************                                   
C     NOTE THAT EIV IS TRANSPOSE!   *                                   
C************************************                                   
      NBFTMP = 0                                                        
      NSHL   = 1                                                        
      NBFTMP = NBFTMP + AOPS(NSHL)                                      
      NTRI   = NBFTMP*(NBFTMP + 1)/2                                    
      IBPT = 1                                                          
      WRITE(ICHECK,10)                                                  
   10 FORMAT(/5X,'IN BACK2')                                            
C*****************************************                              
C     LOOP OVER NEW NUMBER OF BUCKETS    *                              
C*****************************************                              
      DO 90 IBKT2 = 1,NBKT2                                             
        IF(IBKT2.EQ.NBKT2) THEN                                         
          L2 = LGRP2                                                    
        ELSE                                                            
          L2 = NGRPIJ                                                   
        ENDIF                                                           
C*****************************************                              
C     LOOP OVER OLD NUMBER OF BUCKETS    *                              
C*****************************************                              
        DO 40 IBKT = 1,NBKT                                             
          IF(IBKT.EQ.NBKT) THEN                                         
            L1 = LGRPIJ                                                 
          ELSE                                                          
            L1 = NGRPIJ                                                 
          ENDIF                                                         
C*******************************************                            
C     READ IN ANOTHER CHUNK OF FILE95.     *                            
C*******************************************                            
          IPT95 = (IBKT2-1)*NGRPIJ*L1 + (IBKT-1)*NGRPIJ*NIJAO           
          LENGTH = L1*L2                                                
C     WRITE(6,*) 'IT95,IPT95,LENGTH',IT95,IPT95,LENGTH                  
          CALL WREADW(IT95,TPMOLD,INTOWP(LENGTH),INTOWP(IPT95)+1,JUNK)  
C     WRITE(6,*) 'TPMOLD IN IT95'                                       
C     WRITE(6,*) (TPMOLD(I),I=1,NIJAO)                                  
C********************************************************************   
C     SORT FILE 95 INTO ARRAYS OF DIMENSION NIJ, SPACED BY NIJAO    *   
C********************************************************************   
          DO 30 I = 1,L2                                                
            DO 20 J = 1,L1                                              
              TPMNEW( (IBKT-1)*NGRPIJ + (I-1)*NIJAO + J ) =             
     .                             TPMOLD( (I-1)*L1 + J )               
   20       CONTINUE                                                    
   30     CONTINUE                                                      
   40   CONTINUE                                                        
C       WRITE(6,*) 'TPMNEW'                                             
C       WRITE(6,*) TPMNEW                                               
C****************************************************                   
C     TRANSFORM TPMNEW, PUT RESULT INTO TPMNEW.     *                   
C****************************************************                   
        ICNT = 0                                                        
        DO 80 I = 1,L2                                                  
          NIJT = (IBKT2-1)*NGRPIJ + I                                   
          IF(NIJT.GT.NTRI) THEN                                         
            NSHL   = NSHL + 1                                           
            NBFTMP = NBFTMP + AOPS(NSHL)                                
            NTRI   = NBFTMP*(NBFTMP + 1)/2                              
C           WRITE(6,*) 'NSHL',NSHL                                      
C           WRITE(6,*) 'NTRI',NTRI                                      
C           WRITE(6,*) 'NBFTMP',NBFTMP                                  
          ENDIF                                                         
          OFSET = (I-1)*NIJAO + 1                                       
          CALL SQUR(TPMNEW(OFSET),SQ,NBF)                               
C       DO 11 J=1,NBF*NBF                                               
C         IF(SQ(J).GT.1.0D+10) THEN                                     
C           WRITE(6,*) 'TPMNEW TOO LARGE! AT J=',J,'IBKT=',IBKT,'IBKT2='
C    ,      ,IBKT                                                       
C           WRITE(6,*) (SQ(K), K=1,NBF*NBF)                             
C           STOP                                                        
C         ENDIF                                                         
C  11   CONTINUE                                                        
          CALL CTSQC(EIV,SQ,XTP,SQ,NBF,NBFAO,NBFTMP)               
C         CALL MXMB(EIV,1,NBFAO,SQ,1,NBF,XTP,1,NBFAO,NBFAO,NBF,NBF)     
C         IF(I.EQ.1) THEN                                               
C           WRITE(6,*) 'EIV'                                            
C           WRITE(6,*) EIV                                              
C           WRITE(6,*) 'SQ'                                             
C           WRITE(6,*) SQ                                               
C           WRITE(6,*) 'XTP'                                            
C           WRITE(6,*) XTP                                              
C         ENDIF                                                         
          MUNU = 0                                                      
          DO 70 MU = 1,NBFTMP                                           
C           CALL ZERO(SQ,MU)                                            
C           CALL MXMB(XTP,1,NBFAO,EIV(MU,1),NBFAO,1,SQ,1,NBFAO,         
C    .                MU,NBFAO,1)                                       
C         IF(I.EQ.1) THEN                                               
C           WRITE(6,*) 'SQ'                                             
C           WRITE(6,*) SQ                                               
C         ENDIF                                                         
            DO 60 NU = 1,MU                                             
              MUNU = MUNU + 1                                           
              NWMUNU = MAPIJ(MUNU)                                      
              TPMNEW(OFSET + NWMUNU - 1) = SQ(MUNU)                     
   60       CONTINUE                                                    
   70     CONTINUE                                                      
C       WRITE(6,*) 'INTOWP(OFSET),INTOWP(NTRI),IBPT,ILPT'               
C       WRITE(6,*) INTOWP(OFSET),INTOWP(NTRI),IBPT,ILPT                 
C       WRITE(6,*) 'TRANSFORMED TWO PDM'                                
C       WRITE(6,*) (TPMNEW(J),J=OFSET,OFSET + NTRI-1)                   
        CALL WWRITW(IT96,TPMNEW(OFSET),INTOWP(NTRI),IBPT,ILPT)          
        IBPT = ILPT                                                     
   80   CONTINUE                                                        
C***********************************                                    
C     WRITE TPMNEW INTO FILE96     *                                    
C***********************************                                    
C       WRITE(6,*) 'TRANSFORMED TPM WRITTEN TO IT96'                    
C       WRITE(6,*) 'LENGTH,IBPT,ILPT',NIJAO*L2,IBPT,ILPT                
C       WRITE(6,*) TPMNEW                                               
   90 CONTINUE                                                          
      RETURN                                                            
      END                                                               
