      SUBROUTINE GETIT(IOFF,IJADD,IJGRP,KADD,LADD,                      
     1 X,F,H,BUF,ORBSYM,START,ILAST,ILAST2,INTOFF,BKTPNT,               
     2 SRTAR,BKTSZ2,ISRTAR,DUM,ijoff,ijbkt,idfile)                      
C                                                                       
      IMPLICIT INTEGER (A-Z)                                            
      REAL*8 BUF(NMAX),VAL,F(NBF,NBF),SRTAR(IBIG2),X(NBF,NBF)           
      REAL*8 H(NBF,NBF),D8,D4,D2,VAL2                                   
      REAL*8 DUM(56),REP,ECI1                                           
      LOGICAL IJEQ,IJKEQ,JKEQ,CCSD                                      
C                                                                       
      COMMON /BIGD/ IBIG,IBIG2                                          
      COMMON /FILES/ IT40,IT54,CCSD                                     
      COMMON /OFF/ MAXINT                                               
      COMMON /ORD/ IFAC1,IFAC2                                          
      COMMON/TAPES/ IOUT,ICHECK,IT30,IT52,IT53,IT58,IT91,IT93           
      COMMON /SRT/ VAL,IADR,IBKT,NBKT,NGRPIJ,LENR                       
      COMMON /DIMS/ ISYMRB,NUMIJ,NIJ,NBF,NIJAO,NBFAO,NORBS,NMAX,NSYMHF  
      COMMON /NREP/ REP,ECI1                                            
      COMMON /PRNT/ IPRINT                                              
C                                                                       
      DIMENSION IOFF(NIJ),IJGRP(NUMIJ),ORBSYM(NORBS),ISRTAR(IBIG),      
     *          IJADD(NUMIJ),KADD(ISYMRB),LADD(ISYMRB)                  
      DIMENSION START(NBKT),ILAST(NBKT),ILAST2(NBKT),INTOFF(NBKT),      
     *          BKTPNT(NBKT),ijoff(nij),ijbkt(nij)                      
C                                                                       
      WRITE(ICHECK,10)                                                  
   10 FORMAT(/5X,'IN GETIT')                                            
      ECI1=0.0D0                                                        
C************************************                                   
C     DECIDE PERMUTATION FACTORS    *                                   
C************************************                                   
      D8=0.125D0                                                        
      D4=0.25D0                                                         
      D2=0.5D0                                                          
      IF(CCSD) THEN                                                     
        D8=1.0D0                                                        
        D4=1.0D0                                                        
        D2=1.0D0                                                        
      END IF                                                            
C***********************                                                
C     ZERO OUT F,H     *                                                
C***********************                                                
      DO 30 II=1,NBF                                                    
        DO 20 JJ=1,NBF                                                  
          F(II,JJ)=0.0D0                                                
          H(II,JJ)=0.0D0                                                
   20   CONTINUE                                                        
   30 CONTINUE                                                          
C********************************************************************** 
C    READ IN THE ONE AND TWO PARTICLE DENSITY MATRICES IN DRT ORDER   * 
C********************************************************************** 
      CALL SREW(IT53)                                                   
      CALL SRTONE(BKTPNT,START,ILAST,BKTSZ2,SRTAR,ISRTAR,INTOFF,IT91)   
      SRIN=0                                                            
      IGRP=0                                                            
      DO 170 IOR=NORBS,1,-1                                             
        ISM=ORBSYM(IOR)                                                 
        ITYP=1                                                          
        DO 160 JOR=IOR,1,-1                                             
          IF (IGRP .NE. IJGRP(IOFF(IOR)+JOR)) THEN                      
            CALL SREAD(IT53,BUF,INTOWP(NMAX))                           
            SRIN=SRIN+NMAX*2/112                                        
            IGRP=IJGRP(IOFF(IOR)+JOR)                                   
          END IF                                                        
          JSM=ORBSYM(JOR)                                               
          IJEQ=IOR.EQ.JOR                                               
          IJT=ITYP                                                      
          IF (IJEQ) IJT=IJT+2                                           
          IJSYM=IEOR((ISM-1),(JSM-1))                                   
          KOFF=IJSYM*NORBS                                              
          JAD=IJADD(IOFF(IOR)+JOR)                                      
          DO 150 KOR=1,JOR                                              
            KSM=ORBSYM(KOR)                                             
            JKEQ=JOR.EQ.KOR                                             
            IJKT=IJT                                                    
            IF (JKEQ) IJKT=IJKT+1                                       
            IJKEQ=IJEQ.AND.JKEQ                                         
            IJKSYM=IEOR(IJSYM,(KSM-1))                                  
            LOFF=IJKSYM*NORBS                                           
            KAD=JAD+KADD(KOFF+KOR)                                      
            DO 140 LOR=1,KOR                                            
              IF(JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GO TO 130
              LSM=ORBSYM(LOR)                                           
              IF ((LSM-1) .NE. IJKSYM) GO TO 130                        
              IJKLT=IJKT                                                
              IF (IJKEQ .AND. KOR.GT.LOR) IJKLT=IJKLT+1                 
              IF (KOR.EQ.LOR) IJKLT=IJKLT+3                             
              LAD=KAD+LADD(LOFF+LOR)                                    
C************************                                               
C    CHECK FOR ERROR    *                                               
C************************                                               
              IF (IJKLT.GT.7 .OR. IJKLT.LT.1) THEN                      
                WRITE(IOUT,40) IJKLT,IOR,JOR,KOR,LOR,IBF,JBF,KBF,LBF,   
     *                        ISM,JSM,KSM,LSM                           
   40           FORMAT(1X,'  IN GETDM, IJKLT UNKNOWN TYPE OF',I5/       
     *                    '  ORBITALS :',4I5/                           
     *                    '  BASIS FUNCTIONS :',4I5/                    
     *                    '   SYMMETRIES :',4I5)                        
                STOP                                                    
              END IF                                                    
              IF (LAD .GT. NMAX) THEN                                   
                WRITE(IOUT,50) LAD                                      
   50           FORMAT(1X,'  IN GETDM, LAD TOO BIG :',I6)               
                STOP                                                    
              END IF                                                    
              GO TO (60,70,80,90,100,110,120) IJKLT                     
   60         VAL=BUF(LAD+1)*D8                                         
              CALL GORDR(IOR,KOR,JOR,LOR,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast,bktsz2,srtar,isrtar,intoff,it91)
              VAL=BUF(LAD+2)*D8                                         
              CALL GORDR(IOR,JOR,KOR,LOR,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast,bktsz2,srtar,isrtar,intoff,it91)
              VAL=BUF(LAD+3)*D8                                         
              CALL GORDR(IOR,LOR,JOR,KOR,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast,bktsz2,srtar,isrtar,intoff,it91)
              GO TO 130                                                 
   70         VAL=BUF(LAD+1)*D8                                         
              CALL GORDR(IOR,JOR,JOR,LOR,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast,bktsz2,srtar,isrtar,intoff,it91)
              VAL=BUF(LAD+2)*D4                                         
              CALL GORDR(IOR,LOR,JOR,JOR,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast,bktsz2,srtar,isrtar,intoff,it91)
              GO TO 130                                                 
   80         VAL=BUF(LAD+1)*D8                                         
              CALL GORDR(IOR,KOR,IOR,LOR,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast,bktsz2,srtar,isrtar,intoff,it91)
              VAL=BUF(LAD+2)*D4                                         
              CALL GORDR(IOR,IOR,KOR,LOR,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast,bktsz2,srtar,isrtar,intoff,it91)
              GO TO 130                                                 
   90         VAL=BUF(LAD+1)*D8                                         
              CALL GORDR(IOR,LOR,JOR,LOR,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast,bktsz2,srtar,isrtar,intoff,it91)
              VAL=BUF(LAD+2)*D4                                         
              CALL GORDR(IOR,JOR,LOR,LOR,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast,bktsz2,srtar,isrtar,intoff,it91)
              GO TO 130                                                 
  100         VAL=BUF(LAD+1)*D4                                         
              CALL GORDR(IOR,IOR,IOR,LOR,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast,bktsz2,srtar,isrtar,intoff,it91)
              VAL=BUF(LAD+2)*D4                                         
              CALL GORDR(IOR,LOR,LOR,LOR,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast,bktsz2,srtar,isrtar,intoff,it91)
              VAL=BUF(LAD+3)*D2                                         
              F(IOR,LOR)=VAL                                            
              F(LOR,IOR)=VAL                                            
              GO TO 130                                                 
  110         VAL=BUF(LAD+1)*D4                                         
              CALL GORDR(IOR,LOR,IOR,LOR,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast,bktsz2,srtar,isrtar,intoff,it91)
              VAL=BUF(LAD+2)*D2                                         
              CALL GORDR(IOR,IOR,LOR,LOR,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast,bktsz2,srtar,isrtar,intoff,it91)
              GO TO 130                                                 
  120         VAL=BUF(LAD+1)                                            
              CALL GORDR(IOR,IOR,IOR,IOR,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast,bktsz2,srtar,isrtar,intoff,it91)
              VAL=BUF(LAD+2)                                            
              F(IOR,IOR)=VAL                                            
  130         CONTINUE                                                  
  140       CONTINUE                                                    
  150     CONTINUE                                                      
  160   CONTINUE                                                        
  170 CONTINUE                                                          
      CALL SRTLST(BKTPNT,START,ILAST,BKTSZ2,SRTAR,ISRTAR,INTOFF,IT91)   
C***************************                                            
C     GET RID OF FILE53    *                                            
C***************************                                            
ctph
      if(idfile.eq.1) then
        CALL RCLOSE(IT53,4)                                               
      else 
        CALL RCLOSE(IT53,3)                                               
      endif
      WRITE(ICHECK,180)                                                 
      WRITE(6  ,180)                                                    
  180 FORMAT(' GETTING RID OF FILE 53')                                 
C********************************************************************** 
C   POSITION TAPE2 AT THE BEGINNING OF THE MO INTEGRALS, AND PICK UP  * 
C   THE NUCLEAR REPULSION ENERGY IN POSITION 31 OF REAL ARRAY DUM.    * 
C********************************************************************** 
      IF (.NOT.CCSD) CALL SREAD(IT52,DUM,112)                           
      REP=DUM(31)                                                       
C******************************************************************     
C    READ IN THE ONE AND TWO ELECTRON MO INTEGRALS IN DRT ORDER   *     
C******************************************************************     
      WRITE(ICHECK,190)                                                 
  190 FORMAT(1X,'  READING MO INTEGRALS ')                              
C                                                                       
      CALL SRTONE(BKTPNT,START,ILAST2,BKTSZ2,SRTAR,ISRTAR,INTOFF,IT93)  
      SRIN=0                                                            
      IGRP2=0                                                           
      DO 310 IOR=NORBS,1,-1                                             
        ISM=ORBSYM(IOR)                                                 
        ITYP=1                                                          
        DO 300 JOR=IOR,1,-1                                             
          IJ=IOFF(IOR)+JOR                                              
          JJ=IOFF(JOR)+JOR                                              
          IF (IGRP2 .NE. IJGRP(IOFF(IOR)+JOR)) THEN                     
            CALL SREAD(IT52,BUF,INTOWP(NMAX))                           
            SRIN=SRIN+NMAX*2/112                                        
            IGRP2=IJGRP(IOFF(IOR)+JOR)                                  
          END IF                                                        
          JSM=ORBSYM(JOR)                                               
          IJEQ=IOR.EQ.JOR                                               
          IJT=ITYP                                                      
          IF (IJEQ) IJT=IJT+2                                           
          IJSYM=IEOR((ISM-1),(JSM-1))                                   
          KOFF=IJSYM*NORBS                                              
          JAD=IJADD(IOFF(IOR)+JOR)                                      
          DO 290 KOR=1,JOR                                              
            IK=IOFF(IOR)+KOR                                            
            JK=IOFF(JOR)+KOR                                            
            KK=IOFF(KOR)+KOR                                            
            KSM=ORBSYM(KOR)                                             
            JKEQ=JOR.EQ.KOR                                             
            IJKT=IJT                                                    
            IF (JKEQ) IJKT=IJKT+1                                       
            IJKEQ=IJEQ.AND.JKEQ                                         
            IJKSYM=IEOR(IJSYM,(KSM-1))                                  
            LOFF=IJKSYM*NORBS                                           
            KAD=JAD+KADD(KOFF+KOR)                                      
            DO 280 LOR=1,KOR                                            
              IL=IOFF(IOR)+LOR                                          
              JL=IOFF(JOR)+LOR                                          
              KL=IOFF(KOR)+LOR                                          
              LL=IOFF(LOR)+LOR                                          
              IF(JOR.EQ.KOR .AND. KOR.EQ.LOR .AND. LOR.LT.IOR) GO TO 270
              LSM=ORBSYM(LOR)                                           
              IF ((LSM-1) .NE. IJKSYM) GO TO 270                        
              IJKLT=IJKT                                                
              IF (IJKEQ .AND. KOR.GT.LOR) IJKLT=IJKLT+1                 
              IF (KOR.EQ.LOR) IJKLT=IJKLT+3                             
              LAD=KAD+LADD(LOFF+LOR)                                    
C************************                                               
C    CHECK FOR ERROR    *                                               
C************************                                               
              IF (IJKLT.GT.7 .OR. IJKLT.LT.1) THEN                      
                WRITE(IOUT,40) IJKLT,IOR,JOR,KOR,LOR,IBF,JBF,KBF,LBF,   
     *                        ISM,JSM,KSM,LSM                           
                STOP                                                    
              END IF                                                    
              IF (LAD .GT. NMAX) THEN                                   
                WRITE(6,50)LAD                                          
                STOP                                                    
              END IF                                                    
              GO TO (200,210,220,230,240,250,260) IJKLT                 
  200         VAL2=BUF(LAD+1)                                           
              IFAC1=8                                                   
              IFAC2=8                                                   
              CALL LIORDR(IOR,KOR,JOR,LOR,VAL2,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast2,bktsz2,srtar,isrtar,intoff,it93)
              VAL2=BUF(LAD+2)                                           
              CALL LIORDR(IOR,JOR,KOR,LOR,VAL2,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast2,bktsz2,srtar,isrtar,intoff,it93)
              VAL2=BUF(LAD+3)                                           
              CALL LIORDR(IOR,LOR,JOR,KOR,VAL2,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast2,bktsz2,srtar,isrtar,intoff,it93)
              GO TO 270                                                 
  210         VAL2=BUF(LAD+1)                                           
              IFAC1=8                                                   
              IFAC2=8                                                   
              CALL LIORDR(IOR,JOR,JOR,LOR,VAL2,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast2,bktsz2,srtar,isrtar,intoff,it93)
              VAL2=BUF(LAD+2)                                           
              IFAC2=4                                                   
              CALL LIORDR(IOR,LOR,JOR,JOR,VAL2,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast2,bktsz2,srtar,isrtar,intoff,it93)
              GO TO 270                                                 
  220         VAL2=BUF(LAD+1)                                           
              IFAC1=8                                                   
              IFAC2=8                                                   
              CALL LIORDR(IOR,KOR,IOR,LOR,VAL2,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast2,bktsz2,srtar,isrtar,intoff,it93)
              VAL2=BUF(LAD+2)                                           
              IFAC1=4                                                   
              CALL LIORDR(IOR,IOR,KOR,LOR,VAL2,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast2,bktsz2,srtar,isrtar,intoff,it93)
              GO TO 270                                                 
  230         VAL2=BUF(LAD+1)                                           
              IFAC1=8                                                   
              IFAC2=8                                                   
              CALL LIORDR(IOR,LOR,JOR,LOR,VAL2,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast2,bktsz2,srtar,isrtar,intoff,it93)
              VAL2=BUF(LAD+2)                                           
              IFAC2=4                                                   
              CALL LIORDR(IOR,JOR,LOR,LOR,VAL2,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast2,bktsz2,srtar,isrtar,intoff,it93)
              GO TO 270                                                 
  240         VAL2=BUF(LAD+1)                                           
              IFAC1=4                                                   
              IFAC2=8                                                   
              CALL LIORDR(IOR,IOR,IOR,LOR,VAL2,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast2,bktsz2,srtar,isrtar,intoff,it93)
              VAL2=BUF(LAD+2)                                           
              IFAC1=8                                                   
              IFAC2=4                                                   
              CALL LIORDR(IOR,LOR,LOR,LOR,VAL2,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast2,bktsz2,srtar,isrtar,intoff,it93)
              VAL=BUF(LAD+3)                                            
              H(IOR,LOR)=VAL                                            
              H(LOR,IOR)=VAL                                            
              GO TO 270                                                 
  250         VAL2=BUF(LAD+1)                                           
              IFAC1=8                                                   
              IFAC2=8                                                   
              CALL LIORDR(IOR,LOR,IOR,LOR,VAL2,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast2,bktsz2,srtar,isrtar,intoff,it93)
              VAL2=BUF(LAD+2)                                           
              IFAC1=4                                                   
              IFAC2=4                                                   
              CALL LIORDR(IOR,IOR,LOR,LOR,VAL2,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast2,bktsz2,srtar,isrtar,intoff,it93)
              GO TO 270                                                 
  260         VAL2=BUF(LAD+1)                                           
              IFAC1=4                                                   
              IFAC2=4                                                   
              CALL LIORDR(IOR,IOR,IOR,IOR,VAL2,ioff,ijoff,ijbkt,
     *          bktpnt,start,ilast2,bktsz2,srtar,isrtar,intoff,it93)
              VAL=BUF(LAD+2)                                            
              H(IOR,IOR)=VAL                                            
  270         CONTINUE                                                  
  280       CONTINUE                                                    
  290     CONTINUE                                                      
  300   CONTINUE                                                        
  310 CONTINUE                                                          
      CALL SRTLST(BKTPNT,START,ILAST2,BKTSZ2,SRTAR,ISRTAR,INTOFF,IT93)   
C********************************************************************   
C     FINISHED SORT OF TPDM AND MO INTEGRALS INTO IT91 AND IT93     *   
C     AND PUTTING ONE PDM INTO F, ONE ELECTRON MO INTEGRALS IN H    *   
C********************************************************************   
C*****************************************************                  
C     ONE ELECTRON CONTRIBUTION TO THE LAGRANGIAN    *                  
C*****************************************************                  
      CALL ZERO(X,NBF*NBF)                                              
      DO 340 I=1,NORBS                                                  
        DO 330 M=1,NORBS                                                
          DO 320 J=1,NORBS                                              
            X(I,M) = X(I,M) + H(M,J)*F(I,J)*2.0D0                       
  320     CONTINUE                                                      
  330   CONTINUE                                                        
  340 CONTINUE                                                          
      DO 350 I=1,NORBS                                                  
        ECI1 = ECI1 + X(I,I)*0.5D0                                      
  350 CONTINUE                                                          
      IF(IPRINT.GE.2) THEN                                              
        WRITE(6,*) 'ONE ELECTRON LAGRANGIAN AT END OF GETIT'            
        CALL MATOUT(X,NBF,NBF,NBF,NBF,IOUT)                             
        WRITE(6,*) 'H AT END OF GETIT'                                  
        CALL MATOUT(H,NBF,NBF,NBF,NBF,IOUT)                             
        WRITE(6,*) 'F AT END OF GETIT'                                  
        CALL MATOUT(F,NBF,NBF,NBF,NBF,IOUT)                             
      ENDIF                                                             
      RETURN                                                            
      END                                                               
