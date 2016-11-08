      SUBROUTINE FAMAT(EPA,TA,HA,BUF,IOFF,JT1C,IOUT,DROL,LSTEPA)        
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      INTEGER IOFF(NUMIJ),JT1C(1),IOUT(NORBS),DROL(NORBS)               
      DIMENSION EPA(NORBS,NORBS),TA(NUMIJ),HA(NUMIJ),BUF(NUMIJ,NBSET)   
      REAL*8    LSTEPA(NUMIJ)                                           
      COMMON /TAPES/ NOUT,ICHK,IT30,IT42,IT47,IT52,IT58,IT94,IT96,IT98  
      COMMON /DIMSG/ ISYMRB,NUMIJ,NBF,NMAX,NBFORB,NSYM,NDF,NORBS        
      COMMON /POINT/ FOCC(5),IPOINT(5),JPOINT(5),NSORB(5)               
      COMMON /FUNCS/ IOCC,JOCC,KOCC,NTYPES,NBSET,NAT                    
      COMMON /IEXT/ IPRINT,IGRP,ICIDIP                                  
      COMMON /CALCJ/ ICALT,IDERT,ICIT                                   
C                                                                       
C     WRITE(6   ,2222)                                                  
      WRITE(ICHK,2222)                                                  
 2222 FORMAT(/'  IN FAMAT')                                             
      NORBSQ=NORBS*NORBS                                                
      IHA=JT1C(NDF+1)                                                   
      IEA=2*NDF+NBSET*NDF                                               
      ISEA=JT1C(IEA)+INTOWP(NUMIJ)                                      
      DO 1 IDF=1,NDF                                                    
        IEA=IEA+1                                                       
        DO 10 I=1,NORBS                                                 
          DO 10 J=1,NORBS                                               
 10     EPA(I,J)=0.0D+00                                                
C                                                                       
        CALL WREADW(IT94,HA,INTOWP(NUMIJ),IHA,IHA)                      
        IF (IAND(IPRINT,32) .NE. 0) THEN                                
          WRITE(NOUT,1101) IDF,IHA                                      
 1101     FORMAT(/'  HA FOR DEGREE OF FREEDOM ',I5,'   IHA=',I8)        
          CALL PRINT(HA,NUMIJ,NORBS,NOUT)                               
        END IF                                                          
        DO 2 IBSET=1,NBSET                                              
          ITA=2*NDF+(IBSET-1)*NDF+IDF                                   
          CALL WREADW(IT94,TA,INTOWP(NUMIJ),JT1C(ITA),JUNK)             
          IF (IAND(IPRINT,32) .NE. 0) THEN                              
            WRITE(NOUT,1102) IDF,ITA,JT1C(ITA)                          
 1102       FORMAT(/'  TA FOR DEGREE ',I5,'   ITA=',I8,'  JT1C=',I8)    
            CALL PRINT(TA,NUMIJ,NORBS,NOUT)                             
          END IF                                                        
          F=FOCC(IBSET)/(2.0D+00)                                       
          DO 2 J=1,NUMIJ                                                
            BUF(J,IBSET)=HA(J)*F+TA(J)                                  
 2      CONTINUE                                                        
C                                                                       
CGES    WRITE(*,*)'DROL',DROL                                           
        DO 3 ITYP=1,NTYPES                                              
          IMIN=JPOINT(ITYP)                                             
          IMAX=JPOINT(ITYP+1)-1                                         
          DO 3 IYY=IMIN,IMAX                                            
            IBF=DROL(IYY)                                               
            IORB=IOUT(IBF)                                              
            DO 3 JYY=1,NORBS                                            
              JBF=DROL(JYY)                                             
              JORB=IOUT(JBF)                                            
              IJMN=MIN0(IORB,JORB)                                      
              IJMX=MAX0(IORB,JORB)                                      
              IJ=IOFF(IJMX)+IJMN                                        
              EPA(IORB,JORB)=BUF(IJ,ITYP)                               
 3      CONTINUE                                                        
        JT1C(IEA)=ISEA                                                  
        CALL WWRITW(IT94,EPA,INTOWP(NORBS*NORBS),ISEA,ISEA)             
        IF (IAND(IPRINT,32) .NE. 0) THEN                                
          WRITE(NOUT,1103) IDF,IEA,JT1C(IEA)                            
 1103     FORMAT(/'  EPA FOR DEGREE ',I5,'   IEA=',I8,'  JT1C=',        
     *    I8)                                                           
          CALL MATOUT(EPA,NORBS,NORBS,NORBS,NORBS,NOUT)                 
        END IF                                                          
 1    CONTINUE                                                          
C                                                                       
C     CODE FOR NBSET EPA FOR CI FIRST                                   
C                                                                       
      CALL RFILE(IT96)                                                  
      ISEA=1                                                            
      IADHA=JT1C(NDF+1)                                                 
      IADTA=JT1C(2*NDF+(NBSET-1)*NDF+1)                                 
      DO 312 IDF=1,NDF                                                  
      CALL WREADW(IT94,HA,INTOWP(NUMIJ),IADHA,IADHA)                    
      CALL WREADW(IT94,TA,INTOWP(NUMIJ),IADTA,IADTA)                    
C                                                                       
C     THIS MAY NEED TA*0.5D0   (but, then again, maybe not!)            
C                              (see, I thought so!)                     
         DO 311 IJ=1,NUMIJ                                              
           LSTEPA(IJ)=HA(IJ)+TA(IJ)                                     
  311    CONTINUE                                                       
      CALL WWRITW(IT96,LSTEPA,INTOWP(NUMIJ),ISEA,ISEA)                  
        IF (IAND(IPRINT,32) .NE. 0) THEN                                
          WRITE(NOUT,1104) IDF                                          
 1104     FORMAT(/'  LSTEPA FOR DEGREE ',I5)                            
          CALL PRINT(LSTEPA,NUMIJ,NORBS,NOUT)                           
        END IF                                                          
  312 CONTINUE                                                          
      RETURN                                                            
      END                                                               
