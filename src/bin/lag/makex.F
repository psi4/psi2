      SUBROUTINE MAKEX(X,TPDM,INTS,BUF,IBUF,ORBSYM,                     
     1 BKTSZ,INTOFF,ILAST,ILAST2,ECINEW,BKTSZ2)                         
C                                                                       
      IMPLICIT INTEGER (A-Z)                                            
      INTEGER IBUF(BKTSZ2),ORBSYM(NORBS),INTOFF(NBKT),ILAST(NBKT)       
      INTEGER ILAST2(NBKT)                                              
      REAL*8 X(NBF,NBF),TPDM(LENR),VAL                                  
      REAL*8 INTS(LENR),BUF(BKTSZ),ECIOLD,ESCF,REP,ECINEW               
      REAL*8  ECI1,ECI2,ECITOT,Z0                                       
      LOGICAL CCSD                                                      
C                                                                       
      COMMON/TAPES/ IOUT,ICHECK,IT30,IT52,IT53,IT58,IT91,IT93           
      COMMON /FILES/ IT40,IT54,CCSD                                     
      COMMON /TUM/ IJSYM,IJSRT,II,JJ                                    
      COMMON /SRT/ VAL,IADR,IBKT,NBKT,NGRPIJ,LENR                       
      COMMON /LBLS/ ILBL(26),LBLDRT(26)                                 
      COMMON /DIMS/ ISYMRB,NUMIJ,NIJ,NBF,NIJAO,NBFAO,NORBS,NMAX,NSYMHF  
      COMMON /NREP/ REP,ECI1                                            
      COMMON /TRANS/ NMNB,BLEN,IBPT,ILPT                                
      COMMON /PRNT/ IPRINT                                              
      DATA Z0 /0.0D0/                                                   
C                                                                       
ctph  CALL INITUM(ORBSYM,X,TPDM,INTS,LENR)                              
      IF(IPRINT.GE.2) THEN                                              
        WRITE(6,*) 'ONE ELECTRON LAGRANGIAN AFTER INITUM IS CALLED'     
        CALL MATOUT(X,NBF,NBF,NBF,NBF,IOUT)                             
      END IF                                                            
      ECI2 = Z0                                                         
      ECITOT = Z0                                                       
      IBPT=1                                                            
      INTF=INTOFF(1)                                                    
      IIJJ=0                                                            
      IGRP=0                                                            
C********************************************                           
C     WRITE ILAST OUT TO FILE92 FOR NOW     *                           
C********************************************                           
      IT92 = 92                                                         
      CALL RFILE(IT92)                                                  
      CALL WWRITW(IT92,NBKT,1,1,IEND)                                   
      CALL WWRITW(IT92,ILAST,NBKT,IEND,IEND)                            
      CALL WWRITW(IT92,NGRPIJ,1,IEND,IEND)                              
      CALL RCLOSE(IT92,3)                                               
C*********************************                                      
C   START MAIN LOOP OVER MO'S    *                                      
C*********************************                                      
      DO 90 II=1,NORBS                                                  
        ISYM=ORBSYM(II)-1                                               
        DO 80 JJ=1,II                                                   
          IIJJ=IIJJ+1                                                   
          IJBKT = IIJJ/NGRPIJ                                           
          IF(IJBKT*NGRPIJ.NE.IIJJ) IJBKT = IJBKT + 1                    
          IJSRT = (IIJJ - (IJBKT-1)*NGRPIJ - 1)*NIJ                     
C**********************************************                         
C     CHECK TO SEE IF A NEW READ IS NEEDED    *                         
C**********************************************                         
          IF (IGRP .NE. IJBKT) THEN                                     
            IGRP=IJBKT                                                  
C*************************                                              
C     CHECK FOR ERROR    *                                              
C*************************                                              
            IF (IGRP .GT. NBKT) THEN                                    
              WRITE(ICHECK,*)' IN MAKEX II=',II,' JJ=',JJ,' IGRP=',IGRP,
     .                    ' NBKT=',NBKT                                 
              CALL MABORT                                               
            END IF                                                      
            IF(IPRINT.GE.4) THEN                                        
              WRITE(6,*) 'LAGRANGIAN AT THE FIRST PART OF IF'           
              CALL MATOUT(X,NBF,NBF,NBF,NBF,IOUT)                       
            END IF                                                      
            CALL ZERO(TPDM,LENR)                                        
            CALL ZERO(INTS,LENR)                                        
            IADR=ILAST(IGRP)                                            
            IF (IADR .EQ. 0) GO TO 40                                   
            IF (IADR.LT.0) THEN                                         
              PRINT*,'  FOR IGRP=',IGRP,'  IADR=',IADR                  
              CALL MABORT                                               
            END IF                                                      
   20       CONTINUE                                                    
C           WRITE(6,*) 'IADR',IADR                                      
            CALL RREAD(IT91,IBUF,BKTSZ2,IADR)                           
CTEST       CALL WREADW(IT91,IBUF,BKTSZ2,IADR,JUNK)                     
            IADR=IBUF(1)                                                
            NUM=IBUF(2)                                                 
c           WRITE(6,*) 'IADR,NUM',IADR,NUM                              
            DO 30 I=1,NUM                                               
              TPDM(IBUF(I+2))=BUF(I+INTF)                               
   30       CONTINUE                                                    
            IF(IPRINT.GE.5) THEN                                        
              WRITE(6,*) 'TPDM'                                         
              WRITE(6,65) (TPDM(IBUF(I+2)),I=1,NUM)                     
            END IF                                                      
   40       CONTINUE                                                    
            IF (IADR .NE. 0) GO TO 20                                   
            IADR=ILAST2(IGRP)                                           
            IF (IADR.EQ.0) GO TO 70                                     
            IF (IADR.LT.0) THEN                                         
              PRINT*,'  FOR IGRP=',IGRP,'  IADR=',IADR                  
              CALL MABORT                                               
            END IF                                                      
   50       CONTINUE                                                    
C           WRITE(6,*) 'IADR',IADR                                      
            CALL RREAD(IT93,IBUF,BKTSZ2,IADR)                           
CTEST       CALL WREADW(IT93,IBUF,BKTSZ2,IADR,JUNK)                     
            IADR=IBUF(1)                                                
            NUM=IBUF(2)                                                 
C     WRITE(6,*) 'LAGRANGIAN IN IF'                                     
C     WRITE(6,*) X                                                      
            DO 60 I=1,NUM                                               
              INTS(IBUF(I+2))=BUF(I+INTF)                               
   60       CONTINUE                                                    
            IF(IPRINT.GE.5) THEN                                        
              WRITE(6,*) 'INTS'                                         
              WRITE(6,65) (INTS(IBUF(I+2)),I=1,NUM)                     
   65         FORMAT(5F15.10)                                           
            END IF                                                      
            IF (IADR .NE. 0) GO TO 50                                   
   70       CONTINUE                                                    
c           WRITE(6,*) 'FIRST ELEMENTS OF TPDM'                         
c           WRITE(6,*) (TPDM(I),I=1,NIJ)                                
          END IF                                                        
          JSYM=ORBSYM(JJ)-1                                             
          IJSYM=IEOR(ISYM,JSYM)                                         
          CALL TUMULT(orbsym,x,tpdm,ints,tpdm,ints,lenr)                
   80   CONTINUE                                                        
        IF(IPRINT.GE.4) THEN                                            
          WRITE(6,*) 'LAGRANGIAN AFTER 80'                              
          CALL MATOUT(X,NBF,NBF,NBF,NBF,IOUT)                           
        END IF                                                          
   90 CONTINUE                                                          
      IF(IPRINT.GE.3) THEN                                              
        WRITE(6,*) 'LAGRANGIAN AFTER 90'                                
        CALL MATOUT(X,NBF,NBF,NBF,NBF,IOUT)                             
      END IF                                                            
      DO 130 I=1,NORBS                                                  
        ECITOT=ECITOT + X(I,I)                                          
  130 CONTINUE                                                          
      ECI2 = (ECITOT - 2.0D0*ECI1)*0.25D0                               
      IF(.NOT.CCSD) CALL DCDLBL(IT54,ILBL,ILBL,N,N,N,ESCF,ECIOLD,VAL)   
      ECINEW=ECI1+ECI2+REP                                              
      WRITE(6,140) ECI1,ECI2,REP,ECINEW,ECIOLD,ESCF                     
  140 FORMAT(1X,/' ONE ELECTRON CI ENERGY    = ',F25.15/                
     2           ' TWO ELECTRON CI ENERGY    = ',F25.15/                
     7           ' NUCLEAR REPULSION ENERGY  = ',F25.15/                
     3           ' TOTAL CI ENERGY           = ',F25.15/                
     4           ' CI ENERGY FROM TAPE 54    = ',F25.15/                
     5           ' SCF ENERGY FROM TAPE 54   = ',F25.15/)               
C                                                                       
      RETURN                                                            
      END                                                               
