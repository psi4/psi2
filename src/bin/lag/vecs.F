      SUBROUTINE VECS(JOUT,ACRCY,WVEC,EIG,IFLAG,EIV,ia47)               
C***********************************************************************
C     THIS SUBROUTINE READS IN THE SCF VECTOR AND ADJUSTS IT TO THE    *
C     DRT REPRESENTATION.                                              *
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C                                                                       
      COMMON /FILES/ IT40,IT54,CCSD                                     
      COMMON/TAPES/ IOUT,ICHECK,IT30,IT52,IT53,IT58,IT91,IT93           
      COMMON /DIMS/ ISYMRB,NUMIJ,NIJ,NBF,NIJAO,NBFAO,NORBS,NMAX,NSYMHF  
      COMMON /PRNT/ IPRINT                                              
C                                                                       
      LOGICAL CCSD                                                      
      INTEGER JOUT(NBF),NLAMDA(8),NC(8),ITYP(8),IOR(500),ia47(192)      
      REAL*8 WVEC(NBFAO,NBF),ACRCY,EIV(NBFAO,NBF),EIG(NBF)              
C                                                                       
      it47 = 47
      WRITE(6,*)'CCSD =   ', CCSD                                       
C                                                                       
C10-18-90                                                               
C************************************                                   
C     READ IN AO-MO EIGENVALUES     *                                   
C************************************                                   
      CALL MREAD(EIG,10)                                                
C************************************                                   
C     READ IN AO-MO EIGENVECTOR     *                                   
C************************************                                   
      CALL MREAD(WVEC,13)                                               
      IF(IPRINT.GE.1) THEN                                              
        WRITE(6,*) 'AO-MO SCF EIGENVECTOR READ FROM 40'    
        CALL MATOUT(WVEC,NBFAO,NBF,NBFAO,NBF,IOUT)                      
      END IF                                                            
C************************************                                   
C     REORDER FOR CCSD GRADIENTS    *                                   
C************************************                                   
      IF(CCSD)THEN                                                      
        WRITE(*,*)' '                                                   
        WRITE(*,*)'CCSD OPTION, EIGENVALUES AND EIGENVECTORS ARE SORTED'
        WRITE(*,*)' '                                                   
        DO 10 I=1,NBF                                                   
          IOR(I)=I                                                      
   10   CONTINUE                                                        
   20   IFK=0                                                           
        I=0                                                             
   30   I=I+1                                                           
        IF(EIG(I).GT.EIG(I+1))THEN                                      
          VAL=EIG(I)                                                    
          EIG(I)=EIG(I+1)                                               
          EIG(I+1)=VAL                                                  
          IDU=IOR(I)                                                    
          IOR(I)=IOR(I+1)                                               
          IOR(I+1)=IDU                                                  
          IFK=1                                                         
        ENDIF                                                           
        IF(I.LT.(NBF-1))GO TO 30                                        
        IF(IFK.EQ.1) GO TO 20                                           
        DO 50 I=1,NBF                                                   
          NIB=IOR(I)                                                    
          DO 40 J=1,NBF                                                 
            EIV(J,I)=WVEC(J,NIB)                                        
   40     CONTINUE                                                      
   50   CONTINUE                                                        
        DO 60 I=1,NBF                                                   
          DO 60 J=1,NBF                                                 
            WVEC(I,J)=EIV(I,J)                                          
   60   CONTINUE                                                        
      ENDIF                                                             
C************************************                                   
C     REORDER TO DRT ORDERING       *                                   
C************************************                                   
      DO 80 IBF = 1,NBF                                                 
        NIBF = JOUT(IBF)                                                
        DO 70 IAO = 1,NBFAO                                             
          EIV(IAO,NIBF) = WVEC(IAO,IBF)                                 
   70   CONTINUE                                                        
   80 CONTINUE                                                          
C***********************************************                        
C     SWAP DRT ORDERED VECTOR BACK TO WVEC     *                        
C***********************************************                        
      DO 100 IBF = 1,NBF                                                
        DO 90 IAO = 1,NBFAO                                             
c         TEMP=WVEC(IAO,IBF)                                            
          WVEC(IAO,IBF) = EIV(IAO,IBF)                                  
c         EIV(IAO,IBF) = TEMP                                           
   90   CONTINUE                                                        
  100 CONTINUE                                                          
C********************************************************               
C     WRITE DRT ORDERED AO-MO VECTORS TO MASTER FILE    *               
C********************************************************               
      CALL MWRIT(WVEC,42)                                               
ctph
      leiv = ia47(160)
      lseiv = ia47(110)
      call wwritw(it47,wvec,leiv,lseiv,junk)
      IF(IPRINT.GE.1) THEN                                              
        WRITE(6,*) 'DRT ORDERED AO-MO SCF EIGENVECTOR WRITTEN TO 40'    
        CALL MATOUT(WVEC,NBFAO,NBF,NBFAO,NBF,IOUT)                      
      END IF                                                            
C                                                                       
      RETURN                                                            
      END                                                               
