      SUBROUTINE LIORDR(I,J,K,L,VAL2,ioff,ijoff,ijbkt,
     *    bktpnt,start,ilast,bktsz2,srtar,isrtar,intoff,itxxx)
      IMPLICIT INTEGER (A-Z)                                            
      REAL*8 VAL,VAL2,srtar(ibig2)                                      
C                                                                       
      DIMENSION IOFF(NIJ),IJOFF(NIJ),IJBKT(NIJ),isrtar(ibig)            
      DIMENSION bktpnt(nbkt),start(nbkt),ilast(nbkt),intoff(nbkt)       
      COMMON /ORD/ IFAC1,IFAC2                                          
      common /bigd/ ibig,ibig2
      COMMON/TAPES/ IOUT,ICHECK,IT30,IT52,IT53,IT58,IT91,IT93           
      COMMON /DIMS/ ISYMRB,NUMIJ,NIJ,NBF,NIJAO,NBFAO,NORBS,NMAX,NSYMHF  
      COMMON /SRT/ VAL,IADR,IBKT,NBKT,NGRPIJ,LENR                       
      COMMON /PRNT/ IPRINT                                              
      IJ=IOFF(I)+J                                                      
      KL=IOFF(K)+L                                                      
      IBKT=IJBKT(IJ)                                                    
      IADR=IJOFF(IJ)+KL                                                 
      VAL=VAL2*IFAC1                                                    
      IF (IBKT.GT.NBKT .OR. IADR.GT.LENR) THEN                          
        WRITE(IOUT,*)I,J,K,L,IJ,KL,IBKT,IADR                           
        write(iout,*) 'nbkt',nbkt,'lenr',lenr
      END IF                                                            
      CALL SRTOUT(bktpnt,start,ilast,bktsz2,srtar,isrtar,intoff,itxxx)  
      IF (IJ .EQ. KL) RETURN                                            
      IBKT=IJBKT(KL)                                                    
      IADR=IJOFF(KL)+IJ                                                 
      VAL=VAL2*IFAC2                                                    
      IF (IBKT.GT.NBKT .OR. IADR.GT.LENR) THEN                          
        WRITE(IOUT,*)I,J,K,L,IJ,KL,IBKT,IADR                           
        write(iout,*) 'nbkt',nbkt,'lenr',lenr
      END IF                                                            
      CALL SRTOUT(bktpnt,start,ilast,bktsz2,srtar,isrtar,intoff,itxxx)  
C                                                                       
      RETURN                                                            
      END                                                               
