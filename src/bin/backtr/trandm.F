      SUBROUTINE TRANDM(TPDM,BUF,EIV,GRPIJ,XTP,SQ,MAPIJ,ITPDM,IBUF,     
     .                  ILAST)                                          
C                                                                       
      IMPLICIT INTEGER (A-Z)                                            
      INTEGER MAPIJ(NIJAO),ITPDM(2),IBUF(2),ILAST(NBKT)                 
      REAL*8 TPDM(2),BUF(2),EIV(NBF,NBFAO),XTP(NBFAO,NBFAO)             
      REAL*8 SQ(NBFAO*NBFAO),Z0                                         
C                                                                       
      COMMON /DIMS/ NIJ,NBF,NIJAO,NBFAO,LENR,NST,NSHELL,AOSTRT,BKTSZ    
      COMMON/TAPES/ IOUT,ICHECK,IT30,IT55,IT91,IT95,IT96                
      COMMON /SRT/  NBKT,NGRPIJ,LGRPIJ,NBKT2,LGRP2,IBPT,IAPT,IGRP       
      DATA Z0 /0.0D0/                                                   
C*****************************                                          
C     READ IN MO TWO PDM     *                                          
C*****************************                                          
      INTF = ( (BKTSZ*2 - 2) / 6 ) + 1                                  
C     WRITE(6,*) 'INTF',INTF                                            
C     WRITE(6,*) 'IT91,IAPT,LENGTH',IT91,IAPT,INTOWP(NIJ*GRPIJ)         
      CALL ZERO(TPDM,NIJ*GRPIJ)                                         
      IADR = ILAST(IGRP)                                                
C     WRITE(6,*) 'IADR',IADR                                            
      IF(IADR.EQ.0) GO TO 140                                           
      IF(IADR.LT.0) THEN                                                
        PRINT*,'  FOR IGRP=',IGRP,'  IADR=',IADR                        
        CALL MABORT                                                     
      ENDIF                                                             
  120 CONTINUE                                                          
      CALL RREAD(IT91,IBUF,INTOWP(BKTSZ),IADR)                          
      IADR = IBUF(1)                                                    
      NUM = IBUF(2)                                                     
C     WRITE(6,*) 'IADR NUM',IADR,NUM                                    
      DO 130 I=1,NUM                                                    
        TPDM(IBUF(I+2)) = BUF(I+INTF)                                   
  130 CONTINUE                                                          
C     WRITE(6,*) 'IBUF IN IT91'                                         
C     WRITE(6,*) (IBUF(I),I=3,NUM+2)                                    
C     WRITE(6,*) 'TPDM IN IT91'                                         
C     WRITE(6,*) (TPDM(IBUF(I+2)),I=1,NUM)                              
  140 CONTINUE                                                          
      IF(IADR.NE.0) GO TO 120                                           
C     CALL WREADW(IT91,ITPDM,INTOWP(NIJ*GRPIJ),IAPT,IKPT)               
C*********************************                                      
C     HALF TRANSFORM TWO PDM     *                                      
C*********************************                                      
C     WRITE(6,*) 'IADR NUM',IADR,NUM                                    
      IAPT = IKPT                                                       
C     WRITE(6,*) 'UNTRANSFORMED TWO PDM'                                
C     WRITE(6,*) (TPDM(I),I=1,LENR)                                     
      CALL ZERO(BUF,LENR)                                               
      DO 40 IJ = 1,GRPIJ                                                
        OFSET = (IJ-1)*NIJ + 1                                          
        CALL SQUR(TPDM(OFSET),SQ,NBF)                                   
C       DO 11 I=1,NBF*NBF                                               
C         IF(SQ(I).GT.1.0D+10) THEN                                     
C           WRITE(6,*) 'TPDM TOO LARGE! AT IJ=',IJ,'IGRP=',IGRP         
C           WRITE(6,*) (SQ(J), J=1,NBF*NBF)                             
C           STOP                                                        
C         ENDIF                                                         
C  11   CONTINUE                                                        
        CALL CTSQC(EIV,SQ,XTP,SQ,NBF,NBFAO,NBFAO)                 
C       CALL ZERO(XTP,NBFAO*NBFAO)                                      
C       CALL MXMB(EIV,1,NBFAO,SQ,1,NBF,XTP,1,NBFAO,NBFAO,NBF,NBF)       
        MUNU = 0                                                        
        DO 30 MU = 1,NBFAO                                              
C         DO 10 NU = 1,MU                                               
C           SQ(NU) = Z0                                                 
C  10     CONTINUE                                                      
C         CALL MXMB(XTP,1,NBFAO,EIV(MU,1),NBFAO,1,SQ,1,NBFAO,MU,NBFAO,1)
          DO 20 NU = 1,MU                                               
            MUNU = MUNU + 1                                             
            NWMUNU = MAPIJ(MUNU)                                        
            INDEX = GRPIJ*(NWMUNU-1) + IJ                               
            BUF(INDEX) = SQ(MUNU)                                       
   20     CONTINUE                                                      
   30   CONTINUE                                                        
   40 CONTINUE                                                          
      CALL WWRITW(IT95,BUF,INTOWP(NIJAO*GRPIJ),IBPT,ILPT)               
C     WRITE(6,*) 'HALF TRANSFORMED TWO PDM'                             
C     WRITE(6,*) 'LENGTH,IBPT,ILPT',NIJAO*GRPIJ,IBPT,ILPT               
C     WRITE(6,*) (BUF(I),I=1,NIJAO*GRPIJ)                               
      IBPT = ILPT                                                       
C**************************************************************         
C     GO BACK TO BACK1 AND READ IN ANOTHER BUFFER OF TPDM     *         
C**************************************************************         
      RETURN                                                            
      END                                                               
