      Subroutine fentry(A,IA,ACORE)                                              
C********************************************************************** 
C*   NOTICE OF PROGRAM MODIFICATION                                   * 
C********************************************************************** 
c       by: Edward Seidl (cets072091)
c     date: July 20, 1991
c   reason: add option DELETE91 for debugging, also call mmult in
c           subroutine ctsqc
c
C*  By:  Richard Remington                      Search:  CMACHINESIZE * 
C*  Date:  November  17, 1989                                         * 
C*  Reason: Adjust core to run in 7MB on 9370 or 15MB on 3090.        * 
C********************************************************************** 
C*  By:  Tracy P. Hamilton                      Search:  CMACHINESIZE * 
C*  Date:  November  21, 1989                                         * 
C*  Reason: Make SORT2 require less space.  Also, BACK2 has now been  * 
C*          written to not transform most of the redundant elements.  * 
C*          The subroutine DIME figures out the space for SORT2.      * 
C********************************************************************** 
C                                                                       
C********************************************************************** 
C     THIS PROGRAM WILL TRANSFORM THE MO  ONE AND TWO PDM'S IN DRT    * 
C     ORDER TO AO BASIS IN SHELL CANONICAL ORDER.                     * 
C                                                                     * 
C     THIS SPLITS THE BACK TRANSFORMATION FROM THE LAGRANGIAN         * 
C     CONSTRUCTION SO THAT FEWER LARGE FILES ARE NEEDED AT A TIME.    * 
C                                                                     * 
C     THIS IS NOT VERY SOPHISTICATED YET:  BACK2 CAN BE WRITTEN SO    * 
C     THAT NOT ALL OF A MATRIX NIJ ELEMENTS ARE TRANSFORMED.          * 
C     IN THE FUTURE, MAPIJ MAY CHANGE TO A DIFFERENT SHELL ORDER      * 
C     SO THAT ROTATIONAL INVARIANCE COULD POSSIBLY BE APPLIED.        * 
C                                                                     * 
C      BY:  TRACY P. HAMILTON                                         * 
C    DATE:  JUNE 29, 1989.                                            * 
C                                                                     * 
C********************************************************************** 
C                                                                     * 
C     THE FOLLOWING FEATURES HAVE BEEN INCORPORATED:                  * 
C                                                                     * 
C         CORE CHECKING ROUTINES FROM CLJ.                            * 
C         CHANGE FILE55 IO TO INTEGER WORD ADDRESSING.                * 
C         THIS VERSION COMPATIBLE ONLY WITH DERIVTH.                  * 
C         CC SYMM GRADS AND MONGO. EIGENVECTOR READ FROM FILE30.      * 
C         RE-SORTING TO CC ORDERING TO TAKE ADVANTAGE OF SYMMETRY.    * 
C         USE YUKIO'S MASTER FILE40 INSTEAD OF FILE47, SO SAME        * 
C           SCF PROGRAMS CAN BE USED.                                 * 
C         ABILITY TO DO TRANSFORMATION FOR NBFAO NOT EQUAL NBF.       * 
C                                                                     * 
C********************************************************************** 
C                                                                     * 
C DIRECTORY OF FILES                                                  * 
C                                                                     * 
C ICHECK - THE OUTPUT FILE WHERE ERRORS AND MESSAGES ARE SENT.        * 
C IOUT - THE OUTPUT FILE.                                             * 
C IT40 - YUKIO'S MASTER FILE.                                         * 
C IT55 - CONTAINS THE TWOPDM IN THE AO BASIS.                         * 
C IT91 - WORK FILE USED TO STORE THE MO TWOPDM.  MADE IN PROGRAM LAG  * 
C IT92 - info from lag on i/o buffers and ngrpij.                     * 
C IT95 - WORK FILE USED TO STORE THE HALF-TRANSFORMED TWOPDM.         * 
C IT96 - WORK FILE USED TO STORE THE TRANSFORMED TWOPDM BEFORE SORT.  * 
C                                                                     * 
C********************************************************************** 
      IMPLICIT INTEGER (A-Z)                                            
cets030691
#include <error.h>
      integer errcod,frdi,frdboo
C                                                                       
      COMMON /DIMS/ NIJ,NBF,NIJAO,NBFAO,LENR,NST,NSHELL,AOSTRT,BKTSZ    
      COMMON/TAPES/ IOUT,ICHECK,IT30,IT55,IT91,IT95,IT96                
      COMMON /SRT/  NBKT,NGRPIJ,LGRPIJ,NBKT2,LGRP2,IBPT,IAPT,IGRP       
C                                                                       
      INTEGER IA(ACORE*2)                                                  
      REAL*8 A(1)                                                       
C                                                                       
      CALL psinit('APPEND')                                             
      CALL TSTART(6)                                                    
      core = intowp(acore)
      iprint = 0
      errcod = frdi('IPRINT',iprint)
cets072091
      idl91=1
      errcod = frdboo('DELETE91',idl91)
C*************************************************                      
C   INITIALIZE THE FILES USED BY THIS PROGRAM    *                      
C*************************************************                      
      errcod = frdi('BUCKET_SIZE',BKTSZ)
      if (BKTSZ.eq.0) BKTSZ = 4096 
      errcod = frdi('MAX_BUCKET',MXBKT)
      if (MXBKT.eq.0) MXBKT = 500 
      ICHECK= 3                                                         
      INPUT = 5                                                         
      IOUT  = 6                                                         
      IT30  = 30                                                        
      IT55  = 55                                                        
      IT91  = 91                                                        
      IT95  = 95                                                        
      IT96  = 96                                                        
      CALL RFILE(IT30)                                                  
      CALL RFILE(IT55)                                                  
      CALL RFILE(IT91)                                                  
      CALL SREW(IT91)                                                   
      CALL RFILE(IT95)                                                  
      CALL RFILE(IT96)                                                  
C********************                                                   
C   PRINT HEADER    *                                                   
C********************                                                   
      WRITE(IOUT,5)                                                     
    5 FORMAT(//,10X,'  PROGRAM TO TRANSFORM THE TWO PARTICLE  '/        
     1          10X,'     DENSITY MATRIX BACK TO AO BASIS '/            
     2          10X,'          WRITTEN BY T. HAMILTON  '//              
     3          10X,'          VERSION JUNE 29, 1989'/)                 
      WRITE(IOUT,10)                                                    
   10 FORMAT(/,1X,'USING NEW FILE55 IO (INTEGER WORD ADDRESSING)')      
      WRITE(IOUT,11)                                                    
   11 FORMAT(/,1X,'CAPABLE OF TRANSFORMING WITH 5D, 7F ETC. IN BASIS')  
      WRITE(IOUT,15) CORE                                               
   15 FORMAT(/,1X,'USING ',I10,' INTEGER WORDS OF CORE.')               
C************************************                                   
C     INITIALIZE MASTER FILE 40     *                                   
C************************************                                   
      CALL INITMF(1)                                                    
C************************************                                   
C     READ CONSTANTS FROM FILE30    *                                   
C************************************                                   
C     WRITE(6,*) 'IT30,101',IT30,101                                    
      CALL WREADW(IT30,IA(1),200,101,JUNK)                              
      NSHELL = IA(27)                                                   
      NBFAO  = IA(22)                                                   
      NBF    = IA(18)                                                   
      NIJ    = (NBF*(NBF+1))/2                                          
      NIJAO  = NBFAO*(NBFAO+1)/2                                        
      NST    = NSHELL*(NSHELL+1)/2                                      
      WRITE(6,*) ' NBF    =',NBF                                        
      WRITE(6,*) ' NIJ    =',NIJ                                        
      WRITE(6,*) ' NBFAO  =',NBFAO                                      
      WRITE(6,*) ' NIJAO  =',NIJAO                                      
      WRITE(6,*) ' NSHELL =',NSHELL                                     
      WRITE(6,*) ' NST    =',NST                                        
C********************************                                       
C     ZERO OUT BLANK COMMON     *                                       
C********************************                                       
      DO 30 II=1,CORE                                                   
        IA(II)=0                                                        
   30 CONTINUE                                                          
C*********************************                                      
C     ALLOCATE CORE FOR READ30   *                                      
C*********************************                                      
      CALL COINIT(IA,CORE)                                              
      KTYPE  = COALLC(IA,isizof('INT')*NSHELL)                          
      KLOC   = COALLC(IA,isizof('INT')*NSHELL)                          
      KREST  = COALLC(IA,isizof('INT')*200)                             
c     WRITE(6,*) 'KTYPE,KLOC,KREST',KTYPE,KLOC,KREST                    
C**********************                                                 
C     READ TAPE 30    *                                                 
C**********************                                                 
      CALL READ30(IA(KTYPE),IA(KLOC),IA(KREST))                         
C*********************************                                      
C     ALLOCATE CORE FOR SETUP1   *                                      
C*********************************                                      
      CALL COFREE(IA,KREST)                                             
      AOPS   = COALLC(IA,isizof('INT')*NSHELL)                          
      MAPIJ  = COALLC(IA,isizof('INT')*NIJAO)                           
      IOFF   = COALLC(IA,isizof('INT')*NIJAO)                           
      ILAST  = COALLC(IA,isizof('INT')*MXBKT)                           
C*********************************                                      
C   READ IN ILAST FROM FILE92    *                                      
C*********************************                                      
      IT92 = 92                                                         
      CALL RFILE(IT92)                                                  
      CALL WREADW(IT92,NBKT,1,1,JUNK)                                   
      CALL WREADW(IT92,IA(ILAST),NBKT,JUNK,JUNK)                        
      CALL WREADW(IT92,NGRPIJ,1,JUNK,JUNK)                              
      CALL RCLOSE(IT92,3)                                               
c     WRITE(6,*) 'ILAST',(IA(ILAST+I-1),I=1,NBKT)                       
      WRITE(6,*) 'NGRPIJ',NGRPIJ                                        
C*******************************                                        
C     SETUP AOPS AND MAPIJ     *                                        
C*******************************                                        
      CALL SETUP1(IA(AOPS),IA(KTYPE),IA(MAPIJ),IA(IOFF),IA(KLOC))       
C************************************************                       
C     READ DRT ORDERED EIGENVECTOR INTO EIV    *                        
C************************************************                       
      EIV    = COALLC(IA,isizof('REAL')*NBFAO*NBF)                      
      CALL MREAD(IA(EIV),42)                                            
C***********************************                                    
C     ALLOCATE ARRAYS FOR INIT55   *                                    
C***********************************                                    
      IJOFF  = COALLC(IA,isizof('INT')*NIJAO)                           
      KLOFF  = COALLC(IA,isizof('INT')*NIJAO)                           
      LNSR   = COALLC(IA,isizof('INT')*NST)                             
      CALL INIT55(IA(AOPS),IA(LNSR),IA(IJOFF),IA(KLOFF))                
      CALL COFREE(IA,IJOFF)                                             
      CALL COFREE(IA,KLOFF)                                             
      CALL COFREE(IA,LNSR)                                              
      if(iprint.ge.2) then
       CALL CODUMP(IA)                                                   
      endif
C***********************************                                    
C     ALLOCATE ARRAYS FOR BACK1    *                                    
C***********************************                                    
C*********************************************************************  
C     NGRPIJ IS NUMBER OF ARRAYS OF DIMENSION NIJAO CAN BE HELD IN   *  
C     CORE AT ONCE.  LEN IS CORE USED BY THE NGRPIJ ARRAYS.  NBKT    *  
C     IS THE NUMBER OF TIMES YOU HAVE TO FILL CORE WITH NIJAO*NGRPIJ *  
C*********************************************************************  
      SQ     = COALLC(IA,isizof('REAL')*NBFAO*NBFAO)                    
      XTP    = COALLC(IA,isizof('REAL')*NBFAO*NBFAO)                    
      SPACE  = CORE - ( XTP + NBFAO*NBFAO )                             
CTPH  NGRPIJ = SPACE/INTOWP(NIJAO)                                      
CTPH  IF(NGRPIJ.GT.NIJ) NGRPIJ=NIJ                                      
CTPH  IF(NGRPIJ.GT.100) NGRPIJ=100                                      
      LENR   = NGRPIJ*NIJAO                                             
      NBKT   = NIJ/NGRPIJ                                               
      NBKT2  = NIJAO/NGRPIJ                                             
      LGRPIJ = NIJ - (NBKT*NGRPIJ)                                      
      LGRP2  = NIJAO - (NBKT2*NGRPIJ)                                   
      IF(LGRPIJ.NE.0) THEN                                              
        NBKT = NBKT + 1                                                 
      ELSE                                                              
        LGRPIJ = NGRPIJ                                                 
      END IF                                                            
      IF(LGRP2.NE.0) THEN                                               
        NBKT2 = NBKT2 + 1                                               
      ELSE                                                              
        LGRP2 = NGRPIJ                                                  
      END IF                                                            
C                                                                       
      TOLD   = COALLC(IA,isizof('REAL')*NIJ*NGRPIJ)                     
      IF(BKTSZ.GT.LENR) THEN                                            
        TNEW   = COALLC(IA,isizof('REAL')*BKTSZ)                        
      ELSE                                                              
        TNEW   = COALLC(IA,isizof('REAL')*LENR)                         
      ENDIF                                                             
C*********************************************************************  
C     DO THE FIRST HALF TRANSFORMATION:  FILE91 INTO FILE 95         *  
C*********************************************************************  
      CALL EIVT(IA(EIV),IA(XTP),IA(EIV),NBF,NBFAO)                      
      if(iprint.ge.2) then
       CALL CODUMP(IA)                                                   
      endif
      CALL BACK1(IA(TOLD),IA(TNEW),IA(XTP),IA(EIV),IA(SQ),IA(MAPIJ),    
     .           IA(ILAST))                                             
C***********************************************************            
C     WE NO LONGER NEED WORK FILE 91 SO GET RID OF IT.     *            
C***********************************************************            
      if(idl91.eq.1) then
        CALL RCLOSE(IT91,4)                                               
      else
        CALL RCLOSE(IT91,3)                                               
      endif
C*********************************************************************  
C     SECOND HALF OF TRANSFORMATION OF TWO PDM TO AO BASIS,          *  
C     FROM FILE 95 TO 96.                                            *  
C*********************************************************************  
      CALL COFREE(IA,TOLD)                                              
      CALL COFREE(IA,TNEW)                                              
      if(iprint.ge.2) then
       CALL CODUMP(IA)                                                   
      endif
      TOLD   = COALLC(IA,isizof('REAL')*LENR)                           
      TNEW   = COALLC(IA,isizof('REAL')*LENR)                           
C************************************                                   
C     FIRST MAKE EIV TRANSPOSE      *                                   
C************************************                                   
      CALL BACK2(IA(TOLD),IA(TNEW),IA(MAPIJ),IA(SQ),IA(XTP),IA(EIV),    
     .           IA(TOLD),IA(TNEW),IA(AOPS))                            
C***********************************                                    
C     ALLOCATE ARRAYS FOR SORT2    *                                    
C***********************************                                    
      CALL dime(IA(AOPS),mxdim1,mxdim2)                                 
      if(core.lt.nst+nst+40+intowp(mxdim1+mxdim2)) then
       write(6,*) 'ERROR - not enough core for sort2'  
       write(6,*) 'There are',core,'integer words available and' 
       write(6,*) 'There are',ilast + intowp(nst+mxdim1+mxdim2),     
     1            'integer words needed.' 
       call mabort
      endif
C***********************************************************            
C     WE NO LONGER NEED WORK FILE 95 SO GET RID OF IT.     *            
C***********************************************************            
      CALL RCLOSE(IT95,4)                                               
      CALL COFREE(IA,ILAST)                                             
      CALL COFREE(IA,EIV)                                               
      CALL COFREE(IA,SQ)                                                
      CALL COFREE(IA,XTP)                                               
      CALL COFREE(IA,TOLD)                                              
      CALL COFREE(IA,TNEW)                                              
      if(iprint.ge.2) then
       CALL CODUMP(IA)                                                   
      endif
      ISIZE  = COALLC(IA,isizof('INT')*NST)                             
      IREAD  = COALLC(IA,isizof('INT')*NST+1)                           
      TOLD   = COALLC(IA,isizof('REAL')*mxdim1)                         
      TNEW   = COALLC(IA,isizof('REAL')*mxdim2)                         
C***********************************************************            
C     SORT TRANSFORMED TWO PDM IN FILE 96 INTO FILE 55.    *            
C***********************************************************            
C     WRITE(6,*) 'IOFF'                                                 
C     WRITE(6,*) IA(IOFF),IA(IOFF+1),IA(IOFF+2),IA(IOFF+3),IA(IOFF+4)   
      CALL SORT2(IA(TOLD),IA(TNEW),IA(AOPS),IA(ISIZE),IA(IREAD),        
     .           IA(IOFF),IA(MAPIJ),mxdim1,mxdim2)                      
C***********************************************************            
C     WE NO LONGER NEED WORK FILE 96 SO GET RID OF IT.     *            
C***********************************************************            
      CALL RCLOSE(IT96,4)                                               
C******************                                                     
C     THE END     *                                                     
C******************                                                     
      WRITE(IOUT,40)                                                    
  40  FORMAT(//,5X,' >>>> TWO PDM BACK TRANSFORMATION COMPLETED <<<< '/)
      if(iprint.ge.1) then
       CALL CODUMP(IA)                                                   
      endif
C                                                                       
      JCODE = 3                                                         
      CALL RCLOSE(IT30,JCODE)                                           
      CALL RCLOSE(40,JCODE)                                             
      CALL RCLOSE(IT55,JCODE)                                           
C                                                                       
      CALL TSTOP(IOUT)                                                  
      STOP                                                              
      END                                                               
