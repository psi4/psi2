      SUBROUTINE fentry(A,IA,MAXCOR)
C
C  PROGRAM TO FORM 'LAGRANGIAN' MATRIX :
C                    L(I,M)=H(I,J)*P(M,J)+2G(IJKL)*Q(MJKL)
C  *******  BRB VERSION BEGUN ON MARCH 14,1983  ********
C
      IMPLICIT INTEGER (A-Z)
cets030691
#include <error.h>
      integer ip,prcntr,errcod,frdboo,frdi
      INTEGER NSOFF(8),ia47(192)
      INTEGER IA(2*maxcor)
      REAL*8 A(maxcor)
      logical oldinp
      COMMON /BIGD/ IBIG,IBIG2
      COMMON/TAPES/ IOUT,ICHECK,IT30,IT52,IT53,IT58,IT91,IT93
      COMMON /DIMS/ ISYMRB,NUMIJ,NIJ,NBF,NIJAO,NBFAO,NORBS,NMAX,NSYMHF
      COMMON /FILES/ IT40,IT54,CCSD
      COMMON /SRT/ VAL,IADR,IBKT,NBKT,NGRPIJ,LENR         
      COMMON /LBLS/ ILBL(26),LBLDRT(26)
      COMMON /TRANS/ NMNB,BLEN,IBPT,ILPT
      COMMON /PRNT/ IPRINT                
      LOGICAL CCSD
      REAL*8 DRTVER,VAL,ECINEW,ACRCY
      DATA ACRCY /1.0D-20/
C***********************************************************************
C  NOTICE OF PROGRAM MODIFICATIONS                                     *
C***********************************************************************
C  BY:  TRACY P. HAMILTON                                              *
C  DATE:  MAY 8, 1991                                                  *
C  REASON: Put write to file47 back in out1.f and out2.f               *
C          Add in Ed Seidl's new input and replace drum with psinit    *
C***********************************************************************
C  BY:  RICHARD REMINGTON                         SEARCH:  C10-18-90   *
C  DATE:  OCTOBER   18, 1990                                           *
C  REASON: ADD MREAD(EIG,10) AS PER INSTRUCTIONS BY T.P. HAMILTON.     *
C***********************************************************************
C  BY:  TRACY P. HAMILTON                                              *
C  DATE:  APRIL 12, 1990                                               *
C  REASON: PUT IN FENTRY                                               *
C***********************************************************************
C  BY:  TRACY P. HAMILTON                                              *
C  DATE:  JANUARY 21, 1990                                             *
C  REASON: PUT IN CORE CHECKING ROUTINES AND FIX MYSTERIOUS BUG.       *
C          PUT ADRV=OFSET+PNT AS IN ORIGINAL LAGTR.                    *
C***********************************************************************
C  BY:  RICHARD REMINGTON                      SEARCH:  CMACHINESIZE   *
C  DATE:  NOVEMBER  17, 1989                                           *
C  REASON: ADJUST CORE TO RUN IN 7MB ON 9370 OR 15MB ON 3090.          *
C***********************************************************************
C                                                                       
C***********************************************************************
C     THIS PROGRAM WILL CALCULATE THE LAGRANGIAN                       *
C     TO BE VECTORIZED IN THE NEAR FUTURE, AS WELL AS BETTER I/O       *
C     L(I,M)=H(I,J)*P(M,J)+2G(IJKL)*Q(MJKL)                            *
C                                                                      *
C     THIS IS A MAJOR REWRITE OF THE OLD LAGTRTH CODE IN ORDER TO      *
C     HAVE CLEAN CODE AND ELIMINATE AN ELUSIVE BUG THAT SHOWS UP       *
C     ONLY FOR NBF > 60.                                               *
C                                                                      *
C      BY:  TRACY P. HAMILTON                                          *
C    DATE:  JUNE 29, 1989.                                             *
C                                                                      *
C***********************************************************************
C                                                                      *
C     THE FOLLOWING FEATURES HAVE BEEN INCORPORATED:                   *
C                                                                      *
C         THIS VERSION COMPATIBLE ONLY WITH DERIVTH.                   *
C         GET RID OF MO PDMS (FILE53) AFTER FORMING FILE91.            *
C         NOT GETTING RID OF FILE91 AFTER PROGRAM IS DONE.             *
C         CC SYMM GRADS AND MONGO. EIGENVECTOR READ FROM FILE30.       *
C         RE-SORTING TO CC ORDERING TO TAKE ADVANTAGE OF SYMMETRY.     *
C         USE YUKIO'S MASTER FILE40 INSTEAD OF FILE47, SO SAME         *
C           SCF PROGRAMS CAN BE USED.                                  *
C         CURT JANSSEN'S CORE CHECKING ROUTINES                        *
C                                                                      *
C***********************************************************************
C                                                                      *
C DIRECTORY OF FILES                                                   *
C                                                                      *
C ICHECK - THE OUTPUT FILE WHERE ERRORS AND MESSAGES ARE SENT.         *
C IOUT - THE OUTPUT FILE.                                              *
C IT40 - YUKIO'S MASTER FILE.                                          *
C IT52 - CONTAINS THE MO ONE AND TWO ELECTRON INTEGRALS                *
C        IN DRT ORDERING.                                              *
C IT53 - CONTAINS THE ONE AND TWO PDM'S IN THE MO BASIS.               *
C IT58 - CONTAINS THE DRT INFORMATION.                                 *
C IT91 - WORK FILE USED TO STORE THE MO TWOPDM.                        *
C IT93 - WORK FILE USED TO STORE THE MO TWO ELECTRON INTEGRALS.        *
C                                                                      *
C***********************************************************************
C                                                                       
CMachineSize  (For 15Mb)  PARAMETER (CORE = 2800000)                    
CMachineSize  (For  7Mb)  PARAMETER (CORE = 1200000)                    
ctph  PARAMETER (CORE = 1200000)                                        
C                                                                       
      mxcori = maxcor*2
      CORE = MXCORI                                                     
C                                                                       
      CALL psinit('APPEND')                                             
C      CALL TSTART(3)                                                    
      CALL TSTART(6)                                                    
C*************************************************                      
C   INITIALIZE THE FILES USED BY THIS PROGRAM    *                      
C*************************************************                      
      ICHECK= 3                                                         
      INPUT = 5                                                         
      IOUT  = 6                                                         
      IT30  = 30                                                        
      IT40  = 40                                                        
      it47  = 47
      IT52  = 52                                                        
      IT53  = 53                                                        
      IT54  = 54                                                        
      IT58  = 58                                                        
      IT91  = 91                                                        
      IT93  = 93                                                        
      CALL RFILE(IT30)                                                  
      CALL RFILE(IT47)                                                  
      CALL RFILE(IT52)                                                  
      CALL RFILE(IT53)                                                  
      IF(.NOT.CCSD) CALL RFILE(IT54)                                    
      CALL RFILE(IT58)                                                  
      CALL RFILE(IT91)                                                  
      CALL RFILE(IT93)                                                  
C********************                                                   
C   PRINT HEADER    *                                                   
C********************                                                   
      WRITE(IOUT,5)                                                     
    5 FORMAT(//,10X,'  PROGRAM TO CONSTRUCT LAGRANGIAN MATRIX '/        
     1          10X,'      BASED ON THE OLD LAGTR PROGRAM'/             
     2          10X,'          WRITTEN BY T. HAMILTON '//               
     3          10X,'          VERSION JUNE 29, 1989'/)                 
      WRITE(IOUT,10)                                                    
   10 FORMAT(/,1X,'USING NEW FILE55 IO (INTEGER WORD ADDRESSING)')      
      WRITE(IOUT,15) CORE                                               
   15 FORMAT(/,1X,'USING ',I10,' INTEGER WORDS OF CORE.')               
C**************************************                                 
C   LOCATE LAG INPUT IN INPUT FILE    *                                 
C**************************************                                 
      ccsd = .false.
      CALL LOCATE(INPUT,'# LAGTR #####',IERR)                              
      IF (IERR .EQ. 0) THEN                                             
        oldinp = .true.
        READ(INPUT,'(3I5)') IPRINT,IDFILE,ICCSD                         
        IF(ICCSD.EQ.1) CCSD = .TRUE.                                      
      ELSE                                                              
        oldinp = .false.
        IPRINT = 0                                                      
        errcod = frdi('IPRINT',iprint)
cets030291
        idfile = 1
        errcod = frdboo('DELETE53',idfile)
      END IF                                                            
      CALL INITMF(1)                                                    
      CALL WREADW(IT47,ia47,192,1,junk)                                 
C****************************************                               
C     READ DIMENSIONS FROM DRT TAPE     *                               
C****************************************                               
      CALL WREADW(IT58,DRTVER,INTOWP(1),1,END)                          
      CALL WREADW(IT58,LBLDRT,26,END,END)                               
      CALL WREADW(IT58,NBF,1,END,END)                                   
      CALL WREADW(IT58,NSYMHF,1,END,END)                                
      CALL WREADW(IT58,NORBS,1,END,END)                                 
      END = END + 10                                                    
      CALL WREADW(IT58,ISYMRB,1,END,END)                                
      CALL WREADW(IT58,NUMIJ,1,END,END)                                 
      CALL WREADW(IT58,NGROUP,1,END,END)                                
      END = END + 1                                                     
      CALL WREADW(IT58,NMAX,1,END,END)                                  
      END=END+3                                                         
C***********************************************************************
C  WRITE THE DRT VERSION AND OTHER DRT INFORMATION TO THE OUTPUT FILE  *
C***********************************************************************
      WRITE(IOUT,20) DRTVER                                             
   20 FORMAT(/1X,'DRTVER  = ',F10.6)                                    
      WRITE(IOUT,25) NBF,NSYMHF,NORBS,NMAX,NUMIJ,NGROUP,ISYMRB          
   25 FORMAT( 1X,'NBF     = ',I8/                                       
     1       ,1X,'NSYMHF  = ',I8/                                       
     2       ,1X,'NORBS   = ',I8/                                       
     3       ,1X,'NMAX    = ',I8/                                       
     4       ,1X,'NUMIJ   = ',I8/                                       
     5       ,1X,'NGROUP  = ',I8/                                       
     6       ,1X,'ISYMRB  = ',I8/)                                      
C******************************************                             
C     CHECK FOR FROZEN CORE OR VIRTUAL    *                             
C******************************************                             
      IF(NBF.NE.NORBS) THEN                                             
        WRITE(IOUT,*) 'ERROR! CI GRADIENTS WILL NOT RUN WITH FZV OR FZC'
        CALL MABORT                                                     
      END IF                                                            
C**************************************                                 
C     GET INFO FROM TAPE30, PRINT     *                                 
C**************************************                                 
      CALL WREADW(IT30,IA,200,101,JUNK)                                 
      NSHELL = IA(27)                                                   
      NBFAO  = IA(22)                                                   
      NIJ    = (NBF*(NBF+1))/2                                          
      NIJAO  = NBFAO*(NBFAO+1)/2                                        
      WRITE(IOUT,26) NIJ,NIJAO,NSHELL                                   
   26 FORMAT( 1X,'NIJ     = ',I8/                                       
     1       ,1X,'NIJAO   = ',I8/                                       
     2       ,1X,'NSHELL  = ',I8)                                       
C********************************                                       
C     ZERO OUT BLANK COMMON     *                                       
C********************************                                       
      DO 30 II=1,CORE                                                   
        IA(II)=0                                                        
   30 CONTINUE                                                          
      CALL COINIT(IA,CORE)                                              
C********************************************************               
C     ALLOCATE SOME CORE FOR THE DRT INTEGER ARRAYS     *               
C********************************************************               
      IOFF = COALLC(IA,isizof('INT')*NIJ)                               
      KADD = COALLC(IA,isizof('INT')*ISYMRB)                            
      LADD = COALLC(IA,isizof('INT')*ISYMRB)                            
      IJADD = COALLC(IA,isizof('INT')*NUMIJ)                            
      IJGRP = COALLC(IA,isizof('INT')*NUMIJ)                            
      ORBTBF = COALLC(IA,isizof('INT')*NORBS)                           
      JOUT = COALLC(IA,isizof('INT')*NBF)                               
      ORBSYM = COALLC(IA,isizof('INT')*NORBS)                           
C******************************************************************     
C     END OF POINTERS FOR DRT.  READ IN THE DRT ARRAYS NEEDED.    *     
C******************************************************************     
      CALL GETDRT(IA(KADD),IA(LADD),IA(IJADD),IA(IJGRP),                
     .       IA(ORBTBF),IA(JOUT),IA(ORBSYM),END,NGROUP)                 
C********************************************************************** 
C   MXBKT AND BKTSZ ARE ARBITRARY PARAMETERS WHOSE VALUES SHOULD BE   * 
C   OPTIMIZED FOR THE CONDITIONS UNDER WHICH THE PROGRAM IS RUN.      * 
C********************************************************************** 
      MXBKT  = 500                                                      
      BKTSZ  = 4096                                                     
      BKTSZ2 = INTOWP(4096)                                             
C***********************************                                    
C     ALLOCATE ARRAYS FOR INITG    *                                    
C***********************************                                    
      IJOFF = COALLC(IA,isizof('INT')*NIJ)                              
      IJBKT = COALLC(IA,isizof('INT')*NIJ)                              
C***********************************                                    
C     ALLOCATE ARRAYS FOR GETIT    *                                    
C***********************************                                    
      X = COALLC(IA,isizof('REAL')*NBF*NBF)                             
      ILAST = COALLC(IA,isizof('INT')*MXBKT)                            
      ILAST2 = COALLC(IA,isizof('INT')*MXBKT)                           
      INTOFF = COALLC(IA,isizof('INT')*MXBKT)                           
      H = COALLC(IA,isizof('REAL')*NBF*NBF)                             
      F = COALLC(IA,isizof('REAL')*NBF*NBF)                             
      START = COALLC(IA,isizof('INT')*MXBKT)                            
      BKTPNT = COALLC(IA,isizof('INT')*MXBKT)                           
      BUF = COALLC(IA,isizof('REAL')*NMAX)                              
      SPACE  = CORE - F - 2*NBF*NBF - BKTSZ2                            
C*********************************************************************  
C     NGRPIJ IS NUMBER OF ARRAYS OF DIMENSION NIJ THAT CAN BE HELD IN*  
C     CORE AT ONCE.  LENR IS CORE USED BY THE NGRPIJ ARRAYS.  NBKT   *  
C     IS THE NUMBER OF TIMES YOU HAVE TO FILL CORE WITH NIJ*NGRPIJ   *  
C*********************************************************************  
      NGRPIJ = SPACE/(2*INTOWP(NIJ))                                    
      IF(NGRPIJ.GT.NIJ) NGRPIJ=NIJ                                      
      IF(NGRPIJ.GT.100) NGRPIJ=100                                      
      LENR   = NGRPIJ*NIJ                                               
      NBKT   = NIJ/NGRPIJ                                               
      LGRPIJ = NIJ - (NBKT*NGRPIJ)                                      
      IF(LGRPIJ.NE.0) THEN                                              
        NBKT = NBKT + 1                                                 
      ELSE                                                              
        LGRPIJ = NGRPIJ                                                 
      END IF                                                            
      IF(NBKT.GT.MXBKT) THEN                                            
        WRITE(IOUT,40) NBKT,MXBKT                                       
   40   FORMAT(/5X,'NBKT = ',I5,'  MXBKT = ',I5,'INCREASE MXBKT.')      
        CALL MABORT                                                     
      END IF                                                            
C                                                                       
      IBIG   = BKTSZ2*NBKT                                              
      IBIG2  = IBIG/2                                                   
      WRITE(6,*)  'NGRPIJ =',NGRPIJ                                     
      WRITE(6,*)  'LGRPIJ =',LGRPIJ                                     
      WRITE(6,*)  'LENR   =',LENR                                       
      WRITE(6,*)  'IBIG   =',IBIG                                       
      ISRTAR = COALLC(IA,isizof('INT')*IBIG)                            
      DUM = COALLC(IA,isizof('REAL')*56)                                
C************************************                                   
C     CHECK CORE USAGE AND PRINT    *                                   
C************************************                                   
      ITOP1  = DUM + 56                                                 
      IF(IPRINT.GE.1) THEN                                              
        WRITE(6,*) '   '                                                
        WRITE(6,*) 'CORE ALLOCATION: ADDRESSES'                         
        WRITE(6,*) '   '                                                
        WRITE(6,*) 'IOFF   =',IOFF                                      
        WRITE(6,*) 'KADD   =',KADD                                      
        WRITE(6,*) 'LADD   =',LADD                                      
        WRITE(6,*) 'IJADD  =',IJADD                                     
        WRITE(6,*) 'IJGRP  =',IJGRP                                     
        WRITE(6,*) 'ORBTBF =',ORBTBF                                    
        WRITE(6,*) 'JOUT   =',JOUT                                      
        WRITE(6,*) 'ORBSYM =',ORBSYM                                    
        WRITE(6,*) 'IJOFF  =',IJOFF                                     
        WRITE(6,*) 'IJBKT  =',IJBKT                                     
        WRITE(6,*) 'X      =',X                                         
        WRITE(6,*) 'ILAST  =',ILAST                                     
        WRITE(6,*) 'ILAST2 =',ILAST2                                    
        WRITE(6,*) 'INTOFF =',INTOFF                                    
        WRITE(6,*) 'H      =',H                                         
        WRITE(6,*) 'F      =',F                                         
        WRITE(6,*) 'START  =',START                                     
        WRITE(6,*) 'BKTPNT =',BKTPNT                                    
        WRITE(6,*) 'BUF    =',BUF                                       
        WRITE(6,*) 'ISRTAR =',ISRTAR                                    
        WRITE(6,*) 'DUM    =',DUM                                       
        WRITE(6,*) '   '                                                
      ENDIF                                                             
      WRITE(IOUT,45) ITOP1,IBIG2                                        
   45 FORMAT(/,1X,'HIGHEST CORE ADDRESS USED IN ROUTINE GETIT IS ',     
     1        I8,/,1X,'WITH IBIG2 = ',I8,/)                             
C*********************************************************************  
C     MAKE ARRAY IJOFF, WHICH IS THE OFFSET OF THE IJ-TH ARRAY OF    *  
C     DIMENSION NIJ, AND IJBKT WHICH IS WHICH BUCKET THE IJ-TH       *  
C     ARRAY IS IN.                                                   *  
C*********************************************************************  
      CALL INITG(IA(IOFF),IA(IJOFF),IA(IJBKT))                          
C***********************************************************************
C   READ IN THE MO DENSITY MATRICES AND INTEGRALS.  STORE THE ONE      *
C   PARTICLE DENSITY MATRIX IN F AND THE ONE ELECTRON INTEGRALS IN H.  *
C   WRITE THE MO TWO PARTICLE DENSITY MATRIX IN IJ BLOCKS IN DRT ORDER *
C   TO IT91.  WRITE THE MO TWO ELECTRON INTEGRALS IN IJ BLOCKS IN DRT  *
C   ORDERING TO IT93.  THE OFFSETS ARE SEPARATED BY NBATRI, WHEREAS    *
C   THE DIMENSION IS STILL NTRI.                                       *
C***********************************************************************
      CALL GETIT(IA(IOFF),IA(IJADD),IA(IJGRP),IA(KADD),IA(LADD),        
     1 IA(X),IA(F),IA(H),IA(BUF),IA(ORBSYM),                            
     2 IA(START),IA(ILAST),IA(ILAST2),IA(INTOFF),IA(BKTPNT),            
     3 IA(ISRTAR),BKTSZ2,IA(ISRTAR),IA(DUM),ia(ijoff),ia(ijbkt),idfile)
C*********************************************************************  
C   WRITE THE MO LAGRANGIAN AND THE MO ONE PDM IN DRT ORDERING.      *  
C   FORM THE AO LAGRANGIAN AND ONE PDM.  WRITE ALL FOUR MATRICES     *  
C   TO MASTER FILE.                                                  *  
C*********************************************************************  
      CALL OUT1(IA(F),IA(H),ia47)                                       
C**********************************************                         
C     ALLOCATE ARRAYS FOR MAKEX               *                         
C**********************************************                         
      CALL COFREE(IA,START)                                             
      CALL COFREE(IA,BKTPNT)                                            
      CALL COFREE(IA,BUF)                                               
      CALL COFREE(IA,ISRTAR)                                            
      CALL COFREE(IA,DUM)                                               
      BUF2 = COALLC(IA,isizof('REAL')*BKTSZ)                            
      TPDM = COALLC(IA,isizof('REAL')*LENR)                             
      INTS = COALLC(IA,isizof('REAL')*LENR)                             
C************************************                                   
C     CHECK CORE USAGE AND PRINT    *                                   
C************************************                                   
      ITOP2  = INTS + LENR*INTOWP(1)                                    
      IF(IPRINT.GE.1) THEN                                              
        WRITE(6,*) '   '                                                
        WRITE(6,*) 'CORE ALLOCATION: ADDRESSES'                         
        WRITE(6,*) '   '                                                
        WRITE(6,*) 'BUF2   =',BUF2                                      
        WRITE(6,*) 'TPDM   =',TPDM                                      
        WRITE(6,*) 'INTS   =',INTS                                      
        WRITE(6,*) '   '                                                
      ENDIF                                                             
      WRITE(IOUT,46) ITOP2,NBKT,NGRPIJ,LGRPIJ,LENR                      
   46 FORMAT( 1X,'HIGHEST CORE ADDRESS USED IN ROUTINE MAKEX IS ',      
     1        I8,/,1X,'WITH NBKT = ',I4,', NGRPIJ = ',I4,               
     2        ', LGRPIJ = ',I4,' AND LENR = ',I8,/)                     
C*****************************************************                  
C     FORM THE LAGRANGIAN MATRIX IN DRT ORDERING.    *                  
C*****************************************************                  
      CALL MAKEX(IA(X),IA(TPDM),IA(INTS),IA(BUF2),IA(BUF2),             
     1 IA(ORBSYM),BKTSZ,IA(INTOFF),IA(ILAST),                           
     2 IA(ILAST2),ECINEW,BKTSZ2)                                        
C***********************************************************            
C     WE NO LONGER NEED WORK FILE 93 SO GET RID OF IT.     *            
C     GET RID OF MO TWOPDM UNLESS DESIRED AS KEPT.         *            
C***********************************************************            
      CALL RCLOSE(IT93,4)                                               
C****************************************************                   
C     ALLOCATE SOME CORE FOR EIGENVECTOR ARRAYS     *                   
C****************************************************                   
      CALL COFREE(IA,TPDM)                                              
      CALL COFREE(IA,INTS)                                              
      CALL COFREE(IA,BUF2)                                              
      WVEC = COALLC(IA,isizof('REAL')*NBFAO*NBF)                        
      EIV = COALLC(IA,isizof('REAL')*NBFAO*NBF)                         
      EIG = COALLC(IA,isizof('REAL')*NBF)                               
C*********************************************************************  
C     GET EIGENVECTOR IN PITZER ORDERING AND CONVERT TO DRT ORDERING *  
C     EIV IS TEMPORARY ARRAY, WVEC WILL HOLD DRT ORDERED AO-MO VECTOR*  
C*********************************************************************  
      CALL VECS(IA(JOUT),ACRCY,IA(WVEC),IA(EIG),IFLAG,IA(EIV),ia47)     
C***********************************                                    
C     ALLOCATE ARRAYS FOR OUTPUT   *                                    
C***********************************                                    
      CALL COFREE(IA,EIV)                                               
      CALL COFREE(IA,EIG)                                               
      F = COALLC(IA,isizof('REAL')*NBF*NBF)                             
      BUF3 = COALLC(IA,isizof('REAL')*NBFAO*NBF)                        
      XAO = COALLC(IA,isizof('REAL')*NIJAO)                             
C*********************************************************************  
C   WRITE THE MO LAGRANGIAN AND THE MO ONE PDM IN DRT ORDERING.      *  
C   FORM THE AO LAGRANGIAN AND ONE PDM.  WRITE ALL FOUR MATRICES     *  
C   TO MASTER FILE.                                                  *  
C*********************************************************************  
      CALL OUT2(IA(IOFF),IA(F),IA(X),IA(WVEC),IA(ORBTBF),IA(BUF3),      
     1          IA(XAO),ia47)                                           
C******************                                                     
C     THE END     *                                                     
C******************                                                     
      WRITE(IOUT,4)                                                     
   4  FORMAT(//,7X,' >>>>> LAGRANGIAN CALCULATION COMPLETED <<<<< '/)   
      if(iprint.ge.1) then
       CALL CODUMP(IA)                                                   
      endif
      call fioflu(6)
      JCODE = 3                                                         
      CALL RCLOSE(IT30,JCODE)                                           
      CALL RCLOSE(IT40,JCODE)                                           
      CALL RCLOSE(IT52,JCODE)                                           
      IF(.NOT.CCSD) CALL RCLOSE(IT54,JCODE)                             
      CALL RCLOSE(IT58,JCODE)                                           
      CALL RCLOSE(IT91,JCODE)                                           
C                                                                       
C      CALL TSTOP(ICHECK)                                                
      CALL TSTOP(IOUT)                                                  
      stop
      RETURN                                                            
      END                                                               
