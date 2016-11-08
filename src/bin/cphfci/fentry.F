      subroutine fentry(z,iz,maxcor)                               
C********************************************************************** 
C GENERAL COUPLED-PERTURBED HARTREE-FOCK PROGRAM FOR CI AND CC FIRST  * 
C********************************************************************** 
C                                                                     * 
C********************************************************************** 
C*   NOTICE OF PROGRAM MODIFICATION                                   * 
C********************************************************************** 
c       by: Daniel Crawford (tdc101195)
c     date: Oct. 11, 1995
c   reason: Multiple ffile() calls for file11 made things die with the
c           new IBM operating system, AIX4.1
c
c       by: Edward Seidl (cets072091)
c     date: July 20, 1991
c   reason: add option DELETE52 for debugging
c
c       by: Tracy P. Hamilton (ctph)
c     date: Summer 1991 (near as I can tell)
c   reason: 1) changed dimensions of various arrays from norbs*norbs to
c              nbfao*norbs to accomodate pure angular momentum functions
c           2) added common block dimao
c           3) replace read47 with read40
c           4) use new core allocation routines
c
c       by: Edward Seidl (cets030691)
c     date: March 6, 1991
c   reason: change dermat, drmat2, main, and b0mat to read file42
c           using sread so that it will be compatable with derivth
c           (which will soon become deriv?)
c           also write a little extra info to file11
c           and add print control
c
c       by: Edward Seidl (cets030291)
c     date: March 2, 1991
c   reason: add new input and replace drum with psinit
c
c       by: Edward Seidl (cets121890)
c     date: December 18, 1990
c   reason: 1) finish changes to gtints started by tph
c           2) removed all entry points
c           3) changed dimension of zeta from numij*norbs to
c                numij*(ntypes+1) and changed getint, ga* and gb*
c                subroutines
c                (ie throw away ALOT of zeros and save tons of memory)
c           4) added ability to do tcscf-ci (maybe, someone should
c                 check and make sure the theory is valid)
c           5) removed subroutines gtint3 and gtint4, so no more scf
c                second derivatives
c           6) fix bug in dcicof by defining a2 to be 2.0d0
c
C*  BY:  Edward Seidl (cets121890) as added by TPH                    * 
C*  DATE:  December 20, 1990                                          * 
C*  REASON: Added tcscf-ci gradients.                                 * 
C*          Change zeta matrix handling (all 0 for virtual orbitals!) * 
C*          TPH also removed DRMAT2 subroutine (Same as DERMAT!)      * 
C********************************************************************** 
C*  BY:  TRACY HAMILTON                                               * 
C*  DATE:  November 19, 1990                                          * 
C*  REASON: Use iterative method to solve CPHF equation for Z vector. * 
C********************************************************************** 
C*  BY:  TRACY HAMILTON                                               * 
C*  DATE:  November 9, 1990                                           * 
C*  REASON: Remove junk code, put in fentry.                          * 
C*          Use as little core as possible by putting zeta on disk.   * 
C********************************************************************** 
C*  BY:  TRACY HAMILTON                            Search: c07-02-90  * 
C*  DATE:  July 2, 1990                                               * 
C*  REASON: CHANGE SUBROUTINE CALLS TO GA1, GA2, GA3, GA4 and GA6     * 
C*          IN SUBROUTINE GTINTS TO BE THE SAME AS IN GTINT2          * 
C*          SO OPEN-SHELL SINGLET CI GRADIENTS WILL WORK              * 
C********************************************************************** 
C*  BY:  Tracy Hamilton                                               * 
C*  DATE:  April 12, 1990                                             * 
C*  REASON: Make change so that nthcphf is just like ngrcphf as far   * 
C*          as ha + ta is concerned.  This required changing der1th.  * 
C********************************************************************** 
C*  BY:  TRACY HAMILTON                            SEARCH:  CTPH      * 
C*  DATE:  DECEMBER 18, 1988                                          * 
C*  REASON: MAKE CHANGES TO RUN F,G FUNCTIONS  COMMON/GRAD/           * 
C*          COMMON / GRAD/  PUT IN                                    * 
C*          ELIMINATE FILE47 IS THE GRAIL                             * 
C*          DROL IN OLD SCFTFK IS SAME AS LORD IN YY MASTER?          * 
C********************************************************************** 
C*  BY:  RICHARD REMINGTON                     search:  cMachineSize  * 
C*  DATE:  MAY   22,   1988                                           * 
C*  REASON: DECREASE CORE TO RUN IN 7MB ON 9370                       * 
C********************************************************************** 
C      BY:  GES                                                       * 
C    DATE:  OCTOBER 12,1988                                           * 
C  REASON:  GET RID OF ARRAY W IN DYNAMICAL ALLOCATION                * 
C********************************************************************** 
C      BY:  GES                                                       * 
C    DATE:  MARCH 4, 1988                                             * 
C  REASON:  ADD DIPOLE MOMENTS TO FILE11                              * 
C********************************************************************** 
C      BY:  GES                                                       * 
C    DATE:  JANUARY 7, 1988                                           * 
C  REASON:  GET RID OF BLANK LINE IN FILE11                           * 
C********************************************************************** 
C      BY:  ACS                                                       * 
C    DATE:  JANUARY 2, 1988                                           * 
C  REASON:  WRITE CC ENERGY TO FILE 11                                * 
C********************************************************************** 
C   GES VERSION: OCTOBER 6, 1987                                      *\
C   INCLUDES ONLY CLOSED SHELL AND HIGH-SPIN OPEN SHELL CASES.        * 
C********************************************************************** 
C   GES VERSION: JUNE 26, 1987                                        * 
C   INVERTING ONLY A1 -A1 COMPONENT OF A MATRIX!                      * 
C********************************************************************** 
C   GES VERSION: MAY 6, 1987                                          * 
C   CC SYMM GRADS AND MONGO                                           * 
C********************************************************************** 
C      BY:  TJL                                                       * 
C    DATE:  NOVEMBER 1, 1986                                          * 
C  REASON:  CORRECT BUG IN DEPENDENT PAIR CONTRIBUTION TO THE CI      * 
C  REASON:  DIPOLE MOMENT, THIS IS A MINOR CORRECTION AND SHOULD NOT  * 
C  REASON:  AFFECT MOST CASES,  HIGHER EXCITATIONS MAY BE AFFECTED    * 
C  REASON:  (I.E., WHEN THE CPHF CONTRIBUTION SHOULD GO TO ZERO.)     * 
C********************************************************************** 
C      BY:  TJL                                                       * 
C    DATE:  AUGUST 3, 1985                                            * 
C  REASON:  WRITE CI ENERGY TO FILE 11                                * 
C********************************************************************** 
C      BY:  TJL                                                       * 
C    DATE:  APRIL 14, 1985                                            * 
C  REASON:  IMPLEMENT COR'S AND VIR'S FOR OPEN SHELL                  * 
C  REASON:  OPEN SHELL SINGLETS STILL DO NOT WORK; WITH OR WITHOUT C/V* 
C  REASON:  COR'S AND VIR'S FOR DIPOLE MOMENTS HAVE ALSO BEEN IMPL.   * 
C  REASON:  THESE HAVE NOT BEEN TESTED.                               * 
C     UNIT 3 TO CHECK                                                 * 
C  REASON:  CHANGE UNIT NUMBERS !! 2=>52, 82=>56, 21=>55,  JTAP1=>94  * 
C  REASON:  CHANGE UNIT NUMBERS !! 8=>58, 20=>53, 22=>42,  JTAP2=>96  * 
C  REASON:  CHANGE UNIT NUMBERS !!12=>54, 57=>47         , IUTAP=>98  * 
C********************************************************************** 
C      BY:  RICHARD REMINGTON                                         * 
C    DATE:  FEBRUARY 22, 1985                                         * 
C  REASON:  ADD SUBROUTINE PRN11 IN ORDER TO WRITE CI GRADIENTS INTO  * 
C           FILES 11, 30, AND THE OUTPUT (FILE6).                     * 
C********************************************************************** 
      REAL*8 Z(MAXCOR)                                                  
      INTEGER IZ(1)                                                     
      COMMON /TAPES/ NOUT,ICHK,IT30,IT42,IT47,IT52,IT58,IT94,IT96,IT98  
      COMMON /FILES/ IT54                                               
C                                                                       
 1000 FORMAT(///,10X,' >>>>> GENERAL CPHF CALCULATION COMPLETED <<<<< ')
C                                                                       
      ICHK = 3                                                          
      NOUT = 6                                                          
      IT30 = 30                                                         
      IT42 = 42                                                         
      IT47 = 47                                                         
      IT52 = 52                                                         
      IT54 = 54                                                         
      IT58 = 58                                                         
      ITAP86=86                                                         
      ITAP87=87                                                         
      ITAP88=88                                                         
      IT94 = 94                                                         
      IT96 = 96                                                         
C                                                                       
      call psinit('APPEND')                                             
      CALL TSTART(6)                                                    
C                                                                       
      CALL MAININ(IZ,Z,MAXCOR)                                            
C                                                                       
      CALL RCLOSE(IT94,4)                                               
      CALL RCLOSE(IT96,4)                                               
C                                                                       
      WRITE(NOUT,1000)                                                  
C                                                                       
      CALL TSTOP(6)                                                     
C                                                                       
      STOP                                                              
      END                                                               
