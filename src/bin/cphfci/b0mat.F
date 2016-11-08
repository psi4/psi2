      SUBROUTINE B0MAT(ALPHA,BETA,EPS,EPA,SS,ZETA,CC,OCC,               
     1                 INTS,IJORD,IOFF,KADD,LADD,IJGRP,IJADD,ORBSYM,    
     2                 JT1C,LG,BOX,C,EX,RS,IJDEP,EAVF,                  
     3                 LSTEPA,EX2,CCSD,DOFF,DADD,NMA1,IDIMCC,DB,ICASE,  
     4                 ha11,ha22,ha12,e11,e22,e12)
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      INTEGER DOFF(NUMIJ),DADD(NUMIJ)                                   
      REAL*8 det,ha11(ndf),ha12(ndf),ha22(ndf)                         
      REAL*8 INTS(NMAX),LG(NORBS,NORBS),e11(1),e12(1),e22(1)            
      REAL*8 CDIP(3),A2,DEBYE,FAC,A1,A1S2,EAVF(NUMIJ),LSTEPA(NUMIJ)     
      REAL*8 EX2(IJIND),CIDT(3),a4                                      
      LOGICAL IVIR,JVIR,CCSD                                            
      INTEGER FDOC,CORDOC,IJDEP(NUMIJ),FVIR,VUOC                        
      DIMENSION ALPHA(NUMIJ),BETA(NUMIJ),EPS(NUMIJ),EPA(NORBS,NORBS),   
     1 SS(NUMIJ,NDF),DB(*),                                             
cets 2 ZETA(NUMIJ,*),CC(idimcc),OCC(NORBS),                                  
     2 ZETA(NTYPES+1,NUMIJ),CC(idimcc),OCC(NORBS),            
     3 BOX(NDF),C(IJIND,NDF),XCIJA(28),EX(IJIND),RS(NUMDEP,NDF)         
      INTEGER IJORD(NUMIJ),IOFF(NUMIJ),KADD(1),LADD(1),IJGRP(1),        
     1 IJADD(1),ORBSYM(NORBS),                                          
     2 JT1C(1),IAJT(20)                                                 
      COMMON /TAPES/ NOUT,ICHK,IT30,IT42,IT47,IT52,IT58,IT94,IT96,IT98  
      COMMON /DIMAO/ NBFAO,NBATRI                                       
      COMMON /DIMSG/ ISYMRB,NUMIJ,NBF,NMAX,NBFORB,NSYM,NDF,NORBS        
      COMMON /POINT/ FOCC(5),IPOINT(5),JPOINT(5),NSORB(5)               
      COMMON /FUNCS/ IOCC,JOCC,KOCC,NTYPES,NBSET,NAT                    
      COMMON /IEXT/ IPRINT,IGRP,ICIDIP                                  
      COMMON /OCCS/ NSPE,NSA,NSB,NALP,NOP,NDOC,LFERM,NUOC,IJIND,IJIND3  
      COMMON/OCCS2/NUMCOR,NCOR,FDOC,CORDOC,NUMVIR,FVIR,LUOC,NUMDEP,VUOC 
      COMMON /LOC42/ IA42(100),IBF2AT(284)                              
      COMMON /CALCJ/ ICALT,IDERT,ICIT                                   
      COMMON /GRAD/ E1A(150)                                            
      common /tcscf/ derc1,derc2,ci1,ci2,h11,h12,h22,one,itcscf
C                                                                       
      WRITE(ICHK,2222)                                                  
 2222 FORMAT(/'  IN B0MAT')                                             
      NORBSQ=NORBS*NORBS                                                
      A1S2=0.50D+00                                                     
      A4 = 4.0D+00                                                      
      A2 = 2.0D+00                                                      
      A1 = 1.0D+00                                                      
      DEBYE = 2.541765480D+00                                           
      IT98=98                                                           
C                                                                       
C  READ BACK IN S DERIVATIVE INTEGRALS                                  
C                                                                       
      IS=1                                                              
      DO 101 IDF=1,NDF                                                  
        CALL WREADW(IT94,SS(1,IDF),INTOWP(NUMIJ),IS,IS)                 
        IF (IAND(IPRINT,64) .NE. 0) THEN                                
          WRITE(NOUT,1101) IDF                                          
 1101     FORMAT(/'  THE SA MATRIX FOR DEGREE ',I5,'  AS READ IN B0MAT')
          CALL PRINT(SS(1,IDF),NUMIJ,NORBS,NOUT)                        
        END IF                                                          
 101  CONTINUE                                                          
C                                                                       
C  INITIALIZE ROUTINE FOR CALCULATION OF A MATRIX.                      
C                                                                       
      IDIMCC=NMA1*(NMA1+1)                                              
c     IF(ICASE.EQ.1) THEN                                               
c      CALL GETA(CC,EPS,ZETA,ALPHA,BETA,DOFF,DADD,NUMIJ,IDIMCC)       
c      CALL GETB(C ,EPS,ZETA,ALPHA,BETA,SS)                           
c      CALL GETC(RS,EX,EAVF,LG,SS)                                    
c     ELSE IF(ICASE.EQ.2) THEN                                          
c      CALL GETA2(DB,EPS,ZETA,ALPHA,BETA,C)                             
c      CALL GETB2(C ,EPS,ZETA,ALPHA,BETA,SS)                            
c      CALL GETC2(RS,EX,EAVF,LG,SS)                                     
c     ELSE IF(ICASE.EQ.3) THEN                                          
c      CALL GETA3(CC,EPS,ZETA,ALPHA,BETA,DOFF,DADD,NUMIJ,IDIMCC)       
c      CALL GETB3(C ,EPS,ZETA,ALPHA,BETA,SS)                           
c      CALL GETC3(RS,EX,EAVF,LG,SS)                                    
c     ENDIF                                                             
C-------------------------------------------                            
C     READ LAGRANGIAN  AND ZERO OUT BOX    |                            
C-------------------------------------------                            
      CALL ZERO(BOX,NDF)                                                
      CALL ZERO(RS,NUMDEP*NDF)                                          
      CALL ZERO(EX,IJIND)                                               
      CALL ZERO(EX2,IJIND)                                               
      CALL ZERO(C,IJIND*NDF)                                           
      CALL MREAD(LG,33)                                                 
c     CALL MATOUT(LG,NORBS,NORBS,NORBS,NORBS,NOUT)                      
C-----------------------------------------                              
C   MAKE B0 MATRICES, EX, EX2 AND RS.    |                              
C   ALSO A MATRIX FOR CASES 1 AND 3.     |                              
C-----------------------------------------                              
      IF(ICASE.EQ.1) THEN                                               
       CALL GTINTS(C,EPS,ZETA,OCC,ALPHA,BETA,IJORD,KADD,LADD,IJADD,              
     1        IJGRP,ORBSYM,INTS,SS,EPA,JT1C,IOFF,EX,RS,LG,IJDEP,BOX,  
     2        eavf,LSTEPA,cc,doff,dadd,numij,idimcc,e11,e12,e22)
C                                                                       
       CALL GTINT2(C,EPS,ZETA,OCC,ALPHA,BETA,IJORD,KADD,LADD,IJADD,               
     1        IJGRP,ORBSYM,INTS,SS,EPA,JT1C,IOFF,EX,RS,LG,IJDEP,BOX,  
     2        eavf,LSTEPA,EX2,cc,doff,dadd,numij,idimcc,e11,e12,e22,
     3        ha11,ha22,ha12,nma1)
      ENDIF                                                             
c     IF(ICASE.EQ.2) THEN                                               
c      CALL GTINTS2(ZETA,ALPHA,BETA,IJORD,KADD,LADD,IJADD,              
c    1        IJGRP,ORBSYM,INTS,SS,EPA,JT1C,IOFF,RS,LG,IJDEP,BOX,EAVF,  
c    2        LSTEPA,EPS,EX)                                            
C                                                                       
c      CALL GTINT2(C,ZETA,ALPHA,BETA,IJORD,KADD,LADD,IJADD,               
c    1        IJGRP,ORBSYM,INTS,SS,EPA,JT1C,IOFF,RS,LG,IJDEP,BOX,EAVF,  
c    2        LSTEPA,EX2,EPS,EX)                                        
c     ENDIF                                                             
c     WRITE(*,*)'B0 MATRIX AFTER GTINT2'                                
c     CALL MATOUT(C,IJIND,NDF,IJIND,NDF,NOUT)                           
C--------------------------                                             
C     STORE B0 ON IT96    |                                             
C--------------------------                                             
      IC=1                                                              
      CALL WWRITW(IT96,C,INTOWP(IJIND*NDF),IC,JUNK)                     
c     WRITE(ICHK,9011)                                                  
9011  FORMAT(' C WRITTEN TO IT96 ')                                     
C------------------------------------------                             
C     MAKE PART OF BOX, EPSILON ARRAY     |                             
C------------------------------------------                             
      CALL ZERO(C,IJIND)                                                
C                                                                       
      DO 400 IXDF=1,NDF                                                 
c     WRITE(*,*) 'BOX,IXDF',BOX(IXDF),IXDF                              
      DO 410 I=1,NORBS                                                  
      DO 420 R=1,I-1                                                    
      IR=IOFF(I)+R                                                      
      BOX(IXDF)=BOX(IXDF)-SS(IR,IXDF)*(LG(I,R)+LG(R,I))                 
 420  CONTINUE                                                          
      II=IOFF(I)+I                                                      
      BOX(IXDF)=BOX(IXDF)-SS(II,IXDF)*LG(I,I)                           
 410  CONTINUE                                                          
c     WRITE(*,*) 'BOX,IXDF',BOX(IXDF),IXDF                              
 400  CONTINUE                                                          
      DO 120 I=LFERM,NDOC                                               
      IO=IOFF(I)                                                        
      DO 130 J=1,NUOC                                                   
      FAC = OCC(J)-OCC(I)                                               
      IJ=IO+J                                                           
      IJND=IJORD(IJ)                                                    
      C(IJND,1)=((LG(J,I)-LG(I,J)) + FAC*EX(IJND) + EX2(IJND))          
  130 CONTINUE                                                          
  120 CONTINUE                                                          
      IF (NOP.EQ.0) GOTO 200                                            
      DO 140 I=NOP,NALP                                                 
      IO=IOFF(I)                                                        
      DO 150 J=1,NDOC                                                   
      FAC = OCC(J)-OCC(I)                                               
      IJ=IO+J                                                           
      IJND=IJORD(IJ)                                                    
      C(IJND,1)=((LG(J,I)-LG(I,J)) + FAC*EX(IJND) + EX2(IJND))          
  150 CONTINUE                                                          
  140 CONTINUE                                                          
  200 CONTINUE                                                          
      IF (NSPE.EQ.0) GOTO 110                                           
      DO 160 I=NSA,NSB                                                  
      IO=IOFF(I)                                                        
      DO 170 J=1,I-1                                                    
      FAC = OCC(J)-OCC(I)                                               
      IJ=IO+J                                                           
      IJND=IJORD(IJ)                                                    
      C(IJND,1)=((LG(J,I)-LG(I,J)) + FAC*EX(IJND) + EX2(IJND))          
  170 CONTINUE                                                          
  160 CONTINUE                                                          
  110 CONTINUE                                                          
      IF(ICASE.EQ.1.OR.ICASE.EQ.3) THEN                                 
C--------------------------------------------------------               
C     SORT PERTURBATION TO PACKED-SYMMETRY ORDERING     |               
C--------------------------------------------------------               
       IADR=1                                                           
       IADR=IADR+NMA1*NMA1                                              
       DO 918 I=LFERM,NDOC                                              
        ISYM=ORBSYM(I)-1                                                
        DO 916 J=1,NUOC                                                 
         JSYM=ORBSYM(J)-1                                               
         IJSYM=IEOR(ISYM,JSYM)+1                                        
         IF(IJSYM.NE.1)GO TO 915                                        
         IJ=IOFF(I)+J                                                   
         IJND=IJORD(IJ)                                                 
C        WRITE(*,*)'STORING',C(IJND,1),'IN POSITION',IADR               
         CC(IADR)=C(IJND,1)                                             
         IADR=IADR+1                                                    
  915    CONTINUE                                                       
  916   CONTINUE                                                        
  918  CONTINUE                                                         
       IF(NOP.EQ.0)GO TO 929                                            
       DO 928 I=NOP,NALP                                                
        ISYM=ORBSYM(I)-1                                                
        DO 926 J=1,NDOC                                                 
         JSYM=ORBSYM(J)-1                                               
         IJSYM=IEOR(ISYM,JSYM)+1                                        
         IF(IJSYM.NE.1)GO TO 925                                        
         IJ=IOFF(I)+J                                                   
         IJND=IJORD(IJ)                                                 
C        WRITE(*,*)'STORING',C(IJND,1),'IN POSITION',IADR               
         CC(IADR)=C(IJND,1)                                             
         IADR=IADR+1                                                    
  925    CONTINUE                                                       
  926   CONTINUE                                                        
  928  CONTINUE                                                         
  929 CONTINUE                                                          
C                                                                       
       IF(NSPE.EQ.0)GO TO 939                                           
       DO 938 I=NSA,NSB                                                 
        ISYM=ORBSYM(I)-1                                                
        DO 936 J=1,I-1                                                  
         JSYM=ORBSYM(J)-1                                               
         IJSYM=IEOR(ISYM,JSYM)+1                                        
         IF(IJSYM.NE.1)GO TO 935                                        
         IJ=IOFF(I)+J                                                   
         IJND=IJORD(IJ)                                                 
C        WRITE(*,*)'STORING',C(IJND,1),'IN POSITION',IADR               
         CC(IADR)=C(IJND,1)                                             
         IADR=IADR+1                                                    
  935    CONTINUE                                                       
  936   CONTINUE                                                        
  938  CONTINUE                                                         
  939  CONTINUE                                                         
cets
       if(itcscf.eq.1) then
        xxwa1 = derc1
        xxwa2 = derc2
        derc1 = -a4*ci1*xxwa1
        derc2 = -a4*ci2*xxwa2
        cc(idimcc-1) = derc1
        cc(idimcc) = derc2
       endif
C---------------------------------                                      
C     INVERT ONLY IRREP A1*A1    |                                      
C---------------------------------                                      
c     WRITE(*,*)'A MATRIX BEFORE FLINX'                                
c     CALL MATOUT(CC,nma1,nma1+1,nma1,nma1+1,NOUT)                           
      CALL FLINX(CC,NMA1,NMA1,1,DET)                                    
C--------------------------                                             
C     UNPACK SOLUTION     |                                             
C--------------------------                                             
      CALL ZERO(C,IJIND)                                                
      IADR=1                                                            
      IADR=IADR+NMA1*NMA1                                               
      DO 948 I=LFERM,NDOC                                               
      ISYM=ORBSYM(I)-1                                                  
      DO 946 J=1,NUOC                                                   
      JSYM=ORBSYM(J)-1                                                  
      IJSYM=IEOR(ISYM,JSYM)+1                                           
      IF(IJSYM.NE.1)GO TO 945                                           
      IJ=IOFF(I)+J                                                      
      IJND=IJORD(IJ)                                                    
      C(IJND,1)=CC(IADR)                                                
      IADR=IADR+1                                                       
  945 CONTINUE                                                          
  946 CONTINUE                                                          
  948 CONTINUE                                                          
C                                                                       
      IF(NOP.EQ.0)GO TO 959                                             
      DO 958 I=NOP,NALP                                                 
      ISYM=ORBSYM(I)-1                                                  
      DO 956 J=1,NDOC                                                   
      JSYM=ORBSYM(J)-1                                                  
      IJSYM=IEOR(ISYM,JSYM)+1                                           
      IF(IJSYM.NE.1)GO TO 955                                           
      IJ=IOFF(I)+J                                                      
      IJND=IJORD(IJ)                                                    
      C(IJND,1)=CC(IADR)                                                
      IADR=IADR+1                                                       
  955 CONTINUE                                                          
  956 CONTINUE                                                          
  958 CONTINUE                                                          
  959 CONTINUE                                                          
C                                                                       
      IF(NSPE.EQ.0) GO TO 969                                           
      DO 968 I=NSA,NSB                                                  
      ISYM=ORBSYM(I)-1                                                  
      DO 966 J=1,I-1                                                    
      JSYM=ORBSYM(J)-1                                                  
      IJSYM=IEOR(ISYM,JSYM)+1                                           
      IF(IJSYM.NE.1)GO TO 965                                           
      IJ=IOFF(I)+J                                                      
      IJND=IJORD(IJ)                                                    
      C(IJND,1)=CC(IADR)                                                
      IADR=IADR+1                                                       
  965 CONTINUE                                                          
  966 CONTINUE                                                          
  968 CONTINUE                                                          
  969 CONTINUE                                                          
cets
      if(itcscf.eq.1) then
       c(ijind-1,1) = cc(idimcc-1)
       c(ijind,1) = cc(idimcc)
      endif
C                                                                       
      DO 331 I=1,IJIND                                                  
       CC(I)=C(I,1)                                                     
  331 CONTINUE                                                          
c     ELSE IF(ICASE.EQ.2) THEN                                          
C---------------------------------                                      
C     SOLVE CPHF ITERATIVELY     |                                      
C---------------------------------                                      
c      ITER = 1                                                         
c      CALL GTINTS4(ZETA,ALPHA,BETA,IJORD,KADD,LADD,IJADD,              
c    1        IJGRP,ORBSYM,INTS,SS,EPA,JT1C,IOFF,RS,LG,IJDEP,BOX,EAVF,  
c    2        LSTEPA,EPS,EX,CCSD)                                       
C                                                                       
c      CALL GTINT4(ZETA,ALPHA,BETA,IJORD,KADD,LADD,IJADD,               
c    1        IJGRP,ORBSYM,INTS,SS,EPA,JT1C,IOFF,RS,LG,IJDEP,BOX,EAVF,  
c    2        LSTEPA,EX2,EPS,EX)                                        
      ENDIF                                                             
C--------------------                                                   
C     READ IN B0    |                                                   
C--------------------                                                   
c     CALL MATOUT(CC,IJIND,1,IJIND,1,NOUT)                              
      CALL WREADW(IT96,C,INTOWP(IJIND*NDF),IC,JUNK)                     
c     WRITE(NOUT,9021)                                                  
9021  FORMAT(/' B0 MATRIX FROM IT96')                                   
c     CALL MATOUT(C,IJIND,NDF,IJIND,NDF,NOUT)                           
C--------------------                                                   
C     FINISH BOX    |                                                   
C--------------------                                                   
      DO 332 IIJ=1,IJIND                                                
      DO 333 IDF=1,NDF                                                  
      BOX(IDF)=BOX(IDF)+C(IIJ,IDF)*CC(IIJ)                              
  333 CONTINUE                                                          
  332 CONTINUE                                                          
      DO 334 IDF=1,NDF                                                  
      BOX(IDF)=BOX(IDF)*A1S2                                            
      WRITE(ICHK,9023) BOX(IDF),IDF                                     
9023  FORMAT(/' BOX = ',F20.8,' FOR DEGREE OF FREEDOM ',I8)             
  334 CONTINUE                                                          
C                                                                       
      DO 335 IDF=1,NDF                                                  
       BOX(IDF)=BOX(IDF)+E1A(IDF)                                       
       WRITE(ICHK,3000) IDF,BOX(IDF)                                    
3000   FORMAT(/' GRADIENT FOR ',I6,' DEGREE OF FREEDOM = ',F20.10)      
 335  CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
