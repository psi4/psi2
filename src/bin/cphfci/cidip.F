      SUBROUTINE CIDIP(DIPINT,VEC,SCR1,SCR2,Z,IOFF,IJORD,IJDEP,EAVF,LG, 
c     1 OCC,ZZ,IDIMZZ,CIDT)                                              
     1 occ,zz,cidt,orbtbf)
      IMPLICIT INTEGER (A-Z)                                            
C                                                                       
c      REAL*8 ZZ(IDIMZZ)                                                 
      real*8 zz(numij)
      REAL*8 DIPINT(NBATRI),VEC(NBFAO,NORBS),SCR1(NBFAO,NBFAO)          
      REAL*8 SCR2(NBFAO,NBFAO),Z(IJIND),CIDM(3),CIDC(3),CIDT(3)         
      REAL*8 EAVF(NUMIJ),LG(NBF,NBF),OCC(NBF)                           
      REAL*8 DEBYE,FAC,A0,A1,A1S2,CITM,FAC2,NORM,a4,kappa1,kappa2       
      real*8 ci1,ci2,kappa(3),a2,derc1,derc2,h11,h22,h12,one
c      INTEGER IOFF(NUMIJ),IJORD(NUMIJ),IJDEP(NUMIJ)                     
      integer ioff(numij),ijord(numij),ijdep(numij),orbtbf(norbs)
C                                                                       
      COMMON /TAPES/ NOUT,ICHK,IT30,IT42,IT47,IT52,IT58,IT94,IT96,IT98  
      common /tcscf/ derc1,derc2,ci1,ci2,h11,h12,h22,one,itcscf
      COMMON /DIMAO/ NBFAO,NBATRI                                       
      COMMON /DIMSG/ ISYMRB,NUMIJ,NBF,NMAX,NBFORB,NSYM,NDF,NORBS        
      COMMON /FUNCS/ IOCC,JOCC,KOCC,NTYPES,NBSET,NAT                    
      COMMON /IEXT/ IPRINT,IGRP,ICIDIP                                  
      COMMON /OCCS/ NSPE,NSA,NSB,NALP,NOP,NDOC,LFERM,NUOC,IJIND,IJIND3  
      COMMON/OCCS2/NUMCOR,NCOR,FDOC,CORDOC,NUMVIR,FVIR,LUOC,NUMDEP,VUOC 
      COMMON/DEP/IJDC,IFZC,IFZV,IZDIM                                   
C                                                                       
      DATA DEBYE / 2.541765480D+00 /                                    
      DATA A0,A1,A1S2 / 0.0D+00,1.0D+00,0.5D+00 /                       
C                                                                       
    1 FORMAT(5X,'THE CI DIPOLE MOMENT ',10X,3F20.10)                    
    2 FORMAT(5X,'CORRECTION TO CI D.M.',10X,3F20.10)                    
    3 FORMAT(5X,'THE TOTAL CI DIPOLE MOMENT',5X,3F20.10)                
    4 FORMAT(5X,'MAGNITUDE OF THE CI D.M.',7X,F20.10,5X,'DEBYES')       
    5 FORMAT(5X,'KAPPA CONTRIBUTION TO D.M.',5X,3F20.10,5X,'DEBYES')       
C                                                                       
      a2=2.0d+00
      a4=4.0d+00
cedtool020998
      call zero(zz,numij)
      DIPTAP = 59                                                       
      CALL RFILE(DIPTAP)                                                
C                                                                       
      CALL ZERO(CIDC,3)                                                 
C*******************************                                        
C   LOOP OVER X,Y,Z DIPOLES.   *                                        
C*******************************                                        
      WPOINT = 1                                                        
      DO 110 IDIP = 1,3                                                 
C*********************************************                          
C   TRANSFORM DIPOLE INTEGRALS TO MO BASIS   *                          
C*********************************************                          
       CALL WREADW(DIPTAP,DIPINT,INTOWP(NBATRI),WPOINT,WPOINT)          
       CALL SQUR(DIPINT,SCR1,NBFAO)                                     
       CALL ZERO(SCR2,NBFAO*NBFAO)                                      
       CALL MXMB(VEC,NBFAO,1,SCR1,1,NBFAO,SCR2,1,NBFAO,NBF,NBFAO,NBFAO) 
       CALL ZERO(SCR1,NBFAO*NBFAO)                                      
       CALL MXMB(SCR2,1,NBFAO,VEC,1,NBFAO,SCR1,1,NBFAO,NBF,NBFAO,NBF)   
C*******************************                                        
C   DOUBLY OCCUPIED ORBITALS   *                                        
C*******************************                                        
       DO 120 I=LFERM,NDOC                                              
        IO=IOFF(I)
        in = orbtbf(i)                                                      
        DO 130 J=1,NUOC                                                 
         FAC = OCC(I) - OCC(J)                                          
         IJ=IO+J
         jn = orbtbf(j)
         ii = max0(in,jn)
         jj = min0(in,jn)
         ijn = ioff(ii) + jj
         IJND=IJORD(IJ)                                                 
         CIDC(IDIP) = CIDC(IDIP) + SCR1(I,J)*Z(IJND)*FAC                
         zz(ijn) = z(ijnd)*fac
  130   CONTINUE                                                        
  120  CONTINUE                                                         
       IF (NOP.EQ.0) GOTO 200                                           
C**************************                                             
C   OPEN SHELL ORBITALS   *                                             
C**************************                                             
       DO 140 I=NOP,NALP                                                
        IO=IOFF(I)                                                      
        in = orbtbf(i)
        DO 150 J=1,NDOC                                                 
         FAC = OCC(I) - OCC(J)                                          
         IJ=IO+J                                                        
         jn = orbtbf(j)
         ii = max0(in,jn)
         jj = min0(in,jn)
         ijn = ioff(ii) + jj
         IJND=IJORD(IJ)                                                 
         CIDC(IDIP) = CIDC(IDIP) + SCR1(I,J)*Z(IJND)*FAC                
         zz(ijn) = z(ijnd)*fac
  150   CONTINUE                                                        
  140  CONTINUE                                                         
  200  CONTINUE                                                         
       IF (NSPE.EQ.0) GOTO 520                                          
C***********************                                                
C   SPECIAL ORBITALS   *                                                
C***********************                                                
       DO 160 I=NSA,NSB                                                 
        IO=IOFF(I)                                                      
        in = orbtbf(i)
        DO 170 J=1,I-1                                                  
         FAC = OCC(I) - OCC(J)                                          
         IJ=IO+J                                                        
         jn = orbtbf(j)
         ii = max0(in,jn)
         jj = min0(in,jn)
         ijn = ioff(ii) + jj
         IJND=IJORD(IJ)                                                 
         CIDC(IDIP) = CIDC(IDIP) + SCR1(I,J)*Z(IJND)*FAC                
         zz(ijn) = z(ijnd)*fac
  170   CONTINUE                                                        
  160  CONTINUE                                                         
cets121890
       if(itcscf.eq.1) then
        kappa1 = scr1(nsb,nsb)-scr1(nsa,nsa)
        kappa2 = ( z(ijind-1)*ci2 - z(ijind)*ci1 ) *ci1*ci2
        kappa(idip) = -a2*kappa1*kappa2
cedtool020998
        in = orbtbf(nsa)
        jn = orbtbf(nsb)
        ii = ioff(in) + in
        jj = ioff(jn) + jn
        zz(ii) = a2*kappa2
        zz(jj) = -a2*kappa2
cedtool020998 end of changes
        cidc(idip) = cidc(idip) + kappa(idip)
        kappa(idip) = kappa(idip)*debye
       endif
  520  CONTINUE                                                         
C************************************                                   
C   CONTRIBUTION FROM LAGRANGIAN    *                                   
C************************************                                   
       IJ = 0                                                           
       DO 500 I = 1,NBF                                                 
        II = IOFF(I) + I                                                
        in = orbtbf(i)
        DO 510 J = 1,I                                                  
         JJ = IOFF(J) + J                                               
         jn = orbtbf(j)
         k = max0(in,jn)
         l = min0(in,jn)
         ijn = ioff(k) + l
         IJ = IJ + 1                                                    
         IJCV = IJDEP(IJ)                                               
         IF(IJCV.EQ.0) GO TO 510                                        
         FAC = A1/(EAVF(II)-EAVF(JJ))                                   
         CIDC(IDIP) = CIDC(IDIP) + FAC*(LG(J,I)-LG(I,J))*SCR1(I,J)      
         zz(ijn) = fac*(lg(j,i)-lg(i,j))
  510   CONTINUE                                                        
  500  CONTINUE                                                         
  110 CONTINUE                                                          
C**********************************************************             
C   END LOOPS OVER X,Y,Z.   CONVERT TO DEBYE AND PRINT.   *             
C**********************************************************             
      CALL WREADW(DIPTAP,CIDM,INTOWP(3),WPOINT,WPOINT)                  
      DO 210 I=1,3                                                      
       CIDC(I) = CIDC(I)*DEBYE                                          
       CIDT(I) = CIDC(I) + CIDM(I)                                      
  210 CONTINUE                                                          
      CITM = CIDT(1)*CIDT(1) + CIDT(2)*CIDT(2) + CIDT(3)*CIDT(3)        
      CITM = DSQRT(CITM)                                                
      WRITE(*,*)                                                        
      WRITE(*,*)                                                        
      WRITE(*,1) (CIDM(I),I=1,3)                                        
      WRITE(*,2) (CIDC(I),I=1,3)                                        
      WRITE(*,3) (CIDT(I),I=1,3)                                        
      WRITE(*,*)                                                        
      WRITE(*,*)                                                        
      WRITE(*,4) CITM                                                   
      WRITE(*,*)                                                        
      if(itcscf.eq.1) then
       WRITE(*,*)                                                        
       WRITE(*,5) (KAPPA(I),I=1,3)                                        
       WRITE(*,*)                                                        
      endif
C                                                                       
c      IDIMZ=IJIND+IJDC                                                  
c      DO 697 I=1,IJIND                                                  
c       ZZ(I)=Z(I)                                                       
c  697 CONTINUE                                                          
c      IJ = 0                                                            
c      IK = 0                                                            
c      DO 700 I = 1,NBF                                                  
c      II = IOFF(I) + I                                                  
c      DO 710 J = 1,I                                                    
c      JJ = IOFF(J) + J                                                  
c      IJ = IJ + 1                                                       
c      IJCV = IJDEP(IJ)                                                  
c      IF(IJCV.EQ.0) GO TO 710                                           
c      IK = IK + 1                                                       
c      FAC=1.0D0/(EAVF(II)-EAVF(JJ))                                     
c      ZZ(IJIND+IK) =(LG(J,I)-LG(I,J))                                   
c  710 CONTINUE                                                          
c  700 CONTINUE                                                          
C                                                                       
c      ND=NDOC-NUOC                                                      
c      NY=NUOC                                                           
c      if(iprint.gt.1) then                                              
c       WRITE(*,*)'          '                                           
c       WRITE(*,*)'Z VECTOR NO*NV PART'                                  
c       CALL MATOUT(ZZ,NY,ND,NY,ND,6)                                    
c       IF(NOP.NE.0)THEN                                                 
c        WRITE(*,*)'Z VECTOR OPEN SHELLS'                                
c        WRITE(6,6123)(ZZ(ND*NY+I),I=1,ND+NY)                            
c       END IF                                                           
c      IF(IJDC.NE.0)THEN                                                 
c      WRITE(*,*)'Z VECTOR COR-DOC AND VIR-UOC ROTATIONS'                
c      WRITE(6,6123)(ZZ(IJIND+I),I=1,IJDC)                               
c 6123 FORMAT(7X,10F12.7)                                                
c      END IF                                                            
c      endif
C                                                                       
c      CALL LOCATE (5,'# MONGO ###',IERR)                                
c      IF(IERR.NE.0) THEN                                                
c       WRITE (6,6007)                                                   
c 6007  FORMAT(//,2X,'NO INPUT TO MONGO HAS BEEN FOUND.',/)              
c       FLAG=0                                                           
c      ELSE                                                              
C**********************                                                 
C   DO MONGO STUFF    *                                                 
C**********************                                                 
c       READ(5,6004)FLAG                                                 
c 6004  FORMAT(I5)                                                       
       ITAP86=86                                                        
c       ITAP87=87                                                        
c       ITAP88=88                                                        
       CALL RFILE (ITAP86)                                              
c       CALL WWRITW(ITAP86,IDIMZ,1,1     ,IVEC86)                        
c       CALL WWRITW(ITAP86,NDOC-NUOC,1,IVEC86,IVEC86)                    
c       CALL WWRITW(ITAP86,NUOC,1,IVEC86,IVEC86)                         
c       CALL WWRITW(ITAP86,IFZC,1,IVEC86,IVEC86)                         
c       CALL WWRITW(ITAP86,IFZV,1,IVEC86,IVEC86)                         
c       CALL WWRITW(ITAP86,NOP ,1,IVEC86,IVEC86)                         
c       CALL WWRITW(ITAP86,NALP,1,IVEC86,IVEC86)                         
c       CALL WWRITW(ITAP86,IJDC,1,IVEC86,IVEC86)                         
c       CALL WWRITW(ITAP86,IJIND,1,IVEC86,IVEC86)                        
c       CALL WWRITW(ITAP86,ZZ,INTOWP(IDIMZ),IVEC86,IVEC86)               
cedtool020998
       call wwritw(itap86,zz,intowp(numij),1,ivec86)
C     LEAVE PLACE FOR FOCK MATRIX IN THE MO BASIS (CCSD)                
c       IVEC86=IVEC86+INTOWP(NUMIJ)                                      
c       CALL WWRITW(ITAP86,IJDEP,NUMIJ,IVEC86,IVEC86)                    
       CALL RCLOSE(ITAP86,3)                                            
       write(*,*)
       WRITE(*,*)' Z VECTOR HAS BEEN WRITTEN TO FILE86'
       WRITE(*,*)'TO OBTAIN CORRECTED FIRST-ORDER PROPERTIES RUN OEPROP'
       write(*,*)
c       WRITE(*,*)'Z VECTOR HAS BEEN WRITTEN TO FILE86'                  
c       WRITE(*,*)'IJDEP OF DIM ',NUMIJ,' HAS BEEN WRITTEN TO FILE 86'   
c       IF(FLAG.EQ.0)THEN                                                
c        WRITE(*,*)                                                      
c        WRITE(*,*)                                                      
c        WRITE(*,*)                                                      
c        WRITE(*,*)' FLAG EQ ZERO   FILE 87 IS INITIALIZED'              
c        WRITE(*,*)' AND U=Z WRITTEN                      '              
c        CALL RFILE (ITAP87)                                             
c        IDUM=1                                                          
c        CALL WWRITW(ITAP87,IDUM,1,1     ,LOCVEC)                        
c        CALL WWRITW(ITAP87,IDIMZ,1,LOCVEC,LOCVEC)                       
c        CALL WWRITW(ITAP87,ZZ,INTOWP(IDIMZ),LOCVEC,LOCVEC)              
c        CALL RCLOSE(ITAP87,3)                                           
c        WRITE(*,*)' FLAG EQ ZERO   FILE 88 IS INITIALIZED'              
c        CALL RFILE (ITAP88)                                             
c        IDUM=0                                                          
c        CALL WWRITW(ITAP88,IDUM,1,1     ,LOCVEC)                        
c        CALL WWRITW(ITAP88,IDIMZ,1,LOCVEC,LOCVEC)                       
c        CALL RCLOSE(ITAP88,3)                                           
c       ELSE                                                             
c        WRITE(*,*)' FLAG NE ZERO   FILE 87 ALREADY IN DISK'             
c        WRITE(*,*)' FLAG NE ZERO   FILE 88 ALREADY IN DISK'             
c       ENDIF                                                            
c       CALL RFILE (ITAP88)                                              
c       CALL WREADW(ITAP88,NUMBER,1,1,LOCVEC)                            
c       CALL WREADW(ITAP88,IDIM  ,1,LOCVEC,LOCVEC)                       
c       LOCVEC=LOCVEC+NUMBER*INTOWP(IDIM)                                
c       CALL WWRITW(ITAP88,ZZ,INTOWP(IDIMZ),LOCVEC,LOCVEC)               
c       NUMBER=NUMBER+1                                                  
c       CALL WWRITW(ITAP88,NUMBER,1,1,LOCVEC)                            
c       CALL RCLOSE(ITAP88,3)                                            
c       WRITE(*,*)'Z VECTOR HAS BEEN WRITTEN TO FILE 88'                 
c       NORM=0.0D0                                                       
c       DO 1344 I=1,IDIMZ                                                
c        NORM=NORM+ZZ(I)*ZZ(I)                                           
c 1344  CONTINUE                                                         
c       NORM=DSQRT((NORM))                                               
c       WRITE(*,*)                                                       
c       WRITE(*,*) 'NORM OF Z VECTOR',NORM                               
c       WRITE(*,*)                                                       
c       NORM=NORM/IDIMZ                                                  
c       WRITE(*,*) 'STANDARD NORM = ',NORM                               
c      ENDIF                                                             
C                                                                       
      RETURN                                                            
      END                                                               
