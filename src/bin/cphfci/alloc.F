      SUBROUTINE ALLOCA(ORBTBF,JCODE,IJORD,IOFF,IJDEP,                  
     .                 ORBSYM,DOFF,DADD,NMA1,itcscf)       
      IMPLICIT INTEGER (A-Z)                                            
cets030691
#include <error.h>
      integer ip, prcntr
      COMMON /TAPES/ NOUT,ICHK,IT30,IT42,IT47,IT52,IT58,IT94,IT96,IT98  
      COMMON /DIMSG/ SYMORB,NUMIJ,NBF,NMAX,NBFORB,NSYM,NDF,NORBS        
      COMMON /FUNCS/ IOCC,JOCC,KOCC,NTYPES,NBSET,NAT                    
      COMMON /OCCS/ NSPE,NSA,NSB,NALP,NOP,NDOC,LFERM,NUOC,IJIND,IJIND3  
      COMMON/OCCS2/NUMCOR,NCOR,FDOC,CORDOC,NUMVIR,FVIR,LUOC,NUMDEP,VUOC 
      COMMON/DEP/IJDC,IFZC,IFZV,IZDIM                                   
      common /types/ntydoc,ntyuoc
      DIMENSION ORBTBF(NORBS),JCODE(NORBS),IJORD(NUMIJ),IOFF(NUMIJ)     
      DIMENSION IJDEP(NUMIJ)                                            
      DIMENSION ITYP(9)                                                 
      DIMENSION ORBSYM(NORBS),DOFF(NUMIJ),DADD(NUMIJ)                   
C                                                                       
      WRITE(ICHK,2222)                                                  
 2222 FORMAT(/'  IN ALLOC')                                             
      DO 313 I=1,9                                                      
        ITYP(I)=0                                                       
  313 CONTINUE                                                          
C                                                                       
      DO 314 I=1,NUMIJ                                                  
        IOFF(I)=I*(I-1)/2                                               
  314 CONTINUE                                                          
C                                                                       
      DO 315 IORB=1,NORBS                                               
        NTYP=JCODE(IORB)                                                
        IF (NTYP.GT.9 .OR. NTYP.LE.0) THEN                              
          WRITE(NOUT,201) IORB,NTYP                                     
 201      FORMAT(/'   JCODE(',I3,')=',I4/                               
     1    '   I DO NOT RECOGNIZE THIS TYPE.')                           
          CALL ABORTS                                                   
        END IF                                                          
        ITYP(NTYP)=ITYP(NTYP)+1                                         
  315 CONTINUE                                                          
      IF (ITYP(1).NE.0 .OR. ITYP(2).NE.0) THEN                          
        WRITE(NOUT,202)                                                 
 202    FORMAT(/'  THERE ARE SOME FZC OR FZV IN THE DRT'/               
     1  '  I CAN NOT RUN A GRADIENT FROM THIS.')                        
        CALL ABORTS                                                     
      END IF                                                            
      NUOC=ITYP(4)+ITYP(6)                                              
      LFERM=NUOC+1                                                      
      NDOC=NUOC+ITYP(3)+ITYP(5)                                         
      ntydoc = ityp(3)+ityp(5)
      ntyuoc = ityp(4)+ityp(6)
CTJL                                                                    
      NUMVIR = ITYP(4)                                                  
      VUOC = NUMVIR*ITYP(6)                                             
      FVIR = ITYP(6) + 1                                                
      LUOC = FVIR - 1                                                   
      NUMCOR = ITYP(3)                                                  
      NCOR = LFERM + NUMCOR - 1                                         
      FDOC = LFERM + NUMCOR                                             
      IF(NUMCOR.EQ.0) NCOR = -999                                       
      CORDOC = ITYP(5)*NUMCOR                                           
      NUMDEP = VUOC + CORDOC                                            
COUT  WRITE(*,*) ' NUMDEP,CORDOC,VUOC',NUMDEP,CORDOC,VUOC               
C     DO 9240 IJ = 1,NUMIJ                                              
C     IJDEP(IJ) = NUMDEP + 1                                            
C9240 CONTINUE                                                          
      IJDC = 0                                                          
      DO 9200 I = FDOC,NDOC                                             
      DO 9210 J = LFERM,NCOR                                            
      IJ = IOFF(I) + J                                                  
      IJDC = IJDC + 1                                                   
      IJDEP(IJ) = IJDC                                                  
 9210 CONTINUE                                                          
 9200 CONTINUE                                                          
      IFZC  = FDOC-LFERM                                                
cets030691
      ip = prcntr('IS_ON BRIEF')
      if(ip.eq.0) then                                       
       WRITE(*,*) ' THE NUMBER OF FZCOR   ORBITALS IS',IFZC              
      endif
      DO 9220 I = FVIR,NUOC                                             
      DO 9230 J = 1,LUOC                                                
      IJ = IOFF(I) + J                                                  
      IJDC = IJDC + 1                                                   
      IJDEP(IJ) = IJDC                                                  
 9230 CONTINUE                                                          
 9220 CONTINUE                                                          
      IFZV  = NUOC-FVIR +1                                              
cets030691
      ip = prcntr('IS_ON BRIEF')
      if(ip.eq.0) then                                       
       WRITE(*,*) ' THE NUMBER OF FZVIR   ORBITALS IS',IFZV              
       WRITE(*,*) ' THE NUMBER OF EXTRA INDEPENDENT  PAIRS IS',IJDC      
      endif
      IF(NUMDEP.EQ.0) NUMDEP = 1                                        
COUT  WRITE(*,*) ' NUMCOR,NCOR,FDOC,CORDOC',NUMCOR,NCOR,FDOC,CORDOC     
CTJL                                                                    
      NOP=NDOC+1                                                        
      NALP=NDOC+ITYP(7)                                                 
      NSA=NALP+1                                                        
      NSB=NALP+ITYP(9)                                                  
      NSPE=ITYP(9)                                                      
      IF (NSB.NE.NORBS) THEN                                            
        WRITE(NOUT,203) NSB,NORBS                                       
 203    FORMAT(/'   ORBITALS MISCOUNTED IN ALLOC',2I5)                  
        CALL ABORTS                                                     
      END IF                                                            
      IF (ITYP(3)+ITYP(5) .EQ. 0) THEN                                  
        WRITE(NOUT,204)                                                 
 204    FORMAT(/'  THERE ARE NO DOUBLY OCCUPIED ORBITALS.'/             
     1  '  IT IS NECESSARY TO MODIFY THE PROGRAM. ')                    
        CALL ABORTS                                                     
      END IF                                                            
      IF (ITYP(7).EQ.0) THEN                                            
        NOP=0                                                           
        NALP=0                                                          
      ELSE                                                              
        DO 318 IORB=NOP,NALP                                            
          IBF=ORBTBF(IORB)                                              
          IF (JCODE(IBF).NE.7) GOTO 10                                  
  318   CONTINUE                                                        
      END IF                                                            
      IF (NSPE.EQ.0) THEN                                               
        NSA=0                                                           
        NSB=0                                                           
      ELSE                                                              
        DO 319 IORB=NSA,NSB                                             
          IBF=ORBTBF(IORB)                                              
          IF (JCODE(IBF).NE.9) GOTO 10                                  
  319   CONTINUE                                                        
      END IF                                                            
      GOTO 12                                                           
 10   WRITE(NOUT,205) IORB,IBF,JCODE(IBF)                               
 205  FORMAT(/'  THE DRT ORBITALS ARE OUT OF ORDER.'/                   
     1'  IORB=',I5,'  IBF=',I5,'  NTYPE=',I5)                           
      CALL ABORTS                                                       
 12   CONTINUE  
cets030691
      ip = prcntr('IS_ON BRIEF')
      if(ip.eq.0) then                                       
       WRITE(NOUT,206) 1,NUOC,LFERM,NDOC                                 
 206   FORMAT(/'  THE UNOCCUPIED ORBITALS:     ',2I5/                    
     1 '  THE DOUBLY OCCUPIED ORBITALS:',2I5/)                           
       IF (NOP.NE.0) WRITE(NOUT,207) NOP,NALP                            
 207   FORMAT('  THE SINGLY OCCUPIED ORBITALS:',2I5/)                    
       IF (NSPE.NE.0) WRITE(NOUT,208) NSA,NSB                            
 208   FORMAT('  THE SPECIAL ORBITALS:        ',2I5/)                    
      endif
C                                                                       
C  CONSRTUCT DRT TO INDEP PAIR ORDERING ARRAY                           
      DO 320 II=1,NUMIJ                                                 
        IJORD(II)=0                                                     
  320 CONTINUE                                                          
      IJIND=0                                                           
      DO 321 I=LFERM,NDOC                                               
        IO=IOFF(I)                                                      
        DO 322 J=1,NUOC                                                 
          IJ=IO+J                                                       
          IJIND=IJIND+1                                                 
          IJORD(IJ)=IJIND                                               
  322   CONTINUE                                                        
  321 CONTINUE                                                          
      IF (NOP.EQ.0) GOTO 14                                             
      DO 323 I=NOP,NALP                                                 
        IO=IOFF(I)                                                      
        DO 324 J=1,NDOC                                                 
          IJ=IO+J                                                       
          IJIND=IJIND+1                                                 
          IJORD(IJ)=IJIND                                               
  324   CONTINUE                                                        
  323 CONTINUE                                                          
 14   CONTINUE                                                          
      IF (NSPE.EQ.0) GOTO 16                                            
      DO 325 I=NSA,NSB                                                  
        IO=IOFF(I)                                                      
        DO 326 J=1,I-1                                                  
          IJ=IO+J                                                       
          IJIND=IJIND+1                                                 
          IJORD(IJ)=IJIND                                               
  326   CONTINUE                                                        
  325 CONTINUE                                                          
 16   CONTINUE                                                          
cets030691
      ip = prcntr('IS_ON BRIEF')
      if(ip.eq.0) then
       write(NOUT,209) IJIND,IJIND+2                                     
 209   FORMAT('  THE NUMBER OF INDEPENDENT PAIRS IS',I5/                 
     1 '  PLUS THE TWO CI COEFFICIENTS MAKES',I5/)                       
       IEXC=IFZC*(NDOC-NUOC-IFZC)                                        
       WRITE(*,*)'THE NUMBER OF INDEPENDENT COR-DOC ROTATIONS IS', IEXC  
       IEXV=IFZV*(NUOC-IFZV)                                             
       WRITE(*,*)'THE NUMBER OF INDEPENDENT VIR-UOC ROTATIONS IS', IEXV  
      endif
      IZDIM=IJIND+IEXC+IEXV                                             
CGES                                                                    
C     CONSTRUCT SYMMETRY-PACKING ARRAY                                  
C                                                                       
      NOFF=0                                                            
      KTSYM=0                                                           
      DO 561 I=LFERM,NDOC                                               
      ISYM=ORBSYM(I)-1                                                  
      DO 560 J=1,NUOC                                                   
      JSYM=ORBSYM(J)-1                                                  
      IJSYM=IEOR(ISYM,JSYM)+1                                           
      IF(IJSYM.NE.1) GO TO 560                                          
      NADD=0                                                            
      IJ=IOFF(I)+J                                                      
      DOFF(IJ)=NOFF                                                     
      if(itcscf.ne.0) noff = noff+2
      DO 530 K=LFERM,NDOC                                               
      KSYM=ORBSYM(K)-1                                                  
      DO 520 L=1,NUOC                                                   
      LSYM=ORBSYM(L)-1                                                  
      KL=IOFF(K)+L                                                      
      KLSYM=IEOR(KSYM,LSYM)+1                                           
      IF(KLSYM.NE.1)GO TO 510                                           
      NOFF=NOFF+1                                                       
      NADD=NADD+1                                                       
      KTSYM=KTSYM+1                                                     
      DADD(KL)=NADD                                                     
      ITOT=DOFF(IJ)+DADD(KL)                                            
c     WRITE(6,657)I,J,DOFF(IJ),K,L,DADD(KL),ITOT                        
  657 FORMAT(' I=',I3,' J=',I3,' DOFF=',I6,' K=',I3,' L=',I3,' DADD=',I6
     .       ,'  TOT=',I6)                                              
  510 CONTINUE                                                          
  520 CONTINUE                                                          
  530 CONTINUE                                                          
      IF(NOP.EQ.0)GO TO 541                                             
      DO 531 K=NOP,NALP                                                 
      KSYM=ORBSYM(K)-1                                                  
      DO 521 L=1,NDOC                                                   
      LSYM=ORBSYM(L)-1                                                  
      KL=IOFF(K)+L                                                      
      KLSYM=IEOR(KSYM,LSYM)+1                                           
      IF(KLSYM.NE.1)GO TO 511                                           
      NOFF=NOFF+1                                                       
      NADD=NADD+1                                                       
      KTSYM=KTSYM+1                                                     
      DADD(KL)=NADD                                                     
      ITOT=DOFF(IJ)+DADD(KL)                                            
C     WRITE(6,657)I,J,DOFF(IJ),K,L,DADD(KL),ITOT                        
  511 CONTINUE                                                          
  521 CONTINUE                                                          
  531 CONTINUE                                                          
  541 CONTINUE                                                          
      IF(NSPE.EQ.0)GO TO 542                                            
      DO 532 K=NSA,NSB                                                  
      KSYM=ORBSYM(K)-1                                                  
      DO 522 L=1,K-1                                                    
      LSYM=ORBSYM(L)-1                                                  
      KL=IOFF(K)+L                                                      
      KLSYM=IEOR(KSYM,LSYM)+1                                           
      IF(KLSYM.NE.1)GO TO 512                                           
      NOFF=NOFF+1                                                       
      NADD=NADD+1                                                       
      KTSYM=KTSYM+1                                                     
      DADD(KL)=NADD                                                     
      ITOT=DOFF(IJ)+DADD(KL)                                            
C     WRITE(6,657)I,J,DOFF(IJ),K,L,DADD(KL),ITOT                        
  512 CONTINUE                                                          
  522 CONTINUE                                                          
  532 CONTINUE                                                          
  542 CONTINUE                                                          
  560 CONTINUE                                                          
  561 CONTINUE                                                          
C                                                                       
      IF(NOP.EQ.0)GO TO 662                                             
      DO 661 I=NOP,NALP                                                 
      ISYM=ORBSYM(I)-1                                                  
      DO 660 J=1,NDOC                                                   
      JSYM=ORBSYM(J)-1                                                  
      IJSYM=IEOR(ISYM,JSYM)+1                                           
      IF(IJSYM.NE.1) GO TO 660                                          
      NADD=0                                                            
      IJ=IOFF(I)+J                                                      
      DOFF(IJ)=NOFF                                                     
      if(itcscf.ne.0) noff = noff+2
      DO 630 K=LFERM,NDOC                                               
      KSYM=ORBSYM(K)-1                                                  
      DO 620 L=1,NUOC                                                   
      LSYM=ORBSYM(L)-1                                                  
      KL=IOFF(K)+L                                                      
      KLSYM=IEOR(KSYM,LSYM)+1                                           
      IF(KLSYM.NE.1)GO TO 610                                           
      NOFF=NOFF+1                                                       
      NADD=NADD+1                                                       
      KTSYM=KTSYM+1                                                     
      DADD(KL)=NADD                                                     
      ITOT=DOFF(IJ)+DADD(KL)                                            
C     WRITE(6,657)I,J,DOFF(IJ),K,L,DADD(KL),ITOT                        
  610 CONTINUE                                                          
  620 CONTINUE                                                          
  630 CONTINUE                                                          
      IF(NOP.EQ.0)GO TO 641                                             
      DO 631 K=NOP,NALP                                                 
      KSYM=ORBSYM(K)-1                                                  
      DO 621 L=1,NDOC                                                   
      LSYM=ORBSYM(L)-1                                                  
      KL=IOFF(K)+L                                                      
      KLSYM=IEOR(KSYM,LSYM)+1                                           
      IF(KLSYM.NE.1)GO TO 611                                           
      NOFF=NOFF+1                                                       
      NADD=NADD+1                                                       
      KTSYM=KTSYM+1                                                     
      DADD(KL)=NADD                                                     
      ITOT=DOFF(IJ)+DADD(KL)                                            
C     WRITE(6,657)I,J,DOFF(IJ),K,L,DADD(KL),ITOT                        
  611 CONTINUE                                                          
  621 CONTINUE                                                          
  631 CONTINUE                                                          
  641 CONTINUE                                                          
      IF(NSPE.EQ.0)GO TO 642                                            
      DO 632 K=NSA,NSB                                                  
      KSYM=ORBSYM(K)-1                                                  
      DO 622 L=1,K-1                                                    
      LSYM=ORBSYM(L)-1                                                  
      KL=IOFF(K)+L                                                      
      KLSYM=IEOR(KSYM,LSYM)+1                                           
      IF(KLSYM.NE.1)GO TO 612                                           
      NOFF=NOFF+1                                                       
      NADD=NADD+1                                                       
      KTSYM=KTSYM+1                                                     
      DADD(KL)=NADD                                                     
      ITOT=DOFF(IJ)+DADD(KL)                                            
C     WRITE(6,657)I,J,DOFF(IJ),K,L,DADD(KL),ITOT                        
  612 CONTINUE                                                          
  622 CONTINUE                                                          
  632 CONTINUE                                                          
  642 CONTINUE                                                          
  660 CONTINUE                                                          
  661 CONTINUE                                                          
  662 CONTINUE                                                          
C                                                                       
      IF(NSPE.EQ.0) GO TO 762                                           
      DO 761 I=NSA,NSB                                                  
      ISYM=ORBSYM(I)-1                                                  
      DO 760 J=1,I-1                                                    
      JSYM=ORBSYM(J)-1                                                  
      IJSYM=IEOR(ISYM,JSYM)+1                                           
      IF(IJSYM.NE.1) GO TO 760                                          
      NADD=0                                                            
      IJ=IOFF(I)+J                                                      
      DOFF(IJ)=NOFF                                                     
      if(itcscf.ne.0) noff = noff+2
      DO 730 K=LFERM,NDOC                                               
      KSYM=ORBSYM(K)-1                                                  
      DO 720 L=1,NUOC                                                   
      LSYM=ORBSYM(L)-1                                                  
      KL=IOFF(K)+L                                                      
      KLSYM=IEOR(KSYM,LSYM)+1                                           
      IF(KLSYM.NE.1)GO TO 710                                           
      NOFF=NOFF+1                                                       
      NADD=NADD+1                                                       
      KTSYM=KTSYM+1                                                     
      DADD(KL)=NADD                                                     
      ITOT=DOFF(IJ)+DADD(KL)                                            
C     WRITE(6,657)I,J,DOFF(IJ),K,L,DADD(KL),ITOT                        
  710 CONTINUE                                                          
  720 CONTINUE                                                          
  730 CONTINUE                                                          
      IF(NOP.EQ.0)GO TO 741                                             
      DO 731 K=NOP,NALP                                                 
      KSYM=ORBSYM(K)-1                                                  
      DO 721 L=1,NDOC                                                   
      LSYM=ORBSYM(L)-1                                                  
      KL=IOFF(K)+L                                                      
      KLSYM=IEOR(KSYM,LSYM)+1                                           
      IF(KLSYM.NE.1)GO TO 711                                           
      NOFF=NOFF+1                                                       
      NADD=NADD+1                                                       
      KTSYM=KTSYM+1                                                     
      DADD(KL)=NADD                                                     
      ITOT=DOFF(IJ)+DADD(KL)                                            
C     WRITE(6,657)I,J,DOFF(IJ),K,L,DADD(KL),ITOT                        
  711 CONTINUE                                                          
  721 CONTINUE                                                          
  731 CONTINUE                                                          
  741 CONTINUE                                                          
      IF(NSPE.EQ.0)GO TO 742                                            
      DO 732 K=NSA,NSB                                                  
      KSYM=ORBSYM(K)-1                                                  
      DO 722 L=1,K-1                                                    
      LSYM=ORBSYM(L)-1                                                  
      KL=IOFF(K)+L                                                      
      KLSYM=IEOR(KSYM,LSYM)+1                                           
      IF(KLSYM.NE.1)GO TO 712                                           
      NOFF=NOFF+1                                                       
      NADD=NADD+1                                                       
      KTSYM=KTSYM+1                                                     
      DADD(KL)=NADD                                                     
      ITOT=DOFF(IJ)+DADD(KL)                                            
C     WRITE(6,657)I,J,DOFF(IJ),K,L,DADD(KL),ITOT                        
  712 CONTINUE                                                          
  722 CONTINUE                                                          
  732 CONTINUE                                                          
  742 CONTINUE                                                          
  760 CONTINUE                                                          
  761 CONTINUE                                                          
  762 CONTINUE                                                          
C                                                                       
      NMA1=INT(SQRT(KTSYM+0.))                                          
      ip = prcntr('IS_ON BRIEF')
      if(ip.eq.0) then
       WRITE(*,*)                                                        
       WRITE(*,*)'NUMBER OF A TERMS IN A1 SYMMETRY ',KTSYM               
       WRITE(*,*)'NUMBER OF TRANSITIONS NEEDED IS  ',NMA1                
      endif
      RETURN                                                            
      END                                                               
