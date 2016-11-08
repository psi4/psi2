      SUBROUTINE LAGRAN(H,Q,X,G,CC,LBLI,BUFI)
      IMPLICIT REAL*8 (A-H,O-Z)                                         SMM14290
      DIMENSION H(NTRI),Q(NTRI),X(NBASIS,NBASIS)
      DIMENSION G(NOCTOT),CC(NTOT)
      DIMENSION LBLI(MAXBF4),BUFI(MAXBF2)                               SMM14320
      COMMON/BASIS/NBASIS,NTRI,NTOT,NSYM
      COMMON/CIENR/CIELEC,CITOT                                         SMM14340
      COMMON/ENRGY/ESCF,ENUC                                            SMM14350
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          SMM14360
      common/posit/ioemc,iopmc,ilgmc
      COMMON/SIGNS/IOFF(32897),IPRNT  
      COMMON/TAPES/ITAP35,ITAP56
      COMMON/CI101/NOCMO,NOCTRI,NOCTOT
      DATA A00,HALF,TWO / 0.0D+00 , 0.5D+00 , 2.0D+00 /                 SMM14410
    1 FORMAT(//,2X,' ONE ELECTRON MO INTEGRALS in Lagran'/)
    2 FORMAT(//,2X,' ONE PARTICLE DENSITY MATRIX in Lagran'/)
    3 FORMAT(//,2X,' THE LAGRANGIAN MATRIX (ONE ELECTRON PART)'/)       SMM14440
    4 FORMAT(//,2X,' THE LAGRANGIAN MATRIX'/)                           SMM14450
    5 FORMAT(//,2X,' CI ENERGIES IN LAGRAN'/                            SMM14460
     * 2X,' ESCF      = ',F20.10/                                       SMM14470
     * 2X,' ENUC      = ',F20.10/                                       SMM14480
     * 2X,' ECI-ONE   = ',F20.10/                                       SMM14490
     * 2X,' ECI-TWO   = ',F20.10/                                       SMM14500
     * 2X,' ECI-ELEC  = ',F20.10/                                       SMM14510
     * 2X,' ECI-TOTAL = ',F20.10/                                       SMM14520
     * 2X,' ECI-STATE = ',F20.10)                                       SMM14530
C                                                                       SMM14540
C   FORM THE LAGRANGIAN MATRIX                                          SMM14550
C                                                                       SMM14560
      if(iprnt.gt.4) then
        WRITE(6,1)                                            
        CALL PRINT(H,NTRI,NBASIS,6)
        WRITE(6,2)                                       
        CALL PRINT(Q,NTRI,NBASIS,6)                  
      end if
C                                                                       SMM14630
C   READ IN TWO ELECTRON MATRICES                                       SMM14640
      CALL RDINT2(CC,NTOT,ITAP35,NTOT,LBLI,BUFI) 
c
c.....................................................................
cc* The G matrix was already readin in subroutine ``energy''
cc*     CALL RDIN22(G,NOCTOT,ITAP56,NOCTOT,NOCMO,LBLI,BUFI) 
c.....................................................................
C                                                                       SMM14690
C   ONE ELECTRON PART                                                   SMM14700
      CALL ZERO(X,NBASIS*NBASIS)                                        SMM14710
      DO 102 I=1,NBASIS                                                 SMM14720
      DO 102 J=1,NOCMO
      EONE=A00                                                          SMM14740
      DO 101 M=1,NOCMO
      IM=I21D(I,M)                                                      SMM14760
      JM=I21D(J,M)                                                      SMM14770
      EONE=EONE+Q(JM)*H(IM)                                             SMM14780
  101 CONTINUE                                                          SMM14790
      X(I,J)=EONE                                                       SMM14800
  102 CONTINUE                                                          SMM14810
C                                                                       SMM14820
C   CALCULATE THE ONE ELECTRON ENERGY                                   SMM14830
      ECI1=A00                                                          SMM14840
      DO 103 I=1,NBASIS
      ECI1=ECI1+X(I,I)                                                  SMM14860
  103 CONTINUE                                                          SMM14870
      if(iprnt.gt.4) then
        WRITE(6,3)                                          
        CALL MATOUT(X,NBASIS,NBASIS,NBASIS,NBASIS,6)     
      end if
C                                                                       SMM14920
C   TWO ELECTRON PART                                                   SMM14930
      DO 105 I=1,NBASIS                                                 SMM14940
      DO 105 J=1,NOCMO
      ETWO=A00       
      DO 104 M=1,NOCMO
      IM=I21D(I,M)                                                      SMM14980
      JM=I21D(J,M)                                                      SMM14990
      DO 104 K=1,NOCMO
      DO 104 L=1,NOCMO
      KL=I21D(K,L)                                                      SMM15020
      IMKL=I21D(IM,KL)                                                  SMM15030
      JMKL=I21D(JM,KL)                                                  SMM15040
      ETWO=ETWO+G(JMKL)*CC(IMKL)                                        SMM15050
  104 CONTINUE                                                          SMM15060
      X(I,J)=X(I,J)+ETWO*TWO                                            SMM15070
  105 CONTINUE                                                          SMM15080
      if(iprnt.ge.2) then
        WRITE(6,4)                                        
        CALL MATOUT(X,NBASIS,NBASIS,NBASIS,NBASIS,6)   
      end if
C                                                                       SMM15130
C   CALCULATE THE TWO ELECTRON ENERGY                                   SMM15140
      ECIEL=A00                                                         SMM15150
      DO 107 I=1,NBASIS                                                 SMM15160
      ECIEL=ECIEL+X(I,I)                                                SMM15170
  107 CONTINUE                                                          SMM15180
      ECI2=(ECIEL-ECI1)*HALF                                            SMM15190
      ECIEL=ECI1+ECI2                                                   SMM15200
      ECIT=ECIEL+ENUC                                                   SMM15210
      ESTATE=ECIT-ESCF                                                  SMM15220
      CIELEC=ECIEL                                                      SMM15230
      CITOT=ECIT                                                        SMM15240
      WRITE(6,5) ESCF,ENUC,ECI1,ECI2,ECIEL,ECIT,ESTATE                  SMM15260
C                                                                       SMM15270
c   Store the Lagrangian matrix in the master file
      call mwrit(x,ilgmc)
      write(6,*) '  '
      write(6,*) '  The lagrangian matrix is stored in the master file'
c
      RETURN                                                            SMM15290
      END                                                               SMM15300
