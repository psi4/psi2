      SUBROUTINE Energy(H,Q,G,LBLI,BUFI)
      IMPLICIT REAL*8 (A-H,O-Z)                                         SMM13270
      DIMENSION H(NTRI),Q(NTRI),G(NOCTOT)
      DIMENSION LBLI(MAXBF4),BUFI(MAXBF2)                               SMM13290
      COMMON/BASIS/NBASIS,NTRI,NTOT,NSYM
      COMMON/ENRGY/ESCF,ENUC                                            SMM13310
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          SMM13320
      common/posit/ioemc,iopmc,ilgmc
      COMMON/SIGNS/IOFF(32897),IPRNT  
      COMMON/TAPES/ITAP35,ITAP56
      COMMON/CI101/NOCMO,NOCTRI,NOCTOT
      DATA A00,HALF,ONE,TWO,EIGHT / 0.0D+00 , 0.5D+00 , 1.0D+00 ,       SMM13380
     1                              2.0D+00 , 8.0D+00 /                 SMM13390
      data xlimit / 1.0d-10 /
    1 FORMAT(//,2X,' One electron MO integral matrix in Energy'/)
    2 FORMAT(//,2X,' One particle density matrix in Energy'/)
    3 FORMAT(2X,8I5,2F15.8)                                             SMM13420
    4 FORMAT(//,2X,' CI ENERGIES IN ENERGY'/
     * 2X,' ESCF      = ',F20.10/                                       SMM13440
     * 2X,' ENUC      = ',F20.10/                                       SMM13450
     * 2X,' ECI-ONE   = ',F20.10/                                       SMM13460
     * 2X,' ECI-TWO   = ',F20.10/                                       SMM13470
     * 2X,' ECI-ELEC  = ',F20.10/                                       SMM13480
     * 2X,' ECI-TOTAL = ',F20.10/                                       SMM13490
     * 2X,' ECI-STATE = ',F20.10)                                       SMM13500
C                                                                       SMM13510
C   Read in the one-electron integrals and one-PDM
      call mread(h,ioemc)
      call mread(q,iopmc)
      if(iprnt.gt.4) then
        write(6,1)
        call print(h,ntri,nbasis,6)
        write(6,2)
        call print(q,ntri,nbasis,6)
      end if
C                                                                       SMM13570
C   CALCULATE CI ENERGY                                                 SMM13740
C   ONE ELECTRON PART                                                   SMM13750
      ECI1=A00                                                          SMM13780
      IJ=O                                                              SMM13790
      DO 104 I=1,nocmo
      DO 104 J=1,I                                                      SMM13810
      IJ=IJ+1                                                           SMM13820
      FAC=TWO                                                           SMM13830
      IF(I.EQ.J) FAC=ONE                                                SMM13840
      ECI1=ECI1+Q(IJ)*H(IJ)*FAC                                         SMM13850
  104 continue
C                                                                       SMM13860
C   Read in the two-particle density matrix
      call rdin22(g,noctot,itap56,noctot,nocmo,lbli,bufi)
c
C   TWO ELECTRON PART                                                   SMM13870
      CALL RFILE(ITAP35)                                                SMM13880
      ECI2=A00                                                          SMM13890
      INTI=0                                                            SMM13900
      NBLI=0                                                            SMM13910
  200 IBL=0                                                             SMM13920
      NBLI=NBLI+1                                                       SMM13930
      CALL SREAD(ITAP35,LBLI,MAXBF4)                                    SMM13940
  201 IBL=IBL+1                                                         SMM13950
      JA=LBLI(IBL+IBL-1)                                                SMM13960
      IF(JA.EQ.0) GO TO 205                                             SMM13970
      JAA=LBLI(IBL+IBL)                                                 SMM13980
      CALL UNPACK(I,J,K,L,IX,JA,JAA)                                    SMM13990
      IMAX=MAX0(I,J,K,L)                                                SMM14000
      if(imax.gt.nocmo) go to 301
      IJ=IOFF(I)+J                                                      SMM14020
      KL=IOFF(K)+L                                                      SMM14030
      IJKL=IOFF(IJ)+KL                                                  SMM14040
      valm=bufi(maxbuf+ibl)
      valt=g(ijkl)
      if(dabs(valm).lt.xlimit) go to 301
      if(dabs(valt).lt.xlimit) go to 301
      valx=valm*valt
C     IF(IPRNT.LE.5) GO TO 202                                          SMM14060
C     WRITE(6,3) I,J,K,L,IJ,KL,IX,IJKL,VAL,G(IJKL)                      SMM14070
  202 fac=eight*valx
      IF(I.EQ.J) FAC=FAC*HALF                                           SMM14090
      IF(K.EQ.L) FAC=FAC*HALF                                           SMM14100
      IF(I.EQ.K.AND.J.EQ.L) FAC=FAC*HALF                                SMM14110
      eci2=eci2+fac
  301 continue
      IF(IBL.EQ.MAXBUF) GO TO 203                                       SMM14140
      GO TO 201                                                         SMM14150
  203 INTI=INTI+MAXBUF                                                  SMM14160
      GO TO 200                                                         SMM14170
  205 INTI=INTI+IBL                                                     SMM14180
C                                                                       SMM14190
      ECIEL=ECI1+ECI2                                                   SMM14200
      ECIT=ECIEL+ENUC                                                   SMM14210
      ESTATE=ECIT-ESCF                                                  SMM14220
      WRITE(6,4) ESCF,ENUC,ECI1,ECI2,ECIEL,ECIT,ESTATE                  SMM14230
C                                                                       SMM14240
      CALL RCLOSE(ITAP35,3)                                             SMM14250
      RETURN                                                            SMM14260
      END                                                               SMM14270
