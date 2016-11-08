c===================================================================
c   NORMCO.F
C   with the CINFV3_MP06 main program
c===================================================================
c   For all modes
c===================================================================
c     Last updated on December 18, 2006 by Dr. Yukio Yamaguchi
c     Modification for the MOLPRO program
c===================================================================
c     Last updated on July 21, 2005 by Dr. Yukio Yamaguchi
c===================================================================
c     Last updated on February 11, 2003 by Dr. Yukio Yamaguchi
c===================================================================
C
      SUBROUTINE NORMCO(G,F)
      IMPLICIT REAL*8 (A-H,O-Z)                                    
      dimension g(6),f(3,3)
      dimension s(6),ww(3,3),h(3,3)
      dimension aa(3,3),y(3,3),elin(3,3)
      dimension ell(3,3),pot(3,3),eked(3,3)
      dimension ee(3,3),gin(3,3)
      dimension valu(3),root(3),fv1(3),fv2(3)
      dimension rogam(3),freq(3)
      common/com102/iprnt
      COMMON/VIB101/WAVE
      DATA A00,HALF,ONE / 0.0D+00 , 0.5D+00 , 1.0D+00 /                 NOR18230
      DATA RLIMIT / 1.0D-06 /
    1 FORMAT(//,2X,' G MATRIX'/)                                        NOR18240
    2 FORMAT(//,2X,' F MATRIX'/)                                        NOR18250
    3 FORMAT(//,2X,' EIGENVALUES OF G MATRIX')                          NOR18260
    4 FORMAT(13X,8H GAMMA-G,9X,13H ROOT GAMMA-G/)                       NOR18270
    5 FORMAT(2X,I2,2(5X,F12.5))                                         NOR18280
    6 FORMAT(//,2X,' H MATRIX'/)                                        NOR18290
    7 FORMAT(//,2X,' EIGENVALUES OF H MATRIX')                          NOR18300
    8 FORMAT(13X,8H GAMMA-H,6X,10H FREQUENCY/)                          NOR18310
    9 FORMAT(//,2X,' EIGENVECTORS OF H MATRIX'/)                        NOR18320
   10 FORMAT(//,2X,' ZERO-POINT VIBRATIONAL ENERGY'//                   NOR18330
     1          2X,F20.5,' IN CM-1',5X,F20.5,' IN KCAL/MOL'/)           NOR18340
   11 FORMAT(//,2X,' L INVERSE MATRIX'/)                                NOR18350
   12 FORMAT(//,2X,' L MATRIX'/)                                        NOR18360
   13 FORMAT(//,2X,' L X L(T) = G'/)                                    NOR18370
   14 FORMAT(//,2X,' L X L(-1) = E'/)                                   NOR18380
   15 FORMAT(//,2X,' POTENTIAL ENERGY DISTRIBUTION'/)                   NOR18390
   16 FORMAT(//,2X,' KINETIC ENERGY DISTRIBUTION'/)                     NOR18400
C                                                                       NOR18410
      WRITE(6,1)       
      call print(g,6,3,6)
      WRITE(6,2)               
      call matout(f,3,3,3,3,6)
C                                                                       NOR18460
C   THE CALCULATION OF EIGENVALUE AND EIGENVECTOR OF G MATRIX           NOR18530
      call rsp(3,3,6,g,valu,1,ww,fv1,fv2)
      WRITE(6,3)                                                        NOR18550
      WRITE(6,4)                                                        NOR18560
      DO 101 I=1,3       
      ROGAM(I)=A00                                                      NOR18580
      IF(VALU(I).LE.A00) GO TO 201                                      NOR18590
      ROGAM(I)=DSQRT(VALU(I))                                           NOR18600
  201 continue
      WRITE(6,5) I,VALU(I),ROGAM(I)                                     NOR18610
  101 CONTINUE                                                          NOR18620
C                                                                       NOR18630
C   THE CALCULATION OF LGZERO MATRIX                                    NOR18640
      DO 102 I=1,3
      DO 102 J=1,3
      WW(I,J)=WW(I,J)*ROGAM(J)                                          NOR18670
  102 continue
C                                                                       NOR18680
C   THE CALCULATION OF H MATRIX                                         NOR18690
C     H MATRIX CORRESPONDS TO W*F*W                                     NOR18700
      call mtxmpy(ww,3,f,3,h,3,aa,3,3,5)
      WRITE(6,6)                                                        NOR18720
      IJ=0                                                              NOR18730
      DO 103 I=1,3
      DO 103 J=1,I                                                      NOR18750
      IJ=IJ+1                                                           NOR18760
      S(IJ)=H(I,J)                                                      NOR18770
  103 CONTINUE                                                          NOR18780
      call print(s,6,3,6)
C                                                                       NOR18800
C   THE CALCULATION OF FREQUENCIES                                      NOR18810
      call rsp(3,3,6,s,root,1,aa,fv1,fv2)
      WRITE(6,7)                                                        NOR18830
      WRITE(6,8)                                                        NOR18840
      DO 104 I=1,3
      FREQ(I)=DSIGN(ONE,ROOT(I))*WAVE*DSQRT(DABS(ROOT(I)))              NOR18860
      WRITE(6,5) I,ROOT(I),FREQ(I)                                      NOR18870
  104 CONTINUE                                                          NOR18880
      WRITE(6,9)                                                        NOR18890
      call frqout(aa,freq,3,3,3,3,6)
C                                                                       NOR18910
C   THE CALCULATION OF L INVERSE MATRIX                                 NOR19030
      DO 106 I=1,3
      DO 106 J=1,3
      Y(I,J)=AA(J,I)/VALU(J)                                            NOR19060
  106 CONTINUE                                                          NOR19070
      call mtxmpy(y,3,ww,3,elin,3,ee,3,3,3)
      WRITE(6,11)                                                       NOR19090
      call frqout(elin,freq,3,3,3,3,6)
C                                                                       NOR19110
C   THE CALCULATION OF L MATRIX                                         NOR19120
      call mtxmpy(ww,3,aa,3,ell,3,ee,3,3,1) 
      WRITE(6,12)                                                       NOR19140
      call frqout(ell,freq,3,3,3,3,6)
C                                                                       NOR19160
C   THE EXAMINATION OF ABOVE CALCULATION ACCORDING TO G MATRIX          NOR19170
      call mtxmpy(ell,3,ell,3,aa,3,ee,3,3,3)
      WRITE(6,13)                                                       NOR19190
      call matout(aa,3,3,3,3,6)
      call mtxmpy(ell,3,elin,3,aa,3,ee,3,3,1)
      WRITE(6,14)                                                       NOR19220
      call matout(aa,3,3,3,3,6)
C                                                                       NOR19240
C   THE CALCULATION OF POTENTIAL ENERGY DISTRIBUTION                    NOR19250
      DO 110 I=1,3
      DO 110 J=1,3
      POT(I,J)=A00                                                      NOR19280
  110 continue
      DO 111 I=1,3
      DO 111 J=1,3
      if(dabs(root(j)) .lt. rlimit ) go to 111
      POT(I,J)=POT(I,J)+(ELL(I,J)**2*F(I,I))/ROOT(J)                    NOR19310
  111 CONTINUE                                                          NOR19320
      WRITE(6,15)                                                       NOR19330
      call frqout(pot,freq,3,3,3,3,6)
C                                                                       NOR19350
C   THE CALCULATION OF KINETIC ENERGY DISTRIBUTION                      NOR19360
      call mtxmpy(elin,3,elin,3,gin,3,aa,3,3,2)
      DO 112 I=1,3
      DO 112 J=1,3
      EKED(I,J)=A00                                                     NOR19400
  112 continue
      DO 113 I=1,3
      DO 113 J=1,3
      EKED(I,J)=EKED(I,J)+GIN(I,I)*ELL(I,J)**2                          NOR19430
  113 continue
      WRITE(6,16)                                                       NOR19440
      call frqout(eked,freq,3,3,3,3,6)
      write(6,*) ' '
C                                                                       NOR19460
      RETURN                                                            NOR19470
      END                                                               NOR19480
