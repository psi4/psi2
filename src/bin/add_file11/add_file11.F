      subroutine fentry(cc,ia,maxcor)
c     PROGRAM ADD_FILE11
c***********************************************************************
c   Last update on March 03, 2003 by Dr. Yukio Yamaguchi               *
c***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)                     
      dimension cc(maxcor),ia(maxcor*2)
      DIMENSION ATM(18)                            
      integer errcod, frdi, fwkclr, fwkadd
      character*80 TITLE
      character*49 label
      dimension chg(10),x(10),y(10),z(10),w(10)
      dimension chg0(10),x0(10),y0(10),z0(10)
      dimension gx0(10),gy0(10),gz0(10)
      dimension i30(200),coord(3,10)
      COMMON/CALIF/LPARA(1024),APARA(1024)                
      COMMON/MFSEC/MASTER,NSECT  
      COMMON/VIB101/NATOMS
      COMMON/VIB102/R(10,10)   
      COMMON/VIB103/IOFF(150),IPRNT                                     FOR00210
      DATA ATM / 1.007825D+00 , 4.00260D+00  , 7.01600D+00  ,           FOR00250
     *           9.01218D+00  , 11.00931D+00 , 12.00000D+00 ,           FOR00260
     *           14.00307D+00 , 15.99491D+00 , 18.99840D+00 ,           FOR00270
     *           19.99244D+00 , 22.9898D+00  , 23.98504D+00 ,           FOR00280
     *           26.98153D+00 , 27.97693D+00 , 30.97376D+00 ,           FOR00290
     *           31.97207D+00 , 34.96885D+00 , 39.948D+00      /        FOR00300
      DATA A0 / 0.52917706D+00 /
      DATA A00 / 0.0D+00 /                                              FOR00350
    1 FORMAT(//,2X,' The program to add a set of geometry and gradients'
     *        /,2x,'           (Version of 03/03/2003)'/)
    2 FORMAT(2I5)     
    3 FORMAT(A80)
    4 FORMAT(I5,F20.10,F15.10,2F20.10)                                  QPE00190
    5 FORMAT(4F20.10)                                                   QPE00200
    6 FORMAT(20X,3F20.10)                                               QPE00210
    7 FORMAT(//,2X,' PARAMETERS FOR THE CALCULATION'//                  FOR00400
     * 2X,' NATOMS = ',I10/                                             FOR00430
     * 2X,' N3N    = ',I10/                                             FOR00440
     * 2X,' IPOINT = ',I10/                                             FOR00430
     * 2X,' IPRNT  = ',I10/)                                            FOR00430
    8 FORMAT(//,2X,' INITIAL MOLECULAR GEOMETRY (IN BOHR)'/             FOR00470
     1 5X,4H NO.,11X,2H X,18X,2H Y,18X,2H Z,18X,2H W/)                  FOR00480
    9 FORMAT(2X,I7,5X,4F20.10)                                          FOR00490
   10 FORMAT(//,2X,' INITIAL MOLECULAR GEOMETRY (IN A)'/                FOR00500
     1 5X,4H NO.,11X,2H X,18X,2H Y,18X,2H Z,18X,2H W/)                  FOR00510
   11 FORMAT(//,2X,' PRESENT MOLECULAR GEOMETRY (IN BOHR)'/             FOR00470
     1 5X,4H NO.,11X,2H X,18X,2H Y,18X,2H Z,18X,2H W/)                  FOR00480
   12 FORMAT(//,2X,' PRESENT MOLECULAR GEOMETRY (IN A)'/                FOR00500
     1 5X,4H NO.,11X,2H X,18X,2H Y,18X,2H Z,18X,2H W/)                  FOR00510
   13 FORMAT(A49,'PERT = (',I4,')')
C                                                                       FOR00530
      call psinit('APPEND')
      CALL TSTART(6)                                                    FOR00540
      INPUT=5                                                           FOR00640
      ITAPE6=6                                                          FOR00650
      itap11=11
      ITAP30=30                                                         FOR00680
C                                                                       FOR00710
      IOFF(1)=0                                                         FOR00730
      DO 101 I=1,149                                                    FOR00740
      IOFF(I+1)=IOFF(I)+I                                               FOR00750
  101 CONTINUE                                                          FOR00760
C                                                                       FOR00770
C***********************************************************************FOR00780
C*SPECIFY THE TYPE OF CALCULATION                                       FOR00790
C***********************************************************************FOR00800
      call ffile(itap11,'file11',0)
c
      WRITE(6,1)                                                        FOR00720
      write(6,*) '  maxcor = ',maxcor
c
      errcod = fwkclr()
      errcod = fwkadd(5,'ADD_FILE11')
      errcod = fwkadd(5,'DEFAULT')
      iprnt = 0
      errcod = frdi('print',iprnt)
C                                                                       FOR00820
C   READ IN INFORMATION FROM TAPE11                                     QPE00560
      REWIND ITAP11                                                     QPE00570
c   For the reference point
      READ(ITAP11,3) TITLE                                              QPE00580
      WRITE(3,3) TITLE                                                  QPE00590
C*    WRITE(6,3) TITLE                                                  QPE00600
      READ(ITAP11,4) NATOMS,ENERG0,DIPX,DIPY,DIPZ
      n3n=natoms*3
      WRITE(3,4) NATOMS,ENERG0,DIPX,DIPY,DIPZ
C*    WRITE(6,4) NATOMS,ENERG0,DIPX,DIPY,DIPZ
      DO 102 I=1,NATOMS
      READ(ITAP11,5) CHG0(I),X0(I),Y0(I),Z0(I)                          QPE00660
      WRITE(3,5) CHG0(I),X0(I),Y0(I),Z0(I)                              QPE00670
C*    WRITE(6,5) CHG0(I),X0(I),Y0(I),Z0(I)                              QPE00680
  102 CONTINUE                                                          QPE00690
      DO 103 I=1,NATOMS
      READ(ITAP11,6) GX0(I),GY0(I),GZ0(I)                               QPE00710
      WRITE(3,6) GX0(I),GY0(I),GZ0(I)                                   QPE00720
C*    WRITE(6,6) GX0(I),GY0(I),GZ0(I)                                   QPE00730
  103 CONTINUE                                                          QPE00740
C                                                                       QPE00750
      IPOINT=1                                                          QPE00760
  200 CONTINUE                                                          QPE00770
c   For perturbed points
      READ(ITAP11,3,END=201) TITLE                                      QPE00780
      WRITE(3,3) TITLE                                                  QPE00790
C*    WRITE(6,3) TITLE                                                  QPE00800
      READ(ITAP11,4,END=201) NATOMP,ENERGP,DIPXX,DIPYY,DIPZZ            QPE00810
      WRITE(3,4) NATOMP,ENERGP,DIPXX,DIPYY,DIPZZ                        QPE00820
C*    WRITE(6,4) NATOMP,ENERGP,DIPXX,DIPYY,DIPZZ                        QPE00830
      DO 104 I=1,NATOMS
      READ(ITAP11,5,END=201) chg(i),x(i),y(i),z(i)
      WRITE(3,5) chg(i),x(i),y(i),z(i)
c*    WRITE(6,5) chg(i),x(i),y(i),z(i)
  104 CONTINUE         
      DO 105 I=1,NATOMS
      READ(ITAP11,6,END=201) GXX,GYY,GZZ                                QPE00900
      WRITE(3,6) GXX,GYY,GZZ                                            QPE00910
C*    WRITE(6,6) GXX,GYY,GZZ                                            QPE00920
  105 CONTINUE                                                          QPE00930
      IPOINT=IPOINT+1                                                   QPE00940
      GO TO 200          
c
  201 continue
      WRITE(6,7) NATOMS,N3N,IPOINT,IPRNT
C                                                                       CPC01820
C   PRINT INITIAL GEOMETRY                                              FOR01280
      WRITE(6,8)                                                        FOR01290
      DO 106 I=1,NATOMS
      II=CHG0(I)                                               
      W(I)=ATM(II)
      WRITE(6,9) I,X0(I),Y0(I),Z0(I),W(I)    
  106 CONTINUE                                                          FOR01390
      call dist(x0,y0,z0)
      WRITE(6,10)                                                       FOR01410
      DO 107 I=1,NATOMS
      X0(I)=X0(I)*A0        
      Y0(I)=Y0(I)*A0      
      Z0(I)=Z0(I)*A0 
      WRITE(6,9) I,X0(I),Y0(I),Z0(I),W(I)            
  107 CONTINUE                                                          FOR01470
      call dist(x0,y0,z0)
C
C   Find the present geometry
      CALL RFILE(ITAP30)                                                QPE01770
      CALL WREADW(ITAP30,I30,200,101,JUNK)                              QPE01780
      IEND=I30(1)                                                       QPE01790
      MPOINT=I30(2)                                                     QPE01800
      MCONST=I30(3)                                                     QPE01810
      MCALCS=I30(4)                                                     QPE01820
      NCALCS=I30(5)                                                     QPE01830
      JUNK=101+MCONST                                                   QPE01840
      CALL WREADW(ITAP30,IPT,1,JUNK,JUNK)                               QPE01850
      CALL WREADW(ITAP30,CHG,NATOMS*2,IPT,JUNK)
      IGEOP=100+MCONST+MPOINT+NCALCS                                    QPE01870
      CALL WREADW(ITAP30,LOCCAL,1,IGEOP,JUNK)                           QPE01880
      IGEOP=LOCCAL+80                                                   QPE01890
      CALL WREADW(ITAP30,COORD,N3N*2,IGEOP,JUNK)                        QPE01900
      call rclose(itap30,3)
C                                                                       QPE01910
C   PRINT PRESENT GEOMETRY
      WRITE(6,11)                                                       FOR01290
      DO 108 I=1,NATOMS
      II=CHG(I)                                                         FOR01340
      W(I)=ATM(II)
      x(i)=coord(1,i)
      y(i)=coord(2,i)
      z(i)=coord(3,i)
      WRITE(6,9) I,X(I),Y(I),Z(I),W(I)                                  FOR01380
  108 CONTINUE                                                          FOR01390
      call dist(x,y,z)
      WRITE(6,12)                                                       FOR01410
      DO 109 I=1,NATOMS
      X(I)=X(I)*A0        
      Y(I)=Y(I)*A0      
      Z(I)=Z(I)*A0 
      WRITE(6,9) I,X(I),Y(I),Z(I),W(I)                                  FOR01460
  109 CONTINUE                                                          FOR01470
      call dist(x,y,z)
c
c   Read in old geometries and gradients from file11
      rewind itap11
      do 112 iabc=1,ipoint
      READ(ITAP11,3) TITLE     
      READ(ITAP11,4) NATOMP,ENERGP,DIPXX,DIPYY,DIPZZ  
      DO 110 I=1,NATOMS
      READ(ITAP11,5) chg(i),x(i),y(i),z(i)
  110 CONTINUE
      DO 111 I=1,NATOMS
      READ(ITAP11,6) GXX,GYY,GZZ   
  111 continue
  112 continue
c
c   Add a new geometry and dummy gradients to file11
      dumx=a00
c     write(itap11,3) TITLE
      label = TITLE(1:49)
      write(itap11,13) label,ipoint
      write(itap11,4) NATOMS,dumx,dumx,dumx,dumx
      DO 113 I=1,NATOMS   
      write(itap11,5) chg(i),coord(1,i),coord(2,i),coord(3,i)
  113 CONTINUE
      DO 114 I=1,NATOMS
      write(itap11,6) dumx,dumx,dumx
  114 CONTINUE                                                          QPE00930
      rewind itap11
C                                                                       FOR01490
      CALL TSTOP(6)                                                     FOR01810
      STOP                                                              FOR01820
      END                                                               FOR01830
      SUBROUTINE DIST(x,y,z)       
      IMPLICIT REAL*8 (A-H,O-Z)                                         FOR01850
      dimension x(10),y(10),z(10)
      COMMON/VIB101/NATOMS
      COMMON/VIB102/R(10,10)
    1 FORMAT(//,2X,' INTERATOMIC DISTANCE MATRIX'/)                     FOR01890
C                                                                       FOR01900
      CALL ZERO(R,100)   
      DO 101 I=1,NATOMS
      DO 101 J=1,NATOMS
      IF(I.EQ.J) GO TO 101                                              FOR01960
      XD=X(I)-X(J)                                                      FOR01970
      YD=Y(I)-Y(J)                                                      FOR01980
      ZD=Z(I)-Z(J)                                                      FOR01990
      valu=DSQRT(XD*XD+YD*YD+ZD*ZD)   
      R(I,J)=valu
  101 CONTINUE                                                          FOR02010
      WRITE(6,1)                                                        FOR02030
      CALL MATOUT(R,10,10,NATOMS,NATOMS,6)    
c
      RETURN
      END                                                               FOR02060
