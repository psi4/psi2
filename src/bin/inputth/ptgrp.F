      SUBROUTINE PTGRP(T,C,D,INVT)
C***********************************************************************
C     T,C AND D ARE IMPLICITLY EQUIVALENCED.                           *
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /MAXDIM/ MXAT,MXTRN,MXPR,MXSH,MXSYM,MXCF,MXT,MXISO,MXPSYM
      COMMON /PSIZE/  MSU,KAORDS,MCONSU,MGCSU,MSTU,MRU,MCSU,MCTU,MCONU,
     1                MCU,MSFU,MGU,MSFRU,MNSFU,MPRU,MCXU,MCCU,MBLU,MRCRU
      COMMON /OUTPUT/ NPRINT,NGEOM,NORMF,NORMP
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV
      COMMON /INOUT/  IOUT,ICHECK,ITAP30
      COMMON /SYMTRY/ NT
      COMMON /TRANSF/ XSMAL,YSMAL,ZSMAL,XNEW,YNEW,ZNEW,XP,YP,ZP
      COMMON /FRAME/  U1,U2,U3,V1,V2,V3,W1,W2,W3,X0,Y0,Z0
      COMMON /SMORB/  GROUP,INDEX,NAXIS
clj031091
      character*8 symm
      integer errcod, frdc, frdf1, prcntr, fcnt, cout
      integer inold
      common /instyl/ inold
      REAL*8 T(MXT),C(3,3,MXTRN),D(9,MXTRN)
      INTEGER INVT(MXTRN)
      DIMENSION GRP(19),DRC(3)
      DIMENSION TI(3,3),UI(9),VI(3,3),WI(9)
      EQUIVALENCE (TI(1,1),UI(1))
      EQUIVALENCE (VI(1,1),WI(1))
      DATA GRP /5HC1   ,5HCS   ,5HCI   ,5HCN   ,5HS2N  ,5HCNH  ,
     1  5HCNV  ,5HDN   ,5HDNH  ,5HDND  ,5HCINFV,5HDINFH,5HT    ,
     1  5HTH   ,5HTD   ,5HO    ,5HOH   ,5HI    ,5HIH   /
      DATA DRC /8HNORMAL  ,8HPARALLEL,8HTHREE   /
      DATA ZERO,PT5,ONE /0.0D+00,0.5D+00,1.0D+00/
      DATA TOL /1.0D-10/
cges  DATA PI2 /6.2831853071795867D+00/
      DATA PI2 /6.2831853071796D+00/
C***********************************************************************
C  THIS VALUE OF PI2 CAREFULLY SELECTED SO THAT SIN(PI2/2)=0           *
C***********************************************************************
clj031091 old style input
      if (inold.eq.1) then
        READ(IR,9999) GROUP,NAXIS
        IF(NPRINT.EQ.7) WRITE(IP,9999) GROUP,NAXIS
        DO 2 I=1,19
    2   IF(GROUP.EQ.GRP(I)) INDEX=I
      else
        errcod = frdc('SYMMETRY',symm)
        if (errcod.ne.0) then
          write(6,*) 'ERROR: couldn''t read ''SYMMETRY'' from input'
          call pabort
          call qabort
          endif
        write(6,*) 'SYMMETRY = ',symm
        if (symm.eq.'C1') then
          index = 1
          group=grp(index)
        else if (symm.eq.'CS') then
          index = 2
          group=grp(index)
        else if (symm.eq.'CI') then
          index = 3
          group=grp(index)
        else if (symm.eq.'C2') then
          index = 4
          group=grp(index)
          naxis = 2
        else if (symm.eq.'C2H') then
          index = 6
          group=grp(index)
          naxis = 2
        else if (symm.eq.'C2V') then
          index = 7
          group=grp(index)
          naxis = 2
        else if (symm.eq.'D2') then
          index = 8
          group=grp(index)
          naxis = 2
        else if (symm.eq.'D2H') then
          index = 9
          group=grp(index)
          naxis = 2
        else
          write(6,*) 'ERROR: bad value for ''symmetry'''
          call qabort
          endif
        endif
      if (prcntr('IS_ON BRIEF').eq.0) WRITE(IW,9991) GRP(INDEX),NAXIS
clj031091 read the new style input
c get ORIGIN
      if (inold.eq.0) then
        errcod = fcnt('ORIGIN',cout)
        if (errcod.ne.0) then
          x0 = 0.0D0
          y0 = 0.0D0
          z0 = 0.0D0
        else if (cout .ne. 3) then
          write(6,*) 'ERROR: wrong number of elements in ''ORIGIN'''
          call qabort
        else
          errcod = frdf1('ORIGIN',X0,1)
          if (errcod.eq.0) errcod = frdf1('ORIGIN',Y0,2)
          if (errcod.eq.0) errcod = frdf1('ORIGIN',Z0,3)
          if (errcod.ne.0) then
            write(6,*) 'ERROR: problem reading ''ORIGIN'''
            call pabort
            call qabort
            endif
          endif
c get Z_AXIS
        errcod = fcnt('Z_AXIS',cout)
        if (errcod.ne.0) then
          x1 = 0.0D0
          y1 = 0.0D0
          z1 = 1.0D0
        else if (cout .ne. 3) then
          write(6,*) 'ERROR: wrong number of elements in ''Z_AXIS'''
          call qabort
        else
          errcod = frdf1('Z_AXIS',X1,1)
          if (errcod.eq.0) errcod = frdf1('Z_AXIS',Y1,2)
          if (errcod.eq.0) errcod = frdf1('Z_AXIS',Z1,3)
          if (errcod.ne.0) then
            write(6,*) 'ERROR: problem reading ''Z_AXIS'''
            call pabort
            call qabort
            endif
          endif
c get X_AXIS
        errcod = fcnt('X_AXIS',cout)
        if (errcod.ne.0) then
          x2 = 1.0D0
          y2 = 0.0D0
          z2 = 0.0D0
        else if (cout .ne. 3) then
          write(6,*) 'ERROR: wrong number of elements in ''X_AXIS'''
          call qabort
        else
          errcod = frdf1('X_AXIS',X2,1)
          if (errcod.eq.0) errcod = frdf1('X_AXIS',Y2,2)
          if (errcod.eq.0) errcod = frdf1('X_AXIS',Z2,3)
          if (errcod.ne.0) then
            write(6,*) 'ERROR: problem reading ''X_AXIS'''
            call pabort
            call qabort
            endif
          endif
c Write out the coordinates
        write(6,9001) ' ORIGIN = (',x0,y0,z0,')'
        write(6,9001) ' Z_AXIS = (',x1,y1,z1,')'
        write(6,9001) ' X_AXIS = (',x2,y2,z2,')'
 9001   format(1X,A,3(1X,F16.10),A)

c --------------------- Double check the input
c Make sure that the x and z axis are perpendicular
        dprod = x1*x2 + y1*y2 + z1*z2
        if (dprod.gt.tol) then
          write(6,*) 'ERROR: the x and z axes are not perpendicular'
          call qabort
          endif
c Make sure that the x and z axis have a reasonable magnitude.
        if (x1*x1+y1*y1+z1*z1 .lt. 0.01d0) then
          write(6,*) 'ERROR: the magnitude of ''Z_AXIS'' is too small'
          call qabort
          endif
        if (x2*x2+y2*y2+z2*z2 .lt. 0.01d0) then
          write(6,*) 'ERROR: the magnitude of ''X_AXIS'' is too small'
          call qabort
          endif

c ---------------------- Compute the Y axis (point 3)
c If we are using Cs symmetry, then points 0,1,2 must be transformed into
c points on the sigma h plane--the xy plane.
c This will be done by rotating z(point1) about x(point2) by 90 degrees.
        rnorm = 1.D0/dsqrt((x2*x2+y2*y2+z2*z2)*(x1*x1+y1*y1+z1*z1))

        x3 = - (y1*z2 - z1*y2)*rnorm
        y3 = - (x1*z2 - z1*x2)*rnorm
        z3 = - (x1*y2 - y1*x2)*rnorm

        write(6,9001) '(Y_AXIS = (',x3,y3,z3,'))'

c ------------------------ Switch to absolute coordinates
c point one was given relative to the origin, put in absolute coordinates
        x1 = x1 + x0
        y1 = y1 + y0
        z1 = z1 + z0
c point two was given relative to the origin, put in absolute coordinates
        x2 = x2 + x0
        y2 = y2 + y0
        z2 = z2 + z0
c point three was computed relative to the origin, put in absolute coordinates
        x3 = x3 + x0
        y3 = y3 + y0
        z3 = z3 + z0

c ----------------------- Special cases
c For CS symmetry make point 1 the origin
        if (symm.eq.'CS') then
          x1 = x0
          y1 = y0
          z1 = z0
          endif
c ---------------------- Done reading new style symmetry input.
        endif
      IF(INDEX.LE.3) GO TO 200
C***********************************************************************
C    CHECK IF SYMMETRY SPECIFIED IS D2H OR A SUBGROUP.                 *
C***********************************************************************
      IF(INDEX.EQ.11.OR.INDEX.EQ.12) THEN
         write(6,*) 'ERROR: USE DNH OR CNV INSTEAD OF DINFH OR CINFV '
         call qabort
      ENDIF
      IF(INDEX.GT.12) THEN
         write(6,*) 'ERROR: CAN ONLY RUN D2H AND LOWER'
         call qabort
      ENDIF
C***********************************************************************
C     DEFINE LOCAL FRAME                                               *
C     READ IN PRINCIPAL AXIS   ( 1 CARD )                              *
C     READ IN X-LOCAL AXIS   ( 1 CARD )                                *
C     DEFAULT OPTION : LOCAL FRAME IDENTICAL TO MASTER FRAME           *
C***********************************************************************
clj031091 an old style read
      if (inold.eq.1) READ(IR,9997) X0,Y0,Z0,X1,Y1,Z1
      IF(NPRINT.EQ.7) WRITE(IP,9997) X0,Y0,Z0,X1,Y1,Z1
      RHO=DSQRT((X1-X0)**2+(Y1-Y0)**2+(Z1-Z0)**2)
      IF(RHO.GT.TOL) GO TO 6
      X0=ZERO
      Y0=ZERO
      Z0=ZERO
      X1=ZERO
      Y1=ZERO
      Y2=ZERO
      Z2=ZERO
      Z1=ONE
      X2=ONE
      DIRECT=DRC(2)
      RHO=ONE
      GO TO 7
    6 continue
      if (inold.eq.1) READ(IR,9996) X2,Y2,Z2,DIRECT
        IF(NPRINT.EQ.7) WRITE(IP,9996) X2,Y2,Z2,DIRECT
    7 CONTINUE
      W1=(X1-X0)/RHO
      W2=(Y1-Y0)/RHO
      W3=(Z1-Z0)/RHO
      WW=W1*W1+W2*W2+W3*W3
      X02=X2-X0
      Y02=Y2-Y0
      Z02=Z2-Z0
      RHO=(W1*X02+W2*Y02+W3*Z02)/WW
      DUM=RHO*W1
      X0=X0+DUM
      X02=X02-DUM
      DUM=RHO*W2
      Y0=Y0+DUM
      Y02=Y02-DUM
      DUM=RHO*W3
      Z0=Z0+DUM
      Z02=Z02-DUM
      UU=(X02*X02+Y02*Y02+Z02*Z02)
      U=DSQRT(UU)
      U1=X02/U
      U2=Y02/U
      U3=Z02/U
      V3=W1*U2-W2*U1
      V2=W3*U1-W1*U3
      V1=W2*U3-W3*U2
      IF(DIRECT.NE.DRC(1)) GO TO 8
      DUM=U1
      U1=V1
      V1=-DUM
      DUM=U2
      U2=V2
      V2=-DUM
      DUM=U3
      U3=V3
      V3=-DUM
    8 CONTINUE
      IF(NPRINT.NE.1) GO TO 9
      WRITE(IW,9990) X0,Y0,Z0,U1,V1,W1,U2,V2,W2,U3,V3,W3
    9 CONTINUE
      IF(INDEX.GE.13) GO TO 200
C***********************************************************************
C     ROTATION ABOUT PRINCIPAL AXIS                                    *
C***********************************************************************
      NN=0
      N=NAXIS
      ALPHA=ZERO
      ALPH=PI2/DBLE(N)
   10 NN=NN+1
      IF(NN.GT.N) GO TO 20
      COSA=DCOS(ALPHA)
      SINA=DSIN(ALPHA)
      I=9*(NN-1)
      T(I+1)=COSA
      T(I+5)=COSA
      T(I+2)=-SINA
      T(I+4)=SINA
      T(I+3)=ZERO
      T(I+6)=ZERO
      T(I+7)=ZERO
      T(I+8)=ZERO
      T(I+9)=ONE
      ALPHA=ALPHA+ALPH
      GO TO 10
C
C     END OF GROUP 4
C
   20 NT=N
      II=9*NT
      IF(NPRINT.NE.1) GO TO 24
      WRITE(IW,9989)
      N1=1
      N2=NAXIS
      CALL PRINT2(N1,N2,T)
   24 CONTINUE
      IF(INDEX.EQ.4) GO TO 1000
      IF(INDEX.EQ.5) GO TO 500
      IF(INDEX.EQ.7) GO TO 115
      IF(INDEX.NE.6.AND.INDEX.NE.9) GO TO 55
C
C     SIGMA-H PLANE  EQUATION (Z=0) IN LOCAL FRAME
C
      NN=0
   30 NN=NN+1
      IF(NN.GT.NT) GO TO 50
C
C     GROUP 6 0R 9
C
      I=II+9*(NN-1)
      DO 40 J=1,8
   40 T(I+J)=T(I+J-II)
      T(I+9)=-T(I+9-II)
      GO TO 30
   50 NT=NT+NT
      II=9*NT
      IF(NPRINT.NE.1) GO TO 54
      WRITE(IW,9988)
      N1=N2+1
      N2=N2+NAXIS
      CALL PRINT2(N1,N2,T)
   54 CONTINUE
C
C     END OF GROUP 6
C
      IF(INDEX.EQ.6) GO TO 1000
C
C     ONE CP2 AXIS IS THE X-AXIS OF THE LOCAL FRAME
C     GROUP 8 , 9 ,10
C
   55 CONTINUE
      NN=0
   60 NN=NN+1
      IF(NN.GT.NT) GO TO 70
      I=II+9*(NN-1)
      T(I+1)=T(I+1-II)
      T(I+2)=-T(I+2-II)
      T(I+3)=-T(I+3-II)
      T(I+4)=T(I+4-II)
      T(I+5)=-T(I+5-II)
      T(I+6)=-T(I+6-II)
      T(I+7)=T(I+7-II)
      T(I+8)=-T(I+8-II)
      T(I+9)=-T(I+9-II)
      GO TO 60
   70 NT=NT+NT
      II=9*NT
      IF(NPRINT.NE.1) GO TO 99
      WRITE(IW,9987)
      N1=N2+1
      N2=N2+NAXIS
      CALL PRINT2(N1,N2,T)
      IF(INDEX.NE.9) GO TO 99
      WRITE(IW,9981)
      N1=N2+1
      N2=N2+NAXIS
      CALL PRINT2(N1,N2,T)
   99 CONTINUE
C
C     END OF GROUP 8 AND 9
C
      IF(INDEX.EQ.8.OR.INDEX.EQ.9) GO TO 1000
C
C     DND GROUP . EQUATION OF PLANE SIGMA-D IS :
C     DSIN(ALPH/4)*X-DCOS(ALPH/4)*Y=0
C     THE X-AXIS IS THE CP2 AXIS.
C
C     GROUP 10
C
      BETA=PT5*ALPH
      COSA=DCOS(BETA)
      SINA=DSIN(BETA)
      NN=0
  100 NN=NN+1
      IF(NN.GT.NT) GO TO 110
      I=II+9*(NN-1)
      T(I+1)=COSA*T(I+1-II) + SINA*T(I+2-II)
      T(I+2)=SINA*T(I+1-II) - COSA*T(I+2-II)
      T(I+3)=     T(I+3-II)
      T(I+4)=COSA*T(I+4-II) + SINA*T(I+5-II)
      T(I+5)=SINA*T(I+4-II) - COSA*T(I+5-II)
      T(I+6)=     T(I+6-II)
      T(I+7)=COSA*T(I+7-II) + SINA*T(I+8-II)
      T(I+8)=SINA*T(I+7-II) - COSA*T(I+8-II)
      T(I+9)=     T(I+9-II)
      GO TO 100
  110 NT=NT+NT
      II=9*NT
      IF(NPRINT.NE.1) GO TO 114
      WRITE(IW,9986)
      N1=N2+1
      N2=N2+NAXIS
      CALL PRINT2(N1,N2,T)
      WRITE(IW,9980)
      N1=N2+1
      N2=N2+NAXIS
      CALL PRINT2(N1,N2,T)
  114 CONTINUE
C
C     END OF GROUP 10
C
      GO TO 1000
C
C     GROUP 7
C     SIGMA-V IS THE (X-Z) PLANE OF LOCAL FRAME
C
  115 CONTINUE
      NN=0
  120 NN=NN+1
      IF(NN.GT.NT) GO TO 130
      I=II+9*(NN-1)
      T(I+1)=T(I+1-II)
      T(I+2)=-T(I+2-II)
      T(I+3)=T(I+3-II)
      T(I+4)=T(I+4-II)
      T(I+5)=-T(I+5-II)
      T(I+6)=T(I+6-II)
      T(I+7)=T(I+7-II)
      T(I+8)=-T(I+8-II)
      T(I+9)=T(I+9-II)
      GO TO 120
  130 NT=NT+NT
      II=9*NT
C
C     END OF GROUP 7
C
      IF(NPRINT.NE.1) GO TO 1000
      WRITE(IW,9985)
      N1=N2+1
      N2=N2+NAXIS
      CALL PRINT2(N1,N2,T)
      GO TO 1000
  200 CONTINUE
      T(1)=ONE
      T(5)=ONE
      T(9)=ONE
      T(2)=ZERO
      T(3)=ZERO
      T(4)=ZERO
      T(6)=ZERO
      T(7)=ZERO
      T(8)=ZERO
      IF(INDEX.EQ.1) GO TO 210
      IF(INDEX.EQ.2) GO TO 250
      IF(INDEX.EQ.3) GO TO 300
      IF(INDEX.EQ.11.OR.INDEX.EQ.12) GO TO 400
      IF(INDEX.GE.18) GO TO 800
      GO TO 600
  210 NT=1
      X0=ZERO
      Y0=ZERO
      Z0=ZERO
      U1=ONE
      V2=ONE
      W3=ONE
      U2=ZERO
      U3=ZERO
      V1=ZERO
      V3=ZERO
      W1=ZERO
      W2=ZERO
      GO TO 1000
C
C     CS SYMMETRY GROUP
C     THE 3 POINTS 1,2,3 DEFINE SIGMA-H PLANE
C
  250 continue
      if (inold.eq.1) READ(IR,9997) X1,Y1,Z1,X2,Y2,Z2
      RHO=(X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2
      IF(RHO.GT.TOL) GO TO 251
C     DEFAULT OPTION : PLANE IS THE (X,Y) PLANE
      X1=ZERO
      Y1=ZERO
      Z1=ZERO
      Y2=ZERO
      Z2=ZERO
      X3=ZERO
      Z3=ZERO
      X2=ONE
      Y3=ONE
      GO TO 252
  251 continue
      if (inold.eq.1) READ(IR,9997) X3,Y3,Z3
  252 CONTINUE
      NT=2
      W1=(Y2-Y1)*(Z3-Z1)-(Y3-Y1)*(Z2-Z1)
      W2=(Z2-Z1)*(X3-X1)-(Z3-Z1)*(X2-X1)
      W3=(X2-X1)*(Y3-Y1)-(X3-X1)*(Y2-Y1)
      RHO=DSQRT(W1*W1+W2*W2+W3*W3)
      W1=W1/RHO
      W2=W2/RHO
      W3=W3/RHO
      U1=X2-X1
      U2=Y2-Y1
      U3=Z2-Z1
      RHO=DSQRT(U1*U1+U2*U2+U3*U3)
      U1=U1/RHO
      U2=U2/RHO
      U3=U3/RHO
      V1=W2*U3-W3*U2
      V2=W3*U1-W1*U3
      V3=W1*U2-W2*U1
      X0=X1
      Y0=Y1
      Z0=Z1
      T(10)=ONE
      T(14)=ONE
      T(18)=-ONE
      T(11)=ZERO
      T(12)=ZERO
      T(13)=ZERO
      T(15)=ZERO
      T(16)=ZERO
      T(17)=ZERO
      IF(NPRINT.NE.1) GO TO 1000
      WRITE(IW,9983) W1,W2,W3
      WRITE(IW,9982) U1,V1,W1,U2,V2,W2,U3,V3,W3
      GO TO 1000
C
C     CI SYMMETRY GROUP
C     CENTER OF INVERSION IS (X0,Y0,Z0)
C
  300 continue
      if (inold.eq.1) READ(IR,9997) X0,Y0,Z0
      IF(NPRINT.EQ.3) WRITE(IW,9984) X0,Y0,Z0
      T(10)=-ONE
      T(14)=-ONE
      T(18)=-ONE
      T(11)=ZERO
      T(12)=ZERO
      T(13)=ZERO
      T(15)=ZERO
      T(16)=ZERO
      T(17)=ZERO
      NT=2
      U1=ONE
      V2=ONE
      W3=ONE
      U2=ZERO
      U3=ZERO
      V1=ZERO
      V3=ZERO
      W1=ZERO
      W2=ZERO
      GO TO 1000
  400 CONTINUE
      STOP
  500 NN=0
      BETA=PT5*ALPH
      COSB=DCOS(BETA)
      SINB=DSIN(BETA)
  510 NN=NN+1
      IF(NN.GT.NT) GO TO 520
C
C     S2N GROUP
C     THE PLANE OF SYMMETRY FOR THE IMPROPER ROTATION
C     IS THE (X,Y) PLANE OF THE LOCAL FRAME
C
      I=II+9*(NN-1)
      T(I+1)= COSB*T(I+1-II)+SINB*T(I+2-II)
      T(I+2)=-SINB*T(I+1-II)+COSB*T(I+2-II)
      T(I+3)=     -T(I+3-II)
      T(I+4)= COSB*T(I+4-II)+SINB*T(I+5-II)
      T(I+5)=-SINB*T(I+4-II)+COSB*T(I+5-II)
      T(I+6)=     -T(I+6-II)
      T(I+7)= COSB*T(I+7-II)+SINB*T(I+8-II)
      T(I+8)=-SINB*T(I+7-II)+COSB*T(I+8-II)
      T(I+9)=     -T(I+9-II)
      GO TO 510
  520 NT=NT+NT
      II=9*NT
      IF(NPRINT.NE.1) GO TO 1000
      WRITE(IW,9979)
      N1=N2+1
      N2=N2+NAXIS
      CALL PRINT2(N1,N2,T)
      GO TO 1000
C
C     T GROUP AND OTHERS CONTAINING A SUBGROUP T :
C     LOCAL X,Y,Z ARE THE C2 AXES
C
  600 DO 610 I=10,36
  610 T(I)=ZERO
      T(10)=ONE
      T(23)=ONE
      T(36)=ONE
      T(14)=-ONE
      T(18)=-ONE
      T(19)=-ONE
      T(27)=-ONE
      T(28)=-ONE
      T(32)=-ONE
      DO 620 II=5,12
        I=9*(II-1)
        J=9*(II-5)
        T(I+1)=T(J+7)
        T(I+2)=T(J+8)
        T(I+3)=T(J+9)
        T(I+4)=T(J+1)
        T(I+5)=T(J+2)
        T(I+6)=T(J+3)
        T(I+7)=T(J+4)
        T(I+8)=T(J+5)
  620 T(I+9)=T(J+6)
      NT=12
      IF(INDEX.EQ.13) GO TO 1000
      IF(INDEX.EQ.14) GO TO 650
      IF(INDEX.EQ.15) GO TO 680
      GO TO 670
C
C     TH GROUP
C     EXPAND GROUP BY TAKING DIRECT PRODUCT WITH CI
C
  650 I=9*NT
      DO 660 J=1,I
  660 T(J+I)=-T(J)
      NT=NT+NT
      GO TO 1000
C
C     OCTAHEDRAL GROUP IS DIRECT PRODUCT OF T AND A C4 ROTATION
C     ABOUT Z AXIS
C
  670 SIGN=-ONE
      GO TO 685
C
C     TD GROUP IS DIRECT PRODUCT OF T AND A REFLECTION IN A
C     PLANE ( EQUATION OF THE PLANE   X=Y  )
C
  680 SIGN=ONE
  685 DO 690 II=13,24
        I=9*(II-1)
        J=9*(II-13)
        T(I+1)=T(J+4)*SIGN
        T(I+2)=T(J+5)*SIGN
        T(I+3)=T(J+6)*SIGN
        T(I+4)=T(J+1)
        T(I+5)=T(J+2)
        T(I+6)=T(J+3)
        T(I+7)=T(J+7)
        T(I+8)=T(J+8)
  690 T(I+9)=T(J+9)
      NT=24
      IF(INDEX.NE.17) GO TO 1000
C
C     OH GROUP IS DIRECT PRODUCT OF O AND CI
C
      I=9*NT
      DO 700 J=1,I
  700 T(J+I)=-T(J)
      NT=48
      GO TO 1000
  800 CONTINUE
C
C     ICOSAHEDRAL GROUPS I AND IH
C
      DO 801 I=1,27
  801 T(I)=ZERO
      NT=3
      S5=DSQRT(5.0D+00)
      S3=DSQRT(3.0D+00)
      SIN3=S3/2.0
      COS3=-1./2.
      COS5=DCOS(1.2566370614359D+00)
      SIN5=DSIN(1.2566370614359D+00)
      IF(DIRECT.NE.DRC(3))GO TO 880
      SINA=SIN5
      COSA=COS5
      SIN5=SIN3
      COS5=COS3
      SIN3=SINA
      COS3=COSA
  880 A=0.6070619982066862231
      B=DSQRT(1.0-A*A)
      D(1,1)=ONE
      D(5,1)=ONE
      D(9,1)=ONE
      D(1,2)=COS5
      D(2,2)=-SIN5
      D(4,2)=-D(2,2)
      D(5,2)=D(1,2)
      D(9,2)=ONE
      D(1,3)=A*A*(1.0-COS3)+COS3
      D(2,3)=B*SIN3
      D(3,3)=A*B*(1.0-COS3)
      D(4,3)=-D(2,3)
      D(5,3)=COS3
      D(6,3)=A*SIN3
      D(7,3)=D(3,3)
      D(8,3)=-D(6,3)
      D(9,3)=B*B*(1.0-COS3)+COS3
      KMIN=2
  802 IN=2
      NTOLD=NT
  803 DO 804 I=1,9
  804 UI(I)=D(I,IN)
      KN=KMIN
  805 DO 807 I=1,3
        DO 807 K=1,3
          SUM=0.0
          DO 806 J=1,3
  806     SUM=SUM+C(I,J,KN)*TI(J,K)
  807 VI(I,K)=SUM
      J1=1
  808 DO 809 I=1,9
        Q=DABS(D(I,J1)-WI(I))
        IF(Q.GT.TOL) GO TO  810
  809 CONTINUE
      GO TO 813
  810 IF(J1.EQ.NT) GO TO 811
      J1=J1+1
      GO TO 808
  811 NT=NT+1
      DO 812 I=1,9
  812 D(I,NT)=WI(I)
  813 KN=KN+1
      IF(KN.LE.NT)GO TO 805
      IN=IN+1
      IF(IN.LE.3)GO TO 803
      IF(NT.EQ.NTOLD) GO TO 814
      KMIN=NTOLD+1
      GO TO 802
  814 CONTINUE
      IF(INDEX.NE.19) GO TO 850
C
C     IH GROUP IS DIRECT PRODUCT OF I AND CI
C
      DO 815 I=1,540
  815 T(540+I)=-T(I)
      NT=120
  850 IF(NPRINT.NE.1) GO TO 1000
      N1=1
      N2=NT
      CALL PRINT2(N1,N2,T)
C
C     END OF ICOSAHEDRAL GROUPS
C
 1000 CONTINUE
C
C     FIND THE INVERSE TRANSFORMATIONS
C
      DO 1002 ITR=1,NT
        NN=9*(ITR-1)
        DO 1001 IT=1,NT
          II=9*(IT-1)
          TEST=T(NN+1)*T(II+1)+T(NN+2)*T(II+4)+T(NN+3)*T(II+7)
     1        +T(NN+4)*T(II+2)+T(NN+5)*T(II+5)+T(NN+6)*T(II+8)
     2        +T(NN+7)*T(II+3)+T(NN+8)*T(II+6)+T(NN+9)*T(II+9)-3.0D0
          IF(DABS(TEST).GT.TOL) GO TO 1001
          INVT(ITR)=IT
          GO TO 1002
 1001   CONTINUE
 1002 CONTINUE
C
      MGU=NT
C
      RETURN
 9979 FORMAT(' S2N ROTATION FOLLOWED BY ROTATIONS')
 9980 FORMAT(' SIGMA-D FOLLOWED BY C2 FOLLOWED BY ROTATIONS')
 9981 FORMAT(' C2 FOLLOWED BY SIGMA-H FOLLOWED BY ROTATIONS')
 9982 FORMAT(/,10X,3F15.9,/,10X,3F15.9,/,10X,3F15.9)
 9983 FORMAT(/,' PLANE OF SYMMETRY DEFINED BY ITS NORMAL U = ',F10.5,
     1' V = ',F10.5,' W = ',F10.5)
 9984 FORMAT(/,10X,' CENTER OF SYMMETRY AT X = ',F10.5,' Y = ',F10.5,
     1' Z = ',F10.5)
 9985 FORMAT(' SIGMA-V FOLLOWED BY ROTATIONS')
 9986 FORMAT(' SIGMA-D FOLLOWED BY ROTATIONS')
 9987 FORMAT(' C2 FOLLOWED BY ROTATIONS ')
 9988 FORMAT(' SIGMA-H FOLLOWED BY ROTATIONS')
 9989 FORMAT(' ROTATIONS ABOUT PRINCIPAL AXIS')
 9990 FORMAT(/,' THE ORIGIN OF THE LOCAL FRAME IS AT X =  ',F10.5,
     1 ' Y = ',F10.5,' Z = ',F10.5,/,
     1 ' DIRECTOR COSINES OF THE NEW AXES ',/,34X,3(5X,F10.5),/,34X,
     1 3(5X,F10.5),/,34X,3(5X,F10.5))
 9991 FORMAT(/,' THE POINT GROUP OF THE MOLECULE IS ...',A8,/,
     1 ' THE ORDER OF THE PRINCIPAL AXIS IS ...',I5,/)
 9996 FORMAT(3F10.5,A8)
 9997 FORMAT(9F10.5)
 9999 FORMAT (A5,I5)
      END
