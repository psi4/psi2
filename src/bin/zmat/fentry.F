      SUBROUTINE FENTRY(CC,IA,MAXCOR)
c**********************************************************************
c Moved to PSI distribution on disk 020289 - clj                      *
c**********************************************************************
C*   MODIFICATION FOR IMS VERSION                                     *
C*   BY: YUKIO YAMAGUCHI                                              *
C*   DATE: FEBRUARY 21, 1989                                          *
c**********************************************************************
c  Write coordinates to INPUT file for use with GEOMUP$ - bfy 2/17/88 *
c**********************************************************************
C*** LAST UPDATED ON SEPTEMBER 11, 1986 BY R. REMINGTON               *
C* REASON: TO OBTAIN HIGHER ACCURACY BY INCREASING THE                *
C* NUMBER OF DIGITS FOR PI AND CHANGE ATAN TO DATAN #143              *
c**********************************************************************
C***LAST UPDATED ON FEBRUARY 01, 1985 BY YUKIO YAMAGUCHI              *
c**********************************************************************
C   THIS PROGRAM CALCULATES CARTESIAN COORDINATES OF A SYSTEM
C   SPECIFIED BY INTERNAL COORDINATES
C   THE PROGRAM IS ABLE TO HANDLE UP TO 150 INTERNAL COORDINATES
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CC(MAXCOR),IA(1)
      COMMON/ZMATS/INUM(150),IZ(150,4),BL(150),ALP(150),BET(150),
     1             IANZ(150),NZ,input
      DATA A00 / 0.0D+00 /
    1 FORMAT(//,2X,' PROGRAM ZMAT'/)
C
      call drum
      INPUT=5
      iout = 6
      CALL LOCATE(INPUT,'# ZMAT ###',IERR)
      if(ierr.ne.0) then
          write(iout,*) 'No ZMAT input ?'
          stop
      endif
C
      WRITE(6,1)
      READ(5,*) INUM(1),IANZ(1)
      IZ(1,1)=0
      IZ(1,2)=0
      IZ(1,3)=0
      BL(1)=A00
      ALP(1)=A00
      BET(1)=A00
      I=1
  200 I=I+1
      IZ(I,1)=0
      IZ(I,2)=0
      IZ(I,3)=0
      BL(I)=A00
      ALP(I)=A00
      BET(I)=A00
      IF(I.EQ.2) GO TO 201
      IF(I.EQ.3) GO TO 202
      GO TO 203
  201 READ(5,*) INUM(2),IZ(2,1),BL(2),IANZ(2)
      GO TO 204
  202 READ(5,*) INUM(3),IZ(3,1),BL(3),IZ(3,2),ALP(3),IANZ(3)
      GO TO 204
  203 READ(5,*) INUM(I),IZ(I,1),BL(I),IZ(I,2),ALP(I),IZ(I,3),BET(I),
     1          IZ(I,4),IANZ(I)
  204 IF(INUM(I)) 200,205,200
C
  205 NZ=I-1
C
      CALL BUILDZ
C
C2-20-89 STOP
      RETURN
      END
      SUBROUTINE BUILDZ
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/ZMATS/INUM(150),IZ(150,4),BL(150),ALP(150),BET(150),
     1             IANZ(150),NZ,input
      COMMON/COM101/NATOMS
      COMMON/COM201/C(150,3)
      DIMENSION A(150),B(150),CZ(150,3),D(150)
      DIMENSION U1(3),U2(3),U3(3),U4(3),VJ(3),VP(3),V3(3)
      DATA PI / 3.141592653589793D+00 /
      DATA A00,ONE,TWO / 0.0D+00 , 1.0D+00 , 2.0D+00 /
    1 FORMAT(/,2X,' INTERNAL COORDINATES'/)
    2 FORMAT(2X,I3,84X,I5)
    3 FORMAT(2X,I3,I5,F14.7,'(',I3,')',60X,I5)
    4 FORMAT(2X,I3,I5,F14.7,'(',I3,')',I5,F14.7,'(',I3,')',36X,I5)
    5 FORMAT(2X,I3,I5,F14.7,'(',I3,')',I5,F14.7,'(',I3,')',I5,F14.7,
     1 '(',I3,')',I5,7X,I5)
    6 FORMAT(//,25X,' CARTESIAN COORDINATES in Angstrom',
     $      /,14X,'X',24X,'Y',24X,'Z'/)
    7 FORMAT(2X,I3,3(F20.10,5X))
    8 FORMAT(//,25X,' CARTESIAN COORDINATES in bohr',
     $      /,14X,'X',24X,'Y',24X,'Z'/)
    9 format(3f20.10)
C
      AR=PI/180.0D+00
      ANGTOB = 0.52917706D0
      WRITE(6,1)
      WRITE(6,2) INUM(1),IANZ(1)
      IF(NZ.LE.1) GO TO 201
      I1=1
      WRITE(6,3) INUM(2),IZ(2,1),BL(2),I1,IANZ(2)
      IF(NZ.LE.2) GO TO 201
      I1=I1+1
      I2=NZ
      WRITE(6,4) INUM(3),IZ(3,1),BL(3),I1,IZ(3,2),ALP(3),I2,IANZ(3)
      IF(NZ.LE.3) GO TO 201
      I3=2*NZ-3
      DO 101 I=4,NZ
      I1=I1+1
      I2=I2+1
      I3=I3+1
      WRITE(6,5) INUM(I),IZ(I,1),BL(I),I1,IZ(I,2),ALP(I),I2,IZ(I,3),
     1 BET(I),I3,IZ(I,4),IANZ(I)
  101 CONTINUE
  201 DO 102 I=1,NZ
      DO 102 J=1,3
  102 CZ(I,J)=A00
      DO 20 I=1,NZ
      ALP(I)=ALP(I)*AR
   20 BET(I)=BET(I)*AR
      CZ(2,3)=BL(2)
      IF(NZ-3) 260,30,30
   30 CZ(3,1)=BL(3)*DSIN(ALP(3))
      IF(IZ(3,1)-1) 50,40,50
   40 CZ(3,3)=BL(3)*DCOS(ALP(3))
      GO TO 60
   50 CZ(3,3)=CZ(2,3)-BL(3)*DCOS(ALP(3))
   60 DO 80 I=4,NZ
      IF(DABS(CZ(I-1,1))-1.0D-6) 70,90,90
   70 CZ(I,1)=BL(I)*DSIN(ALP(I))
      ITEMP=IZ(I,1)
      JTEMP=IZ(I,2)
   80 CZ(I,3)=CZ(ITEMP,3)-BL(I)*DCOS(ALP(I))*DSIGN(ONE,CZ(ITEMP,3)
     1 -CZ(JTEMP,3))
   90 K=I
      IF(K-NZ) 100,100,260
  100 DO 250 J=K,NZ
      CAJ=DCOS(ALP(J))
      SAJ=DSIN(ALP(J))
      CBJ=DCOS(BET(J))
      SBJ=DSIN(BET(J))
      IF(IZ(J,4)) 135,110,135
  110 CALL VEC(U1,CZ,IZ(J,2),IZ(J,3))
      CALL VEC(U2,CZ,IZ(J,1),IZ(J,2))
      CALL VPROD(VP,U1,U2)
      R=DSQRT(ONE-(U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3))**2)
      DO 120 I=1,3
  120 U3(I)=VP(I)/R
      CALL VPROD(U4,U3,U2)
      DO 130 I=1,3
      VJ(I)=BL(J)*(-U2(I)*CAJ+U4(I)*SAJ*CBJ+U3(I)*SAJ*SBJ)
      ITEMP=IZ(J,1)
  130 CZ(J,I)=VJ(I)+CZ(ITEMP,I)
      GO TO 250
  135 IF(IABS(IZ(J,4)-1)) 210,140,210
  140 CALL VEC(U1,CZ,IZ(J,1),IZ(J,3))
      CALL VEC(U2,CZ,IZ(J,2),IZ(J,1))
      ZETA=-(U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3))
      A(J)=(-CBJ+ZETA*CAJ)/(ONE-ZETA*ZETA)
      B(J)=(CAJ-ZETA*CBJ)/(ONE-ZETA*ZETA)
      R=A00
      GAMMA=PI/TWO
      IF(ZETA) 150,170,160
  150 R=PI
  160 GAMMA=DATAN(DSQRT(ONE-ZETA*ZETA)/ZETA)+R
  170 D(J)=A00
      IF(DABS(GAMMA+ALP(J)+BET(J)-TWO*PI)-1.0D-6) 190,180,180
  180 D(J)=IZ(J,4)*(DSQRT(ONE+A(J)*CBJ-B(J)*CAJ))/DSQRT(ONE-ZETA*ZETA)
  190 CALL VPROD(V3,U1,U2)
      DO 200 I=1,3
      U3(I)=A(J)*U1(I)+B(J)*U2(I)+D(J)*V3(I)
      VJ(I)=BL(J)*U3(I)
      ITEMP=IZ(J,1)
  200 CZ(J,I)=VJ(I)+CZ(ITEMP,I)
      GO TO 250
  210 CALL VEC(U1,CZ,IZ(J,1),IZ(J,3))
      CALL VEC(U2,CZ,IZ(J,2),IZ(J,1))
      ZETA=-(U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3))
      CALL VPROD(V3,U1,U2)
      V3MAG=DSQRT(V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3))
      A(J)=V3MAG*CBJ/(ONE-ZETA*ZETA)
      B(J)=DSQRT((ONE-CAJ*CAJ-A(J)*CBJ*V3MAG)/(ONE-ZETA*ZETA))
      IF(IZ(J,4)-2) 220,230,220
  220 B(J)=-B(J)
  230 D(J)=B(J)*ZETA+CAJ
      DO 240 I=1,3
      U3(I)=B(J)*U1(I)+D(J)*U2(I)+A(J)*V3(I)
      VJ(I)=BL(J)*U3(I)
      ITEMP=IZ(J,1)
  240 CZ(J,I)=VJ(I)+CZ(ITEMP,I)
  250 CONTINUE
C
C   ELIMINATE IMAGINARY ATOM(S)
  260 NATOMS=0
      DO 290 I=1,NZ
      IF(IANZ(I)) 270,290,270
  270 NATOMS=NATOMS+1
      DO 280 J=1,3
      C(NATOMS,J)=CZ(I,J)
  280 CONTINUE
  290 CONTINUE
      WRITE(6,6)
      DO 300 I=1,NATOMS
  300 WRITE(6,7) I,(C(I,J),J=1,3)
cbfy
c  Write cartesian coordinates in bohr to File6 and to INPUT file.
      call locate (input,'# GEOML ##',ierr)
      backspace input
      write(input,'(A10)') '# GEOMUP #'
      WRITE(6,8)
      DO 305 I=1,NATOMS
      WRITE(6,7) I,(C(I,J)/ANGTOB,J=1,3)
      WRITE(input,9) (C(I,J)/ANGTOB,J=1,3)
  305 continue
      write(input,'(A10)') '# GEOML ##'
cbfy
C
C   CALCULATE INTERATOMIC DISTANCES
      CALL DIST
C
C   CALCULATE ANGLES
      CALL ANGLE
C
      RETURN
      END
      SUBROUTINE DIST
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/COM101/NATOM
      COMMON/COM201/X(150,3)
      COMMON/COM202/R(1830)
      DATA A00/ 0.0D+00 /
    1 FORMAT(//,2X,' INTERATOMIC DISTANCE MATRIX'/)
C
      IJ=0
      DO 101 I=1,NATOM
      DO 101 J=1,I
      IJ=IJ+1
      R(IJ)=A00
      IF(I.EQ.J) GO TO 101
      R(IJ)=DSQRT((X(I,1)-X(J,1))**2+(X(I,2)-X(J,2))**2+(X(I,3)-X(J,3))
     1 **2)
  101 CONTINUE
      WRITE(6,1)
      CALL PRINT(R,1830,NATOM,6)
  201 RETURN
      END
      SUBROUTINE VEC(U,C,J,K)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(150,3),R(3),U(3)
      DATA A00 / 0.0D+00 /
C
      R2=A00
      DO 101 I=1,3
      R(I)=C(J,I)-C(K,I)
      R2=R2+R(I)*R(I)
  101 CONTINUE
      R2=DSQRT(R2)
      DO 102 I=1,3
      U(I)=R(I)/R2
  102 CONTINUE
      RETURN
      END
      SUBROUTINE VPROD(VP,X,Y)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION VP(3),X(3),Y(3)
C
      VP(1)=X(2)*Y(3)-X(3)*Y(2)
      VP(2)=X(3)*Y(1)-X(1)*Y(3)
      VP(3)=X(1)*Y(2)-X(2)*Y(1)
      RETURN
      END
      SUBROUTINE ANGLE
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/COM101/NATOM
      COMMON/COM201/X(150,3)
      COMMON/COM202/R(1830)
      DATA PI / 3.141592653589793D+00 /
CRBR  DATA PI / 3.1415926536D+00 /   OLD VALUE
    1 FORMAT(//,2X,' BOND ANGLES'/
     1 2X,' I   J   K',3X,' DEGREES   ',1X,' J   I   K',3X,' DEGREES  ',
     2 1X,' I   K   J',3X,' DEGREES   '/)
    2 FORMAT(3I4,F11.3, 2X,3I4,F11.3,1X,3I4,F11.3)
C
      AR=PI/180.0D+00
      WRITE(6,1)
      DO 101 I=1,NATOM-2
      DO 101 J=I+1,NATOM-1
      IJ=I+J*(J-1)/2
      XIJ=X(I,1)-X(J,1)
      YIJ=X(I,2)-X(J,2)
      ZIJ=X(I,3)-X(J,3)
      EXIJ=-XIJ/R(IJ)
      EYIJ=-YIJ/R(IJ)
      EZIJ=-ZIJ/R(IJ)
      EXJI=-EXIJ
      EYJI=-EYIJ
      EZJI=-EZIJ
      DO 101 K=J+1,NATOM
      IK=I+K*(K-1)/2
      XIK=X(I,1)-X(K,1)
      YIK=X(I,2)-X(K,2)
      ZIK=X(I,3)-X(K,3)
      EXIK=-XIK/R(IK)
      EYIK=-YIK/R(IK)
      EZIK=-ZIK/R(IK)
      EXKI=-EXIK
      EYKI=-EYIK
      EZKI=-EZIK
      JK=J+K*(K-1)/2
      XJK=X(J,1)-X(K,1)
      YJK=X(J,2)-X(K,2)
      ZJK=X(J,3)-X(K,3)
      EXJK=-XJK/R(JK)
      EYJK=-YJK/R(JK)
      EZJK=-ZJK/R(JK)
      EXKJ=-EXJK
      EYKJ=-EYJK
      EZKJ=-EZJK
      AIJK=EXJI*EXJK+EYJI*EYJK+EZJI*EZJK
      AJIK=EXIJ*EXIK+EYIJ*EYIK+EZIJ*EZIK
      AIKJ=EXKI*EXKJ+EYKI*EYKJ+EZKI*EZKJ
      ANGIJK=DACOS(AIJK)/AR
      ANGJIK=DACOS(AJIK)/AR
      ANGIKJ=DACOS(AIKJ)/AR
      WRITE(6,2) I,J,K,ANGIJK,J,I,K,ANGJIK,I,K,J,ANGIKJ
  101 CONTINUE
C
      RETURN
      END
