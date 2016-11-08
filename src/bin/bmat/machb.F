      subroutine machb(bmat,stype,satom,symms,symmc,nsim,nsym,nstre,
     1 nbnd,ntor,nout,nlin1,nlin2,ninv,nad3,xa,ya,za,qq,napdm,iout,
     2 wri,qonly)
C     ********************************
C     TAKEN FROM PULAY'S PROGRAM TEXAS
c     modified by edward seidl may 5, 1991
C     ********************************
      IMPLICIT REAL*8 (A-H,O-Z)
      integer satom(nsim,4),symms(nsym,0:100)
      real*8 symmc(nsym,100),bmat(nad3,nsym),qq(nsym)
      real*8 xa(napdm),ya(napdm),za(napdm)
      character*4 stype(nsim),typ,tipus
      integer nsim,nsym,nstre,nbnd,ntor,nout,nlin1,nlin2,ninv
      integer nad3,napdm,iout
      integer simco
      DOUBLE PRECISION ARCOS
      LOGICAL WRI,QONLY,SFIRST
C
C     nsym is the number of internal coordinates
C     nad3 is 3* the number of atoms (including dummy atoms)
C     bmat is the transpose of the b matrix, dimensioned as
C       b(nad3,nsym) in the calling program
C     xa,ya,and za are the x,y,z coordinates of the atoms in angstroms
C     qonly is true if only internal coordinates are desired,
C       not the b matrix
C     the internal coordinate types can be stretch, bend, torsion, out-
C     of plane bending, linear1, linear2, and inverse bond length
c     the order of the atoms for bend, lin1, and lin2 have been changed
c     so that the center atom is in the second position instead of 3rd
c     the torsion coordinate range is -pi/2 to 3pi/2
c     the simple internal coordinates are printed now
C
      DIMENSION TIPUS(7)
      DIMENSION WORT(3)
      DIMENSION A(4), IA(4), U(3), V(3), W(3), Z(3), X(3), UU(3), VV(3),
     1 WW(3), ZZ(3), UV(12)
      COMMON/DISCON/SSS(200),SFIRST
      EQUIVALENCE (KA,IA(1)), (KB,IA(2)), (KC,   IA(3)), (KD,IA(4))
      EQUIVALENCE (UV(1),UU(1)), (UV(4),VV(1)), (UV(7),WW(1)), (UV(10),Z
     1Z(1))
      PARAMETER(ONE=1.0D0,ZERO=0.0D0,PI=3.14159265358979D0)
      PARAMETER(PIH=1.57079632679490D0,TPI=6.28318530717959D0)
      DATA TIPUS/'STRE','BEND',' OUT','TORS','LIN1','LIN2','INVR'/
      DATA WORT/4HK   ,4H    ,4HF   /
C     WRITE(6,*) 'START OF MACHB'
      RAD=180.0/PI
      ANULL=ONE
      O=ONE
      ISS = 0
      C1=ZERO
      IF (WRI) WRITE(IOUT,330)
      do 10 i=1,nsym
        ns=symms(i,0)
        c1=zero
        qq(i)=zero
        if(.not.qonly) then
          do 11 j=1,nad3
   11     bmat(j,i)=zero
        endif
        do 20 isim=1,ns
          c1=c1+symmc(i,isim)**2
          simco=symms(i,isim)
          c = symmc(i,isim)
          ka=satom(simco,1)
          kb=satom(simco,2)
          kc=satom(simco,3)
          kd=satom(simco,4)
          typ=stype(simco)
          do 130 k=1,7
            if (typ.eq.tipus(k)) go to 140
  130     continue
  140     go to (150,160,180,210,230,250,265), k
C
C     ..... STRETCH
C
  150     CALL VEKTOR (UU,R1,KA,KB,XA,YA,ZA)
          UU(1)=UU(1)*ANULL
          UU(2)=UU(2)*ANULL
          UU(3)=UU(3)*ANULL
          VV(1)=-UU(1)
          VV(2)=-UU(2)
          VV(3)=-UU(3)
          IA(3)=0
          IA(4)=0
          QQ(I) = QQ(I) + R1*C
          IF (WRI) THEN
            IF(QQ(I).GT.900000000.0D0) THEN
              WRITE(IOUT,381) I,TYP,IA,C,R1
            ELSE
              WRITE(IOUT,380) I,TYP,IA,C,R1
            ENDIF
          ENDIF
          GO TO 270
C
C     .....BENDING
C
  160     CALL VEKTOR (U,R1,KA,KB,XA,YA,ZA)
          CALL VEKTOR (V,R2,KC,KB,XA,YA,ZA)
          CO=SCALAR(U,V)
          SI=S2(CO)
          DO 170 L=1,3
            UU(L)=(CO*U(L)-V(L))/(SI*R1)
            WW(L)=(CO*V(L)-U(L))/(SI*R2)
            VV(L)=-UU(L)-WW(L)
  170     CONTINUE
          IA(4)=0
          QQ(I)=QQ(I)+C*ARCOS(CO)
          IF (WRI) WRITE(IOUT,380) I,TYP,IA,C,RAD*ARCOS(CO)
          GO TO 270
C
C     .....OUT OF PLANE
C
  180     CALL VEKTOR (U,R1,KA,KB,XA,YA,ZA)
          CALL VEKTOR (V,R2,KC,KB,XA,YA,ZA)
          CALL VEKTOR (W,R3,KD,KB,XA,YA,ZA)
          CALL NORMAL (V,W,Z)
          STETA=SCALAR(U,Z)
          CTETA=S2(STETA)
          CFI1=SCALAR(V,W)
          SFI1=S2(CFI1)
          CFI2=SCALAR(W,U)
          CFI3=SCALAR(V,U)
          DEN=CTETA*SFI1**2
          ST2=(CFI1*CFI2-CFI3)/(R2*DEN)
          ST3=(CFI1*CFI3-CFI2)/(R3*DEN)
          DO 190 L=1,3
            WW(L)=Z(L)*ST2
            ZZ(L)=Z(L)*ST3
  190     CONTINUE
          CALL NORMAL(Z,U,X)
          CALL NORMAL(U,X,Z)
          DO 200 L=1,3
            UU(L)=Z(L)/R1
            VV(L)=-UU(L)-WW(L)-ZZ(L)
  200     CONTINUE
          CX=-C
          IF (STETA.LT.ZERO) CX=C
          QQ(I)=QQ(I)-CX*ARCOS(CTETA)
          IF (WRI) WRITE(IOUT,380) I,TYP,IA,C,-RAD*CX*ARCOS(CTETA)/C
          GO TO 270
C
C     ..... TORSION
C
  210     CALL VEKTOR (U,R1,KA,KB,XA,YA,ZA)
          CALL VEKTOR (V,R2,KC,KB,XA,YA,ZA)
          CALL VEKTOR (W,R3,KC,KD,XA,YA,ZA)
          CALL NORMAL (U,V,Z)
          CALL NORMAL (W,V,X)
          CO=SCALAR(U,V)
          CO2=SCALAR(V,W)
          SI=S2(CO)
          SI2=S2(CO2)
          DO 220 L=1,3
            UU(L)=Z(L)/(R1*SI)
            ZZ(L)=X(L)/(R3*SI2)
            VV(L)=(R1*CO/R2-ONE)*UU(L)-R3*CO2/R2*ZZ(L)
            WW(L)=-UU(L)-VV(L)-ZZ(L)
  220     CONTINUE
          CO=SCALAR(Z,X)
          U(1)=Z(2)*X(3)-Z(3)*X(2)
          U(2)=Z(3)*X(1)-Z(1)*X(3)
          U(3)=Z(1)*X(2)-Z(2)*X(1)
          CO2=SCALAR(U,V)
          S=ARCOS(-CO)
          IF(CO2.LT.ZERO) S=-S
C     THE FOLLOWING THREE LINES PREVENT DISCONTINUITY DURING THE
C     TRANSFORMATION TO THE NEW CARTESIAN COORDINATES.  SSS CONTAINS
C     THE ORIGINAL VALUE OF THE INTERNAL TORSIONAL COORDINATE.
C     KEEP IN MIND THAT S IS MINUS TAU.
c     by bfy and tph
          ISS = ISS + 1
          IF(SFIRST) THEN
            SSS(ISS) = S
          ENDIF
          IF(SSS(ISS) .GT. PIH) S = S-TPI
  201     CONTINUE
          IF (SSS(ISS)-S.GT.6.0D0) THEN
            S=S+TPI
            GO TO 201
          ELSE IF (SSS(ISS)-S.LT.-6.0D0) THEN
            S=S-TPI
            GO TO 201
          ENDIF
          QQ(I)=QQ(I)-C*S
          IF (WRI) WRITE(IOUT,380) I,TYP,IA,C,-RAD*S
C
C     .... REMEMBER THAT THE RANGE OF THIS COORDINATE IS -PI/2 TO 3*PI/2
C     .... IN ORDER TO SHIFT THE DISCONTINUITY OFF THE PLANAR POSITION
C
      GO TO 270
C
C     .....LINEAR COPLANAR BENDING
C
  230     CALL VEKTOR (U,R1,KA,KB,XA,YA,ZA)
          CALL VEKTOR (V,R2,KD,KB,XA,YA,ZA)
          CALL VEKTOR (X,R2,KC,KB,XA,YA,ZA)
          CO=SCALAR(V,U)
          CO2=SCALAR(X,V)
          QQ(I)=QQ(I)+C*(PI-ARCOS(CO)-ARCOS(CO2))
          CALL NORMAL (V,U,W)
          CALL NORMAL (U,W,Z)
          CALL NORMAL(X,V,W)
          CALL NORMAL (W,X,U)
C
C     ..... COORDINATE POSITIVE IF ATOM A MOVES TOWARDS ATOM D
C
          DO 240 L=1,3
            UU(L)=Z(L)/R1
            WW(L)=U(L)/R2
            VV(L)=-UU(L)-WW(L)
  240     CONTINUE
          IA(4)=0
          IF(WRI)WRITE(IOUT,380)I,TYP,IA,C,RAD*(PI-ARCOS(CO)-ARCOS(CO2))
          GO TO 270
C
C     ..... LINEAR PERPENDICULAR BENDING
C
  250     CALL VEKTOR (U,R1,KA,KB,XA,YA,ZA)
          CALL VEKTOR (V,R2,KD,KB,XA,YA,ZA)
          CALL VEKTOR (Z,R2,KC,KB,XA,YA,ZA)
          CALL NORMAL (V,U,W)
          CALL NORMAL (Z,V,X)
          DO 260 L=1,3
            UU(L)=W(L)/R1
            WW(L)=X(L)/R2
            VV(L)=-UU(L)-WW(L)
  260     CONTINUE
          IA(4)=0
          CO=SCALAR(U,W)
          CO2=SCALAR(Z,W)
          QQ(I)=QQ(I)+C*(PI-ARCOS(CO)-ARCOS(CO2))
          IF(WRI)WRITE(IOUT,380)I,TYP,IA,C,RAD*(PI-ARCOS(CO)-ARCOS(CO2))
          GO TO 270
C
C     ..... INVERSE BOND LENGTH
C
  265     CALL VEKTOR (UU,R1,KA,KB,XA,YA,ZA)
          RM1=ONE/R1
          RM2=RM1*RM1
          UU(1)=-RM2*UU(1)*ANULL
          UU(2)=-RM2*UU(2)*ANULL
          UU(3)=-RM2*UU(3)*ANULL
          VV(1)=-UU(1)
          VV(2)=-UU(2)
          VV(3)=-UU(3)
          IA(3)=0
          IA(4)=0
          QQ(I) = QQ(I) + RM1*C
          IF (WRI) WRITE(IOUT,380) I,TYP,IA,C,RM1
  270     if(.not.qonly) then
            DO 290 J=1,4
              M=IA(J)
              IF (M.LE.0) GO TO 290
              M=M-1
              J1=3*(J-1)
              DO 280 L=1,3
                M1=3*M+L
                L1=J1+L
                BMAT(M1,I)=UV(L1)*C+BMAT(M1,I)
  280         CONTINUE
  290       CONTINUE
        endif
   20   continue
        c1=dsqrt(one/c1)/anull
        qq(i)=qq(i)*c1
        if(.not.qonly) then
          do 40 k=1,nad3
            bmat(k,i)=bmat(k,i)*c1
   40     continue
        endif
        if(wri) write(iout,*)' '
   10 continue
      call fioflu(iout)
      return
c
  330 FORMAT(//,1X,'DEFINITION OF INTERNAL COORDINATES          VALUE (A
     &NGS. OR DEGR.)' ,//)
  380 FORMAT(1X,I3,'.',A8,4I3,F16.12,4X,F18.12)
  381 FORMAT(1X,I3,'.',A8,4I3,F16.12,4X,F18.4)
c
c
      END
