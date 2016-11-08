      SUBROUTINE ECORR(T2,T1,II,XE,UOFF,VADD,IPQ,ORBSYM,FLOV,
     .                 NIT,NO,NV,NT,NT2,NTIN,NDIMT2,NIRRED,FZO,FZV,FOCK,
     .                 DOFF,DADD,NM,ITR,map,spin,ssov)
      IMPLICIT INTEGER (A-Z)
      REAL*8 XE,VAL1,VAL2
      REAL*8 T2(*),T1(NO,NV),FOCK(*),II(*)
      INTEGER ORBSYM(NT),FLOV(8,2,4),IPQ(NT2),UOFF(NO,NO),
     .        VADD(NV,NV),FZO(NO),FZV(NV),DOFF(NM),DADD(NM),ITR(NV),
     .        map(*),spin(*),ssov(2,4)
C
C     CALCULATE THE CORRELATION ENERGY
C
      XE=0.0D0
      DO 81 I=1,NO
         IF(FZO(I).EQ.1)GO TO 81
         ISYM=ORBSYM(I)
         ispn=spin(i)
         xi=map(i)
         DO 80 J=1,i-1
            IF(FZO(J).EQ.1)GO TO 80
            jspn=spin(j)
            xj=map(j)
            JSYM=ORBSYM(J)
            IJSYM=IEOR(ISYM,JSYM)
            fa=ssov(ispn,3)-no
            la=ssov(ispn,4)-no
            DO 78 A=fa,la
               if(fzv(a).eq.1)go to 78
               xa=map(a+no)
               ASYM=ORBSYM(A+NO)
               BSYM=IEOR(IJSYM,ASYM)+1
               FB=FLOV(BSYM,jspn,3)-NO
               LB=FLOV(BSYM,jspn,4)-NO
               if(lb.ge.a)lb=a-1
               xia=ipq(max0(xa,xi))+min0(xa,xi)
               DO 76 B=FB,LB
                  xb=map(b+no)
                  xJB=ipq(max0(xb,xj))+min0(xb,xj)
                  IAJB=DOFF(MAX0(xIA,xJB))+DADD(MIN0(xIA,xJB))
                  IJAB=UOFF(I,J)+VADD(A,B)
                  VAL1=T2(IJAB)+T1(I,A)*T1(J,B)
                  val2=ii(iajb)
                  if(ispn.eq.jspn)then
                  xiB=ipq(max0(xb,xi))+min0(xb,xi)
                  xJa=ipq(max0(xa,xj))+min0(xa,xj)
                  IbJa=DOFF(MAX0(xIb,xJa))+DADD(MIN0(xIb,xJa))
                  val1=val1-t1(i,b)*t1(j,a)
                  val2=val2-ii(ibja)
                  endif  
                  XE=XE+val2*VAL1
   76          CONTINUE
   78       CONTINUE
   80    CONTINUE
   81 CONTINUE
C
C     ADD FOCK MATRIX DEPENDENT TERMS
C
      DO 95 I=1,NO
         IF(FZO(I).EQ.1)GO TO 95
         ispn=spin(i)
         isym=orbsym(i)
         fa=flov(isym+1,ispn,3) -no
         la=flov(isym+1,ispn,4) -no
         DO 93 A=fa,la
            FIA=IPQ(A+NO)+I
            XE=XE+FOCK(FIA)*T1(I,A)
   93    CONTINUE
   95 CONTINUE
C
C     WRITE(6,9000)
      WRITE(6,9010) NIT,XE,XE
      WRITE(4,9010) NIT,XE,XE
 9000 FORMAT(/)
 9010 FORMAT('NIT =',I3,5X,'CORREL.ENERGY =',F16.12,3X,D23.16)
c 
      return
      END
