      SUBROUTINE PM(S,C,EV,
     &              E, TEMP, mxnao,msfru)
CThis subroutine localizes M.O. by Pipek-Mezey method,
C                                 JCP 90, 4916 (1989).
C
C  NC: No. of occupied MO's.  Max of 40.
C  NB: No. of basis set.    Max of mxnao.
C  NA: No. of atoms.        Max of 20.
C  MB(20):  No. of b.f. for all atoms.
C  S:  Overlap integral matrix.
C  C:  Occupied MO coefficient matrix.
C  EV: Orbital energies
C
      IMPLICIT REAL*8 (A-H,O-Z)
      integer maxnc, maxna
c (These must be consistent with the sizes of the arrays that were passed in)
      parameter (maxnc=40,maxna=20)
      DIMENSION C(mxnao,mxnao),S(mxnao,mxnao),EV(mxnao),E(2,mxnao),
     &  U(maxnc,maxnc),
     &  V(maxnc,maxnc),VV(maxnc,maxnc),F(maxnc,maxnc),
     &  TEMP(maxnc,mxnao),IN(maxnc)
      common/funcs/na,n3n,nb,nbatri,nc,mb(maxna)
      integer prcntr
   10 FORMAT(20I5)
   40 FORMAT(/14X,'Pipek-Mezey LOCALIZATION METHOD'/
     &          20X,'*** INPUT DATA ***')
   44 FORMAT(I3,' occupied orbitals'/I3,' basis functions'/I3,' atoms')
   45 FORMAT(I5,' basis functions belong to atom', I3)
   50 FORMAT(5F10.6,2x,5F10.6)
   55 FORMAT(7F10.6)
   60 FORMAT(/1x,'CMO coefficient matrix')
   70 FORMAT(/22X,'*** OUTPUT ***'/)
   80 FORMAT(1X,'iteration',I3,':    optimization functional=', F15.10)
   90 FORMAT(/1x,'LMO coefficient matrix')
   95 FORMAT(/1x,'Lagrangian multiplier matrix')
   99 FORMAT(/1x,'Unitary transformation matrix')

      call fioopn(7,'cmo.dat','USER','WRITE')
      call fioopn(8,'lmo.dat','USER','WRITE')

      PI=3.1415926535D0
      ZERO=0.0D0
      ONE=1.0D0
      TWO=2.0D0
      HALF=0.5D0
C
      WRITE (6,40)
      WRITE (6,44) NC,NB,NA
      DO 110 I=1,NA
110   WRITE (6,45) MB(I),I
C
      if (prcntr('IS_ON BRIEF').eq.0) then
        write(6,*)
        WRITE(6,*) ' Overlap integral matrix in terms of AO'
        call matout(s,mxnao,mxnao,nb,nb,6)
        endif
C
      if (prcntr('IS_ON BRIEF').eq.0) then
        WRITE (6,60)
        call matout(c,mxnao,mxnao,nc,nb,6)
        endif
      DO 130 I=1,msfru
      write (7,'(I5,F15.8)') I,EV(I)
      WRITE (7,55) (C(I,J), J=1,NB)
130   CONTINUE
C
      DO 140 I=1,NC
      DO 140 J=1,NC
140   V(I,J)=ZERO
      DO 150 I=1,NC
150   V(I,I)=ONE
C
      WRITE (6,70)
      IT=0
C
C   calculate the optimum functional P of the population localization.
220   IT=IT+1
      IF (IT.GT.200) GOTO 405
      P=ZERO
      DO 190 I=1,NC
      DO 190 J=1,NA
      PP=ZERO
      IF (J.EQ.1) GOTO 160
      M1=M2+1
      GOTO 170
160   M1=1
170   M2=MB(J)+M1-1
      DO 180 K=M1,M2
      DO 180 L=1,NB
      PP=PP+C(I,K)*C(I,L)*S(K,L)
180   CONTINUE
      P=P+PP*PP
190   CONTINUE
      WRITE (6,80) IT,P
      IF (DABS(P-P0).LT.1.0D-12) GOTO 400
      P0=P
C
C  begin 2 by 2 rotations to calculate angle alpha in eq. (13) in
C  JCP 90, 4916 (1989).
      DO 390 I=1,NC-1
      DO 390 J=I+1,NC
      AST=ZERO
      BST=ZERO
      DO 260 K=1,NA
      PSS=ZERO
      PTT=ZERO
      PST=ZERO
      IF (K.EQ.1) GOTO 230
      M1=M2+1
      GOTO 240
230   M1=1
240   M2=MB(K)+M1-1
      DO 250 L=M1,M2
      DO 250 M=1,NB
      PSS=PSS+C(I,L)*C(I,M)*S(L,M)
      PTT=PTT+C(J,L)*C(J,M)*S(L,M)
      PST=PST+(C(I,L)*C(J,M)+C(J,L)*C(I,M))*S(L,M)*HALF
250   CONTINUE
C  PST stands for the value of eq. (31) in P-M's paper.
C
C  AST stands for the value of eq. (29a).
C  BST stands for the value of eq. (29b).
      PSSTT=PSS-PTT
      AST=AST+PST*PST-PSSTT*PSSTT*0.25D0
      BST=BST+PST*PSSTT
260   CONTINUE
C
C  C4 stands for COS(4*alpha) in eq. (13b).
C  S4 stands for SIN(4*alpha) in eq. (13a).
C  C1 stands for COS(alpha)
C  S1 stands for SIN(alpha)
      AB=DSQRT(AST*AST+BST*BST)
      S4=BST/AB
      C4=-AST/AB
      Q1=HALF*(ONE-C4)
      Q2=DSQRT(ONE-Q1)
      Q3=HALF*(ONE+Q2)
      C1=DSQRT(Q3)
      S1=DSQRT(ONE-Q3)
      SS4=4.0D0*C1*S1*(C1*C1-S1*S1)
      ERR=DABS(S4-SS4)
      IF (ERR.GT.1.0D-5) THEN
      AA=C1
      C1=S1
      S1=AA
      END IF
      DO 330 II=1,NB
      E(1,II)=C(I,II)*C1+C(J,II)*S1
      E(2,II)=-C(I,II)*S1+C(J,II)*C1
330   CONTINUE
      DO 340 II=1,NB
      C(I,II)=E(1,II)
      C(J,II)=E(2,II)
340   CONTINUE
      DO 350 II=1,NC
      DO 350 JJ=1,NC
350   U(II,JJ)=ZERO
      DO 360 II=1,NC
360   U(II,II)=ONE
      U(I,I)=C1
      U(J,J)=C1
      U(I,J)=-S1
      U(J,I)=S1
      DO 370 II=1,NC
      DO 370 JJ=1,NC
      VV(II,JJ)=ZERO
      DO 370 KK=1,NC
      VV(II,JJ)=VV(II,JJ)+V(II,KK)*U(KK,JJ)
370   CONTINUE
      DO 380 II=1,NC
      DO 380 JJ=1,NC
      V(II,JJ)=VV(II,JJ)
380   CONTINUE
390   CONTINUE
C  end of 2 by 2 rotation
C
      GOTO 220
C
C  F is Lagrange multiplier matrix.
405   WRITE (6,*) '****************************************************'
      WRITE (6,*) '*** Warning: No convergence after 200 iterations ***'
      WRITE (6,*) '****************************************************'
400   DO 410 I=1,NC
      DO 410 J=1,NC
      F(I,J)=0.0D0
      DO 410 K=1,NC
      F(I,J)=F(I,J)+V(K,I)*EV(K)*V(K,J)
410   CONTINUE
C
C  Reorder the rows according to the values of F(i,i)
      DO 440 I=1,NC
      FF=1.0D50
      DO 430 J=1,NC
      DO 420 K=1,I
420   IF(J.EQ.IN(K)) GOTO 430
      IF(F(J,J).LT.FF) THEN
      FF=F(J,J)
      IN(I)=J
      END IF
430   CONTINUE
440   CONTINUE
      DO 450 I=1,NC
      DO 450 J=1,NC
      II=IN(I)
      JJ=IN(J)
      TEMP(I,J)=F(II,JJ)
450   CONTINUE
      DO 460 I=1,NC
      DO 460 J=1,NC
      F(I,J)=TEMP(I,J)
460   CONTINUE
      WRITE (6,95)
      DO 470 I=1,NC
      WRITE (6,50) (F(I,J), J=1,NC)
470   CONTINUE
C
C   Reorder the rows of C matrix.
      DO 480 I=1,NC
      DO 480 J=1,NB
      II=IN(I)
      TEMP(I,J)=C(II,J)
480   CONTINUE
      DO 490 I=1,NC
      DO 490 J=1,NB
      C(I,J)=TEMP(I,J)
490   CONTINUE
      if (prcntr('IS_ON BRIEF').eq.0) then
        WRITE (6,90)
        call matout(c,mxnao,mxnao,nb,nb,6)
        endif
      DO 500 I=1,NC
      write (8,'(I5)') I
      WRITE (8,55) (C(I,J), J=1,NB)
500   CONTINUE
C
C   Reorder the rows of V matrix.
      DO 510 I=1,NC
      DO 510 J=1,NC
      II=IN(I)
      JJ=IN(J)
      TEMP(I,J)=V(II,JJ)
510   CONTINUE
      DO 520 I=1,NC
      DO 520 J=1,NC
      V(I,J)=TEMP(I,J)
520   CONTINUE
      if (prcntr('IS_ON BRIEF').eq.0) then
        WRITE (6,99)
        call matout(c,maxnc,maxnc,nc,nc,6)
        endif
C
      RETURN
      END
