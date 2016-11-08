      SUBROUTINE DAVLIX(A,EIG,VEC,G,FV1,FV2,CC,ALPA,B,AB,Q,Z,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NDTRI),EIG(MXVEC),VEC(MXVEC*MXVEC)
      DIMENSION G(MXTRI),FV1(MXVEC),FV2(MXVEC)
      DIMENSION CC(MXTRI),ALPA(MXVEC)
      DIMENSION B(NIND),AB(NIND),Q(NIND),Z(NIND),D(NIND)
      COMMON/CONVS/IQCONV,IDCONV,IECONV
      COMMON/MFSEC/MFILE,NSECT
      COMMON/POSIT/IABP(500)
      COMMON/ROOTS/NROOT,NPAIR
      COMMON/SIGNS/IOFF(32897),IPRNT
      common/stori/istor(50),maxstr
      common/storr/cstor(50),eigsv(50)
      COMMON/TAPES/ITAP35,ITAP91,ITAP92
      COMMON/CI102/NIND,NDEP,NDTRI
      COMMON/CI103/MXITER,MXVEC,MXTRI
      COMMON/CI104/NVEC
      DATA A00,ONE,TEN / 0.0D+00 , 1.0D+00 , 10.0D+00 /
      DATA XLIMIT / 1.0D-06 /
    1 FORMAT(//,2X,' THE SIMULTANEOUS EXPANSION METHOD FOR THE LARGE REA
     *L SYMMETRIC MATRIX'/)
    2 FORMAT(3X,127('*')/)
    3 FORMAT(2X,' ITERATION',12X,' MXQNORM',19X,' MXDNORM',
     1      19X,' MXEDIF',20X,' EIGEN(1)'/)
    4 FORMAT(3X,I5,4(12X,F15.10))
    5 FORMAT(2X,' ')
    6 FORMAT(2X,' NUMBER OF EXPANDED VECTORS = ',I5/)
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C    THE DAVIDSON-LIU ALGORITHM                                        C
C                                                                      C
C    REFERENCES : E.R. DAVIDSON, J. COMP. PHYS. 17,87(1975).           C
C               : B. LIU, PROCEEDING FROM A WORKSHOP OF THE NRCC,      C
C                         BERKELEY, (1978) P.49.                       C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
      WRITE(6,1)
      WRITE(3,1)
C
C   READ BACK THE MO HESSIAN MATRIX
      IF(IPRNT.GT.5) THEN
        WRITE(6,*) '  '
        WRITE(6,*) ' THE A MATRIX'
        CALL PRINT(A,NDTRI,NIND,6)
      END IF
C
C   SET UP AN ARRAY FOR WORKING TAPES
      NIND2=NIND*2
      NINDL=(NIND2-1)/NSECT+1
      DO 101 I=1,MXVEC
      IABP(I)=(I-1)*NINDL+1
      IF(IPRNT.GT.5) THEN
        WRITE(6,*) I,IABP(I)
      END IF
  101 CONTINUE
C
C=======================================================================
C   STEP 1 (LIU)
C   STEP A (DAVIDSON)
C=======================================================================
C
C   FIND THE SMALLEST DIAGONAL ELEMENTS
      call izero(istor,maxstr)
      call zero(cstor,maxstr)
      call zero(eigsv,maxstr)
      DO 104 IPOS=1,NROOT
      CMIN=1.0D+06
      DO 103 I=1,NIND
      II=IOFF(I+1)
      CTEMP=A(II)
      IF(CTEMP.GT.CMIN) GO TO 103
      IF(IPOS.GT.1) THEN
        DO 102 J=1,IPOS-1
        ITEMP=ISTOR(J)
        IF(ITEMP.EQ.I) GO TO 103
  102   CONTINUE
      END IF
      CMIN=CTEMP
      III=I
  103 CONTINUE
      ISTOR(IPOS)=III
      CSTOR(IPOS)=CMIN
  104 CONTINUE
      IF(IPRNT.GT.5) THEN
        WRITE(6,*) ' ISTOR = ',ISTOR
        WRITE(6,*) ' CSTOR = ',CSTOR
      END IF
C
C=======================================================================
C   STEP 2 (LIU)
C   STEP A (DAVIDSON)
C=======================================================================
C
C   SET UP INITIAL EIGENVECTORS
      CALL RFILE(ITAP91)
      CALL RFILE(ITAP92)
C
      IF(IPRNT.GT.5) THEN
        WRITE(6,*) ' INITIAL B  AND AB VECTORS'
      END IF
      DO 106 II=1,NROOT
      CALL ZERO(B,NIND)
      CALL ZERO(AB,NIND)
      IPOS=ISTOR(II)
      B(IPOS)=ONE
      DO 105 I=1,NIND
      IJP=IOFF(MAX0(I,IPOS))+MIN0(I,IPOS)
      AB(I)=A(IJP)
  105 CONTINUE
      CALL RWRIT(ITAP91,B,NIND2,IABP(II))
      CALL RWRIT(ITAP92,AB,NIND2,IABP(II))
CCC     WRITE(6,*) ' 105  B VECTOR, II = ',II
CCC     WRITE(6,*) (B(IV),IV=1,NIND)
CCC     WRITE(6,*) ' 105  AB VECTOR, II = ',II
CCC     WRITE(6,*) (AB(IV),IV=1,NIND)
  106 CONTINUE
C
C   FORM THE INITIAL G MATRIX
      CALL ZERO(G,MXTRI)
      IIJJ=0
      DO 108 II=1,NROOT
      I=ISTOR(II)
      DO 107 JJ=1,II
      IIJJ=IIJJ+1
      J=ISTOR(JJ)
      IJ=IOFF(MAX0(I,J))+MIN0(I,J)
      G(IIJJ)=A(IJ)
  107 CONTINUE
  108 CONTINUE
C
C   SET UP INITIAL VARIABLES
      WRITE(6,2)
      WRITE(6,3)
      WRITE(6,2)
      QLIMIT=TEN**(-IQCONV)
      DLIMIT=TEN**(-IDCONV)
      ELIMIT=TEN**(-IECONV)
      NVEC=NROOT
      NITER=0
      IEND=0
      NADDV=NVEC
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C   ITERATIVE PROCEDURE STARTS HERE
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  300 CONTINUE
      NITER=NITER+1
      IF(NADDV.EQ.0) THEN
        ISTART=IEND
        IEND=ISTART
      ELSE
        ISTART=IEND+1
        IEND=ISTART+NADDV-1
      END IF
      IF(IPRNT.GT.5) THEN
        WRITE(6,*) ' NITER,NVEC,ISTART,IEND = ',NITER,NVEC,ISTART,IEND
        WRITE(3,*) ' NITER,NVEC,ISTART,IEND = ',NITER,NVEC,ISTART,IEND
      END IF
C
      IF(NITER.GT.MXITER) GO TO 299
C
C   FORM THE G MATRIX
      NVTRI=IOFF(NVEC+1)
      IF(NITER.EQ.1) GO TO 301
      DO 109 II=ISTART,IEND
      CALL RREAD(ITAP91,B,NIND2,IABP(II))
      DO 109 JJ=1,IEND
      IIJJ=IOFF(MAX0(II,JJ))+MIN0(II,JJ)
      CALL RREAD(ITAP92,AB,NIND2,IABP(JJ))
      CALL DOTAB(VALU,B,AB,NIND)
      G(IIJJ)=VALU
  109 CONTINUE
  301 CONTINUE
      DO 110 I=1,NVTRI
      CC(I)=G(I)
  110 CONTINUE
      IF(IPRNT.GT.5) THEN
        WRITE(6,*) '   '
        WRITE(6,*) ' G MATRIX'
        CALL PRINT(G,NVTRI,NVEC,6)
      END IF
C
C   DIAGONALIZE THE G MATRIX
      CALL RSP(NVEC,NVEC,NVTRI,CC,EIG,1,VEC,FV1,FV2)
      IF(IPRNT.GT.5) THEN
        WRITE(6,*) '   '
        WRITE(6,*) ' EIGENVECTOR MATRIX OF THE G MATRIX'
        CALL EIVOUT(VEC,EIG,NVEC,NVEC,NVEC,NVEC,6)
        WRITE(6,*) ' EIGENVALUES'
        WRITE(6,*) (EIG(I),I=1,NVEC)
      END IF
C
C=======================================================================
C   STEP 3 (LIU)
C   STEP B (DAVIDSON)
C=======================================================================
C
C   FIND ADDITIONAL VECTORS
      NVOLD=NVEC
      NADDV=0
      QNMAX=A00
      DNMAX=A00
C
C///////////////////////////////////////////////////////////////////////
C   A BIG LOOP STARTS HERE
C///////////////////////////////////////////////////////////////////////
      DO 130 III=1,NROOT
      EPS=EIG(III)
      IADD=(III-1)*NVOLD
      DO 111 I=1,NVOLD
      IPOS=IADD+I
      ALPA(I)=VEC(IPOS)
  111 CONTINUE
      IF(IPRNT.GT.5) THEN
        WRITE(6,*) ' III = ',III,' EPS = ',EPS
        WRITE(6,*) ' ALPA '
        WRITE(6,*) (ALPA(IV),IV=1,NVOLD)
      END IF
C
C   FORM THE Q VECTOR
      CALL ZERO(Q,NIND)
      DO 115 II=1,NVOLD
      CALL RREAD(ITAP91,B,NIND2,IABP(II))
CCC     WRITE(6,*) ' DO 115  B VECTOR, II = ',II
CCC     WRITE(6,*) (B(IV),IV=1,NIND)
      CALL RREAD(ITAP92,AB,NIND2,IABP(II))
CCC     WRITE(6,*) ' DO 115  AB VECTOR, II = ',II
CCC     WRITE(6,*) (AB(IV),IV=1,NIND)
      FACT1=ALPA(II)
      CALL SCALE(Q,AB,FACT1,NIND)
      FACT2=-ALPA(II)*EPS
      CALL SCALE(Q,B,FACT2,NIND)
  115 CONTINUE
      IF(IPRNT.GT.5) THEN
        WRITE(6,*) ' Q VECTOR'
        WRITE(6,*) (Q(IV),IV=1,NIND)
      END IF
C
C=======================================================================
C   STEP C (DAVIDSON)
C=======================================================================
C
C   CHECK CONVERGENCE OF THE Q VECTOR
      CALL DOTAB(QNORM,Q,Q,NIND)
      QNORM=DSQRT(QNORM)
      IF(IPRNT.GT.5) THEN
        WRITE(6,*) ' QNORM = ',QNORM
      END IF
      QNMAX=DMAX1(QNMAX,QNORM)
      IF(QNORM.LT.QLIMIT) GO TO 130
C
C=======================================================================
C   STEP D (DAVIDSON)
C=======================================================================
C
C   FORM THE Z VECTOR
      DO 116 I=1,NIND
      II=IOFF(I+1)
      DENOM=EPS-A(II)
      IF(DABS(DENOM).LT.XLIMIT) GO TO 202
      Z(I)=Q(I)/DENOM
      GO TO 116
  202 CONTINUE
      Z(I)=A00
  116 CONTINUE
      IF(IPRNT.GT.5) THEN
        WRITE(6,*) ' Z VECTOR'
        WRITE(6,*) (Z(IV),IV=1,NIND)
      END IF
C
C=======================================================================
C   STEP E (DAVIDSON)
C=======================================================================
C
C   FORM THE D VECTOR
      DO 118 II=1,NVEC
      CALL RREAD(ITAP91,B,NIND2,IABP(II))
      CALL DOTAB(EBZ,B,Z,NIND)
      FACT=-EBZ
      CALL SCALE(Z,B,FACT,NIND)
  118 CONTINUE
      DO 120 I=1,NIND
      D(I)=Z(I)
  120 CONTINUE
      IF(IPRNT.GT.5) THEN
        WRITE(6,*) ' D VECTOR'
        WRITE(6,*) (D(IV),IV=1,NIND)
      END IF
C
C=======================================================================
C   STEP F (DAVIDSON)
C=======================================================================
C
C   NORMALIZE THE NEW VECTOR
C   AND CHECK CONVERGENCE OF THE D VECTOR
      CALL DOTAB(DNORM,D,D,NIND)
      DNORM=DSQRT(DNORM)
      IF(IPRNT.GT.5) THEN
        WRITE(6,*) ' DNORM = ',DNORM
      END IF
      DNMAX=DMAX1(DNMAX,DNORM)
      IF(DNORM.LT.DLIMIT) GO TO 130
      NADDV=NADDV+1
C
C   FORM THE B VECTOR
      DO 121 I=1,NIND
      B(I)=D(I)/DNORM
  121 CONTINUE
      IF(IPRNT.GT.5) THEN
        WRITE(6,*) ' B VECTOR'
        WRITE(6,*) (B(IV),IV=1,NIND)
      END IF
C
C=======================================================================
C   STEP G (DAVIDSON)
C=======================================================================
C
C   FORM AB VECTOR
      CALL CIMUL(AB,A,B,NIND)
C
C   ADD A NEW VECTOR
      INEW=NVEC+1
      CALL RWRIT(ITAP91,B,NIND2,IABP(INEW))
      CALL RWRIT(ITAP92,AB,NIND2,IABP(INEW))
      NVEC=NVEC+1
      IF(NVEC.GT.MXVEC) GO TO 305
C
  130 CONTINUE
C///////////////////////////////////////////////////////////////////////
C   A BIG LOOP ENDS HERE
C///////////////////////////////////////////////////////////////////////
  305 CONTINUE
      IF(IPRNT.GT.5) THEN
        WRITE(6,*) ' NADDV = ',NADDV
        WRITE(3,*) ' NADDV = ',NADDV
        WRITE(6,*) '  '
      END IF
C
      DEMAX=A00
      DO 131 I=1,NROOT
      DELE=EIG(I)-EIGSV(I)
      EIGSV(I)=EIG(I)
      ABDELE=DABS(DELE)
      IF(ABDELE.GT.DEMAX) DEMAX=ABDELE
      IF(IPRNT.GT.5) THEN
        WRITE(6,*) ' I , DELE , DEMAX ',I,DELE,DEMAX
      END IF
  131 CONTINUE
C
C   PRINT THE PRESENT RESULT
      WRITE(6,4) NITER,QNMAX,DNMAX,DEMAX,EIG(1)
      IF(NVEC.GT.MXVEC) GO TO 399
      IF(NADDV.EQ.0) GO TO 310
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C   ITERATIVE PROCEDURE ENDS HERE
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      GO TO 300
C
  310 CONTINUE
      IF(IPRNT.GT.5) THEN
        WRITE(6,*) '  '
        WRITE(6,*) ' NO MORE EXPANSION OF THE CI VECTORS'
      END IF
      GO TO 400
C
C///////////////////////////////////////////////////////////////////////
C   MAXIMUM ITERATION HAS COMPLETED
C///////////////////////////////////////////////////////////////////////
  299 CONTINUE
      WRITE(6,*) ' MAXIMUM ITERATION HAS BEEN ATTAINED'
      GO TO 400
C///////////////////////////////////////////////////////////////////////
C   MAXIMUM NUMBER OF VECTORS HAVE BEEN EXPANDED
C///////////////////////////////////////////////////////////////////////
  399 CONTINUE
      NVEC=MXVEC
      WRITE(6,*) ' MAXIMUM VECTORS HAVE BEEN OBTAINED --- STOP !'
C
  400 CONTINUE
      WRITE(6,5)
      WRITE(6,2)
      WRITE(6,6) NVEC
C
      CALL RCLOSE(ITAP91,3)
      CALL RCLOSE(ITAP92,3)
      RETURN
      END
