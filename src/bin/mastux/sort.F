      SUBROUTINE SORT(EIG,OCC,ESO,UAS,EAO,EE,EGN,OCN,
     1                KORD,LORD,MORD,NORD)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION EIG(NBASIS),OCC(NBASIS),ESO(NBASIS,NBASIS)
      DIMENSION UAS(NBFAO,NBFSO),EAO(NBFAO,NBFSO)
      DIMENSION EE(NBFAO,NBFAO),EGN(NBASIS),OCN(NBASIS)
      DIMENSION KORD(NBASIS),LORD(NBASIS),MORD(NBASIS),NORD(NBASIS)
      DIMENSION A30(200)
      CHARACTER*8 SCFTYP,CITYP,DERTYP
      dimension itocc(20)
      COMMON/BASIS/NBASIS,NTRI,NST,NSYM,NBFAO,NBFSO,NBATRI
      COMMON/CALIF/LPARA(1024),APARA(1024)
      COMMON/COUPL/ALPC(15),BETC(15),ALPT(15),BETT(15)
      COMMON/ENRGY/ESCF,ENUC
      COMMON/FUNCS/NTYPES,NATOM,N3N
      COMMON/GRSCF/FOCC(10),NSORB(10)
      COMMON/OPEN1/NUNIQ,LOPN(5),MOPN(5,5)
      COMMON/OPEN2/GOCC(10)
      COMMON/SIGNS/IOFF(256),IPRNT
      COMMON/SYMMS/IDEG(10),NDEG(10),LSYM(10)
      COMMON/TEMPS/I30(400)
      COMMON/TYPES/SCFTYP,CITYP,DERTYP
      COMMON/CI101/IOPEN,IOCC,JOCC,KOCC
      EQUIVALENCE (I30,A30)
      DATA A00,ONE,TWO / 0.0D+00 , 1.0D+00 , 2.0D+00 /
    1 FORMAT(//,2X,' IOCC = ',I5,3X,' KOCC = ',I5,3X,' JOCC = ',I5/)
    2 FORMAT(//,2X,' ILLEGAL SCFTYP ',A8/)
    3 FORMAT(//,2X,' ALPC MATRIX'/)
    4 FORMAT(//,2X,' BETC MATRIX'/)
    5 FORMAT(//,2X,' SORTED ESO MATRIX'/)
    6 FORMAT(//,2X,' SORTED EAO MATRIX'/)
C
C   CALCULATE NUMBERS OF OCCUPATIONS
      IOCC=0
      KOCC=0
      DO 101 I=1,NBASIS
      NOCC=INT(OCC(I))
      IF(NOCC-1) 201,202,203
  201 IF(OCC(I).GT.A00) GO TO 202
      GO TO 101
  202 KOCC=KOCC+1
      itocc(kocc)=i
      GO TO 101
  203 IOCC=IOCC+1
  101 CONTINUE
      JOCC=IOCC+KOCC
      WRITE(6,1) IOCC,KOCC,JOCC
c
c   Test ``orbital'' energy ordering for a TCSCF case
      if(SCFTYP.EQ.'TCSCF   ') then
        write(3,*) '  TCSCF wavefunction'
        write(3,*) '  itc1   = ',itc1
        write(3,*) '  itc2   = ',itc2
      end if
C
C   SORT EIGENVALUES AND EIGENVECTORS
cc    write(6,*) ' before ordr'
cc    write(6,*) ' eig  = '
cc    write(6,10) (eig(i),i=1,nbasis)
cc    write(6,*) ' occ  = '
cc    write(6,10) (occ(i),i=1,nbasis)
cc 10 format(2x,10f12.5)
c
      CALL ORDR(EIG,NORD,NBASIS)
c
cc    write(6,*) ' after ordr'
cc    write(6,*) ' nord = '
cc    write(6,11) (nord(i),i=1,nbasis)
cc    write(6,*) ' eig  = '
cc    write(6,10) (eig(i),i=1,nbasis)
cc    write(6,*) ' occ  = '
cc    write(6,10) (occ(i),i=1,nbasis)
cc 11 format(2x,20i5)
cc
      DO 102 II=1,NBASIS
      NN=NORD(II)
      MORD(II)=NN
      OCN(II)=OCC(NN)
  102 CONTINUE
cc    write(6,*) ' before reorder'
cc    write(6,*) ' eig  = '
cc    write(6,10) (eig(i),i=1,nbasis)
cc    write(6,*) ' ocn  = '
cc    write(6,10) (ocn(i),i=1,nbasis)
c
      CALL REORDR(OCN,NORD,NBASIS)
c
cc    write(6,*) ' after reorder'
cc    write(6,*) ' nord = '
cc    write(6,11) (nord(i),i=1,nbasis)
cc    write(6,*) ' eig  = '
cc    write(6,10) (eig(i),i=1,nbasis)
cc    write(6,*) ' ocn  = '
cc    write(6,10) (ocn(i),i=1,nbasis)
C
C I'TH PITZER'S MO # IS GOING TO GFOCK MO # "LORD(I)"
C
      DO 103 II=1,NBASIS
      NN=NORD(II)
      KORD(II)=MORD(NN)
      EGN(II)=EIG(NN)
  103 CONTINUE
      DO 104 I=1,NBASIS
      EIG(I)=EGN(I)
      OCC(I)=OCN(I)
      LORD(I)=KORD(I)
  104 CONTINUE
c
cc    write(6,*) ' before oporder'
cc    write(6,*) ' lord = '
cc    write(6,11) (lord(i),i=1,nbasis)
cc    write(6,*) ' eig  = '
cc    write(6,10) (eig(i),i=1,nbasis)
cc    write(6,*) ' occ  = '
cc    write(6,10) (occ(i),i=1,nbasis)
c
c   Sort eigenvalues among open-shells if necessary
c   This part assumes that the occupation numbers are the same
c   for all open-shells
      IF(SCFTYP.EQ.'TCSCF   ') GO TO 204
c
cc    write(6,*) ' before opordr'
cc    write(6,*) ' eig  = '
cc    write(6,10) (eig(i),i=1,nbasis)
c
      if(kocc.ge.2) then
        do 121 i=1,nbasis
        kord(i)=lord(i)
  121 continue
        call opordr(eig(iocc+1),nord(iocc+1),kocc)
        do 122 i=1,kocc
        nn=nord(iocc+i)
        kord(iocc+i)=lord(iocc+nn)
  122 continue
cc    write(6,*) ' after opordr'
cc    write(6,*) ' nord = '
cc    write(6,11) (nord(iocc+i),i=1,kocc)
cc    write(6,10) (eig(iocc+i),i=1,kocc)
        do 123 i=1,nbasis
        lord(i)=kord(i)
  123 continue
cc    write(6,*) ' after opordr'
cc    write(6,11) (lord(i),i=1,nbasis)
cc    write(6,10) (eig(i),i=1,nbasis)
      end if
C
  204 continue
      IF(SCFTYP.EQ.'CLSCF   ') GO TO 205
      IF(SCFTYP.EQ.'GRSCF   '.AND.IOPEN.EQ.0) GO TO 205
      IF(SCFTYP.EQ.'HSOP    ') GO TO 206
      IF(SCFTYP.EQ.'GRSCF   '.AND.IOPEN.NE.0) GO TO 207
      IF(SCFTYP.EQ.'TCSCF   ') GO TO 208
      WRITE(6,2) SCFTYP
      STOP
C
C   CLOSED SHELL ONLY
  205 CONTINUE
C   TYPE1=CLOSED SHELL
C   TYPE2=VIRTUAL
      NTYPES=1
      FOCC(1)=TWO
      FOCC(2)=A00
      ALPC(1)=A00
      ALPC(2)=A00
      ALPC(3)=A00
      BETC(1)=A00
      BETC(2)=A00
      BETC(3)=A00
      NSORB(1)=IOCC
      NSORB(2)=NBASIS-IOCC
      GO TO 210
C
C   HIGH SPIN OPEN-SHELL
  206 CONTINUE
C   TYPE1=CLOSED SHELL
C   TYPE2=ALPHA SPIN SHELL
C   TYPE3=VIRTUAL
      NTYPES=2
      FOCC(1)=TWO
      FOCC(2)=ONE
      FOCC(3)=A00
      ALPC(1)=A00
      ALPC(2)=A00
      ALPC(3)=A00
      ALPC(4)=A00
      ALPC(5)=A00
      ALPC(6)=A00
      BETC(1)=A00
      BETC(2)=ONE
      BETC(3)=ONE
      BETC(4)=A00
      BETC(5)=A00
      BETC(6)=A00
      NSORB(1)=IOCC
      NSORB(2)=KOCC
      NSORB(3)=NBASIS-JOCC
      GO TO 210
C
C   GENERALIZED SCF OPEN-SHELL
  207 CONTINUE
      NTYPES=NUNIQ+1
      NNP=IOFF(NTYPES+2)
      CALL ZERO(ALPC,NNP)
      CALL ZERO(BETC,NNP)
C   FOR CLOSED SHELL
      FOCC(1)=TWO
      NSORB(1)=IOCC
C   FOR OPEN-SHELL
      II=1
      DO 106 I=1,NUNIQ
      II=II+1
      FOCC(II)=GOCC(I)
      NSORB(II)=LOPN(I)
  106 CONTINUE
      DO 107 I=1,NUNIQ
      II=I+1
      K=MOPN(I,1)
      DO 107 J=1,I
      JJ=J+1
      L=MOPN(J,1)
      IIJJ=IOFF(II)+JJ
      KL=IOFF(MAX0(K,L))+MIN0(K,L)
      ALPC(IIJJ)=ALPT(KL)
      BETC(IIJJ)=BETT(KL)
  107 CONTINUE
C   VIRTUAL
      FOCC(II+1)=A00
      NSORB(II+1)=NBASIS-JOCC
      WRITE(6,3)
      CALL PRINT(ALPC,15,NTYPES+1,6)
      WRITE(6,4)
      CALL PRINT(BETC,15,NTYPES+1,6)
      GO TO 210
C
C   TCSCF
  208 CONTINUE
C   TYPE1=CLOSED SHELL
C   TYPE2=GVA SHELL
C   TYPE3=GVB SHELL
C   TYPE4=VIRTUAL
      CALL ZERO(ALPC,15)
      CALL ZERO(BETC,15)
      NTYPES=3
      FOCC(1)=TWO
      IF(IOPEN.EQ.0) GO TO 209
c
      FOCC(2)=OCC(IOCC+1)
      FOCC(3)=OCC(IOCC+2)
c
      FOCC(4)=A00
      ALPC(1)=A00
      ALPC(2)=A00
      ALPC(3)=ONE-ONE/FOCC(2)
      ALPC(4)=A00
      ALPC(5)=ONE
      ALPC(6)=ONE-ONE/FOCC(3)
      ALPC(7)=A00
      ALPC(8)=A00
      ALPC(9)=A00
      ALPC(10)=A00
      BETC(1)=A00
      BETC(2)=A00
      BETC(3)=ONE
      BETC(4)=A00
      BETC(5)=BETT(2)
      BETC(6)=ONE
      BETC(7)=A00
      BETC(8)=A00
      BETC(9)=A00
      BETC(10)=A00
  209 CONTINUE
      NSORB(1)=IOCC
      NSORB(2)=1
      NSORB(3)=1
      NSORB(4)=NBASIS-IOCC-2
C
  210 CONTINUE
c
      DO 110 II=1,NBASIS
      NN=LORD(II)
      DO 109 I=1,NBASIS
  109 EE(I,II)=ESO(I,NN)
  110 CONTINUE
      DO 111 I=1,NBASIS
      DO 111 J=1,NBASIS
      ESO(I,J)=EE(I,J)
  111 CONTINUE
      IF(IPRNT.LE.3) GO TO 212
      WRITE(6,5)
      CALL EIGOUT(ESO,EIG,OCC,NBFSO,NBFSO,NBASIS,NBASIS,6)
C
C   TRANSFORM SORTED EIGENVECTORS FROM SO TO AO BASIS
  212 CONTINUE
      DO 113 I=1,NBFAO
      DO 113 J=1,NBASIS
      VALU=A00
      DO 112 K=1,NBFSO
      VALU=VALU+UAS(I,K)*ESO(K,J)
  112 CONTINUE
      EAO(I,J)=VALU
  113 CONTINUE
      IF(IPRNT.LE.3) GO TO 213
      WRITE(6,6)
      CALL EIGOUT(EAO,EIG,OCC,NBFAO,NBFSO,NBFAO,NBASIS,6)
C
C   DEFINE REAL CONSTANTS
  213 CONTINUE
      A30(1)=ENUC
      A30(2)=ESCF
      DO 115 I=1,15
      A30(I+10)=ALPC(I)
      A30(I+25)=BETC(I)
  115 CONTINUE
      DO 116 I=1,10
      A30(I+40)=FOCC(I)
  116 CONTINUE
C
C   DEFINE INTEGER CONSTANTS
      LPARA(6)=IOPEN
      LPARA(7)=IOCC
      LPARA(8)=KOCC
      LPARA(9)=JOCC
      LPARA(18)=NTYPES
      IF(SCFTYP.EQ.'CLSCF   ') LPARA(19)=1
      IF(SCFTYP.EQ.'HSOP    ') LPARA(19)=2
      IF(SCFTYP.EQ.'GRSCF   ') LPARA(19)=3
      IF(SCFTYP.EQ.'TCSCF   ') LPARA(19)=4
      IF(SCFTYP.EQ.'MCSCF   ') LPARA(19)=5
      IF(CITYP.EQ.'SCF     ') LPARA(20)=1
      IF(CITYP.EQ.'CI      ') LPARA(20)=2
      IF(CITYP.EQ.'MCSCF   ') LPARA(20)=3
      IF(CITYP.EQ.'MCSCFCI ') LPARA(20)=4
      IF(DERTYP.EQ.'FIRST  ') LPARA(21)=1
      IF(DERTYP.EQ.'SECOND ') LPARA(21)=2
      IF(DERTYP.EQ.'THIRD  ') LPARA(21)=3
      IF(DERTYP.EQ.'FOURTH ') LPARA(21)=4
C
      DO 117 I=1,10
      LPARA(I+30)=NDEG(I)
      LPARA(I+40)=NSORB(I)
  117 CONTINUE
C
      RETURN
      END
