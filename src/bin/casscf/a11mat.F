      SUBROUTINE A11MAT(H,Q,X,NIJ,G,CC,LBLI,BUFI)
      IMPLICIT REAL*8 (A-H,O-Z)                                         SMM16380
      DIMENSION H(NTRI),Q(NTRI),X(NBASIS,NBASIS),NIJ(NTRI,2)
      DIMENSION G(noctot),CC(ntot),LBLI(MAXBF4),BUFI(MAXBF2)
      COMMON/BASIS/NBASIS,NTRI,NTOT,NSYM
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4
      COMMON/SIGNS/IOFF(32897),IPRNT  
      COMMON/TAPES/ITAP35,ITAP56
      common/works/jtape1,jtape2
      COMMON/CI101/NOCMO,NOCTRI,NOCTOT
      COMMON/CI102/NIND,NINTRI
      DATA A00,HALF / 0.0D+00 , 0.5D+00 /                               SMM16500
    1 FORMAT(//,2X,' THE H MATRIX IN A11MAT'/)                          SMM16510
    2 FORMAT(//,2X,' THE Q MATRIX IN A11MAT'/)                          SMM16520
C                                                                       SMM16540
C   FORM THE A11 MATRIX                                                 SMM16550
c......................................................................
c   The G and CC matrices were already readin in preceding subroutines
cc*      call rdin22(g,noctot,itap56,noctot,nocmo,lbli,bufi)
cc*      call rdint2(cc,ntot,itap35,ntot,lbli,bufi)
c......................................................................
      if(iprnt.ge.4) then
        WRITE(6,1)                                         
        CALL PRINT(H,NTRI,NBASIS,6)                    
        WRITE(6,2)                                
        CALL PRINT(Q,NTRI,NBASIS,6)         
      end if
C                                                                       SMM16660
      call rfile(jtape1)
      ibl=0
      DO 102 II=1,NIND                                                  SMM16690
      I=NIJ(II,1)                                                       SMM16700
      J=NIJ(II,2)                                                       SMM16710
      DO 101 JJ=1,II                                                    SMM16720
      K=NIJ(JJ,1)                                                       SMM16730
      L=NIJ(JJ,2)                                                       SMM16740
      HYMAT=Y4MAT(I,J,K,L,H,Q,G,CC)
C                                                                       SMM16780
      HXMAT=A00                                                         SMM16790
      IF(I.NE.L) GO TO 301                                              SMM16800
      HXMAT=HXMAT + (X(J,K) + X(K,J))                                   SMM16810
  301 CONTINUE                                                          SMM16820
      IF(J.NE.L) GO TO 302                                              SMM16830
      HXMAT=HXMAT - (X(I,K) + X(K,I))                                   SMM16840
  302 CONTINUE                                                          SMM16850
      IF(I.NE.K) GO TO 303                                              SMM16860
      HXMAT=HXMAT - (X(J,L) + X(L,J))                                   SMM16870
  303 CONTINUE                                                          SMM16880
      IF(J.NE.K) GO TO 304                                              SMM16890
      HXMAT=HXMAT + (X(I,L) + X(L,I))                                   SMM16900
C                                                                       SMM16910
  304 CONTINUE                                                          SMM16920
      HIJKL=HYMAT + HXMAT*HALF                                          SMM16930
      ibl=ibl+1
      lbli(ibl+ibl-1)=ii
      lbli(ibl+ibl  )=jj
      bufi(maxbuf+ibl)=hijkl
      if(ibl.lt.maxbuf) go to 101
      ibl=0
      call swrit(jtape1,lbli,maxbf4)
  101 CONTINUE                                                          SMM16950
  102 CONTINUE                                                          SMM16960
C                                                                       SMM16970
C   Store the last buffer
      ibl=ibl+1
      lbli(ibl+ibl-1)=0 
      lbli(ibl+ibl  )=0 
      bufi(maxbuf+ibl)=a00
      call swrit(jtape1,lbli,maxbf4)
      write(3,*) ' final ibl in a11mat = ',ibl
      call rclose(jtape1,3)
c
      RETURN                                   
      END                                                   
      FUNCTION Y4MAT(I,J,K,L,H,Q,G,CC)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION H(NTRI),Q(NTRI),G(NOCTOT),CC(NTOT)
      COMMON/BASIS/NBASIS,NTRI,NTOT,NSYM
      COMMON/SIGNS/IOFF(32897),IPRNT
      COMMON/CI101/NOCMO,NOCTRI,NOCTOT
      DATA A00, TWO, FOUR / 0.0D+00 , 2.0D+00 , 4.0D+00 /         
c
c   One-electron contribution
      yone=a00
      jl=i21d(j,l)
      ik=i21d(i,k)
      if(jl.le.noctri) then
        yone=yone+q(jl)*h(ik)
      end if
      if(ik.le.noctri) then
        yone=yone+q(ik)*h(jl)
      end if
      jk=i21d(j,k)
      il=i21d(i,l)
      if(jk.le.noctri) then
        yone=yone-q(jk)*h(il)
      end if
      if(il.le.noctri) then
        yone=yone-q(il)*h(jk)
      end if
c
c   Two-electron contribution
      ytwo2=a00
      ytwo4=a00
      DO 101 M=1,NOCMO
      IM=I21D(I,M)
      JM=I21D(J,M)
      DO 101 N=1,NOCMO
      MN=I21D(M,N)
      KN=I21D(K,N)
      LN=I21D(L,N)
      JLMN=I21D(JL,MN)
      IKMN=I21D(IK,MN)
      JMLN=I21D(JM,LN)
      IMKN=I21D(IM,KN)
      jkmn=i21d(jk,mn)
      ilmn=i21d(il,mn)
      jmkn=i21d(jm,kn)
      imln=i21d(im,ln)
c   Coulomb-type terms
      if(jlmn.le.noctot) then
        ytwo2=ytwo2+g(jlmn)*cc(ikmn)
      end if
      if(ikmn.le.noctot) then
        ytwo2=ytwo2+g(ikmn)*cc(jlmn)
      end if
      if(jkmn.le.noctot) then
        ytwo2=ytwo2-g(jkmn)*cc(ilmn)
      end if
      if(ilmn.le.noctot) then
        ytwo2=ytwo2-g(ilmn)*cc(jkmn)
      end if
c   Exchange-type terms
      if(jmln.le.noctot) then
        ytwo4=ytwo4+g(jmln)*cc(imkn)
      end if
      if(imkn.le.noctot) then
        ytwo4=ytwo4+g(imkn)*cc(jmln)
      end if
      if(jmkn.le.noctot) then
        ytwo4=ytwo4-g(jmkn)*cc(imln)
      end if
      if(imln.le.noctot) then
        ytwo4=ytwo4-g(imln)*cc(jmkn)
      end if
  101 CONTINUE
C
      ytwo= ytwo2*two + ytwo4*four
      y4mat= yone + ytwo
c
      RETURN                                                            SMM17500
      END                          
