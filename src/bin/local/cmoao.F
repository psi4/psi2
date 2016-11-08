      subroutine cmoao(nsymhf,num,msfru,nc,no,nlamda,soao,su,cmo,mxnao,
     & iopen)
c  change MO into the form in terms of AO's, instead of SO's.
      implicit real*8 (a-h,o-z)
      dimension nc(1),no(1),nlamda(1),su(1)
      dimension soao(num,msfru),cmo(mxnao,mxnao)

c First do the occupied orbitals
      kk=0
      do 150 nblc=1,nsymhf
      i1=1
      if (nblc.ne.1) i1=i2+1
      i2=i1-1+nc(nblc)
      k1=1
      if (nblc.ne.1) k1=k2+1
      k2=k1-1+nlamda(nblc)
      jk=kk-nlamda(nblc)
      kk=kk+nlamda(nblc)*nlamda(nblc)
      do 140 i=i1,i2
      jk=jk+nlamda(nblc)
      do 130 j=1,num
      cmo(i,j)=0.0
      ik=jk
      do 120 k=k1,k2
      ik=ik+1
      cmo(i,j)=cmo(i,j)+soao(j,k)*su(ik)
120   continue
130   continue
140   continue
150   continue

      if(iopen.eq.0) goto 160
c Next do the single occupied orbitals
      kk=0
      do 152 nblc=1,nsymhf
      i1=i2+1
      i2=i1-1+no(nblc)
      k1=1
      if (nblc.ne.1) k1=k2+1
      k2=k1-1+nlamda(nblc)
      jk=kk-nlamda(nblc) + nlamda(nblc)*nc(nblc)
      kk=kk+nlamda(nblc)*nlamda(nblc)
      if(no(nblc).eq.0) goto 152
      do 142 i=i1,i2
      jk=jk+nlamda(nblc)
      do 132 j=1,num
      cmo(i,j)=0.0
      ik=jk
      do 122 k=k1,k2
      ik=ik+1
      cmo(i,j)=cmo(i,j)+soao(j,k)*su(ik)
122   continue
132   continue
142   continue
152   continue
c
c Then do the virtual orbitals
160   kk=0
      do 151 nblc=1,nsymhf
      i1=i2+1
      i2=i1-1+nlamda(nblc)-nc(nblc)-no(nblc)
      k1=1
      if (nblc.ne.1) k1=k2+1
      k2=k1-1+nlamda(nblc)
clj   jk=kk-nlamda(nblc)
      jk=kk-nlamda(nblc)+nlamda(nblc)*nc(nblc)+nlamda(nblc)*no(nblc)
      kk=kk+nlamda(nblc)*nlamda(nblc)
      do 141 i=i1,i2
      jk=jk+nlamda(nblc)
      do 131 j=1,num
      cmo(i,j)=0.0
      ik=jk
      do 121 k=k1,k2
      ik=ik+1
      cmo(i,j)=cmo(i,j)+soao(j,k)*su(ik)
121   continue
131   continue
141   continue
151   continue
c
      return
c
      end
