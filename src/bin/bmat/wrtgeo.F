      subroutine wrtgeo(xy,geobf,na,nad3,numit,igmup,iout,ang,
     1 dipno,ndisp,lbls,nsym,ccin,indi)
      implicit integer(a-z)
      integer indi(4,ndisp)
      real*8 xy(nad3),ang,x,y,z,ccin(4,ndisp)
      character*80 geobf(2000),lbls(nsym)
      character*34 geo1
      character*12 geo2

      geo1='% updated geometry from iteration'
      geo2='geometry = ('
      napdm=nad3/3

      rewind(igmup)

      do 10 i=1,2000
  10  read(igmup,100,end=11) geobf(i)

  11  continue

      nlines = i-1
      if(nlines.eq.0) geobf(1)='%%'

      if(nlines.eq.1999) then
        write(iout,*) 'too many lines in geom.dat.'
        write(iout,*) 'some information may be lost'
      endif

      if(ndisp.gt.0) then
        if(dipno.eq.1.and.nlines.eq.0) write(igmup,106)
        write(igmup,105)ccin(1,dipno),lbls(indi(1,dipno)),indi(1,dipno)
  105 format('% displacement by',f8.5,' of the ',a20,' coordinate',
     +  ' (',i5,')')
  106 format('%%')
      else
        rewind(igmup)
        write(igmup,100) geobf(1)
        write(igmup,101) geo1,numit+1
  101 format(a34,i5)
      endif

      write(igmup,102) geo2
  102 format(a12)

      do 20 i=1,na
        j=i+napdm
        j1=j+napdm
        x=xy(i)*ang
        y=xy(j)*ang
        z=xy(j1)*ang
        write(igmup,103) x,y,z
        write(iout,104) i,x,y,z
   20 continue
      write(igmup,*)')'

      if(ndisp.eq.0) then
        do 30 i=2,nlines
   30   write(igmup,100) geobf(i)
      endif
 
      return

  100 format(a80)
  103 format('(',3f20.10,')')
  104 format(1x,i3,3f16.12)
      end
