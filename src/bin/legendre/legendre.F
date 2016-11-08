      subroutine fentry(cc,ia,maxcor)
c**********************************************************************
C     PROGRAM LEGENDRE                                                *
c     This program will test the Legendre polynomials                 *
c**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
c**********************************************************************
C   Last updated on April 18, 2003 by Dr. Yukio Yamaguchi             *
c**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      character*2 kadjst
      data a00,one,two / 0.0d+00 , 1.0d+00 , 2.0d+00 /
      data pi / 3.1415926536d+00 /
c==================================================================
c  For psi 2.0 input
#include <error.h>
      integer frdc, frdi, frdf, fwkclr, fwkadd, frdboo, frdpos
      integer fcnt, frdia1
      integer errcod
c==================================================================
    1 FORMAT(//,2X,' THE LEGENDRE POLYNOMIALS PROGRAM'/
     1          2x,'     (Version of 04/18/2003)'/)
    2 FORMAT(//,2X,' PARAMETERS'/
     1 2X,' IPRNT     = ',I5/ 
     2 2X,' K_POWER   = ',I5/ 
     3 2X,' K_ADJUST  = ',A2/
     4 2X,' P_VALU    = ',f15.10/)
    3 FORMAT(/,2X,' Legendre polynomials'/
     1          2x,'   K   kvalu',
     2         10x,'pvalu',
     3         12x,'pavalu',
     4         11x,'pbvalu'/)
    4 FORMAT(3x,i3,5x,i3,5x,f12.7,5x,f12.7,5x,f12.7)
C
      call psinit('APPEND')
      CALL TSTART(6)                                                    FDG00550
      WRITE(6,1)                                                        FDG00630
      write(6,*) '  maxcor    = ',maxcor
C
      errcod    = fwkclr()
      errcod    = fwkadd(5,'LEGENDRE')
      errcod    = fwkadd(5,'DEFAULT')
      errcod    = frdi('PRINT',iprnt)
      ppk       = a00
      errcod    = frdf('P_VALU',ppk)
      kpower    = 1
      errcod    = frdi('K_POWER',kpower)
      kadjst    = '2K'
      errcod    = frdc('K_ADJUST',kadjst)
      if(kadjst.eq.'2K') then
        k_pref  = 2
      end if
      if(kadjst.eq.'1K') then
        k_pref  = 1
      end if
C
      WRITE(6,2) iprnt,kpower,kadjst,ppk
C
      write(6,3)
      do 101 k = 1, kpower+1
      kvalu  = k_pref*(k-1)
      pvalu  = ppk
      pavalu = pcala(ppk,kvalu)
      pbvalu = pcalb(ppk,kvalu)
      write(6,4) k, kvalu, pvalu, pavalu, pbvalu
  101 continue
      write(6,*) '  '
C
      CALL TSTOP(6)                                                     FDG01730
      STOP                                                              FDG01750
      END                                                               FDG01760
      function pcala(pvalu,lpower)
      implicit real*8 (a-h,o-z)
      dimension q(100)
      data one, two / 1.0d+00 , 2.0d+00 /
    1 format(2x,i5,2f20.10)
C
C  Form the Legendre polynomials
      q(1) = one
      if(lpower.eq.0) go to 200
C
      zvalu = dcos(pvalu)
      q(2)  = zvalu
      if(lpower.eq.1) go to 200
C
      llim = lpower-1
      do 101 l = 1, llim
      rl     = dfloat(l)
      fact   = two*rl + one
      denom  = rl + one
      valu   = fact*zvalu*q(l+1) - rl*q(l)
      q(l+2) = valu / ( denom )
  101 continue
C
  200 continue
      pcala = q(lpower+1)
C
      return
      end
      function pcalb(pvalu,lpower)
      implicit real*8 (a-h,o-z)
      dimension q(100)
      data a00, half / 0.0d+00 , 0.5d+00 /
      data one, two  / 1.0d+00 , 2.0d+00 /
    1 format(2x,i5,2f20.10)
C
C  Form the Legendre polynomials
      go to (201,202,203,204,205,206,207,208,209,210,
     *       211,212,213), lpower+1
      go to 215
C  p_{0} = q(1)
  201 continue
      q(1)   = one
      if(lpower.eq.0) go to 300
C
C  p_{1} = q(2)
  202 continue
      zvalu  = dcos(pvalu)
      q(2)   = zvalu
      if(lpower.eq.1) go to 300
C
C  p_{2} = q(3)
  203 continue
      zvalu  = dcos(pvalu)
      zval2  = zvalu * zvalu
      q(3)   = ( 3.0d+00 * zval2 - one ) * half
      if(lpower.eq.2) go to 300
C
C  p_{3} = q(4)
  204 continue
      zvalu  = dcos(pvalu)
      zval2  = zvalu * zvalu
      zval3  = zval2 * zvalu
      q(4)   = ( 5.0d+00 * zval3 - 3.0d+00 * zvalu ) * half
      if(lpower.eq.3) go to 300
C
C  p_{4} = q(5)
  205 continue
      zvalu  = dcos(pvalu)
      zval2  = zvalu * zvalu
      zval4  = zval2 * zval2
      q(5)   = ( 35.0d+00 * zval4 - 30.0d+00 * zval2
     *         +  3.0d+00 ) * 0.125d+00
      if(lpower.eq.4) go to 300
C
C  p_{5} = q(6)
  206 continue
      zvalu  = dcos(pvalu)
      zval2  = zvalu * zvalu
      zval3  = zval2 * zvalu
      zval5  = zval3 * zval2
      q(6)   = ( 63.0d+00 * zval5 - 70.0d+00 * zval3
     *         + 15.0d+00 * zvalu ) * 0.125d+00
      if(lpower.eq.5) go to 300
C
C  p_{6} = q(7)
  207 continue
      zvalu  = dcos(pvalu)
      zval2  = zvalu * zvalu
      zval4  = zval2 * zval2
      zval6  = zval4 * zval2
      q(7)   = ( 231.0d+00 * zval6 - 315.0d+00 * zval4
     *         + 105.0d+00 * zval2 -   5.0d+00 ) * 0.0625d+00
      if(lpower.eq.6) go to 300
C
C  p_{7} = q(8)
  208 continue
      zvalu  = dcos(pvalu)
      zval2  = zvalu * zvalu
      zval3  = zval2 * zvalu
      zval5  = zval3 * zval2
      zval7  = zval5 * zval2
      q(8)   = ( 429.0d+00 * zval7 - 693.0d+00 * zval5
     *         + 315.0d+00 * zval3 -  35.0d+00 * zvalu ) * 0.0625d+00
      if(lpower.eq.7) go to 300
C
C  p_{8} = q(9)
  209 continue
      zvalu  = dcos(pvalu)
      zval2  = zvalu * zvalu
      zval4  = zval2 * zval2
      zval6  = zval4 * zval2
      zval8  = zval4 * zval4
      q(9)   = ( 6435.0d+00 * zval8 - 12012.0d+00 * zval6
     *         + 6930.0d+00 * zval4 -  1260.0d+00 * zval2
     *         +   35.0d+00 ) * 0.0078125d+00
      if(lpower.eq.8) go to 300
C
C  p_{9} = q(10)
  210 continue
      zvalu  = dcos(pvalu)
      zval2  = zvalu * zvalu
      zval3  = zval2 * zvalu
      zval5  = zval3 * zval2
      zval7  = zval5 * zval2
      zval9  = zval7 * zval2
      q(10)  = ( 12155.0d+00 * zval9 - 25740.0d+00 * zval7
     *         + 18018.0d+00 * zval5 -  4620.0d+00 * zval3
     *         +   315.0d+00 * zvalu ) * 0.0078125d+00
      if(lpower.eq.9) go to 300
C
C  p_{10} = q(11)
  211 continue
      zvalu  = dcos(pvalu)
      zval2  = zvalu * zvalu
      zval4  = zval2 * zval2
      zval6  = zval4 * zval2
      zval8  = zval4 * zval4
      zval10 = zval8 * zval2
      q(11)  = ( 46189.0d+00 * zval10 - 109395.0d+00 * zval8
     *         + 90090.0d+00 * zval6  -  30030.0d+00 * zval4
     *         +  3465.0d+00 * zval2  -     63.0d+00 )
     *         * 0.00390625d+00
      if(lpower.eq.10) go to 300
C
C  p_{11} = q(12)
  212 continue
      zvalu  = dcos(pvalu)
      zval2  = zvalu * zvalu
      zval3  = zval2 * zvalu
      zval5  = zval3 * zval2
      zval7  = zval5 * zval2
      zval9  = zval7 * zval2
      zval11 = zval9 * zval2
      q(12)  = (  88179.0d+00 * zval11 - 230945.0d+00 * zval9
     *         + 218790.0d+00 * zval7  -  90090.0d+00 * zval5
     *         +  15015.0d+00 * zval3  -    693.0d+00 * zvalu )
     *         * 0.00390625d+00
      if(lpower.eq.11) go to 300
C
C  p_{12} = q(13)
  213 continue
      zvalu  = dcos(pvalu)
      zval2  = zvalu * zvalu
      zval4  = zval2 * zval2
      zval6  = zval4 * zval2
      zval8  = zval4 * zval4
      zval10 = zval8 * zval2
      zval12 = zval8 * zval4
      q(13)  = (  676039.0d+00 * zval12 - 1939938.0d+00 * zval10
     *         + 2078505.0d+00 * zval8  - 1021020.0d+00 * zval6
     *         +  225225.0d+00 * zval4  -   18018.0d+00 * zval2   
     *         +     231.0d+00 ) / (1024.0d+00)
      if(lpower.eq.12) go to 300
C
  215 continue
      q(lpower+1) = a00
C
  300 continue
      pcalb = q(lpower+1)
C
      return
      end
