c**********************************************************************
c   SUBROUTINE MK_FORCE.F                                             *
C**********************************************************************
C   For all modes                                                     *
c**********************************************************************
c   Last updated on January 05, 2007 by Dr. Yukio Yamaguchi           *
c**********************************************************************
      SUBROUTINE MK_FORCE(energ,fxx)
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (max_cal=25)
      dimension energ(max_cal)
      dimension fxx(3,3)
      common/vib102/pertr,perta
      data a00,one,two / 0.0d+00 , 1.0d+00, 2.0d+00  /
      data HE / 4.359813653D+00 /
      data A0 / 0.529177249D+00 /
    1 format(/,1x,' freq      = ',f20.10/)
    2 format(/,1x,' freq      = ',f20.10,' i',/)
C
C===========================================
C===========================================
c   Hessian via 5-point formula
C      f(x0+2p),f(x0+1p),f(x0)
C===========================================
C===========================================
c  
      f11 = - energ(1)*(30.0d+00)
     *      + (energ(2) + energ(3))*(16.0d+00)
     *      - (energ(4) + energ(5))
c
      f22 = - energ(1)*(30.0d+00)
     *      + (energ(6) + energ(7))*(16.0d+00)
     *      - (energ(8) + energ(9))
c
      f33 = - energ(1)*(30.0d+00)
     *      + (energ(10) + energ(11))*(16.0d+00)
     *      - (energ(12) + energ(13))
c
      f12 =   energ(14) - energ(15)
     *      - energ(16) + energ(17)
c
      f13 =   energ(18) - energ(19)
     *      - energ(20) + energ(21)
c
      f23 =   energ(22) - energ(23)
     *      - energ(24) + energ(25)
c
      fxx(1,1) = f11
      fxx(2,2) = f22
      fxx(3,3) = f33
      fxx(1,2) = f12
      fxx(2,1) = f12
      fxx(1,3) = f13
      fxx(3,1) = f13
      fxx(2,3) = f23
      fxx(3,2) = f23
c
      write(6,*) '  '
      write(6,*) ' ::::::::::::::::::::::::::::::::::::::::::::::'
      write(6,*) ' :  Hessian with 5-point formula              :'
      write(6,*) ' :  f(x0+2p) f(x0+1p) f(x0)                   :'
      write(6,*) ' ::::::::::::::::::::::::::::::::::::::::::::::'
      write(6,*) '  '
      write(6,*) ' The hessian value with 5-point formula'
      write(6,*) ' [in hartree]'
      write(6,*) '  '
      call matout(fxx,3,3,3,3,6)
c
c   Convert the force constant in atomic units
C
C   For the stretch-stretch force constants
      d12sqh_r = pertr*pertr*(12.0d+00)
      d4sqh_rr = pertr*pertr*(4.0d+00)
      fxx(1,1) = fxx(1,1) / ( d12sqh_r )
      fxx(2,2) = fxx(2,2) / ( d12sqh_r )
      fxx(1,2) = fxx(1,2) / ( d4sqh_rr )
      fxx(2,1) = fxx(2,1) / ( d4sqh_rr )
C
C   For the bend-bend force constant
      d12sqh_a = perta*perta*(12.0d+00)
      fxx(3,3) = fxx(3,3) / ( d12sqh_a )
C
C   For the stretch-bend force constants
      d4sqh_ra = pertr*perta*(4.0d+00)
      fxx(1,3) = fxx(1,3) / ( d4sqh_ra )
      fxx(3,1) = fxx(3,1) / ( d4sqh_ra )
      fxx(2,3) = fxx(2,3) / ( d4sqh_ra )
      fxx(3,2) = fxx(3,2) / ( d4sqh_ra )
      write(6,*) '  '
      write(6,*) ' The hessian value with 5-point formula'
      write(6,*) ' [in hartree/(bohr*bohr)]'
      write(6,*) '  '
      call matout(fxx,3,3,3,3,6)
c
c   Convert the force constant in spectroscopic units
      rausq = one / (A0 * A0)
      rausg = one / (A0)
C
C   For the stretch-stretch force constants
      fxx(1,1) = fxx(1,1) * ( he ) * ( rausq )
      fxx(2,2) = fxx(2,2) * ( he ) * ( rausq )
      fxx(1,2) = fxx(1,2) * ( he ) * ( rausq )
      fxx(2,1) = fxx(2,1) * ( he ) * ( rausq )

C   For the bend-bend force constant
      fxx(3,3) = fxx(3,3) * ( he )

C   For the stretch-bend force constants
      fxx(1,3) = fxx(1,3) * ( he ) * ( rausg )
      fxx(3,1) = fxx(3,1) * ( he ) * ( rausg )
      fxx(2,3) = fxx(2,3) * ( he ) * ( rausg )
      fxx(3,2) = fxx(3,2) * ( he ) * ( rausg )
      write(6,*) '  '
      write(6,*) ' The hessian value with 5-point formula'
      write(6,*) ' [in spectroscopic units]'
      write(6,*) '  '
      call matout(fxx,3,3,3,3,6)
C
      return
      end
