!---------------------------------------------------------------------
!---                             Subroutine Poten                  ---
!---------------------------------------------------------------------

subroutine poten()

use alphasterm ! (ualphas)
use grid     ! (nxyz)
use he4      ! (cp4,cpp4)
use impur    ! (vq,)
use lenard4  ! (wk2,pelj4,fvlj4,delj4,core4,lalphas)
use field    ! (pot4,limp)
use rho      ! (dencg,den,wcgk,)
use work1    ! (sto1,sto2,sto3,sto4)
use constraint

implicit none

real    (kind=8) :: a0,a1
real    (kind=8) :: zcom
integer (kind=4) :: ix,iy,iz

!..............................
! Lennard-Jones contribution  .
!..............................

forall(ix=1:nx/2+1,iy=1:ny,iz=1:nz)
   wk1(ix,iy,iz) = fden(ix,iy,iz)*fvlj4(ix,iy,iz)
end forall 
call fftbk(wk1,delj4) ! Get delj4 -> (   int{ rho_4*V_4 dr'}  )

!........................
! Correlation terms.  ...
!........................

a0 = cp4 /2.d0               ! Auxiliar variable useful for saving operations
a1 = cpp4/3.d0               ! Auxiliar variable useful for saving operations

!   The rest of the correlation contribution is calculated as
!   a convolution product.

forall(ix=1:nx,iy=1:ny,iz=1:nz)
   sto1(ix,iy,iz) = den(ix,iy,iz)*dencg(ix,iy,iz)*(cp4+cpp4*dencg(ix,iy,iz))
end forall 

call fftfw(sto1,wk1)

forall(ix=1:nx/2+1,iy=1:ny,iz=1:nz)
   wk2(ix,iy,iz) = wk1(ix,iy,iz)*wcgk(ix,iy,iz)
end forall 

call fftbk(wk2,sto1)

!..........................
!... Solid penalty term ...
!..........................
If(lsolid)Then
  penalty = dtanh(beta*(den-den_m))
  penalty = C*(1.d0 + penalty + beta*den*(1.d0 - penalty**2) )
Endif

If(Lexcite_state)Then
   If(.Not.Lexcite_state_external)Then
     Call instates()
  Endif
Endif


!.......................
!.. Final 'Potential' ...
!........................

if(lexternal)then
 forall(ix=1:nx,iy=1:ny,iz=1:nz)
    pot4(ix,iy,iz) = delj4(ix,iy,iz) +                                 &   ! Lennard-Jones
                     dencg(ix,iy,iz)**2*(a0+a1*dencg(ix,iy,iz)) +      &   ! Correlation
                     sto1(ix,iy,iz) + uext(ix,iy,iz)                       ! Correlation
 end forall
else
 forall(ix=1:nx,iy=1:ny,iz=1:nz)
    pot4(ix,iy,iz) = delj4(ix,iy,iz) +                                 &   ! Lennard-Jones
                     dencg(ix,iy,iz)**2*(a0+a1*dencg(ix,iy,iz)) +      &   ! Correlation
                     sto1(ix,iy,iz)                                        ! Correlation
 end forall
endif
if(core4.eq.'OTC') then
    call term_alfa()    ! Calculates alpha_s contribution to the field
    forall(ix=1:nx,iy=1:ny,iz=1:nz)
       pot4(ix,iy,iz) = pot4(ix,iy,iz)  +                                 &   ! Mean field...
                        ualphas(ix,iy,iz)                                     ! Alfa_s term
    end forall
end if

if(limp) then
   forall(ix=1:nx/2+1,iy=1:ny,iz=1:nz)
      wk1(ix,iy,iz) = fden(ix,iy,iz) *vq(ix,iy,iz)  ! FFT(den_4)*FFT(V_x)
      wk2(ix,iy,iz) = fdenx(ix,iy,iz)*vq(ix,iy,iz)  ! FFT(den_x)*FFT(V_x)
   end forall 

   call fftbk(wk1,upotx) ! Get  ( int{ rho_4    * V_X dr'} ) Potential due to 4He for the impurity
   call fftbk(wk2,potx4) ! Get  ( int{ Psi_x**2 * V_X dr'} ) Potential due to the impurity for the 4He

  forall(ix=1:nx,iy=1:ny,iz=1:nz)
     pot4(ix,iy,iz) = pot4(ix,iy,iz) + potx4(ix,iy,iz)
  end forall 

end if



If(lsolid)Then
  forall(ix=1:nx,iy=1:ny,iz=1:nz)
    pot4(ix,iy,iz) = pot4(ix,iy,iz) + penalty(ix,iy,iz)
  end forall 
Endif

If(lconstraint)Then

! INTRODUCE CONSTRAINT:

  zcom = 0.d0
  do iz=1,nz
   zcom = zcom + sum(den(:,:,iz))*z(iz)
  enddo
  zcom = zcom*dxyz/n4real
! print*,'zcom=',zcom

  If(limp)Then
    zimp = 0.d0
    do iz=1,nz
      zimp = zimp + sum(denx(:,:,iz))*z(iz)
    enddo
    zimp = zimp*dxyz
  Endif
  If(limp.And.Lconstraint_imp)Then
    forall(ix=1:nx,iy=1:ny,iz=1:nz)
       upotx(ix,iy,iz) = upotx(ix,iy,iz)+Intens*(zimp-zcom-zdist)*z(iz)
    end forall
  Endif

  forall(ix=1:nx,iy=1:ny,iz=1:nz)
     pot4(ix,iy,iz) = pot4(ix,iy,iz)+Intens*(zimp-zcom-zdist)*(-z(iz)/n4real)
  end forall 
Endif

return

end
