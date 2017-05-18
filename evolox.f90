!-----------------------------------------------------------------------------
!--                          Subroutine evolox                             ---
!-----------------------------------------------------------------------------
!
subroutine evolox(deltatx,epsx,epsxerr)

use Para_DerivnD
use deriva ! (icon,npd,dxden,dyden,dzden)
use field  ! (pot4,hpsi,uext,limp)
use grid   ! (nx,ny,nz,nxyz,dxyz)
use gridk  ! (px,py,pz)
use he4    ! (h2o2m4)
use impur
use rho    ! (*psi*)
use util1  ! (vdt,nn,mmx,iw)
use work1  ! (temporal storage)

implicit none

real    (kind=8) :: deltatx,epsx,epsxerr

integer (kind=4) :: ix,iy,iz
real    (kind=8) :: hmean
real    (kind=8) :: cnorm

!........................
!.. Laplacian of Psix ...
!........................

!
!   icon   =  0 ! Take the derivative.
!   icon   = 20 ! Take the derivative. Use Wigner-Seitz conditions.
!   icon   =  8 ! Take the derivative. Use periodic conditions.

!call pderg(2,npd,nn,hx,psix,sto1,3,1,mmx,iw,icon)
!call pderg(2,npd,nn,hy,psix,sto2,3,2,mmx,iw,icon)
!call pderg(2,npd,nn,hz,psix,sto3,3,3,mmx,iw,icon)

  Call derivnD(2,nn,hx,1,psix,sto1,Icon)
  Call derivnD(2,nn,hy,2,psix,sto2,Icon)
  Call derivnD(2,nn,hz,3,psix,sto3,Icon)

forall(ix=1:nx,iy=1:ny,iz=1:nz)
  sto4(ix,iy,iz)  = (sto1(ix,iy,iz)+sto2(ix,iy,iz)+sto3(ix,iy,iz))*h2o2mx
  hpsix(ix,iy,iz) = -sto4(ix,iy,iz)+upotx(ix,iy,iz)*psix(ix,iy,iz)
end forall

hmean   = sum(psix*hpsix)*dxyz     ! Average over all the epsx

forall(ix=1:nx,iy=1:ny,iz=1:nz)
  hpsix(ix,iy,iz) = hpsix(ix,iy,iz)-epsx*psix(ix,iy,iz) ! (T+U)*Psi - eps_x*Psi_x = (H-eps_x*Psi
end forall

epsxerr = abs(1.0d0-epsx/hmean)    ! 'Error' in epsx
epsx    = hmean                    ! New epsilon autovalue

if(ironx.ne.0) then
  do ix=1,ironx
    call ironing(hpsix,nx,ny,nz)       ! Smoothing  (H-epsx)*Psix
  end do
end if

forall(ix=1:nx,iy=1:ny,iz=1:nz)
  psixnw(ix,iy,iz) = psix(ix,iy,iz) - deltatx*hpsix(ix,iy,iz)    &
                    + vdt(1)*(psix(ix,iy,iz)  - psi1x(ix,iy,iz)) &
                    + vdt(2)*(psi1x(ix,iy,iz) - psi2x(ix,iy,iz))
  psi2x(ix,iy,iz)  = psi1x(ix,iy,iz)
  psi1x(ix,iy,iz)  = psix(ix,iy,iz)
end forall

do iz=1,nz
  do iy=1,ny
    do ix=1,nx
      psixnw(ix,iy,iz)= max(psimin,psixnw(ix,iy,iz))
    end do
  end do
end do
cnorm  = sqrt(1.0d0/(sum(psixnw**2)*dxyz))

forall(ix=1:nx,iy=1:ny,iz=1:nz)
  psix(ix,iy,iz) = psixnw(ix,iy,iz)*cnorm
  denx(ix,iy,iz) = psix(ix,iy,iz)**2
end forall

!
return
end
