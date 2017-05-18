!.......................................................................
!...                  Subroutine RHOASIN                             ...
!.......................................................................
subroutine rhoasin0(P,rho,mu,ff)

use he4

real    (kind=8), intent(out) :: p
real    (kind=8), intent(out) :: ff
real    (kind=8), intent(in)  :: rho
real    (kind=8), intent(out) :: mu

  ff  =f(rho)
  p   = df(rho) *rho - ff

  mu  = df(rho)

return

contains
!-------------------------------------------------------------------
!---            Function f 
!-------------------------------------------------------------------
!
!... Free energy for Barcelona-BBAA funtional

function f(rho)

use he4
use lenard4

real (kind=8) ,intent(in) :: rho
real (kind=8)             :: f
real (kind=8)             :: rho2,rho3,rho4,aux1
real (kind=8), parameter  :: c1=0.5d0
real (kind=8), parameter  :: c2=0.5d0
real (kind=8), parameter  :: c3=1.0d0/3.d0
Integer (Kind=4) :: i1

rho2 = rho *rho
rho3 = rho2*rho 
rho4 = rho2*rho2

!b  = bforce(core4,eps4,sigma4,h4)    ! Parameter 'b' calculated from LJ potential
f  = c1*b4*rho2 + c2*cp4*rho3 + c3*cpp4*rho4 

return
end function f
!-------------------------------------------------------------------
!---            Function df
!-------------------------------------------------------------------
!
!... First Derivative of Free energy for Barcelona-BBAA funtional
function df(rho)

use he4
use lenard4

real (kind=8) ,intent(in) :: rho
real (kind=8)             :: df
real (kind=8)             :: rho2,rho3
real (kind=8),parameter   :: c1=3.d0/2.0d0
real (kind=8),parameter   :: c2=4.d0/3.0d0
Integer (Kind=4) :: i1

rho2 = rho *rho
rho3 = rho2*rho

df = b4*rho + c1*cp4*rho2 + c2*cpp4*rho3 

return

end function df
!-------------------------------------------------------------------
!---            Function ddf
!-------------------------------------------------------------------
!
!... Second Derivative of Free energy for Barcelona-BBAA funtional
function ddf(rho)

use he4
use lenard4

real (kind=8) ,intent(in) :: rho
real (kind=8)             :: ddf

ddf = b4 + 3.d0*cp4*rho + 4.d0*cpp4*rho*rho

return

end function ddf


end subroutine rhoasin0
