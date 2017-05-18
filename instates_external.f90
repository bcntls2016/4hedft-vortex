subroutine instates_external()

use grid
use rho
use impur
use energies
use field

!implicit Real*8(A-H,O-Z)
implicit none

complex (kind=8) :: H_D(10,10) ! Matriz para los estados DD
real    (kind=8) :: ar(10,10), ai(10,10), w(10), zr(10,10), zi(10,10), fv1(10), fv2(10), fm1(2,10)
complex (kind=8) :: Mat(10,10), SO(10,10), MatSO(10,10)
REAL (kind=8) :: d0, d1, d2, d3, d4, d5, sd0, sd1, sd2, sd3, sd4, sd5, dV1, dV2
REAL (kind=8) :: A,B,C,Vx,Vy,Vz,Vxy,Vxz,Vyz,dVpi,dVsig,dVdel, dum2,Dv,vls,Ep1,Ep2,Ep3,V0gs,U_ext
REAL (kind=8) :: armonicd(6), rimp(3), pi, twopi, fpi, dx, dy, dz, zx, yx, xx, yy, zz, r, r2, c1osqr2
Integer (kind=4)  :: ierr,i,j,ix,iy,iz,il,nm,np,nd,il1, k, kk, l, ll, njz, ki
Integer (kind=4),save  :: iaux
REAL (kind=8) :: sqr3, sqr3o2, insqr3, in3, difd3d5, sumd3d5, cm1_to_K = 1.439d0, Also2, Uext_new
REAL (kind=8) :: sumVdel(10,10), MP(10,10), MS(10,10), MatD(10,10), U_total, Aux, Over=1.d4
COMPLEX (kind=8) :: SOD(10,10), MatDSO(10,10), eigv(10,10), testv(10,10), invars(10), vlsaux, SO_aux(10)
complex (kind=8) :: auxEVEN, auxODD, auxc
complex   (kind=8) :: xxxr(3)
complex   (kind=8) :: uim = (0.d0,1.d0)

Data nm/10/, np/6/, nd/10/, iaux/0/

pi    = 4.0D0*DATAN(1.0D0)
twopi = 2.d0*pi
fpi   = 4.d0*pi
c1osqr2=1.0d0/Dsqrt(2.0d0)

!iaux = iaux +1

!Write(6,'("1: (",I5,") ","Lstate....:",A)')iaux, Lstate


If (Lstate=='P') then

  Als = Als_P
  Also2 = Als_P*0.5d0

!......................................
!... Dimensionamos algunas matrices ...
!......................................

  rimp(1)=ximp; rimp(2)=yimp; rimp(3)=zimp

  dx   = hx
  dy   = hy
  dz   = hz

  zx = rimp(3)
  yx = rimp(2)
  xx = rimp(1)

  do iz=1,nz
    zz = z(iz) - rimp(3)
    do iy=1,ny
      yy = y(iy) - rimp(2)
      do ix=1,nx
        xx = x(ix) - rimp(1)
        r2=xx**2+yy**2+zz**2
        r=dsqrt(r2)
        If(r.ne.0.d0) then
          dVpi = Vpi(ix,iy,iz)
          dV=Delta(ix,iy,iz)*r2
          Vx  = (dVpi+dV*(xx/r)**2)
          Vy  = (dVpi+dV*(yy/r)**2)
          Vz  = (dVpi+dV*(zz/r)**2)
          Vxy = dV*xx*yy/r2
          Vxz = dV*xx*zz/r2
          Vyz = dV*yy*zz/r2
        Endif
!
!! polinomio caracteristico para encontrar los valores propios:
!! t^3 + C t^2 + B t + A == 0
!
  A = Vxz**2*Vy-2*Vxy*Vxz*Vyz+Vx*Vyz**2+Vxy**2*Vz-Vx*Vy*Vz
  B = -Vxy**2-Vxz**2-Vyz**2+Vx*Vy+Vx*Vz+Vy*Vz
  C = -(Vx+Vy+Vz)
!
!! adicion del termino de spin-orbita
!! Delta = 3/2 Als
  A = A + 0.25d0*(Als**3-C*Als**2)
  B = B - 0.75d0*Als**2
  C = C
!! valores propios:
  call cubic(1.d0,C,B,A,xxxr)
  Ep1=real(xxxr(1))
  Ep2=real(xxxr(2))
  Ep3=real(xxxr(3))
  If(r.Ne.0.0d0)Then
    If(Instate.Eq.0)Uext(ix,iy,iz)= Ep3
    If(Instate.Eq.1)Uext(ix,iy,iz)= Ep2
    If(Instate.Eq.2)Uext(ix,iy,iz)= Ep1
  Else
    Uext(ix,iy,iz)= Over
  Endif
      end do
    end do
  end do
  Eso = 0.0d0
Else
  Write(6,'("De moment instate_external, nomes funciona pels estats P")')
  Stop 'FromInstae_external  0001'
EndIf
return
end

