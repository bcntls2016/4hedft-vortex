      SUBROUTINE STEPRKR(deltat,mu4,mu4err,n4)
!
!      Runge-Kutta-Gill method
!       (Ralston & Wilf Vol I, pag. 117)
!
use Para_derivnD
use deriva ! (icon,npd,dxden,dyden,dzden,pderx,pdery,pderz,llap,xlap)
use field  ! (pot4,hpsi,uext,limp)
use grid   ! (nx,ny,nz,nxyz,dxyz)
use gridk  ! (px,py,pz)
use he4    ! (h2o2m4)
use rho    ! (psi, psiold & hpsiold)
use util1  ! (vdt,nn,mmx,iw)
use work1  ! (temporal storage)
use rkpc   ! (Storage for Steprk & Steppc rutines)

implicit none

real (kind=8) :: arun(4),brun(4),crun(4)

integer (kind=4) :: ix,iy,iz,jrun,n4
real    (kind=8) :: deltat,mu4,mu4err
real    (kind=8) :: hmean,cnorm

arun(1)=0.5d0
arun(2)=1.0d0-1.d0/dsqrt(2.d0)
arun(3)=1.0d0+1.d0/dsqrt(2.d0)
arun(4)=1.d0/6.d0

brun(1)=2.0d0
brun(2)=1.0d0
brun(3)=1.0d0
brun(4)=2.0d0

 crun(1)=0.5d0
 crun(2)=1.0d0-1.d0/dsqrt(2.d0)
 crun(3)=1.0d0+1.d0/dsqrt(2.d0)
 crun(4)=0.5d0



do jrun=1,4

  Call derivnD(2,nn,hx,1,psi,sto1c,Icon)
  Call derivnD(2,nn,hy,2,psi,sto2c,Icon)
  Call derivnD(2,nn,hz,3,psi,sto3c,Icon)

!
!   We compute H·Psi
!

hpsi   = -(sto1c+sto2c+sto3c)*h2o2m4 + pot4*psi
mu4err  = sum(Conjg(psi)*hpsi)*dxyz/n4    ! Average over all the mu4's
If(jrun.Eq.1)Then
  hmean  = mu4err
  timec=deltat/(1.0d0+dexp((Abs(hpsi/psi)-Dmumax)/sdmu))+dtmin
Endif
Sto4c   = -(hpsi-mu4err*psi)

Sto1c = arun(jrun)*(Sto4c - brun(jrun)*q)
q = q + 3.*Sto1c - crun(jrun)*Sto4c

if(jrun.eq.1)then
  hpsiold(:,:,:,2) = hpsiold(:,:,:,1)
  hpsiold(:,:,:,1) = Sto4c
  psiold(:,:,:,3) = psiold(:,:,:,2)
  psiold(:,:,:,2) = psiold(:,:,:,1)
  psiold(:,:,:,1) = psi
endif

psi = psi + timec*Sto1c

!If(psimin.Gt.0.0d0)  psi = Max(psimin,psi )

!psi = psi + deltat*Sto1

cnorm = sqrt(n4/(sum(Conjg(psi)*psi)*dxyz))
psi = psi*cnorm
den = Conjg(psi)*psi

  if(jrun.le.3)then
    call poten()
  endif
enddo

do ix=1,3
  ioldp(ix)=ix
enddo

do ix=1,2
  ioldh(ix)=ix
enddo
mu4err = abs(1.0d0-mu4/hmean)     ! 'Error' in mu4
mu4    = hmean                    ! New chemical potential
return
end
