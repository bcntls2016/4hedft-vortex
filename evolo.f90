!-----------------------------------------------------------------------------
!--                          Subroutine evolo                              ---
!-----------------------------------------------------------------------------
!
subroutine evolo(deltat,mu4,mu4err,n4)

Use Para_DerivnD
use deriva ! (icon,npd,dxden,dyden,dzden,pderx,pdery,pderz,llap,xlap)
use vortex ! (to compute with a vortex)
use field  ! (pot4,hpsi,uext,limp)
use grid   ! (nx,ny,nz,nxyz,dxyz)
use gridk  ! (px,py,pz)
use he4    ! (h2o2m4)
use impur  !
use rho    ! (*psi*)
use util1  ! (vdt,nn,mmx,iw)
use work1  ! (temporal storage)

implicit none

integer (kind=4) :: n4
real    (kind=8) :: deltat,mu,muerr,n4real
real    (kind=8) :: mu4,mu4err

integer (kind=4) :: ix,iy,iz
real    (kind=8) :: hmean,a1
complex    (kind=8) :: caux,ci=(0.d0,1.0d0)
real    (kind=8) :: cnorm,sto,aux,dlaux


dlaux=dlog(2.0d0)


!.......................
!.. Laplacian of Psi ...
!.......................

!
!   icon   =  0 ! Take the derivative.
!   icon   = 20 ! Take the derivative. Use Wigner-Seitz conditions.
!   icon   =  8 ! Take the derivative. Use periodic conditions.
   
  Call derivnD(2,nn,hx,1,psi,sto1c,Icon)
  Call derivnD(2,nn,hy,2,psi,sto2c,Icon)
  Call derivnD(2,nn,hz,3,psi,sto3c,Icon)

   sto4c = -(sto1c+sto2c+sto3c)*h2o2m4 + pot4*psi
!
!     Here we will compute -Omega*|L_axis Psi>
!
      If(Lvortex.And.nv.Gt.1)Then
        caux = (0.d0, 0.d0)      
        If(Vortex_axis.Eq.'Z')Then
          Call derivnD(1,nn,hx,1,psi,sto1c,Icon)
          Call derivnD(1,nn,hy,2,psi,sto2c,Icon)
          Do iz=1, nz
            Do iy=1, ny
              Do ix=1, nx
                sto5c(ix,iy,iz) = Ci*(y(iy)*sto1c(ix,iy,iz) - x(ix)*sto2c(ix,iy,iz)) 
              EndDo
            EndDo
          EndDo
        Endif
        If(Vortex_axis.Eq.'Y')Then
          Call derivnD(1,nn,hx,1,psi,sto1c,Icon)
          Call derivnD(1,nn,hz,3,psi,sto2c,Icon)
          Do iz=1, nz
            Do iy=1, ny
              Do ix=1, nx
                sto5c(ix,iy,iz) = Ci*(x(ix)*sto2c(ix,iy,iz) - z(iz)*sto1c(ix,iy,iz)) 
              EndDo
            EndDo
          EndDo
        Endif
        If(Vortex_axis.Eq.'X')Then
          Call derivnD(1,nn,hy,2,psi,sto1c,Icon)
          Call derivnD(1,nn,hz,3,psi,sto2c,Icon)
          Do iz=1, nz
            Do iy=1, ny
              Do ix=1, nx
                sto5c(ix,iy,iz) = Ci*(z(iz)*sto1c(ix,iy,iz) - y(iy)*sto2c(ix,iy,iz)) 
              EndDo
            EndDo
          EndDo
        Endif
        sto4c = sto4c -omega*sto5c
      Endif        
!.................................................................. hpsi  = (T+U)*Psi_4 (With impurity)

   hpsi = sto4c
   hmean  = sum(conjg(psi)*hpsi)/sum(den)    ! Average over all the mu4's
   mu4err = abs(1.0d0-mu4/hmean)             ! Relative error
  
!   Write(6,'("From evolo: hmean, mu4, mu4err....",1p,3E15.6)')hmean,mu4,mu4err

   a1=1.0d0-deltat*(hmean-mu4)
   a1=deltat/a1

   If(.Not.Lbulk)Then
     mu4    = hmean
   Endif

   hpsi =a1*(hpsi-mu4*psi)

   if(iron.ne.0) then
     do ix=1,iron
       call ironingc(hpsi,nx,ny,nz)       ! Smoothing  (H-mu)*Psi
     end do
   end if

   forall(ix=1:nx,iy=1:ny,iz=1:nz)
     psinw(ix,iy,iz) = psi(ix,iy,iz)-  hpsi(ix,iy,iz)            &
                     + vdt(1) * (psi(ix,iy,iz)-psi1(ix,iy,iz))   &
                     + vdt(2) * (psi1(ix,iy,iz)-psi2(ix,iy,iz))
     psi2(ix,iy,iz)  = psi1(ix,iy,iz)
     psi1(ix,iy,iz)  = psi(ix,iy,iz)
   end forall

   hpsi = sto4c

!   psinw=max(psimin,Abs(psinw))

   If(Lbulk)Then
     cnorm=1.0d0
     n4real=sum(den)*dxyz
     n4=n4real+0.5
   Else
     cnorm = sqrt(n4/(sum(conjg(psinw)*psinw)*dxyz))
   Endif

   forall(ix=1:nx,iy=1:ny,iz=1:nz)
     psi(ix,iy,iz) = psinw(ix,iy,iz)*cnorm
     den(ix,iy,iz) = conjg(psi(ix,iy,iz))*psi(ix,iy,iz)
   end forall

return
end
