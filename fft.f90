!...................................................................
!...                Subroutine fftini                            ...
!...................................................................

subroutine fftini(nx,ny,nz)

use fftmodule  ! fin,fout,fftwplan,pfftfw,pfftbk,nthread,renor,npx,
               ! npy, npz
implicit none

integer   (kind=4) :: nx,ny,nz,iret    ! Size of the grid along X, Y and Z axis

include 'fftw3.f.include'

allocate(fin(nx,ny,nz)  )
allocate(fout(nx/2+1,ny,nz) )

npx   = nx
npy   = ny
npz   = nz
renor = 1.0d0/(nx*ny*nz)

call dfftw_init_threads(iret)
call dfftw_plan_with_nthreads(nthread)

call dfftw_plan_dft_r2c_3d(pfftfw,nx,ny,nz,fin ,fout,fftwplan)
call dfftw_plan_dft_c2r_3d(pfftbk,nx,ny,nz,fout,fin ,fftwplan)


return

end
!...................................................................
!...                Subroutine fftfw                             ...
!...................................................................

subroutine fftfw(a,b)

use fftmodule  ! fin,fout,fftwplan,pfftfw,pfftbk,nthread

implicit none

integer (kind=4) :: ix,iy,iz    ! Working variables.
real    (kind=8) :: a(npx,npy,npz)
complex (kind=8) :: b(npx/2+1,npy,npz)

forall(ix=1:npx, iy=1:npy, iz=1:npz)
  fin(ix,iy,iz) = a(ix,iy,iz)
end forall

call dfftw_execute(pfftfw)

forall(ix=1:npx/2+1, iy=1:npy, iz=1:npz)
  b(ix,iy,iz) = fout(ix,iy,iz)
end forall

return

end
!...................................................................
!...                Subroutine fftbk                             ...
!...................................................................

subroutine fftbk(c,d)

use fftmodule  ! fin,fout,fftwplan,pfftfw,pfftbk,nthread,renor

implicit none

integer (kind=4) :: ix,iy,iz    ! Working variables.
complex (kind=8) :: c(npx/2+1,npy,npz)
real    (kind=8) :: d(npx,npy,npz)



forall(ix=1:npx/2+1, iy=1:npy, iz=1:npz)
  fout(ix,iy,iz)=c(ix,iy,iz)
end forall

call dfftw_execute(pfftbk)

forall(ix=1:npx, iy=1:npy, iz=1:npz)
  d(ix,iy,iz) = fin(ix,iy,iz)*renor
end forall


return

end
