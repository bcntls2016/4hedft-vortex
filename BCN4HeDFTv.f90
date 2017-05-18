program BCN4HeDFT
!--------------------------------------------------------------------------------------
!        
!                This code computes the structure and energetics of pure or doped
!                helium-4 droplets with a vortex within Density Functional Theory.
!
!                The functional can be Orsay-Paris [Dupont-Roc et al, JLTP 81, 31
!                (1990)], Orsay-Trento with alpha_s term but no backflow term -
!                static conditions- [Dalfovo et al, PRB 52, 1193 (1995)] or the solid-
!                like functional by Ancilotto et al, PRB 72, 214522 (2005).
!                
!                The impurity can be treated either as a quantum particle by solving
!                the corresponding Schroedinger equation or as an external field.
!                Constrained calculations on the impurity location are also possible.
!                Handling bulk liquid helium or liquid helium free surfaces under the
!                above conditions is also possible.
!
!                The code is fully 3D and works in Cartesian coordinates
!
!                It can be freely used. Please, refer to it as:
!                BCN-4He-DFT code for helium-4 droplets v.1 (2015)
!                M. Barranco, A. Hernando, D. Mateo, R. Mayol, and M. Pi.
!
!                University of Barcelona.
!
!                Comments and suggestions welcome (marti@ecm.ub.edu)
!
!                         *****************
!                         ** BCN4HeDFT **
!                         *****************
!
!
!--------------------------------------------------------------------------------------
Use For_FT_V_spline
Use seleccio_de_potencial
Use Modificacio_De_Select_Pot
Use Para_DerivnD
use alphasterm
use deriva
use vortex
use energies
use rho
use field
use fftmodule
use grid 
use gridk
use impur
use lenard4
use he4
use util1
use work1
use constraint
use rkpc
implicit none

real (kind=4) :: t0,t1,t2,t3,t4,t5,t6   ! Variables used to printout run time execution.

logical              :: lfilepv         ! T-> print save file when change Paflo parameter.
logical              :: lpaflv=.false.  ! T-> allows change of Paflov coeffient
integer    (kind=4)  :: N_Newton=50     ! Nombre maxim de iterations pel Newton-Rapson
integer    (kind=4)  :: ndmax=2         ! maxima derivada a calcular
integer    (kind=4)  :: naux            ! Auxiliar variable
integer    (kind=4)  :: it,j
integer    (kind=4)  :: nstepp=1        ! Number of 'Paflov parameter'.
integer    (kind=4)  :: ix ,iy ,iz      ! Used as indexs of arrays
integer    (kind=4)  :: ipx,ipy,ipz     ! Used as indexs of arrays
integer    (kind=4)  :: iter,niter      ! Control the number of iterations
integer    (kind=4)  :: pener=50        ! Computes energy each 'pener' iterations
integer    (kind=4)  :: pdespl=500      ! Computes energy each 'pener' iterations
integer    (kind=4)  :: pdenpar=50      ! Writes partial densities each 'pdenpar' iterations
integer    (kind=4)  :: pchem=50        ! Writes partial chemical potential each 'pchem' iter.
integer    (kind=4)  :: ppot=1          ! Computes the potential each 'ppot' iterations
integer    (kind=4)  :: tstgrid,Iterrk
real       (kind=8)  :: prec_Newton=1.d-8       ! Precisio pel Newton-Rapson
real       (kind=8)  :: norma           ! Use for normalization
real       (kind=8)  :: mimpur          ! Impurity mass in uma
real       (kind=8)  :: FT_analitic
real       (kind=8)  :: denxmax         ! Maximum value for density of the impurity
real       (kind=8)  :: r_clust         ! Radius of the helium cluster
real       (kind=8)  :: rimpur=5.5      ! Radius of the impurity
real       (kind=8)  :: p,px2,py2,pz2   ! Temporary variables for momentum values
real       (kind=8)  :: mu4=10.d0,mu4err      ! Value of Chemical potential and associated error
real       (kind=8)  :: epsx,epsxerr    ! Value of autovalue and associated error
real       (kind=8)  :: errmu4          ! Relative change betwen iteration for chemical potential
real       (kind=8)  :: Select_pot 
real       (kind=8)  :: deltat0,deltat  ! Step of time unvariable part and 'real step'
real       (kind=8)  :: deltat0x,deltatx ! Step of time unvariable part and 'real step'
real       (kind=8)  :: pafl=0.10       ! Default value for Paflov parameter
real       (kind=8)  :: precie=1.d-8    ! Cuando se alcanza esta precision en la energia, el programa se para
real       (kind=8)  :: eold            ! Auxiliar variables
complex    (kind=8)  :: ci=(0.d0,1.d0)  ! Imaginary unit
complex    (kind=8)  :: caux            ! Auxiliar variables
real       (kind=8)  :: aux,aux1        ! Auxiliar variables
integer    (kind=4),allocatable  :: nitera(:) ! In wich iteration change to the
real       (kind=8),allocatable  :: paflv(:)  !    corresponging Paflov coeffiecient
real       (kind=8)  :: cnorm           !    Normalization constant
real       (kind=8)  :: pmod1        !    Work variable for controlling p-values
real       (kind=8)  :: pmod2        !    Work variable for controlling p-values
real       (kind=8)  :: xcm4,ycm4,zcm4,xcmx,ycmx,zcmx ! Center of mass Drop and Impurity
real       (kind=8)  :: distx,disty,distz  ! Distance between center of masses
real       (kind=8)  :: rz,rx,ry,r,rmax,r2,qmax,Hq,FT_V_spline

integer    (kind=4)  :: n4=300          ! Number of helium_4 atoms
integer    (kind=4)  :: mode=0          ! Way to start the program (see readen subroutine)

character  (len=40)  :: title     = 'Helium4 - 3dim.  '
character  (len=60)  :: fileout   = 'DFT4He3d.res'
character  (len=60)  :: filedenin = 'he4.input.density'
character  (len=60)  :: filedenout= 'he4.output.density'
character  (len=60)  :: fileimpin = 'X.input.wf'
character  (len=60)  :: fileimpout= 'X.output.wf'
character  (len=60)  :: namefile,namefile1
!character  (len=1)   :: Vortex_axis='Z'  ! Are defined in the vortex module
!Logical              :: Lvortex=.false.
Logical              :: Lprint=.false.
Logical              :: LHe_frozen=.false., LFT_Test_numerica=.false.
real       (kind=8)  ::  sigma_x,sigma_y,sigma_z,zcom, V_Pi, V_Sigma, V_Delta

interface
  double precision function v_alka(d,elem)
       character (len=3),  intent(in) :: elem
       real      (kind=8), intent(in) :: d
  end function v_alka
end interface



!....................Variables read in a NAMELIST statement ..............................

namelist /vortex_z/xv,yv
namelist /vortex_y/xv,zv
namelist /vortex_x/yv,zv
namelist /input/title,fftwplan,nthread,nsfiles,                         &
                fileout,                                                &
                filedenin,filedenout,                                   &
                fileimpin,fileimpout,                                   &
                n4,mode,                                                &
                nx,ny,nz,xmax,ymax,zmax,                                &
                xc,yc,zc,afermi,                                        &
                eps4,sigma4,core4,l,                                    &
                cp4,cpp4,den4c,alphas,h2o2m4,                           &
                denmin,psimin,npd,ndmax,icon,icon_local,                    &
                niter,pafl,printpot,pchem,irespar,                          &
                pdenpar,pener,ppot,iron,ironx,lsolid,lden_max,solid_denmax, &
                limp,ximp,yimp,zimp,rimpur,mimpur,lrandom,r_shell,den_m,    &
                gwf,lpaflv,nstepp,leepot,nq,tol,rmin,vxpot,rinfi,lexternal, &
                Intens,zdist,precie,limp_despl,pdespl,pas_imp,pas_imp_max,  &
                N_newton, Prec_Newton, lmillorar_singularitats,lconstraint, &
                Lconstraint_imp,vdt,LHe_frozen,LFT_Test_numerica,           &
                dmu4,underpsil,Lbulk,mu4,rfermi,afermi,                     &
                selec_gs,r_cutoff_gs,umax_gs,Lvortex,Vortex_axis,nv,omega,  &
                selec_pi,r_cutoff_pi,umax_pi,Dtmin,Sdmu,Dmumax,             &
                selec_sigma,r_cutoff_sigma,umax_sigma,                      &
                selec_delta,r_cutoff_delta,umax_delta,Lrkpc,Lprint,         &
                Lexcite_state,instate,Als_P,Als_D,Xnew,defx,defy,defz,      &
                Lexcite_state_fix,Lparabolic_fit_FT_V_spline,               &
                Quita_C4_Ba_plus_gs_fix_C4,                                 &
                Quita_C4_Ba_plus_pi_fix_C4,                                 &
                Quita_C4_Ba_plus_sigma_fix_C4,Lstate,Ldiag_jz,Ljz,          &
                Laverage_P_value,Lexcite_state_external,                    &
                L_anell, r_anell, a_anell,                                  &
                L_esfera, r_esfera, a_esfera,                               &
                Exciplex, Lexciplex_state_fix,r_exc                         
!................................ Start main Program ..............................
call timer(t0)

!.............................................
!... Inicializate some numerical constants ...
!.............................................

pi     = 4.0d0*datan(1.0d0) ! Initialization of pi
twopi  = 2.0d0*pi
fourpi = 4.0d0*pi
piq    = pi*pi


!...............................
!... Read  master-input file ...
!...............................

read(5,nml=input,end=999)
open(10,file="DFT4He3d.namelist.read")
write(10,nml=input)
If(mode.Ne.0.And.mode.Ne.4.And.Lvortex)Then
  If(Vortex_Axis.Eq.'Z')Then
    Allocate(xv(nv)); Allocate(yv(nv))
    Open(unit=1,File='vortex.input')    
    Read(1,vortex_z)
    Write(6,vortex_z)
    Close(Unit=1)
  EndIf  
  If(Vortex_Axis.Eq.'Y')Then
    Allocate(xv(nv)); Allocate(zv(nv))
    Open(unit=1,File='vortex.input')    
    Read(1,vortex_y)
    Write(6,vortex_y)
    Close(Unit=1)
  EndIf  
  If(Vortex_Axis.Eq.'X')Then
    Allocate(yv(nv)); Allocate(zv(nv))
    Open(unit=1,File='vortex.input')    
    Read(1,vortex_x)
    Write(6,vortex_x)
    Close(Unit=1)
  EndIf  
EndIf        
call flush(10)

If(Lstate.Eq.'P')Then
  ninvar=6
Else
  ninvar=10
EndIf

If(lbulk)then
 Call rhoasin0(aux,densat4,mu4err,aux1)
 Write(6,'(" Potencial quimic a la saturacio......:",1p,E18.9)')mu4err
 If(mu4.Gt.0.0d0)mu4=mu4err
Endif

!
!  Aqui modifiquem els coeficients C4 dels potencials: 
!      Ba_plus_gs_fixC4, Ba_plus_pi_fixC4 i Ba_plus_sigma_fixC4
!
Pg_Ba_plus_gs_fixC4(1)    = Pg_Ba_plus_gs_fixC4(1)   *(1.0d0-Quita_C4_Ba_plus_gs_fix_C4   )
Pg_Ba_plus_pi_fixC4(1)    = Pg_Ba_plus_pi_fixC4(1)   *(1.0d0-Quita_C4_Ba_plus_pi_fix_C4   )
Pg_Ba_plus_sigma_fixC4(1) = Pg_Ba_plus_sigma_fixC4(1)*(1.0d0-Quita_C4_Ba_plus_sigma_fix_C4)

nn(1)  = nx ; nn(2)  = ny ; nn(3)  = nz;                ! Initialize things for PDERG
mmx(1) = nx ; mmx(2) = ny ; mmx(3) = nx ; mmx(4) = ny   ! (NO NOT MOVE THAT!!!!!!!)

!.............................................................
!.. Check if the size of the grid is among the valid values ..
!.............................................................

if(tstgrid(nx).ne.0 ) stop 'SEVERE ERROR: NX is not correct'
if(tstgrid(ny).ne.0 ) stop 'SEVERE ERROR: NY is not correct'
if(tstgrid(nz).ne.0 ) stop 'SEVERE ERROR: NZ is not correct'

!...................................................
!.. Controls Paflov parameters (read and storage ...
!...................................................

if(nstepp.lt.0) nstepp=1   ! Prevents possible errors

allocate(nitera(nstepp))   ! Allocate arrays with Paflov things...
allocate(paflv(nstepp))

!
    Call Init_deriv_p(npd,ndmax,nthread)
!

if(lpaflv) then   !........................ Change of Paflov allowed
  if(nstepp.eq.1) then    !(Only one factor...)
     nitera(1) = niter
     paflv(1)  = pafl
  else                    !(multiple factors... prepare in which iteration changes..)
     do ix=1,nstepp
       read(5,*)   nitera(ix),paflv(ix)
     end do
     naux = 0
     do ix=1,nstepp-1
       naux       = naux+nitera(ix)
       nitera(ix) = naux
     end do
     nitera(nstepp) = niter
  end if
else  !............................ If Change is not allowed. fix Paflov
  nstepp    = 1
  nitera(1) = niter
  paflv(1)  = pafl
end if

close(5)
close(10)

!................................................
!.. Some consistency check on input variables ...
!................................................

nthread=abs(nthread)

if(mode.eq.2) then
  if(.not.limp) then
     write(6,*) ' '
     write(6,*) ' Inconsistency input error.'
     write(6,*) ' '
     write(6,*) ' If mode=2. MUST be limp=.true.'
     write(6,*) ' '
     write(6,*) ' ACTION: ABORT EXECUTION'
     write(6,*) ' '
     STOP
  end if
end if

hx    = 2.0d0*abs(xmax)/(nx)  ! Step in x-grid
hy    = 2.0d0*abs(ymax)/(ny)  ! Step in y-grid
hz    = 2.0d0*abs(zmax)/(nz)  ! Step in z-grid

dxyz  = hx*hy*hz              ! Element of volum in real space
nxyz  = nx*ny*nz              ! Total number of points

hpx   = 1.0d0/(nx*hx)         ! Step in px-grid
hpy   = 1.0d0/(ny*hy)         ! Step in py-grid
hpz   = 1.0d0/(nz*hz)         ! Step in pz-grid

pmaxx = 1.0d0/(2.0d0*hx)      ! Maximum 'frequency' in X-grid
pmaxy = 1.0d0/(2.0d0*hy)      ! Maximum 'frequency' in Y-grid
pmaxz = 1.0d0/(2.0d0*hz)      ! Maximum 'frequency' in Z-grid


!...............................
!.. Dimensionate main ARRAYS ...
!...............................

call dimen()

!................................
!... Build grid in real space ...
!................................

do ix=1,nx  !.................... Grid X
 x(ix) = -xmax+hx*(ix-1)
end do
do iy=1,ny  !.................... Grid Y
 y(iy) = -ymax+hy*(iy-1)
end do
do iz=1,nz  !.................... Grid  Z
 z(iz) = -zmax+hz*(iz-1)
end do

!....................................
!... Build grid in momentum space ...
!....................................

!.... Build p-grid. In order to use FFTW the grid must
!     start from frequency zero to the maximum and then continue from
!     (-maximum) to zero (negative).

!............................................ grid Px
do ipx=1,nx/2+1
   px(ipx) =        hpx*(ipx-1)
end do
do ipx=nx/2+2,nx
   px(ipx) = -pmaxx+hpx*(ipx-(nx/2+1))
end do

!............................................ grid Py
do ipy=1,ny/2+1
   py(ipy) =        hpy*(ipy-1)
end do
do ipy=ny/2+2,ny
   py(ipy) = -pmaxy+hpy*(ipy-(ny/2+1))
end do

!............................................ grid Pz
do ipz=1,nz/2+1
   pz(ipz) =        hpz*(ipz-1)
end do
do ipz=nz/2+2,nz
   pz(ipz) = -pmaxz+hpz*(ipz-(nz/2+1))
end do

!............................................ Compule modulus of p
do ipz=1,nz
  pz2=pz(ipz)**2
  do ipy=1,ny
    py2=py(ipy)**2
    do ipx=1,nx/2+1
      px2               = px(ipx)**2
      pmod(ipx,ipy,ipz) = sqrt(px2+py2+pz2)
    end do
  end do
end do

pmod1=maxval(pmod)
pmod2=sqrt(pmaxx**2+pmaxy**2+pmaxz**2)

!................................
!... read density or build-it ...
!................................

n4real = dfloat(n4)



call readen(n4,densat4,filedenin,fileimpin,mode,rimpur,r_clust)

!===============================================================================
!===============================================================================
!===============================================================================
!
!  We imprint a vortex if we do not continue a calculation (if Lvortex=.true.)
!
If(mode.Ne.0.And.mode.Ne.4.And.Lvortex)Then
 If(Vortex_axis.Eq.'Z')Then
   Do iz = 1, nz
     Do iy = 1, ny
       Do ix = 1, nx
         Caux = (1.0d0, 0.d0)
         Do j=1, nv
           aux = DSqrt((x(ix)-xv(j))**2+(y(iy)-yv(j))**2)
           If(aux.eq.0.0d0)then
             Caux = (0.d0, 0.d0)
             Exit
           Else
             Caux = Caux * Cmplx((x(ix)-xv(j)),(y(iy)-yv(j)))/aux
           Endif
         EndDo        
         Psi(ix,iy,iz) = Sqrt(den(ix,iy,iz))*Caux       
       EndDo
     EndDo
   EndDo
 ElseIf(Vortex_axis.Eq.'Y')Then
   Do iz = 1, nz
     Do iy = 1, ny
       Do ix = 1, nx
         Caux = (1.0d0, 0.d0)
         Do j=1, nv
           aux = DSqrt((x(ix)-xv(j))**2+(z(iz)-zv(j))**2)
           If(aux.eq.0.0d0)then
             Caux = (0.d0, 0.d0)
             Exit
           Else
             Caux = Caux * Cmplx((x(ix)-xv(j)),(z(iz)-zv(j)))/aux
           Endif
         EndDo        
         Psi(ix,iy,iz) = Sqrt(den(ix,iy,iz))*Caux       
       EndDo
     EndDo
   EndDo
 ElseIf(Vortex_axis.Eq.'X')Then
   Do iz = 1, nz
     Do iy = 1, ny
       Do ix = 1, nx
         Caux = (1.0d0, 0.d0)
         Do j=1, nv
           aux = DSqrt((y(iy)-yv(j))**2+(z(iz)-zv(j))**2)
           If(aux.eq.0.0d0)then
             Caux = (0.d0, 0.d0)
             Exit
           Else
             Caux = Caux * Cmplx((y(iy)-yv(j)),(z(iz)-zv(j)))/aux
           Endif
         EndDo        
         Psi(ix,iy,iz) = Sqrt(den(ix,iy,iz))*Caux       
       EndDo
     EndDo
   EndDo
 Else
    Stop '009'
 Endif
 den=conjg(psi)*psi
Endif
!===============================================================================
!===============================================================================
!===============================================================================
!..............................................................................
     Open(Unit=1,file='den0-x.dat')
     Do ix=1,nx
       Write(1,'(1p,E15.6,5E19.11)')x(ix),den(ix,ny/2+1,nz/2+1)
     Enddo
     Close(Unit=1)

     Open(Unit=1,file='den0-y.dat')
     Do iy=1,ny
       Write(1,'(1p,E15.6,5E19.11)')y(iy),den(nx/2+1,iy,nz/2+1)
     Enddo
     Close(Unit=1)

     Open(Unit=1,file='den0-z.dat')
     Do iz=1,nz
       Write(1,'(1p,E15.6,5E19.11)')z(iz),den(nx/2+1,ny/2+1,iz)
     Enddo
     Close(Unit=1)

If(limp)Then

     aux = sum(denx)*dxyz
     Write(6,'(" Normalitzacio de la impureça...:",1p,E15.6)')aux

     Open(Unit=1,file='den0x-x.dat')
     Do ix=1,nx
       Write(1,'(1p,E15.6,5E19.11)')x(ix),denx(ix,ny/2+1,nz/2+1)
     Enddo
     Close(Unit=1)

     Open(Unit=1,file='den0x-y.dat')
     Do iy=1,ny
       Write(1,'(1p,E15.6,5E19.11)')y(iy),denx(nx/2+1,iy,nz/2+1)
     Enddo
     Close(Unit=1)

     Open(Unit=1,file='den0x-z.dat')
     Do iz=1,nz
       Write(1,'(1p,E15.6,5E19.11)')z(iz),denx(nx/2+1,ny/2+1,iz)
     Enddo
     Close(Unit=1)

Endif


  If(lbulk)Then
    n4real=sum(den)*dxyz
    n4=n4+0.5
    Write(6,'(" Nombre de partícules...:",1p,E15.6)')n4real
  Else 
    Write(6,'(" Nombre de partícules...:",1p,E15.6)')(sum(den)*dxyz)
  Endif


if(lexternal)then
   If(Lexcite_state.Or.Lexcite_state_external)Then
     Write(6,'(" Potencial exterior per:",A80,/,"r_cutoff_pi,umax_pi..:",1p,2E15.6)')selec_pi,r_cutoff_pi,umax_pi
     Write(6,'(" Potencial exterior per:",A80,/,"r_cutoff_sigma,umax_sigma..:",1p,2E15.6)')selec_sigma,r_cutoff_sigma,umax_sigma
     Write(6,'(" Potencial exterior per:",A80,/,"r_cutoff_delta,umax_delta..:",1p,2E15.6)')selec_delta,r_cutoff_delta,umax_delta
     Write(6,'("Lstate....:",A)')LState
     do iz=1,nz
       rz = (z(iz)-zimp)**2
       do iy=1,ny
         ry = (y(iy)-yimp)**2+rz
         do ix=1,nx
           r2 = (x(ix)-ximp)**2+ry
           r = dsqrt(r2)
           If (LState=='P') then
             aux=V_pi(r)
             Vpi(ix,iy,iz) = aux
             If(r2.Ne.0.0d0)Then
               Delta(ix,iy,iz) = (V_Sigma(r) - aux)/r2
             Else
               Delta(ix,iy,iz) = umax_Sigma
             Endif
           Else
             aux=V_Delta(r)
             Delta(ix,iy,iz)=aux
             Pi_Del(ix,iy,iz)=V_Pi(r)-aux
             Sig_Del(ix,iy,iz)=V_Sigma(r)-aux
           End if
         enddo
       enddo
     enddo
     If(Lexcite_state)Then
       If(Lexcite_state_external)Then
          Lprint_Invar=.true.
          Call Instates_external()
       Else
          Lprint_Invar=.true.
          Call Instates()
          Lprint_Invar=.false.
       EndIf
     EndIf
   Else
     Write(6,'(" Potencial exterior per:",A80,/,"r_cutoff_gs,umax_gs..:",1p,2E15.6)')selec_gs,r_cutoff_gs,umax_gs
     do iz=1,nz
       rz = (z(iz)-zimp)**2
       do iy=1,ny
         ry = (y(iy)-yimp)**2+rz
         do ix=1,nx
           r = dsqrt((x(ix)-ximp)**2+ry)
           uext(ix,iy,iz) = Select_Pot(selec_gs,r,r_cutoff_gs,umax_gs)
         enddo
       enddo
     enddo
   Endif

!
!  Construim la matriu factt, responsable de la variació local del temps imaginary 
!

  Do ix=1,nx
    Do iy=1,ny
      Do iz=1,nz
        aux=Uext(ix,iy,iz)
        factt(ix,iy,iz)=1.0d0
        If(aux.Ge.0.0d0)factt(ix,iy,iz)=1.d-1/(1.0d0+aux)
      Enddo
    Enddo
  Enddo

!  call generate_paflmap(uext,paflmap)
!  call respar(x,y,z,nx,ny,nz,2,'uext','paflmap',uext,paflmap)

!  call pderg(2,npd,nn,hx,uext,dx2uext,3,1,mmx,iw,icon)

!  call pderg(2,npd,nn,hx,uext,dx2uext,3,1,mmx,iw,icon)

!  call pderg(2,npd,nn,hx,uext,dx2uext,3,1,mmx,iw,icon)

!  call pderg(2,npd,nn,hx,uext,dx2uext,3,1,mmx,iw,icon)
!  call pderg(2,npd,nn,hy,uext,dy2uext,3,2,mmx,iw,icon)
!  call pderg(2,npd,nn,hz,uext,dz2uext,3,3,mmx,iw,icon)

!  Call derivnD(2,nn,hx,1,uext,dx2uext,Icon)
!  Call derivnD(2,nn,hy,2,uext,dy2uext,Icon)
!  Call derivnD(2,nn,hz,3,uext,dz2uext,Icon)

   h2o2mx = h2o2m4*(4.002602d0/mimpur)
Endif


if(irespar.ne.0) then
   if(limp) then
      call respar(x,y,z,nx,ny,nz,2,'hedenini','impurdenini',den,denx)
   else
      call respar(x,y,z,nx,ny,nz,1,'hedenini','hedenini',den,den)
   end if
end if

!....................................
!.. Print-out initial quantities  ...
!....................................

!open(6,file=fileout)

write(6,6010) title

select case(mode)
   case(0,4,5) !................................... Continue a calculation
      if(limp) then
         write(6,6111) filedenin,fileimpin,filedenout,fileimpout
      else
         write(6,6011) filedenin,filedenout
      end if
   case(1) !................................... Start from scratch a pure drop
      write(6,6012) filedenout
   case(2,3) !................................... Start from scratch one drop with impurity
      write(6,6013) filedenout,fileimpout
   case default !.............................. Start still not programed.
      write(6,*) ' '
      write(6,*) ' The variable mode has no acceptable value.'
      write(6,*) ' '
      stop
end select

!...............................................................
!................................ Write the input parameters ...
!...............................................................

write(6,6018) nthread,niter
if(mode.ne.0) then
  write(6,6020) n4,r_clust
else
  write(6,6025) n4
end if
write(6,6030) nx,ny,nz,hx, hy, hz, x(1) ,y(1) ,z(1) ,x(nx),y(ny),z(nz)
write(6,6035)          hpx,hpz,hpz,px(1),py(1),pz(1),pmaxx,pmaxy,pmaxz,&
                       pmod1,pmod2
write(6,6037) cp4,cpp4,den4c,alphas,l,den0s,h2o2m4


!...............
!.. Impurity ...
!...............

if(limp) then
   if(mimpur.lt.1.d-10) then
     print *,'********************'
     print *,'*** SEVERE ERROR ***'
     print *,'********************'
     print *,' '
     print *,'If LIMP=.TRUE. Mimpur cannor be zero...'
     print *,'  '
     print *,'PROGRAM CANCELLED'
   end if
   h2o2mx = h2o2m4*(4.002602d0/mimpur)
   write(6,6138) h2o2mx
   write(6,*) ' '
   write(6,*) '    Calculation with the impurity: ',elem
   write(6,fmt='(''         Atomic mass '',F10.3)') mimpur
   write(6,*) ' '
   write(6,*) '    Initial position of the impurity:'
   write(6,fmt='(''            X_imp = '',F10.5,'' A'')') ximp
   write(6,fmt='(''            Y_imp = '',F10.5,'' A'')') yimp
   write(6,fmt='(''            Z_imp = '',F10.5,'' A'')') zimp
   write(6,*) ' '
else
   write(6,*) '    Calculation without impurities.'
end if

!............................ Print everything  about Paflov parameters

if(lpaflv) then
  if(nstepp.eq.1) then
    write(6,6150) pafl
  else
    write(6,6038)
    write(6,6039) 1,nitera(1),nitera(1),paflv(1)
    do ix=2,nstepp
       write(6,6039) nitera(ix-1)+1,nitera(ix),nitera(ix)-nitera(ix-1)+1,paflv(ix)
    end do
  end if
else
  write(6,6150) pafl
end if

!...................................................................
!... Compute the FFT of Lennard-Jones                            ...
!... Prepara \alpha_s term in case of Orsay-Trento Interaction.  ...
!...................................................................

!
!  Per calcular amb la funcional sólida s'ha de treure el terme amb Alpha_s
!

If(Lsolid)core4='OT '

select case(core4)
   case('OP ')
     h4=h4op
     write(6,*) '    Use Orsay-Paris-Barcelona Interaction.'
     write(6,6040) core4,h4,eps4,sigma4,b4
   case('OT ')
     h4=h4ot
     write(6,*) '    Use Orsay-Trento Interaction. (ONLY CORE)'
     write(6,*) '    Do not calculate Alpha_s in the field neither energy.'
     If(Lsolid)write(6,'("    We will use the solid functional (see: PRB72, 214522(2005)) ")')
     write(6,6040) core4,h4,eps4,sigma4,b4
   case('OTE')
     h4=h4ot
     write(6,*) '    Use Orsay-Trento Interaction. '
     write(6,*) '    Calculate Alpha_s contribution ONLY in the energy..'
     write(6,6040) core4,h4,eps4,sigma4,b4
     allocate( denalf(nx  ,ny,nz))                                                            
     allocate(  falfs(nx  ,ny,nz))
     allocate(kalfs(nx/2+1,ny,nz))
     allocate(intxalf(nx  ,ny,nz))
     allocate(intyalf(nx  ,ny,nz))
     allocate(intzalf(nx  ,ny,nz))
   case('OTC')
     h4=h4ot
     write(6,*) '    Use Orsay-Trento Interaction.'
     write(6,*) '    Full Orsay-Trento calculation. (Field and Energy)'
     write(6,6040) core4,h4,eps4,sigma4,b4
     allocate( denalf(nx  ,ny,nz))                                                            
     allocate(  falfs(nx  ,ny,nz))
     allocate(kalfs(nx/2+1,ny,nz))
     allocate(intxalf(nx  ,ny,nz))
     allocate(intyalf(nx  ,ny,nz))
     allocate(intzalf(nx  ,ny,nz))
     allocate(ualphas(nx  ,ny,nz))
   case default
     print *,' ***************** WARNING ************************'
     print *,' '
     print *,' I do not know how to work with the core: ',core4
     print *,' '
     print *,' **************************************************'
     print *,' '
     STOP ' ... The program stops due to a severe error.'
end select

!...............................
!... Prepare plans for FFTWs ...
!...............................

write(6,*) '    Initialize Plans for FFT.'
call fftini(nx,ny,nz)

!...........................................................
!... Form  factor for Lennard-Jones and for the impurity ...
!...........................................................

write(6,*) '    Compute the FFT of the kernel of Lennard-Jones integrals.'

call fforma(core4,b4,eps4,sigma4,h4,nx,ny,nz,pmod,fvlj4)

!...........................................................
!... Form factor for Impurity potential and prepare      ...
!... arrays for the inclusion of the impurity            ...
!...                                                     ...
!... If (mode=2), start new calculation with impurity    ...
!... a hole un the helium density must be done....       ...
!...........................................................


if(limp) then
!   if(irespar.ne.0) then
!      open(31,file='PATIL.dat')
!      do iz=1,nz
!        aux = abs(z(iz))
!        write(31,3100) z(iz),v_alka(aux,elem)
!      end do
!      close(31)
!   end if
   qmax = Sqrt(pmaxx**2+pmaxy**2+pmaxz**2)
   Hq   = qmax/(Nq-1)
   If(Trim(selec_gs).Eq.'LJ_OT')Then
      Write(6,"('Usamos la interacción de OT para la impuerza')")
   Endif        
   Do iz=1, nz
     Do iy=1, ny
       Do ix=1, nx/2 + 1
         p=pmod(ix,iy,iz)
         If(Trim(selec_gs).Eq.'LJ_OT')Then
           Vq(ix,iy,iz) = fvlj4(ix,iy,iz)
         Else        
           Vq(ix,iy,iz) = cmplx(FT_V_spline(p,qmax,Hq,Selec_gs,r_cutoff_gs,umax_gs))
         Endif
       EndDo
     EndDo
   EndDo

   open(32,file='Interpolated.Impurity.Potential.new.dat')

do iz=nz/2+2,nz
   write(32,32000) pz(iz),real(vq(1,1,iz))
end do
do iz=1,nz/2+1
   write(32,32000) pz(iz),real(vq(1,1,iz))
end do

32000 format(5x,0P,F10.5,3x,1P,E13.5)
close(32)



!   selec=selec_gs
!   umax=umax_gs
!   call potimpx()
   forall(ix=1:nx,iy=1:ny,iz=1:nz)
      psi1x(ix,iy,iz)  = psix(ix,iy,iz)      ! Array initialized for imaginary step-time
      psi2x(ix,iy,iz)  = psix(ix,iy,iz)      ! Array initialized for imaginary step-time
   end forall
   call fftfw(den ,fden)
   call fftfw(denx,fdenx)
   forall(ix=1:nx/2+1,iy=1:ny,iz=1:nz)
     wk1(ix,iy,iz) = fden(ix,iy,iz)*vq(ix,iy,iz)
     wk2(ix,iy,iz) = fdenx(ix,iy,iz)*vq(ix,iy,iz)
   end forall
   call fftbk(wk1,upotx)                  ! Integral{denx*V_x}
   call fftbk(wk2,potx4)                  ! Integral{|Psi_x|**2 V_X}
!   if(irespar.ne.0) call respar(x,y,z,nx,ny,nz,2,'potfol','upotx',potx4,upotx)
!   do iz=1,nz
!     do iy=1,ny
!       do ix=1,nx
!         if(potx4(ix,iy,iz).ge.umax) then
!            den(ix,iy,iz) = denmin
!         end if
!       end do
!     end do
!   end do
   cnorm = n4/(sum(den)*dxyz)
   forall(ix=1:nx,iy=1:ny,iz=1:nz)
      den(ix,iy,iz) = cnorm*den(ix,iy,iz)
      psi(ix,iy,iz) = sqrt(den(ix,iy,iz))
   end forall
endif

!...............................................
!... Calculate of Fourier Transforms for 4He ...
!...............................................

call fftfw(den,fden)
!call fftfw(psi,fpsi)  ! Get FFT(Psi)

!........................................
!.. Initialize coarse-graining kernel ...
!........................................



If(Lsolid)Then

  aux=h4*1.065d0
  write(6,'("    Initialize Coarse-graining kernel, for Solid DF, h_cg=",1p, E15.6)')aux
  call initcg(aux,wcgk)

Else

  write(6,'("    Initialize Coarse-graining kernel, h_cg=",1p, E15.6)')h4
  call initcg(h4,wcgk)

Endif

!....................................
!.. First coarse-graining density ...
!.. First alfa_s density          ...
!....................................

forall(ix=1:nx/2+1,iy=1:ny,iz=1:nz)
   wk1(ix,iy,iz) = fden(ix,iy,iz)*wcgk(ix,iy,iz)
end forall
call fftbk(wk1,dencg)  ! get Coarse graining density

if(core4.eq.'OTE'.or.core4.eq.'OTC') then
   forall(ix=1:nx/2+1,iy=1:ny,iz=1:nz)
      kalfs(ix,iy,iz) = exp(-(pi*l*pmod(ix,iy,iz))**2)
      wk1(ix,iy,iz)   = fden(ix,iy,iz)*kalfs(ix,iy,iz)
   end forall
   call fftbk(wk1,denalf) ! Get Alfa_s density
end if

!...............................
!... Calculate Psi=sqrt(den) ...
!...............................

forall(ix=1:nx,iy=1:ny,iz=1:nz)
!   psi(ix,iy,iz)   = sqrt(den(ix,iy,iz))
   psi1(ix,iy,iz)  = psi(ix,iy,iz)           ! Array initialized for imaginary step-time
   psi2(ix,iy,iz)  = psi(ix,iy,iz)           ! Array initialized for imaginary step-time
end forall

!.................................
!.. First call to total energy ...
!.................................

call derden() ! Calculate derivatives of the density

if(core4.eq.'OTE'.or.core4.eq.'OTC')  then
  call term_alfa()
end if

call poten()              ! First Potential  (for Lagrange Equation)

If(Limp)Then

!  Write(6,'("Control(1): Normalitzacio, den ,denx..",1p,2E15.6)')(sum(den)*dxyz), (sum(denx)*dxyz)

  Open(Unit=11,File="upotx-x.dat")
  do ix=1,nx
    Write(11,'(1p,3E15.6)')x(ix),upotx(ix,(ny/2+1),(nz/2+1)),potx4(ix,(ny/2+1),(nz/2+1))
  EndDo
  Close(Unit=11)

  Open(Unit=11,File="upotx-y.dat")
  do iy=1,ny
    Write(11,'(1p,3E15.6)')y(iy),upotx((nx/2+1),iy,(nz/2+1)),potx4((nx/2+1),iy,(nz/2+1))
  EndDo
  Close(Unit=11)

  Open(Unit=11,File="upotx-z.dat")
  do iz=1,nz
    Write(11,'(1p,3E15.6)')z(iz),upotx((nx/2+1),(ny/2+1),iz),potx4((nx/2+1),(ny/2+1),iz)
  EndDo
  Close(Unit=11)

Endif




If(lexternal)Then
!
!   Impresio del potencial: uext & pot4
!
   Open(Unit=1,File='uext-xy.dat')
   Write(1,'("#")')
   Write(1,'("# Uext(x,y,0) & pot4(x,y,0) at begining of the process")')
   Write(1,'("#")')
   Do ix=1,nx
     Do iy=1,ny
       Write(1,'(1p,6E18.10)')x(ix),y(iy),uext(ix,iy,(nz/2+1)),pot4(ix,iy,(nz/2+1))
     EndDo
     Write(1,*)
   EndDo
   Close(Unit=1)

   Open(Unit=1,File='uext-yz.dat')
   Write(1,'("#")')
   Write(1,'("# Uext(0,y,z) & pot4(0,y,z) at begining of the process")')
   Write(1,'("#")')
   Do iy=1,ny
     Do iz=1,nz
       Write(1,'(1p,6E18.10)')y(iy),z(iz),uext((nx/2+1),iy,iz),pot4((nx/2+1),iy,iz)
     EndDo
     Write(1,*)
   EndDo
   Close(Unit=1)

   Open(Unit=1,File='uext-xz.dat')
   Write(1,'("#")')
   Write(1,'("# Uext(x,0,z) & pot4(x,0,z) at begining of the process")')
   Write(1,'("#")')
   Do ix=1,nx
     Do iz=1,nz
       Write(1,'(1p,6E18.10)')x(ix),z(iz),uext(ix,(ny/2+1),iz),pot4(ix,(ny/2+1),iz)
     EndDo
     Write(1,*)
   EndDo
   Close(Unit=1)

   Open(Unit=1,File='uext-x.dat')
   Write(1,'("#")')
   Write(1,'("# Uext(x,0,0) & pot4(x,0,0) at begining of the process")')
   Write(1,'("#")')
   Do ix=1,nx
       Write(1,'(1p,6E18.10)')x(ix),uext(ix,(ny/2+1),(nz/2+1)),pot4(ix,(ny/2+1),(nz/2+1))
   EndDo
   Close(Unit=1)

   Open(Unit=1,File='uext-y.dat')
   Write(1,'("#")')
   Write(1,'("# Uext(0,y,0) & pot4(0,y,0) at begining of the process")')
   Write(1,'("#")')
   Do iy=1,ny
     Write(1,'(1p,6E18.10)')y(iy),uext((nx/2+1),iy,(nz/2+1)),pot4((nx/2+1),iy,(nz/2+1))
   EndDo
   Close(Unit=1)

   Open(Unit=1,File='uext-z.dat')
   Write(1,'("#")')
   Write(1,'("# Uext(0,0,z) & pot4(0,0,z) at begining of the process")')
   Write(1,'("#")')
   Do iz=1,nz
     Write(1,'(1p,6E18.10)')z(iz),uext((nx/2+1),(ny/2+1),iz),pot4((nx/2+1),(ny/2+1),iz)
   EndDo
   Close(Unit=1)
Endif

if(limp) then
!  write(6,6060) eimpu,ekinx,etot
!  irimp = maxloc(denx)
!  ximp = x(irimp(1))
!  yimp = y(irimp(2))
!  zimp = z(irimp(3))
  call r_cm(denx,1,xcmx,ycmx,zcmx)    ! Center of mass of Impurity
  ximp = xcmx; yimp = ycmx; zimp = zcmx
!  write(6,6065) ximp,yimp,zimp
endif

call energy()             ! Calculate energies

If(lconstraint)Then

! ..............................................
zcom = 0.d0
do iz=1,nz
 zcom = zcom + sum(den(:,:,iz))*z(iz)
enddo
zcom = zcom*dxyz/n4real
print*,'---- Center of Mass z: ',zcom
print*,'---- CONSTRAINT ENERGY: ',enercons
! ..............................................

Endif

write(6,6050) etot4,etot/n4,ekin4,elj4,ealphas,esolid,ecor4

If(lexternal)Then
  Write(6,'(T5,"Spin-Orbit energy.............:",8x,1p,E15.6," K")')eso
  Write(6,'(T5,"Impurity energy...............:",8x,1p,E15.6," K")')eimpu
  Write(6,'(T5,"Total energy (He+X)...........:",8x,1p,E15.6," K")')etot
end if

!
!     Here we will compute <L_axis>
!
      If(Lvortex)Then
        caux = (0.d0, 0.d0)      
        If(Vortex_axis.Eq.'Z')Then
          Call derivnD(1,nn,hx,1,psi,sto1c,Icon)
          Call derivnD(1,nn,hy,2,psi,sto2c,Icon)
          Do iz=1, nz
            Do iy=1, ny
              Do ix=1, nx
                caux = caux + Ci*Conjg(Psi(ix,iy,iz))*                  &
                (y(iy)*sto1c(ix,iy,iz) - x(ix)*sto2c(ix,iy,iz)) 
              EndDo
            EndDo
          EndDo
          caux = caux*dxyz
          Write(6,'(T5,"<L_z>.........................:",8x,1p,2E15.6," \hbar")')caux 
          IF(nv.gt.1) then
            Write(6,'(T5,"-Omega*<L_z>..................:",8x,1p,2E15.6," K")')-omega*caux 
            Write(6,'(T5,"<H>-Omega*<L_z>...............:",8x,1p,2E15.6," K")')(etot4-omega*caux)
          Endif 
        Endif
        If(Vortex_axis.Eq.'Y')Then
          Call derivnD(1,nn,hx,1,psi,sto1c,Icon)
          Call derivnD(1,nn,hz,3,psi,sto2c,Icon)
          Do iz=1, nz
            Do iy=1, ny
              Do ix=1, nx
                caux = caux + Ci*Conjg(Psi(ix,iy,iz))*                  &
                (x(ix)*sto2c(ix,iy,iz) - z(iz)*sto1c(ix,iy,iz)) 
              EndDo
            EndDo
          EndDo
          caux = caux*dxyz
          Write(6,'(T5,"<L_y>.........................:",8x,1p,2E15.6," \hbar")')caux 
          IF(nv.gt.1) then
            Write(6,'(T5,"-Omega*<L_y>..................:",8x,1p,2E15.6," K")')-omega*caux 
            Write(6,'(T5,"<H>-Omega*<L_y>...............:",8x,1p,2E15.6," K")')(etot4-omega*caux)
          Endif 
        Endif
        If(Vortex_axis.Eq.'X')Then
          Call derivnD(1,nn,hy,2,psi,sto1c,Icon)
          Call derivnD(1,nn,hz,3,psi,sto2c,Icon)
          Do iz=1, nz
            Do iy=1, ny
              Do ix=1, nx
                caux = caux + Ci*Conjg(Psi(ix,iy,iz))*                  &
                (z(iz)*sto1c(ix,iy,iz) - y(iy)*sto2c(ix,iy,iz)) 
              EndDo
            EndDo
          EndDo
          caux = caux*dxyz
          Write(6,'(T5,"<L_x>.........................:",8x,1p,2E15.6," \hbar")')caux 
          IF(nv.gt.1) then
            Write(6,'(T5,"-Omega*<L_x>..................:",8x,1p,2E15.6," K")')-omega*caux 
            Write(6,'(T5,"<H>-Omega*<L_x>...............:",8x,1p,2E15.6," K")')(etot4-omega*caux)
          Endif 
        Endif
      Endif        
if(limp) then
  write(6,6060) eimpu,ekinx,etot
  write(6,6065) ximp,yimp,zimp
endif

eold = etot
call flush(6)

!-------------------------------------------------------------------------------
!---                            Iterative procedure                           --
!-------------------------------------------------------------------------------

call evolo(0.0d0, mu4, mu4err,n4)               ! First Chemical potential
if(limp) call evolox(0.0d0,epsx,epsxerr)        ! First autovalue for impurity

deltat0  = (min(hx,hy,hz))**2/(4.0d0*h2o2m4)    ! 'piece for the Step-time
deltat0x = (min(hx,hy,hz))**2/(4.0d0*h2o2mx)    ! 'piece for the Step-time

if(limp) then
   deltat0  = (min(hx,hy,hz))**2/(4.0d0*h2o2m4) ! 'piece for the Step-time
   deltat0x = (min(hx,hy,hz))**2/(4.0d0*h2o2mx) ! 'piece for the Step-time
   deltat0  = min(deltat0,deltat0x)
   deltatx  = deltat0                           !<--   OJO!!!!!!!!  PAY ATTENTION!!!!!!
   deltatx  = pafl*deltatx                      ! Step-Time (For impurity)
   deltat   = pafl*deltat0                      ! Step-Time (For Helium)
else
   deltat0  = (min(hx,hy,hz))**2/(4.0d0*h2o2m4) ! 'piece for the Step-time
   deltat   = pafl*deltat0                      ! Step-Time (For Helium)
end if

lfilepv  = .false.

! write(6,7000)
call timer(t5)


      if(limp) then
        write(6,7001)
      else
        write(6,7000)
      end if

Iterrk=0

do iter=1,niter       ! <--------------------------------- Iterative procedure starts here.

!
!.......................................... if the change of the Paflov parameter
!.......................................... calculate the deltat adequate.

   if(lpaflv) then                        
     do ix=1,nstepp
       naux   = iter-nitera(ix)
       deltat = paflv(ix)*deltat0
       if(limp) deltatx= paflv(ix)*deltat0x
       if(naux.lt.0) then
          lfilepv=.false.
          exit
       elseif(naux.eq.0) then
          lfilepv=.true.
          Iterrk=0
          exit
       else
          continue
       end if
     end do
   end if

   Iterrk=Iterrk+1

!...........................................
 If(.Not.LHe_frozen)Then
   If(.Not.Lrkpc)Then
       call evolo(deltat,mu4,mu4err,n4)              ! Get new (density, chemical potential...)
   Else
     If(Iterrk.le.4)Then
       call steprkr(deltat,mu4,mu4err,n4)
     Else
       call steppcr(deltat,mu4,mu4err,n4)
       If(lprint)write(6,'("Error(He) (From Steppc)...",1p,e15.6)')ErrPC
     Endif
   Endif
   call derden()                                 ! Derivatives of the density
!
!  Desplaçem la impureça per tal que vagi cap el minim de potencial
!
   call fftfw(den,fden)   ! FFT of den
 EndIf   
   if(limp) then 
       call evolox(deltatx,epsx,epsxerr)         ! Get new (density, chemical potential...)
       call fftfw(denx,fdenx)                    ! FFT of den
   end if
 If(.Not.LHe_frozen)Then
   if(core4.eq.'OTC'.or.core4.eq.'OTE') then     
      forall(ix=1:nx/2+1,iy=1:ny,iz=1:nz)
         wk1(ix,iy,iz) = fden(ix,iy,iz)*wcgk(ix,iy,iz)
         wk2(ix,iy,iz) = fden(ix,iy,iz)*kalfs(ix,iy,iz)
      end forall
      call fftbk(wk1,dencg)                      ! Get Den of Coarse graining
      call fftbk(wk2,denalf)                     ! Get Alfa_s density
   else
      forall(ix=1:nx/2+1,iy=1:ny,iz=1:nz)
         wk1(ix,iy,iz) = fden(ix,iy,iz)*wcgk(ix,iy,iz)
      end forall
      call fftbk(wk1,dencg)                      ! Get Den of Coarse graining
   end if
 EndIf
!..............................................................................

   if(mod(iter,ppot).eq.0) then           ! Compute New Potential
      call poten()                       
   end if
!
!  Desplacem la impureça si ho creiem convenient
!
   If(Mod(iter,pdespl).Eq.0.And.limp_despl.And..Not.limp)Then
     
     Do it=1,N_Newton
       Call derivnD(1,nn,hz,3,den,dzden,Icon)
       Call derivnD(2,nn,hz,3,den,sto1,Icon)
       aux=-pas_imp*Sum(dzden*uext)/Sum(sto1*uext)
       If(N_Newton.Eq.it.Or.Dabs(aux).lt.prec_Newton)Then
         If(Dabs(aux).Gt.Pas_imp_max)Then
           aux=aux/Dabs(aux)*pas_imp_max
         Endif
         Write(6,'("Desplaçament de la impureça i nova posició:",1p,2E18.10)')aux,zimp
       Endif
       Write(6,'("it,aux",i5,1p,E15.6)')it,aux
       zimp=zimp+pas_imp_max*aux
       If(dabs(aux).Lt.Prec_Newton)Exit
       if(lexternal)then
!         If(elem.Eq.'AG ')Then
!           Write(6,'(" Potencial exterior de la plata")')
!           call generatepoten(uext)
!         ElseIf(elem.Eq.'BA+')Then
!!           Write(6,'(" Potencial exterior del Ba+")')
!!      call generatepoten_Ba_plus(uext)
!           call generatepoten_Ba_plus_new(uext)
!         Else

!           Write(6,'(" Potencial exterior per:",A4," (Patil)")')elem
           do iz=1,nz
             rz = (z(iz)-zimp)**2
             do iy=1,ny
               ry = (y(iy)-yimp)**2+rz
               do ix=1,nx
                 r = dsqrt((x(ix)-ximp)**2+ry)
                 uext(ix,iy,iz) = Select_Pot(selec_gs,r,r_cutoff_gs,umax_gs)
!                 if(r.gt.rmax)then
!                   uext(ix,iy,iz) = min(v_alka(r,elem),umax)
!                 Else
!                   uext(ix,iy,iz) = umax
!                 Endif
               enddo
             enddo
           enddo
!         Endif
       Endif
     EndDo
     If(N_Newton.Gt.1.and.Dabs(aux).Gt.Prec_newton)Write(6,'(" El process de Newton no ha convergit",1p,E15.6)')aux
   Endif

!..............................................................................

   if(mod(iter,pener).eq.0) then          ! Compute New energy and max of density
!      call fftfw(psi,fpsi)  ! Get FFT(Psi)


     call energy()

If(lconstraint)Then

! ..............................................
zcom = 0.d0
do iz=1,nz
 zcom = zcom + sum(den(:,:,iz))*z(iz)
enddo
zcom = zcom*dxyz/n4real
print*,'---- Center of Mass z: ',zcom
print*,'---- CONSTRAINT ENERGY: ',enercons
! ..............................................

Endif
 
      write(6,7010) etot4,(etot-eold),etot/n4,ekin4,elj4,ealphas,esolid,ecor4
      call varmu(n4,mu4,errmu4)           ! Error in mu4
      write(6,7017) mu4,errmu4

      If(Lbulk)Then
        n4real=sum(den)*dxyz
        Write(6,'("    Number of particles........ ",1p,E18.9)')n4real
        n4=n4real+0.5
      Endif

      if(limp) then
         call vareps(epsx,epsxerr)           ! Error in mu4
         write(6,7018) epsx,epsxerr
         write(6,7015) eimpu,ekinx,etot

!         irimp = maxloc(denx)
!         ximp = x(irimp(1))
!         yimp = y(irimp(2))
!         zimp = z(irimp(3))
         call r_cm(denx,1,xcmx,ycmx,zcmx)    ! Center of mass of Impurity
         ximp = xcmx; yimp = ycmx; zimp = zcmx
         write(6,7016) ximp,yimp,zimp
      elseif(lexternal)then
         Write(6,'(T5,"Spin-Orbit energy..........",8x,1p,E15.6," K")')eso
         Write(6,'(T5,"Impurity energy............",8x,1p,E15.6," K")')eimpu
         Write(6,'(T5,"Total energy (He+X)........",8x,1p,E15.6," K")')etot
      end if
!
!     Here we will compute <L_axis>
!
      If(Lvortex)Then
        caux = (0.d0, 0.d0)      
        If(Vortex_axis.Eq.'Z')Then
          Call derivnD(1,nn,hx,1,psi,sto1c,Icon)
          Call derivnD(1,nn,hy,2,psi,sto2c,Icon)
          Do iz=1, nz
            Do iy=1, ny
              Do ix=1, nx
                caux = caux + Ci*Conjg(Psi(ix,iy,iz))*                  &
                (y(iy)*sto1c(ix,iy,iz) - x(ix)*sto2c(ix,iy,iz)) 
              EndDo
            EndDo
          EndDo
          caux = caux*dxyz
          Write(6,'(T5,"<L_z>.....................",8x,1p,2E15.6," \hbar")')caux 
          IF(nv.gt.1) then
            Write(6,'(T5,"-Omega*<L_z>..............",8x,1p,2E15.6," K")')-omega*caux 
            Write(6,'(T5,"<H>-Omega*<L_z>...........",8x,1p,2E15.6," K")')(etot4-omega*caux)
          Endif 
        Endif
        If(Vortex_axis.Eq.'Y')Then
          Call derivnD(1,nn,hx,1,psi,sto1c,Icon)
          Call derivnD(1,nn,hz,3,psi,sto2c,Icon)
          Do iz=1, nz
            Do iy=1, ny
              Do ix=1, nx
                caux = caux + Ci*Conjg(Psi(ix,iy,iz))*                  &
                (x(ix)*sto2c(ix,iy,iz) - z(iz)*sto1c(ix,iy,iz)) 
              EndDo
            EndDo
          EndDo
          caux = caux*dxyz
          Write(6,'(T5,"<L_z>.....................",8x,1p,2E15.6," \hbar")')caux 
          IF(nv.gt.1) then
            Write(6,'(T5,"-Omega*<L_z>..............",8x,1p,2E15.6," K")')-omega*caux 
            Write(6,'(T5,"<H>-Omega*<L_z>...........",8x,1p,2E15.6," K")')(etot4-omega*caux)
          Endif 
        Endif
        If(Vortex_axis.Eq.'X')Then
          Call derivnD(1,nn,hy,2,psi,sto1c,Icon)
          Call derivnD(1,nn,hz,3,psi,sto2c,Icon)
          Do iz=1, nz
            Do iy=1, ny
              Do ix=1, nx
                caux = caux + Ci*Conjg(Psi(ix,iy,iz))*                  &
                (z(iz)*sto1c(ix,iy,iz) - y(iy)*sto2c(ix,iy,iz)) 
              EndDo
            EndDo
          EndDo
          caux = caux*dxyz
          Write(6,'(T5,"<L_x>.....................",8x,1p,2E15.6," \hbar")')caux 
          IF(nv.gt.1) then
            Write(6,'(T5,"-Omega*<L_x>..............",8x,1p,2E15.6," K")')-omega*caux 
            Write(6,'(T5,"<H>-Omega*<L_x>...........",8x,1p,2E15.6," K")')(etot4-omega*caux)
          Endif 
        Endif
      Endif        
      If(Abs(eold-etot).Le.precie)Exit
      eold = etot

      if(limp) then
        call r_cm(den,n4,xcm4,ycm4,zcm4)    ! Center of mass of 4He Drop
        call r_cm(denx,1,xcmx,ycmx,zcmx)    ! Center of mass of Impurity
        distx=abs(xcmx-xcm4)
        disty=abs(ycmx-ycm4)
        distz=abs(zcmx-zcm4)
        write(6,7110) xcm4,ycm4,zcm4,xcmx,ycmx,zcmx,distx,disty,distz
      else
        call r_cm(den,n4,xcm4,ycm4,zcm4)    ! Center of mass of 4He Drop
        write(6,7100) xcm4,ycm4,zcm4
      end if
      if(limp) then
        write(6,7001)
      else
        write(6,7000)
      end if
   end if

!..............................................................................

   if(mod(iter,pdenpar).eq.0) then        ! Save partial densities
     nsfaux = nsfaux+1
     if(nsfaux.gt.nsfiles) nsfaux=mod(nsfaux,nsfiles)
     select case(nsfaux)
         case(1:9)
           write(namefile, 8010) nsfaux
           write(namefile1,8015) nsfaux
         case(10:99)
           write(namefile, 8020) nsfaux
           write(namefile1,8025) nsfaux
         case(100:999)
           write(namefile, 8030) nsfaux
           write(namefile1,8035) nsfaux
     end select

     if(limp) then
       call respar(x,y,z,nx,ny,nz,2,'den','denx',den,denx)
       if(irespar.ne.0) call respar(x,y,z,nx,ny,nz,2,'folding','foldingx',potx4,upotx)
       call printout(3,namefile,namefile1,psi,elem,nx,ny,nz,hx,hy,hz,limp, &
                    xmax,ymax,zmax,ximp,yimp,zimp,psix,         &
                    deltat/deltat0,iter)
     else
       call respar(x,y,z,nx,ny,nz,1,'den','den',den,den)
       call printout(3,namefile,namefile1,psi,elem,nx,ny,nz,hx,hy,hz,limp, &
                    xmax,ymax,zmax,ximp,yimp,zimp,den,         &
                    deltat/deltat0,iter)
     end if

     ! ................................................................... !
     ! Let's compute the width of the impurity wavefunction to first order
     ! ................................................................... !
     if(lexternal)then
       sigma_x = ( 0.5d0*h2o2mx / (sum(den*dx2uext)*dxyz) )**0.25d0      
       sigma_y = ( 0.5d0*h2o2mx / (sum(den*dy2uext)*dxyz) )**0.25d0      
       sigma_z = ( 0.5d0*h2o2mx / (sum(den*dz2uext)*dxyz) )**0.25d0      
       write(*,*)'::::... Width of the impurity wavefunction ...:::'
       write(*,*)'Sigma(x/y/z)= ',sigma_x,sigma_y,sigma_z
       write(*,*)'::::..........................................:::'
     endif

!..............................................................................
     Open(Unit=1,file='den-x.dat')
     Do ix=1,nx
          Write(1,'(1p,E15.6,5E19.11)')x(ix),den(ix,ny/2+1,nz/2+1),(hpsi(ix,ny/2+1,nz/2+1)/psi(ix,ny/2+1,nz/2+1)-mu4)
     Enddo
     Close(Unit=1)

     Open(Unit=1,file='den-y.dat')
     Do iy=1,ny
          Write(1,'(1p,E15.6,5E19.11)')y(iy),den(nx/2+1,iy,nz/2+1),(hpsi(nx/2+1,iy,nz/2+1)/psi(nx/2+1,iy,nz/2+1)-mu4)
     Enddo
     Close(Unit=1)

     Open(Unit=1,file='den-z.dat')
     Do iz=1,nz
          Write(1,'(1p,E15.6,5E19.11)')z(iz),den(nx/2+1,ny/2+1,iz),(hpsi(nx/2+1,ny/2+1,iz)/psi(nx/2+1,ny/2+1,iz)-mu4)
     Enddo
     Close(Unit=1)

   end if
!..............................................................................

   if(lfilepv) then
     select case(iter)
         case(1:9)
           write(namefile ,5010) iter
           write(namefile1,5015) iter
         case(10:99)
           write(namefile ,5020) iter
           write(namefile1,5025) iter
         case(100:999)
           write(namefile ,5030) iter
           write(namefile1,5035) iter
         case(1000:9999)
           write(namefile ,5040) iter
           write(namefile1,5045) iter
         case(10000:99999)
           write(namefile ,5050) iter
           write(namefile1,5055) iter
         case(100000:999999)
           write(namefile ,5060) iter
           write(namefile1,5065) iter
     end select


     if(limp) then
        call printout(2,namefile,namefile1,psi,elem,nx,ny,nz,hx,hy,hz,limp, &
                      xmax,ymax,zmax,ximp,yimp,zimp,psix,                   &
                      deltat/deltat0,iter)
     else
        call printout(2,namefile,namefile1,psi,elem,nx,ny,nz,hx,hy,hz,limp, &
                      xmax,ymax,zmax,ximp,yimp,zimp,den,                   &
                      deltat/deltat0,iter)
     end if


   end if
!..............................................................................

   call timer(t6)                         ! Compute use time

   if(mod(iter,pchem).eq.0) then          ! Save partial densities
     if(limp) then
!       write(6,7001)
       write(6,7035) iter,mu4,mu4err,epsx,epsxerr,t6,(t6-t5)
     else
!       write(6,7000)
       write(6,7030) iter,mu4,mu4err,t6,(t6-t5)
     end if
     If((lrkpc).And..Not.Lprint)write(6,'("Error(He) (From Steppc)...",1p,e15.6)')ErrPC
   end if
   t5=t6
end do
if(limp) then
   call printout(1,filedenout,fileimpout,psi,elem,nx,ny,nz,hx,hy,hz,limp, &
                 xmax,ymax,zmax,ximp,yimp,zimp,psix,         &
                 deltat/deltat0,iter)
else
   call printout(1,filedenout,fileimpout,psi,elem,nx,ny,nz,hx,hy,hz,limp, &
                 xmax,ymax,zmax,ximp,yimp,zimp,den,         &
                 deltat/deltat0,iter)
end if

!..............................................................................
     Open(Unit=1,file='den-x.dat')
     Do ix=1,nx
          Write(1,'(1p,E15.6,5E19.11)')x(ix),den(ix,ny/2+1,nz/2+1),(hpsi(ix,ny/2+1,nz/2+1)/psi(ix,ny/2+1,nz/2+1)-mu4)
     Enddo
     Close(Unit=1)

     Open(Unit=1,file='den-y.dat')
     Do iy=1,ny
          Write(1,'(1p,E15.6,5E19.11)')y(iy),den(nx/2+1,iy,nz/2+1),(hpsi(nx/2+1,iy,nz/2+1)/psi(nx/2+1,iy,nz/2+1)-mu4)
     Enddo
     Close(Unit=1)

     Open(Unit=1,file='den-z.dat')
     Do iz=1,nz
          Write(1,'(1p,E15.6,5E19.11)')z(iz),den(nx/2+1,ny/2+1,iz),(hpsi(nx/2+1,ny/2+1,iz)/psi(nx/2+1,ny/2+1,iz)-mu4)
     Enddo
     Close(Unit=1)
!..............................................................................

call timer(t4)
print *,' Total  ',t4-t0

stop
999 stop 'DFT3He3d. Error in input master file. Too short'

!...............
!... Formats ...
!...............

3100 format(3x,0P,f9.4,2x,1P,E13.5)

6010 format(//,&
T10,'   ######  ####### ####### #       #     #          #####          ',/,  &
T10,'   #     # #          #    #    #  #     #  ###### #     #  #####  ',/,  &
T10,'   #     # #          #    #    #  #     #  #            #  #    # ',/,  &
T10,'   #     # #####      #    #    #  #######  #####   #####   #    # ',/,  &
T10,'   #     # #          #    ####### #     #  #            #  #    # ',/,  &
T10,'   #     # #          #         #  #     #  #      #     #  #    # ',/,  &
T10,'   ######  #          #         #  #     #  ######  #####   #####  ',//, &
T6,'Title of the run: ',A)

6011 format(//,T6,'CONTINUE a calculation:',//,&
               T6,'Input  densitity file: ',A,/,&
               T6,'Output densitity file: ',A)

6111 format(//,T6,'CONTINUE a calculation:',//,&
               T6,'Input file with helium densitity    : ',A,/,&
               T6,'Input file with impurity wave func. : ',A,/,&
               T6,'Output file with helium densitity   : ',A,/,&
               T6,'Output file with impurity wave func.: ',A,/,' ')

6012 format(//,T6,'Start a new calculation:',//,&
               T6,'Output densitity file: ',A)
6013 format(//,T6,'Start a new calculation with an impurity:',//,&
               T6,'Output file for Helium density: ',A,/,        &
               T6,'Output file for the impurity wave function: ',A)
6018 format(//,T6,'Number of threads:    ',I6,/,&
               T6,'Number of iterations: ',i16)
6020 format(//,T6,'Number of particles:    ',0P,I12,/,&
               T6,'Radius of the cluster : ',F10.3,' A')
6025 format(//,T6,'Number of particles:    ',0P,I12,/,' ')
6030 format(//,T6,'+-------------------+----------------------------------------+',/,&
               T6,'| REAL GRID         |     X-grid       Y-grid       Z-grid   |',/,&
               T6,'+-------------------+----------------------------------------+',/,&
               T6,'| Number of points  |',0P,T32,I4,T45,I4,T58,I4,T66,' |',/,&
               T6,'| Step              |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'| Min value         |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'| Max value         |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'+-------------------+----------------------------------------+')
6035 format(//,T6,'+-------------------+----------------------------------------+',/,&
               T6,'| MOMEMTUM GRID     |    Px-grid      Py-grid      Pz-grid   |',/,&
               T6,'+-------------------+----------------------------------------+',/,&
               T6,'| Step              |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'| Min value         |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'| Max value         |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'+-------------------+----------------------------------------+' ,//,&
               T6,'Maximum modulus of p ',0P,F11.6,3x,F11.6,/,' ')
6037 format(//,T6,'Parameters for the funcional:',//,          &
               T8,'cp4 ...... ',1P,E13.5 ,' K \AA**6'      ,/, &
               T8,'cpp4 ..... ',1P,E13.5 ,' K \AA**9'      ,/, &
               T8,'den4c .... ',0P,F13.3 ,' \AA**{-3}'     ,/, &
               T8,'Alphas ... ',0P,F13.3 ,' K ^-1 \AA**3'  ,/, &
               T8,'L ........ ',0P,F13.2 ,' \AA'           ,/, &
               T8,'den0s .... ',0P,F13.2 ,' \AA**-3'       ,/, &
               T8,'h2o2m4 ... ',0P,F14.11,' hbar**2 / (2 m_4)' )
6138 format(' ',T8,'h2o2mx ... ',0P,F14.11,' hbar**2 / (2 m_x)' )
6038 format(//,T6,'Change of Paflov parameter allowed: ',//, &
     T18,'From     to      iter     Factor',/, &
     T19,'------  ------  ------  -----------')

6039 format(1x,0P,T17,i6,T25,I6,T33,i6,t42,f11.7)

6150 format(//,T6,'Pavlov parameter fixed for all the run to: ',F8.4)

6040 format( /,T6,'Lennard-Jones parameters:',//,&
               T10,'Core    ',A3,/,&
               T10,'h     ',F11.7,' A',/,&
               T10,'eps   ',F11.7,' K',/,&
               T10,'sigma ',F11.7,' A',/,&
               T10,'b     ',F11.3,' K A**3 '//,' ')
6050 format(//,T5,'FIRST ENERGY BALANCE: ',                    //     &
              ,T5,'TOTAL   energy (He) ..........: ',F18.6,' K',/,    &
               T5,'Energy per particle (He+X) ...: ',F18.6,' K',/,    &
               T5,'Kinetic energy (He) ..........: ',F18.6,' K',/,    &
               T5,'Lennard-Jones energy (He) ....: ',F18.6,' K',/,    &
               T5,'Alpha_s term  energy (He) ....: ',F18.6,' K',/,    &
               T5,'Solid energy (He)  ...........: ',F18.6,' K',/,    &
               T5,'Correlation energy   (He) ....: ',F18.6,' K')
6060 format(1x,T5,'Impurity energy (X-He) .......: ',F18.6,' K',/,    &
               T5,'Kinetic energy (X) ...........: ',F18.6,' K',/,    &
               T5,'TOTAL ENERGY (He+X) ..........: ',F18.6,' K',/)
6065 format(1x,T5,'Impurity location  (x-axix) ..: ',F18.6,' K',/,    &
               T5,'                   (y-axis) ..: ',F18.6,' K',/,    &
               T5,'                   (z-axis) ..: ',F18.6,' K',/)
7000 format(//,1x,T2,'Iter             Mu(K)      Err(Mu)    Ttime  / Lap Time',/,&
 '----------------------------------------------------------')
7001 format(//,1x,T2, &
     'Iter             Mu(K)      Err(Mu)    Autovalue(K)   err(K)   ETtime  / Lap Time',&
     /,87('-'))

7010 format(//,T5,'ITERATIVE PROCEDURE ',                                    //  &
              ,T5,'Total Energy (He).......... ',0P,F18.6,' K +- ',1P,e12.4,' K',&
             /,T5,'Energy per particle (He+X). ',0P,F18.6,' K',/,                &
             /,T5,'Kinetic Energy (He)........ ',0P,F18.6,' K',                  &
             /,T5,'Lennard-Jones Energy (He).. ',0P,F18.6,' K',                  &
             /,T5,'Alpha_s term  energy (He).. ',0P,F18.6,' K',                  &
             /,T5,'Solid energy (He) ......... ',0P,F18.6,' K',                  &
             /,T5,'Correlation Energy  (He)... ',0P,F18.6,' K')
7015 format(   T5,'Impurity energy (X->He) ... ',0P,F18.6,' K',/,&
               T5,'Kinetic energy (X) ........ ',0P,F18.6,' K',/,&
               T5,'TOTAL energy (He+X) ....... ',0P,F18.6,' K',/,' ')

7016 format(1x,T5,'Impurity location  (x-axix) ..: ',F18.6,' K',/,    &
               T5,'                   (y-axis) ..: ',F18.6,' K',/,    &
               T5,'                   (z-axis) ..: ',F18.6,' K',/,' ')

7017 format(   T5,'Chemical Potential ........ ',0P,F18.6,' K +- ',1P,e12.4,'K')
7018 format(   T5,'Autovalue (impurity) ...... ',0P,F18.6,' K +- ',1P,e12.4,'K')

7100 format(/1x,T5,'Center of Mass of the Helium ...(', &
                    0P,F9.4,',',F9.4,',',F9.4,') A')
7110 format(/1x,T5,'Center of Mass of the Helium ........(', &
                    0P,F9.4,',',F9.4,',',F9.4,') A',/,       &
             1x,T5,'Center of Mass of the Impurity ......(', &
                    0P,F9.4,',',F9.4,',',F9.4,') A',/,       &
             1x,T5,'Distances between centers of mass ...(', &
                    0P,F9.4,',',F9.4,',',F9.4,') A',/,' ')

7020 format(1x,1P,2E14.5)
7030 format(0P,I10,T12,1p,E15.7,T29,1P,E9.2,T39,0P,F8.0,'/',F8.2)
7035 format(0P,I10,T12,1p,E15.7,T29,1P,E9.2,T39,E15.7,T56,1P,E9.2,T66,0P,F8.0,'/',F8.2)


8010 format('partial.density' ,i1)
8015 format('partial.densityx',i1)
8020 format('partial.density' ,i2)
8025 format('partial.densityx',i2)
8030 format('partial.density' ,i3)
8035 format('partial.densityx',i3)

5010 format('density.',SS,i1,'.out')
5020 format('density.',SS,i2,'.out')
5030 format('density.',SS,i3,'.out')
5040 format('density.',SS,i4,'.out')
5050 format('density.',SS,i5,'.out')
5060 format('density.',SS,i6,'.out')

5015 format('densityx.',SS,i1,'.out')
5025 format('densityx.',SS,i2,'.out')
5035 format('densityx.',SS,i3,'.out')
5045 format('densityx.',SS,i4,'.out')
5055 format('densityx.',SS,i5,'.out')
5065 format('densityx.',SS,i6,'.out')
!         1         2         3         4         5         6         7         8
!|2345678901234567890123456789012345678901234567890123456789012345678901234567890

end program 
!--------------------------------------------------------------------
!---                  Subroutine FForma                           ---
!--------------------------------------------------------------------

! Gives the FFT of lennard-Jones potential.
! Thre are two cores possible: 
!    OP -> Orsay-Paris  core.
!    OT -> Orsay-Trento core.


Double Precision Function FT_analitic(p,core,b,eps,sigma,hc)
implicit none

real      (kind=8) :: eps,b,sigma,hc,p
character  (len=2) :: core

real      (kind=8)  :: pi,twopi
real      (kind=8)  :: aux1,aux2,aux3,sigma6,v0
real      (kind=8)  :: t1,t2,t3

!........................
!... Useful constants ...
!........................

pi    = 4.0d0*datan(1.0d0)
twopi = 2.0d0*pi

!...............................................................................
!... calculo de la transformada del potencial en p
!...............................................................................

sigma6 = sigma**6
aux1   = 8.d0*eps
aux3   = sigma6/hc**6

       if(abs(p).lt.1.d-10) then
          FT_analitic  = bforce(core,eps,sigma,hc)
       else
          aux2 = twopi*p
          t1   =  sinx11(aux2,hc)
          t2   =   sinx5(aux2,hc)
          if(core.eq.'OP') then    !........................ Core Orsay-Paris
             v0                = aux1*aux3*(aux3-1.0d0)/hc**4
             t3                = intcore(aux2,hc)
             FT_analitic = (aux1*sigma6*(sigma6*t1-t2)+v0*t3)/p
          else                     !........................ Core Orsay-Trento
             FT_analitic = (aux1*sigma6*(sigma6*t1-t2))/p
          end if
       end if
return

!
!.. Subrutinas internas del calculo de la transformada del Lennard-Jones
!   (sinxn, cosxn, sinint, cosint, factorial)

contains

!--------------------------------------------------------------------
!                  Function sinx11
!--------------------------------------------------------------------

function sinx11(a,x) result(res)
 
! Esta funcion calcula la integral \int{ \frac{\sin(ax)}{x^11} }
! entre x e infinito. La expresion fue calculada con Mathematica
 
implicit none

real    (kind=8)             :: a,x
real    (kind=8)             :: b(5),c(5),d(5)
real    (kind=8)             :: res,res1,res2
real    (kind=8)             :: y,y2
integer (kind=4)             :: ifail


y = a*x
y2= y*y

b(1) = 40320.d0 ; c(1) = 362880.d0  ; d(1) =   1.d0   
b(2) =  -720.d0 ; c(2) =  -5040.d0  ; d(2) =   y2
b(3) =    24.d0 ; c(3) =    120.d0  ; d(3) =   y2**2
b(4) =    -2.d0 ; c(4) =     -6.d0  ; d(4) =   y2**3
b(5) =     1.d0 ; c(5) =      1.d0  ; d(5) =   y2**4

res1 = (a*dot_product(b,d)*cos(y) + dot_product(c,d)*sin(y)/x) &
       /(3628800*x**9)
res2 = -(a**10*sinint(a,x))/3628800.d0

res = res1+res2
return
end function sinx11

!--------------------------------------------------------------------
!                  Function sinx5
!--------------------------------------------------------------------

function sinx5(a,x) result(res)
 
! Esta funcion calcula la integral \int{ \frac{\sin(ax)}{x^11} }
! entre x e infinito. La expresion fue calculada con Mathematica
 
implicit none

real    (kind=8), intent(in) :: a,x
real    (kind=8)             :: res1,res2,res
real    (kind=8)             :: y,y2
real    (kind=8)             :: b(2),c(2),d(2)
integer (kind=4) :: ifail


y  = a*x
y2 = y*y

b(1) = -2.d0   ; c(1) = -6.d0 ; d(1) = 1.d0
b(2) =  1.d0   ; c(2) =  1.d0 ; d(2) = y2


res1 = -(a*dot_product(b,d)*cos(y) + dot_product(c,d)*sin(y)/x) &
       /(24*x**3)
res2 =  (a**4*sinint(a,x))/24.d0


res  =  res1+res2

return
end function sinx5

!--------------------------------------------------------------------
!                  Function sinint
!--------------------------------------------------------------------

function sinint(a,x) result(res)

! Funcion seno integral definido como Schaum (viejo) 1.285 p230

implicit none

real    (kind=8)             :: res
real    (kind=8)             :: a,x
real    (kind=8)             :: y,y2,yact
real    (kind=8)             :: term,signo
real    (kind=8), parameter  :: aux0 = -1
real    (kind=8)             :: halfpi,pi
integer (kind=4), parameter  :: maxn   = 500! 'Solo 51 terminos'
integer (kind=4)             :: n,aux1
integer (kind=4)             :: ifail 

interface
  double precision function s13adf(x,ifail)
         real (kind=8) :: x
         integer (kind=4) :: ifail
  end function s13adf
end interface

pi     = 4.0d0*datan(1.0d0)
halfpi = pi*0.5d0

y     = (a*x)


!y2    = y*y
!yact  = y
!res   = 0.0d0
 
!signo = -1.0d0
!do n=1,maxn
! signo  = signo*aux0
!  aux1   = (2*n-1)
!  term   = signo*yact/(aux1*factorial(aux1))
!  res    = res+term
!  if(abs(term/res).lt.1.d-16) exit
!  yact   = yact*y2
!end do


ifail=0
res = s13adf(y,ifail)

if(a.gt.0.0d0) then
  res = halfpi-res
else
  res = halfpi+res
end if
return
end function sinint

!--------------------------------------------------------------------
!                  Function factorial
!--------------------------------------------------------------------

! Ojo no tiene limitaciones. para 100! el orden es 10**157

recursive function factorial(n) result(res)
implicit none
integer, intent(in) :: n
real    (kind=8)    :: res
if(n==1) then
  res = 1.0d0
else
  res = n*factorial(n-1)
end if
return
end function factorial


!--------------------------------------------------------------------
!                  Function intcore
!--------------------------------------------------------------------

function intcore(a,x) result(res)
 
! Esta funcion calcula la integral la transformada de fourier del 
! core Orsay-Paris del potencial de Lennard-Jones
 
implicit none

real (kind=8) :: a,x
real (kind=8) :: b(3),c(3),d(3),e(3)
real (kind=8) :: y,y2
real (kind=8) :: res

y  = a*x
y2 = y*y

b(1) = -120.d0 ; c(1) = 120.d0  ;  d(1) = x/a**5   ; e(1) = 1.d0/a**6
b(2) =   20.d0 ; c(2) = -60.d0  ;  d(2) = d(1)*y2  ; e(2) = e(1)*y2
b(3) =   -1.d0 ; c(3) =   5.d0  ;  d(3) = d(2)*y2  ; e(3) = e(2)*y2


res  = (dot_product(b,d)*cos(y) + dot_product(c,e)*sin(y) )

return

end function intcore

!-----------------------------------------------------------------
!--               Functio bforce                               ---
!-----------------------------------------------------------------
!
!.. Esta funcion integral el Potencial de Lennard-Jones con Core
!   Orsay Paris y con el core de Orsay-Trento
!   en coordenadas esfericas a todo el espacio.
!   (sirve para comprobar el valor de 'b')

function bforce(core,eps,sigma,hc) result(res)
implicit none
character (len=2) :: core
real (kind=8) :: eps4,xx,v0,s6,s12,bb1,bb2,bb3,pi,pi4,res
real (kind=8) :: eps,hc,sigma

pi      = 4.0d0*datan(1.0d0)
pi4     = 4.0d0*pi
eps4    = 4.0d0*eps
xx      = (sigma/hc)**6
v0      = xx*(xx-1.0d0)/hc**4
s6      = sigma**6
s12     = s6*s6
if(core.eq.'OP') then
  bb1     =  pi4*eps4*(v0*hc**7)/7.d0
  bb2     =  pi4*eps4*s12/(9.d0*hc**9)
  bb3     = -pi4*eps4*s6 /(3.d0*hc**3)
  res     = bb1+bb2+bb3
else
  bb2     =  pi4*eps4*s12/(9.d0*hc**9)
  bb3     = -pi4*eps4*s6 /(3.d0*hc**3)
  res     = bb2+bb3
end if
return
end function bforce

end 
