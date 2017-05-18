!------------------------------------------------------------------
!---                    MODULES                                 ---
!------------------------------------------------------------------
module vortex
logical                          :: Lvortex=.false.     ! Flag to imprint vortex
character (len =1)               :: Vortex_Axis='Z'     ! Could be: Z, Y or X
real      (Kind=8)               :: Omega=0.d0          ! Angular velocity in the corotating frame
integer   (kind=4)               :: nv=1                ! Number of vortex
real      (kind=8), allocatable  :: xv(:), yv(:), zv(:) ! Position of the initiali mprinting vortex
end module vortex
module alphasterm
real    (kind=8) , allocatable  :: denalf(:,:,:)  ! Gaussian weighted density
real    (kind=8) , allocatable  :: falfs(:,:,:)   ! f(r) function for alpha_s term
real    (kind=8) , allocatable  :: intxalf(:,:,:) ! Intermediate integral in alpha_s terms
real    (kind=8) , allocatable  :: intyalf(:,:,:) !    "
real    (kind=8) , allocatable  :: intzalf(:,:,:) !    " 
real    (kind=8) , allocatable  :: kalfs(:,:,:)   ! FFT-kernel for the gaussian.
real    (kind=8) , allocatable  :: ualphas(:,:,:)  ! Piece of field due to the alpha_s term.
end module alphasterm
!------------------------------------------------------------------
module deriva
integer (kind=4)              :: npd=7             ! Number of points for derivatives
real    (kind=8), allocatable :: dxden(:,:,:)      ! Partial derivative in X for den in Real-space
real    (kind=8), allocatable :: dyden(:,:,:)      ! Partial derivative in X for den in Real-space
real    (kind=8), allocatable :: dzden(:,:,:)      ! Partial derivative in X for den in Real-space
integer (kind=4)              :: icon=8            ! Bounday conditions for derivatives
end module deriva
!------------------------------------------------------------------
module energies
real    (kind=8)              :: etot             ! Total energy (system)
real    (kind=8)              :: etot4            ! Total energy (helium)
real    (kind=8)              :: ekin4            ! Kinetic energy (helium)
real    (kind=8)              :: ekinx            ! Kinetic energy (impurity)
real    (kind=8)              :: elj4             ! Lennard-Jones energy
real    (kind=8)              :: ealphas          ! Alpha_s term
real    (kind=8)              :: ecor4            ! Correlation contribution
real    (kind=8)              :: eimpu            ! Impurity energy
real    (kind=8)              :: esolid           ! Solid term
real    (kind=8)              :: eso=0.0d0        ! Spin-Orbit contribution
end module energies
!------------------------------------------------------------------
module field
complex (kind=8), allocatable :: hpsi(:,:,:)     ! H Psi
real    (kind=8), allocatable :: hpsix(:,:,:)     ! H Psi
real    (kind=8), allocatable :: pot4(:,:,:)     ! 'Potential' of lagrange equation
real    (kind=8), allocatable :: uext(:,:,:)     ! 'Potential due to impurities
real    (kind=8), allocatable :: penalty(:,:,:)     ! 'Potential due to impurities
real    (kind=8), allocatable :: factt(:,:,:)    ! Local factor for imginary time evolution
real    (kind=8)  :: dmu4=0.5d0    ! Valor que sumem a la mu4 per determinar si estem en la zona clasica o no
real    (kind=8), allocatable :: paflmap(:,:,:)     ! Mapping of the paflov parameter to correct large uext
real    (kind=8), allocatable :: dx2uext(:,:,:)     ! 'Potential due to impurities
real    (kind=8), allocatable :: dy2uext(:,:,:)     ! 'Potential due to impurities
real    (kind=8), allocatable :: dz2uext(:,:,:)     ! 'Potential due to impurities
integer (kind=4)              :: iron=1          ! Number of Smoothings
integer (kind=4)              :: ironx=1          ! Number of Smoothings
end module field
!------------------------------------------------------------------
module fftmodule

character (len=15)            :: fftwplan="FFTW_ESTIMATE"
real    (kind=8),target ,allocatable :: fin(:,:,:)  ! Work Array for FFT
complex (kind=8),target ,allocatable :: fout(:,:,:) ! Work Array for FFT
integer (kind=8)              :: pfftfw      ! Pointer for FFT forward
integer (kind=8)              :: pfftbk      ! Pointer for FFT bakward
integer (kind=4)              :: nthread=1   ! Number of threads
integer (kind=4)              :: npx,npy,npz ! Number of of points (axis)
real    (kind=8)              :: renor       ! Inverse of (nx*ny*nz)

end module fftmodule
!------------------------------------------------------------------
module grid
real    (kind=8), allocatable :: x(:)            ! Values in X
real    (kind=8), allocatable :: y(:)            !   Id.     Y
real    (kind=8), allocatable :: z(:)            !   Id.     Z
real    (kind=8)              :: xmax=64.0d0     ! Last  point of the grid
real    (kind=8)              :: ymax=64.0d0     ! Last  point of the grid
real    (kind=8)              :: zmax=64.0d0     ! Last  point of the grid
real    (kind=8)              :: xc=0.0d0        ! Center of the cluster. X
real    (kind=8)              :: yc=0.0d0        ! Center of the cluster. Y
real    (kind=8)              :: zc=0.0d0        ! Center of the cluster. Z
real    (kind=8)              :: hx,hy,hz        ! Steps in X, Y and Z
integer (kind=4)              :: nx=64           ! Number of points in X
integer (kind=4)              :: ny=64           ! Number of points in Y
integer (kind=4)              :: nz=64           ! Number of points in Z
integer (kind=4)              :: nxyz            ! nxyz=nx*ny*nz
real    (kind=8)              :: dxyz            ! dxyz=hx*hy*hz
end module grid
!------------------------------------------------------------------
module gridk
real    (kind=8), allocatable :: px(:)          ! Values in PX
real    (kind=8), allocatable :: py(:)          !   Id.     PY
real    (kind=8), allocatable :: pz(:)          !   Id.     PZ
real    (kind=8), allocatable :: pmod(:,:,:)    ! Module of p
real    (kind=8)              :: pmaxx,pmaxy,pmaxz
real    (kind=8)              :: hpx,hpy,hpz
end module gridk
!------------------------------------------------------------------
module he3
real    (kind=8)              :: densat3=0.0163 ! Density of saturation for 3He
!
!.......... Density functional parameter  (Orsay-Trento)
!           (From Barranco et al. PRB 56(1997)8997-9003
!
!..Helium-3  (Refs. 1 and 2)
real      (kind=8) ::    cp3 =1.588790d6   ! K \AA**{3+\gamma_3}
real      (kind=8) ::   cpp3 =    -3.5d4   ! K \AA**{3+\gamma_3}
real      (kind=8) ::   gam3 =  2.1251d0
real      (kind=8) ::  den3c =  0.0406d0   ! \AA**{-3}
real      (kind=8) :: h2o2m3 =  8.041775d0 ! \hbar**2 / (2 m_3)
real      (kind=8) ::   beta =5.55555555555556d-02 !  = 1/18 See Ref.2 p-5026
real      (kind=8) ::  betap =3.33333333333333d-01 !  = 1/3  See Ref.2 p-5026
!real     (kind=8) ::  betap =     0.0d0
end module he3
!------------------------------------------------------------------
module he4

real    (kind=8)              :: densat4=0.02184 ! Density of saturation for 4He
!
!.......... Density functional parameter  (Orsay-Trento)
!           (From Barranco et al. PRB 56(1997)8997-9003
!
!..Helium-4  (Refs. 1 and 2)
real      (kind=8)               ::    cp4=-2.41186d4       ! K \AA**6
real      (kind=8)               ::   cpp4= 1.85850d6       ! K \AA**9
real      (kind=8)               ::  den4c= 0.062d0         ! \AA**{-3}
real      (kind=8)               :: alphas= 54.31d0         ! K ^-1 \AA**3
real      (kind=8)               ::      l= 1.0d0           ! \AA
real      (kind=8)               ::  den0s= 0.04d0          ! \AA**-3
real      (kind=8)               :: h2o2m4= 6.05969638298d0 ! \hbar**2 / (2 m_4)
logical                          :: lsolid=.false. ! Sumem els termes de la funcional solida
logical                          :: lden_max=.false. ! Controlem que la funcional solida no s'en vagi de mare
real      (kind=8)               :: solid_denmax = 0.5d0
real      (kind=8)               :: C = 3.1577504d4
real      (kind=8)               :: beta= 40.d0
real      (kind=8)               :: den_m = 0.37d0
real      (kind=8)               :: Chempo = -7.15d0 ! Potencial químico del sistema homogeneo
logical                          :: lbulk=.false.    ! Etiqueta para imponer el cálculo del sistema homogeneo
!
end module he4
!------------------------------------------------------------------
module impur
real    (kind=8)                 :: Als, Als_P=0.0d0, Als_D=0.0d0    ! Spin-orbit spliting
real    (kind=8)                 :: ximp=0.0d0   ! Position of the impurity  X-ccordinate.
real    (kind=8)                 :: yimp=0.0d0   !  "        "  "   "        Y-Coordinate
real    (kind=8)                 :: zimp=0.0d0   !  "        "  "   "        Z-Coordinate
integer (kind=4)                 :: irimp(3)      ! Coordinates of the impurity.
logical                          :: limp_despl=.false. ! Desplaçament de la impureça durant la mimització
real      (kind=8)               :: pas_imp=1.0d0! Pas pel Newton-Rapson del desplaçament de l'impureça
real      (kind=8)               :: pas_imp_max=1.0d-1! Pas maxim pel Newton-Rapson
logical                          :: lrandom=.false. ! We add a random process on the start density to init localization
logical                          :: limp=.false. ! T-> Impurity, F->Pure
logical                          :: lexternal = .false. ! T-> uext included, F-> not true
real    (kind=8)                 :: umax=400.d0   ! Maximum value for the potential
real    (kind=8)                 :: r_cutoff=0.d0 ! Maximum value for the potential
character (Len=80)               :: selec=''      ! Aqui seleccionem el potencial
integer   (kind=8)               :: nq=1000      ! Number of q-values for Patil
real      (kind=8)               :: tol=1.d-7    ! Tolerance for Romberg
real      (kind=8)               :: rmin=3.d0    ! Core For Patil
complex   (kind=8),allocatable   :: vq(:,:,:)    ! Fourier of Patil Potential
real      (kind=8),allocatable   :: psix(:,:,:)  ! Wave function for the impurity
real      (kind=8),allocatable   :: psixnw(:,:,:)  ! Wave function for the impurity
real      (kind=8),target,allocatable   :: denx(:,:,:)  ! Density for the impurity
real      (kind=8),allocatable   :: psi1x(:,:,:) ! Wave function for the impurity
real      (kind=8),allocatable   :: psi2x(:,:,:) ! Wave function for the impurity
real      (kind=8),allocatable   :: upotx(:,:,:) ! Mean field for the impurity
real      (kind=8),allocatable   :: potx4(:,:,:) ! 
real      (kind=8)               :: h2o2mx       ! Term hbar**2/(2*m_x)
real      (kind=8)               :: rinfi=500        ! Infinite for Patil
real      (kind=8)               :: gwf=5.69      ! Parameter for starting gaussian
real      (kind=8)               :: r_shell=0.0d0 ! To fix a spherical w.f. at radius r_shell
!complex (kind=8), allocatable :: fpsix(:,:,:)     ! FFT of Psix
complex (kind=8),target, allocatable :: fdenx(:,:,:)     ! FFT of Psix
Logical                          :: Lprint_invar=.false.
Logical                          :: Lexcite_state=.false.
Logical                          :: Lexcite_state_external=.false.
Logical                          :: Lexcite_state_fix=.false.
Logical                          :: Lexciplex_state_fix=.false.
character (len=6)                :: Exciplex='Ring' ! Or 'Linear'
real      (kind=8)               :: r_exc=3.9d0    ! To fix the radial position of lineal or ring shape exciplex
character (len=3)                :: elem        ! Chemical Symbol for Alkali
character (len=3)                :: leepot='NO '! YES/NO Read VXQ external
character (len=40)               :: vxpot       ! Name of External file with VXQ
complex (kind=8), allocatable :: invar(:), invar0(:)
integer (kind=4)              :: ninvar =10 
integer (kind=4)              :: instate
real      (kind=8) , allocatable :: Vpi(:,:,:),Pi_Del(:,:,:),Sig_Del(:,:,:)
real      (kind=8) , allocatable :: Delta(:,:,:)
real      (kind=8)               :: Xnew=1.d0     ! Parameter to mix old & new Uext excited potentials (see instates)
Logical                          :: Lfirst=.true.
character (len = 1)              :: Lstate='P'
Logical                          :: Laverage_P_value=.false.
Logical                          :: Ldiag_jz=.false.
character (len = 4)  :: Ljz='' ! Per  Lstate='P',  -3/2, -1/2, +1/2, +3/2; Si posem '', per instate=0 | Correspond a j=1/2, Ljz ='', Ljz='-1/2', '+1/2'
                               !             per instate=1 ---> Ljz='-3/2'  | correspond a j=3/2
                               !             per instate=2 ---> Ljz='-1/2'  |
                               !
                               ! Per Lstate='D',  -5/2, -3/2, +1/2, +5/2, +3/2, +1/2; Si posem Ljz='', per instate=0 ---> Ljz='-3/2'  |  Correspond a j=3/2
                               !                                                                       per instate=1 ---> Ljz='-1/2'  |
                               !                                                                   per instate=2 ---> Ljz='-5/2'  |
                               !                                                                   per instate=3 ---> Ljz='-3/2'  |  Correspond a j=5/2
                               !                                                                   per instate=4 ---> Ljz='-1/2'  |
end module impur
!------------------------------------------------------------------
module lenard3
real      (kind=8), allocatable  ::   fvlj3(:,:,:)
real      (kind=8)               ::   h3   ! \AA

real      (kind=8)               ::     h3op =2.356415d0   ! \AA
real      (kind=8)               ::     h3ot =2.356415d0   ! \AA
real      (kind=8)               ::   eps3   =   10.22d0   ! K
real      (kind=8)               :: sigma3   =   2.556d0   ! \AA
real      (kind=8)               ::     b3   =-684.676d0   ! K \AA**3
character (len=2)                ::  core3   ='OT'

real      (kind=8), allocatable  :: delj3(:,:,:) ! Density of energy-lennard-Jones
end module lenard3
!------------------------------------------------------------------
module lenard4
real      (kind=8), allocatable  ::   fvlj4(:,:,:)
real      (kind=8)               ::   h4  

real      (kind=8)               ::     h4op =2.359665d0   ! \AA
real      (kind=8)               ::     h4ot =2.190323d0   ! \AA
real      (kind=8)               ::     eps4 =   10.22d0   ! K
real      (kind=8)               ::   sigma4 =   2.556d0   ! \AA
real      (kind=8)               ::       b4 =-718.99d0    ! K \AA**3
character (len=3)                :: core4 ='OT '
real      (kind=8), allocatable  :: delj4(:,:,:) ! Density of energy-lennard-Jones

end module lenard4
!------------------------------------------------------------------
module rho
real    (kind=8),target, allocatable :: den(:,:,:)      ! Density in Real-space
complex (kind=8), allocatable :: Psi(:,:,:)      ! Psi=sqrt(density) in real space
complex (kind=8), allocatable :: Psinw(:,:,:)    ! New value of Psi (temporal)
complex (kind=8), allocatable :: Psi1(:,:,:)     ! Old Psi 1-step  backward
complex (kind=8), allocatable :: Psi2(:,:,:)     ! Old Psi 1-step  backward
complex (kind=8),target, allocatable :: fden(:,:,:)     ! Density in K-space
!complex (kind=8), allocatable :: fpsi(:,:,:)     ! FFT of Psi
real    (kind=8), allocatable :: dencg(:,:,:)    ! Coarse-Graining density
real    (kind=8), allocatable :: wcgk(:,:,:)     ! Kernel of coarse graining
                                                 ! in fourier space
real    (kind=8)              :: denmin=1.d-99          ! Minimum value for densities.
real    (kind=8)              :: psimin=3.162277661d-50    ! Minimum value for w.f.
real    (kind=8)              :: defx=1.0d0, defy=1.0d0, defz=1.0d0  ! Deform factors for the initial density

end module rho
!------------------------------------------------------------------
module util1
character (len=1)  :: cchar="#"  ! Used in routine 'Titols'
real      (kind=8) :: pi         ! Pi=3.141592.... value
real      (kind=8) :: twopi      ! twopi  = 2*pi
real      (kind=8) :: fourpi     ! fourpi = 4*pi
real      (kind=8) :: piq        ! piq    = pi*pi
real      (kind=8) :: afermi     ! Parameter for the initial fermi-distribution
real      (kind=8) :: rfermi     ! Parameter for the initial fermi-distribution
real      (kind=8), dimension (2) :: vdt=(/0.6d0, 0.1d0/)   ! Speeds for imaginary step-time method.
integer   (kind=4) :: nn(3)      ! Auxiliar array for pderg
integer   (kind=4) :: mmx(4)     ! Auxiliar array for pderg
integer   (kind=4) :: iw(11)     ! Auxiliar array for pderg
integer   (kind=4) :: nsfiles=10 ! Number-of-save-files
integer   (kind=4) :: nsfaux=0   ! Actual generation of backup file
integer   (kind=4) :: irespar=0  ! Does not write partial auxiliar plot files...
logical            :: printpot=.false. ! Prints (no) impurity potential
logical            :: lmillorar_singularitats=.false. ! Intent de millorar les baixes densitats
Logical            :: L_anell=.false., L_esfera=.false.
real      (kind=8) :: r_anell=3.d0, r_esfera=3.0d0  ! Radi del anell al voltant del eix z
real      (kind=8) :: a_anell=0.5d0, a_esfera=0.5d0 ! Parametre per la gausiana per construir l'anell (veure readen)

end module util1
!------------------------------------------------------------------
module work1
complex (kind=8),target, allocatable :: wk1(:,:,:)
complex (kind=8),target, allocatable :: wk2(:,:,:)
complex (kind=8),target, allocatable :: wk3(:,:,:)
real    (kind=8),target, allocatable :: sto1(:,:,:)
real    (kind=8),target, allocatable :: sto2(:,:,:)
real    (kind=8),target, allocatable :: sto3(:,:,:)
real    (kind=8), allocatable :: sto4(:,:,:)
real    (kind=8), allocatable :: sto5(:,:,:)
real    (kind=8), allocatable :: sto6(:,:,:)
complex (kind=8), allocatable :: sto1c(:,:,:)
complex (kind=8), allocatable :: sto2c(:,:,:)
complex (kind=8), allocatable :: sto3c(:,:,:)
complex (kind=8), allocatable :: sto4c(:,:,:)
complex (kind=8), allocatable :: sto5c(:,:,:)
complex (kind=8), allocatable :: sto6c(:,:,:)
end module work1
!------------------------------------------------------------------
module constraint
logical       :: lconstraint=.false.
logical       :: lconstraint_imp=.false.
real (kind=8) :: Intens = 3000.d0 ! Intensity of the constraint.
real (kind=8) :: n4real ! Intensity of the constraint.
real (kind=8) :: zdist  = 0.d0    ! distance in z axis between the He com and Ag atom
real (kind=8) :: enercons
end module constraint
module rkpc
logical          :: lrkpc=.false.  ! Control to start de Steprkr & Steppcr
integer (kind=4) :: ioldp(3),ioldh(2),Icon_local=0                         ! Auxiliar arrays for steppcr
complex (kind=8), allocatable ::  q(:,:,:), pc(:,:,:)                      ! Auxiliar arrays for Steprkr & Steppcr
complex (kind=8), allocatable ::  Psiold(:,:,:,:), hpsiold(:,:,:,:)        ! Auxiliar arrays for Steppcr
real    (kind=8), allocatable ::  timec(:,:,:)                             ! Auxiliar arrays for Steppcr & Steprkr
real    (kind=8) ::  ErrPC=0.0d0, underpsil=-700.0d0                        ! Auxiliar to compute the error in Steppcr
real    (kind=8) ::  Dtmin=1.0d-10, Sdmu=1.d-1, Dmumax=1.d2                ! Auxiliar to compute the efective Dt
end module rkpc
!------------------------------------------------------------------

