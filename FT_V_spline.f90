Module For_FT_V_spline
  Parameter(N_intervals=4)                      ! Intervals que fem servir per separar les zones del spline
  Integer (Kind=4), Save :: N_p = 0
  Real    (Kind=8), Save :: S1 = 0.d0, Sn = 0.d0
  Real    (Kind=8), Save, Allocatable ::  P_p(:), Vpp(:), A_Vpp(:), B_Vpp(:), C_Vpp(:), D_Vpp(:)
  Integer (Kind=4) :: N_FT_V_spline(N_intervals), N_it_max_FT_VSpline=100
  Real    (Kind=8) :: Fi(N_intervals), Eps=1.d-10, Eps2=1.d-8, dr_FT_V_Spline=5.0d-1
  Character (Len=80), Save  :: Selec_FT_V_spline
  Logical :: Lprint_FT_V_spline=.false., Lparabolic_fit_FT_V_spline=.false.
  Data  N_FT_V_spline/4001,4001,4001,40001/ ! Valors dels punts en els N_intervals, definits per fer l'spline
  Data  f_value/-1.d-4/                     ! Valor del potencial per trobar el tercer punt de tall
  Data  r_max/2000.0d0/                     ! Valor supossat com infinit
End Module For_FT_V_spline
!
!  Funció que ens dona la T.F. d'una funció tipus L.J. (mirar el directori test/test_FFT_Real_to_Real pel seu funcionament practic)
!--------------------------------------------------------------------------!
!--------------------------------------------------------------------------!

Double Precision Function FT_V_spline(Pp,Pmax,H_p,Selec,r_cutoff,f_cutoff) !

!--------------------------------------------------------------------------!
!--------------------------------------------------------------------------!
!                                                          
!       Pp: Moment on volem saber la T.F.
!     Pmax: Màxim moment on voldrem que es calculi la T.F.
!      H_P: Pas amb que calcularem la T.F. en el espai de moments
!    Seelc: Nom amb el que selecionarem la función en el espai de coordenades
! r_cutoff: Punt de tall per sota del qual soposarem que función es constant i val f_cutoff
! f_cutoff: Valor de la funcio per r=r_cutoff
!
Use For_FT_V_spline
Implicit Real*8(A-H,O-Z)
Integer (Kind=4) :: N(N_intervals)
Real  (Kind=8) :: Pp, Pmax, H_p, r_cutoff, f_cutoff
Character (Len=80)  :: Selec 
Intent (In) Pp, Pmax, H_p, Selec, r_cutoff, f_cutoff

Real (Kind=8), Allocatable :: X(:,:), Y(:,:), A(:,:), B(:,:), C(:,:), D(:,:)

Data Pi/3.141592653589793238462643D0/
!
!  Distribucio de punts en els diferents intervals
!
N = N_FT_V_spline
!
!  Maxim nombre d'iteracions per trobar els punts de tall
!
N_It_max=N_it_max_FT_VSpline

If(Lprint_FT_V_spline)Then
  Open(Unit=10,File='FT_V_spline.out.p')
Endif

100 Continue

   N_p_Aux = (Pmax/H_p + 1.5)

If(N_p_Aux.Eq.N_p.And.selec.Eq.Selec_FT_V_spline)Then
!
!  Executem l'spline realitzat anteriorment (veure just per sota del Else)
!
  Call Findi(P_p,pp,N_p,k)
  FT_V_spline = A_Vpp(k) + pp*( B_Vpp(k) + pp*( C_Vpp(k) + pp*D_Vpp(k) ) )
  Return
Else
  Selec_FT_V_spline = selec
  N_P = N_p_Aux
!
!  Calculem els valors de Vpp, per realitzar l'spline posteriorment
!
  Allocate (P_p(N_p)); Allocate(Vpp(N_p)); Allocate(A_Vpp(N_p)); Allocate(B_Vpp(N_p)); Allocate(C_Vpp(N_p)); Allocate(D_Vpp(N_p));
  N_max = 0
  Do k = 1, N_intervals
    N_max=Max(N_max,N(k))
  EndDo

Allocate(X(N_intervals,N_max))  ! Matrius definides pels spline
Allocate(Y(N_intervals,N_max))  !              "
Allocate(A(N_intervals,N_max))  !              "
Allocate(B(N_intervals,N_max))  !              "
Allocate(C(N_intervals,N_max))  !              "
Allocate(D(N_intervals,N_max))  !              "

!
!  Calculem el zero del potencial
!

h=1.d-3
h2=h*h

Aux = 10.*h
x_0=r_cutoff + Aux

Step=1.0d0
i = 0
Do while(Abs(Step).Gt.Eps)
  i = i + 1
  xx = x_0 - h
  fm1=Select_Pot(selec,xx,r_cutoff,f_cutoff)
  xx = x_0 + h
  fp1=Select_Pot(selec,xx,r_cutoff,f_cutoff)
  xx = x_0
  f=Select_Pot(selec,xx,r_cutoff,f_cutoff)
  df=(fp1-fm1)/(2.*h)
  step=f/df
  x_0 = x_0 - step
  if(i.Gt.N_it_max)Then
    If(Lprint_FT_V_spline)Then
      Write(6,'("From FT_V_spline: i, x_0, f, Step",I5,1p,2E15.6)')i,x_0,f,Step
    Endif
    Stop 'From FT_V_spline 001'
  Endif
EndDo
!
!
If(Lprint_FT_V_spline)Then
  xx = x_0
  f=Select_Pot(selec,xx,r_cutoff,f_cutoff)
  Write(6,'("From FT_V_spline: Valor al primer punt de tall del potencial...:",1p,2E15.6)')x_0,f
  Call Flush(6)
Endif
!
!  Calculem el mínim del potencial
!
Aux = 10.*h
xmin = x_0 + Aux
Step=1.0d0
i = 0
Do while(Abs(Step).Gt.Eps)
  i = i + 1
  xx = xmin - h
  fm1=Select_Pot(selec,xx,r_cutoff,f_cutoff)
  xx = xmin + h
  fp1=Select_Pot(selec,xx,r_cutoff,f_cutoff)
  xx = xmin
  f=Select_Pot(selec,xx,r_cutoff,f_cutoff)
  df=(fp1-fm1)/(2.*h)
  ddf=(fp1-2.*f+fm1)/h2
  step=df/ddf
  xmin = xmin - step
  if(i.Gt.N_it_max)Then
    If(Lprint_FT_V_spline)Then
      Write(6,'("From FT_V_spline: i, xmin, df, Step",I5,1p,2E15.6)')i,xmin,df,Step
    Endif
    Stop 'From FT_V_spline 002'
  Endif
EndDo
!
!
If(Lprint_FT_V_spline)Then
  xx=xmin
  f=Select_Pot(selec,xx,r_cutoff,f_cutoff)
  Write(6,'("From FT_V_spline: Mínim del potencial...:",1p,2E15.6)')xx,f
  Call Flush(6)
Endif
!
!  Calculem el punt on el potencial val f_value
!
x_value = 2.*xmin
Step=1.0d0
i = 0
Do while(Abs(Step).Gt.Eps)
  i = i + 1
  xx = x_value - h
  fm1=Select_Pot(selec,xx,r_cutoff,f_cutoff)
  xx = x_value + h
  fp1=Select_Pot(selec,xx,r_cutoff,f_cutoff)
  xx = x_value
  f=Select_Pot(selec,xx,r_cutoff,f_cutoff)
  f = f - f_value
  df=(fp1-fm1)/(2.*h)
  step=f/df
  x_value = x_value - step
  if(i.Gt.N_it_max)Then
    If(Lprint_FT_V_spline)Then
      Write(6,'("From FT_V_spline: i, x_value, f, Step",I5,1p,2E15.6)')i,x_value,f,Step
    Endif
    Stop 'From FT_V_spline 003'
  Endif
EndDo
!
!
If(Lprint_FT_V_spline)Then
  xx=x_value
  f=Select_Pot(selec,xx,r_cutoff,f_cutoff)
  Write(6,'("From FT_V_spline: Valor al segon punt de tall del potencial...:",1p,2E15.6)')x_value,f
  Call Flush(6)
Endif
!
!  Definicio dels diferents intervals dels splines
!
X(1,1)    = R_cutoff + Eps2
X(2,1)    = x_0
X(3,1)    = xmin
X(4,1)    = x_value
X(1,N(1)) = x_0
X(2,N(2)) = xmin
X(3,N(3)) = x_value
X(4,N(4)) = r_max
!
!  Calculem els splines de la funció, per a posteriorment fer l'integració analítica
!
!  Calculem la derivada segona per millorar l'spline
!

r0 = X(1,1)
r1 = r0 + dr_FT_V_Spline
r2 = r1 + dr_FT_V_Spline
f0 = Select_Pot(selec,r0,r_cutoff,f_cutoff)
f1 = Select_Pot(selec,r1,r_cutoff,f_cutoff)
f2 = Select_Pot(selec,r2,r_cutoff,f_cutoff)

S1 = (f0 - 2.*f1 + f2)/dr_FT_V_Spline**2

!
!  Cas especial del LJ de la interacció de Orsay-Trento
!
If(Trim(selec).Eq.'LJ_OT')Then
  Lparabolic_fit_FT_V_spline=.false.
Endif

Do k=1, N_intervals
  r0 = X(k,N(k))
  rm1 = r0 - dr_FT_V_Spline
  rp1 = r0 + dr_FT_V_Spline
  fm1 = Select_Pot(selec,rm1,r_cutoff,f_cutoff)
  f0  = Select_Pot(selec,r0,r_cutoff,f_cutoff)
  fp1 = Select_Pot(selec,rp1,r_cutoff,f_cutoff)

  Sn = (fm1 - 2.*f0 + fp1)/dr_FT_V_Spline**2

  h  = (X(k,N(k))-X(k,1))/(N(k)-1)
  Do i=1, N(k)
    r = (i-1)*h + X(k,1)
    x(k,i)=r
    y(k,i)=Select_Pot(selec,r,r_cutoff,f_cutoff)
  EndDo
  Call Spline(X(k,:),Y(k,:),A(k,:),B(k,:),C(K,:),D(k,:),S1,Sn,N(k))
  S1 = Sn
EndDo
Qpi=4.0d0*Pi
Pio2= Pi*0.5d0
r0  = r_cutoff
f0  = f_cutoff
!
!  Per si volem ajustar a un paràbola invertida fins l'origen
!
If(Lparabolic_fit_FT_V_spline)Then
  r0 = r_cutoff
  f0 = f_cutoff
  r1 = r0 + dr_FT_V_Spline
  f1 = Select_Pot(selec,r1,r_cutoff,f_cutoff)
  aa = (f0-f1)/(r0**2-r1**2)
  cc = (r0**2*f1-r1**2*f0)/(r0**2-r1**2)
!
!  si la paràbola no és invertida hem de parar
!
  If(aa.Gt.0d0)Then
    Write(6,'("From FT_V_Spline: the parabola must be inverted..(a*r^2+b^r+c)..,a,b,c....:",/,1p,3E18.10)')aa,bb,cc
    Stop 'FromFT_V_Spline  004'
  Endif
Endif
rf  = x(4,N(4))
!
!  Aqui començema calcular la T.F. analítica amb els coeficients dels diferents splines
!
Do i=1,N_p
  p  = (i-1)*H_p
  P_p(i) = p
  q  = 2.0d0*Pi*P
  q2 = q*q
  q3 = q2*q
  c1oq=1.0d0/q
  c1oq2=c1oq*c1oq
  c1oq3=c1oq2*c1oq
  c1oq4=c1oq2*c1oq2
  c1oq5=c1oq4*c1oq
!
!  Cas especial per p=0
!
  If(P.eq.0.0d0)Then
    If(Lparabolic_fit_FT_V_spline)Then
      Vp =  aa*r0**5*0.2d0 + cc*r0**3/3.d0
    Else
      Vp = r0**3*f0/3.0d0
    Endif
    Fi =0.d0
!
! Integrem analiticament els splines
!
    Do k=1,N_intervals
      Do j=1,N(k)-1
        Fi(k) = Fi(k) + a(k,j)/3.0d0*(x(k,(j+1))**3-x(k,j)**3)          &
                      + b(k,j)/4.0d0*(x(k,(j+1))**4-x(k,j)**4)          &
                      + c(k,j)/5.0d0*(x(k,(j+1))**5-x(k,j)**5)          &
                      + d(k,j)/6.0d0*(x(k,(j+1))**6-x(k,j)**6)
      EndDo
      Fi(k) = Fi(k)*Qpi
    EndDo
    Sto1 = Qpi*Vp
    Sto2 =  Fi(1) + Fi(2) + Fi(3) + Fi(4)
    Vp = Sto1 +  Sto2
    If(Lprint_FT_V_spline)Then
      Write(10,'(F10.3,1p,7E15.6)')p, Vp, Sto1, (Fi(j), j=1, N_intervals)
    Endif
    Vpp(i) = Vp
  Else
!
!  Casos per p#0
!
    If(Lparabolic_fit_FT_V_spline)Then
      Vp = c1oq*(Cos(q*r0)*(6.*aa*r0*c1oq2 - aa*r0**3 - cc*r0) +                                 &
                 Sin(q*r0)*((3.*aa*r0**2+cc)*c1oq -6.*aa*c1oq3))
    Else
      Vp = f0*c1oq*(Sin(q*r0)*c1oq - r0*Cos(q*r0))
    Endif
    Fi = 0.0d0
!
! Integrem analiticament els splines
!
    Do k=1, N_intervals
      Do j=1, N(k)-1
        rj = x(k,j)
        rj1 = x(k,j+1)
        sqrj = Sin(q*rj)
        cqrj = Cos(q*rj)
        sqrj1 = Sin(q*rj1)
        cqrj1 = Cos(q*rj1)
        Fi(k) = Fi(k) + (A(k,j) + (B(k,j) + (C(k,j) + D(k,j)*rj)*rj)*rj)*rj*c1oq*cqrj            &
                      - (A(k,j) + (B(k,j) + (C(k,j) + D(k,j)*rj1)*rj1)*rj1)*rj1*c1oq*cqrj1       &
                      + (A(k,j) + (2.*B(k,j) + (3.*C(k,j) + 4.*D(k,j)*rj1)*rj1)*rj1)*c1oq2*sqrj1 &
                      - (A(k,j) + (2.*B(k,j) + (3.*C(k,j) + 4.*D(k,j)*rj)*rj)*rj)*c1oq2*sqrj     &
                      + (2.*B(k,j) + (6.*C(k,j) + 12.*D(k,j)*rj1)*rj1)*c1oq3*cqrj1               &
                      - (2.*B(k,j) + (6.*C(k,j) + 12.*D(k,j)*rj)*rj)*c1oq3*cqrj                  &
                      + (6.*C(k,j) + 24.*D(k,j)*rj)*C1oq4*sqrj                                   &
                      - (6.*C(k,j) + 24.*D(k,j)*rj1)*C1oq4*sqrj1                                 &
                      + 24.*D(k,j)*c1oq5*(cqrj - cqrj1)
      EndDo
      Fi(k) = Fi(k)*C1oq*Qpi
    EndDo
    Sto   = Vp*QPi*c1oq
    Vp    = Sto + Fi(1) + Fi(2) + Fi(3) + Fi(4)
    If(Lprint_FT_V_spline)Then
      Write(10,'(F10.3,1p,7E15.6)')p, Vp, Sto, (Fi(j), j=1, N_intervals)
    Endif
    VPP(i) = Vp
  Endif
EndDo
!
!  Fem l'spline en l'espài de moments per poder-lo emprar posteriorment (veure les instruccions a partir del 100 Continue)
!
  Call Spline(P_p,Vpp,A_Vpp,B_Vpp,C_Vpp,D_Vpp,S1,Sn,N_p)

!
! Anem al començament per fer l'spline
!
  Go to 100
Endif
Return
End
