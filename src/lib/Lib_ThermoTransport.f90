
module FLINT_Lib_Thermodynamic
  implicit none
  public:: E0, H0
  public:: prim2cons, cons2prim
  !-----------------------------------------------------------------------------------------------------------------------------------

  character(len=128)      :: FLINT_phase_prefix=''
  integer, parameter      :: s_str_len = 20
  real(kind=8), parameter :: Runiv = 8314.51d0  ! Universal gas constant [J/(kmol*K)]

  integer                 :: ns                 ! Number of species
  integer                 :: ne                 ! Number of atomic elements
  character(len=s_str_len), allocatable :: species_names(:)
  character(len=s_str_len), allocatable :: elements_names(:)

  character(len=s_str_len)              :: real_species_names(1)

  real(kind=8), dimension(:,:), allocatable :: species_composition

  ! Ideal gas properties tables limits
  integer                 :: Tmin, Tmax

  ! Real fluid properties tables limits and step
  real(kind=8)            :: pmin, deltap, hmin, deltah, Tmin2, deltaT

  ! Ideal gas properties tables (T, specie)
  real(kind=8), dimension(:), allocatable   :: wm_tab, Ri_tab
  real(kind=8), dimension(:,:), allocatable :: cp_tab, dcpi_tab, h_tab, s_tab
  real(kind=8), dimension(:,:), allocatable :: mi_tab, k_tab

  ! Real fluid properties tables (p, h). One species only
  real(kind=8), dimension(:,:), allocatable :: rho_tab, T_tab, dT_tab, rh_tab, hT_tab, rp_tab, sound_tab, s_tab2D, h_tab2D
  real(kind=8), dimension(:,:), allocatable :: mi_tab2D, k_tab2D

  ! Diffusion coefficients tables (T, specie-specie interaction)
  real(kind=8), dimension(:,:), allocatable :: dij_tab

contains

  !> Function for computing the total specific entalpy (per unit of mass)
  pure function H0(p,r,u) result(entalpy)
    implicit none
    real(kind=8), intent(IN):: p       !< Pressure.
    real(kind=8), intent(IN):: r(:)    !< Density (\f$\rho\f$).
    real(kind=8), intent(IN):: u       !< Module of velocity vector.
    real(kind=8)            :: entalpy !< Total specific entalpy (per unit of mass).
    entalpy = H(p,r)+0.5d0*u*u
  endfunction H0

  !> Function for computing the static specific entalpy (per unit of mass)
  pure function H(p,r) result(entalpy)
    implicit none
    real(kind=8), intent(IN):: p       !< Pressure.
    real(kind=8), intent(IN):: r(:)    !< Density (\f$\rho\f$).
    real(kind=8)            :: entalpy !< Static specific entalpy (per unit of mass).
    integer  :: s
    real(kind=8) :: rho, h_, Rgas
    entalpy = 0.d0
    rho = sum(r)
    Rgas = f_Rtot(r)
    do s = 1, ns
      h_ = f_tab(p,rho,Rgas,s,h_tab)
      entalpy = entalpy+h_*r(s)/rho
    enddo
  endfunction H


  !> Function for computing the total specific energy (per unit of mass)
  pure function E0(p,r,u) result(energy)
    real(kind=8), intent(IN):: p      !< Pressure.
    real(kind=8), intent(IN):: r(:)   !< Density (\f$\rho\f$).
    real(kind=8), intent(IN):: u      !< Module of velocity vector.
    real(kind=8)               energy !< Total specific energy (per unit of mass).
    energy = E(p,r)+0.5d0*u*u
  endfunction E0


  !> Function for computing the static specific energy (per unit of mass)
  pure function E(p,r) result(energy)
    implicit none
    real(kind=8), intent(IN):: p      !< Pressure.
    real(kind=8), intent(IN):: r(:)   !< Density (\f$\rho\f$).
    real(kind=8)               energy !< Static specific energy (per unit of mass).
    energy = H(p,r)-p/sum(r)
  endfunction E


  pure function S0i(i,T) result(s)
    implicit none
    integer, intent(IN)     :: i
    real(kind=8), intent(IN):: T
    real(kind=8)            :: s

    s = f_tabT(T,i,s_tab)

  end function S0i


  pure function H0i(i,T) result(h)
    implicit none
    integer, intent(IN)     :: i
    real(kind=8), intent(IN):: T
    real(kind=8)            :: h

    h = f_tabT(T,i,h_tab)

  end function H0i


  !> Computes the values of C_P_0, H_0 and S_0
  pure subroutine COMP_THERMO_QUANTS(temp,N_reac,C_P_0,H_0,S_0)
    use FLINT_CEA_data, only: Rr
    real(8), intent(in)  :: temp
    integer, intent(in)  :: N_reac
    real(8), intent(out) :: C_P_0(N_reac), H_0(N_reac), S_0(N_reac)
    integer :: i_reac

    do i_reac = 1, N_reac
      C_P_0(i_reac) = cpi(i_reac,temp) * wm_tab(i_reac) / (Rr)
      H_0(i_reac) = H0i(i_reac,temp) * wm_tab(i_reac) / (Rr*temp)
      S_0(i_reac) = S0i(i_reac,temp) * wm_tab(i_reac) / (Rr)
    end do

  end subroutine COMP_THERMO_QUANTS


  !> Function for the perfect gas equation of state
  pure function EOS(p,rho,R) result(T)
    implicit none
    real(kind=8), intent(IN):: p      !< Pressure.
    real(kind=8), intent(IN):: rho    !< Density (\f$\rho\f$).
    real(kind=8), intent(IN):: R      !< Gas constant.
    real(kind=8)            :: T      !< Temperature.
    T = p/(rho*R)
  endfunction EOS


  !> Check for non-physical states (negative density) and correct them if necessary
  pure subroutine check_gas_state ( rhoi, p, error )
    implicit none
    real(kind=8), intent(inout) :: rhoi(ns)
    real(kind=8), intent(in) :: p
    integer, intent(out) :: error
    integer :: s
    real(kind=8) :: T

    error = 0
    
    ! Set a minimum value for density to avoid numerical issues in the flux computation.
    do s = 1, ns
      rhoi(s) = max(rhoi(s),1d-20)
    end do

    ! Check for negative pressure and temperature out of bounds. Set error flag if non-physical state is detected.
    if (p < 0d0) error = 1
    T = EOS(p,sum(rhoi),f_Rtot(rhoi))
    if (T < Tmin .and. T>Tmax) error = 2

  end subroutine check_gas_state


  !> Function to compute the conservative variables starting from the primitive ones.
  pure function cons2prim(cons,T) result(prim)
    implicit none
    real(kind=8), intent(IN) :: cons(:)
    real(kind=8), intent(IN) :: T
    real(kind=8)             :: prim(size(cons))
    real(kind=8) :: Rgas, rho, rhoi(1:ns)
    real(kind=8) :: T_iter ,p_iter, cv_iter, energy_iter, diff_en
    real(kind=8) :: energy, vel
    integer :: iter_T, np, nui, nuf

    np = size(cons)
    nui = ns + 1
    nuf = np - 1

    Rgas = f_Rtot(cons(1:ns))
    rhoi = cons(1:ns)
    rho = sum(rhoi)
    prim(1:ns) = rhoi                       ! density
    prim(nui:nuf) = cons(nui:nuf)/rho        ! velocity

    ! Newton Raphson
    vel = norm2(prim(nui:nuf))
    energy = cons(np)/rho-0.5d0*vel*vel
    diff_en = 1.0d+3; iter_T = 1; T_iter = T
    do while ( (abs(diff_en/energy) >= 1d-6) .and. (iter_T <= 20) )
      p_iter = rho*Rgas*T_iter
      cv_iter = f_cp(rhoi,T_iter,rho)-Rgas
      energy_iter = E(p_iter,rhoi)
      diff_en = abs(energy_iter-energy)
      T_iter = T_iter+(energy-energy_iter)/(cv_iter)
      iter_T = iter_T+1
    enddo
    !if (iter_T==21) write(*,*)'Warning: max NR iter reached'
    prim(np) = rho*Rgas*T_iter               ! pressure

  end function cons2prim


  !> Subroutine to compute the primitive variables starting from the conserveative ones 
  pure function prim2cons(prim) result(cons)
    implicit none
    real(kind=8), intent(IN) :: prim(:)
    real(kind=8)             :: cons(size(prim))
    real(kind=8) :: Rgas, rho, vel
    integer      :: nui, nuf

    nui = ns + 1
    nuf = size(prim) - 1

    Rgas = f_Rtot(prim(1:ns))
    rho = sum(prim(1:ns))
    vel = norm2(prim(nui:nuf))
    cons(1:ns)   = prim(1:ns)                           ! density
    cons(nui:nuf) = rho*prim(nui:nuf)                     ! momentum
    cons(nuf+1)   = rho*E0(prim(nuf+1),prim(1:ns),vel)   ! energy

  end function prim2cons


  pure subroutine co_rotot_Rtot ( rhoi, rho, R )
    implicit none
    real(8), intent(in)  :: rhoi(ns)
    real(8), intent(out) :: rho, R
  
    rho = sum(rhoi)
    R =f_Rtot(rhoi)
  
  endsubroutine co_rotot_Rtot


  pure function f_molecularWeight(rhoi) result(result)
    implicit none
    real(8), intent(in)  :: rhoi(ns)
    real(8) :: result
    integer :: s
    real(8) :: rho

    rho = sum(rhoi)
    result = 0.d0
    do s = 1, ns
      result = result+(rhoi(s)/rho/Wm_tab(s))
    enddo
    result = 1.d0/result

  endfunction f_molecularWeight


  pure function f_Rtot(rhoi) result(result)
    implicit none
    real(8), intent(in)  :: rhoi(ns)
    real(8) :: result
    integer :: s
    real(8) :: rho

    rho = sum(rhoi)
    result = 0.d0
    do s = 1, ns
      result = result+Ri_tab(s)*rhoi(s)/rho
    enddo

  endfunction f_Rtot


  pure function cpi(i,T) result(result)
    implicit none
    integer, intent(in) :: i
    real(8), intent(in) :: T
    real(8) :: Tdiff
    integer :: T_i, Tint(2)
    real(8) :: result

    T_i = idint(T)
    Tdiff  = T-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    result = f_tabT_expr(i,cp_tab,Tint,Tdiff)

  endfunction cpi


  pure function f_cp(rhoi,T,rho) result(result)
    implicit none
    real(8), intent(in) :: rhoi(ns)
    real(8), intent(in) :: T, rho
    real(8) :: cpi(ns)
    real(8) :: Tdiff
    integer :: s, T_i, Tint(2)
    real(8) :: result

    T_i = idint(T)
    Tdiff  = T-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    result = 0.d0
    do s = 1, ns
      cpi(s) = f_tabT_expr(s,cp_tab,Tint,Tdiff)
      result = result+rhoi(s)/rho*cpi(s)
    enddo

  endfunction f_cp


  pure function f_cp_expr(rhoi,Tint,Tdiff,rho) result(result)
    implicit none
    real(8), intent(in) :: rhoi(ns)
    real(8), intent(in) :: Tdiff, rho
    integer, intent(in) :: Tint(2)
    real(8) :: cpi(ns)
    integer :: s
    real(8) :: result

    result = 0.d0
    do s = 1, ns
      cpi(s) = f_tabT_expr(s,cp_tab,Tint,Tdiff)
      result = result+rhoi(s)/rho*cpi(s)
    enddo

  endfunction f_cp_expr


  pure function f_ss(rhoi,p,rho,Rtot) result(result)
    implicit none
    real(8), intent(in) :: rhoi(ns)
    real(8), intent(in) :: p, Rtot, rho
    real(8) :: T, cp, gam
    real(8) :: result
    
    T = p/(Rtot*rho)
    cp = f_cp(rhoi,T,rho)
    gam = cp/(cp-Rtot)
    result = dsqrt(gam*Rtot*T)

  endfunction f_ss


  pure function f_gamma(rhoi,p,rho,Rtot) result(gam)
    implicit none
    real(8), intent(in)  :: rhoi(ns), p, rho, Rtot
    real(8) :: gam
    real(8) :: T, cp

    T = p/(Rtot*rho)
    cp = f_cp(rhoi,T,rho)
    gam = cp/(cp-Rtot)

  endfunction f_gamma


  pure function mass2molar(Y_species) result(x_species)
    implicit none
    real(8), intent(in) :: Y_species(ns)
    real(8)             :: x_species(ns)
    integer  :: j
    real(8)  :: n(ns)

    do j = 1, ns
       n(j) = Y_species(j) / wm_tab(j)
    end do
    x_species(:) = n(:) / sum(n)

  end function mass2molar


  pure function molar2mass(x_species) result(y_species)
    implicit none
    real(8), intent(in) :: x_species(ns)
    real(8)             :: y_species(ns)
    integer  :: j
    real(8)  :: n(ns)

    do j = 1, ns
       n(j) = x_species(j) * wm_tab(j)
    end do
    y_species(:) = n(:) / sum(n)

  end function molar2mass


  pure subroutine co_fiij(fi,mi)
    implicit none
    real(8), intent(in)  :: mi(ns)
    real(8), intent(out) :: fi(ns,ns)
    real(8) :: Mi_Mj
    integer :: i, j

    do j = 1, ns; do i = 1, ns
        Mi_Mj = Wm_tab(i)/Wm_tab(j)
        fi(i,j)= (1+dsqrt(mi(i)/(mi(j)+1d-20))/(Mi_Mj**0.25d0))**2/dsqrt(8*(1+Mi_Mj))
    enddo; enddo

  endsubroutine co_fiij


  pure function f_laminarViscosity(rhoi,p,rho,Rtot) result(mil)
    implicit none
    real(8), intent(in)     :: rhoi(ns), p, rho, Rtot
    real(8)                 :: mil
    real(8), dimension(ns)  :: mil_den, Xi, mi_fiij
    real(8)                 :: fi(ns,ns)
    real(8)                 :: Wmtot, T, Tdiff
    integer                 :: s ,i ,j ,T_i, Tint(2)

    T = p/(Rtot*rho)
    T_i = idint(T)
    Tdiff = T-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    Wmtot = f_molecularWeight(rhoi)

    !calcolo delle frazioni molari Xi, viscosità laminare da tabella, mi(s)
    do s = 1, ns
      Xi(s) = (rhoi(s)*Wmtot)/(rho*Wm_tab(s))
      mi_fiij(s) = f_tabT_expr(s,mi_tab,Tint,Tdiff)
    enddo

    ! calcolo del denominatore della legge di Wilke
    call co_fiij(fi,mi_fiij)

    do i = 1, ns
      mil_den(i)=0.d0
      do j = 1, ns
        mil_den(i) = mil_den(i)+Xi(j)*fi(i,j)
      enddo
    enddo
    
    ! calcolo della vicosità laminare
    mil = 0.d0
    do s = 1, ns
      mil = mil+Xi(s)*mi_fiij(s)/mil_den(s)
    enddo

  endfunction f_laminarViscosity


  pure subroutine co_k_mi_lam_Wilke(rhoi,rho,T,milam,klam)
    implicit none
    real(8), intent(in)  :: rhoi(ns), rho, T
    real(8), intent(out) :: milam, klam
    real(8) :: Tdiff
    integer :: T_i, Tint(2)

    T_i = idint(T)
    Tdiff  = T-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    call co_k_mi_lam_Wilke_expr(rhoi,rho,Tint,Tdiff,milam,klam)

  endsubroutine co_k_mi_lam_Wilke


  pure subroutine co_k_mi_lam_Wilke_expr(rhoi,rho,Tint,Tdiff,milam,klam)
    implicit none
    integer, intent(in)  :: Tint(2)
    real(8), intent(in)  :: rhoi(ns), rho, Tdiff
    real(8), intent(out) :: milam, klam
    real(8), dimension(ns) :: lam_den, Xi
    real(8)                 :: fi(ns,ns), klam_i(ns), milam_i(ns)
    real(8)                 :: Wmtot, inv_lam_den
    integer                 :: s ,i ,j

    Wmtot = f_molecularWeight(rhoi)

    ! calcolo delle frazioni molari Xi, viscosità laminare da tabella, mi(s)
    do s = 1, ns
      Xi(s) = rhoi(s)*Wmtot/(rho*Wm_tab(s))
      milam_i(s) = f_tabT_expr(s,mi_tab,Tint,Tdiff)
      klam_i(s)  = f_tabT_expr(s,k_tab,Tint,Tdiff)
    enddo

    ! calcolo del denominatore della legge di Wilke (same for k and mi computations)
    call co_fiij(fi,milam_i)

    do i=1,ns
      lam_den(i)=0.d0
      do j=1,ns
        lam_den(i)=lam_den(i)+Xi(j)*fi(i,j)
      enddo
    enddo

    ! calcolo della vicosità laminare
    milam = 0.d0
    klam = 0.d0
    do s = 1, ns
      inv_lam_den = 1/lam_den(s)
      milam = milam+Xi(s)*milam_i(s)*inv_lam_den
      klam = klam+Xi(s)*klam_i(s)*inv_lam_den
    enddo

  endsubroutine co_k_mi_lam_Wilke_expr


  pure subroutine co_DS(rhoi,p,rho,Rtot,Dm)
    implicit none
    real(8), intent(in)  :: rhoi(ns), p, rho, Rtot
    real(8), intent(out) :: Dm(ns)
    real(8), dimension(ns)     :: Xi, Dm_den
    real(8), dimension(ns,ns) :: Dij
    real(8) :: Wmtot
    integer :: s, i, j, iter

    Wmtot = f_molecularWeight(rhoi)

    ! calcolo delle frazioni molari Xi
    do s = 1, ns
      Xi(s) = rhoi(s)*Wmtot/(rho*Wm_tab(s))
    enddo

    ! creazione matrice di coefficienti di diffusione
    iter = 0
    do j = 1, ns; do i = 1, ns
        if (i/=j) then
          iter = iter+1
          Dij(i,j) = f_tab(p,rho,Rtot,iter,dij_tab)
        endif
    enddo; enddo

    !c calcolo del coefficiente di miscela
    do s = 1, ns
      do j=1,ns
        if (j/=s) Dm_den(s)=Dm_den(s)+Xi(j)/Dij(s,j)
      enddo
      Dm(s) = (1-Xi(s))/Dm_den(s)
    enddo

  endsubroutine co_DS


  pure function f_tab(p,rho,Rtot,s,tab) result(result)
    implicit none
    real(8), intent(in) :: rho, p, Rtot, tab(Tmin:Tmax,ns)
    integer, intent(in) :: s
    real(8) :: result
    real(8) :: T
    T = p/(Rtot*rho)
    ! chiamata alla funzione di interpolazione delle proprietà di tabella
    ! le proprietà per ora dipendono solo dalla temperatura, ma in generale
    ! possono anche dipendere dalla densità, motivo per cui in input alla funzione
    ! è anche data la variabile "ro"
    result = f_tabT(T,s,tab)
  endfunction f_tab


  pure function f_tabT(T,sp,tab) result(result)
    implicit none
    real(8), intent(in) :: T, tab(Tmin:Tmax,ns)
    integer, intent(in) :: sp
    real(8) :: result
    integer :: Tint(2)
    real(8) :: Tdiff, Tc

    ! Clamp to table range to protect against overshoots
    Tc = max(dble(Tmin), min(T, dble(Tmax-1)))
    Tint(1) = idint(Tc); Tint(2) = Tint(1)+1
    Tdiff = Tc-Tint(1)
    result = f_tabT_expr(sp,tab,Tint,Tdiff)

  endfunction f_tabT


  pure function f_tabT_expr(sp,tab,Tint,Tdiff) result(result)
    implicit none
    real(8), intent(in) :: Tdiff, tab(Tmin:Tmax,ns)
    integer, intent(in) :: sp, Tint(2)
    real(8) :: result
    real(8) :: Vij,Viij
    integer :: Ti1, Ti2

    ! Clamp indices to valid range
    Ti1 = max(Tmin, min(Tint(1), Tmax))
    Ti2 = max(Tmin, min(Tint(2), Tmax))
    Vij=tab(Ti1,sp)
    Viij=tab(Ti2,sp)
    result = Vij+(Viij-Vij)*Tdiff

  endfunction f_tabT_expr


  ! Real fluid properties interpolation function.
  ! The input variables are pressure and enthalpy, the output variable is the one stored in the table.
  pure function ph2vars ( p, h, tab2D ) result(result)
    implicit none

    ! pressure-entalpy tabulation to have wanted variables 
    real(8), intent(in) :: p, h
    real(8),intent(in) :: tab2D(0:, 0:)
    real(8) :: result
    ! Local variables
    integer :: i,j
    real(8) :: Vij, Viij, Vijj, Viijj
    real(8) :: A, B, C, D
    real(8) :: dist_p, dist_h

    i = idint( (p - pmin) / deltap )
    j = idint( (h - hmin) / deltah )

    Vij   = tab2D ( i   , j   )
    Viij  = tab2D ( i+1 , j   )
    Vijj  = tab2D ( i   , j+1 )
    Viijj = tab2D ( i+1 , j+1 )

    A = Vij
    B = ( Viij - Vij ) / deltap
    C = ( Vijj - Vij ) / deltah 
    D = ( Vij + Viijj - Viij - Vijj ) / ( deltap * deltah )

    dist_p = p - pmin - i*deltap
    dist_h = h - hmin - j*deltah

    result = A + B*dist_p + C*dist_h + D*dist_p*dist_h

  endfunction ph2vars


  ! Real fluid properties interpolation function.
  ! The input variables are pressure and temperature, the output variable is enthalpy stored in the table.
  pure function pT2h ( p, T ) result(result)
    implicit none

    ! pressure-entalpy tabulation to have wanted variables 
    real(8), intent(in) :: p, T
    real(8) :: result
    ! Local variables
    integer :: i,j
    real(8) :: Vij, Viij, Vijj, Viijj
    real(8) :: A, B, C, D
    real(8) :: dist_p, dist_T

    i = idint( (p - pmin) / deltap )
    j = idint( (T - Tmin2) / deltaT )

    Vij   = h_tab2D ( i   , j   )
    Viij  = h_tab2D ( i+1 , j   )
    Vijj  = h_tab2D ( i   , j+1 )
    Viijj = h_tab2D ( i+1 , j+1 )

    A = Vij
    B = ( Viij - Vij ) / deltap
    C = ( Vijj - Vij ) / deltaT
    D = ( Vij + Viijj - Viij - Vijj ) / ( deltap * deltaT )

    dist_p = p - pmin - i*deltap
    dist_T = T - Tmin2 - j*deltaT

    result = A + B*dist_p + C*dist_T + D*dist_p*dist_T

  endfunction pT2h

endmodule FLINT_Lib_Thermodynamic