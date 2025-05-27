!> @ingroup Library
!> @{
!> @defgroup Lib_Thermodynamic_LawsLibrary Lib_Thermodynamic_Laws
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Lib_Thermodynamic_LawsPublicProcedure Lib_Thermodynamic_Laws
!> @}

!> This module contains the definition of Thermodynamic Laws for ideal thermally perfect gas.
!> @ingroup Lib_Thermodynamic_LawsLibrary
module MOSE_Lib_Thermodynamic_Laws
  use MOSE_Global_m
  implicit none
  public:: E0, H0
  public:: prim2cons, cons2prim
  !-----------------------------------------------------------------------------------------------------------------------------------

  integer                 :: Tmin,Tmax    ! Extreme temperatures in tables
  real(kind=8), parameter :: Runiv = 8314.51d0

  real(kind=8), dimension(:), allocatable   :: wm_tab, Ri_tab
  real(kind=8), dimension(:,:), allocatable :: cp_tab, dcpi_tab, h_tab, s_tab, mi_tab, k_tab
  real(kind=8), dimension(:,:), allocatable :: dij_tab  !> coefficiente di diffusione binaria (T,interazione)
  integer                                   :: inter    !> numero di interazioni tra le specie in miscela

contains

  !> Function for computing the total specific entalpy (per unit of mass)
  function H0(p,r,u) result(entalpy)
    !---------------------------------------------------------------------------------------------------------------------------------
    implicit none
    real(kind=8), intent(IN):: p       !< Pressure.
    real(kind=8), intent(IN):: r(:)    !< Density (\f$\rho\f$).
    real(kind=8), intent(IN):: u       !< Module of velocity vector.
    real(kind=8)            :: entalpy !< Total specific entalpy (per unit of mass).
    !---------------------------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------------------------------------------
    entalpy = H(p,r)+0.5d0*u*u
    !---------------------------------------------------------------------------------------------------------------------------------
  endfunction H0

  !> Function for computing the static specific entalpy (per unit of mass)
  function H(p,r) result(entalpy)
    !---------------------------------------------------------------------------------------------------------------------------------
    implicit none
    real(kind=8), intent(IN):: p       !< Pressure.
    real(kind=8), intent(IN):: r(:)    !< Density (\f$\rho\f$).
    real(kind=8)            :: entalpy !< Static specific entalpy (per unit of mass).
    integer  :: s
    real(kind=8) :: rho, h_, Rgas
    !---------------------------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------------------------------------------
    entalpy = 0.d0
    rho = sum(r)
    Rgas = co_Rtot(r)
    do s = 1, nsc
      h_ = comp_ms_tab(p,rho,Rgas,s,h_tab)
      entalpy = entalpy+h_*r(s)/rho
    enddo
    !---------------------------------------------------------------------------------------------------------------------------------
  endfunction H


  !> Function for computing the total specific energy (per unit of mass)
  function E0(p,r,u) result(energy)
    !---------------------------------------------------------------------------------------------------------------------------------
    implicit none
    real(kind=8), intent(IN):: p      !< Pressure.
    real(kind=8), intent(IN):: r(:)   !< Density (\f$\rho\f$).
    real(kind=8), intent(IN):: u      !< Module of velocity vector.
    real(kind=8)               energy !< Total specific energy (per unit of mass).
    !---------------------------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------------------------------------------
    energy = E(p,r)+0.5d0*u*u
    !---------------------------------------------------------------------------------------------------------------------------------
  endfunction E0


  !> Function for computing the static specific energy (per unit of mass)
  function E(p,r) result(energy)
    !---------------------------------------------------------------------------------------------------------------------------------
    implicit none
    real(kind=8), intent(IN):: p      !< Pressure.
    real(kind=8), intent(IN):: r(:)   !< Density (\f$\rho\f$).
    real(kind=8)               energy !< Static specific energy (per unit of mass).
    !---------------------------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------------------------------------------
    energy = H(p,r)-p/sum(r)
    !---------------------------------------------------------------------------------------------------------------------------------
  endfunction E


  !> Function for the perfect gas equation of state
  function EOS(p,ro,R) result(T)
    !---------------------------------------------------------------------------------------------------------------------------------
    implicit none
    real(kind=8), intent(IN):: p      !< Pressure.
    real(kind=8), intent(IN):: ro     !< Density (\f$\rho\f$).
    real(kind=8), intent(IN):: R      !< Gas constant.
    real(kind=8)            :: T      !< Temperature.

    T = p/(ro*R)
    
  endfunction EOS


  !> Function to compute the conservative variables starting from the primitive ones.
  function cons2prim(cons,T) result(prim)
    !---------------------------------------------------------------------------------------------------------------------------------
    implicit none
    real(kind=8), intent(IN) :: cons(:)
    real(kind=8), intent(IN) :: T
    real(kind=8)             :: prim(size(cons))
    real(kind=8) :: Rgas, rho, rhoi(1:nsc)
    real(kind=8) :: T_iter ,p_iter, cv_iter, energy_iter, diff_en
    real(kind=8) :: energy, vel!, gamma=1.39997848512154
    integer :: iter_T
    !---------------------------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------------------------------------------
    call co_Rtot(cons(1:nsc),Rgas)
    rhoi = cons(1:nsc)
    rho = sum(rhoi)
    prim(1:nsc) = rhoi                       ! density
    prim(nu:nw) = cons(nu:nw)/rho            ! velocity

    ! Newton Raphson
    vel = norm2(prim(nu:nw))
    energy = cons(np)/rho-0.5d0*vel*vel
    diff_en = 1.0d+3; iter_T = 1; T_iter = T
    do while ( (abs(diff_en/energy) >= 1d-6) .and. (iter_T <= 20) )
      p_iter = rho*Rgas*T_iter
      cv_iter = co_cp(rhoi,T_iter,rho)-Rgas
      energy_iter = E(p_iter,rhoi,vel,Rgas)-0.5*vel*vel
      diff_en = abs(energy_iter-energy)
      T_iter = T_iter+(energy-energy_iter)/(cv_iter)
      iter_T = iter_T+1
    enddo
    if (iter_T==21) write(*,*)'Warning: max NR iter reached'
    prim(np) = rho*Rgas*T_iter               ! pressure
    !prim(np) = (cons(np)-0.5*rho*prim(nu)**2)*(gamma-1)
    !---------------------------------------------------------------------------------------------------------------------------------
  end function cons2prim


  !> Subroutine to compute the primitive variables starting from the conserveative ones 
  function prim2cons(prim) result(cons)
    !---------------------------------------------------------------------------------------------------------------------------------
    implicit none
    real(kind=8), intent(IN) :: prim(:)
    real(kind=8)             :: cons(size(prim))
    real(kind=8) :: Rgas, rho, vel
    !---------------------------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------------------------------------------
    call co_Rtot(prim(1:nsc),Rgas)
    rho = sum(prim(1:nsc))
    vel = norm2(prim(nu:nw))
    cons(1:nsc) = prim(1:nsc)                             ! density
    cons(nu:nw) = rho*prim(nu:nw)                         ! momentum
    cons(np)    = rho*E(prim(np),prim(1:nsc),vel,Rgas)    ! energy
    !---------------------------------------------------------------------------------------------------------------------------------
  end function prim2cons


  function molecularWeight(rhoi) result(result)
    implicit none
    real(8), intent(in)  :: rhoi(nsc)
    real(8) :: result
    integer :: s
    real(8) :: rho
    rho = sum(rhoi)
    result = 0.d0
    do s = 1, nsc
      result = result+(rhoi(s)/rho/Wm_tab(s))
    enddo
    result = 1.d0/result
  endfunction molecularWeight


  subroutine co_Rtot(rhoi,Rtot)
    implicit none
    real(8), intent(in)  :: rhoi(nsc)
    real(8), intent(out) :: Rtot
    integer :: s
    real(8) :: rho

    rho = sum(rhoi)
    Rtot = 0.d0
    do s = 1, nsc
      Rtot = Rtot+Ri_tab(s)*rhoi(s)/rho
    enddo

  endsubroutine co_Rtot


  subroutine co_rotot_Rtot(rhoi,rho,Rtot)
    implicit none
    real(8), intent(in)  :: rhoi(nsc)
    real(8), intent(out) :: rho, Rtot
    integer :: s

    rho = sum(rhoi)
    Rtot = 0.d0
    do s = 1, nsc
      Rtot = Rtot+Ri_tab(s)*rhoi(s)/rho
    enddo

  endsubroutine co_rotot_Rtot


  function co_cp(rhoi,T,rho) result(result)
    implicit none
    real(8), intent(in) :: rhoi(nsc)
    real(8), intent(in) :: T, rho
    real(8) :: cpi(nsc)
    real(8) :: Tdiff
    integer :: s, T_i, Tint(2)
    real(8) :: result

    T_i = idint(T)
    Tdiff  = T-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    result = 0.d0
    do s = 1, nsc
      cpi(s) = comp_ms_tabT_expr(s,cp_tab,Tint,Tdiff)
      result = result+rhoi(s)/rho*cpi(s)
    enddo

  endfunction co_cp


  function co_cp_expr(rhoi,Tint,Tdiff,rho) result(result)
    implicit none
    real(8), intent(in) :: rhoi(nsc)
    real(8), intent(in) :: Tdiff, rho
    integer, intent(in) :: Tint(2)
    real(8) :: cpi(nsc)
    integer :: s
    real(8) :: result

    result = 0.d0
    do s = 1, nsc
      cpi(s) = comp_ms_tabT_expr(s,cp_tab,Tint,Tdiff)
      result = result+rhoi(s)/rho*cpi(s)
    enddo

  endfunction co_cp_expr


  function co_a(rhoi,p,rho,Rtot) result(result)
    implicit none
    real(8), intent(in) :: rhoi(nsc)
    real(8), intent(in) :: p, Rtot, rho
    real(8) :: T, cp, gam
    real(8) :: result
    
    T = p/(Rtot*rho)
    cp = co_cp(rhoi,T,rho)
    gam = cp/(cp-Rtot)
    result = dsqrt(gam*Rtot*T)

  endfunction co_a


  subroutine co_gamma(rhoi,p,rho,Rtot,gam)
    implicit none
    real(8), intent(in)  :: rhoi(nsc), p, rho, Rtot
    real(8), intent(out) :: gam
    real(8) :: T, cp

    T = p/(Rtot*rho)
    cp = co_cp(rhoi,T,rho)
    gam = cp/(cp-Rtot)

  endsubroutine co_gamma


  subroutine co_fiij(fi,mi)
    implicit none
    real(8), intent(in)  :: mi(nsc)
    real(8), intent(out) :: fi(nsc,nsc)
    real(8) :: Mi_Mj
    integer :: i, j

    do j = 1, nsc; do i = 1, nsc
        Mi_Mj = Wm_tab(i)/Wm_tab(j)
        fi(i,j)= (1+dsqrt(mi(i)/(mi(j)+1d-20))/(Mi_Mj**0.25d0))**2/dsqrt(8*(1+Mi_Mj))
    enddo; enddo

  endsubroutine co_fiij


  subroutine co_lam(rhoi,p,rho,Rtot,tab,lam)
    implicit none
    real(8), intent(in)  :: rhoi(nsc), p, rho, Rtot
    real(8), intent(out) :: lam
    real(8), dimension(nsc) :: lam_i, lam_den, Xi, mi_fiij
    real(8)                 :: fi(nsc,nsc)
    real(8)                 :: Wmtot, T, Tdiff
    real(8)                 :: tab(Tmin:Tmax,1:nsc)
    integer                 :: s ,i ,j ,T_i, Tint(2)

    T = p/(Rtot*rho)
    T_i = idint(T)
    Tdiff = T-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    Wmtot = molecularWeight(rhoi)

    !calcolo delle frazioni molari Xi, viscosità laminare da tabella, mi(s)
    do s = 1, nsc
      Xi(s) = (rhoi(s)*Wmtot)/(rho*Wm_tab(s))
      lam_i(s) = comp_ms_tabT_expr(s,tab,Tint,Tdiff)
      mi_fiij(s) = comp_ms_tabT_expr(s,mi_tab,Tint,Tdiff)
    enddo

    ! calcolo del denominatore della legge di Wilke
    call co_fiij(fi,mi_fiij)

    do i = 1, nsc
      lam_den(i)=0.d0
      do j = 1, nsc
        lam_den(i) = lam_den(i)+Xi(j)*fi(i,j)
      enddo
    enddo
    
    ! calcolo della vicosità laminare
    lam = 0.d0
    do s = 1, nsc
      lam = lam+Xi(s)*lam_i(s)/lam_den(s)
    enddo

  endsubroutine co_lam


  subroutine co_k_mi_lam_Wilke(rhoi,rho,T,milam,klam)
    implicit none
    real(8), intent(in)  :: rhoi(nsc), rho, T
    real(8), intent(out) :: milam, klam
    real(8) :: Tdiff
    integer :: T_i, Tint(2)

    T_i = idint(T)
    Tdiff  = T-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    call co_k_mi_lam_Wilke_expr(rhoi,rho,Tint,Tdiff,milam,klam)

  endsubroutine co_k_mi_lam_Wilke


  subroutine co_k_mi_lam_Wilke_expr(rhoi,rho,Tint,Tdiff,milam,klam)
    implicit none
    integer, intent(in)  :: Tint(2)
    real(8), intent(in)  :: rhoi(nsc), rho, Tdiff
    real(8), intent(out) :: milam, klam
    real(8), dimension(nsc) :: lam_den, Xi
    real(8)                 :: fi(nsc,nsc), klam_i(nsc), milam_i(nsc)
    real(8)                 :: Wmtot, inv_lam_den
    integer                 :: s ,i ,j

    Wmtot = molecularWeight(rhoi)

    ! calcolo delle frazioni molari Xi, viscosità laminare da tabella, mi(s)
    do s = 1, nsc
      Xi(s) = rhoi(s)*Wmtot/(rho*Wm_tab(s))
      milam_i(s) = comp_ms_tabT_expr(s,mi_tab,Tint,Tdiff)
      klam_i(s)  = comp_ms_tabT_expr(s,k_tab,Tint,Tdiff)
    enddo

    ! calcolo del denominatore della legge di Wilke (same for k and mi computations)
    call co_fiij(fi,milam_i)

    do i=1,nsc
      lam_den(i)=0.d0
      do j=1,nsc
        lam_den(i)=lam_den(i)+Xi(j)*fi(i,j)
      enddo
    enddo

    ! calcolo della vicosità laminare
    milam = 0.d0
    klam = 0.d0
    do s = 1, nsc
      inv_lam_den = 1/lam_den(s)
      milam = milam+Xi(s)*milam_i(s)*inv_lam_den
      klam = klam+Xi(s)*klam_i(s)*inv_lam_den
    enddo

  endsubroutine co_k_mi_lam_Wilke_expr


  subroutine co_DS(rhoi,p,rho,Rtot,Dm)
    implicit none
    real(8), intent(in)  :: rhoi(nsc), p, rho, Rtot
    real(8), intent(out) :: Dm(nsc)
    real(8), dimension(nsc)     :: Xi, Dm_den
    real(8), dimension(nsc,nsc) :: Dij
    real(8) :: Wmtot
    integer :: s, i, j, iter

    Wmtot = molecularWeight(rhoi)

    ! calcolo delle frazioni molari Xi
    do s = 1, nsc
      Xi(s) = rhoi(s)*Wmtot/(rho*Wm_tab(s))
    enddo

    ! creazione matrice di coefficienti di diffusione
    iter = 0
    do j = 1, nsc; do i = 1, nsc
        if (i/=j) then
          iter = iter+1
          Dij(i,j) = comp_ms_tab(p,rho,Rtot,iter,dij_tab)
        endif
    enddo; enddo

    !c calcolo del coefficiente di miscela
    do s = 1, nsc
      do j=1,nsc
        if (j/=s) Dm_den(s)=Dm_den(s)+Xi(j)/Dij(s,j)
      enddo
      Dm(s) = (1-Xi(s))/Dm_den(s)
    enddo

  endsubroutine co_DS


  function comp_ms_tab(p,rho,Rtot,s,tab) result(result)
    implicit none
    real(8), intent(in) :: rho, p, Rtot, tab(Tmin:Tmax,nsc)
    integer, intent(in) :: s
    real(8) :: result
    real(8) :: T
    T = p/(Rtot*rho)
    ! chiamata alla funzione di interpolazione delle proprietà di tabella
    ! le proprietà per ora dipendono solo dalla temperatura, ma in generale
    ! possono anche dipendere dalla densità, motivo per cui in input alla funzione
    ! è anche data la variabile "ro"
    result = comp_ms_tabT(T,s,tab)
  endfunction comp_ms_tab


  function comp_ms_tabT(T,sp,tab) result(result)
    implicit none
    real(8), intent(in) :: T, tab(Tmin:Tmax,nsc)
    integer, intent(in) :: sp
    real(8) :: result
    integer :: Tint(2)
    real(8) :: Tdiff
    
    Tint(1) = idint(T); Tint(2) = Tint(1)+1
    Tdiff = T-Tint(1)
    result = comp_ms_tabT_expr(sp,tab,Tint,Tdiff)

  endfunction comp_ms_tabT


  function comp_ms_tabT_expr(sp,tab,Tint,Tdiff) result(result)
    implicit none
    real(8), intent(in) :: Tdiff, tab(Tmin:Tmax,nsc)
    integer, intent(in) :: sp, Tint(2)
    real(8) :: result
    real(8) :: Vij,Viij

    Vij=tab(Tint(1),sp)       ! int(T)   <- Tint(1)
    Viij=tab(Tint(2),sp)      ! int(T)+1 <- Tint(2)
    result = Vij+(Viij-Vij)*Tdiff

  endfunction comp_ms_tabT_expr

!> @}
endmodule MOSE_Lib_Thermodynamic_Laws