! test-pyrolysis-gas: the pyrolysis-gas tables of FLINT_Lib_PyrolysisGas (fixtures made by
! test/thermo/pyrolysis/make_pyrolysis.py; run in test/thermo):
!  - good.dat, three tables, read onto a solver list with species the tables do not have (N2, AR, O2) and in
!    another order: the columns land on their species by name, the others are 0;
!  - table bilin (3 pressures x 4 temperatures): its columns are bilinear in (P, Q), so the interpolation gives the
!    function back: exactly at the nodes and at points where P and Q are multiples of 1/4 (dyadic values: no
!    rounding), and within 4e-15 at a point where the weights round (the bound of the evaluation, 8 roundings of
!    values below 1, 1.8e-15, with a margin);
!  - out of the table the value at the nearest point of its boundary, on each of the four sides and at a corner;
!    p = pmax and T = Tmax exactly read the last cell, not a row past it;
!  - table single (one pressure): linear in T, the same at any pressure; its species CH4, not a species of the
!    solver, has a column of zeros and is accepted;
!  - table cond (condensed 1): x_condensed interpolated as the species;
!  - pyrolysis_table_info: name, treatment, condensed column and range of each table;
!  - a file that does not exist (ios = 1) and eleven files that must be refused (ios = 3 to 11, D18: 10 and 11), each
!    with no table left loaded.
! Exit code 1 on failure.
program test
  use FLINT_Lib_PyrolysisGas
  implicit none
  integer, parameter :: NS = 7
  character(len=20), parameter :: SOLVER(NS) = [character(len=20) :: 'N2', 'H2O', 'H2', 'AR', 'C2H2,acetylene', 'CO', 'O2']
  ! the columns of table bilin, a + b P + c Q + d P Q, on the order CO, H2, C2H2,acetylene (keep in step with
  ! make_pyrolysis.py); H2O = 1 - the three
  real(8), parameter :: COEF(4,3) = reshape([0.40625d0,  1d0/64,  -1d0/128,  1d0/512,  &
                                             0.125d0,   -1d0/128,  1d0/256,  -1d0/1024, &
                                             0.25d0,     1d0/256,  1d0/64,   -1d0/512], [4, 3])
  integer, parameter :: ICO = 6, IH2 = 3, IC2H2 = 5, IH2O = 2      ! their places in SOLVER
  integer, parameter :: ICOL(3) = [ICO, IH2, IC2H2]
  real(8), parameter :: PANY(3) = [1d3, 2d5, 1d7]
  real(8), parameter :: PG(3) = [1d5, 3d5, 7d5], TG(4) = [1000d0, 1500d0, 2500d0, 3000d0]
  integer :: err, nfail, i, j
  real(8) :: y(NS), want(NS), xc, err_max, pmin, pmax, Tmin, Tmax
  character(len=64) :: name
  character(len=16) :: treatment
  logical :: cond, exact

  nfail = 0

  err = read_pyrolysis_tables('pyrolysis/good.dat', SOLVER)
  call verdict('good.dat: three tables read', err == 0 .and. pyrolysis_table_count() == 3)
  if (err /= 0) error stop 1

  ! info
  call pyrolysis_table_info(1, name, treatment, cond, pmin, pmax, Tmin, Tmax)
  call verdict('table 1: bilin, gas-only, no condensed column, 1e5..7e5 Pa, 1000..3000 K', name == 'bilin' .and. &
               treatment == 'gas-only' .and. .not. cond .and. pmin == 1d5 .and. pmax == 7d5 .and. Tmin == 1d3 .and. &
               Tmax == 3d3)
  call pyrolysis_table_info(2, name, treatment, cond, pmin, pmax, Tmin, Tmax)
  call verdict('table 2: single, strip, one pressure (pmin = pmax = 2e5 Pa)', name == 'single' .and. &
               treatment == 'strip' .and. .not. cond .and. pmin == 2d5 .and. pmax == 2d5)
  call pyrolysis_table_info(3, name, treatment, cond, pmin, pmax, Tmin, Tmax)
  call verdict('table 3: cond, condensed, with the condensed column', name == 'cond' .and. treatment == 'condensed' &
               .and. cond)

  ! bilin at its 12 nodes and at points where P and Q are multiples of 1/4: exact
  exact = .true.
  do j = 1, 3
    do i = 1, 4
      exact = exact .and. same(PG(j), TG(i), PG(j), TG(i), 0d0)
    end do
  end do
  call verdict('bilin: the 12 nodes, exactly', exact)
  exact = same(2d5, 1250d0, 2d5, 1250d0, 0d0) .and. same(5d5, 2000d0, 5d5, 2000d0, 0d0) .and. &
          same(4d5, 2750d0, 4d5, 2750d0, 0d0) .and. same(1.5d5, 2875d0, 1.5d5, 2875d0, 0d0)
  call verdict('bilin: (2e5, 1250), (5e5, 2000), (4e5, 2750), (1.5e5, 2875): the bilinear function, exactly', exact)
  call verdict('bilin: (4.2e5, 2718.28), weights that round: within 4e-15', same(4.2d5, 2718.28d0, 4.2d5, 2718.28d0, 4d-15))
  write(*,'(A,ES10.3)') '        largest difference there: ', err_max
  ! the species of the solver that the table does not have are 0, exactly
  call pyrolysis_gas_composition(1, 4.2d5, 2718.28d0, y)
  call verdict('bilin: N2, AR and O2, not in the table, are 0', y(1) == 0d0 .and. y(4) == 0d0 .and. y(7) == 0d0)

  ! clamped: each side, a corner, and p = pmax, T = Tmax exactly
  call verdict('bilin: p below the table takes p = 1e5', same(5d4, 2000d0, 1d5, 2000d0, 0d0))
  call verdict('bilin: p above the table takes p = 7e5', same(1d6, 2000d0, 7d5, 2000d0, 0d0))
  call verdict('bilin: T below the table takes T = 1000', same(2d5, 300d0, 2d5, 1000d0, 0d0))
  call verdict('bilin: T above the table takes T = 3000', same(2d5, 4000d0, 2d5, 3000d0, 0d0))
  call verdict('bilin: a corner beyond both (1e3 Pa, 5000 K) takes (1e5, 3000)', same(1d3, 5d3, 1d5, 3d3, 0d0))
  call verdict('bilin: p = pmax and T = Tmax exactly', same(7d5, 3000d0, 7d5, 3000d0, 0d0) .and. &
               same(7d5, 1250d0, 7d5, 1250d0, 0d0))

  ! single: one pressure, linear in T, any pressure
  exact = .true.
  do i = 1, 3
    call pyrolysis_gas_composition(2, PANY(i), 1500d0, y, xc)
    exact = exact .and. y(ICO) == 0.71875d0 .and. y(IH2) == 0.28125d0 .and. count(y /= 0d0) == 2 .and. xc == 0d0
  end do
  call verdict('single: (any p, 1500 K) = CO 0.71875, H2 0.28125, the rest 0, x_condensed 0', exact)
  call pyrolysis_gas_composition(2, 2d5, 500d0, y)
  exact = y(ICO) == 0.75d0
  call pyrolysis_gas_composition(2, 2d5, 3500d0, y)
  call verdict('single: clamped in T (500 K -> CO 0.75, 3500 K -> CO 0.625)', exact .and. y(ICO) == 0.625d0)

  ! cond: x_condensed interpolated
  call pyrolysis_gas_composition(3, 1.5d5, 2500d0, y, xc)
  call verdict('cond: (1.5e5, 2500) = CO 0.40625, H2 0.15625, x_condensed 0.4375', y(ICO) == 0.40625d0 .and. &
               y(IH2) == 0.15625d0 .and. xc == 0.4375d0)

  ! refusals: each code, and no table left loaded
  call refused('a file that does not exist', 'pyrolysis/none.dat', 1)
  call refused('format 2', 'pyrolysis/bad-format.dat', 3)
  call refused('species XYZ, not of the solver, with a column not 0', 'pyrolysis/bad-species.dat', 5)
  call refused('temperatures 1000 1500 1400 3000', 'pyrolysis/bad-grid.dat', 6)
  call refused('row 5 summing to 1 + 1e-6', 'pyrolysis/bad-sum.dat', 7)
  call refused('one data row short', 'pyrolysis/bad-short.dat', 8)
  call refused('a fifth value on every row, with condensed 0', 'pyrolysis/bad-extra.dat', 8)
  call refused('no end table', 'pyrolysis/bad-noend.dat', 9)
  call refused('condensed 1 with treatment gas-only', 'pyrolysis/bad-condflag.dat', 4)
  call refused('elemental mass fractions summing to 0.99', 'pyrolysis/bad-makeup.dat', 4)
  call refused('row 5 with H2O -0.04296875, summing to 1', 'pyrolysis/bad-negative.dat', 10)
  call refused('two tables named bilin', 'pyrolysis/bad-dupname.dat', 11)

  if (nfail > 0) then
    write(*,'(A,I0,A)') ' Verdict -> fail (', nfail, ' checks)'
    stop 1
  end if
  write(*,'(A)') ' Verdict -> pass'

contains

  !> The composition of table bilin at (p, T) against the bilinear function at (pf, Tf), every species within tol;
  !> err_max the largest difference.
  logical function same ( p, T, pf, Tf, tol )
    real(8), intent(in) :: p, T, pf, Tf, tol
    real(8) :: P1, Q1
    integer :: k
    call pyrolysis_gas_composition(1, p, T, y)
    P1 = (pf - 1d5) / 2d5
    Q1 = (Tf - 1000d0) / 500d0
    want = 0d0
    do k = 1, 3
      want(ICOL(k)) = COEF(1,k) + COEF(2,k) * P1 + COEF(3,k) * Q1 + COEF(4,k) * P1 * Q1
    end do
    want(IH2O) = 1d0 - want(ICO) - want(IH2) - want(IC2H2)
    err_max = maxval(abs(y - want))
    same = err_max <= tol
  end function same

  subroutine refused ( what, file, code )
    character(*), intent(in) :: what, file
    integer, intent(in)      :: code
    integer :: e
    character(len=8) :: c
    e = read_pyrolysis_tables(file, SOLVER)
    write(c, '(I0)') code
    call verdict(what//': refused, ios = '//trim(c), e == code .and. pyrolysis_table_count() == 0)
  end subroutine refused

  subroutine verdict ( what, good )
    character(*), intent(in) :: what
    logical, intent(in)      :: good
    if (good) then
      write(*,'(A)') ' [ok]   '//what
    else
      write(*,'(A)') ' [FAIL] '//what
      nfail = nfail + 1
    end if
  end subroutine verdict

end program test
