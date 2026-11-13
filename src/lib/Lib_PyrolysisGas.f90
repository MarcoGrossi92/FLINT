! read_pyrolysis_tables ios map:
! 0 -> no error
! 1 -> file not found, or it cannot be opened
! 2 -> read error, or the file ends inside a table
! 3 -> the first line is not 'format 1'
! 4 -> a line of the header is not the one expected, or holds a value out of range (table count, treatment, condensed
!      flag, counts, elemental makeup that does not sum to 1, a line after the last table)
! 5 -> a species of a table twice, or absent from the solver's list with a column that is not identically 0
! 6 -> pressures or temperatures not positive and strictly ascending
! 7 -> a row whose mass fractions (and condensed fraction) do not sum to 1 within 1e-8
! 8 -> a data row that does not hold the number of values of its table (nsp, or nsp + 1 with 'condensed 1')
! 9 -> no 'end table' after the last data row of a table
! 10 -> a data row with a negative value (a mass fraction, or X_C with 'condensed 1')
! 11 -> a table with the name of a table before it (every table has its own name)


!> Pyrolysis-gas composition tables y_pg(p, T): the equilibrium composition of the gas that a pyrolysing material
!> blows through its char, tabulated on a rectilinear grid of pressures and temperatures, one table per material.
!>
!> The file (INPUT/<prefix>pyrolysis.dat, written by ATLAS GPB) is plain text: a '#' at the start of a line or after a
!> blank starts a comment, blank lines are skipped, tokens are separated by blanks, a species name is one token (commas
!> allowed):
!>
!>     format 1
!>     tables <ntables>
!>     table <name>                                  (one block per table, ntables blocks)
!>     treatment <gas-only|strip|condensed>
!>     elements <nel> <El1> <massf1> ... <Elnel> <massfnel>
!>     condensed <0|1>
!>     species <nsp>
!>     <name_1>                                      (one name per line)
!>     ...
!>     pressures <np>                                (Pa, strictly ascending; the values may run over several lines)
!>     <p_1> ... <p_np>
!>     temperatures <nT>                             (K, strictly ascending; nT >= 2)
!>     <T_1> ... <T_nT>
!>     data                                          (np blocks of nT rows, pressure-major: row (iT, ip) is line
!>     <y_1> ... <y_nsp> [X_C]                        (ip-1)*nT + iT; one row per line)
!>     end table
!>
!> 'treatment gas-only' and 'strip' have 'condensed 0'; 'treatment condensed' has 'condensed 1', and every row ends
!> with X_C, the mass fraction of the condensed carbon. Every row sums to 1 within 1e-8 (y_1 + ... + y_nsp + X_C), and
!> no value of a row is negative. Every table has its own name (a boundary condition names its table by it).
!> The columns are mapped by name onto the species list of the solver that loads the tables: a solver species absent
!> from a table gets 0; a table species absent from the solver list is an error, unless its column is identically 0.
!>
!> The tables are read once (read_pyrolysis_tables) and only read afterwards: pyrolysis_gas_composition is pure, and
!> may be called from any thread.
module FLINT_Lib_PyrolysisGas
  use iso_fortran_env, only: error_unit
  implicit none
  private
  public :: read_pyrolysis_tables, pyrolysis_table_count, pyrolysis_table_info, pyrolysis_gas_composition

  integer, parameter      :: PG_NAME_LEN = 64     ! longest table or species name kept
  real(kind=8), parameter :: PG_SUM_TOL  = 1.d-8  ! |sum of a row - 1|, and |sum of the elemental makeup - 1|

  type :: pyrolysis_table_t
    character(len=PG_NAME_LEN) :: name = ''
    character(len=16)          :: treatment = ''
    logical                    :: has_condensed = .false.
    integer                    :: np = 0, nT = 0
    real(kind=8), allocatable  :: p(:), T(:)        ! [Pa], [K], strictly ascending
    real(kind=8), allocatable  :: y(:,:,:)          ! (solver species, nT, np): the rows on the solver's species order
    real(kind=8), allocatable  :: xc(:,:)           ! (nT, np): condensed-carbon mass fraction, 0 without the column
  end type pyrolysis_table_t

  type(pyrolysis_table_t), allocatable :: tables(:)
  integer                              :: nsolver = 0   ! species of the solver's list the tables are mapped onto

contains

  !> Read every table of the file filename and map its columns onto species_names, the species list of the solver.
  !> Tables read before are dropped first. ios: see the map at the top of this file; on an error the message, with the
  !> line of the file, is already printed on the standard output and on error_unit, and no table is loaded.
  function read_pyrolysis_tables ( filename, species_names ) result ( ios )
    implicit none
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: species_names(:)
    integer                      :: ios
    ! Local
    integer :: unitfile, nline, ntab, itab, k, s, ncol, nsp, nel, kp, kT, ntok, stat
    character(len=:), allocatable :: line
    character(len=PG_NAME_LEN)    :: tok(3)
    character(len=PG_NAME_LEN), allocatable :: tspec(:)
    integer, allocatable          :: map(:)
    real(kind=8), allocatable     :: row(:), makeup(:)
    real(kind=8)                  :: total
    character(len=256)            :: msg

    if (allocated(tables)) deallocate(tables)
    nsolver = size(species_names)

    open(newunit=unitfile, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      call report('cannot open the file')
      ios = 1
      return
    end if
    nline = 0

    ! format 1
    call next_line(unitfile, line, nline, ios)
    if (ios /= 0) then; call fail(2, 'the file holds no line'); return; end if
    call split(line, tok, ntok)
    if (ntok /= 2 .or. tok(1) /= 'format' .or. tok(2) /= '1') then
      call fail(3, 'the first line is "'//trim(line)//'", not "format 1": this reader knows format 1 only'); return
    end if

    ! tables <ntables>
    call next_line(unitfile, line, nline, ios)
    if (ios /= 0) then; call fail(2, 'the file ends before "tables"'); return; end if
    call split(line, tok, ntok)
    ntab = 0
    if (ntok == 2 .and. tok(1) == 'tables') read(tok(2), *, iostat=stat) ntab
    if (ntok /= 2 .or. tok(1) /= 'tables' .or. ntab < 1) then
      call fail(4, 'expected "tables <n>", n >= 1, found "'//trim(line)//'"'); return
    end if
    allocate(tables(ntab))

    do itab = 1, ntab
      associate ( tb => tables(itab) )

        ! table <name>
        call next_line(unitfile, line, nline, ios)
        if (ios /= 0) then; call fail(2, 'the file ends before table '//itoa(itab)//' of '//itoa(ntab)); return; end if
        call split(line, tok, ntok)
        if (ntok /= 2 .or. tok(1) /= 'table') then
          call fail(4, 'expected "table <name>" (table '//itoa(itab)//' of '//itoa(ntab)//'), found "'//trim(line)//'"')
          return
        end if
        tb%name = tok(2)
        if (any(tables(1:itab-1)%name == tb%name)) then
          call fail(11, 'table '//itoa(itab)//' of '//itoa(ntab)//' is named '//trim(tb%name)//', as a table before'// &
                    ' it: every table has its own name (a boundary condition names its table by it)'); return
        end if

        ! treatment <gas-only|strip|condensed>
        call header('treatment')
        if (ios /= 0) return
        if (ntok /= 2 .or. (tok(2) /= 'gas-only' .and. tok(2) /= 'strip' .and. tok(2) /= 'condensed')) then
          call fail(4, 'table '//trim(tb%name)//': expected "treatment gas-only|strip|condensed", found "'// &
                    trim(line)//'"'); return
        end if
        tb%treatment = tok(2)

        ! elements <nel> <El1> <massf1> ...: the elemental makeup, mass fractions that sum to 1
        call header('elements')
        if (ios /= 0) return
        nel = 0
        if (ntok >= 2) read(tok(2), *, iostat=stat) nel
        if (ntok < 2 .or. nel < 1) then
          call fail(4, 'table '//trim(tb%name)//': expected "elements <nel> <El> <massf> ...", nel >= 1, found "'// &
                    trim(line)//'"'); return
        end if
        allocate(makeup(nel))
        if (.not. values_after(line, 2, 2*nel, makeup, pairs=.true.)) then
          call fail(4, 'table '//trim(tb%name)//': "elements '//itoa(nel)//'" needs '//itoa(nel)// &
                    ' pairs of an element and its mass fraction, found "'//trim(line)//'"'); deallocate(makeup); return
        end if
        total = sum(makeup)
        deallocate(makeup)
        if (.not. (abs(total - 1.d0) <= PG_SUM_TOL)) then
          write(msg, '(ES24.16)') total
          call fail(4, 'table '//trim(tb%name)//': the elemental mass fractions sum to '//trim(adjustl(msg))// &
                    ', not 1 within 1e-8'); return
        end if

        ! condensed <0|1>, consistent with the treatment
        call header('condensed')
        if (ios /= 0) return
        if (ntok /= 2 .or. (tok(2) /= '0' .and. tok(2) /= '1')) then
          call fail(4, 'table '//trim(tb%name)//': expected "condensed 0" or "condensed 1", found "'//trim(line)//'"')
          return
        end if
        tb%has_condensed = tok(2) == '1'
        if (tb%has_condensed .neqv. (tb%treatment == 'condensed')) then
          call fail(4, 'table '//trim(tb%name)//': treatment '//trim(tb%treatment)//' with "'//trim(line)// &
                    '": treatment condensed has condensed 1, gas-only and strip condensed 0'); return
        end if

        ! species <nsp>, then one name per line, mapped onto the solver's list
        call header('species')
        if (ios /= 0) return
        nsp = 0
        if (ntok == 2) read(tok(2), *, iostat=stat) nsp
        if (ntok /= 2 .or. nsp < 1) then
          call fail(4, 'table '//trim(tb%name)//': expected "species <nsp>", nsp >= 1, found "'//trim(line)//'"')
          return
        end if
        allocate(tspec(nsp), map(nsp))
        do k = 1, nsp
          call next_line(unitfile, line, nline, ios)
          if (ios /= 0) then; call fail(2, 'table '//trim(tb%name)//': the file ends inside the species list'); return
          end if
          call split(line, tok, ntok)
          if (ntok /= 1) then
            call fail(4, 'table '//trim(tb%name)//': expected species '//itoa(k)//' of '//itoa(nsp)// &
                      ', one name, found "'//trim(line)//'"'); return
          end if
          tspec(k) = tok(1)
          if (any(tspec(1:k-1) == tspec(k))) then
            call fail(5, 'table '//trim(tb%name)//': species '//trim(tspec(k))//' given twice'); return
          end if
          map(k) = 0
          do s = 1, nsolver
            if (trim(adjustl(species_names(s))) == trim(tspec(k))) then
              map(k) = s
              exit
            end if
          end do
        end do

        ! pressures <np>, temperatures <nT>: positive, strictly ascending
        call grid('pressures', 1, tb%p, tb%np)
        if (ios /= 0) return
        call grid('temperatures', 2, tb%T, tb%nT)
        if (ios /= 0) return

        ! data: np blocks of nT rows, one row per line
        call header('data')
        if (ios /= 0) return
        if (ntok /= 1) then
          call fail(4, 'table '//trim(tb%name)//': expected "data", found "'//trim(line)//'"'); return
        end if
        ncol = nsp
        if (tb%has_condensed) ncol = nsp + 1
        allocate(row(ncol))
        allocate(tb%y(nsolver, tb%nT, tb%np), tb%xc(tb%nT, tb%np))
        tb%y  = 0.d0
        tb%xc = 0.d0
        do kp = 1, tb%np
          do kT = 1, tb%nT
            call next_line(unitfile, line, nline, ios)
            if (ios /= 0) then
              call fail(2, 'table '//trim(tb%name)//': the file ends at data row '//itoa((kp-1)*tb%nT + kT)//' of '// &
                        itoa(tb%np*tb%nT)); return
            end if
            if (count_tokens(line) /= ncol) then
              call fail(8, 'table '//trim(tb%name)//': data row '//itoa((kp-1)*tb%nT + kT)//' of '// &
                        itoa(tb%np*tb%nT)//' holds '//itoa(count_tokens(line))//' values, not '//itoa(ncol)// &
                        ' (the species, and X_C with condensed 1)'); return
            end if
            read(line, *, iostat=stat) row
            if (stat /= 0) then
              call fail(8, 'table '//trim(tb%name)//': data row '//itoa((kp-1)*tb%nT + kT)//' holds a value that is'// &
                        ' not a number: "'//trim(line)//'"'); return
            end if
            if (.not. (abs(sum(row) - 1.d0) <= PG_SUM_TOL)) then
              write(msg, '(ES24.16)') sum(row)
              call fail(7, 'table '//trim(tb%name)//': data row '//itoa((kp-1)*tb%nT + kT)//' (p = '// &
                        rtoa(tb%p(kp))//' Pa, T = '//rtoa(tb%T(kT))//' K) sums to '//trim(adjustl(msg))// &
                        ', not 1 within 1e-8'); return
            end if
            if (any(row < 0.d0)) then
              k = minloc(row, 1)
              msg = 'X_C'
              if (k <= nsp) msg = tspec(k)
              call fail(10, 'table '//trim(tb%name)//': data row '//itoa((kp-1)*tb%nT + kT)//' (p = '// &
                        rtoa(tb%p(kp))//' Pa, T = '//rtoa(tb%T(kT))//' K) holds '//rtoa(row(k))//' for '//trim(msg)// &
                        ': a mass fraction is 0 or positive'); return
            end if
            do k = 1, nsp
              if (map(k) > 0) then
                tb%y(map(k), kT, kp) = row(k)
              else if (row(k) /= 0.d0) then
                call fail(5, 'table '//trim(tb%name)//': species '//trim(tspec(k))//' is not a species of the '// &
                          'solver, and its mass fraction is not 0 (data row '//itoa((kp-1)*tb%nT + kT)//')'); return
              end if
            end do
            if (tb%has_condensed) tb%xc(kT, kp) = row(ncol)
          end do
        end do
        deallocate(row, tspec, map)

        ! end table
        call next_line(unitfile, line, nline, ios)
        call split(line, tok, ntok)
        if (ios /= 0 .or. ntok /= 2 .or. tok(1) /= 'end' .or. tok(2) /= 'table') then
          if (ios /= 0) line = 'the end of the file'
          call fail(9, 'table '//trim(tb%name)//': expected "end table" after its '//itoa(tb%np*tb%nT)// &
                    ' data rows, found "'//trim(line)//'"'); return
        end if

      end associate
    end do

    ! nothing after the last table
    call next_line(unitfile, line, nline, ios)
    if (ios == 0) then
      call fail(4, '"tables '//itoa(ntab)//'", but the file goes on after table '//itoa(ntab)//': "'//trim(line)//'"')
      return
    end if
    ios = 0
    close(unitfile)

  contains

    !> Stop the reading with code: print the message with the line of the file, close it and drop the tables.
    subroutine fail ( code, what )
      integer, intent(in)          :: code
      character(len=*), intent(in) :: what
      character(len=24)            :: where
      write(where, '(A,I0)') ' line ', nline
      call report(trim(where)//': '//what)
      close(unitfile)
      if (allocated(tables)) deallocate(tables)
      ios = code
    end subroutine fail

    subroutine report ( what )
      character(len=*), intent(in) :: what
      write(*,'(A)')          '[ERROR] FLINT read_pyrolysis_tables: '//trim(filename)//what_sep(what)//what
      write(error_unit,'(A)') '[ERROR] FLINT read_pyrolysis_tables: '//trim(filename)//what_sep(what)//what
    end subroutine report

    !> The next line (line, its tokens tok and ntok), which must start with keyword; ios /= 0 when it does not.
    subroutine header ( keyword )
      character(len=*), intent(in) :: keyword
      call next_line(unitfile, line, nline, ios)
      if (ios /= 0) then
        call fail(2, 'table '//trim(tables(itab)%name)//': the file ends before "'//keyword//'"'); return
      end if
      call split(line, tok, ntok)
      if (tok(1) /= keyword) then
        call fail(4, 'table '//trim(tables(itab)%name)//': expected "'//keyword//'", found "'//trim(line)//'"')
      end if
    end subroutine header

    !> "<keyword> <n>", then n values over one or more lines: positive and strictly ascending, n >= nmin.
    subroutine grid ( keyword, nmin, v, n )
      character(len=*), intent(in)           :: keyword
      integer, intent(in)                    :: nmin
      real(kind=8), allocatable, intent(out) :: v(:)
      integer, intent(out)                   :: n
      integer :: got, m, st
      n = 0
      call header(keyword)
      if (ios /= 0) return
      if (ntok == 2) read(tok(2), *, iostat=st) n
      if (ntok /= 2 .or. n < nmin) then
        call fail(4, 'table '//trim(tables(itab)%name)//': expected "'//keyword//' <n>", n >= '//itoa(nmin)// &
                  ', found "'//trim(line)//'"'); return
      end if
      allocate(v(n))
      got = 0
      do while (got < n)
        call next_line(unitfile, line, nline, ios)
        if (ios /= 0) then
          call fail(2, 'table '//trim(tables(itab)%name)//': the file ends inside the '//keyword); return
        end if
        m = min(count_tokens(line), n - got)
        read(line, *, iostat=st) v(got+1:got+m)
        if (st /= 0 .or. m < count_tokens(line)) then
          call fail(4, 'table '//trim(tables(itab)%name)//': the '//keyword//' are '//itoa(n)// &
                    ' numbers, found "'//trim(line)//'"'); return
        end if
        got = got + m
      end do
      if (.not. all(v > 0.d0)) then
        call fail(6, 'table '//trim(tables(itab)%name)//': the '//keyword//' must be positive'); return
      end if
      if (n > 1) then
        if (.not. all(v(2:n) > v(1:n-1))) then
          call fail(6, 'table '//trim(tables(itab)%name)//': the '//keyword//' are not strictly ascending'); return
        end if
      end if
    end subroutine grid

  end function read_pyrolysis_tables


  !> The number of tables loaded: 0 before read_pyrolysis_tables, or after it failed.
  integer function pyrolysis_table_count ( )
    implicit none
    pyrolysis_table_count = 0
    if (allocated(tables)) pyrolysis_table_count = size(tables)
  end function pyrolysis_table_count


  !> Name, treatment, condensed column and range of table id (1 <= id <= pyrolysis_table_count()).
  subroutine pyrolysis_table_info ( id, name, treatment, has_condensed, pmin, pmax, Tmin, Tmax )
    implicit none
    integer, intent(in)           :: id
    character(len=*), intent(out) :: name, treatment
    logical, intent(out)          :: has_condensed
    real(kind=8), intent(out)     :: pmin, pmax, Tmin, Tmax

    if (id < 1 .or. id > pyrolysis_table_count()) error stop 'pyrolysis_table_info: no pyrolysis-gas table of that id'
    name          = tables(id)%name
    treatment     = tables(id)%treatment
    has_condensed = tables(id)%has_condensed
    pmin = tables(id)%p(1)
    pmax = tables(id)%p(tables(id)%np)
    Tmin = tables(id)%T(1)
    Tmax = tables(id)%T(tables(id)%nT)

  end subroutine pyrolysis_table_info


  !> Composition of the pyrolysis gas of table id at pressure p [Pa] and temperature T [K], on the species order of the
  !> solver (y of the size of the species list given to read_pyrolysis_tables), and the condensed-carbon mass fraction
  !> (0 for a table without that column). Bilinear in (p, T) on the cell of the grid that holds the point; a point out
  !> of the table takes the value at the nearest point of its boundary (clamped). A table of one pressure is linear
  !> in T.
  pure subroutine pyrolysis_gas_composition ( id, p, T, y, x_condensed )
    implicit none
    integer, intent(in)                 :: id
    real(kind=8), intent(in)            :: p, T
    real(kind=8), intent(out)           :: y(:)
    real(kind=8), intent(out), optional :: x_condensed
    ! Local
    integer      :: i, j
    real(kind=8) :: Tc, pc, wT, wp, xc

    if (.not. allocated(tables)) error stop 'pyrolysis_gas_composition: no pyrolysis-gas table loaded'
    if (id < 1 .or. id > size(tables)) error stop 'pyrolysis_gas_composition: no pyrolysis-gas table of that id'
    if (size(y) /= nsolver) error stop 'pyrolysis_gas_composition: y is not of the size of the species list'

    associate ( tb => tables(id) )
      Tc = min(max(T, tb%T(1)), tb%T(tb%nT))
      i  = cell(Tc, tb%T)
      wT = (Tc - tb%T(i)) / (tb%T(i+1) - tb%T(i))
      if (tb%np == 1) then
        y  = (1.d0 - wT) * tb%y(:,i,1) + wT * tb%y(:,i+1,1)
        xc = (1.d0 - wT) * tb%xc(i,1)  + wT * tb%xc(i+1,1)
      else
        pc = min(max(p, tb%p(1)), tb%p(tb%np))
        j  = cell(pc, tb%p)
        wp = (pc - tb%p(j)) / (tb%p(j+1) - tb%p(j))
        y  = (1.d0 - wp) * ( (1.d0 - wT) * tb%y(:,i,j)   + wT * tb%y(:,i+1,j)   ) &
           +         wp  * ( (1.d0 - wT) * tb%y(:,i,j+1) + wT * tb%y(:,i+1,j+1) )
        xc = (1.d0 - wp) * ( (1.d0 - wT) * tb%xc(i,j)    + wT * tb%xc(i+1,j)    ) &
           +         wp  * ( (1.d0 - wT) * tb%xc(i,j+1)  + wT * tb%xc(i+1,j+1)  )
      end if
    end associate
    if (present(x_condensed)) x_condensed = xc

  end subroutine pyrolysis_gas_composition


  !> The cell of the strictly ascending grid g (size >= 2) that holds x, g(1) <= x <= g(n): the i, 1 <= i <= n-1, with
  !> g(i) <= x <= g(i+1) (bisection; x = g(n) gives n-1, and a NaN gives 1).
  pure integer function cell ( x, g )
    implicit none
    real(kind=8), intent(in) :: x, g(:)
    integer :: lo, hi, mid
    lo = 1
    hi = size(g)
    do while (hi - lo > 1)
      mid = (lo + hi) / 2
      if (x >= g(mid)) then
        lo = mid
      else
        hi = mid
      end if
    end do
    cell = lo
  end function cell


  !> The next line of unit that is not blank once its comment is removed, the comment removed: a '#' that starts the
  !> line or follows a blank starts a comment. nline counts the lines read. ios /= 0 at the end of the file.
  subroutine next_line ( unit, line, nline, ios )
    implicit none
    integer, intent(in)                        :: unit
    character(len=:), allocatable, intent(out) :: line
    integer, intent(inout)                     :: nline
    integer, intent(out)                       :: ios
    integer :: k
    do
      call read_line(unit, line, ios)
      if (ios /= 0) return
      nline = nline + 1
      do k = 1, len(line)
        if (line(k:k) == char(9)) line(k:k) = ' '
      end do
      do k = 1, len(line)
        if (line(k:k) == '#') then
          if (k == 1) then
            line = ''
          else if (line(k-1:k-1) == ' ') then
            line = line(1:k-1)
          else
            cycle
          end if
          exit
        end if
      end do
      if (len_trim(line) > 0) return
    end do
  end subroutine next_line


  !> One line of unit, of any length.
  subroutine read_line ( unit, line, ios )
    implicit none
    integer, intent(in)                        :: unit
    character(len=:), allocatable, intent(out) :: line
    integer, intent(out)                       :: ios
    character(len=1024) :: buf
    integer :: n
    line = ''
    do
      read(unit, '(A)', advance='no', size=n, iostat=ios) buf
      line = line // buf(1:n)
      if (ios /= 0) exit
    end do
    if (is_iostat_eor(ios)) ios = 0
  end subroutine read_line


  !> The first size(tok) blank-separated tokens of line (each cut to the length of tok) and their total number, ntok.
  pure subroutine split ( line, tok, ntok )
    implicit none
    character(len=*), intent(in)  :: line
    character(len=*), intent(out) :: tok(:)
    integer, intent(out)          :: ntok
    integer :: k, b
    tok  = ''
    ntok = 0
    k = 1
    do
      do while (k <= len(line))
        if (line(k:k) /= ' ') exit
        k = k + 1
      end do
      if (k > len(line)) exit
      b = k
      do while (k <= len(line))
        if (line(k:k) == ' ') exit
        k = k + 1
      end do
      ntok = ntok + 1
      if (ntok <= size(tok)) tok(ntok) = line(b:k-1)
    end do
  end subroutine split


  pure integer function count_tokens ( line )
    implicit none
    character(len=*), intent(in) :: line
    character(len=1) :: t(1)
    call split(line, t, count_tokens)
  end function count_tokens


  !> The n numbers of line after its first skip tokens (every second one when pairs: element, mass fraction, ...).
  logical function values_after ( line, skip, n, v, pairs ) result ( ok )
    implicit none
    character(len=*), intent(in) :: line
    integer, intent(in)          :: skip, n
    real(kind=8), intent(out)    :: v(:)
    logical, intent(in)          :: pairs
    character(len=PG_NAME_LEN), allocatable :: t(:)
    integer :: nt, k, st
    allocate(t(skip + n))
    call split(line, t, nt)
    ok = nt == skip + n
    if (.not. ok) return
    do k = 1, size(v)
      if (pairs) then
        read(t(skip + 2*k), *, iostat=st) v(k)
      else
        read(t(skip + k), *, iostat=st) v(k)
      end if
      if (st /= 0) then
        ok = .false.
        return
      end if
    end do
  end function values_after


  pure function itoa ( i ) result ( s )
    implicit none
    integer, intent(in)           :: i
    character(len=:), allocatable :: s
    character(len=16) :: b
    write(b, '(I0)') i
    s = trim(b)
  end function itoa


  pure function rtoa ( x ) result ( s )
    implicit none
    real(kind=8), intent(in)      :: x
    character(len=:), allocatable :: s
    character(len=24) :: b
    write(b, '(ES12.5)') x
    s = trim(adjustl(b))
  end function rtoa


  pure function what_sep ( what ) result ( s )
    implicit none
    character(len=*), intent(in)  :: what
    character(len=:), allocatable :: s
    if (len(what) > 0) then
      if (what(1:1) == ' ') then
        s = ''
        return
      end if
    end if
    s = ': '
  end function what_sep

end module FLINT_Lib_PyrolysisGas
