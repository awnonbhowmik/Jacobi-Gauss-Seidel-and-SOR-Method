module iterative_methods
  use iso_fortran_env, only : real64, error_unit
  implicit none
  real(real64), parameter :: eps = 1.0e-15_real64
contains
  pure function relative_error(new_values, old_values) result(err)
    real(real64), intent(in) :: new_values(:), old_values(:)
    real(real64) :: err, denom
    integer :: i

    err = 0.0_real64
    do i = 1, size(new_values)
      denom = max(abs(new_values(i)), eps)
      err = max(err, abs(new_values(i) - old_values(i)) / denom)
    end do
  end function relative_error

  subroutine write_iteration(unit_out, iter, values)
    integer, intent(in) :: unit_out, iter
    real(real64), intent(in) :: values(:)

    write(unit_out, '(2X,I4,1X,*(1X,F12.6))') iter, values
  end subroutine write_iteration

  subroutine jacobi_method(a, b, x_start, max_iter, tol, unit_out)
    real(real64), intent(in) :: a(:, :), b(:), x_start(:)
    integer, intent(in) :: max_iter, unit_out
    real(real64), intent(in) :: tol

    integer :: n, i, j, iter
    real(real64) :: sum_val, err
    real(real64), allocatable :: x(:), x_old(:)

    n = size(b)
    call validate_dimensions(a, n)
    call validate_iteration_input(max_iter, tol)

    allocate(x(n), x_old(n))
    x_old = x_start

    write(unit_out, '(A)') 'Jacobi Method Result'
    do iter = 1, max_iter
      do i = 1, n
        sum_val = 0.0_real64
        do j = 1, n
          if (j == i) cycle
          sum_val = sum_val + a(i, j) * x_old(j)
        end do
        call ensure_nonzero_diagonal(a(i, i), 'Jacobi', i)
        x(i) = (b(i) - sum_val) / a(i, i)
      end do

      call write_iteration(unit_out, iter, x)

      err = relative_error(x, x_old)
      if (err < tol) exit
      x_old = x
    end do

    write(unit_out, '(A)') ''
  end subroutine jacobi_method

  subroutine gauss_seidel_method(a, b, x_start, max_iter, tol, unit_out)
    real(real64), intent(in) :: a(:, :), b(:), x_start(:)
    integer, intent(in) :: max_iter, unit_out
    real(real64), intent(in) :: tol

    integer :: n, i, j, iter
    real(real64) :: sum_val, err
    real(real64), allocatable :: x(:), x_old(:)

    n = size(b)
    call validate_dimensions(a, n)
    call validate_iteration_input(max_iter, tol)

    allocate(x(n), x_old(n))
    x = x_start

    write(unit_out, '(A)') 'Gauss-Seidel Method Result'
    do iter = 1, max_iter
      x_old = x
      do i = 1, n
        sum_val = 0.0_real64
        do j = 1, n
          if (j < i) then
            sum_val = sum_val + a(i, j) * x(j)
          else if (j > i) then
            sum_val = sum_val + a(i, j) * x_old(j)
          end if
        end do
        call ensure_nonzero_diagonal(a(i, i), 'Gauss-Seidel', i)
        x(i) = (b(i) - sum_val) / a(i, i)
      end do

      call write_iteration(unit_out, iter, x)

      err = relative_error(x, x_old)
      if (err < tol) exit
    end do

    write(unit_out, '(A)') ''
  end subroutine gauss_seidel_method

  subroutine sor_method(a, b, x_start, max_iter, tol, omega, unit_out)
    real(real64), intent(in) :: a(:, :), b(:), x_start(:)
    integer, intent(in) :: max_iter, unit_out
    real(real64), intent(in) :: tol, omega

    integer :: n, i, j, iter
    real(real64) :: sum_val, err
    real(real64), allocatable :: x(:), x_old(:)

    if (omega <= 0.0_real64 .or. omega >= 2.0_real64) then
      error stop 'SOR relaxation factor (w) must be in the interval (0, 2).'
    end if

    n = size(b)
    call validate_dimensions(a, n)
    call validate_iteration_input(max_iter, tol)

    allocate(x(n), x_old(n))
    x = x_start

    write(unit_out, '(A)') 'Successive Over-Relaxation (SOR) Result'
    do iter = 1, max_iter
      x_old = x
      do i = 1, n
        sum_val = 0.0_real64
        do j = 1, n
          if (j < i) then
            sum_val = sum_val + a(i, j) * x(j)
          else if (j > i) then
            sum_val = sum_val + a(i, j) * x_old(j)
          end if
        end do
        call ensure_nonzero_diagonal(a(i, i), 'SOR', i)
        x(i) = (1.0_real64 - omega) * x_old(i) + omega * (b(i) - sum_val) / a(i, i)
      end do

      call write_iteration(unit_out, iter, x)

      err = relative_error(x, x_old)
      if (err < tol) exit
    end do

    write(unit_out, '(A)') ''
  end subroutine sor_method

  subroutine validate_dimensions(a, n)
    real(real64), intent(in) :: a(:, :)
    integer, intent(in) :: n
    if (size(a, 1) /= n .or. size(a, 2) /= n) then
      error stop 'Matrix A must be square and match the length of vector b.'
    end if
  end subroutine validate_dimensions

  subroutine validate_iteration_input(max_iter, tol)
    integer, intent(in) :: max_iter
    real(real64), intent(in) :: tol
    if (max_iter <= 0) error stop 'Maximum iterations must be positive.'
    if (tol <= 0.0_real64) error stop 'Tolerance must be positive.'
  end subroutine validate_iteration_input

  subroutine ensure_nonzero_diagonal(value, method, row)
    real(real64), intent(in) :: value
    character(len=*), intent(in) :: method
    integer, intent(in) :: row
    if (abs(value) < eps) then
      write(error_unit, '(A,1X,I0)') 'Zero diagonal encountered in '//trim(method)//' at row', row
      error stop 'Diagonal entries must be non-zero.'
    end if
  end subroutine ensure_nonzero_diagonal
end module iterative_methods
