program iterative_solvers
  use iso_fortran_env, only : real64, error_unit
  use iterative_methods
  implicit none

  character(len=:), allocatable :: input_path, output_path
  integer :: n, max_iter, ios, input_unit, output_unit
  real(real64) :: tol, omega
  real(real64), allocatable :: a(:, :), b(:), x0(:)

  call resolve_paths(input_path, output_path)

  call open_file_read(input_path, input_unit)
  call read_problem(input_unit, n, max_iter, tol, a, b, x0, omega)
  close(input_unit)

  open(newunit=output_unit, file=output_path, status='replace', action='write', iostat=ios)
  if (ios /= 0) then
    write(error_unit, '(A)') 'Could not open output file: '//trim(output_path)
    error stop 1
  end if

  call jacobi_method(a, b, x0, max_iter, tol, output_unit)
  call gauss_seidel_method(a, b, x0, max_iter, tol, output_unit)
  call sor_method(a, b, x0, max_iter, tol, omega, output_unit)

  close(output_unit)
contains
  subroutine resolve_paths(input_path, output_path)
    character(len=:), allocatable, intent(out) :: input_path, output_path
    call get_argument(1, input_path, 'data/input.txt')
    call get_argument(2, output_path, 'data/output_fortran.txt')
  end subroutine resolve_paths

  subroutine get_argument(index, value, default)
    integer, intent(in) :: index
    character(len=:), allocatable, intent(out) :: value
    character(len=*), intent(in) :: default
    integer :: len_arg

    call get_command_argument(index, length=len_arg)
    if (len_arg > 0) then
      allocate(character(len=len_arg) :: value)
      call get_command_argument(index, value)
    else
      allocate(character(len=len(default)) :: value)
      value = default
    end if
  end subroutine get_argument

  subroutine open_file_read(path, unit_number)
    character(len=*), intent(in) :: path
    integer, intent(out) :: unit_number
    integer :: ios

    open(newunit=unit_number, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(error_unit, '(A)') 'Could not open input file: '//trim(path)
      error stop 1
    end if
  end subroutine open_file_read

  subroutine read_problem(unit_number, n, max_iter, tol, a, b, x0, omega)
    integer, intent(in) :: unit_number
    integer, intent(out) :: n, max_iter
    real(real64), intent(out) :: tol, omega
    real(real64), allocatable, intent(out) :: a(:, :), b(:), x0(:)

    integer :: i

    call read_integer(unit_number, n)
    call read_integer(unit_number, max_iter)
    call read_real(unit_number, tol)

    allocate(a(n, n), b(n), x0(n))

    do i = 1, n
      call read_real_array(unit_number, a(i, :))
    end do
    call read_real_array(unit_number, b)
    call read_real_array(unit_number, x0)
    call read_real(unit_number, omega)
  end subroutine read_problem

  subroutine read_integer(unit_number, value)
    integer, intent(in) :: unit_number
    integer, intent(out) :: value
    character(len=:), allocatable :: line
    integer :: ios

    call read_nonempty_line(unit_number, line, ios)
    if (ios /= 0) error stop 'Unexpected end of file while reading integer.'
    read(line, *, iostat=ios) value
    if (ios /= 0) error stop 'Invalid integer in input file.'
  end subroutine read_integer

  subroutine read_real(unit_number, value)
    integer, intent(in) :: unit_number
    real(real64), intent(out) :: value
    character(len=:), allocatable :: line
    integer :: ios

    call read_nonempty_line(unit_number, line, ios)
    if (ios /= 0) error stop 'Unexpected end of file while reading real value.'
    read(line, *, iostat=ios) value
    if (ios /= 0) error stop 'Invalid real number in input file.'
  end subroutine read_real

  subroutine read_real_array(unit_number, values)
    integer, intent(in) :: unit_number
    real(real64), intent(out) :: values(:)
    character(len=:), allocatable :: line
    integer :: ios

    call read_nonempty_line(unit_number, line, ios)
    if (ios /= 0) error stop 'Unexpected end of file while reading array.'
    read(line, *, iostat=ios) values
    if (ios /= 0) error stop 'Invalid array values in input file.'
  end subroutine read_real_array

  subroutine read_nonempty_line(unit_number, line, ios)
    integer, intent(in) :: unit_number
    character(len=:), allocatable, intent(out) :: line
    integer, intent(out) :: ios
    character(len=512) :: buffer

    do
      read(unit_number, '(A)', iostat=ios) buffer
      if (ios /= 0) then
        line = ''
        return
      end if

      if (len_trim(buffer) == 0) cycle

      line = trim(buffer)
      return
    end do
  end subroutine read_nonempty_line
end program iterative_solvers
