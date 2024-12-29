module interpolate_mod
implicit none
public

! Public subroutines and functions
public :: min_value, interpolation, find_y, read_data_file

contains

! Function to calculate the minimum dx value
double precision function min_value(x, n)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)
    double precision :: min_x(n-1)
    integer :: ii

    do ii = 1, n-1
        min_x(ii) = abs(x(ii+1) - x(ii))
    end do
    min_value = minval(min_x)
end function min_value

! Subroutine for interpolation
subroutine interpolation(x, y, n, x_final, y_final, n_final)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: x(n), y(n)
    double precision, allocatable, intent(out) :: x_final(:), y_final(:)
    integer, intent(out) :: n_final
    double precision :: mx, dx
    integer :: ii, idx
    double precision, allocatable :: x_new(:), y_new(:)

    mx = min_value(x, n)
    dx = mx / 100.d0

    allocate(x_final(0))
    allocate(y_final(0))

    do ii = 1, n-1
        allocate(x_new(ceiling((x(ii+1) - x(ii)) / dx)))
        allocate(y_new(size(x_new)))

        x_new = [(x(ii) + (idx - 1) * dx, idx = 1, size(x_new))]

        if (ii == n-1) then
            y_new = [(y(ii) + (idx - 1) * (y(ii+1) - y(ii)) / (size(x_new) - 1), idx = 1, size(x_new))]
        else
            y_new = [(y(ii) + (idx - 1) * (y(ii+1) - y(ii)) / size(x_new), idx = 1, size(x_new))]
        end if

        x_final = [x_final, x_new]
        y_final = [y_final, y_new]

        deallocate(x_new, y_new)
    end do

    n_final = size(x_final)
end subroutine interpolation

! Function to find y value for a given radius
function find_y(x, r_new, em_new, n) result(y_new)
    implicit none
    double precision, intent(in) :: x
    double precision, intent(in) :: r_new(n), em_new(n)
    integer, intent(in) :: n
    double precision :: y_new, dx
    integer :: idx

    dx = 0.0005d0
    idx = maxloc((r_new >= x - dx) .and. (r_new <= x + dx), dim=1)
    if (size(idx) > 0) then
        y_new = sum(em_new(idx)) / size(idx)
    else
        y_new = 0.d0
    end if
end function find_y

! Subroutine to read data from a file
subroutine read_data_file(filename, x, y, n_points)
    implicit none
    character(len=*), intent(in) :: filename
    double precision, allocatable, intent(out) :: x(:), y(:)
    integer, intent(out) :: n_points
    integer :: unit_num, ios, i
    double precision :: temp_x, temp_y
    character(len=256) :: line

    unit_num = 20  ! File unit number
    open(unit=unit_num, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, "Error: Unable to open file", trim(filename)
        stop
    end if

    ! Count lines to determine the number of data points
    n_points = 0
    do
        read(unit_num, '(A)', iostat=ios) line
        if (ios /= 0) exit
        n_points = n_points + 1
    end do
    rewind(unit_num)

    ! Allocate arrays based on the number of lines
    allocate(x(n_points), y(n_points))

    ! Read data into arrays
    i = 1
    do
        read(unit_num, *, iostat=ios) temp_x, temp_y
        if (ios /= 0) exit
        x(i) = temp_x
        y(i) = temp_y
        i = i + 1
    end do

    close(unit_num)
end subroutine read_data_file

end module interpolate_mod

