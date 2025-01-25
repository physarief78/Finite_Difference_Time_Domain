program double_slit_fdtd
  implicit none

  ! Parameters
  integer, parameter :: nx = 200, ny = 200
  real(kind=8), parameter :: dx = 1.0d0, dy = 1.0d0
  real(kind=8), parameter :: dt = 0.5d0
  integer, parameter :: nt = 1000

  ! Constants
  real(kind=8), parameter :: c = 1.0d0
  real(kind=8), parameter :: eps0 = 1.0d0
  real(kind=8), parameter :: mu0 = 1.0d0

  ! Derived parameters
  integer, parameter :: source_position_x = 2
  real(kind=8), parameter :: frequency = 0.025d0
  real(kind=8), parameter :: wavelength = c / frequency
  real(kind=8), parameter :: k = 2.0d0 * 3.141592653589793d0 / wavelength

  ! Slit parameters
  integer, parameter :: slit_center_x = nx / 2, slit_width = 5, slit_gap = 30
  integer :: slit_start1, slit_end1, slit_start2, slit_end2

  ! Field arrays
  real(kind=8) :: Ez(nx, ny), Hx(nx, ny-1), Hy(nx-1, ny)
  real(kind=8) :: Ez_prev(nx, ny), slit_mask(nx, ny)

  ! Update coefficients
  real(kind=8), parameter :: cezh = dt / (eps0 * dx)
  real(kind=8), parameter :: chxe = dt / (mu0 * dy)
  real(kind=8), parameter :: chye = dt / (mu0 * dx)

  ! Variables
  integer :: i, j, t
  real(kind=8) :: start_time, end_time, elapsed_time

  ! Initialize fields
  Ez = 0.0d0
  Hx = 0.0d0
  Hy = 0.0d0
  Ez_prev = 0.0d0
  slit_mask = 1.0d0

  ! Define slits
  slit_mask(slit_center_x, :) = 0.0d0
  slit_start1 = ny / 2 - slit_gap / 2 - slit_width
  slit_end1 = slit_start1 + slit_width - 1
  slit_start2 = ny / 2 + slit_gap / 2
  slit_end2 = slit_start2 + slit_width - 1

  slit_mask(slit_center_x, slit_start1:slit_end1) = 1.0d0
  slit_mask(slit_center_x, slit_start2:slit_end2) = 1.0d0

  ! Open file for writing all time steps
  open(unit=10, file='fdtd_output.dat', status='unknown')

  ! Write the slit mask to the file
  write(10, '("Slit Mask:")')
  do i = 1, nx
    write(10, *) (slit_mask(i, j), j = 1, ny)
  end do
  write(10, *) ! Empty line for separation

  ! Record the start time
  call cpu_time(start_time)

  ! Time-stepping loop
  do t = 1, nt

    ! Save previous electric field
    Ez_prev = Ez

    ! Update magnetic fields
    do i = 1, nx
      do j = 1, ny-1
        Hx(i, j) = Hx(i, j) - chye * (Ez(i, j+1) - Ez(i, j))
      end do
    end do

    do i = 1, nx-1
      do j = 1, ny
        Hy(i, j) = Hy(i, j) + chxe * (Ez(i+1, j) - Ez(i, j))
      end do
    end do

    ! Update electric field
    do i = 2, nx-1
      do j = 2, ny-1
        Ez(i, j) = Ez(i, j) + cezh * ((Hy(i, j) - Hy(i-1, j)) - (Hx(i, j) - Hx(i, j-1)))
      end do
    end do

    ! Apply the slit mask
    Ez(slit_center_x, :) = Ez(slit_center_x, :) * slit_mask(slit_center_x, :)

    ! Add sinusoidal wave source
    do j = 1, ny
      Ez(source_position_x, j) = sin(2.0d0 * 3.141592653589793d0 * frequency * t - k * source_position_x)
    end do

    ! Absorbing boundary conditions (Mur's first-order ABC)
    do j = 2, ny-1
      Ez(1, j) = Ez_prev(2, j) + (c * dt - dx) / (c * dt + dx) * (Ez(2, j) - Ez_prev(1, j))
      Ez(nx, j) = Ez_prev(nx-1, j) + (c * dt - dx) / (c * dt + dx) * (Ez(nx-1, j) - Ez_prev(nx, j))
    end do

    do i = 2, nx-1
      Ez(i, 1) = Ez_prev(i, 2) + (c * dt - dy) / (c * dt + dy) * (Ez(i, 2) - Ez_prev(i, 1))
      Ez(i, ny) = Ez_prev(i, ny-1) + (c * dt - dy) / (c * dt + dy) * (Ez(i, ny-1) - Ez_prev(i, ny))
    end do

    ! Write time step data to file
    write(10, '("Time step: ", i5)') t
    do i = 1, nx
      write(10, *) (Ez(i, j), j = 1, ny)
    end do

  end do

  ! Record the end time
  call cpu_time(end_time)

  ! Calculate elapsed time
  elapsed_time = end_time - start_time
  print *, "Execution Time (seconds): ", elapsed_time

  ! Close the file
  close(10)

end program double_slit_fdtd
