program fdtd_bending_waveguide_simulation_optimized
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  !-----------------------------------------------------------------
  ! Parameters and variable declarations.
  !-----------------------------------------------------------------
  integer, parameter :: output_interval = 1  ! Frame output frequency.
  integer :: nx, ny, nsteps, i, j, n, i_source, j_source
  real(dp) :: c, dx, dy, dt, t, cfl, Cx, Cy, source_val, omega
  real(dp) :: x_min, x_max, y_min, y_max
  real(dp) :: channel_halfwidth
  real(dp), allocatable, target :: field_old(:,:), field(:,:), field_new(:,:)
  real(dp), pointer :: Ez_old(:,:), Ez(:,:), Ez_new(:,:)
  real(dp), allocatable :: x(:), y(:)
  real(dp) :: pi, alpha
  real(dp) :: start_time, end_time
  integer :: file_unit
  real(dp), allocatable :: center(:), bottom_bound(:), top_bound(:)
  real(dp), pointer :: tmp(:,:)  ! Temporary pointer for swapping

  !-----------------------------------------------------------------
  ! Define physical constants and grid parameters.
  !-----------------------------------------------------------------
  pi = acos(-1.0_dp)
  c = 3.0e8_dp           ! Speed of light in m/s

  ! Domain: x from 0 to 100 m; y from -50 to 50 m.
  x_min = 0.0_dp
  x_max = 100.0_dp
  y_min = -50.0_dp
  y_max = 50.0_dp

  ! Grid spacing and grid dimensions.
  dx = 0.1_dp
  dy = 0.1_dp
  nx = int((x_max - x_min) / dx + 0.5_dp) + 1
  ny = int((y_max - y_min) / dy + 0.5_dp) + 1

  ! Total number of time steps.
  nsteps = 2000

  ! CFL condition with safety factor.
  cfl = 0.999999999999999_dp
  dt = cfl * dx / (c * sqrt(2.0_dp))
  print *, "Adaptive dt based on CFL condition: ", dt, " s"

  !-----------------------------------------------------------------
  ! Allocate field arrays and associate pointers for pointer swapping.
  ! The target attribute is needed for pointer association.
  !-----------------------------------------------------------------
  allocate(field_old(nx, ny), field(nx, ny), field_new(nx, ny))
  Ez_old => field_old
  Ez     => field
  Ez_new => field_new

  ! Allocate coordinate arrays.
  allocate(x(nx), y(ny))
  do i = 1, nx
     x(i) = x_min + (i - 1) * dx
  end do
  do j = 1, ny
     y(j) = y_min + (j - 1) * dy
  end do

  ! Initialize all fields to zero.
  Ez_old = 0.0_dp
  Ez     = 0.0_dp
  Ez_new = 0.0_dp

  ! Finite-difference coefficients.
  Cx = (c * dt / dx)**2
  Cy = (c * dt / dy)**2

  !-----------------------------------------------------------------
  ! Define the bending waveguide geometry.
  !-----------------------------------------------------------------
  channel_halfwidth = 10.0_dp  ! Half–width of the channel (in meters).

  ! Precompute the channel center and PEC boundaries for each x.
  allocate(center(nx), bottom_bound(nx), top_bound(nx))
  do i = 1, nx
     center(i)      = 5.0_dp * sin(2.0_dp * pi * x(i) / 50.0_dp)
     bottom_bound(i) = center(i) - channel_halfwidth
     top_bound(i)    = center(i) + channel_halfwidth
  end do

  !-----------------------------------------------------------------
  ! Set up the source.
  ! Place the source near x = 5 m at the local channel center.
  !-----------------------------------------------------------------
  i_source = int((5.0_dp - x_min) / dx + 0.5_dp) + 1
  j_source = int((center(i_source) - y_min) / dy + 0.5_dp) + 1

  ! Source frequency: 10 MHz (angular frequency).
  omega = 2.0_dp * pi * 1.0e7_dp

  !-----------------------------------------------------------------
  ! Open a file to save all frame data.
  !-----------------------------------------------------------------
  open(newunit = file_unit, file = "bending_waveguide_fdtd.dat", status = "replace", &
       action = "write", form = "formatted")

  !-----------------------------------------------------------------
  ! Start timing the simulation.
  !-----------------------------------------------------------------
  call CPU_TIME(start_time)

  !-----------------------------------------------------------------
  ! Main time-stepping loop.
  !-----------------------------------------------------------------
  do n = 1, nsteps
     t = n * dt

     !--------------------------------------------------------------
     ! Finite-difference update for interior points.
     ! The loops are re-ordered (j outer, i inner) to improve memory
     ! access (Fortran arrays are column-major) and are parallelized.
     !--------------------------------------------------------------
!$omp parallel do collapse(2) default(shared) private(i,j)
     do j = 2, ny - 1
        do i = 2, nx - 1
           Ez_new(i, j) = 2.0_dp * Ez(i, j) - Ez_old(i, j) + &
                           Cx * (Ez(i + 1, j) - 2.0_dp * Ez(i, j) + Ez(i - 1, j)) + &
                           Cy * (Ez(i, j + 1) - 2.0_dp * Ez(i, j) + Ez(i, j - 1))
        end do
     end do
!$omp end parallel do

     !--------------------------------------------------------------
     ! Add the sinusoidal source term.
     !--------------------------------------------------------------
     source_val = sin(omega * t)
     Ez_new(i_source, j_source) = Ez_new(i_source, j_source) + source_val

     !--------------------------------------------------------------
     ! Enforce PEC (perfect electric conductor) conditions outside
     ! the bending waveguide.
     ! Since the channel boundaries (bottom_bound and top_bound) were
     ! precomputed, we simply test each point.
     !--------------------------------------------------------------
!$omp parallel do default(shared) private(j)
     do i = 1, nx
        do j = 1, ny
           if (y(j) < bottom_bound(i) .or. y(j) > top_bound(i)) then
              Ez_new(i, j) = 0.0_dp
           end if
        end do
     end do
!$omp end parallel do

     !--------------------------------------------------------------
     ! Apply Mur absorbing boundary conditions on the left and right.
     !--------------------------------------------------------------
!$omp parallel do default(shared) private(j)
     alpha = (c * dt - dx) / (c * dt + dx)
     do j = 2, ny - 1
        Ez_new(1, j)  = Ez(2, j) + alpha * (Ez_new(2, j) - Ez(1, j))
        Ez_new(nx, j) = Ez(nx - 1, j) + alpha * (Ez_new(nx - 1, j) - Ez(nx, j))
     end do
!$omp end parallel do

     !--------------------------------------------------------------
     ! Save the current frame to file.
     ! (Here we write every output_interval time step.)
     !--------------------------------------------------------------
     if (mod(n, output_interval) == 0) then
        write(file_unit, *) "Frame ", n, " Time ", t, " NX ", nx, " NY ", ny
        do j = 1, ny
           write(file_unit, *) (Ez_new(i, j), i = 1, nx)
        end do
        write(file_unit, *)  ! Blank line to separate frames.
     end if

     !--------------------------------------------------------------
     ! Swap pointers (rotate the roles of old, current, and new arrays)
     ! instead of copying arrays. This avoids a full-array copy.
     !--------------------------------------------------------------
     tmp     => Ez_old
     Ez_old  => Ez
     Ez      => Ez_new
     Ez_new  => tmp
  end do

  !-----------------------------------------------------------------
  ! End timing the simulation.
  !-----------------------------------------------------------------
  call CPU_TIME(end_time)
  close(file_unit)
  
  print *, "Simulation completed. All frames saved in bending_waveguide_fdtd.dat."
  print *, "Elapsed CPU time: ", end_time - start_time, " seconds."
  
end program fdtd_bending_waveguide_simulation_optimized
