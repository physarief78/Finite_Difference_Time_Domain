program fdtd_simulation
  implicit none
  ! Define double precision kind.
  integer, parameter :: dp = selected_real_kind(15, 307)
  
  ! Declaration of variables.
  integer :: nsteps, n
  integer :: Nx, Ny
  integer :: i, j
  integer :: source_i, source_j
  integer :: out_unit, ios

  ! Physical constants.
  real(dp), parameter :: c       = 3.0d8               ! Speed of light [m/s]
  real(dp), parameter :: epsilon0 = 8.854d-12           ! Permittivity [F/m]
  real(dp), parameter :: pi      = 3.141592653589793d0
  real(dp), parameter :: mu0     = 4.0d0 * pi * 1.0d-7   ! Permeability [H/m]
  real(dp) :: dt

  ! Domain parameters and grid spacing.
  real(dp) :: xmin, xmax, ymin, ymax, dx, dy

  ! FDTD coefficients.
  real(dp) :: ch, ce

  ! Geometry for the reflector.
  real(dp) :: f_param, reflector_y_min, reflector_y_max, reflector_x_offset

  ! Source pulse parameters.
  real(dp) :: t0, spread, t1, spread2

  ! Arrays for grid, fields, and reflector curve/mask.
  real(dp), allocatable :: xgrid(:), ygrid(:)
  real(dp), allocatable :: Ez(:,:), Ez_old(:,:)
  real(dp), allocatable :: Hx(:,:), Hy(:,:)
  real(dp), allocatable :: y_curve(:), x_curve(:)
  logical, allocatable :: reflector_mask(:,:)

  !--------------------------------------------------------------------
  ! Set up grid and simulation parameters.
  !--------------------------------------------------------------------
  Nx = 400
  Ny = 300
  xmin = -10.0d0
  xmax = 80.0d0
  ymin = -25.0d0
  ymax = 25.0d0
  dx = (xmax - xmin) / (Nx - 1)
  dy = (ymax - ymin) / (Ny - 1)

  dt = 0.99d0 / ( c * sqrt((1.0d0/dx**2) + (1.0d0/dy**2)) )
  nsteps = 1000

  ch = dt / mu0
  ce = dt / epsilon0

  allocate(xgrid(Nx))
  allocate(ygrid(Ny))
  allocate(Ez(Nx,Ny))
  allocate(Ez_old(Nx,Ny))
  allocate(Hx(Nx,Ny))
  allocate(Hy(Nx,Ny))
  allocate(reflector_mask(Nx,Ny))
  
  ! Set up spatial grid.
  do i = 1, Nx
     xgrid(i) = xmin + (i - 1) * dx
  end do
  do j = 1, Ny
     ygrid(j) = ymin + (j - 1) * dy
  end do

  ! Initialize fields to zero.
  Ez     = 0.0d0
  Ez_old = 0.0d0
  Hx     = 0.0d0
  Hy     = 0.0d0

  !--------------------------------------------------------------------
  ! Define the parabolic reflector.
  !--------------------------------------------------------------------
  f_param            = 10.0d0    ! Focal length (focal point will be at (f_param,0))
  reflector_y_min    = -15.0d0   ! y-range for reflector definition
  reflector_y_max    = 15.0d0
  reflector_x_offset = 0.0d0

  ! Allocate arrays to store the reflector curve.
  allocate(y_curve(Ny))
  allocate(x_curve(Ny))
  do j = 1, Ny
     y_curve(j) = reflector_y_min + (j - 1) * (reflector_y_max - reflector_y_min) / (Ny - 1)
     x_curve(j) = parabolic_curve(y_curve(j), f_param, reflector_x_offset)
  end do

  ! Build a Boolean mask to mark PEC cells near the reflector.
  reflector_mask = .false.
  do j = 1, Ny
     do i = 1, Nx
        if ( xgrid(i) < parabolic_curve(ygrid(j), f_param, reflector_x_offset) .and. &
             xgrid(i) > parabolic_curve(ygrid(j), f_param, reflector_x_offset) - dx ) then
           reflector_mask(i,j) = .true.
        end if
     end do
  end do

  !--------------------------------------------------------------------
  ! Define the source at the focal point.
  !--------------------------------------------------------------------
  source_i = 1
  source_j = 1
  do i = 1, Nx
     if ( abs(xgrid(i) - f_param) < abs(xgrid(source_i) - f_param) ) then
        source_i = i
     end if
  end do
  do j = 1, Ny
     if ( abs(ygrid(j) - 0.0d0) < abs(ygrid(source_j) - 0.0d0) ) then
        source_j = j
     end if
  end do

  !--------------------------------------------------------------------
  ! Set up the source pulse parameters.
  !--------------------------------------------------------------------
  t0     = 20.0d0   ! Delay for the first Gaussian pulse peak.
  spread = 6.0d0    ! Width of the first Gaussian pulse.
  t1     = 500.0d0  ! Delay for the second Gaussian pulse peak.
  spread2= 6.0d0    ! Width of the second Gaussian pulse.

  !--------------------------------------------------------------------
  ! Open one output file to save the entire time-series.
  !--------------------------------------------------------------------
  out_unit = 10
  open(out_unit, file="Ez_time_series.dat", status='replace', action='write', iostat=ios)
  if (ios /= 0) then
     print *, "Error opening file Ez_time_series.dat"
     stop
  end if

  ! Write header information.
  write(out_unit,*) "# FDTD Simulation Time Series Data"
  write(out_unit,*) "# Grid dimensions: Nx =", Nx, " Ny =", Ny
  write(out_unit,*) "# x grid: xmin =", xgrid(1), " xmax =", xgrid(Nx)
  write(out_unit,*) "# y grid: ymin =", ygrid(1), " ymax =", ygrid(Ny)
  write(out_unit,*) "# Each block corresponds to one time step (n) and includes the Ez field values."
  write(out_unit,*) "# Format: Step_number, time, then for each row (j) the Ez values for all x (i)."
  write(out_unit,*) ""

  !--------------------------------------------------------------------
  ! Main FDTD Simulation Loop
  !--------------------------------------------------------------------
  do n = 1, nsteps

     ! --- Update Magnetic Fields ---
     do i = 1, Nx
        do j = 2, Ny
           Hx(i,j) = Hx(i,j) - ch * ( Ez(i,j) - Ez(i,j-1) ) / dy
        end do
     end do

     do i = 2, Nx
        do j = 1, Ny
           Hy(i,j) = Hy(i,j) + ch * ( Ez(i,j) - Ez(i-1,j) ) / dx
        end do
     end do

     ! --- Update Electric Field (interior points) ---
     do i = 2, Nx-1
        do j = 2, Ny-1
           Ez(i,j) = Ez(i,j) + ce * ( (Hy(i+1,j) - Hy(i,j))/dx - (Hx(i,j+1) - Hx(i,j))/dy )
        end do
     end do

     ! --- Add Gaussian Pulse Sources at the focal point ---
     Ez(source_i, source_j) = Ez(source_i, source_j) + exp( -((n - t0)**2)/(spread**2) )
     Ez(source_i, source_j) = Ez(source_i, source_j) + exp( -((n - t1)**2)/(spread2**2) )

     ! --- Impose PEC on the reflector cells ---
     do i = 1, Nx
        do j = 1, Ny
           if ( reflector_mask(i,j) ) then
              Ez(i,j) = 0.0d0
           end if
        end do
     end do

     ! --- Mur Absorbing Boundary Conditions ---
     do i = 2, Nx-1
        Ez(i,1)  = Ez_old(i,2)   + ((c*dt - dy)/(c*dt + dy)) * ( Ez(i,2)   - Ez_old(i,1) )
        Ez(i,Ny) = Ez_old(i,Ny-1)+ ((c*dt - dy)/(c*dt + dy)) * ( Ez(i,Ny-1)- Ez_old(i,Ny) )
     end do

     do j = 2, Ny-1
        Ez(1,j)  = Ez_old(2,j)   + ((c*dt - dx)/(c*dt + dx)) * ( Ez(2,j)   - Ez_old(1,j) )
        Ez(Nx,j) = Ez_old(Nx-1,j)+ ((c*dt - dx)/(c*dt + dx)) * ( Ez(Nx-1,j)- Ez_old(Nx,j) )
     end do

     ! --- Update Ez_old for the next time step ---
     Ez_old = Ez

     ! --- Write current time step data to the output file ---
     call write_Ez_step(out_unit, n, dt, Nx, Ny, Ez)

  end do

  close(out_unit)
  print *, "Simulation complete. Data saved in Ez_time_series.dat"

contains

  !--------------------------------------------------------------------
  ! Elemental function for the parabolic reflector curve.
  ! Computes: x = y^2/(4*f) + x_offset.
  !--------------------------------------------------------------------
  elemental function parabolic_curve(y, f, x_offset) result(x)
    implicit none
    real(dp), intent(in) :: y, f, x_offset
    real(dp) :: x
    x = (y**2) / (4.0d0 * f) + x_offset
  end function parabolic_curve

  !--------------------------------------------------------------------
  ! Subroutine to write one time step of Ez data to an open file.
  ! It writes the step number, simulation time, and the Ez array.
  !--------------------------------------------------------------------
  subroutine write_Ez_step(out_unit, step, dt, Nx, Ny, Ez)
    implicit none
    integer, intent(in) :: out_unit, step, Nx, Ny
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: Ez(Nx,Ny)
    integer :: i, j
    real(dp) :: t

    t = step * dt
    write(out_unit,*) "Step:", step, " Time:", t
    do j = 1, Ny
       write(out_unit,*) (Ez(i,j), i = 1, Nx)
    end do
    write(out_unit,*) ""  ! Blank line separator between time steps.
  end subroutine write_Ez_step

end program fdtd_simulation
