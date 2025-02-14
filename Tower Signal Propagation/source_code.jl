using Plots

# ---------------------------
# Simulation parameters (SI units)
# ---------------------------
nx = 500                # Number of grid points in x
ny = 500                # Number of grid points in y
nt = 1000               # Number of time steps (increase this if needed)

dx = 1.0                # Spatial step in x (meters)
dy = 1.0                # Spatial step in y (meters)

# ---------------------------
# Physical constants (SI units)
# ---------------------------
ϵ0 = 8.854187817e-12    # Permittivity of free space (F/m)
μ0 = 4π * 1e-7          # Permeability of free space (H/m)
c = 1/sqrt(ϵ0*μ0)       # Speed of light (~2.9979e8 m/s)

# ---------------------------
# Adaptive CFL condition (2D)
# ---------------------------
dt_max = 1 / (c * sqrt((1/dx^2) + (1/dy^2)))
cfl = 0.99              # Safety factor (<= 1)
dt = cfl * dt_max       # Time step (seconds)
println("Adaptive time step dt = $(dt) s  (", dt*1e9, " ns)")

# ---------------------------
# Field arrays (initialized to zero)
# ---------------------------
Ez = zeros(Float64, nx, ny)  # Out-of-plane electric field (V/m)
Hx = zeros(Float64, nx, ny)  # Magnetic field (x–component)
Hy = zeros(Float64, nx, ny)  # Magnetic field (y–component)

# ---------------------------
# Dipole source parameters
# ---------------------------
# Place the dipole at the center of the grid.
src_i = div(nx+1, 2)
src_j = div(ny+1, 2)

# Choose a source frequency (Hz) appropriate for the simulation.
f = 1e7                   # Frequency in Hz (1 MHz)
ω = 2π * f                # Angular frequency (rad/s)

# ---------------------------
# Define coordinate arrays so that the source is at (0,0)
# ---------------------------
xcoords = collect(1:nx) .- src_i  # x coordinates (meters)
ycoords = collect(1:ny) .- src_j  # y coordinates (meters)

# ---------------------------
# Environment: fixed point & slopes overlay
# ---------------------------
# Sloping lines (environment overlay):
x_left  = range(-20, stop=0, length=500)
y_left  = range(-250, stop=0, length=500)
x_right = range(0, stop=20, length=500)
y_right = range(0, stop=-250, length=500)

# ---------------------------
# Animation of the simulation with Mur boundaries
# ---------------------------
anim = @animate for n in 1:nt
    t = n * dt          # simulation time in seconds
    t_ns = t * 1e9      # simulation time in nanoseconds

    # Save a copy of the old Ez for Mur boundary calculations.
    Ez_old = copy(Ez)
    
    # --- Update Magnetic Fields ---
    for i in 1:nx
        for j in 1:ny-1
            Hx[i,j] -= (dt/(μ0 * dy)) * (Ez[i,j+1] - Ez[i,j])
        end
    end

    for i in 1:nx-1
        for j in 1:ny
            Hy[i,j] += (dt/(μ0 * dx)) * (Ez[i+1,j] - Ez[i,j])
        end
    end

    # --- Update Electric Field at interior points (avoid boundaries) ---
    for i in 2:nx-1
        for j in 2:ny-1
            Ez[i,j] += (dt/ϵ0) * ( (Hy[i,j] - Hy[i-1,j])/dx - (Hx[i,j] - Hx[i,j-1])/dy )
        end
    end

    # --- Inject the Dipole Source (pure sinusoidal) ---
    src_val = sin(ω * t)
    Ez[src_i,   src_j]   += src_val
    Ez[src_i, src_j+1]   -= src_val

    # --- Apply Mur Absorbing Boundary Conditions ---
    mur_x = (c*dt - dx) / (c*dt + dx)
    mur_y = (c*dt - dy) / (c*dt + dy)
    
    # Left & Right boundaries (for j = 2:ny-1)
    for j in 2:ny-1
        Ez[1,j]  = Ez_old[2,j] + mur_x*(Ez[2,j] - Ez_old[1,j])
        Ez[nx,j] = Ez_old[nx-1,j] + mur_x*(Ez[nx-1,j] - Ez_old[nx,j])
    end

    # Bottom & Top boundaries (for i = 2:nx-1)
    for i in 2:nx-1
        Ez[i,1]  = Ez_old[i,2] + mur_y*(Ez[i,2] - Ez_old[i,1])
        Ez[i,ny] = Ez_old[i,ny-1] + mur_y*(Ez[i,ny-1] - Ez_old[i,ny])
    end

    # Corner boundaries (average of adjacent boundaries)
    Ez[1,1]   = 0.5*(Ez[2,1]   + Ez[1,2])
    Ez[nx,1]  = 0.5*(Ez[nx-1,1]  + Ez[nx,2])
    Ez[1,ny]  = 0.5*(Ez[2,ny]  + Ez[1,ny-1])
    Ez[nx,ny] = 0.5*(Ez[nx-1,ny] + Ez[nx,ny-1])

    # --- Plot the Current Ez Field with Heatmap and Environment ---
    heatmap(xcoords, ycoords, Ez,
        color = :RdBu, clims = (-0.0025, 0.0025), aspect_ratio = 1,
        xlims = (minimum(xcoords), maximum(xcoords)),
        ylims = (minimum(ycoords), maximum(ycoords)),
        title = "Tower Signal Propagation: Ez at t = $(round(t_ns, digits=2)) ns",
        xlabel = "x (m)", ylabel = "y (m)", legend = true)
    
    scatter!([0], [0], markersize = 8, label = false)
    plot!(x_left, y_left, label = false, lw = 2)
    plot!(x_right, y_right, label = false, lw = 2)
end

# Save the animation as a GIF.
gif(anim, "tower_signal_propagation.mp4", fps = 30)
