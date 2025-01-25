using Plots

# Parameters
nx, ny = 200, 200                # Grid size
dx, dy = 1.0, 1.0                # Grid spacing
dt = 0.5                         # Time step (must satisfy CFL condition)
nt = 1000                        # Number of time steps

# Constants
c = 1.0                          # Speed of light in vacuum
eps0 = 1.0                       # Permittivity of free space
mu0 = 1.0                        # Permeability of free space

# Fields
Ez = zeros(nx, ny)               # z-component of electric field
Hx = zeros(nx, ny-1)             # x-component of magnetic field (staggered in y)
Hy = zeros(nx-1, ny)             # y-component of magnetic field (staggered in x)

# Previous time-step fields (for Mur's ABC)
Ez_prev = zeros(nx, ny)

# Update coefficients
cezh = dt / (eps0 * dx)
chxe = dt / (mu0 * dy)
chye = dt / (mu0 * dx)

# Source parameters
source_position_x = 2           # Source position along the x-axis
frequency = 0.025                # Frequency of the wave
wavelength = c / frequency       # Wavelength
k = 2 * π / wavelength           # Wave number

# Slits parameters
slit_center_x = nx ÷ 2          # Center of the slits in the x-direction
slit_width = 5                  # Width of each slit
slit_gap = 30                   # Gap between the two slits

# Create a mask for the slits
slit_mask = ones(nx, ny)        # Default: allow propagation everywhere
slit_mask[slit_center_x, :] .= 0.0  # Block the entire center row initially

# Define the slits
slit_start1 = ny ÷ 2 - slit_gap ÷ 2 - slit_width
slit_end1 = slit_start1 + slit_width
slit_start2 = ny ÷ 2 + slit_gap ÷ 2
slit_end2 = slit_start2 + slit_width

# Open up the slits
slit_mask[slit_center_x, slit_start1:slit_end1] .= 1.0
slit_mask[slit_center_x, slit_start2:slit_end2] .= 1.0

println("Simulation started...")
start_time = time()  # Start timing

# Animation loop
anim = @animate for t in 1:nt
    # Save the previous electric field for Mur's ABC
    Ez_prev[:, :] .= Ez[:, :]

    # Update magnetic fields (staggered in y and x)
    Hx[:, :] .= Hx[:, :] - chye .* (Ez[:, 2:end] - Ez[:, 1:end-1])
    Hy[:, :] .= Hy[:, :] + chxe .* (Ez[2:end, :] - Ez[1:end-1, :])

    # Update electric field (centered)
    Ez[2:end-1, 2:end-1] .= Ez[2:end-1, 2:end-1] + cezh .* (
        (Hy[2:end, 2:end-1] - Hy[1:end-1, 2:end-1]) -
        (Hx[2:end-1, 2:end] - Hx[2:end-1, 1:end-1])
    )

    # Apply the slit mask to the electric field
    Ez[slit_center_x, :] .= Ez[slit_center_x, :] .* slit_mask[slit_center_x, :]

    # Add sinusoidal wave source
    Ez[source_position_x, :] .= sin(2 * π * frequency * t - k * source_position_x)

    # Absorbing boundary conditions (Mur's first-order ABC)
    Ez[1, 2:end-1] .= Ez_prev[2, 2:end-1] + (c * dt - dx) / (c * dt + dx) * (Ez[2, 2:end-1] - Ez_prev[1, 2:end-1])
    Ez[end, 2:end-1] .= Ez_prev[end-1, 2:end-1] + (c * dt - dx) / (c * dt + dx) * (Ez[end-1, 2:end-1] - Ez_prev[end, 2:end-1])
    Ez[2:end-1, 1] .= Ez_prev[2:end-1, 2] + (c * dt - dy) / (c * dt + dy) * (Ez[2:end-1, 2] - Ez_prev[2:end-1, 1])
    Ez[2:end-1, end] .= Ez_prev[2:end-1, end-1] + (c * dt - dy) / (c * dt + dy) * (Ez[2:end-1, end-1] - Ez_prev[2:end-1, end])

    Ez_visual = copy(Ez)
    Ez_visual[slit_mask .== 0] .= NaN

    # Plot the field
    heatmap(Ez_visual', color=cgrad([:white, :black, :gold]), title="2D EM Wave with Double Slits", clims=(-0.4, 0.4), nan_color=:white)
end

# Save the animation and stop timing
gif(anim, "2d_fdtd_em_wave_double_slits.gif", fps=30)

end_time = time()  # End timing
execution_time = end_time - start_time
println("Simulation finished.")
println("Total execution time: $(round(execution_time, digits=2)) seconds.")
