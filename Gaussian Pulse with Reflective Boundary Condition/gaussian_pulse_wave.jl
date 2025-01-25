using Plots

# Parameters
nx, ny = 200, 200                # Grid size
dx, dy = 1.0, 1.0                # Grid spacing
dt = 0.5                         # Time step (must satisfy CFL condition)
nt = 1000                         # Number of time steps

# Constants
c = 1.0                          # Speed of light in vacuum
eps0 = 1.0                       # Permittivity of free space
mu0 = 1.0                        # Permeability of free space

# Fields
Ez = zeros(nx, ny)               # z-component of electric field
Hx = zeros(nx, ny-1)             # x-component of magnetic field (staggered in y)
Hy = zeros(nx-1, ny)             # y-component of magnetic field (staggered in x)

# Update coefficients
cezh = dt / (eps0 * dx)
chxe = dt / (mu0 * dy)
chye = dt / (mu0 * dx)

# Source parameters
source_position = (nx ÷ 2, ny ÷ 2)

# Visualization setup
anim = @animate for t in 1:nt
    # Update magnetic fields (staggered in y and x)
    Hx[:, :] .= Hx[:, :] - chye .* (Ez[:, 2:end] - Ez[:, 1:end-1])
    Hy[:, :] .= Hy[:, :] + chxe .* (Ez[2:end, :] - Ez[1:end-1, :])

    # Update electric field (centered)
    Ez[2:end-1, 2:end-1] .= Ez[2:end-1, 2:end-1] + cezh .* (
        (Hy[2:end, 2:end-1] - Hy[1:end-1, 2:end-1]) -
        (Hx[2:end-1, 2:end] - Hx[2:end-1, 1:end-1])
    )

    # Add source
    Ez[source_position...] += exp(-((t - 30) / 10)^2)

    # Plot the field
    heatmap(Ez', color=cgrad([:blue, :white, :red]), title="2D Electromagnetic Wave", clims=(-0.05, 0.05))
end

# Save or display the animation
gif(anim, "gaussian_pulse_source.gif", fps=30)
