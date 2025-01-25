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
source_position = (nx ÷ 2, ny ÷ 2)  # Set x=1 and y=middle of the grid
frequency = 0.025               # Frequency of the wave

println("Simulation started...")
start_time = time()  # Start timing

# Animation loop
anim = @animate for t in 1:nt
    # Update magnetic fields (staggered in y and x)
    Hx[:, :] .= Hx[:, :] - chye .* (Ez[:, 2:end] - Ez[:, 1:end-1])
    Hy[:, :] .= Hy[:, :] + chxe .* (Ez[2:end, :] - Ez[1:end-1, :])

    # Update electric field (centered)
    Ez[2:end-1, 2:end-1] .= Ez[2:end-1, 2:end-1] + cezh .* (
        (Hy[2:end, 2:end-1] - Hy[1:end-1, 2:end-1]) -
        (Hx[2:end-1, 2:end] - Hx[2:end-1, 1:end-1])
    )

    # Add sinusoidal-cosinusoidal wave source
    Ez[source_position...] += sin(2 * π * frequency * t) + cos(2 * π * frequency * t)

    # Plot the field
    heatmap(Ez', color=cgrad([:red, :white, :blue]), title="2D Electromagnetic Wave", clims=(-0.1, 0.1))
end

# Save the animation and stop timing
gif(anim, "em_wave_sin_cos_source.gif", fps=30)

end_time = time()  # End timing
execution_time = end_time - start_time
println("Simulation finished.")
println("Total execution time: $(round(execution_time, digits=2)) seconds.")
