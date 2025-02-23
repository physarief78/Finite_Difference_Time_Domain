using Plots

function run_fdtd_waveguide_simulation()
    # ---------------------------
    # Simulation and Grid Setup
    # ---------------------------
    c = 1.0                         # Wave speed (normalized)
    # Domain: x from 0 to 50; y between the waveguide walls: -2 and 2
    x_min, x_max = 0.0, 100.0
    y_min, y_max = -5.0, 5.0

    # Grid spacing (note: dx ≠ dy here)
    dx = 0.1
    dy = 0.1
    nx = Int(round((x_max - x_min) / dx)) + 1
    ny = Int(round((y_max - y_min) / dy)) + 1

    # Create coordinate arrays for plotting
    x = range(x_min, x_max, length=nx)
    y = range(y_min, y_max, length=ny)

    # ---------------------------
    # CFL Condition: Adaptive dt
    # ---------------------------
    # For simplicity we use dx as the limiting length:
    cfl = 0.9999999  # CFL safety factor (≤ 1.0)
    dt = cfl * dx / (c * sqrt(2))
    println("Adaptive dt based on CFL condition: ", dt)

    nsteps = 1500  # Total number of time steps for the animation

    # ---------------------------
    # Field Initialization
    # ---------------------------
    # We simulate the 2nd-order-in-time scalar wave equation:
    #    ∂²Ez/∂t² = c² (∂²Ez/∂x² + ∂²Ez/∂y²)
    # Using a finite-difference scheme, we require three time-level arrays.
    Ez_old = zeros(Float64, nx, ny)  # Field at time step n-1
    Ez     = zeros(Float64, nx, ny)  # Field at time step n
    Ez_new = zeros(Float64, nx, ny)  # Field at time step n+1

    # Since dx ≠ dy, compute separate finite-difference coefficients:
    Cx = (c * dt / dx)^2
    Cy = (c * dt / dy)^2

    # ---------------------------
    # Source Parameters
    # ---------------------------
    # Inject a sinusoidal wave source near x = 5, centered in y.
    i_source = Int(round((5.0 - x_min) / dx)) + 1
    j_source = Int(round((0.0 - y_min) / dy)) + 1

    # For the sinusoidal source we only need an angular frequency.
    omega = 2π * 0.1    # Angular frequency

    # ---------------------------
    # Animation Loop
    # ---------------------------
    anim = @animate for n = 1:nsteps
        t = n * dt

        # Update interior grid points (finite-difference update)
        @inbounds for i in 2:nx-1
            for j in 2:ny-1
                Ez_new[i, j] = 2 * Ez[i, j] - Ez_old[i, j] +
                    Cx * (Ez[i+1, j] - 2 * Ez[i, j] + Ez[i-1, j]) +
                    Cy * (Ez[i, j+1] - 2 * Ez[i, j] + Ez[i, j-1])
            end
        end

        # Add the sinusoidal source term at the designated grid point.
        source = sin(omega * t)
        Ez_new[i_source, j_source] += source

        # ---------------------------
        # Boundary Conditions for y (waveguide walls)
        # ---------------------------
        # Enforce PEC conditions: the field is zero at y = y_min and y = y_max.
        for i in 1:nx
            Ez_new[i, 1]  = 0.0    # Bottom wall (y = y_min)
            Ez_new[i, ny] = 0.0    # Top wall (y = y_max)
        end

        # ---------------------------
        # Mur Boundary Condition for x (left and right edges)
        # ---------------------------
        # Compute the Mur parameter
        alpha = (c * dt - dx) / (c * dt + dx)
        # Apply Mur BC for j = 2:ny-1 (avoid overriding the PEC boundaries at y=1 and y=ny)
        for j in 2:ny-1
            # Left boundary (i = 1)
            Ez_new[1, j] = Ez[2, j] + alpha * (Ez_new[2, j] - Ez[1, j])
            # Right boundary (i = nx)
            Ez_new[nx, j] = Ez[nx-1, j] + alpha * (Ez_new[nx-1, j] - Ez[nx, j])
        end

        # ---------------------------
        # Plot the Field
        # ---------------------------
        # (Note: transpose the field so that rows correspond to y.)
        p = heatmap(x, y, Ez', aspect_ratio=1,
            clims=(-0.1, 0.1), color=:RdBu,
            xlabel="x", ylabel="y",
            title="2D FDTD Waveguide Simulation, t = $(round(t, digits=2))",
            xlims=(x_min, x_max), ylims=(-30, 30))
        # Overlay the waveguide walls (red lines) at y = -2 and y = 2.
        plot!(p, [x_min, x_max], [y_min, y_min], lw=2, color=:red, label="Wall")
        plot!(p, [x_min, x_max], [y_max, y_max], lw=2, color=:red, label="")

        # ---------------------------
        # Advance to Next Time Step
        # ---------------------------
        # Update the previous time levels for the next iteration.
        Ez_old .= Ez
        Ez     .= Ez_new

        p  # Return the plot for the current frame.
    end

    # Save the animation as a GIF (adjust fps as desired)
    gif(anim, "fdtd_waveguide.gif", fps=30)
end

# Run the simulation
run_fdtd_waveguide_simulation()
