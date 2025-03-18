using Plots

function run_fdtd_dielectric_waveguide_coupler_gap_shifted()
    # ---------------------------
    # Simulation and Grid Setup
    # ---------------------------
    c = 1.0                         # Speed in vacuum (normalized)
    x_min, x_max = 0.0, 200.0       # Extended domain in x
    y_min, y_max = -40.0, 40.0      # Domain in y

    dx = 0.1
    dy = 0.1
    nx = Int(round((x_max - x_min) / dx)) + 1
    ny = Int(round((y_max - y_min) / dy)) + 1

    x = range(x_min, x_max, length=nx)
    y = range(y_min, y_max, length=ny)

    # ---------------------------
    # CFL Condition: Adaptive dt
    # ---------------------------
    cfl = 0.9999999
    dt = cfl * dx / (c * sqrt(2))
    println("Adaptive dt based on CFL condition: ", dt)

    nsteps = 6000  # Total time steps

    # ---------------------------
    # Define Piecewise Geometry
    # ---------------------------
    # Regions along x:
    #   1)  0 <= x < 20: waveguides straight at ±12  (farther from center)
    #   2) 20 <= x < 50: transition from ±12 to ±6
    #   3) 50 <= x < 90: waveguides remain at ±6     (coupler region, narrow gap)
    #   4) 90 <= x < 120: transition from ±6 back to ±12
    #   5) 120 <= x <=200: waveguides remain at ±12
    function y_upper(x_val)
        if x_val < 20
            return 12.0
        elseif x_val < 50
            # from +12 down to +6 (linear transition)
            return 12.0 - 4.0 * (x_val - 20) / 20.0
        elseif x_val < 90
            return 6.0
        elseif x_val < 120
            # from +6 back up to +12 (linear transition)
            return 6.0 + 4.0 * (x_val - 90) / 20.0
        else
            return 12.0
        end
    end

    function y_lower(x_val)
        if x_val < 20
            return -12.0
        elseif x_val < 50
            # from -12 up to -6
            return -12.0 + 4.0 * (x_val - 20) / 20.0
        elseif x_val < 90
            return -6.0
        elseif x_val < 120
            # from -6 back down to -12
            return -6.0 - 4.0 * (x_val - 90) / 20.0
        else
            return -12.0
        end
    end

    # Waveguide parameters
    core_half_width = 5.0
    core_eps = 3.8
    clad_eps = 1.0

    # ---------------------------
    # Define Dielectric Profile
    # ---------------------------
    eps = fill(clad_eps, nx, ny)
    for i in 1:nx
        for j in 1:ny
            yu = y_upper(x[i])
            yl = y_lower(x[i])
            if (y[j] >= yu - core_half_width) && (y[j] <= yu + core_half_width)
                eps[i, j] = core_eps
            elseif (y[j] >= yl - core_half_width) && (y[j] <= yl + core_half_width)
                eps[i, j] = core_eps
            end
        end
    end

    # ---------------------------
    # Field Initialization
    # ---------------------------
    Ez_old = zeros(Float64, nx, ny)
    Ez     = zeros(Float64, nx, ny)
    Ez_new = zeros(Float64, nx, ny)

    # Precompute finite-difference coefficients
    Cx = (c * dt / dx)^2
    Cy = (c * dt / dy)^2

    # ---------------------------
    # Source: Gaussian Beam Injection at x_min
    # ---------------------------
    omega = 2π * 0.1
    sigma = 2.0  # Gaussian beam width

    # ---------------------------
    # Prepare Animation
    # ---------------------------
    anim = Animation()

    # ---------------------------
    # Simulation Loop
    # ---------------------------
    for n in 1:nsteps
        t = n * dt

        # Update interior grid points using FDTD scheme
        @inbounds for i in 2:nx-1
            for j in 2:ny-1
                Ez_new[i, j] = 2 * Ez[i, j] - Ez_old[i, j] +
                               (Cx/eps[i, j]) * (Ez[i+1, j] - 2 * Ez[i, j] + Ez[i-1, j]) +
                               (Cy/eps[i, j]) * (Ez[i, j+1] - 2 * Ez[i, j] + Ez[i, j-1])
            end
        end

        # Left boundary (Mur + source)
        alpha = (c * dt - dx) / (c * dt + dx)
        for j in 2:ny-1
            # Inject a Gaussian beam centered at y_upper(x_min)
            Ez_new[1, j] = sin(omega * t) * exp(-((y[j] - y_upper(x_min))^2) / (2*sigma^2))
        end

        # Right boundary (Mur)
        for j in 2:ny-1
            Ez_new[nx, j] = Ez[nx-1, j] + alpha * (Ez_new[nx-1, j] - Ez[nx, j])
        end

        # Top/Bottom boundaries (Mur)
        beta = (c * dt - dy) / (c * dt + dy)
        for i in 2:nx-1
            Ez_new[i, 1]  = Ez[i, 2] + beta * (Ez_new[i, 2] - Ez[i, 1])
            Ez_new[i, ny] = Ez[i, ny-1] + beta * (Ez_new[i, ny-1] - Ez[i, ny])
        end

        # Advance time
        Ez_old .= Ez
        Ez .= Ez_new

        # Save frame every 5 steps
        if mod(n, 5) == 0
            p = heatmap(
                x, y, Ez',
                aspect_ratio = 1,
                clims = (-0.1, 0.1),
                color = :RdBu,
                xlabel = "x",
                ylabel = "y",
                title = "Coupler Waveguides Simulation \n with Gaussian Beam Source (t = $(round(t, digits=2)))",
                xlims = (x_min, x_max),
                ylims = (y_min, y_max)
            )

            # Plot waveguide boundaries
            upper_top = [y_upper(xv) + core_half_width for xv in x]
            upper_bot = [y_upper(xv) - core_half_width for xv in x]
            lower_top = [y_lower(xv) + core_half_width for xv in x]
            lower_bot = [y_lower(xv) - core_half_width for xv in x]

            plot!(p, x, upper_top, lw=2, color=:red, label="Upper Core Boundaries")
            plot!(p, x, upper_bot, lw=2, color=:red, label="")
            plot!(p, x, lower_top, lw=2, color=:red, label="Lower Core Boundaries")
            plot!(p, x, lower_bot, lw=2, color=:red, label="")

            frame(anim, p)
        end
    end

    gif(anim, "fdtd_directional_coupler_gaussian_beam.gif", fps=30)
end

# Run the simulation
run_fdtd_dielectric_waveguide_coupler_gap_shifted()
