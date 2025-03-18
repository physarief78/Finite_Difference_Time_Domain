using Plots

"""
Smooth transition function S(x) that goes from 0 to 1 as x goes from x_start to x_end.
Here we use a half-cosine shape for smoothness, but you can replace it with a simpler
linear interpolation if you prefer.
"""
function smooth_transition(x, x_start, x_end)
    if x <= x_start
        return 0.0
    elseif x >= x_end
        return 1.0
    else
        # Half-cosine interpolation from 0 to 1
        ξ = (x - x_start) / (x_end - x_start)
        return 0.5 * (1.0 - cos(π * ξ))
    end
end

"""
Compute the center of the top waveguide as x moves from x_start to x_end,
going from +5 (at x_start) to +15 (at x_end).
"""
function center_top(x, x_start, x_end)
    S = smooth_transition(x, x_start, x_end)
    return 5.0 + 10.0 * S   # transitions from +5 to +15
end

"""
Compute the center of the bottom waveguide as x moves from x_start to x_end,
going from -5 (at x_start) to -15 (at x_end).
"""
function center_bot(x, x_start, x_end)
    S = smooth_transition(x, x_start, x_end)
    return -5.0 - 10.0 * S  # transitions from -5 to -15
end

function run_fdtd_y_splitter_sbend_simulation()
    # ---------------------------
    # Simulation and Grid Setup
    # ---------------------------
    c = 1.0
    x_min, x_max = 0.0, 100.0       # Expanded x-range
    y_min, y_max = -30.0, 30.0      # Expanded y-range

    dx = 0.1
    dy = 0.1
    nx = Int(round((x_max - x_min) / dx)) + 1
    ny = Int(round((y_max - y_min) / dy)) + 1

    x = range(x_min, x_max, length=nx)
    y = range(y_min, y_max, length=ny)

    # ---------------------------
    # CFL Condition & Time Steps
    # ---------------------------
    cfl = 0.9999999
    dt = cfl * dx / (c * sqrt(2))
    println("Adaptive dt based on CFL condition: ", dt)

    nsteps = 3000

    # ---------------------------
    # Define S-bend Waveguide Geometry
    # ---------------------------
    # Four regions:
    # 1) For x < x_transition_start (x < 15): single waveguide at y in [-5,5]
    # 2) For x in [x_transition_start, x_start) (15 ≤ x < 30): transition region where
    #    the waveguide edges linearly interpolate from [-5,5] at x=15 to [-10,10] at x=30.
    # 3) For x in [x_start, x_end] (30 ≤ x ≤ 70): S-bend region with two waveguides whose centers
    #    smoothly shift from ±5 to ±15.
    # 4) For x > x_end (x > 70): two parallel waveguides at y in [10,20] and [-20,-10].
    x_transition_start = 15.0
    x_start = 30.0
    x_end   = 70.0
    half_width = 5.0     # Each waveguide is 10 wide

    core_eps = 2.56
    eps = ones(Float64, nx, ny)

    # Fill eps in a piecewise manner
    for i in 1:nx
        xi = x[i]
        for j in 1:ny
            yj = y[j]

            if xi < x_transition_start
                # Region 1: Single waveguide from y = -5 to y = +5
                if (yj >= -5.0) && (yj <= 5.0)
                    eps[i, j] = core_eps
                end

            elseif xi < x_start
                # Region 2: Transition region (x from 15 to 30)
                # Linearly interpolate boundaries:
                # Top edge: from 5.0 at x=15 to 10.0 at x=30
                # Bottom edge: from -5.0 at x=15 to -10.0 at x=30
                top_boundary = 5.0 + (xi - x_transition_start) * (10.0 - 5.0) / (x_start - x_transition_start)
                bottom_boundary = -5.0 - (xi - x_transition_start) * (10.0 - 5.0) / (x_start - x_transition_start)
                if (yj >= bottom_boundary) && (yj <= top_boundary)
                    eps[i, j] = core_eps
                end

            elseif xi <= x_end
                # Region 3: S-bend region with smoothly shifting centers
                y_top = center_top(xi, x_start, x_end)
                y_bot = center_bot(xi, x_start, x_end)
                # Top waveguide region
                if (yj >= (y_top - half_width)) && (yj <= (y_top + half_width))
                    eps[i, j] = core_eps
                end
                # Bottom waveguide region
                if (yj >= (y_bot - half_width)) && (yj <= (y_bot + half_width))
                    eps[i, j] = core_eps
                end

            else
                # Region 4: x > x_end, two parallel waveguides:
                # Top waveguide from y=10 to y=20 and bottom from y=-20 to y=-10.
                if (yj >= 10) && (yj <= 20)
                    eps[i, j] = core_eps
                elseif (yj >= -20) && (yj <= -10)
                    eps[i, j] = core_eps
                end
            end
        end
    end

    # ---------------------------
    # Field Initialization
    # ---------------------------
    Ez_old = zeros(Float64, nx, ny)
    Ez     = zeros(Float64, nx, ny)
    Ez_new = zeros(Float64, nx, ny)

    Cx = (c * dt / dx)^2
    Cy = (c * dt / dy)^2

    # ---------------------------
    # Source Parameters
    # ---------------------------
    omega = 2π * 0.1

    # ---------------------------
    # Prepare Animation
    # ---------------------------
    anim = Animation()

    # ---------------------------
    # Simulation Loop
    # ---------------------------
    for n in 1:nsteps
        t = n * dt

        # Update interior points
        @inbounds for i in 2:nx-1
            for j in 2:ny-1
                Ez_new[i, j] = 2 * Ez[i, j] - Ez_old[i, j] +
                    (Cx / eps[i, j]) * (Ez[i+1, j] - 2*Ez[i, j] + Ez[i-1, j]) +
                    (Cy / eps[i, j]) * (Ez[i, j+1] - 2*Ez[i, j] + Ez[i, j-1])
            end
        end

        # Left boundary condition (inject into the single waveguide region)
        alpha = (c * dt - dx) / (c * dt + dx)
        for j in 2:ny-1
            # Injection is done in the region where the waveguide is originally defined
            if abs(y[j]) <= 5
                Ez_new[1, j] = sin(omega * t) * cos((π/2) * (y[j] / 5))
            else
                Ez_new[1, j] = Ez[2, j] + alpha * (Ez_new[2, j] - Ez[1, j])
            end
        end

        # Right boundary (Mur)
        for j in 2:ny-1
            Ez_new[nx, j] = Ez[nx-1, j] + alpha * (Ez_new[nx-1, j] - Ez[nx, j])
        end

        # Top/Bottom Mur boundaries
        beta = (c * dt - dy) / (c * dt + dy)
        for i in 2:nx-1
            # Bottom
            Ez_new[i, 1] = Ez[i, 2] + beta * (Ez_new[i, 2] - Ez[i, 1])
            # Top
            Ez_new[i, ny] = Ez[i, ny-1] + beta * (Ez_new[i, ny-1] - Ez[i, ny])
        end

        # Advance time
        Ez_old .= Ez
        Ez .= Ez_new

        # Plot every few steps
        if mod(n, 3) == 0
            p = heatmap(x, y, Ez', aspect_ratio=1,
                        clims=(-0.75, 0.75), color=:RdBu,
                        xlabel="x", ylabel="y",
                        title="2D FDTD Y-Splitter S-bend (t = $(round(t, digits=2)))",
                        xlims=(x_min, x_max), ylims=(y_min, y_max))

            # --------------------------------------------------
            # Overlay the waveguide boundaries for visualization
            # --------------------------------------------------
            # Region 1: x in [0,15]
            plot!(p, [0, 15], [5.0, 5.0], lw=2, color=:red, label="")
            plot!(p, [0, 15], [-5.0, -5.0], lw=2, color=:red, label="Waveguide")
            # Region 2: x in [15,30] transition
            plot!(p, [15, 30], [5.0, 10.0], lw=2, color=:red, label="")
            plot!(p, [15, 30], [-5.0, -10.0], lw=2, color=:red, label="")
            # Region 3: x in [30,70] S-bend transition
            Nx_bend = 50  # number of points for the bend
            xvals = range(x_start, x_end, length=Nx_bend)
            top_plus = [center_top(xx, x_start, x_end) + half_width for xx in xvals]
            top_minus = [center_top(xx, x_start, x_end) - half_width for xx in xvals]
            bot_plus = [center_bot(xx, x_start, x_end) + half_width for xx in xvals]
            bot_minus = [center_bot(xx, x_start, x_end) - half_width for xx in xvals]

            plot!(p, xvals, top_plus, lw=2, color=:red, label="")
            plot!(p, xvals, top_minus, lw=2, color=:red, label="")
            plot!(p, xvals, bot_plus, lw=2, color=:red, label="")
            plot!(p, xvals, bot_minus, lw=2, color=:red, label="")
            # Region 4: x in [70,100]
            plot!(p, [x_end, x_max], [20, 20], lw=2, color=:red, label="")
            plot!(p, [x_end, x_max], [10, 10], lw=2, color=:red, label="")
            plot!(p, [x_end, x_max], [-10, -10], lw=2, color=:red, label="")
            plot!(p, [x_end, x_max], [-20, -20], lw=2, color=:red, label="")

            frame(anim, p)
        end
    end

    # Save the animation
    gif(anim, "fdtd_y_splitter_sbend.gif", fps=30)
end

# Run the Y-splitter simulation
run_fdtd_y_splitter_sbend_simulation()
