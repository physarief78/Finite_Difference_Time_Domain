    using Plots

    function run_fdtd_bending_waveguide_simulation()
        # ---------------------------
        # Simulation and Grid Setup
        # ---------------------------
        # Use the real speed of light (in m/s)
        c = 3.0e8                        # Wave speed (m/s)
        # Domain: x from 0 to 100 m; y from -10 to 10 m (wider domain to accommodate the bend)
        x_min, x_max = 0.0, 100.0
        y_min, y_max = -50.0, 50.0

        # Grid spacing (in meters)
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
        # With physical c, the time step becomes very small.
        cfl = 0.999999999999999  # CFL safety factor (≤ 1.0)
        dt = cfl * dx / (c * sqrt(2))
        println("Adaptive dt based on CFL condition: ", dt, " s")

        nsteps = 2000  # Total number of time steps

        # ---------------------------
        # Field Initialization
        # ---------------------------
        Ez_old = zeros(Float64, nx, ny)  # Field at time step n-1
        Ez     = zeros(Float64, nx, ny)  # Field at time step n
        Ez_new = zeros(Float64, nx, ny)  # Field at time step n+1

        # Finite-difference coefficients (for second derivatives):
        Cx = (c * dt / dx)^2
        Cy = (c * dt / dy)^2

        # ---------------------------
        # Define Bending Waveguide Geometry
        # ---------------------------
        # Define the channel half–width and a bending centerline.
        channel_halfwidth = 10.0  # in meters
        center_line(x_val) = 5.0 * sin(2π * x_val / 50)  # centerline oscillates with amplitude 2 m

        # ---------------------------
        # Source Parameters
        # ---------------------------
        # Place the source near x = 5 m, at the channel center.
        source_x = 5.0
        source_center = center_line(source_x)
        i_source = Int(round((source_x - x_min) / dx)) + 1
        j_source = Int(round((source_center - y_min) / dy)) + 1

        # Set the source frequency to 10 MHz (so that the wavelength ~ 30 m)
        # (You may adjust this along with grid resolution in a real simulation.)
        omega = 2π * 1.0e7  # Angular frequency in rad/s

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
            source_val = sin(omega * t)
            Ez_new[i_source, j_source] += source_val

            # ---------------------------
            # Enforce PEC Conditions Outside the Bending Waveguide
            # ---------------------------
            # For each x, compute the local channel boundaries and zero the field outside.
            for i in 1:nx
                c_val = center_line(x[i])
                bottom_bound = c_val - channel_halfwidth
                top_bound = c_val + channel_halfwidth
                for j in 1:ny
                    if y[j] < bottom_bound || y[j] > top_bound
                        Ez_new[i, j] = 0.0
                    end
                end
            end

            # ---------------------------
            # Mur Boundary Condition for x (left and right edges)
            # ---------------------------
            alpha = (c * dt - dx) / (c * dt + dx)
            for j in 2:ny-1
                # Left boundary (i = 1)
                Ez_new[1, j] = Ez[2, j] + alpha * (Ez_new[2, j] - Ez[1, j])
                # Right boundary (i = nx)
                Ez_new[nx, j] = Ez[nx-1, j] + alpha * (Ez_new[nx-1, j] - Ez[nx, j])
            end

            # ---------------------------
            # Plot the Field
            # ---------------------------
            # Transpose the field so that rows correspond to y.
            p = heatmap(x, y, Ez', aspect_ratio=1,
                clims=(-0.1, 0.1), color=:RdBu,
                xlabel="x (m)", ylabel="y (m)",
                title="2D FDTD Bending Waveguide Simulation, t = $(round(t, sigdigits=3)) s",
                xlims=(x_min, x_max), ylims=(y_min, y_max))
            # Overlay the bending waveguide walls
            bottom_curve = [center_line(x_val) - channel_halfwidth for x_val in x]
            top_curve = [center_line(x_val) + channel_halfwidth for x_val in x]
            plot!(p, x, bottom_curve, lw=2, color=:red, label="Wall")
            plot!(p, x, top_curve, lw=2, color=:red, label="")

            # ---------------------------
            # Advance to Next Time Step
            # ---------------------------
            Ez_old .= Ez
            Ez     .= Ez_new

            p  # Return the plot for the current frame.
        end

        # Save the animation as a GIF (adjust fps as desired)
        gif(anim, "fdtd_bending_waveguide.gif", fps=30)
    end

    # Run the simulation
    run_fdtd_bending_waveguide_simulation()
