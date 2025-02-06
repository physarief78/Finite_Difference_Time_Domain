using Plots

# Define a function that computes the x-position of the parabolic curve
# given a vector (or scalar) of y values, the focal length f, and an optional x_offset.
function parabolic_curve(y, f; x_offset=0.0)
    return (y.^2) ./ (4*f) .+ x_offset
end

function run_fdtd_simulation()
    # ----------------------------
    # 1. Define Physical Constants
    # ----------------------------
    c = 3e8                     # Speed of light [m/s]
    ε0 = 8.854e-12              # Permittivity of free space [F/m]
    μ0 = 4π*1e-7                # Permeability of free space [H/m]

    # ----------------------------
    # 2. Set up the simulation domain
    # ----------------------------
    xmin, xmax = -10.0, 80.0
    ymin, ymax = -25.0, 25.0

    Nx = 400                    # Number of grid points in x
    Ny = 300                    # Number of grid points in y

    dx = (xmax - xmin) / (Nx-1)
    dy = (ymax - ymin) / (Ny-1)

    # Create coordinate arrays
    xgrid = collect(range(xmin, xmax, length=Nx))
    ygrid = collect(range(ymin, ymax, length=Ny))

    # ----------------------------
    # 3. Time stepping parameters
    # ----------------------------
    # Courant stability condition for 2D
    dt = 0.99/(c*sqrt((1/dx^2) + (1/dy^2)))
    nsteps = 1000                # Total number of time steps

    # ----------------------------
    # 4. Allocate field arrays
    # ----------------------------
    Ez = zeros(Nx, Ny)          # z-directed electric field
    Hx = zeros(Nx, Ny)          # x-directed magnetic field
    Hy = zeros(Nx, Ny)          # y-directed magnetic field

    # ----------------------------
    # 5. Define the adjustable parabolic reflector geometry (PEC boundary)
    # ----------------------------
    # Adjust these parameters as desired:
    f_param = 10.0              # Focal length for the parabola
    reflector_y_min = -15       # Minimum y-value over which to define the reflector
    reflector_y_max = 15        # Maximum y-value over which to define the reflector
    reflector_x_offset = 0.0    # Horizontal shift of the parabola (default 0)

    # Create a vector of y values to generate the curve (for plotting later)
    y_curve = collect(range(reflector_y_min, reflector_y_max, length=Ny))
    # Compute the corresponding x values (adjustable via f_param and x_offset)
    x_curve = parabolic_curve(y_curve, f_param; x_offset=reflector_x_offset)

    # Build a Boolean mask to flag PEC cells near the reflector curve.
    reflector_mask = falses(Nx, Ny)
    for j in 1:Ny
        y_val = ygrid[j]
        # For each y value, compute the corresponding x position on the reflector.
        x_ref = parabolic_curve(y_val, f_param; x_offset=reflector_x_offset)
        for i in 1:Nx
            # Mark cells within one grid spacing dx of the reflector curve as PEC.
            if xgrid[i] < x_ref && xgrid[i] > (x_ref - dx)
                reflector_mask[i,j] = true
            end
        end
    end

    # ----------------------------
    # 6. Define the source location (at the focal point)
    # ----------------------------
    # With this definition of the parabola, the focus is at (f_param, 0)
    focal_x = f_param          # Focal point x-coordinate
    focal_y = 0.0              # Focal point y-coordinate

    # Find nearest grid indices for the source location.
    source_i = findmin(abs.(xgrid .- focal_x))[2]
    source_j = findmin(abs.(ygrid .- focal_y))[2]

    # ----------------------------
    # 7. FDTD Coefficients
    # ----------------------------
    ch = dt/μ0
    ce = dt/ε0

    # ----------------------------
    # 8. Initialize storage for Ez from the previous time step (for Mur ABC)
    # ----------------------------
    Ez_old = copy(Ez)

    # ----------------------------
    # 9. FDTD Simulation Loop with Animation
    # ----------------------------
    anim = @animate for n in 1:nsteps

        # --- Update Magnetic Fields ---
        for i in 1:Nx
            for j in 2:Ny
                Hx[i,j] -= ch * (Ez[i,j] - Ez[i,j-1]) / dy
            end
        end
        for i in 2:Nx
            for j in 1:Ny
                Hy[i,j] += ch * (Ez[i,j] - Ez[i-1,j]) / dx
            end
        end

        # --- Update Electric Field (interior points) ---
        for i in 2:(Nx-1)
            for j in 2:(Ny-1)
                Ez[i,j] += ce * ((Hy[i+1,j] - Hy[i,j])/dx - (Hx[i,j+1] - Hx[i,j])/dy)
            end
        end

        # --- Add first Gaussian Pulse Source at the Focal Point ---
        t0 = 20.0             # Delay for the first Gaussian pulse peak
        spread = 6.0          # Width of the first Gaussian pulse
        Ez[source_i, source_j] += exp(-((n - t0)^2) / (spread^2))

        # --- Add second Gaussian Pulse Source (after some delay) ---
        t1 = 500.0            # Delay for the second Gaussian pulse peak
        spread2 = 6.0         # Width of the second Gaussian pulse
        Ez[source_i, source_j] += exp(-((n - t1)^2) / (spread2^2))

        # --- Impose PEC (Perfect Electric Conductor) on the reflector ---
        for i in 1:Nx
            for j in 1:Ny
                if reflector_mask[i,j]
                    Ez[i,j] = 0.0
                end
            end
        end

        # --- Mur Absorbing Boundary Conditions ---
        # Update the top and bottom (y) boundaries:
        for i in 2:(Nx-1)
            # Bottom boundary (j = 1)
            Ez[i,1] = Ez_old[i,2] + ((c*dt - dy)/(c*dt + dy))*(Ez[i,2] - Ez_old[i,1])
            # Top boundary (j = Ny)
            Ez[i,Ny] = Ez_old[i,Ny-1] + ((c*dt - dy)/(c*dt + dy))*(Ez[i,Ny-1] - Ez_old[i,Ny])
        end
        # Update the left and right (x) boundaries:
        for j in 2:(Ny-1)
            # Left boundary (i = 1)
            Ez[1,j] = Ez_old[2,j] + ((c*dt - dx)/(c*dt + dx))*(Ez[2,j] - Ez_old[1,j])
            # Right boundary (i = Nx)
            Ez[Nx,j] = Ez_old[Nx-1,j] + ((c*dt - dx)/(c*dt + dx))*(Ez[Nx-1,j] - Ez_old[Nx,j])
        end

        # --- Update Ez_old for the next time step ---
        Ez_old = copy(Ez)

        # --- Plot the instantaneous field ---
        p = contourf(xgrid, ygrid, Ez', levels=50, c=:balance,
                     xlabel="x (m)", ylabel="y (m)",
                     title="Radio Wave transmit simulation at (Step = $n)",
                     aspect_ratio=1)
        # Overlay the adjustable parabolic reflector curve.
        plot!(p, x_curve, y_curve, lw=3, color=:black, label="Reflector")
        # Add a straight line at y=0 from x=0 to x=-3 without legend entry.
        plot!(p, [-3, 0], [0, 0], lw=3, color=:black, label=nothing)
        # Add two sloped lines approaching ymin without legend entry.
        plot!(p, [-6, -3], [ymin, 0], lw=3, color=:black, label=nothing)
        plot!(p, [0, -3], [ymin, 0], lw=3, color=:black, label=nothing)
        # Mark the source location.
        scatter!(p, [focal_x], [focal_y], color=:red, label="Source")
        p  # Return the plot to be captured by the animation macro
    end

    # Save the animation as a GIF file
    gif(anim, "fdtd_transmitted_light_wave_animation.gif", fps=30)
    println("Animation saved as fdtd_transmitted_light_wave_animation.gif")
end

# Run the simulation (runs silently and only saves the animation)
run_fdtd_simulation()
