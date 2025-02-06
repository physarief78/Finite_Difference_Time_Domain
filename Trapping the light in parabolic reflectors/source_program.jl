using Plots

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
    nsteps = 2000                # Total number of time steps

    # ----------------------------
    # 4. Allocate field arrays
    # ----------------------------
    Ez = zeros(Nx, Ny)          # z-directed electric field
    Hx = zeros(Nx, Ny)          # x-directed magnetic field
    Hy = zeros(Nx, Ny)          # y-directed magnetic field

    # ----------------------------
    # 5. Define the parabolic reflector geometries (PEC boundaries)
    # ----------------------------
    # Left reflector parameters
    f_left = 10.0             # Focal length for left reflector
    # Equation: x = y^2/(4*f_left) = y^2/40, with vertex at (0,0) and focus at (10,0)
    reflector_mask_left = falses(Nx, Ny)

    # Right reflector parameters
    f_right = 10.0            # Focal length for right reflector
    vertex_right = 65.0       # Vertex of the right reflector is at (80,0)
    # Equation for a parabola opening left: x = vertex_right - y^2/(4*f_right) = 80 - y^2/80.
    # Its focus is at (vertex_right - f_right, 0) = (60,0).
    reflector_mask_right = falses(Nx, Ny)

    # Loop over grid points in y to flag PEC cells near the reflector curves.
    for j in 1:Ny
        y_val = ygrid[j]
        # Left reflector curve: x = y^2/40.
        x_left = (y_val^2) / (4*f_left)
        for i in 1:Nx
            # Mark cells within one grid spacing dx of the left reflector curve.
            if xgrid[i] < x_left && xgrid[i] > (x_left - dx)
                reflector_mask_left[i,j] = true
            end
        end

        # Right reflector curve: x = 80 - y^2/80.
        x_right = vertex_right - (y_val^2) / (4*f_right)
        for i in 1:Nx
            # Mark cells within one dx of the right reflector curve.
            if xgrid[i] > x_right && xgrid[i] < (x_right + dx)
                reflector_mask_right[i,j] = true
            end
        end
    end

    # Combine the two reflector masks into one overall PEC mask.
    reflector_mask = reflector_mask_left .| reflector_mask_right

    # ----------------------------
    # 6. Define the source location 
    # ----------------------------
    # We inject the Gaussian pulse at the left reflector’s focus.
    focal_x = f_left         # (10,0)
    focal_y = 0.0

    # Find nearest grid indices for the source location.
    source_i = findmin(abs.(xgrid .- focal_x))[2]
    source_j = findmin(abs.(ygrid .- focal_y))[2]

    # ----------------------------
    # 7. FDTD Coefficients
    # ----------------------------
    ch = dt/μ0
    ce = dt/ε0

    # ----------------------------
    # 8. FDTD Simulation Loop with Animation
    # ----------------------------
    anim = @animate for n in 1:nsteps

        # --- Update Magnetic Fields ---
        # Update Hx (avoid j=1 because of j-1 index)
        for i in 1:Nx
            for j in 2:Ny
                Hx[i,j] -= ch * (Ez[i,j] - Ez[i,j-1]) / dy
            end
        end
        # Update Hy (avoid i=1 because of i-1 index)
        for i in 2:Nx
            for j in 1:Ny
                Hy[i,j] += ch * (Ez[i,j] - Ez[i-1,j]) / dx
            end
        end

        # --- Update Electric Field ---
        # Avoid the boundaries (i=1 and i=Nx, j=1 and j=Ny) for finite differences.
        for i in 2:(Nx-1)
            for j in 2:(Ny-1)
                Ez[i,j] += ce * ((Hy[i+1,j] - Hy[i,j])/dx - (Hx[i,j+1] - Hx[i,j])/dy)
            end
        end

        # --- Add Gaussian Pulse Source at the left focus ---
        t0 = 20.0             # Delay for the Gaussian pulse peak
        spread = 6.0          # Width of the Gaussian pulse
        Ez[source_i, source_j] += exp(-((n - t0)^2) / (spread^2))

        # --- Impose PEC (Perfect Electric Conductor) on the reflectors ---
        for i in 1:Nx
            for j in 1:Ny
                if reflector_mask[i,j]
                    Ez[i,j] = 0.0
                end
            end
        end

        # --- Simple Absorbing Boundary Conditions ---
        # Force Ez=0 at the edges to reduce reflections.
        Ez[1, :] .= 0.0
        Ez[end, :] .= 0.0
        Ez[:, 1] .= 0.0
        Ez[:, end] .= 0.0

        # --- Plot the instantaneous field ---
        p = contourf(xgrid, ygrid, Ez', levels=50, c=:balance,
                     xlabel="x (um)", ylabel="y (um)",
                     title="FDTD Light Wave of Ez Simulation at (Step = $n)",
                     aspect_ratio=1)
        # Overlay the left reflector curve.
        y_curve = collect(ymin:dy:ymax)
        x_curve_left = [(y^2)/(4*f_left) for y in y_curve]
        plot!(p, x_curve_left, y_curve, lw=3, color=:black, label="Reflector")
        # Overlay the right reflector curve.
        x_curve_right = [vertex_right - (y^2)/(4*f_right) for y in y_curve]
        plot!(p, x_curve_right, y_curve, lw=3, color=:black, label=nothing)
        # Mark the source location.
        scatter!(p, [focal_x], [focal_y], color=:red, label="Source")
        p  # Return the plot for the animation frame
    end

    # Save the animation as a GIF file
    gif(anim, "trapping_the_light_wave_animation.gif", fps=30)
    println("Animation saved as fdtd_tez_animation.gif")
end

# Run the simulation (runs silently and only saves the animation)
run_fdtd_simulation()
