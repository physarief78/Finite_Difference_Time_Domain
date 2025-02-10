using Plots

# ------------------------------------------------------------
# 1. Helper Function: Distance from a Point to a Line Segment
# ------------------------------------------------------------
function distance_to_segment(x, y, x1, y1, x2, y2)
    dx = x2 - x1
    dy = y2 - y1
    if dx == 0 && dy == 0
        return sqrt((x - x1)^2 + (y - y1)^2)
    end
    t = ((x - x1)*dx + (y - y1)*dy) / (dx^2 + dy^2)
    if t < 0
        return sqrt((x - x1)^2 + (y - y1)^2)
    elseif t > 1
        return sqrt((x - x2)^2 + (y - y2)^2)
    else
        proj_x = x1 + t * dx
        proj_y = y1 + t * dy
        return sqrt((x - proj_x)^2 + (y - proj_y)^2)
    end
end

# ------------------------------------------------------------
# 2. Centerline Distance Function (Radius Method) for Downward Turn
# ------------------------------------------------------------
function centerline_distance(x, y)
    # Horizontal segment: from (0,0) to (40,0)
    d1 = distance_to_segment(x, y, 0.0, 0.0, 40.0, 0.0)
    
    # Circular arc: quarter–circle from (40,0) to (60,-20)
    # New center is (40,-20) with radius = 20.
    # For a point on the arc, relative coordinates are: (x - 40, y - (-20)) = (x - 40, y + 20)
    dx_arc = x - 40.0
    dy_arc = y + 20.0
    d_center = sqrt(dx_arc^2 + dy_arc^2)
    # For points on the arc, the angle (in radians) should lie between 0 and π/2.
    angle = atan(dy_arc, dx_arc)
    if angle >= 0 && angle <= π/2
        d2 = abs(d_center - 20.0)
    else
        d2 = min(sqrt((x - 40.0)^2 + (y - 0.0)^2),
                 sqrt((x - 60.0)^2 + (y - (-20.0))^2))
    end
    
    # Vertical segment: from (60,-20) to (60,-50)
    d3 = distance_to_segment(x, y, 60.0, -20.0, 60.0, -50.0)
    
    return min(d1, min(d2, d3))
end

# ------------------------------------------------------------
# 3. Compute the Wall Curves for Visualization (Downward Turn)
# ------------------------------------------------------------
# Channel half–width (in meters)
channel_halfwidth = 10.0

# -- Horizontal segment (centerline: (x,0) for x in [0,40]) --
x_h = collect(range(0, 40, length=100))
# Centerline: y = 0.
# Walls are offset vertically by ± channel_halfwidth.
y_h_inner = fill(channel_halfwidth, length(x_h))   # Upper wall of horizontal segment
y_h_outer = fill(-channel_halfwidth, length(x_h))    # Lower wall of horizontal segment

# -- Circular arc (centerline: (40+20cosθ, -20+20sinθ) for θ in [π/2, 0]) --
θ_arc = collect(range(π/2, 0, length=100))
# Centerline arc:
x_arc_center = [40.0 + 20*cos(θi) for θi in θ_arc]
y_arc_center = [-20.0 + 20*sin(θi) for θi in θ_arc]
# Offset arcs: for a circular arc, the walls are concentric arcs.
# Inner wall: radius = 20 - channel_halfwidth = 10.
x_arc_inner = [40.0 + (20 + channel_halfwidth)*cos(θi) for θi in θ_arc]
y_arc_inner = [-20.0 + (20 + channel_halfwidth)*sin(θi) for θi in θ_arc]
# Outer wall: radius = 20 + channel_halfwidth = 30.
x_arc_outer = [40.0 + (20 - channel_halfwidth)*cos(θi) for θi in θ_arc]
y_arc_outer = [ -20.0 + (20 - channel_halfwidth)*sin(θi) for θi in θ_arc]

# -- Vertical segment (centerline: (60,y) for y in [-20,-50]) --
y_v = collect(range(-20, -50, length=100))
x_v_center = fill(60.0, length(y_v))
# Walls are offset horizontally by ± channel_halfwidth.
x_v_inner = fill(60.0 + channel_halfwidth, length(y_v))  # Left wall of vertical segment (x = 50)
x_v_outer = fill(60.0 - channel_halfwidth, length(y_v))  # Right wall of vertical segment (x = 70)

# Concatenate the pieces to form complete wall curves.
# Inner wall (continuous curve):
inner_wall_x = vcat(x_h, x_arc_inner, x_v_inner)
inner_wall_y = vcat(y_h_inner, y_arc_inner, y_v)
# Outer wall (continuous curve):
outer_wall_x = vcat(x_h, x_arc_outer, x_v_outer)
outer_wall_y = vcat(y_h_outer, y_arc_outer, y_v)

# ------------------------------------------------------------
# 4. FDTD Simulation Using the L–Shaped Waveguide (Radius Method, Downward Turn)
# ------------------------------------------------------------
function run_fdtd_lshaped_waveguide_simulation()
    # ---------------------------
    # Simulation and Grid Setup
    # ---------------------------
    c = 3.0e8                        # Speed of light (m/s)
    x_min, x_max = 0.0, 100.0          # Domain in x (m)
    y_min, y_max = -50.0, 50.0         # Domain in y (m)
    dx = 0.1                         # Grid spacing (m)
    dy = 0.1
    nx = Int(round((x_max - x_min)/dx)) + 1
    ny = Int(round((y_max - y_min)/dy)) + 1
    x = range(x_min, x_max, length=nx)
    y = range(y_min, y_max, length=ny)
    
    # ---------------------------
    # CFL Condition: Adaptive dt
    # ---------------------------
    cfl = 0.999999999999999          # CFL safety factor
    dt = cfl * dx / (c * sqrt(2))
    println("Adaptive dt: ", dt, " s")
    nsteps = 2000                  # Total number of time steps
    
    # ---------------------------
    # Field Initialization
    # ---------------------------
    Ez_old = zeros(Float64, nx, ny)  # Field at time step n-1
    Ez     = zeros(Float64, nx, ny)  # Field at time step n
    Ez_new = zeros(Float64, nx, ny)  # Field at time step n+1
    
    # Finite-difference coefficients (for second derivatives)
    Cx = (c * dt / dx)^2
    Cy = (c * dt / dy)^2
    
    # ---------------------------
    # Source Parameters
    # ---------------------------
    # Place the source inside the waveguide. For example, on the horizontal segment:
    source_x = 5.0
    source_y = 0.0
    i_source = Int(round((source_x - x_min)/dx)) + 1
    j_source = Int(round((source_y - y_min)/dy)) + 1
    omega = 2π * 1.0e7   # 10 MHz
    
    # ---------------------------
    # Animation Loop
    # ---------------------------
    anim = @animate for n = 1:nsteps
        t = n * dt
        
        # Update interior grid points with the standard FDTD scheme.
        @inbounds for i in 2:nx-1
            for j in 2:ny-1
                Ez_new[i,j] = 2 * Ez[i,j] - Ez_old[i,j] +
                              Cx * (Ez[i+1,j] - 2 * Ez[i,j] + Ez[i-1,j]) +
                              Cy * (Ez[i,j+1] - 2 * Ez[i,j] + Ez[i,j-1])
            end
        end
        
        # Inject the sinusoidal source.
        Ez_new[i_source, j_source] += sin(omega * t)
        
        # ---------------------------
        # Enforce PEC Outside the Waveguide (Radius Method)
        # ---------------------------
        # For every grid point (x[i], y[j]), if its distance from the centerline
        # exceeds the channel halfwidth, set the field to zero.
        for i in 1:nx
            for j in 1:ny
                if centerline_distance(x[i], y[j]) > channel_halfwidth
                    Ez_new[i,j] = 0.0
                end
            end
        end
        
        # ---------------------------
        # Apply Mur Boundary Conditions (Left & Right Edges)
        # ---------------------------
        alpha = (c*dt - dx) / (c*dt + dx)
        for j in 2:ny-1
            Ez_new[1,j]  = Ez[2,j]  + alpha*(Ez_new[2,j] - Ez[1,j])
            Ez_new[nx,j] = Ez[nx-1,j] + alpha*(Ez_new[nx-1,j] - Ez[nx,j])
        end
        
        # ---------------------------
        # Plot the Field and Overlay the Waveguide Walls
        # ---------------------------
        p = heatmap(x, y, Ez', aspect_ratio=1,
                    clims=(-0.1,0.1), color=:RdBu,
                    xlabel="x (m)", ylabel="y (m)",
                    title="2D FDTD L–Shaped Waveguide (Downward Turn), t = $(round(t, sigdigits=3)) s",
                    xlims=(x_min, x_max), ylims=(y_min, y_max))
        # Overlay the centerline segments (optional)
        plot!(p, [0, 40], [0, 0], lw=2, color=:black, label="Centerline")
        θs = collect(range(π/2, 0, length=100))
        arc_center_x = [40.0 + 20*cos(θi) for θi in θs]
        arc_center_y = [ -20.0 + 20*sin(θi) for θi in θs]
        plot!(p, arc_center_x, arc_center_y, lw=2, color=:black, label="")
        plot!(p, [60, 60], [-20, -50], lw=2, color=:black, label="")
        
        # Overlay the two waveguide walls.
        plot!(p, inner_wall_x, inner_wall_y, lw=2, color=:green, label="Inner Wall")
        plot!(p, outer_wall_x, outer_wall_y, lw=2, color=:orange, label="Outer Wall")
        
        # Advance the time steps.
        Ez_old .= Ez
        Ez     .= Ez_new
        
        p  # Return the current frame.
    end
    
    # Save the animation as a GIF.
    gif(anim, "fdtd_lshaped_waveguide_radius_downward.gif", fps=30)
end

# Run the simulation.
run_fdtd_lshaped_waveguide_simulation()
