using Plots

# Define the parabolic function
function parabolic_curve(y, f)
    return y.^2 / (4 * f)
end

# Set parameters
f = 10.0  # Focal length
y_range = -15.0:0.1:15.0  # Range of y values

# Calculate the parabolic curve
x_values = parabolic_curve(y_range, f)

# Plot the parabolic antenna
plot(x_values, y_range, label="Parabolic Antenna", xlabel="x", ylabel="y", title="2D Parabolic Antenna", aspect_ratio=:equal)

# Add the focal point
scatter!([f], [0], label="Focal Point", color=:red)

# Add a horizontal line segment from x = 0 to x = -3 at y = 0
plot!([0, -3], [0, 0], label="Line Segment", color=:blue, linewidth=2)

# Add a triangle approaching -ylims! from x = -3
triangle_x = [-3, -3, -3, -3]  # x-coordinates of the triangle
triangle_y = [0, -25, -25, 0]  # y-coordinates of the triangle (extending to y = -25)
plot!(triangle_x, triangle_y, label="Tower", color=:green, linewidth=2, fill=(0, :green))

# Adjust x and y limits of the plot
xlims!(-10, 80)  # Set x-axis limits
ylims!(-25, 25)  # Set y-axis limits