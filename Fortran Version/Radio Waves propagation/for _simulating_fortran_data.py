import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

# Simulation grid and reflector parameters
Nx, Ny = 400, 300
xmin, xmax, ymin, ymax = -10.0, 80.0, -25.0, 25.0
xgrid, ygrid = np.linspace(xmin, xmax, Nx), np.linspace(ymin, ymax, Ny)
f_param, source_y = 10.0, 0.0
source_x = f_param
y_curve = np.linspace(-15.0, 15.0, 300)
x_curve = y_curve**2 / (4.0 * f_param)

# Load time-series data
def load_ez_time_series(filename):
    blocks = []
    with open(filename, 'r') as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith("Step:"):
            step, time_val = int(lines[i].split()[1]), float(lines[i].split()[3])
            i += 1
            ez_data = [np.array(list(map(float, lines[j].split()))) for j in range(i, i + Ny)]
            blocks.append({'step': step, 'time': time_val, 'Ez': np.array(ez_data)})
            i += Ny
        else:
            i += 1
    return blocks

time_series = load_ez_time_series("Ez_time_series.dat")
if not time_series:
    raise RuntimeError("No data blocks found in the file.")

# Set up the figure and axis
fig, ax = plt.subplots(figsize=(8, 6))
contour = ax.contourf(xgrid, ygrid, time_series[0]['Ez'], levels=50, cmap='RdBu_r')
ax.set_title(f"Step: {time_series[0]['step']}, Time: {time_series[0]['time']:.3e}")
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
cbar = fig.colorbar(contour, ax=ax, label="Ez")

# Animation update function
def update_frame(frame):
    ax.clear()
    block = time_series[frame]
    ax.contourf(xgrid, ygrid, block['Ez'], levels=50, cmap='RdBu_r')
    ax.set_title(f"Step: {block['step']}, Time: {block['time']:.3e}")
    ax.plot(x_curve, y_curve, 'k', linewidth=3, label='Reflector')
    ax.plot([-3, 0], [0, 0], 'k', linewidth=3)
    ax.plot([-6, -3], [ymin, 0], 'k', linewidth=3)
    ax.plot([0, -3], [ymin, 0], 'k', linewidth=3)
    ax.scatter([source_x], [source_y], color='red', s=100, zorder=5, label='Source')
    ax.legend(loc='upper right')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

# Create and show animation
anim = animation.FuncAnimation(fig, update_frame, frames=len(time_series), interval=50, blit=False)
plt.show()

# Save animation (optional)
anim.save("Ez_animation_overlay.gif", writer='imagemagick', fps=30)