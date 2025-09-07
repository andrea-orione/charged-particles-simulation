import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

N_TIMESTAMPS = 10000 # (keep under 500)
N_PARTICLES = 4

def animate(i):
    for j,body in enumerate(bodies):
        body.set_data([data[j*3 + 1][i]],[data[j*3 + 2][i]])

# IMPORT DATA
data = np.genfromtxt("../output/main.csv", delimiter=",", unpack=True)

# DRAW
# plt.style.use('dark_background')
fig, ax = plt.subplots()

bodies = []
ax.set_title("Motion")
ax.set_aspect('equal')
ax.set_xlim((-50,50))
ax.set_ylim((-50.,50.))

for i in range(N_PARTICLES):
    bodies.append(ax.plot([], [], 'b.', markersize=5)[0])

ax.grid(True, lw=0.3)

indices = range(N_TIMESTAMPS)
anim = FuncAnimation(fig, animate, frames=indices, interval=0.01, repeat=True)
plt.show()
