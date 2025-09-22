import numpy as np
import matplotlib.pyplot as plt

#aspen data for n-hexane + cycloheptane at 353.15 K source: Aspen Plus done by the man the myth the ledgend Rylin Singh
xy_data = {
    'x': [0.00,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0],
    'y': [0.0, 0.1824197, 0.3218958, 0.430985, 0.5181176, 0.5890586, 0.6478281, 0.6972854, 0.7395021,0.7760066,0.8079447, 0.8361889, 0.8614127
          , 0.8841429, 0.9047961, 0.9237058, 0.9411415, 0.9573229, 0.9724306, 0.9866146, 1.0]
          }
#separating the data into arrays
x_pts= np.array(xy_data['x'])
y_pts = np.array(xy_data['y'])
#Design specs wrt n-hexane...most volatile component
xD = 0.95
xB = 0.05
xF = 0.5
R = 2.0



def y_equil(x): return np.interp(x, x_pts, y_pts)
def x_from_y(y): return np.interp(y, y_pts, x_pts)

# Rectifying line
def y_rect(x, R=R):
    return (R/(R+1.0))*x + xD/(R+1.0)

yR_at_xF = y_rect(xF, R)

# Correct stripping line: pass through (xB, xB) and (xF, yR_at_xF)
slope_strip = (yR_at_xF - xB) / (xF - xB)
def y_strip(x):
    return slope_strip * (x - xB) + xB

# Stepping (top-down)
stair_x = [xD]
stair_y = [xD]
current_y = xD
max_steps = 300
tol = 1e-6
feed_stage_index = None

for step in range(max_steps):
    # horizontal to equilibrium curve
    x_new = x_from_y(current_y)
    stair_x.append(x_new)
    stair_y.append(current_y)
    if x_new <= xB + tol:
        break
    # which operating line?
    if x_new > xF + 1e-12:
        next_y = y_rect(x_new, R)
    else:
        if feed_stage_index is None:
            feed_stage_index = len(stair_x)-1
        next_y = y_strip(x_new)
    stair_x.append(x_new)
    stair_y.append(next_y)
    current_y = next_y

num_stage_steps = (len(stair_x)-1)//2 + ((len(stair_x)-1)%2)
# Plot
fig, ax = plt.subplots(figsize=(8,8))
ax.plot(x_pts, y_pts, label='Equilibrium curve', linewidth=2)
ax.plot([0,1],[0,1],'--',label='Diagonal y=x')
xx = np.linspace(0,1,201)
ax.plot(xx, y_rect(xx, R), label=f'Rectifying (R={R})')
ax.plot(xx, y_strip(xx), label='Stripping (anchored at diagonal xB)')
ax.vlines(xF, 0, 1, linestyles=':', label='q-line (q=1)')

# staircase
for i in range(0, len(stair_x)-1, 2):
    x_h0 = stair_x[i]; y_h0 = stair_y[i]
    x_h1 = stair_x[i+1]; y_h1 = stair_y[i+1]
    ax.plot([x_h0, x_h1],[y_h0,y_h1], color='orange')  # horizontal
    if i+2 < len(stair_x):
        x_v = stair_x[i+1]; y_v0 = stair_y[i+1]; y_v1 = stair_y[i+2]
        ax.plot([x_v,x_v],[y_v0,y_v1], color='orange')  # vertical

ax.scatter([xD],[xD], marker='D', color='green', label='xD')
ax.scatter([xB],[xB], marker='s', color='purple', label='xB (diagonal)')
ax.scatter([xF],[yR_at_xF], marker='o', color='red', label='feed intersection (rect line)')

ax.set_xlim(-0.02,1.02); ax.set_ylim(-0.02,1.02)
ax.set_xlabel('x (liquid mole fraction n-hexane)'); ax.set_ylabel('y (vapor mole fraction n-hexane)')
ax.set_title('Corrected McCabe-Thiele (stripping line anchored at diagonal xB)')
ax.grid(True); ax.legend(loc='lower right')
plt.gca().set_aspect('equal', adjustable='box')

# summary text
txt = f"R = {R}\nStages (est): {num_stage_steps}\nFeed stage index: {feed_stage_index}"
ax.text(0.02, 0.95, txt, transform=ax.transAxes, va='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.6))
plt.show()

print("Rectifying y at xF: ", yR_at_xF)
print("Stripping line slope: ", slope_strip)
print("Estimated theoretical stages: ", num_stage_steps)
print("Feed stage index (approx): ", feed_stage_index)

