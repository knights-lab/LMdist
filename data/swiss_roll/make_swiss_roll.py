import matplotlib.pyplot as plt
from sklearn import manifold, datasets
from matplotlib import cm
from matplotlib.colors import Normalize, to_hex
import numpy as np

# create dataset
sr_pts, sr_color = datasets.make_swiss_roll(300, random_state=0)

# visualize dataset (3D)
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection="3d")
fig.add_axes(ax)
ax.scatter( sr_pts[:,0], sr_pts[:,1], sr_pts[:,2], c=sr_color, s=50, alpha=0.9, cmap="plasma")
ax.set_title("Swiss Roll 3D", size=16, fontweight="bold", fontname="Arial", x=0.25, y=0.9)
ax.view_init(azim=-75, elev=10)
_ = ax.text2D(0.15, 0.71, s="n samples = 300", transform=ax.transAxes, fontname="Arial", size=12)
plt.savefig('swiss_roll_3d.png', dpi=1200)
plt.savefig('swiss_roll_3d_lowres.png', dpi=350)
plt.show()


# export dataset (3D)
np.savetxt("swiss_roll.txt", sr_pts)

# export colors & manifold numbers
norm = Normalize(sr_color.min(), sr_color.max())
cmap = cm.plasma
sr_color_hex = []
sr_color_num = []
for col in sr_color:
	sr_color_hex.append( to_hex( cmap(norm(col)) ) )
	sr_color_num.append(str(norm(col)))
with open("swiss_roll_colors.txt", "w") as color_file:
	color_file.write("\n".join(sr_color_hex))
with open("swiss_roll_manifold.txt", "w") as color_file2:
	color_file2.write("\n".join(sr_color_num))
