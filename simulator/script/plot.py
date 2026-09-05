from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Locate the build directory relative to this script file
script_dir = Path(__file__).resolve().parent
build_dir = script_dir.parent / 'build'

t_path = build_dir / 'T.dat'
phi_path = build_dir / 'phi.dat'

# Load 2D matrix data
T = np.loadtxt(t_path)
phi = np.loadtxt(phi_path)

# Create a 1x2 subplot figure
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Plot T matrix
im1 = axes[0].imshow(T, cmap='viridis', origin='lower', aspect='auto')
axes[0].set_title('Temperature ($T$)')
axes[0].set_xlabel('X')
axes[0].set_ylabel('Y')
fig.colorbar(im1, ax=axes[0], fraction=0.046, pad=0.04, label='T value')

# Plot phi matrix
im2 = axes[1].imshow(phi, cmap='plasma', origin='lower', aspect='auto')
axes[1].set_title(r'Field ($\phi$)')
axes[1].set_xlabel('X')
axes[1].set_ylabel('Y')
fig.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.04, label=r'$\phi$ value')

plt.tight_layout()

# Save combined figure to build directory and display
plt.savefig(build_dir / 'combined_plot.png', dpi=300, bbox_inches='tight')
plt.show()