import matplotlib.pyplot as plt
import numpy as np

# File names (update these with your actual filenames if needed)
files = ['comp-batch-cantera.dat', 'comp-batch-canteraFor.dat', 'comp-batch-explicit.dat', 'comp-batch-general.dat']
labels = ['Cantera (CXX)', 'hydra-Cantera', 'hydra-explicit', 'hydra-general']
data = {}

# Step 1: Load data
for i, file in enumerate(files):
    with open(file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            name, time_str = parts
            time = float(time_str)
            if name not in data:
                data[name] = [0.0] * len(files)
            data[name][i] = time

# Step 2: Normalize to File 1
normalized_data = {}
for name, times in data.items():
    base = times[0]
    if base == 0.0:
        # Optional: skip normalization for this name
        normalized_data[name] = [0.0] * len(files)
    else:
        normalized_data[name] = [t / base for t in times]

# Step 3: Plot
names = sorted(normalized_data.keys())
x = np.arange(len(names))
width = 0.8 / len(files)

# Use a color-blind friendly palette
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

fig, ax = plt.subplots(figsize=(14, 6))

for i in range(len(files)):
    times = [normalized_data[name][i] for name in names]
    ax.bar(
        x + i * width,
        times,
        width,
        label=labels[i],
        color=colors[i % len(colors)],
        edgecolor='black',
        linewidth=1.0
    )

# Formatting
ax.set_ylabel('Normalized to Cantera-CXX Time', fontsize=12)
#ax.set_xlabel('Name', fontsize=12)
#ax.set_title('Normalized Execution Times per Sample', fontsize=14, pad=15)

ax.set_xticks(x + width * (len(files) - 1) / 2)
ax.set_xticklabels(names, rotation=0, fontsize=10)
ax.tick_params(axis='y', labelsize=10)

ax.grid(axis='y', linestyle='--', alpha=0.5)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=len(files), fontsize=10, frameon=False)

plt.tight_layout()
plt.show()