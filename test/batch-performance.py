import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

mpl.rcParams.update({
    "font.size": 14,          # base font size
    "axes.titlesize": 18,     # title
    "axes.labelsize": 16,     # axis labels
    "xtick.labelsize": 14,    # x tick labels
    "ytick.labelsize": 14,    # y tick labels
    "legend.fontsize": 14,    # legend text
})

# --- Make background transparent ---
plt.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",  # keep text as text
})

# File names (update these with your actual filenames if needed)
files = ['comp-batch-cantera.dat', 'comp-batch-canteraFor.dat', 'comp-batch-explicit.dat', 'comp-batch-general.dat']
labels = ['Cantera (C++)', 'FLINT-Cantera', 'FLINT-explicit', 'FLINT-general']
data = {}

# Step 1: Load data
for i, file in enumerate(files):
    with open(file, 'r') as f:
        for line in f:
            if 'WD' in line or 'Pelucchi' in line:
                continue
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

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

fig, ax = plt.subplots(figsize=(14, 6), facecolor="none")
ax.set_facecolor("none")

for i in range(len(files)):
    times = [normalized_data[name][i] for name in names]
    ax.bar(
        x + i * width,
        times,
        width,
        label=labels[i],
        color=colors[i % len(colors)],
        edgecolor='none'
    )

ax.set_ylabel('Normalized to Cantera Time')

ax.set_xticks(x + width * (len(files) - 1) / 2)
ax.set_xticklabels(names)

ax.grid(axis='y', linestyle='--', alpha=0.3)

ax.legend(
    loc='upper center',
    bbox_to_anchor=(0.5, 1.15),
    ncol=len(files),
    frameon=False
)

plt.tight_layout()
plt.savefig("../docs/user/images/barplot.svg", format="svg", transparent=True)
plt.show()
plt.close()