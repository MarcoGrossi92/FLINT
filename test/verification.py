import os
import matplotlib.pyplot as plt

# Root directory containing the folders
root_dir = "./"  # Change this to your root folder path

# Collect all subfolders
folders = [os.path.join(root_dir, d) for d in os.listdir(root_dir) if os.path.isdir(os.path.join(root_dir, d))]

# Optionally define ranges for each subplot (per folder)
# Format: {"folder_name": {"xlim": (xmin, xmax), "ylim": (ymin, ymax)}}
custom_ranges = {
    # "folder1": {"xlim": (0, 100), "ylim": (20, 80)},
    # "folder2": {"xlim": (10, 200)}
}

# Define styles for up to 4 files per folder
styles = [
    {"color": "orange", "linestyle": "--", "marker": None},   # 1st file
    {"color": "green", "linestyle": "-", "marker": None},         # 2nd file
    {"color": "black", "linestyle": "None", "marker": "D", "markevery": 50, "markersize": 4},  # 3rd file
    {"color": "red", "linestyle": "-.", "marker": None},          # 4th file
]

# Create tiled plots
n_folders = len(folders)

# Automatically choose a near-square layout for clarity
cols = int(n_folders**0.5)
if cols * cols < n_folders:
    cols += 1
rows = (n_folders + cols - 1) // cols

# Scale figure size based on rows and cols
fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 4*rows), sharex=False, sharey=False)
axes = axes.flatten()

for idx, folder_ in enumerate(folders):
    folder = folder_ + "/OUTPUT/"
    ax = axes[idx]
    files = [f for f in os.listdir(folder) if 'batch' in f]

    for i, file in enumerate(files):
        filepath = os.path.join(folder, file)
        # Ensure the file path exists before opening
        if not os.path.isfile(filepath):
            print(f"Warning: file not found {filepath}")
            continue

        # Read data manually without pandas
        time, temp = [], []
        with open(filepath, "r") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) == 2:
                    try:
                        t, val = float(parts[0]), float(parts[1])
                        time.append(t)
                        temp.append(val)
                    except ValueError:
                        continue

        # Choose style based on file index
        style = styles[i % len(styles)]
        ax.plot(
            time,
            temp,
            label=os.path.splitext(file)[0],
            color=style["color"],
            linestyle=style.get("linestyle", "-"),
            marker=style.get("marker", None),
            markevery=style.get("markevery", None),
            markersize=style.get("markersize", None)
        )

    # Apply custom ranges if specified for this folder
    folder_name = os.path.basename(folder_)
    if folder_name in custom_ranges:
        if "xlim" in custom_ranges[folder_name]:
            ax.set_xlim(custom_ranges[folder_name]["xlim"])
        if "ylim" in custom_ranges[folder_name]:
            ax.set_ylim(custom_ranges[folder_name]["ylim"])

    ax.set_title(folder_name)
    ax.set_xlabel("Time")
    ax.set_ylabel("Temperature")
    ax.legend(fontsize=8)

# Remove unused subplots if any
for j in range(n_folders, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()
