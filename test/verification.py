import os
import matplotlib.pyplot as plt

plt.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "axes.edgecolor": "#666666",
    "axes.labelcolor": "#444444",
    "xtick.color": "#444444",
    "ytick.color": "#444444",
    "text.color": "#444444",
    "grid.color": "#888888",
    "grid.alpha": 0.3,
})

# Root directory containing the folders
root_dir = "./"  # Change this to your root folder path
output_dir = "../docs/images/"

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
    {"color": "black", "linestyle": "None", "marker": "D", "markevery": 50, "markersize": 4},
    {"color": "orange", "linestyle": "--", "marker": None},
    {"color": "green", "linestyle": "-", "marker": None},
    {"color": "red", "linestyle": "-.", "marker": None},
]

for idx, folder_ in enumerate(folders):
    folder = folder_ + "/OUTPUT/"

    fig, ax = plt.subplots(figsize=(5, 4))

    order = [0, 1, 2, 3]
    files = [f for f in os.listdir(folder) if 'batch' in f]
    files.sort()
    files = [files[i] for i in order if i < len(files)]

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
        legend_labels = [
            "Cantera",
            "FLINT Cantera",
            "FLINT Explicit",
            "FLINT General",
        ]
        label = legend_labels[i] if i < len(legend_labels) else f"Series {i+1}"

        ax.plot(
            time,
            temp,
            label=label,
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

        style = styles[i % len(styles)]
        ax.plot(time, temp, label=file, **style)

    ax.set_title(folder_name)
    ax.set_xlabel("Time")
    ax.set_ylabel("Temperature")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    out_file = os.path.join(output_dir, f"{folder_name}.svg")
    plt.savefig(out_file, bbox_inches="tight", transparent=True)
    plt.close(fig)

    print(f"Saved {out_file}")
