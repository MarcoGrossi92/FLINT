import os
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib as mpl

mpl.rcParams.update({
    "font.size": 12,          # base font size
    "axes.labelsize": 14,     # axis labels
    "xtick.labelsize": 12,    # x tick labels
    "ytick.labelsize": 12,    # y tick labels
    "legend.fontsize": 12,    # legend text
})

# --- Make background transparent ---
plt.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
})

# Root directory containing the folders
root_dir = "./"  # Change this to your root folder path
output_dir = "../../docs/examples/images/"

# The case folders of the batch reactor (first column of cases.txt); run from test/batch
with open(os.path.join(root_dir, "cases.txt")) as f:
    cases = [line.split()[0] for line in f if line.strip() and not line.startswith("#")]
folders = [os.path.join(root_dir, c) for c in cases if os.path.isdir(os.path.join(root_dir, c))]

# Optionally define ranges for each subplot (per folder)
# Format: {"folder_name": {"xlim": (xmin, xmax), "ylim": (ymin, ymax)}}
custom_ranges = {
    # "folder1": {"xlim": (0, 100), "ylim": (20, 80)},
    # "folder2": {"xlim": (10, 200)}
}

# The datasets of a case folder, in drawing order: file written by the drivers, legend label, style
datasets = [
    ("batch-CXX.dat",      "Cantera",        {"color": "gray", "linestyle": "None", "marker": "D", "markevery": 50, "markersize": 4}),
    ("batch-cantera.dat",  "FLINT Cantera",  {"color": "orange", "linestyle": "--"}),
    ("batch-explicit.dat", "FLINT Explicit", {"color": "green", "linestyle": "-"}),
    ("batch-general.dat",  "FLINT General",  {"color": "red", "linestyle": "-."}),
]

for idx, folder_ in enumerate(folders):
    folder = folder_ + "/"

    fig, ax = plt.subplots(figsize=(5, 4), facecolor="none")
    ax.set_facecolor("none")

    for file, label, style in datasets:
        filepath = os.path.join(folder, file)
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

        ax.plot(time, temp, label=label, **style)

    # Apply custom ranges if specified for this folder
    folder_name = os.path.basename(folder_)
    if folder_name in custom_ranges:
        if "xlim" in custom_ranges[folder_name]:
            ax.set_xlim(custom_ranges[folder_name]["xlim"])
        if "ylim" in custom_ranges[folder_name]:
            ax.set_ylim(custom_ranges[folder_name]["ylim"])

    ax.set_xlabel("Time")
    ax.set_ylabel("Temperature")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))

    ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    ax.tick_params(axis='both', labelsize=9)

    out_file = os.path.join(output_dir, f"{folder_name}.svg")
    plt.savefig(out_file, bbox_inches="tight", transparent=True)
    plt.close(fig)

    print(f"Saved {out_file}")
