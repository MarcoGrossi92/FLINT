import os
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib as mpl

mpl.rcParams.update({
    "font.size": 10,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
})

plt.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
})

root_dir = "./"
output_dir = "../docs/examples/images/"

folders = [
    os.path.join(root_dir, d)
    for d in os.listdir(root_dir)
    if os.path.isdir(os.path.join(root_dir, d))
]


# -------------------------------------------------------
# Function to read two-column data file
# -------------------------------------------------------
def read_two_column_file(filepath):
    x, y = [], []
    with open(filepath, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                try:
                    x.append(float(parts[0]))
                    y.append(float(parts[1]))
                except ValueError:
                    continue
    return x, y


# -------------------------------------------------------

for folder_ in folders:
    folder_name = os.path.basename(folder_)

    folder = os.path.join(folder_, "OUTPUT")
    if not os.path.exists(folder):
        continue

    # -------- Load reference file automatically --------
    ref_file = f"{folder_name}/eq-ref.txt"
    if not os.path.exists(ref_file):
        print(f"Reference file {ref_file} not found")
        continue

    ref_time, ref_temp = read_two_column_file(ref_file)

    # -------- Find ONE FLINT file --------
    flint_file = None
    for f in os.listdir(folder):
        if "CEA" in f:
            flint_file = os.path.join(folder, f)
            break

    if flint_file is None:
        continue

    flint_time, flint_temp = read_two_column_file(flint_file)

    # -------- Plot --------
    fig, ax = plt.subplots(figsize=(5, 4), facecolor="none")
    ax.set_facecolor("none")

    # Cantera (reference)
    ax.plot(
        ref_time,
        ref_temp,
        label="Cantera",
        color="gray",
        linestyle="None",
        marker="D",
        markersize=4
    )

    # FLINT
    ax.plot(
        flint_time,
        flint_temp,
        label="FLINT",
        color="green",
        linestyle="-"
    )

    ax.set_xscale("log")
    ax.set_xlabel("Mixture Ratio")
    ax.set_ylabel("Temperature")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    out_file = os.path.join(output_dir, f"{folder_name}-eq.svg")

    # Uncomment to save
    plt.savefig(out_file, bbox_inches="tight", transparent=True)
    plt.close(fig)

    plt.show()

    print(f"Processed {folder_name}")
