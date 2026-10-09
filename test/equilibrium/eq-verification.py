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

# Run from test/equilibrium, after test-CEA: reference/<mechanism>.dat (Cantera, test-equilCXX) against
# <mechanism>/FLINT-CEA.txt (test-CEA)
root_dir = "./"
output_dir = "../../docs/examples/images/"

mechanisms = sorted(f[:-4] for f in os.listdir(os.path.join(root_dir, "reference")) if f.endswith(".dat"))

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

for folder_name in mechanisms:
    ref_file = os.path.join(root_dir, "reference", f"{folder_name}.dat")
    flint_file = os.path.join(root_dir, folder_name, "FLINT-CEA.txt")
    if not os.path.exists(flint_file):
        print(f"FLINT file {flint_file} not found (run test-CEA)")
        continue

    ref_time, ref_temp = read_two_column_file(ref_file)
    flint_time, flint_temp = read_two_column_file(flint_file)
    with open(ref_file) as f:
        pressure_sweep = any(line.startswith("# p [Pa]") for line in f)

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
        markevery=50,
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
    ax.set_xlabel("Pressure [Pa]" if pressure_sweep else "Mixture Ratio")
    ax.set_ylabel("Temperature")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    out_file = os.path.join(output_dir, f"{folder_name}-eq.svg")

    # Uncomment to save
    plt.savefig(out_file, bbox_inches="tight", transparent=True)
    plt.close(fig)


    print(f"Processed {folder_name}")
