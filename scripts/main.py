import sys
import os
import numpy as np
import matplotlib.pyplot as plt

def main():
    # The dirs come from CLI (included in the makefile ofc)
    if len(sys.argv) != 3:
        print("Usage: python3 main.py <data_dir> <plot_dir>")
        sys.exit(1)

    data_dir = sys.argv[1]
    plot_dir = sys.argv[2]

    # Create output directory if non-existent
    os.makedirs(plot_dir, exist_ok=True)
    print(f"Reading simulation data from {data_dir}...")

    # Energy Evolution
    energy_file = os.path.join(data_dir, "energy.dat")
    if os.path.exists(energy_file):
        data = np.loadtxt(energy_file)
        steps = data[:, 0]
        energy = data[:, 1]
        
        plt.figure(figsize = (10, 6))
        plt.plot(steps, energy, label = "Total Energy", color = "#1f77b4", linewidth = 1.5)
        plt.xlabel("Monte Carlo Step")
        plt.ylabel("Energy (Reduced Units)")
        plt.title("System Energy Evolution")
        plt.grid(True, linestyle = "--", alpha = 0.6)
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, "energy_evolution.png"), dpi = 300)
        plt.close()
        print(" > Generated energy_evolution.png")

    # Torsion angle distribution
    tors_file = os.path.join(data_dir, "torsions.dat")
    if os.path.exists(tors_file):
        data = np.loadtxt(tors_file)
        # index 1: angles in radians
        angles_rad = data[:, 1]
        # Convert to degrees
        angles_deg = np.degrees(angles_rad)
        
        plt.figure(figsize = (10, 6))
        plt.hist(angles_deg, bins = 90, density = True, color = "#2ca02c", alpha = 0.75, edgecolor = "black", linewidth = 0.5)
        plt.xlabel("Torsion Angle (Degrees)")
        plt.ylabel("Probability Density")
        plt.title("Torsion Angle Distribution")
        plt.xlim(-180, 180)
        plt.grid(True, linestyle = "--", alpha = 0.6)
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, "torsion_distribution.png"), dpi = 300)
        plt.close()
        print(" > Generated torsion_distribution.png")

    # Structural properties (gyration and end-to-end distance)
    struc_file = os.path.join(data_dir, "structure.dat")
    if os.path.exists(struc_file):
        data = np.loadtxt(struc_file)
        steps = data[:, 0]
        ree = data[:, 1]
        rg = data[:, 2]
        
        fig, ax1 = plt.subplots(figsize = (10, 6))

        color1 = "#ff7f0e"
        ax1.set_xlabel("Monte Carlo Step")
        ax1.set_ylabel("End-to-End Distance ($R_{ee}$)", color = color1)
        ax1.plot(steps, ree, color = color1, alpha = 0.8, label = "$R_{ee}$")
        ax1.tick_params(axis = "y", labelcolor = color1)
        ax1.grid(True, linestyle = "--", alpha = 0.6)

        # Second y-axis that shares the same x-axis
        ax2 = ax1.twinx()  
        color2 = "#9467bd"
        ax2.set_ylabel("Radius of Gyration ($R_g$)", color=color2)
        ax2.plot(steps, rg, color = color2, alpha = 0.8, label = "$R_g$")
        ax2.tick_params(axis = "y", labelcolor = color2)

        plt.title("Evolution of Polymer Dimensions")
        fig.tight_layout()
        plt.savefig(os.path.join(plot_dir, "structure_evolution.png"), dpi = 300)
        plt.close()
        print(" > Generated structure_evolution.png")

if __name__ == "__main__":
    main()