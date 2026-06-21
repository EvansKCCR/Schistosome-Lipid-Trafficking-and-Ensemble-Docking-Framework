import os
import glob
import subprocess
import numpy as np
import pandas as pd
import MDAnalysis as mda
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

# ==============================
# USER SETTINGS
# ==============================

BASE_DIR = "./"
PROTEINS = ["RP1", "RP3", "RP4"]
LIGANDS = ["CE18_1", "CE18_2", "CE18_3", "CE20_4", "CE20_5","CE16_0"]

CUTOFF = 4.5
FREQ_THRESHOLD = 0.5   # 🔥 KEY FILTER

# ==============================
# STEP 1: GENERATE POSES
# ==============================

def generate_pose_files(work_dir):

    pdbqt_files = glob.glob(os.path.join(work_dir, "*.pdbqt"))

    for f in pdbqt_files:
        base = os.path.splitext(os.path.basename(f))[0]

        cmd = f"obabel {f} -O {work_dir}/poses.pdb -m"
        subprocess.run(cmd, shell=True,
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)

        split_files = sorted(glob.glob(os.path.join(work_dir, "poses*.pdb")))

        for i, sf in enumerate(split_files):
            new_name = os.path.join(work_dir, f"{base}_pose{i+1}.pdb")
            os.rename(sf, new_name)

# ==============================
# VALIDATION
# ==============================

def is_valid_pdb(pose):
    try:
        if os.path.getsize(pose) < 200:
            return False
        with open(pose) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    return True
        return False
    except:
        return False

# ==============================
# CONTACT ANALYSIS
# ==============================

def compute_contacts(receptor_file, pose_files):

    u = mda.Universe(receptor_file)
    protein = u.select_atoms("protein")

    contact_counts = defaultdict(int)
    total_frames = len(pose_files)

    print(f"Valid poses: {total_frames}")

    for pose in pose_files:
        ligand = mda.Universe(pose).atoms

        for res in protein.residues:
            distances = mda.lib.distances.distance_array(
                res.atoms.positions,
                ligand.positions
            )

            if np.min(distances) < CUTOFF:
                contact_counts[f"{res.resname}-{res.resid}"] += 1

    if total_frames == 0:
        return pd.DataFrame(), pd.DataFrame()

    df = pd.DataFrame(
        [(k, v / total_frames) for k, v in contact_counts.items()],
        columns=["Residue", "Frequency"]
    ).sort_values("Frequency", ascending=False)

    # 🔥 FILTERED DATA
    df_filtered = df[df["Frequency"] > FREQ_THRESHOLD]

    return df, df_filtered

# ==============================
# STATISTICS
# ==============================

def compute_statistics(df, protein, ligand):

    if df.empty:
        return None

    mean_freq = df["Frequency"].mean()
    high_contacts = (df["Frequency"] >= 0.8).sum()

    freqs = df["Frequency"].values
    freqs = freqs / freqs.sum()
    entropy = -np.sum(freqs * np.log(freqs + 1e-9))

    return {
        "Protein": protein,
        "Ligand": ligand,
        "MeanFrequency": mean_freq,
        "HighFreqResidues": high_contacts,
        "Entropy": entropy
    }

# ==============================
# INDIVIDUAL PLOTS
# ==============================

def plot_individual(df, prefix):

    if df.empty:
        return

    df = df[df["Frequency"] > FREQ_THRESHOLD]
    if df.empty:
        return

    # Bar plot
    plt.figure()
    top = df.head(20)
    plt.barh(top["Residue"], top["Frequency"])
    plt.gca().invert_yaxis()
    plt.xlabel("Frequency")
    plt.title(prefix + " (Filtered)")
    plt.tight_layout()
    plt.savefig(prefix + "_bar_filtered.png", dpi=600)
    plt.close()

    # Histogram (keep FULL distribution)
    plt.figure()
    plt.hist(df["Frequency"], bins=20)
    plt.xlabel("Frequency")
    plt.ylabel("Count")
    plt.title(prefix + " Distribution")
    plt.tight_layout()
    plt.savefig(prefix + "_hist.png", dpi=600)
    plt.close()

# ==============================
# SPLIT PLOTS BY PROTEIN
# ==============================

def plot_split_by_protein(combined):

    for protein in combined["Protein"].unique():

        df_p = combined[combined["Protein"] == protein]

        # FILTERED
        df_pf = df_p[df_p["Frequency"] > FREQ_THRESHOLD]

        # =========================
        # BOXPLOT (FILTERED)
        # =========================
        if not df_pf.empty:
            plt.figure()
            sns.boxplot(data=df_pf, x="Ligand", y="Frequency")
            plt.title(f"{protein} Contact Frequency (>0.5)")
            plt.tight_layout()
            plt.savefig(f"{protein}_boxplot_filtered.png", dpi=600)
            plt.close()

        # =========================
        # HISTOGRAM (FULL)
        # =========================
        plt.figure()
        sns.histplot(data=df_p, x="Frequency", hue="Ligand", bins=20, kde=True)
        plt.title(f"{protein} Frequency Distribution")
        plt.tight_layout()
        plt.savefig(f"{protein}_histogram.png", dpi=600)
        plt.close()

        # =========================
        # HEATMAP (FILTERED)
        # =========================
        if not df_pf.empty:
            pivot = df_pf.pivot_table(
                index="Residue",
                columns="Ligand",
                values="Frequency",
                fill_value=0
            )

            if not pivot.empty:
                plt.figure(figsize=(10, 12))
                sns.heatmap(pivot, cmap="viridis")
                plt.title(f"{protein} Residue Contact Heatmap (>0.5)")
                plt.tight_layout()
                plt.savefig(f"{protein}_heatmap_filtered.png", dpi=600)
                plt.close()

# ==============================
# SIDE-BY-SIDE COMPARISON
# ==============================

def plot_side_by_side(combined):

    combined_f = combined[combined["Frequency"] > FREQ_THRESHOLD]

    if combined_f.empty:
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    for i, protein in enumerate(["RP1", "RP3"]):
        df_p = combined_f[combined_f["Protein"] == protein]
        sns.boxplot(data=df_p, x="Ligand", y="Frequency", ax=axes[i])
        axes[i].set_title(protein)

    plt.tight_layout()
    plt.savefig("RP1_vs_RP3_filtered.png", dpi=600)
    plt.close()

# ==============================
# MAIN PIPELINE
# ==============================

all_data = []
stats_data = []

for protein in PROTEINS:

    receptor = os.path.join(BASE_DIR, protein, "receptor.pdb")

    for ligand in LIGANDS:

        work_dir = os.path.join(BASE_DIR, protein, ligand)

        if not os.path.exists(work_dir):
            continue

        print(f"\nProcessing {protein} - {ligand}")

        generate_pose_files(work_dir)

        pose_files = glob.glob(os.path.join(work_dir, "*_pose*.pdb"))
        pose_files = [p for p in pose_files if is_valid_pdb(p)]

        if not pose_files:
            print("No valid poses")
            continue

        df, df_filtered = compute_contacts(receptor, pose_files)

        prefix = os.path.join(work_dir, f"{protein}_{ligand}")

        # SAVE BOTH
        df.to_csv(prefix + "_contacts_all.csv", index=False)
        df_filtered.to_csv(prefix + "_contacts_filtered.csv", index=False)

        plot_individual(df, prefix)

        df["Protein"] = protein
        df["Ligand"] = ligand
        all_data.append(df)

        stats = compute_statistics(df, protein, ligand)
        if stats:
            stats_data.append(stats)

# ==============================
# GLOBAL ANALYSIS
# ==============================

if all_data:

    combined = pd.concat(all_data)
    combined.to_csv("ALL_contact_summary.csv", index=False)

    plot_split_by_protein(combined)
    plot_side_by_side(combined)

if stats_data:
    stats_df = pd.DataFrame(stats_data)
    stats_df.to_csv("binding_statistics.csv", index=False)

    print("\n=== STATISTICS ===")
    print(stats_df)

print("\nPipeline complete.")