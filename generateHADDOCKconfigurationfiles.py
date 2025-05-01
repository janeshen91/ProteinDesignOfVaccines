import csv
from pathlib import Path

def generate_configs_from_csv(csv_path, output_dir, sampling=400, tolerance=5):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fixed_molecule2 = "/output/data/vx22Web/protein2.pdb"

    with open(csv_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)

        for row in reader:
            run_dir = row["run_dir"]
            mol1 = row["molecule_path1"]
            mol2 = fixed_molecule2
            ambig_fname = row["ambig_fname"]
            unambig_fname = row["unambig_fname"]

            run_id = Path(run_dir).name
            output_path = output_dir / f"{run_id}.toml"

            config = f"""\
# ====================================================================
# Protein-protein docking example with NMR-derived ambiguous interaction restraints

# directory in which the scoring will be done
run_dir = "{run_dir}"

# molecules to be docked
molecules =  [
    "{mol1}",
    "{mol2}"
    ]

# ====================================================================
# Parameters for each stage are defined below, prefer full paths
# ====================================================================
[topoaa]

[rigidbody]
tolerance = {tolerance}
##got this based on reading the VX22 paper on epitope interaction with FAB
ambig_fname = "{ambig_fname}"
# got it from haddock-restraint restraint
unambig_fname = "{unambig_fname}"
#  sampling
sampling = {sampling}

[caprieval]



# ====================================================================
"""
            with open(output_path, 'w') as f:
                f.write(config)
            print(f"Generated config for {run_id}: {output_path}")

# Example usage:
# generate_configs_from_csv("input_table.csv", "generated_configs")
