from openbabel import pybel
import os

def sdf_to_smiles(input_sdf):
    # Create output file name with .smiles extension
    base = os.path.splitext(input_sdf)[0]
    output_smi = base + ".smiles"

    with open(output_smi, "w") as out:
        for mol in pybel.readfile("sdf", input_sdf):
            smiles = mol.write("smi").strip()
            out.write(smiles + "\n")

    print(f"Conversion complete! SMILES written to {output_smi}")

if __name__ == "__main__":
    sdf_to_smiles("1a0q_ligand.sdf")
