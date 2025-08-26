from openbabel import pybel
import sys

def sdf_to_smiles(input_sdf, output_smi):
    # Open output file
    with open(output_smi, "w") as out:
        # Read molecules from SDF
        for mol in pybel.readfile("sdf", input_sdf):
            smiles = mol.write("smi").strip()
            out.write(smiles + "\n")
    print(f"Conversion complete! SMILES written to {output_smi}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python sdf_to_smiles.py input.sdf output.smi")
    else:
        sdf_to_smiles(sys.argv[1], sys.argv[2])
