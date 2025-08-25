#!/usr/bin/env python3
import os
from Bio.PDB import PDBParser, PPBuilder
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# === CONFIGURATION ===
main_folder = "/path/to/main_database"   # change this to your database path
output_base = "/path/to/output_sequences"  # where sequences will be stored

parser = PDBParser(QUIET=True)
ppb = PPBuilder()

# Make sure output base exists
os.makedirs(output_base, exist_ok=True)

# Walk through all subfolders
for protein_folder in os.listdir(main_folder):
    protein_path = os.path.join(main_folder, protein_folder)

    # Only process if it's a folder
    if not os.path.isdir(protein_path):
        continue

    # Find a .pdb file inside this protein folder
    pdb_file = None
    for file in os.listdir(protein_path):
        if file.lower().endswith(".pdb"):
            pdb_file = os.path.join(protein_path, file)
            break

    if pdb_file is None:
        print(f"❌ No PDB file found in {protein_folder}")
        continue

    print(f"✅ Processing {pdb_file}")

    # Parse structure
    structure = parser.get_structure(protein_folder, pdb_file)

    # Extract sequence(s)
    sequences = []
    for i, pp in enumerate(ppb.build_peptides(structure)):
        seq = pp.get_sequence()
        record = SeqRecord(seq, id=f"{protein_folder}_chain{i+1}", description="")
        sequences.append(record)

    # Prepare output folder for this protein
    output_folder = os.path.join(output_base, protein_folder + "_sequence")
    os.makedirs(output_folder, exist_ok=True)

    # Save sequence(s) to FASTA
    output_file = os.path.join(output_folder, protein_folder + ".fasta")
    SeqIO.write(sequences, output_file, "fasta")

    print(f"   → Sequences saved to {output_file}")

