from Bio import SeqIO
from Bio.PDB import PDBParser, PPBuilder

# Load PDB file
parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", "1a3n.pdb")

# Extract sequence from structure
ppb = PPBuilder()
for pp in ppb.build_peptides(structure):
    print(pp.get_sequence())
