from Bio import SeqIO
from Bio.PDB import PDBParser, PPBuilder

# Load PDB file
parser = PDBParser(QUIET=True)
x=input("Enter protein number")
structure = parser.get_structure("protein", x+".pdb")

# Extract sequence from structure
ppb = PPBuilder()
for pp in ppb.build_peptides(structure):
    print(pp.get_sequence())
