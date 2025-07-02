#!/usr/bin/env python3
"""
gmx2pdb_strip.py

Description:
    Convert GROMACS .gro structures to PDB, then filter the resulting PDB to include only:
      - Standard amino acid residues
      - Palmitoylated cysteine backbone atoms (CYSP)
      - Zinc ions (ZN)
    Produces an alignment-ready PDB suitable for input to the RCSB RMSD alignment tool.

Usage:
    python gmx2pdb_select.py -i input.gro -o output.pdb

    Arguments:
        -i, --input    Path to the input GROMACS .gro file
        -o, --output   Path for the filtered PDB output (default: align-ready.pdb)

Requirements:
    - GROMACS CLI (`gmx`), with either `gmx2pdb` or `editconf` available
    - Biopython (install via `pip install biopython`)

Examples:
    python gmx2pdb_select.py -i pre_prod.gro -o pre_align_ready.pdb
    python gmx2pdb_select.py -i post_prod.gro -o post_align_ready.pdb
"""

import os
import subprocess
import argparse
from Bio.PDB import PDBParser, PDBIO, Select

class ProteinSelect(Select):
    # Standard amino acids
    standard_aa = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}

    def accept_atom(self, atom):
        residue = atom.get_parent()
        res_id = residue.get_id()
        # Protein residues
        if res_id[0] == " ":
            if residue.resname in self.standard_aa:
                return True
            if residue.resname == "CYSP":
                return atom.name in ["N", "CA", "C", "O", "CB", "SG"]
            return False
        # Heteroatoms: include only zinc ions
        if res_id[0].startswith("H_"):
            return residue.resname == "ZN" and atom.name == "ZN"
        return False


def run_gmx_convert(input_gro, output_pdb):
    """
    Converts a GROMACS .gro file to PDB using gmx editconf (or gmx2pdb if available).
    """
    # Try using gmx2pdb first
    cmd = ["gmx", "gmx2pdb", "-f", input_gro, "-o", output_pdb]
    try:
        subprocess.run(cmd, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        # Fallback to editconf
        cmd = ["gmx", "editconf", "-f", input_gro, "-o", output_pdb]
        subprocess.run(cmd, check=True)


def filter_structure(input_pdb, output_pdb):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(os.path.basename(input_pdb), input_pdb)
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, ProteinSelect())


def main():
    parser = argparse.ArgumentParser(description="Convert .gro to PDB and filter for RCSB alignment.")
    parser.add_argument("-i", "--input", required=True, help="Input .gro file")
    parser.add_argument("-o", "--output", default="align-ready.pdb", help="Filtered PDB output")
    args = parser.parse_args()

    gro_file = args.input
    temp_pdb = "temp_conv.pdb"

    print(f"Converting {gro_file} -> {temp_pdb}...")
    run_gmx_convert(gro_file, temp_pdb)

    print(f"Filtering atoms and writing to {args.output}...")
    filter_structure(temp_pdb, args.output)

    # Clean up
    os.remove(temp_pdb)
    print("Done.")

if __name__ == "__main__":
    main()
