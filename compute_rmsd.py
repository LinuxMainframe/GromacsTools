#!/usr/bin/env python3
"""
compute_rmsd.py

Fetches a reference PDB from RCSB, loads your converted PDB, aligns them,
and reports RMSD for several atom selections: all‐atom, backbone, loops,
extracellular (EC) and intracellular (IC) regions.

Usage:
    python compute_rmsd.py \
      --converted myprot_conv.pdb \
      --pdb-id 6A94 \
      --chain A \
      --loops 45-55,120-130 \
      --ec 1-30 \
      --ic 200-232 \
      --out rmsd.csv

Requirements:
    - Biopython
    - Internet access to fetch PDB
"""

import argparse
import csv
from Bio.PDB import PDBParser, PDBList, Superimposer, is_aa

def parse_ranges(rangestr):
    """Convert '45-55,120-130' → [(45,55),(120,130)]."""
    regions = []
    for part in rangestr.split(','):
        start,end = part.split('-')
        regions.append((int(start), int(end)))
    return regions

def atom_groups(chain, loop_regs, ec_regs, ic_regs):
    """
    Given a Chain object, return dict of lists of Atom objects:
      - 'all': every atom in standard residues
      - 'backbone': N, CA, C
      - 'loop': any atom in residues of loop_regs
      - 'ec': any atom in residues of ec_regs
      - 'ic': any atom in residues of ic_regs
    """
    groups = {k: [] for k in ['all','backbone','loop','ec','ic']}
    for res in chain:
        if not is_aa(res, standard=True):
            continue
        resnum = res.get_id()[1]
        for atom in res:
            groups['all'].append(atom)
            if atom.get_name() in ('N','CA','C'):
                groups['backbone'].append(atom)
            # check each specialized region
            for name, regs in [('loop', loop_regs), ('ec', ec_regs), ('ic', ic_regs)]:
                if any(start <= resnum <= end for start,end in regs):
                    groups[name].append(atom)
    return groups

def extract_atoms(structure, chain_id):
    """Return the first Model → specified Chain."""
    model = structure[0]
    return model[chain_id]

def main():
    p = argparse.ArgumentParser(description="Compute several RMSDs between a local PDB and an RCSB PDB")
    p.add_argument('--converted', '-c', required=True,
                   help="Your .gro→.pdb converted file")
    p.add_argument('--pdb-id',  '-p', required=True,
                   help="RCSB PDB ID (e.g. 6A94)")
    p.add_argument('--chain',   '-C', default='A',
                   help="Chain identifier (both local & RCSB)")
    p.add_argument('--loops',   default='',
                   help="Comma‐sep ranges for loop RMSD, e.g. 45-55,120-130")
    p.add_argument('--ec',      default='',
                   help="Comma‐sep ranges for extracellular region")
    p.add_argument('--ic',      default='',
                   help="Comma‐sep ranges for intracellular region")
    p.add_argument('--out',     '-o', required=True,
                   help="CSV file to write RMSD results")
    args = p.parse_args()

    # Parse region definitions
    loops = parse_ranges(args.loops) if args.loops else []
    ec    = parse_ranges(args.ec)    if args.ec else []
    ic    = parse_ranges(args.ic)    if args.ic else []

    # 1) Fetch reference from RCSB
    pdb_fetcher = PDBList()
    ref_fn = pdb_fetcher.retrieve_pdb_file(args.pdb_id, pdir='.', file_format='pdb')
    # Bio.PDBList names file like 'pdbXXXX.ent'
    
    # 2) Parse both structures
    parser = PDBParser(QUIET=True)
    ref_struct = parser.get_structure('ref', ref_fn)
    conv_struct = parser.get_structure('conv', args.converted)

    # 3) Extract the chain of interest
    ref_chain = extract_atoms(ref_struct, args.chain)
    cv_chain  = extract_atoms(conv_struct, args.chain)

    # 4) Build atom‐groups
    ref_atoms = atom_groups(ref_chain, loops, ec, ic)
    cv_atoms  = atom_groups(cv_chain, loops, ec, ic)

    # 5) For each group, find matching atom pairs by (resnum, atom name)
    results = {}
    for grp in ref_atoms:
        # map (resid, atomname) → Atom
        ref_map = { (atom.get_parent().get_id()[1], atom.get_name()): atom
                    for atom in ref_atoms[grp] }
        cv_map  = { (atom.get_parent().get_id()[1], atom.get_name()): atom
                    for atom in cv_atoms[grp] }

        # intersection of keys
        keys = set(ref_map) & set(cv_map)
        if not keys:
            results[grp] = None
            continue

        # lists of coords
        ref_coords = [ ref_map[k].get_coord() for k in keys ]
        cv_coords  = [ cv_map[k].get_coord()  for k in keys  ]

        # superimpose and compute RMSD
        sup = Superimposer()
        sup.set_atoms([ ref_map[k] for k in keys ],
                      [ cv_map[k]  for k in keys ])
        sup.apply(cv_chain.get_atoms())  # align in‐place
        results[grp] = sup.rms

    # 6) Write out
    with open(args.out, 'w', newline='') as csvf:
        w = csv.writer(csvf)
        w.writerow(['group','n_atoms','rmsd'])
        for grp,rms in results.items():
            n = len(set((a.get_parent().get_id()[1],a.get_name()) for a in ref_atoms[grp]))
            w.writerow([grp, n, f"{rms:.3f}" if rms is not None else 'NA'])

    print(f"Written RMSD breakdown to {args.out}")

if __name__ == '__main__':
    main()

