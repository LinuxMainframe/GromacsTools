# GromacsTools

**GromacsTools** is an open-source collection of Python scripts designed to streamline and document common GROMACS-based workflows in structural biology and molecular dynamics research. Each script encapsulates a specific step—conversion, selection, analysis, or reporting—allowing users to adapt, extend, or repurpose code for their own pipelines. While GromacsTools is not intended as a commercial product or “finished” software package, it provides transparent, reusable building blocks and preserves the provenance of your research steps for reproducibility and publication.

---

## Table of Contents

1. [Highlights](#highlights)  
2. [Installation](#installation)  
3. [Project Structure](#project-structure)  
4. [Getting Started](#getting-started)  
5. [Usage Examples](#usage-examples)  
6. [Contributing & Extensibility](#contributing--extensibility)  
7. [Citation & Acknowledgments](#citation--acknowledgments)  
8. [About the Author](#about-the-author)  
9. [License](#license)  

---

## Highlights

- **Modular Scripts**: Each tool lives in its own Python script, handling one discrete task (e.g., .gro → .pdb conversion, atom selection, RMSD analysis).  
- **Transparency**: Clear inline documentation and headers explain every step—ideal for journal supplementary materials.  
- **Customizable**: Easily adapt residue selections, I/O formats, or add new analysis modules.  
- **Academic Focus**: Preserves a record of methods and software citations to satisfy rigorous publication requirements.

---

## Installation

1. **Clone the repository**  
   ```bash
   git clone https://github.com/LinuxMainframe/GromacsTools.git
   cd GromacsTools
   ```

2. **Set up a Python environment** (recommended):  
   ```bash
   python3 -m venv venv
   source venv/bin/activate
   ```

3. **Install dependencies**  
   ```bash
   pip install --upgrade pip
   pip install -r requirements.txt
   ```

> **Dependencies include**:  
> - [BioPython](https://biopython.org/) (for PDB parsing/selection)  
> - [GROMACS](http://www.gromacs.org/) command-line tools (e.g., `gmx2pdb`, `editconf`)  
> - Others as specified per script (e.g., MDAnalysis, NumPy, Matplotlib).

---

## Project Structure

```text
GromacsTools/
├── README.md
├── LICENSE
├── requirements.txt
├── convert_gro2pdb.py         # .gro → .pdb conversion + selection
├── select_protein_atoms.py    # Atom selection filter for standard AA / palmitoylation
├── compute_rmsd.py            # Automate RMSD calculations using BioPython / MDAnalysis
├── analyze_rmsf.py            # RMSF analysis & plotting
└── docs/
    └── citations.md           # Tool citations & versions
```

---

## Getting Started

1. **Convert a GROMACS `.gro` file to a filtered PDB**  
   ```bash
   python convert_gro2pdb.py      --input molecule.gro      --output molecule_ready.pdb
   ```

2. **Run RMSD calculation**  
   ```bash
   python compute_rmsd.py      --ref reference.pdb      --mobile molecule_ready.pdb      --out rmsd_results.csv
   ```

3. **Plot RMSF profile**  
   ```bash
   python analyze_rmsf.py      --traj production.xtc      --topology molecule_ready.pdb      --out rmsf_plot.png
   ```

---

## Usage Examples

```bash
# 1. Convert pre- and post-production trajectories
python convert_gro2pdb.py -i pre_prod.gro  -o pre_ready.pdb
python convert_gro2pdb.py -i post_prod.gro -o post_ready.pdb

# 2. Compute and compare RMSD
python compute_rmsd.py -r pre_ready.pdb -m post_ready.pdb -o rmsd.csv

# 3. Visualize deviations
python analyze_rmsf.py -t trajectory.xtc -p post_ready.pdb -o rmsf.png
```

---

## Contributing & Extensibility

We welcome improvements, bug fixes, or new analysis modules. To contribute:

1. Fork the repository.  
2. Create a feature branch (`git checkout -b feature/YourFeature`).  
3. Add tests and update `requirements.txt` if needed.  
4. Submit a pull request with a clear description of changes.

Please follow PEP 8 guidelines and include in-code documentation for each function or class.

---

## Citation & Acknowledgments

When using GromacsTools in publications, please cite:

- **BioPython**:  
  Cock, P.J.A. *et al.* Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* **25**, 1422–1423 (2009).  
- **GROMACS**:  
  Abraham, M.J. *et al.* GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers. *SoftwareX* **1–2**, 19–25 (2015).

Additional citations and version details are maintained in `docs/citations.md`.

---

## About the Author

**Aidan A. Bradley**  
- University of Connecticut, Department of Chemistry & Biology  
- GitHub: https://github.com/LinuxMainframe/   
- LinkedIn: https://www.linkedin.com/in/aidanbradleyresearch/

---

## License

Distributed under the [MIT License](LICENSE).  
© 2025 Aidan A. Bradley. All rights reserved—provided “as is” without warranty.  

---

*This file will evolve as new scripts and capabilities are added to GromacsTools.*  
