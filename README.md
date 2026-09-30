# Gut Microbiome-Targeted Peptide Design Pipeline

A complete computational pipeline for identifying species-specific surface proteins from human gut
bacteria and designing targeting peptide binders, from reference proteome download through
molecular-dynamics validation.

## Pipeline Overview

```text
32 gut bacterial species (NCBI RefSeq)
        |
        v
[Stage 1] Target screening
   1. Download 32 reference proteomes
   2. All-vs-32 BLASTp comparison (whole-sequence identity)
   3. Grade proteins S / A / B / C by cross-species identity
   4. Final validation against NCBI nr (online BLASTp)
        |
        v
Selected species-specific target proteins
        |
        v
[Stage 2] Binder design and validation
   1. RFdiffusion backbone generation + ProteinMPNN sequence design + AlphaFold2 prediction
   2. Rosetta scoring and InterfaceAnalyzer interface evaluation
   3. GROMACS GPU molecular-dynamics simulation
   4. Trajectory post-processing (nojump, centering, fitting)
   5. MD analysis (RMSD, Rg, RMSF, min distance, H-bonds, temperature, pressure, potential)
   6. Batch plotting of results
        |
        v
Validated targeting peptide binders
```

## Directory Structure

```
.
├── README.md                 # This file
├── LICENSE                   # MIT License
├── .gitignore
├── .gitattributes
├── requirements-plot.txt     # Python dependencies for the plotting script
├── 01_target_screening/      # Stage 1: species-specific target protein screening
│   ├── README.md             # Detailed documentation for screening scripts
│   ├── 01_download_proteomes.sh
│   ├── 02_blast_specificity_v2.py
│   ├── 03_batch_blast_ncbi.py
│   └── 04_rank_specificity_graded.py
└── 02_binder_design/         # Stage 2: binder design, screening and MD validation
    ├── README.md             # Detailed documentation for binder design scripts
    ├── 01_rfdiffusion_binder_design.sh
    ├── 02_rosetta_binder_screening.sh
    ├── 03_gromacs_md_simulation_gpu.sh
    ├── 04_gromacs_trajectory_postprocess.sh
    ├── 05_gromacs_md_analysis.sh
    └── 06_plot_xvg_results.py
```

## Stage 1 — Target Screening (`01_target_screening/`)

| Script | Function |
|---|---|
| `01_download_proteomes.sh` | Download 32 RefSeq proteomes (fixed assembly accessions, resumable) |
| `02_blast_specificity_v2.py` | Each strain vs the combined other-31 BLASTp library; ranked by identity |
| `04_rank_specificity_graded.py` | Grade proteins: S (no hit) / A (<30%) / B (30–40%) / C (>40%, rejected) |
| `03_batch_blast_ncbi.py` | Submit candidates to NCBI nr online BLASTp for final validation (resumable) |

Grading criteria: S — no cross-species hit; A — identity <30%; B — identity 30–40%;
C — identity >40% (rejected by default). See `01_target_screening/README.md` for full usage.

## Stage 2 — Binder Design (`02_binder_design/`)

| Script | Function |
|---|---|
| `01_rfdiffusion_binder_design.sh` | RFdiffusion (2000 designs) + ProteinMPNN + AlphaFold2 initial guess |
| `02_rosetta_binder_screening.sh` | Extract AF2 candidates; Rosetta score + InterfaceAnalyzer |
| `03_gromacs_md_simulation_gpu.sh` | Solvation, ionization, minimization, NVT/NPT, production MD |
| `04_gromacs_trajectory_postprocess.sh` | Nojump, centering, rotational/translational fitting |
| `05_gromacs_md_analysis.sh` | RMSD, Rg, RMSF, min distance, H-bonds and thermodynamic analysis |
| `06_plot_xvg_results.py` | Batch plotting of 8 `.xvg` result files to PNG |

Key parameters: CHARMM27 force field, TIP3P water, dodecahedral box (-d 3.0), 50,000 ps production.
See `02_binder_design/README.md` for the full software environment and output description.

## Requirements

- **Stage 1:** NCBI Datasets CLI (v14+), BLAST+ (2.10+), Python 3.6+ (standard library only)
- **Stage 2:** RFdiffusion v1.1.0, ProteinMPNN, dl_binder_design, Rosetta 2024.09,
  AlphaFold2, GROMACS 2025.0; Python 3 with NumPy/Matplotlib for plotting (`requirements-plot.txt`)

All scripts run on Linux / WSL / macOS. Every long-running script supports checkpoint/resume.

## License

Released under the MIT License. See [LICENSE](LICENSE).
