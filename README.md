# gtedemedtest
=======
This repository contains data and script for reproducing results in the follow manuscript:

*The Effect of Genome Graph Expressiveness onthe Discrepancy Between Genome Graph Distanceand String Set Distance*

`plot_GTED_EMED.ipynb` contains necessary scripts to analyze results presented in the manuscript
`scripts/` contains necessary scripts to
1. Sample flow decomposition from input graphs
2. Compute FGTED between graphs
3. Construct MSA/dBG graphs with conforming formats to FGTED solver

## How to compute FGTED / sample diameters from different data
- See `pipeline_tcr.sh`, `pipeline_hbv.sh` for scripts to produce pair-wise FGTED.
- See `pipeline_diameter.sh` for scripts to produce sampled diameters.

## Data
### T-Cell Receptor
- Data: MSA graph and dBG4 constructed on 50 string sets are contained in `TCR_graphs`
    - The instruction to simulate TCR repertoires is included in `data/`
    - `TCR_all_EMED.csv` -- contains all pair-wise EMED
    - `TCR_sampled_diameters_msa.csv` -- contains all sampled diameters on MSA graphs
    - `TCR_sampled_diameters_dbg4.csv` -- contains all sampled diameters on dBG4 graphs
    - `TCR_fgted_logs` -- contains all output from gurobi for solving FGTED
    - `TCR_summary.csv` -- contains summary stats for all pairs of MSA graphs
    - `TCR_FGTED.csv` -- contains FGTED between pairs of dBG4 graphs

## Hepatitis B Virus
- Data: MSA graphs constructed on 9 string sets are contained in `HBV_graphs`
    - `HBV_all_EMED.csv` -- contains all pair-wise EMED
    - `HBV_sampled_diameters.txt` -- contains all sampled diameters
    - `HBV_fgted_logs` -- contains all output from gurobi for solving FGTED
    - `HBV_summary.csv` -- contains all summary stats for all pairs
- Scripts: see `pipeline_hbv.sh` for scripts to produce pair-wise FGTED
    - see `GTED_MSA_HBV_solver.py` for details
