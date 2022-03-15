This repository contains data and script for reproducing results in the follow manuscript:

*The Effect of Genome Graph Expressiveness onthe Discrepancy Between Genome Graph Distanceand String Set Distance*

`plot_GTED_EMED.ipynb` contains necessary scripts to analyze results presented in the manuscript

## T-Cell Receptor
- Data: MSA graph constructed on 50 string sets are contained in `TCR_graphs`
    - `TCR_all_EMED.csv` -- contains all pair-wise EMED
    - `TCR_sampled_diameters.txt` -- contains all sampled diameters
    - `TCR_fgted_logs` -- contains all output from gurobi for solving FGTED
    - `TCR_summary.csv` -- contains summary stats for all pairs
- Scripts: see `pipeline_tcr.sh` for scripts to produce pair-wise FGTED
    - see `GTED_MSA_TCR_solver.py` for details

## Hepatitis B Virus
- Data: MSA graphs constructed on 9 string sets are contained in `HBV_graphs`
    - `HBV_all_EMED.csv` -- contains all pair-wise EMED
    - `HBV_sampled_diameters.txt` -- contains all sampled diameters
    - `HBV_fgted_logs` -- contains all output from gurobi for solving FGTED
    - `HBV_summary.csv` -- contains all summary stats for all pairs
- Scripts: see `pipeline_hbv.sh` for scripts to produce pair-wise FGTED
    - see `GTED_MSA_HBV_solver.py` for details

## To compute diameters
- See pipeline_diameter.sh for how to run `./scripts/diameter_flow_decomp.py`