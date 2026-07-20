# Simulation Workflow

This folder contains the simulation workflow used to compare MiXcan, S-MiXcan,
and PrediXcan.

The current workflow uses binary disease outcomes and the real haplotype pool
stored in Dropbox. Continuous-outcome simulation scripts and older sensitivity
scripts were moved to the archive because they are not part of the current
figure workflow.

## Input

The simulation uses this haplotype pool:

```text
/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Data/Simulation/X_pool_filtered.csv
```

By default, the scripts read this file through:

```r
PAPER_SMIXCAN_DIR/Data/Simulation/X_pool_filtered.csv
```

You can override it with:

```bash
SMIXCAN_SIM_X_POOL=/path/to/X_pool_filtered.csv
```

## Current Simulation Setting

The simulated model is:

```text
y1 = b0 + X b1 + e1
y2 =      X b2 + e2
y  = pi * y1 + (1 - pi) * y2
logit P(D = 1) = eta0 + eta1 * y1 + eta2 * y2
```

The cell-type proportion is generated as:

```text
pi ~ Beta(2, 3)
```

So `y1` is the minor cell type on average and `y2` is the major cell type on
average.

The current regularization setting is:

```text
fixed reg_scale = 0.05
```

For the all-panel simulation:

```text
Type I error scenario: 2000 replicates
Power scenarios:       200 replicates
```

## Scripts

### Step 0: Haplotype Sampling Helper

```text
0_hap_generation.R
```

Helper functions for sampling haplotypes and generating genotype matrices from
the real haplotype pool.

### Step 1: Binary Data Generation

```text
1_data_generation_binary.R
```

Generates simulated training and application data for binary disease outcomes.

### Step 2: Shared Simulation Runner

```text
2_run_sim.R
```

Trains prediction models and runs association testing for:

```text
MiXcan
S-MiXcan
PrediXcan
```

### Step 3: Run All-Panel Simulation

```bash
Rscript Analysis_code/Simulation/3_run_all_panel_simulation.R
```

Default output:

```text
Results/simulation/simulation_figure_with_predixcan_latest_reps/
```

Main outputs:

```text
condition_results/*.csv
simulation_full_results_with_predixcan.csv
simulation_summary_with_predixcan.csv
```

This step is expensive. The full all-panel simulation is approximately:

```text
6 panels * 5 x-values * (2000 + 4 * 200) = 84,000 replicate simulations
```

Expected runtime on the current machine is roughly several hours.

### Step 4: Plot All-Panel Figure

```bash
Rscript Analysis_code/Simulation/4_plot_all_panel_simulation.R
```

Default input:

```text
Results/simulation/simulation_figure_with_predixcan_latest_reps/
```

Default output:

```text
simulation_figure_with_predixcan.pdf
simulation_figure_with_predixcan.png
```

To plot separate association p-values:

```bash
SMIXCAN_SIM_PVALUE_MODE=sep Rscript Analysis_code/Simulation/4_plot_all_panel_simulation.R
```

### Step 5: Run Type I Error for ABC Figures

```bash
Rscript Analysis_code/Simulation/5_run_type1_error_2000rep.R
```

Default setting:

```text
b0 = 1
b1 = 0.5
b2 = 1
group = heter xy
replicates = 2000
```

This script is used when regenerating the Type I error component for the ABC
summary figures.

### Step 6: Run Power for ABC Figures

```bash
Rscript Analysis_code/Simulation/6_run_power_eta4_200rep.R
```

Default setting:

```text
b0 = 1
b1 = 0.5
b2 = 1
group = heter xy
replicates = 200 per disease scenario
```

Power scenarios:

```text
eta1 = 0.2,  eta2 = 0.2
eta1 = 0.2,  eta2 = 0
eta1 = 0,    eta2 = 0.2
eta1 = -0.2, eta2 = 0.2
```

### Step 7: Plot ABC Summary Figures

Homogeneous SNP-Exp setting:

```bash
Rscript Analysis_code/Simulation/7_plot_homogeneous_abc_summary.R
```

Output:

```text
Results/simulation/final_bar_summary_homogeneous/
```

Heterogeneous SNP-Exp setting:

```bash
Rscript Analysis_code/Simulation/7_plot_heterogeneous_abc_summary.R
```

Output:

```text
Results/simulation/final_bar_summary_heterogeneous/
```

Each final folder contains:

```text
*.pdf
*.png
*_bar_summary_values.csv
replicate_results/
```

The `replicate_results/` folder stores the full replicate-level tables used to
produce the final ABC figures.

## Final Results Currently Kept

Only the latest homogeneous and heterogeneous ABC results are kept in the main
simulation results folder:

```text
Results/simulation/final_bar_summary_homogeneous/
Results/simulation/final_bar_summary_heterogeneous/
```

Older simulation outputs were moved to:

```text
Archive/Simulation_results_unused_20260720/
```

Older unused simulation scripts were moved to:

```text
Archive/Github/Simulation_unused_20260720/
```

## Useful Overrides

Run a smaller smoke test:

```bash
SMIXCAN_SIM_PANELS=a \
SMIXCAN_SIM_SCENARIOS=S1,S2 \
SMIXCAN_SIM_BATCH_SIZE=10 \
SMIXCAN_SIM_ITERATIONS=1 \
SMIXCAN_SIM_TYPE1_ITERATIONS=1 \
Rscript Analysis_code/Simulation/3_run_all_panel_simulation.R
```

Use LD from the application/test genotypes:

```bash
SMIXCAN_SIM_LD_REFERENCE=test Rscript Analysis_code/Simulation/3_run_all_panel_simulation.R
```

Use more workers:

```bash
SMIXCAN_SIM_WORKERS=4 Rscript Analysis_code/Simulation/3_run_all_panel_simulation.R
```

