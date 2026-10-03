# Generalized Causal Regularization

Code, data, and results for the computational experiments in
*Generalized Causal Regularization*.

<!--
# Generalized Causal Regularization

Code and data for reproducing the computational experiments in
*Generalized Causal Regularization*.

## Requirements

R 4.6.1 and its standard packages. No additional R packages are required.
All input data are included.

## Reproduce the experiments

Open a terminal in the repository folder and run:

```sh
Rscript --vanilla code/run_all.R
```

This command runs the simulations, RNA analysis, and Causal Chamber analysis,
and saves figures and numerical results in `output/`. 
Each analysis runs in a separate R session.

To run one part:

```sh
Rscript --vanilla code/run_all.R simulations
Rscript --vanilla code/run_all.R rna
Rscript --vanilla code/run_all.R chamber
```

Individual entry scripts can also be run with `Rscript --vanilla` from the
repository folder. Execution logs are saved locally in `output/logs/`.

## Repository structure

- `code/`: experiment scripts and shared functions.
- `data/`: input datasets, source information, and variable descriptions.
- `output/`: figures and numerical results.

## Figures and table

Output paths below are relative to `output/`.

| Result | Entry script | Output |
|---|---|---|
| Figure 3 | `code/simulations/example_A.R` | `simulations/figure_03_a.pdf`, `simulations/figure_03_b.pdf` |
| Figure 5 | `code/simulations/example_B.R` | `simulations/figure_05_a.pdf`, `simulations/figure_05_b.pdf` |
| Figure 7 | `code/simulations/example_C.R` | `simulations/figure_07_a.pdf`, `simulations/figure_07_b.pdf` |
| Figure 9 | `code/simulations/example_D.R` | `simulations/figure_09_a.pdf`, `simulations/figure_09_b.pdf` |
| Figure 10 | `code/simulations/anchor_comparison.R` | `simulations/figure_10.pdf` |
| Figure 11(a) | `code/RNA/paper_RNA.R` | `rna/rpe_subsample_OLS_CV_scatterplot_log_full_sample_alternative_n2.pdf` |
| Figure 11(b) | `code/RNA/paper_RNA.R` | `rna/rpe_subsample_OLS_CV_scatterplot_log_n100_alternative_n2.pdf` |
| Figure 12(a) | `code/RNA/paper_RNA.R` | `rna/rpe_OLS_minus_CR_MSE_histogram_full_sample_alternative_n2.pdf` |
| Figure 12(b) | `code/RNA/paper_RNA.R` | `rna/rpe_OLS_minus_CR_MSE_histogram_n100_alternative_n2.pdf` |
| Table 1 | `code/chamber/paper_chamber.R` | `table_1.csv` |
| Tests accompanying Table 1 | `code/chamber/paper_chamber.R` | `table_1_tests.csv` |

Figures 1, 2, 4, 6, 8, and 13 are conceptual diagrams rather than computed experiments.
CSV files and RDS objects contain coefficient paths, cross-validation results,
and intermediate fitted objects. Figure 10 uses the same observations and
causal regularization path as Figure 9.

## Data sources

RNA data are the processed RPE1 measurements used by Shen, Bühlmann, and Taeb in
[Causality-oriented robustness: exploiting general noise interventions](https://arxiv.org/abs/2307.10299).
Chamber data are the `lt_interventions_standard_v1` dataset from
[Causal Chamber](https://github.com/juangamella/causal-chamber).

See [data/README.md](data/README.md) for source versions, preprocessing,
and download instructions, and
[data/data_dictionary.md](data/data_dictionary.md) for variable descriptions.
-->
