# Data dictionary

## RNA: `dataset_rpe.csv`

Each row represents one RPE1 cell. The ten gene columns are numeric expression
values on the published preprocessed scale. `interventions` is a categorical
text label identifying the targeted gene, `non-targeting` (control), or `excluded`.

| Variable | Role |
|---|---|
| `ENSG00000173812` | Response gene |
| `ENSG00000187514` | Predictor gene |
| `ENSG00000075624` | Predictor gene |
| `ENSG00000147604` | Predictor gene |
| `ENSG00000110700` | Predictor gene |
| `ENSG00000172757` | Predictor gene |
| `ENSG00000133112` | Predictor gene |
| `ENSG00000067225` | Predictor gene |
| `ENSG00000108518` | Predictor gene |
| `ENSG00000125691` | Predictor gene |
| `interventions` | Experimental environment |

The two training groups are `non-targeting` and pooled perturbations targeting
`ENSG00000133112`, `ENSG00000075624`, `ENSG00000187514`, and `ENSG00000067225`.
The analysis excludes the `excluded` label and environments with fewer than
100 cells. Expression values are centered using the training data.

## Chamber: `lt_interventions_standard_v1/`

Each row is one measurement; the CSV filename identifies the environment.
All variables listed below are numeric. Other CSV columns are not used by
this analysis. Original definitions are given in
[Gamella, Peters, and Bühlmann, Appendix II](https://arxiv.org/html/2404.11341v2).

| Variables | Meaning | Role | Units / scale in the CSV |
|---|---|---|---|
| `ir_2` | Sensor-2 infrared intensity | Response | Uncalibrated sensor counts, 0–65535 |
| `red`, `green`, `blue` | Main light-source brightness settings | Predictors | Controller levels, 0–255 |
| `current` | Light-source current measurement | Predictor | Uncalibrated sensor counts, 0–1023 |
| `pol_1`, `pol_2` | Polarizer motor angle settings | Predictors | Degrees |
| `ir_1`, `ir_3` | Infrared intensity measurements | Predictors | Uncalibrated sensor counts, 0–65535 |
| `vis_1`, `vis_2`, `vis_3` | Visible-light intensity measurements | Predictors | Uncalibrated sensor counts, 0–65535 |
| `l_11`, `l_12`, `l_21`, `l_22`, `l_31`, `l_32` | Local LED brightness settings | Predictors | Controller levels, 0–255 |

The classification reference marks `red`, `green`, `blue`, `l_21`, and `l_22`
as direct causes of `ir_2`; the other twelve predictors are non-direct causes.
Predictors and spline columns are standardized, and the response is centered.
