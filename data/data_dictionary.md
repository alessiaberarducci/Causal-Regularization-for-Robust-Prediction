# Analysis data dictionary

## RNA table

One row is one RPE1 cell. The numeric gene columns are inherited prepared
expression values on the published table's scale, rather than physical SI units
or raw integer sequencing counts. All ten columns are numeric.

| Variable | Role |
|---|---|
| `ENSG00000173812` | Response gene |
| `ENSG00000187514` | Predictor gene |
| `ENSG00000075624` | Predictor gene; pooled intervention training group |
| `ENSG00000147604` | Predictor gene |
| `ENSG00000110700` | Predictor gene |
| `ENSG00000172757` | Predictor gene |
| `ENSG00000133112` | Predictor gene; pooled intervention training group |
| `ENSG00000067225` | Predictor gene; pooled intervention training group |
| `ENSG00000108518` | Predictor gene |
| `ENSG00000125691` | Predictor gene |
| `interventions` | Character environment label: targeted Ensembl gene ID, `non-targeting`, or `excluded` |

The pooled intervention training group also includes `ENSG00000187514`.
The comparison group is `non-targeting`. Cells are nested in intervention
environments; tests are computed separately for each eligible unseen label.

## Chamber tables

One row is one measurement; the CSV filename identifies the experimental
environment. All analyzed columns are numeric. Values are retained in native
instrument/control units. `chamber_variables_upstream.csv` defines all raw
columns, including metadata and hardware diagnostics excluded from analysis.
The source paper's Appendix II describes calibration and measurement units.

| Variables | Meaning | Analysis role / scale |
|---|---|---|
| `ir_2` | Infrared sensor-2 intensity | Response, native sensor reading |
| `red`, `green`, `blue` | Main light-source color settings | Predictors and direct causes; controller intensity levels |
| `current` | Light-source current measurement | Predictor, native sensor reading |
| `pol_1`, `pol_2` | Polarizer angle settings | Predictors, degrees |
| `ir_1`, `ir_3` | Infrared sensor intensities | Predictors, native sensor readings |
| `vis_1`, `vis_2`, `vis_3` | Visible-light sensor intensities | Predictors, native sensor readings |
| `l_11`, `l_12`, `l_21`, `l_22`, `l_31`, `l_32` | Local LED brightness settings | Predictors, controller intensity levels |

The known direct causes of `ir_2` are `red`, `green`, `blue`, `l_21`, `l_22`.
The remaining 12 predictors are treated as non-direct causes for the
classification metrics. Standardized predictor and spline-basis values are
dimensionless. Each predictor yields three quadratic-spline coefficients at
the single knot zero. Predictions are on the centered response scale.
