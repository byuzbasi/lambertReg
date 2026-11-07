# Run the Lambert penalty study

These workflows accompany *The Lambert Penalty: Logarithmic Shrinkage for
Sparse Regression* ([preprint](https://arxiv.org/abs/2610.09627)).
The directory contains run code and the licensed real-data inputs, with hashes
and attribution in [DATA.md](DATA.md). Simulation samples are generated at run
time. Outputs are written locally and ignored by Git.

Clone the repository and run the commands from the indicated directories.
Use R ≥ 4.1 and an R-compatible C++ compiler. Install digest, jsonlite, knitr,
Rcpp, RcppArmadillo and readxl. The benchmark paths use glmnet 5.0 and
ncvreg 3.16.0. The simulation requires mltools 0.3.5. Install lambertReg 0.1.11
from the repository for the maize run.

## Simulation

Run on macOS or Linux with a worker count that fits available memory.
Each scenario/replication is validated and committed as it finishes.

```sh
cd lambertReg/reproduction/simulation
Rscript --vanilla run.R design
Rscript --vanilla run.R smoke runs/smoke 2
Rscript --vanilla run.R start runs/manuscript 6
# If interrupted:
Rscript --vanilla run.R resume runs/manuscript 6
Rscript --vanilla run.R verify runs/manuscript
```

The full study generates 15,600 datasets in 52 scenarios (300 replications
each). The CV results are in `runs/manuscript/results.csv`. The paired
16-scenario subset uses the same datasets and full-training paths for
BIC, HBIC and EBIC selection. The known-support least-squares reference
is calculated for every scenario. The runner writes `-bic.csv`,
`-information_criteria.csv` and `-oracle.csv` beside the fit directory.
The fit and selection engines follow the archived study, with a portable
run wrapper and output layout.

## QSAR aquatic toxicity

```sh
cd lambertReg/reproduction/qsar
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
Rscript --vanilla scripts/prepare_qsar_aquatic_v1.R validate
Rscript --vanilla scripts/lambert_fixed_c1_realdata_v2.R --dataset qsar --mode smoke --out runs/smoke
Rscript --vanilla scripts/qsar_information_criteria_v1.R smoke
Rscript --vanilla scripts/lambert_fixed_c1_realdata_v2.R --dataset qsar --mode production --out runs/cv
Rscript --vanilla scripts/qsar_information_criteria_v1.R run
```

The study evaluates 100 repeated outer five-fold partitions and tunes within
the corresponding training data. BIC and other information criteria reuse
the same partitions and verify the completed CV input before running. To
resume CV, add `--resume` to its production command. To resume information
criteria, replace its final `run` argument with `resume`. The output folders
are `runs/cv` and `runs/information-criteria`.

## MAGIC maize plant height

```sh
cd lambertReg
R CMD INSTALL .
cd reproduction/magic-maize
Rscript --vanilla run.R design
Rscript --vanilla run.R smoke runs/smoke
Rscript --vanilla run.R start runs/manuscript
# If interrupted:
Rscript --vanilla run.R resume runs/manuscript
```

The study uses 50 frozen random 423/106 training/test splits and five-fold
CV within each training sample. SNP filtering, imputation and scaling are
learned in each training role. The portable runner processes one split at
a time and writes a resumable split checkpoint, held-out predictions and
mean/SD summaries by method and tuning criterion. The smoke run uses a
small subset of rows and markers for an execution check.

The full studies are computationally intensive. Runtime depends on the
machine and numerical libraries. The package gallery is a separate synthetic
demonstration.
