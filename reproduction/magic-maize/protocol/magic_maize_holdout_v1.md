# MAGIC maize: repeated holdout design

Retain all 529 source PH responses and 54,234 original dosages. Each split
samples 423 training rows without replacement and leaves 106 test rows.
For split r=1,...,50, training seed is 2026101801+1000*r; inner-fold seed is
2026101802+1000*r. Use the unchanged sch_sample Mersenne-Twister/Inversion/
Rejection helper and source-sorted training IDs. Inner labels are shuffled
balanced 1:5. The pilot seeds are excluded, so the pilot is separate from the
50-result summary. Persist seeds, IDs and inner labels before any fitting.
Legacy field outer_fold is always 1 and denotes the one held-out test group;
there is no outer K-fold loop. Do not pool the 50 test samples as independent
people or interpret split SD as a population standard error.

Use mmp_fit unchanged on each full training sample and its five inner-training
subsets. All methods share deterministic transformations for each role:
call rate .95, MAF .05, observed-pattern duplicate removal, LD r-squared .4
within 500 kb, observed-training-mean filling and training RMS scaling. No
whole-cohort reference matrix, genotype filling from test/validation, response
transformation, new exclusions or different tolerance is introduced. SNP sets
may differ by role/split. Cache at most the six roles of one split per worker;
clear the cache on a split change. Learn contexts on demand from raw dosages.

Keep Lambert c=1, MCP=3, SCAD=3.7; LASSO, Ridge and the existing ENET alpha
grid {0,.25,.5,.75,1}. Keep 60 lambda positions, normalized-fraction CV,
tie/start rules, early_stop=FALSE, 10,000-sweep rescue and numerical gates:
normalized KKT 1e-7, support 1e-8, objective increase 1e-9, curvature 1e-8.
Reuse lambertReg 0.1.11 and existing solver/curvature/IC functions. CV and BIC cover six methods;
EBIC(.5), EBIC(1), HBIC cover five sparse methods; Ridge stays CV/BIC only.

