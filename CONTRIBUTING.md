# Contributing to lambertReg

Thank you for helping improve lambertReg. Start with a GitHub issue describing
the proposed change or a small reproducible example of the problem.

## A useful report

Include the package and R versions, operating system, function call, controls,
and the numerical status of the fit. Use a synthetic or openly licensed example
and omit confidential data. Include a seed and `foldid` when the issue involves
cross-validation. For plotting issues, include the plotting call and device size.

## Scope

Documentation, examples, accessibility, and interface improvements are welcome.
Changes to the estimator, preprocessing, scalar update, coordinate order, or
acceptance criteria require an explicit methodological discussion. A faster
implementation must preserve the relevant coefficients, CV scores, and diagnostics
against the reference implementation; runtime alone is not sufficient evidence.

## Local checks

From a clone's parent directory, run:

```sh
R CMD build lambertReg
R CMD check lambertReg_0.1.10.tar.gz
```

Use the tests in `tests/` and small examples. Do not attach large datasets or
launch a production study as part of a routine contribution. Package examples
use dense finite predictor matrices and Gaussian responses; preserve the current
function signatures and defaults unless the change is agreed in advance.

## Documentation site

The optional site build uses pkgdown, knitr, rmarkdown, and Pandoc. In R,
from the repository directory:

```r
pkgdown::build_site(preview = FALSE)
```

The local site is written to `docs/`. Inspect the introduction, reference pages,
and example plots on both desktop and a narrow screen before proposing changes.
