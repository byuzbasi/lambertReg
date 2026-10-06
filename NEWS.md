# lambertReg 0.1.10

- Use the previously validated column-contiguous C++ kernel, reusing Lambert
  roots for penalty and derivative evaluations with unchanged scalar arithmetic,
  coordinate order and convergence criteria.
- Remove an unused singular-value decomposition from training preprocessing.
- Standardize each validation fold once per CV call rather than at every
  penalty position, using the same training-only statistics.
- Reuse verified warm-start fits only when the warm and zero initial vectors
  are exactly identical. Both logical attempt records remain available;
  `attempts$executed` and `attempts$reused` distinguish actual solver calls.
- Preserve public fitting arguments, the complete coefficient path and CV plot
  data. Early stopping is not enabled.

# lambertReg 0.1.9

* Added an executable introductory vignette covering preprocessing, cross-validation,
  prediction, diagnostics, CV plots, coefficient paths and penalty evaluation.
* Expanded help for returned objects, dimensions, status codes and reproducibility.
* Added private GitHub installation, citation and support information.
* Added full vignette builds to the local package-validation workflow.
* Numerical estimation code, defaults and stored-study results are unchanged.

# lambertReg 0.1.8

* Set CV plot labels, titles and annotation text to black.

# lambertReg 0.1.7

* Introduced blue CV curves, orange selected candidates and compact annotations.

# lambertReg 0.1.6

* Aligned native connected axes with explicit tick limits.

# lambertReg 0.1.4

* Added coefficient-path plots and CV-to-path plot dispatch.

# lambertReg 0.1.1

* Completed the full-training path for CV variable-count displays.

# lambertReg 0.1.0

* Initial fixed-shape Gaussian Lambert regression implementation with
  regularization paths, fold-specific CV, prediction and numerical diagnostics.
