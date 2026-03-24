# Predict method for DR-based Cox component models

This function provides prediction methods for the rich objects returned
by `coxplsDR(..., allres = TRUE)`, `coxsplsDR(..., allres = TRUE)` and
`coxDKsplsDR(..., allres = TRUE)`.

## Usage

``` r
# S3 method for class 'coxplsDRmodel'
predict(
  object,
  newdata,
  comps = ncol(object$tt_plsDR),
  type = c("lp", "risk", "expected", "terms", "scores"),
  se.fit = FALSE,
  reference = c("strata", "sample", "zero"),
  y = NULL,
  weights = NULL,
  verbose = TRUE,
  ...
)

# S3 method for class 'coxsplsDRmodel'
predict(
  object,
  newdata,
  comps = ncol(object$tt_splsDR),
  type = c("lp", "risk", "expected", "terms", "scores"),
  se.fit = FALSE,
  reference = c("strata", "sample", "zero"),
  y = NULL,
  weights = NULL,
  verbose = TRUE,
  ...
)

# S3 method for class 'coxDKsplsDRmodel'
predict(
  object,
  newdata,
  comps = ncol(object$tt_DKsplsDR),
  type = c("lp", "risk", "expected", "terms", "scores"),
  se.fit = FALSE,
  reference = c("strata", "sample", "zero"),
  y = NULL,
  weights = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  An object returned by one of the DR-based fitting functions with
  `allres = TRUE`.

- newdata:

  An optional data frame or matrix containing original covariates. If
  omitted, predictions are computed on the training data.

- comps:

  Number of latent components to use for prediction.

- type:

  Type of predicted value. Choices are the linear predictor (`"lp"`),
  the risk score exp(lp) (`"risk"`), the expected number of events
  (`"expected"`), the terms of the linear predictor (`"terms"`) or the
  latent component scores (`"scores"`).

- se.fit:

  If `TRUE`, pointwise standard errors are produced by the underlying
  Cox model.

- reference:

  Reference level used to center relative predictions. This is passed to
  [`predict.coxph`](https://rdrr.io/pkg/survival/man/predict.coxph.html)
  and affects `type = "lp"`, `"risk"` and `"terms"`.

- y:

  Optional [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html) response
  to use with `type = "expected"` when predicting on genuinely new
  observations.

- weights:

  Optional case weights used when rebuilding a model matrix for
  formula-based fits.

- verbose:

  Should some details be displayed?

- ...:

  Additional arguments passed to
  [`predict.coxph`](https://rdrr.io/pkg/survival/man/predict.coxph.html).

## Value

A vector, matrix or list of predictions depending on `type` and
`se.fit`.

## See also

[`predict.coxph`](https://rdrr.io/pkg/survival/man/predict.coxph.html),
[`coxplsDR`](https://fbertran.github.io/plsRcox/reference/coxplsDR.md),
[`coxsplsDR`](https://fbertran.github.io/plsRcox/reference/coxsplsDR.md),
[`coxDKsplsDR`](https://fbertran.github.io/plsRcox/reference/coxDKsplsDR.md)
