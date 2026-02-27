# Fitting a PLSR model on the (Deviance) Residuals

This function computes the PLSR model with the Residuals of a Cox-Model
fitted with an intercept as the only explanatory variable as the
response and Xplan as explanatory variables. Default behaviour uses the
Deviance residuals. It uses the package `pls`.

## Usage

``` r
coxpls2DR(Xplan, ...)

# Default S3 method
coxpls2DR(
  Xplan,
  time,
  time2,
  event,
  type,
  origin,
  typeres = "deviance",
  collapse,
  weighted,
  scaleX = TRUE,
  scaleY = TRUE,
  ncomp = min(7, ncol(Xplan)),
  methodpls = "kernelpls",
  validation = "CV",
  plot = FALSE,
  allres = FALSE,
  ...
)

# S3 method for class 'formula'
coxpls2DR(
  Xplan,
  time,
  time2,
  event,
  type,
  origin,
  typeres = "deviance",
  collapse,
  weighted,
  scaleX = TRUE,
  scaleY = TRUE,
  ncomp = min(7, ncol(Xplan)),
  methodpls = "kernelpls",
  validation = "CV",
  plot = FALSE,
  allres = FALSE,
  dataXplan = NULL,
  subset,
  weights,
  model_frame = FALSE,
  model_matrix = FALSE,
  contrasts.arg = NULL,
  ...
)
```

## Arguments

- Xplan:

  a formula or a matrix with the eXplanatory variables (training)
  dataset

- ...:

  Arguments to be passed on to
  [`survival::coxph`](https://rdrr.io/pkg/survival/man/coxph.html).

- time:

  for right censored data, this is the follow up time. For interval
  data, the first argument is the starting time for the interval.

- time2:

  The status indicator, normally 0=alive, 1=dead. Other choices are
  `TRUE/FALSE` (`TRUE` = death) or 1/2 (2=death). For interval censored
  data, the status indicator is 0=right censored, 1=event at `time`,
  2=left censored, 3=interval censored. Although unusual, the event
  indicator can be omitted, in which case all subjects are assumed to
  have an event.

- event:

  ending time of the interval for interval censored or counting process
  data only. Intervals are assumed to be open on the left and closed on
  the right, `(start, end]`. For counting process data, event indicates
  whether an event occurred at the end of the interval.

- type:

  character string specifying the type of censoring. Possible values are
  `"right"`, `"left"`, `"counting"`, `"interval"`, or `"interval2"`. The
  default is `"right"` or `"counting"` depending on whether the `time2`
  argument is absent or present, respectively.

- origin:

  for counting process data, the hazard function origin. This option was
  intended to be used in conjunction with a model containing time
  dependent strata in order to align the subjects properly when they
  cross over from one strata to another, but it has rarely proven
  useful.

- typeres:

  character string indicating the type of residual desired. Possible
  values are `"martingale"`, `"deviance"`, `"score"`, `"schoenfeld"`,
  `"dfbeta"`, `"dfbetas"`, and `"scaledsch"`. Only enough of the string
  to determine a unique match is required.

- collapse:

  vector indicating which rows to collapse (sum) over. In time-dependent
  models more than one row data can pertain to a single individual. If
  there were 4 individuals represented by 3, 1, 2 and 4 rows of data
  respectively, then `collapse=c(1,1,1,2,3,3,4,4,4,4)` could be used to
  obtain per subject rather than per observation residuals.

- weighted:

  if `TRUE` and the model was fit with case weights, then the weighted
  residuals are returned.

- scaleX:

  Should the `Xplan` columns be standardized ?

- scaleY:

  Should the `time` values be standardized ?

- ncomp:

  The number of components to include in the model. The number of
  components to fit is specified with the argument ncomp. It this is not
  supplied, the maximal number of components is used (taking account of
  any cross-validation).

- methodpls:

  The multivariate regression method to be used. See
  [`mvrCv`](https://rdrr.io/pkg/pls/man/mvrCv.html) for details.

- validation:

  character. What kind of (internal) validation to use. If
  `validation = "CV"`, cross-validation is performed. The number and
  type of cross-validation segments are specified with the arguments
  `segments` and `segment.type`. See
  [`mvrCv`](https://rdrr.io/pkg/pls/man/mvrCv.html) for details. If
  `validation = "LOO"`, leave-one-out cross-validation is performed. It
  is an error to specify the segments when `validation = "LOO"` is
  specified.

- plot:

  Should the survival function be plotted ?)

- allres:

  FALSE to return only the Cox model and TRUE for additionnal results.
  See details. Defaults to FALSE.

- dataXplan:

  an optional data frame, list or environment (or object coercible by
  [`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html) to a data
  frame) containing the variables in the model. If not found in
  `dataXplan`, the variables are taken from `environment(Xplan)`,
  typically the environment from which `coxpls2DR` is called.

- subset:

  an optional vector specifying a subset of observations to be used in
  the fitting process.

- weights:

  an optional vector of 'prior weights' to be used in the fitting
  process. Should be `NULL` or a numeric vector.

- model_frame:

  If `TRUE`, the model frame is returned.

- model_matrix:

  If `TRUE`, the model matrix is returned.

- contrasts.arg:

  a list, whose entries are values (numeric matrices, functions or
  character strings naming functions) to be used as replacement values
  for the contrasts replacement function and whose names are the names
  of columns of data containing factors.

## Value

If `allres=FALSE` :

- cox_pls2DR:

  Final Cox-model.

If `allres=TRUE` :

- tt_pls2DR:

  PLSR components.

- cox_pls2DR:

  Final Cox-model.

- pls2DR_mod:

  The PLSR model.

## Details

If `allres=FALSE` returns only the final Cox-model. If `allres=TRUE`
returns a list with the PLS components, the final Cox-model and the PLSR
model. `allres=TRUE` is useful for evluating model prediction accuracy
on a test sample.

## References

plsRcox, Cox-Models in a high dimensional setting in R, Frederic
Bertrand, Philippe Bastien, Nicolas Meyer and Myriam Maumy-Bertrand
(2014). Proceedings of User2014!, Los Angeles, page 152.  

Deviance residuals-based sparse PLS and sparse kernel PLS regression for
censored data, Philippe Bastien, Frederic Bertrand, Nicolas Meyer and
Myriam Maumy-Bertrand (2015), Bioinformatics, 31(3):397-404,
doi:10.1093/bioinformatics/btu660.

## See also

[`coxph`](https://rdrr.io/pkg/survival/man/coxph.html),
[`plsr`](https://rdrr.io/pkg/pls/man/mvr.html)

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
data(micro.censure)
data(Xmicro.censure_compl_imp)

X_train_micro <- apply((as.matrix(Xmicro.censure_compl_imp)),FUN="as.numeric",MARGIN=2)[1:80,]
X_train_micro_df <- data.frame(X_train_micro)
Y_train_micro <- micro.censure$survyear[1:80]
C_train_micro <- micro.censure$DC[1:80]

(cox_pls2DR_fit=coxpls2DR(X_train_micro,Y_train_micro,C_train_micro,ncomp=6,validation="none"))
#> Call:
#> coxph(formula = YCsurv ~ ., data = tt_pls2DR)
#> 
#>          coef exp(coef) se(coef)     z        p
#> Comp.1 0.7784    2.1781   0.1987 3.917 8.96e-05
#> Comp.2 0.9626    2.6186   0.2982 3.228  0.00125
#> Comp.3 0.9110    2.4868   0.4075 2.236  0.02536
#> Comp.4 0.9022    2.4650   0.4004 2.253  0.02424
#> Comp.5 0.1844    1.2026   0.2664 0.692  0.48865
#> Comp.6 0.7448    2.1059   0.4228 1.761  0.07819
#> 
#> Likelihood ratio test=54.95  on 6 df, p=4.745e-10
#> n= 80, number of events= 17 
(cox_pls2DR_fit2=coxpls2DR(~X_train_micro,Y_train_micro,C_train_micro,ncomp=6,validation="none"))
#> Call:
#> coxph(formula = YCsurv ~ ., data = tt_pls2DR)
#> 
#>          coef exp(coef) se(coef)     z        p
#> Comp.1 0.7784    2.1781   0.1987 3.917 8.96e-05
#> Comp.2 0.9626    2.6186   0.2982 3.228  0.00125
#> Comp.3 0.9110    2.4868   0.4075 2.236  0.02536
#> Comp.4 0.9022    2.4650   0.4004 2.253  0.02424
#> Comp.5 0.1844    1.2026   0.2664 0.692  0.48865
#> Comp.6 0.7448    2.1059   0.4228 1.761  0.07819
#> 
#> Likelihood ratio test=54.95  on 6 df, p=4.745e-10
#> n= 80, number of events= 17 
(cox_pls2DR_fit3=coxpls2DR(~.,Y_train_micro,C_train_micro,ncomp=6,validation="none",
dataXplan=X_train_micro_df))
#> Call:
#> coxph(formula = YCsurv ~ ., data = tt_pls2DR)
#> 
#>          coef exp(coef) se(coef)     z        p
#> Comp.1 0.7784    2.1781   0.1987 3.917 8.96e-05
#> Comp.2 0.9626    2.6186   0.2982 3.228  0.00125
#> Comp.3 0.9110    2.4868   0.4075 2.236  0.02536
#> Comp.4 0.9022    2.4650   0.4004 2.253  0.02424
#> Comp.5 0.1844    1.2026   0.2664 0.692  0.48865
#> Comp.6 0.7448    2.1059   0.4228 1.761  0.07819
#> 
#> Likelihood ratio test=54.95  on 6 df, p=4.745e-10
#> n= 80, number of events= 17 

rm(X_train_micro,Y_train_micro,C_train_micro,cox_pls2DR_fit,cox_pls2DR_fit2,cox_pls2DR_fit3)
```
