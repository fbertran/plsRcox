# Fitting a sPLSR model on the (Deviance) Residuals

This function computes the Cox Model based on sPLSR components computed
model with

- as the response: the Residuals of a Cox-Model fitted with no covariate

- as explanatory variables: Xplan.

It uses the package `spls` to perform the first step in SPLSR then
`mixOmics` to perform PLSR step fit.

## Usage

``` r
coxsplsDR(Xplan, ...)

# Default S3 method
coxsplsDR(
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
  modepls = "regression",
  plot = FALSE,
  allres = FALSE,
  eta = 0.5,
  trace = FALSE,
  ...
)

# S3 method for class 'formula'
coxsplsDR(
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
  modepls = "regression",
  plot = FALSE,
  allres = FALSE,
  dataXplan = NULL,
  subset,
  weights,
  model_frame = FALSE,
  eta = 0.5,
  trace = FALSE,
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
  supplied, the maximal number of components is used.

- modepls:

  character string. What type of algorithm to use, (partially) matching
  one of "regression", "canonical", "invariant" or "classic". See
  [`pls`](https://rdrr.io/pkg/mixOmics/man/pls.html) for details

- plot:

  Should the survival function be plotted ?)

- allres:

  FALSE to return only the Cox model and TRUE for additionnal results.
  See details. Defaults to FALSE.

- eta:

  Thresholding parameter. `eta` should be between 0 and 1.

- trace:

  Print out the progress of variable selection?

- dataXplan:

  an optional data frame, list or environment (or object coercible by
  [`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html) to a data
  frame) containing the variables in the model. If not found in
  `dataXplan`, the variables are taken from `environment(Xplan)`,
  typically the environment from which `coxsplsDR` is called.

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

- cox_splsDR:

  Final Cox-model.

If `allres=TRUE` :

- tt_splsDR:

  sPLSR components.

- cox_splsDR:

  Final Cox-model.

- splsDR_mod:

  The sPLSR model.

## Details

If `allres=FALSE` returns only the final Cox-model. If `allres=TRUE`
returns a list with the sPLS components, the final Cox-model and the
sPLSR model. `allres=TRUE` is useful for evluating model prediction
accuracy on a test sample.

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

(cox_splsDR_fit=coxsplsDR(X_train_micro,Y_train_micro,C_train_micro,ncomp=6,eta=.5))
#> Call:
#> coxph(formula = YCsurv ~ ., data = tt_splsDR)
#> 
#>         coef exp(coef) se(coef)     z        p
#> dim.1 0.8093    2.2462   0.2029 3.989 6.63e-05
#> dim.2 0.9295    2.5333   0.2939 3.163  0.00156
#> dim.3 0.9968    2.7096   0.4190 2.379  0.01736
#> dim.4 0.9705    2.6391   0.3793 2.558  0.01052
#> dim.5 0.2162    1.2413   0.2811 0.769  0.44192
#> dim.6 0.4380    1.5496   0.3608 1.214  0.22473
#> 
#> Likelihood ratio test=55.06  on 6 df, p=4.51e-10
#> n= 80, number of events= 17 
(cox_splsDR_fit2=coxsplsDR(~X_train_micro,Y_train_micro,C_train_micro,ncomp=6,eta=.5,trace=TRUE))
#> The variables that join the set of selected variables at each step:
#> - 1th step (K=1):
#> X_train_microD20S107 X_train_microD5S346 X_train_microD1S225 X_train_microD3S1282 X_train_microD15S127 X_train_microD1S207 X_train_microD2S138 X_train_microD10S191 X_train_microD14S65 X_train_microD4S414
#> X_train_microD16S408 X_train_microT X_train_microN X_train_microSTADE
#> - 2th step (K=2):
#> X_train_microD22S928 X_train_microD16S422 X_train_microD3S1283 X_train_microAgediag X_train_microM
#> - 3th step (K=3):
#> X_train_microD1S305 X_train_microD8S283 X_train_microD10S192 X_train_microsexe X_train_microSiege
#> - 4th step (K=4):
#> X_train_microD17S794 X_train_microD13S173 X_train_microTP53 X_train_microD6S264 X_train_microD2S159 X_train_microD6S275
#> - 5th step (K=5):
#> X_train_microD18S61 X_train_microD9S171 X_train_microD8S264 X_train_microD18S53 X_train_microD4S394 X_train_microD11S916
#> - 6th step (K=6):
#> X_train_microD17S790
#> Call:
#> coxph(formula = YCsurv ~ ., data = tt_splsDR)
#> 
#>         coef exp(coef) se(coef)     z        p
#> dim.1 0.8093    2.2462   0.2029 3.989 6.63e-05
#> dim.2 0.9295    2.5333   0.2939 3.163  0.00156
#> dim.3 0.9968    2.7096   0.4190 2.379  0.01736
#> dim.4 0.9705    2.6391   0.3793 2.558  0.01052
#> dim.5 0.2162    1.2413   0.2811 0.769  0.44192
#> dim.6 0.4380    1.5496   0.3608 1.214  0.22473
#> 
#> Likelihood ratio test=55.06  on 6 df, p=4.51e-10
#> n= 80, number of events= 17 
(cox_splsDR_fit3=coxsplsDR(~.,Y_train_micro,C_train_micro,ncomp=6,
dataXplan=X_train_micro_df,eta=.5))
#> Call:
#> coxph(formula = YCsurv ~ ., data = tt_splsDR)
#> 
#>         coef exp(coef) se(coef)     z        p
#> dim.1 0.8093    2.2462   0.2029 3.989 6.63e-05
#> dim.2 0.9295    2.5333   0.2939 3.163  0.00156
#> dim.3 0.9968    2.7096   0.4190 2.379  0.01736
#> dim.4 0.9705    2.6391   0.3793 2.558  0.01052
#> dim.5 0.2162    1.2413   0.2811 0.769  0.44192
#> dim.6 0.4380    1.5496   0.3608 1.214  0.22473
#> 
#> Likelihood ratio test=55.06  on 6 df, p=4.51e-10
#> n= 80, number of events= 17 

rm(X_train_micro,Y_train_micro,C_train_micro,cox_splsDR_fit,cox_splsDR_fit2,cox_splsDR_fit3)
```
