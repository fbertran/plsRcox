# Fitting a PLSR model on the (Deviance) Residuals

This function computes the PLSR model with the Residuals of a Cox-Model
fitted with an intercept as the only explanatory variable as the
response and Xplan as explanatory variables. Default behaviour uses the
Deviance residuals. It uses the package `plsRglm`.

## Usage

``` r
coxpls3DR(Xplan, ...)

# Default S3 method
coxpls3DR(
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
  nt = min(7, ncol(Xplan)),
  typeVC = "none",
  plot = FALSE,
  allres = FALSE,
  sparse = FALSE,
  sparseStop = TRUE,
  ...
)

# S3 method for class 'formula'
coxpls3DR(
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
  nt = min(7, ncol(Xplan)),
  typeVC = "none",
  plot = FALSE,
  allres = FALSE,
  dataXplan = NULL,
  subset,
  weights,
  model_frame = FALSE,
  sparse = FALSE,
  sparseStop = TRUE,
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
  [`survival::coxph`](https://rdrr.io/pkg/survival/man/coxph.html) and
  to
  [`plsRglm::PLS_lm`](https://fbertran.github.io/plsRglm/reference/plsR.html).

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

- nt:

  Number of PLSR components to fit.

- typeVC:

  type of leave one out crossed validation. Several procedures are
  available and may be forced.

  list("none")

  :   no crossed validation

  list("standard")

  :   as in SIMCA for datasets without missing values and with all
      values predicted as those with missing values for datasets with
      any missing values

  list("missingdata")

  :   all values predicted as those with missing values for datasets
      with any missing values

  list("adaptative")

  :   predict a response value for an x with any missing value as those
      with missing values and for an x without any missing value as
      those without missing values.

- plot:

  Should the survival function be plotted ?)

- allres:

  FALSE to return only the Cox model and TRUE for additionnal results.
  See details. Defaults to FALSE.

- sparse:

  should the coefficients of non-significant predictors
  (\<`alpha.pvals.expli`) be set to 0

- sparseStop:

  should component extraction stop when no significant predictors
  (\<`alpha.pvals.expli`) are found

- dataXplan:

  an optional data frame, list or environment (or object coercible by
  [`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html) to a data
  frame) containing the variables in the model. If not found in
  `dataXplan`, the variables are taken from `environment(Xplan)`,
  typically the environment from which `coxpls3DR` is called.

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

- cox_pls3DR:

  Final Cox-model.

If `allres=TRUE` :

- tt_pls3DR:

  PLSR components.

- cox_pls3DR:

  Final Cox-model.

- pls3DR_mod:

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
[`PLS_lm`](https://fbertran.github.io/plsRglm/reference/plsR.html)

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

(cox_pls3DR_fit <- coxpls3DR(X_train_micro,Y_train_micro,C_train_micro,nt=7))
#> ____************************************************____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
#> Call:
#> coxph(formula = YCsurv ~ ., data = tt_pls3DR)
#> 
#>          coef exp(coef) se(coef)     z        p
#> Comp_1 0.8534    2.3476   0.2209 3.864 0.000111
#> Comp_2 1.0599    2.8862   0.3260 3.251 0.001150
#> Comp_3 0.9169    2.5016   0.4353 2.106 0.035166
#> Comp_4 0.9094    2.4827   0.4183 2.174 0.029697
#> Comp_5 0.1663    1.1809   0.2790 0.596 0.551259
#> Comp_6 0.8359    2.3068   0.4452 1.877 0.060455
#> Comp_7 0.4368    1.5478   0.2469 1.769 0.076872
#> 
#> Likelihood ratio test=57.92  on 7 df, p=3.919e-10
#> n= 80, number of events= 17 
(cox_pls3DR_fit2 <- coxpls3DR(~X_train_micro,Y_train_micro,C_train_micro,nt=7))
#> ____************************************************____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
#> Call:
#> coxph(formula = YCsurv ~ ., data = tt_pls3DR)
#> 
#>          coef exp(coef) se(coef)     z        p
#> Comp_1 0.8534    2.3476   0.2209 3.864 0.000111
#> Comp_2 1.0599    2.8862   0.3260 3.251 0.001150
#> Comp_3 0.9169    2.5016   0.4353 2.106 0.035166
#> Comp_4 0.9094    2.4827   0.4183 2.174 0.029697
#> Comp_5 0.1663    1.1809   0.2790 0.596 0.551259
#> Comp_6 0.8359    2.3068   0.4452 1.877 0.060455
#> Comp_7 0.4368    1.5478   0.2469 1.769 0.076872
#> 
#> Likelihood ratio test=57.92  on 7 df, p=3.919e-10
#> n= 80, number of events= 17 
(cox_pls3DR_fit3 <- coxpls3DR(~.,Y_train_micro,C_train_micro,nt=7,dataXplan=X_train_micro_df))
#> ____************************************************____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Component____ 6 ____
#> ____Component____ 7 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
#> Call:
#> coxph(formula = YCsurv ~ ., data = tt_pls3DR)
#> 
#>          coef exp(coef) se(coef)     z        p
#> Comp_1 0.8534    2.3476   0.2209 3.864 0.000111
#> Comp_2 1.0599    2.8862   0.3260 3.251 0.001150
#> Comp_3 0.9169    2.5016   0.4353 2.106 0.035166
#> Comp_4 0.9094    2.4827   0.4183 2.174 0.029697
#> Comp_5 0.1663    1.1809   0.2790 0.596 0.551259
#> Comp_6 0.8359    2.3068   0.4452 1.877 0.060455
#> Comp_7 0.4368    1.5478   0.2469 1.769 0.076872
#> 
#> Likelihood ratio test=57.92  on 7 df, p=3.919e-10
#> n= 80, number of events= 17 
(cox_pls3DR_fit4 <- coxpls3DR(~.,Y_train_micro,C_train_micro,nt=7,typeVC="none",
data=X_train_micro_df,sparse=TRUE))
#> ____************************************************____
#> No significant predictors (<0.05) found
#> Warning only one standard component (without sparse option) was thus extracted
#> ____Component____ 1 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
#> Call:
#> coxph(formula = YCsurv ~ ., data = tt_pls3DR)
#> 
#>          coef exp(coef) se(coef)     z        p
#> Comp_1 0.8922    2.4405   0.1859 4.799 1.59e-06
#> 
#> Likelihood ratio test=29.74  on 1 df, p=4.936e-08
#> n= 80, number of events= 17 
(cox_pls3DR_fit5 <- coxpls3DR(~.,Y_train_micro,C_train_micro,nt=7,typeVC="none",
data=X_train_micro_df,sparse=TRUE,sparseStop=FALSE))
#> ____************************************************____
#> No significant predictors (<0.05) found
#> Warning only one standard component (without sparse option) was thus extracted
#> ____Component____ 1 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
#> Call:
#> coxph(formula = YCsurv ~ ., data = tt_pls3DR)
#> 
#>          coef exp(coef) se(coef)     z        p
#> Comp_1 0.8922    2.4405   0.1859 4.799 1.59e-06
#> 
#> Likelihood ratio test=29.74  on 1 df, p=4.936e-08
#> n= 80, number of events= 17 

rm(X_train_micro,Y_train_micro,C_train_micro,cox_pls3DR_fit,cox_pls3DR_fit2,
cox_pls3DR_fit3,cox_pls3DR_fit4,cox_pls3DR_fit5)

```
