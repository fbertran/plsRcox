# Partial least squares Regression generalized linear models

This function implements an extension of Partial least squares
Regression to Cox Models.

## Usage

``` r
DKplsRcox(Xplan, ...)

DKplsRcoxmodel(Xplan, ...)

# Default S3 method
DKplsRcoxmodel(
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
  nt = min(2, ncol(Xplan)),
  limQ2set = 0.0975,
  dataPredictY = Xplan,
  pvals.expli = FALSE,
  alpha.pvals.expli = 0.05,
  tol_Xi = 10^(-12),
  weights,
  control,
  sparse = FALSE,
  sparseStop = TRUE,
  plot = FALSE,
  allres = FALSE,
  kernel = "rbfdot",
  hyperkernel,
  verbose = TRUE,
  ...
)

# S3 method for class 'formula'
DKplsRcoxmodel(
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
  scaleY = NULL,
  dataXplan = NULL,
  nt = min(2, ncol(Xplan)),
  limQ2set = 0.0975,
  dataPredictY = Xplan,
  pvals.expli = FALSE,
  model_frame = FALSE,
  alpha.pvals.expli = 0.05,
  tol_Xi = 10^(-12),
  weights,
  subset,
  control,
  sparse = FALSE,
  sparseStop = TRUE,
  plot = FALSE,
  allres = FALSE,
  kernel = "rbfdot",
  hyperkernel,
  verbose = TRUE,
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

  arguments to pass to `plsRmodel.default` or to `plsRmodel.formula`

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

  number of components to be extracted

- limQ2set:

  limit value for the Q2

- dataPredictY:

  predictor(s) (testing) dataset

- pvals.expli:

  should individual p-values be reported to tune model selection ?

- alpha.pvals.expli:

  level of significance for predictors when pvals.expli=TRUE

- tol_Xi:

  minimal value for Norm2(Xi) and \\\mathrm{det}(pp' \times pp)\\ if
  there is any missing value in the `dataX`. It defaults to \\10^{-12}\\

- weights:

  an optional vector of 'prior weights' to be used in the fitting
  process. Should be `NULL` or a numeric vector.

- control:

  a list of parameters for controlling the fitting process. For
  `glm.fit` this is passed to
  [`glm.control`](https://rdrr.io/r/stats/glm.control.html).

- sparse:

  should the coefficients of non-significant predictors
  (\<`alpha.pvals.expli`) be set to 0

- sparseStop:

  should component extraction stop when no significant predictors
  (\<`alpha.pvals.expli`) are found

- plot:

  Should the survival function be plotted ?)

- allres:

  FALSE to return only the Cox model and TRUE for additionnal results.
  See details. Defaults to FALSE.

- kernel:

  the kernel function used in training and predicting. This parameter
  can be set to any function, of class kernel, which computes the inner
  product in feature space between two vector arguments (see
  [kernels](https://rdrr.io/pkg/kernlab/man/dots.html)). The `kernlab`
  package provides the most popular kernel functions which can be used
  by setting the kernel parameter to the following strings:

  list("rbfdot")

  :   Radial Basis kernel "Gaussian"

  list("polydot")

  :   Polynomial kernel

  list("vanilladot")

  :   Linear kernel

  list("tanhdot")

  :   Hyperbolic tangent kernel

  list("laplacedot")

  :   Laplacian kernel

  list("besseldot")

  :   Bessel kernel

  list("anovadot")

  :   ANOVA RBF kernel

  list("splinedot")

  :   Spline kernel

- hyperkernel:

  the list of hyper-parameters (kernel parameters). This is a list which
  contains the parameters to be used with the kernel function. For valid
  parameters for existing kernels are :

  - `sigma`, inverse kernel width for the Radial Basis kernel function
    "rbfdot" and the Laplacian kernel "laplacedot".

  - `degree`, `scale`, `offset` for the Polynomial kernel "polydot".

  - `scale`, offset for the Hyperbolic tangent kernel function
    "tanhdot".

  - `sigma`, `order`, `degree` for the Bessel kernel "besseldot".

  - `sigma`, `degree` for the ANOVA kernel "anovadot".

  In the case of a Radial Basis kernel function (Gaussian) or Laplacian
  kernel, if `hyperkernel` is missing, the heuristics in sigest are used
  to calculate a good sigma value from the data.

- verbose:

  Should some details be displayed ?

- dataXplan:

  an optional data frame, list or environment (or object coercible by
  [`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html) to a data
  frame) containing the variables in the model. If not found in
  `dataXplan`, the variables are taken from `environment(Xplan)`,
  typically the environment from which `coxDKplsDR` is called.

- model_frame:

  If `TRUE`, the model frame is returned.

- subset:

  an optional vector specifying a subset of observations to be used in
  the fitting process.

- model_matrix:

  If `TRUE`, the model matrix is returned.

- contrasts.arg:

  a list, whose entries are values (numeric matrices, functions or
  character strings naming functions) to be used as replacement values
  for the contrasts replacement function and whose names are the names
  of columns of data containing factors.

- method:

  the method to be used in fitting the model. The default method
  `"glm.fit"` uses iteratively reweighted least squares (IWLS).
  User-supplied fitting functions can be supplied either as a function
  or a character string naming a function, with a function which takes
  the same arguments as `glm.fit`.

## Value

Depends on the model that was used to fit the model.

## Details

A typical predictor has the form response ~ terms where response is the
(numeric) response vector and terms is a series of terms which specifies
a linear predictor for response. A terms specification of the form
first + second indicates all the terms in first together with all the
terms in second with any duplicates removed.

A specification of the form first:second indicates the the set of terms
obtained by taking the interactions of all terms in first with all terms
in second. The specification first\*second indicates the cross of first
and second. This is the same as first + second + first:second.

The terms in the formula will be re-ordered so that main effects come
first, followed by the interactions, all second-order, all third-order
and so on: to avoid this pass a terms object as the formula.

Non-NULL weights can be used to indicate that different observations
have different dispersions (with the values in weights being inversely
proportional to the dispersions); or equivalently, when the elements of
weights are positive integers w_i, that each response y_i is the mean of
w_i unit-weight observations.

## References

plsRcox, Cox-Models in a high dimensional setting in R, Frederic
Bertrand, Philippe Bastien, Nicolas Meyer and Myriam Maumy-Bertrand
(2014). Proceedings of User2014!, Los Angeles, page 152.  

Deviance residuals-based sparse PLS and sparse kernel PLS regression for
censored data, Philippe Bastien, Frederic Bertrand, Nicolas Meyer and
Myriam Maumy-Bertrand (2015), Bioinformatics, 31(3):397-404,
doi:10.1093/bioinformatics/btu660.

## See also

[`plsR`](https://fbertran.github.io/plsRglm/reference/plsR.html) and
[`plsRglm`](https://fbertran.github.io/plsRglm/reference/plsRglm.html)

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

# \donttest{
DKplsRcox(X_train_micro,time=Y_train_micro,event=C_train_micro,nt=5)
#> Kernel :  rbfdot 
#> Estimated_sigma  0.01201981 
#> ____************************************************____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
#> Call:
#> coxph(formula = YCsurv ~ ., data = tt_DKplsRcox)
#> 
#>              coef  exp(coef)   se(coef)      z        p
#> Comp_1 -9.226e+00  9.843e-05  3.242e+00 -2.846 0.004430
#> Comp_2  2.706e+01  5.631e+11  7.135e+00  3.792 0.000149
#> Comp_3  2.084e+01  1.127e+09  5.748e+00  3.626 0.000288
#> Comp_4  8.992e+00  8.037e+03  2.990e+00  3.007 0.002638
#> Comp_5  8.378e+00  4.350e+03  2.612e+00  3.208 0.001338
#> 
#> Likelihood ratio test=69.52  on 5 df, p=1.288e-13
#> n= 80, number of events= 17 
DKplsRcox(~X_train_micro,time=Y_train_micro,event=C_train_micro,nt=5)
#> Kernel :  rbfdot 
#> Estimated_sigma  0.01273187 
#> ____************************************************____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
#> Call:
#> coxph(formula = YCsurv ~ ., data = tt_DKplsRcox)
#> 
#>              coef  exp(coef)   se(coef)      z        p
#> Comp_1 -7.537e+00  5.332e-04  2.895e+00 -2.604 0.009221
#> Comp_2  2.660e+01  3.566e+11  7.108e+00  3.743 0.000182
#> Comp_3  2.130e+01  1.788e+09  5.929e+00  3.593 0.000327
#> Comp_4  9.790e+00  1.786e+04  3.195e+00  3.064 0.002185
#> Comp_5  8.626e+00  5.577e+03  2.669e+00  3.232 0.001230
#> 
#> Likelihood ratio test=71.01  on 5 df, p=6.301e-14
#> n= 80, number of events= 17 
# }

# \donttest{
DKplsRcox(Xplan=X_train_micro,time=Y_train_micro,event=C_train_micro,nt=5,sparse=TRUE, 
alpha.pvals.expli=.15)
#> Kernel :  rbfdot 
#> Estimated_sigma  0.01358702 
#> ____************************************************____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
#> Call:
#> coxph(formula = YCsurv ~ ., data = tt_DKplsRcox)
#> 
#>              coef  exp(coef)   se(coef)      z        p
#> Comp_1 -6.158e+01  1.795e-27  1.906e+01 -3.231 0.001235
#> Comp_2  9.489e+01  1.619e+41  2.860e+01  3.317 0.000909
#> Comp_3  4.526e+01  4.544e+19  1.381e+01  3.276 0.001051
#> Comp_4  4.668e+01  1.867e+20  1.506e+01  3.100 0.001933
#> Comp_5  9.479e+00  1.309e+04  3.505e+00  2.705 0.006838
#> 
#> Likelihood ratio test=82.52  on 5 df, p=2.492e-16
#> n= 80, number of events= 17 
DKplsRcox(Xplan=~X_train_micro,time=Y_train_micro,event=C_train_micro,nt=5,sparse=TRUE,
alpha.pvals.expli=.15)
#> Kernel :  rbfdot 
#> Estimated_sigma  0.01185157 
#> ____************************************************____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
#> Call:
#> coxph(formula = YCsurv ~ ., data = tt_DKplsRcox)
#> 
#>              coef  exp(coef)   se(coef)      z       p
#> Comp_1 -6.081e+01  3.907e-27  1.910e+01 -3.183 0.00146
#> Comp_2  8.404e+01  3.154e+36  2.583e+01  3.254 0.00114
#> Comp_3  3.885e+01  7.421e+16  1.211e+01  3.207 0.00134
#> Comp_4  3.647e+01  6.894e+15  1.207e+01  3.020 0.00252
#> Comp_5  7.230e+00  1.381e+03  2.550e+00  2.835 0.00458
#> 
#> Likelihood ratio test=77.29  on 5 df, p=3.09e-15
#> n= 80, number of events= 17 
# }
```
