# Partial least squares Regression generalized linear models

This function implements an extension of Partial least squares
Regression to Cox Models.

## Usage

``` r
plsRcox(Xplan, ...)

plsRcoxmodel(Xplan, ...)

# Default S3 method
plsRcoxmodel(
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
  allres = TRUE,
  verbose = TRUE,
  ...
)

# S3 method for class 'formula'
plsRcoxmodel(
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
  allres = TRUE,
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

- allres:

  FALSE to return only the Cox model and TRUE for additionnal results.
  See details. Defaults to FALSE.

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

plsRcox(X_train_micro,time=Y_train_micro,event=C_train_micro,nt=5)
#> ____************************************************____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
#> Number of required components:
#> [1] 5
#> Number of successfully computed components:
#> [1] 5
#> Coefficients:
#>                [,1]
#> D18S61   0.68964919
#> D17S794 -1.14362392
#> D13S173  1.37632457
#> D20S107  4.96128745
#> TP53     1.68453950
#> D9S171  -1.46691252
#> D8S264   0.66710776
#> D5S346  -4.61338196
#> D22S928 -1.82005524
#> D18S53   0.79853646
#> D1S225  -1.46234986
#> D3S1282 -1.67925042
#> D15S127  3.92225537
#> D1S305  -2.29680161
#> D1S207   2.02539691
#> D2S138  -3.48975878
#> D16S422 -2.92189625
#> D9S179  -0.59484679
#> D10S191 -1.30136747
#> D4S394   1.34265359
#> D1S197  -0.75014044
#> D6S264   1.32746604
#> D14S65  -3.20882866
#> D17S790  0.55427680
#> D5S430   3.40654627
#> D3S1283  2.12510239
#> D4S414   2.73619967
#> D8S283   0.71955323
#> D11S916  1.45026508
#> D2S159   0.90293134
#> D16S408 -0.59719901
#> D6S275  -1.02204186
#> D10S192  1.14220367
#> sexe     0.67314561
#> Agediag  0.04908478
#> Siege   -0.41985924
#> T        2.70581463
#> N        2.47039973
#> M       -4.53213922
#> STADE    0.48221697
#> Information criteria and Fit statistics:
#>                 AIC       BIC
#> Nb_Comp_0 112.87990 112.87990
#> Nb_Comp_1  85.11075  87.49278
#> Nb_Comp_2  75.49537  80.25942
#> Nb_Comp_3  68.45852  75.60460
#> Nb_Comp_4  63.09284  72.62094
#> Nb_Comp_5  55.30567  67.21581
plsRcox(~X_train_micro,time=Y_train_micro,event=C_train_micro,nt=5)
#> ____************************************************____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Component____ 4 ____
#> ____Component____ 5 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
#> Number of required components:
#> [1] 5
#> Number of successfully computed components:
#> [1] 5
#> Coefficients:
#>                             [,1]
#> X_train_microD18S61   0.68964919
#> X_train_microD17S794 -1.14362392
#> X_train_microD13S173  1.37632457
#> X_train_microD20S107  4.96128745
#> X_train_microTP53     1.68453950
#> X_train_microD9S171  -1.46691252
#> X_train_microD8S264   0.66710776
#> X_train_microD5S346  -4.61338196
#> X_train_microD22S928 -1.82005524
#> X_train_microD18S53   0.79853646
#> X_train_microD1S225  -1.46234986
#> X_train_microD3S1282 -1.67925042
#> X_train_microD15S127  3.92225537
#> X_train_microD1S305  -2.29680161
#> X_train_microD1S207   2.02539691
#> X_train_microD2S138  -3.48975878
#> X_train_microD16S422 -2.92189625
#> X_train_microD9S179  -0.59484679
#> X_train_microD10S191 -1.30136747
#> X_train_microD4S394   1.34265359
#> X_train_microD1S197  -0.75014044
#> X_train_microD6S264   1.32746604
#> X_train_microD14S65  -3.20882866
#> X_train_microD17S790  0.55427680
#> X_train_microD5S430   3.40654627
#> X_train_microD3S1283  2.12510239
#> X_train_microD4S414   2.73619967
#> X_train_microD8S283   0.71955323
#> X_train_microD11S916  1.45026508
#> X_train_microD2S159   0.90293134
#> X_train_microD16S408 -0.59719901
#> X_train_microD6S275  -1.02204186
#> X_train_microD10S192  1.14220367
#> X_train_microsexe     0.67314561
#> X_train_microAgediag  0.04908478
#> X_train_microSiege   -0.41985924
#> X_train_microT        2.70581463
#> X_train_microN        2.47039973
#> X_train_microM       -4.53213922
#> X_train_microSTADE    0.48221697
#> Information criteria and Fit statistics:
#>                 AIC       BIC
#> Nb_Comp_0 112.87990 112.87990
#> Nb_Comp_1  85.11075  87.49278
#> Nb_Comp_2  75.49537  80.25942
#> Nb_Comp_3  68.45852  75.60460
#> Nb_Comp_4  63.09284  72.62094
#> Nb_Comp_5  55.30567  67.21581

plsRcox(Xplan=X_train_micro,time=Y_train_micro,event=C_train_micro,nt=5,sparse=TRUE,
alpha.pvals.expli=.15)
#> ____************************************************____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> Warning : 25 < 10^{-12}
#> Warning only 3 components could thus be extracted
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
#> Number of required components:
#> [1] 5
#> Number of successfully computed components:
#> [1] 3
#> Coefficients:
#>                [,1]
#> D18S61   0.00000000
#> D17S794  0.00000000
#> D13S173  0.00000000
#> D20S107  2.22871454
#> TP53     0.00000000
#> D9S171   0.00000000
#> D8S264   0.00000000
#> D5S346  -1.20298526
#> D22S928  0.00000000
#> D18S53   0.00000000
#> D1S225  -1.29459798
#> D3S1282 -1.99426291
#> D15S127  1.39645601
#> D1S305   0.00000000
#> D1S207   1.25164327
#> D2S138  -1.65740160
#> D16S422  0.00000000
#> D9S179   0.00000000
#> D10S191 -1.25360805
#> D4S394   0.00000000
#> D1S197   0.00000000
#> D6S264   0.00000000
#> D14S65  -1.33587373
#> D17S790  0.00000000
#> D5S430   1.72799213
#> D3S1283  0.00000000
#> D4S414   1.03558702
#> D8S283   0.00000000
#> D11S916  0.00000000
#> D2S159   0.00000000
#> D16S408 -1.75748257
#> D6S275   0.00000000
#> D10S192  0.00000000
#> sexe     0.00000000
#> Agediag  0.05075304
#> Siege    0.00000000
#> T        1.36569407
#> N        1.27485618
#> M       -1.17682617
#> STADE   -0.65106093
#> Information criteria and Fit statistics:
#>                 AIC       BIC
#> Nb_Comp_0 112.87990 112.87990
#> Nb_Comp_1  85.54313  87.92516
#> Nb_Comp_2  75.16125  79.92530
#> Nb_Comp_3  73.63097  80.77705
plsRcox(Xplan=~X_train_micro,time=Y_train_micro,event=C_train_micro,nt=5,sparse=TRUE,
alpha.pvals.expli=.15)
#> ____************************************************____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> Warning : 25 < 10^{-12}
#> Warning only 3 components could thus be extracted
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
#> Number of required components:
#> [1] 5
#> Number of successfully computed components:
#> [1] 3
#> Coefficients:
#>                             [,1]
#> X_train_microD18S61   0.00000000
#> X_train_microD17S794  0.00000000
#> X_train_microD13S173  0.00000000
#> X_train_microD20S107  2.22871454
#> X_train_microTP53     0.00000000
#> X_train_microD9S171   0.00000000
#> X_train_microD8S264   0.00000000
#> X_train_microD5S346  -1.20298526
#> X_train_microD22S928  0.00000000
#> X_train_microD18S53   0.00000000
#> X_train_microD1S225  -1.29459798
#> X_train_microD3S1282 -1.99426291
#> X_train_microD15S127  1.39645601
#> X_train_microD1S305   0.00000000
#> X_train_microD1S207   1.25164327
#> X_train_microD2S138  -1.65740160
#> X_train_microD16S422  0.00000000
#> X_train_microD9S179   0.00000000
#> X_train_microD10S191 -1.25360805
#> X_train_microD4S394   0.00000000
#> X_train_microD1S197   0.00000000
#> X_train_microD6S264   0.00000000
#> X_train_microD14S65  -1.33587373
#> X_train_microD17S790  0.00000000
#> X_train_microD5S430   1.72799213
#> X_train_microD3S1283  0.00000000
#> X_train_microD4S414   1.03558702
#> X_train_microD8S283   0.00000000
#> X_train_microD11S916  0.00000000
#> X_train_microD2S159   0.00000000
#> X_train_microD16S408 -1.75748257
#> X_train_microD6S275   0.00000000
#> X_train_microD10S192  0.00000000
#> X_train_microsexe     0.00000000
#> X_train_microAgediag  0.05075304
#> X_train_microSiege    0.00000000
#> X_train_microT        1.36569407
#> X_train_microN        1.27485618
#> X_train_microM       -1.17682617
#> X_train_microSTADE   -0.65106093
#> Information criteria and Fit statistics:
#>                 AIC       BIC
#> Nb_Comp_0 112.87990 112.87990
#> Nb_Comp_1  85.54313  87.92516
#> Nb_Comp_2  75.16125  79.92530
#> Nb_Comp_3  73.63097  80.77705
```
