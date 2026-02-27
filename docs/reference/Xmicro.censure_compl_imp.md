# Imputed Microsat features

This dataset provides imputed microsat specifications. Imputations were
computed using Multivariate Imputation by Chained Equations (MICE) using
predictive mean matching for the numeric columns, logistic regression
imputation for the binary data or the factors with 2 levels and
polytomous regression imputation for categorical data i.e. factors with
three or more levels.

## Format

A data frame with 117 observations on the following 40 variables.

- D18S61:

  a numeric vector

- D17S794:

  a numeric vector

- D13S173:

  a numeric vector

- D20S107:

  a numeric vector

- TP53:

  a numeric vector

- D9S171:

  a numeric vector

- D8S264:

  a numeric vector

- D5S346:

  a numeric vector

- D22S928:

  a numeric vector

- D18S53:

  a numeric vector

- D1S225:

  a numeric vector

- D3S1282:

  a numeric vector

- D15S127:

  a numeric vector

- D1S305:

  a numeric vector

- D1S207:

  a numeric vector

- D2S138:

  a numeric vector

- D16S422:

  a numeric vector

- D9S179:

  a numeric vector

- D10S191:

  a numeric vector

- D4S394:

  a numeric vector

- D1S197:

  a numeric vector

- D6S264:

  a numeric vector

- D14S65:

  a numeric vector

- D17S790:

  a numeric vector

- D5S430:

  a numeric vector

- D3S1283:

  a numeric vector

- D4S414:

  a numeric vector

- D8S283:

  a numeric vector

- D11S916:

  a numeric vector

- D2S159:

  a numeric vector

- D16S408:

  a numeric vector

- D6S275:

  a numeric vector

- D10S192:

  a numeric vector

- sexe:

  a numeric vector

- Agediag:

  a numeric vector

- Siege:

  a numeric vector

- T:

  a numeric vector

- N:

  a numeric vector

- M:

  a numeric vector

- STADE:

  a factor with levels `0` `1` `2` `3` `4`

## Source

Allelotyping identification of genomic alterations in rectal
chromosomally unstable tumors without preoperative treatment, Benoît
Romain, Agnès Neuville, Nicolas Meyer, Cécile Brigand, Serge Rohr, Anne
Schneider, Marie-Pierre Gaub and Dominique Guenot, *BMC Cancer 2010*,
10:561, doi:10.1186/1471-2407-10-561.

## References

plsRcox, Cox-Models in a high dimensional setting in R, Frederic
Bertrand, Philippe Bastien, Nicolas Meyer and Myriam Maumy-Bertrand
(2014). Proceedings of User2014!, Los Angeles, page 152.  

Deviance residuals-based sparse PLS and sparse kernel PLS regression for
censored data, Philippe Bastien, Frederic Bertrand, Nicolas Meyer and
Myriam Maumy-Bertrand (2015), Bioinformatics, 31(3):397-404,
doi:10.1093/bioinformatics/btu660.

## Examples

``` r
# \donttest{
data(Xmicro.censure_compl_imp)
X_train_micro <- Xmicro.censure_compl_imp[1:80,]
X_test_micro <- Xmicro.censure_compl_imp[81:117,]
rm(X_train_micro,X_test_micro)
# }
```
