# plsRcox: Partial Least Squares Regression for Cox Models and Related Techniques

Provides Partial least squares Regression and various regular, sparse or
kernel, techniques for fitting Cox models in high dimensional settings
[doi:10.1093/bioinformatics/btu660](https://doi.org/10.1093/bioinformatics/btu660)
, Bastien, P., Bertrand, F., Meyer N., Maumy-Bertrand, M. (2015),
Deviance residuals-based sparse PLS and sparse kernel PLS regression for
censored data, Bioinformatics, 31(3):397-404. Cross validation criteria
were studied in
[doi:10.48550/arXiv.1810.02962](https://doi.org/10.48550/arXiv.1810.02962)
, Bertrand, F., Bastien, Ph. and Maumy-Bertrand, M. (2018), Cross
validating extensions of kernel, sparse or regular partial least squares
regression models to censored data.

## References

Bastien, P., Bertrand, F., Meyer N., Maumy-Bertrand, M. (2015), Deviance
residuals-based sparse PLS and sparse kernel PLS regression for censored
data, Bioinformatics, 31(3):397-404.
\<doi:10.1093/bioinformatics/btu660\>.

Cross validation criteria were studied in \<arXiv:1810.02962\>,
Bertrand, F., Bastien, Ph. and Maumy-Bertrand, M. (2018), Cross
validating extensions of kernel, sparse or regular partial least squares
regression models to censored data.

## See also

Useful links:

- <https://fbertran.github.io/plsRcox/>

- <https://github.com/fbertran/plsRcox>

- Report bugs at <https://github.com/fbertran/plsRcox/issues>

## Author

This package has been written by Frédéric Bertrand and Myriam
Maumy-Bertrand. Maintainer: Frédéric Bertrand
\<frederic.bertrand@lecnam.net\>

## Examples

``` r
# The original allelotyping dataset

library(plsRcox)
data(micro.censure)
Y_train_micro <- micro.censure$survyear[1:80]
C_train_micro <- micro.censure$DC[1:80]
Y_test_micro <- micro.censure$survyear[81:117]
C_test_micro <- micro.censure$DC[81:117]

data(Xmicro.censure_compl_imp)
X_train_micro <- apply((as.matrix(Xmicro.censure_compl_imp)),
FUN="as.numeric",MARGIN=2)[1:80,]
X_train_micro_df <- data.frame(X_train_micro)

# coxsplsDR
cox_splsDR_fit=coxsplsDR(X_train_micro,Y_train_micro,C_train_micro,
ncomp=6,eta=.5)
cox_splsDR_fit
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
cox_splsDR_fit2=coxsplsDR(~X_train_micro,Y_train_micro,C_train_micro,
ncomp=6,eta=.5,trace=TRUE)
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
cox_splsDR_fit2
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
cox_splsDR_fit3=coxsplsDR(~.,Y_train_micro,C_train_micro,ncomp=6,
dataXplan=X_train_micro_df,eta=.5)
cox_splsDR_fit3
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
rm(cox_splsDR_fit,cox_splsDR_fit2,cox_splsDR_fit3)
```
