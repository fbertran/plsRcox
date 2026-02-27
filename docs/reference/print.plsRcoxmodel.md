# Print method for plsRcox models

This function provides a print method for the class `"plsRcoxmodel"`

## Usage

``` r
# S3 method for class 'plsRcoxmodel'
print(x, ...)
```

## Arguments

- x:

  an object of the class `"plsRcoxmodel"`

- ...:

  not used

## Value

`NULL`

## References

plsRcox, Cox-Models in a high dimensional setting in R, Frederic
Bertrand, Philippe Bastien, Nicolas Meyer and Myriam Maumy-Bertrand
(2014). Proceedings of User2014!, Los Angeles, page 152.  

Deviance residuals-based sparse PLS and sparse kernel PLS regression for
censored data, Philippe Bastien, Frederic Bertrand, Nicolas Meyer and
Myriam Maumy-Bertrand (2015), Bioinformatics, 31(3):397-404,
doi:10.1093/bioinformatics/btu660.

## See also

[`print`](https://rdrr.io/r/base/print.html)

## Author

Frédéric Bertrand  
<frederic.bertrand@lecnam.net>  
<https://fbertran.github.io/homepage/>

## Examples

``` r
data(micro.censure)
data(Xmicro.censure_compl_imp)

X_train_micro <- apply((as.matrix(Xmicro.censure_compl_imp)),FUN="as.numeric",MARGIN=2)[1:80,]
Y_train_micro <- micro.censure$survyear[1:80]
C_train_micro <- micro.censure$DC[1:80]

modpls <- plsRcox(X_train_micro,time=Y_train_micro,event=C_train_micro,nt=3)
#> ____************************************************____
#> ____Component____ 1 ____
#> ____Component____ 2 ____
#> ____Component____ 3 ____
#> ____Predicting X without NA neither in X nor in Y____
#> ****________________________________________________****
#> 
print(modpls)
#> Number of required components:
#> [1] 3
#> Number of successfully computed components:
#> [1] 3
#> Coefficients:
#>                 [,1]
#> D18S61   0.657688859
#> D17S794 -0.265485544
#> D13S173  0.532071747
#> D20S107  2.764628048
#> TP53     0.635427658
#> D9S171   0.008139129
#> D8S264  -0.346586438
#> D5S346  -1.628707075
#> D22S928 -1.199432030
#> D18S53   0.550835752
#> D1S225  -1.098480981
#> D3S1282 -1.784482327
#> D15S127  1.905056253
#> D1S305  -1.028283057
#> D1S207   1.202494887
#> D2S138  -1.610961966
#> D16S422 -0.970535096
#> D9S179  -0.209672191
#> D10S191 -1.143815474
#> D4S394   0.239525569
#> D1S197   0.087674404
#> D6S264   0.289838007
#> D14S65  -1.281410428
#> D17S790 -0.335500453
#> D5S430   0.789195774
#> D3S1283  0.453349027
#> D4S414   1.313974219
#> D8S283  -0.179467540
#> D11S916  0.457823141
#> D2S159   0.719452513
#> D16S408 -1.343339387
#> D6S275  -0.568676682
#> D10S192 -0.011708963
#> sexe    -0.080266201
#> Agediag  0.051845736
#> Siege   -0.157190141
#> T        0.865566178
#> N        0.903857312
#> M       -0.883429770
#> STADE   -0.079069085
#> Information criteria and Fit statistics:
#>                 AIC       BIC
#> Nb_Comp_0 112.87990 112.87990
#> Nb_Comp_1  85.11075  87.49278
#> Nb_Comp_2  75.49537  80.25942
#> Nb_Comp_3  68.45852  75.60460
```
