# Microsat features and survival times

This dataset provides Microsat specifications and survival times.

## Format

A data frame with 117 observations on the following 43 variables.

- numpat:

  a factor with levels `B1006` `B1017` `B1028` `B1031` `B1046` `B1059`
  `B1068` `B1071` `B1102` `B1115` `B1124` `B1139` `B1157` `B1161`
  `B1164` `B1188` `B1190` `B1192` `B1203` `B1211` `B1221` `B1225`
  `B1226` `B1227` `B1237` `B1251` `B1258` `B1266` `B1271` `B1282`
  `B1284` `B1285` `B1286` `B1287` `B1290` `B1292` `B1298` `B1302`
  `B1304` `B1310` `B1319` `B1327` `B1353` `B1357` `B1363` `B1368`
  `B1372` `B1373` `B1379` `B1388` `B1392` `B1397` `B1403` `B1418`
  `B1421t1` `B1421t2` `B1448` `B1451` `B1455` `B1460` `B1462` `B1466`
  `B1469` `B1493` `B1500` `B1502` `B1519` `B1523` `B1529` `B1530`
  `B1544` `B1548` `B500` `B532` `B550` `B558` `B563` `B582` `B605`
  `B609` `B634` `B652` `B667` `B679` `B701` `B722` `B728` `B731` `B736`
  `B739` `B744` `B766` `B771` `B777` `B788` `B800` `B836` `B838` `B841`
  `B848` `B871` `B873` `B883` `B889` `B912` `B924` `B925` `B927` `B938`
  `B952` `B954` `B955` `B968` `B972` `B976` `B982` `B984`

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

- survyear:

  a numeric vector

- DC:

  a numeric vector

## Source

Allelotyping identification of genomic alterations in rectal
chromosomally unstable tumors without preoperative treatment, \#' Benoît
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
data(micro.censure)
Y_train_micro <- micro.censure$survyear[1:80]
C_train_micro <- micro.censure$DC[1:80]
Y_test_micro <- micro.censure$survyear[81:117]
C_test_micro <- micro.censure$DC[81:117]
rm(Y_train_micro,C_train_micro,Y_test_micro,C_test_micro)
# }
```
