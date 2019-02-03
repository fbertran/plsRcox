\name{Xmicro.censure_compl_imp}
\alias{Xmicro.censure_compl_imp}
\docType{data}
\title{Imputed Microsat features}
\description{
This dataset provides imputed microsat specifications. Imputations were computed using Multivariate Imputation by Chained Equations (MICE) using predictive mean matching for the numeric columns, logistic regression imputation for the binary data or the factors with 2 levels and polytomous regression imputation for categorical data i.e. factors with three or more levels.
}
\format{
  A data frame with 117 observations on the following 40 variables.
  \describe{
    \item{\code{D18S61}}{a numeric vector}
    \item{\code{D17S794}}{a numeric vector}
    \item{\code{D13S173}}{a numeric vector}
    \item{\code{D20S107}}{a numeric vector}
    \item{\code{TP53}}{a numeric vector}
    \item{\code{D9S171}}{a numeric vector}
    \item{\code{D8S264}}{a numeric vector}
    \item{\code{D5S346}}{a numeric vector}
    \item{\code{D22S928}}{a numeric vector}
    \item{\code{D18S53}}{a numeric vector}
    \item{\code{D1S225}}{a numeric vector}
    \item{\code{D3S1282}}{a numeric vector}
    \item{\code{D15S127}}{a numeric vector}
    \item{\code{D1S305}}{a numeric vector}
    \item{\code{D1S207}}{a numeric vector}
    \item{\code{D2S138}}{a numeric vector}
    \item{\code{D16S422}}{a numeric vector}
    \item{\code{D9S179}}{a numeric vector}
    \item{\code{D10S191}}{a numeric vector}
    \item{\code{D4S394}}{a numeric vector}
    \item{\code{D1S197}}{a numeric vector}
    \item{\code{D6S264}}{a numeric vector}
    \item{\code{D14S65}}{a numeric vector}
    \item{\code{D17S790}}{a numeric vector}
    \item{\code{D5S430}}{a numeric vector}
    \item{\code{D3S1283}}{a numeric vector}
    \item{\code{D4S414}}{a numeric vector}
    \item{\code{D8S283}}{a numeric vector}
    \item{\code{D11S916}}{a numeric vector}
    \item{\code{D2S159}}{a numeric vector}
    \item{\code{D16S408}}{a numeric vector}
    \item{\code{D6S275}}{a numeric vector}
    \item{\code{D10S192}}{a numeric vector}
    \item{\code{sexe}}{a numeric vector}
    \item{\code{Agediag}}{a numeric vector}
    \item{\code{Siege}}{a numeric vector}
    \item{\code{T}}{a numeric vector}
    \item{\code{N}}{a numeric vector}
    \item{\code{M}}{a numeric vector}
    \item{\code{STADE}}{a factor with levels \code{0} \code{1} \code{2} \code{3} \code{4}}
  }
}
\source{Allelotyping identification of genomic alterations in rectal chromosomally unstable tumors without preoperative treatment, \enc{Benoît}{Beno\^it} Romain, \enc{Agnès}{Agn\`es} Neuville, Nicolas Meyer, \enc{Cécile}{C\'ecile} Brigand, Serge Rohr, Anne Schneider, Marie-Pierre Gaub and Dominique Guenot, \emph{BMC Cancer 2010}, 10:561, doi:10.1186/1471-2407-10-561.}
\references{
plsRcox, Cox-Models in a high dimensional setting in R, Frederic Bertrand, Philippe Bastien, Nicolas Meyer and Myriam Maumy-Bertrand (2014). Proceedings of User2014!, Los Angeles, page 152.\cr
  
Deviance residuals-based sparse PLS and sparse kernel PLS regression for censored data, Philippe Bastien, Frederic Bertrand, Nicolas Meyer and Myriam Maumy-Bertrand (2015), Bioinformatics, 31(3):397-404, doi:10.1093/bioinformatics/btu660.
}
\examples{
\donttest{
data(Xmicro.censure_compl_imp)
X_train_micro <- Xmicro.censure_compl_imp[1:80,]
X_test_micro <- Xmicro.censure_compl_imp[81:117,]
rm(X_train_micro,X_test_micro)
}
}
\keyword{datasets}

