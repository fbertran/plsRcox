% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plsRcox-package.R
\docType{package}
\name{plsRcox-package}
\alias{plsRcox-package}
\title{plsRcox-package: Partial Least Squares Regression for Cox Models and Related Techniques}
\description{
Provides Partial least squares Regression and various regular, sparse or kernel, techniques for fitting Cox models in high dimensional settings <doi:10.1093/bioinformatics/btu660>, Bastien, P., Bertrand, F., Meyer N., Maumy-Bertrand, M. (2015), Deviance residuals-based sparse PLS and sparse kernel PLS regression for censored data, Bioinformatics, 31(3):397-404. Cross validation criteria were studied in <arXiv:1810.02962>, Bertrand, F., Bastien, Ph. and Maumy-Bertrand, M. (2018), Cross validating extensions of kernel, sparse or regular partial least squares regression models to censored data.
}
\examples{
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
cox_splsDR_fit2=coxsplsDR(~X_train_micro,Y_train_micro,C_train_micro,
ncomp=6,eta=.5,trace=TRUE)
cox_splsDR_fit2
cox_splsDR_fit3=coxsplsDR(~.,Y_train_micro,C_train_micro,ncomp=6,
dataXplan=X_train_micro_df,eta=.5)
cox_splsDR_fit3
rm(cox_splsDR_fit,cox_splsDR_fit2,cox_splsDR_fit3)

}
\references{
Bastien, P., Bertrand, F., Meyer N., Maumy-Bertrand, M. (2015), Deviance residuals-based sparse PLS and sparse kernel PLS regression for censored data, Bioinformatics, 31(3):397-404. <doi:10.1093/bioinformatics/btu660>. Cross validation criteria were studied in <arXiv:1810.02962>, Bertrand, F., Bastien, Ph. and Maumy-Bertrand, M. (2018), Cross validating extensions of kernel, sparse or regular partial least squares regression models to censored data.
}
