\name{print.summary.plsRcoxmodel}
\alias{print.summary.plsRcoxmodel}
\title{Print method for summaries of plsRcox models}
\description{
This function provides a print method for the class \code{"summary.plsRcoxmodel"}
}
\usage{
\method{print}{summary.plsRcoxmodel}(x, \dots)
}
\arguments{
  \item{x}{an object of the class \code{"summary.plsRcoxmodel"}}
  \item{\dots}{not used}
}
%\details{}
\value{
  \item{language}{call of the model}
}
\references{\enc{Fr?d?ric}{Fr\'ed\'eric} Bertrand, Myriam Maumy-Bertrand et Nicolas Meyer (2011). \enc{R?gression}{R\'egression} \enc{B?ta}{B\^eta} PLS. \emph{Preprint}.
}
\author{\enc{Fr?d?ric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{\code{\link{print}} and \code{\link{summary}}}
\examples{
data(micro.censure)
data(Xmicro.censure_compl_imp)

X_train_micro <- apply((as.matrix(Xmicro.censure_compl_imp)),FUN="as.numeric",MARGIN=2)[1:80,]
Y_train_micro <- micro.censure$survyear[1:80]
C_train_micro <- micro.censure$DC[1:80]

modpls <- plsRcox(X_train_micro,time=Y_train_micro,event=C_train_micro,nt=3)
print(summary(modpls))
}
\keyword{methods}
\keyword{print}
