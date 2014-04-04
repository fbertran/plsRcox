\name{print.plsRcoxmodel}
\alias{print.plsRcoxmodel}
\title{Print method for plsRcox models}
\description{
This function provides a print method for the class \code{"plsRcoxmodel"}
}
\usage{
\method{print}{plsRcoxmodel}(x, \dots)
}
\arguments{
  \item{x}{an object of the class \code{"plsRcoxmodel"}}
  \item{\dots}{not used}
}
%\details{}
\value{
\code{NULL}
}
\references{\enc{Frédéric}{Fr\'ed\'eric} Bertrand, Myriam Maumy-Bertrand et Nicolas Meyer (2011). \enc{Régression}{R\'egression} \enc{Bêta}{B\^eta} PLS. \emph{Preprint}.
}
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{\code{\link{print}}}
\examples{
\donttest{
data(micro.censure)
data(Xmicro.censure_compl_imp)

X_train_micro <- apply((as.matrix(Xmicro.censure_compl_imp)),FUN="as.numeric",MARGIN=2)[1:80,]
Y_train_micro <- micro.censure$survyear[1:80]
C_train_micro <- micro.censure$DC[1:80]

modpls <- plsRcox(X_train_micro,time=Y_train_micro,event=C_train_micro,nt=3)
print(modpls)
}
}
\keyword{methods}
\keyword{print}
