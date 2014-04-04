\name{summary.plsRcoxmodel}
\alias{summary.plsRcoxmodel}
\title{Summary method for plsRcox models}
\description{
This function provides a summary method for the class \code{"plsRcoxmodel"}
}
\usage{
\method{summary}{plsRcoxmodel}(object, \dots)
}
\arguments{
  \item{object}{an object of the class \code{"plsRcoxmodel"}}
  \item{\dots}{further arguments to be passed to or from methods.}  
}
%\details{}
\value{
  \item{call }{function call of plsR beta models}
}
\references{\enc{Frédéric}{Fr\'ed\'eric} Bertrand, Myriam Maumy-Bertrand et Nicolas Meyer (2011). \enc{Régression}{R\'egression} \enc{Bêta}{B\^eta} PLS. \emph{Preprint}.
}
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{\code{\link{summary}}}
\examples{
\donttest{
data(micro.censure)
data(Xmicro.censure_compl_imp)

X_train_micro <- apply((as.matrix(Xmicro.censure_compl_imp)),FUN="as.numeric",MARGIN=2)[1:80,]
Y_train_micro <- micro.censure$survyear[1:80]
C_train_micro <- micro.censure$DC[1:80]

modpls <- plsRcox(X_train_micro,time=Y_train_micro,event=C_train_micro,nt=3)
summary(modpls)
}
}
\keyword{methods}
\keyword{print}
