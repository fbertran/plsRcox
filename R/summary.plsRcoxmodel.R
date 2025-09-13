#' Summary method for plsRcox models
#' 
#' This function provides a summary method for the class \code{"plsRcoxmodel"}
#' 
#' 
#' @param object an object of the class \code{"plsRcoxmodel"}
#' @param \dots further arguments to be passed to or from methods.
#' @return \item{call }{function call of plsRcox models}
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@lecnam.net}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link{summary}}
#' @references plsRcox, Cox-Models in a high dimensional setting in R, Frederic
#' Bertrand, Philippe Bastien, Nicolas Meyer and Myriam Maumy-Bertrand (2014).
#' Proceedings of User2014!, Los Angeles, page 152.\cr
#' 
#' Deviance residuals-based sparse PLS and sparse kernel PLS regression for
#' censored data, Philippe Bastien, Frederic Bertrand, Nicolas Meyer and Myriam
#' Maumy-Bertrand (2015), Bioinformatics, 31(3):397-404,
#' doi:10.1093/bioinformatics/btu660.
#' @keywords methods print
#' @examples
#' 
#' data(micro.censure)
#' data(Xmicro.censure_compl_imp)
#' 
#' X_train_micro <- apply((as.matrix(Xmicro.censure_compl_imp)),FUN="as.numeric",MARGIN=2)[1:80,]
#' Y_train_micro <- micro.censure$survyear[1:80]
#' C_train_micro <- micro.censure$DC[1:80]
#' 
#' modpls <- plsRcox(X_train_micro,time=Y_train_micro,event=C_train_micro,nt=3)
#' summary(modpls)
#' 
#' @export 
summary.plsRcoxmodel <- function(object, ...)
{
res <- list(call=object$call)
class(res) <- "summary.plsRcoxmodel"
res
}
