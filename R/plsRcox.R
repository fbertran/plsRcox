#' Partial least squares Regression generalized linear models
#' 
#' This function implements an extension of Partial least squares Regression to
#' Cox Models.
#' 
#' 
#' A typical predictor has the form response ~ terms where response is the
#' (numeric) response vector and terms is a series of terms which specifies a
#' linear predictor for response. A terms specification of the form first +
#' second indicates all the terms in first together with all the terms in
#' second with any duplicates removed.
#' 
#' A specification of the form first:second indicates the the set of terms
#' obtained by taking the interactions of all terms in first with all terms in
#' second. The specification first*second indicates the cross of first and
#' second. This is the same as first + second + first:second.
#' 
#' The terms in the formula will be re-ordered so that main effects come first,
#' followed by the interactions, all second-order, all third-order and so on:
#' to avoid this pass a terms object as the formula.
#' 
#' Non-NULL weights can be used to indicate that different observations have
#' different dispersions (with the values in weights being inversely
#' proportional to the dispersions); or equivalently, when the elements of
#' weights are positive integers w_i, that each response y_i is the mean of w_i
#' unit-weight observations.
#' 
#' @aliases plsRcox plsRcoxmodel.default plsRcoxmodel.formula
#' @param Xplan a formula or a matrix with the eXplanatory variables (training)
#' dataset
#' @param time for right censored data, this is the follow up time. For
#' interval data, the first argument is the starting time for the interval.
#' @param time2 The status indicator, normally 0=alive, 1=dead. Other choices
#' are \code{TRUE/FALSE} (\code{TRUE} = death) or 1/2 (2=death). For interval
#' censored data, the status indicator is 0=right censored, 1=event at
#' \code{time}, 2=left censored, 3=interval censored. Although unusual, the
#' event indicator can be omitted, in which case all subjects are assumed to
#' have an event.
#' @param event ending time of the interval for interval censored or counting
#' process data only. Intervals are assumed to be open on the left and closed
#' on the right, \code{(start, end]}. For counting process data, event
#' indicates whether an event occurred at the end of the interval.
#' @param type character string specifying the type of censoring. Possible
#' values are \code{"right"}, \code{"left"}, \code{"counting"},
#' \code{"interval"}, or \code{"interval2"}. The default is \code{"right"} or
#' \code{"counting"} depending on whether the \code{time2} argument is absent
#' or present, respectively.
#' @param origin for counting process data, the hazard function origin. This
#' option was intended to be used in conjunction with a model containing time
#' dependent strata in order to align the subjects properly when they cross
#' over from one strata to another, but it has rarely proven useful.
#' @param typeres character string indicating the type of residual desired.
#' Possible values are \code{"martingale"}, \code{"deviance"}, \code{"score"},
#' \code{"schoenfeld"}, \code{"dfbeta"}, \code{"dfbetas"}, and
#' \code{"scaledsch"}. Only enough of the string to determine a unique match is
#' required.
#' @param collapse vector indicating which rows to collapse (sum) over. In
#' time-dependent models more than one row data can pertain to a single
#' individual. If there were 4 individuals represented by 3, 1, 2 and 4 rows of
#' data respectively, then \code{collapse=c(1,1,1,2,3,3,4,4,4,4)} could be used
#' to obtain per subject rather than per observation residuals.
#' @param weighted if \code{TRUE} and the model was fit with case weights, then
#' the weighted residuals are returned.
#' @param scaleX Should the \code{Xplan} columns be standardized ?
#' @param scaleY Should the \code{time} values be standardized ?
#' @param nt number of components to be extracted
#' @param limQ2set limit value for the Q2
#' @param dataPredictY predictor(s) (testing) dataset
#' @param pvals.expli should individual p-values be reported to tune model
#' selection ?
#' @param alpha.pvals.expli level of significance for predictors when
#' pvals.expli=TRUE
#' @param tol_Xi minimal value for Norm2(Xi) and \eqn{\mathrm{det}(pp' \times
#' pp)}{det(pp'*pp)} if there is any missing value in the \code{dataX}. It
#' defaults to \eqn{10^{-12}}{10^{-12}}
#' @param weights an optional vector of 'prior weights' to be used in the
#' fitting process. Should be \code{NULL} or a numeric vector.
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param allres FALSE to return only the Cox model and TRUE for additionnal
#' results. See details. Defaults to FALSE.
#' @param dataXplan an optional data frame, list or environment (or object
#' coercible by \code{\link{as.data.frame}} to a data frame) containing the
#' variables in the model. If not found in \code{dataXplan}, the variables are
#' taken from \code{environment(Xplan)}, typically the environment from which
#' \code{coxDKplsDR} is called.
#' @param model_frame If \code{TRUE}, the model frame is returned.
#' @param method the method to be used in fitting the model. The default method
#' \code{"glm.fit"} uses iteratively reweighted least squares (IWLS).
#' User-supplied fitting functions can be supplied either as a function or a
#' character string naming a function, with a function which takes the same
#' arguments as \code{glm.fit}.
#' @param control a list of parameters for controlling the fitting process. For
#' \code{glm.fit} this is passed to \code{\link{glm.control}}.
#' @param sparse should the coefficients of non-significant predictors
#' (<\code{alpha.pvals.expli}) be set to 0
#' @param sparseStop should component extraction stop when no significant
#' predictors (<\code{alpha.pvals.expli}) are found
#' @param verbose Should some details be displayed ?
#' @param model_matrix If \code{TRUE}, the model matrix is returned.
#' @param contrasts.arg a list, whose entries are values (numeric matrices, 
#' functions or character strings naming functions) to be used as replacement 
#' values for the contrasts replacement function and whose names are the names 
#' of columns of data containing factors.
#' @param \dots arguments to pass to \code{plsRmodel.default} or to
#' \code{plsRmodel.formula}
#' @return Depends on the model that was used to fit the model.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@lecnam.net}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link[plsRglm]{plsR}} and \code{\link[plsRglm]{plsRglm}}
#' @references plsRcox, Cox-Models in a high dimensional setting in R, Frederic
#' Bertrand, Philippe Bastien, Nicolas Meyer and Myriam Maumy-Bertrand (2014).
#' Proceedings of User2014!, Los Angeles, page 152.\cr
#' 
#' Deviance residuals-based sparse PLS and sparse kernel PLS regression for
#' censored data, Philippe Bastien, Frederic Bertrand, Nicolas Meyer and Myriam
#' Maumy-Bertrand (2015), Bioinformatics, 31(3):397-404,
#' doi:10.1093/bioinformatics/btu660.
#' @keywords models regression
#' @examples
#' 
#' data(micro.censure)
#' data(Xmicro.censure_compl_imp)
#' 
#' X_train_micro <- apply((as.matrix(Xmicro.censure_compl_imp)),FUN="as.numeric",MARGIN=2)[1:80,]
#' X_train_micro_df <- data.frame(X_train_micro)
#' Y_train_micro <- micro.censure$survyear[1:80]
#' C_train_micro <- micro.censure$DC[1:80]
#' 
#' plsRcox(X_train_micro,time=Y_train_micro,event=C_train_micro,nt=5)
#' plsRcox(~X_train_micro,time=Y_train_micro,event=C_train_micro,nt=5)
#' 
#' plsRcox(Xplan=X_train_micro,time=Y_train_micro,event=C_train_micro,nt=5,sparse=TRUE,
#' alpha.pvals.expli=.15)
#' plsRcox(Xplan=~X_train_micro,time=Y_train_micro,event=C_train_micro,nt=5,sparse=TRUE,
#' alpha.pvals.expli=.15)
#' 
#' @export plsRcox
plsRcox <- function (Xplan, ...) UseMethod("plsRcoxmodel")

#' @rdname plsRcox
#' @aliases plsRcox
#' @export plsRcoxmodel
plsRcoxmodel <- plsRcox
