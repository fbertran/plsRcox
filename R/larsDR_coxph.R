#' Fitting a LASSO/LARS model on the (Deviance) Residuals
#' 
#' This function computes the Cox Model based on lars variables computed model
#' with \itemize{ \item as the response: the Residuals of a Cox-Model fitted
#' with no covariate \item as explanatory variables: Xplan.  } It uses the
#' package \code{lars} to perform PLSR fit.
#' 
#' This function computes the LASSO/LARS model with the Residuals of a
#' Cox-Model fitted with an intercept as the only explanatory variable as the
#' response and Xplan as explanatory variables. Default behaviour uses the
#' Deviance residuals.
#' 
#' If \code{allres=FALSE} returns only the final Cox-model. If
#' \code{allres=TRUE} returns a list with the (Deviance) Residuals, the
#' LASSO/LARS model fitted to the (Deviance) Residuals, the eXplanatory
#' variables and the final Cox-model. \code{allres=TRUE} is useful for
#' evluating model prediction accuracy on a test sample.
#' 
#' @aliases larsDR_coxph larsDR_coxph.default larsDR_coxph.formula
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
#' @param plot Should the survival function be plotted ?)
#' @param typelars One of \code{"lasso"}, \code{"lar"},
#' \code{"forward.stagewise"} or \code{"stepwise"}. The names can be
#' abbreviated to any unique substring. Default is \code{"lasso"}.
#' @param normalize If TRUE, each variable is standardized to have unit L2
#' norm, otherwise it is left alone. Default is TRUE.
#' @param max.steps Limit the number of steps taken; the default is \code{8 *
#' min(m, n-intercept)}, with m the number of variables, and n the number of
#' samples. For \code{type="lar"} or \code{type="stepwise"}, the maximum number
#' of steps is \code{min(m,n-intercept)}. For \code{type="lasso"} and
#' especially \code{type="forward.stagewise"}, there can be many more terms,
#' because although no more than min(m,n-intercept) variables can be active
#' during any step, variables are frequently droppped and added as the
#' algorithm proceeds. Although the default usually guarantees that the
#' algorithm has proceeded to the saturated fit, users should check.
#' @param use.Gram When the number m of variables is very large, i.e. larger
#' than N, then you may not want LARS to precompute the Gram matrix. Default is
#' \code{use.Gram=TRUE}
#' @param allres FALSE to return only the Cox model and TRUE for additionnal
#' results. See details. Defaults to FALSE.
#' @param dataXplan an optional data frame, list or environment (or object
#' coercible by \code{\link{as.data.frame}} to a data frame) containing the
#' variables in the model. If not found in \code{dataXplan}, the variables are
#' taken from \code{environment(Xplan)}, typically the environment from which
#' \code{plscox} is called.
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param weights an optional vector of 'prior weights' to be used in the
#' fitting process. Should be \code{NULL} or a numeric vector.
#' @param model_frame If \code{TRUE}, the model frame is returned.
#' @param model_matrix If \code{TRUE}, the "unweighted" model matrix is
#' returned.
#' @param verbose Should some details be displayed ?
#' @param model_matrix If \code{TRUE}, the model matrix is returned.
#' @param contrasts.arg a list, whose entries are values (numeric matrices, 
#' functions or character strings naming functions) to be used as replacement 
#' values for the contrasts replacement function and whose names are the names 
#' of columns of data containing factors.
#' @param \dots Arguments to be passed on to \code{survival::coxph} or to
#' \code{lars::lars}.
#' @return If \code{allres=FALSE} : \item{cox_larsDR}{Final Cox-model.} If
#' \code{allres=TRUE} : \item{DR_coxph}{The (Deviance) Residuals.}
#' \item{larsDR}{The LASSO/LARS model fitted to the (Deviance) Residuals.}
#' \item{X_larsDR}{The eXplanatory variables.} \item{cox_larsDR}{Final
#' Cox-model.}
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@lecnam.net}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link[survival]{coxph}}, \code{\link[lars]{lars}}
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
#' (cox_larsDR_fit <- larsDR_coxph(X_train_micro,Y_train_micro,C_train_micro,max.steps=6,
#' use.Gram=FALSE,scaleX=TRUE))
#' (cox_larsDR_fit <- larsDR_coxph(~X_train_micro,Y_train_micro,C_train_micro,max.steps=6,
#' use.Gram=FALSE,scaleX=TRUE))
#' (cox_larsDR_fit <- larsDR_coxph(~.,Y_train_micro,C_train_micro,max.steps=6,
#' use.Gram=FALSE,scaleX=TRUE,dataXplan=X_train_micro_df))
#' 
#' larsDR_coxph(~X_train_micro,Y_train_micro,C_train_micro,max.steps=6,use.Gram=FALSE)
#' larsDR_coxph(~X_train_micro,Y_train_micro,C_train_micro,max.steps=6,use.Gram=FALSE,scaleX=FALSE)
#' larsDR_coxph(~X_train_micro,Y_train_micro,C_train_micro,max.steps=6,use.Gram=FALSE,
#' scaleX=TRUE,allres=TRUE)
#' 
#' rm(X_train_micro,Y_train_micro,C_train_micro,cox_larsDR_fit)
#' 
#' @export larsDR_coxph
larsDR_coxph <- function (Xplan, ...) UseMethod("larsDR_coxph")
