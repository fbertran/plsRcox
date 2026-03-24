#' Print method for plsRcox models
#' 
#' This function provides a predict method for the class \code{"plsRcoxmodel"}
#' 
#' 
#' @param object An object of the class \code{"plsRcoxmodel"}.
#' @param newdata An optional data frame in which to look for variables with
#' which to predict. If omitted, the fitted values are used. For
#' \code{type = "expected"} on genuinely new observations, supply the
#' associated survival response with \code{y}; when \code{newdata} is a subset
#' of the training rows and keeps compatible row names, the training response
#' is matched automatically.
#' @param comps A value with a single value of component to use for prediction.
#' @param type Type of predicted value. Choices are the linear predictor
#' ("\code{lp}"), the risk score exp(lp) ("\code{risk}"), the expected number
#' of events given the covariates and follow-up time ("\code{expected}"), the
#' terms of the linear predictor ("\code{terms}") or the scores
#' ("\code{scores}").
#' @param se.fit If TRUE, pointwise standard errors are produced for the
#' predictions using the Cox model.
#' @param reference Reference level used to center relative predictions. This
#' is passed to \code{\link[survival]{predict.coxph}} and affects
#' \code{type = "lp"}, \code{"risk"} and \code{"terms"}.
#' @param y Optional \code{\link[survival]{Surv}} response to use with
#' \code{type = "expected"} when predicting on new observations. The time
#' values are interpreted on the original scale used to fit the model and are
#' rescaled internally when needed.
#' @param weights Vector of case weights. If \code{weights} is a vector of
#' integers, then the estimated coefficients are equivalent to estimating the
#' model from data with the individual \code{cases} replicated as many times as
#' indicated by \code{weights}.
#' @param methodNA Selects the way of predicting the response or the scores of
#' the new data. For complete rows, without any missing value, there are two
#' different ways of computing the prediction. As a consequence, for mixed
#' datasets, with complete and incomplete rows, there are two ways of computing
#' prediction : either predicts any row as if there were missing values in it
#' (\code{missingdata}) or selects the prediction method accordingly to the
#' completeness of the row (\code{adaptative}).
#' @param verbose Should some details be displayed ?
#' @param \dots Additional arguments passed on to
#' \code{\link[survival]{predict.coxph}}.
#' @return When type is "\code{response}", a matrix of predicted response
#' values is returned.\cr When type is "\code{scores}", a score matrix is
#' returned.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@lecnam.net}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link[survival]{predict.coxph}}
#' @references plsRcox, Cox-Models in a high dimensional setting in R, Frederic
#' Bertrand, Philippe Bastien, Nicolas Meyer and Myriam Maumy-Bertrand (2014).
#' Proceedings of User2014!, Los Angeles, page 152.\cr
#' 
#' Deviance residuals-based sparse PLS and sparse kernel PLS regression for
#' censored data, Philippe Bastien, Frederic Bertrand, Nicolas Meyer and Myriam
#' Maumy-Bertrand (2015), Bioinformatics, 31(3):397-404,
#' doi:10.1093/bioinformatics/btu660.
#' @keywords methods predict
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
#' 
#' predict(modpls)    
#' #Identical to predict(modpls,type="lp")    
#' 
#' predict(modpls,type="risk")    
#' predict(modpls,type="expected")    
#' predict(modpls,type="terms")    
#' predict(modpls,type="scores")    
#' 
#' predict(modpls,se.fit=TRUE)    
#' #Identical to predict(modpls,type="lp")    
#' predict(modpls,type="risk",se.fit=TRUE)    
#' predict(modpls,type="expected",se.fit=TRUE)    
#' predict(modpls,type="terms",se.fit=TRUE)    
#' predict(modpls,type="scores",se.fit=TRUE)    
#' 
#' 
#' #Identical to predict(modpls,type="lp")    
#' predict(modpls,newdata=X_train_micro[1:5,],type="risk")    
#' predict(modpls,newdata=X_train_micro[1:5,],type="expected")    
#' predict(modpls,newdata=X_train_micro[1:5,],type="terms")    
#' predict(modpls,newdata=X_train_micro[1:5,],type="scores")    
#' 
#' #Identical to predict(modpls,type="lp")    
#' predict(modpls,newdata=X_train_micro[1:5,],type="risk",se.fit=TRUE)    
#' predict(modpls,newdata=X_train_micro[1:5,],type="expected",se.fit=TRUE)    
#' predict(modpls,newdata=X_train_micro[1:5,],type="terms",se.fit=TRUE)    
#' predict(modpls,newdata=X_train_micro[1:5,],type="scores")    
#' 
#' newY_micro <- survival::Surv(Y_train_micro[1:5], C_train_micro[1:5])
#' predict(modpls,newdata=unname(X_train_micro[1:5,]),type="expected",y=newY_micro)    
#' 
#' predict(modpls,newdata=X_train_micro[1:5,],type="risk",comps=1)    
#' predict(modpls,newdata=X_train_micro[1:5,],type="risk",comps=2)    
#' predict(modpls,newdata=X_train_micro[1:5,],type="risk",comps=3)    
#' try(predict(modpls,newdata=X_train_micro[1:5,],type="risk",comps=4))
#' 
#' predict(modpls,newdata=X_train_micro[1:5,],type="terms",comps=1)    
#' predict(modpls,newdata=X_train_micro[1:5,],type="terms",comps=2)    
#' predict(modpls,newdata=X_train_micro[1:5,],type="terms",comps=3)    
#' try(predict(modpls,newdata=X_train_micro[1:5,],type="terms",comps=4))
#' 
#' predict(modpls,newdata=X_train_micro[1:5,],type="scores",comps=1)    
#' predict(modpls,newdata=X_train_micro[1:5,],type="scores",comps=2)    
#' predict(modpls,newdata=X_train_micro[1:5,],type="scores",comps=3)    
#' try(predict(modpls,newdata=X_train_micro[1:5,],type="scores",comps=4))
#' 
##' @method predict plsRcoxmodel
##' @export
predict.plsRcoxmodel <- function(object,newdata,comps=object$computed_nt,type=c("lp", "risk", "expected", "terms", "scores"),se.fit=FALSE,reference=c("strata", "sample", "zero"),y=NULL,weights,methodNA="adaptative", verbose=TRUE,...)
{
  if (!inherits(object, "plsRcoxmodel")) 
    stop("Primary argument much be a plsRcoxmodel object")
  if(missing(type)){type="lp"}
  if (!(type %in% c("lp", "risk", "expected", "terms", "scores")))
    stop("Invalid type specification")
  if (comps>object$computed_nt){stop("Cannot predict using more components than extracted.")}
  type <- match.arg(type)
  reference <- match.arg(reference)

  build_prediction_frame <- function(scores, response = NULL) {
    scores <- as.matrix(scores)
    padding <- object$computed_nt - ncol(scores)
    if (padding > 0) {
      scores <- cbind(scores, matrix(0, nrow = nrow(scores), ncol = padding))
    }
    if (ncol(scores) > object$computed_nt) {
      scores <- scores[, 1:object$computed_nt, drop = FALSE]
    }
    colnames(scores) <- NULL
    if (is.null(response)) {
      data.frame(tt = scores)
    } else {
      data.frame(YwotNA = response, tt = scores)
    }
  }

  scale_expected_response <- function(response) {
    center <- attr(object$RepY, "scaled:center")
    scale <- attr(object$RepY, "scaled:scale")
    if (is.null(center) || is.null(scale)) {
      return(response)
    }
    if (isTRUE(all.equal(center, 0)) && isTRUE(all.equal(scale, 1))) {
      return(response)
    }

    response_type <- attr(response, "type")
    time_columns <- switch(response_type,
                           counting = 1:2,
                           interval2 = 1:2,
                           1)
    response[, time_columns] <- sweep(sweep(response[, time_columns, drop = FALSE], 2, center), 2, scale, "/")
    response
  }

  match_training_rows <- function(candidate) {
    training_data <- as.matrix(object$dataX)
    candidate <- as.matrix(candidate)

    if (ncol(candidate) != ncol(training_data)) {
      return(NULL)
    }
    if (!identical(colnames(candidate), colnames(training_data))) {
      return(NULL)
    }

    row_key <- function(x) {
      apply(x, 1, function(row_values) {
        paste(ifelse(is.na(row_values), "NA", as.character(row_values)), collapse = "\r")
      })
    }

    training_keys <- row_key(training_data)
    candidate_keys <- row_key(candidate)
    candidate_matches <- match(candidate_keys, training_keys)

    if (anyNA(candidate_matches)) {
      return(NULL)
    }

    duplicate_matches <- vapply(candidate_keys, function(row_value) sum(training_keys == row_value), integer(1))
    if (any(duplicate_matches > 1L)) {
      stop("type = 'expected' could not match newdata rows uniquely to the training response. Supply 'y = survival::Surv(...)' or keep unique training row names.")
    }

    candidate_matches
  }

  resolve_expected_response <- function(newdata_rows = NULL, score_rows = NULL, newdata_matrix = NULL) {
    if (!is.null(y)) {
      if (!inherits(y, "Surv")) {
        stop("'y' must inherit from 'Surv' when type = 'expected'.")
      }
      if (!is.null(score_rows) && NROW(y) != score_rows) {
        stop("'y' must have the same number of rows as the predicted newdata.")
      }
      return(scale_expected_response(y))
    }

    if (is.null(newdata_rows)) {
      return(object$RepY)
    }

    train_rows <- rownames(object$ExpliX)
    if (!is.null(train_rows) && !is.null(newdata_rows)) {
      row_index <- match(newdata_rows, train_rows)
      if (!anyNA(row_index)) {
        return(object$RepY[row_index])
      }
    }

    if (!is.null(newdata_matrix)) {
      row_index <- match_training_rows(newdata_matrix)
      if (!is.null(row_index)) {
        return(object$RepY[row_index])
      }
    }

    stop("type = 'expected' requires follow-up information for newdata. Supply 'y = survival::Surv(...)' or keep row names that match the training rows.")
  }

  if (missing(newdata) || is.null(newdata)) {
    scores <- object$tt[,1:comps,drop=FALSE]
    if (type=="scores"){return(scores)}
    if (type=="expected") {
      ttpred <- build_prediction_frame(scores, response = resolve_expected_response())
    } else {
      ttpred <- build_prediction_frame(scores)
    }
    return(predict(object$FinalModel,newdata=ttpred,type = type,se.fit=se.fit,reference=reference,...))
  } else {
    nrnd <- nrow(newdata)
    if(any(apply(is.na(newdata),MARGIN=1,"all"))){return(vector("list",0)); cat("One of the rows of newdata is completely filled with missing data\n"); stop()}
    if(any(is.na(newdata))) {na.miss.newdata <- TRUE} else {na.miss.newdata <- FALSE}
    if(!is.null(object$XplanFormula)){
    mf0 <- match.call(expand.dots = FALSE)
    m0 <- match(c("weights"), names(mf0), 0L)
    mf0 <- mf0[c(1L, m0)]
    mf0$data <- newdata
    mf0$formula <- object$XplanFormula
    mf0$drop.unused.levels <- TRUE
    mf0[[1L]] <- as.name("model.frame")
    mf0 <- eval(mf0, parent.frame())
    mt0 <- attr(mf0, "terms")
    attr(mt0,"intercept")<-0L
    newdata.frame <- if (!is.empty.model(mt0)) model.matrix(mt0, mf0, contrasts)[,-1]
    else matrix(, nrnd, 0L)
    weights <- as.vector(model.weights(mf0))
    } else {newdata.frame <- newdata}
    newdata.scaled <- sweep(sweep(newdata.frame, 2, attr(object$ExpliX,"scaled:center")), 2 ,attr(object$ExpliX,"scaled:scale"), "/")
  newdataNA <- !is.na(newdata)
  newdata.scaledwotNA <- as.matrix(newdata.scaled)
  newdata.scaledwotNA[!newdataNA] <- 0
  ttpredY <- NULL
  if (methodNA=="adaptative") {
    for(ii in 1:nrnd){
      if (all(newdataNA[ii,])){
        ttpredY <- rbind(ttpredY, c(newdata.scaledwotNA[ii,]%*%object$wwetoile[,1:comps],rep(0,object$computed_nt-comps))) 
      }
      else {
        if(verbose){cat("Missing value in row ",ii,".\n")}
        ttpredY <- rbind(ttpredY, c(t(solve(t(object$pp[newdataNA[ii,],,drop=FALSE])%*%object$pp[newdataNA[ii,],,drop=FALSE])%*%t(object$pp[newdataNA[ii,],,drop=FALSE])%*%(newdata.scaledwotNA[ii,])[newdataNA[ii,]])[1:comps],rep(0,object$computed_nt-comps)))
      }}}
  if(methodNA=="missingdata") {
    if(verbose){cat("Prediction as if missing values in every row.\n")}
    for (ii in 1:nrnd) {  
      ttpredY <- rbind(ttpredY, c(t(solve(t(object$pp[newdataNA[ii,],,drop=FALSE])%*%object$pp[newdataNA[ii,],,drop=FALSE])%*%t(object$pp[newdataNA[ii,],,drop=FALSE])%*%(newdata.scaledwotNA[ii,])[newdataNA[ii,]])[1:comps],rep(0,object$computed_nt-comps)))
    }
  }
    if (type=="scores"){
      colnames(ttpredY) <- paste("Comp_",1:object$computed_nt,sep="")
      return(ttpredY[,1:comps,drop=FALSE])
    }
    if (type=="expected") {
      ttpred <- build_prediction_frame(ttpredY, response = resolve_expected_response(newdata_rows = rownames(newdata.frame), score_rows = nrnd, newdata_matrix = newdata.frame))
    } else {
      ttpred <- build_prediction_frame(ttpredY)
    }
    return(predict(object$FinalModel,newdata=ttpred,type = type,se.fit=se.fit,reference=reference,...))
}
}
