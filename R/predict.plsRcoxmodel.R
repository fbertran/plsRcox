#' Print method for plsRcox models
#' 
#' This function provides a predict method for the class \code{"plsRcoxmodel"}
#' 
#' 
#' @param object An object of the class \code{"plsRcoxmodel"}.
#' @param newdata An optional data frame in which to look for variables with
#' which to predict. If omitted, the fitted values are used.
#' @param comps A value with a single value of component to use for prediction.
#' @param type Type of predicted value. Choices are the linear predictor
#' ("\code{lp}"), the risk score exp(lp) ("\code{risk}"), the expected number
#' of events given the covariates and follow-up time ("\code{expected}"), the
#' terms of the linear predictor ("\code{terms}") or the scores
#' ("\code{scores}").
#' @param se.fit If TRUE, pointwise standard errors are produced for the
#' predictions using the Cox model.
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
#' @param \dots Arguments to be passed on to \code{survival::coxph} and to
#' \code{plsRglm::PLS_lm}.
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
#' #predict(modpls,newdata=X_train_micro[1:5,],type="expected")    
#' predict(modpls,newdata=X_train_micro[1:5,],type="terms")    
#' predict(modpls,newdata=X_train_micro[1:5,],type="scores")    
#' 
#' #Identical to predict(modpls,type="lp")    
#' predict(modpls,newdata=X_train_micro[1:5,],type="risk",se.fit=TRUE)    
#' #predict(modpls,newdata=X_train_micro[1:5,],type="expected",se.fit=TRUE)    
#' predict(modpls,newdata=X_train_micro[1:5,],type="terms",se.fit=TRUE)    
#' predict(modpls,newdata=X_train_micro[1:5,],type="scores")    
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
predict.plsRcoxmodel <- function(object,newdata,comps=object$computed_nt,type=c("lp", "risk", "expected", "terms", "scores"),se.fit=FALSE,weights,methodNA="adaptative", verbose=TRUE,...)
{
  if (!inherits(object, "plsRcoxmodel")) 
    stop("Primary argument much be a plsRcoxmodel object")
  if(missing(type)){type="lp"}
  if (!(type %in% c("lp", "risk", "expected", "terms", "scores")))
    stop("Invalid type specification")
  if (comps>object$computed_nt){stop("Cannot predict using more components than extracted.")}
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata)) {
    nrtt <- nrow(object$tt)
    if (type=="lp"){ttpredY<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(ttpredY) <- NULL
    ttpred<-data.frame(tt=ttpredY);return(predict(object$FinalModel,newdata=ttpred,type = "lp",se.fit=se.fit))
    }
    if (type=="risk"){ttpredY<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(ttpredY) <- NULL
    ttpred<-data.frame(tt=ttpredY);return(predict(object$FinalModel,newdata=ttpred,type = "risk",se.fit=se.fit))
    }
    if (type=="expected"){ttpredY<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(ttpredY) <- NULL
    ttpred<-data.frame(tt=ttpredY);return(predict(object$FinalModel,newdata=ttpred,type = "expected",se.fit=se.fit))
    }
    if (type=="terms"){ttpredY<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(ttpredY) <- NULL
    ttpred<-data.frame(tt=ttpredY);return(predict(object$FinalModel,newdata=ttpred,type = "terms",se.fit=se.fit))
    }  
    if (type=="scores"){return(object$tt[,1:comps])}   
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
colnames(ttpredY) <- NULL
ttpred<-data.frame(tt=ttpredY)
    if (type=="lp"){return(predict(object$FinalModel,newdata=ttpred,type = "lp",se.fit=se.fit))
    }
    if (type=="risk"){return(predict(object$FinalModel,newdata=ttpred,type = "risk",se.fit=se.fit))
    }
    if (type=="expected"){return(predict(object$FinalModel,newdata=ttpred,type = "expected",se.fit=se.fit))
    }
    if (type=="terms"){return(predict(object$FinalModel,newdata=ttpred,type = "terms",se.fit=se.fit))
    }  

    if (type=="scores"){colnames(ttpredY) <- paste("Comp_",1:object$computed_nt,sep="");
                      return(ttpredY[,1:comps,drop=FALSE])}
}
}

