predict.plsRcoxmodel <- function(object,newdata,comps=object$computed_nt,type=c("lp", "risk", "expected", "terms", "scores"),se.fit=FALSE,weights,methodNA="adaptative",...)
{
    if (!inherits(object, "plsRcoxmodel")) 
        stop("Primary argument much be a plsRcoxmodel object")
    if (comps>object$computed_nt){stop("Cannot predict using more components than extracted.")}
    type <- match.arg(type)
    if (missing(newdata) || is.null(newdata)) {
    nrtt <- nrow(object$tt)
    if (type=="lp"){tt<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(tt) <- paste("Comp_",1:object$computed_nt,sep="")
    return(predict(object$FinalModel,newdata=tt,type = "lp",se.fit=se.fit))
    }
    if (type=="risk"){tt<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(tt) <- paste("Comp_",1:object$computed_nt,sep="")
    return(predict(object$FinalModel,newdata=tt,type = "risk",se.fit=se.fit))
    }
    if (type=="expected"){tt<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(tt) <- paste("Comp_",1:object$computed_nt,sep="")
    return(predict(object$FinalModel,newdata=tt,type = "expected",se.fit=se.fit))
    }
    if (type=="terms"){tt<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(tt) <- paste("Comp_",1:object$computed_nt,sep="")
    return(predict(object$FinalModel,newdata=tt,type = "terms",se.fit=se.fit))
    }  
    if (type=="scores"){return(object$tt[,1:comps])}   
    } else {
    nrnd <- nrow(newdata)
    if(!is.null(object$XplanFormula)){
    datanewdata <- environment(newdata)
    mf0 <- match.call(expand.dots = FALSE)
    m0 <- match(c("weights"), names(mf0), 0L)
    mf0 <- mf0[c(1L, m0)]
    mf0$data <- datanewdata
    mf0$formula <- object$XplanFormula
    mf0$drop.unused.levels <- TRUE
    mf0[[1L]] <- as.name("model.frame")
    mf0 <- eval(mf0, parent.frame())
    mt0 <- attr(mf0, "terms")
    Y <- model.response(mf0, "any")
    if (length(dim(Y)) == 1L) {
      nm <- rownames(Y)
      dim(Y) <- NULL
      if (!is.null(nm)) 
        names(Y) <- nm
    }
    newdata.frame <- if (!is.empty.model(mt0)) model.matrix(mt0, mf0, contrasts)[,-1]
    else matrix(, NROW(Y), 0L)
    weights <- as.vector(model.weights(mf0))
    } else {newdata.frame <- newdata}
    newdata.scaled <- sweep(sweep(newdata.frame, 2, attr(object$ExpliX,"scaled:center")), 2 ,attr(object$ExpliX,"scaled:scale"), "/")}
    newdataNA <- is.na(newdata)    
    newdata.scaledwotNA <- newdata.scaled
    newdata.scaledwotNA[newdataNA] <- 0
    tt <- NULL
    if (methodNA=="adaptative") {
    for(ii in 1:nrnd){
    if (!any(newdataNA[ii,])){
    tt <- rbind(tt, c(newdata.scaledwotNA[ii,]%*%object$wwetoile[,1:comps],rep(0,object$computed_nt-comps))) 
    }
    else {
    cat("Missing value in row ",ii,".")
    tt <- rbind(tt, c(t(solve(t(object$pp[newdataNA[ii,],,drop=FALSE])%*%object$pp[newdataNA[ii,],,drop=FALSE])%*%t(object$pp[newdataNA[ii,],,drop=FALSE])%*%(newdata.scaledwotNA[ii,])[newdataNA[ii,]])[1:comps],rep(0,object$computed_nt-comps)))
    }}}
    if(methodNA=="missingdata") {
    cat("Prediction as if missing values in every row.")
    for (ii in 1:nrnd) {  
    tt <- rbind(tt, c(t(solve(t(object$pp[newdataNA[ii,],,drop=FALSE])%*%object$pp[newdataNA[ii,],,drop=FALSE])%*%t(object$pp[newdataNA[ii,],,drop=FALSE])%*%(newdata.scaledwotNA[ii,])[newdataNA[ii,]])[1:comps],rep(0,object$computed_nt-comps)))
    }
    }
    colnames(tt) <- paste("Comp_",1:object$computed_nt,sep="")
    if (type=="lp"){return(predict(object$FinalModel,newdata=data.frame(tt),type = "lp",se.fit=se.fit))
    }
    if (type=="risk"){return(predict(object$FinalModel,newdata=data.frame(tt),type = "risk",se.fit=se.fit))
    }
    if (type=="expected"){return(predict(object$FinalModel,newdata=data.frame(tt),type = "expected",se.fit=se.fit))
    }
    if (type=="terms"){return(predict(object$FinalModel,newdata=data.frame(tt),type = "terms",se.fit=se.fit))
    }  
    if (type=="scores"){return(tt[,1:comps,drop=FALSE])}      
}
 
