\name{predict.plsRcoxmodel}
\alias{predict.plsRcoxmodel}
\title{Print method for plsRcox models}
\description{
This function provides a predict method for the class \code{"plsRcoxmodel"}
}
\usage{
\method{predict}{plsRcoxmodel}(object,newdata,comps=object$computed_nt,
type=c("lp", "risk", "expected", "terms", "scores"),se.fit=FALSE,
weights,methodNA="adaptative",verbose=TRUE,...)
}
\arguments{
  \item{object}{An object of the class \code{"plsRcoxmodel"}.}
  \item{newdata}{An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.}
  \item{comps}{A value with a single value of component to use for prediction.}
  \item{type}{Type of predicted value. Choices are the linear predictor ("\code{lp}"), the risk score exp(lp) ("\code{risk}"), the expected number of events given the covariates and follow-up time ("\code{expected}"), the terms of the linear predictor ("\code{terms}") or the scores ("\code{scores}").}
  \item{se.fit}{If TRUE, pointwise standard errors are produced for the predictions using the Cox model.}
  \item{weights}{Vector of case weights. If \code{weights} is a vector of integers, then the estimated coefficients are equivalent to estimating the model from data with the individual \code{cases} replicated as many times as indicated by \code{weights}.}
  \item{methodNA}{Selects the way of predicting the response or the scores of the new data. For complete rows, without any missing value, there are two different ways of computing the prediction. As a consequence, for mixed datasets, with complete and incomplete rows, there are two ways of computing prediction : either predicts any row as if there were missing values in it (\code{missingdata}) or selects the prediction method accordingly to the completeness of the row (\code{adaptative}).}
  \item{verbose}{Should some details be displayed ?}
  \item{\dots}{Arguments to be passed on to \code{survival::coxph} and to \code{plsRglm::PLS_lm}.}
}
\value{
When type is "\code{response}", a matrix of predicted response values is returned.\cr
When type is "\code{scores}", a score matrix is returned. 
}
\references{\enc{Fr?d?ric}{Fr\'ed\'eric} Bertrand, Myriam Maumy-Bertrand et Nicolas Meyer (2011). \enc{R?gression}{R\'egression} \enc{B?ta}{B\^eta} PLS. \emph{Preprint}.
}
\author{\enc{Fr?d?ric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{\code{\link[survival]{predict.coxph}}}
\examples{
data(micro.censure)
data(Xmicro.censure_compl_imp)

X_train_micro <- apply((as.matrix(Xmicro.censure_compl_imp)),FUN="as.numeric",MARGIN=2)[1:80,]
Y_train_micro <- micro.censure$survyear[1:80]
C_train_micro <- micro.censure$DC[1:80]

modpls <- plsRcox(X_train_micro,time=Y_train_micro,event=C_train_micro,nt=3)

predict(modpls)    
#Identical to predict(modpls,type="lp")    

predict(modpls,type="risk")    
predict(modpls,type="expected")    
predict(modpls,type="terms")    
predict(modpls,type="scores")    

predict(modpls,se.fit=TRUE)    
#Identical to predict(modpls,type="lp")    
predict(modpls,type="risk",se.fit=TRUE)    
predict(modpls,type="expected",se.fit=TRUE)    
predict(modpls,type="terms",se.fit=TRUE)    
predict(modpls,type="scores",se.fit=TRUE)    


#Identical to predict(modpls,type="lp")    
predict(modpls,newdata=X_train_micro[1:5,],type="risk")    
#predict(modpls,newdata=X_train_micro[1:5,],type="expected")    
predict(modpls,newdata=X_train_micro[1:5,],type="terms")    
predict(modpls,newdata=X_train_micro[1:5,],type="scores")    

#Identical to predict(modpls,type="lp")    
predict(modpls,newdata=X_train_micro[1:5,],type="risk",se.fit=TRUE)    
#predict(modpls,newdata=X_train_micro[1:5,],type="expected",se.fit=TRUE)    
predict(modpls,newdata=X_train_micro[1:5,],type="terms",se.fit=TRUE)    
predict(modpls,newdata=X_train_micro[1:5,],type="scores")    

predict(modpls,newdata=X_train_micro[1:5,],type="risk",comps=1)    
predict(modpls,newdata=X_train_micro[1:5,],type="risk",comps=2)    
predict(modpls,newdata=X_train_micro[1:5,],type="risk",comps=3)    
try(predict(modpls,newdata=X_train_micro[1:5,],type="risk",comps=4))

predict(modpls,newdata=X_train_micro[1:5,],type="terms",comps=1)    
predict(modpls,newdata=X_train_micro[1:5,],type="terms",comps=2)    
predict(modpls,newdata=X_train_micro[1:5,],type="terms",comps=3)    
try(predict(modpls,newdata=X_train_micro[1:5,],type="terms",comps=4))

predict(modpls,newdata=X_train_micro[1:5,],type="scores",comps=1)    
predict(modpls,newdata=X_train_micro[1:5,],type="scores",comps=2)    
predict(modpls,newdata=X_train_micro[1:5,],type="scores",comps=3)    
try(predict(modpls,newdata=X_train_micro[1:5,],type="scores",comps=4))
}
\keyword{methods}
\keyword{predict}

