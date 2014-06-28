\name{coxpls3}
\alias{coxpls3}
\alias{coxpls3.default}
\alias{coxpls3.formula}
\title{Fitting a Cox-Model on PLSR components}
\description{
This function computes the the Cox-Model with PLSR components as the explanatory variables. It uses the package \code{plsRglm}.
}
\usage{
coxpls3(Xplan, \dots)
\method{coxpls3}{default}(Xplan,time,time2,event,type,
origin,typeres="deviance", collapse, weighted, scaleX=TRUE, 
scaleY=TRUE, nt=min(7,ncol(Xplan)), typeVC="none", 
plot=FALSE, allres=FALSE,sparse=FALSE,sparseStop=TRUE,...)
\method{coxpls3}{formula}(Xplan,time,time2,event,type,
origin,typeres="deviance", collapse, weighted, scaleX=TRUE, 
scaleY=TRUE, nt=min(7,ncol(Xplan)), typeVC="none", 
plot=FALSE, allres=FALSE,dataXplan=NULL,subset,weights,
model_frame=FALSE,sparse=FALSE,sparseStop=TRUE,...)
}
\arguments{
  \item{Xplan}{a formula or a matrix with the eXplanatory variables (training) dataset}
  \item{time}{for right censored data, this is the follow up time. For interval data, the first argument is the starting time for the interval.}
  \item{time2}{The status indicator, normally 0=alive, 1=dead. Other choices are \code{TRUE/FALSE} (\code{TRUE} = death) or 1/2 (2=death). For interval censored data, the status indicator is 0=right censored, 1=event at \code{time}, 2=left censored, 3=interval censored. Although unusual, the event indicator can be omitted, in which case all subjects are assumed to have an event.}
  \item{event}{ending time of the interval for interval censored or counting process data only. Intervals are assumed to be open on the left and closed on the right, \code{(start, end]}. For counting process data, event indicates whether an event occurred at the end of the interval.}
  \item{type}{character string specifying the type of censoring. Possible values are \code{"right"}, \code{"left"}, \code{"counting"}, \code{"interval"}, or \code{"interval2"}. The default is \code{"right"} or \code{"counting"} depending on whether the \code{time2} argument is absent or present, respectively.}
  \item{origin}{for counting process data, the hazard function origin. This option was intended to be used in conjunction with a model containing time dependent strata in order to align the subjects properly when they cross over from one strata to another, but it has rarely proven useful.}
  \item{typeres}{character string indicating the type of residual desired. Possible values are \code{"martingale"}, \code{"deviance"}, \code{"score"}, \code{"schoenfeld"}, \code{"dfbeta"}, \code{"dfbetas"}, and \code{"scaledsch"}. Only enough of the string to determine a unique match is required.}
  \item{collapse}{vector indicating which rows to collapse (sum) over. In time-dependent models more than one row data can pertain to a single individual. If there were 4 individuals represented by 3, 1, 2 and 4 rows of data respectively, then \code{collapse=c(1,1,1,2,3,3,4,4,4,4)} could be used to obtain per subject rather than per observation residuals.}
  \item{weighted}{if \code{TRUE} and the model was fit with case weights, then the weighted residuals are returned.}
  \item{scaleX}{Should the \code{Xplan} columns be standardized ?}
  \item{scaleY}{Should the \code{time} values be standardized ?}
  \item{nt}{Number of PLSR components to fit.}
  \item{typeVC}{type of leave one out crossed validation. Several procedures are available and may be forced.
    \describe{
      \item{\code{none}}{no crossed validation}
      \item{\code{standard}}{as in SIMCA for datasets without missing values and with all values predicted as those with missing values for datasets with any missing values}
      \item{\code{missingdata}}{all values predicted as those with missing values for datasets with any missing values}
      \item{\code{adaptative}}{predict a response value for an x with any missing value as those with missing values and for an x without any missing value as those without missing values.}
      }
  }
  \item{plot}{Should the survival function be plotted ?)}
  \item{allres}{FALSE to return only the Cox model and TRUE for additionnal results. See details. Defaults to FALSE.}
  \item{dataXplan}{an optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model. If not found in \code{dataXplan}, the variables are taken from \code{environment(Xplan)}, typically the environment from which \code{coxpls3} is called.}  
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
  \item{weights}{an optional vector of 'prior weights' to be used in the fitting process. Should be \code{NULL} or a numeric vector.}
  \item{model_frame}{If \code{TRUE}, the model frame is returned.}
  \item{sparse}{should the coefficients of non-significant predictors (<\code{alpha.pvals.expli}) be set to 0}
  \item{sparseStop}{should component extraction stop when no significant predictors (<\code{alpha.pvals.expli}) are found}
  \item{\dots}{Arguments to be passed on to \code{survival::coxph} and to \code{plsRglm::PLS_lm}.}
}
\details{
If \code{allres=FALSE} returns only the final Cox-model.
If \code{allres=TRUE} returns a list with the PLS components, the final Cox-model and the PLSR model.
\code{allres=TRUE} is useful for evluating model prediction accuracy on a test sample.
}
\value{If \code{allres=FALSE} : 
\item{cox_pls3}{Final Cox-model.}
If \code{allres=TRUE} :
\item{tt_pls3}{PLSR components.}
\item{cox_pls3}{Final Cox-model.}
\item{pls3_mod}{The PLSR model.}
}
\references{plsRcox : \enc{mod?les}{mod\`eles} de Cox en \enc{pr?sence}{pr\'esence} d'un grand nombre de variables explicatives, \enc{Fr?d?ric}{Fr\'ed\'eric} Bertrand, Myriam Maumy-Bertrand, Marie-Pierre Gaub, Nicolas Meyer, \enc{Chimiom?trie}{Chimiom\'etrie} 2010, Paris, 2010.}
\author{\enc{Fr?d?ric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{\code{\link[survival]{coxph}}, \code{\link[plsRglm]{PLS_lm}}}
\examples{
data(micro.censure)
data(Xmicro.censure_compl_imp)

X_train_micro <- apply((as.matrix(Xmicro.censure_compl_imp)),FUN="as.numeric",MARGIN=2)[1:80,]
X_train_micro_df <- data.frame(X_train_micro)
Y_train_micro <- micro.censure$survyear[1:80]
C_train_micro <- micro.censure$DC[1:80]

(cox_pls3_fit <- coxpls3(X_train_micro,Y_train_micro,C_train_micro,nt=7,typeVC="none"))
(cox_pls3_fit2 <- coxpls3(~X_train_micro,Y_train_micro,C_train_micro,nt=7,typeVC="none"))
(cox_pls3_fit3 <- coxpls3(~.,Y_train_micro,C_train_micro,nt=7,typeVC="none",data=X_train_micro_df))
(cox_pls3_fit4 <- coxpls3(~.,Y_train_micro,C_train_micro,nt=7,typeVC="none",
data=X_train_micro_df,sparse=TRUE))
(cox_pls3_fit5 <- coxpls3(~.,Y_train_micro,C_train_micro,nt=7,typeVC="none",
data=X_train_micro_df,sparse=FALSE,sparseStop=TRUE))

rm(X_train_micro,Y_train_micro,C_train_micro,cox_pls3_fit,cox_pls3_fit2,
cox_pls3_fit3,cox_pls3_fit4,cox_pls3_fit5)
}
\keyword{models}
\keyword{regression}
