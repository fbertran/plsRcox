\name{DKplsRcox}
\alias{DKplsRcox}
\alias{DKplsRcoxmodel.default}
\alias{DKplsRcoxmodel.formula}
\title{Partial least squares Regression generalized linear models}
\description{
This function implements an extension of Partial least squares Regression to Cox Models.
}
\usage{
DKplsRcox(Xplan, \dots)
\method{DKplsRcoxmodel}{default}(Xplan,time,time2,event,type,
origin,typeres="deviance", collapse, weighted, scaleX=TRUE, 
scaleY=TRUE, nt=min(2,ncol(Xplan)),limQ2set=.0975, 
dataPredictY=Xplan, pvals.expli=FALSE,alpha.pvals.expli=.05,
tol_Xi=10^(-12),weights,control, sparse=FALSE,
sparseStop=TRUE,plot=FALSE,allres=FALSE, kernel = "rbfdot", 
hyperkernel, verbose=TRUE,...)
\method{DKplsRcoxmodel}{formula}(Xplan,time,time2,event,type,
origin,typeres="deviance", collapse, weighted,scaleX=TRUE,
scaleY=NULL,dataXplan=NULL, nt=min(2,ncol(Xplan)),
limQ2set=.0975, dataPredictY=Xplan, pvals.expli=FALSE, 
model_frame=FALSE, alpha.pvals.expli=.05,tol_Xi=10^(-12),
weights,subset,control,sparse=FALSE,sparseStop=TRUE,
plot=FALSE,allres=FALSE, kernel = "rbfdot", hyperkernel, 
verbose=TRUE,...)
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
  \item{nt}{number of components to be extracted}
  \item{limQ2set}{limit value for the Q2}
  \item{dataPredictY}{predictor(s) (testing) dataset}
  \item{pvals.expli}{should individual p-values be reported to tune model selection ?}
  \item{alpha.pvals.expli}{level of significance for predictors when pvals.expli=TRUE}
  \item{tol_Xi}{minimal value for Norm2(Xi) and \eqn{\mathrm{det}(pp' \times pp)}{det(pp'*pp)} if there is any missing value in the \code{dataX}. It defaults to \eqn{10^{-12}}{10^{-12}}}
  \item{weights}{an optional vector of 'prior weights' to be used in the fitting process. Should be \code{NULL} or a numeric vector.}
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
  \item{plot}{Should the survival function be plotted ?)}  
  \item{allres}{FALSE to return only the Cox model and TRUE for additionnal results. See details. Defaults to FALSE.}
  \item{dataXplan}{an optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model. If not found in \code{dataXplan}, the variables are taken from \code{environment(Xplan)}, typically the environment from which \code{coxDKplsDR} is called.}  
  \item{model_frame}{If \code{TRUE}, the model frame is returned.}
  \item{method}{the method to be used in fitting the model. The default method \code{"glm.fit"} uses iteratively reweighted least squares (IWLS). User-supplied fitting functions can be supplied either as a function or a character string naming a function, with a function which takes the same arguments as \code{glm.fit}.}
  \item{control}{a list of parameters for controlling the fitting process. For \code{glm.fit} this is passed to \code{\link{glm.control}}.}
  \item{sparse}{should the coefficients of non-significant predictors (<\code{alpha.pvals.expli}) be set to 0}
  \item{sparseStop}{should component extraction stop when no significant predictors (<\code{alpha.pvals.expli}) are found}
  \item{kernel}{the kernel function used in training and predicting. This parameter can be set to any function, of class kernel, which computes the inner product in feature space between two vector arguments (see \link[kernlab]{kernels}). The \code{kernlab} package provides the most popular kernel functions which can be used by setting the kernel parameter to the following strings: 
  \describe{
\item{\code{rbfdot}}{Radial Basis kernel "Gaussian"}
\item{\code{polydot}}{Polynomial kernel}
\item{\code{vanilladot}}{Linear kernel} 
\item{\code{tanhdot}}{Hyperbolic tangent kernel} 
\item{\code{laplacedot}}{Laplacian kernel} 
\item{\code{besseldot}}{Bessel kernel} 
\item{\code{anovadot}}{ANOVA RBF kernel} 
\item{\code{splinedot}}{Spline kernel} 
  }}
  \item{hyperkernel}{the list of hyper-parameters (kernel parameters). This is a list which contains the parameters to be used with the kernel function. For valid parameters for existing kernels are : 
  \itemize{
\item \code{sigma}, inverse kernel width for the Radial Basis kernel function "rbfdot" and the Laplacian kernel "laplacedot". 
\item \code{degree}, \code{scale}, \code{offset} for the Polynomial kernel "polydot". 
\item \code{scale}, offset for the Hyperbolic tangent kernel function "tanhdot". 
\item \code{sigma}, \code{order}, \code{degree} for the Bessel kernel "besseldot". 
\item \code{sigma}, \code{degree} for the ANOVA kernel "anovadot". 
}
 In the case of a Radial Basis kernel function (Gaussian) or Laplacian kernel, if \code{hyperkernel} is missing, the heuristics in sigest are used to calculate a good sigma value from the data.
} 
  \item{verbose}{Should some details be displayed ?}
  \item{\dots}{arguments to pass to \code{plsRmodel.default} or to \code{plsRmodel.formula}}
}
\details{

A typical predictor has the form response ~ terms where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response. A terms specification of the form first + second indicates all the terms in first together with all the terms in second with any duplicates removed. 

A specification of the form first:second indicates the the set of terms obtained by taking the interactions of all terms in first with all terms in second. The specification first*second indicates the cross of first and second. This is the same as first + second + first:second. 

The terms in the formula will be re-ordered so that main effects come first, followed by the interactions, all second-order, all third-order and so on: to avoid this pass a terms object as the formula. 

Non-NULL weights can be used to indicate that different observations have different dispersions (with the values in weights being inversely proportional to the dispersions); or equivalently, when the elements of weights are positive integers w_i, that each response y_i is the mean of w_i unit-weight observations. 
}
\value{
Depends on the model that was used to fit the model.
}
\references{\enc{Frederic}{Fr\'ed\'eric} Bertrand, Myriam Maumy-Bertrand et Nicolas Meyer (2011). \enc{R?gression}{R\'egression} \enc{B?ta}{B\^eta} PLS. \emph{Preprint}.
}
\author{\enc{Frederic}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{\code{\link[plsRglm]{plsR}} and \code{\link[plsRglm]{plsRglm}}}
\examples{
data(micro.censure)
data(Xmicro.censure_compl_imp)

X_train_micro <- apply((as.matrix(Xmicro.censure_compl_imp)),FUN="as.numeric",MARGIN=2)[1:80,]
X_train_micro_df <- data.frame(X_train_micro)
Y_train_micro <- micro.censure$survyear[1:80]
C_train_micro <- micro.censure$DC[1:80]

\donttest{
DKplsRcox(X_train_micro,time=Y_train_micro,event=C_train_micro,nt=5)
DKplsRcox(~X_train_micro,time=Y_train_micro,event=C_train_micro,nt=5)
}

\donttest{
DKplsRcox(Xplan=X_train_micro,time=Y_train_micro,event=C_train_micro,nt=5,sparse=TRUE, 
alpha.pvals.expli=.15)
DKplsRcox(Xplan=~X_train_micro,time=Y_train_micro,event=C_train_micro,nt=5,sparse=TRUE,
alpha.pvals.expli=.15)
}
}\keyword{models}
\keyword{regression}
