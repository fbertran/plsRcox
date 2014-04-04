\name{DR_coxph}
\alias{DR_coxph}
\title{(Deviance) Residuals Computation}
\description{
This function computes the Residuals for a Cox-Model fitted with an intercept as the only explanatory variable. Default behaviour gives the Deviance residuals.
}
\usage{
DR_coxph(time, time2, event, type, origin, typeres = "deviance",
collapse, weighted, scaleY = TRUE, plot = FALSE, ...)
}
\arguments{
  \item{time}{for right censored data, this is the follow up time. For interval data, the first argument is the starting time for the interval.}
  \item{time2}{The status indicator, normally 0=alive, 1=dead. Other choices are \code{TRUE/FALSE} (\code{TRUE} = death) or 1/2 (2=death). For interval censored data, the status indicator is 0=right censored, 1=event at \code{time}, 2=left censored, 3=interval censored. Although unusual, the event indicator can be omitted, in which case all subjects are assumed to have an event.}
  \item{event}{ending time of the interval for interval censored or counting process data only. Intervals are assumed to be open on the left and closed on the right, \code{(start, end]}. For counting process data, event indicates whether an event occurred at the end of the interval.}
  \item{type}{character string specifying the type of censoring. Possible values are \code{"right"}, \code{"left"}, \code{"counting"}, \code{"interval"}, or \code{"interval2"}. The default is \code{"right"} or \code{"counting"} depending on whether the \code{time2} argument is absent or present, respectively.}
  \item{origin}{for counting process data, the hazard function origin. This option was intended to be used in conjunction with a model containing time dependent strata in order to align the subjects properly when they cross over from one strata to another, but it has rarely proven useful.}
  \item{typeres}{character string indicating the type of residual desired. Possible values are \code{"martingale"}, \code{"deviance"}, \code{"score"}, \code{"schoenfeld"}, \code{"dfbeta"}, \code{"dfbetas"}, and \code{"scaledsch"}. Only enough of the string to determine a unique match is required.}
  \item{collapse}{vector indicating which rows to collapse (sum) over. In time-dependent models more than one row data can pertain to a single individual. If there were 4 individuals represented by 3, 1, 2 and 4 rows of data respectively, then \code{collapse=c(1,1,1,2,3,3,4,4,4,4)} could be used to obtain per subject rather than per observation residuals.}
  \item{weighted}{if \code{TRUE} and the model was fit with case weights, then the weighted residuals are returned.}
  \item{scaleY}{Should the \code{time} values be standardized ?}
  \item{plot}{Should the survival function be plotted ?)}
  \item{\dots}{Arguments to be passed on to \code{survival::coxph}.}
}
\value{\item{Named num}{Vector of the residual values.}
}
\references{plsRcox : \enc{mod?les}{mod\`eles} de Cox en \enc{pr?sence}{pr\'esence} d'un grand nombre de variables explicatives, \enc{Fr?d?ric}{Fr\'ed\'eric} Bertrand, Myriam Maumy-Bertrand, Marie-Pierre Gaub, Nicolas Meyer, \enc{Chimiom?trie}{Chimiom\'etrie} 2010, Paris, 2010.}
\author{\enc{Fr?d?ric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{\code{\link[survival]{coxph}}}
\examples{
data(micro.censure)
Y_train_micro <- micro.censure$survyear[1:80]
C_train_micro <- micro.censure$DC[1:80]

DR_coxph(Y_train_micro,C_train_micro,plot=TRUE)
DR_coxph(Y_train_micro,C_train_micro,scaleY=FALSE,plot=TRUE)
DR_coxph(Y_train_micro,C_train_micro,scaleY=TRUE,plot=TRUE)

rm(Y_train_micro,C_train_micro)
}
\keyword{models}
\keyword{regression}
