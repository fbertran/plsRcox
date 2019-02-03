\name{internal-plsRcox}
\alias{getIndic}
\alias{getIndicCV}
\alias{getIndicCViAUCSH}
\alias{getIndicCViAUCSurvROCTest}
\alias{correctp.cox}
\alias{spls.cox}
\alias{ust}
\alias{spls.dv}
\alias{pls.cox}
\alias{predict.pls.cox}
\title{
Internal plsRcox functions
}
\description{
Internal plsRcox functions.
}
\usage{
getIndic(lp,lpnew,Surv.rsp,Surv.rsp.new,times.auc=seq(10,1000,10),
times.prederr=1:500,train.fit,train.fit.cph,tmax.train=365,
tmax.test=365,TR,TE,plot.it=TRUE)
getIndicCV(lp,lpnew,Surv.rsp,Surv.rsp.new,times.auc=seq(10,1000,10),
times.prederr=1:500,train.fit,plot.it=FALSE,
tmax.train=max(Surv.rsp[,"time"][ object$Surv.rsp[,"status"] == 1 ]),
tmax.test=max(Surv.rsp.new[,"time"][ object$Surv.rsp.new[,"status"] == 1 ]))
getIndicCViAUCSH(lp,lpnew,Surv.rsp,Surv.rsp.new,times.auc=seq(10,1000,10),
times.prederr=1:500,train.fit,plot.it=FALSE,
tmax.train=max(Surv.rsp[,"time"][ object$Surv.rsp[,"status"] == 1 ]),
tmax.test=max(Surv.rsp.new[,"time"][ object$Surv.rsp.new[,"status"] == 1 ]))
getIndicCViAUCSurvROCTest(lp,lpnew,Surv.rsp,Surv.rsp.new,times.auc=seq(10,1000,10),
times.prederr=1:500,train.fit,plot.it=FALSE,
tmax.train=max(Surv.rsp[,"time"][ object$Surv.rsp[,"status"] == 1 ]),
tmax.test=max(Surv.rsp.new[,"time"][ object$Surv.rsp.new[,"status"] == 1 ]))
correctp.cox(x, y, eta, K, kappa, select, fit, verbose=FALSE) 
spls.cox(x, y, K, eta, kappa = 0.5, select = "pls2", 
fit = "regression", scale.x = TRUE, scale.y = FALSE, eps = 1e-04, 
maxstep = 100, trace = FALSE) 
ust(b, eta) 
spls.dv(Z, eta, kappa, eps, maxstep) 
pls.cox(X, Y, ncomp = 2, mode = c("regression", "canonical", 
"invariant", "classic"), max.iter = 500, tol = 1e-06, scale.X=TRUE, scale.Y=TRUE, ...) 
predict.pls.cox(object, newdata, scale.X=TRUE, scale.Y=TRUE,...) 

}
%- maybe also 'usage' for other objects documented here.
\details{
These are not to be called by the user.
}
\references{
plsRcox, Cox-Models in a high dimensional setting in R, Frederic Bertrand, Philippe Bastien, Nicolas Meyer and Myriam Maumy-Bertrand (2014). Proceedings of User2014!, Los Angeles, page 152.\cr
  
Deviance residuals-based sparse PLS and sparse kernel PLS regression for censored data, Philippe Bastien, Frederic Bertrand, Nicolas Meyer and Myriam Maumy-Bertrand (2015), Bioinformatics, 31(3):397-404, doi:10.1093/bioinformatics/btu660.
}
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
