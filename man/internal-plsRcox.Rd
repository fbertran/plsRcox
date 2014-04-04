\name{internal-plsRcox}
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
correctp.cox(x, y, eta, K, kappa, select, fit) 
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
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
