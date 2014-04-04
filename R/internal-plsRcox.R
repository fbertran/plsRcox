correctp.cox=function (x, y, eta, K, kappa, select, fit) 
{
  force(K)
  if (min(eta) < 0 | max(eta) >= 1) {
    if (max(eta) == 1) {
      stop("eta should be strictly less than 1!")
    }
    if (length(eta) == 1) {
      stop("eta should be between 0 and 1!")
    }
    else {
      stop("eta should be between 0 and 1! \n  Choose appropriate range of eta!")
    }
  }
  if (max(K) > ncol(x)) {
    stop("K cannot exceed the number of predictors! Pick up smaller K!")
  }
  if (max(K) >= nrow(x)) {
    stop("K cannot exceed the sample size! Pick up smaller K!")
  }
  if (min(K) <= 0 | !all(K%%1 == 0)) {
    if (length(K) == 1) {
      stop("K should be a positive integer!")
    }
    else {
      stop("K should be a positive integer! \n  Choose appropriate range of K!")
    }
  }
  if (kappa > 0.5 | kappa < 0) {
    cat("kappa should be between 0 and 0.5! kappa=0.5 is used. \n\n")
    kappa <- 0.5
  }
  if (select != "pls2" & select != "simpls") {
    cat("Invalid PLS algorithm for variable selection.\n")
    cat("pls2 algorithm is used. \n\n")
    select <- "pls2"
  }
  fits <- c("regression", "canonical", "invariant", "classic")
  if (!any(fit == fits)) {
    cat("Invalid PLS algorithm for model fitting\n")
    cat("regression algorithm is used. \n\n")
    fit <- "regression"
  }
  list(K = K, eta = eta, kappa = kappa, select = select, fit = fit)
}

spls.cox=function (x, y, K, eta, kappa = 0.5, select = "pls2", fit = "regression", 
                   scale.x = TRUE, scale.y = FALSE, eps = 1e-04, maxstep = 100, 
                   trace = FALSE) 
{
  force(K)
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  ip <- c(1:p)
  y <- as.matrix(y)
  q <- ncol(y)
  one <- matrix(1, 1, n)
  mu <- one %*% y/n
  y <- scale(y, drop(mu), FALSE)
  meanx <- drop(one %*% x)/n
  x <- scale(x, meanx, FALSE)
  if (scale.x) {
    normx <- sqrt(drop(one %*% (x^2))/(n - 1))
    if (any(normx < .Machine$double.eps)) {
      stop("Some of the columns of the predictor matrix have zero variance.")
    }
    x <- scale(x, FALSE, normx)
  }
  else {
    normx <- rep(1, p)
  }
  if (scale.y) {
    normy <- sqrt(drop(one %*% (y^2))/(n - 1))
    if (any(normy < .Machine$double.eps)) {
      stop("Some of the columns of the response matrix have zero variance.")
    }
    y <- scale(y, FALSE, normy)
  }
  else {
    normy <- rep(1, q)
  }
  betahat <- matrix(0, p, q)
  betamat <- list()
  x1 <- x
  y1 <- y
  type <- correctp.cox(x, y, eta, K, kappa, select, fit)
  eta <- type$eta
  K <- type$K
  kappa <- type$kappa
  select <- type$select
  fit <- type$fit
  if (is.null(colnames(x))) {
    xnames <- c(1:p)
  }
  else {
    xnames <- colnames(x)
  }
  new2As <- list()
  if (trace) {
    cat("The variables that join the set of selected variables at each step:\n")
  }
  for (k in 1:K) {
    Z <- t(x1) %*% y1
    what <- spls.dv(Z, eta, kappa, eps, maxstep)
    A <- unique(ip[what != 0 | betahat[, 1] != 0])
    new2A <- ip[what != 0 & betahat[, 1] == 0]
    xA <- x[, A, drop = FALSE]
    plsfit <- pls.cox(X=xA, Y=y, ncomp = min(k, length(A)), mode = fit, 
                      scale.X = FALSE, scale.Y=FALSE)
    predplsfit <- predict.pls.cox(plsfit,newdata=xA,scale.X = FALSE, scale.Y=FALSE)
    betahat <- matrix(0, p, q)
    betahat[A, ] <- matrix(predplsfit$B.hat[,,plsfit$ncomp], length(A), q)
    betamat[[k]] <- betahat
    #        pj <- plsfit$projection
    if (select == "pls2") {
      y1 <- y - predplsfit$predict[,,plsfit$ncomp]
    }
    #        if (select == "simpls") {
    #            pw <- pj %*% solve(t(pj) %*% pj) %*% t(pj)
    #            x1 <- x
    #            x1[, A] <- x[, A, drop = FALSE] - x[, A, drop = FALSE] %*% 
    #                pw
    #        }
    new2As[[k]] <- new2A
    if (trace) {
      if (length(new2A) <= 10) {
        cat(paste("- ", k, "th step (K=", k, "):\n", 
                  sep = ""))
        cat(xnames[new2A])
        cat("\n")
      }
      else {
        cat(paste("- ", k, "th step (K=", k, "):\n", 
                  sep = ""))
        nlines <- ceiling(length(new2A)/10)
        for (i in 0:(nlines - 2)) {
          cat(xnames[new2A[(10 * i + 1):(10 * (i + 1))]])
          cat("\n")
        }
        cat(xnames[new2A[(10 * (nlines - 1) + 1):length(new2A)]])
        cat("\n")
      }
    }
  }
  if (!is.null(colnames(x))) {
    rownames(betahat) <- colnames(x)
  }
  if (q > 1 & !is.null(colnames(y))) {
    colnames(betahat) <- colnames(y)
  }
  object <- list(x = x, y = y, betahat = betahat, A = A, betamat = betamat, 
                 new2As = new2As, mu = mu, meanx = meanx, normx = normx, 
                 normy = normy, eta = eta, K = K, kappa = kappa, select = select, 
                 fit = fit, projection = NA, plsmod=plsfit)
  class(object) <- "spls"
  object
}

ust=function (b, eta) 
{
  b.ust <- matrix(0, length(b), 1)
  if (eta < 1) {
    valb <- abs(b) - eta * max(abs(b))
    b.ust[valb >= 0] <- valb[valb >= 0] * (sign(b))[valb >= 
                                                      0]
  }
  return(b.ust)
}

spls.dv <- function (Z, eta, kappa, eps, maxstep) 
{
  p <- nrow(Z)
  q <- ncol(Z)
  Znorm1 <- median(abs(Z))
  Z <- Z/Znorm1
  if (q == 1) {
    c <- ust(Z, eta)
  }
  if (q > 1) {
    M <- Z %*% t(Z)
    dis <- 10
    i <- 1
    if (kappa == 0.5) {
      c <- matrix(10, p, 1)
      c.old <- c
      while (dis > eps & i <= maxstep) {
        mcsvd <- svd(M %*% c)
        a <- mcsvd$u %*% t(mcsvd$v)
        c <- ust(M %*% a, eta)
        dis <- max(abs(c - c.old))
        c.old <- c
        i <- i + 1
      }
    }
    if (kappa > 0 & kappa < 0.5) {
      kappa2 <- (1 - kappa)/(1 - 2 * kappa)
      c <- matrix(10, p, 1)
      c.old <- c
      h <- function(lambda) {
        alpha <- solve(M + lambda * diag(p)) %*% M %*% 
          c
        obj <- t(alpha) %*% alpha - 1/kappa2^2
        return(obj)
      }
      if (h(eps) * h(1e+30) > 0) {
        while (h(eps) <= 1e+05) {
          M <- 2 * M
          c <- 2 * c
        }
      }
      while (dis > eps & i <= maxstep) {
        if (h(eps) * h(1e+30) > 0) {
          while (h(eps) <= 1e+05) {
            M <- 2 * M
            c <- 2 * c
          }
        }
        lambdas <- uniroot(h, c(eps, 1e+30))$root
        a <- kappa2 * solve(M + lambdas * diag(p)) %*% 
          M %*% c
        c <- ust(M %*% a, eta)
        dis <- max(abs(c - c.old))
        c.old <- c
        i <- i + 1
      }
    }
  }
  return(c)
}

pls.cox=function(X, Y, ncomp = 2, mode = c("regression", "canonical", "invariant", "classic"), max.iter = 500, tol = 1e-06, scale.X=TRUE, scale.Y=TRUE, ...) 
{
  force(ncomp)
  if (length(dim(X)) != 2) 
    stop("'X' must be a numeric matrix.")
  X = as.matrix(X)
  Y = as.matrix(Y)
  if (!is.numeric(X) || !is.numeric(Y)) 
    stop("'X' and/or 'Y' must be a numeric matrix.")
  n = nrow(X)
  q = ncol(Y)
  if ((n != nrow(Y))) 
    stop("unequal number of rows in 'X' and 'Y'.")
  if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0) 
    stop("invalid number of variates, 'ncomp'.")
  nzv = mixOmics::nearZeroVar(X, ...)
  if (length(nzv$Position > 0)) {
    warning("Zero- or near-zero variance predictors. \n  Reset predictors matrix to not near-zero variance predictors.\n  See $nzv for problematic predictors.")
    X = X[, -nzv$Position]
  }
  p = ncol(X)
  ncomp = round(ncomp)
  if (ncomp > p) {
    warning("Reset maximum number of variates 'ncomp' to ncol(X) = ", 
            p, ".")
    ncomp = p
  }
  mode = match.arg(mode)
  X.names = dimnames(X)[[2]]
  if (is.null(X.names)) 
    X.names = paste("X", 1:p, sep = "")
  if (dim(Y)[2] == 1) 
    Y.names = "Y"
  else {
    Y.names = dimnames(Y)[[2]]
    if (is.null(Y.names)) 
      Y.names = paste("Y", 1:q, sep = "")
  }
  ind.names = dimnames(X)[[1]]
  if (is.null(ind.names)) {
    ind.names = dimnames(Y)[[1]]
    rownames(X) = ind.names
  }
  if (is.null(ind.names)) {
    ind.names = 1:n
    rownames(X) = rownames(Y) = ind.names
  }
  if(scale.X){X = scale(X, center = TRUE, scale = TRUE)}
  if(scale.Y){Y = scale(Y, center = TRUE, scale = TRUE)}
  X.temp = X
  Y.temp = Y
  mat.t = matrix(nrow = n, ncol = ncomp)
  mat.u = matrix(nrow = n, ncol = ncomp)
  mat.a = matrix(nrow = p, ncol = ncomp)
  mat.b = matrix(nrow = q, ncol = ncomp)
  mat.c = matrix(nrow = p, ncol = ncomp)
  mat.d = matrix(nrow = q, ncol = ncomp)
  mat.e = matrix(nrow = q, ncol = ncomp)
  n.ones = rep(1, n)
  p.ones = rep(1, p)
  q.ones = rep(1, q)
  na.X = FALSE
  na.Y = FALSE
  is.na.X = is.na(X)
  is.na.Y = is.na(Y)
  if (any(is.na.X)) 
    na.X = TRUE
  if (any(is.na.Y)) 
    na.Y = TRUE
  for (h in 1:ncomp) {
    u = Y.temp[, 1]
    if (any(is.na(u))) 
      u[is.na(u)] = 0
    a.old = 0
    b.old = 0
    iter = 1
    if (na.X) {
      X.aux = X.temp
      X.aux[is.na.X] = 0
    }
    if (na.Y) {
      Y.aux = Y.temp
      Y.aux[is.na.Y] = 0
    }
    repeat {
      if (na.X) {
        a = crossprod(X.aux, u)
        U = drop(u) %o% p.ones
        U[is.na.X] = 0
        u.norm = crossprod(U)
        a = a/diag(u.norm)
        a = a/drop(sqrt(crossprod(a)))
        t = X.aux %*% a
        A = drop(a) %o% n.ones
        A[t(is.na.X)] = 0
        a.norm = crossprod(A)
        t = t/diag(a.norm)
      }
      else {
        a = crossprod(X.temp, u)/drop(crossprod(u))
        a = a/drop(sqrt(crossprod(a)))
        t = X.temp %*% a/drop(crossprod(a))
      }
      if (na.Y) {
        b = crossprod(Y.aux, t)
        T = drop(t) %o% q.ones
        T[is.na.Y] = 0
        t.norm = crossprod(T)
        b = b/diag(t.norm)
        u = Y.aux %*% b
        B = drop(b) %o% n.ones
        B[t(is.na.Y)] = 0
        b.norm = crossprod(B)
        u = u/diag(b.norm)
      }
      else {
        b = crossprod(Y.temp, t)/drop(crossprod(t))
        u = Y.temp %*% b/drop(crossprod(b))
      }
      if (crossprod(a - a.old) < tol) 
        break
      if (iter == max.iter) {
        warning(paste("Maximum number of iterations reached for dimension", 
                      h), call. = FALSE)
        break
      }
      a.old = a
      b.old = b
      iter = iter + 1
    }
    if (na.X) {
      X.aux = X.temp
      X.aux[is.na.X] = 0
      c = crossprod(X.aux, t)
      T = drop(t) %o% p.ones
      T[is.na.X] = 0
      t.norm = crossprod(T)
      c = c/diag(t.norm)
    }
    else {
      c = crossprod(X.temp, t)/drop(crossprod(t))
    }
    X.temp = X.temp - t %*% t(c)
    if (mode == "canonical") {
      if (na.Y) {
        Y.aux = Y.temp
        Y.aux[is.na.Y] = 0
        e = crossprod(Y.aux, u)
        U = drop(u) %o% q.ones
        U[is.na.Y] = 0
        u.norm = crossprod(U)
        e = e/diag(u.norm)
      }
      else {
        e = crossprod(Y.temp, u)/drop(crossprod(u))
      }
      Y.temp = Y.temp - u %*% t(e)
    }
    if (mode == "classic") 
      Y.temp = Y.temp - t %*% t(b)
    if (mode == "regression") {
      if (na.Y) {
        Y.aux = Y.temp
        Y.aux[is.na.Y] = 0
        d = crossprod(Y.aux, t)
        T = drop(t) %o% q.ones
        T[is.na.Y] = 0
        t.norm = crossprod(T)
        d = d/diag(t.norm)
      }
      else {
        d = crossprod(Y.temp, t)/drop(crossprod(t))
      }
      Y.temp = Y.temp - t %*% t(d)
    }
    if (mode == "invariant") 
      Y.temp = Y
    mat.t[, h] = t
    mat.u[, h] = u
    mat.a[, h] = a
    mat.b[, h] = b
    mat.c[, h] = c
    if (mode == "regression") 
      mat.d[, h] = d
    if (mode == "canonical") 
      mat.e[, h] = e
  }
  rownames(mat.a) = rownames(mat.c) = X.names
  rownames(mat.b) = Y.names
  rownames(mat.t) = rownames(mat.u) = ind.names
  comp = paste("comp", 1:ncomp)
  colnames(mat.t) = colnames(mat.u) = comp
  colnames(mat.a) = colnames(mat.b) = colnames(mat.c) = comp
  cl = match.call()
  cl[[1]] = as.name("pls")
  result = list(call = cl, X = X, Y = Y, ncomp = ncomp, mode = mode, 
                mat.c = mat.c, mat.d = mat.d, mat.e = mat.e, variates = list(X = mat.t, 
                                                                             Y = mat.u), loadings = list(X = mat.a, Y = mat.b), 
                names = list(X = X.names, Y = Y.names, indiv = ind.names))
  if (length(nzv$Position > 0)) 
    result$nzv = nzv
  class(result) = "pls"
  return(invisible(result))
}


predict.pls.cox=function(object, newdata, scale.X=TRUE, scale.Y=TRUE,...) 
{
  if (missing(newdata)) 
    stop("No new data available.")
  X = object$X
  Y = object$Y
  q = ncol(Y)
  p = ncol(X)
  if (length(dim(newdata)) == 2) {
    if (ncol(newdata) != p) 
      stop("'newdata' must be a numeric matrix with ncol = ", 
           p, " or a vector of length = ", p, ".")
  }
  if (length(dim(newdata)) == 0) {
    if (length(newdata) != p) 
      stop("'newdata' must be a numeric matrix with ncol = ", 
           p, " or a vector of length = ", p, ".")
    dim(newdata) = c(1, p)
  }
  ncomp = object$ncomp
  a = object$loadings$X
  b = object$loadings$Y
  c = object$mat.c
  if(scale.X){means.X = attr(X, "scaled:center")}
  if(scale.Y){means.Y = attr(Y, "scaled:center")}
  if(scale.X){sigma.X = attr(X, "scaled:scale")}
  if(scale.Y){sigma.Y = attr(Y, "scaled:scale")}
  newdata = as.matrix(newdata)
  ones = matrix(rep(1, nrow(newdata)), ncol = 1)
  B.hat = array(0, dim = c(p, q, ncomp))
  Y.hat = array(0, dim = c(nrow(newdata), q, ncomp))
  t.pred = array(0, dim = c(nrow(newdata), ncomp))
  for (h in 1:ncomp) {
    W = a[, 1:h] %*% solve(t(c[, 1:h]) %*% a[, 1:h])
    B = W %*% drop(t(b[, 1:h]))
    if(scale.Y){B = scale(B, center = FALSE, scale = 1/sigma.Y)}
    if(scale.X){B = as.matrix(scale(t(B), center = FALSE, scale = sigma.X))}
    if(!scale.X){B = as.matrix(t(B))}
    if(scale.X&scale.Y){intercept = -scale(B, center = FALSE, scale = 1/means.X)
                        intercept = matrix(apply(intercept, 1, sum) + means.Y, 
                                           nrow = 1)
                        Y.hat[, , h] = newdata %*% t(B) + ones %*% intercept}
    if(scale.X&!scale.Y){intercept = -scale(B, center = FALSE, scale = 1/means.X)
                         intercept = matrix(apply(intercept, 1, sum), nrow = 1)
                         Y.hat[, , h] = newdata %*% t(B) + ones %*% intercept}
    if(!scale.X&scale.Y){intercept = -B
                         intercept = matrix(apply(intercept, 1, sum) + means.Y, 
                                            nrow = 1)
                         Y.hat[, , h] = newdata %*% t(B) + ones %*% intercept}
    if(!scale.X&!scale.Y){Y.hat[, , h] = newdata %*% t(B)}
    if(!scale.X){t.pred[, h] = newdata %*% W[, h]}
    if(scale.X){t.pred[, h] = scale(newdata, center = means.X, scale = sigma.X) %*% W[, h]}
    B.hat[, , h] = B
  }
  rownames(t.pred) = rownames(newdata)
  colnames(t.pred) = paste("dim", c(1:ncomp), sep = " ")
  rownames(Y.hat) = rownames(newdata)
  colnames(Y.hat) = colnames(Y)
  return(invisible(list(predict = Y.hat, variates = t.pred, 
                        B.hat = B.hat)))
}

