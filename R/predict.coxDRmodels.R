#' Predict method for DR-based Cox component models
#'
#' This function provides prediction methods for the rich objects returned by
#' \code{coxplsDR(..., allres = TRUE)}, \code{coxsplsDR(..., allres = TRUE)}
#' and \code{coxDKsplsDR(..., allres = TRUE)}.
#'
#' @param object An object returned by one of the DR-based fitting functions
#' with \code{allres = TRUE}.
#' @param newdata An optional data frame or matrix containing original
#' covariates. If omitted, predictions are computed on the training data.
#' @param comps Number of latent components to use for prediction.
#' @param type Type of predicted value. Choices are the linear predictor
#' (\code{"lp"}), the risk score exp(lp) (\code{"risk"}), the expected number
#' of events (\code{"expected"}), the terms of the linear predictor
#' (\code{"terms"}) or the latent component scores (\code{"scores"}).
#' @param se.fit If \code{TRUE}, pointwise standard errors are produced by the
#' underlying Cox model.
#' @param reference Reference level used to center relative predictions. This
#' is passed to \code{\link[survival]{predict.coxph}} and affects
#' \code{type = "lp"}, \code{"risk"} and \code{"terms"}.
#' @param y Optional \code{\link[survival]{Surv}} response to use with
#' \code{type = "expected"} when predicting on genuinely new observations.
#' @param weights Optional case weights used when rebuilding a model matrix for
#' formula-based fits.
#' @param verbose Should some details be displayed?
#' @param \dots Additional arguments passed to \code{\link[survival]{predict.coxph}}.
#' @return A vector, matrix or list of predictions depending on \code{type} and
#' \code{se.fit}.
#' @seealso \code{\link[survival]{predict.coxph}}, \code{\link{coxplsDR}},
#' \code{\link{coxsplsDR}}, \code{\link{coxDKsplsDR}}
#' @name predict.coxDRmodel
NULL

.plsRcox_prepare_component_newdata <- function(object, newdata, weights = NULL) {
  nrnd <- nrow(newdata)
  if (!is.null(object$XplanFormula)) {
    mf0 <- match.call(expand.dots = FALSE)
    m0 <- match(c("weights"), names(mf0), 0L)
    mf0 <- mf0[c(1L, m0)]
    mf0$data <- newdata
    mf0$formula <- object$XplanFormula
    mf0$drop.unused.levels <- TRUE
    mf0[[1L]] <- as.name("model.frame")
    mf0 <- eval(mf0, parent.frame())
    mt0 <- attr(mf0, "terms")
    attr(mt0, "intercept") <- 0L
    newdata.frame <- if (!is.empty.model(mt0)) {
      model.matrix(mt0, mf0, contrasts)[, -1, drop = FALSE]
    } else {
      matrix(, nrnd, 0L)
    }
    weights <- as.vector(model.weights(mf0))
  } else {
    newdata.frame <- newdata
  }

  list(newdata.frame = newdata.frame, weights = weights)
}

.plsRcox_scale_expected_response_shared <- function(response, rep_y) {
  center <- attr(rep_y, "scaled:center")
  scale <- attr(rep_y, "scaled:scale")
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
  response[, time_columns] <- sweep(
    sweep(response[, time_columns, drop = FALSE], 2, center),
    2,
    scale,
    "/"
  )
  response
}

.plsRcox_match_training_rows_shared <- function(candidate, training_data) {
  candidate <- as.matrix(candidate)
  training_data <- as.matrix(training_data)

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

.plsRcox_resolve_expected_response_shared <- function(object, y, newdata_rows = NULL, score_rows = NULL, newdata_matrix = NULL) {
  if (!is.null(y)) {
    if (!inherits(y, "Surv")) {
      stop("'y' must inherit from 'Surv' when type = 'expected'.")
    }
    if (!is.null(score_rows) && NROW(y) != score_rows) {
      stop("'y' must have the same number of rows as the predicted newdata.")
    }
    return(.plsRcox_scale_expected_response_shared(y, object$RepY))
  }

  if (is.null(newdata_rows) && is.null(newdata_matrix)) {
    return(object$RepY)
  }

  train_rows <- rownames(object$dataX)
  if (!is.null(train_rows) && !is.null(newdata_rows)) {
    row_index <- match(newdata_rows, train_rows)
    if (!anyNA(row_index)) {
      return(object$RepY[row_index])
    }
  }

  if (!is.null(newdata_matrix)) {
    row_index <- .plsRcox_match_training_rows_shared(newdata_matrix, object$dataX)
    if (!is.null(row_index)) {
      return(object$RepY[row_index])
    }
  }

  stop("type = 'expected' requires follow-up information for newdata. Supply 'y = survival::Surv(...)' or keep row names that match the training rows.")
}

.plsRcox_build_component_prediction_frame <- function(score_names, response_name, scores, response = NULL) {
  scores <- as.matrix(scores)
  padding <- length(score_names) - ncol(scores)
  if (padding > 0) {
    scores <- cbind(scores, matrix(0, nrow = nrow(scores), ncol = padding))
  }
  if (ncol(scores) > length(score_names)) {
    scores <- scores[, seq_len(length(score_names)), drop = FALSE]
  }

  score_df <- as.data.frame(scores)
  if (length(score_names) > 0) {
    names(score_df) <- score_names
  }

  if (is.null(response)) {
    return(score_df)
  }

  score_df[[response_name]] <- response
  score_df[c(response_name, score_names)]
}

.plsRcox_scores_coxplsDR <- function(object, newdata.frame) {
  if (ncol(object$tt_plsDR) == 0L) {
    return(matrix(, nrow = nrow(newdata.frame), ncol = 0L))
  }
  predict.pls.cox(
    object$plsDR_mod,
    newdata = scale(newdata.frame, object$XplanCent, object$XplanScal),
    scale.X = FALSE,
    scale.Y = FALSE
  )$variates
}

.plsRcox_scores_coxsplsDR <- function(object, newdata.frame) {
  if (ncol(object$tt_splsDR) == 0L) {
    return(matrix(, nrow = nrow(newdata.frame), ncol = 0L))
  }
  avalues <- object$splsDR_mod$A
  predict.pls.cox(
    object$splsDR_modplsr,
    newdata = scale(newdata.frame[, avalues, drop = FALSE], object$XplanCent[avalues], object$XplanScal[avalues]),
    scale.X = FALSE,
    scale.Y = FALSE
  )$variates
}

.plsRcox_scores_coxDKsplsDR <- function(object, newdata.frame) {
  if (ncol(object$tt_DKsplsDR) == 0L) {
    return(matrix(, nrow = nrow(newdata.frame), ncol = 0L))
  }
  train_scaled <- scale(object$dataX, object$XplanCent, object$XplanScal)
  new_scaled <- scale(newdata.frame, object$XplanCent, object$XplanScal)
  kernel_new <- kernelMatrix(object$kernDKsplsDR_mod, as.matrix(new_scaled), as.matrix(train_scaled))
  avalues <- object$DKsplsDR_mod$A
  predict.pls.cox(
    object$DKsplsDR_modplsr,
    newdata = as.data.frame(kernel_new)[, avalues, drop = FALSE],
    scale.X = FALSE,
    scale.Y = FALSE
  )$variates
}

.plsRcox_predict_component_model <- function(object, newdata, comps, type, se.fit, reference, y, weights, verbose, score_slot, cox_slot, score_fun, ...) {
  if (missing(type)) {
    type <- "lp"
  }
  if (!(type %in% c("lp", "risk", "expected", "terms", "scores"))) {
    stop("Invalid type specification")
  }
  type <- match.arg(type, c("lp", "risk", "expected", "terms", "scores"))
  reference <- match.arg(reference, c("strata", "sample", "zero"))

  score_template <- object[[score_slot]]
  cox_model <- object[[cox_slot]]
  max_comps <- ncol(score_template)

  if (missing(comps)) {
    comps <- max_comps
  }
  if (comps > max_comps) {
    stop("Cannot predict using more components than extracted.")
  }
  if (comps < 0) {
    stop("'comps' must be non-negative.")
  }

  score_names <- names(score_template)
  response_name <- all.vars(stats::formula(cox_model))[1]

  if (missing(newdata) || is.null(newdata)) {
    scores <- as.matrix(score_template)
    if (comps < max_comps) {
      scores <- scores[, seq_len(comps), drop = FALSE]
    }
    if (type == "scores") {
      return(scores)
    }
    if (type == "expected") {
      ttpred <- .plsRcox_build_component_prediction_frame(score_names, response_name, scores, .plsRcox_resolve_expected_response_shared(object, y))
    } else {
      ttpred <- .plsRcox_build_component_prediction_frame(score_names, response_name, scores)
    }
    return(predict(cox_model, newdata = ttpred, type = type, se.fit = se.fit, reference = reference, ...))
  }

  prepared <- .plsRcox_prepare_component_newdata(object, newdata, weights)
  newdata.frame <- prepared$newdata.frame
  row_index <- NULL
  train_rows <- rownames(object$dataX)
  if (!is.null(train_rows) && !is.null(rownames(newdata.frame))) {
    matched_rows <- match(rownames(newdata.frame), train_rows)
    if (!anyNA(matched_rows)) {
      row_index <- matched_rows
    }
  }
  if (is.null(row_index)) {
    row_index <- .plsRcox_match_training_rows_shared(newdata.frame, object$dataX)
  }

  if (!is.null(row_index)) {
    scores <- as.matrix(score_template[row_index, , drop = FALSE])
  } else {
    scores <- as.matrix(score_fun(object, newdata.frame))
  }
  if (comps < max_comps) {
    scores <- scores[, seq_len(comps), drop = FALSE]
  }

  if (type == "scores") {
    colnames(scores) <- score_names[seq_len(comps)]
    return(scores)
  }

  if (type == "expected") {
    ttpred <- .plsRcox_build_component_prediction_frame(
      score_names,
      response_name,
      scores,
      .plsRcox_resolve_expected_response_shared(
        object,
        y,
        newdata_rows = rownames(newdata.frame),
        score_rows = nrow(newdata.frame),
        newdata_matrix = newdata.frame
      )
    )
  } else {
    ttpred <- .plsRcox_build_component_prediction_frame(score_names, response_name, scores)
  }

  predict(cox_model, newdata = ttpred, type = type, se.fit = se.fit, reference = reference, ...)
}

#' @rdname predict.coxDRmodel
#' @method predict coxplsDRmodel
#' @export
predict.coxplsDRmodel <- function(object, newdata, comps = ncol(object$tt_plsDR), type = c("lp", "risk", "expected", "terms", "scores"), se.fit = FALSE, reference = c("strata", "sample", "zero"), y = NULL, weights = NULL, verbose = TRUE, ...) {
  .plsRcox_predict_component_model(
    object = object,
    newdata = newdata,
    comps = comps,
    type = type,
    se.fit = se.fit,
    reference = reference,
    y = y,
    weights = weights,
    verbose = verbose,
    score_slot = "tt_plsDR",
    cox_slot = "cox_plsDR",
    score_fun = .plsRcox_scores_coxplsDR,
    ...
  )
}

#' @rdname predict.coxDRmodel
#' @method predict coxsplsDRmodel
#' @export
predict.coxsplsDRmodel <- function(object, newdata, comps = ncol(object$tt_splsDR), type = c("lp", "risk", "expected", "terms", "scores"), se.fit = FALSE, reference = c("strata", "sample", "zero"), y = NULL, weights = NULL, verbose = TRUE, ...) {
  .plsRcox_predict_component_model(
    object = object,
    newdata = newdata,
    comps = comps,
    type = type,
    se.fit = se.fit,
    reference = reference,
    y = y,
    weights = weights,
    verbose = verbose,
    score_slot = "tt_splsDR",
    cox_slot = "cox_splsDR",
    score_fun = .plsRcox_scores_coxsplsDR,
    ...
  )
}

#' @rdname predict.coxDRmodel
#' @method predict coxDKsplsDRmodel
#' @export
predict.coxDKsplsDRmodel <- function(object, newdata, comps = ncol(object$tt_DKsplsDR), type = c("lp", "risk", "expected", "terms", "scores"), se.fit = FALSE, reference = c("strata", "sample", "zero"), y = NULL, weights = NULL, verbose = TRUE, ...) {
  .plsRcox_predict_component_model(
    object = object,
    newdata = newdata,
    comps = comps,
    type = type,
    se.fit = se.fit,
    reference = reference,
    y = y,
    weights = weights,
    verbose = verbose,
    score_slot = "tt_DKsplsDR",
    cox_slot = "cox_DKsplsDR",
    score_fun = .plsRcox_scores_coxDKsplsDR,
    ...
  )
}
