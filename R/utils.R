# Internal helpers - not exported.
#
# `balanced.folds.vec()` is the workhorse used throughout CASPEN to
# generate stratified cross-validation folds for binary, multi-class,
# and survival outcomes. It returns an integer vector of fold
# assignments of length `length(y)`.

#' @keywords internal
#' @noRd
balanced.folds <- function(y, nfolds = min(min(table(y)), 10)) {
  y[is.na(y)] <- resample(y[!is.na(y)],
                          size = sum(is.na(y)),
                          replace = TRUE)
  totals <- table(y)

  if (length(totals) < 2) {
    return(cv.folds(length(y), nfolds))
  } else {
    fmax   <- max(totals)
    nfolds <- min(nfolds, fmax)
    nfolds <- max(nfolds, 2)
    yids   <- split(seq(y), y)
    bigmat <- matrix(NA, ceiling(fmax / nfolds) * nfolds, length(totals))

    for (i in seq(totals)) {
      if (length(yids[[i]]) > 1) {
        bigmat[seq(totals[i]), i] <- sample(yids[[i]])
      }
      if (length(yids[[i]]) == 1) {
        bigmat[seq(totals[i]), i] <- yids[[i]]
      }
    }

    smallmat <- matrix(bigmat, nrow = nfolds)
    smallmat <- permute.rows(t(smallmat))
    res <- vector("list", nfolds)

    for (j in seq_len(nfolds)) {
      jj <- !is.na(smallmat[, j])
      res[[j]] <- smallmat[jj, j]
    }
    return(res)
  }
}

#' @keywords internal
#' @noRd
cv.folds <- function(n, nfolds) {
  ids <- sample(seq_len(n))
  split(ids, rep(seq_len(nfolds), length.out = n))
}

#' @keywords internal
#' @noRd
resample <- function(x, size, ...) {
  if (length(x) <= 1) {
    if (!missing(size) && size == 0) x[FALSE] else x
  } else {
    sample(x, size, ...)
  }
}

#' @keywords internal
#' @noRd
permute.rows <- function(A) {
  B <- t(A)
  row.list <- split(B, rep(1:ncol(B), each = nrow(B)))
  t(sapply(row.list, permute.vector))
}

#' @keywords internal
#' @noRd
permute.vector <- function(x) {
  x[sample(1:length(x))]
}

#' @keywords internal
#' @noRd
svm_binary_probability <- function(fit, newdata, positive = "1") {
  pred <- predict(fit, newdata = as.data.frame(newdata), probability = TRUE)
  probs <- attr(pred, "probabilities")
  if (is.null(probs)) {
    stop("SVM did not return class probabilities. Fit the model with probability = TRUE.")
  }
  probs <- as.matrix(probs)
  class_names <- colnames(probs)
  positive <- as.character(positive)
  pos_col <- match(positive, class_names)
  if (is.na(pos_col)) {
    pos_col <- match(paste0("Class", positive), class_names)
  }
  if (is.na(pos_col)) {
    pos_col <- match(make.names(positive), make.names(class_names))
  }
  if (is.na(pos_col)) {
    pos_col <- if (ncol(probs) >= 2) 2 else 1
  }
  as.numeric(probs[, pos_col])
}

#' @keywords internal
#' @noRd
balanced.folds.vec <- function(y, k.folds) {
  bfold <- balanced.folds(y = y, nfolds = k.folds)
  bfold.vec <- unlist(bfold)
  folds <- unlist(lapply(seq_len(k.folds),
                         function(x) rep(x, length(bfold[[x]]))))
  folds[order(bfold.vec)]
}

#' @keywords internal
#' @noRd
survival_cindex <- function(time, status, marker) {
  keep <- stats::complete.cases(time, status, marker) &
    is.finite(time) & is.finite(status) & is.finite(marker)

  if (sum(keep) < 2 || length(unique(status[keep])) < 2 ||
      length(unique(marker[keep])) < 2) {
    return(c(cindex = NA_real_, ci.low = NA_real_, ci.up = NA_real_))
  }

  fit <- tryCatch(
    survival::concordance(
      survival::Surv(time[keep], status[keep]) ~ I(-marker[keep])
    ),
    error = function(e) NULL
  )

  if (is.null(fit) || is.null(fit$concordance)) {
    return(c(cindex = NA_real_, ci.low = NA_real_, ci.up = NA_real_))
  }

  cindex <- as.numeric(fit$concordance)
  se <- if (!is.null(fit$std.err)) {
    as.numeric(fit$std.err)
  } else if (!is.null(fit$var)) {
    sqrt(as.numeric(fit$var))
  } else {
    NA_real_
  }
  if (!is.finite(se)) {
    return(c(cindex = cindex, ci.low = NA_real_, ci.up = NA_real_))
  }

  ci <- cindex + c(-1, 1) * stats::qnorm(0.975) * se
  c(cindex = cindex, ci.low = max(0, ci[1]), ci.up = min(1, ci[2]))
}

#' @keywords internal
#' @noRd
survival_knn_predict <- function(xtrain, time, status, xtest = NULL, k = 5,
                                 leave_one_out = FALSE) {
  xtrain <- as.matrix(xtrain)
  xtest <- if (is.null(xtest)) xtrain else as.matrix(xtest)
  storage.mode(xtrain) <- "double"
  storage.mode(xtest) <- "double"

  k <- max(1, as.integer(k))
  m <- colMeans(xtrain, na.rm = TRUE)
  s <- apply(xtrain, 2, stats::sd, na.rm = TRUE)
  s[!is.finite(s) | s < 1e-8] <- 1

  xtrain <- scale(xtrain, center = m, scale = s)
  xtest <- scale(xtest, center = m, scale = s)
  xtrain[!is.finite(xtrain)] <- 0
  xtest[!is.finite(xtest)] <- 0

  neighbor_risk <- as.numeric(status) / pmax(as.numeric(time), .Machine$double.eps)
  neighbor_risk[!is.finite(neighbor_risk)] <- 0
  out <- numeric(nrow(xtest))

  for (i in seq_len(nrow(xtest))) {
    d <- rowSums((t(t(xtrain) - xtest[i, ]))^2)
    if (leave_one_out && nrow(xtest) == nrow(xtrain)) d[i] <- Inf
    keep_k <- min(k, sum(is.finite(d)))
    if (keep_k == 0) {
      out[i] <- NA_real_
    } else {
      idx <- order(d)[seq_len(keep_k)]
      out[i] <- mean(neighbor_risk[idx], na.rm = TRUE)
    }
  }
  out
}

#' @keywords internal
#' @noRd
binary_mb_predict <- function(xtrain, outcome, xtest = NULL, mstop = 300, ...) {
  xtrain <- as.data.frame(xtrain)
  xtest <- if (is.null(xtest)) xtrain else as.data.frame(xtest)
  colnames(xtrain) <- make.names(colnames(xtrain), unique = TRUE)
  colnames(xtest) <- colnames(xtrain)

  keep_cols <- vapply(xtrain, function(x) {
    v <- stats::var(as.numeric(x), na.rm = TRUE)
    is.finite(v) && v > 1e-12
  }, logical(1))
  if (!any(keep_cols)) {
    return(rep(mean(as.numeric(outcome), na.rm = TRUE), nrow(xtest)))
  }
  xtrain <- xtrain[, keep_cols, drop = FALSE]
  xtest <- xtest[, keep_cols, drop = FALSE]

  train_df <- data.frame(y = factor(outcome, levels = c(0, 1)), xtrain)
  control_object <- mboost::boost_control(mstop = mstop)
  fit <- tryCatch(
    mboost::glmboost(
      y ~ .,
      data = train_df,
      family = mboost::Binomial(),
      control = control_object,
      ...
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(rep(mean(as.numeric(outcome), na.rm = TRUE), nrow(xtest)))
  }

  as.numeric(stats::predict(fit, newdata = xtest, type = "response"))
}

#' @keywords internal
#' @noRd
normalize_probability_rows <- function(prob_mat) {
  prob_mat <- as.matrix(prob_mat)
  prob_mat[!is.finite(prob_mat)] <- 0
  prob_mat[prob_mat < 0] <- 0
  row_sums <- rowSums(prob_mat)
  bad <- !is.finite(row_sums) | row_sums <= 0
  if (any(bad)) {
    prob_mat[bad, ] <- 1 / ncol(prob_mat)
    row_sums <- rowSums(prob_mat)
  }
  prob_mat / row_sums
}

#' @keywords internal
#' @noRd
multiclass_simple_average <- function(prob_mat, num_classes) {
  prob_mat <- as.matrix(prob_mat)
  if (ncol(prob_mat) == num_classes) {
    return(normalize_probability_rows(prob_mat))
  }
  if (ncol(prob_mat) %% num_classes != 0) {
    stop("Probability matrix column count is not divisible by num_classes.")
  }
  num_models <- ncol(prob_mat) / num_classes
  arr <- array(as.numeric(prob_mat), dim = c(nrow(prob_mat), num_classes, num_models))
  normalize_probability_rows(apply(arr, c(1, 2), mean, na.rm = TRUE))
}

#' @keywords internal
#' @noRd
extract_pathway_genes <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    if (!is.null(names(x)) && any(nzchar(names(x)))) {
      return(unique(names(x)))
    }
    return(unique(unlist(x, use.names = FALSE)))
  }
  unique(as.character(x))
}

#' @keywords internal
#' @noRd
pool_pathway_features <- function(features, literature.features = NULL,
                                  redundancy.thresh = NULL) {
  if (is.null(literature.features) || length(literature.features) == 0) {
    return(lapply(features, extract_pathway_genes))
  }

  out <- lapply(features, extract_pathway_genes)
  if (is.null(names(out))) names(out) <- paste0("Pathway_", seq_along(out))

  lit <- lapply(literature.features, extract_pathway_genes)
  if (is.null(names(lit))) names(lit) <- paste0("Literature_", seq_along(lit))

  existing_sets <- lapply(out, function(x) unique(as.character(x)))
  for (nm in names(lit)) {
    genes <- unique(as.character(lit[[nm]]))
    genes <- genes[nzchar(genes)]
    if (!length(genes)) next

    if (!is.null(redundancy.thresh) && length(existing_sets)) {
      jmax <- max(vapply(existing_sets, function(ref) {
        union_n <- length(union(genes, ref))
        if (union_n == 0) 0 else length(intersect(genes, ref)) / union_n
      }, numeric(1)), na.rm = TRUE)
      if (is.finite(jmax) && jmax > redundancy.thresh) next
    }

    out_name <- nm
    if (out_name %in% names(out)) out_name <- paste0(out_name, "_literature")
    out[[out_name]] <- genes
    existing_sets[[out_name]] <- genes
  }
  out
}

#' @keywords internal
#' @noRd
normalize_celltype_token <- function(x) {
  tolower(gsub("[^[:alnum:]]+", "", x))
}

#' @keywords internal
#' @noRd
validate_celltype_aliases <- function(celltype.aliases, celltypes) {
  if (is.null(celltype.aliases)) return(NULL)
  if (!is.list(celltype.aliases) || is.null(names(celltype.aliases)) ||
      any(!nzchar(names(celltype.aliases)))) {
    stop("celltype.aliases must be a named list keyed by entries in celltypes.")
  }

  unknown <- setdiff(names(celltype.aliases), celltypes)
  if (length(unknown)) {
    stop(
      "celltype.aliases contains names not present in celltypes: ",
      paste(unknown, collapse = ", ")
    )
  }

  lapply(celltype.aliases, function(x) {
    x <- unique(as.character(unlist(x, use.names = FALSE)))
    x[nzchar(x)]
  })
}

#' @keywords internal
#' @noRd
pathway_celltype_matches <- function(pathway.name, celltypes,
                                     celltype.aliases = NULL) {
  if (is.null(celltypes) || !length(celltypes)) return(character(0))
  path_norm <- normalize_celltype_token(pathway.name)

  vapply(celltypes, function(ct) {
    aliases <- ct
    if (!is.null(celltype.aliases) && ct %in% names(celltype.aliases)) {
      aliases <- unique(c(aliases, celltype.aliases[[ct]]))
    }
    any(vapply(aliases, function(alias) {
      grepl(normalize_celltype_token(alias), path_norm, fixed = TRUE)
    }, logical(1)))
  }, logical(1)) |>
    which() |>
    {\(idx) celltypes[idx]}()
}

#' @keywords internal
#' @noRd
celltype_gene_columns <- function(genes, celltype, data.cols,
                                  celltype.sep = c("_", ".", "-", ":", "|"),
                                  celltype.position = c("auto", "suffix", "prefix")) {
  celltype.position <- match.arg(celltype.position)
  genes <- unique(as.character(genes))
  genes <- genes[nzchar(genes)]
  if (!length(genes)) return(character(0))

  candidates <- character(0)
  for (sep in celltype.sep) {
    if (celltype.position %in% c("auto", "suffix")) {
      candidates <- c(candidates, paste0(genes, sep, celltype))
    }
    if (celltype.position %in% c("auto", "prefix")) {
      candidates <- c(candidates, paste0(celltype, sep, genes))
    }
  }
  candidates <- unique(c(candidates, make.names(candidates, unique = FALSE)))
  unique(data.cols[data.cols %in% candidates])
}

#' @keywords internal
#' @noRd
expand_celltype_pathways <- function(features, data.cols, celltypes = NULL,
                                     celltype.sep = c("_", ".", "-", ":", "|"),
                                     celltype.position = c("auto", "suffix", "prefix"),
                                     celltype.aliases = NULL) {
  celltype.position <- match.arg(celltype.position)
  if (is.null(celltypes) || !length(celltypes)) return(features)
  celltype.aliases <- validate_celltype_aliases(celltype.aliases, celltypes)
  if (is.null(names(features))) {
    names(features) <- paste0("Pathway_", seq_along(features))
  }

  out <- list()
  for (nm in names(features)) {
    genes <- extract_pathway_genes(features[[nm]])
    matched_ct <- pathway_celltype_matches(nm, celltypes, celltype.aliases)
    target_ct <- if (length(matched_ct)) matched_ct else celltypes

    for (ct in target_ct) {
      cols <- celltype_gene_columns(genes, ct, data.cols, celltype.sep,
                                    celltype.position)
      if (!length(cols)) next
      out_name <- if (length(matched_ct) == 1 && length(target_ct) == 1) {
        nm
      } else {
        paste0(nm, "_", ct)
      }
      out[[out_name]] <- cols
    }
  }

  if (!length(out)) {
    stop("No pathway genes could be mapped to celltype-specific expression columns.")
  }
  out
}

#' @keywords internal
#' @noRd
align_pathway_prior <- function(pathway.names, pathway.prior = NULL,
                                prior.default = NA_real_,
                                prior.clip = c(0, 1)) {
  prior <- rep(prior.default, length(pathway.names))
  names(prior) <- pathway.names
  if (is.null(pathway.prior)) return(prior)

  pathway.prior <- normalize_pathway_prior(pathway.prior)
  if (!length(pathway.prior)) return(prior)
  if (is.null(names(pathway.prior))) {
    stop("pathway.prior must be a named vector/list or a data.frame with pathway and prior columns.")
  }

  matched <- intersect(pathway.names, names(pathway.prior))
  prior[matched] <- as.numeric(pathway.prior[matched])
  unmatched <- setdiff(pathway.names, matched)
  if (length(unmatched)) {
    prior.names <- names(pathway.prior)
    clean.prior.names <- clean_signature_name(prior.names)
    for (nm in unmatched) {
      clean.nm <- clean_signature_name(nm)
      clean.hit <- match(clean.nm, clean.prior.names)
      if (!is.na(clean.hit)) {
        prior[nm] <- as.numeric(pathway.prior[[clean.hit]])
        next
      }
      hits <- prior.names[vapply(prior.names, function(pr) {
        startsWith(nm, paste0(pr, "_"))
      }, logical(1))]
      if (length(hits)) {
        best <- hits[which.max(nchar(hits))]
        prior[nm] <- as.numeric(pathway.prior[[best]])
      }
    }
  }
  prior[!is.finite(prior)] <- prior.default
  if (!is.null(prior.clip) && length(prior.clip) == 2) {
    prior <- pmin(max(prior.clip), pmax(min(prior.clip), prior))
  }
  names(prior) <- pathway.names
  prior
}

#' @keywords internal
#' @noRd
normalize_pathway_prior <- function(pathway.prior = NULL) {
  if (is.null(pathway.prior)) return(numeric(0))
  if (is.data.frame(pathway.prior)) {
    nm_col <- intersect(c("pathway", "feature", "geneset", "name"), tolower(colnames(pathway.prior)))[1]
    pr_col <- intersect(c("prior", "weight", "pi", "prior_weight"), tolower(colnames(pathway.prior)))[1]
    if (!is.na(nm_col) && !is.na(pr_col)) {
      names(pathway.prior) <- tolower(names(pathway.prior))
      vals <- pathway.prior[[pr_col]]
      names(vals) <- pathway.prior[[nm_col]]
      pathway.prior <- vals
    }
  }
  if (is.list(pathway.prior) && !is.data.frame(pathway.prior)) {
    pathway.prior <- unlist(pathway.prior, use.names = TRUE)
  }
  pathway.prior
}

#' @keywords internal
#' @noRd
merge_pathway_prior <- function(llm.prior = NULL, user.prior = NULL) {
  llm <- normalize_pathway_prior(llm.prior)
  usr <- normalize_pathway_prior(user.prior)
  if (!length(llm)) return(usr)
  if (!length(usr)) return(llm)
  out <- llm
  out[names(usr)] <- as.numeric(usr)
  out
}

#' @keywords internal
#' @noRd
combine_selection_score <- function(metric, prior, prior.weight = 0,
                                    prior.mode = c("none", "multiply", "weighted_sum")) {
  prior.mode <- match.arg(prior.mode)
  metric.names <- names(metric)
  metric <- as.numeric(metric)
  prior <- as.numeric(prior)
  prior.weight <- max(0, min(1, prior.weight))
  if (prior.mode == "none" || prior.weight == 0) {
    names(metric) <- metric.names
    return(metric)
  }

  if (prior.mode == "multiply") {
    out <- metric * ((1 - prior.weight) + prior.weight * prior)
    names(out) <- metric.names
    return(out)
  }

  scale01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (!all(is.finite(rng)) || diff(rng) < 1e-12) return(rep(0.5, length(x)))
    (x - rng[1]) / diff(rng)
  }
  out <- (1 - prior.weight) * scale01(metric) + prior.weight * scale01(prior)
  names(out) <- metric.names
  out
}

#' @keywords internal
#' @noRd
split_iteration_chunks <- function(iter, workers) {
  workers <- max(1, min(as.integer(workers), as.integer(iter)))
  counts <- rep(iter %/% workers, workers)
  counts[seq_len(iter %% workers)] <- counts[seq_len(iter %% workers)] + 1
  offsets <- if (length(counts) == 1) 0 else cumsum(c(0, counts[-length(counts)]))
  data.frame(iter = counts, offset = offsets)
}

#' @keywords internal
#' @noRd
combine_iteration_results <- function(parts, ensemble = FALSE) {
  individual_names <- unique(unlist(lapply(parts, function(x) names(x$individual))))
  individual <- stats::setNames(lapply(individual_names, function(model) {
    unlist(lapply(parts, function(x) x$individual[[model]]), use.names = FALSE)
  }), individual_names)

  ensemble_out <- NULL
  if (ensemble) {
    ensemble_names <- unique(unlist(lapply(parts, function(x) names(x$ensemble))))
    ensemble_out <- stats::setNames(lapply(ensemble_names, function(model) {
      unlist(lapply(parts, function(x) x$ensemble[[model]]), use.names = FALSE)
    }), ensemble_names)
  }

  list(
    individual = individual,
    ensemble = ensemble_out,
    y.pred.avg.pool.all = unlist(lapply(parts, function(x) x$y.pred.avg.pool.all), use.names = FALSE),
    y.pred.avg.pool1.all = unlist(lapply(parts, function(x) x$y.pred.avg.pool1.all), use.names = FALSE)
  )
}

#' @keywords internal
#' @noRd
run_model_iterations <- function(model.fun, iter,
                                 workers = future::availableCores(),
                                 parallel.iter = FALSE, future.seed = NULL,
                                 future.strategy = "multisession",
                                 ensemble = FALSE, ...) {
  if (!isTRUE(parallel.iter) || iter <= 1) {
    return(model.fun(iter = iter, ensemble = ensemble, ..., iter.offset = 0))
  }

  if (is.null(workers)) workers <- future::availableCores()
  workers <- max(1, min(as.integer(workers), as.integer(iter)))
  if (workers <= 1) {
    return(model.fun(iter = iter, ensemble = ensemble, ..., iter.offset = 0))
  }

  chunks <- split_iteration_chunks(iter, workers)
  message("Parallel iterations: ", iter, " total iterations split across ",
          nrow(chunks), " future workers.")

  if (!is.null(future.strategy)) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    do.call(future::plan, c(list(future.strategy), list(workers = nrow(chunks))))
  }

  args <- list(...)
  parts <- future.apply::future_lapply(
    seq_len(nrow(chunks)),
    function(k) {
      do.call(model.fun, c(
        args,
        list(iter = chunks$iter[k], ensemble = ensemble,
             iter.offset = chunks$offset[k])
      ))
    },
    future.seed = future.seed
  )

  combine_iteration_results(parts, ensemble = ensemble)
}

#' @keywords internal
#' @noRd
combine_train_performance_chunks <- function(parts, chunk.iter) {
  if (length(parts) == 1) return(parts[[1]])
  models <- unique(unlist(lapply(parts, function(x) as.character(x$model))))
  out <- lapply(models, function(model) {
    rows <- lapply(parts, function(x) x[as.character(x$model) == model, , drop = FALSE])
    keep <- vapply(rows, nrow, integer(1)) > 0
    rows <- rows[keep]
    weights <- chunk.iter[keep]
    template <- rows[[1]][1, , drop = FALSE]
    numeric.cols <- names(template)[vapply(template, is.numeric, logical(1))]
    for (col in numeric.cols) {
      vals <- vapply(rows, function(r) as.numeric(r[[col]][1]), numeric(1))
      ok <- is.finite(vals)
      template[[col]] <- if (any(ok)) stats::weighted.mean(vals[ok], weights[ok]) else NA_real_
    }
    template
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

#' @keywords internal
#' @noRd
run_train_performance_iterations <- function(perf.fun, iter,
                                             workers = future::availableCores(),
                                             parallel.iter = FALSE,
                                             future.seed = NULL,
                                             future.strategy = "multisession",
                                             ...) {
  if (!isTRUE(parallel.iter) || iter <= 1) {
    return(perf.fun(iter = iter, ..., iter.offset = 0))
  }

  if (is.null(workers)) workers <- future::availableCores()
  workers <- max(1, min(as.integer(workers), as.integer(iter)))
  if (workers <= 1) {
    return(perf.fun(iter = iter, ..., iter.offset = 0))
  }

  chunks <- split_iteration_chunks(iter, workers)
  message("Parallel Train_perform_CV iterations: ", iter,
          " total iterations split across ", nrow(chunks),
          " future workers.")

  if (!is.null(future.strategy)) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    do.call(future::plan, c(list(future.strategy), list(workers = nrow(chunks))))
  }

  args <- list(...)
  parts <- future.apply::future_lapply(
    seq_len(nrow(chunks)),
    function(k) {
      do.call(perf.fun, c(
        args,
        list(iter = chunks$iter[k], iter.offset = chunks$offset[k])
      ))
    },
    future.seed = future.seed
  )
  combine_train_performance_chunks(parts, chunks$iter)
}

#' @keywords internal
#' @noRd
importance_percentile <- function(importance) {
  original_names <- names(importance) %||% paste0("V", seq_along(importance))
  importance <- abs(as.numeric(importance))
  names(importance) <- original_names
  importance[!is.finite(importance)] <- 0
  if (length(importance) == 1) {
    out <- 1
    names(out) <- names(importance)
    return(out)
  }
  if (diff(range(importance, na.rm = TRUE)) < 1e-12) {
    out <- rep(0.5, length(importance))
    names(out) <- names(importance)
    return(out)
  }
  ranks <- rank(importance, ties.method = "average")
  out <- (ranks - 1) / (length(importance) - 1)
  names(out) <- names(importance)
  out
}

#' @keywords internal
#' @noRd
gene_importance_summary <- function(raw, genes, stability.quantile = 0.75) {
  rows <- lapply(genes, function(gene) {
    vals <- raw$percentile[raw$gene == gene]
    if (!length(vals)) vals <- NA_real_
    data.frame(
      gene = gene,
      median_percentile = stats::median(vals, na.rm = TRUE),
      mean_percentile = mean(vals, na.rm = TRUE),
      stability = mean(vals >= stability.quantile, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)
  tab$median_percentile[!is.finite(tab$median_percentile)] <- 0
  tab$mean_percentile[!is.finite(tab$mean_percentile)] <- 0
  tab$stability[!is.finite(tab$stability)] <- 0
  tab
}

#' @keywords internal
#' @noRd
native_gene_importance <- function(outcome.type, model, xtrain, ytrain,
                                   survdays = NULL, params = list(),
                                   suppress.warnings = TRUE) {
  original_names <- colnames(xtrain)
  xtrain <- as.data.frame(xtrain)
  colnames(xtrain) <- make.names(original_names, unique = TRUE)
  name_map <- setNames(original_names, colnames(xtrain))
  keep_cols <- vapply(xtrain, function(x) {
    v <- stats::var(as.numeric(x), na.rm = TRUE)
    is.finite(v) && v > 1e-12
  }, logical(1))
  out <- rep(0, ncol(xtrain))
  names(out) <- colnames(xtrain)
  if (!any(keep_cols)) {
    names(out) <- name_map[names(out)]
    return(out)
  }
  xtrain <- xtrain[, keep_cols, drop = FALSE]

  fit_call <- function() {
    switch(
      toupper(model),
      XG = native_importance_xg(outcome.type, xtrain, ytrain, survdays, params),
      RF = native_importance_rf(outcome.type, xtrain, ytrain, survdays, params),
      GB = native_importance_gb(outcome.type, xtrain, ytrain, survdays, params),
      EN = native_importance_en(outcome.type, xtrain, ytrain, survdays, params),
      MB = native_importance_mb(outcome.type, xtrain, ytrain, survdays, params),
      DCT = native_importance_dct(outcome.type, xtrain, ytrain, survdays, params),
      NULL
    )
  }
  imp <- tryCatch(
    if (isTRUE(suppress.warnings)) suppressWarnings(fit_call()) else fit_call(),
    error = function(e) NULL
  )
  if (is.null(imp)) return(NULL)
  out[names(imp)] <- abs(as.numeric(imp))
  names(out) <- name_map[names(out)]
  out
}

#' @keywords internal
#' @noRd
native_importance_xg <- function(outcome.type, xtrain, ytrain, survdays, params) {
  nrounds <- params$nrounds %||% 50
  params$nrounds <- NULL
  if (outcome.type == "binary") {
    pars <- utils::modifyList(
      list(objective = "binary:logistic", max_depth = 2, eta = 0.1,
           verbosity = 0, nthread = 1),
      params
    )
    dtrain <- xgboost::xgb.DMatrix(as.matrix(xtrain), label = as.numeric(ytrain))
  } else if (outcome.type == "survival") {
    pars <- utils::modifyList(
      list(objective = "survival:cox", eta = 0.1, verbosity = 0,
           nthread = 1),
      params
    )
    label <- ifelse(as.numeric(ytrain) == 1, survdays, -survdays)
    dtrain <- xgboost::xgb.DMatrix(as.matrix(xtrain), label = label)
  } else {
    cls <- as.integer(as.factor(ytrain)) - 1
    pars <- utils::modifyList(
      list(objective = "multi:softprob", num_class = length(unique(cls)),
           max_depth = 2, eta = 0.1, verbosity = 0, nthread = 1),
      params
    )
    dtrain <- xgboost::xgb.DMatrix(as.matrix(xtrain), label = cls)
  }
  fit <- xgboost::xgb.train(params = pars, data = dtrain, nrounds = nrounds)
  tab <- xgboost::xgb.importance(model = fit)
  imp <- setNames(rep(0, ncol(xtrain)), colnames(xtrain))
  if (nrow(tab)) imp[tab$Feature] <- tab$Gain
  imp
}

#' @keywords internal
#' @noRd
native_importance_rf <- function(outcome.type, xtrain, ytrain, survdays, params) {
  if (outcome.type == "survival") {
    dat <- data.frame(time = survdays, status = ytrain, xtrain)
    params <- utils::modifyList(list(ntree = 300, importance = TRUE), params)
    fit <- do.call(
      randomForestSRC::rfsrc,
      c(list(survival::Surv(time, status) ~ ., data = dat), params)
    )
    imp <- fit$importance
    return(setNames(abs(as.numeric(imp[colnames(xtrain)])), colnames(xtrain)))
  }
  params <- utils::modifyList(list(ntree = 300, importance = TRUE), params)
  fit <- do.call(
    randomForest::randomForest,
    c(list(x = xtrain, y = as.factor(ytrain)), params)
  )
  imp <- randomForest::importance(fit)
  if (is.matrix(imp)) imp <- rowMeans(abs(imp), na.rm = TRUE)
  setNames(abs(as.numeric(imp[colnames(xtrain)])), colnames(xtrain))
}

#' @keywords internal
#' @noRd
native_importance_gb <- function(outcome.type, xtrain, ytrain, survdays, params) {
  dat <- data.frame(y = ytrain, time = survdays, xtrain)
  if (outcome.type == "binary") {
    params <- utils::modifyList(
      list(distribution = "bernoulli", n.trees = 100,
           interaction.depth = 2, shrinkage = 0.05, n.cores = 1,
           verbose = FALSE),
      params
    )
    fit <- do.call(gbm::gbm, c(list(formula = y ~ ., data = dat), params))
  } else if (outcome.type == "survival") {
    params <- utils::modifyList(
      list(distribution = "coxph", n.trees = 100,
           interaction.depth = 2, shrinkage = 0.05, n.cores = 1,
           verbose = FALSE),
      params
    )
    fit <- do.call(
      gbm::gbm,
      c(list(formula = survival::Surv(time, y) ~ .,
             data = dat[, c("time", "y", colnames(xtrain)), drop = FALSE]),
        params)
    )
  } else {
    dat$y <- as.factor(dat$y)
    params <- utils::modifyList(
      list(distribution = "multinomial", n.trees = 100,
           interaction.depth = 2, shrinkage = 0.05, n.cores = 1,
           verbose = FALSE),
      params
    )
    fit <- do.call(gbm::gbm, c(list(formula = y ~ ., data = dat), params))
  }
  tab <- summary(fit, plotit = FALSE)
  imp <- setNames(rep(0, ncol(xtrain)), colnames(xtrain))
  if (nrow(tab)) imp[as.character(tab$var)] <- tab$rel.inf
  imp
}

#' @keywords internal
#' @noRd
native_importance_en <- function(outcome.type, xtrain, ytrain, survdays, params) {
  s.val <- params$s %||% "lambda.1se"
  params$s <- NULL
  if (outcome.type == "binary") {
    params <- utils::modifyList(list(alpha = 0.2, family = "binomial", nfolds = 3), params)
    fit <- do.call(glmnet::cv.glmnet, c(list(x = as.matrix(xtrain), y = ytrain), params))
    cf <- as.matrix(stats::coef(fit, s = s.val))
    imp <- abs(cf[setdiff(rownames(cf), "(Intercept)"), 1])
  } else if (outcome.type == "survival") {
    params <- utils::modifyList(list(alpha = 0.2, family = "cox", nfolds = 3), params)
    ysurv <- survival::Surv(survdays, ytrain)
    fit <- do.call(glmnet::cv.glmnet, c(list(x = as.matrix(xtrain), y = ysurv), params))
    cf <- as.matrix(stats::coef(fit, s = params$s %||% "lambda.min"))
    imp <- abs(cf[, 1])
  } else {
    params <- utils::modifyList(list(alpha = 0.2, family = "multinomial", nfolds = 3), params)
    fit <- do.call(glmnet::cv.glmnet, c(list(x = as.matrix(xtrain), y = as.factor(ytrain)), params))
    cf <- stats::coef(fit, s = s.val)
    mat <- do.call(cbind, lapply(cf, function(z) {
      z <- as.matrix(z)
      abs(z[setdiff(rownames(z), "(Intercept)"), 1])
    }))
    imp <- rowMeans(mat, na.rm = TRUE)
  }
  setNames(as.numeric(imp[colnames(xtrain)]), colnames(xtrain))
}

#' @keywords internal
#' @noRd
native_importance_mb <- function(outcome.type, xtrain, ytrain, survdays, params) {
  mstop_val <- params$mstop %||% 100
  params$mstop <- NULL
  ctrl <- mboost::boost_control(mstop = mstop_val)
  if (outcome.type == "binary") {
    dat <- data.frame(y = factor(ytrain, levels = c(0, 1)), xtrain)
    fit <- do.call(
      mboost::glmboost,
      c(list(formula = y ~ ., data = dat, family = mboost::Binomial(),
             control = ctrl), params)
    )
    return(mb_coef_importance(stats::coef(fit, off2int = TRUE), colnames(xtrain)))
  }
  if (outcome.type == "survival") {
    dat <- data.frame(time = survdays, status = ytrain, xtrain)
    fit <- do.call(
      mboost::glmboost,
      c(list(formula = survival::Surv(time, status) ~ ., data = dat,
             family = mboost::CoxPH(), control = ctrl), params)
    )
    return(mb_coef_importance(stats::coef(fit, off2int = TRUE), colnames(xtrain)))
  }
  classes <- sort(unique(ytrain))
  imps <- lapply(classes, function(cls) {
    y01 <- as.integer(ytrain == cls)
    if (length(unique(y01)) < 2) {
      return(setNames(rep(0, ncol(xtrain)), colnames(xtrain)))
    }
    dat <- data.frame(y = factor(y01, levels = c(0, 1)), xtrain)
    fit <- do.call(
      mboost::glmboost,
      c(list(formula = y ~ ., data = dat, family = mboost::Binomial(),
             control = ctrl), params)
    )
    mb_coef_importance(stats::coef(fit, off2int = TRUE), colnames(xtrain))
  })
  Reduce("+", imps) / length(imps)
}

#' @keywords internal
#' @noRd
mb_coef_importance <- function(coefs, genes) {
  imp <- setNames(rep(0, length(genes)), genes)
  if (!length(coefs)) return(imp)
  nm <- names(coefs)
  for (g in genes) {
    hit <- nm == g | grepl(paste0("^", g, "\\b"), nm)
    if (any(hit)) imp[g] <- sum(abs(as.numeric(coefs[hit])), na.rm = TRUE)
  }
  imp
}

#' @keywords internal
#' @noRd
native_importance_dct <- function(outcome.type, xtrain, ytrain, survdays, params) {
  if (outcome.type == "survival") {
    dat <- data.frame(time = survdays, status = ytrain, xtrain)
    control_args <- utils::modifyList(list(cp = 0.001, minsplit = 10), params)
    fit <- rpart::rpart(
      survival::Surv(time, status) ~ .,
      data = dat,
      method = "exp",
      control = do.call(rpart::rpart.control, control_args)
    )
  } else {
    dat <- data.frame(y = as.factor(ytrain), xtrain)
    method <- params$method %||% "class"
    params$method <- NULL
    control_args <- utils::modifyList(list(cp = 0.001), params)
    fit <- rpart::rpart(
      y ~ .,
      data = dat,
      method = method,
      control = do.call(rpart::rpart.control, control_args)
    )
  }
  imp <- setNames(rep(0, ncol(xtrain)), colnames(xtrain))
  if (!is.null(fit$variable.importance)) {
    vi <- fit$variable.importance
    imp[names(vi)] <- as.numeric(vi)
  }
  imp
}

#' @keywords internal
#' @noRd
multiclass_adb_predict <- function(xtrain, outcome, xtest = NULL, classes = NULL,
                                   tree_depth = 3, n_rounds = 200,
                                   verbose = FALSE, ...) {
  xtrain <- as.data.frame(xtrain)
  xtest <- if (is.null(xtest)) xtrain else as.data.frame(xtest)
  colnames(xtrain) <- make.names(colnames(xtrain), unique = TRUE)
  colnames(xtest) <- colnames(xtrain)

  classes <- classes %||% sort(unique(outcome))
  out <- matrix(NA_real_, nrow = nrow(xtest), ncol = length(classes))
  colnames(out) <- paste0("C", classes)

  keep_cols <- vapply(xtrain, function(x) {
    v <- stats::var(as.numeric(x), na.rm = TRUE)
    is.finite(v) && v > 1e-12
  }, logical(1))
  if (!any(keep_cols)) {
    out[,] <- 1 / length(classes)
    return(out)
  }
  xtrain <- xtrain[, keep_cols, drop = FALSE]
  xtest <- xtest[, keep_cols, drop = FALSE]

  for (k in seq_along(classes)) {
    y01 <- as.integer(outcome == classes[k])
    if (length(unique(y01)) < 2) {
      out[, k] <- mean(y01, na.rm = TRUE)
      next
    }
    fit <- tryCatch(
      JOUSBoost::adaboost(
        X = as.matrix(xtrain),
        y = ifelse(y01 == 1, 1, -1),
        tree_depth = tree_depth,
        n_rounds = n_rounds,
        verbose = verbose,
        ...
      ),
      error = function(e) NULL
    )
    out[, k] <- if (is.null(fit)) {
      mean(y01, na.rm = TRUE)
    } else {
      as.numeric(stats::predict(fit, as.matrix(xtest), type = "prob"))
    }
  }

  normalize_probability_rows(out)
}

#' @keywords internal
#' @noRd
multiclass_mb_predict <- function(xtrain, outcome, xtest = NULL, classes = NULL,
                                  mstop = 300, ...) {
  xtrain <- as.data.frame(xtrain)
  xtest <- if (is.null(xtest)) xtrain else as.data.frame(xtest)
  colnames(xtrain) <- make.names(colnames(xtrain), unique = TRUE)
  colnames(xtest) <- colnames(xtrain)

  classes <- classes %||% sort(unique(outcome))
  out <- matrix(NA_real_, nrow = nrow(xtest), ncol = length(classes))
  colnames(out) <- paste0("C", classes)

  keep_cols <- vapply(xtrain, function(x) {
    v <- stats::var(as.numeric(x), na.rm = TRUE)
    is.finite(v) && v > 1e-12
  }, logical(1))
  if (!any(keep_cols)) {
    out[,] <- 1 / length(classes)
    return(out)
  }
  xtrain <- xtrain[, keep_cols, drop = FALSE]
  xtest <- xtest[, keep_cols, drop = FALSE]

  for (k in seq_along(classes)) {
    y01 <- as.integer(outcome == classes[k])
    if (length(unique(y01)) < 2) {
      out[, k] <- mean(y01, na.rm = TRUE)
      next
    }
    train_df <- data.frame(y = factor(y01, levels = c(0, 1)), xtrain)
    fit <- tryCatch(
      mboost::glmboost(
        y ~ .,
        data = train_df,
        family = mboost::Binomial(),
        control = mboost::boost_control(mstop = mstop),
        ...
      ),
      error = function(e) NULL
    )
    out[, k] <- if (is.null(fit)) {
      mean(y01, na.rm = TRUE)
    } else {
      as.numeric(stats::predict(fit, newdata = xtest, type = "response"))
    }
  }

  normalize_probability_rows(out)
}

#' @keywords internal
#' @noRd
survival_adb_predict <- function(xtrain, time, status, time_point = NULL,
                                 xtest = NULL, tree_depth = 3,
                                 n_rounds = 200, verbose = FALSE, ...) {
  xtrain_all <- as.data.frame(xtrain)
  xtest <- if (is.null(xtest)) xtrain_all else as.data.frame(xtest)
  colnames(xtrain_all) <- make.names(colnames(xtrain_all), unique = TRUE)
  colnames(xtest) <- colnames(xtrain_all)

  if (is.null(time_point) || !is.finite(time_point)) {
    time_point <- stats::median(time[status == 1], na.rm = TRUE)
    if (!is.finite(time_point)) time_point <- stats::median(time, na.rm = TRUE)
  }

  usable <- (status == 1 & time <= time_point) | (time > time_point)
  if (sum(usable, na.rm = TRUE) < 2) {
    return(rep(NA_real_, nrow(xtest)))
  }

  xtrain <- xtrain_all[usable, , drop = FALSE]
  y01 <- ifelse(status[usable] == 1 & time[usable] <= time_point, 1, 0)
  if (length(unique(y01)) < 2) {
    return(rep(mean(y01, na.rm = TRUE), nrow(xtest)))
  }

  keep_cols <- vapply(xtrain, function(x) stats::var(as.numeric(x), na.rm = TRUE) > 1e-12, logical(1))
  if (!any(keep_cols)) {
    return(rep(mean(y01, na.rm = TRUE), nrow(xtest)))
  }
  xtrain <- xtrain[, keep_cols, drop = FALSE]
  xtest <- xtest[, keep_cols, drop = FALSE]

  y_adb <- ifelse(y01 == 1, 1, -1)
  fit <- tryCatch(
    JOUSBoost::adaboost(
      X = as.matrix(xtrain),
      y = y_adb,
      tree_depth = tree_depth,
      n_rounds = n_rounds,
      verbose = verbose,
      ...
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(rep(mean(y01, na.rm = TRUE), nrow(xtest)))
  }

  as.numeric(stats::predict(fit, as.matrix(xtest), type = "prob"))
}

#' @keywords internal
#' @noRd
survival_nn_predict <- function(xtrain, time, status, xtest = NULL, size = 5,
                                decay = 1e-4, maxit = 200, trace = FALSE,
                                MaxNWts = 10000, ...) {
  xtrain <- as.matrix(xtrain)
  xtest <- if (is.null(xtest)) xtrain else as.matrix(xtest)
  storage.mode(xtrain) <- "double"
  storage.mode(xtest) <- "double"

  keep_cols <- apply(xtrain, 2, function(x) stats::var(x, na.rm = TRUE) > 1e-12)
  if (!any(keep_cols)) {
    return(rep(NA_real_, nrow(xtest)))
  }
  xtrain <- xtrain[, keep_cols, drop = FALSE]
  xtest <- xtest[, keep_cols, drop = FALSE]

  m <- colMeans(xtrain, na.rm = TRUE)
  s <- apply(xtrain, 2, stats::sd, na.rm = TRUE)
  s[!is.finite(s) | s < 1e-8] <- 1

  xtrain <- scale(xtrain, center = m, scale = s)
  xtest <- scale(xtest, center = m, scale = s)
  xtrain[!is.finite(xtrain)] <- 0
  xtest[!is.finite(xtest)] <- 0

  risk <- as.numeric(status) / pmax(as.numeric(time), .Machine$double.eps)
  risk[!is.finite(risk)] <- 0
  y_center <- mean(risk, na.rm = TRUE)
  y_scale <- stats::sd(risk, na.rm = TRUE)
  if (!is.finite(y_scale) || y_scale < 1e-8 || length(unique(risk)) < 2) {
    return(rep(y_center, nrow(xtest)))
  }
  y_scaled <- as.numeric(scale(risk, center = y_center, scale = y_scale))

  fit <- tryCatch(
    nnet::nnet(
      x = xtrain,
      y = y_scaled,
      size = size,
      linout = TRUE,
      decay = decay,
      maxit = maxit,
      trace = trace,
      MaxNWts = MaxNWts,
      ...
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(rep(NA_real_, nrow(xtest)))
  }

  pred <- as.numeric(stats::predict(fit, xtest))
  pred * y_scale + y_center
}

#' @keywords internal
#' @noRd
binary_sp_ci <- function(outcome, marker, specificity = 0.95, nboot = 200, seed = 123) {
  keep <- stats::complete.cases(outcome, marker) &
    is.finite(outcome) & is.finite(marker)
  outcome <- outcome[keep]
  marker <- marker[keep]
  if (length(outcome) < 2 || length(unique(outcome)) < 2 ||
      length(unique(marker)) < 2) {
    return(c(ci.low = NA_real_, ci.up = NA_real_))
  }

  set.seed(seed)
  vals <- rep(NA_real_, nboot)
  n <- length(outcome)
  for (i in seq_len(nboot)) {
    idx <- sample.int(n, n, replace = TRUE)
    if (length(unique(outcome[idx])) < 2 || length(unique(marker[idx])) < 2) next
    roc_obj <- tryCatch(
      pROC::roc(outcome[idx], marker[idx], levels = c(0, 1), direction = "<", quiet = TRUE),
      error = function(e) NULL
    )
    if (!is.null(roc_obj)) {
      vals[i] <- as.numeric(pROC::coords(
        roc_obj, x = specificity, input = "specificity",
        ret = "sensitivity", transpose = FALSE
      ))
    }
  }
  if (all(is.na(vals))) return(c(ci.low = NA_real_, ci.up = NA_real_))
  ci <- stats::quantile(vals, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  c(ci.low = max(0, ci[1]), ci.up = min(1, ci[2]))
}

#' @keywords internal
#' @noRd
survival_sp_ci <- function(time, status, marker, time_point, specificity = 0.95,
                           nboot = 200, seed = 123) {
  keep <- stats::complete.cases(time, status, marker) &
    is.finite(time) & is.finite(status) & is.finite(marker)
  time <- time[keep]
  status <- status[keep]
  marker <- marker[keep]
  if (length(time) < 2 || length(unique(status)) < 2 || length(unique(marker)) < 2) {
    return(c(ci.low = NA_real_, ci.up = NA_real_))
  }

  set.seed(seed)
  vals <- rep(NA_real_, nboot)
  n <- length(time)
  for (i in seq_len(nboot)) {
    idx <- sample.int(n, n, replace = TRUE)
    if (length(unique(status[idx])) < 2 || length(unique(marker[idx])) < 2) next
    roc_obj <- tryCatch(
      survivalROC::survivalROC(
        Stime = time[idx], status = status[idx], marker = marker[idx],
        predict.time = time_point, method = "KM"
      ),
      error = function(e) NULL
    )
    if (!is.null(roc_obj) && length(roc_obj$TP) > 0) {
      vals[i] <- roc_obj$TP[which.min(abs(1 - roc_obj$FP - specificity))]
    }
  }
  if (all(is.na(vals))) return(c(ci.low = NA_real_, ci.up = NA_real_))
  ci <- stats::quantile(vals, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  c(ci.low = max(0, ci[1]), ci.up = min(1, ci[2]))
}

#' @keywords internal
#' @noRd
survival_auc_ci <- function(time, status, marker, time_point,
                            nboot = 200, seed = 123) {
  keep <- stats::complete.cases(time, status, marker) &
    is.finite(time) & is.finite(status) & is.finite(marker)
  time <- time[keep]
  status <- status[keep]
  marker <- marker[keep]
  if (length(time) < 2 || length(unique(status)) < 2 || length(unique(marker)) < 2) {
    return(c(ci.low = NA_real_, ci.up = NA_real_))
  }

  set.seed(seed)
  vals <- rep(NA_real_, nboot)
  n <- length(time)
  for (i in seq_len(nboot)) {
    idx <- sample.int(n, n, replace = TRUE)
    if (length(unique(status[idx])) < 2 || length(unique(marker[idx])) < 2) next
    roc_obj <- tryCatch(
      survivalROC::survivalROC(
        Stime = time[idx], status = status[idx], marker = marker[idx],
        predict.time = time_point, method = "KM"
      ),
      error = function(e) NULL
    )
    if (!is.null(roc_obj) && is.finite(roc_obj$AUC)) {
      vals[i] <- pmin(1, pmax(0, ifelse(roc_obj$AUC < 0.5, 1 - roc_obj$AUC, roc_obj$AUC)))
    }
  }
  if (all(is.na(vals))) return(c(ci.low = NA_real_, ci.up = NA_real_))
  ci <- stats::quantile(vals, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  c(ci.low = max(0, ci[1]), ci.up = min(1, ci[2]))
}

#' @keywords internal
#' @noRd
multiclass_sp95 <- function(outcome, pred_matrix, specificity = 0.95) {
  classes <- levels(as.factor(outcome))
  if (is.null(colnames(pred_matrix))) colnames(pred_matrix) <- classes
  vals <- vapply(classes, function(cls) {
    y_bin <- as.integer(outcome == cls)
    score <- pred_matrix[, cls]
    keep <- stats::complete.cases(y_bin, score)
    if (length(unique(y_bin[keep])) < 2 || length(unique(score[keep])) < 2) return(NA_real_)
    roc_obj <- tryCatch(
      pROC::roc(y_bin[keep], score[keep], levels = c(0, 1), direction = "<", quiet = TRUE),
      error = function(e) NULL
    )
    if (is.null(roc_obj)) return(NA_real_)
    as.numeric(pROC::coords(
      roc_obj, x = specificity, input = "specificity",
      ret = "sensitivity", transpose = FALSE
    ))
  }, numeric(1))
  mean(vals, na.rm = TRUE)
}

#' @keywords internal
#' @noRd
categorical_sp_ci <- function(outcome, pred_matrix, specificity = 0.95,
                              nboot = 200, seed = 123) {
  keep <- stats::complete.cases(outcome, pred_matrix)
  outcome <- outcome[keep]
  pred_matrix <- pred_matrix[keep, , drop = FALSE]
  if (length(outcome) < 2 || length(unique(outcome)) < 2) {
    return(c(ci.low = NA_real_, ci.up = NA_real_))
  }

  set.seed(seed)
  vals <- rep(NA_real_, nboot)
  n <- length(outcome)
  for (i in seq_len(nboot)) {
    idx <- sample.int(n, n, replace = TRUE)
    if (length(unique(outcome[idx])) < 2) next
    vals[i] <- multiclass_sp95(outcome[idx], pred_matrix[idx, , drop = FALSE], specificity)
  }
  if (all(is.na(vals))) return(c(ci.low = NA_real_, ci.up = NA_real_))
  ci <- stats::quantile(vals, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  c(ci.low = ci[1], ci.up = ci[2])
}
