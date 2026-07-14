# Internal Functions: Not exported
fix_auc <- function(x) {
  pmin(1, pmax(0, ifelse(x < 0.5, 1 - x, x)))
}
categorical_roc <- function(result, outcome, models.indiv, iter, ensemble = FALSE, models.ens = NULL) {

  safe_num <- function(x) {
    if (is.numeric(x)) return(x)
    if (is.list(x)) return(as.numeric(unlist(x)))
    return(as.numeric(x))
  }
  get_avg_pool <- function(model_name, type = c("individual", "ensemble")) {
    type <- match.arg(type)
    pool <- lapply(1:iter, function(i) {
      start <- (i-1) * length(outcome) * length(unique(outcome)) + 1
      end   <- i * length(outcome) * length(unique(outcome))
      matrix(result[[type]][[model_name]][start:end],
             nrow = length(outcome),
             ncol = length(unique(outcome)),
             byrow = TRUE)
    })
    apply(simplify2array(pool), c(1,2), mean)
  }

  # Initialize prediction pools for individual models
  y.pred <- list()
  for (m in models.indiv) y.pred[[m]] <- get_avg_pool(m, type = "individual")

  # Ensemble prediction pools
  y.ens.pred <- list()
  if (ensemble && !is.null(models.ens)) {
    for (m in models.ens) y.ens.pred[[m]] <- get_avg_pool(m, type = "ensemble")
  }

  # SimpleAvg and WeightedAvg predictions
  y.pred.avg.pool.all.W <- apply(simplify2array(lapply(1:iter, function(i) {
    start <- (i-1) * length(outcome) * length(unique(outcome)) + 1
    end   <- i * length(outcome) * length(unique(outcome))
    matrix(result$y.pred.avg.pool.all[start:end],
           nrow = length(outcome),
           ncol = length(unique(outcome)), byrow = TRUE)
  })), c(1,2), mean)

  y.pred.avg.pool1.all.S <- apply(simplify2array(lapply(1:iter, function(i) {
    start <- (i-1) * length(outcome) * length(unique(outcome)) + 1
    end   <- i * length(outcome) * length(unique(outcome))
    matrix(result$y.pred.avg.pool1.all[start:end],
           nrow = length(outcome),
           ncol = length(unique(outcome)), byrow = TRUE)
  })), c(1,2), mean)

  # Compute ROC for a prediction matrix
  compute_roc <- function(pred_matrix) {
    mat <- matrix(pred_matrix, ncol = length(unique(outcome)), byrow = TRUE)
    colnames(mat) <- levels(outcome)
    roc_obj <- multiclass.roc(outcome, mat)
    roc_obj$ci <- get_multiclass_auc_ci(outcome, mat, nboot = 200)
    return(roc_obj)
  }

  # Compute ROC for all individual models
  roc_list <- lapply(names(y.pred), function(m) compute_roc(y.pred[[m]]))
  names(roc_list) <- names(y.pred)

  # Compute ROC for ensemble models
  roc_ens_list <- list()
  if (ensemble && length(y.ens.pred) > 0) {
    roc_ens_list <- lapply(names(y.ens.pred), function(m) compute_roc(y.ens.pred[[m]]))
    names(roc_ens_list) <- names(y.ens.pred)
  }

  # SimpleAvg and WeightedAvg
  roc_list$SimpleAvg   <- compute_roc(y.pred.avg.pool1.all.S)
  roc_list$WeightedAvg <- compute_roc(y.pred.avg.pool.all.W)

  # Prepare model names
  all_models <- c(names(roc_list), if (ensemble) paste0(names(roc_ens_list), ".ens"))

  # Initialize empty named vectors
  roc.tab    <- numeric()
  ci.low.tab <- numeric()
  ci.up.tab  <- numeric()
  sp.98.tab  <- numeric()
  sp.95.tab  <- numeric()
  sp.95.ci.low.tab <- numeric()
  sp.95.ci.up.tab <- numeric()

  append_metrics <- function(roc_obj, model_name) {
    auc_val <- safe_num(roc_obj$auc)[1]
    ci_vals <- safe_num(roc_obj$ci)

    roc.tab[model_name]    <<- fix_auc(auc_val)
    ci.low.tab[model_name] <<- fix_auc(ci_vals[1])
    ci.up.tab[model_name]  <<- fix_auc(ci_vals[3])

    sens_98 <- extract_coords(roc_obj, 0.98)
    sens_95 <- extract_coords(roc_obj, 0.95)

    sp.98.tab[model_name] <<- sens_98
    sp.95.tab[model_name] <<- sens_95
  }

  # Apply metrics to all models
  for (m in names(roc_list)) append_metrics(roc_list[[m]], m)
  if (ensemble) for (m in names(roc_ens_list)) append_metrics(roc_ens_list[[m]], paste0(m, ".ens"))

  # Assign results
  AUC.mat    <- roc.tab
  CI.low.mat <- ci.low.tab
  CI.up.mat  <- ci.up.tab
  SP.98.mat  <- sp.98.tab
  SP.95.mat  <- sp.95.tab
  pred_mats <- c(
    y.pred,
    list(SimpleAvg = y.pred.avg.pool1.all.S, WeightedAvg = y.pred.avg.pool.all.W),
    if (ensemble) stats::setNames(y.ens.pred, paste0(names(y.ens.pred), ".ens")) else list()
  )
  sp95_ci <- t(vapply(
    names(SP.95.mat),
    function(model_name) categorical_sp_ci(outcome, pred_mats[[model_name]]),
    numeric(2)
  ))

  ret_list <- list(
    auc    = AUC.mat,
    ci.low = CI.low.mat,
    ci.up  = CI.up.mat,
    sp.98  = SP.98.mat,
    sp.95  = SP.95.mat,
    sp.95.ci.low = sp95_ci[, "ci.low"],
    sp.95.ci.up = sp95_ci[, "ci.up"],
    roc = c(roc_list, if (ensemble) stats::setNames(roc_ens_list, paste0(names(roc_ens_list), ".ens")) else list()),
    predictions = pred_mats
  )

  return(ret_list)
}
binary_roc <- function(outcome, result, models.indiv, iter, ensemble = FALSE, models.ens = NULL) {
  # ----- Compute average predictions for individual models -----
  y.pred.avg.pool <- list()
  for (m in intersect(models.indiv, names(result$individual))) {
    mat <- result$individual[[m]]
    if (!is.null(mat)) {
      y.pred.avg.pool[[m]] <- apply(matrix(mat, length(outcome), iter), 1, mean, na.rm = TRUE)
    }
  }

  # ----- Compute average predictions for ensemble models if needed -----
  y.ens.avg.pool <- list()
  if (ensemble && !is.null(models.ens)) {
    for (m in intersect(models.ens, names(result$ensemble))) {
      mat <- result$ensemble[[m]]
      if (!is.null(mat)) {
        y.ens.avg.pool[[m]] <- apply(matrix(mat, length(outcome), iter), 1, mean, na.rm = TRUE)
      }
    }
  }

  # ----- Compute predictions for all-model averages -----
  y.pred.avg.pool.all.W <- apply(matrix(result$y.pred.avg.pool.all, length(outcome), iter), 1, mean, na.rm = TRUE)
  y.pred.avg.pool1.all.S <- apply(matrix(result$y.pred.avg.pool1.all, length(outcome), iter), 1, mean, na.rm = TRUE)

  # ----- Compute ROC curves for individual and ensemble models -----
  roc.list <- list()

  for (m in names(y.pred.avg.pool)) {
    roc.list[[m]] <- roc(outcome, y.pred.avg.pool[[m]], ci = TRUE, direction = "<")
  }

  for (m in names(y.ens.avg.pool)) {
    roc.list[[paste0(m, ".ens")]] <- roc(outcome, y.ens.avg.pool[[m]], ci = TRUE, direction = "<")
  }

  # Add the overall predictions
  roc.list[["SimpleAvg"]] <- roc(outcome, y.pred.avg.pool1.all.S, ci = TRUE, direction = "<")
  roc.list[["WeightedAvg"]] <- roc(outcome, y.pred.avg.pool.all.W, ci = TRUE, direction = "<")

  # ----- Build tables -----

  roc.tab <- sapply(roc.list, function(x) fix_auc(as.numeric(x$auc)))

  ci.low.tab <- sapply(roc.list, function(x) fix_auc(x$ci[1]))
  ci.up.tab  <- sapply(roc.list, function(x) fix_auc(x$ci[3]))
  sp.98.tab <- sapply(roc.list, function(x) as.numeric(coords(x, x = 0.98, input = "specificity", ret = "sensitivity")))
  sp.95.tab <- sapply(roc.list, function(x) as.numeric(coords(x, x = 0.95, input = "specificity", ret = "sensitivity")))
  pred_list <- c(y.pred.avg.pool, y.ens.avg.pool, list(SimpleAvg = y.pred.avg.pool1.all.S, WeightedAvg = y.pred.avg.pool.all.W))
  names(pred_list) <- names(roc.list)
  sp95_ci <- t(vapply(
    names(roc.list),
    function(model_name) binary_sp_ci(outcome, pred_list[[model_name]]),
    numeric(2)
  ))

  # ----- Return results -----
  ret_list <- list(
    auc    = roc.tab,
    ci.low = ci.low.tab,
    ci.up  = ci.up.tab,
    sp.98  = sp.98.tab,
    sp.95  = sp.95.tab,
    sp.95.ci.low = sp95_ci[, "ci.low"],
    sp.95.ci.up = sp95_ci[, "ci.up"],
    roc = roc.list,
    predictions = pred_list
  )

  return(ret_list)
}
survival_roc <- function(outcome, survdays, time_point, result, models.indiv, iter, ensemble = FALSE, models.ens = NULL) {

  # Helper to compute mean predictions
  get_avg_pred <- function(model_name, source) {
    apply(matrix(source[[model_name]], length(outcome), iter), 1, mean, na.rm = TRUE)
  }

  # Compute average predictions for each model
  preds <- list()
  for (m in models.indiv) {
    if (m %in% names(result$individual)) {
      preds[[m]] <- get_avg_pred(m, result$individual)
    }
  }

  # Compute average predictions for ensemble models
  if (ensemble) {
    for (m in models.ens) {
      if (m %in% names(result$ensemble)) {
        preds[[paste0(m, ".ens")]] <- get_avg_pred(m, result$ensemble)
      }
    }
  }

  # Simple and weighted averages
  preds$WeightedAvg <- apply(matrix(result$y.pred.avg.pool.all, length(outcome), iter), 1, mean, na.rm = TRUE)
  preds$SimpleAvg   <- apply(matrix(result$y.pred.avg.pool1.all, length(outcome), iter), 1, mean, na.rm = TRUE)

  # Compute ROC for each prediction
  rocs <- lapply(preds, function(marker) {
    survivalROC(Stime = survdays, status = outcome, marker = marker,
                predict.time = time_point, method = "KM")
  })

  # Extract AUC, sensitivity at 98% and 95% specificity
  auc_tab <- sapply(rocs, function(x) fix_auc(x$AUC))
  sp98_tab  <- sapply(rocs, function(x) x$TP[which.min(abs(1 - x$FP - 0.98))])
  sp95_tab  <- sapply(rocs, function(x) x$TP[which.min(abs(1 - x$FP - 0.95))])
  cindex_stats <- t(vapply(
    preds,
    function(marker) survival_cindex(survdays, outcome, marker),
    numeric(3)
  ))
  sp95_ci <- t(vapply(
    preds,
    function(marker) survival_sp_ci(survdays, outcome, marker, time_point),
    numeric(2)
  ))
  auc_ci <- t(vapply(
    preds,
    function(marker) survival_auc_ci(survdays, outcome, marker, time_point),
    numeric(2)
  ))

  return(list(
    auc    = auc_tab,
    ci.low = auc_ci[, "ci.low"],
    ci.up  = auc_ci[, "ci.up"],
    sp.98  = sp98_tab,
    sp.95  = sp95_tab,
    sp.95.ci.low = sp95_ci[, "ci.low"],
    sp.95.ci.up = sp95_ci[, "ci.up"],
    cindex = cindex_stats[, "cindex"],
    cindex.ci.low = cindex_stats[, "ci.low"],
    cindex.ci.up = cindex_stats[, "ci.up"],
    roc = rocs,
    predictions = preds
  ))
}
