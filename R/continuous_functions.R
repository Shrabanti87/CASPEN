# Internal functions for continuous/regression outcomes

continuous_folds <- function(y, k.folds) {
  k.folds <- max(2L, min(as.integer(k.folds %||% 3L), length(y)))
  sample(rep(seq_len(k.folds), length.out = length(y)))
}

continuous_metric_values <- function(y, pred) {
  y <- as.numeric(y)
  pred <- as.numeric(pred)
  ok <- is.finite(y) & is.finite(pred)
  if (sum(ok) < 2) {
    return(c(r2 = NA_real_, rmse = NA_real_, mae = NA_real_,
             cindex = NA_real_))
  }
  sst <- sum((y[ok] - mean(y[ok]))^2)
  sse <- sum((y[ok] - pred[ok])^2)
  r2 <- if (is.finite(sst) && sst > 0) 1 - sse / sst else NA_real_
  rmse <- sqrt(mean((y[ok] - pred[ok])^2))
  mae <- mean(abs(y[ok] - pred[ok]))
  cindex <- continuous_cindex(y[ok], pred[ok])
  c(r2 = r2, rmse = rmse, mae = mae, cindex = cindex)
}

continuous_cindex <- function(y, pred) {
  y <- as.numeric(y)
  pred <- as.numeric(pred)
  ok <- is.finite(y) & is.finite(pred)
  y <- y[ok]
  pred <- pred[ok]
  n <- length(y)
  if (n < 2) return(NA_real_)
  total <- 0
  conc <- 0
  for (i in seq_len(n - 1L)) {
    dy <- y[i] - y[(i + 1L):n]
    dp <- pred[i] - pred[(i + 1L):n]
    keep <- dy != 0
    if (!any(keep)) next
    total <- total + sum(keep)
    prod <- dy[keep] * dp[keep]
    conc <- conc + sum(prod > 0) + 0.5 * sum(prod == 0)
  }
  if (total == 0) NA_real_ else conc / total
}

continuous_scale_data <- function(xtrain, xtest = NULL) {
  xtrain <- as.data.frame(xtrain, check.names = FALSE)
  m <- vapply(xtrain, mean, numeric(1), na.rm = TRUE)
  s <- vapply(xtrain, stats::sd, numeric(1), na.rm = TRUE)
  s[!is.finite(s) | s == 0] <- 1
  xtrain.sd <- as.data.frame(scale(xtrain, center = m, scale = s),
                             check.names = FALSE)
  if (is.null(xtest)) return(list(train = xtrain.sd, test = NULL))
  xtest <- as.data.frame(xtest, check.names = FALSE)
  xtest.sd <- as.data.frame(scale(xtest, center = m, scale = s),
                            check.names = FALSE)
  list(train = xtrain.sd, test = xtest.sd)
}

continuous_fit_predict <- function(model_name, xtrain, ytrain, xtest,
                                   param = NULL, defaults = NULL) {
  model_name <- toupper(model_name)
  param <- modifyList(defaults[[model_name]] %||% list(), param %||% list())
  xtrain <- as.data.frame(xtrain, check.names = FALSE)
  xtest <- as.data.frame(xtest, check.names = FALSE)
  ytrain <- as.numeric(ytrain)

  if (!ncol(xtrain)) {
    return(list(train = rep(mean(ytrain, na.rm = TRUE), length(ytrain)),
                test = rep(mean(ytrain, na.rm = TRUE), nrow(xtest))))
  }

  if (model_name == "XG") {
    nrounds <- param$nrounds %||% 100
    param$nrounds <- NULL
    param$objective <- param$objective %||% "reg:squarederror"
    param$eval_metric <- param$eval_metric %||% "rmse"
    param$verbosity <- param$verbosity %||% 0
    if (!is.null(param$eta)) {
      param$learning_rate <- param$eta
      param$eta <- NULL
    }
    dtrain <- xgboost::xgb.DMatrix(as.matrix(xtrain), label = ytrain)
    dtest <- xgboost::xgb.DMatrix(as.matrix(xtest))
    fit <- quiet(xgboost::xgb.train(params = param, data = dtrain,
                                    nrounds = nrounds))
    return(list(train = as.numeric(predict(fit, dtrain)),
                test = as.numeric(predict(fit, dtest))))
  }

  if (model_name == "RF") {
    fit <- quiet(do.call(randomForest::randomForest,
                         c(list(x = xtrain, y = ytrain), param)))
    return(list(train = as.numeric(predict(fit, newdata = xtrain)),
                test = as.numeric(predict(fit, newdata = xtest))))
  }

  if (model_name == "EN") {
    s.val <- param$s %||% "lambda.1se"
    param$s <- NULL
    param$family <- param$family %||% "gaussian"
    fit <- quiet(do.call(glmnet::cv.glmnet,
                         c(list(x = as.matrix(xtrain), y = ytrain), param)))
    return(list(train = as.numeric(predict(fit, newx = as.matrix(xtrain),
                                           s = s.val)),
                test = as.numeric(predict(fit, newx = as.matrix(xtest),
                                          s = s.val))))
  }

  if (model_name == "GB" || model_name == "ADB") {
    if (model_name == "ADB") {
      param$interaction.depth <- param$interaction.depth %||%
        param$tree_depth %||% 2
      param$tree_depth <- NULL
    }
    param$distribution <- param$distribution %||% "gaussian"
    param$n.trees <- param$n.trees %||% 200
    fit <- quiet(do.call(gbm::gbm,
                         c(list(formula = y ~ ., data = data.frame(y = ytrain,
                                                                    xtrain)),
                           param)))
    return(list(train = as.numeric(predict(fit, data.frame(xtrain),
                                           n.trees = fit$n.trees)),
                test = as.numeric(predict(fit, data.frame(xtest),
                                          n.trees = fit$n.trees))))
  }

  if (model_name == "KNN") {
    k.grid <- param$tuneGrid %||% expand.grid(k = seq(3, min(25, max(3, nrow(xtrain) - 1)), by = 2))
    param$tuneGrid <- NULL
    ctrl <- param$trControl %||% caret::trainControl(method = "cv", number = 3)
    param$trControl <- NULL
    fit <- quiet(do.call(caret::train,
                         c(list(x = xtrain, y = ytrain, method = "knn",
                                tuneGrid = k.grid, trControl = ctrl,
                                preProcess = c("center", "scale")),
                           param)))
    return(list(train = as.numeric(predict(fit, newdata = xtrain)),
                test = as.numeric(predict(fit, newdata = xtest))))
  }

  if (model_name == "SVM") {
    fit <- quiet(do.call(e1071::svm,
                         c(list(x = xtrain, y = ytrain, type = "eps-regression"),
                           param)))
    return(list(train = as.numeric(predict(fit, newdata = xtrain)),
                test = as.numeric(predict(fit, newdata = xtest))))
  }

  if (model_name == "DCT") {
    fit <- quiet(do.call(rpart::rpart,
                         c(list(formula = y ~ ., data = data.frame(y = ytrain,
                                                                    xtrain),
                                method = "anova"),
                           param)))
    return(list(train = as.numeric(predict(fit, newdata = data.frame(xtrain))),
                test = as.numeric(predict(fit, newdata = data.frame(xtest)))))
  }

  if (model_name == "NN") {
    scaled <- continuous_scale_data(xtrain, xtest)
    param$linear.output <- param$linear.output %||% TRUE
    param$hidden <- param$hidden %||% c(8, 4)
    param$rep <- param$rep %||% 1
    fit <- quiet(do.call(neuralnet::neuralnet,
                         c(list(formula = y ~ .,
                                data = data.frame(y = ytrain, scaled$train)),
                           param)))
    return(list(train = as.numeric(predict(fit, newdata = scaled$train)),
                test = as.numeric(predict(fit, newdata = scaled$test))))
  }

  if (model_name == "MB") {
    param$mstop <- param$mstop %||% 300
    fit <- quiet(do.call(mboost::glmboost,
                         c(list(formula = y ~ .,
                                data = data.frame(y = ytrain, xtrain),
                                family = mboost::Gaussian()),
                           param)))
    return(list(train = as.numeric(predict(fit, newdata = data.frame(xtrain))),
                test = as.numeric(predict(fit, newdata = data.frame(xtest)))))
  }

  if (model_name == "NB") {
    warning("NB is classification-only and is not fit for continuous outcomes.")
    return(list(train = rep(NA_real_, length(ytrain)),
                test = rep(NA_real_, nrow(xtest))))
  }

  warning("Unknown continuous model: ", model_name)
  list(train = rep(NA_real_, length(ytrain)),
       test = rep(NA_real_, nrow(xtest)))
}

continuous_defaults <- function() {
  list(
    XG = list(max_depth = 2, eta = 0.05, nrounds = 100,
              objective = "reg:squarederror", eval_metric = "rmse",
              nthread = 2, verbosity = 0),
    RF = list(ntree = 500),
    EN = list(alpha = 0.2, family = "gaussian", s = "lambda.1se"),
    ADB = list(n.trees = 200, interaction.depth = 2, shrinkage = 0.05,
               n.minobsinnode = 2, distribution = "gaussian"),
    MB = list(mstop = 300),
    GB = list(n.trees = 300, interaction.depth = 3, shrinkage = 0.05,
              n.minobsinnode = 2, distribution = "gaussian", n.cores = 1),
    KNN = list(),
    SVM = list(kernel = "radial"),
    NN = list(hidden = c(8, 4), rep = 1, linear.output = TRUE),
    NB = list(),
    DCT = list(cp = 0.001)
  )
}

continuous_model <- function(outcome, x, models.indiv, iter, num.folds = NULL,
                             ensemble = FALSE, models.ens = NULL,
                             param.indiv = NULL, param.ens = NULL,
                             iter.offset = 0) {
  outcome <- as.numeric(outcome)
  x <- as.data.frame(x, check.names = FALSE)
  param.indiv <- param.indiv %||% list()
  param.ens <- param.ens %||% list()
  models.indiv <- toupper(models.indiv)
  models.ens <- if (isTRUE(ensemble)) toupper(models.ens %||% models.indiv) else character(0)
  nfold <- max(2L, as.integer(num.folds %||% 3L))
  defaults <- continuous_defaults()

  pred_pools <- stats::setNames(vector("list", length(models.indiv)), models.indiv)
  ens_pools <- if (ensemble) stats::setNames(vector("list", length(models.ens)), models.ens) else NULL
  y.pred.avg.pool.all <- NULL
  y.pred.avg.pool1.all <- NULL

  for (j in iter.offset + seq_len(iter)) {
    message("Iteration: ", j)
    set.seed(j)
    cv.fold <- continuous_folds(outcome, nfold)
    y.pred <- stats::setNames(lapply(models.indiv, function(z) rep(NA_real_, length(outcome))), models.indiv)
    y.train <- stats::setNames(lapply(models.indiv, function(z) rep(NA_real_, length(outcome))), models.indiv)
    y.ens <- if (ensemble) stats::setNames(lapply(models.ens, function(z) rep(NA_real_, length(outcome))), models.ens) else list()
    y.pred.avg.all <- rep(NA_real_, length(outcome))

    for (i in seq_len(nfold)) {
      cur.train <- which(cv.fold != i)
      cur.test <- which(cv.fold == i)
      xtrain <- x[cur.train, , drop = FALSE]
      xtest <- x[cur.test, , drop = FALSE]

      for (model_name in models.indiv) {
        message(model_name, " regression running...")
        pred <- continuous_fit_predict(
          model_name, xtrain, outcome[cur.train], xtest,
          param = param.indiv[[model_name]], defaults = defaults
        )
        y.train[[model_name]][cur.train] <- pred$train
        y.pred[[model_name]][cur.test] <- pred$test
      }

      xtrain.ens <- as.data.frame(lapply(y.train[models.indiv], `[`, cur.train))
      xtest.ens <- as.data.frame(lapply(y.pred[models.indiv], `[`, cur.test))
      colnames(xtrain.ens) <- make.names(models.indiv, unique = TRUE)
      colnames(xtest.ens) <- colnames(xtrain.ens)
      keep.cols <- vapply(xtrain.ens, function(z) any(is.finite(z)), logical(1))
      xtrain.ens <- xtrain.ens[, keep.cols, drop = FALSE]
      xtest.ens <- xtest.ens[, keep.cols, drop = FALSE]
      for (cc in colnames(xtrain.ens)) {
        fill <- mean(xtrain.ens[[cc]], na.rm = TRUE)
        if (!is.finite(fill)) fill <- 0
        xtrain.ens[[cc]][!is.finite(xtrain.ens[[cc]])] <- fill
        xtest.ens[[cc]][!is.finite(xtest.ens[[cc]])] <- fill
      }

      if (ensemble && ncol(xtrain.ens)) {
        for (model_name in models.ens) {
          message(model_name, " regression ensemble...")
          pred <- continuous_fit_predict(
            model_name, xtrain.ens, outcome[cur.train], xtest.ens,
            param = param.ens[[model_name]], defaults = defaults
          )
          y.ens[[model_name]][cur.test] <- pred$test
        }
      }

      if (ncol(xtrain.ens)) {
        fit.avg <- quiet(stats::lm(y ~ ., data = data.frame(y = outcome[cur.train],
                                                            xtrain.ens)))
        y.pred.avg.all[cur.test] <- quiet(as.numeric(predict(fit.avg,
                                                             newdata = xtest.ens)))
      } else {
        y.pred.avg.all[cur.test] <- mean(outcome[cur.train], na.rm = TRUE)
      }
    }

    y.pred.mat <- as.data.frame(y.pred)
    y.pred.avg1.all <- rowMeans(y.pred.mat, na.rm = TRUE)
    for (model in models.indiv) pred_pools[[model]] <- c(pred_pools[[model]], y.pred[[model]])
    if (ensemble) {
      for (model in models.ens) ens_pools[[model]] <- c(ens_pools[[model]], y.ens[[model]])
    }
    y.pred.avg.pool.all <- c(y.pred.avg.pool.all, y.pred.avg.all)
    y.pred.avg.pool1.all <- c(y.pred.avg.pool1.all, y.pred.avg1.all)
  }

  list(
    individual = pred_pools,
    ensemble = if (ensemble) ens_pools else NULL,
    y.pred.avg.pool.all = y.pred.avg.pool.all,
    y.pred.avg.pool1.all = y.pred.avg.pool1.all
  )
}

continuous_metrics <- function(outcome, result, models.indiv, iter,
                               ensemble = FALSE, models.ens = NULL) {
  outcome <- as.numeric(outcome)
  get_avg <- function(source, model_name) {
    vals <- source[[model_name]]
    if (is.null(vals)) return(rep(NA_real_, length(outcome)))
    rowMeans(matrix(vals, nrow = length(outcome), ncol = iter), na.rm = TRUE)
  }
  preds <- list()
  for (m in intersect(toupper(models.indiv), names(result$individual))) {
    preds[[m]] <- get_avg(result$individual, m)
  }
  if (ensemble && !is.null(result$ensemble)) {
    for (m in intersect(toupper(models.ens), names(result$ensemble))) {
      preds[[paste0(m, ".ens")]] <- get_avg(result$ensemble, m)
    }
  }
  preds$SimpleAvg <- rowMeans(matrix(result$y.pred.avg.pool1.all,
                                     nrow = length(outcome), ncol = iter),
                              na.rm = TRUE)
  preds$WeightedAvg <- rowMeans(matrix(result$y.pred.avg.pool.all,
                                       nrow = length(outcome), ncol = iter),
                                na.rm = TRUE)

  metric.mat <- t(vapply(preds, function(pred) continuous_metric_values(outcome, pred),
                         numeric(4)))
  r2 <- metric.mat[, "r2"]
  rmse <- metric.mat[, "rmse"]
  mae <- metric.mat[, "mae"]
  cindex <- metric.mat[, "cindex"]
  names(r2) <- names(rmse) <- names(mae) <- names(cindex) <- names(preds)
  na_vec <- stats::setNames(rep(NA_real_, length(preds)), names(preds))
  list(
    auc = r2,
    ci.low = na_vec,
    ci.up = na_vec,
    sp.98 = na_vec,
    sp.95 = na_vec,
    sp.95.ci.low = na_vec,
    sp.95.ci.up = na_vec,
    r2 = r2,
    rmse = rmse,
    mae = mae,
    cindex = cindex,
    predictions = preds
  )
}

continuous_train_performance <- function(outcome, train_data, features,
                                         models.indiv, iter, num.folds = NULL,
                                         ensemble = FALSE, models.ens = NULL,
                                         param.indiv = NULL, param.ens = NULL,
                                         iter.offset = 0) {
  genes <- unique(unlist(lapply(features, extract_pathway_genes), use.names = FALSE))
  genes <- intersect(genes, colnames(train_data))
  if (!length(genes)) stop("No selected continuous-outcome features found in train_data.")
  res <- continuous_model(
    outcome = outcome,
    x = train_data[, genes, drop = FALSE],
    models.indiv = models.indiv,
    iter = iter,
    num.folds = num.folds,
    ensemble = ensemble,
    models.ens = models.ens,
    param.indiv = param.indiv,
    param.ens = param.ens,
    iter.offset = iter.offset
  )
  met <- continuous_metrics(outcome, res, models.indiv, iter, ensemble, models.ens)
  list(
    R2 = met$r2,
    r2 = met$r2,
    RMSE = met$rmse,
    rmse = met$rmse,
    MAE = met$mae,
    mae = met$mae,
    C.index = met$cindex,
    cindex = met$cindex,
    auc = met$auc,
    SP.95 = met$sp.95,
    predictions = met$predictions
  )
}

continuous_test_performance <- function(outcome, train_data, test_data,
                                        test_outcome, features, models.indiv,
                                        ensemble = FALSE, models.ens = NULL,
                                        param.indiv = NULL, param.ens = NULL) {
  outcome <- as.numeric(outcome)
  test_outcome <- as.numeric(test_outcome)
  models.indiv <- toupper(models.indiv)
  models.ens <- if (isTRUE(ensemble)) toupper(models.ens %||% models.indiv) else character(0)
  genes <- unique(unlist(lapply(features, extract_pathway_genes), use.names = FALSE))
  genes <- intersect(genes, intersect(colnames(train_data), colnames(test_data)))
  if (!length(genes)) stop("No shared selected continuous-outcome features found.")
  defaults <- continuous_defaults()
  xtrain <- as.data.frame(train_data[, genes, drop = FALSE], check.names = FALSE)
  xtest <- as.data.frame(test_data[, genes, drop = FALSE], check.names = FALSE)

  train.preds <- test.preds <- list()
  for (model_name in models.indiv) {
    pred <- continuous_fit_predict(model_name, xtrain, outcome, xtest,
                                   param = param.indiv[[model_name]],
                                   defaults = defaults)
    train.preds[[model_name]] <- pred$train
    test.preds[[model_name]] <- pred$test
  }

  xtrain.ens <- as.data.frame(train.preds)
  xtest.ens <- as.data.frame(test.preds)
  keep.cols <- vapply(xtrain.ens, function(z) any(is.finite(z)), logical(1))
  xtrain.ens <- xtrain.ens[, keep.cols, drop = FALSE]
  xtest.ens <- xtest.ens[, keep.cols, drop = FALSE]
  for (cc in colnames(xtrain.ens)) {
    fill <- mean(xtrain.ens[[cc]], na.rm = TRUE)
    if (!is.finite(fill)) fill <- 0
    xtrain.ens[[cc]][!is.finite(xtrain.ens[[cc]])] <- fill
    xtest.ens[[cc]][!is.finite(xtest.ens[[cc]])] <- fill
  }

  if (ensemble && ncol(xtrain.ens)) {
    for (model_name in models.ens) {
      pred <- continuous_fit_predict(model_name, xtrain.ens, outcome, xtest.ens,
                                     param = param.ens[[model_name]],
                                     defaults = defaults)
      test.preds[[paste0(model_name, ".ens")]] <- pred$test
    }
  }
  test.preds$SimpleAvg <- rowMeans(as.data.frame(test.preds[models.indiv]),
                                   na.rm = TRUE)
  if (ncol(xtrain.ens)) {
    fit.avg <- quiet(stats::lm(y ~ ., data = data.frame(y = outcome, xtrain.ens)))
    test.preds$WeightedAvg <- quiet(as.numeric(predict(fit.avg, newdata = xtest.ens)))
  } else {
    test.preds$WeightedAvg <- rep(mean(outcome, na.rm = TRUE), nrow(xtest))
  }

  metric.mat <- t(vapply(test.preds, function(pred) continuous_metric_values(test_outcome, pred),
                         numeric(4)))
  r2 <- metric.mat[, "r2"]
  rmse <- metric.mat[, "rmse"]
  mae <- metric.mat[, "mae"]
  cindex <- metric.mat[, "cindex"]
  list(
    R2 = r2,
    r2 = r2,
    RMSE = rmse,
    rmse = rmse,
    MAE = mae,
    mae = mae,
    C.index = cindex,
    cindex = cindex,
    auc = r2,
    SP.95 = stats::setNames(rep(NA_real_, length(r2)), names(r2)),
    predictions = test.preds
  )
}
