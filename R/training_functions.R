# Internal functions: not exported
"%!in%" <- Negate("%in%")
`%||%` <- function(a, b) if (!is.null(a)) a else b
extract_coords <- function(multi_roc, spec_val) {
  sapply(names(multi_roc$rocs), function(model_name) {
    r <- multi_roc$rocs[[model_name]]
    vals <- sapply(r, function(single_class_roc) {
      coords(single_class_roc,
             x = spec_val,
             input = "specificity",
             ret = "sensitivity",
             transpose = FALSE,
             as.matrix = FALSE,
             best.method = "closest") |> as.numeric()
    })
    mean(vals, na.rm = TRUE)
  })
}
quiet <- function(expr) {
  suppressWarnings(suppressMessages(expr))
}


categorical_model <- function(outcome, x, models.indiv, iter, num.folds = NULL, ensemble = FALSE, models.ens = NULL, param.indiv = NULL, param.ens = NULL, iter.offset = 0) {
  nfold <- if (!is.null(num.folds)) num.folds else 3
  indiv.defaults <- list(
    XG = list(learning_rate = 0.1, objective = "multi:softprob", nrounds = 1000, verbosity = 0),
    RF = list(ntree = 500),
    EN = list(nfolds = 3, s = "lambda.1se"),
    ADB = list(tree_depth = 3, n_rounds = 200, verbose = FALSE),
    MB = list(mstop = 300),
    GB = list(n.trees = 50, interaction.depth = 6, shrinkage = 0.1),
    KNN = list(trControl = trainControl(method = "cv", number = 3, classProbs = TRUE)),
    SVM = list(kernel = "radial", cost = 50, gamma = 0.75),
    NB = list(laplace = 10, threshold = 0.001),
    DCT = list(cp = 0.001),
    NN = list(hidden = c(8,4), rep = 1)
  )

  ens.defaults <- list(
    XG = list(learning_rate = 0.1, objective = "multi:softprob", nrounds = 50, verbosity = 0),
    RF = list(ntree = 500),
    EN = list(nfolds = 3, s = "lambda.1se"),
    ADB = list(tree_depth = 3, n_rounds = 200, verbose = FALSE),
    MB = list(mstop = 300),
    GB = list(n.trees = 50, interaction.depth = 6, shrinkage = 0.1),
    KNN = list(trControl = trainControl(method = "cv", number = 3, classProbs = TRUE)),
    SVM = list(kernel = "radial", cost = 50, gamma = 0.75),
    NB = list(laplace = 10, threshold = 0.001),
    DCT = list(cp = 0.001),
    NN = list(hidden = c(8,4), rep = 1)
  )

  # Initialize pooled results as lists
  pred_pools <- list()
  for (model in models.indiv) pred_pools[[model]] <- NULL

  ens_pools <- list()
  if (ensemble) {
    for (model in models.ens) ens_pools[[model]] <- NULL
  }

  y.pred.avg.pool.all <- NULL
  y.pred.avg.pool1.all <- NULL

  for (j in iter.offset + seq_len(iter)) {
    message(j)
    set.seed(j)

    cv.fold <-balanced.folds.vec(outcome,nfold)

    # Initialize predictions for each iteration
    y.pred <- list()
    y.train <- list()
    y.ens <- list()

    for (m in models.indiv) {
      y.pred[[m]] <- matrix(NA, nrow = length(outcome), ncol = length(unique(outcome)))
      y.train[[m]] <- matrix(NA, nrow = length(outcome), ncol = length(unique(outcome)))
    }
    if (ensemble) {
      for (m in models.ens) y.ens[[m]] <- matrix(NA, nrow = length(outcome), ncol = length(unique(outcome)))
    }

    y.pred.avg.all <- matrix(NA, nrow = length(outcome), ncol = length(unique(outcome)))

    # Cross-validation loop

    for (i in 1:nfold) {
      cur.train <- which(cv.fold != i)
      cur.test <- which(cv.fold == i)

      xtrain <- as.matrix(x[cur.train, ])
      xtest  <- as.matrix(x[cur.test, ])

      xtrain_cl <- as.data.frame(x[cur.train, ])
      xtrain_cl$cluster <- factor(paste0("C", outcome[cur.train]))

      xtest_cl <- as.data.frame(x[cur.test, ])
      xtest_cl$cluster <- factor(paste0("C", outcome[cur.test]))

      ytrain <- outcome[cur.train]
      ytest <- outcome[cur.test]
      num_classes <- length(unique(outcome))

      # Run individual models
      for (model_name in models.indiv) {
        if (toupper(model_name) == "XG") {
          message("XgBoost running...")
          params_user <- modifyList(indiv.defaults$XG,param.indiv[["XG"]] %||% list())

          temp <- quiet(do.call(
            xgboost::xgboost,
            c(
              list(
                x = as.matrix(xtrain),
                y = as.factor(ytrain)
              ),
              params_user
            )
          ))
          pred_probs <- quiet(predict(temp, newdata = as.matrix(xtest)))
          pred_matrix <- matrix(pred_probs, nrow = nrow(xtest), ncol = num_classes)
          y.pred[[model_name]][cur.test, ] <- pred_matrix
          pred_probs.t <- quiet(predict(temp, newdata = as.matrix(xtrain)))
          pred_matrix.t <- matrix(pred_probs.t, nrow = nrow(xtrain), ncol = num_classes)
          y.train[[model_name]][cur.train, 1:num_classes] <- pred_matrix.t

        }
        if (toupper(model_name) == "RF") {
          message("RandomForest running...")
          params <- modifyList(indiv.defaults$RF, param.indiv[["RF"]] %||% list())
          temp <- quiet(do.call(randomForest, c(list(x = xtrain, y = as.factor(ytrain)), params)))
          y.pred[[model_name]][cur.test,] <- quiet(predict(temp, newdata = xtest, type = "prob"))
          y.train[[model_name]][cur.train,1:length(unique(outcome))] <- quiet(predict(temp, newdata = xtrain, type = "prob"))
        }
        if (toupper(model_name) == "EN") {
          message("ElasticNet running...")
          params <- modifyList(indiv.defaults$EN, param.indiv[["EN"]] %||% list())
          train.params <- params[setdiff(names(params), "s")]
          temp <- quiet(do.call(cv.glmnet, c(list(xtrain, as.factor(ytrain), family = "multinomial"), train.params)))
          s.val <- param.indiv[["EN"]]$s %||% "lambda.1se"
          y.pred[[model_name]][cur.test,] <- quiet(predict(temp, newx = xtest, s = s.val, type = "response"))
          y.train[[model_name]][cur.train,] <- quiet(predict(temp, newx = xtrain, s = s.val, type = "response"))
        }
        if (toupper(model_name) == "ADB") {
          message("AdaBoost one-vs-rest running...")
          params <- modifyList(indiv.defaults$ADB, param.indiv[["ADB"]] %||% list())
          y.pred[[model_name]][cur.test, ] <- quiet(do.call(
            multiclass_adb_predict,
            c(list(xtrain = xtrain, outcome = ytrain, xtest = xtest,
                   classes = sort(unique(outcome))), params)
          ))
          y.train[[model_name]][cur.train, ] <- quiet(do.call(
            multiclass_adb_predict,
            c(list(xtrain = xtrain, outcome = ytrain,
                   classes = sort(unique(outcome))), params)
          ))
        }
        if (toupper(model_name) == "MB") {
          message("MBoost one-vs-rest running...")
          params <- modifyList(indiv.defaults$MB, param.indiv[["MB"]] %||% list())
          y.pred[[model_name]][cur.test, ] <- quiet(do.call(
            multiclass_mb_predict,
            c(list(xtrain = xtrain, outcome = ytrain, xtest = xtest,
                   classes = sort(unique(outcome))), params)
          ))
          y.train[[model_name]][cur.train, ] <- quiet(do.call(
            multiclass_mb_predict,
            c(list(xtrain = xtrain, outcome = ytrain,
                   classes = sort(unique(outcome))), params)
          ))
        }
        if (toupper(model_name) == "GB") {
          message("GradientBoost running...")
          params <- modifyList(indiv.defaults$GB, param.indiv[["GB"]] %||% list())
          temp <- quiet(do.call(gbm.fit, c(list(xtrain, as.factor(ytrain), distribution = "multinomial"), params)))
          y.pred[[model_name]][cur.test,] <- quiet(predict.gbm(temp, data.frame(xtest), type = "response"))
          y.train[[model_name]][cur.train,] <- quiet(predict.gbm(temp, data.frame(xtrain), type = "response"))
        }
        if (toupper(model_name) == "KNN") {
          message("KNN running...")
          trc <- param.indiv[["KNN"]]$trControl %||% trainControl(method = "cv", number = 3, classProbs = TRUE)
          params <- modifyList(indiv.defaults$KNN, param.indiv[["KNN"]] %||% list())
          params$trControl <- NULL  # remove so not passed twice

          temp <- quiet(do.call(caret::train, c(
            list(
              cluster ~ .,
              data = xtrain_cl,
              method = "knn",
              trControl = trc,
              preProcess = c("center", "scale")
            ),
            params
          )))
          y.pred[[model_name]][cur.test, 1:length(unique(outcome))] <- quiet(as.matrix(predict(temp, newdata = xtest_cl, type = "prob")))
          y.train[[model_name]][cur.train, 1:length(unique(outcome))] <- quiet(as.matrix(predict(temp, newdata = xtrain_cl, type = "prob")))
        }
        if (toupper(model_name) == "SVM") {
          message("SVM running...")
          params <- modifyList(indiv.defaults$SVM, param.indiv[["SVM"]] %||% list())
          temp <- quiet(do.call(svm, c(list(cluster ~ ., data = xtrain_cl, probability = TRUE), params)))

          pred <- quiet(predict(temp, newdata= xtest_cl, probability = TRUE))
          probs <- attr(pred, "probabilities")
          class_order <- paste0("C", seq_len(length(unique(outcome))))
          probs_reordered <- probs[, class_order, drop = FALSE]
          y.pred[[model_name]][cur.test, ] <- probs_reordered

          pred <- quiet(predict(temp, newdata= xtrain_cl, probability = TRUE))
          probs <- attr(pred, "probabilities")
          probs_reordered <- probs[, class_order, drop = FALSE]
          y.train[[model_name]][cur.train, ] <- probs_reordered
        }
        if (toupper(model_name) == "NB") {
          message("NB running...")
          params <- modifyList(indiv.defaults$NB, param.indiv[["NB"]] %||% list())
          temp <- quiet(do.call(naiveBayes, c(list(cluster ~., data = xtrain_cl), params)))
          y.pred[[model_name]][cur.test, 1:length(unique(outcome))] <- quiet(predict(temp, xtest_cl, type = "raw"))
          y.train[[model_name]][cur.train, 1:length(unique(outcome))] <- quiet(predict(temp, xtrain_cl, type = "raw"))
        }
        if (toupper(model_name) == "DCT") {
          message("DecisionTree running...")
          params <- modifyList(indiv.defaults$DCT, param.indiv[["DCT"]] %||% list())
          temp <- quiet(do.call(rpart, c(list(cluster ~., data = xtrain_cl, method = "class"), params)))
          y.pred[[model_name]][cur.test, 1:length(unique(outcome))] <- quiet(as.matrix(predict(temp, newdata= as.data.frame(xtest_cl), type="prob")))
          y.train[[model_name]][cur.train, 1:length(unique(outcome))] <- quiet(as.matrix(predict(temp, newdata= as.data.frame(xtrain_cl), type="prob")))
        }
        if (toupper(model_name) == "NN") {
          message("NeuralNet running...")
          m <- colMeans(as.data.frame(xtrain))
          s <- apply(as.data.frame(xtrain), 2, sd)
          s[s == 0] <- 1
          xtrain.sd <- scale(xtrain, center = m, scale = s)
          xtest.sd  <- scale(xtest,  center = m, scale = s)

          y.mat  <- model.matrix(~ ytrain - 1)  # one-hot encoding
          colnames(y.mat) <- paste0("class", 1:length(unique(outcome)))
          train.df <- data.frame(xtrain.sd, y.mat)

          f <- as.formula(paste(paste(colnames(y.mat), collapse = " + "), "~ ."))
          params <- modifyList(indiv.defaults$NN, param.indiv[["NN"]] %||% list())
          temp.nn <- quiet(do.call(neuralnet, c(list(f, data = train.df, linear.output = FALSE), params)))

          y.pred[[model_name]][cur.test, 1:length(unique(outcome))] <- quiet(predict(temp.nn, newdata = xtest.sd))
          y.train[[model_name]][cur.train, 1:length(unique(outcome))] <- quiet(predict(temp.nn, newdata = xtrain.sd))
        }
      }

      # Ensemble models
      if (ensemble) {
        message("Ensemble running...")
        models_to_use <- intersect(names(y.train), models.ens)

        xtrain.ens <- do.call(cbind, lapply(y.train[models_to_use], function(mat) mat[cur.train, ]))
        xtest.ens  <- do.call(cbind, lapply(y.pred[models_to_use], function(mat) mat[cur.test, ]))

        xtrain.ens_cl <- as.data.frame(xtrain.ens)
        xtrain.ens_cl$cluster <- factor(paste0("C", outcome[cur.train]))

        xtest.ens_cl <- as.data.frame(xtest.ens)
        xtest.ens_cl$cluster <- factor(paste0("C", outcome[cur.test]))

        for (model_name in models.ens) {

          if (toupper(model_name) == "XG") {
            message("XgBoost running...")
            params <- modifyList(ens.defaults$XG, param.ens[["XG"]] %||% list())

            temp <- quiet(do.call(
              xgboost::xgboost,
              c(
                list(
                  x = as.matrix(xtrain.ens),
                  y = as.factor(ytrain)
                ),
                params
              )
            ))
            pred_probs <- quiet(predict(temp, newdata = as.matrix(xtest.ens)))
            pred_matrix <- matrix(pred_probs, nrow = nrow(xtest.ens), ncol = num_classes)
            y.ens[[model_name]][cur.test, ] <- pred_matrix
          }
          if (toupper(model_name) == "RF") {
            message("RF ensemble...")
            params <- modifyList(ens.defaults$RF, param.ens[["RF"]] %||% list())
            temp <- quiet(do.call(randomForest, c(list(x = xtrain.ens, y = as.factor(ytrain)), params)))
            y.ens[[model_name]][cur.test,] <- quiet(predict(temp, newdata = xtest.ens, type = "prob"))
          }
          if (toupper(model_name) == "EN") {
            message("EN ensemble...")
            params <- modifyList(ens.defaults$EN, param.ens[["EN"]] %||% list())
            train.params <- params[setdiff(names(params), "s")]
            temp <- quiet(do.call(cv.glmnet, c(list(xtrain.ens, as.factor(ytrain), family = "multinomial"), train.params)))
            s.val <- param.ens[["EN"]]$s %||% "lambda.1se"
            y.ens[[model_name]][cur.test,] <- quiet(predict(temp, newx = xtest.ens, s = s.val, type = "response"))
          }
          if (toupper(model_name) == "ADB") {
            message("ADB one-vs-rest ensemble...")
            params <- modifyList(ens.defaults$ADB, param.ens[["ADB"]] %||% list())
            y.ens[[model_name]][cur.test, ] <- quiet(do.call(
              multiclass_adb_predict,
              c(list(xtrain = xtrain.ens, outcome = ytrain, xtest = xtest.ens,
                     classes = sort(unique(outcome))), params)
            ))
          }
          if (toupper(model_name) == "MB") {
            message("MB one-vs-rest ensemble...")
            params <- modifyList(ens.defaults$MB, param.ens[["MB"]] %||% list())
            y.ens[[model_name]][cur.test, ] <- quiet(do.call(
              multiclass_mb_predict,
              c(list(xtrain = xtrain.ens, outcome = ytrain, xtest = xtest.ens,
                     classes = sort(unique(outcome))), params)
            ))
          }
          if (toupper(model_name) == "NN") {
            message("NN ensemble...")
            m <- colMeans(as.data.frame(xtrain.ens))
            s <- apply(as.data.frame(xtrain.ens), 2, sd)
            xtrain.ens.sd <- scale(xtrain.ens, center = m, scale = s)
            xtest.ens.sd  <- scale(xtest.ens,  center = m, scale = s)

            y.mat  <- model.matrix(~ ytrain - 1)  # one-hot encoding
            colnames(y.mat) <- paste0("class", 1:length(unique(outcome)))
            train.df <- data.frame(xtrain.ens.sd, y.mat)

            f <- as.formula(paste(paste(colnames(y.mat), collapse = " + "), "~ ."))
            params <- modifyList(ens.defaults$NN, param.ens[["NN"]] %||% list())
            temp.nn <- quiet(do.call(neuralnet, c(list(f, data = train.df, linear.output = FALSE), params)))
            y.ens[[model_name]][cur.test, 1:length(unique(outcome))] <- quiet(predict(temp.nn, newdata = xtest.ens.sd))
          }
          if (toupper(model_name) == "GB") {
            message("GB ensemble...")

            xtrain.mat <- as.matrix(xtrain.ens)
            xtest.mat <- as.matrix(xtest.ens)

            # 1. Store the expected number of features (columns in the ensemble matrix)
            num_features_ens <- ncol(xtrain.mat)

            feature_names <- paste0("EnsVar", 1:num_features_ens)

            colnames(xtrain.mat) <- feature_names
            colnames(xtest.mat) <- feature_names

            if (ncol(xtest.mat) != num_features_ens) {
              stop("FATAL ERROR: Test data columns do not match expected number in GB ensemble.")
            }

            params <- modifyList(ens.defaults$GB, param.ens[["GB"]] %||% list())
            message(class(ytrain))
            message(length(ytrain))

            args_list <- c(
              list(
                x = xtrain.mat,
                y = ytrain,
                distribution = "multinomial"
              ),
              params
            )
            temp <- quiet(do.call(gbm.fit, args_list))

            n_trees_final <- temp$n.trees

            y.pred.array <- quiet(predict.gbm(
              temp,
              newdata = as.data.frame(xtest.mat),
              n.trees = n_trees_final, # Pass total trees built
              type = "response"
            ))

            # 3. Use the LAST dimension index from the *prediction array* itself.
            #    We use dim(y.pred.array)[3] as the highest valid index.
            #    This is the most robust way to get the true last dimension.
            last_valid_index <- dim(y.pred.array)[3]

            # Slice the 3D array to get the 2D matrix of final probabilities
            y.ens[[model_name]][cur.test,] <- y.pred.array[,,last_valid_index]

          }
          if (toupper(model_name) == "KNN") {
            message("KNN ensemble...")
            trc <- param.ens[["KNN"]]$trControl %||% trainControl(method = "cv", number = 3, classProbs = TRUE)
            params <- modifyList(ens.defaults$KNN, param.ens[["KNN"]] %||% list())
            params$trControl <- NULL  # remove so not passed twice

            temp <- quiet(do.call(caret::train, c(
              list(
                cluster ~ .,
                data = xtrain.ens_cl,
                method = "knn",
                trControl = trc,
                preProcess = c("center", "scale")
              ),
              params
            )))
            y.ens[[model_name]][cur.test, 1:length(unique(outcome))] <- quiet(as.matrix(predict(temp, newdata = xtest.ens_cl, type = "prob")))
          }
          if (toupper(model_name) == "SVM") {
            message("SVM ensemble...")
            params <- modifyList(ens.defaults$SVM, param.ens[["SVM"]] %||% list())
            temp <- quiet(do.call(svm, c(list(cluster ~ ., data = data.frame(xtrain.ens_cl), probability = TRUE), params)))

            pred <- quiet(predict(temp, newdata= xtest.ens_cl, probability = TRUE))
            probs <- attr(pred, "probabilities")
            class_order <- levels(xtrain.ens_cl$cluster)
            y.ens[[model_name]][cur.test, ] <- probs[, class_order, drop = FALSE]
          }
          if (toupper(model_name) == "NB") {
            message("NB ensemble...")
            params <- modifyList(ens.defaults$NB, param.ens[["NB"]] %||% list())
            temp <- quiet(do.call(naiveBayes, c(list(cluster ~., data = xtrain.ens_cl), params)))
            y.ens[[model_name]][cur.test, 1:length(unique(outcome))] <- quiet(predict(temp, xtest.ens_cl, type = "raw"))
          }
          if (toupper(model_name) == "DCT") {
            message("DCT ensemble...")
            params <- modifyList(ens.defaults$DCT, param.ens[["DCT"]] %||% list())
            temp <- quiet(do.call(rpart, c(list(cluster ~., data = xtrain.ens_cl, method = "class"), params)))
            y.ens[[model_name]][cur.test, 1:length(unique(outcome))] <- quiet(as.matrix(predict(temp, newdata= as.data.frame(xtest.ens_cl), type="prob")))
          }
        }
      }

      # Weighted average
      message("Weighted avg running...")
      # Train matrix
      y.pred.train.all <- quiet(do.call(cbind, lapply(models.indiv, function(m) {
        mat <- y.train[[m]][cur.train, ]
        colnames(mat) <- paste0(m, "_", 1:ncol(mat))
        mat
      })))
      y.pred.train.all <- as.data.frame(y.pred.train.all)
      y.pred.train.all$outcome <- factor(outcome[cur.train])

      y.pred.mat.all <- quiet(do.call(cbind, lapply(models.indiv, function(m) {
        mat <- y.pred[[m]][cur.test, ]
        colnames(mat) <- paste0(m, "_", 1:ncol(mat))
        mat
      })))
      y.pred.mat.all <- as.data.frame(y.pred.mat.all)

      # Match column names exactly
      temp.avg.all <- quiet(multinom(outcome ~ ., data = y.pred.train.all, trace = FALSE))
      y.pred.avg.all[cur.test, 1:length(unique(outcome))] <- quiet(predict(temp.avg.all, newdata = y.pred.mat.all, type = "probs"))
    }

    # Simple average
    message("Simple average running...")
    pred_list <- lapply(models.indiv, function(m) y.pred[[m]])
    names(pred_list) <- models.indiv
    if (length(pred_list) > 0) {
      pred_array <- array(
        unlist(pred_list, use.names = FALSE),
        dim = c(length(outcome), length(unique(outcome)), length(pred_list))
      )
      y.pred.avg1.all <- apply(pred_array, c(1, 2), mean)
    }

    # Pooling
    message("Pooling iterations.")
    for (model in models.indiv) pred_pools[[model]] <- c(pred_pools[[model]], y.pred[[model]])
    if (ensemble) {
      for (model in models.ens) ens_pools[[model]] <- c(ens_pools[[model]], y.ens[[model]])
    }
    y.pred.avg.pool.all  <- c(y.pred.avg.pool.all, y.pred.avg.all)
    y.pred.avg.pool1.all <- c(y.pred.avg.pool1.all, y.pred.avg1.all)
  }

  # Return all pooled results
  ret_list <- list(
    individual = pred_pools,
    ensemble   = if (ensemble) ens_pools else NULL,
    y.pred.avg.pool.all  = y.pred.avg.pool.all,
    y.pred.avg.pool1.all = y.pred.avg.pool1.all
  )

  return(ret_list)
}
binary_model <- function(outcome, x, models.indiv, iter, num.folds = NULL, ensemble = FALSE, models.ens = NULL, param.indiv = NULL, param.ens = NULL, iter.offset = 0) {
  param.indiv <- param.indiv %||% list()
  param.ens <- param.ens %||% list()
  models.indiv <- toupper(models.indiv)
  models.ens <- if (isTRUE(ensemble)) toupper(models.ens %||% models.indiv) else character(0)
  nfold <- if (!is.null(num.folds)) num.folds else 3
  indiv.defaults <- list(
    XG = list(max_depth = 2, eta = 0.1, nthread = 2, nrounds = 500, objective = "binary:logistic", verbose = 0),
    RF = list(ntree = 500),
    EN = list(alpha = 0.2, family = "binomial", s = "lambda.1se"),
    ADB = list(tree_depth = 3, n_rounds = 200, verbose = FALSE),
    MB = list(mstop = 300),
    GB = list(cv.folds = 10, shrinkage = 0.01, n.minobsinnode = 2, n.trees = 500, distribution = "bernoulli", n.cores = 1),
    KNN = list(k = 5),
    SVM = list(),
    NN = list(hidden = c(8,4), rep = 1, linear.output = FALSE),
    NB = list(usekernel = TRUE),
    DCT = list(method = "class")
  )

  ens.defaults <- list(
    XG = list(
      booster = "gblinear", eta = 0.05, nrounds = 200,
      objective = "binary:logistic", eval_metric = "logloss",
      lambda = 1, alpha = 0, nthread = 2, verbosity = 0
    ),
    RF = list(ntree = 500),
    EN = list(alpha = 0.2, family = "binomial", s = "lambda.1se"),
    ADB = list(tree_depth = 3, n_rounds = 200, verbose = FALSE),
    MB = list(mstop = 300),
    GB = list(cv.folds = 10, shrinkage = 0.01, n.minobsinnode = 2, n.trees = 500, distribution = "bernoulli", n.cores = 1),
    KNN = list(k = 5),
    SVM = list(),
    NN = list(hidden = c(8,4), rep = 1, linear.output = FALSE),
    NB = list(usekernel = TRUE),
    DCT = list(method = "class")
  )


  # Initialize pooled results as lists
  pred_pools <- list()
  for (model in models.indiv) pred_pools[[model]] <- NULL

  ens_pools <- list()
  if (ensemble) {
    for (model in models.ens) ens_pools[[model]] <- NULL
  }

  y.pred.avg.pool.all <- NULL
  y.pred.avg.pool1.all <- NULL

  for (j in iter.offset + seq_len(iter)) {
    message(paste("Iteration:", j))
    set.seed(j)

    cv.fold <-balanced.folds.vec(outcome,nfold)

    # Initialize predictions for each iteration
    y.pred <- list()
    y.train <- list()
    y.ens <- list()

    for (m in models.indiv) {
      y.pred[[m]] <- rep(NA, length(outcome))
      y.train[[m]] <- rep(NA, length(outcome))
    }
    if (ensemble) {
      for (m in models.ens) y.ens[[m]] <- rep(NA, length(outcome))
    }
    y.pred.avg.all<- rep(NA, length(outcome))

    # Cross-validation loop
    for (i in 1:nfold) {
      cur.train <- which(cv.fold != i)
      cur.test <- which(cv.fold == i)

      xtrain <- x[cur.train, ]
      xtest  <- x[cur.test, ]

      ctrl.ens <- trainControl(
        method = "cv",
        number = num.folds %||% 3,
        classProbs = TRUE,
        summaryFunction = twoClassSummary
      )

      # Run individual models
      for (model_name in models.indiv) {
        if (toupper(model_name) == "XG") {
          message("XgBoost running...")
          params <- modifyList(
            indiv.defaults$XG,
            param.indiv[["XG"]] %||% list()
          )
          params$data  <- NULL
          params$label <- NULL
          if (!is.null(params$eta)) {
            params$learning_rate <- params$eta
            params$eta <- NULL
          }
          params$verbose <- NULL
          params$verbosity <- 0
          params$objective <- "binary:logistic"
          dtrain <- xgboost::xgb.DMatrix(data  = as.matrix(xtrain), label = as.numeric(outcome[cur.train]))
          dtest <- xgboost::xgb.DMatrix(data = as.matrix(xtest))
          temp <- quiet(xgboost::xgb.train(params = params, data = dtrain, nrounds = params$nrounds %||% 500))
          y.pred[[model_name]][cur.test] <- quiet(predict(temp, dtest))
          y.train[[model_name]][cur.train] <-quiet(predict(temp, dtrain))

        }
        if (toupper(model_name) == "RF") {
          message("RandomForest running...")
          params <- modifyList(indiv.defaults$RF, param.indiv[["RF"]] %||% list())
          y <- factor(outcome[cur.train], levels = c(0,1))
          params$x <- as.data.frame(xtrain)
          params$y <- y
          temp <- quiet(do.call(randomForest, params))
          y.pred[[model_name]][cur.test]=quiet(predict(temp, newdata=xtest, type="prob")[,2])
          y.train[[model_name]][cur.train]=quiet(predict(temp, newdata=xtrain, type="prob")[,2])
        }
        if (toupper(model_name) == "EN") {
          message("ElasticNet running...")
          en.params <- modifyList(indiv.defaults$EN, param.indiv[["EN"]] %||% list())
          train.params <- en.params[setdiff(names(en.params), "s")]
          temp <- quiet(do.call(cv.glmnet, c(list(x = as.matrix(xtrain), y = outcome[cur.train]), train.params)))
          s.val <- param.indiv[["EN"]]$s %||% "lambda.1se"   # user-specified or default
          y.pred[[model_name]][cur.test] <- quiet(predict(temp, newx = as.matrix(xtest), type = "response", s = s.val))
          y.train[[model_name]][cur.train] <- quiet(predict(temp, newx = as.matrix(xtrain), type = "response", s = s.val))
        }
        if (toupper(model_name) == "ADB") {
          message("AdaBoost running...")
          params <- modifyList(indiv.defaults$ADB, param.indiv[["ADB"]] %||% list())
          y_factor <- factor(outcome[cur.train], levels = c(0,1))
          ada.sr <- ifelse(outcome == 0, -1, 1)
          params$formula <- NULL
          params$x <- NULL
          params$y <- NULL
          temp <- quiet(do.call(adaboost, c(list(xtrain, ada.sr[cur.train]), params)))
          y.pred[[model_name]][cur.test]=quiet(predict(temp, xtest, type = "prob"))
          y.train[[model_name]][cur.train]=quiet(predict(temp, xtrain, type = "prob"))
        }
        if (toupper(model_name) == "MB") {
          message("MBoost running...")
          params <- modifyList(indiv.defaults$MB, param.indiv[["MB"]] %||% list())
          y.train[[model_name]][cur.train] <- quiet(do.call(
            binary_mb_predict,
            c(list(xtrain = xtrain, outcome = outcome[cur.train]), params)
          ))
          y.pred[[model_name]][cur.test] <- quiet(do.call(
            binary_mb_predict,
            c(list(xtrain = xtrain, outcome = outcome[cur.train], xtest = xtest), params)
          ))
        }
        if (toupper(model_name) == "GB") {
          message("GradientBoost running...")
          params <- modifyList(indiv.defaults$GB, param.indiv[["GB"]] %||% list())
          train.df <- data.frame(y = outcome[cur.train], xtrain)
          params$formula <- y ~.
          params$data <- train.df
          temp <- quiet(do.call(gbm::gbm, params))
          params$x <- NULL
          params$y <- NULL
          y.pred[[model_name]][cur.test]=quiet(predict.gbm(temp, data.frame(xtest), n.trees = temp$n.trees, type = "response"))
          y.train[[model_name]][cur.train]=quiet(predict(temp, data.frame(xtrain), n.trees = temp$n.trees, type = "response"))
        }
        if (toupper(model_name) == "KNN") {
          y.fold <- outcome[cur.train]
          message("KNN running...")

          knn.params <- modifyList(indiv.defaults$KNN, param.ens[["KNN"]] %||% list())

          y.train.factor <- factor(
            ifelse(outcome[cur.train] == 1, "Class1", "Class0"),
            levels = c("Class0","Class1")
          )

          # Grid of k values to try
          knn.grid <- expand.grid(k = seq(3, 25, by = 2))

          # Cross-validation control
          ctrl <- trainControl(method = "cv", number = 5)

          temp <- quiet(caret::train(
            x = xtrain,
            y = y.train.factor,
            method = "knn",
            tuneGrid = knn.grid,
            metric = "Accuracy",
            trControl = ctrl,
            preProcess = c("center","scale")
          ))

          # Extract best k
          best_k <- temp$bestTune$k
          # Predictions
          y.train[[model_name]][cur.train] <- quiet(predict(temp, newdata = xtrain, type = "prob")[,"Class1"])

          y.pred[[model_name]][cur.test] <- quiet(predict(temp, newdata = xtest, type = "prob")[,"Class1"])
        }
        if (toupper(model_name) == "SVM") {
          message("SVM running...")
          params <- modifyList(indiv.defaults$SVM %||% list(), param.indiv[["SVM"]] %||% list())
          params$formula <- NULL
          params$data <- NULL
          params$x <- NULL
          params$y <- NULL
          train.df <- data.frame(y = factor(outcome[cur.train], levels = c(0, 1)), xtrain)
          temp <- quiet(do.call(svm, c(list(y ~ ., data = train.df, probability = TRUE), params)))
          y.train[[model_name]][cur.train] <- quiet(svm_binary_probability(temp, xtrain, positive = "1"))
          y.pred[[model_name]][cur.test] <- quiet(svm_binary_probability(temp, xtest, positive = "1"))
        }
        if (toupper(model_name) == "NB") {
          message("NB running...")
          params <- modifyList(indiv.defaults$NB %||% list(), param.indiv[["NB"]] %||% list())
          usekernel_flag <- params$usekernel %||% TRUE
          params$usekernel <- NULL
          y1<-as.factor(as.character(as.numeric(outcome)))
          dat.nb<- cbind.data.frame("class" = y1[cur.train], xtrain)
          params$formula <- class ~.
          params$data <- dat.nb
          params$x <- NULL
          params$y <- NULL
          temp <- quiet(do.call(naive_bayes, c(params, list(usekernel = usekernel_flag))))
          y.pred[[model_name]][cur.test] <- quiet(predict(temp, data.frame(xtest), type="prob")[,2])
          y.train[[model_name]][cur.train] <- quiet(predict(temp, data.frame(xtrain), type="prob")[,2])
        }
        if (toupper(model_name) == "DCT") {
          message("DecisionTree running...")
          y_outcome <- factor(outcome[cur.train], levels = c(0,1))
          train.df <- data.frame(y = y_outcome, xtrain)
          params <- modifyList(indiv.defaults$DCT %||% list(), param.indiv[["DCT"]] %||% list())
          params$x <- NULL
          params$y <- NULL
          params$formula <- y ~.
          params$data <- train.df
          params$method <- "class"
          temp <- quiet(do.call(rpart, params))
          y.pred[[model_name]][cur.test]=quiet(predict(temp, newdata=data.frame(xtest), type="prob")[,2])
          y.train[[model_name]][cur.train]=quiet(predict(temp, newdata=data.frame(xtrain), type="prob")[,2])
        }
        if (toupper(model_name) == "NN") {
          message("NeuralNet running...")
          params <- modifyList(indiv.defaults$NN, param.indiv[["NN"]] %||% list())
          params$x <- NULL
          params$y <- NULL

          m <- colMeans(as.data.frame(xtrain))
          s <- apply(as.data.frame(xtrain), 2, sd)
          xtrain.sd <- as.data.frame(scale(xtrain, center = m, scale = s))
          s[s==0] <- 1
          xtest.sd <- as.data.frame(scale(xtest, center = m, scale = s))
          colnames(xtrain.sd) <- colnames(xtrain)
          colnames(xtest.sd) <- colnames(xtest)
          train_df <- cbind(y = as.numeric(outcome[cur.train]), xtrain.sd)

          params$formula <- y ~.
          params$data <- train_df
          temp.nn <- quiet(do.call(neuralnet::neuralnet, params))
          y.pred[[model_name]][cur.test] <- quiet(as.vector(predict(temp.nn, newdata = xtest.sd)))
          y.train[[model_name]][cur.train] <- quiet(as.vector(predict(temp.nn, newdata = xtrain.sd)))
        }
      }

      if (ensemble) {
        message("Ensemble running...")
        xtrain.ens <- as.data.frame(lapply(y.train[models.indiv], `[`, cur.train))
        xtest.ens  <- as.data.frame(lapply(y.pred[models.indiv], `[`, cur.test))
        colnames(xtrain.ens) <- make.names(models.indiv, unique = TRUE)
        colnames(xtest.ens) <- colnames(xtrain.ens)
        xtest.ens <- xtest.ens[, colnames(xtrain.ens), drop = FALSE]
        xtrain.ens[xtrain.ens < 0] <- 0
        xtrain.ens[xtrain.ens > 1] <- 1
        xtest.ens[xtest.ens < 0] <- 0
        xtest.ens[xtest.ens > 1] <- 1

        for (model_name in models.ens) {
          if (toupper(model_name) == "XG") {
            message("Xg ensemble...")
            params <- modifyList(ens.defaults$XG, param.ens[["XG"]] %||% list())
            params$data  <- NULL
            params$label <- NULL
            if (!is.null(params$eta)) {
              params$learning_rate <- params$eta
              params$eta <- NULL
            }
            params$verbose <- NULL
            params$objective <- params$objective %||% "binary:logistic"
            params$verbosity <- params$verbosity %||% 0
            nrounds <- params$nrounds %||% 300
            params$nrounds <- NULL
            dtrain.ens <- xgboost::xgb.DMatrix(data  = as.matrix(xtrain.ens), label = ifelse(outcome[cur.train] == 1, 1, 0))
            dtest.ens <- xgboost::xgb.DMatrix(data = as.matrix(xtest.ens))
            temp <- quiet(xgboost::xgb.train(params = params, data = dtrain.ens, nrounds = nrounds))
            y.ens[[model_name]][cur.test] <- quiet(predict(temp, dtest.ens))
          }
          if (toupper(model_name) == "RF") {
            message("RF ensemble...")
            params <- modifyList(ens.defaults$RF, param.ens[["RF"]] %||% list())
            y <- factor(outcome[cur.train], levels = c(0,1))
            params$x <- as.data.frame(xtrain.ens)
            params$y <- y
            temp <- quiet(do.call(randomForest, params))
            y.ens[[model_name]][cur.test] = quiet(predict(temp, newdata=xtest.ens, type="prob")[,2])
          }
          if (toupper(model_name) == "EN") {
            message("EN ensemble...")
            en.params <- modifyList(ens.defaults$EN, param.ens[["EN"]] %||% list())
            train.params <- en.params[setdiff(names(en.params), "s")]
            temp <- quiet(do.call(cv.glmnet, c(list(x = as.matrix(xtrain.ens), y = outcome[cur.train]), train.params)))
            table(outcome[cur.train])
            s.val <- param.ens[["EN"]]$s %||% "lambda.1se" # user-specified or default
            y.ens[[model_name]][cur.test] <- quiet(predict(temp, newx = as.matrix(xtest.ens), type = "response", s = s.val))
          }
          if (toupper(model_name) == "NN") {
            message("NN ensemble...")
            zero.var <- which(apply(xtrain.ens, 2, sd) == 0)
            if (length(zero.var) > 0) {
              xtrain.ens <- xtrain.ens[, -zero.var, drop = FALSE]
              xtest.ens  <- xtest.ens[, -zero.var, drop = FALSE]
            }

            params <- modifyList(ens.defaults$NN, param.ens[["NN"]] %||% list())
            params$x <- NULL
            params$y <- NULL

            m <- colMeans(as.data.frame(xtrain.ens))
            s <- apply(as.data.frame(xtrain.ens), 2, sd)
            xtrain.sd <- as.data.frame(scale(xtrain.ens, center = m, scale = s))
            xtest.sd <- as.data.frame(scale(xtest.ens, center = m, scale = s))
            colnames(xtrain.sd) <- colnames(xtrain.ens)
            colnames(xtest.sd) <- colnames(xtest.ens)
            train_df <- cbind(y = outcome[cur.train], xtrain.sd)
            params$formula <- y ~.
            params$data <- train_df

            if (length(unique(train_df$y)) > 1) {
              temp.nn <- quiet(do.call(neuralnet::neuralnet, params))
              y.ens[[model_name]][cur.test] <- quiet(as.vector(predict(temp.nn, newdata = xtest.sd)))
            } else {
              warning("Only one class in this fold - skipping NN ensemble training")
              y.ens[[model_name]][cur.test] <- quiet(mean(train_df$y))  # fallback
            }
          }
          if (toupper(model_name) == "ADB") {
            message("ADB ensemble...")
            params <- modifyList(ens.defaults$ADB, param.ens[["ADB"]] %||% list())
            ada.sr <- ifelse(outcome == 0, -1, 1)
            params$formula <- NULL
            params$x <- NULL
            params$y <- NULL
            temp <- quiet(do.call(adaboost, c(list(as.matrix(xtrain.ens), ada.sr[cur.train]), params)))
            y.ens[[model_name]][cur.test]=quiet(predict(temp, as.matrix(xtest.ens), type = "prob"))
          }
          if (toupper(model_name) == "MB") {
            message("MB ensemble...")
            params <- modifyList(ens.defaults$MB, param.ens[["MB"]] %||% list())
            y.ens[[model_name]][cur.test] <- quiet(do.call(
              binary_mb_predict,
              c(list(xtrain = xtrain.ens, outcome = outcome[cur.train], xtest = xtest.ens), params)
            ))
          }
          if (toupper(model_name) == "GB") {
            message("GB ensemble...")
            params <- modifyList(ens.defaults$GB, param.ens[["GB"]] %||% list())
            params$formula <- y ~.
            params$data <- data.frame(y = outcome[cur.train], xtrain.ens)
            params$distribution <- "bernoulli"
            temp <- quiet(do.call(gbm, params))
            params$x <- NULL
            params$y <- NULL
            y.ens[[model_name]][cur.test]=quiet(predict.gbm(temp, data.frame(xtest.ens), n.trees = params$n.trees, type = "response"))
          }
          if (toupper(model_name) == "KNN") {
            y.fold <- outcome[cur.train]

            if (length(unique(y.fold)) < 2) {
              warning("Skipped fold: only one class")
              y.ens[[model_name]][cur.test] <- mean(y.fold)
              next
            }

            message("KNN ensemble...")

            knn.params <- modifyList(indiv.defaults$KNN, param.ens[["KNN"]] %||% list())

            y.train.factor <- factor(
              ifelse(outcome[cur.train] == 1, "Class1", "Class0"),
              levels = c("Class0","Class1")
            )

            # Adaptive k range (based on training size)
            max_k <- floor(sqrt(nrow(xtrain.ens)))
            k_seq <- seq(3, max_k, by = 2)
            if (length(k_seq) == 0) k_seq <- 3  # safeguard

            knn.grid <- expand.grid(k = k_seq)

            # Cross-validation instead of "none"
            ctrl <- trainControl(method = "cv", number = 5)

            temp <- quiet(caret::train(
              x = xtrain.ens,
              y = y.train.factor,
              method = "knn",
              tuneGrid = knn.grid,
              metric = "Accuracy",
              trControl = ctrl,
              preProcess = c("center","scale")
            ))

            # Optional: inspect best k
            best_k <- temp$bestTune$k
            message(paste("Best k (ensemble):", best_k))

            # Predict
            y.ens[[model_name]][cur.test] <- quiet(
              predict(temp, newdata = xtest.ens, type = "prob")[,"Class1"]
            )
          }
          if (toupper(model_name) == "SVM") {
            message("SVM ensemble...")
            params <- modifyList(ens.defaults$SVM %||% list(), param.ens[["SVM"]] %||% list())
            params$formula <- NULL
            params$data <- NULL
            params$x <- NULL
            params$y <- NULL
            train.df <- data.frame(y = factor(outcome[cur.train], levels = c(0, 1)), xtrain.ens)
            temp <- quiet(do.call(svm, c(list(y ~ ., data = train.df, probability = TRUE), params)))
            y.ens[[model_name]][cur.test] <- quiet(svm_binary_probability(temp, xtest.ens, positive = "1"))
          }
          if (toupper(model_name) == "NB") {
            message("NB ensemble...")
            params <- modifyList(ens.defaults$NB %||% list(), param.ens[["NB"]] %||% list())
            usekernel_flag <- params$usekernel %||% TRUE
            params$usekernel <- NULL
            y1<-as.factor(as.character(as.numeric(outcome)))
            dat.nb<- cbind.data.frame("class" = y1[cur.train], xtrain.ens)
            params$formula <- class ~.
            params$data <- dat.nb
            params$x <- NULL
            params$y <- NULL
            temp <- quiet(do.call(naive_bayes, c(params, list(usekernel = usekernel_flag))))
            y.ens[[model_name]][cur.test] <- quiet(predict(temp, data.frame(xtest.ens), type="prob")[,2])
          }
          if (toupper(model_name) == "DCT") {
            message("DCT ensemble...")
            y_outcome <- factor(outcome[cur.train], levels = c(0,1))
            train.df <- data.frame(y = y_outcome, xtrain.ens)
            params <- modifyList(ens.defaults$DCT %||% list(), param.ens[["DCT"]] %||% list())
            params$x <- NULL
            params$y <- NULL
            params$formula <- y ~.
            params$data <- train.df
            params$method <- "class"
            temp <- quiet(do.call(rpart, params))
            y.ens[[model_name]][cur.test]=quiet(predict(temp, newdata=data.frame(xtest.ens), type="prob")[,2])
          }
        }
      }

      # Weighted average logistic regression
      y.pred.train.all <- as.data.frame(lapply(y.train, `[`, cur.train))
      y.pred.mat.all   <- as.data.frame(lapply(y.pred, `[`, cur.test))
      y.pred.train.all[y.pred.train.all < 0] <- 0
      y.pred.train.all[y.pred.train.all > 1] <- 1
      y.pred.mat.all[y.pred.mat.all < 0] <- 0
      y.pred.mat.all[y.pred.mat.all > 1] <- 1

      formula_str <- paste("outcome ~", paste(names(y.pred.train.all), collapse = " + "))
      formula_obj <- as.formula(formula_str)
      temp.avg.all <- quiet(glm(formula_obj, family = binomial(link = "logit"),
                                data = cbind(outcome = outcome[cur.train], y.pred.train.all)))
      y.pred.avg.all[cur.test] <- quiet(predict(temp.avg.all, newdata = y.pred.mat.all, type = "response"))
    }

    # Simple average
    y.pred.mat1.all <- as.data.frame(y.pred)
    y.pred.mat1.all[y.pred.mat1.all < 0] <- 0
    y.pred.mat1.all[y.pred.mat1.all > 1] <- 1
    y.pred.avg1.all <- quiet(apply(y.pred.mat1.all, 1, mean, na.rm = TRUE))

    # Pooling
    for (model in models.indiv) pred_pools[[model]] <- c(pred_pools[[model]], y.pred[[model]])
    if (ensemble) {
      for (model in models.ens) ens_pools[[model]] <- c(ens_pools[[model]], y.ens[[model]])
    }
    y.pred.avg.pool.all  <- c(y.pred.avg.pool.all, y.pred.avg.all)
    y.pred.avg.pool1.all <- c(y.pred.avg.pool1.all, y.pred.avg1.all)
  }

  # Return all pooled results
  ret_list <- list(
    individual = pred_pools,
    ensemble   = if (ensemble) ens_pools else NULL,
    y.pred.avg.pool.all  = y.pred.avg.pool.all,
    y.pred.avg.pool1.all = y.pred.avg.pool1.all
  )

  return(ret_list)
}
survival_model <- function(outcome, survdays, x, models.indiv, iter, num.folds = NULL, ensemble = FALSE, models.ens = NULL, param.indiv = NULL, param.ens = NULL, time_point = NULL, iter.offset = 0) {
  nfold <- if (!is.null(num.folds)) num.folds else 3
  indiv.defaults <- list(
    XG = list(eta = 0.1, objective = "survival:cox", nrounds = 1000, verbosity = 0),
    RF = list(ntree = 2000, max.depth = 10, importance = TRUE),
    EN = list(alpha = 0.2, s = "lambda.min"),
    MB = list(mstop = 300),
    ADB = list(tree_depth = 3, n_rounds = 200, verbose = FALSE),
    GB = list(n.trees = 1000, interaction.depth = 10, shrinkage = 0.005, bag.fraction = 0.5),
    KNN = list(k = 5),
    NN = list(size = 5, decay = 1e-4, maxit = 200, trace = FALSE),
    SVM = list(gamma.mu = 0.01, opt.meth = "quadprog", kernel = "add_kernel"),
    DCT = list(minsplit = 10, maxdepth = 10, cp = 0.001)
  )


  # Initialize pooled results as lists
  pred_pools <- list()
  for (model in models.indiv) pred_pools[[model]] <- NULL


  ens_pools <- list()
  if (ensemble) {
    for (model in models.ens) ens_pools[[model]] <- NULL
  }


  y.pred.avg.pool.all <- NULL
  y.pred.avg.pool1.all <- NULL


  for (j in iter.offset + seq_len(iter)) {
    message(j)
    set.seed(j)


    cv.fold <-balanced.folds.vec(outcome,nfold)

    # Initialize predictions for each iteration
    y.pred <- list()
    y.train <- list()
    y.ens <- list()

    for (m in models.indiv) {
      y.pred[[m]] <- rep(NA, length(outcome))
      y.train[[m]] <- rep(NA, length(outcome))
    }
    if (ensemble) {
      for (m in models.ens) y.ens[[m]] <- rep(NA, length(outcome))
    }

    y.pred.avg.iter <- rep(NA, length(outcome))  # full length of your data

    # Cross-validation loop
    for (i in 1:nfold) {
      cur.train <- which(cv.fold != i)
      cur.test <- which(cv.fold == i)

      xtrain <- x[cur.train, ]
      xtest  <- x[cur.test, ]

      xtrain_time <- survdays[cur.train]
      xtrain_status <- outcome[cur.train]
      xtest_time <- survdays[cur.test]
      xtest_status <- outcome[cur.test]

      # Run individual models
      for (model_name in models.indiv) {
        if (toupper(model_name) == "XG") {
          message("XgBoost running...")
          params_user <- modifyList(indiv.defaults$XG, param.indiv[["XG"]] %||% list())
          # Cox label encoding
          label_train <- ifelse(xtrain_status == 1,
                                xtrain_time,
                                -xtrain_time)
          # validation split
          val_ind <- sample(seq_len(nrow(xtrain)),
                            0.1 * nrow(xtrain),
                            replace = FALSE)
          dtrain <- xgboost::xgb.DMatrix(
            data = as.matrix(xtrain[-val_ind,]),
            label = label_train[-val_ind]
          )
          dval <- xgboost::xgb.DMatrix(
            data = as.matrix(xtrain[val_ind,]),
            label = label_train[val_ind]
          )
          fixed_params <- list(
            objective = "survival:cox",
            eval_metric = "cox-nloglik"
          )
          final_params <- modifyList(fixed_params, params_user)
          nrounds_val <- final_params$nrounds %||% 1000
          final_params$nrounds <- NULL
          temp <- quiet(xgboost::xgb.train(
            params = final_params,
            data = dtrain,
            nrounds = nrounds_val,
            evals = list(validation = dval),
            verbose = 0
          ))

          # risk score predictions
          y.train[[model_name]][cur.train] <- quiet(predict(temp, newdata = xgboost::xgb.DMatrix(as.matrix(xtrain))))
          y.pred[[model_name]][cur.test] <- quiet(predict(temp, newdata = xgboost::xgb.DMatrix(as.matrix(xtest))))
        }
        if (toupper(model_name) == "RF") {
          message("RandomForest running...")
          params <- modifyList(indiv.defaults$RF, param.indiv[["RF"]] %||% list())
          surv_train_data <- data.frame(time = xtrain_time, status = xtrain_status, xtrain)
          surv_test_data <- data.frame(time = xtest_time, status = xtest_status, xtest)
          formula_args <- list(formula = Surv(time, status) ~., data = surv_train_data)
          temp <- quiet(do.call(randomForestSRC::rfsrc, c(formula_args, params)))
          predictions_train <- quiet(predict(temp, newdata = surv_train_data, type = "response"))
          predictions_test <- quiet(predict(temp, newdata = data.frame(xtest), type = "response"))
          y.train[[model_name]][cur.train] <- predictions_train$predicted
          y.pred[[model_name]][cur.test] <- predictions_test$predicted
        }
        if (toupper(model_name) == "EN") {

          message("ElasticNet running...")

          params <- modifyList(indiv.defaults$EN, param.indiv[["EN"]] %||% list())

          surv_train_object <- Surv(time = xtrain_time, event = xtrain_status)

          train.params <- params[setdiff(names(params), "s")]

          cv_glmnet <- quiet(do.call(glmnet::cv.glmnet, c(
            list(
              x = as.matrix(xtrain),
              y = surv_train_object,
              family = "cox",
              nfolds = 5
            ),
            train.params
          )))

          # choose lambda AFTER fitting
          if (is.null(params$s)) {
            lambda_use <- cv_glmnet$lambda.min
          } else if (is.character(params$s)) {
            lambda_use <- cv_glmnet[[params$s]]
          } else {
            lambda_use <- params$s
          }

          risk_scores_train <- quiet(predict(
            cv_glmnet,
            newx = as.matrix(xtrain),
            s = lambda_use,
            type = "link"
          ))

          risk_scores_test <- quiet(predict(
            cv_glmnet,
            newx = as.matrix(xtest),
            s = lambda_use,
            type = "link"
          ))
          y.pred[[model_name]][cur.test] <- as.numeric(risk_scores_test)
          y.train[[model_name]][cur.train] <- as.numeric(risk_scores_train)

        }
        if (toupper(model_name) == "MB") {
          message("MBoost running...")
          params <- modifyList(indiv.defaults$MB, param.indiv[["MB"]] %||% list())
          mstop_val <- params$mstop %||% 300
          params$mstop <- NULL

          control_object <- mboost::boost_control(mstop = mstop_val)

          x_matrix <- as.matrix(xtrain)
          y_surv <- Surv(time = xtrain_time, event = xtrain_status)

          fixed_args <- list(
            x = x_matrix,
            y = y_surv,
            family = mboost::CoxPH(), # Fixed survival family
            control = control_object
          )

          temp <- quiet(do.call(mboost::glmboost, c(fixed_args, params)))
          risk_scores_train <- quiet(predict(temp, newdata = as.matrix(xtrain), type = "link"))
          risk_scores_test <- quiet(predict(temp, newdata = as.matrix(xtest), type = "link"))
          y.pred[[model_name]][cur.test] <- risk_scores_test
          y.train[[model_name]][cur.train] <- risk_scores_train
        }
        if (toupper(model_name) == "ADB") {
          message("AdaBoost survival running...")
          params <- modifyList(indiv.defaults$ADB, param.indiv[["ADB"]] %||% list())
          y.train[[model_name]][cur.train] <- quiet(do.call(
            survival_adb_predict,
            c(list(xtrain = xtrain, time = xtrain_time, status = xtrain_status,
                   time_point = time_point), params)
          ))
          y.pred[[model_name]][cur.test] <- quiet(do.call(
            survival_adb_predict,
            c(list(xtrain = xtrain, time = xtrain_time, status = xtrain_status,
                   time_point = time_point, xtest = xtest), params)
          ))
        }
        if (toupper(model_name) == "GB") {
          message("GradientBoost running...")
          params <- modifyList(indiv.defaults$GB, param.indiv[["GB"]] %||% list())
          surv_train_data <- data.frame(time = xtrain_time, status = xtrain_status, xtrain)
          surv_test_data <- data.frame(time = xtest_time, status = xtest_status, xtest)
          fixed_args <- list(
            formula = Surv(time, status) ~.,
            data = surv_train_data,
            distribution = "coxph" # Fixed survival distribution
          )
          n.trees_pred <- params$n.trees %||% 500
          temp <- quiet(do.call(gbm::gbm, c(fixed_args, params)))
          risk_scores_train <- quiet(predict(temp, newdata = surv_train_data, n.trees = n.trees_pred, type = "link"))
          risk_scores_test <- quiet(predict(temp, newdata = surv_test_data, n.trees = n.trees_pred, type = "link"))

          y.pred[[model_name]][cur.test] <- risk_scores_test
          y.train[[model_name]][cur.train] <- risk_scores_train
        }
        if (toupper(model_name) == "KNN") {
          message("KNN survival running...")
          params <- modifyList(indiv.defaults$KNN, param.indiv[["KNN"]] %||% list())
          k_val <- params$k %||% 5
          y.train[[model_name]][cur.train] <- survival_knn_predict(
            xtrain, xtrain_time, xtrain_status, k = k_val, leave_one_out = TRUE
          )
          y.pred[[model_name]][cur.test] <- survival_knn_predict(
            xtrain, xtrain_time, xtrain_status, xtest, k = k_val
          )
        }
        if (toupper(model_name) == "NN") {
          message("NN survival running...")
          params <- modifyList(indiv.defaults$NN, param.indiv[["NN"]] %||% list())
          y.train[[model_name]][cur.train] <- do.call(
            survival_nn_predict,
            c(list(xtrain = xtrain, time = xtrain_time, status = xtrain_status), params)
          )
          y.pred[[model_name]][cur.test] <- do.call(
            survival_nn_predict,
            c(list(xtrain = xtrain, time = xtrain_time, status = xtrain_status,
                   xtest = xtest), params)
          )
        }
        if (toupper(model_name) == "SVM") {
          message("SVM running...")

          result <- tryCatch({
            params <- modifyList(indiv.defaults$SVM, param.indiv[["SVM"]] %||% list())
            surv_train_data <- data.frame(time = xtrain_time, status = xtrain_status, xtrain)

            formula_args <- list(
              formula = Surv(time, status) ~ .,
              data = surv_train_data,
              type = "regression"
            )

            temp <- quiet(do.call(survivalsvm::survivalsvm, c(formula_args, params)))

            risk_scores_train <- quiet(predict(object = temp, newdata = data.frame(xtrain)))
            risk_scores_test  <- quiet(predict(object = temp, newdata = data.frame(xtest)))

            list(
              train = as.numeric(risk_scores_train$predicted),
              test  = as.numeric(risk_scores_test$predicted)
            )

          }, error = function(e) {
            message("SVM failed: ", e$message)
            return(NULL)
          })

          # Only assign if it worked
          if (!is.null(result)) {
            y.pred[[model_name]][cur.test]  <- result$test
            y.train[[model_name]][cur.train] <- result$train
          }
        }
        if (toupper(model_name) == "DCT") {
          message("DecisionTree running...")
          params <- modifyList(indiv.defaults$DCT, param.indiv[["DCT"]] %||% list())
          surv_train_data <- data.frame(time = xtrain_time, status = xtrain_status, xtrain)
          formula_args <- list(
            formula = Surv(time, status) ~ .,
            data = surv_train_data,
            method = "exp"
          )

          # The user/default parameters (minsplit, maxdepth, cp) go into rpart.control
          control_params <- params[names(params)]

          # Execute training with control() wrapper
          temp <- quiet(do.call(rpart::rpart, c(formula_args, list(control = do.call(rpart::rpart.control, control_params)))))

          # Prediction
          risk_scores_train <- quiet(predict(temp, newdata = data.frame(xtrain)))
          risk_scores_test <- quiet(predict(temp, newdata = data.frame(xtest)))

          y.pred[[model_name]][cur.test] <- risk_scores_test
          y.train[[model_name]][cur.train] <- risk_scores_train
        }
      }

      # Ensemble models
      if (ensemble) {
        message("Ensemble running...")
        xtrain.ens <- as.data.frame(lapply(y.train[models.ens], `[`, cur.train))
        xtest.ens  <- as.data.frame(lapply(y.pred[models.ens], `[`, cur.test))

        max.train <- apply(xtrain.ens, 2, max)
        xtrain.ens[, which(max.train > 1e30)] <- log(xtrain.ens[, which(max.train > 1e30)])

        max.test <- apply(xtest.ens, 2, max)
        xtest.ens[, which(max.test > 1e30)] <- log(xtest.ens[, which(max.test > 1e30)])

        surv_ensemble_train_data <- data.frame(time = xtrain_time, status = xtrain_status, xtrain.ens)
        surv_ensemble_test_data <- data.frame(time = xtest_time, status = xtest_status, xtest.ens)
        safe_names <- c("time", "status", paste0(models.ens, "_pred"))
        colnames(surv_ensemble_train_data) <- safe_names
        colnames(surv_ensemble_test_data) <- safe_names

        for (model_name in models.ens) {
          if (toupper(model_name) == "XG") {
            message("Xg ensemble...")

            params_xg <- modifyList(
              list(
                objective = "survival:cox",
                eval_metric = "cox-nloglik",
                eta = 0.1,
                tree_method = "hist",
                nrounds = 50
              ),
              param.ens[[model_name]] %||% list()
            )

            # extract nrounds (REQUIRED for >=3)
            nrounds_val <- params_xg$nrounds %||% 50
            params_xg$nrounds <- NULL

            # validation split
            val_ind <- sample(
              seq_len(nrow(xtrain.ens)),
              floor(0.1 * nrow(xtrain.ens))
            )

            # Cox label encoding
            label_train <- ifelse(
              xtrain_status == 1,
              xtrain_time,
              -xtrain_time
            )

            # training matrices
            dtrain <- xgboost::xgb.DMatrix(
              data = as.matrix(xtrain.ens[-val_ind, ]),
              label = label_train[-val_ind]
            )

            dval <- xgboost::xgb.DMatrix(
              data = as.matrix(xtrain.ens[val_ind, ]),
              label = label_train[val_ind]
            )

            temp <- quiet(xgboost::xgb.train(
              params = params_xg,
              data = dtrain,
              nrounds = nrounds_val,
              evals = list(validation = dval),
              verbose = 0
            ))

            # prediction
            dtest <- xgboost::xgb.DMatrix(
              data = as.matrix(xtest.ens)
            )

            y.ens[[model_name]][cur.test] <- quiet(predict(temp, dtest))
          }
          if (toupper(model_name) == "RF") {
            message("RF ensemble...")
            params <- modifyList(list(ntree = 2000, max.depth = 10, importance = TRUE),
                                 param.ens[[model_name]] %||% list())
            temp <- quiet(do.call(rfsrc, c(list(formula = Surv(time,status) ~ ., data = surv_ensemble_train_data), params)))
            predictions_test <- predict(temp, newdata = surv_ensemble_test_data, type = "response")
            y.ens[[model_name]][cur.test] <- quiet(predictions_test$predicted)
          }
          if (toupper(model_name) == "EN") {
            message("EN ensemble...")
            params <- modifyList(list(alpha = 0.1), param.ens[[model_name]] %||% list())
            surv_train_object <- Surv(time = xtrain_time, event = xtrain_status)
            cv_glmnet <- do.call(cv.glmnet, c(list(x = as.matrix(xtrain.ens), y = surv_train_object, family = "cox", nfolds = 5), params))
            temp <- quiet(glmnet(as.matrix(xtrain.ens), surv_train_object, family = "cox", lambda = cv_glmnet$lambda.min, alpha = params$alpha %||% 0.1))
            y.ens[[model_name]][cur.test] <- quiet(predict(temp, newx = as.matrix(xtest.ens), s = cv_glmnet$lambda.min, type = "response"))
          }
          if (toupper(model_name) == "MB") {
            message("MBoost ensemble...")
            params <- modifyList(list(mstop = 300), param.ens[[model_name]] %||% list())
            control_object <- boost_control(mstop = params$mstop %||% 300)
            temp <- quiet(glmboost(Surv(time, status) ~ ., data = surv_ensemble_train_data, family = CoxPH(), control = control_object))
            y.ens[[model_name]][cur.test] <- quiet(predict(temp, newdata = surv_ensemble_test_data, type = "link"))
          }
          if (toupper(model_name) == "ADB") {
            message("ADB survival ensemble...")
            params <- modifyList(indiv.defaults$ADB, param.ens[[model_name]] %||% list())
            y.ens[[model_name]][cur.test] <- quiet(do.call(
              survival_adb_predict,
              c(list(xtrain = xtrain.ens, time = xtrain_time, status = xtrain_status,
                     time_point = time_point, xtest = xtest.ens), params)
            ))
          }
          if (toupper(model_name) == "GB") {
            message("GB ensemble...")
            summary(surv_ensemble_train_data$XG_pred)
            params <- modifyList(list(n.trees = 1000, interaction.depth = 10, shrinkage = 0.005, bag.fraction = 0.5), param.ens[[model_name]] %||% list())
            temp <- quiet(do.call(gbm, c(list(formula = Surv(time, status) ~ ., data = surv_ensemble_train_data, distribution = "coxph"), params)))
            y.ens[[model_name]][cur.test] <- quiet(predict(temp, newdata = surv_ensemble_test_data, n.trees = params$n.trees %||% 1000, type = "link"))
          }
          if (toupper(model_name) == "KNN") {
            message("KNN ensemble...")
            params <- modifyList(list(k = 5), param.ens[[model_name]] %||% list())
            y.ens[[model_name]][cur.test] <- survival_knn_predict(
              xtrain.ens, xtrain_time, xtrain_status, xtest.ens, k = params$k %||% 5
            )
          }
          if (toupper(model_name) == "NN") {
            message("NN ensemble...")
            params <- modifyList(indiv.defaults$NN, param.ens[[model_name]] %||% list())
            y.ens[[model_name]][cur.test] <- do.call(
              survival_nn_predict,
              c(list(xtrain = xtrain.ens, time = xtrain_time, status = xtrain_status,
                     xtest = xtest.ens), params)
            )
          }
          if (toupper(model_name) == "SVM") {
            message("SVM ensemble...")

            result <- tryCatch({

              params <- modifyList(
                list(gamma.mu = 0.5, opt.meth = "quadprog", kernel = "add_kernel"),
                param.ens[[model_name]] %||% list()
              )

              temp <- quiet(do.call(
                survivalsvm,
                c(
                  list(
                    formula = Surv(time, status) ~ .,
                    data = surv_ensemble_train_data,
                    type = "regression"
                  ),
                  params
                )
              ))

              risk_scores_test <- quiet(predict(temp, newdata = surv_ensemble_test_data))

              as.numeric(risk_scores_test$predicted)

            }, error = function(e) {
              message("SVM ensemble failed: ", e$message)
              return(NULL)
            })

            # Only assign if it worked
            if (!is.null(result)) {
              y.ens[[model_name]][cur.test] <- result
            }
          }
          if (toupper(model_name) == "DCT") {
            message("DCT ensemble...")
            params <- modifyList(list(minsplit = 20, maxdepth = 10, cp = 0.01), param.ens[[model_name]] %||% list())
            temp <- quiet(rpart(Surv(time, status) ~ ., data = surv_ensemble_train_data, method = "exp", control = do.call(rpart.control, params)))
            risk_scores_test <- quiet(predict(temp, newdata = surv_ensemble_test_data))
            y.ens[[model_name]][cur.test] <- unname(risk_scores_test)
          }
        }
      }

      message("Weighted average (Cox stacking)...")

      # --------------------------------------------------
      # 1. Stack base-model predictions
      # --------------------------------------------------
      y.pred.train.all <- as.data.frame(lapply(y.train, `[`, cur.train))
      y.pred.test.all  <- as.data.frame(lapply(y.pred,  `[`, cur.test))

      time.train   <- survdays[cur.train]
      status.train <- outcome[cur.train]

      # --------------------------------------------------
      # 2. Remove rows containing NA / Inf
      # --------------------------------------------------
      valid_rows <- apply(
        y.pred.train.all,
        1,
        function(x) all(is.finite(x))
      )

      y.pred.train.all <- y.pred.train.all[valid_rows, , drop = FALSE]
      time.train       <- time.train[valid_rows]
      status.train     <- status.train[valid_rows]

      # --------------------------------------------------
      # 3. Remove zero-variance predictors
      # --------------------------------------------------
      zero_var <- sapply(y.pred.train.all, function(x)
        var(x, na.rm = TRUE) < 1e-12
      )

      if (any(zero_var)) {
        y.pred.train.all <- y.pred.train.all[, !zero_var, drop = FALSE]
        y.pred.test.all  <- y.pred.test.all[,  !zero_var, drop = FALSE]
      }

      # --------------------------------------------------
      # 4. Standardize predictors safely
      # --------------------------------------------------
      mu  <- colMeans(y.pred.train.all)
      sdv <- apply(y.pred.train.all, 2, sd)

      # prevent divide-by-zero
      sdv[sdv < 1e-8] <- 1

      y.pred.train.all <- scale(y.pred.train.all, center = mu, scale = sdv)
      y.pred.test.all  <- scale(y.pred.test.all,  center = mu, scale = sdv)

      y.pred.train.all <- as.data.frame(y.pred.train.all)
      y.pred.test.all  <- as.data.frame(y.pred.test.all)

      # --------------------------------------------------
      # 5. Add survival outcome
      # --------------------------------------------------
      y.pred.train.all$time   <- time.train
      y.pred.train.all$status <- status.train

      # --------------------------------------------------
      # 6. Remove bad columns AGAIN (after scaling)
      # --------------------------------------------------
      good_cols <- sapply(y.pred.train.all, function(x)
        all(is.finite(x)) && sd(x) > 1e-8
      )

      good_cols[c("time","status")] <- TRUE

      y.pred.train.all <- y.pred.train.all[, good_cols, drop = FALSE]

      predictors <- setdiff(colnames(y.pred.train.all), c("time","status"))

      # align test data EXACTLY
      y.pred.test.all <- y.pred.test.all[, predictors, drop = FALSE]

      formula_obj <- as.formula(
        paste("Surv(time, status) ~",
              paste(predictors, collapse = " + "))
      )

      temp.avg.all <- quiet(survival::coxph(
        formula_obj,
        data = y.pred.train.all,
        ties = "breslow",
        x = TRUE,
        model = TRUE
      ))

      y.pred.avg <- quiet(as.numeric(
        predict(
          temp.avg.all,
          newdata = y.pred.test.all,
          type = "lp"
        )
      ))

      y.pred.avg.iter[cur.test] <- y.pred.avg
    }

    # Simple average
    message("Simple average...")
    y.pred.mat1.all <- as.data.frame(y.pred)
    y.pred.mat1.all[y.pred.mat1.all < 0] <- 0
    y.pred.mat1.all[y.pred.mat1.all > 1] <- 1
    y.pred.avg1.all <- apply(y.pred.mat1.all, 1, mean)

    # Pooling
    message("pooling...")
    for (model in models.indiv) pred_pools[[model]] <- c(pred_pools[[model]], y.pred[[model]])
    if (ensemble) {
      for (model in models.ens) ens_pools[[model]] <- c(ens_pools[[model]], y.ens[[model]])
    }
    y.pred.avg.pool.all <- c(y.pred.avg.pool.all, y.pred.avg.iter)
    y.pred.avg.pool1.all <- c(y.pred.avg.pool1.all, y.pred.avg1.all)
  }


  # Return all pooled results
  ret_list <- list(
    individual = pred_pools,
    ensemble   = if (ensemble) ens_pools else NULL,
    y.pred.avg.pool.all  = y.pred.avg.pool.all,
    y.pred.avg.pool1.all = y.pred.avg.pool1.all
  )

  return(ret_list)
}
