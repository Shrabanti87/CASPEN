# Internal functions - not exported

#=======BINARY (Works) ========
binary_train_performance <- function(outcome, train_data, features, models.indiv, iter, num.folds = NULL, ensemble = FALSE, models.ens = NULL, param.indiv = NULL, param.ens = NULL, iter.offset = 0) {
  param.indiv <- param.indiv %||% list()
  param.ens <- param.ens %||% list()
  models.indiv <- toupper(models.indiv)
  models.ens <- if (isTRUE(ensemble)) toupper(models.ens %||% models.indiv) else character(0)
  nfold <- max(2L, as.integer(num.folds %||% 3L))

  indiv.defaults <- list(
    XG = list(
      max_depth = 1, eta = 0.03, nrounds = 80,
      objective = "binary:logistic", eval_metric = "logloss",
      min_child_weight = 5, subsample = 0.8, colsample_bytree = 1,
      lambda = 10, alpha = 1, nthread = 2, verbosity = 0
    ),
    RF = list(ntree = 500),
    EN = list(alpha = 0.2, family = "binomial", s = "lambda.1se"),
    ADB = list(tree_depth = 3, n_rounds = 200, verbose = FALSE),
    MB = list(mstop = 300),
    GB = list(cv.folds = 10, shrinkage = 0.01, n.minobsinnode = 2, n.trees = 500, distribution = "bernoulli", n.cores = 1),
    KNN = list(k = 5),
    SVM = list(),
    NN = list(hidden = c(8,4), rep = 1),
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
    NN = list(hidden = c(8,4), rep = 1),
    NB = list(usekernel = TRUE),
    DCT = list(method = "class")
  )
  if (ensemble == TRUE) {
    for (model_name in models.ens) {
      if (grepl("KNN", model_name, ignore.case = TRUE)) {
        y.ens.KNN.pool <- NULL
      }
      if (grepl("XG", model_name, ignore.case = TRUE)) {
        y.ens.XG.pool <- NULL
      }
      if (grepl("RF", model_name, ignore.case = TRUE)) {
        y.ens.RF.pool <- NULL
      }
      if (grepl("EN", model_name, ignore.case = TRUE)) {
        y.ens.EN.pool <- NULL
      }
      if (grepl("ADB", model_name, ignore.case = TRUE)) {
        y.ens.ADB.pool <- NULL
      }
      if (grepl("MB", model_name, ignore.case = TRUE)) {
        y.ens.MB.pool <- NULL
      }
      if (grepl("NN", model_name, ignore.case = TRUE)) {
        y.ens.NN.pool <- NULL
      }
      if (grepl("GB", model_name, ignore.case = TRUE)) {
        y.ens.GB.pool <- NULL
      }
      if (grepl("SVM", model_name, ignore.case = TRUE)) {
        y.ens.SVM.pool <- NULL
      }
      if (grepl("DCT", model_name, ignore.case = TRUE)) {
        y.ens.DCT.pool <- NULL
      }
      if (grepl("NB", model_name, ignore.case = TRUE)) {
        y.ens.NB.pool <- NULL
      }
    }
  }

  y.pred.avg.pool.all <- NULL
  y.pred.avg.pool1.all <- NULL

  for(j in 1:iter){
    iter.id <- j + iter.offset
    message(iter.id)
    set.seed(iter.id)

    cv.fold <- balanced.folds.vec(outcome, nfold)

    if (ensemble == TRUE) {
      for (model_name in models.ens) {
        if (toupper(model_name) == "KNN") {
          y.ens.KNN <- rep(NA, length(outcome))
        }
        if (toupper(model_name) == "XG") {
          y.ens.XG=rep(NA, length(outcome))
        }
        if (toupper(model_name) == "RF") {
          y.ens.RF=rep(NA, length(outcome))
        }
        if (toupper(model_name) == "EN") {
          y.ens.EN=rep(NA, length(outcome))
        }
        if (toupper(model_name) == "ADB") {
          y.ens.ADB<-rep(NA, length(outcome))
        }
        if (toupper(model_name) == "MB") {
          y.ens.MB<-rep(NA, length(outcome))
        }
        if (toupper(model_name) == "NN") {
          y.ens.NN <- rep(NA, length(outcome))
        }
        if (toupper(model_name) == "GB") {
          y.ens.GB<- rep(NA, length(outcome))
        }
        if (toupper(model_name) == "SVM") {
          y.ens.SVM <- rep(NA, length(outcome))
        }
        if (toupper(model_name) == "DCT") {
          y.ens.DCT <- rep(NA, length(outcome))
        }
        if (toupper(model_name) == "NB") {
          y.ens.NB <- rep(NA, length(outcome))
        }
      }
    }

    y.pred.avg.all<- rep(NA, length(outcome))
    mat.test.full<- NULL
    for(i in seq_len(nfold))
    {
      cur.train=which(cv.fold!=i)
      cur.test=which(cv.fold==i)

      mat.train<- NULL
      mat.test<- NULL
      y.pred.train.all<- NULL
      y.pred.mat.all<- NULL

      for(f in seq_along(features))
      {
        message(paste("Pathway Number:", f))
        com<-unique(features[[f]])
        data_train <- train_data[, colnames(train_data) %in% com, drop = FALSE]
        colnames(data_train) <- make.names(colnames(data_train), unique = TRUE)

        xtrain<-(data_train)[cur.train,]
        xtest<-(data_train)[cur.test,]

        for (model_name in models.indiv) {
          if (toupper(model_name) == "XG") {
            assign(paste0("y.pred.XG_",f), rep(NA, length(cur.test)))
            assign(paste0("y.train_tr.XG_",f), rep(NA, length(cur.train)))
            assign(paste0("y.train_ts.XG_",f), rep(NA, length(cur.test)))
          }
          if (toupper(model_name) == "EN") {
            assign(paste0("y.pred.EN_",f), rep(NA, length(cur.test)))
            assign(paste0("y.train_tr.EN_",f), rep(NA, length(cur.train)))
            assign(paste0("y.train_ts.EN_",f), rep(NA, length(cur.test)))
          }
          if (toupper(model_name) == "RF") {
            assign(paste0("y.pred.RF_",f), rep(NA, length(cur.test)))
            assign(paste0("y.train_tr.RF_",f), rep(NA, length(cur.train)))
            assign(paste0("y.train_ts.RF_",f), rep(NA, length(cur.test)))
          }
          if (toupper(model_name) == "ADB") {
            assign(paste0("y.pred.ADB_",f), rep(NA, length(cur.test)))
            assign(paste0("y.train_tr.ADB_",f), rep(NA, length(cur.train)))
            assign(paste0("y.train_ts.ADB_",f), rep(NA, length(cur.test)))
          }
          if (toupper(model_name) == "MB") {
            assign(paste0("y.pred.MB_",f), rep(NA, length(cur.test)))
            assign(paste0("y.train_tr.MB_",f), rep(NA, length(cur.train)))
            assign(paste0("y.train_ts.MB_",f), rep(NA, length(cur.test)))
          }
          if (toupper(model_name) == "NN") {
            assign(paste0("y.pred.NN_",f), rep(NA, length(cur.test)))
            assign(paste0("y.train_tr.NN_",f), rep(NA, length(cur.train)))
            assign(paste0("y.train_ts.NN_",f), rep(NA, length(cur.test)))
          }
          if (toupper(model_name) == "KNN") {
            assign(paste0("y.pred.KNN_",f), rep(NA, length(cur.test)))
            assign(paste0("y.train_tr.KNN_",f), rep(NA, length(cur.train)))
            assign(paste0("y.train_ts.KNN_",f), rep(NA, length(cur.test)))
          }
          if (toupper(model_name) == "GB") {
            assign(paste0("y.pred.GB_",f), rep(NA, length(cur.test)))
            assign(paste0("y.train_tr.GB_",f), rep(NA, length(cur.train)))
            assign(paste0("y.train_ts.GB_",f), rep(NA, length(cur.test)))
          }
          if (toupper(model_name) == "SVM") {
            assign(paste0("y.pred.SVM_",f), rep(NA, length(cur.test)))
            assign(paste0("y.train_tr.SVM_",f), rep(NA, length(cur.train)))
            assign(paste0("y.train_ts.SVM_",f), rep(NA, length(cur.test)))
          }
          if (toupper(model_name) == "NB") {
            assign(paste0("y.pred.NB_",f), rep(NA, length(cur.test)))
            assign(paste0("y.train_tr.NB_",f), rep(NA, length(cur.train)))
            assign(paste0("y.train_ts.NB_",f), rep(NA, length(cur.test)))
          }
          if (toupper(model_name) == "DCT") {
            assign(paste0("y.pred.DCT_",f), rep(NA, length(cur.test)))
            assign(paste0("y.train_tr.DCT_",f), rep(NA, length(cur.train)))
            assign(paste0("y.train_ts.DCT_",f), rep(NA, length(cur.test)))
          }
        }
        for (model_name in models.indiv) {
          if (toupper(model_name) == "XG") {
            message("XgBoost running...")
            params_xg <- modifyList(indiv.defaults$XG, param.indiv[["XG"]] %||% list())
            dtrain <- quiet(xgboost::xgb.DMatrix(data = as.matrix(xtrain), label = as.numeric(outcome[cur.train])))
            dtest <- xgboost::xgb.DMatrix(data = as.matrix(xtest))
            params_xg$data <- NULL
            params_xg$label <- NULL
            if (!is.null(params_xg$eta)) {
              params_xg$learning_rate <- params_xg$eta
              params_xg$eta <- NULL
            }
            params_xg$verbose <- NULL
            params_xg$verbosity <- 0
            params_xg$x <- NULL
            params_xg$y <- NULL
            temp <- quiet(xgboost::xgb.train(params = params_xg, data = dtrain, nrounds = params_xg$nrounds %||% 500))
            assign(paste0("y.pred.XG_",f), quiet(predict(temp, dtest)))
            assign(paste0("y.train_tr.XG_",f), quiet(predict(temp, dtrain)))
            assign(paste0("y.train_ts.XG_",f), quiet(predict(temp, dtest)))
          }
          if (toupper(model_name) == "RF") {
            message("RandomForest running...")
            params_rf <- modifyList(indiv.defaults$RF, param.indiv[["RF"]] %||% list())
            y <- factor(outcome[cur.train], levels = c(0,1))
            params_rf$x <- as.data.frame(xtrain)
            params_rf$y <- y
            temp <- quiet(do.call(randomForest, params_rf))
            assign(paste0("y.pred.RF_",f), quiet(predict(temp, newdata = xtest, type = "prob")[,2]))
            assign(paste0("y.train_tr.RF_",f), quiet(predict(temp, newdata=xtrain, type="prob")[,2]))
            assign(paste0("y.train_ts.RF_",f), quiet(predict(temp, newdata=xtest, type="prob")[,2]))
          }
          if (toupper(model_name) == "EN") {
            message("ElasticNet running...")
            en.params <- modifyList(indiv.defaults$EN, param.indiv[["EN"]] %||% list())
            train.params <- en.params[setdiff(names(en.params), "s")]
            temp <- do.call(cv.glmnet, c(list(x = as.matrix(xtrain), y = outcome[cur.train]), train.params))
            s.val <- param.indiv[["EN"]]$s %||% "lambda.1se"   # user-specified or default
            assign(paste0("y.pred.EN_",f), quiet(predict(temp, newx=xtest, s=s.val, type = "response")))
            assign(paste0("y.train_tr.EN_",f), quiet(predict(temp, newx=xtrain, s=s.val, type="response")))
            assign(paste0("y.train_ts.EN_",f), quiet(predict(temp, newx=xtest, s=s.val, type="response")))
          }
          if (toupper(model_name) == "ADB") {
            message("AdaBoost running...")
            params_adb <- modifyList(indiv.defaults$ADB, param.indiv[["ADB"]] %||% list())
            params_adb$formula <- NULL
            params_adb$x <- NULL
            params_adb$y <- NULL
            ada.sr <- ifelse(outcome == 0, -1, 1)
            args_list <- c(list(xtrain, ada.sr[cur.train]), params_adb)
            temp <- quiet(do.call(adaboost, args_list))
            assign(paste0("y.pred.ADB_",f), quiet(predict(temp, xtest, type = "prob")))
            assign(paste0("y.train_tr.ADB_",f), quiet(predict(temp, xtrain, type = "prob")))
            assign(paste0("y.train_ts.ADB_",f), quiet(predict(temp, xtest, type = "prob")))
          }
          if (toupper(model_name) == "MB") {
            message("MBoost running...")
            params_mb <- modifyList(indiv.defaults$MB, param.indiv[["MB"]] %||% list())
            assign(paste0("y.pred.MB_",f), quiet(do.call(
              binary_mb_predict,
              c(list(xtrain = xtrain, outcome = outcome[cur.train], xtest = xtest), params_mb)
            )))
            assign(paste0("y.train_tr.MB_",f), quiet(do.call(
              binary_mb_predict,
              c(list(xtrain = xtrain, outcome = outcome[cur.train]), params_mb)
            )))
            assign(paste0("y.train_ts.MB_",f), get(paste0("y.pred.MB_",f)))
          }
          if (toupper(model_name) == "GB") {
            message("GradientBoost running...")
            params_gb <- modifyList(indiv.defaults$GB, param.indiv[["GB"]] %||% list())
            train.df <- data.frame(y = outcome[cur.train], xtrain)
            params_gb$formula <- y ~.
            params_gb$data <- train.df
            params_gb$x <- NULL
            params_gb$y <- NULL
            temp <- quiet(do.call(gbm::gbm, params_gb))
            assign(paste0("y.pred.GB_",f), quiet(predict(temp, data.frame(xtest), n.trees = params_gb$n.trees, type = "response")))
            assign(paste0("y.train_tr.GB_",f), quiet(predict(temp, data.frame(xtrain), n.trees = params_gb$n.trees, type = "response")))
            assign(paste0("y.train_ts.GB_",f), quiet(predict(temp, data.frame(xtest), n.trees = params_gb$n.trees, type = "response")))
          }
          if (toupper(model_name) == "KNN") {
            message("KNN running...")
            train.df <- data.frame(y = factor(outcome[cur.train], levels = c(0,1)), xtrain)
            y.train.factor <- factor(outcome[cur.train], levels = c(0,1))
            test.df  <- data.frame(xtest)
            colnames(test.df) <- colnames(xtrain)  # ensure names match
            knn.params <- modifyList(indiv.defaults$KNN, param.indiv[["KNN"]] %||% list())
            knn.params$x <- as.matrix(xtrain)
            knn.params$y <- y.train.factor
            temp <- quiet(do.call(knn3, knn.params))
            assign(paste0("y.pred.KNN_",f), quiet(predict(temp, data.frame(xtest), type="prob")[,2]))
            assign(paste0("y.train_tr.KNN_",f), quiet(predict(temp, data.frame(xtrain),type = "prob")[,2]))
            assign(paste0("y.train_ts.KNN_",f), quiet(predict(temp, data.frame(xtest),type = "prob")[,2]))
          }
          if (toupper(model_name) == "SVM") {
            message("SVM running...")
            params_svm <- modifyList(indiv.defaults$SVM %||% list(), param.indiv[["SVM"]] %||% list())
            params_svm$formula <- NULL
            params_svm$data <- NULL
            params_svm$x <- NULL
            params_svm$y <- NULL
            train.df <- data.frame(y = factor(outcome[cur.train], levels = c(0, 1)), xtrain)
            temp <- quiet(do.call(svm, c(list(y ~ ., data = train.df, probability = TRUE), params_svm)))
            assign(paste0("y.pred.SVM_",f), svm_binary_probability(temp, xtest, positive = "1"))
            assign(paste0("y.train_tr.SVM_",f), svm_binary_probability(temp, xtrain, positive = "1"))
            assign(paste0("y.train_ts.SVM_",f), get(paste0("y.pred.SVM_",f)))
          }
          if (toupper(model_name) == "NB") {
            message("NB running...")
            params <- modifyList(indiv.defaults$NB %||% list(), param.indiv[["NB"]] %||% list())
            usekernel_flag <- params$usekernel %||% TRUE
            params$usekernel <- NULL
            y1<-as.factor(as.character(as.numeric(outcome[cur.train])))
            dat.nb<- cbind.data.frame("class" = y1, xtrain)
            params$formula <- class ~.
            params$data <- dat.nb
            params$x <- NULL
            params$y <- NULL
            temp <- quiet(do.call(naive_bayes, c(params, list(usekernel = usekernel_flag))))
            assign(paste0("y.pred.NB_",f), quiet(predict(temp, data.frame(xtest), type="prob")[,2]))
            assign(paste0("y.train_tr.NB_",f), quiet(predict(temp, data.frame(xtrain), type="prob")[,2]))
            assign(paste0("y.train_ts.NB_",f), quiet(predict(temp, data.frame(xtest), type="prob")[,2]))
          }
          if (toupper(model_name) == "DCT") {
            message("DecisionTree running...")
            y_outcome <- factor(outcome[cur.train], levels = c(0,1))
            train.df <- data.frame(y = y_outcome, xtrain)
            params_dct <- modifyList(indiv.defaults$DCT %||% list(), param.indiv[["DCT"]] %||% list())
            params_dct$x <- NULL
            params_dct$y <- NULL
            params_dct$formula <- y ~.
            params_dct$data <- train.df
            params_dct$method <- "class"
            temp <- quiet(do.call(rpart, params_dct))
            assign(paste0("y.pred.DCT_",f), quiet(predict(temp, newdata = data.frame(xtest), type = "prob")[,2]))
            assign(paste0("y.train_tr.DCT_",f), quiet(predict(temp, newdata=data.frame(xtrain), type="prob")[,2]))
            assign(paste0("y.train_ts.DCT_",f), quiet(predict(temp, newdata=data.frame(xtest), type="prob")[,2]))

          }
          if (toupper(model_name) == "NN") {
            # Drop zero-variance columns safely
            vars.ok <- apply(xtrain, 2, sd) > 0
            if(sum(vars.ok) == 0) stop("No variables left after dropping zero-variance columns!")
            xtrain <- xtrain[, vars.ok, drop = FALSE]
            xtest  <- xtest[, vars.ok, drop = FALSE]

            # Keep only common columns
            common <- intersect(colnames(xtrain), colnames(xtest))
            if(length(common) == 0) stop("No common features between train and test!")
            xtrain <- xtrain[, common, drop = FALSE]
            xtest  <- xtest[, common, drop = FALSE]

            # Scale
            m <- colMeans(xtrain)
            s <- apply(xtrain, 2, sd)
            xtrain.sd <- as.data.frame(scale(xtrain, center = m, scale = s))
            xtest.sd  <- as.data.frame(scale(xtest,  center = m, scale = s))

            # Prepare training data
            train_df <- cbind(y = as.numeric(outcome[cur.train]), xtrain.sd)
            params_nn <- param.indiv[["NN"]] %||% indiv.defaults$NN
            params_nn$formula <- y ~ .
            params_nn$data    <- train_df
            temp.nn <- quiet(do.call(neuralnet::neuralnet, params_nn))

            # Predict
            assign(paste0("y.pred.NN_", f), as.vector(predict(temp.nn, newdata = xtest.sd)))
            assign(paste0("y.train_tr.NN_", f), as.vector(predict(temp.nn, newdata = xtrain.sd)))
            assign(paste0("y.train_ts.NN_", f), as.vector(predict(temp.nn, newdata = xtest.sd)))

          }
        }
        ############## ensemble of 10 methods #############
        foo <- do.call(cbind, lapply(models.indiv, function(m) {
          get(paste0("y.train_tr.", m, "_", f))
        }))

        foo[foo < 0] <- 0
        foo[foo > 1] <- 1

        mat.train<- cbind(mat.train, foo)

        foo2 <- do.call(cbind, lapply(models.indiv, function(m) {
          get(paste0("y.train_ts.", m, "_", f))
        }))

        foo2[foo2 < 0] <- 0
        foo2[foo2 > 1] <- 1

        mat.test<- cbind(mat.test, foo2)

        ################# weighted avg ###############
        foo.w <- do.call(cbind, lapply(models.indiv, function(m) {
          get(paste0("y.train_tr.", m, "_", f))
        }))

        foo.w[foo.w < 0] <- 0
        foo.w[foo.w > 1] <- 1
        y.pred.train.all<- cbind(y.pred.train.all, foo.w)

        foo2.w <- do.call(cbind, lapply(models.indiv, function(m) {
          get(paste0("y.train_ts.", m, "_", f))
        }))

        foo2.w[foo2.w < 0] <- 0
        foo2.w[foo2.w > 1] <- 1

        y.pred.mat.all<- cbind(y.pred.mat.all, foo2.w)
      }

      colnames(mat.train)<- colnames(mat.test)<- colnames(y.pred.train.all)<- colnames(y.pred.mat.all) <- paste0("M",1:(length(models.indiv)*length(features)))

      for (model_name in models.ens) {
        mat.test.df <- data.frame(mat.test, check.names = FALSE)
        mat.train.df <- data.frame(mat.train, check.names = FALSE)
        colnames(mat.train.df) <- colnames(mat.train)
        colnames(mat.test.df) <- colnames(mat.train)
        if (toupper(model_name) == "XG") {
          message("XgBoost ensemble...")
          params_xg <- modifyList(ens.defaults$XG, param.ens[["XG"]] %||% list())
          params_xg$data  <- NULL
          params_xg$label <- NULL
          if (!is.null(params_xg$eta)) {
            params_xg$learning_rate <- params_xg$eta
            params_xg$eta <- NULL
          }
          dtrain.ens <- xgboost::xgb.DMatrix(data = as.matrix(mat.train), label = as.numeric(outcome[cur.train]))
          params_xg$data <- NULL
          params_xg$label <- NULL
          params_xg$verbose <- NULL
          params_xg$objective <- params_xg$objective %||% "binary:logistic"
          params_xg$verbosity <- params_xg$verbosity %||% 0
          nrounds <- params_xg$nrounds %||% 80
          params_xg$nrounds <- NULL
          temp.XG <- xgboost::xgb.train(params = params_xg, data = dtrain.ens, nrounds = nrounds)
          y.ens.XG[cur.test] = predict(temp.XG, newdata = xgboost::xgb.DMatrix(as.matrix(mat.test)))
        }
        if (toupper(model_name) == "RF") {
          message("RandomForest running...")
          params_rf <- modifyList(ens.defaults$RF, param.ens[["RF"]] %||% list())
          y <- factor(outcome[cur.train], levels = c(0,1))
          params_rf$x <- mat.train.df
          params_rf$y <- y
          temp.RF <- do.call(randomForest, params_rf)
          y.ens.RF[cur.test] = predict(temp.RF, newdata = mat.test.df)
        }
        if (toupper(model_name) == "EN") {
          message("ElasticNet running...")
          en.params <- modifyList(ens.defaults$EN, param.ens[["EN"]] %||% list())
          s.val <- param.ens[["EN"]]$s %||% "lambda.1se"   # user-specified or default
          train.params <- en.params[setdiff(names(en.params), "s")]
          temp.EN <- do.call(cv.glmnet, c(list(x = as.matrix(mat.train), y = outcome[cur.train]), train.params))
          y.ens.EN[cur.test] = predict(temp.EN, newx=mat.test, type = "response", s=s.val)
        }
        if (toupper(model_name) == "ADB") {
          message("AdaBoost running...")
          params_adb <- modifyList(ens.defaults$ADB, param.indiv[["ADB"]] %||% list())
          params_adb$formula <- NULL
          params_adb$x <- NULL
          params_adb$y <- NULL
          ada.sr <- ifelse(outcome[cur.train] == 0, -1, 1)
          args_list <- c(list(X = as.matrix(mat.train), y = ada.sr), params_adb)
          temp.ADB <- do.call(adaboost, args_list)
          y.ens.ADB[cur.test] = predict(temp.ADB, mat.test, type = "prob")
        }
        if (toupper(model_name) == "MB") {
          message("MBoost running...")
          params_mb <- modifyList(ens.defaults$MB, param.ens[["MB"]] %||% list())
          y.ens.MB[cur.test] <- do.call(
            binary_mb_predict,
            c(list(xtrain = mat.train.df, outcome = outcome[cur.train],
                   xtest = mat.test.df), params_mb)
          )
        }
        if (toupper(model_name) == "GB") {
          message("GradientBoost running...")
          params_gb <- modifyList(ens.defaults$GB, param.ens[["GB"]] %||% list())
          train.df <- data.frame(y = outcome[cur.train], mat.train)
          params_gb$formula <- y ~.
          params_gb$data <- train.df
          params_gb$x <- NULL
          params_gb$y <- NULL
          temp.GB <- do.call(gbm::gbm, params_gb)
          y.ens.GB[cur.test] = predict(temp.GB, mat.test.df, n.trees = params_gb$n.trees, type = "response")
        }
        if (toupper(model_name) == "KNN") {
          message("KNN running...")
          train.df <- data.frame(y = factor(outcome[cur.train], levels = c(0,1)), mat.train)
          y.train.factor <- factor(outcome[cur.train], levels = c(0,1))
          knn.params <- modifyList(ens.defaults$KNN, param.ens[["KNN"]] %||% list())
          knn.params$x <- mat.train
          knn.params$y <- y.train.factor
          temp.KNN <- do.call(knn3, knn.params)
          y.ens.KNN[cur.test] = predict(temp.KNN, mat.test.df, type = "prob")[,2]
        }
        if (toupper(model_name) == "SVM") {
          message("SVM running...")
          params_svm <- modifyList(ens.defaults$SVM %||% list(), param.ens[["SVM"]] %||% list())
          params_svm$formula <- NULL
          params_svm$data <- NULL
          params_svm$x <- NULL
          params_svm$y <- NULL
          train.df <- data.frame(y = factor(outcome[cur.train], levels = c(0, 1)), mat.train.df)
          temp.SVM <- do.call(svm, c(list(y ~ ., data = train.df, probability = TRUE), params_svm))
          y.ens.SVM[cur.test] <- svm_binary_probability(temp.SVM, mat.test.df, positive = "1")
        }
        if (toupper(model_name) == "NB") {
          message("NB running...")
          params_nb <- modifyList(ens.defaults$NB %||% list(), param.ens[["NB"]] %||% list())
          usekernel_flag <- params_nb$usekernel %||% TRUE
          params_nb$usekernel <- NULL
          y1<-as.factor(as.character(as.numeric(outcome[cur.train])))
          dat.nb<- cbind.data.frame("class" = y1, mat.train)
          params_nb$formula <- class ~.
          params_nb$data <- dat.nb
          params_nb$x <- NULL
          params_nb$y <- NULL
          temp.NB <- do.call(naive_bayes, c(params_nb, list(usekernel = usekernel_flag)))
          y.ens.NB[cur.test] <- predict(temp.NB, mat.test.df, threshold = 0.01, type="prob")[,2]
        }
        if (toupper(model_name) == "DCT") {
          message("DecisionTree running...")
          y_outcome <- factor(outcome[cur.train], levels = c(0,1))
          train.df <- data.frame(y = y_outcome, mat.train.df)
          params_dct <- modifyList(ens.defaults$DCT %||% list(), param.ens[["DCT"]] %||% list())
          params_dct$x <- NULL
          params_dct$y <- NULL
          params_dct$formula <- y ~.
          params_dct$data <- train.df
          params_dct$method <- "class"
          temp.DCT <- do.call(rpart, params_dct)
          y.ens.DCT[cur.test]=predict(temp.DCT, newdata=mat.test.df, type="matrix")[,6]
        }
        if (toupper(model_name) == "NN") {
          message("NeuralNet running...")
          params_nn <- modifyList(ens.defaults$NN, param.ens[["NN"]] %||% list())
          params_nn$x <- NULL
          params_nn$y <- NULL
          vars.ok <- apply(mat.train.df, 2, sd) > 0
          mat.train.nn <- mat.train.df[, vars.ok, drop = FALSE]
          m <- colMeans(mat.train.nn)
          s <- apply(mat.train.nn, 2, sd)
          mat.xtrain.sd <- as.data.frame(scale(mat.train.nn, center = m, scale = s))
          mat.train_nn <- cbind(y = as.numeric(outcome[cur.train]), mat.train.nn)
          params_nn$formula <- y ~ .
          params_nn$data <- mat.train_nn
          temp.NN <- do.call(neuralnet::neuralnet, params_nn)
          y.ens.NN[cur.test]=predict(temp.NN, newdata=mat.test, type="prob")
        }
      }

      ########### weighted avg ###############
      ind<- which(apply(y.pred.train.all,2,sd) == 0)
      if(length(ind) > 0){
        y.pred.train.all <- y.pred.train.all[,-ind, drop = FALSE]
        y.pred.mat.all <- y.pred.mat.all[,-ind, drop = FALSE]
      }

      if (ncol(y.pred.train.all) > 0) {
        y.pred.train.all.df <- as.data.frame(y.pred.train.all, check.names = FALSE)
        y.pred.mat.all.df <- as.data.frame(y.pred.mat.all, check.names = FALSE)
        colnames(y.pred.mat.all.df) <- colnames(y.pred.train.all.df)

        weighted_train_df <- data.frame(outcome = outcome[cur.train], y.pred.train.all.df,
                                        check.names = FALSE)
        temp.avg.all <- glm(outcome ~ ., family = binomial(link = "logit"),
                            data = weighted_train_df)

        y.pred.avg.all[cur.test] <- predict(temp.avg.all, newdata = y.pred.mat.all.df,
                                            type = "response")
      } else {
        y.pred.avg.all[cur.test] <- rowMeans(mat.test, na.rm = TRUE)
      }

      mat.test.full <- rbind(mat.test.full, mat.test)

    }
    ########### simple avg ###############

    mat.test.full[mat.test.full < 0] <- 0
    mat.test.full[mat.test.full > 1] <- 1

    y.pred.avg1.all<-apply(mat.test.full,1,mean)

    ############# pooling ensemble of methods ####################

    for (model_name in models.ens) {
      if (toupper(model_name) == "XG") {
        y.ens.XG.pool<- c(y.ens.XG.pool,y.ens.XG)
      }
      if (toupper(model_name) == "RF") {
        y.ens.RF.pool<- c(y.ens.RF.pool,y.ens.RF)
      }
      if (toupper(model_name) == "EN") {
        y.ens.EN.pool<- c(y.ens.EN.pool,y.ens.EN)
      }
      if (toupper(model_name) == "NN") {
        y.ens.NN.pool<- c(y.ens.NN.pool,y.ens.NN)
      }
      if (toupper(model_name) == "ADB") {
        y.ens.ADB.pool<- c(y.ens.ADB.pool,y.ens.ADB)
      }
      if (toupper(model_name) == "MB") {
        y.ens.MB.pool<- c(y.ens.MB.pool,y.ens.MB)
      }
      if (toupper(model_name) == "GB") {
        y.ens.GB.pool<- c(y.ens.GB.pool,y.ens.GB)
      }
      if (toupper(model_name) == "KNN") {
        y.ens.KNN.pool<- c(y.ens.KNN.pool,y.ens.KNN)
      }
      if (toupper(model_name) == "SVM") {
        y.ens.SVM.pool<- c(y.ens.SVM.pool,y.ens.SVM)
      }
      if (toupper(model_name) == "NB") {
        y.ens.NB.pool<- c(y.ens.NB.pool,y.ens.NB)
      }
      if (toupper(model_name) == "DCT") {
        y.ens.DCT.pool<- c(y.ens.DCT.pool,y.ens.DCT)
      }
    }

    ############# pooling simple and weighted avg of methods ####################
    y.pred.avg.pool.all<- c(y.pred.avg.pool.all,y.pred.avg.all)

    y.pred.avg.pool1.all<- c(y.pred.avg.pool1.all,y.pred.avg1.all)
  }

  ############# avg of iterations ####################

  for (model_name in models.ens) {
    if (toupper(model_name) == "XG") {
      y.ens.avg.pool.XG<-apply(matrix(y.ens.XG.pool,length(outcome),iter),1,function(x) mean(x, na.rm=T))
    }
    if (toupper(model_name) == "RF") {
      y.ens.avg.pool.RF<-apply(matrix(y.ens.RF.pool,length(outcome),iter),1,function(x) mean(x, na.rm=T))
    }
    if (toupper(model_name) == "EN") {
      y.ens.avg.pool.EN<-apply(matrix(y.ens.EN.pool,length(outcome),iter),1,function(x) mean(x, na.rm=T))
    }
    if (toupper(model_name) == "NN") {
      y.ens.avg.pool.NN<-apply(matrix(y.ens.NN.pool,length(outcome),iter),1,function(x) mean(x, na.rm=T))
    }
    if (toupper(model_name) == "ADB") {
      y.ens.avg.pool.ADB<-apply(matrix(y.ens.ADB.pool,length(outcome),iter),1,function(x) mean(x, na.rm=T))
    }
    if (toupper(model_name) == "MB") {
      y.ens.avg.pool.MB<-apply(matrix(y.ens.MB.pool,length(outcome),iter),1,function(x) mean(x, na.rm=T))
    }
    if (toupper(model_name) == "GB") {
      y.ens.avg.pool.GB<-apply(matrix(y.ens.GB.pool,length(outcome),iter),1,function(x) mean(x, na.rm=T))
    }
    if (toupper(model_name) == "KNN") {
      y.ens.avg.pool.KNN<-apply(matrix(y.ens.KNN.pool,length(outcome),iter),1,function(x) mean(x, na.rm=T))
    }
    if (toupper(model_name) == "SVM") {
      y.ens.avg.pool.SVM<-apply(matrix(y.ens.SVM.pool,length(outcome),iter),1,function(x) mean(x, na.rm=T))
    }
    if (toupper(model_name) == "NB") {
      y.ens.avg.pool.NB<-apply(matrix(y.ens.NB.pool,length(outcome),iter),1,function(x) mean(x, na.rm=T))
    }
    if (toupper(model_name) == "DCT") {
      y.ens.avg.pool.DCT<-apply(matrix(y.ens.DCT.pool,length(outcome),iter),1,function(x) mean(x, na.rm=T))
    }
  }

  y.pred.avg.pool.all.W<-apply(matrix(y.pred.avg.pool.all,length(outcome),iter),1,function(x) mean(x, na.rm=T))

  y.pred.avg.pool1.all.S<-apply(matrix(y.pred.avg.pool1.all,length(outcome),iter),1,function(x) mean(x, na.rm=T))

  outcome.pool <- outcome
  model.names <- c()
  auc.tab <- c()
  ci.low.tab <- c()
  ci.up.tab <- c()
  sp.98.tab <- c()
  sp.95.tab <- c()

  for (model_name in models.ens) {
    if (toupper(model_name) == "XG") {
      roc.XG.ens <- roc(outcome.pool, y.ens.avg.pool.XG, ci = TRUE, direction = "<")
      auc.tab <- c(auc.tab, roc.XG.ens$auc)
      ci.low.tab <- c(ci.low.tab, roc.XG.ens$ci[1])
      ci.up.tab <- c(ci.up.tab, roc.XG.ens$ci[3])
      sp.98.tab <- c(sp.98.tab, as.numeric(coords(roc.XG.ens,x=0.98,input="specificity",ret="sensitivity")))
      sp.95.tab <- c(sp.95.tab, as.numeric(coords(roc.XG.ens,x=0.95,input="specificity",ret="sensitivity")))
      model.names <- c(model.names, "XG.ens")
    }
    if (toupper(model_name) == "RF") {
      roc.RF.ens <- roc(outcome.pool, y.ens.avg.pool.RF, ci = TRUE, direction = "<")
      auc.tab <- c(auc.tab, roc.RF.ens$auc)
      ci.low.tab <- c(ci.low.tab, roc.RF.ens$ci[1])
      ci.up.tab <- c(ci.up.tab, roc.RF.ens$ci[3])
      sp.98.tab <- c(sp.98.tab, as.numeric(coords(roc.RF.ens,x=0.98,input="specificity",ret="sensitivity")))
      sp.95.tab <- c(sp.95.tab, as.numeric(coords(roc.RF.ens,x=0.95,input="specificity",ret="sensitivity")))
      model.names <- c(model.names, "RF.ens")
    }
    if (toupper(model_name) == "EN") {
      roc.EN.ens <- roc(outcome.pool, y.ens.avg.pool.EN, ci = TRUE, direction = "<")
      auc.tab <- c(auc.tab, roc.EN.ens$auc)
      ci.low.tab <- c(ci.low.tab, roc.EN.ens$ci[1])
      ci.up.tab <- c(ci.up.tab, roc.EN.ens$ci[3])
      sp.98.tab <- c(sp.98.tab, as.numeric(coords(roc.EN.ens,x=0.98,input="specificity",ret="sensitivity")))
      sp.95.tab <- c(sp.95.tab, as.numeric(coords(roc.EN.ens,x=0.95,input="specificity",ret="sensitivity")))
      model.names <- c(model.names, "EN.ens")
    }
    if (toupper(model_name) == "NN") {
      roc.NN.ens <- roc(outcome.pool, y.ens.avg.pool.NN, ci = TRUE, direction = "<")
      auc.tab <- c(auc.tab, roc.NN.ens$auc)
      ci.low.tab <- c(ci.low.tab, roc.NN.ens$ci[1])
      ci.up.tab <- c(ci.up.tab, roc.NN.ens$ci[3])
      sp.98.tab <- c(sp.98.tab, as.numeric(coords(roc.NN.ens,x=0.98,input="specificity",ret="sensitivity")))
      sp.95.tab <- c(sp.95.tab, as.numeric(coords(roc.NN.ens,x=0.95,input="specificity",ret="sensitivity")))
      model.names <- c(model.names, "NN.ens")
    }
    if (toupper(model_name) == "ADB") {
      roc.ADB.ens <- roc(outcome.pool, y.ens.avg.pool.ADB, ci = TRUE, direction = "<")
      auc.tab <- c(auc.tab, roc.ADB.ens$auc)
      ci.low.tab <- c(ci.low.tab, roc.ADB.ens$ci[1])
      ci.up.tab <- c(ci.up.tab, roc.ADB.ens$ci[3])
      sp.98.tab <- c(sp.98.tab, as.numeric(coords(roc.ADB.ens,x=0.98,input="specificity",ret="sensitivity")))
      sp.95.tab <- c(sp.95.tab, as.numeric(coords(roc.ADB.ens,x=0.95,input="specificity",ret="sensitivity")))
      model.names <- c(model.names, "ADB.ens")
    }
    if (toupper(model_name) == "MB") {
      roc.MB.ens <- roc(outcome.pool, y.ens.avg.pool.MB, ci = TRUE, direction = "<")
      auc.tab <- c(auc.tab, roc.MB.ens$auc)
      ci.low.tab <- c(ci.low.tab, roc.MB.ens$ci[1])
      ci.up.tab <- c(ci.up.tab, roc.MB.ens$ci[3])
      sp.98.tab <- c(sp.98.tab, as.numeric(coords(roc.MB.ens,x=0.98,input="specificity",ret="sensitivity")))
      sp.95.tab <- c(sp.95.tab, as.numeric(coords(roc.MB.ens,x=0.95,input="specificity",ret="sensitivity")))
      model.names <- c(model.names, "MB.ens")
    }
    if (toupper(model_name) == "GB") {
      roc.GB.ens <- roc(outcome.pool, y.ens.avg.pool.GB, ci = TRUE, direction = "<")
      auc.tab <- c(auc.tab, roc.GB.ens$auc)
      ci.low.tab <- c(ci.low.tab, roc.GB.ens$ci[1])
      ci.up.tab <- c(ci.up.tab, roc.GB.ens$ci[3])
      sp.98.tab <- c(sp.98.tab, as.numeric(coords(roc.GB.ens,x=0.98,input="specificity",ret="sensitivity")))
      sp.95.tab <- c(sp.95.tab, as.numeric(coords(roc.GB.ens,x=0.95,input="specificity",ret="sensitivity")))
      model.names <- c(model.names, "GB.ens")
    }
    if (toupper(model_name) == "KNN") {
      roc.KNN.ens <- roc(outcome.pool, y.ens.avg.pool.KNN, ci = TRUE, direction = "<")
      auc.tab <- c(auc.tab, roc.KNN.ens$auc)
      ci.low.tab <- c(ci.low.tab, roc.KNN.ens$ci[1])
      ci.up.tab <- c(ci.up.tab, roc.KNN.ens$ci[3])
      sp.98.tab <- c(sp.98.tab, as.numeric(coords(roc.KNN.ens,x=0.98,input="specificity",ret="sensitivity")))
      sp.95.tab <- c(sp.95.tab, as.numeric(coords(roc.KNN.ens,x=0.95,input="specificity",ret="sensitivity")))
      model.names <- c(model.names, "KNN.ens")
    }
    if (toupper(model_name) == "SVM") {
      roc.SVM.ens <- roc(outcome.pool, y.ens.avg.pool.SVM, ci = TRUE, direction = "<")
      auc.tab <- c(auc.tab, roc.SVM.ens$auc)
      ci.low.tab <- c(ci.low.tab, roc.SVM.ens$ci[1])
      ci.up.tab <- c(ci.up.tab, roc.SVM.ens$ci[3])
      sp.98.tab <- c(sp.98.tab, as.numeric(coords(roc.SVM.ens,x=0.98,input="specificity",ret="sensitivity")))
      sp.95.tab <- c(sp.95.tab, as.numeric(coords(roc.SVM.ens,x=0.95,input="specificity",ret="sensitivity")))
      model.names <- c(model.names, "SVM.ens")
    }
    if (toupper(model_name) == "NB") {
      roc.NB.ens <- roc(outcome.pool, y.ens.avg.pool.NB, ci = TRUE, direction = "<")
      auc.tab <- c(auc.tab, roc.NB.ens$auc)
      ci.low.tab <- c(ci.low.tab, roc.NB.ens$ci[1])
      ci.up.tab <- c(ci.up.tab, roc.NB.ens$ci[3])
      sp.98.tab <- c(sp.98.tab, as.numeric(coords(roc.NB.ens,x=0.98,input="specificity",ret="sensitivity")))
      sp.95.tab <- c(sp.95.tab, as.numeric(coords(roc.NB.ens,x=0.95,input="specificity",ret="sensitivity")))
      model.names <- c(model.names, "NB.ens")
    }
    if (toupper(model_name) == "DCT") {
      roc.DCT.ens <- roc(outcome.pool, y.ens.avg.pool.DCT, ci = TRUE, direction = "<")
      auc.tab <- c(auc.tab, roc.DCT.ens$auc)
      ci.low.tab <- c(ci.low.tab, roc.DCT.ens$ci[1])
      ci.up.tab <- c(ci.up.tab, roc.DCT.ens$ci[3])
      sp.98.tab <- c(sp.98.tab, as.numeric(coords(roc.DCT.ens,x=0.98,input="specificity",ret="sensitivity")))
      sp.95.tab <- c(sp.95.tab, as.numeric(coords(roc.DCT.ens,x=0.95,input="specificity",ret="sensitivity")))
      model.names <- c(model.names, "DCT.ens")
    }
  }

  roc.all <- roc(outcome.pool, y.pred.avg.pool.all.W, ci = TRUE, direction="<")
  auc.tab <- c(auc.tab, roc.all$auc)
  ci.low.tab <- c(ci.low.tab, roc.all$ci[1])
  ci.up.tab <- c(ci.up.tab, roc.all$ci[3])
  sp.98.tab <- c(sp.98.tab, as.numeric(coords(roc.all,x=0.98,input="specificity",ret="sensitivity")))
  sp.95.tab <- c(sp.95.tab, as.numeric(coords(roc.all,x=0.95,input="specificity",ret="sensitivity")))
  model.names <- c(model.names, "Weighted")

  roc.all1 <- roc(outcome.pool, y.pred.avg.pool1.all.S, ci = TRUE, direction="<")
  auc.tab <- c(auc.tab, roc.all1$auc)
  ci.low.tab <- c(ci.low.tab, roc.all1$ci[1])
  ci.up.tab <- c(ci.up.tab, roc.all1$ci[3])
  sp.98.tab <- c(sp.98.tab, as.numeric(coords(roc.all1,x=0.98,input="specificity",ret="sensitivity")))
  sp.95.tab <- c(sp.95.tab, as.numeric(coords(roc.all1,x=0.95,input="specificity",ret="sensitivity")))
  model.names <- c(model.names, "Simple")

  get_marker <- function(model_name) {
    if (identical(model_name, "Weighted")) return(y.pred.avg.pool.all.W)
    if (identical(model_name, "Simple")) return(y.pred.avg.pool1.all.S)
    base_name <- sub("\\.ens$", "", model_name)
    get(paste0("y.ens.avg.pool.", base_name), inherits = TRUE)
  }
  sp95_ci <- t(vapply(
    model.names,
    function(model_name) binary_sp_ci(outcome.pool, get_marker(model_name)),
    numeric(2)
  ))

  df <- data.frame(
    model = model.names,
    auc = auc.tab,
    ci.low = ci.low.tab,
    ci.up = ci.up.tab,
    sp95 = sp.95.tab,
    sp95.ci.low = sp95_ci[, "ci.low"],
    sp95.ci.up = sp95_ci[, "ci.up"],
    sp98 = sp.98.tab
  )

  return(df)
}

#=======Survival (Works) #####
surv_train_performance <- function(outcome, survdays, time_point, train_data, features, models.indiv, iter, num.folds = NULL, ensemble = TRUE, models.ens = NULL, param.indiv = NULL, param.ens = NULL, iter.offset = 0) {
  nfold <- max(2L, as.integer(num.folds %||% 3L))
  indiv.defaults <- list(
    XG = list(eta = 0.1, nrounds = 1000),
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

  ens.defaults <- list(
    XG = list(eta = 0.1, nrounds = 1000),
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

  if (ensemble == TRUE) {
    for (model_name in models.ens) {
      if (grepl("XG", model_name, ignore.case = TRUE)) {
        y.ens.XG.pool <- NULL
      }
      if (grepl("RF", model_name, ignore.case = TRUE)) {
        y.ens.RF.pool <- NULL
      }
      if (grepl("EN", model_name, ignore.case = TRUE)) {
        y.ens.EN.pool <- NULL
      }
      if (grepl("MB", model_name, ignore.case = TRUE)) {
        y.ens.MB.pool <- NULL
      }
      if (grepl("ADB", model_name, ignore.case = TRUE)) {
        y.ens.ADB.pool <- NULL
      }
      if (grepl("GB", model_name, ignore.case = TRUE)) {
        y.ens.GB.pool <- NULL
      }
      if (grepl("KNN", model_name, ignore.case = TRUE)) {
        y.ens.KNN.pool <- NULL
      }
      if (grepl("NN", model_name, ignore.case = TRUE)) {
        y.ens.NN.pool <- NULL
      }
      if (grepl("SVM", model_name, ignore.case = TRUE)) {
        y.ens.SVM.pool <- NULL
      }
      if (grepl("DCT", model_name, ignore.case = TRUE)) {
        y.ens.DCT.pool <- NULL
      }
    }
  }

  y.pred.avg.pool.all <- NULL
  y.pred.avg.pool1.all <- NULL

  for (j in 1:iter) {
    iter.id <- j + iter.offset
    message(iter.id)
    set.seed(iter.id)

    cv.fold <- balanced.folds.vec(outcome, nfold)
    if (ensemble == TRUE) {
      for (model_name in models.ens) {
        if (grepl("XG", model_name, ignore.case = TRUE)) {
          y.ens.XG=rep(NA, length(outcome))
        }
        if (grepl("RF", model_name, ignore.case = TRUE)) {
          y.ens.RF=rep(NA, length(outcome))
        }
        if (grepl("EN", model_name, ignore.case = TRUE)) {
          y.ens.EN=rep(NA, length(outcome))
        }
        if (grepl("MB", model_name, ignore.case = TRUE)) {
          y.ens.MB <- rep(NA, length(outcome))
        }
        if (grepl("ADB", model_name, ignore.case = TRUE)) {
          y.ens.ADB <- rep(NA, length(outcome))
        }
        if (grepl("GB", model_name, ignore.case = TRUE)) {
          y.ens.GB<- rep(NA, length(outcome))
        }
        if (grepl("KNN", model_name, ignore.case = TRUE)) {
          y.ens.KNN <- rep(NA, length(outcome))
        }
        if (grepl("NN", model_name, ignore.case = TRUE)) {
          y.ens.NN <- rep(NA, length(outcome))
        }
        if (grepl("SVM", model_name, ignore.case = TRUE)) {
          y.ens.SVM <- rep(NA, length(outcome))
        }
        if (grepl("DCT", model_name, ignore.case = TRUE)) {
          y.ens.DCT <- rep(NA, length(outcome))
        }
      }
    }

    y.pred.avg.all<- rep(NA, length(outcome))
    mat.test.full<- NULL

    for (i in seq_len(nfold)) {
      cur.train=which(cv.fold!=i)
      cur.test=which(cv.fold==i)

      mat.train<- NULL
      mat.test<- NULL
      y.pred.train.all<- NULL
      y.pred.mat.all<- NULL

      for (f in 1:(length(features))) {
        message(f)
        com<-unique(features[[f]])
        data_train <- train_data[, colnames(train_data) %in% com, drop = FALSE]
        colnames(data_train) <- make.names(colnames(data_train), unique = TRUE)

        xtrain<-(data_train)[cur.train,]
        xtest<-(data_train)[cur.test,]

        xtrain_time <- survdays[cur.train]
        xtrain_status <- outcome[cur.train]
        xtest_time <- survdays[cur.test]
        xtest_status <- outcome[cur.test]

        for (model_name in models.indiv) {
          if (grepl("XG", model_name, ignore.case = TRUE)) {
            assign(paste0("y.pred.XG_",f), rep(NA, length(cur.test)))
          }
          if (grepl("EN", model_name, ignore.case = TRUE)) {
            assign(paste0("y.pred.EN_",f), rep(NA, length(cur.test)))
          }
          if (grepl("RF", model_name, ignore.case = TRUE)) {
            assign(paste0("y.pred.RF_",f), rep(NA, length(cur.test)))
          }
          if (grepl("MB", model_name, ignore.case = TRUE)) {
            assign(paste0("y.pred.MB_",f), rep(NA, length(cur.test)))
          }
          if (grepl("ADB", model_name, ignore.case = TRUE)) {
            assign(paste0("y.pred.ADB_",f), rep(NA, length(cur.test)))
          }
          if (grepl("GB", model_name, ignore.case = TRUE)) {
            assign(paste0("y.pred.GB_",f), rep(NA, length(cur.test)))
          }
          if (grepl("KNN", model_name, ignore.case = TRUE)) {
            assign(paste0("y.pred.KNN_",f), rep(NA, length(cur.test)))
          }
          if (grepl("NN", model_name, ignore.case = TRUE)) {
            assign(paste0("y.pred.NN_",f), rep(NA, length(cur.test)))
          }
          if (grepl("SVM", model_name, ignore.case = TRUE)) {
            assign(paste0("y.pred.SVM_",f), rep(NA, length(cur.test)))
          }
          if (grepl("DCT", model_name, ignore.case = TRUE)) {
            assign(paste0("y.pred.DCT_",f), rep(NA, length(cur.test)))
          }
        }
        for (model_name in models.indiv) {
          if (grepl("XG", model_name, ignore.case = TRUE)) {
            assign(paste0("y.train_tr.XG_",f), rep(NA, length(cur.train)))
          }
          if (grepl("EN", model_name, ignore.case = TRUE)) {
            assign(paste0("y.train_tr.EN_",f), rep(NA, length(cur.train)))
          }
          if (grepl("RF", model_name, ignore.case = TRUE)) {
            assign(paste0("y.train_tr.RF_",f), rep(NA, length(cur.train)))
          }
          if (grepl("MB", model_name, ignore.case = TRUE)) {
            assign(paste0("y.train_tr.MB_",f), rep(NA, length(cur.train)))
          }
          if (grepl("ADB", model_name, ignore.case = TRUE)) {
            assign(paste0("y.train_tr.ADB_",f), rep(NA, length(cur.train)))
          }
          if (grepl("GB", model_name, ignore.case = TRUE)) {
            assign(paste0("y.train_tr.GB_",f), rep(NA, length(cur.train)))
          }
          if (grepl("KNN", model_name, ignore.case = TRUE)) {
            assign(paste0("y.train_tr.KNN_",f), rep(NA, length(cur.train)))
          }
          if (grepl("NN", model_name, ignore.case = TRUE)) {
            assign(paste0("y.train_tr.NN_",f), rep(NA, length(cur.train)))
          }
          if (grepl("SVM", model_name, ignore.case = TRUE)) {
            assign(paste0("y.train_tr.SVM_",f), rep(NA, length(cur.train)))
          }
          if (grepl("DCT", model_name, ignore.case = TRUE)) {
            assign(paste0("y.train_tr.DCT_",f), rep(NA, length(cur.train)))
          }
        }
        for (model_name in models.indiv) {
          if (grepl("XG", model_name, ignore.case = TRUE)) {
            assign(paste0("y.train_ts.XG_",f), rep(NA, length(cur.test)))
          }
          if (grepl("EN", model_name, ignore.case = TRUE)) {
            assign(paste0("y.train_ts.EN_",f), rep(NA, length(cur.test)))
          }
          if (grepl("RF", model_name, ignore.case = TRUE)) {
            assign(paste0("y.train_ts.RF_",f), rep(NA, length(cur.test)))
          }
          if (grepl("MB", model_name, ignore.case = TRUE)) {
            assign(paste0("y.train_ts.MB_",f), rep(NA, length(cur.test)))
          }
          if (grepl("ADB", model_name, ignore.case = TRUE)) {
            assign(paste0("y.train_ts.ADB_",f), rep(NA, length(cur.test)))
          }
          if (grepl("GB", model_name, ignore.case = TRUE)) {
            assign(paste0("y.train_ts.GB_",f), rep(NA, length(cur.test)))
          }
          if (grepl("KNN", model_name, ignore.case = TRUE)) {
            assign(paste0("y.train_ts.KNN_",f), rep(NA, length(cur.test)))
          }
          if (grepl("NN", model_name, ignore.case = TRUE)) {
            assign(paste0("y.train_ts.NN_",f), rep(NA, length(cur.test)))
          }
          if (grepl("SVM", model_name, ignore.case = TRUE)) {
            assign(paste0("y.train_ts.SVM_",f), rep(NA, length(cur.test)))
          }
          if (grepl("DCT", model_name, ignore.case = TRUE)) {
            assign(paste0("y.train_ts.DCT_",f), rep(NA, length(cur.test)))
          }
        }

        for (model_name in models.indiv) {
          if (grepl("XG", model_name, ignore.case = TRUE)) {
            message("XgBoost running...")
            params_user <- modifyList(indiv.defaults$XG, param.indiv[["XG"]] %||% list())
            params_user$nrounds <- NULL

            label_train <- as.numeric(ifelse(xtrain_status == 1, xtrain_time, -xtrain_time))

            val_ind <- sample(1:nrow(xtrain), 0.1 * nrow(xtrain))

            dtrain <- xgboost::xgb.DMatrix(
              data = as.matrix(xtrain[-val_ind, ]),
              label = label_train[-val_ind]
            )

            dval <- xgboost::xgb.DMatrix(
              data = as.matrix(xtrain[val_ind, ]),
              label = label_train[val_ind]
            )

            final_params <- modifyList(params_user, list(
              objective = "survival:cox",
              eval_metric = "cox-nloglik",
              disable_default_eval_metric = 1
            ))

            nrounds_val <- param.indiv[["XG"]]$nrounds %||% 1000

            temp <- xgboost::xgb.train(
              params = final_params,
              data = dtrain,
              nrounds = nrounds_val,
              evals = list(val = dval),
              verbose = 0
            )

            assign(paste0("y.pred.XG_",f), quiet(predict(temp, newdata = xgboost::xgb.DMatrix(as.matrix(xtest)))))
            assign(paste0("y.train_tr.XG_",f), quiet(predict(temp, newdata = xgboost::xgb.DMatrix(as.matrix(xtrain)))))
            assign(paste0("y.train_ts.XG_",f), quiet(predict(temp, newdata = xgboost::xgb.DMatrix(as.matrix(xtest)))))
          }
          if (grepl("RF", model_name, ignore.case = TRUE)) {
            message("RandomForest running...")
            params_rf <- modifyList(indiv.defaults$RF, param.indiv[["RF"]] %||% list())
            surv_train_data <- data.frame(time = xtrain_time, status = xtrain_status, xtrain)
            surv_test_data <- data.frame(time = xtest_time, status = xtest_status, xtest)
            formula_args <- list(formula = Surv(time, status) ~., data = surv_train_data)
            temp <- quiet(do.call(randomForestSRC::rfsrc, c(formula_args, params_rf)))
            predictions_train <- quiet(predict(temp, newdata = surv_train_data, type = "response"))
            predictions_test <- quiet(predict(temp, newdata = data.frame(xtest), type = "response"))

            assign(paste0("y.pred.RF_",f), predictions_test$predicted)
            assign(paste0("y.train_tr.RF_",f), predictions_train$predicted)
            assign(paste0("y.train_ts.RF_",f), predictions_test$predicted)

          }
          if (grepl("EN", model_name, ignore.case = TRUE)) {
            message("ElasticNet running...")
            params_en <- modifyList(indiv.defaults$EN, param.indiv[["EN"]] %||% list())
            surv_train_object <- Surv(time = xtrain_time, event = xtrain_status)
            s.val <- params_en$s %||% "lambda.min"
            train.params <- params_en[setdiff(names(params_en), "s")]
            cv_glmnet <- quiet(do.call(glmnet::cv.glmnet, c(
              list(x = as.matrix(xtrain), y = surv_train_object, family = "cox", nfolds = 10), train.params)))
            risk_scores_train <- quiet(predict(cv_glmnet, newx = as.matrix(xtrain), s = s.val, type = "response"))
            risk_scores_test <- quiet(predict(cv_glmnet, newx = as.matrix(xtest), s = s.val, type = "response"))

            assign(paste0("y.pred.EN_",f), risk_scores_test)
            assign(paste0("y.train_tr.EN_",f), risk_scores_train)
            assign(paste0("y.train_ts.EN_",f), risk_scores_test)
          }
          if (grepl("MB", model_name, ignore.case = TRUE)) {
            message("MBoost running...")
            params_mb <- modifyList(indiv.defaults$MB, param.indiv[["MB"]] %||% list())
            mstop_val <- params_mb$mstop %||% 300
            params_mb$mstop <- NULL

            control_object <- mboost::boost_control(mstop = mstop_val)

            x_matrix <- as.matrix(xtrain)
            y_surv <- Surv(time = xtrain_time, event = xtrain_status)

            fixed_args <- list(
              x = x_matrix,
              y = y_surv,
              family = mboost::CoxPH(), # Fixed survival family
              control = control_object
            )

            temp <- quiet(do.call(mboost::glmboost, c(fixed_args, params_mb)))
            risk_scores_train <- quiet(predict(temp, newdata = as.matrix(xtrain), type = "link"))
            risk_scores_test <- quiet(predict(temp, newdata = as.matrix(xtest), type = "link"))

            assign(paste0("y.pred.MB_",f), risk_scores_test)
            assign(paste0("y.train_tr.MB_",f), risk_scores_train)
            assign(paste0("y.train_ts.MB_",f), risk_scores_test)
          }
          if (grepl("ADB", model_name, ignore.case = TRUE)) {
            message("AdaBoost survival running...")
            params_adb <- modifyList(indiv.defaults$ADB, param.indiv[["ADB"]] %||% list())
            risk_scores_train <- do.call(
              survival_adb_predict,
              c(list(xtrain = xtrain, time = xtrain_time, status = xtrain_status,
                     time_point = time_point), params_adb)
            )
            risk_scores_test <- do.call(
              survival_adb_predict,
              c(list(xtrain = xtrain, time = xtrain_time, status = xtrain_status,
                     time_point = time_point, xtest = xtest), params_adb)
            )
            assign(paste0("y.pred.ADB_",f), risk_scores_test)
            assign(paste0("y.train_tr.ADB_",f), risk_scores_train)
            assign(paste0("y.train_ts.ADB_",f), risk_scores_test)
          }
          if (grepl("GB", model_name, ignore.case = TRUE)) {
            message("GradientBoost running...")
            params_gb <- modifyList(indiv.defaults$GB, param.indiv[["GB"]] %||% list())
            surv_train_data <- data.frame(time = xtrain_time, status = xtrain_status, xtrain)
            surv_test_data <- data.frame(time = xtest_time, status = xtest_status, xtest)
            fixed_args <- list(
              formula = Surv(time, status) ~.,
              data = surv_train_data,
              distribution = "coxph" # Fixed survival distribution
            )
            n.trees_pred <- params_gb$n.trees %||% 500
            temp <- quiet(do.call(gbm::gbm, c(fixed_args, params_gb)))

            risk_scores_train <- quiet(predict(temp, newdata = surv_train_data, n.trees = n.trees_pred, type = "link"))
            risk_scores_test <- quiet(predict(temp, newdata = surv_test_data, n.trees = n.trees_pred, type = "link"))
            assign(paste0("y.pred.GB_",f), risk_scores_test)
            assign(paste0("y.train_tr.GB_",f), risk_scores_train)
            assign(paste0("y.train_ts.GB_",f), risk_scores_test)
          }
          if (grepl("KNN", model_name, ignore.case = TRUE)) {
            message("KNN survival running...")
            params_knn <- modifyList(indiv.defaults$KNN, param.indiv[["KNN"]] %||% list())
            k_val <- params_knn$k %||% 5
            risk_scores_train <- survival_knn_predict(
              xtrain, xtrain_time, xtrain_status, k = k_val, leave_one_out = TRUE
            )
            risk_scores_test <- survival_knn_predict(
              xtrain, xtrain_time, xtrain_status, xtest, k = k_val
            )
            assign(paste0("y.pred.KNN_",f), risk_scores_test)
            assign(paste0("y.train_tr.KNN_",f), risk_scores_train)
            assign(paste0("y.train_ts.KNN_",f), risk_scores_test)
          }
          if (grepl("NN", model_name, ignore.case = TRUE)) {
            message("NN survival running...")
            params_nn <- modifyList(indiv.defaults$NN, param.indiv[["NN"]] %||% list())
            risk_scores_train <- do.call(
              survival_nn_predict,
              c(list(xtrain = xtrain, time = xtrain_time, status = xtrain_status), params_nn)
            )
            risk_scores_test <- do.call(
              survival_nn_predict,
              c(list(xtrain = xtrain, time = xtrain_time, status = xtrain_status,
                     xtest = xtest), params_nn)
            )
            assign(paste0("y.pred.NN_",f), risk_scores_test)
            assign(paste0("y.train_tr.NN_",f), risk_scores_train)
            assign(paste0("y.train_ts.NN_",f), risk_scores_test)
          }
          if (grepl("SVM", model_name, ignore.case = TRUE)) {
            message("SVM running...")
            params_svm <- modifyList(indiv.defaults$SVM, param.indiv[["SVM"]] %||% list())
            surv_train_data <- data.frame(time = xtrain_time, status = xtrain_status, xtrain)
            formula_args <- list(
              formula = Surv(time, status) ~.,
              data = surv_train_data,
              type = "regression"
            )
            temp <- quiet(do.call(survivalsvm::survivalsvm, c(formula_args, params_svm)))
            risk_scores_train <- quiet(predict(object = temp, newdata = data.frame(xtrain)))
            risk_scores_test <- quiet(predict(object = temp, newdata = data.frame(xtest)))

            assign(paste0("y.pred.SVM_",f), matrix(risk_scores_test$predicted, ncol = 1))
            assign(paste0("y.train_tr.SVM_",f), matrix(risk_scores_train$predicted, ncol = 1))
            assign(paste0("y.train_ts.SVM_",f), matrix(risk_scores_test$predicted, ncol = 1))
          }
          if (grepl("DCT", model_name, ignore.case = TRUE)) {
            message("DecisionTree running...")
            params_dct <- modifyList(indiv.defaults$DCT, param.indiv[["DCT"]] %||% list())
            surv_train_data <- data.frame(time = xtrain_time, status = xtrain_status, xtrain)
            formula_args <- list(
              formula = Surv(time, status) ~ .,
              data = surv_train_data,
              method = "exp"
            )

            # The user/default parameters (minsplit, maxdepth, cp) go into rpart.control
            control_params <- params_dct[names(params_dct)]

            temp <- quiet(do.call(rpart::rpart, c(formula_args, list(control = do.call(rpart::rpart.control, control_params)))))

            risk_scores_train <- quiet(predict(temp, newdata = data.frame(xtrain)))
            risk_scores_test <- quiet(predict(temp, newdata = data.frame(xtest)))

            assign(paste0("y.pred.DCT_",f), risk_scores_test)
            assign(paste0("y.train_tr.DCT_",f), risk_scores_train)
            assign(paste0("y.train_ts.DCT_",f), risk_scores_test)
          }
        }
        ############## ensemble #############
        foo <- do.call(cbind, lapply(models.indiv, function(m) {
          get(paste0("y.train_tr.", m, "_", f))
        }))

        foo[foo < 0] <- 0
        foo[foo > 1] <- 1

        mat.train<- cbind(mat.train, foo)

        foo2 <- do.call(cbind, lapply(models.indiv, function(m) {
          get(paste0("y.train_ts.", m, "_", f))
        }))

        foo2[foo2 < 0] <- 0
        foo2[foo2 > 1] <- 1

        mat.test<- cbind(mat.test, foo2)



        ################# weighted avg ###############
        foo.w <- do.call(cbind, lapply(models.indiv, function(m) {
          get(paste0("y.train_tr.", m, "_", f))
        }))

        foo.w[foo.w < 0] <- 0
        foo.w[foo.w > 1] <- 1
        y.pred.train.all<- cbind(y.pred.train.all, foo.w)

        foo2.w <- do.call(cbind, lapply(models.indiv, function(m) {
          get(paste0("y.train_ts.", m, "_", f))
        }))

        foo2.w[foo2.w < 0] <- 0
        foo2.w[foo2.w > 1] <- 1

        y.pred.mat.all<- cbind(y.pred.mat.all, foo2.w)

      }

      colnames(mat.test)<- colnames(mat.train)<-colnames(y.pred.mat.all)<- colnames(y.pred.train.all)<- paste0("M",1:(length(models.indiv)*length(features)))

      K <- length(models.indiv) * length(features)
      ensemble_colnames <- c("time", "status",paste0("M", seq_len(K)))
      surv_ensemble_train_data <- data.frame(time   = xtrain_time, status = xtrain_status, mat.train, check.names = FALSE)
      surv_ensemble_test_data <- data.frame(time   = xtest_time,status = xtest_status, mat.test, check.names = FALSE)
      colnames(surv_ensemble_train_data) <- ensemble_colnames
      colnames(surv_ensemble_test_data)  <- ensemble_colnames

      for (model_name in models.ens) {
        if (grepl("XG", model_name, ignore.case = TRUE)) {
          message("Xg ensemble...")
          params_user <- modifyList(list(objective = "survival:cox", eval_metric = "cox-nloglik", eta = 0.1), param.ens[[model_name]] %||% list())
          val_ind <- sample(1:nrow(mat.train), 0.1 * nrow(surv_ensemble_train_data))
          label_train <- ifelse(xtrain_status == 1, xtrain_time, -xtrain_time)
          label_train_test <- label_train[-val_ind]

          xtrain_xgboost <- as.matrix(mat.train[-val_ind,])
          val_data <- as.matrix(mat.train[val_ind,])
          x_val <- xgb.DMatrix(data = val_data, label = label_train[val_ind])

          nrounds_val <- params_user$nrounds %||% 1000
          params_user$nrounds <- NULL

          dtrain <- xgboost::xgb.DMatrix(
            data = as.matrix(mat.train[-val_ind, ]),
            label = as.numeric(label_train[-val_ind])
          )

          dval <- xgboost::xgb.DMatrix(
            data = as.matrix(mat.train[val_ind, ]),
            label = as.numeric(label_train[val_ind])
          )

          params_user$nrounds <- NULL

          params <- modifyList(params_user, list(
            objective = "survival:cox",
            eval_metric = "cox-nloglik"
          ))

          xgboost_fit <- xgboost::xgb.train(
            params = params,
            data = dtrain,
            nrounds = nrounds_val,
            watchlist = list(val = dval),
            verbose = 0
          )

          # prediction (no type needed)
          y.ens.XG[cur.test] <- predict(xgboost_fit, as.matrix(mat.test))
        }
        if (grepl("RF", model_name, ignore.case = TRUE)) {
          message("RF ensemble...")
          params_rf <- modifyList(list(ntree = 2000, max.depth = 10, importance = TRUE),param.ens[[model_name]] %||% list())
          temp <- quiet(do.call(rfsrc, c(list(formula = Surv(time,status) ~ ., data = surv_ensemble_train_data), params_rf)))
          predictions_test <- quiet(predict(temp, newdata = surv_ensemble_test_data, type = "response"))
          y.ens.RF[cur.test] <- predictions_test$predicted
        }
        if (grepl("EN", model_name, ignore.case = TRUE)) {
          message("EN ensemble...")
          params_en <- modifyList(list(alpha = 0.1), param.ens[[model_name]] %||% list())
          surv_train_object <- Surv(time = xtrain_time, event = xtrain_status)
          cv_glmnet <- do.call(cv.glmnet, c(list(x = as.matrix(mat.train), y = surv_train_object, family = "cox", nfolds = 10), params_en))
          s.val <- param.ens[["EN"]]$s %||% "lambda.1se"   # user-specified or default
          lambda_val <- switch(
            s.val,
            lambda.min = cv_glmnet$lambda.min,
            lambda.1se = cv_glmnet$lambda.1se,
            as.numeric(s.val)   # allows numeric override
          )
          temp <- quiet(glmnet(as.matrix(mat.train), surv_train_object, family = "cox", lambda = lambda_val, alpha = params_en$alpha %||% 0.1))
          risk_scores_test <- quiet(predict(temp, newx = as.matrix(mat.test), s = lambda_val, type = "response"))
          y.ens.EN[cur.test] <- risk_scores_test
        }
        if (grepl("MB", model_name, ignore.case = TRUE)) {
          message("MBoost ensemble...")
          params_mb <- modifyList(list(mstop = 300), param.ens[[model_name]] %||% list())
          control_object <- boost_control(mstop = params_mb$mstop %||% 300)
          temp <- quiet(glmboost(Surv(time, status) ~ ., data = surv_ensemble_train_data, family = CoxPH(), control = control_object))
          risk_scores_test <- quiet(predict(temp, newdata = surv_ensemble_test_data, type = "link"))
          y.ens.MB[cur.test] <- risk_scores_test
        }
        if (grepl("ADB", model_name, ignore.case = TRUE)) {
          message("ADB survival ensemble...")
          params_adb <- modifyList(ens.defaults$ADB, param.ens[[model_name]] %||% list())
          y.ens.ADB[cur.test] <- do.call(
            survival_adb_predict,
            c(list(xtrain = mat.train, time = xtrain_time, status = xtrain_status,
                   time_point = time_point, xtest = mat.test), params_adb)
          )
        }
        if (grepl("GB", model_name, ignore.case = TRUE)) {
          message("GB ensemble...")
          params_gb <- modifyList(list(n.trees = 1000, interaction.depth = 10, shrinkage = 0.005, bag.fraction = 0.5), param.ens[[model_name]] %||% list())
          temp <- quiet(do.call(gbm, c(list(formula = Surv(time, status) ~ ., data = surv_ensemble_train_data, distribution = "coxph"), params_gb)))
          y.ens.GB[cur.test] <- quiet(predict(temp, newdata = surv_ensemble_test_data, n.trees = params_gb$n.trees %||% 1000, type = "link"))
        }
        if (grepl("KNN", model_name, ignore.case = TRUE)) {
          message("KNN ensemble...")
          params_knn <- modifyList(list(k = 5), param.ens[[model_name]] %||% list())
          y.ens.KNN[cur.test] <- survival_knn_predict(
            mat.train, xtrain_time, xtrain_status, mat.test, k = params_knn$k %||% 5
          )
        }
        if (grepl("NN", model_name, ignore.case = TRUE)) {
          message("NN ensemble...")
          params_nn <- modifyList(ens.defaults$NN, param.ens[[model_name]] %||% list())
          y.ens.NN[cur.test] <- do.call(
            survival_nn_predict,
            c(list(xtrain = mat.train, time = xtrain_time, status = xtrain_status,
                   xtest = mat.test), params_nn)
          )
        }
        if (grepl("SVM", model_name, ignore.case = TRUE)) {
          message("SVM ensemble...")
          params_svm <- modifyList(list(gamma.mu = 0.5, opt.meth = "quadprog", kernel = "add_kernel"), param.ens[[model_name]] %||% list())
          temp <- quiet(do.call(survivalsvm, c(list(formula = Surv(time, status) ~ ., data = surv_ensemble_train_data, type = "regression"), params_svm)))
          risk_scores_test <- quiet(predict(temp, newdata = surv_ensemble_test_data))
          y.ens.SVM[cur.test] <- as.numeric(risk_scores_test$predicted)
        }
        if (grepl("DCT", model_name, ignore.case = TRUE)) {
          message("DCT ensemble...")
          params_dct <- modifyList(list(minsplit = 20, maxdepth = 10, cp = 0.01), param.ens[[model_name]] %||% list())
          temp <- quiet(rpart(Surv(time, status) ~ ., data = surv_ensemble_train_data, method = "exp", control = do.call(rpart.control, params_dct)))
          risk_scores_test <- quiet(predict(temp, newdata = surv_ensemble_test_data))
          y.ens.DCT[cur.test] <- unname(risk_scores_test)
        }
      }

      ########### weighted avg ###############
      message("Weighted average (Cox stacking)...")

      # --------------------------------------------------
      # 1. Stack base-model predictions (your naming)
      # --------------------------------------------------
      y.pred.train.all <- as.data.frame(y.pred.train.all)
      y.pred.mat.all   <- as.data.frame(y.pred.mat.all)

      time.train   <- survdays[cur.train]
      status.train <- outcome[cur.train]

      # --------------------------------------------------
      # 2. Remove NA / Inf rows
      # --------------------------------------------------
      valid_rows <- apply(
        y.pred.train.all,
        1,
        function(x) all(is.finite(x))
      )

      y.pred.train.all <- y.pred.train.all[valid_rows, , drop = FALSE]

      time.train   <- time.train[valid_rows]
      status.train <- status.train[valid_rows]

      # --------------------------------------------------
      # 3. Remove zero-variance predictors
      # --------------------------------------------------
      zero_var <- sapply(y.pred.train.all, function(x)
        var(x, na.rm = TRUE) < 1e-12
      )

      if (any(zero_var)) {
        y.pred.train.all <- y.pred.train.all[, !zero_var, drop = FALSE]
        y.pred.mat.all   <- y.pred.mat.all[,  !zero_var, drop = FALSE]
      }

      # --------------------------------------------------
      # 4. Standardize predictors safely
      # --------------------------------------------------
      mu  <- colMeans(y.pred.train.all)
      sdv <- apply(y.pred.train.all, 2, sd)

      sdv[sdv < 1e-8] <- 1

      y.pred.train.all <- scale(y.pred.train.all, center = mu, scale = sdv)
      y.pred.mat.all   <- scale(y.pred.mat.all,   center = mu, scale = sdv)

      y.pred.train.all <- as.data.frame(y.pred.train.all)
      y.pred.mat.all   <- as.data.frame(y.pred.mat.all)

      # --------------------------------------------------
      # 5. Add survival outcome
      # --------------------------------------------------
      y.pred.train.all$time   <- time.train
      y.pred.train.all$status <- status.train

      # --------------------------------------------------
      # 6. Clean again after scaling
      # --------------------------------------------------
      good_cols <- sapply(y.pred.train.all, function(x)
        all(is.finite(x)) && sd(x) > 1e-8
      )

      good_cols[c("time","status")] <- TRUE

      y.pred.train.all <- y.pred.train.all[, good_cols, drop = FALSE]

      predictors <- setdiff(colnames(y.pred.train.all), c("time","status"))

      # align test data EXACTLY
      y.pred.mat.all <- y.pred.mat.all[, predictors, drop = FALSE]

      # --------------------------------------------------
      # 7. Cox model
      # --------------------------------------------------
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

      # --------------------------------------------------
      # 8. Prediction
      # --------------------------------------------------
      y.pred.avg <- quiet(as.numeric(
        predict(
          temp.avg.all,
          newdata = y.pred.mat.all,
          type = "lp"
        )
      ))

      y.pred.avg.all[cur.test] <- y.pred.avg

      mat.test.full <- rbind(mat.test.full, mat.test)
    }
    ########### simple avg ###############

    mat.test.full[mat.test.full < 0] <- 0
    mat.test.full[mat.test.full > 1] <- 1

    y.pred.avg1.all<-apply(mat.test.full,1,mean)

    ############# pooling ensemble of methods ####################

    for (model_name in models.ens) {
      if (grepl("XG", model_name, ignore.case = TRUE)) {
        y.ens.XG.pool<- c(y.ens.XG.pool,y.ens.XG)
      }
      if (grepl("RF", model_name, ignore.case = TRUE)) {
        y.ens.RF.pool<- c(y.ens.RF.pool,y.ens.RF)
      }
      if (grepl("EN", model_name, ignore.case = TRUE)) {
        y.ens.EN.pool<- c(y.ens.EN.pool,y.ens.EN)
      }
      if (grepl("MB", model_name, ignore.case = TRUE)) {
        y.ens.MB.pool<- c(y.ens.MB.pool,y.ens.MB)
      }
      if (grepl("ADB", model_name, ignore.case = TRUE)) {
        y.ens.ADB.pool<- c(y.ens.ADB.pool,y.ens.ADB)
      }
      if (grepl("GB", model_name, ignore.case = TRUE)) {
        y.ens.GB.pool<- c(y.ens.GB.pool,y.ens.GB)
      }
      if (grepl("KNN", model_name, ignore.case = TRUE)) {
        y.ens.KNN.pool<- c(y.ens.KNN.pool,y.ens.KNN)
      }
      if (grepl("NN", model_name, ignore.case = TRUE)) {
        y.ens.NN.pool<- c(y.ens.NN.pool,y.ens.NN)
      }
      if (grepl("SVM", model_name, ignore.case = TRUE)) {
        y.ens.SVM.pool<- c(y.ens.SVM.pool,y.ens.SVM)
      }
      if (grepl("DCT", model_name, ignore.case = TRUE)) {
        y.ens.DCT.pool<- c(y.ens.DCT.pool,y.ens.DCT)
      }
    }

    ############# pooling simple and weighted avg of methods ####################
    y.pred.avg.pool.all<- c(y.pred.avg.pool.all,y.pred.avg.all)

    y.pred.avg.pool1.all<- c(y.pred.avg.pool1.all,y.pred.avg1.all)

  }

  ############# avg of iterations ####################

  for (model_name in models.ens) {
    if (grepl("XG", model_name, ignore.case = TRUE)) {
      y.ens.avg.pool.XG<-apply(matrix(y.ens.XG.pool,length(outcome),iter),1,function(x) mean(x, na.rm=T))
    }
    if (grepl("RF", model_name, ignore.case = TRUE)) {
      y.ens.avg.pool.RF<-apply(matrix(y.ens.RF.pool,length(outcome),iter),1,function(x) mean(x, na.rm=T))
    }
    if (grepl("EN", model_name, ignore.case = TRUE)) {
      y.ens.avg.pool.EN<-apply(matrix(y.ens.EN.pool,length(outcome),iter),1,function(x) mean(x, na.rm=T))
    }
    if (grepl("MB", model_name, ignore.case = TRUE)) {
      y.ens.avg.pool.MB<-apply(matrix(y.ens.MB.pool,length(outcome),iter),1,function(x) mean(x, na.rm=T))
    }
    if (grepl("ADB", model_name, ignore.case = TRUE)) {
      y.ens.avg.pool.ADB<-apply(matrix(y.ens.ADB.pool,length(outcome),iter),1,function(x) mean(x, na.rm=T))
    }
    if (grepl("GB", model_name, ignore.case = TRUE)) {
      y.ens.avg.pool.GB<-apply(matrix(y.ens.GB.pool,length(outcome),iter),1,function(x) mean(x, na.rm=T))
    }
    if (grepl("KNN", model_name, ignore.case = TRUE)) {
      y.ens.avg.pool.KNN<-apply(matrix(y.ens.KNN.pool,length(outcome),iter),1,function(x) mean(x, na.rm=T))
    }
    if (grepl("NN", model_name, ignore.case = TRUE)) {
      y.ens.avg.pool.NN<-apply(matrix(y.ens.NN.pool,length(outcome),iter),1,function(x) mean(x, na.rm=T))
    }
    if (grepl("SVM", model_name, ignore.case = TRUE)) {
      y.ens.avg.pool.SVM<-apply(matrix(y.ens.SVM.pool,length(outcome),iter),1,function(x) mean(x, na.rm=T))
    }
    if (grepl("DCT", model_name, ignore.case = TRUE)) {
      y.ens.avg.pool.DCT<-apply(matrix(y.ens.DCT.pool,length(outcome),iter),1,function(x) mean(x, na.rm=T))
    }
  }

  y.pred.avg.pool.all.W<-apply(matrix(y.pred.avg.pool.all,length(outcome),iter),1,function(x) mean(x, na.rm=T))

  y.pred.avg.pool1.all.S<-apply(matrix(y.pred.avg.pool1.all,length(outcome),iter),1,function(x) mean(x, na.rm=T))

  outcome.pool <- outcome
  model.names <- c()
  auc.tab <- c()
  ci.low.tab <- c()
  ci.up.tab <- c()
  sp.98.tab <- c()
  sp.95.tab <- c()
  for (model_name in models.ens) {
    if (grepl("XG", model_name, ignore.case = TRUE)) {
      roc.XG.ens <- survivalROC(Stime = survdays, status = outcome, marker = y.ens.avg.pool.XG, predict.time = time_point, method = "KM")
      auc.tab <- c(auc.tab, roc.XG.ens$AUC)
      #ci.low.tab <- c(ci.low.tab, roc.XG.ens$ci[1])
      #ci.up.tab <- c(ci.up.tab, roc.XG.ens$ci[3])
      sp.98.tab <- c(sp.98.tab, roc.XG.ens$TP[roc.XG.ens$FP <= 0.02][1])
      sp.95.tab <- c(sp.95.tab, roc.XG.ens$TP[roc.XG.ens$FP <= 0.05][1])
      model.names <- c(model.names, "XG.ens")
    }
    if (grepl("RF", model_name, ignore.case = TRUE)) {
      roc.RF.ens <- survivalROC(Stime = survdays, status = outcome, marker = y.ens.avg.pool.RF, predict.time = time_point, method = "KM")
      auc.tab <- c(auc.tab, roc.RF.ens$AUC)
      #ci.low.tab <- c(ci.low.tab, roc.RF.ens$ci[1])
      #ci.up.tab <- c(ci.up.tab, roc.RF.ens$ci[3])
      sp.98.tab <- c(sp.98.tab, roc.RF.ens$TP[roc.RF.ens$FP <= 0.02][1])
      sp.95.tab <- c(sp.95.tab, roc.RF.ens$TP[roc.RF.ens$FP <= 0.05][1])
      model.names <- c(model.names, "RF.ens")
    }
    if (grepl("EN", model_name, ignore.case = TRUE)) {
      roc.EN.ens <- survivalROC(Stime = survdays, status = outcome, marker = y.ens.avg.pool.EN, predict.time = time_point, method = "KM")
      auc.tab <- c(auc.tab, roc.EN.ens$AUC)
      #ci.low.tab <- c(ci.low.tab, roc.EN.ens$ci[1])
      #ci.up.tab <- c(ci.up.tab, roc.EN.ens$ci[3])
      sp.98.tab <- c(sp.98.tab, roc.EN.ens$TP[roc.EN.ens$FP <= 0.02][1])
      sp.95.tab <- c(sp.95.tab, roc.EN.ens$TP[roc.EN.ens$FP <= 0.05][1])
      model.names <- c(model.names, "EN.ens")
    }
    if (grepl("MB", model_name, ignore.case = TRUE)) {
      roc.MB.ens <- survivalROC(Stime = survdays, status = outcome, marker = y.ens.avg.pool.MB, predict.time = time_point, method = "KM")
      auc.tab <- c(auc.tab, roc.MB.ens$AUC)
      #ci.low.tab <- c(ci.low.tab, roc.NN.ens$ci[1])
      #ci.up.tab <- c(ci.up.tab, roc.NN.ens$ci[3])
      sp.98.tab <- c(sp.98.tab, roc.MB.ens$TP[roc.MB.ens$FP <= 0.02][1])
      sp.95.tab <- c(sp.95.tab, roc.MB.ens$TP[roc.MB.ens$FP <= 0.05][1])
      model.names <- c(model.names, "MB.ens")
    }
    if (grepl("ADB", model_name, ignore.case = TRUE)) {
      roc.ADB.ens <- survivalROC(Stime = survdays, status = outcome, marker = y.ens.avg.pool.ADB, predict.time = time_point, method = "KM")
      auc.tab <- c(auc.tab, roc.ADB.ens$AUC)
      sp.98.tab <- c(sp.98.tab, roc.ADB.ens$TP[roc.ADB.ens$FP <= 0.02][1])
      sp.95.tab <- c(sp.95.tab, roc.ADB.ens$TP[roc.ADB.ens$FP <= 0.05][1])
      model.names <- c(model.names, "ADB.ens")
    }
    if (grepl("GB", model_name, ignore.case = TRUE)) {
      roc.GB.ens <- survivalROC(Stime = survdays, status = outcome, marker = y.ens.avg.pool.GB, predict.time = time_point, method = "KM")
      auc.tab <- c(auc.tab, roc.GB.ens$AUC)
      #ci.low.tab <- c(ci.low.tab, roc.GB.ens$ci[1])
      #ci.up.tab <- c(ci.up.tab, roc.GB.ens$ci[3])
      sp.98.tab <- c(sp.98.tab, roc.GB.ens$TP[roc.GB.ens$FP <= 0.02][1])
      sp.95.tab <- c(sp.95.tab, roc.GB.ens$TP[roc.GB.ens$FP <= 0.05][1])
      model.names <- c(model.names, "GB.ens")
    }
    if (grepl("KNN", model_name, ignore.case = TRUE)) {
      roc.KNN.ens <- survivalROC(Stime = survdays, status = outcome, marker = y.ens.avg.pool.KNN, predict.time = time_point, method = "KM")
      auc.tab <- c(auc.tab, roc.KNN.ens$AUC)
      sp.98.tab <- c(sp.98.tab, roc.KNN.ens$TP[roc.KNN.ens$FP <= 0.02][1])
      sp.95.tab <- c(sp.95.tab, roc.KNN.ens$TP[roc.KNN.ens$FP <= 0.05][1])
      model.names <- c(model.names, "KNN.ens")
    }
    if (grepl("NN", model_name, ignore.case = TRUE)) {
      roc.NN.ens <- survivalROC(Stime = survdays, status = outcome, marker = y.ens.avg.pool.NN, predict.time = time_point, method = "KM")
      auc.tab <- c(auc.tab, roc.NN.ens$AUC)
      sp.98.tab <- c(sp.98.tab, roc.NN.ens$TP[roc.NN.ens$FP <= 0.02][1])
      sp.95.tab <- c(sp.95.tab, roc.NN.ens$TP[roc.NN.ens$FP <= 0.05][1])
      model.names <- c(model.names, "NN.ens")
    }
    if (grepl("SVM", model_name, ignore.case = TRUE)) {
      roc.SVM.ens <- survivalROC(Stime = survdays, status = outcome, marker = y.ens.avg.pool.SVM, predict.time = time_point, method = "KM")
      auc.tab <- c(auc.tab, roc.SVM.ens$AUC)
      #ci.low.tab <- c(ci.low.tab, roc.SVM.ens$ci[1])
      #ci.up.tab <- c(ci.up.tab, roc.SVM.ens$ci[3])
      sp.98.tab <- c(sp.98.tab, roc.SVM.ens$TP[roc.SVM.ens$FP <= 0.02][1])
      sp.95.tab <- c(sp.95.tab, roc.SVM.ens$TP[roc.SVM.ens$FP <= 0.05][1])
      model.names <- c(model.names, "SVM.ens")
    }
    if (grepl("DCT", model_name, ignore.case = TRUE)) {
      roc.DCT.ens <- survivalROC(Stime = survdays, status = outcome, marker = y.ens.avg.pool.DCT, predict.time = time_point, method = "KM")
      auc.tab <- c(auc.tab, roc.DCT.ens$AUC)
      #ci.low.tab <- c(ci.low.tab, roc.DCT.ens$ci[1])
      #ci.up.tab <- c(ci.up.tab, roc.DCT.ens$ci[3])
      sp.98.tab <- c(sp.98.tab, roc.DCT.ens$TP[roc.DCT.ens$FP <= 0.02][1])
      sp.95.tab <- c(sp.95.tab, roc.DCT.ens$TP[roc.DCT.ens$FP <= 0.05][1])
      model.names <- c(model.names, "DCT.ens")
    }
  }
  roc.all <- survivalROC(Stime = survdays, status = outcome, marker = y.pred.avg.pool.all.W, predict.time = time_point, method = "KM")
  auc.tab <- c(auc.tab, roc.all$AUC)

  #ci.low.tab <- c(ci.low.tab, roc.all$ci[1])
  #ci.up.tab <- c(ci.up.tab, roc.all$ci[3])
  sp.98.tab <- c(sp.98.tab, roc.all$TP[roc.all$FP <= 0.02][1])
  sp.95.tab <- c(sp.95.tab, roc.all$TP[roc.all$FP <= 0.05][1])
  model.names <- c(model.names, "Weighted")

  roc.all1 <- survivalROC(Stime = survdays, status = outcome, marker = y.pred.avg.pool1.all.S, predict.time = time_point, method = "KM")
  auc.tab <- c(auc.tab, roc.all1$AUC)
  #ci.low.tab <- c(ci.low.tab, roc.all1$ci[1])
  #ci.up.tab <- c(ci.up.tab, roc.all1$ci[3])
  sp.98.tab <- c(sp.98.tab, roc.all1$TP[roc.all1$FP <= 0.02][1])
  sp.95.tab <- c(sp.95.tab, roc.all1$TP[roc.all1$FP <= 0.05][1])
  model.names <- c(model.names, "Simple")

  get_marker <- function(model_name) {
    if (identical(model_name, "Weighted")) return(y.pred.avg.pool.all.W)
    if (identical(model_name, "Simple")) return(y.pred.avg.pool1.all.S)
    base_name <- sub("\\.ens$", "", model_name)
    get(paste0("y.ens.avg.pool.", base_name), inherits = TRUE)
  }
  cindex_stats <- t(vapply(
    model.names,
    function(model_name) survival_cindex(survdays, outcome, get_marker(model_name)),
    numeric(3)
  ))
  sp95_ci <- t(vapply(
    model.names,
    function(model_name) survival_sp_ci(survdays, outcome, get_marker(model_name), time_point),
    numeric(2)
  ))
  auc_ci <- t(vapply(
    model.names,
    function(model_name) survival_auc_ci(survdays, outcome, get_marker(model_name), time_point),
    numeric(2)
  ))

  df <- data.frame(
    model = model.names,
    auc = fix_auc(auc.tab),
    auc.ci.low = auc_ci[, "ci.low"],
    auc.ci.up = auc_ci[, "ci.up"],
    cindex = cindex_stats[, "cindex"],
    cindex.ci.low = cindex_stats[, "ci.low"],
    cindex.ci.up = cindex_stats[, "ci.up"],
    sp95 = sp.95.tab,
    sp95.ci.low = sp95_ci[, "ci.low"],
    sp95.ci.up = sp95_ci[, "ci.up"],
    sp98 = sp.98.tab
  )

  return(df)
}

#=======CATEGORICAL (Works) ======
get_multiclass_auc_ci <- function(y_true, y_pred_mat, nboot = 200, seed = 123) {
  set.seed(seed)
  aucs <- numeric(nboot)

  n <- length(y_true)
  num_classes <- ncol(y_pred_mat)

  for (i in 1:nboot) {
    # sample indices with replacement
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    y_true_boot <- y_true[idx]
    y_pred_boot <- y_pred_mat[idx, , drop = FALSE] #subset

    # compute multiclass AUC
    roc_boot <- pROC::multiclass.roc(y_true_boot, y_pred_boot)
    aucs[i] <- pROC::auc(roc_boot)
  }

  # 95% confidence interval
  ci <- quantile(aucs, probs = c(0.025, 0.975)) #bootstrap 95% CI

  return(list(
    auc_mean = mean(aucs),
    ci_lower = ci[1],
    ci_upper = ci[2],
    all_bootstrap_aucs = aucs
  ))
}
cat_train_performance <- function(outcome, train_data, features, models.indiv, iter, num.folds = NULL, ensemble = FALSE, models.ens = NULL, param.indiv = NULL, param.ens = NULL, iter.offset = 0) {

  num_classes <- length(unique(outcome))
  nfold <- max(2L, as.integer(num.folds %||% 3L))
  indiv.defaults <- list(
    XG = list(max_depth = 2, eta = 0.1, nrounds = 50, verbose = FALSE),
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
    XG = list(max_depth = 2, eta = 0.1, nrounds = 50),
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

  # Initialize pooling stuff

  if (ensemble == TRUE) {
    for (model_name in models.ens) {
      if (toupper(model_name) == "XG") {
        y.ens.XG.pool <- list()
      }
      if (toupper(model_name) == "RF") {
        y.ens.RF.pool <- list()
      }
      if (toupper(model_name) == "EN") {
        y.ens.EN.pool <- list()
      }
      if (toupper(model_name) == "ADB") {
        y.ens.ADB.pool <- list()
      }
      if (toupper(model_name) == "MB") {
        y.ens.MB.pool <- list()
      }
      if (toupper(model_name) == "GB") {
        y.ens.GB.pool <- list()
      }
      if (toupper(model_name) == "SVM") {
        y.ens.SVM.pool <- list()
      }
      if (toupper(model_name) == "KNN") {
        y.ens.KNN.pool <- list()
      }
      if (toupper(model_name) == "NN") {
        y.ens.NN.pool <- list()
      }
      if (toupper(model_name) == "DCT") {
        y.ens.DCT.pool <- list()
      }
      if (toupper(model_name) == "NB") {
        y.ens.NB.pool <- list()
      }
    }
  }
  y.pred.avg.pool1.all <- list()

  for(j in 1:iter){
    iter.id <- j + iter.offset
    message(iter.id)
    set.seed(iter.id)

    cv.fold <- balanced.folds.vec(outcome, nfold)

    if (ensemble == TRUE) {
      for (model_name in models.ens) {
        if (toupper(model_name) == "KNN") {
          y.ens.KNN <- matrix(NA, nrow = length(outcome), ncol = length(unique(outcome)))
        }
        if (toupper(model_name) == "XG") {
          y.ens.XG <- matrix(NA, nrow = length(outcome), ncol = length(unique(outcome)))
        }
        if (toupper(model_name) == "RF") {
          y.ens.RF <- matrix(NA, nrow = length(outcome), ncol = length(unique(outcome)))
        }
        if (toupper(model_name) == "EN") {
          y.ens.EN <- matrix(NA, nrow = length(outcome), ncol = length(unique(outcome)))
        }
        if (toupper(model_name) == "ADB") {
          y.ens.ADB <- matrix(NA, nrow = length(outcome), ncol = length(unique(outcome)))
        }
        if (toupper(model_name) == "MB") {
          y.ens.MB <- matrix(NA, nrow = length(outcome), ncol = length(unique(outcome)))
        }
        if (toupper(model_name) == "NN") {
          y.ens.NN <- matrix(NA, nrow = length(outcome), ncol = length(unique(outcome)))
        }
        if (toupper(model_name) == "GB") {
          y.ens.GB <- matrix(NA, nrow = length(outcome), ncol = length(unique(outcome)))
        }
        if (toupper(model_name) == "SVM") {
          y.ens.SVM <- matrix(NA, nrow = length(outcome), ncol = length(unique(outcome)))
        }
        if (toupper(model_name) == "DCT") {
          y.ens.DCT <- matrix(NA, nrow = length(outcome), ncol = length(unique(outcome)))
        }
        if (toupper(model_name) == "NB") {
          y.ens.NB <- matrix(NA, nrow = length(outcome), ncol = length(unique(outcome)))
        }
      }
    }

    #y.pred.avg.all <- matrix(NA, nrow = length(outcome), ncol = length(unique(outcome))) # No weighted avg

    mat.test.full<- NULL
    for(i in seq_len(nfold))
    {
      cur.train=which(cv.fold!=i)
      cur.test=which(cv.fold==i)

      mat.train <- matrix(nrow = length(cur.train), ncol = 0)
      mat.test  <- matrix(nrow = length(cur.test),  ncol = 0)
      y.pred.train.all<- matrix(nrow = length(cur.train), ncol = 0)
      y.pred.mat.all<- matrix(nrow = length(cur.test),  ncol = 0)

      for(f in seq_along(features))
      {
        message(paste("Pathway Number:", f))
        com<-unique(features[[f]])
        data_train <- train_data[, colnames(train_data) %in% com, drop = FALSE]
        colnames(data_train) <- make.names(colnames(data_train), unique = TRUE)

        xtrain<-(data_train)[cur.train,]
        xtest<-(data_train)[cur.test,]

        xtrain_cl <- as.data.frame(xtrain)
        xtrain_cl$cluster <- factor(paste0("C", outcome[cur.train]))
        ytrain <- outcome[cur.train]

        xtest_cl <- as.data.frame(xtest)
        xtest_cl$cluster <- factor(paste0("C", outcome[cur.test]))

        for (model_name in models.indiv) {
          if (toupper(model_name) == "XG") {
            y.pred.XG <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
            y.train_tr.XG <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
            y.train_ts.XG <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
          }
          if (toupper(model_name) == "EN") {
            y.pred.EN <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
            y.train_tr.EN <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
            y.train_ts.EN <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
          }
          if (toupper(model_name) == "RF") {
            y.pred.RF <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
            y.train_tr.RF <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
            y.train_ts.RF <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
          }
          if (toupper(model_name) == "ADB") {
            y.pred.ADB <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
            y.train_tr.ADB <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
            y.train_ts.ADB <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
          }
          if (toupper(model_name) == "MB") {
            y.pred.MB <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
            y.train_tr.MB <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
            y.train_ts.MB <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
          }
          if (toupper(model_name) == "GB") {
            y.pred.GB <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
            y.train_tr.GB <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
            y.train_ts.GB <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
          }
          if (toupper(model_name) == "SVM") {
            y.pred.SVM <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
            y.train_tr.SVM <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
            y.train_ts.SVM <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
          }
          if (toupper(model_name) == "NB") {
            y.pred.NB <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
            y.train_tr.NB <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
            y.train_ts.NB <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
          }
          if (toupper(model_name) == "KNN") {
            y.pred.KNN <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
            y.train_tr.KNN <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
            y.train_ts.KNN <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
          }
          if (toupper(model_name) == "DCT") {
            y.pred.DCT <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
            y.train_tr.DCT <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
            y.train_ts.DCT <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
          }
          if (toupper(model_name) == "NN") {
            y.pred.NN <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
            y.train_tr.NN <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
            y.train_ts.NN <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
          }
        }

        for (model_name in models.indiv) {
          num_classes <- length(unique(outcome))
          if (toupper(model_name) == "XG") {
            message("XgBoost running...")
            params_xg <- modifyList(indiv.defaults$XG,param.indiv[["XG"]] %||% list())
            num_classes <- length(levels(outcome[cur.train]))
            temp <- quiet(do.call(
              xgboost::xgboost,
              c(
                list(
                  x = as.matrix(xtrain),
                  y = as.factor(ytrain)
                ),
                params_xg
              )
            ))
            pred_probs <- quiet(predict(temp, newdata = as.matrix(xtest)))
            pred_matrix <- matrix(pred_probs, nrow = nrow(xtest), ncol = num_classes)
            y.pred.XG <- pred_matrix
            y.train_ts.XG <- pred_matrix

            pred_probs <- quiet(predict(temp, newdata = as.matrix(xtrain)))
            pred_matrix <- matrix(pred_probs, nrow = nrow(xtrain), ncol = num_classes)
            y.train_tr.XG <- pred_matrix
          }
          if (toupper(model_name) == "RF") {
            message("RandomForest running...")
            params_rf <- modifyList(indiv.defaults$RF, param.indiv[["RF"]] %||% list())
            temp <- quiet(do.call(randomForest, c(list(x = xtrain, y = as.factor(outcome[cur.train])), params_rf)))
            y.pred.RF <- quiet(predict(temp, newdata = xtest, type = "prob"))
            #y.train_ts.RF[cur.test,1:num_classes] <- predict(temp, newdata = xtest, type = "prob")
            y.train_ts.RF <- quiet(predict(temp, newdata = xtest, type = "prob"))
            y.train_tr.RF <- quiet(predict(temp, newdata = xtrain, type = "prob"))
          }
          if (toupper(model_name) == "EN") {
            message("ElasticNet running...")
            params_en <- modifyList(indiv.defaults$EN, param.indiv[["EN"]] %||% list())
            train.params <- params_en[setdiff(names(params_en), "s")]
            temp <- quiet(do.call(cv.glmnet, c(list(xtrain, as.factor(outcome[cur.train]), family = "multinomial"), train.params)))
            s.val <- param.indiv[["EN"]]$s %||% "lambda.1se"
            y.pred.EN <- quiet(predict(temp, newx = xtest, s = s.val, type = "response"))
            y.train_tr.EN <- quiet(predict(temp, newx = xtrain, s = s.val, type = "response"))
            y.train_ts.EN <- quiet(predict(temp, newx = xtest, s = s.val, type = "response"))
          }
          if (toupper(model_name) == "ADB") {
            message("AdaBoost one-vs-rest running...")
            params_adb <- modifyList(indiv.defaults$ADB, param.indiv[["ADB"]] %||% list())
            y.pred.ADB <- quiet(do.call(
              multiclass_adb_predict,
              c(list(xtrain = xtrain, outcome = outcome[cur.train], xtest = xtest,
                     classes = sort(unique(outcome))), params_adb)
            ))
            y.train_tr.ADB <- quiet(do.call(
              multiclass_adb_predict,
              c(list(xtrain = xtrain, outcome = outcome[cur.train],
                     classes = sort(unique(outcome))), params_adb)
            ))
            y.train_ts.ADB <- y.pred.ADB
          }
          if (toupper(model_name) == "MB") {
            message("MBoost one-vs-rest running...")
            params_mb <- modifyList(indiv.defaults$MB, param.indiv[["MB"]] %||% list())
            y.pred.MB <- quiet(do.call(
              multiclass_mb_predict,
              c(list(xtrain = xtrain, outcome = outcome[cur.train], xtest = xtest,
                     classes = sort(unique(outcome))), params_mb)
            ))
            y.train_tr.MB <- quiet(do.call(
              multiclass_mb_predict,
              c(list(xtrain = xtrain, outcome = outcome[cur.train],
                     classes = sort(unique(outcome))), params_mb)
            ))
            y.train_ts.MB <- y.pred.MB
          }
          if (toupper(model_name) == "GB") {
            message("GradientBoost running...")
            params_gb <- modifyList(indiv.defaults$GB, param.indiv[["GB"]] %||% list())
            temp <- quiet(do.call(gbm.fit, c(list(xtrain, as.factor(outcome[cur.train]), distribution = "multinomial"), params_gb)))
            y.pred.GB <- quiet(predict.gbm(temp, data.frame(xtest), type = "response"))
            y.train_tr.GB <- quiet(predict.gbm(temp, data.frame(xtrain), type = "response"))
            y.train_ts.GB <- quiet(predict.gbm(temp, data.frame(xtest), type = "response"))
          }
          if (toupper(model_name) == "KNN") {
            message("KNN running...")
            params_knn <- modifyList(indiv.defaults$KNN, param.indiv[["KNN"]] %||% list())
            params_knn$trControl <- NULL  # remove so not passed twice

            temp <- quiet(do.call(caret::train, c(
              list(
                cluster ~ .,
                data = xtrain_cl,
                method = "knn",
                trControl = trc,
                preProcess = c("center", "scale")
              ),
              params_knn
            )))
            y.pred.KNN <- quiet(as.matrix(predict(temp, newdata = xtest_cl, type = "prob")))
            y.train_tr.KNN <- quiet(as.matrix(predict(temp, newdata = xtrain_cl, type = "prob")))
            y.train_ts.KNN <- quiet(as.matrix(predict(temp, newdata = xtest_cl, type = "prob")))
          }
          if (toupper(model_name) == "SVM") {
            message("SVM running...")
            params_svm <- modifyList(indiv.defaults$SVM, param.indiv[["SVM"]] %||% list())
            temp <- quiet(do.call(svm, c(list(cluster ~ ., data = xtrain_cl, probability = TRUE), params_svm)))

            pred <- quiet(predict(temp, newdata= xtest_cl, probability = TRUE))
            probs <- attr(pred, "probabilities")
            class_order <- paste0("C", seq_len(num_classes))
            probs_reordered <- probs[, class_order, drop = FALSE]
            y.pred.SVM <- probs_reordered
            y.train_ts.SVM <- probs_reordered

            pred <- quiet(predict(temp, newdata= xtrain_cl, probability = TRUE))
            probs <- attr(pred, "probabilities")
            probs_reordered <- probs[, class_order, drop = FALSE]
            y.train_tr.SVM <- probs_reordered
          }
          if (toupper(model_name) == "NB") {
            message("NB running...")
            params_nb <- modifyList(indiv.defaults$NB, param.indiv[["NB"]] %||% list())
            temp <- quiet(do.call(naiveBayes, c(list(cluster ~., data = xtrain_cl), params_nb)))
            y.pred.NB <- quiet(predict(temp, xtest_cl, type = "raw"))
            y.train_tr.NB <- quiet(predict(temp, xtrain_cl, type = "raw"))
            y.train_ts.NB <- quiet(predict(temp, xtest_cl, type = "raw"))
          }
          if (toupper(model_name) == "DCT") {
            message("DecisionTree running...")
            params_dct <- modifyList(indiv.defaults$DCT, param.indiv[["DCT"]] %||% list())
            temp <- quiet(do.call(rpart, c(list(cluster ~., data = xtrain_cl, method = "class"), params_dct)))
            y.pred.DCT <- quiet(as.matrix(predict(temp, newdata= as.data.frame(xtest_cl), type="prob")))
            y.train_tr.DCT <- quiet(as.matrix(predict(temp, newdata= as.data.frame(xtrain_cl), type="prob")))
            y.train_ts.DCT <- quiet(as.matrix(predict(temp, newdata= as.data.frame(xtest_cl), type="prob")))
          }
          if (toupper(model_name) == "NN") {
            message("NeuralNet running...")
            m <- colMeans(as.data.frame(xtrain))
            s <- apply(as.data.frame(xtrain), 2, sd)
            xtrain.sd <- scale(xtrain, center = m, scale = s)
            xtest.sd  <- scale(xtest,  center = m, scale = s)

            y.mat  <- model.matrix(~ outcome[cur.train] - 1)  # one-hot encoding
            colnames(y.mat) <- paste0("class", 1:num_classes)
            train.df <- data.frame(xtrain.sd, y.mat)

            f <- as.formula(paste(paste(colnames(y.mat), collapse = " + "), "~ ."))
            params <- modifyList(indiv.defaults$NN, param.indiv[["NN"]] %||% list())
            temp.nn <- quiet(do.call(neuralnet, c(list(f, data = train.df, linear.output = FALSE), params)))

            y.pred.NN <- quiet(predict(temp.nn, newdata = xtest.sd))
            y.train_tr.NN <- quiet(predict(temp.nn, newdata = xtrain.sd))
            y.train_ts.NN <- quiet(predict(temp.nn, newdata = xtest.sd))
          }
        }

        for (model_name in models.ens) {
          if (toupper(model_name) == "XG") {
            mat.train <- cbind(mat.train, y.train_tr.XG)
            mat.test <- cbind(mat.test, y.train_ts.XG)
            y.pred.train.all <- cbind(y.pred.train.all, y.train_tr.XG)
            y.pred.mat.all <- cbind(y.pred.mat.all, y.train_ts.XG)
          }
          if (toupper(model_name) == "RF") {
            mat.train <- cbind(mat.train, y.train_tr.RF)
            mat.test <- cbind(mat.test, y.train_ts.RF)
            y.pred.train.all <- cbind(y.pred.train.all, y.train_tr.RF)
            y.pred.mat.all <- cbind(y.pred.mat.all, y.train_ts.RF)
          }
          if (toupper(model_name) == "EN") {
            mat.train <- cbind(mat.train, y.train_tr.EN)
            mat.test <- cbind(mat.test, y.train_ts.EN)
            y.pred.train.all <- cbind(y.pred.train.all, y.train_tr.EN)
            y.pred.mat.all <- cbind(y.pred.mat.all, y.train_ts.EN)
          }
          if (toupper(model_name) == "ADB") {
            mat.train <- cbind(mat.train, y.train_tr.ADB)
            mat.test <- cbind(mat.test, y.train_ts.ADB)
            y.pred.train.all <- cbind(y.pred.train.all, y.train_tr.ADB)
            y.pred.mat.all <- cbind(y.pred.mat.all, y.train_ts.ADB)
          }
          if (toupper(model_name) == "MB") {
            mat.train <- cbind(mat.train, y.train_tr.MB)
            mat.test <- cbind(mat.test, y.train_ts.MB)
            y.pred.train.all <- cbind(y.pred.train.all, y.train_tr.MB)
            y.pred.mat.all <- cbind(y.pred.mat.all, y.train_ts.MB)
          }
          if (toupper(model_name) == "GB") {
            mat.train <- cbind(mat.train, y.train_tr.GB)
            mat.test <- cbind(mat.test, y.train_ts.GB)
            y.pred.train.all <- cbind(y.pred.train.all, y.train_tr.GB)
            y.pred.mat.all <- cbind(y.pred.mat.all, y.train_ts.GB)
          }
          if (toupper(model_name) == "KNN") {
            mat.train <- cbind(mat.train, y.train_tr.KNN)
            mat.test <- cbind(mat.test, y.train_ts.KNN)
            y.pred.train.all <- cbind(y.pred.train.all, y.train_tr.KNN)
            y.pred.mat.all <- cbind(y.pred.mat.all, y.train_ts.KNN)
          }
          if (toupper(model_name) == "SVM") {
            mat.train <- cbind(mat.train, y.train_tr.SVM)
            mat.test <- cbind(mat.test, y.train_ts.SVM)
            y.pred.train.all <- cbind(y.pred.train.all, y.train_tr.SVM)
            y.pred.mat.all <- cbind(y.pred.mat.all, y.train_ts.SVM)
          }
          if (toupper(model_name) == "NB") {
            mat.train <- cbind(mat.train, y.train_tr.NB)
            mat.test <- cbind(mat.test, y.train_ts.NB)
            y.pred.train.all <- cbind(y.pred.train.all, y.train_tr.NB)
            y.pred.mat.all <- cbind(y.pred.mat.all, y.train_ts.NB)
          }
          if (toupper(model_name) == "DCT") {
            mat.train <- cbind(mat.train, y.train_tr.DCT)
            mat.test <- cbind(mat.test, y.train_ts.DCT)
            y.pred.train.all <- cbind(y.pred.train.all, y.train_tr.DCT)
            y.pred.mat.all <- cbind(y.pred.mat.all, y.train_ts.DCT)
          }
          if (toupper(model_name) == "NN") {
            mat.train <- cbind(mat.train, y.train_tr.NN)
            mat.test <- cbind(mat.test, y.train_ts.NN)
            y.pred.train.all <- cbind(y.pred.train.all, y.train_tr.NN)
            y.pred.mat.all <- cbind(y.pred.mat.all, y.train_ts.NN)
          }
        }
      }

      colnames(mat.train) <- NULL
      colnames(mat.test)  <- NULL
      colnames(y.pred.train.all) <- NULL
      colnames(y.pred.mat.all) <- NULL

      for (model_name in models.ens) {
        mat.train.ens_cl <- as.data.frame(mat.train)
        mat.train.ens_cl$cluster <- factor(paste0("C", outcome[cur.train]))

        mat.test.ens_cl <- as.data.frame(mat.test)
        mat.test.ens_cl$cluster <- factor(paste0("C", outcome[cur.test]))

        if (grepl("XG", model_name, ignore.case = TRUE)) {
          params_xg <- modifyList(ens.defaults$XG, param.ens[["XG"]] %||% list())
          temp <- quiet(do.call(
            xgboost::xgboost,
            c(
              list(
                x = as.matrix(mat.train),
                y = as.factor(ytrain)
              ),
              params_xg
            )
          ))
          pred_probs <- quiet(predict(temp, newdata = as.matrix(mat.test)))
          pred_matrix <- matrix(pred_probs, nrow = nrow(mat.test), ncol = num_classes)
          y.ens.XG[cur.test, ] <- pred_matrix
        }
        if (grepl("RF", model_name, ignore.case = TRUE)) {
          params_rf <- modifyList(ens.defaults$RF, param.ens[["RF"]] %||% list())
          temp <- quiet(do.call(randomForest, c(list(x = mat.train, y = as.factor(outcome[cur.train])), params_rf)))
          y.ens.RF[cur.test,] <- quiet(predict(temp, newdata = mat.test, type = "prob"))
        }
        if (grepl("EN", model_name, ignore.case = TRUE)) {
          params_en <- modifyList(ens.defaults$EN, param.ens[["EN"]] %||% list())
          train.params <- params_en[setdiff(names(params_en), "s")]
          temp <- quiet(do.call(cv.glmnet, c(list(mat.train, as.factor(outcome[cur.train]), family = "multinomial"), train.params)))
          s.val <- param.ens[["EN"]]$s %||% "lambda.1se"
          y.ens.EN[cur.test,] <- quiet(predict(temp, newx = mat.test, s = s.val, type = "response"))
        }
        if (grepl("ADB", model_name, ignore.case = TRUE)) {
          params_adb <- modifyList(ens.defaults$ADB, param.ens[["ADB"]] %||% list())
          y.ens.ADB[cur.test, ] <- quiet(do.call(
            multiclass_adb_predict,
            c(list(xtrain = mat.train, outcome = outcome[cur.train],
                   xtest = mat.test, classes = sort(unique(outcome))), params_adb)
          ))
        }
        if (grepl("MB", model_name, ignore.case = TRUE)) {
          params_mb <- modifyList(ens.defaults$MB, param.ens[["MB"]] %||% list())
          y.ens.MB[cur.test, ] <- quiet(do.call(
            multiclass_mb_predict,
            c(list(xtrain = mat.train, outcome = outcome[cur.train],
                   xtest = mat.test, classes = sort(unique(outcome))), params_mb)
          ))
        }
        if (grepl("GB", model_name, ignore.case = TRUE)) {
          xtrain.mat <- as.matrix(mat.train)
          xtest.mat <- as.matrix(mat.test)

          num_features_ens <- ncol(xtrain.mat)

          #    This prevents gbm from breaking if names are NULL or complex.
          feature_names <- paste0("EnsVar", 1:num_features_ens)

          colnames(xtrain.mat) <- feature_names
          colnames(xtest.mat) <- feature_names

          # Safety check (optional, but good practice):
          if (ncol(xtest.mat) != num_features_ens) {
            stop("FATAL ERROR: Test data columns do not match expected number in GB ensemble.")
          }

          params_gb <- modifyList(ens.defaults$GB, param.ens[["GB"]] %||% list())

          args_list <- c(
            list(
              x = xtrain.mat,
              y = outcome[cur.train],
              distribution = "multinomial"
            ),
            params_gb
          )
          temp <- quiet(do.call(gbm.fit, args_list))

          n_trees_final <- temp$n.trees

          # 2. Predict, explicitly passing the final tree count.
          y.pred.array <- quiet(predict.gbm(
            temp,
            newdata = as.data.frame(xtest.mat),
            n.trees = n_trees_final, # Pass total trees built
            type = "response"
          ))

          last_valid_index <- dim(y.pred.array)[3]

          # Slice the 3D array to get the 2D matrix of final probabilities
          final_probs_matrix <- y.pred.array[,,last_valid_index]
          y.ens.GB[cur.test,] <- final_probs_matrix
        }
        if (grepl("KNN", model_name, ignore.case = TRUE)) {
          trc <- param.ens[["KNN"]]$trControl %||% trainControl(method = "cv", number = 3, classProbs = TRUE)
          params_knn <- modifyList(ens.defaults$KNN, param.ens[["KNN"]] %||% list())
          params_knn$trControl <- NULL  # remove so not passed twice

          temp <- quiet(do.call(caret::train, c(
            list(
              cluster ~ .,
              data = mat.train.ens_cl,
              method = "knn",
              trControl = trc,
              preProcess = c("center", "scale")
            ),
            params_knn
          )))
          y.ens.KNN[cur.test, 1:length(unique(outcome))] <- quiet(as.matrix(predict(temp, newdata = mat.test.ens_cl, type = "prob")))
        }
        if (grepl("SVM", model_name, ignore.case = TRUE)) {
          params_svm <- modifyList(ens.defaults$SVM, param.ens[["SVM"]] %||% list())
          temp <- quiet(do.call(svm, c(list(cluster ~ ., data = data.frame(mat.train.ens_cl), probability = TRUE), params_svm)))

          pred <- quiet(predict(temp, newdata= mat.test.ens_cl, probability = TRUE))
          probs <- attr(pred, "probabilities")
          class_order <- levels(mat.train.ens_cl$cluster)
          probs_reordered <- probs[, class_order, drop = FALSE]
          y.ens.SVM[cur.test, ] <- probs_reordered
        }
        if (grepl("NB", model_name, ignore.case = TRUE)) {
          params_nb <- modifyList(ens.defaults$NB, param.ens[["NB"]] %||% list())
          temp <- do.call(naiveBayes, c(list(cluster ~., data = mat.train.ens_cl), params_nb))
          y.ens.NB[cur.test, 1:length(unique(outcome))] <- predict(temp, mat.test.ens_cl, type = "raw")
        }
        if (grepl("DCT", model_name, ignore.case = TRUE)) {
          params_dct <- modifyList(ens.defaults$DCT, param.ens[["DCT"]] %||% list())
          temp <- do.call(rpart, c(list(cluster ~., data = mat.train.ens_cl, method = "class"), params_dct))
          y.ens.DCT[cur.test, 1:length(unique(outcome))] <- as.matrix(predict(temp, newdata= as.data.frame(mat.test.ens_cl), type="prob"))
        }
        if (grepl("NN", model_name, ignore.case = TRUE)) {
          m <- colMeans(as.data.frame(mat.train))
          s <- apply(as.data.frame(mat.train), 2, sd)
          xtrain.ens.sd <- scale(mat.train, center = m, scale = s)
          xtest.ens.sd  <- scale(mat.test,  center = m, scale = s)

          y.mat  <- model.matrix(~ outcome[cur.train] - 1)  # one-hot encoding
          colnames(y.mat) <- paste0("class", 1:length(unique(outcome)))
          train.df <- data.frame(xtrain.ens.sd, y.mat)

          f <- as.formula(paste(paste(colnames(y.mat), collapse = " + "), "~ ."))
          params_nn <- modifyList(ens.defaults$NN, param.ens[["NN"]] %||% list())
          temp.nn <- do.call(neuralnet, c(list(f, data = train.df, linear.output = FALSE), params_nn))
          y.ens.NN[cur.test, ] <- predict(temp.nn, newdata = xtest.ens.sd)
        }
      }

      ########### weighted avg ###############
      #y.pred.train.all <- cbind.data.frame(outcome = outcome[cur.train], y.pred.train.all)

      #temp.avg.all <- multinom(outcome[cur.train] ~ ., data = y.pred.train.all, trace = FALSE)
      #y.pred.avg.all[cur.test, ] <- predict(temp.avg.all, newdata = y.pred.mat.all, type = "probs")

      mat.test.full <- rbind(mat.test.full, mat.test)

    }
    ########### simple avg ###############

    y.pred.avg1.all <- multiclass_simple_average(mat.test.full, num_classes)

    ############# pooling ensemble of methods ####################

    for (model_name in models.ens) {
      if (grepl("XG", model_name, ignore.case = TRUE)) {
        y.ens.XG.pool[[j]] <- y.ens.XG
      }
      if (grepl("RF", model_name, ignore.case = TRUE)) {
        y.ens.RF.pool[[j]] <- y.ens.RF
      }
      if (grepl("EN", model_name, ignore.case = TRUE)) {
        y.ens.EN.pool[[j]] <- y.ens.EN
      }
      if (grepl("ADB", model_name, ignore.case = TRUE)) {
        y.ens.ADB.pool[[j]] <- y.ens.ADB
      }
      if (grepl("MB", model_name, ignore.case = TRUE)) {
        y.ens.MB.pool[[j]] <- y.ens.MB
      }
      if (grepl("NN", model_name, ignore.case = TRUE)) {
        y.ens.NN.pool[[j]] <- y.ens.NN
      }
      if (grepl("GB", model_name, ignore.case = TRUE)) {
        y.ens.GB.pool[[j]] <- y.ens.GB
      }
      if (grepl("KNN", model_name, ignore.case = TRUE)) {
        y.ens.KNN.pool[[j]] <- y.ens.KNN
      }
      if (grepl("SVM", model_name, ignore.case = TRUE)) {
        y.ens.SVM.pool[[j]] <- y.ens.SVM
      }
      if (grepl("NB", model_name, ignore.case = TRUE)) {
        y.ens.NB.pool[[j]] <- y.ens.NB
      }
      if (grepl("DCT", model_name, ignore.case = TRUE)) {
        y.ens.DCT.pool[[j]] <- y.ens.DCT
      }
    }

    ############# pooling simple (no weighted) avg of methods ####################
    y.pred.avg.pool1.all[[j]] <- y.pred.avg1.all
  }

  ############# avg of iterations ####################

  for (model_name in models.ens) {
    if (grepl("XG", model_name, ignore.case = TRUE)) {
      y.ens.avg.pool.XG <- apply(simplify2array(y.ens.XG.pool), c(1,2), mean)
    }
    if (grepl("RF", model_name, ignore.case = TRUE)) {
      y.ens.avg.pool.RF <- apply(simplify2array(y.ens.RF.pool), c(1,2), mean)
    }
    if (grepl("EN", model_name, ignore.case = TRUE)) {
      y.ens.avg.pool.EN <- apply(simplify2array(y.ens.EN.pool), c(1,2), mean)
    }
    if (grepl("ADB", model_name, ignore.case = TRUE)) {
      y.ens.avg.pool.ADB <- apply(simplify2array(y.ens.ADB.pool), c(1,2), mean)
    }
    if (grepl("MB", model_name, ignore.case = TRUE)) {
      y.ens.avg.pool.MB <- apply(simplify2array(y.ens.MB.pool), c(1,2), mean)
    }
    if (grepl("NN", model_name, ignore.case = TRUE)) {
      y.ens.avg.pool.NN <- apply(simplify2array(y.ens.NN.pool), c(1,2), mean)
    }
    if (grepl("GB", model_name, ignore.case = TRUE)) {
      y.ens.avg.pool.GB <- apply(simplify2array(y.ens.GB.pool), c(1,2), mean)
    }
    if (grepl("KNN", model_name, ignore.case = TRUE)) {
      y.ens.avg.pool.KNN <- apply(simplify2array(y.ens.KNN.pool), c(1,2), mean)
    }
    if (grepl("SVM", model_name, ignore.case = TRUE)) {
      y.ens.avg.pool.SVM <- apply(simplify2array(y.ens.SVM.pool), c(1,2), mean)
    }
    if (grepl("NB", model_name, ignore.case = TRUE)) {
      y.ens.avg.pool.NB <- apply(simplify2array(y.ens.NB.pool), c(1,2), mean)
    }
    if (grepl("DCT", model_name, ignore.case = TRUE)) {
      y.ens.avg.pool.DCT <- apply(simplify2array(y.ens.DCT.pool), c(1,2), mean)
    }
  }

  #y.pred.avg.pool.all.W<-apply(matrix(y.pred.avg.pool.all,length(outcome),iter),1,function(x) mean(x, na.rm=T))

  y.pred.avg.pool1.all.S<-apply(simplify2array(y.pred.avg.pool1.all), c(1,2), mean)

  outcome.pool <- outcome
  model.names <- c()
  auc.tab <- c()
  ci.low.tab <- c()
  ci.up.tab <- c()
  #sp.98.tab <- c()
  #sp.95.tab <- c()
  for (model_name in models.ens) {
    num_classes <- length(unique(outcome))
    if (grepl("XG", model_name, ignore.case = TRUE)) {
      XG_mat <- matrix(y.ens.avg.pool.XG, ncol = num_classes, byrow = TRUE)
      colnames(XG_mat) <- levels(outcome)
      roc.XG <- multiclass.roc(outcome.pool, XG_mat)
      roc.XG$ci <- get_multiclass_auc_ci(outcome.pool, XG_mat, nboot = 200)
      model.names <- c(model.names, "XG.ens")
      ci.low.tab <- c(ci.low.tab, roc.XG$ci$ci_lower)
      ci.up.tab <- c(ci.up.tab, roc.XG$ci$ci_upper)
      auc.tab <- c(auc.tab, roc.XG$auc)
    }
    if (grepl("RF", model_name, ignore.case = TRUE)) {
      RF_mat <- matrix(y.ens.avg.pool.RF, ncol = num_classes, byrow = TRUE)
      colnames(RF_mat) <- levels(outcome)
      roc.RF <- multiclass.roc(outcome.pool, RF_mat)
      roc.RF$ci <- get_multiclass_auc_ci(outcome.pool, RF_mat, nboot = 200)
      model.names <- c(model.names, "RF.ens")
      ci.low.tab <- c(ci.low.tab, roc.RF$ci$ci_lower)
      ci.up.tab <- c(ci.up.tab, roc.RF$ci$ci_upper)
      auc.tab <- c(auc.tab, roc.RF$auc)
    }
    if (grepl("EN", model_name, ignore.case = TRUE)) {
      EN_mat <- matrix(y.ens.avg.pool.EN, ncol = num_classes, byrow = TRUE)
      colnames(EN_mat) <- levels(outcome)
      roc.EN <- multiclass.roc(outcome.pool, EN_mat)
      roc.EN$ci <- get_multiclass_auc_ci(outcome.pool, EN_mat, nboot = 200)
      model.names <- c(model.names, "EN.ens")
      ci.low.tab <- c(ci.low.tab, roc.EN$ci$ci_lower)
      ci.up.tab <- c(ci.up.tab, roc.EN$ci$ci_upper)
      auc.tab <- c(auc.tab, roc.EN$auc)
    }
    if (grepl("ADB", model_name, ignore.case = TRUE)) {
      ADB_mat <- matrix(y.ens.avg.pool.ADB, ncol = num_classes, byrow = TRUE)
      colnames(ADB_mat) <- levels(outcome)
      roc.ADB <- multiclass.roc(outcome.pool, ADB_mat)
      roc.ADB$ci <- get_multiclass_auc_ci(outcome.pool, ADB_mat, nboot = 200)
      model.names <- c(model.names, "ADB.ens")
      ci.low.tab <- c(ci.low.tab, roc.ADB$ci$ci_lower)
      ci.up.tab <- c(ci.up.tab, roc.ADB$ci$ci_upper)
      auc.tab <- c(auc.tab, roc.ADB$auc)
    }
    if (grepl("MB", model_name, ignore.case = TRUE)) {
      MB_mat <- matrix(y.ens.avg.pool.MB, ncol = num_classes, byrow = TRUE)
      colnames(MB_mat) <- levels(outcome)
      roc.MB <- multiclass.roc(outcome.pool, MB_mat)
      roc.MB$ci <- get_multiclass_auc_ci(outcome.pool, MB_mat, nboot = 200)
      model.names <- c(model.names, "MB.ens")
      ci.low.tab <- c(ci.low.tab, roc.MB$ci$ci_lower)
      ci.up.tab <- c(ci.up.tab, roc.MB$ci$ci_upper)
      auc.tab <- c(auc.tab, roc.MB$auc)
    }
    if (grepl("NN", model_name, ignore.case = TRUE)) {
      NN_mat <- matrix(y.ens.avg.pool.NN, ncol = num_classes, byrow = TRUE)
      colnames(NN_mat) <- levels(outcome)
      roc.NN <- multiclass.roc(outcome.pool, NN_mat)
      roc.NN$ci <- get_multiclass_auc_ci(outcome.pool, NN_mat, nboot = 200)
      model.names <- c(model.names, "NN.ens")
      ci.low.tab <- c(ci.low.tab, roc.NN$ci$ci_lower)
      ci.up.tab <- c(ci.up.tab, roc.NN$ci$ci_upper)
      auc.tab <- c(auc.tab, roc.NN$auc)
    }
    if (grepl("GB", model_name, ignore.case = TRUE)) {
      GB_mat <- matrix(y.ens.avg.pool.GB, ncol = num_classes, byrow = TRUE)
      colnames(GB_mat) <- levels(outcome)
      roc.GB <- multiclass.roc(outcome.pool, GB_mat)
      roc.GB$ci <- get_multiclass_auc_ci(outcome.pool, GB_mat, nboot = 200)
      model.names <- c(model.names, "GB.ens")
      ci.low.tab <- c(ci.low.tab, roc.GB$ci$ci_lower)
      ci.up.tab <- c(ci.up.tab, roc.GB$ci$ci_upper)
      auc.tab <- c(auc.tab, roc.GB$auc)
    }
    if (grepl("KNN", model_name, ignore.case = TRUE)) {
      KNN_mat <- matrix(y.ens.avg.pool.KNN, ncol = num_classes, byrow = TRUE)
      colnames(KNN_mat) <- levels(outcome)
      roc.KNN <- multiclass.roc(outcome.pool, KNN_mat)
      roc.KNN$ci <- get_multiclass_auc_ci(outcome.pool, KNN_mat, nboot = 200)
      model.names <- c(model.names, "KNN.ens")
      ci.low.tab <- c(ci.low.tab, roc.KNN$ci$ci_lower)
      ci.up.tab <- c(ci.up.tab, roc.KNN$ci$ci_upper)
      auc.tab <- c(auc.tab, roc.KNN$auc)
    }
    if (grepl("SVM", model_name, ignore.case = TRUE)) {
      SVM_mat <- matrix(y.ens.avg.pool.SVM, ncol = num_classes, byrow = TRUE)
      colnames(SVM_mat) <- levels(outcome)
      roc.SVM <- multiclass.roc(outcome.pool, SVM_mat)
      roc.SVM$ci <- get_multiclass_auc_ci(outcome.pool, SVM_mat, nboot = 200)
      model.names <- c(model.names, "SVM.ens")
      ci.low.tab <- c(ci.low.tab, roc.SVM$ci$ci_lower)
      ci.up.tab <- c(ci.up.tab, roc.SVM$ci$ci_upper)
      auc.tab <- c(auc.tab, roc.SVM$auc)
    }
    if (grepl("NB", model_name, ignore.case = TRUE)) {
      NB_mat <- matrix(y.ens.avg.pool.NB, ncol = num_classes, byrow = TRUE)
      colnames(NB_mat) <- levels(outcome)
      roc.NB <- multiclass.roc(outcome.pool, NB_mat)
      roc.NB$ci <- get_multiclass_auc_ci(outcome.pool, NB_mat, nboot = 200)
      model.names <- c(model.names, "NB.ens")
      ci.low.tab <- c(ci.low.tab, roc.NB$ci$ci_lower)
      ci.up.tab <- c(ci.up.tab, roc.NB$ci$ci_upper)
      auc.tab <- c(auc.tab, roc.NB$auc)
    }
    if (grepl("DCT", model_name, ignore.case = TRUE)) {
      DCT_mat <- matrix(y.ens.avg.pool.DCT, ncol = num_classes, byrow = TRUE)
      colnames(DCT_mat) <- levels(outcome)
      roc.DCT <- multiclass.roc(outcome.pool, DCT_mat)
      roc.DCT$ci <- get_multiclass_auc_ci(outcome.pool, DCT_mat, nboot = 200)
      model.names <- c(model.names, "DCT.ens")
      ci.low.tab <- c(ci.low.tab, roc.DCT$ci$ci_lower)
      ci.up.tab <- c(ci.up.tab, roc.DCT$ci$ci_upper)
      auc.tab <- c(auc.tab, roc.DCT$auc)
    }
  }

  num_samples <- nrow(y.pred.avg.pool1.all.S)
  num_models  <- length(models.ens)
  num_features <- length(features)

  # reshape to 3D: samples x classes x (features*models)
  arr <- array(y.pred.avg.pool1.all.S, dim = c(num_samples, num_classes, num_features*num_models))

  # average per class across features and models
  s_mat <- apply(arr, c(1,2), mean)

  # now s_mat is n_samples x n_classes
  colnames(s_mat) <- levels(outcome)
  roc.all1 <- multiclass.roc(outcome.pool, s_mat)
  roc.all1$ci <- get_multiclass_auc_ci(outcome.pool, s_mat, nboot = 200)
  model.names <- c(model.names, "simple")
  ci.low.tab <- c(ci.low.tab, roc.all1$ci$ci_lower)
  ci.up.tab <- c(ci.up.tab, roc.all1$ci$ci_upper)
  auc.tab <- c(auc.tab, roc.all1$auc)

  get_marker <- function(model_name) {
    if (identical(model_name, "simple")) return(s_mat)
    base_name <- sub("\\.ens$", "", model_name)
    get(paste0(base_name, "_mat"), inherits = TRUE)
  }
  sp95_tab <- vapply(
    model.names,
    function(model_name) multiclass_sp95(outcome.pool, get_marker(model_name)),
    numeric(1)
  )
  sp95_ci <- t(vapply(
    model.names,
    function(model_name) categorical_sp_ci(outcome.pool, get_marker(model_name)),
    numeric(2)
  ))

  df <- data.frame(
    model = model.names,
    auc = auc.tab,
    ci.low = ci.low.tab,
    ci.up = ci.up.tab,
    sp95 = sp95_tab,
    sp95.ci.low = sp95_ci[, "ci.low"],
    sp95.ci.up = sp95_ci[, "ci.up"]
  )

  return(df)
}
