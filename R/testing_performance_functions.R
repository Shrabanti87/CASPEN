# Internal Function - Not exported
`%||%` <- function(a, b) if (!is.null(a)) a else b
######### BINARY ###########
binary_test_performance <- function(outcome, train_data, test_data, test_outcome, features, models.indiv, ensemble = FALSE, models.ens = NULL, param.indiv = NULL, param.ens = NULL) {

  ##### Parameters ####
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
    NN = list(hidden = c(8,4), rep = 1),
    NB = list(usekernel = TRUE),
    DCT = list(method = "class")
  )

  mat.train<- NULL; mat.test<- NULL
  y.pred.train.all<- NULL; y.pred.mat.all<- NULL

  for (f in 1:length(features)) {
    message(paste("Pathway number:", f))
    com<-unique(features[[f]])

    common_features <- intersect(com, intersect(colnames(train_data), colnames(test_data)))
    if (!length(common_features)) {
      stop("No shared training/test features for pathway ", f, ".")
    }
    data_train <- train_data[, common_features, drop = FALSE]
    data_test  <- test_data[,  common_features, drop = FALSE]

    xtrain <- data_train
    xtest <- data_test

    for (model_name in models.indiv) {
      if (toupper(model_name) == "XG") {
        assign(paste0("y.train_tr.XG_",f), rep(NA, length(outcome)))
      }
      if (toupper(model_name) == "EN") {
        assign(paste0("y.train_tr.EN_",f), rep(NA, length(outcome)))
      }
      if (toupper(model_name) == "RF") {
        assign(paste0("y.train_tr.RF_",f), rep(NA, length(outcome)))
      }
      if (toupper(model_name) == "ADB") {
        assign(paste0("y.train_tr.ADB_",f), rep(NA, length(outcome)))
      }
      if (toupper(model_name) == "MB") {
        assign(paste0("y.train_tr.MB_",f), rep(NA, length(outcome)))
      }
      if (toupper(model_name) == "NN") {
        assign(paste0("y.train_tr.NN_",f), rep(NA, length(outcome)))
      }
      if (toupper(model_name) == "KNN") {
        assign(paste0("y.train_tr.KNN_",f), rep(NA, length(outcome)))
      }
      if (toupper(model_name) == "GB") {
        assign(paste0("y.train_tr.GB_",f), rep(NA, length(outcome)))
      }
      if (toupper(model_name) == "SVM") {
        assign(paste0("y.train_tr.SVM_",f), rep(NA, length(outcome)))
      }
      if (toupper(model_name) == "NB") {
        assign(paste0("y.train_tr.NB_",f), rep(NA, length(outcome)))
      }
      if (toupper(model_name) == "DCT") {
        assign(paste0("y.train_tr.DCT_",f), rep(NA, length(outcome)))
      }
    }
    for (model_name in models.indiv) {
      if (toupper(model_name) == "XG") {
        assign(paste0("y.train_ts.XG_",f), rep(NA, length(outcome)))
      }
      if (toupper(model_name) == "EN") {
        assign(paste0("y.train_ts.EN_",f), rep(NA, length(outcome)))
      }
      if (toupper(model_name) == "RF") {
        assign(paste0("y.train_ts.RF_",f), rep(NA, length(outcome)))
      }
      if (toupper(model_name) == "ADB") {
        assign(paste0("y.train_ts.ADB_",f), rep(NA, length(outcome)))
      }
      if (toupper(model_name) == "MB") {
        assign(paste0("y.train_ts.MB_",f), rep(NA, nrow(test_data)))
      }
      if (toupper(model_name) == "NN") {
        assign(paste0("y.train_ts.NN_",f), rep(NA, length(outcome)))
      }
      if (toupper(model_name) == "KNN") {
        assign(paste0("y.train_ts.KNN_",f), rep(NA, length(outcome)))
      }
      if (toupper(model_name) == "GB") {
        assign(paste0("y.train_ts.GB_",f), rep(NA, length(outcome)))
      }
      if (toupper(model_name) == "SVM") {
        assign(paste0("y.train_ts.SVM_",f), rep(NA, length(outcome)))
      }
      if (toupper(model_name) == "NB") {
        assign(paste0("y.train_ts.NB_",f), rep(NA, length(outcome)))
      }
      if (toupper(model_name) == "DCT") {
        assign(paste0("y.train_ts.DCT_",f), rep(NA, length(outcome)))
      }
    }
    for (model_name in models.indiv) {
      if (toupper(model_name) == "XG") {
        message("XgBoost running...")
        params <- modifyList(indiv.defaults$XG, param.indiv[["XG"]] %||% list())
        dtrain <- xgboost::xgb.DMatrix(data = as.matrix(xtrain), label = as.numeric(outcome))
        params$data <- NULL
        params$label <- NULL
        params$objective <- "binary:logistic"
        temp <- do.call(
          xgboost::xgb.train,
          list(
            data    = dtrain,
            params  = params,
            nrounds = params$nrounds %||% 500
          )
        )
        params$x <- NULL
        params$y <- NULL
        assign(paste0("y.train_tr.XG_",f), predict(temp, dtrain))
        assign(paste0("y.train_ts.XG_",f), predict(temp, newdata=xtest))
      }
      if (toupper(model_name) == "RF") {
        message("RandomForest running...")
        params <- modifyList(indiv.defaults$RF, param.indiv[["RF"]] %||% list())
        y <- factor(outcome, levels = c(0,1))
        params$x <- as.data.frame(xtrain)
        params$y <- y
        temp <- do.call(randomForest, params)
        assign(paste0("y.train_tr.RF_",f), predict(temp, newdata=xtrain, type="prob")[,2])
        assign(paste0("y.train_ts.RF_",f), predict(temp, newdata=xtest, type="prob")[,2])
      }
      if (toupper(model_name) == "EN") {
        message("ElasticNet running...")
        en.params <- modifyList(indiv.defaults$EN, param.indiv[["EN"]] %||% list())
        train.params <- en.params[setdiff(names(en.params), "s")]
        temp <- do.call(cv.glmnet, c(list(x = as.matrix(xtrain), y = outcome), train.params))
        s.val <- param.indiv[["EN"]]$s %||% "lambda.1se"   # user-specified or default
        assign(paste0("y.train_tr.EN_",f),predict(temp, newx=as.matrix(xtrain), s=s.val, type="response"))
        assign(paste0("y.train_ts.EN_",f),predict(temp, newx=as.matrix(xtest), s=s.val, type="response"))
      }
      if (toupper(model_name) == "ADB") {
        message("AdaBoost running...")
        params <- modifyList(indiv.defaults$ADB, param.indiv[["ADB"]] %||% list())
        y_factor <- factor(outcome, levels = c(0,1))
        ada.sr <- ifelse(outcome == 0, -1, 1)
        params$formula <- NULL
        params$x <- NULL
        params$y <- NULL
        temp <- do.call(adaboost, c(list(xtrain, ada.sr), params))
        assign(paste0("y.train_tr.ADB_",f), predict(temp, xtrain, type = "prob"))
        assign(paste0("y.train_ts.ADB_",f), predict(temp, xtest, type = "prob"))
      }
      if (toupper(model_name) == "MB") {
        message("MBoost running...")
        params <- modifyList(indiv.defaults$MB, param.indiv[["MB"]] %||% list())
        assign(paste0("y.train_tr.MB_",f), do.call(
          binary_mb_predict,
          c(list(xtrain = xtrain, outcome = outcome), params)
        ))
        assign(paste0("y.train_ts.MB_",f), do.call(
          binary_mb_predict,
          c(list(xtrain = xtrain, outcome = outcome, xtest = xtest), params)
        ))
      }
      if (toupper(model_name) == "GB") {
        message("GradientBoost running...")
        params <- modifyList(indiv.defaults$GB, param.indiv[["GB"]] %||% list())
        train.df <- data.frame(y = outcome, xtrain)
        params$formula <- y ~.
        params$data <- train.df
        temp <- do.call(gbm::gbm, params)
        params$x <- NULL
        params$y <- NULL
        assign(paste0("y.train_tr.GB_",f), predict(temp, data.frame(xtrain), n.trees = 1000, type = "response"))
        assign(paste0("y.train_ts.GB_",f), predict(temp, data.frame(xtest), n.trees = 1000, type = "response"))
      }
      if (toupper(model_name) == "KNN") {
        message("KNN running...")
        train.df <- data.frame(y = factor(outcome, levels = c(0,1)), xtrain)
        y.train.factor <- factor(outcome, levels = c(0,1))
        test.df  <- data.frame(xtest)
        colnames(test.df) <- colnames(xtrain)  # ensure names match
        knn.params <- modifyList(indiv.defaults$KNN, param.indiv[["KNN"]] %||% list())
        knn.params$x <- as.matrix(xtrain)
        knn.params$y <- y.train.factor
        temp <- do.call(knn3, knn.params)
        assign(paste0("y.train_tr.KNN_",f), predict(temp, data.frame(xtrain),type = "prob")[,2])
        assign(paste0("y.train_ts.KNN_",f), predict(temp, data.frame(xtest),type = "prob")[,2])
      }
      if (toupper(model_name) == "SVM") {
        message("SVM running...")
        temp<-svm(outcome ~ ., data=xtrain)
        assign(paste0("y.train_tr.SVM_",f), predict(temp, newdata=xtrain, probability = TRUE))
        assign(paste0("y.train_ts.SVM_",f), predict(temp, newdata=xtest, probability = TRUE))
      }
      if (toupper(model_name) == "NB") {
        message("NB running...")
        params <- modifyList(indiv.defaults$NB %||% list(), param.indiv[["NB"]] %||% list())
        usekernel_flag <- params$usekernel %||% TRUE
        params$usekernel <- NULL
        y1<-as.factor(as.character(as.numeric(outcome)))
        dat.nb<- cbind.data.frame("class" = y1, xtrain)
        params$formula <- class ~.
        params$data <- dat.nb
        params$x <- NULL
        params$y <- NULL
        temp <- do.call(naive_bayes, c(params, list(usekernel = usekernel_flag)))
        assign(paste0("y.train_tr.NB_",f), predict(temp, data.frame(xtrain), type="prob")[,2])
        assign(paste0("y.train_ts.NB_",f), predict(temp, data.frame(xtest), type="prob")[,2])
      }
      if (toupper(model_name) == "DCT") {
        message("DecisionTree running...")
        y_outcome <- factor(outcome, levels = c(0,1))
        train.df <- data.frame(y = y_outcome, xtrain)
        params <- modifyList(indiv.defaults$DCT %||% list(), param.indiv[["DCT"]] %||% list())
        params$x <- NULL
        params$y <- NULL
        params$formula <- y ~.
        params$data <- train.df
        params$method <- "class"
        temp <- do.call(rpart, params)
        assign(paste0("y.train_tr.DCT_",f), predict(temp, newdata=data.frame(xtrain), type="prob")[,2])
        assign(paste0("y.train_ts.DCT_",f), predict(temp, newdata=data.frame(xtest), type="prob")[,2])
      }
      if (toupper(model_name) == "NN") {
        message("NeuralNet running...")
        nn.defaults <- indiv.defaults$NN
        nn.overrides <- param.indiv[["NN"]] %||% list()

        params <- modifyList(nn.defaults, nn.overrides)

        # HARD reset anything neuralnet might cache
        params$x <- NULL
        params$y <- NULL
        params$data <- NULL
        params$formula <- NULL
        params$covariate <- NULL
        params$response <- NULL
        # Force data frames
        xtrain <- as.data.frame(xtrain)
        xtest  <- as.data.frame(xtest)

        # Keep only common features
        common <- intersect(colnames(xtrain), colnames(xtest))
        xtrain <- xtrain[, common, drop = FALSE]
        xtest  <- xtest[,  common, drop = FALSE]

        # Drop zero-variance columns (important!)
        vars.ok <- apply(xtrain, 2, sd) > 0
        xtrain <- xtrain[, vars.ok, drop = FALSE]
        xtest  <- xtest[,  vars.ok, drop = FALSE]

        # Scale
        m <- colMeans(xtrain)
        s <- apply(xtrain, 2, sd)
        xtrain.sd <- as.data.frame(scale(xtrain, center = m, scale = s))
        xtest.sd  <- as.data.frame(scale(xtest,  center = m, scale = s))

        # Train data frame
        train_df <- cbind(y = as.numeric(outcome), xtrain.sd)
        colnames(train_df) <- make.names(colnames(train_df), unique = TRUE)

        # Fit neural network
        params$formula <- y ~ .
        params$data <- train_df
        temp.nn <- do.call(neuralnet::neuralnet, params)

        # Predict
        assign(paste0("y.train_tr.NN_", f), predict(temp.nn, newdata = xtrain.sd))
        assign(paste0("y.train_ts.NN_", f), predict(temp.nn, newdata = xtest.sd))

      }
    }

    foo <- do.call(cbind, lapply(models.indiv, function(m) {
      get(paste0("y.train_tr.", m, "_", f))
    }))

    foo[foo < 0] <- 0
    foo[foo > 1] <- 1

    foo2 <- do.call(cbind, lapply(models.indiv, function(m) {
      get(paste0("y.train_ts.", m, "_", f))
    }))

    foo2[foo2 < 0] <- 0
    foo2[foo2 > 1] <- 1

    mat.train<- cbind(mat.train, foo)

    mat.test<- cbind(mat.test, foo2)

    # Weighted
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

  ########## ensemble of methods ###########
  colnames(mat.train)<- colnames(mat.test)<- colnames(y.pred.train.all)<- colnames(y.pred.mat.all)<- paste0("M",1:(length(models.indiv)*length(features)))

  for (model_name in models.ens) {
    if (toupper(model_name) == "XG") {
      message("XgBoost ensemble..")
      params <- modifyList(ens.defaults$XG, param.ens[["XG"]] %||% list())
      dtrain <- xgboost::xgb.DMatrix(data = as.matrix(mat.train), label = as.numeric(outcome))
      params$data <- NULL
      params$label <- NULL
      if (!is.null(params$eta)) {
        params$learning_rate <- params$eta
        params$eta <- NULL
      }
      params$verbose <- NULL
      params$objective <- params$objective %||% "binary:logistic"
      params$verbosity <- params$verbosity %||% 0
      nrounds <- params$nrounds %||% 80
      params$nrounds <- NULL
      temp.XG <- xgboost::xgb.train(params = params, data = dtrain, nrounds = nrounds)
    }
    if (toupper(model_name) == "RF") {
      message("RandomForest running...")
      params <- modifyList(indiv.defaults$RF, param.indiv[["RF"]] %||% list())
      y <- factor(outcome, levels = c(0,1))
      params$x <- as.data.frame(mat.train)
      params$y <- y
      temp.RF <- do.call(randomForest, params)
    }
    if (toupper(model_name) == "EN") {
      message("ElasticNet running...")
      en.params <- modifyList(indiv.defaults$EN, param.indiv[["EN"]] %||% list())
      train.params <- en.params[setdiff(names(en.params), "s")]
      temp.EN <- do.call(cv.glmnet, c(list(x = as.matrix(mat.train), y = outcome), train.params))
    }
    if (toupper(model_name) == "ADB") {
      message("AdaBoost running...")
      params <- modifyList(indiv.defaults$ADB, param.indiv[["ADB"]] %||% list())
      y_factor <- factor(outcome, levels = c(0,1))
      ada.sr <- ifelse(outcome == 0, -1, 1)
      params$formula <- NULL
      params$x <- NULL
      params$y <- NULL
      temp.ADB <- do.call(adaboost, c(list(mat.train, ada.sr), params))
    }
    if (toupper(model_name) == "MB") {
      message("MBoost running...")
      params <- modifyList(ens.defaults$MB, param.ens[["MB"]] %||% list())
      temp.MB.params <- params
    }
    if (toupper(model_name) == "GB") {
      message("GradientBoost running...")
      params <- modifyList(indiv.defaults$GB, param.indiv[["GB"]] %||% list())
      train.df <- data.frame(y = outcome, mat.train)
      params$formula <- y ~.
      params$data <- train.df
      temp.GB <- do.call(gbm::gbm, params)
      params$x <- NULL
      params$y <- NULL
    }
    if (toupper(model_name) == "KNN") {
      message("KNN running...")
      train.df <- data.frame(y = factor(outcome, levels = c(0,1)), mat.train)
      y.train.factor <- factor(outcome, levels = c(0,1))
      knn.params <- modifyList(indiv.defaults$KNN, param.indiv[["KNN"]] %||% list())
      knn.params$x <- as.matrix(mat.train)
      knn.params$y <- y.train.factor
      temp.KNN <- do.call(knn3, knn.params)
    }
    if (toupper(model_name) == "SVM") {
      message("SVM running...")
      temp.SVM<-svm(outcome ~ ., data=mat.train)
    }
    if (toupper(model_name) == "NB") {
      message("NB running...")
      params <- modifyList(indiv.defaults$NB %||% list(), param.indiv[["NB"]] %||% list())
      usekernel_flag <- params$usekernel %||% TRUE
      params$usekernel <- NULL
      y1<-as.factor(as.character(as.numeric(outcome)))
      dat.nb<- cbind.data.frame("class" = y1, mat.train)
      params$formula <- class ~.
      params$data <- dat.nb
      params$x <- NULL
      params$y <- NULL
      temp.NB <- do.call(naive_bayes, c(params, list(usekernel = usekernel_flag)))
    }
    if (toupper(model_name) == "DCT") {
      message("DecisionTree running...")
      y_outcome <- factor(outcome, levels = c(0,1))
      train.df <- data.frame(y = y_outcome, mat.train)
      params <- modifyList(indiv.defaults$DCT %||% list(), param.indiv[["DCT"]] %||% list())
      params$x <- NULL
      params$y <- NULL
      params$formula <- y ~.
      params$data <- train.df
      params$method <- "class"
      temp.DCT <- do.call(rpart, params)
    }
    if (toupper(model_name) == "NN") {
      message("NeuralNet running...")
      params <- modifyList(indiv.defaults$NN, param.indiv[["NN"]] %||% list())
      params$x <- NULL
      params$y <- NULL
      mat.train.nn <- as.data.frame(mat.train)
      vars.ok <- apply(mat.train.nn, 2, sd) > 0
      mat.train.nn <- mat.train.nn[, vars.ok, drop = FALSE]
      m <- colMeans(mat.train.nn)
      s <- apply(mat.train.nn, 2, sd)
      xtrain.sd <- as.data.frame(scale(mat.train.nn, center = m, scale = s))
      train_df <- cbind(y = as.numeric(outcome), xtrain.sd)
      colnames(train_df) <- make.names(colnames(train_df), unique = TRUE)
      params$formula <- y ~ .
      params$data <- train_df
      temp.NN <- do.call(neuralnet::neuralnet, params)
    }
  }

  ########### weighted avg ###############

  ind<- which(apply(y.pred.train.all,2,sd) == 0)
  if(length(ind) > 0){
    y.pred.train.all <- y.pred.train.all[,-ind]
    y.pred.mat.all <- y.pred.mat.all[,-ind]
  }

  y.pred.train.all<- cbind.data.frame(outcome = outcome, y.pred.train.all)
  temp.avg.all<- glm(outcome ~  ., family=binomial(link = "logit"), data = y.pred.train.all)

  ########### simple avg ###############
  y.pred.avg1.all<-apply(mat.test,1,mean)

  ########## predict on Frozen using ensemble ###########
  if (ensemble == TRUE) {
    for (model_name in models.ens) {
      if (toupper(model_name) == "XG") {
        dtest <- xgboost::xgb.DMatrix(data = as.matrix(mat.test))
        y.ens.XG=predict(temp.XG, dtest)
      }
      if (toupper(model_name) == "RF") {
        message("RandomForest running...")
        y.ens.RF=predict(temp.RF, newdata=mat.test, type="prob")[,2]
      }
      if (toupper(model_name) == "EN") {
        message("ElasticNet running...")
        s.val <- param.indiv[["EN"]]$s %||% "lambda.1se"   # user-specified or default
        y.ens.EN=predict(temp.EN, newx=mat.test, type="response",s=s.val)
      }
      if (toupper(model_name) == "ADB") {
        message("AdaBoost running...")
        y.ens.ADB=predict(temp.ADB, mat.test, type="prob")
      }
      if (toupper(model_name) == "MB") {
        message("MBoost running...")
        y.ens.MB <- do.call(
          binary_mb_predict,
          c(list(xtrain = mat.train, outcome = outcome, xtest = mat.test), temp.MB.params)
        )
      }
      if (toupper(model_name) == "GB") {
        message("GradientBoost running...")
        y.ens.GB=predict.gbm(temp.GB, data.frame(mat.test), type = "response")
      }
      if (toupper(model_name) == "KNN") {
        message("KNN running...")
        y.ens.KNN=predict(temp.KNN, data.frame(mat.test),type = "prob")[,2]
      }
      if (toupper(model_name) == "SVM") {
        message("SVM running...")
        y.ens.SVM=predict(temp.SVM, newdata=mat.test, probability = TRUE)
      }
      if (toupper(model_name) == "NB") {
        message("NB running...")
        y.ens.NB <- predict(temp.NB, data.frame(mat.test), type="prob")[,2]
      }
      if (toupper(model_name) == "DCT") {
        message("DecisionTree running...")
        y.ens.DCT <- predict(temp.DCT, newdata=data.frame(mat.test), type="prob")[,2]
      }
      if (toupper(model_name) == "NN") {
        message("NeuralNet running...")
        m <- colMeans(as.data.frame(mat.train))
        s <- apply(as.data.frame(mat.train), 2, sd)
        xtest.sd <- as.data.frame(scale(mat.test, center = m, scale = s))
        colnames(xtest.sd) <- colnames(mat.test)
        y.ens.NN=predict(temp.NN, newdata=mat.test, type="prob")
      }
    }
  }
  y.pred.avg.all<- predict(temp.avg.all, newdata=as.data.frame(mat.test), type="response")

  # Create storage vectors
  auc.tab <- c()
  ci.low.tab <- c()
  ci.up.tab <- c()
  sp.95.tab <- c()
  sp.98.tab <- c()
  model.names <- c()

  # Loop through ensemble models only
  if (ensemble == TRUE) {
    for (m in models.ens) {
      pred <- as.numeric(get(paste0("y.ens.", m)))
      roc.obj <- roc(test_outcome, pred, ci = TRUE, direction = "<")

      auc.tab <- c(auc.tab, as.numeric(roc.obj$auc))
      ci.low.tab <- c(ci.low.tab, roc.obj$ci[1])
      ci.up.tab <- c(ci.up.tab, roc.obj$ci[3])

      sp.95.tab <- c(sp.95.tab, as.numeric(coords(roc.obj, x = 0.95, input = "specificity", ret = "sensitivity")))
      sp.98.tab <- c(sp.98.tab, as.numeric(coords(roc.obj, x = 0.98, input = "specificity", ret = "sensitivity")))

      model.names <- c(model.names, paste0(m, ".ens"))
    }
  }

  # Add ensemble-level models (weighted and simple averages)
  roc.all  <- roc(test_outcome, y.pred.avg.all,  ci = TRUE, direction = "<")
  roc.all1 <- roc(test_outcome, y.pred.avg1.all, ci = TRUE, direction = "<")

  for (name in c("WeightedAvg", "SimpleAvg")) {
    roc.obj <- if (name == "WeightedAvg") roc.all else roc.all1
    auc.tab <- c(auc.tab, as.numeric(roc.obj$auc))
    ci.low.tab <- c(ci.low.tab, roc.obj$ci[1])
    ci.up.tab <- c(ci.up.tab, roc.obj$ci[3])

    sp.95.tab <- c(sp.95.tab, as.numeric(coords(roc.obj, x = 0.95, input = "specificity", ret = "sensitivity")))
    sp.98.tab <- c(sp.98.tab, as.numeric(coords(roc.obj, x = 0.98, input = "specificity", ret = "sensitivity")))

    model.names <- c(model.names, name)
  }

  get_marker <- function(model_name) {
    if (identical(model_name, "WeightedAvg")) return(y.pred.avg.all)
    if (identical(model_name, "SimpleAvg")) return(y.pred.avg1.all)
    base_name <- sub("\\.ens$", "", model_name)
    get(paste0("y.ens.", base_name), inherits = TRUE)
  }
  sp95_ci <- t(vapply(
    model.names,
    function(model_name) binary_sp_ci(test_outcome, get_marker(model_name)),
    numeric(2)
  ))

  # Combine results
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

####### Categorical ######
cat_test_performance <- function(outcome, train_data, test_data, test_outcome, features, models.indiv, ensemble = FALSE, models.ens = NULL, param.indiv = NULL, param.ens = NULL) {

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

  mat.train <- matrix(nrow = length(outcome), ncol = 0)
  mat.test  <- matrix(nrow = nrow(test_data),  ncol = 0)
  y.pred.train.all<- matrix(nrow = length(outcome), ncol = 0)
  y.pred.mat.all<- matrix(nrow = nrow(test_data),  ncol = 0)

  num_classes <- length(unique(outcome))

  for (f in seq_along(features)) {
    message(paste("Pathway number:", f))
    com<-unique(features[[f]])

    common_features <- intersect(com, intersect(colnames(train_data), colnames(test_data)))
    if (!length(common_features)) {
      stop("No shared training/test features for pathway ", f, ".")
    }
    data_train <- train_data[, common_features, drop = FALSE]
    data_test  <- test_data[,  common_features, drop = FALSE]

    xtrain <- data_train
    xtest <- data_test

    for (model_name in models.indiv) {
      if (toupper(model_name) == "XG") {
        y.train_tr.XG <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
        y.train_ts.XG <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
      }
      if (toupper(model_name) == "EN") {
        y.train_tr.EN <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
        y.train_ts.EN <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
      }
      if (toupper(model_name) == "RF") {
        y.train_tr.RF <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
        y.train_ts.RF <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
      }
      if (toupper(model_name) == "ADB") {
        y.train_tr.ADB <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
        y.train_ts.ADB <- matrix(NA_real_, nrow = nrow(test_data), ncol = length(unique(outcome)))
      }
      if (toupper(model_name) == "MB") {
        y.train_tr.MB <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
        y.train_ts.MB <- matrix(NA_real_, nrow = nrow(test_data), ncol = length(unique(outcome)))
      }
      if (toupper(model_name) == "GB") {
        y.train_tr.GB <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
        y.train_ts.GB <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
      }
      if (toupper(model_name) == "SVM") {
        y.train_tr.SVM <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
        y.train_ts.SVM <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
      }
      if (toupper(model_name) == "NB") {
        y.train_tr.NB <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
        y.train_ts.NB <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
      }
      if (toupper(model_name) == "KNN") {
        y.train_tr.KNN <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
        y.train_ts.KNN <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
      }
      if (toupper(model_name) == "DCT") {
        y.train_tr.DCT <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
        y.train_ts.DCT <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
      }
      if (toupper(model_name) == "NN") {
        y.train_tr.NN <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
        y.train_ts.NN <- matrix(NA_real_, nrow = length(outcome), ncol = length(unique(outcome)))
      }
    }

    xtrain_cl <- as.data.frame(xtrain)
    xtrain_cl$cluster <- factor(paste0("C", outcome))

    ytrain <- outcome

    xtest_cl <- as.data.frame(xtest)

    for (model_name in models.indiv) {
      num_classes <- length(unique(outcome))
      if (toupper(model_name) == "XG") {
        message("XgBoost running...")
        params_xg <- modifyList(indiv.defaults$XG, param.indiv[["XG"]] %||% list())
        num_classes <- length(levels(outcome))
        temp <- do.call(xgboost::xgboost, c(
          list(
            x = as.matrix(xtrain),
            y = as.factor(ytrain)
          ),
          params_xg
        ))
        pred_probs <- predict(temp, newdata = as.matrix(xtest))
        pred_matrix <- matrix(pred_probs, nrow = nrow(xtest), ncol = num_classes)
        y.train_ts.XG <- pred_matrix

        pred_probs <- predict(temp, newdata = as.matrix(xtrain))
        pred_matrix <- matrix(pred_probs, nrow = nrow(xtrain), ncol = num_classes)
        y.train_tr.XG <- pred_matrix
      }
      if (toupper(model_name) == "RF") {
        message("RandomForest running...")
        params_rf <- modifyList(indiv.defaults$RF, param.indiv[["RF"]] %||% list())
        temp <- do.call(randomForest, c(list(x = xtrain, y = as.factor(outcome)), params_rf))
        y.train_ts.RF <- predict(temp, newdata = xtest, type = "prob")
        y.train_tr.RF <- predict(temp, newdata = xtrain, type = "prob")
      }
      if (toupper(model_name) == "EN") {
        message("ElasticNet running...")
        params_en <- modifyList(indiv.defaults$EN, param.indiv[["EN"]] %||% list())
        train.params <- params[setdiff(names(params_en), "s")]
        temp <- do.call(cv.glmnet, c(list(xtrain, as.factor(outcome), family = "multinomial"), train.params))
        s.val <- param.indiv[["EN"]]$s %||% "lambda.1se"
        y.train_tr.EN <- (predict(temp, newx = xtrain, s = s.val, type = "response"))
        y.train_tr.EN <- as.matrix(y.train_tr.EN[,,1])
        y.train_ts.EN <- predict(temp, newx = xtest, s = s.val, type = "response")
        y.train_ts.EN <- as.matrix(y.train_ts.EN[,,1])
      }
      if (toupper(model_name) == "ADB") {
        message("AdaBoost one-vs-rest running...")
        params_adb <- modifyList(indiv.defaults$ADB, param.indiv[["ADB"]] %||% list())
        y.train_tr.ADB <- do.call(
          multiclass_adb_predict,
          c(list(xtrain = xtrain, outcome = outcome,
                 classes = sort(unique(outcome))), params_adb)
        )
        y.train_ts.ADB <- do.call(
          multiclass_adb_predict,
          c(list(xtrain = xtrain, outcome = outcome, xtest = xtest,
                 classes = sort(unique(outcome))), params_adb)
        )
      }
      if (toupper(model_name) == "MB") {
        message("MBoost one-vs-rest running...")
        params_mb <- modifyList(indiv.defaults$MB, param.indiv[["MB"]] %||% list())
        y.train_tr.MB <- do.call(
          multiclass_mb_predict,
          c(list(xtrain = xtrain, outcome = outcome,
                 classes = sort(unique(outcome))), params_mb)
        )
        y.train_ts.MB <- do.call(
          multiclass_mb_predict,
          c(list(xtrain = xtrain, outcome = outcome, xtest = xtest,
                 classes = sort(unique(outcome))), params_mb)
        )
      }
      if (toupper(model_name) == "GB") {
        message("GradientBoost running...")
        params_gb <- modifyList(indiv.defaults$GB, param.indiv[["GB"]] %||% list())
        temp <- do.call(gbm.fit, c(list(xtrain, as.factor(outcome), distribution = "multinomial"), params_gb))
        y.train_tr.GB <- predict.gbm(temp, data.frame(xtrain), type = "response")
        y.train_tr.GB <- y.train_tr.GB[, , 1]
        y.train_ts.GB <- predict.gbm(temp, data.frame(xtest), type = "response")
        y.train_ts.GB <- y.train_ts.GB[, , 1]
      }
      if (toupper(model_name) == "KNN") {
        message("KNN running...")
        params_knn <- modifyList(indiv.defaults$KNN, param.indiv[["KNN"]] %||% list())
        params_knn$trControl <- NULL  # remove so not passed twice

        temp <- do.call(caret::train, c(
          list(
            cluster ~ .,
            data = xtrain_cl,
            method = "knn",
            trControl = trc,
            preProcess = c("center", "scale")
          ),
          params_knn
        ))
        y.train_tr.KNN <- as.matrix(predict(temp, newdata = xtrain_cl, type = "prob"))
        y.train_ts.KNN <- as.matrix(predict(temp, newdata = xtest_cl, type = "prob"))
      }
      if (toupper(model_name) == "SVM") {
        message("SVM running...")
        params_svm <- modifyList(indiv.defaults$SVM, param.indiv[["SVM"]] %||% list())
        temp <- do.call(svm, c(list(cluster ~ ., data = xtrain_cl, probability = TRUE), params_svm))

        pred <- predict(temp, newdata= xtest_cl, probability = TRUE)
        probs <- attr(pred, "probabilities")
        class_order <- paste0("C", seq_len(num_classes))
        y.train_ts.SVM <- probs[, class_order, drop = FALSE]

        pred <- predict(temp, newdata= xtrain_cl, probability = TRUE)
        y.train_tr.SVM <- attr(pred, "probabilities")[, class_order, drop = FALSE]
      }
      if (toupper(model_name) == "NB") {
        message("NB running...")
        params_nb <- modifyList(indiv.defaults$NB, param.indiv[["NB"]] %||% list())
        temp <- do.call(naiveBayes, c(list(cluster ~., data = xtrain_cl), params_nb))
        y.train_tr.NB <- predict(temp, xtrain_cl, type = "raw")
        y.train_ts.NB <- predict(temp, xtest_cl, type = "raw")
      }
      if (toupper(model_name) == "DCT") {
        message("DecisionTree running...")
        params_dct <- modifyList(indiv.defaults$DCT, param.indiv[["DCT"]] %||% list())
        temp <- do.call(rpart, c(list(cluster ~., data = xtrain_cl, method = "class"), params_dct))
        y.train_tr.DCT <- as.matrix(predict(temp, newdata= xtrain_cl, type="prob"))
        y.train_ts.DCT <- as.matrix(predict(temp, newdata= xtest_cl, type="prob"))
      }
      if (toupper(model_name) == "NN") {
        message("NeuralNet running...")
        m <- colMeans(as.data.frame(xtrain))
        s <- apply(as.data.frame(xtrain), 2, sd)
        xtrain.sd <- scale(xtrain, center = m, scale = s)
        xtest.sd  <- scale(xtest,  center = m, scale = s)

        y.mat  <- model.matrix(~ outcome - 1)  # one-hot encoding
        colnames(y.mat) <- paste0("class", 1:num_classes)
        train.df <- data.frame(xtrain.sd, y.mat)

        f <- as.formula(paste(paste(colnames(y.mat), collapse = " + "), "~ ."))
        params <- modifyList(indiv.defaults$NN, param.indiv[["NN"]] %||% list())
        temp.nn <- do.call(neuralnet, c(list(f, data = train.df, linear.output = FALSE), params))

        y.train_tr.NN <- predict(temp.nn, newdata = xtrain.sd)
        y.train_ts.NN <- predict(temp.nn, newdata = xtest.sd)
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

  ########## ensemble of methods ###########

  colnames(mat.train)<- colnames(mat.test)<- colnames(y.pred.train.all)<- colnames(y.pred.mat.all)<- paste0("M",1:(length(models.ens)*length(features)*num_classes))

  mat.train.ens_cl <- as.data.frame(mat.train)
  mat.train.ens_cl$cluster <- factor(paste0("C", outcome))

  mat.test.ens_cl <- as.data.frame(mat.test)

  if (any(grepl("XG", models.ens, ignore.case = TRUE))) {
    params_xg <- modifyList(ens.defaults$XG, param.ens[["XG"]] %||% list())
    temp.XG <- do.call(xgboost::xgboost, c(
      list(x = as.matrix(mat.train), y = as.factor(ytrain)),
      params_xg
    ))
  }
  if (any(grepl("RF", models.ens, ignore.case = TRUE))) {
    params_rf <- modifyList(ens.defaults$RF, param.ens[["RF"]] %||% list())
    temp.RF <- do.call(randomForest, c(list(x = mat.train, y = as.factor(outcome)), params_rf))
  }
  if (any(grepl("EN", models.ens, ignore.case = TRUE))) {
    params_en <- modifyList(ens.defaults$EN, param.ens[["EN"]] %||% list())
    train.params <- params_en[setdiff(names(params_en), "s")]
    temp.EN <- do.call(cv.glmnet, c(list(mat.train, as.factor(outcome), family = "multinomial"), train.params))
  }
  if (any(grepl("ADB", models.ens, ignore.case = TRUE))) {
    params_adb <- modifyList(ens.defaults$ADB, param.ens[["ADB"]] %||% list())
    y.ens.ADB <- do.call(
      multiclass_adb_predict,
      c(list(xtrain = mat.train, outcome = outcome, xtest = mat.test,
             classes = sort(unique(outcome))), params_adb)
    )
  }
  if (any(grepl("MB", models.ens, ignore.case = TRUE))) {
    params_mb <- modifyList(ens.defaults$MB, param.ens[["MB"]] %||% list())
    y.ens.MB <- do.call(
      multiclass_mb_predict,
      c(list(xtrain = mat.train, outcome = outcome, xtest = mat.test,
             classes = sort(unique(outcome))), params_mb)
    )
  }
  if (any(grepl("GB", models.ens, ignore.case = TRUE))) {
    xtrain.mat <- as.matrix(mat.train)
    xtest.mat <- as.matrix(mat.test)
    num_features_ens <- ncol(xtrain.mat)
    feature_names <- paste0("EnsVar", 1:num_features_ens)
    colnames(xtrain.mat) <- feature_names
    colnames(xtest.mat) <- feature_names
    params_gb <- modifyList(ens.defaults$GB, param.ens[["GB"]] %||% list())
    temp.GB <- do.call(gbm.fit, c(list(x = xtrain.mat, y = outcome, distribution = "multinomial"), params_gb))
    n_trees_final <- temp.GB$n.trees
  }
  if (any(grepl("KNN", models.ens, ignore.case = TRUE))) {
    trc <- param.ens[["KNN"]]$trControl %||% trainControl(method = "cv", number = 3, classProbs = TRUE)
    params_knn <- modifyList(ens.defaults$KNN, param.ens[["KNN"]] %||% list())
    params_knn$trControl <- NULL
    temp.KNN <- do.call(caret::train, c(
      list(cluster ~ ., data = mat.train.ens_cl, method = "knn",
           trControl = trc, preProcess = c("center", "scale")),
      params_knn
    ))
  }
  if (any(grepl("SVM", models.ens, ignore.case = TRUE))) {
    params_svm <- modifyList(ens.defaults$SVM, param.ens[["SVM"]] %||% list())
    temp.SVM <- do.call(svm, c(list(cluster ~ ., data = data.frame(mat.train.ens_cl), probability = TRUE), params_svm))
  }
  if (any(grepl("NB", models.ens, ignore.case = TRUE))) {
    params_nb <- modifyList(ens.defaults$NB, param.ens[["NB"]] %||% list())
    temp.NB <- do.call(naiveBayes, c(list(cluster ~., data = mat.train.ens_cl), params_nb))
  }
  if (any(grepl("DCT", models.ens, ignore.case = TRUE))) {
    params_dct <- modifyList(ens.defaults$DCT, param.ens[["DCT"]] %||% list())
    temp.DCT <- do.call(rpart, c(list(cluster ~., data = mat.train.ens_cl, method = "class"), params_dct))
  }
  if (any(grepl("NN", models.ens, ignore.case = TRUE))) {
    m <- colMeans(as.data.frame(mat.train))
    s <- apply(as.data.frame(mat.train), 2, sd)
    xtrain.ens.sd <- scale(mat.train, center = m, scale = s)
    xtest.ens.sd  <- scale(mat.test,  center = m, scale = s)
    y.mat  <- model.matrix(~ outcome - 1)
    colnames(y.mat) <- paste0("class", 1:length(unique(outcome)))
    train.df <- data.frame(xtrain.ens.sd, y.mat)
    f <- as.formula(paste(paste(colnames(y.mat), collapse = " + "), "~ ."))
    params_nn <- modifyList(ens.defaults$NN, param.ens[["NN"]] %||% list())
    temp.NN <- do.call(neuralnet, c(list(f, data = train.df, linear.output = FALSE), params_nn))
  }

  ########### simple avg ###############
  N <- nrow(mat.test)
  K <- length(unique(outcome))
  M <- length(models.ens)

  y.pred.avg1.all <- multiclass_simple_average(mat.test, K)
  ########## predict on test set using ensemble ###########

  if (exists("temp.XG")) {
    pred_probs <- predict(temp.XG, newdata = as.matrix(mat.test))
    y.ens.XG <- matrix(pred_probs, nrow = nrow(mat.test), ncol = num_classes, byrow = TRUE)
  }
  if (exists("temp.RF")) y.ens.RF <- predict(temp.RF, newdata = mat.test, type = "prob")
  if (exists("temp.EN")) {
    s.val <- param.ens[["EN"]]$s %||% "lambda.1se"
    y.ens.EN <- predict(temp.EN, newx = mat.test, s = s.val, type = "response")
  }
  if (exists("temp.GB")) {
    y.pred.array <- predict.gbm(temp.GB, newdata = as.data.frame(xtest.mat), n.trees = n_trees_final, type = "response")
    y.ens.GB <- y.pred.array[, , dim(y.pred.array)[3]]
  }
  if (exists("temp.KNN")) y.ens.KNN <- as.matrix(predict(temp.KNN, newdata = mat.test.ens_cl, type = "prob"))
  if (exists("temp.SVM")) {
    pred <- predict(temp.SVM, newdata= mat.test.ens_cl, probability = TRUE)
    probs <- attr(pred, "probabilities")
    class_order <- levels(mat.train.ens_cl$cluster)
    y.ens.SVM <- probs[, class_order, drop = FALSE]
  }
  if (exists("temp.NB")) y.ens.NB <- predict(temp.NB, mat.test.ens_cl, type = "raw")
  if (exists("temp.DCT")) y.ens.DCT <- as.matrix(predict(temp.DCT, newdata= as.data.frame(mat.test.ens_cl), type="prob"))
  if (exists("temp.NN")) y.ens.NN <- predict(temp.NN, newdata = xtest.ens.sd)
  ############### roc #############
  outcome.pool <- test_outcome
  model.names <- c()
  auc.tab <- c()
  ci.low.tab <- c()
  ci.up.tab <- c()
  #sp.98.tab <- c()
  #sp.95.tab <- c()

  for (model_name in models.ens) {
    num_classes <- length(unique(outcome))
    if (grepl("XG", model_name, ignore.case = TRUE)) {
      XG_mat <- matrix(y.ens.XG, ncol = num_classes, byrow = TRUE)
      colnames(XG_mat) <- levels(outcome)
      roc.XG <- multiclass.roc(outcome.pool, XG_mat)
      roc.XG$ci <- get_multiclass_auc_ci(outcome.pool, XG_mat, nboot = 200)
      model.names <- c(model.names, "XG.ens")
      ci.low.tab <- c(ci.low.tab, roc.XG$ci$ci_lower)
      ci.up.tab <- c(ci.up.tab, roc.XG$ci$ci_upper)
      auc.tab <- c(auc.tab, roc.XG$auc)
    }
    if (grepl("RF", model_name, ignore.case = TRUE)) {
      RF_mat <- matrix(y.ens.RF, ncol = num_classes, byrow = TRUE)
      colnames(RF_mat) <- levels(outcome)
      roc.RF <- multiclass.roc(outcome.pool, RF_mat)
      roc.RF$ci <- get_multiclass_auc_ci(outcome.pool, RF_mat, nboot = 200)
      model.names <- c(model.names, "RF.ens")
      ci.low.tab <- c(ci.low.tab, roc.RF$ci$ci_lower)
      ci.up.tab <- c(ci.up.tab, roc.RF$ci$ci_upper)
      auc.tab <- c(auc.tab, roc.RF$auc)
    }
    if (grepl("EN", model_name, ignore.case = TRUE)) {
      EN_mat <- matrix(y.ens.EN, ncol = num_classes, byrow = TRUE)
      colnames(EN_mat) <- levels(outcome)
      roc.EN <- multiclass.roc(outcome.pool, EN_mat)
      roc.EN$ci <- get_multiclass_auc_ci(outcome.pool, EN_mat, nboot = 200)
      model.names <- c(model.names, "EN.ens")
      ci.low.tab <- c(ci.low.tab, roc.EN$ci$ci_lower)
      ci.up.tab <- c(ci.up.tab, roc.EN$ci$ci_upper)
      auc.tab <- c(auc.tab, roc.EN$auc)
    }
    if (grepl("ADB", model_name, ignore.case = TRUE)) {
      ADB_mat <- matrix(y.ens.ADB, ncol = num_classes, byrow = TRUE)
      colnames(ADB_mat) <- levels(outcome)
      roc.ADB <- multiclass.roc(outcome.pool, ADB_mat)
      roc.ADB$ci <- get_multiclass_auc_ci(outcome.pool, ADB_mat, nboot = 200)
      model.names <- c(model.names, "ADB.ens")
      ci.low.tab <- c(ci.low.tab, roc.ADB$ci$ci_lower)
      ci.up.tab <- c(ci.up.tab, roc.ADB$ci$ci_upper)
      auc.tab <- c(auc.tab, roc.ADB$auc)
    }
    if (grepl("MB", model_name, ignore.case = TRUE)) {
      MB_mat <- matrix(y.ens.MB, ncol = num_classes, byrow = TRUE)
      colnames(MB_mat) <- levels(outcome)
      roc.MB <- multiclass.roc(outcome.pool, MB_mat)
      roc.MB$ci <- get_multiclass_auc_ci(outcome.pool, MB_mat, nboot = 200)
      model.names <- c(model.names, "MB.ens")
      ci.low.tab <- c(ci.low.tab, roc.MB$ci$ci_lower)
      ci.up.tab <- c(ci.up.tab, roc.MB$ci$ci_upper)
      auc.tab <- c(auc.tab, roc.MB$auc)
    }
    if (grepl("NN", model_name, ignore.case = TRUE)) {
      NN_mat <- matrix(y.ens.NN, ncol = num_classes, byrow = TRUE)
      colnames(NN_mat) <- levels(outcome)
      roc.NN <- multiclass.roc(outcome.pool, NN_mat)
      roc.NN$ci <- get_multiclass_auc_ci(outcome.pool, NN_mat, nboot = 200)
      model.names <- c(model.names, "NN.ens")
      ci.low.tab <- c(ci.low.tab, roc.NN$ci$ci_lower)
      ci.up.tab <- c(ci.up.tab, roc.NN$ci$ci_upper)
      auc.tab <- c(auc.tab, roc.NN$auc)
    }
    if (grepl("GB", model_name, ignore.case = TRUE)) {
      GB_mat <- matrix(y.ens.GB, ncol = num_classes, byrow = TRUE)
      colnames(GB_mat) <- levels(outcome)
      roc.GB <- multiclass.roc(outcome.pool, GB_mat)
      roc.GB$ci <- get_multiclass_auc_ci(outcome.pool, GB_mat, nboot = 200)
      model.names <- c(model.names, "GB.ens")
      ci.low.tab <- c(ci.low.tab, roc.GB$ci$ci_lower)
      ci.up.tab <- c(ci.up.tab, roc.GB$ci$ci_upper)
      auc.tab <- c(auc.tab, roc.GB$auc)
    }
    if (grepl("KNN", model_name, ignore.case = TRUE)) {
      KNN_mat <- matrix(y.ens.KNN, ncol = num_classes, byrow = TRUE)
      colnames(KNN_mat) <- levels(outcome)
      roc.KNN <- multiclass.roc(outcome.pool, KNN_mat)
      roc.KNN$ci <- get_multiclass_auc_ci(outcome.pool, KNN_mat, nboot = 200)
      model.names <- c(model.names, "KNN.ens")
      ci.low.tab <- c(ci.low.tab, roc.KNN$ci$ci_lower)
      ci.up.tab <- c(ci.up.tab, roc.KNN$ci$ci_upper)
      auc.tab <- c(auc.tab, roc.KNN$auc)
    }
    if (grepl("SVM", model_name, ignore.case = TRUE)) {
      SVM_mat <- matrix(y.ens.SVM, ncol = num_classes, byrow = TRUE)
      colnames(SVM_mat) <- levels(outcome)
      roc.SVM <- multiclass.roc(outcome.pool, SVM_mat)
      roc.SVM$ci <- get_multiclass_auc_ci(outcome.pool, SVM_mat, nboot = 200)
      model.names <- c(model.names, "SVM.ens")
      ci.low.tab <- c(ci.low.tab, roc.SVM$ci$ci_lower)
      ci.up.tab <- c(ci.up.tab, roc.SVM$ci$ci_upper)
      auc.tab <- c(auc.tab, roc.SVM$auc)
    }
    if (grepl("NB", model_name, ignore.case = TRUE)) {
      NB_mat <- matrix(y.ens.NB, ncol = num_classes, byrow = TRUE)
      colnames(NB_mat) <- levels(outcome)
      roc.NB <- multiclass.roc(outcome.pool, NB_mat)
      roc.NB$ci <- get_multiclass_auc_ci(outcome.pool, NB_mat, nboot = 200)
      model.names <- c(model.names, "NB.ens")
      ci.low.tab <- c(ci.low.tab, roc.NB$ci$ci_lower)
      ci.up.tab <- c(ci.up.tab, roc.NB$ci$ci_upper)
      auc.tab <- c(auc.tab, roc.NB$auc)
    }
    if (grepl("DCT", model_name, ignore.case = TRUE)) {
      DCT_mat <- matrix(y.ens.DCT, ncol = num_classes, byrow = TRUE)
      colnames(DCT_mat) <- levels(outcome)
      roc.DCT <- multiclass.roc(outcome.pool, DCT_mat)
      roc.DCT$ci <- get_multiclass_auc_ci(outcome.pool, DCT_mat, nboot = 200)
      model.names <- c(model.names, "DCT.ens")
      ci.low.tab <- c(ci.low.tab, roc.DCT$ci$ci_lower)
      ci.up.tab <- c(ci.up.tab, roc.DCT$ci$ci_upper)
      auc.tab <- c(auc.tab, roc.DCT$auc)
    }
  }

  all1.mat <- matrix(y.pred.avg1.all, ncol = num_classes, byrow = TRUE)
  colnames(all1.mat) <- levels(outcome)
  roc.all1 <- multiclass.roc(outcome.pool, all1.mat)
  roc.all1$ci <- get_multiclass_auc_ci(outcome.pool, all1.mat, nboot = 200)
  model.names <- c(model.names, "simple")
  ci.low.tab <- c(ci.low.tab, roc.all1$ci$ci_lower)
  ci.up.tab <- c(ci.up.tab, roc.all1$ci$ci_upper)
  auc.tab <- c(auc.tab, roc.all1$auc)

  get_marker <- function(model_name) {
    if (identical(model_name, "simple")) return(all1.mat)
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

############# Survival #############
surv_test_performance <- function(xtrain_time, outcome, xtest_time, xtest_status, time_point, train_data, test_data, features, models.indiv, ensemble = FALSE, models.ens = NULL, param.indiv = NULL, param.ens = NULL) {
  mat.train<- NULL
  mat.test<- NULL
  y.pred.train.all<- NULL
  y.pred.mat.all<- NULL
  xtrain_status <- outcome

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

  for (f in 1:length(features)) {
    com<-unique(features[[f]])
    #train_data = t(data1[,cur.train])
    #test_data = t(data1[,cur.test])
    common_features <- intersect(com, intersect(colnames(train_data), colnames(test_data)))
    if (!length(common_features)) {
      stop("No shared training/test features for pathway ", f, ".")
    }
    data_train <- train_data[, common_features, drop = FALSE]
    data_test  <- test_data[,  common_features, drop = FALSE]

    xtrain <- data_train
    xtest <- data_test

    for (model_name in models.indiv) {
      if (grepl("XG", model_name, ignore.case = TRUE)) {
        assign(paste0("y.train_tr.XG_",f), rep(NA, length(outcome)))
      }
      if (grepl("EN", model_name, ignore.case = TRUE)) {
        assign(paste0("y.train_tr.EN_",f), rep(NA, length(outcome)))
      }
      if (grepl("RF", model_name, ignore.case = TRUE)) {
        assign(paste0("y.train_tr.RF_",f), rep(NA, length(outcome)))
      }
      if (grepl("MB", model_name, ignore.case = TRUE)) {
        assign(paste0("y.train_tr.MB_",f), rep(NA, length(outcome)))
      }
      if (grepl("ADB", model_name, ignore.case = TRUE)) {
        assign(paste0("y.train_tr.ADB_",f), rep(NA, length(outcome)))
      }
      if (grepl("GB", model_name, ignore.case = TRUE)) {
        assign(paste0("y.train_tr.GB_",f), rep(NA, length(outcome)))
      }
      if (grepl("KNN", model_name, ignore.case = TRUE)) {
        assign(paste0("y.train_tr.KNN_",f), rep(NA, length(outcome)))
      }
      if (grepl("NN", model_name, ignore.case = TRUE)) {
        assign(paste0("y.train_tr.NN_",f), rep(NA, length(outcome)))
      }
      if (grepl("SVM", model_name, ignore.case = TRUE)) {
        assign(paste0("y.train_tr.SVM_",f), rep(NA, length(outcome)))
      }
      if (grepl("DCT", model_name, ignore.case = TRUE)) {
        assign(paste0("y.train_tr.DCT_",f), rep(NA, length(outcome)))
      }
    }
    for (model_name in models.indiv) {
      if (grepl("XG", model_name, ignore.case = TRUE)) {
        assign(paste0("y.train_ts.XG_",f), rep(NA, length(outcome)))
      }
      if (grepl("EN", model_name, ignore.case = TRUE)) {
        assign(paste0("y.train_ts.EN_",f), rep(NA, length(outcome)))
      }
      if (grepl("RF", model_name, ignore.case = TRUE)) {
        assign(paste0("y.train_ts.RF_",f), rep(NA, length(outcome)))
      }
      if (grepl("MB", model_name, ignore.case = TRUE)) {
        assign(paste0("y.train_ts.MB_",f), rep(NA, length(outcome)))
      }
      if (grepl("ADB", model_name, ignore.case = TRUE)) {
        assign(paste0("y.train_ts.ADB_",f), rep(NA, nrow(test_data)))
      }
      if (grepl("GB", model_name, ignore.case = TRUE)) {
        assign(paste0("y.train_ts.GB_",f), rep(NA, length(outcome)))
      }
      if (grepl("KNN", model_name, ignore.case = TRUE)) {
        assign(paste0("y.train_ts.KNN_",f), rep(NA, nrow(test_data)))
      }
      if (grepl("NN", model_name, ignore.case = TRUE)) {
        assign(paste0("y.train_ts.NN_",f), rep(NA, nrow(test_data)))
      }
      if (grepl("SVM", model_name, ignore.case = TRUE)) {
        assign(paste0("y.train_ts.SVM_",f), rep(NA, length(outcome)))
      }
      if (grepl("DCT", model_name, ignore.case = TRUE)) {
        assign(paste0("y.train_ts.DCT_",f), rep(NA, length(outcome)))
      }
    }

    for (model_name in models.indiv) {
      if (grepl("XG", model_name, ignore.case = TRUE)) {
        if (grepl("XG", model_name, ignore.case = TRUE)) {
          message("XgBoost running...")

          params_user <- modifyList(
            indiv.defaults$XG,
            param.indiv[["XG"]] %||% list()
          )

          # survival label
          label_train <- ifelse(xtrain_status == 1, xtrain_time, -xtrain_time)

          # validation split (within fold train)
          val_ind <- sample(seq_len(nrow(xtrain)), 0.1 * nrow(xtrain))

          x_train_mat <- as.matrix(xtrain[-val_ind, , drop = FALSE])
          x_val_mat   <- as.matrix(xtrain[val_ind, , drop = FALSE])

          y_train_sub <- label_train[-val_ind]
          y_val_sub   <- label_train[val_ind]

          # split parameters correctly
          nrounds_val <- params_user$nrounds %||% 1000
          params_user$nrounds <- NULL

          fixed_params <- list(
            objective = "survival:cox",
            eval_metric = "cox-nloglik"
          )

          final_params <- modifyList(fixed_params, params_user)

          # IMPORTANT: survXgboost expects matrices (NOT DMatrix)
          train_args <- list(
            data = x_train_mat,
            label = y_train_sub,
            watchlist = list(val = x_val_mat),
            params = final_params,
            nrounds = nrounds_val,
            verbose = 0
          )

          temp <- do.call(survXgboost::xgb.train.surv, train_args)
          assign(
            paste0("y.train_tr.XG_", f),
            predict(temp, newdata = as.matrix(xtrain), type = "risk")
          )

          assign(
            paste0("y.train_ts.XG_", f),
            predict(temp, newdata = as.matrix(xtest), type = "risk")
          )
        }
      }
      if (grepl("RF", model_name, ignore.case = TRUE)) {
        message("RandomForest running...")
        params_rf <- modifyList(indiv.defaults$RF, param.indiv[["RF"]] %||% list())
        surv_train_data <- data.frame(time = xtrain_time, status = xtrain_status, xtrain)
        surv_test_data <- data.frame(time = xtest_time, status = xtest_status, xtest)
        formula_args <- list(formula = Surv(time, status) ~., data = surv_train_data)
        temp <- do.call(randomForestSRC::rfsrc, c(formula_args, params_rf))
        predictions_train <- predict(temp, newdata = surv_train_data, type = "response")
        predictions_test <- predict(temp, newdata = data.frame(xtest), type = "response")

        assign(paste0("y.train_tr.RF_",f), predictions_train$predicted)
        assign(paste0("y.train_ts.RF_",f), predictions_test$predicted)

      }
      if (grepl("EN", model_name, ignore.case = TRUE)) {
        message("ElasticNet running...")
        params_en <- modifyList(indiv.defaults$EN, param.indiv[["EN"]] %||% list())
        surv_train_object <- Surv(time = xtrain_time, event = xtrain_status)
        s.val <- params_en$s %||% "lambda.min"
        train.params <- params_en[setdiff(names(params_en), "s")]
        cv_glmnet <- do.call(glmnet::cv.glmnet, c(
          list(x = as.matrix(xtrain), y = surv_train_object, family = "cox", nfolds = 10), train.params))
        risk_scores_train <- predict(cv_glmnet, newx = as.matrix(xtrain), s = s.val, type = "response")
        risk_scores_test <- predict(cv_glmnet, newx = as.matrix(xtest), s = s.val, type = "response")

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

        temp <- do.call(mboost::glmboost, c(fixed_args, params_mb))
        risk_scores_train <- predict(temp, newdata = as.matrix(xtrain), type = "link")
        risk_scores_test <- predict(temp, newdata = as.matrix(xtest), type = "link")

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
        temp <- do.call(gbm::gbm, c(fixed_args, params_gb))

        risk_scores_train <- predict(temp, newdata = surv_train_data, n.trees = n.trees_pred, type = "link")
        risk_scores_test <- predict(temp, newdata = surv_test_data, n.trees = n.trees_pred, type = "link")
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
        temp <- do.call(survivalsvm::survivalsvm, c(formula_args, params_svm))
        risk_scores_train <- predict(object = temp, newdata = data.frame(xtrain))
        risk_scores_test <- predict(object = temp, newdata = data.frame(xtest))

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

        temp <- do.call(rpart::rpart, c(formula_args, list(control = do.call(rpart::rpart.control, control_params))))

        risk_scores_train <- predict(temp, newdata = data.frame(xtrain))
        risk_scores_test <- predict(temp, newdata = data.frame(xtest))

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

  ########## ensemble of methods ###########
  colnames(mat.test)<- colnames(mat.train)<-colnames(y.pred.mat.all)<- colnames(y.pred.train.all)<- paste0("M",1:(length(models.indiv)*length(features)))

  K <- length(models.indiv) * length(features)
  ensemble_colnames <- c("time", "status",paste0("M", seq_len(K)))
  surv_ensemble_train_data <- data.frame(time = xtrain_time, status = outcome, mat.train, check.names = FALSE)
  surv_ensemble_test_data <- data.frame(time = xtest_time,status = xtest_status, mat.test, check.names = FALSE)
  colnames(surv_ensemble_train_data) <- ensemble_colnames
  colnames(surv_ensemble_test_data)  <- ensemble_colnames

  for (model_name in models.ens) {
    if (grepl("XG", model_name, ignore.case = TRUE)) {

      params_user <- modifyList(
        list(objective = "survival:cox",
             eval_metric = "cox-nloglik",
             eta = 0.1),
        param.ens[[model_name]] %||% list()
      )

      val_ind <- sample(seq_len(nrow(mat.train)), 0.1 * nrow(mat.train))

      label_train <- ifelse(xtrain_status == 1, xtrain_time, -xtrain_time)

      # split matrices
      xtrain_xgboost <- as.matrix(mat.train[-val_ind, , drop = FALSE])
      x_val_mat      <- as.matrix(mat.train[val_ind, , drop = FALSE])

      label_train_test <- label_train[-val_ind]
      label_val        <- label_train[val_ind]

      # nrounds handling (must be top-level, not in params)
      nrounds_val <- params_user$nrounds %||% 1000
      params_user$nrounds <- NULL

      xgboost_fit <- do.call(
        survXgboost::xgb.train.surv,
        list(
          data = xtrain_xgboost,
          label = label_train_test,
          watchlist = list(val = x_val_mat),
          params = params_user,
          nrounds = nrounds_val,
          verbose = 0
        )
      )
    }
    if (grepl("RF", model_name, ignore.case = TRUE)) {
      params_rf <- modifyList(list(ntree = 2000, max.depth = 10, importance = TRUE),param.ens[[model_name]] %||% list())
      temp.RF <- do.call(rfsrc, c(list(formula = Surv(time,status) ~ ., data = surv_ensemble_train_data), params_rf))
    }
    if (grepl("EN", model_name, ignore.case = TRUE)) {
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
      temp.EN <- glmnet(as.matrix(mat.train), surv_train_object, family = "cox", lambda = lambda_val, alpha = params_en$alpha %||% 0.1)
    }
    if (grepl("MB", model_name, ignore.case = TRUE)) {
      params_mb <- modifyList(list(mstop = 300), param.ens[[model_name]] %||% list())
      control_object <- boost_control(mstop = params_mb$mstop %||% 300)
      temp.MB <- glmboost(Surv(time, status) ~ ., data = surv_ensemble_train_data, family = CoxPH(), control = control_object)
    }
    if (grepl("ADB", model_name, ignore.case = TRUE)) {
      params_adb <- modifyList(ens.defaults$ADB, param.ens[[model_name]] %||% list())
    }
    if (grepl("GB", model_name, ignore.case = TRUE)) {
      params_gb <- modifyList(list(n.trees = 1000, interaction.depth = 10, shrinkage = 0.005, bag.fraction = 0.5), param.ens[[model_name]] %||% list())
      temp.GB <- do.call(gbm, c(list(formula = Surv(time, status) ~ ., data = surv_ensemble_train_data, distribution = "coxph"), params_gb))
    }
    if (grepl("KNN", model_name, ignore.case = TRUE)) {
      params_knn <- modifyList(list(k = 5), param.ens[[model_name]] %||% list())
    }
    if (grepl("NN", model_name, ignore.case = TRUE)) {
      params_nn <- modifyList(ens.defaults$NN, param.ens[[model_name]] %||% list())
    }
    if (grepl("SVM", model_name, ignore.case = TRUE)) {
      params_svm <- modifyList(list(gamma.mu = 0.5, opt.meth = "quadprog", kernel = "add_kernel"), param.ens[[model_name]] %||% list())
      temp.SVM <- do.call(survivalsvm, c(list(formula = Surv(time, status) ~ ., data = surv_ensemble_train_data, type = "regression"), params_svm))
    }
    if (grepl("DCT", model_name, ignore.case = TRUE)) {
      params_dct <- modifyList(list(minsplit = 20, maxdepth = 10, cp = 0.01), param.ens[[model_name]] %||% list())
      temp.DCT <- rpart(Surv(time, status) ~ ., data = surv_ensemble_train_data, method = "exp", control = do.call(rpart.control, params_dct))
    }
  }

  ########### weighted avg ###############
  train_stack <- as.data.frame(y.pred.train.all)
  test_stack <- as.data.frame(y.pred.mat.all)
  valid_rows <- apply(train_stack, 1, function(x) all(is.finite(x)))
  train_stack <- train_stack[valid_rows, , drop = FALSE]
  time.train <- xtrain_time[valid_rows]
  status.train <- outcome[valid_rows]

  zero_var <- sapply(train_stack, function(x) var(x, na.rm = TRUE) < 1e-12)
  if (any(zero_var)) {
    train_stack <- train_stack[, !zero_var, drop = FALSE]
    test_stack <- test_stack[, !zero_var, drop = FALSE]
  }

  if (ncol(train_stack) > 0) {
    mu <- colMeans(train_stack)
    sdv <- apply(train_stack, 2, sd)
    sdv[!is.finite(sdv) | sdv < 1e-8] <- 1
    train_stack <- as.data.frame(scale(train_stack, center = mu, scale = sdv))
    test_stack <- as.data.frame(scale(test_stack, center = mu, scale = sdv))
    train_stack$time <- time.train
    train_stack$status <- status.train
    predictors <- setdiff(colnames(train_stack), c("time", "status"))
    temp.avg.all <- survival::coxph(
      as.formula(paste("Surv(time, status) ~", paste(predictors, collapse = " + "))),
      data = train_stack,
      ties = "breslow",
      x = TRUE,
      model = TRUE
    )
    y.pred.avg.all <- as.numeric(predict(temp.avg.all, newdata = test_stack, type = "lp"))
  } else {
    y.pred.avg.all <- rep(NA_real_, nrow(test_data))
  }

  ########### simple avg ###############
  y.pred.avg1.all<-apply(mat.test,1,mean)

  ########## predict on test set using ensemble ###########
  for (model_name in models.ens) {
    if (grepl("XG", model_name, ignore.case = TRUE)) {
      y.ens.XG <- predict(xgboost_fit, newdata = as.matrix(mat.test), type = "risk")
    }
    if (grepl("RF", model_name, ignore.case = TRUE)) {
      y.ens.RF <- predict(temp.RF, newdata = surv_ensemble_test_data, type = "response")$predicted
    }
    if (grepl("EN", model_name, ignore.case = TRUE)) {
      y.ens.EN <- predict(temp.EN, newx = as.matrix(mat.test), s = lambda_val, type = "response")
    }
    if (grepl("MB", model_name, ignore.case = TRUE)) {
      y.ens.MB <- predict(temp.MB, newdata = surv_ensemble_test_data, type = "link")
    }
    if (grepl("ADB", model_name, ignore.case = TRUE)) {
      y.ens.ADB <- do.call(
        survival_adb_predict,
        c(list(xtrain = mat.train, time = xtrain_time, status = xtrain_status,
               time_point = time_point, xtest = mat.test), params_adb)
      )
    }
    if (grepl("GB", model_name, ignore.case = TRUE)) {
      y.ens.GB <- predict(temp.GB, newdata = surv_ensemble_test_data, n.trees = params_gb$n.trees %||% 1000, type = "link")
    }
    if (grepl("KNN", model_name, ignore.case = TRUE)) {
      y.ens.KNN <- survival_knn_predict(
        mat.train, xtrain_time, xtrain_status, mat.test, k = params_knn$k %||% 5
      )
    }
    if (grepl("NN", model_name, ignore.case = TRUE)) {
      y.ens.NN <- do.call(
        survival_nn_predict,
        c(list(xtrain = mat.train, time = xtrain_time, status = xtrain_status,
               xtest = mat.test), params_nn)
      )
    }
    if (grepl("SVM", model_name, ignore.case = TRUE)) {
      y.ens.SVM <- as.numeric(predict(temp.SVM, newdata = surv_ensemble_test_data)$predicted)
    }
    if (grepl("DCT", model_name, ignore.case = TRUE)) {
      y.ens.DCT <- unname(predict(temp.DCT, newdata = surv_ensemble_test_data))
    }
  }

  ############### roc #############
  auc.tab <- c()
  #ci.low.tab <- c()
  #ci.up.tab <- c()
  sp.95.tab <- c()
  sp.98.tab <- c()
  model.names <- c()

  for (model_name in models.ens) {
    if (grepl("XG", model_name, ignore.case = TRUE)) {
      roc.XG.ens <- survivalROC(Stime = xtest_time, status = xtest_status, marker = y.ens.XG, predict.time = time_point, method = "KM")
      auc.tab <- c(auc.tab, roc.XG.ens$AUC)
      sp.98.tab <- c(sp.98.tab, roc.XG.ens$TP[roc.XG.ens$FP <= 0.02][1])
      sp.95.tab <- c(sp.95.tab, roc.XG.ens$TP[roc.XG.ens$FP <= 0.05][1])
      model.names <- c(model.names, "XG")
    }
    if (grepl("RF", model_name, ignore.case = TRUE)) {
      roc.RF.ens <- survivalROC(Stime = xtest_time, status = xtest_status, marker = y.ens.RF, predict.time = time_point, method = "KM")
      auc.tab <- c(auc.tab, roc.RF.ens$AUC)
      sp.98.tab <- c(sp.98.tab, roc.RF.ens$TP[roc.RF.ens$FP <= 0.02][1])
      sp.95.tab <- c(sp.95.tab, roc.RF.ens$TP[roc.RF.ens$FP <= 0.05][1])
      model.names <- c(model.names, "RF")
    }
    if (grepl("EN", model_name, ignore.case = TRUE)) {
      roc.EN.ens <- survivalROC(Stime = xtest_time, status = xtest_status, marker = y.ens.EN, predict.time = time_point, method = "KM")
      auc.tab <- c(auc.tab, roc.EN.ens$AUC)
      sp.98.tab <- c(sp.98.tab, roc.EN.ens$TP[roc.EN.ens$FP <= 0.02][1])
      sp.95.tab <- c(sp.95.tab, roc.EN.ens$TP[roc.EN.ens$FP <= 0.05][1])
      model.names <- c(model.names, "EN")
    }
    if (grepl("MB", model_name, ignore.case = TRUE)) {
      roc.MB.ens <- survivalROC(Stime = xtest_time, status = xtest_status, marker = y.ens.MB, predict.time = time_point, method = "KM")
      auc.tab <- c(auc.tab, roc.MB.ens$AUC)
      sp.98.tab <- c(sp.98.tab, roc.MB.ens$TP[roc.MB.ens$FP <= 0.02][1])
      sp.95.tab <- c(sp.95.tab, roc.MB.ens$TP[roc.MB.ens$FP <= 0.05][1])
      model.names <- c(model.names, "MB")
    }
    if (grepl("ADB", model_name, ignore.case = TRUE)) {
      roc.ADB.ens <- survivalROC(Stime = xtest_time, status = xtest_status, marker = y.ens.ADB, predict.time = time_point, method = "KM")
      auc.tab <- c(auc.tab, roc.ADB.ens$AUC)
      sp.98.tab <- c(sp.98.tab, roc.ADB.ens$TP[roc.ADB.ens$FP <= 0.02][1])
      sp.95.tab <- c(sp.95.tab, roc.ADB.ens$TP[roc.ADB.ens$FP <= 0.05][1])
      model.names <- c(model.names, "ADB")
    }
    if (grepl("GB", model_name, ignore.case = TRUE)) {
      roc.GB.ens <- survivalROC(Stime = xtest_time, status = xtest_status, marker = y.ens.GB, predict.time = time_point, method = "KM")
      auc.tab <- c(auc.tab, roc.GB.ens$AUC)
      sp.98.tab <- c(sp.98.tab, roc.GB.ens$TP[roc.GB.ens$FP <= 0.02][1])
      sp.95.tab <- c(sp.95.tab, roc.GB.ens$TP[roc.GB.ens$FP <= 0.05][1])
      model.names <- c(model.names, "GB")
    }
    if (grepl("KNN", model_name, ignore.case = TRUE)) {
      roc.KNN.ens <- survivalROC(Stime = xtest_time, status = xtest_status, marker = y.ens.KNN, predict.time = time_point, method = "KM")
      auc.tab <- c(auc.tab, roc.KNN.ens$AUC)
      sp.98.tab <- c(sp.98.tab, roc.KNN.ens$TP[roc.KNN.ens$FP <= 0.02][1])
      sp.95.tab <- c(sp.95.tab, roc.KNN.ens$TP[roc.KNN.ens$FP <= 0.05][1])
      model.names <- c(model.names, "KNN")
    }
    if (grepl("NN", model_name, ignore.case = TRUE)) {
      roc.NN.ens <- survivalROC(Stime = xtest_time, status = xtest_status, marker = y.ens.NN, predict.time = time_point, method = "KM")
      auc.tab <- c(auc.tab, roc.NN.ens$AUC)
      sp.98.tab <- c(sp.98.tab, roc.NN.ens$TP[roc.NN.ens$FP <= 0.02][1])
      sp.95.tab <- c(sp.95.tab, roc.NN.ens$TP[roc.NN.ens$FP <= 0.05][1])
      model.names <- c(model.names, "NN")
    }
    if (grepl("SVM", model_name, ignore.case = TRUE)) {
      roc.SVM.ens <- survivalROC(Stime = xtest_time, status = xtest_status, marker = y.ens.SVM, predict.time = time_point, method = "KM")
      auc.tab <- c(auc.tab, roc.SVM.ens$AUC)
      sp.98.tab <- c(sp.98.tab, roc.SVM.ens$TP[roc.SVM.ens$FP <= 0.02][1])
      sp.95.tab <- c(sp.95.tab, roc.SVM.ens$TP[roc.SVM.ens$FP <= 0.05][1])
      model.names <- c(model.names, "SVM")
    }
    if (grepl("DCT", model_name, ignore.case = TRUE)) {
      roc.DCT.ens <- survivalROC(Stime = xtest_time, status = xtest_status, marker = y.ens.DCT, predict.time = time_point, method = "KM")
      auc.tab <- c(auc.tab, roc.DCT.ens$AUC)
      sp.98.tab <- c(sp.98.tab, roc.DCT.ens$TP[roc.DCT.ens$FP <= 0.02][1])
      sp.95.tab <- c(sp.95.tab, roc.DCT.ens$TP[roc.DCT.ens$FP <= 0.05][1])
      model.names <- c(model.names, "DCT")
    }
  }

  roc.all <- survivalROC(Stime = xtest_time, status = xtest_status, marker = y.pred.avg.all, predict.time = time_point, method = "KM")
  auc.tab <- c(auc.tab, roc.all$AUC)
  sp.98.tab <- c(sp.98.tab, roc.all$TP[roc.all$FP <= 0.02][1])
  sp.95.tab <- c(sp.95.tab, roc.all$TP[roc.all$FP <= 0.05][1])
  model.names <- c(model.names, "Weighted")

  roc.all1 <- survivalROC(Stime = xtest_time, status = xtest_status, marker = y.pred.avg1.all, predict.time = time_point, method = "KM")
  auc.tab <- c(auc.tab, roc.all1$AUC)
  sp.98.tab <- c(sp.98.tab, roc.all1$TP[roc.all1$FP <= 0.02][1])
  sp.95.tab <- c(sp.95.tab, roc.all1$TP[roc.all1$FP <= 0.05][1])
  model.names <- c(model.names, "Simple")

  get_marker <- function(model_name) {
    if (identical(model_name, "Weighted")) return(y.pred.avg.all)
    if (identical(model_name, "Simple")) return(y.pred.avg1.all)
    get(paste0("y.ens.", model_name), inherits = TRUE)
  }
  cindex_stats <- t(vapply(
    model.names,
    function(model_name) survival_cindex(xtest_time, xtest_status, get_marker(model_name)),
    numeric(3)
  ))
  sp95_ci <- t(vapply(
    model.names,
    function(model_name) survival_sp_ci(xtest_time, xtest_status, get_marker(model_name), time_point),
    numeric(2)
  ))
  auc_ci <- t(vapply(
    model.names,
    function(model_name) survival_auc_ci(xtest_time, xtest_status, get_marker(model_name), time_point),
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
