#' Evaluate model performance on the training data via cross-validation
#'
#' This function evaluates the predictive performance of trained models
#' on the training dataset after pathway or feature selection has been
#' completed. Performance metrics are computed for each specified base
#' model and, if enabled, the ensemble model. The results provide an
#' estimate of model fit and allow comparison between individual
#' algorithms used during the training stage.
#' @param outcome.type Character string specifying the prediction task. Must be one of `"survival"`, `"binary"`, or `"categorical"` corresponding to survival analysis, binary classification, and multi-class classification.
#' @param outcome Outcome vector for the training data. For `"binary"` and `"survival"` outcomes, this should be a binary vector (0/1). For `"categorical"` outcomes, this should be a numeric or factor vector representing class labels.
#' @param survdays Numeric vector of survival times required when
#' `outcome.type = "survival"`.
#' @param time_point Numeric value specifying the time point at which survival predictions are evaluated.
#' @param train_data Training data matrix or data frame.
#' @param features Feature specification. For pathway-based modeling this should be a list of pathways where each element is a character vector of genes belonging to that pathway. For protein-level modeling this should be a character vector of proteins used as features.
#' @param models.indiv Character vector specifying base machine learning or deep learning models used for individual model training. Supported abbreviations include: `"XG"` (xgboost), `"RF"` (random forest), `"EN"` (elastic net), `"GB"` (gradient boosting), `"KNN"` (k-nearest neighbors), `"SVM"` (support vector machine), `"NB"` (naive Bayes), `"DCT"` (decision tree), `"NN"` (neural network), `"ADB"` (AdaBoost), and `"MB"` (model-based boosting).
#' @param iter Integer specifying the number of iterations to run during model training.
#' @param num.folds Integer number of cross-validation folds. Values below 2
#' are treated as 2 internally. Default is 3.
#' @param ensemble Logical. If `TRUE`, ensemble (stacking) prediction modeling is performed.
#' @param models.ens Character vector specifying models used to construct the stacking/ensemble model. Model abbreviations follow the same convention as `models.indiv`.
#' @param param.indiv Named list of user-specified parameters for the individual base models listed in `models.indiv`. If `NULL`, default parameters appropriate for the specified outcome type are used.
#' @param param.ens Named list of user-specified parameters for ensemble model training. If `NULL`, default parameters are used.
#' @param parallel.iter Logical; if `TRUE`, repeated cross-validation
#' iterations are split across a temporary `future` backend before chunk-level
#' metrics are combined.
#' @param workers Integer number of future workers to use when
#' `parallel.iter = TRUE`. Defaults to `future::availableCores()` and is capped
#' at `iter`.
#' @param future.seed Passed to `future.apply::future_lapply()`.
#' @param future.strategy Future backend strategy used when `parallel.iter =
#' TRUE`. Default is `"multisession"`. Use `NULL` to use the caller's existing
#' future plan.
#'
#' @examples
#' \dontrun{
#' # -----------------------------------------------------------------------
#' # 1) Select pathways first
#' # -----------------------------------------------------------------------
#' data(caspen_example)
#'
#' path_res <- pathway_select(
#'   outcome.type = "binary",
#'   outcome = caspen_example$outcome,
#'   data = caspen_example$data_to_pass,
#'   features = caspen_example$pathways[1:10],
#'   models.indiv = c("RF", "EN", "SVM"),
#'   iter = 5,
#'   num.folds = 3,
#'   fs.AUC.cut = 0.70,
#'   fs.sp.cut = 0.10
#' )
#'
#' selected_features <- path_res$features
#'
#' # -----------------------------------------------------------------------
#' # 2) Estimate cross-validated training performance on selected features
#' # -----------------------------------------------------------------------
#' train_res <- Train_perform_CV(
#'   outcome.type = "binary",
#'   outcome = caspen_example$outcome,
#'   train_data = caspen_example$data_to_pass,
#'   features = selected_features,
#'   models.indiv = c("RF", "EN", "SVM"),
#'   iter = 10,
#'   ensemble = TRUE,
#'   models.ens = c("RF", "EN", "SVM"),
#'   parallel.iter = TRUE
#' )
#'
#' train_res$auc
#' train_res$SP.95
#'
#' # -----------------------------------------------------------------------
#' # Survival example
#' # -----------------------------------------------------------------------
#' set.seed(1)
#' survival_time <- sample(120:1800, nrow(caspen_example$data_to_pass),
#'                         replace = TRUE)
#' survival_status <- caspen_example$outcome
#'
#' surv_train <- Train_perform_CV(
#'   outcome.type = "survival",
#'   outcome = survival_status,
#'   survdays = survival_time,
#'   time_point = 365,
#'   train_data = caspen_example$data_to_pass,
#'   features = selected_features,
#'   models.indiv = c("RF", "EN", "DCT"),
#'   iter = 10,
#'   ensemble = TRUE,
#'   models.ens = c("RF", "EN")
#' )
#'
#' surv_train$C.index
#' surv_train$SP.95
#' }
#'
#' @export

Train_perform_CV <- function(outcome.type, outcome, survdays = NULL, time_point = NULL,
                             train_data, features, models.indiv, iter,
                             num.folds = NULL,
                             ensemble = FALSE, models.ens = NULL,
                             param.indiv = NULL, param.ens = NULL,
                             auto.tune = FALSE,
                             tune.method = c("random", "grid", "successive_halving",
                                             "hyperband", "bayes"),
                             tune.n = 5,
                             tune.folds = 3,
                             tune.iter = 1,
                             tune.models = NULL,
                             parallel.iter = FALSE,
                             workers = future::availableCores(),
                             future.seed = NULL,
                             future.strategy = "multisession") {
  type_out <- tolower(outcome.type)
  tune.method <- match.arg(tune.method)
  if (!type_out %in% c("binary", "categorical", "survival", "continuous")) {
    stop("Invalid outcome.type: Must be 'binary', 'categorical', 'survival', or 'continuous'.")
  }
  tuning.table <- NULL
  if (isTRUE(auto.tune)) {
    tune.genes <- unique(unlist(lapply(features, extract_pathway_genes)))
    tune.genes <- tune.genes[tune.genes %in% colnames(train_data)]
    if (length(tune.genes)) {
      tuned <- caspen_auto_tune_params(
        outcome.type = type_out,
        outcome = outcome,
        x = train_data[, tune.genes, drop = FALSE],
        models = models.indiv,
        param.indiv = param.indiv,
        tune.method = tune.method,
        tune.n = tune.n,
        tune.folds = tune.folds,
        tune.iter = tune.iter,
        tune.models = tune.models,
        survdays = survdays,
        time_point = time_point,
        seed = 1
      )
      param.indiv <- tuned$param.indiv
      tuning.table <- tuned$tuning.table
    }
  }
  if (type_out == "binary") {
    res <- run_train_performance_iterations(
      binary_train_performance,
      iter         = iter,
      parallel.iter = parallel.iter,
      workers      = workers,
      future.seed  = future.seed,
      future.strategy = future.strategy,
      outcome      = outcome,
      train_data   = train_data,
      features     = features,
      models.indiv = models.indiv,
      num.folds    = num.folds,
      ensemble     = ensemble,
      models.ens   = models.ens,
      param.indiv  = param.indiv,
      param.ens    = param.ens
    )
  } else if (type_out == "survival") {
    if (is.null(survdays) || is.null(time_point)) {
      stop("survdays and time_point are required for survival outcomes.")
    }
    res <- run_train_performance_iterations(
      surv_train_performance,
      iter         = iter,
      parallel.iter = parallel.iter,
      workers      = workers,
      future.seed  = future.seed,
      future.strategy = future.strategy,
      outcome      = outcome,
      survdays     = survdays,
      time_point   = time_point,
      train_data   = train_data,
      features     = features,
      models.indiv = models.indiv,
      num.folds    = num.folds,
      ensemble     = ensemble,
      models.ens   = models.ens,
      param.indiv  = param.indiv,
      param.ens    = param.ens
    )
  } else if (type_out == "categorical") {
    res <- run_train_performance_iterations(
      cat_train_performance,
      iter         = iter,
      parallel.iter = parallel.iter,
      workers      = workers,
      future.seed  = future.seed,
      future.strategy = future.strategy,
      outcome      = outcome,
      train_data   = train_data,
      features     = features,
      models.indiv = models.indiv,
      num.folds    = num.folds,
      ensemble     = ensemble,
      models.ens   = models.ens,
      param.indiv  = param.indiv,
      param.ens    = param.ens
    )
  } else if (type_out == "continuous") {
    res <- run_train_performance_iterations(
      continuous_train_performance,
      iter         = iter,
      parallel.iter = parallel.iter,
      workers      = workers,
      future.seed  = future.seed,
      future.strategy = future.strategy,
      outcome      = outcome,
      train_data   = train_data,
      features     = features,
      models.indiv = models.indiv,
      num.folds    = num.folds,
      ensemble     = ensemble,
      models.ens   = models.ens,
      param.indiv  = param.indiv,
      param.ens    = param.ens
    )
  }

  if (isTRUE(auto.tune)) attr(res, "tuning") <- tuning.table
  return(res)
}
