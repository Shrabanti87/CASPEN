#' Evaluate gene-level model performance on training data by cross-validation
#'
#' `Gene_Train_perform_CV()` fits CASPEN models directly on a gene/protein
#' feature matrix instead of looping over pathways. In each fold, base learners
#' are fit on the selected genes and their predictions become the inputs to the
#' stacked ensemble learners. Predictions are aggregated over cross-validation
#' folds and repeated iterations. With `ensemble = TRUE`, returned rows are the
#' requested ensemble learners plus simple and weighted average ensembles; with
#' `ensemble = FALSE`, returned rows are the individual base learners.
#'
#' @param outcome.type Character string specifying the prediction task. Must be
#' one of `"binary"`, `"categorical"`, `"survival"`, or `"continuous"`.
#' @param outcome Outcome vector for the training data.
#' @param survdays Numeric vector of survival times, required when
#' `outcome.type = "survival"`.
#' @param time_point Numeric value specifying the time point at which survival
#' time-dependent AUC and SP95 are evaluated.
#' @param train_data Training data matrix or data frame with samples in rows and
#' genes/proteins/features in columns.
#' @param genes Optional character vector of genes/proteins/features to use.
#' If `NULL`, all columns of `train_data` are used.
#' @param models.indiv Character vector specifying base models. If `NULL`, all
#' models supported by the selected outcome type are used.
#' @param iter Integer number of repeated cross-validation iterations.
#' @param num.folds Integer number of cross-validation folds. Default is 3.
#' @param ensemble Logical. If `TRUE`, fit stacked ensemble models. Defaults to
#' `TRUE` for gene-level performance because there is no pathway block and the
#' aggregation is over base learners.
#' @param models.ens Character vector specifying ensemble learners.
#' @param param.indiv Named list of model parameters for base models.
#' @param param.ens Named list of model parameters for ensemble models.
#' @param auto.tune Logical. If `TRUE`, tune base-model parameters on the gene
#' matrix before performance evaluation.
#' @param tune.method,tune.n,tune.folds,tune.iter,tune.models Tuning controls.
#' @param parallel.iter Logical; if `TRUE`, repeated iterations are split across
#' a temporary `future` backend.
#' @param workers Integer number of future workers. Defaults to available cores.
#' @param future.seed Passed to `future.apply::future_lapply()`.
#' @param future.strategy Future backend strategy used when `parallel.iter =
#' TRUE`.
#'
#' @return A data frame of model-level performance metrics. For binary and
#' categorical outcomes, AUC and SP95 are returned. For survival outcomes,
#' C-index, time-dependent AUC, and SP95 are returned. For continuous outcomes,
#' C-index, R2, MAE, and RMSE are returned when available.
#'
#' @examples
#' \dontrun{
#' data(caspen_example)
#'
#' gene_train <- Gene_Train_perform_CV(
#'   outcome.type = "binary",
#'   outcome = caspen_example$outcome,
#'   train_data = caspen_example$data_to_pass,
#'   genes = colnames(caspen_example$data_to_pass)[1:40],
#'   models.indiv = c("RF", "EN", "XG"),
#'   iter = 5,
#'   num.folds = 3,
#'   ensemble = TRUE,
#'   models.ens = c("RF", "EN", "XG")
#' )
#'
#' gene_train$auc
#' gene_train$SP.95
#' }
#'
#' @export
Gene_Train_perform_CV <- function(outcome.type, outcome, survdays = NULL,
                                  time_point = NULL, train_data, genes = NULL,
                                  models.indiv = NULL, iter, num.folds = NULL,
                                  ensemble = TRUE, models.ens = NULL,
                                  param.indiv = NULL, param.ens = NULL,
                                  auto.tune = FALSE,
                                  tune.method = c("random", "grid",
                                                  "successive_halving",
                                                  "hyperband", "bayes"),
                                  tune.n = 5,
                                  tune.folds = 3,
                                  tune.iter = 1,
                                  tune.models = NULL,
                                  parallel.iter = FALSE,
                                  workers = future::availableCores(),
                                  future.seed = NULL,
                                  future.strategy = "multisession") {
  models.indiv <- caspen_default_models_for_outcome(outcome.type, models.indiv)
  if (isTRUE(ensemble) && is.null(models.ens)) models.ens <- models.indiv
  features <- caspen_gene_level_features(train_data = train_data, genes = genes)
  res <- Train_perform_CV(
    outcome.type = outcome.type,
    outcome = outcome,
    survdays = survdays,
    time_point = time_point,
    train_data = train_data,
    features = features,
    models.indiv = models.indiv,
    iter = iter,
    num.folds = num.folds,
    ensemble = ensemble,
    models.ens = models.ens,
    param.indiv = param.indiv,
    param.ens = param.ens,
    auto.tune = auto.tune,
    tune.method = tune.method,
    tune.n = tune.n,
    tune.folds = tune.folds,
    tune.iter = tune.iter,
    tune.models = tune.models,
    parallel.iter = parallel.iter,
    workers = workers,
    future.seed = future.seed,
    future.strategy = future.strategy
  )
  attr(res, "analysis.level") <- "gene"
  attr(res, "genes") <- features[[1]]
  res
}

#' Evaluate gene-level model performance on independent testing data
#'
#' `Gene_Test_perform()` fits CASPEN models directly on a gene/protein feature
#' matrix and evaluates them on an independent test set. There is no pathway
#' loop. Base learners are fit on the selected genes, and their fitted
#' prediction columns are used by the stacked ensemble learners. With `ensemble
#' = TRUE`, returned rows are the requested ensemble learners plus simple and
#' weighted average ensembles; with `ensemble = FALSE`, returned rows are the
#' individual base learners.
#'
#' @param outcome.type Character string specifying the prediction task. Must be
#' one of `"binary"`, `"categorical"`, `"survival"`, or `"continuous"`.
#' @param outcome Outcome vector for the training data.
#' @param survdays Numeric vector of training survival times, required when
#' `outcome.type = "survival"`.
#' @param time_point Numeric value specifying the time point at which survival
#' time-dependent AUC and SP95 are evaluated.
#' @param train_data Training data matrix or data frame.
#' @param test_data Independent test data matrix or data frame.
#' @param test_outcome Outcome vector for the test data.
#' @param test_time Survival times for the test data, required when
#' `outcome.type = "survival"`.
#' @param genes Optional character vector of genes/proteins/features to use.
#' If `NULL`, all shared columns between `train_data` and `test_data` are used.
#' @param models.indiv Character vector specifying base models. If `NULL`, all
#' models supported by the selected outcome type are used.
#' @param ensemble Logical. If `TRUE`, fit stacked ensemble models. Defaults to
#' `TRUE` for gene-level performance because there is no pathway block and the
#' aggregation is over base learners.
#' @param models.ens Character vector specifying ensemble learners.
#' @param param.indiv Named list of model parameters for base models.
#' @param param.ens Named list of model parameters for ensemble models.
#' @param auto.tune Logical. If `TRUE`, tune base-model parameters on the
#' training gene matrix before independent test evaluation.
#' @param tune.method,tune.n,tune.folds,tune.iter,tune.models Tuning controls.
#'
#' @return A data frame of model-level independent test performance metrics.
#'
#' @examples
#' \dontrun{
#' data(caspen_example)
#' data(caspen_test_example)
#'
#' genes <- colnames(caspen_example$data_to_pass)[1:40]
#' gene_test <- Gene_Test_perform(
#'   outcome.type = "binary",
#'   outcome = caspen_example$outcome,
#'   train_data = caspen_example$data_to_pass,
#'   test_data = caspen_test_example$test_data,
#'   test_outcome = caspen_test_example$test_outcome,
#'   genes = genes,
#'   models.indiv = c("RF", "EN", "XG"),
#'   ensemble = TRUE,
#'   models.ens = c("RF", "EN", "XG")
#' )
#'
#' gene_test$auc
#' gene_test$SP.95
#' }
#'
#' @export
Gene_Test_perform <- function(outcome.type, outcome, survdays = NULL,
                              time_point = NULL, train_data, test_data,
                              test_outcome, test_time = NULL, genes = NULL,
                              models.indiv = NULL, ensemble = TRUE,
                              models.ens = NULL, param.indiv = NULL,
                              param.ens = NULL, auto.tune = FALSE,
                              tune.method = c("random", "grid",
                                              "successive_halving",
                                              "hyperband", "bayes"),
                              tune.n = 5,
                              tune.folds = 3,
                              tune.iter = 1,
                              tune.models = NULL) {
  models.indiv <- caspen_default_models_for_outcome(outcome.type, models.indiv)
  if (isTRUE(ensemble) && is.null(models.ens)) models.ens <- models.indiv
  features <- caspen_gene_level_features(
    train_data = train_data,
    test_data = test_data,
    genes = genes
  )
  res <- Test_perform(
    outcome.type = outcome.type,
    outcome = outcome,
    survdays = survdays,
    time_point = time_point,
    train_data = train_data,
    test_data = test_data,
    test_outcome = test_outcome,
    test_time = test_time,
    features = features,
    models.indiv = models.indiv,
    ensemble = ensemble,
    models.ens = models.ens,
    param.indiv = param.indiv,
    param.ens = param.ens,
    auto.tune = auto.tune,
    tune.method = tune.method,
    tune.n = tune.n,
    tune.folds = tune.folds,
    tune.iter = tune.iter,
    tune.models = tune.models
  )
  attr(res, "analysis.level") <- "gene"
  attr(res, "genes") <- features[[1]]
  res
}

#' @keywords internal
#' @noRd
caspen_gene_level_features <- function(train_data, test_data = NULL,
                                       genes = NULL,
                                       feature.name = "All_genes") {
  train_cols <- colnames(train_data)
  if (is.null(train_cols)) stop("train_data must have column names.")
  if (is.null(genes)) {
    genes <- train_cols
  } else {
    genes <- unique(as.character(genes))
  }

  if (!is.null(test_data)) {
    test_cols <- colnames(test_data)
    if (is.null(test_cols)) stop("test_data must have column names.")
    genes <- intersect(genes, intersect(train_cols, test_cols))
  } else {
    genes <- intersect(genes, train_cols)
  }

  if (!length(genes)) {
    stop("No requested genes/features were found in the supplied data.")
  }
  stats::setNames(list(genes), feature.name)
}

#' @keywords internal
#' @noRd
caspen_default_models_for_outcome <- function(outcome.type, models = NULL) {
  if (!is.null(models)) return(toupper(models))
  type_out <- tolower(outcome.type)
  switch(
    type_out,
    binary = c("XG", "RF", "EN", "GB", "KNN", "SVM", "NB", "DCT",
               "NN", "ADB", "MB"),
    categorical = c("XG", "RF", "EN", "GB", "KNN", "SVM", "NB", "DCT",
                    "NN", "ADB", "MB"),
    survival = c("XG", "RF", "EN", "GB", "KNN", "SVM", "DCT",
                 "NN", "ADB", "MB"),
    continuous = c("XG", "RF", "EN", "GB", "KNN", "SVM", "DCT",
                   "NN", "ADB", "MB"),
    stop("Invalid outcome.type: Must be 'binary', 'categorical', 'survival', or 'continuous'.")
  )
}
