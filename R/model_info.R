#' CASPEN Model Backend Information
#'
#' Returns the model abbreviations used by CASPEN, the main R package used for
#' each backend, and the installed package version on the user's machine. This
#' is useful for reproducibility sections in reports, vignettes, and manuscript
#' supplements.
#'
#' @return A data frame with model abbreviation, model family, backend package,
#' installed package version, and notes.
#'
#' @examples
#' caspen_model_info()
#'
#' @export
caspen_model_info <- function() {
  info <- data.frame(
    model = c("XG", "RF", "EN", "GB", "KNN", "SVM", "NB", "DCT",
              "NN", "ADB", "MB"),
    family = c("gradient boosted trees", "random forest",
               "elastic net / Cox elastic net", "gradient boosting",
               "k-nearest neighbors", "support vector machine",
               "naive Bayes", "decision tree", "neural network",
               "AdaBoost", "model-based boosting"),
    package = c("xgboost", "randomForest / randomForestSRC", "glmnet",
                "gbm", "caret", "e1071", "naivebayes / e1071",
                "rpart", "neuralnet / nnet", "JOUSBoost", "mboost"),
    notes = c(
      "Used for binary, categorical, and Cox-style survival risk models.",
      "randomForest for classification; randomForestSRC for survival.",
      "cv.glmnet/glmnet for binary, multiclass, and Cox models.",
      "gbm/gmb.fit for classification and Cox proportional hazards boosting.",
      "caret KNN helpers for classification; CASPEN KNN risk helper for survival.",
      "e1071 SVM for classification; survivalsvm for survival.",
      "naivebayes/e1071 implementations depending on task.",
      "rpart decision trees for classification and survival-style trees.",
      "neuralnet/nnet depending on outcome type and helper function.",
      "JOUSBoost AdaBoost helper for classification and time-point survival risk.",
      "mboost glmboost for binary, multiclass one-vs-rest, and Cox boosting."
    ),
    stringsAsFactors = FALSE
  )

  info$installed_version <- vapply(info$package, function(pkg_string) {
    pkgs <- unlist(strsplit(pkg_string, "\\s*/\\s*"))
    versions <- vapply(pkgs, function(pkg) {
      if (requireNamespace(pkg, quietly = TRUE)) {
        as.character(utils::packageVersion(pkg))
      } else {
        NA_character_
      }
    }, character(1))
    paste(ifelse(is.na(versions), "not installed", versions), collapse = " / ")
  }, character(1))

  info[, c("model", "family", "package", "installed_version", "notes")]
}
