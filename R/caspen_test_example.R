#' Independent Test Dataset For CASPEN Examples
#'
#' An independent binary-outcome protein test dataset used to demonstrate
#' `Test_perform()` after training or pathway selection on `caspen_example`.
#' The original CSV files are also installed under `inst/extdata`.
#'
#' @format A list containing:
#' \describe{
#'
#'   \item{test_data}{Protein expression matrix for independent testing
#'   (samples x genes/features).}
#'
#'   \item{test_outcome}{Binary outcome vector for the independent test
#'   samples.}
#'
#'   \item{sample_id}{Sample identifiers for `test_data` and
#'   `test_outcome`.}
#'
#'   \item{raw_data_file}{Name of the raw packaged test data CSV.}
#'
#'   \item{raw_outcome_file}{Name of the raw packaged test outcome CSV.}
#' }
#'
#' @source Derived from internal study datasets
"caspen_test_example"
