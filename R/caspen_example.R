#' Example dataset for CASPEN
#'
#' An example dataset used to demonstrate the workflow of CASPEN
#' binary pathway-based and protein-based prediction modeling.
#'
#' @format A list containing:
#' \describe{
#'
#'   \item{data_to_pass}{Protein expression matrix used for pathway-based
#'   modeling (samples x genes).}
#'
#'   \item{outcome}{Binary outcome vector for classification tasks.}
#'
#'   \item{pathways}{List of pathways where each element is a character
#'   vector of gene symbols.}
#' }
#'
#' @source Derived from internal study datasets
"caspen_example"
