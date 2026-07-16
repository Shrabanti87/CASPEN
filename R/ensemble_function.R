#' Train ML Models with Optional Ensemble
#'
#' This function trains models based on survival analysis, binary, or multi-class classification. You can train individual ML/DL models or an ensemble of models. The function supports specifying hyper-parameters for individual models or the ensemble, and feature selection AUC, specificity cutoff.
#'
#' @param outcome.type Character string specifying the prediction task. Must be one of `"survival"`, `"binary"`, or `"categorical"` corresponding to survival analysis, binary classification, and multi-class classification.
#' @param outcome Outcome vector for the training data. For `"binary"` and `"survival"` outcomes, this should be a binary vector (0/1). For `"categorical"` outcomes, this should be a numeric or factor vector representing class labels.
#' @param data Numeric matrix of dimension n x p, where n is the number of samples and p is the number of genes or proteins.
#' @param iter Integer specifying the number of iterations to run during model training.
#' @param num.folds Integer specifying the number of folds to perform cross validation on. Default is 3.
#' @param features Feature specification. For pathway-based modeling this should be a list of pathways where each element is a character vector of genes belonging to that pathway. For protein-level modeling this should be a character vector of proteins used as features.
#' @param ensemble Logical indicating whether ensemble (stacking) prediction modeling should be performed.
#' @param models.indiv Character vector specifying base machine learning or deep learning models used for individual model training. Supported abbreviations include: `"XG"` (xgboost), `"RF"` (random forest), `"EN"` (elastic net), `"GB"` (gradient boosting), `"KNN"` (k-nearest neighbors), `"SVM"` (support vector machine), `"NB"` (naive Bayes), `"DCT"` (decision tree), `"NN"` (neural network), `"ADB"` (AdaBoost), and `"MB"` (model-based boosting).
#' @param models.ens Character vector specifying models used to construct the stacking/ensemble model. Model abbreviations follow the same convention as `models.indiv`.
#' @param param.indiv Named list of user-specified parameters for the individual base models listed in `models.indiv`. If `NULL`, default parameters appropriate for the specified outcome type are used.
#' @param param.ens Named list of user-specified parameters for ensemble model training. If `NULL`, default parameters are used.
#' @param fs.AUC.cut Numeric value specifying the AUC threshold used for pathway selection during feature selection.
#' For survival outcomes, this threshold is applied to the C-index when
#' available.
#' @param fs.sp.cut Numeric value specifying the sensitivity at 95% specificity threshold used for pathway selection.
#' @param num.sel.path Integer specifying the maximum number of pathways to select. If provided, the final number of selected pathways is the minimum of `num.sel.path` and the number selected using `fs.AUC.cut` and `fs.sp.cut`.
#' @param literature.features Optional list of literature-curated pathways
#' or gene sets to pool with `features`, for example from an external LLM
#' curation workflow. Each element can be a character vector of genes, or
#' a named list where names are genes and values store citation metadata.
#' @param pathway.prior Optional prior weights for pathways, supplied as a
#' named numeric vector/list or a data frame with pathway/name and
#' prior/weight columns. Weights are aligned to the pooled pathway names.
#' @param prior.weight Deprecated; retained for backward compatibility. Prior
#' weights are reported separately and are not blended into the data metric.
#' @param prior.mode Deprecated; retained for backward compatibility. Prior
#' weights are reported separately and are not blended into the data metric.
#' @param prior.cut Optional minimum prior weight required for pathway
#' selection.
#' @param prior.default Default prior weight assigned to pathways missing
#' from `pathway.prior`. Defaults to `NA_real_` so missing priors are visible
#' in the returned `Prior.weight` matrix.
#' @param prior.clip Numeric length-two range used to clip supplied prior
#' weights.
#' @param redundancy.thresh Optional Jaccard similarity threshold used to
#' skip literature gene sets that are redundant with already supplied
#' pathways.
#' @param llm.mode Character string controlling live LLM/literature use. One of
#' `"none"`, `"curate"`, or `"prior"`. `"curate"` asks the LLM for
#' disease/outcome-specific biological concepts, scores priors, maps those
#' concepts to the user-provided pathway universe, and models only matched
#' pathways. `"prior"` scores existing pathway names without curating concepts.
#' No LLM call is made when `"none"`.
#' @param llm.curate.action How curated LLM concepts are used. `"match_pathways"`
#' maps concepts to the supplied pathway list and runs only matched pathways;
#' `"add_signatures"` adds LLM-curated gene sets as new signatures.
#' @param disease Disease context used for LLM curation, required when
#' `llm.mode != "none"`.
#' @param outcome.context Outcome context used for LLM curation, required when
#' `llm.mode != "none"`.
#' @param llm.provider LLM provider passed to `llm_literature_signatures()`.
#' Currently `"openai"` or `"custom"`.
#' @param llm.model OpenAI model name. Defaults to `CASPEN_OPENAI_MODEL` or
#' `"gpt-4.1-mini"`.
#' @param llm.api.key OpenAI API key. Defaults to `OPENAI_API_KEY`.
#' @param llm.fn Optional custom LLM function accepting a prompt and returning
#' text. If supplied, no OpenAI request is made.
#' @param llm.n.signatures Maximum number of new signatures to request.
#' @param llm.max.genes Maximum number of genes per LLM-curated signature.
#' @param llm.min.genes Minimum number of curated genes that must map to the
#' input data for a new LLM-curated signature to be retained. Default is 5.
#' @param llm.pubmed Logical; if `TRUE`, retrieves PubMed abstracts as
#' grounding context.
#' @param llm.pubmed.max Maximum number of PubMed records to retrieve.
#' @param llm.email Optional email passed to NCBI E-utilities.
#' @param llm.valid.genes Optional valid gene-symbol universe for filtering LLM
#' gene symbols. If `NULL`, CASPEN infers this from `features`.
#' @param llm.prior.method Character string specifying how LLM/literature prior
#' weights are computed. `"llm"` uses the LLM-assigned score,
#' `"pubmed_count"` uses normalized PubMed paper counts for disease + outcome +
#' pathway, `"pubtator"` uses PubTator3 gene-annotation evidence, and
#' `"hybrid"` combines component scores.
#' @param llm.hybrid.weights Named numeric vector used when
#' `llm.prior.method = "hybrid"`. Supported names are `"pubmed_count"`,
#' `"pubtator"`, `"llm"`, `"celltype"`, and `"evidence"`.
#' @param llm.temperature OpenAI sampling temperature.
#' @param celltypes Optional character vector of cell type names present in
#' celltype-specific expression data. When supplied, bulk pathway gene lists
#' are expanded into pathway-celltype features using expression columns that
#' encode both gene and cell type.
#' @param celltype.sep Character vector of separators used to match
#' celltype-specific expression columns, for example `"GENE_Celltype"` or
#' `"Celltype_GENE"`.
#' @param celltype.position Character string specifying whether cell type
#' appears as a column-name suffix, prefix, or should be inferred from both.
#' One of `"auto"`, `"suffix"`, or `"prefix"`.
#' @param celltype.aliases Optional named list of aliases used when matching
#' pathway names to cell types. Names should be entries of `celltypes`; values
#' are character vectors of alternative labels. For example,
#' `list(CD4_Tcells = c("immune", "T cells", "CD4"),
#' CD8_Tcells = c("immune", "T cells", "CD8"),
#' Macrophages = c("immune", "macrophage"))`.
#' @param parallel.iter Logical; if `TRUE`, repeated iterations for each
#' pathway are split across a temporary `future` multisession backend.
#' @param parallel.pathways Logical; if `TRUE`, pathway or pathway-celltype
#' tasks are split across a temporary `future` backend. This is usually the
#' preferred parallelization mode when many pathways/celltypes are evaluated.
#' If both `parallel.pathways` and `parallel.iter` are `TRUE`, CASPEN
#' parallelizes pathways and runs iterations serially within each pathway to
#' avoid nested worker oversubscription.
#' @param workers Integer number of future workers to use when
#' `parallel.iter = TRUE` or `parallel.pathways = TRUE`. Defaults to
#' `future::availableCores()`. The actual worker count is capped at the number
#' of iterations or pathway tasks being split.
#' @param future.seed Passed to `future.apply::future_lapply()`. Default is
#' `NULL` because CASPEN already calls `set.seed()` by iteration number inside
#' each worker, preserving serial-compatible fold generation without additional
#' future RNG streams.
#' @param future.strategy Future backend strategy used when `parallel.iter =
#' TRUE`. Default is `"multisession"`. Use `"multicore"` on platforms where
#' forked workers are preferred, or `NULL` to use the caller's existing future
#' plan.
#'
#' @param survdays Numeric vector of survival times required when
#' `outcome.type = "survival"`.
#' @param time_point Numeric value specifying the time point at which survival predictions are evaluated.
#' @return A list containing separate matrices for the primary data metric
#' (`data.metric`; AUC for binary/categorical outcomes and C-index for survival),
#' sensitivity at 95% specificity (`SP.95`), prior weights (`Prior.weight`),
#' selected pathway labels, and a selection summary table.
#'
#' @details
#' Default hyperparameters are provided for each supported model depending
#' on the prediction task (`"binary"`, `"categorical"`, or `"survival"`).
#' Users may override these defaults by supplying custom parameter lists
#' through `param.indiv` or `param.ens`. If these arguments are set to
#' `NULL`, the defaults described below are used.
#'
#' **Multi-class classification defaults**
#'
#' - `"XG"` (xgboost): `max_depth = 2`, `eta = 0.1`, `nrounds = 50`
#' - `"RF"` (random forest): `ntree = 500`
#' - `"EN"` (elastic net): `nfolds = 3`, `s = "lambda.1se"`
#' - `"ADB"` (one-vs-rest AdaBoost): `tree_depth = 3`, `n_rounds = 200`
#' - `"MB"` (one-vs-rest model-based boosting): `mstop = 300`
#' - `"GB"` (gradient boosting): `n.trees = 50`, `interaction.depth = 6`,
#'   `shrinkage = 0.1`
#' - `"KNN"` (k-nearest neighbors): cross-validation with
#'   `trainControl(method = "cv", number = 3, classProbs = TRUE)`
#' - `"SVM"` (support vector machine): `kernel = "radial"`, `cost = 50`,
#'   `gamma = 0.75`
#' - `"NB"` (naive Bayes): `laplace = 10`, `threshold = 0.001`
#' - `"DCT"` (decision tree): `cp = 0.001`
#' - `"NN"` (neural network): `hidden = c(8,4)`, `rep = 1`
#'
#' **Binary classification defaults**
#'
#' - `"XG"` (xgboost): `max_depth = 2`, `eta = 0.1`, `nthread = 2`,
#'   `nrounds = 500`, `objective = "binary:logistic"`
#' - `"RF"` (random forest): `ntree = 500`
#' - `"EN"` (elastic net): `alpha = 0.2`, `family = "binomial"`,
#'   `s = "lambda.1se"`
#' - `"ADB"` (AdaBoost): `tree_depth = 3`, `n_rounds = 200`
#' - `"MB"` (model-based boosting): `mstop = 300`
#' - `"GB"` (gradient boosting): `cv.folds = 10`, `shrinkage = 0.01`,
#'   `n.minobsinnode = 2`, `n.trees = 500`, `distribution = "bernoulli"`
#' - `"KNN"` (k-nearest neighbors): `k = 5`
#' - `"NN"` (neural network): `hidden = c(8,4)`, `rep = 1`
#' - `"NB"` (naive Bayes): `usekernel = TRUE`
#' - `"DCT"` (decision tree): `method = "class"`
#'
#' **Survival analysis defaults**
#'
#' - `"XG"` (survival xgboost): `eta = 0.1`, `nrounds = 1000`
#' - `"RF"` (random survival forest): `ntree = 2000`, `max.depth = 10`,
#'   `importance = TRUE`
#' - `"EN"` (elastic net Cox model): `alpha = 0.2`, `s = "lambda.min"`
#' - `"MB"` (model-based boosting): `mstop = 300`
#' - `"ADB"` (time-point AdaBoost survival risk): `tree_depth = 3`,
#'   `n_rounds = 200`
#' - `"GB"` (gradient boosting): `n.trees = 1000`, `interaction.depth = 10`,
#'   `shrinkage = 0.005`, `bag.fraction = 0.5`
#' - `"SVM"` (survival SVM): `gamma.mu = 0.01`, `opt.meth = "quadprog"`,
#'   `kernel = "add_kernel"`
#' - `"DCT"` (decision tree): `minsplit = 10`, `maxdepth = 10`,
#'   `cp = 0.001`
#'
#' Ensemble model defaults follow the same structure but may use reduced
#' parameter settings to avoid overfitting when stacking predictions.
#'
#' @examples
#' # -------------------------------------------------------------------------
#' # Minimal binary pathway-selection example
#' # -------------------------------------------------------------------------
#' data(caspen_example)
#' data_to_pass <- caspen_example$data_to_pass
#' output <- caspen_example$outcome
#' dq <- caspen_example$pathways
#'
#' AUC.mat <- pathway_select(
#'   outcome.type = "binary",
#'   data = data_to_pass,
#'   outcome = output,
#'   iter = 2,
#'   features = dq[1:5],
#'   models.indiv = c("XG","RF","EN","KNN","GB"),
#'   ensemble = FALSE,
#'   num.sel.path = 3,
#'   fs.AUC.cut = 0.75,
#'   fs.sp.cut = 0.05
#' )
#'
#' # Primary outputs are returned as separate matrices/tables:
#' AUC.mat$data.metric
#' AUC.mat$SP.95
#' AUC.mat$Prior.weight
#' AUC.mat$selection.table
#'
#' \dontrun{
#' # -----------------------------------------------------------------------
#' # Celltype-aware pathway selection
#' # -----------------------------------------------------------------------
#' # If data columns contain celltype-specific features such as
#' # "TP53_CD8_Tcells" or "CD8_Tcells_TP53", CASPEN maps each bulk pathway
#' # gene list to the matching celltype-specific columns. Pathway names that
#' # contain a celltype alias are evaluated only for matching celltypes.
#' celltypes <- c("Tumor", "CD4_Tcells", "CD8_Tcells", "Macrophages")
#' celltype.aliases <- list(
#'   CD4_Tcells = c("CD4", "T cells", "Tcell", "immune"),
#'   CD8_Tcells = c("CD8", "T cells", "Tcell", "immune"),
#'   Macrophages = c("macrophage", "macrophages", "myeloid", "immune"),
#'   Tumor = c("tumor", "cancer", "epithelial")
#' )
#'
#' ct_res <- pathway_select(
#'   outcome.type = "binary",
#'   outcome = output,
#'   data = data_to_pass,
#'   features = dq[1:5],
#'   models.indiv = c("RF", "EN", "SVM"),
#'   iter = 10,
#'   num.folds = 3,
#'   celltypes = celltypes,
#'   celltype.aliases = celltype.aliases,
#'   parallel.iter = TRUE
#' )
#'
#' # -----------------------------------------------------------------------
#' # Add user-supplied prior weights
#' # -----------------------------------------------------------------------
#' manual_prior <- c(
#'   AKT = 0.80,
#'   DNA_repair = 0.90,
#'   Cell_cycle = 0.70
#' )
#'
#' prior_res <- pathway_select(
#'   outcome.type = "binary",
#'   outcome = output,
#'   data = data_to_pass,
#'   features = dq[1:5],
#'   models.indiv = c("RF", "EN"),
#'   iter = 5,
#'   num.folds = 3,
#'   pathway.prior = manual_prior,
#'   prior.cut = 0.6
#' )
#'
#' prior_res$data.metric      # AUC for binary/categorical, C-index for survival
#' prior_res$SP.95            # sensitivity at 95 percent specificity
#' prior_res$Prior.weight     # literature/user prior weights
#' prior_res$selection.table  # one row per pathway-celltype/model summary
#'
#' # -----------------------------------------------------------------------
#' # Ask CASPEN to score pathway priors from literature
#' # -----------------------------------------------------------------------
#' # No LLM key is needed for PubMed-count priors. CASPEN searches PubMed for
#' # disease + outcome + pathway and normalizes log paper counts to 0-1.
#' pubmed_prior_res <- pathway_select(
#'   outcome.type = "binary",
#'   outcome = output,
#'   data = data_to_pass,
#'   features = dq[1:5],
#'   models.indiv = c("RF", "EN"),
#'   iter = 5,
#'   num.folds = 3,
#'   llm.mode = "prior",
#'   llm.prior.method = "pubmed_count",
#'   disease = "high-grade serous ovarian cancer",
#'   outcome.context = "platinum resistance",
#'   llm.email = "your.email@example.com"
#' )
#'
#' pubmed_prior_res$Prior.weight
#' pubmed_prior_res$llm$evidence[, c("pathway", "pubmed.count",
#'                                   "pubmed.score", "prior")]
#'
#' # PubTator3 priors use PubMed IDs plus PubTator gene annotations.
#' pubtator_prior_res <- pathway_select(
#'   outcome.type = "binary",
#'   outcome = output,
#'   data = data_to_pass,
#'   features = dq[1:5],
#'   models.indiv = c("RF", "EN"),
#'   iter = 5,
#'   num.folds = 3,
#'   llm.mode = "prior",
#'   llm.prior.method = "pubtator",
#'   disease = "high-grade serous ovarian cancer",
#'   outcome.context = "platinum resistance",
#'   llm.pubmed.max = 20,
#'   llm.email = "your.email@example.com"
#' )
#'
#' pubtator_prior_res$llm$evidence[, c("pathway", "pubtator.count",
#'                                     "pubtator.score", "prior")]
#'
#' # OpenAI or a local LLM can also be used for LLM or hybrid priors.
#' Sys.setenv(OPENAI_API_KEY = "sk-proj-your_real_key_here")
#' hybrid_prior_res <- pathway_select(
#'   outcome.type = "binary",
#'   outcome = output,
#'   data = data_to_pass,
#'   features = dq[1:5],
#'   models.indiv = c("RF", "EN"),
#'   iter = 5,
#'   num.folds = 3,
#'   llm.mode = "prior",
#'   llm.prior.method = "hybrid",
#'   llm.hybrid.weights = c(pubmed_count = 0.25, pubtator = 0.25,
#'                          llm = 0.4, evidence = 0.1),
#'   disease = "high-grade serous ovarian cancer",
#'   outcome.context = "platinum resistance",
#'   llm.pubmed = FALSE,
#'   llm.model = "gpt-4o-mini",
#'   llm.temperature = NULL
#' )
#'
#' hybrid_prior_res$llm$evidence[, c("pathway", "pubmed.score",
#'                                   "pubtator.score", "llm.score",
#'                                   "evidence.score", "prior")]
#'
#' # -----------------------------------------------------------------------
#' # Survival and categorical outcomes
#' # -----------------------------------------------------------------------
#' # For survival, outcome is the event/status vector and survdays is the
#' # survival time vector. The primary data metric is C-index.
#' set.seed(1)
#' survival_time <- sample(120:1800, nrow(data_to_pass), replace = TRUE)
#' survival_status <- output
#'
#' surv_res <- pathway_select(
#'   outcome.type = "survival",
#'   outcome = survival_status,
#'   survdays = survival_time,
#'   data = data_to_pass,
#'   features = dq[1:5],
#'   models.indiv = c("RF", "EN", "DCT"),
#'   iter = 5,
#'   num.folds = 3
#' )
#' surv_res$data.metric
#'
#' # For categorical outcomes, outcome should be a factor or multiclass label.
#' category <- cut(
#'   rowMeans(data_to_pass[, 1:10]),
#'   breaks = stats::quantile(rowMeans(data_to_pass[, 1:10]),
#'                            probs = c(0, 1/3, 2/3, 1)),
#'   include.lowest = TRUE,
#'   labels = c("low", "medium", "high")
#' )
#'
#' cat_res <- pathway_select(
#'   outcome.type = "categorical",
#'   outcome = category,
#'   data = data_to_pass,
#'   features = dq[1:5],
#'   models.indiv = c("RF", "EN", "DCT"),
#'   iter = 5,
#'   num.folds = 3
#' )
#' cat_res$data.metric
#' }
#'
#' @export
pathway_select <- function(outcome.type, outcome, features, data, models.indiv, iter,
                           num.folds = NULL, ensemble = FALSE,
                           models.ens = NULL, survdays = NULL,
                           time_point = NULL, param.indiv = NULL,
                           param.ens = NULL, num.sel.path = NULL,
                           fs.AUC.cut = NULL, fs.sp.cut = NULL,
                           literature.features = NULL,
                           pathway.prior = NULL, prior.weight = 0,
                           prior.mode = c("none", "multiply", "weighted_sum"),
                           prior.cut = NULL, prior.default = NA_real_,
                           prior.clip = c(0, 1),
                           redundancy.thresh = NULL,
                           llm.mode = c("none", "curate", "prior"),
                           llm.curate.action = c("match_pathways", "add_signatures"),
                           disease = NULL,
                           outcome.context = NULL,
                           llm.provider = c("openai", "custom"),
                           llm.model = Sys.getenv("CASPEN_OPENAI_MODEL", "gpt-4.1-mini"),
                           llm.api.key = Sys.getenv("OPENAI_API_KEY"),
                           llm.fn = NULL,
                           llm.n.signatures = 10,
                           llm.max.genes = 80,
                           llm.min.genes = 5,
                           llm.pubmed = TRUE,
                           llm.pubmed.max = 8,
                           llm.email = NULL,
                           llm.valid.genes = NULL,
                           llm.prior.method = c("llm", "pubmed_count", "pubtator", "hybrid"),
                           llm.hybrid.weights = c(pubmed_count = 0.5, llm = 0.5),
                           llm.temperature = 0,
                           celltypes = NULL,
                           celltype.sep = c("_", ".", "-", ":", "|"),
                           celltype.position = c("auto", "suffix", "prefix"),
                           celltype.aliases = NULL,
                           continuous.metric = c("r2", "cindex", "mae", "rmse"),
                           auto.tune = FALSE,
                           tune.method = c("random", "grid", "successive_halving",
                                           "hyperband", "bayes"),
                           tune.n = 5,
                           tune.folds = 3,
                           tune.iter = 1,
                           tune.models = NULL,
                           parallel.iter = FALSE,
                           parallel.pathways = FALSE,
                           workers = future::availableCores(),
                           future.seed = NULL,
                           future.strategy = "multisession") {

  AUC.mat <- CI.low.mat <- CI.up.mat <- SP.98.mat <- SP.95.mat <- NULL
  SP.95.CI.low.mat <- SP.95.CI.up.mat <- NULL
  CINDEX.mat <- CINDEX.CI.low.mat <- CINDEX.CI.up.mat <- NULL
  R2.mat <- RMSE.mat <- MAE.mat <- CONT.CINDEX.mat <- NULL
  prior.mode <- match.arg(prior.mode)
  if (identical(llm.mode, "both")) {
    warning("llm.mode = 'both' is deprecated. Use llm.mode = 'curate'; ",
            "'curate' now maps concepts and scores priors.")
    llm.mode <- "curate"
  }
  llm.mode <- match.arg(llm.mode)
  if (identical(llm.curate.action, "both")) {
    warning("llm.curate.action = 'both' is deprecated. Using 'match_pathways'.")
    llm.curate.action <- "match_pathways"
  }
  llm.curate.action <- match.arg(llm.curate.action)
  llm.provider <- match.arg(llm.provider)
  llm.prior.method <- match.arg(llm.prior.method)
  celltype.position <- match.arg(celltype.position)
  continuous.metric <- match.arg(continuous.metric)
  tune.method <- match.arg(tune.method)

  type_out <- tolower(outcome.type)
  if (!type_out %in% c("binary", "categorical", "survival", "continuous")) {
    stop("Invalid outcome.type.")
  }

  llm.result <- NULL
  llm.concept.map <- NULL
  if (llm.mode != "none") {
    llm.result <- llm_literature_signatures(
      disease = disease,
      outcome.context = outcome.context,
      features = features,
      celltypes = celltypes,
      mode = llm.mode,
      n.signatures = llm.n.signatures,
      max.genes = llm.max.genes,
      min.genes = llm.min.genes,
      provider = llm.provider,
      model = llm.model,
      api.key = llm.api.key,
      llm.fn = llm.fn,
      pubmed = llm.pubmed,
      pubmed.max = llm.pubmed.max,
      email = llm.email,
      valid.genes = llm.valid.genes,
      prior.method = llm.prior.method,
      hybrid.weights = llm.hybrid.weights,
      prior.clip = prior.clip,
      temperature = llm.temperature
    )
    if (llm.mode == "curate" &&
        llm.curate.action == "match_pathways") {
      matched <- caspen_match_curated_concepts(features, llm.result)
      llm.concept.map <- matched$map
      if (!is.null(llm.concept.map) && nrow(llm.concept.map)) {
        features <- matched$features
        concept.prior <- normalize_pathway_prior(llm.result$pathway.prior)
        if (length(concept.prior)) {
          clean.prior <- clean_signature_name(names(concept.prior))
          mapped.prior <- stats::setNames(rep(NA_real_, length(features)),
                                          names(features))
          for (pw in names(features)) {
            concepts <- unique(llm.concept.map$llm.concept[llm.concept.map$pathway == pw])
            hit <- match(clean_signature_name(concepts), clean.prior)
            vals <- as.numeric(concept.prior[hit[!is.na(hit)]])
            vals <- vals[is.finite(vals)]
            if (length(vals)) mapped.prior[[pw]] <- max(vals, na.rm = TRUE)
          }
          pathway.prior <- merge_pathway_prior(mapped.prior, pathway.prior)
        }
        message("LLM curated concepts matched ", length(features),
                " user-provided pathway(s) for modeling.")
      } else {
        stop("No user-provided pathways matched the LLM-curated concepts. ",
             "Use llm.curate.action = 'add_signatures' to model LLM gene sets directly, ",
             "or use llm.mode = 'prior' to score all supplied pathways.")
      }
    }
    if (llm.mode == "curate" &&
        llm.curate.action == "add_signatures" &&
        length(llm.result$literature.features)) {
      literature.features <- c(literature.features %||% list(),
                               llm.result$literature.features)
    }
    if (llm.mode == "prior" &&
        length(llm.result$pathway.prior)) {
      pathway.prior <- merge_pathway_prior(llm.result$pathway.prior,
                                           pathway.prior)
    }
  }

  features <- pool_pathway_features(features, literature.features,
                                    redundancy.thresh)
  features <- expand_celltype_pathways(features, colnames(data), celltypes,
                                       celltype.sep, celltype.position,
                                       celltype.aliases)
  celltype.map <- attr(features, "celltype.map")
  path.names <- names(features)
  if (is.null(path.names)) path.names <- paste0("Pathway_", seq_along(features))
  names(features) <- path.names
  pathway.prior.vec <- align_pathway_prior(path.names, pathway.prior,
                                           prior.default, prior.clip)

  run_one_pathway <- function(f) {
    message("Pathway Number: ", f)
    com <- extract_pathway_genes(features[[f]])
    data_to_pass <- data[, colnames(data) %in% com, drop = FALSE]
    colnames(data_to_pass) <- make.names(colnames(data_to_pass), unique = TRUE)
    inner_parallel_iter <- isTRUE(parallel.iter) && !isTRUE(parallel.pathways)
    tuned.param.indiv <- param.indiv
    tuning.table <- NULL
    if (isTRUE(auto.tune) && ncol(data_to_pass) > 0) {
      tuned <- caspen_auto_tune_params(
        outcome.type = type_out,
        outcome = outcome,
        x = data_to_pass,
        models = models.indiv,
        param.indiv = param.indiv,
        tune.method = tune.method,
        tune.n = tune.n,
        tune.folds = tune.folds,
        tune.iter = tune.iter,
        tune.models = tune.models,
        survdays = survdays,
        time_point = time_point,
        continuous.metric = continuous.metric,
        seed = f
      )
      tuned.param.indiv <- tuned$param.indiv
      tuning.table <- tuned$tuning.table
      if (!is.null(tuning.table)) tuning.table$pathway <- path.names[[f]]
    }

    if (type_out == "binary") {
      res <- run_model_iterations(
        binary_model, iter = iter, workers = workers,
        parallel.iter = inner_parallel_iter, future.seed = future.seed,
        future.strategy = future.strategy,
        ensemble = ensemble,
        outcome = outcome, x = data_to_pass, models.indiv = models.indiv,
        num.folds = num.folds, models.ens = models.ens,
        param.indiv = tuned.param.indiv, param.ens = param.ens
      )
    } else if (type_out == "survival") {
      if (is.null(survdays) || is.null(time_point))
        stop("survdays and time_point required")
      res <- run_model_iterations(
        survival_model, iter = iter, workers = workers,
        parallel.iter = inner_parallel_iter, future.seed = future.seed,
        future.strategy = future.strategy,
        ensemble = ensemble,
        outcome = outcome, survdays = survdays, x = data_to_pass,
        models.indiv = models.indiv, num.folds = num.folds,
        models.ens = models.ens, param.indiv = tuned.param.indiv,
        param.ens = param.ens, time_point = time_point
      )
    } else if (type_out == "continuous") {
      res <- run_model_iterations(
        continuous_model, iter = iter, workers = workers,
        parallel.iter = inner_parallel_iter, future.seed = future.seed,
        future.strategy = future.strategy,
        ensemble = ensemble,
        outcome = outcome, x = data_to_pass, models.indiv = models.indiv,
        num.folds = num.folds, models.ens = models.ens,
        param.indiv = tuned.param.indiv, param.ens = param.ens
      )
    } else {
      res <- run_model_iterations(
        categorical_model, iter = iter, workers = workers,
        parallel.iter = inner_parallel_iter, future.seed = future.seed,
        future.strategy = future.strategy,
        ensemble = ensemble,
        outcome = outcome, x = data_to_pass, models.indiv = models.indiv,
        num.folds = num.folds, models.ens = models.ens,
        param.indiv = tuned.param.indiv, param.ens = param.ens
      )
    }
    message("Done pathway training")

    result <- if (type_out == "binary") {
      binary_roc(outcome, res, models.indiv, iter, ensemble, models.ens)
    } else if (type_out == "survival") {
      survival_roc(outcome, survdays, time_point, res, models.indiv, iter,
                   ensemble, models.ens)
    } else if (type_out == "continuous") {
      continuous_metrics(outcome, res, models.indiv, iter, ensemble, models.ens)
    } else {
      categorical_roc(res, outcome, models.indiv, iter, ensemble, models.ens)
    }
    message("ROC done")
    attr(result, "tuning") <- tuning.table
    result
  }

  if (isTRUE(parallel.pathways) && length(features) > 1) {
    path.workers <- max(1, min(as.integer(workers %||% future::availableCores()),
                               length(features)))
    message("Parallel pathways: ", length(features),
            " pathway tasks split across ", path.workers,
            " future workers.")
    if (!is.null(future.strategy)) {
      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      do.call(future::plan, c(list(future.strategy), list(workers = path.workers)))
    }
    pathway.results <- future.apply::future_lapply(
      seq_along(features),
      run_one_pathway,
      future.seed = future.seed
    )
  } else {
    pathway.results <- lapply(seq_along(features), run_one_pathway)
  }
  tuning.results <- lapply(pathway.results, attr, "tuning")
  names(tuning.results) <- path.names

  bind_result_metric <- function(metric.name) {
    mat <- do.call(rbind, lapply(pathway.results, `[[`, metric.name))
    if (is.null(dim(mat))) mat <- matrix(mat, nrow = length(pathway.results), byrow = TRUE)
    mat
  }

  AUC.mat <- bind_result_metric("auc")
  CI.low.mat <- bind_result_metric("ci.low")
  CI.up.mat <- bind_result_metric("ci.up")
  SP.98.mat <- bind_result_metric("sp.98")
  SP.95.mat <- bind_result_metric("sp.95")
  SP.95.CI.low.mat <- bind_result_metric("sp.95.ci.low")
  SP.95.CI.up.mat <- bind_result_metric("sp.95.ci.up")
  if (type_out == "survival") {
    CINDEX.mat <- bind_result_metric("cindex")
    CINDEX.CI.low.mat <- bind_result_metric("cindex.ci.low")
    CINDEX.CI.up.mat <- bind_result_metric("cindex.ci.up")
  } else if (type_out == "continuous") {
    R2.mat <- bind_result_metric("r2")
    RMSE.mat <- bind_result_metric("rmse")
    MAE.mat <- bind_result_metric("mae")
    CONT.CINDEX.mat <- bind_result_metric("cindex")
  }

  # Compute mean metrics
  mean.auc  <- rowMeans(AUC.mat, na.rm = TRUE)
  mean.sp95 <- rowMeans(SP.95.mat, na.rm = TRUE)
  primary.metric <- mean.auc
  if (type_out == "survival" && !is.null(CINDEX.mat)) {
    primary.metric <- rowMeans(CINDEX.mat, na.rm = TRUE)
  } else if (type_out == "continuous" && !is.null(R2.mat)) {
    cont.metric.mat <- switch(
      continuous.metric,
      r2 = R2.mat,
      cindex = CONT.CINDEX.mat,
      mae = -MAE.mat,
      rmse = -RMSE.mat
    )
    primary.metric <- rowMeans(cont.metric.mat, na.rm = TRUE)
  }
  selection.score <- primary.metric

  if (!is.null(num.sel.path) || !is.null(fs.AUC.cut) ||
      !is.null(fs.sp.cut) || !is.null(prior.cut)) {
    pass <- seq_along(primary.metric)
    if (!is.null(fs.AUC.cut)) {
      pass <- pass[primary.metric[pass] > fs.AUC.cut]
    }
    if (!is.null(fs.sp.cut)) {
      pass <- pass[mean.sp95[pass] > fs.sp.cut]
    }
    if (!is.null(prior.cut)) {
      prior.pass <- is.finite(pathway.prior.vec[pass]) &
        pathway.prior.vec[pass] >= prior.cut
      pass <- pass[prior.pass]
    }
    if (length(pass) > 0) {
      topk <- order(primary.metric[pass], decreasing = TRUE, na.last = TRUE)
      pass <- pass[topk]
      if (!is.null(num.sel.path)) {
        pass <- pass[seq_len(min(num.sel.path, length(pass)))]
      }
      keep <- pass
    } else {
      keep <- integer(0)
    }
  } else {
    # no selection → keep all pathways
    keep <- seq_len(nrow(AUC.mat))
  }

  # subset matrices
  AUC.mat    <- AUC.mat[keep, , drop = FALSE]
  CI.low.mat <- CI.low.mat[keep, , drop = FALSE]
  CI.up.mat  <- CI.up.mat[keep, , drop = FALSE]
  SP.98.mat  <- SP.98.mat[keep, , drop = FALSE]
  SP.95.mat  <- SP.95.mat[keep, , drop = FALSE]
  SP.95.CI.low.mat <- SP.95.CI.low.mat[keep, , drop = FALSE]
  SP.95.CI.up.mat <- SP.95.CI.up.mat[keep, , drop = FALSE]
  if (type_out == "survival") {
    CINDEX.mat <- CINDEX.mat[keep, , drop = FALSE]
    CINDEX.CI.low.mat <- CINDEX.CI.low.mat[keep, , drop = FALSE]
    CINDEX.CI.up.mat <- CINDEX.CI.up.mat[keep, , drop = FALSE]
  } else if (type_out == "continuous") {
    R2.mat <- R2.mat[keep, , drop = FALSE]
    RMSE.mat <- RMSE.mat[keep, , drop = FALSE]
    MAE.mat <- MAE.mat[keep, , drop = FALSE]
    CONT.CINDEX.mat <- CONT.CINDEX.mat[keep, , drop = FALSE]
  }

  labels <- path.names[keep]
  prior.keep <- pathway.prior.vec[keep]
  primary.keep <- primary.metric[keep]
  selection.keep <- selection.score[keep]
  mean.auc.keep <- mean.auc[keep]
  mean.sp95.keep <- mean.sp95[keep]
  prior.mat <- matrix(prior.keep, ncol = 1,
                      dimnames = list(labels, "prior.weight"))
  data.metric.mat <- if (type_out == "survival") {
    CINDEX.mat
  } else if (type_out == "continuous") {
    R2.mat
  } else {
    AUC.mat
  }

  rownames(AUC.mat)    <- labels
  rownames(CI.low.mat) <- labels
  rownames(CI.up.mat)  <- labels
  rownames(SP.98.mat)  <- labels
  rownames(SP.95.mat)  <- labels
  rownames(SP.95.CI.low.mat) <- labels
  rownames(SP.95.CI.up.mat) <- labels
  if (type_out == "survival") {
    rownames(CINDEX.mat) <- labels
    rownames(CINDEX.CI.low.mat) <- labels
    rownames(CINDEX.CI.up.mat) <- labels
    rownames(data.metric.mat) <- labels
  } else if (type_out == "continuous") {
    rownames(R2.mat) <- labels
    rownames(RMSE.mat) <- labels
    rownames(MAE.mat) <- labels
    rownames(CONT.CINDEX.mat) <- labels
    rownames(data.metric.mat) <- labels
  }

  all_models <- c(models.indiv,
                  if (ensemble) paste0(models.ens, ".ens"),
                  "SimpleAvg", "WeightedAvg")
  colnames(AUC.mat)    <- all_models
  colnames(CI.low.mat) <- all_models
  colnames(CI.up.mat)  <- all_models
  colnames(SP.98.mat)  <- all_models
  colnames(SP.95.mat)  <- all_models
  colnames(SP.95.CI.low.mat) <- all_models
  colnames(SP.95.CI.up.mat) <- all_models
  if (type_out == "survival") {
    colnames(CINDEX.mat) <- all_models
    colnames(CINDEX.CI.low.mat) <- all_models
    colnames(CINDEX.CI.up.mat) <- all_models
    colnames(data.metric.mat) <- all_models
  } else if (type_out == "continuous") {
    colnames(R2.mat) <- all_models
    colnames(RMSE.mat) <- all_models
    colnames(MAE.mat) <- all_models
    colnames(CONT.CINDEX.mat) <- all_models
    colnames(data.metric.mat) <- all_models
  }
  if (!type_out %in% c("survival", "continuous")) {
    data.metric.mat <- AUC.mat
  }
  celltype.map.keep <- NULL
  if (!is.null(celltype.map) && nrow(celltype.map)) {
    celltype.map.keep <- celltype.map[match(labels, celltype.map$label), ,
                                      drop = FALSE]
  }
  specificity.metric.mat <- data.metric.mat
  specificity.chance <- 0.5
  specificity.lower <- FALSE
  if (type_out == "continuous") {
    specificity.metric.mat <- switch(
      continuous.metric,
      r2 = R2.mat,
      cindex = CONT.CINDEX.mat,
      mae = MAE.mat,
      rmse = RMSE.mat
    )
    specificity.chance <- if (continuous.metric %in% c("r2", "mae", "rmse")) 0 else 0.5
    specificity.lower <- continuous.metric %in% c("mae", "rmse")
  }
  celltype.specificity <- if (!is.null(celltype.map.keep) && nrow(celltype.map.keep)) {
    celltype_specificity_summary(
      metric.mat = specificity.metric.mat,
      labels = labels,
      celltype.map = celltype.map.keep,
      chance = specificity.chance,
      lower.is.better = specificity.lower
    )
  } else {
    list(summary = NULL, by.model = NULL)
  }

  ret <- list(
    data.metric = data.metric.mat,
    auc    = AUC.mat,
    CI.low = CI.low.mat,
    CI.up  = CI.up.mat,
    SP.98  = SP.98.mat,
    SP.95  = SP.95.mat,
    SP.95.CI.low = SP.95.CI.low.mat,
    SP.95.CI.up = SP.95.CI.up.mat,
    Prior.weight = prior.mat,
    prior = prior.keep,
    primary.metric = primary.keep,
    selection.score = selection.keep,
    selection.table = data.frame(
      pathway = labels,
      primary.metric = primary.keep,
      primary.metric.name = if (type_out == "continuous") continuous.metric else
        if (type_out == "survival") "cindex" else "auc",
      mean.auc = mean.auc.keep,
      mean.sp95 = mean.sp95.keep,
      parent.pathway = if (!is.null(celltype.map.keep)) celltype.map.keep$pathway else labels,
      celltype = if (!is.null(celltype.map.keep)) celltype.map.keep$celltype else NA_character_,
      prior = prior.keep,
      row.names = labels,
      check.names = FALSE
    ),
    features = features[keep],
    label  = labels,
    celltype.map = celltype.map.keep,
    celltype.specificity = celltype.specificity$summary,
    celltype.specificity.by.model = celltype.specificity$by.model
  )
  if (isTRUE(auto.tune)) {
    ret$tuning <- tuning.results[keep]
    ret$tuning.table <- do.call(rbind, tuning.results[keep])
  }
  if (!is.null(llm.result)) ret$llm <- llm.result
  if (!is.null(llm.concept.map)) ret$llm.concept.map <- llm.concept.map
  if (type_out == "survival") {
    ret$C.index <- CINDEX.mat
    ret$C.index.CI.low <- CINDEX.CI.low.mat
    ret$C.index.CI.up <- CINDEX.CI.up.mat
  } else if (type_out == "continuous") {
    ret$R2 <- R2.mat
    ret$RMSE <- RMSE.mat
    ret$MAE <- MAE.mat
    ret$C.index <- CONT.CINDEX.mat
  }
  ret
}
