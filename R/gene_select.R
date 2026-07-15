#' Select Genes Within Pathways Using Native Model Importance
#'
#' This function ranks genes within selected pathways using native importance
#' from models that expose feature-level importance. Importance values are
#' converted to within-model percentiles before aggregation so scores from
#' different model families are not averaged on incompatible raw scales. Genes
#' are selected by user-specified stability cutoffs and then ranked by median
#' percentile.
#'
#' @param outcome.type Character string specifying the prediction task. Must be
#' one of `"survival"`, `"binary"`, `"categorical"`, or `"continuous"`.
#' @param outcome Outcome vector. For survival, this is the event/status vector.
#' @param data Numeric matrix or data frame with samples in rows and genes or
#' celltype-specific gene features in columns.
#' @param features Named list of selected pathways. Each element should contain
#' gene or feature column names. If `celltypes` is supplied, bulk gene symbols
#' are expanded to celltype-specific columns using the same mapping rules as
#' `pathway_select()`.
#' @param models Character vector of native-importance models to use. Supported
#' values are `"XG"`, `"RF"`, `"GB"`, `"EN"`, `"MB"`, `"DCT"`, and `"NN"`.
#' @param iter Integer number of repeated fold splits.
#' @param num.folds Integer number of folds per iteration.
#' @param survdays Numeric survival time vector required when
#' `outcome.type = "survival"`.
#' @param min.pathway.size Minimum pathway size required to run gene selection.
#' Pathways with `min.pathway.size` genes or fewer are returned unchanged.
#' @param gene.stability.quantile Percentile threshold used to define whether a
#' gene is selected in one model/fold/iteration run. Default is 0.75.
#' @param gene.min.stability Minimum proportion of model/fold/iteration runs in
#' which a gene must meet `gene.stability.quantile`. Default is 0.25.
#' @param max.genes Maximum number of genes to keep for pathways where gene
#' selection is run. Defaults to `min.pathway.size`.
#' @param param.models Optional named list of model-specific parameters.
#' @param celltypes Optional character vector of cell type names. When supplied,
#' bulk pathway genes are expanded into pathway-celltype feature columns.
#' @param celltype.sep Character vector of separators used to match
#' celltype-specific expression columns.
#' @param celltype.position One of `"auto"`, `"suffix"`, or `"prefix"`.
#' @param celltype.aliases Optional named list keyed by entries in `celltypes`;
#' values are aliases used to match pathway names to cell types.
#' @param suppress.model.warnings Logical; if `TRUE`, suppresses warnings from
#' individual native-importance fits. This is useful for small folds where
#' elastic-net models may warn that one class has fewer than eight observations.
#' @param gene.prior.method Character string specifying optional literature
#' prior scoring for genes. `"none"` skips gene priors, `"pubmed_count"` uses
#' normalized PubMed counts for disease + outcome + gene, `"pubtator"` uses
#' PubTator3 gene-annotation evidence, `"llm"` uses LLM-assigned gene
#' relevance, and `"hybrid"` combines available components.
#' Gene priors are reported as additional columns and do not change the default
#' stability-based gene selection.
#' @param gene.hybrid.weights Named numeric vector used when
#' `gene.prior.method = "hybrid"`. Supported names are `"pubmed_count"`,
#' `"pubtator"`, `"llm"`, `"evidence"`, `"model_importance"`, and
#' `"stability"`.
#' @param disease Disease context used for gene-level PubMed/LLM priors.
#' @param outcome.context Outcome context used for gene-level PubMed/LLM priors.
#' @param llm.provider LLM provider for gene priors. Currently `"openai"` or
#' `"custom"`.
#' @param llm.model OpenAI model name.
#' @param llm.api.key OpenAI API key. Defaults to `OPENAI_API_KEY`.
#' @param llm.fn Optional custom LLM function accepting a prompt and returning
#' text, for example an Ollama wrapper.
#' @param llm.pubmed Logical; if `TRUE`, retrieves PubMed abstracts as LLM
#' grounding context.
#' @param llm.email Optional email passed to NCBI E-utilities.
#' @param llm.temperature OpenAI sampling temperature.
#' @param verbose Logical; if `TRUE`, prints pathway, iteration, fold, and model.
#'
#' @return A list with selected pathway gene lists, a gene ranking table, and
#' raw percentile records used to compute stability.
#'
#' @examples
#' # -------------------------------------------------------------------------
#' # Basic gene selection inside selected pathways
#' # -------------------------------------------------------------------------
#' data(caspen_example)
#' gs <- gene_select(
#'   outcome.type = "binary",
#'   outcome = caspen_example$outcome,
#'   data = caspen_example$data_to_pass,
#'   features = caspen_example$pathways[1],
#'   models = c("RF", "EN", "DCT"),
#'   iter = 1,
#'   num.folds = 2,
#'   min.pathway.size = 5,
#'   max.genes = 5
#' )
#'
#' gs$selected.features
#' gs$ranking
#' gs$settings
#'
#' \dontrun{
#' # -----------------------------------------------------------------------
#' # Typical workflow: pathway_select() first, then gene_select()
#' # -----------------------------------------------------------------------
#' path_res <- pathway_select(
#'   outcome.type = "binary",
#'   outcome = caspen_example$outcome,
#'   data = caspen_example$data_to_pass,
#'   features = caspen_example$pathways[1:10],
#'   models.indiv = c("RF", "EN", "DCT"),
#'   iter = 10,
#'   num.folds = 3,
#'   fs.AUC.cut = 0.70,
#'   fs.sp.cut = 0.10
#' )
#'
#' selected_pathways <- path_res$features
#'
#' gene_res <- gene_select(
#'   outcome.type = "binary",
#'   outcome = caspen_example$outcome,
#'   data = caspen_example$data_to_pass,
#'   features = selected_pathways,
#'   models = c("RF", "GB", "EN", "DCT"),
#'   iter = 10,
#'   num.folds = 3,
#'   min.pathway.size = 30,
#'   gene.stability.quantile = 0.75,
#'   gene.min.stability = 0.25,
#'   max.genes = 30
#' )
#'
#' gene_res$selected.features
#' head(gene_res$ranking)
#'
#' # -----------------------------------------------------------------------
#' # Add literature priors for genes using PubMed counts
#' # -----------------------------------------------------------------------
#' # This does not change the default selection rule. It appends
#' # pubmed.count, pubmed.score, and hybrid/priority columns so the user can
#' # decide their own cutoffs after seeing model stability and literature
#' # support side by side.
#' gene_pubmed <- gene_select(
#'   outcome.type = "binary",
#'   outcome = caspen_example$outcome,
#'   data = caspen_example$data_to_pass,
#'   features = caspen_example$pathways[1:2],
#'   models = c("RF", "EN", "DCT"),
#'   iter = 5,
#'   num.folds = 3,
#'   min.pathway.size = 30,
#'   max.genes = 30,
#'   gene.prior.method = "pubmed_count",
#'   disease = "high-grade serous ovarian cancer",
#'   outcome.context = "platinum resistance",
#'   llm.email = "your.email@example.com"
#' )
#'
#' gene_pubmed$ranking[, c("pathway", "gene", "median_percentile",
#'                         "stability", "pubmed.count", "pubmed.score")]
#'
#' # PubTator3 gene priors use PubMed IDs plus PubTator gene annotations.
#' gene_pubtator <- gene_select(
#'   outcome.type = "binary",
#'   outcome = caspen_example$outcome,
#'   data = caspen_example$data_to_pass,
#'   features = caspen_example$pathways[1:2],
#'   models = c("RF", "EN", "DCT"),
#'   iter = 5,
#'   num.folds = 3,
#'   min.pathway.size = 30,
#'   max.genes = 30,
#'   gene.prior.method = "pubtator",
#'   disease = "high-grade serous ovarian cancer",
#'   outcome.context = "platinum resistance",
#'   llm.email = "your.email@example.com"
#' )
#'
#' gene_pubtator$ranking[, c("pathway", "gene", "pubtator.count",
#'                           "pubtator.score", "hybrid.score")]
#'
#' # -----------------------------------------------------------------------
#' # Hybrid gene ranking with model importance, stability, PubMed, and LLM
#' # -----------------------------------------------------------------------
#' Sys.setenv(OPENAI_API_KEY = "sk-proj-your_real_key_here")
#'
#' gene_hybrid <- gene_select(
#'   outcome.type = "binary",
#'   outcome = caspen_example$outcome,
#'   data = caspen_example$data_to_pass,
#'   features = caspen_example$pathways[1:2],
#'   models = c("RF", "GB", "EN", "DCT"),
#'   iter = 5,
#'   num.folds = 3,
#'   min.pathway.size = 30,
#'   max.genes = 30,
#'   gene.prior.method = "hybrid",
#'   gene.hybrid.weights = c(
#'     model_importance = 0.5,
#'     stability = 0.2,
#'     pubmed_count = 0.1,
#'     pubtator = 0.1,
#'     llm = 0.1
#'   ),
#'   disease = "high-grade serous ovarian cancer",
#'   outcome.context = "platinum resistance",
#'   llm.pubmed = FALSE,
#'   llm.model = "gpt-4o-mini",
#'   llm.temperature = NULL
#' )
#'
#' # The user controls the final cutoff. For example:
#' top_genes <- subset(
#'   gene_hybrid$ranking,
#'   selected & stability >= 0.25 & hybrid.score >= 0.6
#' )
#' top_genes[order(top_genes$pathway, -top_genes$hybrid.score), ]
#'
#' # -----------------------------------------------------------------------
#' # Local free LLM option through Ollama
#' # -----------------------------------------------------------------------
#' # Terminal:
#' #   ollama pull qwen2.5:7b
#' # Keep Ollama running while R calls http://localhost:11434.
#' ollama_llm <- function(prompt) {
#'   req <- httr2::request("http://localhost:11434/api/generate") |>
#'     httr2::req_body_json(list(
#'       model = "qwen2.5:7b",
#'       prompt = prompt,
#'       stream = FALSE
#'     ))
#'   res <- httr2::req_perform(req)
#'   out <- jsonlite::fromJSON(httr2::resp_body_string(res))
#'   out$response
#' }
#'
#' gene_ollama <- gene_select(
#'   outcome.type = "binary",
#'   outcome = caspen_example$outcome,
#'   data = caspen_example$data_to_pass,
#'   features = caspen_example$pathways[1:2],
#'   models = c("RF", "EN", "DCT"),
#'   iter = 5,
#'   num.folds = 3,
#'   min.pathway.size = 30,
#'   gene.prior.method = "llm",
#'   disease = "high-grade serous ovarian cancer",
#'   outcome.context = "platinum resistance",
#'   llm.provider = "custom",
#'   llm.fn = ollama_llm
#' )
#' }
#'
#' @export
gene_select <- function(outcome.type, outcome, data, features,
                        models = c("XG", "RF", "GB", "EN", "MB", "DCT", "NN"),
                        iter = 1, num.folds = NULL, survdays = NULL,
                        min.pathway.size = 30,
                        gene.stability.quantile = 0.75,
                        gene.min.stability = 0.25,
                        max.genes = min.pathway.size,
                        param.models = NULL,
                        auto.tune = FALSE,
                        tune.method = c("random", "grid", "successive_halving",
                                        "hyperband", "bayes"),
                        tune.n = 5,
                        tune.folds = 3,
                        tune.iter = 1,
                        tune.models = NULL,
                        celltypes = NULL,
                        celltype.sep = c("_", ".", "-", ":", "|"),
                        celltype.position = c("auto", "suffix", "prefix"),
                        celltype.aliases = NULL,
                        suppress.model.warnings = TRUE,
                        gene.prior.method = c("none", "pubmed_count", "pubtator", "llm", "hybrid"),
                        gene.hybrid.weights = c(model_importance = 0.5,
                                                stability = 0.2,
                                                pubmed_count = 0.1,
                                                pubtator = 0.1,
                                                llm = 0.1),
                        disease = NULL,
                        outcome.context = NULL,
                        llm.provider = c("openai", "custom"),
                        llm.model = Sys.getenv("CASPEN_OPENAI_MODEL", "gpt-4.1-mini"),
                        llm.api.key = Sys.getenv("OPENAI_API_KEY"),
                        llm.fn = NULL,
                        llm.pubmed = FALSE,
                        llm.email = NULL,
                        llm.temperature = 0,
                        verbose = TRUE) {
  type_out <- tolower(outcome.type)
  if (!type_out %in% c("binary", "categorical", "survival", "continuous")) {
    stop("Invalid outcome.type.")
  }
  if (type_out == "survival" && is.null(survdays)) {
    stop("survdays is required when outcome.type = 'survival'.")
  }
  data <- as.data.frame(data)
  models <- unique(toupper(models))
  supported <- c("XG", "RF", "GB", "EN", "MB", "DCT", "NN")
  models <- intersect(models, supported)
  if (!length(models)) stop("No supported native-importance models supplied.")
  if (is.null(num.folds)) num.folds <- 3
  if (!is.list(features)) stop("features must be a list of pathways.")
  tune.method <- match.arg(tune.method)
  gene.prior.method <- match.arg(gene.prior.method)
  llm.provider <- match.arg(llm.provider)
  if (gene.prior.method != "none") {
    if (!nzchar(disease %||% "")) {
      stop("disease is required when gene.prior.method != 'none'.")
    }
    if (!nzchar(outcome.context %||% "")) {
      stop("outcome.context is required when gene.prior.method != 'none'.")
    }
  }

  features <- lapply(features, extract_pathway_genes)
  if (is.null(names(features))) {
    names(features) <- paste0("Pathway_", seq_along(features))
  }
  features <- expand_celltype_pathways(
    features, colnames(data), celltypes, celltype.sep,
    match.arg(celltype.position), celltype.aliases
  )

  selected <- list()
  ranking_tables <- list()
  raw_records <- list()
  tuning_records <- list()

  for (pathway_name in names(features)) {
    genes <- unique(features[[pathway_name]])
    genes <- genes[genes %in% colnames(data)]
    if (!length(genes)) next
    pathway.params <- param.models %||% list()
    if (isTRUE(auto.tune)) {
      tuned <- caspen_auto_tune_params(
        outcome.type = type_out,
        outcome = outcome,
        x = data[, genes, drop = FALSE],
        models = models,
        param.indiv = pathway.params,
        tune.method = tune.method,
        tune.n = tune.n,
        tune.folds = tune.folds,
        tune.iter = tune.iter,
        tune.models = tune.models,
        survdays = if (type_out == "survival") survdays else NULL,
        seed = match(pathway_name, names(features))
      )
      pathway.params <- tuned$param.indiv
      if (!is.null(tuned$tuning.table)) {
        tuned$tuning.table$pathway <- pathway_name
        tuning_records[[pathway_name]] <- tuned$tuning.table
      }
    }

    if (length(genes) <= min.pathway.size) {
      tab <- data.frame(
        pathway = pathway_name,
        gene = genes,
        median_percentile = NA_real_,
        mean_percentile = NA_real_,
        stability = NA_real_,
        rank = seq_along(genes),
        selected = TRUE,
        skipped = TRUE,
        stringsAsFactors = FALSE
      )
      tab <- add_gene_prior_scores(
        tab = tab,
        gene.prior.method = gene.prior.method,
        gene.hybrid.weights = gene.hybrid.weights,
        disease = disease,
        outcome.context = outcome.context,
        llm.provider = llm.provider,
        llm.model = llm.model,
        llm.api.key = llm.api.key,
        llm.fn = llm.fn,
        llm.pubmed = llm.pubmed,
        llm.email = llm.email,
        llm.temperature = llm.temperature
      )
      selected[[pathway_name]] <- genes
      ranking_tables[[pathway_name]] <- tab
      next
    }

    records <- list()
    rec_i <- 1
    for (it in seq_len(iter)) {
      fold_y <- if (type_out == "survival") outcome else outcome
      folds <- balanced.folds.vec(fold_y, num.folds)
      for (fold in seq_len(num.folds)) {
        train_idx <- which(folds != fold)
        xtrain <- data[train_idx, genes, drop = FALSE]
        ytrain <- outcome[train_idx]
        ttrain <- if (type_out == "survival") survdays[train_idx] else NULL

        for (model in models) {
          if (verbose) {
            message("Pathway: ", pathway_name, " | Iteration: ", it,
                    " | Fold: ", fold, " | Model: ", model)
          }
          imp <- native_gene_importance(
            type_out, model, xtrain, ytrain, ttrain,
            pathway.params[[model]] %||% list(),
            suppress.warnings = suppress.model.warnings
          )
          if (is.null(imp)) next
          imp <- imp[genes]
          imp[is.na(imp)] <- 0
          pct <- importance_percentile(imp)
          records[[rec_i]] <- data.frame(
            pathway = pathway_name,
            gene = names(pct),
            model = model,
            iteration = it,
            fold = fold,
            percentile = as.numeric(pct),
            stringsAsFactors = FALSE
          )
          rec_i <- rec_i + 1
        }
      }
    }

    if (!length(records)) {
      tab <- data.frame(
        pathway = pathway_name,
        gene = genes,
        median_percentile = NA_real_,
        mean_percentile = NA_real_,
        stability = NA_real_,
        rank = seq_along(genes),
        selected = TRUE,
        skipped = TRUE,
        stringsAsFactors = FALSE
      )
      tab <- add_gene_prior_scores(
        tab = tab,
        gene.prior.method = gene.prior.method,
        gene.hybrid.weights = gene.hybrid.weights,
        disease = disease,
        outcome.context = outcome.context,
        llm.provider = llm.provider,
        llm.model = llm.model,
        llm.api.key = llm.api.key,
        llm.fn = llm.fn,
        llm.pubmed = llm.pubmed,
        llm.email = llm.email,
        llm.temperature = llm.temperature
      )
      selected[[pathway_name]] <- genes
      ranking_tables[[pathway_name]] <- tab
      next
    }

    raw <- do.call(rbind, records)
    raw_records[[pathway_name]] <- raw
    tab <- gene_importance_summary(raw, genes, gene.stability.quantile)
    tab$pathway <- pathway_name
    tab$rank <- rank(-tab$median_percentile, ties.method = "first")
    tab <- tab[order(tab$rank), c(
      "pathway", "gene", "median_percentile", "mean_percentile",
      "stability", "rank"
    )]

    pass <- which(tab$stability >= gene.min.stability)
    keep_genes <- character(0)
    if (length(pass)) {
      keep_n <- min(max.genes, length(pass))
      keep_genes <- tab$gene[pass][seq_len(keep_n)]
    }
    tab$selected <- tab$gene %in% keep_genes
    tab$skipped <- FALSE
    tab <- add_gene_prior_scores(
      tab = tab,
      gene.prior.method = gene.prior.method,
      gene.hybrid.weights = gene.hybrid.weights,
      disease = disease,
      outcome.context = outcome.context,
      llm.provider = llm.provider,
      llm.model = llm.model,
      llm.api.key = llm.api.key,
      llm.fn = llm.fn,
      llm.pubmed = llm.pubmed,
      llm.email = llm.email,
      llm.temperature = llm.temperature
    )

    selected[[pathway_name]] <- keep_genes
    ranking_tables[[pathway_name]] <- tab
  }

  list(
    selected.features = selected,
    ranking = do.call(rbind, ranking_tables),
    raw.percentiles = raw_records,
    settings = list(
      models = models,
      iter = iter,
      num.folds = num.folds,
      min.pathway.size = min.pathway.size,
      gene.stability.quantile = gene.stability.quantile,
      gene.min.stability = gene.min.stability,
      max.genes = max.genes,
      gene.prior.method = gene.prior.method,
      gene.hybrid.weights = gene.hybrid.weights,
      auto.tune = auto.tune,
      tune.method = tune.method,
      tune.n = tune.n,
      tune.folds = tune.folds,
      tune.iter = tune.iter,
      suppress.model.warnings = suppress.model.warnings
    ),
    tuning = if (length(tuning_records)) do.call(rbind, tuning_records) else NULL
  )
}

#' @keywords internal
#' @noRd
add_gene_prior_scores <- function(tab, gene.prior.method = "none",
                                  gene.hybrid.weights = c(model_importance = 0.5,
                                                          stability = 0.2,
                                                          pubmed_count = 0.1,
                                                          pubtator = 0.1,
                                                          llm = 0.1),
                                  disease = NULL, outcome.context = NULL,
                                  llm.provider = c("openai", "custom"),
                                  llm.model = Sys.getenv("CASPEN_OPENAI_MODEL", "gpt-4.1-mini"),
                                  llm.api.key = Sys.getenv("OPENAI_API_KEY"),
                                  llm.fn = NULL,
                                  llm.pubmed = FALSE,
                                  llm.email = NULL,
                                  llm.temperature = 0) {
  gene.prior.method <- match.arg(gene.prior.method,
                                 c("none", "pubmed_count", "pubtator", "llm", "hybrid"))
  if (gene.prior.method == "none" || !nrow(tab)) return(tab)

  genes <- unique(as.character(tab$gene))
  genes <- genes[nzchar(genes)]
  need.pubmed <- gene.prior.method == "pubmed_count" ||
    (gene.prior.method == "hybrid" &&
       "pubmed_count" %in% names(caspen_normalize_hybrid_weights(gene.hybrid.weights)))
  need.pubtator <- gene.prior.method == "pubtator" ||
    (gene.prior.method == "hybrid" &&
       "pubtator" %in% names(caspen_normalize_hybrid_weights(gene.hybrid.weights)))
  need.llm <- gene.prior.method %in% c("llm", "hybrid")

  pubmed.df <- data.frame(pathway = genes, pubmed.count = NA_real_,
                          pubmed.score = NA_real_, stringsAsFactors = FALSE)
  if (isTRUE(need.pubmed)) {
    pubmed.df <- caspen_pubmed_prior_scores(
      disease = disease,
      outcome.context = outcome.context,
      pathway.names = genes,
      email = llm.email
    )
  }

  pubtator.df <- data.frame(pathway = genes, pubtator.count = NA_real_,
                            pubtator.score = NA_real_, stringsAsFactors = FALSE)
  if (isTRUE(need.pubtator)) {
    gene.features <- stats::setNames(as.list(genes), genes)
    pubtator.df <- caspen_pubtator_prior_scores(
      disease = disease,
      outcome.context = outcome.context,
      features = gene.features,
      pathway.names = genes,
      email = llm.email
    )
  }

  llm.df <- data.frame(pathway = genes, llm.score = NA_real_,
                       evidence.score = NA_real_, stringsAsFactors = FALSE)
  if (isTRUE(need.llm)) {
    gene.features <- stats::setNames(as.list(genes), genes)
    llm.res <- llm_literature_signatures(
      disease = disease,
      outcome.context = outcome.context,
      features = gene.features,
      mode = "prior",
      provider = match.arg(llm.provider),
      model = llm.model,
      api.key = llm.api.key,
      llm.fn = llm.fn,
      pubmed = llm.pubmed,
      email = llm.email,
      prior.method = "llm",
      temperature = llm.temperature
    )
    if (nrow(llm.res$evidence)) {
      llm.df <- llm.res$evidence[, c("pathway", "llm.score", "evidence.score"),
                                 drop = FALSE]
      names(llm.df)[1] <- "pathway"
    }
  }

  tab$pubmed.count <- pubmed.df$pubmed.count[match(tab$gene, pubmed.df$pathway)]
  tab$pubmed.score <- pubmed.df$pubmed.score[match(tab$gene, pubmed.df$pathway)]
  tab$pubtator.count <- pubtator.df$pubtator.count[match(tab$gene, pubtator.df$pathway)]
  tab$pubtator.score <- pubtator.df$pubtator.score[match(tab$gene, pubtator.df$pathway)]
  tab$llm.score <- llm.df$llm.score[match(tab$gene, llm.df$pathway)]
  tab$evidence.score <- llm.df$evidence.score[match(tab$gene, llm.df$pathway)]
  tab$model_importance.score <- tab$median_percentile
  tab$hybrid.score <- NA_real_
  tab$gene.prior.method <- gene.prior.method

  for (i in seq_len(nrow(tab))) {
    components <- c(
      pubmed_count = unname(tab$pubmed.score[i]),
      pubtator = unname(tab$pubtator.score[i]),
      llm = unname(tab$llm.score[i]),
      evidence = unname(tab$evidence.score[i]),
      model_importance = unname(tab$model_importance.score[i]),
      stability = unname(tab$stability[i])
    )
    tab$hybrid.score[i] <- caspen_combine_prior_components(
      prior.method = if (gene.prior.method == "hybrid") "hybrid" else gene.prior.method,
      components = components,
      hybrid.weights = gene.hybrid.weights
    )
  }
  tab
}
