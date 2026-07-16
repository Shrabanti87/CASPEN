#' Curate Literature Signatures and Priors With an LLM
#'
#' Queries an LLM for disease/outcome-specific literature signatures and/or
#' pathway prior weights. This function is used by `pathway_select()` when
#' `llm.mode` is enabled, but can also be called directly so users can inspect
#' and save the curated signatures or prior table before model fitting.
#'
#' @details
#' CASPEN supports four prior-weight strategies. With `prior.method = "llm"`,
#' the prior is the LLM-assigned 0-1 evidence score. With `prior.method =
#' "pubmed_count"`, CASPEN queries PubMed for each disease + outcome + pathway
#' combination and computes `pubmed.score = log1p(count) / max(log1p(count))`
#' across the submitted pathways. With `prior.method = "pubtator"`, CASPEN
#' retrieves PubMed IDs and checks PubTator3 gene annotations for submitted
#' pathway genes. With `prior.method = "hybrid"`, CASPEN forms a weighted
#' average of available component scores:
#'
#' `prior = sum(component_score * hybrid.weight) / sum(available weights)`.
#'
#' Supported hybrid components are `pubmed_count`, `pubtator`, `llm`,
#' `celltype`, and `evidence`. Components that are unavailable or `NA` are
#' omitted and the remaining weights are re-normalized.
#'
#' @param disease Character disease context, for example `"high-grade serous
#' ovarian cancer"`.
#' @param outcome.context Character outcome context, for example `"platinum
#' resistance"` or `"overall survival"`.
#' @param features Optional user pathway list. In prior-only mode, the LLM
#' scores these existing pathway names.
#' @param celltypes Optional cell type names. These are included in the prompt
#' so the LLM can prioritize celltype-relevant biology.
#' @param mode One of `"curate"` or `"prior"`. `"curate"` returns
#' disease/outcome-relevant concepts or signatures with prior weights;
#' `"prior"` returns prior weights for existing user pathways only.
#' @param n.signatures Maximum number of new literature signatures to request.
#' @param max.genes Maximum number of genes per LLM-curated signature.
#' @param min.genes Minimum number of curated genes that must map to the input
#' data for a new LLM-curated signature to be retained. Default is 5.
#' @param provider LLM provider. Currently `"openai"` or `"custom"`.
#' @param model OpenAI model name. Defaults to `CASPEN_OPENAI_MODEL` or
#' `"gpt-4.1-mini"`.
#' @param api.key OpenAI API key. Defaults to `OPENAI_API_KEY`.
#' @param llm.fn Optional custom function accepting a prompt and returning text.
#' If supplied, `provider` and OpenAI settings are ignored.
#' @param pubmed Logical; if `TRUE`, retrieves PubMed abstracts and provides
#' them as grounding context.
#' @param pubmed.max Maximum number of PubMed records to retrieve.
#' @param email Optional email passed to NCBI E-utilities.
#' @param valid.genes Optional character vector of valid gene symbols. If
#' supplied, LLM genes outside this set are dropped.
#' @param prior.method Character string specifying how prior weights are
#' computed. `"llm"` uses the LLM-assigned score, `"pubmed_count"` uses a
#' normalized PubMed paper-count score for disease + outcome + pathway,
#' `"pubtator"` uses PubTator3 gene-annotation evidence, and `"hybrid"`
#' combines multiple components using `hybrid.weights`.
#' @param hybrid.weights Named numeric vector used when `prior.method =
#' "hybrid"`. Supported names are `"pubmed_count"`, `"pubtator"`, `"llm"`,
#' `"celltype"`, and `"evidence"`. Weights are normalized over available
#' finite components.
#' @param prior.clip Numeric length-two range used to clip prior weights.
#' @param temperature OpenAI sampling temperature.
#'
#' @return A list with `literature.features`, `pathway.prior`, `evidence`,
#' `raw.response`, and `pubmed.records`. In `mode = "prior"`, the `n.genes`
#' column in `evidence` reports the true number of genes in the user-supplied
#' pathway when the LLM pathway name matches the input feature name. The
#' `evidence` table also includes component scores used to form the final
#' prior.
#'
#' @examples
#' \dontrun{
#' # 1) Load example data and pathway definitions.
#' library(CASPEN)
#' data(caspen_example)
#'
#' pathways <- caspen_example$pathways[1:5]
#'
#' # 2) Score user-provided pathways with OpenAI.
#' # Set this once per R session, or store it in ~/.Renviron.
#' Sys.setenv(OPENAI_API_KEY = "sk-proj-your_real_key_here")
#'
#' llm_res <- llm_literature_signatures(
#'   disease = "high-grade serous ovarian cancer",
#'   outcome.context = "platinum resistance",
#'   features = pathways,
#'   mode = "prior",
#'   pubmed = FALSE,
#'   model = "gpt-4o-mini",
#'   temperature = NULL
#' )
#'
#' llm_res$pathway.prior
#' llm_res$evidence
#'
#' # The evidence table separates each score component:
#' # pubmed.count  = number of PubMed hits for disease + outcome + pathway
#' # pubmed.score  = normalized log PubMed count, scaled to 0-1
#' # pubtator.count= PubTator3 gene-annotation evidence count
#' # pubtator.score= normalized PubTator3 evidence, scaled to 0-1
#' # llm.score     = LLM relevance score, scaled to 0-1
#' # evidence.score= score based on PMID/citation support returned by the LLM
#' # prior         = final prior according to prior.method
#' llm_res$evidence[, c("pathway", "n.genes", "llm.score", "prior")]
#'
#' # 2b) Compute priors without any LLM call using normalized PubMed counts.
#' # This queries PubMed for disease + outcome + each pathway name.
#' pubmed_res <- llm_literature_signatures(
#'   disease = "high-grade serous ovarian cancer",
#'   outcome.context = "platinum resistance",
#'   features = pathways,
#'   mode = "prior",
#'   prior.method = "pubmed_count",
#'   email = "your.email@example.com"
#' )
#'
#' pubmed_res$evidence[, c("pathway", "pubmed.count", "pubmed.score", "prior")]
#'
#' # 2c) Compute PubTator3 priors without an LLM call.
#' pubtator_res <- llm_literature_signatures(
#'   disease = "high-grade serous ovarian cancer",
#'   outcome.context = "platinum resistance",
#'   features = pathways,
#'   mode = "prior",
#'   prior.method = "pubtator",
#'   pubmed.max = 20,
#'   email = "your.email@example.com"
#' )
#'
#' pubtator_res$evidence[, c("pathway", "pubtator.count",
#'                           "pubtator.score", "prior")]
#'
#' # 2d) Compute a hybrid prior. Component weights are user-controlled.
#' hybrid_res <- llm_literature_signatures(
#'   disease = "high-grade serous ovarian cancer",
#'   outcome.context = "platinum resistance",
#'   features = pathways,
#'   mode = "prior",
#'   pubmed = FALSE,
#'   model = "gpt-4o-mini",
#'   temperature = NULL,
#'   prior.method = "hybrid",
#'   hybrid.weights = c(pubmed_count = 0.25, pubtator = 0.25,
#'                      llm = 0.4, evidence = 0.1)
#' )
#'
#' hybrid_res$evidence[, c("pathway", "pubmed.score", "pubtator.score",
#'                         "llm.score", "evidence.score", "prior")]
#'
#' # 2e) Curate disease/outcome-specific literature concepts/signatures.
#' # In pathway_select(), curated concepts can be mapped back to user pathways.
#' curated_res <- llm_literature_signatures(
#'   disease = "high-grade serous ovarian cancer",
#'   outcome.context = "platinum resistance",
#'   features = pathways,
#'   mode = "curate",
#'   n.signatures = 5,
#'   max.genes = 40,
#'   prior.method = "hybrid",
#'   hybrid.weights = c(pubmed_count = 0.5, llm = 0.5),
#'   pubmed = TRUE,
#'   pubmed.max = 5,
#'   email = "your.email@example.com",
#'   model = "gpt-4o-mini",
#'   temperature = NULL
#' )
#'
#' names(curated_res$literature.features)
#' curated_res$pathway.prior
#' curated_res$evidence
#'
#' # 3) Use the prior weights during pathway selection.
#' path_res <- pathway_select(
#'   outcome.type = "binary",
#'   outcome = caspen_example$outcome,
#'   data = caspen_example$data_to_pass,
#'   features = pathways,
#'   models.indiv = c("RF", "EN"),
#'   iter = 5,
#'   num.folds = 3,
#'   pathway.prior = llm_res$pathway.prior,
#'   llm.mode = "none"
#' )
#'
#' path_res$data.metric
#' path_res$SP.95
#' path_res$Prior.weight
#' path_res$selection.table
#'
#' # 4) Free local alternative: Ollama.
#' # Install Ollama from https://ollama.com/download, then in Terminal run:
#' # ollama pull qwen2.5:7b
#' # Keep the Ollama app/service running while R calls localhost:11434.
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
#' ollama_res <- llm_literature_signatures(
#'   disease = "high-grade serous ovarian cancer",
#'   outcome.context = "platinum resistance",
#'   features = pathways,
#'   mode = "prior",
#'   pubmed = FALSE,
#'   provider = "custom",
#'   llm.fn = ollama_llm
#' )
#'
#' ollama_res$pathway.prior
#' ollama_res$evidence
#' }
#'
#' @export
llm_literature_signatures <- function(disease, outcome.context,
                                      features = NULL, celltypes = NULL,
                                      mode = c("curate", "prior"),
                                      n.signatures = 10, max.genes = 80,
                                      min.genes = 5,
                                      provider = c("openai", "custom"),
                                      model = Sys.getenv("CASPEN_OPENAI_MODEL", "gpt-4.1-mini"),
                                      api.key = Sys.getenv("OPENAI_API_KEY"),
                                      llm.fn = NULL,
                                      pubmed = TRUE, pubmed.max = 8,
                                      email = NULL,
                                      valid.genes = NULL,
                                      prior.method = c("llm", "pubmed_count", "pubtator", "hybrid"),
                                      hybrid.weights = c(pubmed_count = 0.5, llm = 0.5),
                                      prior.clip = c(0, 1),
                                      temperature = 0) {
  if (identical(mode, "both")) {
    warning("mode = 'both' is deprecated. Use mode = 'curate'; ",
            "'curate' now curates concepts and scores priors.")
    mode <- "curate"
  }
  mode <- match.arg(mode)
  provider <- match.arg(provider)
  prior.method <- match.arg(prior.method)
  min.genes <- max(1L, as.integer(min.genes %||% 5L))
  if (!nzchar(disease %||% "")) stop("disease is required for LLM curation.")
  if (!nzchar(outcome.context %||% "")) stop("outcome.context is required for LLM curation.")
  if (mode == "prior" && (is.null(features) || !length(features))) {
    stop("features is required when mode = 'prior'.")
  }

  feature.names <- if (is.null(features)) character(0) else names(features)
  if (!length(feature.names) && length(features)) {
    feature.names <- paste0("Pathway_", seq_along(features))
  }
  input.gene.counts <- caspen_feature_gene_counts(features, feature.names)
  needs.llm <- mode == "curate" || prior.method %in% c("llm", "hybrid")
  hybrid.weights.normalized <- caspen_normalize_hybrid_weights(hybrid.weights)
  pubmed.weight <- hybrid.weights.normalized[match("pubmed_count",
                                                   names(hybrid.weights.normalized))]
  if (!length(pubmed.weight)) pubmed.weight <- NA_real_
  pubtator.weight <- hybrid.weights.normalized[match("pubtator",
                                                     names(hybrid.weights.normalized))]
  if (!length(pubtator.weight)) pubtator.weight <- NA_real_
  needs.prior.components <- mode %in% c("prior", "curate")
  needs.pubmed.count <- isTRUE(needs.prior.components) &&
    (prior.method == "pubmed_count" ||
    (prior.method == "hybrid" &&
       is.finite(pubmed.weight) && pubmed.weight > 0)
  )
  needs.pubtator <- isTRUE(needs.prior.components) &&
    (prior.method == "pubtator" ||
    (prior.method == "hybrid" &&
       is.finite(pubtator.weight) && pubtator.weight > 0)
  )
  records <- if (isTRUE(pubmed) && isTRUE(needs.llm)) {
    caspen_pubmed_retrieve(
      query = paste(disease, outcome.context, "gene signature pathway"),
      k = pubmed.max,
      email = email
    )
  } else {
    data.frame(pmid = character(0), abstract = character(0), stringsAsFactors = FALSE)
  }

  if (needs.llm) {
    prompt <- caspen_llm_prompt(
      disease = disease,
      outcome.context = outcome.context,
      mode = mode,
      feature.names = feature.names,
      celltypes = celltypes,
      records = records,
      n.signatures = n.signatures,
      max.genes = max.genes,
      min.genes = min.genes
    )
    raw <- if (is.function(llm.fn)) {
      llm.fn(prompt)
    } else {
      if (provider == "custom") stop("llm.fn must be supplied when provider = 'custom'.")
      caspen_openai_response(prompt, model = model, api.key = api.key,
                             temperature = temperature)
    }
    parsed <- caspen_parse_llm_json(raw)
    signatures <- parsed$signatures
  } else {
    raw <- NA_character_
    signatures <- list()
  }
  raw.signatures.all <- caspen_raw_signature_table(signatures, min.genes = 1)
  raw.signatures <- caspen_raw_signature_table(signatures, min.genes = min.genes)

  if (mode == "prior" && prior.method %in% c("pubmed_count", "pubtator")) {
    signatures <- lapply(feature.names, function(nm) {
      list(name = nm, genes = character(0), prior = NA_real_,
           pmids = character(0),
           rationale = paste("Prior computed from", prior.method, "evidence."))
    })
  }

  signature.names <- unique(vapply(signatures, function(sig) {
    clean_signature_name(sig$name)
  }, character(1)))
  signature.names <- signature.names[nzchar(signature.names)]
  pubmed.prior <- if (isTRUE(needs.pubmed.count) && length(signature.names)) {
    caspen_pubmed_prior_scores(
      disease = disease,
      outcome.context = outcome.context,
      pathway.names = signature.names,
      email = email
    )
  } else {
    data.frame(pathway = character(0), pubmed.count = numeric(0),
               pubmed.score = numeric(0), stringsAsFactors = FALSE)
  }
  pubtator.prior <- if (isTRUE(needs.pubtator) && length(signature.names)) {
    caspen_pubtator_prior_scores(
      disease = disease,
      outcome.context = outcome.context,
      features = features,
      pathway.names = signature.names,
      email = email,
      pubmed.max = pubmed.max
    )
  } else {
    data.frame(pathway = character(0), pubtator.count = numeric(0),
               pubtator.score = numeric(0), stringsAsFactors = FALSE)
  }

  literature.features <- list()
  pathway.prior <- numeric(0)
  evidence_rows <- list()
  if (length(signatures)) {
    if (is.null(valid.genes)) valid.genes <- infer_valid_gene_symbols(features)
    valid.upper <- if (is.null(valid.genes)) NULL else clean_gene_symbols(valid.genes)
    existing.clean <- clean_signature_name(feature.names)
    existing.canon <- canonical_signature_name(feature.names)

    for (sig in signatures) {
      nm <- clean_signature_name(sig$name)
      if (!nzchar(nm)) next
      nm.canon <- canonical_signature_name(nm)
      existing.hit <- match(nm.canon, existing.canon)
      if (is.na(existing.hit)) existing.hit <- match(nm, existing.clean)
      nm.out <- if (!is.na(existing.hit)) existing.clean[[existing.hit]] else nm
      llm.genes <- clean_gene_symbols(sig$genes %||% character(0))
      genes <- llm.genes
      if (!is.null(valid.upper)) genes <- genes[genes %in% valid.upper]
      is.existing.pathway <- !is.na(existing.hit)
      is.new.curated.signature <- mode == "curate" && !is.existing.pathway
      if (is.new.curated.signature && length(genes) >= min.genes) {
        literature.features[[nm.out]] <- genes
      }
      llm.score <- suppressWarnings(as.numeric(sig$prior %||% NA_real_))
      pubmed.row <- pubmed.prior[match(nm.out, pubmed.prior$pathway), , drop = FALSE]
      pubmed.count <- if (nrow(pubmed.row)) pubmed.row$pubmed.count else NA_real_
      pubmed.score <- if (nrow(pubmed.row)) pubmed.row$pubmed.score else NA_real_
      pubtator.row <- pubtator.prior[match(nm.out, pubtator.prior$pathway), , drop = FALSE]
      pubtator.count <- if (nrow(pubtator.row)) pubtator.row$pubtator.count else NA_real_
      pubtator.score <- if (nrow(pubtator.row)) pubtator.row$pubtator.score else NA_real_
      celltype.score <- caspen_signature_score(sig, c("celltype_score",
                                                      "celltype_relevance"))
      evidence.score <- caspen_signature_score(sig, c("evidence_score",
                                                      "evidence_quality",
                                                      "citation_strength"))
      if (!is.finite(evidence.score)) {
        evidence.score <- caspen_pmid_evidence_score(sig$pmids %||% character(0))
      }
      prior.out <- if (mode %in% c("prior", "curate")) {
        caspen_combine_prior_components(
          prior.method = prior.method,
          components = c(
            pubmed_count = unname(pubmed.score),
            pubtator = unname(pubtator.score),
            llm = unname(llm.score),
            celltype = unname(celltype.score),
            evidence = unname(evidence.score)
          ),
          hybrid.weights = hybrid.weights.normalized
        )
      } else NA_real_
      if (mode %in% c("prior", "curate")) {
        pathway.prior <- c(pathway.prior, stats::setNames(prior.out, nm.out))
      }
      evidence_rows[[length(evidence_rows) + 1]] <- data.frame(
        pathway = nm.out,
        llm.name = nm,
        duplicate.of.user.pathway = is.existing.pathway,
        prior = prior.out,
        n.genes = length(genes),
        n.llm.genes = length(llm.genes),
        n.mapped.genes = length(genes),
        mapped.genes = paste(genes, collapse = ", "),
        pubmed.count = pubmed.count,
        pubmed.score = pubmed.score,
        pubtator.count = pubtator.count,
        pubtator.score = pubtator.score,
        llm.score = llm.score,
        celltype.score = celltype.score,
        evidence.score = evidence.score,
        prior.method = prior.method,
        pmids = paste(as.character(sig$pmids %||% character(0)), collapse = ";"),
        rationale = as.character(sig$rationale %||% ""),
        stringsAsFactors = FALSE
      )
    }
  }

  evidence <- if (length(evidence_rows)) do.call(rbind, evidence_rows) else data.frame()

  if (mode == "prior" && length(feature.names)) {
    completed <- caspen_complete_prior_output(
      evidence = evidence,
      pathway.prior = pathway.prior,
      feature.names = feature.names,
      input.gene.counts = input.gene.counts,
      prior.method = prior.method
    )
    evidence <- completed$evidence
    pathway.prior <- completed$pathway.prior
  }

  if (length(pathway.prior)) {
    prior.names <- names(pathway.prior)
    pathway.prior <- pmin(max(prior.clip), pmax(min(prior.clip), pathway.prior))
    names(pathway.prior) <- prior.names
  }

  list(
    literature.features = literature.features,
    pathway.prior = pathway.prior,
    evidence = evidence,
    raw.signatures.all = raw.signatures.all,
    raw.signatures = raw.signatures,
    raw.response = raw,
    pubmed.records = records
  )
}

#' @keywords internal
#' @noRd
caspen_pubmed_retrieve <- function(query, k = 8, email = NULL) {
  base <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  search_url <- paste0(base, "esearch.fcgi")
  req <- httr2::request(search_url) |>
    httr2::req_url_query(
      db = "pubmed", term = query, retmax = k, retmode = "json",
      email = email %||% NULL
    )
  ids <- tryCatch({
    res <- httr2::req_perform(req)
    txt <- httr2::resp_body_string(res)
    jsonlite::fromJSON(txt)$esearchresult$idlist
  }, error = function(e) character(0))
  ids <- ids[nzchar(ids)]
  if (!length(ids)) {
    return(data.frame(pmid = character(0), abstract = character(0), stringsAsFactors = FALSE))
  }

  fetch_url <- paste0(base, "efetch.fcgi")
  req2 <- httr2::request(fetch_url) |>
    httr2::req_url_query(
      db = "pubmed", id = paste(ids, collapse = ","), rettype = "abstract",
      retmode = "text", email = email %||% NULL
    )
  txt <- tryCatch(httr2::resp_body_string(httr2::req_perform(req2)),
                  error = function(e) "")
  chunks <- strsplit(txt, "\n\n\n+", perl = TRUE)[[1]]
  chunks <- chunks[nzchar(trimws(chunks))]
  n <- min(length(ids), length(chunks))
  if (!n) {
    chunks <- rep("", length(ids))
    n <- length(ids)
  }
  data.frame(
    pmid = ids[seq_len(n)],
    abstract = chunks[seq_len(n)],
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
caspen_pubmed_count <- function(query, email = NULL) {
  base <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  req <- httr2::request(paste0(base, "esearch.fcgi")) |>
    httr2::req_url_query(
      db = "pubmed", term = query, retmax = 0, retmode = "json",
      email = email %||% NULL
    )
  tryCatch({
    res <- httr2::req_perform(req)
    txt <- httr2::resp_body_string(res)
    out <- jsonlite::fromJSON(txt)
    as.numeric(out$esearchresult$count %||% 0)
  }, error = function(e) NA_real_)
}

#' @keywords internal
#' @noRd
caspen_pubmed_search_ids <- function(query, k = 20, email = NULL) {
  base <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  req <- httr2::request(paste0(base, "esearch.fcgi")) |>
    httr2::req_url_query(
      db = "pubmed", term = query, retmax = k, retmode = "json",
      email = email %||% NULL
    )
  ids <- tryCatch({
    res <- httr2::req_perform(req)
    txt <- httr2::resp_body_string(res)
    jsonlite::fromJSON(txt)$esearchresult$idlist
  }, error = function(e) character(0))
  ids <- unique(as.character(ids))
  ids[nzchar(ids)]
}

#' @keywords internal
#' @noRd
caspen_pubmed_prior_scores <- function(disease, outcome.context, pathway.names,
                                       email = NULL) {
  pathway.names <- unique(as.character(pathway.names))
  pathway.names <- pathway.names[nzchar(pathway.names)]
  if (!length(pathway.names)) {
    return(data.frame(pathway = character(0), pubmed.count = numeric(0),
                      pubmed.score = numeric(0), stringsAsFactors = FALSE))
  }
  counts <- vapply(pathway.names, function(pathway) {
    pathway.term <- gsub("_+", " ", pathway)
    query <- paste(disease, outcome.context, pathway.term)
    caspen_pubmed_count(query, email = email)
  }, numeric(1))
  scores <- caspen_normalize_counts(counts)
  data.frame(
    pathway = pathway.names,
    pubmed.count = as.numeric(counts),
    pubmed.score = as.numeric(scores),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
caspen_pubtator_export <- function(pmids) {
  pmids <- unique(as.character(pmids %||% character(0)))
  pmids <- pmids[nzchar(pmids)]
  if (!length(pmids)) return(NULL)
  url <- "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/biocjson"
  req <- httr2::request(url) |>
    httr2::req_url_query(pmids = paste(pmids, collapse = ","))
  tryCatch({
    res <- httr2::req_perform(req)
    jsonlite::fromJSON(httr2::resp_body_string(res), simplifyVector = FALSE)
  }, error = function(e) NULL)
}

#' @keywords internal
#' @noRd
caspen_pubtator_annotations <- function(x, pmid = NA_character_) {
  rows <- list()
  walk <- function(node, current.pmid = pmid) {
    if (is.null(node)) return(NULL)
    if (is.list(node)) {
      if (!is.null(node$id) && length(node$id) == 1 && nzchar(as.character(node$id))) {
        current.pmid <- as.character(node$id)
      }
      if (!is.null(node$infons) && !is.null(node$text)) {
        rows[[length(rows) + 1]] <<- data.frame(
          pmid = current.pmid,
          type = as.character(node$infons$type %||% node$infons$Type %||% ""),
          text = as.character(node$text %||% ""),
          identifier = as.character(node$infons$identifier %||%
                                      node$infons$Identifier %||% ""),
          stringsAsFactors = FALSE
        )
      }
      invisible(lapply(node, walk, current.pmid = current.pmid))
    }
    NULL
  }
  walk(x)
  if (!length(rows)) {
    return(data.frame(pmid = character(0), type = character(0),
                      text = character(0), identifier = character(0),
                      stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

#' @keywords internal
#' @noRd
caspen_pubtator_query_count <- function(disease, outcome.context, feature.name,
                                        genes = NULL, email = NULL,
                                        pubmed.max = 20) {
  feature.term <- gsub("_+", " ", feature.name)
  query <- paste(disease, outcome.context, feature.term)
  pmids <- caspen_pubmed_search_ids(query, k = pubmed.max, email = email)
  if (!length(pmids)) return(0)
  ann <- caspen_pubtator_annotations(caspen_pubtator_export(pmids))
  if (!nrow(ann)) return(0)
  ann <- ann[tolower(ann$type) %in% c("gene", "genes", "protein"), , drop = FALSE]
  if (!nrow(ann)) return(0)
  genes <- unique(toupper(as.character(genes %||% character(0))))
  genes <- genes[nzchar(genes)]
  if (length(genes)) {
    hit <- toupper(ann$text) %in% genes
    hit <- hit | toupper(ann$identifier) %in% genes
    ann <- ann[hit, , drop = FALSE]
    if (!nrow(ann)) return(0)
    return(length(unique(ann$pmid)) + length(unique(toupper(ann$text))))
  }
  length(unique(ann$pmid)) + nrow(ann)
}

#' @keywords internal
#' @noRd
caspen_pubtator_prior_scores <- function(disease, outcome.context, features,
                                         pathway.names, email = NULL,
                                         pubmed.max = 20) {
  pathway.names <- unique(as.character(pathway.names))
  pathway.names <- pathway.names[nzchar(pathway.names)]
  if (!length(pathway.names)) {
    return(data.frame(pathway = character(0), pubtator.count = numeric(0),
                      pubtator.score = numeric(0), stringsAsFactors = FALSE))
  }
  feature.genes <- caspen_feature_gene_map(features, pathway.names)
  counts <- vapply(pathway.names, function(pathway) {
    caspen_pubtator_query_count(
      disease = disease,
      outcome.context = outcome.context,
      feature.name = pathway,
      genes = feature.genes[[pathway]] %||% character(0),
      email = email,
      pubmed.max = pubmed.max
    )
  }, numeric(1))
  scores <- caspen_normalize_counts(counts)
  data.frame(
    pathway = pathway.names,
    pubtator.count = as.numeric(counts),
    pubtator.score = as.numeric(scores),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
caspen_normalize_counts <- function(counts) {
  counts <- as.numeric(counts)
  out <- rep(NA_real_, length(counts))
  ok <- is.finite(counts) & counts >= 0
  if (!any(ok)) return(out)
  scaled <- log1p(counts[ok])
  denom <- max(scaled, na.rm = TRUE)
  out[ok] <- if (is.finite(denom) && denom > 0) scaled / denom else 0
  out
}

#' @keywords internal
#' @noRd
caspen_llm_prompt <- function(disease, outcome.context, mode, feature.names,
                              celltypes, records, n.signatures, max.genes,
                              min.genes = 5) {
  literature <- if (nrow(records)) {
    paste(sprintf("PMID %s: %s", records$pmid, substr(records$abstract, 1, 1400)),
          collapse = "\n\n")
  } else {
    "No retrieved abstracts were available. Use general biomedical knowledge, but mark rationale as uncited."
  }
  existing <- if (length(feature.names)) paste(feature.names, collapse = ", ") else "none"
  ct <- if (length(celltypes)) paste(celltypes, collapse = ", ") else "not supplied"
  task <- switch(
    mode,
    curate = paste0("Propose up to ", n.signatures,
                    " disease/outcome-relevant biological concepts or gene signatures and assign prior weights to each."),
    prior = "Assign prior weights to the existing pathway names only. Do not curate new signatures."
  )
  paste(
    "You are helping CASPEN, a celltype-aware pathway-guided ensemble prediction method.",
    "Use the literature context when possible. Prefer human gene symbols.",
    task,
    paste("Disease:", disease),
    paste("Outcome:", outcome.context),
    paste("Cell types of interest:", ct),
    paste("Existing user pathway names:", existing),
    paste("Minimum genes per new curated signature:", min.genes),
    paste("Maximum genes per new signature:", max.genes),
    "Return strict JSON only, with this schema:",
    '{"signatures":[{"name":"AKT_activation","genes":["GENE1","GENE2"],"prior":0.2,"celltype_score":0.2,"evidence_score":0.2,"pmids":["123"],"rationale":"brief evidence"}]}',
    "Use descriptive biological names such as AKT_activation, DNA_repair, TGF_beta_stromal_signaling, or macrophage_inflammation; do not use generic names such as signature_1.",
    "For new curated signatures, merge highly related biological concepts when needed so each returned signature has at least the minimum gene count.",
    "Prior must be numeric from 0 to 1 where 1 is strongest direct evidence.",
    "If possible, celltype_score should quantify cell-type-specific relevance from 0 to 1, and evidence_score should quantify literature/citation strength from 0 to 1.",
    "For prior-only mode, include existing pathway names with empty genes if gene membership is not being curated.",
    "Do not include markdown.",
    "Literature context:",
    literature,
    sep = "\n"
  )
}

#' @keywords internal
#' @noRd
caspen_openai_response <- function(prompt, model, api.key, temperature = 0) {
  if (!nzchar(api.key %||% "")) {
    stop("OPENAI_API_KEY is not set. Supply api.key or llm.fn, or set llm.mode = 'none'.")
  }
  body <- list(
    model = model,
    input = list(list(
      role = "user",
      content = list(list(type = "input_text", text = prompt))
    ))
  )
  if (!is.null(temperature)) body$temperature <- temperature
  req <- httr2::request("https://api.openai.com/v1/responses") |>
    httr2::req_headers(Authorization = paste("Bearer", api.key)) |>
    httr2::req_body_json(body)
  res <- tryCatch(
    httr2::req_perform(req),
    error = function(e) {
      msg <- conditionMessage(e)
      if (!is.null(e$response)) {
        body_txt <- tryCatch(httr2::resp_body_string(e$response),
                             error = function(err) "")
        if (nzchar(body_txt)) {
          msg <- paste0(msg, "\nOpenAI response body:\n", body_txt)
        }
      }
      stop(msg, call. = FALSE)
    }
  )
  parsed <- jsonlite::fromJSON(httr2::resp_body_string(res), simplifyVector = FALSE)
  txt <- parsed$output_text %||% NULL
  if (!is.null(txt) && nzchar(txt)) return(txt)
  parts <- unlist(lapply(parsed$output %||% list(), function(item) {
    unlist(lapply(item$content %||% list(), function(z) z$text %||% z$annotations %||% ""))
  }), use.names = FALSE)
  paste(parts[nzchar(parts)], collapse = "\n")
}

#' @keywords internal
#' @noRd
caspen_parse_llm_json <- function(raw) {
  raw <- trimws(as.character(raw %||% ""))
  raw <- sub("^```json\\s*", "", raw)
  raw <- sub("^```\\s*", "", raw)
  raw <- sub("\\s*```$", "", raw)
  parsed <- tryCatch(jsonlite::fromJSON(raw, simplifyVector = FALSE),
                     error = function(e) NULL)
  if (is.null(parsed)) {
    start <- regexpr("\\[", raw)[1]
    end <- max(gregexpr("\\]", raw)[[1]])
    if (is.finite(start) && start > 0 && is.finite(end) && end > start) {
      parsed <- tryCatch(jsonlite::fromJSON(substr(raw, start, end),
                                            simplifyVector = FALSE),
                         error = function(e) NULL)
    }
  }
  if (is.null(parsed)) {
    start <- regexpr("\\{", raw)[1]
    end <- max(gregexpr("\\}", raw)[[1]])
    if (is.finite(start) && start > 0 && is.finite(end) && end > start) {
      parsed <- tryCatch(jsonlite::fromJSON(substr(raw, start, end),
                                            simplifyVector = FALSE),
                         error = function(e) NULL)
    }
  }
  if (is.null(parsed)) stop("Could not parse LLM response as JSON.")
  parsed.names <- names(parsed)
  array.like <- is.null(parsed.names) ||
    !length(parsed.names) ||
    all(!nzchar(parsed.names)) ||
    identical(parsed.names, as.character(seq_along(parsed)))
  if (is.list(parsed) && isTRUE(array.like) &&
      all(vapply(parsed, is.list, logical(1)))) {
    parsed <- list(signatures = parsed)
  }
  if (is.null(parsed$signatures)) parsed$signatures <- list()
  parsed
}

#' @keywords internal
#' @noRd
caspen_raw_signature_table <- function(signatures, min.genes = 1) {
  if (!length(signatures)) {
    return(data.frame(
      signature.name = character(0),
      genes = character(0),
      n.genes = integer(0),
      llm.weight = numeric(0),
      celltype.score = numeric(0),
      evidence.score = numeric(0),
      pmids = character(0),
      rationale = character(0),
      stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(signatures, function(sig) {
    genes <- clean_gene_symbols(sig$genes %||% character(0))
    data.frame(
      signature.name = caspen_display_signature_name(sig$name),
      genes = paste(genes, collapse = ", "),
      n.genes = length(genes),
      llm.weight = suppressWarnings(as.numeric(sig$prior %||% NA_real_)),
      celltype.score = caspen_signature_score(sig, c("celltype_score",
                                                     "celltype_relevance")),
      evidence.score = caspen_signature_score(sig, c("evidence_score",
                                                     "evidence_quality",
                                                     "citation_strength")),
      pmids = paste(as.character(sig$pmids %||% character(0)), collapse = ";"),
      rationale = as.character(sig$rationale %||% ""),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  min.genes <- max(1L, as.integer(min.genes %||% 1L))
  out[out$n.genes >= min.genes, , drop = FALSE]
}

#' @keywords internal
#' @noRd
clean_signature_name <- function(x) {
  x <- gsub("[^[:alnum:]_ .-]+", "", as.character(x %||% ""))
  x <- gsub("[[:space:].-]+", "_", trimws(x))
  x
}

#' @keywords internal
#' @noRd
caspen_display_signature_name <- function(x) {
  x <- clean_signature_name(x)
  x <- gsub("^(signature|signatures|pathway|pathways|geneset|genesets)_+", "",
            x, ignore.case = TRUE)
  x <- gsub("_+(signature|signatures|pathway|pathways|geneset|genesets)$", "",
            x, ignore.case = TRUE)
  parts <- strsplit(x, "_", fixed = TRUE)[[1]]
  acronyms <- c("akt", "pi3k", "mtor", "egfr", "erbb", "mapk", "erk", "jak",
                "stat", "dna", "ddr", "emt", "tgf", "ifn", "nfkb", "p53",
                "hif", "vegf", "caf", "nk")
  parts <- vapply(parts, function(part) {
    low <- tolower(part)
    if (low %in% acronyms) toupper(part) else part
  }, character(1))
  paste(parts[nzchar(parts)], collapse = "_")
}

#' @keywords internal
#' @noRd
canonical_signature_name <- function(x) {
  x <- tolower(clean_signature_name(x))
  x <- gsub("^(signature|signatures|pathway|pathways|geneset|genesets)_+", "", x)
  x <- gsub("_+(signature|signatures|pathway|pathways|geneset|genesets)$", "", x)
  x
}

#' @keywords internal
#' @noRd
infer_valid_gene_symbols <- function(features) {
  if (is.null(features) || !length(features)) return(NULL)
  genes <- unique(unlist(lapply(features, extract_pathway_genes), use.names = FALSE))
  genes <- genes[nzchar(genes)]
  if (!length(genes)) NULL else genes
}

#' @keywords internal
#' @noRd
clean_gene_symbols <- function(x) {
  x <- trimws(as.character(unlist(x %||% character(0), use.names = FALSE)))
  paren <- regmatches(x, regexpr("\\([A-Za-z0-9_.-]+\\)", x))
  has.paren <- nzchar(paren)
  x[has.paren] <- gsub("^\\(|\\)$", "", paren[has.paren])
  x <- toupper(trimws(x))
  x <- x[!is.na(x) & nzchar(x)]
  x <- x[!x %in% c("NA", "NULL", "NAN")]
  aliases <- c(
    LC3B = "MAP1LC3B",
    VEGF = "VEGFA",
    P53 = "TP53",
    MTOR = "MTOR",
    CYCLINB1 = "CCNB1",
    BECLIN1 = "BECN1",
    SNAIL = "SNAI1"
  )
  hits <- match(x, names(aliases))
  x[!is.na(hits)] <- unname(aliases[hits[!is.na(hits)]])
  x <- x[grepl("^[A-Z0-9][A-Z0-9_.-]*$", x)]
  unique(x)
}

#' @keywords internal
#' @noRd
caspen_signature_score <- function(sig, fields) {
  for (field in fields) {
    value <- sig[[field]]
    if (!is.null(value)) {
      score <- suppressWarnings(as.numeric(value))
      if (is.finite(score)) return(max(0, min(1, score)))
    }
  }
  NA_real_
}

#' @keywords internal
#' @noRd
caspen_pmid_evidence_score <- function(pmids) {
  pmids <- unique(as.character(pmids %||% character(0)))
  pmids <- pmids[nzchar(pmids)]
  if (!length(pmids)) return(NA_real_)
  min(1, length(pmids) / 5)
}

#' @keywords internal
#' @noRd
caspen_complete_prior_output <- function(evidence, pathway.prior, feature.names,
                                         input.gene.counts = numeric(0),
                                         prior.method = NA_character_) {
  feature.names <- unique(as.character(feature.names %||% character(0)))
  feature.names <- feature.names[nzchar(feature.names)]
  if (!length(feature.names)) {
    return(list(evidence = evidence, pathway.prior = pathway.prior))
  }

  clean.features <- clean_signature_name(feature.names)
  prior.out <- stats::setNames(rep(NA_real_, length(feature.names)), feature.names)
  if (length(pathway.prior)) {
    prior.names <- names(pathway.prior)
    clean.prior.names <- clean_signature_name(prior.names)
    for (i in seq_along(feature.names)) {
      hit <- match(feature.names[i], prior.names)
      if (is.na(hit)) hit <- match(clean.features[i], clean.prior.names)
      if (!is.na(hit)) prior.out[i] <- as.numeric(pathway.prior[[hit]])
    }

    extra <- pathway.prior[!clean.prior.names %in% clean.features &
                             !prior.names %in% feature.names]
    pathway.prior <- c(prior.out, extra)
  } else {
    pathway.prior <- prior.out
  }

  if (is.null(evidence) || !nrow(evidence)) {
    evidence <- data.frame(
      pathway = character(0), prior = numeric(0), n.genes = integer(0),
      llm.name = character(0), duplicate.of.user.pathway = logical(0),
      n.llm.genes = integer(0), n.mapped.genes = integer(0),
      mapped.genes = character(0),
      pubmed.count = numeric(0), pubmed.score = numeric(0),
      pubtator.count = numeric(0), pubtator.score = numeric(0),
      llm.score = numeric(0), celltype.score = numeric(0),
      evidence.score = numeric(0), prior.method = character(0),
      pmids = character(0), rationale = character(0),
      stringsAsFactors = FALSE
    )
  }

  clean.evidence <- clean_signature_name(evidence$pathway)
  missing <- feature.names[!feature.names %in% evidence$pathway &
                             !clean.features %in% clean.evidence]
  if (length(missing)) {
    missing.rows <- lapply(missing, function(nm) {
      data.frame(
        pathway = nm,
        llm.name = "",
        duplicate.of.user.pathway = TRUE,
        prior = unname(pathway.prior[[nm]] %||% NA_real_),
        n.genes = caspen_evidence_gene_count(
          nm, character(0), input.gene.counts, use.input = TRUE
        ),
        n.llm.genes = NA_integer_,
        n.mapped.genes = caspen_evidence_gene_count(
          nm, character(0), input.gene.counts, use.input = TRUE
        ),
        mapped.genes = "",
        pubmed.count = NA_real_,
        pubmed.score = NA_real_,
        pubtator.count = NA_real_,
        pubtator.score = NA_real_,
        llm.score = NA_real_,
        celltype.score = NA_real_,
        evidence.score = NA_real_,
        prior.method = as.character(prior.method %||% NA_character_),
        pmids = "",
        rationale = "No prior evidence was returned for this selected pathway.",
        stringsAsFactors = FALSE
      )
    })
    evidence <- rbind(evidence, do.call(rbind, missing.rows))
  }

  list(evidence = evidence, pathway.prior = pathway.prior)
}

#' @keywords internal
#' @noRd
caspen_combine_prior_components <- function(prior.method, components,
                                            hybrid.weights = c(pubmed_count = 0.5,
                                                               llm = 0.5)) {
  component.names <- names(components)
  components <- as.numeric(components)
  names(components) <- component.names %||% character(length(components))
  components <- pmin(1, pmax(0, components))
  names(components) <- component.names %||% character(length(components))

  get_component <- function(name) {
    idx <- match(name, names(components))
    if (is.na(idx)) NA_real_ else components[[idx]]
  }

  if (prior.method == "llm") return(get_component("llm"))
  if (prior.method == "pubmed_count") {
    return(get_component("pubmed_count"))
  }
  if (prior.method == "pubtator") {
    return(get_component("pubtator"))
  }

  weights <- caspen_normalize_hybrid_weights(hybrid.weights)
  common <- intersect(names(weights), names(components))
  common <- common[is.finite(components[common]) & is.finite(weights[common]) &
                     weights[common] > 0]
  if (!length(common)) return(NA_real_)
  sum(components[common] * weights[common]) / sum(weights[common])
}

#' @keywords internal
#' @noRd
caspen_normalize_hybrid_weights <- function(weights) {
  if (is.null(weights) || !length(weights)) {
    weights <- c(pubmed_count = 0.5, llm = 0.5)
  }
  nms <- names(weights)
  if (is.null(nms) || any(!nzchar(nms))) {
    stop("hybrid.weights must be a named numeric vector.")
  }
  weights <- as.numeric(weights)
  nms <- tolower(gsub("[^[:alnum:]]+", "_", nms))
  nms[nms %in% c("pubmed", "pubmed_score", "pubmedcount")] <- "pubmed_count"
  nms[nms %in% c("pubtator3", "pubtator_score", "pubtator_count")] <- "pubtator"
  nms[nms %in% c("llm_score", "literature", "literature_score")] <- "llm"
  nms[nms %in% c("celltype_score", "cell_type", "cell_type_score")] <- "celltype"
  nms[nms %in% c("evidence_score", "evidence_quality",
                 "citation", "citation_strength")] <- "evidence"
  names(weights) <- nms
  weights <- stats::aggregate(weights, by = list(names(weights)), FUN = sum)
  stats::setNames(as.numeric(weights$x), weights$Group.1)
}

#' @keywords internal
#' @noRd
caspen_feature_gene_counts <- function(features, feature.names) {
  if (is.null(features) || !length(features)) return(numeric(0))
  if (!length(feature.names)) feature.names <- paste0("Pathway_", seq_along(features))
  counts <- vapply(features, function(x) {
    length(extract_pathway_genes(x))
  }, numeric(1))
  names(counts) <- feature.names
  clean.names <- clean_signature_name(feature.names)
  extra <- counts
  names(extra) <- clean.names
  c(counts, extra[!names(extra) %in% names(counts)])
}

#' @keywords internal
#' @noRd
caspen_feature_gene_map <- function(features, feature.names) {
  out <- stats::setNames(vector("list", length(feature.names)), feature.names)
  if (is.null(features) || !length(features)) return(out)
  original.names <- names(features)
  if (is.null(original.names) || !length(original.names)) {
    original.names <- paste0("Pathway_", seq_along(features))
  }
  clean.names <- clean_signature_name(original.names)
  for (i in seq_along(features)) {
    genes <- unique(toupper(extract_pathway_genes(features[[i]])))
    genes <- genes[nzchar(genes)]
    hits <- feature.names %in% c(original.names[i], clean.names[i])
    if (any(hits)) {
      out[hits] <- list(genes)
    }
  }
  out
}

#' @keywords internal
#' @noRd
caspen_evidence_gene_count <- function(pathway.name, genes, input.gene.counts,
                                       use.input = FALSE) {
  if (isTRUE(use.input) && length(input.gene.counts)) {
    hit <- if (pathway.name %in% names(input.gene.counts)) {
      input.gene.counts[[pathway.name]]
    } else NULL
    if (!is.null(hit) && is.finite(hit)) return(as.integer(hit))
    clean.name <- clean_signature_name(pathway.name)
    clean.hit <- if (clean.name %in% names(input.gene.counts)) {
      input.gene.counts[[clean.name]]
    } else NULL
    if (!is.null(clean.hit) && is.finite(clean.hit)) return(as.integer(clean.hit))
  }
  length(genes)
}
