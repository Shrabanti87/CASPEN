library(shiny)
library(DT)
library(CASPEN)

`%||%` <- function(a, b) if (!is.null(a)) a else b

data(caspen_example, package = "CASPEN")
demo_data <- caspen_example$data_to_pass
demo_outcome <- caspen_example$outcome
demo_pathways <- caspen_example$pathways
default_test_data_file <- system.file("extdata", "protein_data_Test.csv", package = "CASPEN")
default_test_outcome_file <- system.file("extdata", "outcome_binary_Test.csv", package = "CASPEN")

available_models <- c("XG", "RF", "EN", "GB", "KNN", "SVM", "NB", "DCT", "NN", "ADB", "MB")

code_block <- function(...) paste(..., collapse = "\n")

read_demo_csv <- function(path, reference_cols = NULL, outcome_samples = NULL) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(x) || !ncol(x)) return(as.data.frame(x, check.names = FALSE))

  first_col <- as.character(x[[1]])
  row_feature_overlap <- if (!is.null(reference_cols)) {
    sum(first_col %in% reference_cols, na.rm = TRUE)
  } else 0
  col_feature_overlap <- if (!is.null(reference_cols)) {
    sum(colnames(x) %in% reference_cols, na.rm = TRUE)
  } else 0

  sample_cols <- colnames(x)[-1]
  sample_overlap <- if (!is.null(outcome_samples)) {
    sum(sample_cols %in% outcome_samples, na.rm = TRUE)
  } else 0

  feature_by_sample <- ncol(x) > 1 &&
    row_feature_overlap > col_feature_overlap &&
    (is.null(outcome_samples) || sample_overlap > 0)

  if (feature_by_sample) {
    feature_names <- make.unique(first_col)
    mat <- as.data.frame(t(as.matrix(x[-1])), check.names = FALSE)
    colnames(mat) <- feature_names
    rownames(mat) <- sample_cols
    return(mat)
  }

  if (ncol(x) > 1 && !is.numeric(x[[1]])) {
    rownames(x) <- make.unique(first_col)
    x[[1]] <- NULL
  }
  as.data.frame(x, check.names = FALSE)
}

read_demo_outcome <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(x) == 1) return(as.numeric(x[[1]]))
  outcome_col <- grep("outcome|status|response|class|label", names(x),
                      ignore.case = TRUE, value = TRUE)[1]
  if (is.na(outcome_col)) outcome_col <- names(x)[ncol(x)]
  as.numeric(x[[outcome_col]])
}

read_demo_outcome_table <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  outcome_col <- grep("outcome|status|response|class|label", names(x),
                      ignore.case = TRUE, value = TRUE)[1]
  if (is.na(outcome_col)) outcome_col <- names(x)[ncol(x)]
  sample_col <- grep("sample|sample_id|id", names(x),
                     ignore.case = TRUE, value = TRUE)[1]
  if (is.na(sample_col) && ncol(x) > 1) sample_col <- names(x)[1]
  list(
    outcome = as.numeric(x[[outcome_col]]),
    sample = if (!is.na(sample_col)) as.character(x[[sample_col]]) else NULL
  )
}

default_test_available <- nzchar(default_test_data_file) && nzchar(default_test_outcome_file)

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body { background: #f7f8fa; }
      .title-panel { padding: 18px 0 8px; }
      .caspen-card {
        background: white;
        border: 1px solid #e5e7eb;
        border-radius: 8px;
        padding: 16px;
        margin-bottom: 14px;
      }
      .metric-note { color: #4b5563; font-size: 13px; }
      pre { background: #111827; color: #f9fafb; border: 0; }
    "))
  ),
  div(
    class = "title-panel",
    h2("CASPEN Demo"),
    p("Celltype-aware stacked pathway-guided ensemble prediction on a bundled toy dataset.")
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Demo Controls"),
      selectInput(
        "pathways",
        "Pathways",
        choices = names(demo_pathways),
        selected = names(demo_pathways)[seq_len(min(5, length(demo_pathways)))],
        multiple = TRUE
      ),
      selectInput(
        "gene_pathway",
        "Gene Selection Pathway",
        choices = names(demo_pathways)[seq_len(min(5, length(demo_pathways)))],
        selected = names(demo_pathways)[1],
        multiple = FALSE
      ),
      selectInput(
        "models",
        "Base Models",
        choices = available_models,
        selected = c("DCT"),
        multiple = TRUE
      ),
      numericInput("iter", "Iterations", value = 1, min = 1, max = 5, step = 1),
      numericInput("folds", "CV Folds", value = 2, min = 2, max = 5, step = 1),
      checkboxInput("ensemble", "Fit Ensemble Models", value = FALSE),
      checkboxInput("parallel_pathways", "Parallelize Pathways", value = FALSE),
      checkboxInput("parallel_iter", "Parallelize Iterations", value = FALSE),
      numericInput("workers", "Workers", value = min(2, future::availableCores()),
                   min = 1, max = max(1, future::availableCores()), step = 1),
      hr(),
      h4("Literature Priors"),
      selectInput(
        "llm_mode",
        "Pathway LLM Mode",
        choices = c("none", "prior", "curate", "both"),
        selected = "none"
      ),
      selectInput(
        "gene_prior_method",
        "Gene Prior Method",
        choices = c("none", "pubmed_count", "pubtator", "llm", "hybrid"),
        selected = "none"
      ),
      selectInput(
        "llm_prior_method",
        "Pathway Prior Method",
        choices = c("pubmed_count", "pubtator", "llm", "hybrid"),
        selected = "pubmed_count"
      ),
      textInput("disease", "Disease", value = "high-grade serous ovarian cancer"),
      textInput("outcome_context", "Outcome Context", value = "platinum resistance"),
      selectInput(
        "llm_provider_ui",
        "LLM Provider",
        choices = c("OpenAI" = "openai", "Ollama/local" = "ollama"),
        selected = "ollama"
      ),
      textInput("llm_model", "LLM Model", value = "gpt-4o-mini"),
      textInput("ollama_model", "Ollama Model", value = "qwen2.5:7b"),
      textInput("ollama_url", "Ollama URL", value = "http://localhost:11434"),
      checkboxInput("llm_pubmed", "Use PubMed Abstract Context", value = FALSE),
      numericInput("llm_pubmed_max", "PubMed Records", value = 5, min = 1, max = 20, step = 1),
      numericInput("llm_min_genes", "Min Curated Genes", value = 5, min = 1, max = 50, step = 1),
      textInput("llm_email", "NCBI Email", value = ""),
      actionButton("run_llm_demo", "Run LLM Signature Demo"),
      fileInput("test_data_file", "Replace Test Data CSV", accept = c(".csv")),
      fileInput("test_outcome_file", "Replace Test Outcome CSV", accept = c(".csv")),
      actionButton("run_pathway", "Run Pathway Selection", class = "btn-primary"),
      actionButton("run_gene", "Run Gene Selection"),
      actionButton("run_train", "Run Train CV"),
      actionButton("run_test", "Run Test Performance"),
      hr(),
      p(class = "metric-note",
        "Defaults are intentionally tiny for live demos. Increase iterations in R for real analyses.")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel(
          "Overview",
          div(
            class = "caspen-card",
            h4("Included Example"),
            verbatimTextOutput("data_summary")
          ),
          div(
            class = "caspen-card",
            h4("Workflow"),
            tags$ol(
              tags$li("Select a few pathways and base learners."),
              tags$li("Run pathway_select() to rank pathway models by AUC and SP95."),
              tags$li("Optionally run gene_select() inside the selected pathways."),
              tags$li("Run Train_perform_CV() on the selected pathway features."),
              tags$li("Run Test_perform() on the bundled independent test cohort, or replace it with uploaded CSV files.")
            )
          )
        ),
        tabPanel(
          "Pathway Selection",
          div(class = "caspen-card", h4("AUC / Primary Metric"), DTOutput("path_auc")),
          div(class = "caspen-card", h4("Sensitivity At 95% Specificity"), DTOutput("path_sp95")),
          div(class = "caspen-card", h4("Pathway Scatter"), plotOutput("path_plot", height = 560))
        ),
        tabPanel(
          "Gene Selection",
          div(class = "caspen-card", h4("Selected Genes"), DTOutput("gene_selected")),
          div(class = "caspen-card", h4("Gene Ranking"), DTOutput("gene_ranking")),
          div(class = "caspen-card", h4("Gene Stability Plot"), plotOutput("gene_plot", height = 560))
        ),
        tabPanel(
          "Literature Priors",
          div(class = "caspen-card", h4("Direct llm_literature_signatures() Evidence"), DTOutput("llm_demo_evidence")),
          div(class = "caspen-card", h4("Direct llm_literature_signatures() Priors"), DTOutput("llm_demo_prior")),
          div(class = "caspen-card", h4("Direct llm_literature_signatures() Curated Signatures"), DTOutput("llm_demo_signatures")),
          div(class = "caspen-card", h4("Pathway LLM / Literature Evidence"), DTOutput("llm_evidence")),
          div(class = "caspen-card", h4("Pathway Prior Weights"), DTOutput("llm_prior")),
          div(class = "caspen-card", h4("Gene Prior Columns"), DTOutput("gene_prior_table"))
        ),
        tabPanel(
          "Train CV",
          div(class = "caspen-card", h4("Cross-Validated Training Performance"), DTOutput("train_table"))
        ),
        tabPanel(
          "Test Performance",
          div(class = "caspen-card", h4("External Test Data"), verbatimTextOutput("test_summary")),
          div(class = "caspen-card", h4("Independent Test Performance"), DTOutput("test_table"))
        ),
        tabPanel(
          "Generated R Code",
          div(class = "caspen-card", h4("LLM Signature Demo Code"), verbatimTextOutput("llm_demo_code")),
          div(class = "caspen-card", h4("Pathway Selection Code"), verbatimTextOutput("path_code")),
          div(class = "caspen-card", h4("Gene Selection Code"), verbatimTextOutput("gene_code")),
          div(class = "caspen-card", h4("Train CV Code"), verbatimTextOutput("train_code")),
          div(class = "caspen-card", h4("Test Performance Code"), verbatimTextOutput("test_code"))
        )
      )
    )
  )
)

server <- function(input, output, session) {
  selected_features <- reactive({
    req(input$pathways)
    demo_pathways[input$pathways]
  })

  observeEvent(input$pathways, {
    choices <- input$pathways
    if (!length(choices)) choices <- names(demo_pathways)[1]
    selected <- input$gene_pathway
    if (is.null(selected) || !selected %in% choices) selected <- choices[1]
    updateSelectInput(session, "gene_pathway", choices = choices, selected = selected)
  }, ignoreInit = FALSE)

  selected_models <- reactive({
    models <- input$models
    if (!length(models)) "DCT" else models
  })

  ens_models <- reactive({
    if (isTRUE(input$ensemble)) selected_models() else NULL
  })

  llm_settings <- reactive({
    ollama_fn <- NULL
    provider <- input$llm_provider_ui %||% "openai"
    model <- input$llm_model %||% "gpt-4o-mini"
    if (identical(provider, "ollama")) {
      provider <- "custom"
      model <- input$ollama_model %||% "qwen2.5:7b"
      base_url <- sub("/+$", "", input$ollama_url %||% "http://localhost:11434")
      base_url <- sub("/v1$", "", base_url)
      endpoint <- if (grepl("/api/generate$", base_url)) {
        base_url
      } else {
        paste0(base_url, "/api/generate")
      }
      ollama_fn <- function(prompt) {
        req <- httr2::request(endpoint) |>
          httr2::req_body_json(list(model = model, prompt = prompt, stream = FALSE)) |>
          httr2::req_error(is_error = function(resp) FALSE)
        resp <- tryCatch(
          httr2::req_perform(req),
          error = function(e) {
            msg <- paste0(
              "Ollama request failed at ", endpoint, ". ",
              "Check that Ollama is running and reachable from R."
            )
            stop(msg, call. = FALSE)
          }
        )
        status <- httr2::resp_status(resp)
        if (status >= 400) {
          body_txt <- tryCatch(httr2::resp_body_string(resp), error = function(err) "")
          msg <- paste0(
            "Ollama request failed at ", endpoint, " with HTTP ", status, ". ",
            "Use Ollama URL http://localhost:11434 and a locally installed model ",
            "(installed here: qwen2.5:7b)."
          )
          if (nzchar(body_txt)) msg <- paste0(msg, "\nOllama response body:\n", body_txt)
          stop(msg, call. = FALSE)
        }
        parsed <- httr2::resp_body_json(resp, simplifyVector = TRUE)
        parsed$response %||% ""
      }
    }
    list(
      provider = provider,
      model = model,
      fn = ollama_fn,
      disease = input$disease,
      outcome.context = input$outcome_context,
      pubmed = isTRUE(input$llm_pubmed),
      pubmed.max = input$llm_pubmed_max,
      email = if (nzchar(input$llm_email %||% "")) input$llm_email else NULL
    )
  })

  gene_models <- reactive({
    models <- intersect(selected_models(), c("RF", "EN", "GB", "DCT"))
    if (!length(models)) "DCT" else models
  })

  llm_mode_for_direct_demo <- reactive({
    mode <- input$llm_mode %||% "prior"
    if (identical(mode, "none")) "prior" else mode
  })

  test_files <- reactive({
    data_path <- if (!is.null(input$test_data_file)) input$test_data_file$datapath else default_test_data_file
    outcome_path <- if (!is.null(input$test_outcome_file)) input$test_outcome_file$datapath else default_test_outcome_file
    list(data = data_path, outcome = outcome_path,
         bundled = is.null(input$test_data_file) && is.null(input$test_outcome_file))
  })

  external_test <- reactive({
    files <- test_files()
    validate(
      need(nzchar(files$data) && file.exists(files$data),
           "No test data file is available. Upload an external test data CSV."),
      need(nzchar(files$outcome) && file.exists(files$outcome),
           "No test outcome file is available. Upload an external test outcome CSV.")
    )
    outcome_info <- read_demo_outcome_table(files$outcome)
    test_x <- read_demo_csv(
      files$data,
      reference_cols = colnames(demo_data),
      outcome_samples = outcome_info$sample
    )
    test_y <- outcome_info$outcome
    if (!is.null(outcome_info$sample) && all(outcome_info$sample %in% rownames(test_x))) {
      test_x <- test_x[outcome_info$sample, , drop = FALSE]
    }
    validate(
      need(nrow(test_x) == length(test_y),
           "External test data rows must match the length of the external test outcome."),
      need(length(intersect(colnames(test_x), colnames(demo_data))) > 0,
           "External test data must share feature/gene column names with the training data.")
    )
    list(data = test_x, outcome = test_y, sample = outcome_info$sample, bundled = files$bundled)
  })

  output$data_summary <- renderPrint({
    cat("Samples:", nrow(demo_data), "\n")
    cat("Features:", ncol(demo_data), "\n")
    cat("Available pathways:", length(demo_pathways), "\n")
    cat("Binary outcome counts:\n")
    print(table(demo_outcome))
  })

  pathway_result <- eventReactive(input$run_pathway, {
    withProgress(message = "Running pathway_select()", value = 0.2, {
      CASPEN::pathway_select(
        outcome.type = "binary",
        outcome = demo_outcome,
        data = demo_data,
        features = selected_features(),
        models.indiv = selected_models(),
        iter = min(input$iter, 5),
        num.folds = input$folds,
        ensemble = isTRUE(input$ensemble),
        models.ens = ens_models(),
        parallel.pathways = isTRUE(input$parallel_pathways),
        parallel.iter = isTRUE(input$parallel_iter),
        workers = input$workers,
        future.seed = TRUE,
        llm.mode = input$llm_mode,
        disease = llm_settings()$disease,
        outcome.context = llm_settings()$outcome.context,
        llm.provider = llm_settings()$provider,
        llm.model = llm_settings()$model,
        llm.fn = llm_settings()$fn,
        llm.pubmed = llm_settings()$pubmed,
        llm.pubmed.max = llm_settings()$pubmed.max,
        llm.email = llm_settings()$email,
        llm.valid.genes = colnames(demo_data),
        llm.min.genes = input$llm_min_genes,
        llm.prior.method = input$llm_prior_method,
        llm.temperature = NULL
      )
    })
  })

  output$path_auc <- renderDT({
    req(pathway_result())
    datatable(round(as.data.frame(pathway_result()$data.metric), 3),
              options = list(pageLength = 8, scrollX = TRUE))
  })

  output$path_sp95 <- renderDT({
    req(pathway_result())
    datatable(round(as.data.frame(pathway_result()$SP.95), 3),
              options = list(pageLength = 8, scrollX = TRUE))
  })

  output$path_plot <- renderPlot({
    req(pathway_result())
    plot_pathway_selection(pathway_result(), outcome.type = "binary", label.top = 6)
  })

  llm_demo_result <- eventReactive(input$run_llm_demo, {
    settings <- llm_settings()
    withProgress(message = "Running llm_literature_signatures()", value = 0.2, {
      CASPEN::llm_literature_signatures(
        disease = settings$disease,
        outcome.context = settings$outcome.context,
        features = selected_features(),
        mode = llm_mode_for_direct_demo(),
        provider = settings$provider,
        model = settings$model,
        llm.fn = settings$fn,
        pubmed = settings$pubmed,
        pubmed.max = settings$pubmed.max,
        email = settings$email,
        valid.genes = colnames(demo_data),
        min.genes = input$llm_min_genes,
        prior.method = input$llm_prior_method,
        temperature = NULL
      )
    })
  })

  output$llm_demo_evidence <- renderDT({
    req(llm_demo_result())
    ev <- llm_demo_result()$evidence
    if (is.null(ev) || !nrow(ev)) {
      ev <- llm_demo_result()$raw.signatures
      if (!is.null(ev) && nrow(ev)) {
        ev$note <- paste0("Shown from raw LLM output; filtered results may be empty because Min Curated Genes = ",
                          input$llm_min_genes, ".")
      } else {
        ev <- llm_demo_result()$raw.signatures.all
        if (!is.null(ev) && nrow(ev)) {
          ev$status <- paste0("Below Min Curated Genes = ", input$llm_min_genes,
                              "; not printed as a curated pathway and not used.")
        } else {
          ev <- data.frame(
            note = "No parseable LLM signatures were returned.",
            stringsAsFactors = FALSE
          )
        }
      }
    }
    numeric_cols <- vapply(ev, is.numeric, logical(1))
    ev[numeric_cols] <- lapply(ev[numeric_cols], round, 3)
    datatable(ev, options = list(pageLength = 8, scrollX = TRUE))
  })

  output$llm_demo_prior <- renderDT({
    req(llm_demo_result())
    if (identical(llm_mode_for_direct_demo(), "curate")) {
      return(datatable(
        data.frame(note = "Curate mode returns new literature signatures only; it does not print or score prior pathways."),
        options = list(dom = "t")
      ))
    }
    prior <- llm_demo_result()$pathway.prior
    tab <- if (length(prior) > 0) {
      data.frame(pathway = names(prior), prior = as.numeric(prior),
                 stringsAsFactors = FALSE)
    } else {
      data.frame(pathway = character(0), prior = numeric(0),
                 stringsAsFactors = FALSE)
    }
    tab$prior <- round(tab$prior, 3)
    datatable(tab, options = list(pageLength = 8, scrollX = TRUE))
  })

  output$llm_demo_signatures <- renderDT({
    req(llm_demo_result())
    sigs <- llm_demo_result()$literature.features
    tab <- if (length(sigs) > 0) {
      data.frame(
        signature = names(sigs),
        n.genes = lengths(sigs),
        genes = vapply(sigs, paste, collapse = ", ", FUN.VALUE = character(1)),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(signature = character(0), n.genes = integer(0),
                 genes = character(0), stringsAsFactors = FALSE)
    }
    datatable(tab, options = list(pageLength = 8, scrollX = TRUE))
  })

  gene_result <- eventReactive(input$run_gene, {
    features <- tryCatch(pathway_result()$features, error = function(e) NULL)
    if (is.null(features)) features <- selected_features()
    gene_pathway <- input$gene_pathway %||% names(features)[1]
    if (!gene_pathway %in% names(features)) {
      features <- selected_features()
    }
    features <- features[intersect(gene_pathway, names(features))]
    validate(need(length(features) == 1, "Choose one pathway for gene selection."))
    withProgress(message = "Running gene_select()", value = 0.2, {
      gene_select(
        outcome.type = "binary",
        outcome = demo_outcome,
        data = demo_data,
        features = features,
        models = gene_models(),
        iter = min(input$iter, 3),
        num.folds = input$folds,
        min.pathway.size = 10,
        max.genes = 10,
        gene.prior.method = input$gene_prior_method,
        disease = llm_settings()$disease,
        outcome.context = llm_settings()$outcome.context,
        llm.provider = llm_settings()$provider,
        llm.model = llm_settings()$model,
        llm.fn = llm_settings()$fn,
        llm.pubmed = llm_settings()$pubmed,
        llm.email = llm_settings()$email,
        llm.temperature = NULL,
        verbose = FALSE
      )
    })
  })

  output$gene_selected <- renderDT({
    req(gene_result())
    selected <- gene_result()$selected.features
    tab <- data.frame(
      pathway = names(selected),
      genes = vapply(selected, paste, collapse = ", ", FUN.VALUE = character(1)),
      stringsAsFactors = FALSE
    )
    datatable(tab, options = list(pageLength = 8, scrollX = TRUE))
  })

  output$gene_ranking <- renderDT({
    req(gene_result())
    datatable(gene_result()$ranking, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$gene_plot <- renderPlot({
    req(gene_result())
    try(plot_gene_selection(gene_result(), label.top = 8), silent = TRUE)
  })

  output$llm_evidence <- renderDT({
    req(pathway_result())
    ev <- pathway_result()$llm$evidence
    if (is.null(ev) || !nrow(ev)) {
      ev <- pathway_result()$llm$raw.signatures
      if (!is.null(ev) && nrow(ev)) {
        ev$note <- paste0("Shown from raw LLM output; filtered results may be empty because Min Curated Genes = ",
                          input$llm_min_genes, ".")
      } else {
        ev <- pathway_result()$llm$raw.signatures.all
        if (!is.null(ev) && nrow(ev)) {
          ev$status <- paste0("Below Min Curated Genes = ", input$llm_min_genes,
                              "; not printed as a curated pathway and not used.")
        } else {
          ev <- data.frame(
            note = "No parseable LLM signatures were returned.",
            stringsAsFactors = FALSE
          )
        }
      }
    }
    numeric_cols <- vapply(ev, is.numeric, logical(1))
    ev[numeric_cols] <- lapply(ev[numeric_cols], round, 3)
    datatable(ev, options = list(pageLength = 8, scrollX = TRUE))
  })

  output$llm_prior <- renderDT({
    req(pathway_result())
    if (identical(input$llm_mode, "curate")) {
      return(datatable(
        data.frame(note = "Curate mode returns new literature signatures only; pathway prior weights are not printed."),
        options = list(dom = "t")
      ))
    }
    prior <- pathway_result()$Prior.weight
    tab <- if (!is.null(prior)) {
      round(as.data.frame(prior), 3)
    } else {
      data.frame(prior.weight = rep(NA_real_, length(input$pathways)),
                 row.names = input$pathways)
    }
    datatable(tab,
              options = list(pageLength = 8, scrollX = TRUE))
  })

  output$gene_prior_table <- renderDT({
    req(gene_result())
    tab <- gene_result()$ranking
    cols <- intersect(
      c("pathway", "gene", "pubmed.count", "pubmed.score", "pubtator.count",
        "pubtator.score", "llm.score", "evidence.score", "hybrid.score",
        "gene.prior.method"),
      names(tab)
    )
    validate(need(length(cols) > 2,
                  "Run gene_select() with Gene Prior Method other than none."))
    out <- tab[, cols, drop = FALSE]
    numeric_cols <- vapply(out, is.numeric, logical(1))
    out[numeric_cols] <- lapply(out[numeric_cols], round, 3)
    datatable(out, options = list(pageLength = 10, scrollX = TRUE))
  })

  train_result <- eventReactive(input$run_train, {
    features <- tryCatch(pathway_result()$features, error = function(e) NULL)
    if (is.null(features)) features <- selected_features()
    withProgress(message = "Running Train_perform_CV()", value = 0.2, {
      Train_perform_CV(
        outcome.type = "binary",
        outcome = demo_outcome,
        train_data = demo_data,
        features = features,
        models.indiv = selected_models(),
        iter = min(input$iter, 5),
        num.folds = input$folds,
        ensemble = isTRUE(input$ensemble),
        models.ens = ens_models(),
        parallel.iter = isTRUE(input$parallel_iter),
        workers = input$workers,
        future.seed = TRUE
      )
    })
  })

  output$train_table <- renderDT({
    req(train_result())
    datatable(train_result(), options = list(pageLength = 10, scrollX = TRUE))
  })

  output$test_summary <- renderPrint({
    test_obj <- external_test()
    cat("Training data: bundled CASPEN example\n")
    cat("Training samples:", nrow(demo_data), "\n")
    cat("Training outcome counts:\n")
    print(table(demo_outcome))
    cat("\nTest data source:", if (isTRUE(test_obj$bundled)) "bundled protein_data_Test.csv" else "uploaded CSV", "\n")
    cat("External test samples:", nrow(test_obj$data), "\n")
    cat("External test features:", ncol(test_obj$data), "\n")
    cat("Shared training/test features:", length(intersect(colnames(test_obj$data), colnames(demo_data))), "\n")
    cat("External test outcome counts:\n")
    print(table(test_obj$outcome))
  })

  test_result <- eventReactive(input$run_test, {
    features <- tryCatch(pathway_result()$features, error = function(e) NULL)
    if (is.null(features)) features <- selected_features()
    test_obj <- external_test()
    withProgress(message = "Running Test_perform()", value = 0.2, {
      Test_perform(
        outcome.type = "binary",
        outcome = demo_outcome,
        train_data = demo_data,
        test_data = test_obj$data,
        test_outcome = test_obj$outcome,
        features = features,
        models.indiv = selected_models(),
        ensemble = isTRUE(input$ensemble),
        models.ens = ens_models()
      )
    })
  })

  output$test_table <- renderDT({
    req(test_result())
    datatable(test_result(), options = list(pageLength = 10, scrollX = TRUE))
  })

  output$path_code <- renderText({
    code_block(
      "library(CASPEN)",
      "data(caspen_example)",
      "path_res <- pathway_select(",
      "  outcome.type = 'binary',",
      "  outcome = caspen_example$outcome,",
      "  data = caspen_example$data_to_pass,",
      paste0("  features = caspen_example$pathways[c(",
             paste(shQuote(input$pathways), collapse = ", "), ")],"),
      paste0("  models.indiv = c(", paste(shQuote(selected_models()), collapse = ", "), "),"),
      paste0("  iter = ", min(input$iter, 5), ","),
      paste0("  num.folds = ", input$folds, ","),
      paste0("  ensemble = ", isTRUE(input$ensemble), ","),
      paste0("  parallel.pathways = ", isTRUE(input$parallel_pathways), ","),
      paste0("  parallel.iter = ", isTRUE(input$parallel_iter), ","),
      paste0("  workers = ", input$workers, ","),
      paste0("  llm.mode = ", shQuote(input$llm_mode), ","),
      paste0("  llm.prior.method = ", shQuote(input$llm_prior_method), ","),
      paste0("  llm.min.genes = ", input$llm_min_genes, ","),
      paste0("  disease = ", shQuote(input$disease), ","),
      paste0("  outcome.context = ", shQuote(input$outcome_context)),
      ")"
    )
  })

  output$llm_demo_code <- renderText({
    tail_lines <- if (identical(llm_mode_for_direct_demo(), "curate")) {
      c("llm_res$evidence", "llm_res$literature.features")
    } else {
      c("llm_res$evidence", "llm_res$pathway.prior", "llm_res$literature.features")
    }
    code_block(
      "llm_res <- llm_literature_signatures(",
      paste0("  disease = ", shQuote(input$disease), ","),
      paste0("  outcome.context = ", shQuote(input$outcome_context), ","),
      paste0("  features = caspen_example$pathways[c(",
             paste(shQuote(input$pathways), collapse = ", "), ")],"),
      paste0("  mode = ", shQuote(llm_mode_for_direct_demo()), ","),
      paste0("  prior.method = ", shQuote(input$llm_prior_method), ","),
      paste0("  min.genes = ", input$llm_min_genes, ","),
      paste0("  pubmed = ", isTRUE(input$llm_pubmed), ","),
      paste0("  pubmed.max = ", input$llm_pubmed_max),
      ")",
      "",
      tail_lines
    )
  })

  output$gene_code <- renderText({
    code_block(
      "gene_res <- gene_select(",
      "  outcome.type = 'binary',",
      "  outcome = caspen_example$outcome,",
      "  data = caspen_example$data_to_pass,",
      paste0("  features = path_res$features[", shQuote(input$gene_pathway), "],"),
      paste0("  models = c(", paste(shQuote(gene_models()), collapse = ", "), "),"),
      paste0("  iter = ", min(input$iter, 3), ","),
      paste0("  num.folds = ", input$folds, ","),
      "  min.pathway.size = 10,",
      "  max.genes = 10,",
      paste0("  gene.prior.method = ", shQuote(input$gene_prior_method), ","),
      paste0("  disease = ", shQuote(input$disease), ","),
      paste0("  outcome.context = ", shQuote(input$outcome_context)),
      ")"
    )
  })

  output$train_code <- renderText({
    code_block(
      "train_res <- Train_perform_CV(",
      "  outcome.type = 'binary',",
      "  outcome = caspen_example$outcome,",
      "  train_data = caspen_example$data_to_pass,",
      "  features = path_res$features,",
      paste0("  models.indiv = c(", paste(shQuote(selected_models()), collapse = ", "), "),"),
      paste0("  iter = ", min(input$iter, 5), ","),
      paste0("  num.folds = ", input$folds, ","),
      paste0("  ensemble = ", isTRUE(input$ensemble), ","),
      paste0("  parallel.iter = ", isTRUE(input$parallel_iter), ","),
      paste0("  workers = ", input$workers),
      ")"
    )
  })

  output$test_code <- renderText({
    code_block(
      "test_raw <- read.csv('protein_data_Test.csv', check.names = FALSE)",
      "test_out <- read.csv('outcome_binary_Test.csv', check.names = FALSE)",
      "test_x <- as.data.frame(t(as.matrix(test_raw[-1])), check.names = FALSE)",
      "colnames(test_x) <- make.unique(as.character(test_raw[[1]]))",
      "rownames(test_x) <- colnames(test_raw)[-1]",
      "test_x <- test_x[test_out$sample, , drop = FALSE]",
      "test_y <- as.numeric(test_out$outcome)",
      "",
      "test_res <- Test_perform(",
      "  outcome.type = 'binary',",
      "  outcome = caspen_example$outcome,",
      "  train_data = caspen_example$data_to_pass,",
      "  test_data = test_x,",
      "  test_outcome = test_y,",
      "  features = path_res$features,",
      paste0("  models.indiv = c(", paste(shQuote(selected_models()), collapse = ", "), "),"),
      paste0("  ensemble = ", isTRUE(input$ensemble), ","),
      paste0("  models.ens = c(", paste(shQuote(selected_models()), collapse = ", "), ")"),
      ")"
    )
  })
}

shinyApp(ui, server)
