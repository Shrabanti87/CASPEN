# CASPEN

**CASPEN** (**c**elltype **a**ware **s**tacked **p**athway-guided
**en**semble) is an R package for building interpretable pathway-level
prediction models from bulk or celltype-specific omics data.

CASPEN supports **binary**, **survival**, and **categorical** outcomes. It can
evaluate individual models and stacked ensemble models, report AUC or C-index,
calculate sensitivity at 95% specificity (SP95), select pathways, select genes
within pathways, and optionally use literature/LLM priors.

## What CASPEN Is For

CASPEN is designed for analyses where the user has:

- an omics matrix, such as RNA, protein, or celltype-specific expression;
- pathway or gene-set definitions, for example MSigDB/Reactome/KEGG/user lists;
- an outcome, such as response/non-response, survival time/status, or subtype;
- optional celltype names and aliases;
- optional literature prior information.

The main idea is:

```text
omics data + pathway definitions + outcome
        -> pathway-celltype models
        -> base learners and ensembles
        -> AUC/C-index, SP95, pathway selection, gene selection
```

## Important LLM Design Choice

CASPEN should **not** rely on an LLM to invent complete gene sets from scratch.
LLMs often return very small marker-style lists, such as 2-3 genes, which are
not ideal for pathway-based prediction.

The recommended CASPEN workflow is:

1. Provide curated pathways/gene sets, such as MSigDB signatures.
2. Let the LLM curate disease/outcome-relevant **biological concepts** and/or
   score prior evidence.
3. CASPEN maps those concepts back to the user-provided pathway list and uses
   the curated pathway genes for modeling.

For example:

```text
LLM concept: AKT_activation
Matched user/MSigDB pathway: HALLMARK_PI3K_AKT_MTOR_SIGNALING
Genes used for modeling: genes from the supplied MSigDB pathway
```

This keeps CASPEN literature-aware while keeping gene membership curated,
interpretable, and large enough for prediction.

## Installation

Install from GitHub:

```r
install.packages("remotes")
remotes::install_github("Shrabanti87/CASPEN")
```

Install from a local checkout:

```r
install.packages("remotes")
remotes::install_local(".")
```

Load the package:

```r
library(CASPEN)
```

Check model backend versions:

```r
caspen_model_info()
```

## Main Functions

- `pathway_select()`: select predictive pathways or pathway-celltypes.
- `gene_select()`: select genes within selected pathways.
- `Train_perform_CV()`: estimate training performance by cross-validation.
- `Test_perform()`: evaluate independent test-set performance.
- `llm_literature_signatures()`: curate/score literature priors using PubMed,
  PubTator3, OpenAI, or a custom local LLM such as Ollama.
- `plot_pathway_selection()`: plot AUC/C-index by SP95.
- `plot_gene_selection()`: plot gene-selection ranking/stability.
- `plot_roc_curve()`: plot ROC curves.
- `launch_caspen_demo()`: open the bundled Shiny demo.

## Supported Model Abbreviations

CASPEN can use these base learners, depending on the outcome type:

- `XG`: XGBoost
- `RF`: random forest
- `EN`: elastic net
- `GB`: gradient boosting
- `KNN`: k-nearest neighbors
- `SVM`: support vector machine
- `NB`: naive Bayes
- `DCT`: decision tree
- `NN`: neural network
- `ADB`: AdaBoost
- `MB`: model-based boosting

Ensemble outputs use names such as `RF.ens`, `SVM.ens`, `XG.ens`,
`SimpleAvg`, and `WeightedAvg`, depending on the models requested.

## Example Data

CASPEN includes small example datasets for package demos:

```r
data(caspen_example)
data(caspen_test_example)

x <- caspen_example$data_to_pass
y <- caspen_example$outcome
pathways <- caspen_example$pathways

test_x <- caspen_test_example$data_to_pass
test_y <- caspen_test_example$outcome
```

Expected input format:

- `data`: samples as rows and features/genes as columns.
- `outcome`: vector, factor, or survival outcome object/data frame expected by
  the chosen function.
- `features`: named list where each element is a character vector of genes.

Example:

```r
pathways <- list(
  AKT = c("AKT1", "PIK3CA", "PTEN", "MTOR"),
  DNA_repair = c("BRCA1", "BRCA2", "RAD51", "ERCC1", "XRCC3")
)
```

## 1. Pathway Selection For Binary Outcomes

```r
set.seed(1)

path_res <- pathway_select(
  outcome.type = "binary",
  outcome = y,
  data = x,
  features = pathways[1:5],
  models.indiv = c("RF", "EN", "DCT"),
  iter = 5,
  num.folds = 3,
  ensemble = TRUE,
  models.ens = c("RF", "EN"),
  num.sel.path = 3,
  fs.AUC.cut = 0.60,
  fs.sp.cut = 0.05
)

path_res$auc
path_res$SP.95
path_res$selection.table
```

Plot pathway selection:

```r
plot_pathway_selection(path_res)
```

## 2. Survival Outcomes

For survival outcomes, CASPEN uses C-index as the primary metric when
available, while also supporting time-dependent AUC/SP95 where implemented.

```r
surv_res <- pathway_select(
  outcome.type = "survival",
  outcome = survival_outcome,
  survdays = "time",
  data = x,
  features = pathways[1:5],
  models.indiv = c("RF", "EN", "GB"),
  iter = 5,
  num.folds = 3,
  ensemble = TRUE,
  models.ens = c("RF", "EN")
)

surv_res$cindex
surv_res$SP.95
```

## 3. Categorical Outcomes

```r
cat_res <- pathway_select(
  outcome.type = "categorical",
  outcome = subtype,
  data = x,
  features = pathways[1:5],
  models.indiv = c("RF", "EN", "DCT"),
  iter = 5,
  num.folds = 3
)

cat_res$auc
cat_res$SP.95
```

## 4. Celltype-Aware Pathway Mapping

If the data columns encode gene and celltype-specific expression, CASPEN can
evaluate pathway-celltype tasks. Provide the celltypes and optional aliases.

```r
celltypes <- c("CAF", "CD4_Tcells", "CD8_Tcells", "Macrophages")

celltype_aliases <- list(
  CAF = c("caf", "fibroblast", "stromal"),
  CD4_Tcells = c("cd4", "tcell", "t cells", "immune"),
  CD8_Tcells = c("cd8", "tcell", "t cells", "immune"),
  Macrophages = c("macrophage", "myeloid", "immune")
)

ct_res <- pathway_select(
  outcome.type = "binary",
  outcome = y,
  data = celltype_specific_x,
  features = pathways,
  celltypes = celltypes,
  celltype.aliases = celltype_aliases,
  models.indiv = c("RF", "EN"),
  iter = 5,
  num.folds = 3
)
```

Rules:

- If a pathway name does not contain a celltype/alias, CASPEN can map it across
  the supplied celltypes.
- If a pathway name contains a celltype/alias, CASPEN restricts that pathway to
  matching celltypes.
- Users control mapping through `celltype.aliases`.

## 5. Literature Priors With Existing MSigDB/User Pathways

If your pathway list already contains MSigDB signatures, use LLM/literature
support mainly for **prior weights** and **concept matching**, not for final
gene membership.

### PubMed count prior, no LLM

```r
pubmed_res <- llm_literature_signatures(
  disease = "high-grade serous ovarian cancer",
  outcome.context = "platinum resistance",
  features = pathways,
  mode = "prior",
  prior.method = "pubmed_count",
  email = "your.email@example.com"
)

pubmed_res$evidence[, c("pathway", "pubmed.count", "pubmed.score", "prior")]
```

### PubTator3 prior, no LLM

PubTator3 is an online NCBI service; no local PubTator installation is needed.

```r
pubtator_res <- llm_literature_signatures(
  disease = "high-grade serous ovarian cancer",
  outcome.context = "platinum resistance",
  features = pathways,
  mode = "prior",
  prior.method = "pubtator",
  pubmed.max = 20,
  email = "your.email@example.com"
)

pubtator_res$evidence[, c("pathway", "pubtator.count",
                          "pubtator.score", "prior")]
```

### Hybrid PubMed + PubTator + LLM prior

```r
hybrid_res <- llm_literature_signatures(
  disease = "high-grade serous ovarian cancer",
  outcome.context = "platinum resistance",
  features = pathways,
  mode = "prior",
  prior.method = "hybrid",
  hybrid.weights = c(pubmed_count = 0.25, pubtator = 0.25, llm = 0.50),
  model = "gpt-4o-mini",
  pubmed = FALSE,
  temperature = NULL
)

hybrid_res$evidence
hybrid_res$pathway.prior
```

### OpenAI setup

Create an API key from the OpenAI dashboard, then set it in R:

```r
Sys.setenv(OPENAI_API_KEY = "sk-proj-your_real_key_here")
```

Then call:

```r
llm_res <- llm_literature_signatures(
  disease = "high-grade serous ovarian cancer",
  outcome.context = "platinum resistance",
  features = pathways,
  mode = "prior",
  model = "gpt-4o-mini",
  pubmed = FALSE,
  temperature = NULL
)
```

### Ollama/local LLM setup

Install Ollama from <https://ollama.com/download>, then in Terminal:

```bash
ollama pull qwen2.5:7b
ollama serve
```

In R:

```r
ollama_llm <- function(prompt) {
  req <- httr2::request("http://localhost:11434/api/generate") |>
    httr2::req_body_json(list(
      model = "qwen2.5:7b",
      prompt = prompt,
      stream = FALSE
    ))
  res <- httr2::req_perform(req)
  out <- jsonlite::fromJSON(httr2::resp_body_string(res))
  out$response
}

ollama_res <- llm_literature_signatures(
  disease = "high-grade serous ovarian cancer",
  outcome.context = "platinum resistance",
  features = pathways,
  mode = "prior",
  provider = "custom",
  llm.fn = ollama_llm,
  pubmed = FALSE
)
```

## LLM Modes

- `llm.mode = "none"`: no LLM/literature curation.
- `llm.mode = "prior"`: score user-provided pathways only.
- `llm.mode = "curate"`: curate disease/outcome-relevant concepts. If user
  pathways/MSigDB lists are supplied, use those as the trusted source of genes.
- `llm.mode = "both"`: score user pathways and curate additional concepts.

CASPEN reports literature support separately from data performance:

```r
path_res$data.metric    # AUC or C-index matrix
path_res$SP.95          # SP95 matrix
path_res$Prior.weight   # literature prior matrix/vector
path_res$llm$evidence   # evidence table
```

The user can then apply study-specific cutoffs for data metric, SP95, and prior
weight.

## 6. Gene Selection Within Pathways

`gene_select()` ranks genes within selected pathways using native model
importance where available. Importance scores from different models are
converted to within-run percentiles, then summarized with stability.

```r
gene_res <- gene_select(
  outcome.type = "binary",
  outcome = y,
  data = x,
  features = path_res$features,
  models = c("RF", "EN", "DCT", "GB"),
  iter = 10,
  num.folds = 3,
  min.pathway.size = 30,
  gene.stability.quantile = 0.75,
  gene.min.stability = 0.25,
  max.genes = 30
)

head(gene_res$ranking)
gene_res$selected.features
plot_gene_selection(gene_res)
```

Gene-level literature priors can also be requested:

```r
gene_res <- gene_select(
  outcome.type = "binary",
  outcome = y,
  data = x,
  features = path_res$features,
  models = c("RF", "EN"),
  iter = 10,
  num.folds = 3,
  gene.prior.method = "hybrid",
  disease = "high-grade serous ovarian cancer",
  outcome.context = "platinum resistance"
)
```

## 7. Training CV And Independent Test Performance

```r
train_res <- Train_perform_CV(
  outcome.type = "binary",
  outcome = y,
  train_data = x,
  features = path_res$features,
  models.indiv = c("RF", "EN", "DCT"),
  iter = 10,
  num.folds = 3,
  ensemble = TRUE,
  models.ens = c("RF", "EN")
)

train_res$auc
train_res$SP.95
```

Independent test set:

```r
test_res <- Test_perform(
  outcome.type = "binary",
  outcome = y,
  train_data = x,
  test_data = test_x,
  test_outcome = test_y,
  features = path_res$features,
  models.indiv = c("RF", "EN", "DCT"),
  ensemble = TRUE,
  models.ens = c("RF", "EN")
)

test_res$auc
test_res$SP.95
```

## 8. Parallel Computing

CASPEN uses the `future` backend.

Pathway-level parallelism:

```r
path_res <- pathway_select(
  outcome.type = "binary",
  outcome = y,
  data = x,
  features = pathways,
  models.indiv = c("RF", "EN"),
  iter = 100,
  num.folds = 3,
  parallel.pathways = TRUE,
  workers = future::availableCores()
)
```

Iteration-level parallelism:

```r
train_res <- Train_perform_CV(
  outcome.type = "binary",
  outcome = y,
  train_data = x,
  features = path_res$features,
  models.indiv = c("RF", "EN"),
  iter = 100,
  num.folds = 3,
  parallel.iter = TRUE,
  workers = future::availableCores()
)
```

## 9. Shiny Demo

CASPEN includes a lightweight Shiny app for demos and collaborator
walkthroughs:

```r
library(CASPEN)
launch_caspen_demo()
```

The app demonstrates pathway selection, gene selection, training/test
performance, plots, LLM/PubMed/PubTator priors, and the R code used for each
run. It is intended for small examples. Use R scripts or HPC/server jobs for
large analyses.


