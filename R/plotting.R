#' Plot Pathway Selection Performance
#'
#' Draws a scatter plot of pathway-level prediction performance against
#' sensitivity at 95 percent specificity. For survival results, the y-axis uses
#' C-index when available; otherwise it uses AUC.
#'
#' @param x A result object from `pathway_select()`, `Train_perform_CV()`, or
#' `Test_perform()` containing `auc` or `C.index`/`cindex` and `SP.95`/`sp.95`.
#' @param outcome.type One of `"auto"`, `"binary"`, `"categorical"`, or
#' `"survival"`. `"auto"` uses C-index if present.
#' @param model Optional model column/name to plot. If `NULL`, row means across
#' all available models are used.
#' @param label.top Number of top-performing points to label.
#' @param main Optional plot title.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param selected.color Point color.
#' @param ... Additional arguments passed to `plot()`.
#'
#' @return Invisibly returns the plotted data frame.
#'
#' @examples
#' \dontrun{
#' data(caspen_example)
#'
#' path_res <- pathway_select(
#'   outcome.type = "binary",
#'   outcome = caspen_example$outcome,
#'   data = caspen_example$data_to_pass,
#'   features = caspen_example$pathways[1:5],
#'   models.indiv = c("RF", "EN"),
#'   iter = 5,
#'   num.folds = 3
#' )
#'
#' plot_pathway_selection(path_res, outcome.type = "binary", label.top = 5)
#'
#' # For survival pathway_select() results, the y-axis becomes C-index.
#' plot_pathway_selection(path_res, outcome.type = "survival", model = "RF")
#' }
#'
#' @export
plot_pathway_selection <- function(x, outcome.type = c("auto", "binary", "categorical", "survival"),
                                   model = NULL, label.top = 8, main = NULL,
                                   xlab = "Sensitivity at 95% specificity",
                                   ylab = NULL, selected.color = "#1f7a3d", ...) {
  outcome.type <- match.arg(outcome.type)
  sp95 <- x$SP.95 %||% x$sp.95
  auc <- x$auc %||% x$AUC
  cindex <- x$C.index %||% x$cindex

  if (is.null(sp95)) stop("x must contain SP.95 or sp.95.")
  if (outcome.type == "auto") outcome.type <- if (!is.null(cindex)) "survival" else "binary"
  metric <- if (outcome.type == "survival" && !is.null(cindex)) cindex else auc
  if (is.null(metric)) stop("x must contain auc/AUC, or C.index/cindex for survival.")

  sp <- .metric_vector(sp95, model)
  y <- .metric_vector(metric, model)
  labels <- names(y) %||% names(sp) %||% rownames(as.matrix(metric)) %||% paste0("Feature_", seq_along(y))
  labels <- labels[seq_along(y)]
  keep <- is.finite(sp) & is.finite(y)
  if (!any(keep)) stop("No finite pathway metrics available to plot.")

  sp <- sp[keep]
  y <- y[keep]
  labels <- labels[keep]
  if (is.null(ylab)) ylab <- if (outcome.type == "survival" && !is.null(cindex)) "C-index" else "AUC"
  if (is.null(main)) main <- paste(ylab, "by SP95")

  old.par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old.par), add = TRUE)
  graphics::par(mar = c(5.2, 5.2, 4.2, 2.4), cex = 1.2, cex.lab = 1.25,
                cex.axis = 1.1, cex.main = 1.2, xpd = NA)
  xlim <- .padded_range(sp, lower = 0, upper = 1)
  ylim <- .padded_range(y, lower = if (min(y, na.rm = TRUE) >= 0) 0 else NULL,
                        upper = if (max(y, na.rm = TRUE) <= 1) 1 else NULL)

  graphics::plot(
    sp, y,
    pch = 21, bg = selected.color, col = "white",
    xlab = xlab, ylab = ylab, main = main,
    xlim = xlim, ylim = ylim,
    cex = 1.35,
    ...
  )
  graphics::abline(h = 0.5, lty = 2, col = "grey65")

  if (label.top > 0) {
    ord <- order(y, decreasing = TRUE, na.last = NA)
    ord <- ord[seq_len(min(label.top, length(ord)))]
    graphics::text(sp[ord], y[ord], labels = labels[ord], pos = 3, cex = 0.9, offset = 0.7)
  }

  invisible(data.frame(feature = labels, sp95 = sp, metric = y, stringsAsFactors = FALSE))
}

#' Plot Gene Selection Stability
#'
#' Draws the gene-level selection scatter from a `gene_select()` result.
#'
#' @param x A result object from `gene_select()`.
#' @param pathway Optional pathway name(s) to plot.
#' @param label.top Number of selected/high-ranking genes to label.
#' @param main Optional plot title.
#' @param selected.color Point color for selected genes.
#' @param other.color Point color for unselected genes.
#' @param ... Additional arguments passed to `plot()`.
#'
#' @return Invisibly returns the plotted ranking table.
#'
#' @examples
#' \dontrun{
#' data(caspen_example)
#'
#' gene_res <- gene_select(
#'   outcome.type = "binary",
#'   outcome = caspen_example$outcome,
#'   data = caspen_example$data_to_pass,
#'   features = caspen_example$pathways[1],
#'   models = c("RF", "EN", "DCT"),
#'   iter = 3,
#'   num.folds = 3,
#'   min.pathway.size = 10,
#'   max.genes = 10
#' )
#'
#' plot_gene_selection(gene_res, label.top = 10)
#' }
#'
#' @export
plot_gene_selection <- function(x, pathway = NULL, label.top = 12, main = NULL,
                                selected.color = "#1f7a3d",
                                other.color = "#9ca3af", ...) {
  tab <- x$ranking
  if (is.null(tab)) stop("x must be a gene_select() result with a ranking table.")
  if (!is.null(pathway)) tab <- tab[tab$pathway %in% pathway, , drop = FALSE]
  tab <- tab[is.finite(tab$stability) & is.finite(tab$median_percentile), , drop = FALSE]
  if (!nrow(tab)) stop("No finite gene-selection metrics available to plot.")
  if (is.null(main)) main <- "Gene selection: median importance percentile by stability"

  cols <- ifelse(tab$selected, selected.color, other.color)
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old.par), add = TRUE)
  graphics::par(mar = c(5.2, 5.2, 4.2, 2.4), cex = 1.2, cex.lab = 1.25,
                cex.axis = 1.1, cex.main = 1.2, xpd = NA)
  xlim <- .padded_range(tab$stability, lower = 0, upper = 1)
  ylim <- .padded_range(tab$median_percentile, lower = 0, upper = 1)

  graphics::plot(
    tab$stability, tab$median_percentile,
    pch = 21, bg = cols, col = "white",
    xlab = "Stability",
    ylab = "Median importance percentile",
    main = main,
    xlim = xlim, ylim = ylim,
    cex = 1.35,
    ...
  )
  if (!is.null(x$settings$gene.min.stability)) {
    graphics::abline(v = x$settings$gene.min.stability, lty = 2, col = "grey65")
  }
  if (!is.null(x$settings$gene.stability.quantile)) {
    graphics::abline(h = x$settings$gene.stability.quantile, lty = 3, col = "grey65")
  }
  graphics::legend(
    "bottomright", legend = c("Selected", "Not selected"),
    pt.bg = c(selected.color, other.color), pch = 21, col = "white", bty = "n",
    cex = 1.05
  )

  if (label.top > 0) {
    ord <- order(tab$selected, tab$median_percentile, decreasing = TRUE, na.last = NA)
    ord <- ord[seq_len(min(label.top, length(ord)))]
    graphics::text(tab$stability[ord], tab$median_percentile[ord],
                   labels = tab$gene[ord], pos = 3, cex = 0.9, offset = 0.7)
  }
  invisible(tab)
}

#' Plot ROC Curves
#'
#' Plots binary, multiclass, or survival ROC curves from CASPEN ROC objects or
#' prediction scores. `Train_perform_CV()` and `Test_perform()` summary tables
#' include AUC/SP95, while lower-level ROC helper results may include `roc` and
#' `predictions` elements that can be passed directly to this function.
#'
#' @param x A pROC ROC object, survivalROC object, list of ROC objects, or a
#' numeric prediction vector/matrix.
#' @param outcome Binary/categorical outcome vector. Required when `x` contains
#' predictions rather than ROC objects.
#' @param outcome.type One of `"binary"`, `"categorical"`, or `"survival"`.
#' @param survdays Survival time vector for `outcome.type = "survival"`.
#' @param time.point Prediction time for survival ROC.
#' @param model Optional model name when `x` is a named list.
#' @param main Optional plot title.
#' @param ... Additional plotting arguments.
#'
#' @return Invisibly returns ROC objects used for plotting.
#'
#' @examples
#' \dontrun{
#' data(caspen_example)
#'
#' # Plot from a numeric risk/probability vector.
#' marker <- caspen_example$data_to_pass[, 1]
#' plot_roc_curve(
#'   x = marker,
#'   outcome = caspen_example$outcome,
#'   outcome.type = "binary",
#'   main = "Example binary ROC"
#' )
#'
#' # Survival ROC at one year from a numeric risk score.
#' set.seed(1)
#' survival_time <- sample(120:1800, nrow(caspen_example$data_to_pass),
#'                         replace = TRUE)
#' survival_status <- caspen_example$outcome
#'
#' plot_roc_curve(
#'   x = marker,
#'   outcome = survival_status,
#'   outcome.type = "survival",
#'   survdays = survival_time,
#'   time.point = 365,
#'   main = "Example one-year survival ROC"
#' )
#' }
#'
#' @export
plot_roc_curve <- function(x, outcome = NULL,
                           outcome.type = c("binary", "categorical", "survival"),
                           survdays = NULL, time.point = NULL,
                           model = NULL, main = NULL, ...) {
  outcome.type <- match.arg(outcome.type)
  if (is.list(x) && !is.null(model)) x <- x[[model]]
  rocs <- .coerce_roc_objects(x, outcome, outcome.type, survdays, time.point)
  if (!length(rocs)) stop("No ROC curves could be constructed.")
  if (is.null(main)) main <- "ROC curve"

  cols <- c("#7c3aed", "#1f7a3d", "#2563eb", "#dc2626", "#f59e0b", "#111827")
  cols <- rep(cols, length.out = length(rocs))
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old.par), add = TRUE)
  graphics::par(mar = c(5.2, 5.2, 4.2, 2.4), cex = 1.2, cex.lab = 1.25,
                cex.axis = 1.1, cex.main = 1.2, xpd = NA)
  first <- TRUE
  for (i in seq_along(rocs)) {
    r <- rocs[[i]]
    if (inherits(r, "roc")) {
      if (first) {
        graphics::plot(1 - r$specificities, r$sensitivities, type = "l",
                       col = cols[i], lwd = 2, xlab = "1 - Specificity",
                       ylab = "Sensitivity", main = main,
                       xlim = c(0, 1), ylim = c(0, 1), ...)
        first <- FALSE
      } else {
        graphics::lines(1 - r$specificities, r$sensitivities, col = cols[i], lwd = 2)
      }
    } else if (!is.null(r$FP) && !is.null(r$TP)) {
      if (first) {
        graphics::plot(r$FP, r$TP, type = "l", col = cols[i], lwd = 2,
                       xlab = "1 - Specificity", ylab = "Sensitivity",
                       main = main, xlim = c(0, 1), ylim = c(0, 1), ...)
        first <- FALSE
      } else {
        graphics::lines(r$FP, r$TP, col = cols[i], lwd = 2)
      }
    }
  }
  graphics::abline(0, 1, lty = 2, col = "grey65")
  if (!is.null(names(rocs))) {
    graphics::legend("bottomright", legend = names(rocs), col = cols,
                     lwd = 2, bty = "n", cex = 1.05)
  }
  invisible(rocs)
}

.metric_vector <- function(x, model = NULL) {
  if (is.null(dim(x))) return(as.numeric(stats::setNames(x, names(x))))
  mat <- as.matrix(x)
  if (!is.null(model)) {
    if (!model %in% colnames(mat)) stop("model not found in metric matrix.")
    out <- mat[, model]
  } else {
    out <- rowMeans(mat, na.rm = TRUE)
  }
  names(out) <- rownames(mat) %||% names(out)
  out
}

.padded_range <- function(x, pad = 0.12, lower = NULL, upper = NULL) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(0, 1))
  rng <- range(x, na.rm = TRUE)
  span <- diff(rng)
  if (!is.finite(span) || span < 1e-12) span <- max(abs(rng[1]), 1) * 0.1
  out <- c(rng[1] - span * pad, rng[2] + span * pad)
  if (!is.null(lower)) out[1] <- max(lower, out[1])
  if (!is.null(upper)) out[2] <- min(upper, out[2])
  if (diff(out) < 1e-12) out <- out + c(-0.05, 0.05)
  out
}

.coerce_roc_objects <- function(x, outcome, outcome.type, survdays, time.point) {
  if (inherits(x, "roc") || (!is.null(x$FP) && !is.null(x$TP))) return(list(ROC = x))
  if (is.list(x) && !is.null(x$roc)) x <- x$roc
  if (is.list(x) && all(vapply(x, function(z) inherits(z, "roc") || (!is.null(z$FP) && !is.null(z$TP)), logical(1)))) {
    return(x)
  }
  if (is.list(x) && !is.null(x$predictions)) x <- x$predictions
  if (is.null(outcome)) stop("outcome is required when x contains prediction scores.")

  if (outcome.type == "binary") {
    if (is.null(dim(x))) {
      return(list(ROC = pROC::roc(outcome, as.numeric(x), ci = TRUE, direction = "<", quiet = TRUE)))
    }
    mat <- as.matrix(x)
    out <- lapply(seq_len(ncol(mat)), function(j) {
      pROC::roc(outcome, mat[, j], ci = TRUE, direction = "<", quiet = TRUE)
    })
    names(out) <- colnames(mat) %||% paste0("Model_", seq_len(ncol(mat)))
    return(out)
  }

  if (outcome.type == "categorical") {
    mat <- as.matrix(x)
    classes <- colnames(mat) %||% sort(unique(as.character(outcome)))
    if (ncol(mat) != length(classes)) {
      stop("Categorical ROC requires one prediction column per class.")
    }
    out <- lapply(seq_len(ncol(mat)), function(j) {
      pROC::roc(as.integer(as.character(outcome) == classes[j]), mat[, j],
                ci = TRUE, direction = "<", quiet = TRUE)
    })
    names(out) <- paste0(classes, " vs rest")
    return(out)
  }

  if (is.null(survdays) || is.null(time.point)) {
    stop("survdays and time.point are required for survival ROC prediction scores.")
  }
  if (is.null(dim(x))) {
    return(list(ROC = survivalROC::survivalROC(
      Stime = survdays, status = outcome, marker = as.numeric(x),
      predict.time = time.point, method = "KM"
    )))
  }
  mat <- as.matrix(x)
  out <- lapply(seq_len(ncol(mat)), function(j) {
    survivalROC::survivalROC(
      Stime = survdays, status = outcome, marker = mat[, j],
      predict.time = time.point, method = "KM"
    )
  })
  names(out) <- colnames(mat) %||% paste0("Model_", seq_len(ncol(mat)))
  out
}
