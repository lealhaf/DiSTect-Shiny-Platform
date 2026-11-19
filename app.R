# COMP 401 — DiSTect Spatial Transcriptomics Explorer, Version 3.4
# Observed-data exploration + read-only saved Stan fits.
# No Stan inference is executed by this application.

required_packages <- c("shiny", "ggplot2", "data.table")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages)) {
  stop("Install required packages: install.packages(c(",
       paste(sprintf('"%s"', missing_packages), collapse = ", "), "))")
}
library(shiny)
library(ggplot2)
library(data.table)
options(shiny.maxRequestSize = 100 * 1024^2)

# Works with shiny::runApp('app.R') from the project directory.
data_path <- "all_cleaned 1.tsv"
if (!file.exists(data_path)) stop("Dataset not found at: ", file.path(getwd(), data_path))
requested_genes <- c("ERBB2", "MUCL1", "TMSB10", "PTPRF", "MIEN1", "CISD3",
  "SCGB2A2", "PSMB3", "IGKC", "TPT1", "GAPDH", "SERF2", "GRB7", "ACTB",
  "RACK1", "FASN", "PHB", "EEF2", "MUC1", "KRT19")
metadata <- c("slice_id", "x", "y", "label")
all_columns <- names(fread(data_path, nrows = 0, check.names = FALSE))
if (!all(metadata %in% all_columns)) stop("Dataset is missing required metadata columns.")
genes <- intersect(requested_genes, all_columns)
if (!length(genes)) stop("No genes from the working panel were found.")
dat <- fread(data_path, select = c(metadata, genes), check.names = FALSE, showProgress = TRUE)
dat[, slice_id := as.character(slice_id)]
dat[, patient := substr(slice_id, 1L, 1L)]
if (anyNA(dat[, ..metadata]) || !all(dat$label %in% c(0L, 1L))) {
  stop("Missing metadata or unexpected cancer labels detected.")
}
default_dat <- copy(dat)
patients <- sort(unique(dat$patient))
patient_slices <- function(patient_id) sort(unique(dat[dat$patient == patient_id, slice_id]))

fit_catalog <- data.frame(
  file = c("DiSTect_A1_mean_VI.rds", "DiSTect_A1_mean_VI_strict.rds",
           "DiSTect_A1_scaled_VI.rds", "DiSTect_A1_NUTS_3genes.rds",
           "DiSTect_A1_NUTS_3genes_adapted.rds",
           "DiSTect_A1_NUTS_wide_spike.rds",
           "DiSTect_A1_simple_prior_NUTS.rds",
           "DiSTect_A1_verified_NUTS.rds"),
  method = c(rep("VI", 3), rep("NUTS", 5)),
  genes = c(rep("20-gene panel", 3), rep("ERBB2, MUCL1, TMSB10", 5)),
  prior = c(rep("Original spike-and-slab", 5),
            "Wider spike (diagnostic experiment)",
            "Normal(0,1) (diagnostic experiment)",
             "Original spike-and-slab"),
  preprocessing = c("Raw counts", "Raw counts", "Standardized genes",
                    rep("Standardized genes", 4), "log1p + z-score"),
  reported_global_pareto_k = c(6.47, 6.01, NA, rep(NA_real_, 5)),
  stringsAsFactors = FALSE
)

plot_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", colour = "#17324D", size = 15),
          plot.subtitle = element_text(colour = "#64748B"),
          panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold", colour = "#17324D"),
          legend.position = "bottom", plot.background = element_rect(fill = "white", colour = NA))
}
plot_cancer <- function(df) {
  df <- copy(df)
  df[, status := factor(label, levels = c(0, 1), labels = c("Non-cancer", "Cancer"))]
  ggplot(df, aes(x, y, color = status)) +
    geom_point(size = 2.5, alpha = .92) + facet_wrap(~slice_id) + coord_fixed() +
    scale_color_manual(values = c("Non-cancer" = "#3B82B8", "Cancer" = "#D65A60"), drop = FALSE) +
    labs(title = "Observed cancer status", subtitle = "Each point represents one tissue spot",
         x = "Spatial x", y = "Spatial y", color = "Spot label") + plot_theme()
}
plot_gene <- function(df, gene, transform) {
  df <- copy(df)
  values <- as.numeric(df[[gene]])
  if (transform == "log1p") values <- log1p(values)
  df[, expression := values]
  ggplot(df, aes(x, y, color = expression)) +
    geom_point(size = 2.5, alpha = .94) + facet_wrap(~slice_id) + coord_fixed() +
    scale_color_gradient(low = "#DDEEF1", high = "#006D77") +
    labs(title = paste("Spatial gene expression ·", gene),
         subtitle = if (transform == "log1p") "Observed log(count + 1)" else "Observed raw counts",
         x = "Spatial x", y = "Spatial y", color = "Expression") + plot_theme()
}

load_saved_fit <- function(filename) {
  if (!filename %in% fit_catalog$file) stop("Unknown fit selected.")
  path <- file.path("results", filename)
  if (!file.exists(path)) stop("Fit not found: ", path)
  if (!requireNamespace("rstan", quietly = TRUE)) stop("Install rstan to inspect saved Stan fits.")
  obj <- readRDS(path)
  if (inherits(obj, "stanfit")) {
    return(list(fit = obj, genes = NULL, provenance = "Gene order not verified"))
  }
  if (is.list(obj) && inherits(obj$fit, "stanfit")) {
    genes <- obj$genes
    n_beta <- as.integer(obj$fit@par_dims$beta[1])
    if (!is.character(genes) || length(genes) != n_beta ||
        anyNA(genes) || anyDuplicated(genes) || any(!nzchar(genes))) {
      stop("Recorded gene names do not match the fit's beta dimension.")
    }
    return(list(fit = obj$fit, genes = genes,
                provenance = "Gene order recorded with this experiment",
                record = obj))
  }
  stop("Unsupported saved-fit format.")
}
posterior_table <- function(fit, gene_names) {
  draws <- rstan::extract(fit, pars = "beta", permuted = TRUE, inc_warmup = FALSE)$beta
  if (is.null(draws)) stop("No beta draws available.")
  b <- as.matrix(draws)
  if (length(gene_names) != ncol(b)) gene_names <- paste0("beta[", seq_len(ncol(b)), "]")
  q <- apply(b, 2, quantile, probs = c(.025, .5, .975), na.rm = TRUE)
  data.frame(Gene = gene_names, Mean = colMeans(b, na.rm = TRUE),
             Lower95 = q[1, ], Median = q[2, ], Upper95 = q[3, ], check.names = FALSE)
}
model_diagnostics <- function(fit, method, global_k = NA_real_) {
  notes <- character()
  if (method == "VI") {
    notes <- c(notes, "VI approximates the posterior; NUTS chain diagnostics do not apply.")
    if (is.finite(global_k)) {
      notes <- c(notes, sprintf("Reported global VI Pareto-k: %.2f (historical console log; not extracted from RDS). This indicates an unreliable approximation.", global_k))
    } else {
      notes <- c(notes, "Global VI Pareto-k is unavailable for this fit. Parameter-level khat is not interchangeable with global Pareto-k.")
    }
  } else {
    pars <- tryCatch(rstan::get_sampler_params(fit, inc_warmup = FALSE), error = function(e) NULL)
    if (is.null(pars) || !length(pars) || !all(vapply(pars, function(z) "divergent__" %in% colnames(z), logical(1)))) {
      notes <- c(notes, "NUTS divergence diagnostics unavailable.")
    } else {
      div <- sum(vapply(pars, function(z) sum(z[, "divergent__"]), numeric(1)))
      notes <- c(notes, paste("NUTS divergent transitions:", div))
      if (div > 0) notes <- c(notes, "Divergences detected: posterior estimates require review.")
      if (all(vapply(pars, function(z) "treedepth__" %in% colnames(z), logical(1)))) {
        td <- vapply(pars, function(z) max(z[, "treedepth__"]), numeric(1))
        notes <- c(notes, paste("Maximum observed tree depth by chain:", paste(td, collapse = ", ")))
      }
    }
    summ <- tryCatch(rstan::summary(fit)$summary, error = function(e) NULL)
    if (!is.null(summ) && "Rhat" %in% colnames(summ)) {
      rh <- summ[, "Rhat"]
      rh <- rh[is.finite(rh)]
      if (length(rh)) {
        notes <- c(notes, sprintf("Maximum finite R-hat: %.3f", max(rh)))
        if (any(rh > 1.01)) notes <- c(notes, "Some R-hat values exceed 1.01.")
      } else notes <- c(notes, "Finite R-hat values unavailable.")
    }
  }
  paste(c(notes, "RESEARCH WARNING: This fit is not scientifically validated. Gene rankings and credible intervals are preliminary."), collapse = "\n")
}

metric_card <- function(label, value, note = NULL) {
  div(class = "metric-card", div(class = "metric-label", label),
      div(class = "metric-value", value), if (!is.null(note)) div(class = "metric-note", note))
}
section_head <- function(title, subtitle = NULL) {
  div(class = "section-heading", h2(title), if (!is.null(subtitle)) p(subtitle))
}
warning_box <- function(...) div(class = "warning-box", strong("Research status · preliminary inference"), p(...))

css <- "
body {background:#F5F8FC; color:#18324B; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;}
.container-fluid {max-width:1600px; padding:0 30px 36px;}
.app-header {background:linear-gradient(115deg,#113B55,#176B79); color:white; padding:26px 30px; margin:0 -30px 24px; border-radius:0 0 20px 20px; box-shadow:0 8px 24px rgba(17,59,85,.12);}
.app-header h1 {margin:0 0 6px; font-size:29px; font-weight:750; letter-spacing:-.5px;}
.app-header p {margin:0; color:#D5EBF1; font-size:14px;}
.header-tag {float:right; border:1px solid rgba(255,255,255,.4); border-radius:20px; padding:7px 13px; color:white; font-size:12px;}
.control-card,.panel-card,.metric-card {background:white; border:1px solid #E1E9F1; border-radius:16px; box-shadow:0 3px 14px rgba(18,48,72,.045);}
.control-card {padding:22px; position:sticky; top:16px;}
.control-card h3 {font-size:17px; margin:0 0 18px; font-weight:750;}
.panel-card {padding:24px; margin-bottom:18px;}
.metric-grid {display:grid; grid-template-columns:repeat(4,minmax(0,1fr)); gap:14px; margin-bottom:22px;}
.metric-card {padding:20px 18px;}
.metric-label {font-size:12px; text-transform:uppercase; letter-spacing:.8px; color:#64788C; font-weight:700;}
.metric-value {font-size:30px; color:#164F62; font-weight:750; margin-top:6px;}
.metric-note {font-size:12px; color:#7A8DA0; margin-top:4px;}
.section-heading h2 {font-size:23px; font-weight:750; margin:4px 0 7px; color:#17324D;}
.section-heading p {color:#61778B; margin-bottom:20px;}
.nav-tabs {border-bottom:1px solid #DCE6EE; margin-bottom:20px; display:flex; flex-wrap:wrap; gap:5px;}
.nav-tabs>li>a {border:0!important; border-radius:10px 10px 0 0!important; color:#536B7F; font-weight:650; padding:13px 16px;}
.nav-tabs>li.active>a,.nav-tabs>li.active>a:focus,.nav-tabs>li.active>a:hover {color:#006D77!important; background:#E8F5F5!important; border-bottom:3px solid #008C95!important;}
.form-control,.selectize-input {border-radius:9px!important; border-color:#CFDDE7!important; box-shadow:none!important;}
.control-label {color:#324C62; font-size:13px; font-weight:700;}
.btn {border-radius:9px; font-weight:650;}
.btn-default {background:#EDF5F7; border-color:#C7E0E5; color:#14586B;}
.btn-default:hover {background:#DDEFF1; border-color:#8EC3CA; color:#0B5261;}
.warning-box {border-left:4px solid #D99735; background:#FFF9EE; border-radius:10px; padding:15px 18px; margin:14px 0 20px; color:#72511F;}
.warning-box p {margin:6px 0 0;}
.info-box {background:#EDF7F8; border:1px solid #D3EAED; padding:17px; border-radius:12px; margin-bottom:16px; color:#245669;}
.muted {color:#6D8192; font-size:13px;}
hr {border-color:#E6EDF3;}
pre {background:#F5F8FB; border:1px solid #E5ECF3; border-radius:10px; color:#355168; white-space:pre-wrap; font-size:12px;}
.table>thead>tr>th {background:#F3F8FA; color:#36576B; border-bottom:1px solid #DDE9EF;}
.shiny-plot-output {max-width:100%;}
@media(max-width:1050px) {.metric-grid{grid-template-columns:repeat(2,1fr);} .control-card{position:static;}}
@media(max-width:650px) {.container-fluid{padding:0 12px 20px;} .app-header{margin:0 -12px 16px;padding:20px;} .header-tag{float:none;display:inline-block;margin-bottom:12px;} .metric-grid{grid-template-columns:1fr 1fr;} .metric-value{font-size:24px;} .panel-card{padding:15px;}}
"

ui <- fluidPage(
  tags$head(tags$style(HTML(paste0(css, "
.selection-metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(130px,1fr));gap:12px;margin:18px 0 26px;} .selection-metric{background:#F4F9FB;border:1px solid #DDEBF0;border-radius:12px;padding:16px;display:flex;flex-direction:column;gap:8px;} .selection-metric span{font-size:12px;color:#597287;font-weight:600;} .selection-metric strong{font-size:25px;color:#135B70;} .overview-visuals{display:flex;flex-wrap:wrap;gap:24px;}")))),
  div(class = "app-header",
      span(class = "header-tag", "COMP 401 · Research prototype"),
      h1("DiSTect", tags$span(style = "font-weight:350;opacity:.85;", " / Spatial Transcriptomics Explorer")),
      p("HER2+ breast cancer · Tissue-level spatial exploration and experimental Bayesian analysis")),
  fluidRow(
    column(3,
      div(class = "control-card",
        h3("Explore tissue data"),
        radioButtons("data_source", "Data source", choices = c("Built-in HER2+ dataset" = "built_in", "Upload my CSV/TSV" = "upload"), selected = "built_in"),
        conditionalPanel("input.data_source == 'upload'",
          fileInput("user_data", "Upload tissue data", accept = c(".csv", ".tsv", ".txt")),
          p(class = "muted", "Required: slice_id, x, y, label (0/1), plus at least one numeric gene column. CSV or tab-separated; max 100 MB."),
          verbatimTextOutput("upload_status")
        ),
        selectInput("patient", "Patient", choices = patients, selected = if ("A" %in% patients) "A" else patients[1]),
        selectizeInput("slices", "Tissue slices", choices = patient_slices(if ("A" %in% patients) "A" else patients[1]),
                       selected = head(patient_slices(if ("A" %in% patients) "A" else patients[1]), 2), multiple = TRUE),
        selectizeInput("gene", "Gene", choices = genes, selected = if ("ERBB2" %in% genes) "ERBB2" else genes[1]),
        radioButtons("transform", "Expression scale",
                     choices = c("log(count + 1)" = "log1p", "Raw counts" = "raw"), selected = "log1p"),
        hr(), strong("Current selection"), br(), br(), textOutput("selection_summary"),
        br(), downloadButton("export_csv", "Download spot data", class = "btn-block"),
        br(), div(class = "muted", "These controls affect observed-data plots only. Saved model fits are selected separately."))
    ),
    column(9,
      uiOutput("dataset_metrics"),
      tabsetPanel(id = "main_tabs",
        tabPanel("Overview", value = "overview",
          section_head("Research dashboard", "Explore observed spatial patterns and inspect saved experimental DiSTect fits."),
          div(class = "panel-card",
            h3("About this project"),
            p("DiSTect investigates associations between gene expression and spatially distributed disease labels. This application presents observed HER2+ breast cancer tissue data and separately displays posterior summaries from saved validation experiments."),
            div(class = "info-box", strong("What you can explore"),
                p("Choose a patient and multiple slices, compare observed cancer labels and gene-expression patterns, and inspect saved A1 Bayesian results. Selecting a gene does not refit the model.")),
            h4("Selected tissue overview"),
            uiOutput("overview_metrics"),
            div(class = "overview-visuals",
              div(style = "flex:1;min-width:260px;", h4("Cancer-label composition"), plotOutput("overview_composition", height = "245px")),
              div(style = "flex:1;min-width:260px;", h4("Cancer proportion by slice"), plotOutput("overview_proportions", height = "245px"))
            )
          ),
          warning_box("Saved Bayesian results are exploratory and have not passed an agreed scientific validation standard. No model fitting occurs in the dashboard.")
        ),
        tabPanel("Spatial Explorer", value = "explorer",
          section_head("Spatial tissue maps", "Switch between observed cancer status and expression. Click a point to inspect its spot data."),
          div(class = "panel-card",
            radioButtons("map_type", "Map layer", choices = c("Cancer labels" = "cancer", "Gene expression" = "gene"),
                         selected = "cancer", inline = TRUE),
            conditionalPanel("input.map_type == 'cancer'",
              plotOutput("cancer_plot", height = "490px", click = "cancer_click")),
            conditionalPanel("input.map_type == 'gene'",
              plotOutput("gene_plot", height = "490px", click = "gene_click")),
            div(class = "info-box", strong("Selected spot"), tableOutput("spot_details")),
            downloadButton("export_cancer", "Cancer map PNG"),
            downloadButton("export_gene", "Gene map PNG")
          )
        ),
        tabPanel("Slice Comparison", value = "comparison",
          section_head("Compare tissue slices", "Summary statistics use the current patient, slices, gene and expression scale."),
          div(class = "panel-card", h4("Slice-level summary"), tableOutput("comparison")),
          div(class = "panel-card", h4("Observed gene expression by slice"),
              plotOutput("comparison_plot", height = "330px"),
              p(class = "muted", "Mean expression uses the selected raw or log-transformed scale; values are descriptive, not fitted effects.")),
          div(class = "panel-card", h4("Interactive expression distributions"),
              radioButtons("distribution_type", "Display", choices = c("Violin and boxplot" = "violin", "Boxplot" = "box", "Jittered spots" = "jitter"), selected = "violin", inline = TRUE),
              checkboxInput("distribution_by_label", "Separate cancer and non-cancer spots", value = FALSE),
              plotOutput("distribution_plot", height = "420px", brush = brushOpts("distribution_brush")),
              verbatimTextOutput("distribution_selection"),
              downloadButton("export_distribution", "Download distribution PNG"),
              p(class = "muted", "Observed values only. Brush vertically to count spots within an expression range across selected slices; no Bayesian inference is run."))
        ),
        tabPanel("DiSTect Analysis", value = "analysis",
          section_head("Saved Bayesian experiments", "Inspect posterior summaries without executing Stan."),
          warning_box("A fit with zero divergences under an altered prior does not validate the original DiSTect model. Interpret all estimates as preliminary."),
          div(class = "panel-card",
            selectInput("fit_choice", "Saved A1 validation experiment",
                        choices = setNames(fit_catalog$file, paste(fit_catalog$method, fit_catalog$file, sep = " · ")),
                        selected = "DiSTect_A1_verified_NUTS.rds"),
            div(class = "info-box", strong("Fit provenance"), verbatimTextOutput("fit_metadata")),
            h4("Posterior gene effects · 95% credible intervals"),
            fluidRow(
              column(6, selectInput("coef_filter", "Show genes", choices = c("All genes" = "all", "Intervals excluding zero" = "exclude", "Intervals including zero" = "include"), selected = "all")),
              column(6, selectInput("coef_order", "Sort by", choices = c("Absolute effect" = "absolute", "Positive to negative" = "descending", "Negative to positive" = "ascending", "Gene name" = "gene"), selected = "absolute"))
            ),
            p(class = "muted", "Filtering describes the saved experimental posterior only; intervals excluding zero do not establish validated discoveries."),
            plotOutput("coef_plot", height = "400px"), tableOutput("coef_table"),
            downloadButton("export_coefficients", "Export coefficients CSV"),
            hr(), h4("Spatial effect and mixture weight"), verbatimTextOutput("spatial_summary")
          )
        ),
        tabPanel("Diagnostics", value = "diagnostics",
          section_head("Inference quality & limitations", "Diagnostics refer to the selected saved fit, not to the sidebar's observed-data selection."),
          warning_box("No saved fit has been declared scientifically validated. Unavailable diagnostics are not evidence of convergence."),
          div(class = "panel-card", h4("Selected fit diagnostics"), verbatimTextOutput("fit_diagnostics")),
          div(class = "panel-card", h4("Display validation checks"),
              p("These checks verify observed-data arithmetic and selected posterior-summary calculations; they cannot establish that the Bayesian model is scientifically correct."),
              verbatimTextOutput("validation_checks")),
          div(class = "panel-card", h4("How to interpret these checks"),
            tags$ul(
              tags$li(strong("VI Pareto-k:"), " large values can indicate an unreliable variational approximation; global and parameter-level diagnostics must not be confused."),
              tags$li(strong("NUTS divergences:"), " indicate potential problems exploring the posterior distribution."),
              tags$li(strong("R-hat:"), " values above approximately 1.01 warrant further review; R-hat alone cannot establish correctness."),
              tags$li(strong("Prior sensitivity:"), " improving diagnostics by changing the prior changes the statistical model.")
            )
          )
        ),
        tabPanel("Gene Network", value = "network",
          section_head("Observed gene co-expression network", "Explore pairwise expression correlations across the selected tissue spots."),
          div(class = "info-box", strong("Descriptive analysis, not inferred biological interactions"),
              p("Nodes are genes and edges represent Pearson or Spearman correlations between observed expression values. The network is independent of the saved Bayesian fits. Correlation does not establish causation, regulation, or direct molecular interaction.")),
          div(class = "panel-card",
              h4("Network settings"),
              selectizeInput("network_genes", "Genes to include (2–20)", choices = genes,
                             selected = head(genes, 8), multiple = TRUE),
              fluidRow(
                column(6, selectInput("network_method", "Correlation method",
                                      choices = c("Spearman (rank-based)" = "spearman", "Pearson (linear)" = "pearson"), selected = "spearman")),
                column(6, sliderInput("network_threshold", "Minimum absolute correlation |r|", min = 0, max = 1, value = 0.5, step = 0.05))
              ),
              checkboxInput("network_log", "Use log(count + 1) expression", value = TRUE),
              textOutput("network_status"),
              plotOutput("network_plot", height = "510px"),
              downloadButton("export_network_png", "Download network PNG"),
              downloadButton("export_network_csv", "Download edge table CSV")),
          div(class = "panel-card", h4("Gene-pair correlations"),
              tableOutput("network_table"),
              p(class = "muted", "Correlations are calculated across the current patient and selected slices. Pooled spots may be spatially dependent; this is an exploratory visualization without significance testing or multiple-testing claims."))
        ),
        tabPanel("Methods & About", value = "methods",
          section_head("Methods & definitions", "A concise guide to the data, model parameters and current limitations."),
          div(class = "panel-card",
            h4("Bring your own data"),
            p("Select Upload my CSV/TSV in the sidebar to explore another spatial dataset. Required columns: slice_id, x, y, label (binary 0/1), and one or more numeric gene columns. The app derives a patient ID from the first character of slice_id; for other naming conventions, patient grouping may need customization. Uploaded data affect observed plots only, not the saved A1 fits."),
            h4("Observed dataset"),
            p("HER2+ breast cancer spatial transcriptomics data with tissue coordinates (x, y), binary observed disease labels and gene-expression measurements. The application loads only a 20-gene working panel for responsiveness."),
            h4("Expression transformation"),
            p("The spatial explorer can display raw gene counts or log(count + 1). Changing this display setting does not alter saved Bayesian model fits."),
            h4("Spatial neighborhood"),
            p("The validation implementation uses tissue spots within Euclidean distance ≤ 1, excluding the spot itself. In the multi-slice code, neighbors are restricted to the same slice."),
            h4("Model parameters"),
            tags$ul(tags$li(strong("β:"), " gene-expression coefficients."),
                    tags$li(strong("η:"), " spatial-neighbor association coefficient."),
                    tags$li(strong("w:"), " spike-and-slab mixture weight."),
                    tags$li(strong("U:"), " experimental slice-specific random effects in the multi-slice model.")),
            h4("Important limitations"),
            p("All currently catalogued fits are A1-only experiments. Multi-slice heterogeneity, spot-level fitted probabilities and Bayesian interaction networks are not implemented or validated in this version. The co-expression network is a separate, descriptive analysis of observed gene correlations. The current model's use of observed neighboring labels should not be confused with out-of-sample predictive validation."),
            p(class = "muted", "COMP 401 research prototype · Read-only saved Stan results")
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  uploaded_data <- reactive({
    req(input$user_data)
    f <- input$user_data
    if (is.na(f$size) || f$size > 100 * 1024^2) stop("File exceeds 100 MB limit.")
    ext <- tolower(tools::file_ext(f$name))
    if (!ext %in% c("csv", "tsv", "txt")) stop("Use a CSV or TSV file.")
    sep <- if (ext == "csv") "," else "\t"
    dt <- data.table::fread(f$datapath, sep = sep, check.names = FALSE, showProgress = FALSE)
    if (!all(metadata %in% names(dt))) stop("Missing required columns: ", paste(setdiff(metadata, names(dt)), collapse = ", "))
    if (nrow(dt) == 0L) stop("The uploaded dataset is empty.")
    if (anyDuplicated(names(dt))) stop("Duplicate column names are not supported.")
    dt[, slice_id := as.character(slice_id)]
    if (anyNA(dt[, ..metadata]) || any(!nzchar(dt$slice_id))) stop("Missing slice IDs, coordinates, or labels.")
    if (!is.numeric(dt$x) || !is.numeric(dt$y) || any(!is.finite(dt$x)) || any(!is.finite(dt$y))) stop("Coordinates x and y must be finite numeric values.")
    if (!is.numeric(dt$label) || !all(dt$label %in% c(0, 1))) stop("Cancer label must be numeric 0 or 1.")
    dt[, label := as.integer(label)]
    gene_cols <- setdiff(names(dt), c(metadata, "patient"))
    gene_cols <- gene_cols[vapply(dt[, ..gene_cols], is.numeric, logical(1))]
    gene_cols <- gene_cols[vapply(dt[, ..gene_cols], function(v) all(is.finite(v)) && all(v >= 0), logical(1))]
    if (!length(gene_cols)) stop("At least one finite, nonnegative numeric gene column is required.")
    dt <- dt[, c(metadata, gene_cols), with = FALSE]
    dt[, patient := substr(slice_id, 1L, 1L)]
    dt
  })
  active_data <- reactive({
    if (identical(input$data_source, "upload")) {
      tryCatch(uploaded_data(), error = function(e) e)
    } else default_dat
  })
  output$upload_status <- renderText({
    req(input$data_source == "upload")
    if (is.null(input$user_data)) return("Choose a file to begin.")
    d <- active_data()
    if (inherits(d, "error")) return(paste("Upload error:", conditionMessage(d)))
    sprintf("Loaded %s spots, %s slices and %s genes (observed data only).", format(nrow(d), big.mark = ","), length(unique(d$slice_id)), ncol(d) - 5L)
  })
  observeEvent(active_data(), {
    d <- active_data()
    if (inherits(d, "error")) return()
    ps <- sort(unique(d$patient))
    gs <- setdiff(names(d), c(metadata, "patient"))
    updateSelectInput(session, "patient", choices = ps, selected = ps[1])
    updateSelectizeInput(session, "slices", choices = sort(unique(d[d$patient == ps[1], slice_id])), selected = head(sort(unique(d[d$patient == ps[1], slice_id])), 2), server = FALSE)
    updateSelectizeInput(session, "gene", choices = gs, selected = if ("ERBB2" %in% gs) "ERBB2" else gs[1], server = TRUE)
    updateSelectizeInput(session, "network_genes", choices = gs, selected = head(gs, min(8L, length(gs))), server = TRUE)
  })
  output$dataset_metrics <- renderUI({
    d <- active_data()
    if (inherits(d, "error")) return(div(class = "info-box", paste("Upload error:", conditionMessage(d))))
    div(class = "metric-grid",
        metric_card("Patients", length(unique(d$patient)), "Derived from slice IDs"),
        metric_card("Tissue slices", length(unique(d$slice_id)), "In selected dataset"),
        metric_card("Spatial spots", format(nrow(d), big.mark = ","), "Observed tissue locations"),
        metric_card("Available genes", ncol(d) - 5L, "Numeric gene columns"))
  })
  active_slices <- function(patient_id, d) sort(unique(d[d$patient == patient_id, slice_id]))
  observeEvent(input$patient, {
    d <- active_data()
    if (inherits(d, "error")) return()
    choices <- active_slices(input$patient, d)
    updateSelectizeInput(session, "slices", choices = choices, selected = head(choices, 2), server = FALSE)
  }, ignoreInit = TRUE)

  selected <- reactive({
    req(input$patient, input$slices)
    d <- active_data()
    validate(need(!inherits(d, "error"), if (inherits(d, "error")) conditionMessage(d) else ""))
    valid <- intersect(input$slices, active_slices(input$patient, d))
    validate(need(length(valid) > 0, "Select at least one tissue slice."))
    d[d$patient == input$patient & d$slice_id %in% valid]
  })
  output$selection_summary <- renderText({
    df <- selected()
    sprintf("%s · %s · %s spots", input$patient,
            paste(sort(unique(df$slice_id)), collapse = ", "), format(nrow(df), big.mark = ","))
  })
  output$overview_metrics <- renderUI({
    df <- selected()
    pct <- round(100 * mean(df$label == 1L), 1)
    div(class = "selection-metrics",
      div(class = "selection-metric", span("Selected slices"), strong(length(unique(df$slice_id)))),
      div(class = "selection-metric", span("Observed spots"), strong(format(nrow(df), big.mark = ","))),
      div(class = "selection-metric", span("Cancer spots"), strong(format(sum(df$label == 1L), big.mark = ","))),
      div(class = "selection-metric", span("Non-cancer spots"), strong(format(sum(df$label == 0L), big.mark = ","))),
      div(class = "selection-metric", span("Cancer proportion"), strong(paste0(pct, "%")))
    )
  })
  output$overview_composition <- renderPlot({
    df <- selected()
    counts <- data.frame(Status = factor(c("Cancer", "Non-cancer"), levels = c("Cancer", "Non-cancer")),
      Spots = c(sum(df$label == 1L), sum(df$label == 0L)))
    ggplot(counts, aes(x = "", y = Spots, fill = Status)) +
      geom_col(width = 1, colour = "white", linewidth = 1) +
      coord_polar(theta = "y") +
      scale_fill_manual(values = c("Cancer" = "#D65A60", "Non-cancer" = "#3B82B8")) +
      theme_void(base_size = 12) + theme(legend.position = "bottom") + labs(fill = NULL)
  })
  output$overview_proportions <- renderPlot({
    df <- selected()[, .(Cancer_percent = 100 * mean(label == 1L)), by = slice_id]
    ggplot(df, aes(x = slice_id, y = Cancer_percent)) +
      geom_col(fill = "#148796", width = .62) +
      geom_text(aes(label = sprintf("%.1f%%", Cancer_percent)), vjust = -0.5, size = 3.5) +
      scale_y_continuous(limits = c(0, max(100, max(df$Cancer_percent) * 1.15)), expand = expansion(mult = c(0, .02))) +
      labs(x = "Tissue slice", y = "Cancer-labelled spots (%)") + plot_theme()
  })
  output$overview_table <- renderTable({
    df <- selected()
    data.frame(Measure = c("Selected slices", "Observed spots", "Cancer spots", "Non-cancer spots", "Cancer proportion"),
               Value = c(length(unique(df$slice_id)), nrow(df), sum(df$label == 1L),
                         sum(df$label == 0L), paste0(round(100 * mean(df$label == 1L), 1), "%")))
  }, striped = TRUE)

  output$cancer_plot <- renderPlot(plot_cancer(selected()))
  output$gene_plot <- renderPlot({req(input$gene); validate(need(input$gene %in% names(selected()), "Select a valid gene.")); plot_gene(selected(), input$gene, input$transform)})
  observeEvent(list(input$patient, input$slices, input$gene, input$transform, input$map_type), {
    # Invalidate the previous spot after the displayed data or map changes.
    spot_click_version(spot_click_version() + 1L)
  }, ignoreInit = TRUE)
  spot_click_version <- reactiveVal(0L)
  clicked_spot <- reactive({
    spot_click_version()
    req(input$map_type)
    click <- if (input$map_type == "cancer") input$cancer_click else input$gene_click
    req(click)
    # nearPoints handles ggplot facets; show only one nearest spot.
    df <- selected()
    nearPoints(df, click, xvar = "x", yvar = "y", threshold = 12, maxpoints = 1, addDist = FALSE)
  })
  output$spot_details <- renderTable({
    df <- clicked_spot()
    validate(need(nrow(df) > 0, "Click a tissue spot to see its coordinates, label and gene expression."))
    data.frame(Slice = df$slice_id, X = df$x, Y = df$y,
               Label = ifelse(df$label == 1L, "Cancer", "Non-cancer"),
               Gene = input$gene,
               Expression = round(if (input$transform == "log1p") log1p(df[[input$gene]]) else df[[input$gene]], 3),
               check.names = FALSE)
  })

  comparison_data <- reactive({
    req(input$gene)
    validate(need(input$gene %in% names(selected()), "Select a valid gene."))
    df <- copy(selected())
    values <- as.numeric(df[[input$gene]])
    if (input$transform == "log1p") values <- log1p(values)
    df[, expression := values]
    df[, .(Spots = .N, Cancer = sum(label == 1L), Non_cancer = sum(label == 0L),
           Cancer_percent = round(100 * mean(label == 1L), 1),
           Mean_selected_gene = round(mean(expression), 3)), by = .(Slice = slice_id)][order(Slice)]
  })
  output$comparison <- renderTable(comparison_data(), striped = TRUE)
  output$comparison_plot <- renderPlot({
    df <- comparison_data()
    ggplot(df, aes(x = Slice, y = Mean_selected_gene, fill = Slice)) +
      geom_col(width = .6, show.legend = FALSE) +
      scale_fill_manual(values = rep(c("#0F8290", "#74B8C0", "#33728E", "#A5CED3"), length.out = nrow(df))) +
      labs(x = "Tissue slice", y = "Mean selected-gene expression", title = input$gene) + plot_theme()
  })

  distribution_data <- reactive({
    req(input$gene)
    d <- copy(selected())
    validate(need(input$gene %in% names(d), "Choose a gene available in the selected dataset."))
    d[, expression := as.numeric(get(input$gene))]
    if (identical(input$transform, "log1p")) d[, expression := log1p(expression)]
    d[, status := factor(label, levels = c(0, 1), labels = c("Non-cancer", "Cancer"))]
    d[, slice_id := factor(slice_id, levels = sort(unique(slice_id)))]
    d
  })
  distribution_figure <- reactive({
    d <- distribution_data()
    p <- ggplot(d, aes(x = slice_id, y = expression))
    if (isTRUE(input$distribution_by_label)) p <- p + aes(fill = status, colour = status)
    if (identical(input$distribution_type, "box")) {
      p <- p + geom_boxplot(alpha = .65, outlier.alpha = .25)
    } else if (identical(input$distribution_type, "jitter")) {
      p <- p + geom_jitter(width = .16, height = 0, alpha = .35, size = 1.5)
    } else {
      p <- p + geom_violin(alpha = .3, trim = FALSE) + geom_boxplot(width = .16, alpha = .65, outlier.alpha = .15)
    }
    p + scale_fill_manual(values = c("Non-cancer" = "#3B82B8", "Cancer" = "#D65A60")) +
      scale_colour_manual(values = c("Non-cancer" = "#3B82B8", "Cancer" = "#D65A60")) +
      labs(title = paste("Distribution of", input$gene), x = "Tissue slice",
           y = if (input$transform == "log1p") "log(count + 1)" else "Raw count",
           fill = "Observed label", colour = "Observed label") + plot_theme()
  })
  output$distribution_plot <- renderPlot(distribution_figure())
  output$distribution_selection <- renderText({
    d <- distribution_data()
    brush <- input$distribution_brush
    if (is.null(brush)) return(sprintf("%s observed spots across %s selected slices. Brush the chart to inspect a region.", nrow(d), length(unique(d$slice_id))))
    within_range <- d$expression >= min(brush$ymin, brush$ymax) & d$expression <= max(brush$ymin, brush$ymax)
    sprintf("Spots within brushed expression range: %s (of %s total; all selected slices).", sum(within_range), nrow(d))
  })
  output$export_distribution <- downloadHandler(
    filename = function() paste0("DiSTect_distribution_", input$gene, ".png"),
    content = function(file) ggsave(file, distribution_figure(), width = 10, height = 6, dpi = 200))
  output$validation_checks <- renderText({
    d <- selected()
    lines <- c(
      sprintf("PASS: %s cancer + %s non-cancer = %s observed spots.", sum(d$label == 1), sum(d$label == 0), nrow(d)),
      sprintf("PASS: displayed cancer proportion calculated as %.2f%%.", 100 * mean(d$label == 1)),
      if (all(is.finite(d$x), is.finite(d$y))) "PASS: selected spatial coordinates are finite." else "FAIL: invalid coordinates.",
      "GENE PROVENANCE: The new A1 verified-order record stores its gene mapping; older fits remain unverified.",
      "LIMITATION: Original model/preprocessing provenance is incomplete for historical RDS fits.",
      "NOT VALIDATED: posterior convergence or scientific interpretation."
    )
    paste(lines, collapse = "\n")
  })

  selected_fit_info <- reactive({
    req(input$fit_choice)
    fit_catalog[match(input$fit_choice, fit_catalog$file), , drop = FALSE]
  })
  selected_fit <- reactive({
    req(input$fit_choice)
    tryCatch(load_saved_fit(input$fit_choice), error = function(e) e)
  })
  fit_coefficients <- reactive({
    fit <- selected_fit()
    validate(need(!inherits(fit, "error"), if (inherits(fit, "error")) conditionMessage(fit) else ""))
    info <- selected_fit_info()
    gene_names <- if (!is.null(fit$genes)) fit$genes else
      paste0("beta[", seq_len(as.integer(fit$fit@par_dims$beta[1])), "]")
    tryCatch(posterior_table(fit$fit, gene_names), error = function(e) e)
  })
  output$fit_metadata <- renderText({
    info <- selected_fit_info()
    loaded <- selected_fit()
    provenance <- if (inherits(loaded, "error")) conditionMessage(loaded) else loaded$provenance
    paste0("Tissue: A1 only\nMethod: ", info$method, "\nGenes: ", info$genes,
           "\nGene provenance: ", provenance,
           "\nPrior: ", info$prior, "\nPreprocessing: ", info$preprocessing,
           if (!file.exists(file.path("results", info$file))) paste0("\nMissing file: results/", info$file) else "",
           "\n\nSidebar selections and uploaded datasets do not change this saved fit.")
  })
  output$fit_diagnostics <- renderText({
    fit <- selected_fit()
    if (inherits(fit, "error")) return(conditionMessage(fit))
    info <- selected_fit_info()
    tryCatch(model_diagnostics(fit$fit, info$method, info$reported_global_pareto_k),
             error = function(e) paste("Could not extract diagnostics:", conditionMessage(e)))
  })
  filtered_coefficients <- reactive({
    tab <- fit_coefficients()
    validate(need(!inherits(tab, "error"), if (inherits(tab, "error")) conditionMessage(tab) else ""))
    if (identical(input$coef_filter, "exclude")) tab <- tab[tab$Lower95 > 0 | tab$Upper95 < 0, , drop = FALSE]
    if (identical(input$coef_filter, "include")) tab <- tab[tab$Lower95 <= 0 & tab$Upper95 >= 0, , drop = FALSE]
    switch(input$coef_order,
      descending = tab <- tab[order(-tab$Mean), , drop = FALSE],
      ascending = tab <- tab[order(tab$Mean), , drop = FALSE],
      gene = tab <- tab[order(tab$Gene), , drop = FALSE],
      tab <- tab[order(-abs(tab$Mean)), , drop = FALSE])
    tab
  })
  output$coef_table <- renderTable({
    tab <- filtered_coefficients()
    validate(need(nrow(tab) > 0, "No genes match this filter for the selected saved fit."))
    tab[, -1] <- lapply(tab[, -1, drop = FALSE], round, 3)
    tab
  }, striped = TRUE)
  output$coef_plot <- renderPlot({
    tab <- filtered_coefficients()
    validate(need(nrow(tab) > 0, "No genes match this filter for the selected saved fit."))
    ggplot(tab, aes(x = reorder(Gene, Mean), y = Mean)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "#94A3B8") +
      geom_pointrange(aes(ymin = Lower95, ymax = Upper95), colour = "#13798B") +
      coord_flip() + labs(x = "Gene", y = "Posterior coefficient (95% credible interval)",
                          title = "A1 posterior gene effects · preliminary") + plot_theme()
  })
  output$spatial_summary <- renderText({
    fit <- selected_fit()
    if (inherits(fit, "error")) return(conditionMessage(fit))
    draws <- tryCatch(rstan::extract(fit$fit, pars = c("eta", "w"), permuted = TRUE), error = function(e) list())
    result <- character()
    for (nm in c("eta", "w")) {
      if (!is.null(draws[[nm]])) {
        v <- as.numeric(draws[[nm]])
        result <- c(result, sprintf("%s: mean = %.3f; 95%% credible interval = [%.3f, %.3f]",
                                    nm, mean(v), quantile(v, .025), quantile(v, .975)))
      }
    }
    paste(c(result, "These estimates are conditional on the selected experimental fit."), collapse = "\n")
  })
  # Descriptive co-expression network: uses only the selected observed spots.
  network_result <- reactive({
    df <- selected()
    gs <- intersect(input$network_genes, setdiff(names(df), c(metadata, "patient")))
    validate(need(length(gs) >= 2L, "Select at least two genes."),
             need(length(gs) <= 20L, "Select no more than 20 genes."),
             need(nrow(df) >= 3L, "At least three tissue spots are required."))
    mat <- as.matrix(df[, ..gs])
    storage.mode(mat) <- "double"
    validate(need(all(is.finite(mat)), "Selected genes contain missing or non-finite values."),
             need(all(mat >= 0), "Expression values must be nonnegative for log transformation."))
    if (isTRUE(input$network_log)) mat <- log1p(mat)
    variable <- apply(mat, 2, function(z) length(unique(z)) > 1L)
    # Constant genes cannot have defined correlations; retain them as isolated nodes.
    correlations <- matrix(NA_real_, length(gs), length(gs), dimnames = list(gs, gs))
    if (sum(variable) >= 2L) {
      correlations[variable, variable] <- stats::cor(mat[, variable, drop = FALSE],
        method = input$network_method, use = "everything")
    }
    pairs <- which(upper.tri(correlations) & is.finite(correlations) &
                     abs(correlations) >= input$network_threshold, arr.ind = TRUE)
    edges <- data.frame(Gene1 = gs[pairs[, 1]], Gene2 = gs[pairs[, 2]],
                        Correlation = correlations[pairs], stringsAsFactors = FALSE)
    if (nrow(edges)) edges <- edges[order(-abs(edges$Correlation)), , drop = FALSE]
    list(genes = gs, edges = edges, spots = nrow(df), constant = gs[!variable])
  })
  output$network_status <- renderText({
    net <- network_result()
    paste0(length(net$genes), " genes · ", net$spots, " spots · ", nrow(net$edges),
           " connections above threshold",
           if (length(net$constant)) paste0(" · constant genes (no defined correlation): ", paste(net$constant, collapse = ", ")) else "")
  })
  network_figure <- reactive({
    net <- network_result()
    n <- length(net$genes)
    angle <- seq(0, 2 * pi, length.out = n + 1L)[-(n + 1L)] + pi / 2
    nodes <- data.frame(Gene = net$genes, x = cos(angle), y = sin(angle))
    p <- ggplot() + coord_fixed(xlim = c(-1.35, 1.35), ylim = c(-1.28, 1.28), clip = "off") +
      theme_void(base_size = 12) + theme(plot.margin = margin(18, 24, 18, 24), legend.position = "bottom")
    if (nrow(net$edges)) {
      ed <- merge(net$edges, nodes, by.x = "Gene1", by.y = "Gene", sort = FALSE)
      ed <- merge(ed, nodes, by.x = "Gene2", by.y = "Gene", suffixes = c("1", "2"), sort = FALSE)
      ed$Direction <- ifelse(ed$Correlation >= 0, "Positive", "Negative")
      p <- p + geom_segment(data = ed,
        aes(x = x1, y = y1, xend = x2, yend = y2,
            linewidth = abs(Correlation), colour = Direction), alpha = 0.75) +
        scale_colour_manual(values = c("Positive" = "#138C91", "Negative" = "#D65A60")) +
        scale_linewidth(range = c(0.5, 3), guide = "none")
    }
    p + geom_point(data = nodes, aes(x = x, y = y), size = 12, colour = "#176B79", fill = "#E4F4F4", shape = 21, stroke = 1.6) +
      geom_text(data = nodes, aes(x = x, y = y, label = Gene), size = 3.2, fontface = "bold", colour = "#17324D") +
      labs(title = "Observed gene co-expression", subtitle = paste0(net$spots, " spots · ", input$network_method, " correlation · |r| ≥ ", input$network_threshold)) +
      theme(plot.title = element_text(face = "bold", colour = "#17324D"),
            plot.subtitle = element_text(colour = "#64748B"))
  })
  output$network_plot <- renderPlot(network_figure())
  output$network_table <- renderTable({
    edges <- network_result()$edges
    validate(need(nrow(edges) > 0, "No gene pairs meet the selected threshold. Try lowering it."))
    edges$Correlation <- round(edges$Correlation, 3)
    edges
  }, striped = TRUE)
  output$export_network_csv <- downloadHandler(
    filename = function() "DiSTect_coexpression_edges.csv",
    content = function(file) data.table::fwrite(network_result()$edges, file))
  output$export_network_png <- downloadHandler(
    filename = function() "DiSTect_coexpression_network.png",
    content = function(file) ggplot2::ggsave(file, network_figure(), width = 10, height = 7, dpi = 200))

  output$export_csv <- downloadHandler(
    filename = function() paste0("DiSTect_", input$patient, "_spots.csv"),
    content = function(file) fwrite(selected(), file))
  output$export_cancer <- downloadHandler(
    filename = function() paste0("DiSTect_cancer_", input$patient, ".png"),
    content = function(file) ggsave(file, plot_cancer(selected()), width = 10, height = 6, dpi = 200))
  output$export_gene <- downloadHandler(
    filename = function() paste0("DiSTect_", input$gene, "_", input$patient, ".png"),
    content = function(file) ggsave(file, plot_gene(selected(), input$gene, input$transform), width = 10, height = 6, dpi = 200))
  output$export_coefficients <- downloadHandler(
    filename = function() paste0("DiSTect_", sub("\\.rds$", "", input$fit_choice), "_coefficients.csv"),
    content = function(file) {
      tab <- filtered_coefficients()
      fwrite(tab, file)
    })
}
shinyApp(ui, server)
