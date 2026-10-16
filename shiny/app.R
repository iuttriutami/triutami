# =============================================================================
# SHINY APP: METAHEURISTIC-OPTIMIZED K-MEDOIDS CLUSTERING
# Methods: GA, PSO, CSA (Crow Search), GWO (Grey Wolf Optimizer)
# =============================================================================

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cluster)
library(clusterCrit)
library(GA)
library(pso)

# ─────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────
decode_medoid <- function(sol, k, n) {
  idx <- as.integer(round(sol))
  idx <- pmax(1L, pmin(n, idx))
  idx <- unique(idx)
  if (length(idx) < k) {
    extra <- setdiff(seq_len(n), idx)
    idx <- c(idx, sample(extra, k - length(idx)))
  }
  idx[seq_len(k)]
}

labels_of_medoids <- function(DMAT, medoid_idx)
  apply(DMAT[, medoid_idx, drop = FALSE], 1, which.min)

sil_of_medoids <- function(DMAT, DIST, medoid_idx) {
  labels <- labels_of_medoids(DMAT, medoid_idx)
  if (length(unique(labels)) < 2) return(list(sil = -1, labels = labels))
  sw <- silhouette(labels, DIST)
  list(sil = mean(sw[, 3]), labels = labels)
}

extra_indices <- function(XMAT, labels) {
  if (length(unique(labels)) < 2) return(c(DB = NA, CH = NA))
  tryCatch({
    cr <- intCriteria(XMAT, as.integer(labels),
                      c("Davies_Bouldin", "Calinski_Harabasz"))
    c(DB = cr$davies_bouldin, CH = cr$calinski_harabasz)
  }, error = function(e) c(DB = NA, CH = NA))
}

# Every candidate is decoded exactly once, and the medoid set is stored together with the
# Silhouette it produced. The decoder repairs duplicates at random, so decoding the final
# continuous solution a second time could return a different medoid set from the one that
# was evaluated; the tracker below prevents that (same rule as the benchmarking script).
make_tracker <- function(k, n, DMAT, DIST) {
  env <- new.env()
  env$best_sil <- -Inf; env$best_med <- NULL; env$n_eval <- 0L
  env$eval <- function(sol) {
    med <- decode_medoid(sol, k, n)
    s   <- sil_of_medoids(DMAT, DIST, med)$sil
    env$n_eval <- env$n_eval + 1L
    if (s > env$best_sil) { env$best_sil <- s; env$best_med <- med }
    s
  }
  env
}

# Check that the reported Silhouette is reproduced by the returned medoids.
finish_run <- function(tr, DMAT, DIST, conv = NULL, algo = "") {
  chk <- sil_of_medoids(DMAT, DIST, tr$best_med)$sil
  if (abs(chk - tr$best_sil) > 1e-9)
    stop(sprintf("%s: returned medoids give Silhouette %.10f, reported %.10f", algo, chk, tr$best_sil))
  list(med = tr$best_med, sil = tr$best_sil, conv = conv, n_eval = tr$n_eval)
}

vif_calc <- function(X) {
  X <- as.data.frame(X)
  sapply(seq_along(X), function(j) {
    r2 <- summary(lm(X[[j]] ~ ., data = X[, -j, drop = FALSE]))$r.squared
    1 / (1 - r2)
  })
}

# ─────────────────────────────────────────────────────────────────────────────
# METAHEURISTIC FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────
run_GA <- function(k, n, DMAT, DIST, seed = 42, popSize = 50, maxiter = 100,
                   trace = FALSE, progress_fn = NULL) {
  set.seed(seed)
  tr <- make_tracker(k, n, DMAT, DIST)
  g <- ga(
    type    = "real-valued",
    fitness = function(sol) tr$eval(sol),
    lower   = rep(1, k), upper = rep(n, k),
    popSize = popSize, maxiter = maxiter,
    pcrossover = 0.8, pmutation = 0.1, elitism = 2,
    monitor = if (!is.null(progress_fn)) {
      function(obj) {
        iter <- obj@iter
        progress_fn(iter, maxiter, paste0("GA generation ", iter, "/", maxiter))
      }
    } else FALSE,
    seed = seed
  )
  finish_run(tr, DMAT, DIST, conv = if (trace) cummax(as.numeric(g@summary[, "max"])) else NULL, algo = "GA")
}

run_PSO <- function(k, n, DMAT, DIST, seed = 42, s = 30, maxit = 100,
                    trace = FALSE, progress_fn = NULL) {
  set.seed(seed)
  tr <- make_tracker(k, n, DMAT, DIST)
  conv <- numeric(0)
  ctrl <- list(maxit = maxit, s = s, w = 0.729, c.p = 1.494, c.g = 1.494, trace = 0)
  psoptim(
    par    = sample(seq_len(n), k),
    fn     = function(sol) {
      v <- -tr$eval(sol)
      if (tr$n_eval %% s == 0) {                      # one swarm sweep finished
        conv <<- c(conv, tr$best_sil)
        it <- length(conv)
        if (!is.null(progress_fn)) progress_fn(min(it, maxit), maxit, paste0("PSO iteration ", min(it, maxit), "/", maxit))
      }
      v
    },
    lower  = rep(1, k), upper = rep(n, k),
    control = ctrl
  )
  finish_run(tr, DMAT, DIST, conv = if (trace) conv else NULL, algo = "PSO")
}

run_CSA <- function(k, n, DMAT, DIST, seed = 42, flock = 30, num_iter = 100,
                    fl = 2.0, AP = 0.1, trace = FALSE, progress_fn = NULL) {
  set.seed(seed)
  tr <- make_tracker(k, n, DMAT, DIST)
  memory <- pos <- vector("list", flock); sil_mem <- numeric(flock)
  for (i in seq_len(flock)) {                          # initial flock; memory = evaluated set
    m <- decode_medoid(sample(n, k), k, n)
    sil_mem[i] <- tr$eval(m); memory[[i]] <- m; pos[[i]] <- m
  }
  conv <- numeric(0)
  for (it in seq_len(num_iter)) {
    if (!is.null(progress_fn)) progress_fn(it, num_iter, paste0("CSA iteration ", it, "/", num_iter))
    for (i in seq_len(flock)) {
      j <- sample(setdiff(seq_len(flock), i), 1)
      if (runif(1) >= AP) {
        new_med <- integer(k)
        for (d in seq_len(k))
          new_med[d] <- if (runif(1) < fl / (fl + 1)) memory[[j]][sample(k, 1)] else pos[[i]][d]
      } else new_med <- sample(n, k)
      new_med <- decode_medoid(new_med, k, n)          # decoded once; tr$eval keeps it unchanged
      s <- tr$eval(new_med)
      pos[[i]] <- new_med
      if (s > sil_mem[i]) { memory[[i]] <- new_med; sil_mem[i] <- s }
    }
    conv <- c(conv, tr$best_sil)
  }
  finish_run(tr, DMAT, DIST, conv = if (trace) conv else NULL, algo = "CSA")
}

run_GWO <- function(k, n, DMAT, DIST, seed = 42, num_wolves = 30, num_iter = 100,
                    trace = FALSE, progress_fn = NULL) {
  set.seed(seed)
  tr  <- make_tracker(k, n, DMAT, DIST)
  pos <- matrix(runif(num_wolves * k, 1, n), nrow = num_wolves)
  fit <- apply(pos, 1, tr$eval)
  ord <- order(fit, decreasing = TRUE)
  alpha <- pos[ord[1], ]; beta <- pos[ord[2], ]; delta <- pos[ord[3], ]
  conv <- numeric(0)
  for (it in seq_len(num_iter)) {
    if (!is.null(progress_fn)) progress_fn(it, num_iter, paste0("GWO iteration ", it, "/", num_iter))
    a <- 2 - it * (2 / num_iter)
    for (i in seq_len(num_wolves)) {
      for (d in seq_len(k)) {
        X <- numeric(3)
        for (l in 1:3) {
          leader <- switch(l, alpha, beta, delta)
          A <- 2 * a * runif(1) - a; C <- 2 * runif(1)
          X[l] <- leader[d] - A * abs(C * leader[d] - pos[i, d])
        }
        pos[i, d] <- mean(X)
      }
      pos[i, ] <- pmin(pmax(pos[i, ], 1), n)
    }
    fit <- apply(pos, 1, tr$eval)
    ord <- order(fit, decreasing = TRUE)
    alpha <- pos[ord[1], ]; beta <- pos[ord[2], ]; delta <- pos[ord[3], ]
    conv <- c(conv, tr$best_sil)
  }
  finish_run(tr, DMAT, DIST, conv = if (trace) conv else NULL, algo = "GWO")
}

run_PAM <- function(k, data_scaled, DMAT, DIST) {
  set.seed(42)
  pam_res <- pam(data_scaled, k = k)
  med     <- pam_res$id.med
  list(med = med, sil = sil_of_medoids(DMAT, DIST, med)$sil)
}

# ─────────────────────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue", title = "Metaheuristic K-Medoids",
  dashboardHeader(title = span(icon("leaf"), " Metaheuristic K-Medoids"), titleWidth = 330),

  dashboardSidebar(width = 330,
    sidebarMenu(
      menuItem("📂 Upload & Preprocessing", tabName = "preproc", icon = icon("database")),
      menuItem("⚙️ Configure & Run",  tabName = "run",     icon = icon("cogs")),
      menuItem("📊 Clustering Results",       tabName = "results", icon = icon("chart-bar")),
      menuItem("📈 Convergence Curves",        tabName = "conv",    icon = icon("line-chart")),
      menuItem("🔍 Validation Analysis",        tabName = "valid",   icon = icon("check-circle")),
      menuItem("📋 Full Tables",                tabName = "tables",  icon = icon("table"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .skin-blue .main-header .logo { background-color: #1a5276; font-weight: bold; }
      .skin-blue .main-header .navbar { background-color: #1a5276; }
      .skin-blue .main-sidebar { background-color: #1c2833; }
      .skin-blue .sidebar-menu > li.active > a { background-color: #2e86c1; border-left-color: #f39c12; }
      .box-primary { border-top-color: #2e86c1; }
      .info-box-icon { min-height: 70px; line-height: 70px; }
      .info-box { min-height: 70px; }
      .info-box-content { padding-top: 10px; }
      body { font-family: 'Segoe UI', sans-serif; }
      .nav-tabs-custom > .nav-tabs > li.active > a { color: #2e86c1; border-top: 3px solid #2e86c1; }
      .progress-text { font-size: 13px; color: #555; margin-top: 4px; }
      .cluster-badge { display:inline-block; padding:3px 8px; border-radius:12px;
                       color:#fff; font-size:11px; font-weight:600; }
    "))),

    tabItems(

      # ── TAB 1: Upload & Preprocessing ─────────────────────────────────────
      tabItem(tabName = "preproc",
        fluidRow(
          box(title = "Upload Data", width = 6, status = "primary", solidHeader = TRUE,
            fileInput("file_csv", "Choose CSV file:", accept = ".csv"),
            checkboxInput("header_cb", "First row as header", value = TRUE),
            selectInput("sep_cb", "Separator:", choices = c("Comma ,"=",","Semicolon ;"=";","Tab"="\t"), selected = ","),
            uiOutput("col_name_ui"),
            uiOutput("col_drop_ui"),
            actionButton("btn_preprocess", "Process Data", icon = icon("play"), class = "btn-primary btn-lg")
          ),
          box(title = "Data Format Guide", width = 6, status = "info", solidHeader = TRUE,
            tags$ul(
              tags$li("The first column should be the ", tags$b("observation name"), " (province, city, etc.)"),
              tags$li("Numeric columns are used as clustering variables"),
              tags$li("Non-numeric columns are ignored automatically"),
              tags$li("Missing values (NA) are detected and reported"),
              tags$li("Data are standardized ", tags$b("z-score"), " before clustering")
            ),
            tags$hr(),
            tags$p("Example CSV format (EQI by province, Indonesia, 2023):"),
            tags$code("Province,AQI,WQI,LCQI,LQI"),
            tags$br(),
            tags$code("Aceh,72.3,68.1,80.2,74.5"),
            tags$br(), tags$code("...etc.")
          )
        ),
        fluidRow(
          box(title = "Raw Data Preview", width = 12, status = "primary",
            DT::dataTableOutput("tbl_raw"))
        ),
        fluidRow(
          box(title = "Mahalanobis Distance (Outlier Detection)", width = 6, status = "warning", solidHeader = TRUE,
            DT::dataTableOutput("tbl_mahal")),
          box(title = "VIF – Feature Redundancy", width = 6, status = "warning", solidHeader = TRUE,
            DT::dataTableOutput("tbl_vif"),
            tags$p(tags$b("Interpretation:"), " VIF < 5 = no serious redundancy", class = "text-muted")
          )
        )
      ),

      # ── TAB 2: Configure & Run ─────────────────────────────────────
      tabItem(tabName = "run",
        fluidRow(
          box(title = "Algorithm Configuration", width = 5, status = "primary", solidHeader = TRUE,
            sliderInput("k_val", "Number of clusters (k):", min = 2, max = 10, value = 3, step = 1),
            checkboxGroupInput("algo_sel", "Select algorithm:",
              choices  = c("PAM (baseline)"="PAM", "GA"="GA", "PSO"="PSO", "CSA"="CSA", "GWO"="GWO"),
              selected = c("PAM","GA","PSO","CSA","GWO")),
            numericInput("seed_val", "Random Seed:", value = 42, min = 1),
            tags$hr(),
            tags$h5("GA parameters"),
            sliderInput("ga_pop",  "Population size:", min = 10, max = 200, value = 50, step = 10),
            sliderInput("ga_iter", "Max. generations:",  min = 50, max = 500, value = 100, step = 50),
            tags$hr(),
            tags$h5("PSO/CSA/GWO parameters"),
            sliderInput("swarm", "Swarm/flock size:", min = 10, max = 100, value = 30, step = 10),
            sliderInput("gen_iter", "Max. iterations:",   min = 50, max = 500, value = 100, step = 50)
          ),
          box(title = "Analysis Mode", width = 4, status = "info", solidHeader = TRUE,
            radioButtons("mode_sel", "Select mode:",
              choices = c("Single k (fast)" = "single",
                          "Sweep k = 2..10 (comprehensive)" = "sweep"),
              selected = "single"),
            tags$p("Sweep mode runs GA for every k from 2 to 10.",
                   class = "text-muted"),
            tags$hr(),
            tags$h5("Medoid labels"),
            tags$p("The name column (first column of the data) is used to label the medoids.", class = "text-muted")
          ),
          box(title = "Run", width = 3, status = "success", solidHeader = TRUE,
            actionButton("btn_run", "▶ Start Clustering", icon = icon("rocket"),
                         class = "btn-success btn-lg", style = "width:100%; margin-bottom:10px"),
            tags$hr(),
            uiOutput("progress_ui"),
            verbatimTextOutput("log_out", placeholder = TRUE)
          )
        )
      ),

      # ── TAB 3: Clustering Results ────────────────────────────────────────────
      tabItem(tabName = "results",
        fluidRow(
          uiOutput("info_boxes")
        ),
        fluidRow(
          box(title = "Silhouette Plot", width = 6, status = "primary", solidHeader = TRUE,
            plotOutput("plot_sil", height = "350px")),
          box(title = "Silhouette Score Comparison", width = 6, status = "primary", solidHeader = TRUE,
            plotOutput("plot_compare", height = "350px"))
        ),
        fluidRow(
          box(title = "Cluster Membership (Best Algorithm)", width = 12, status = "success", solidHeader = TRUE,
            DT::dataTableOutput("tbl_member"))
        ),
        fluidRow(
          box(title = "Mean Profile per Cluster", width = 12, status = "warning", solidHeader = TRUE,
            plotOutput("plot_profile", height = "350px"))
        )
      ),

      # ── TAB 4: Convergence Curves ───────────────────────────────────────────
      tabItem(tabName = "conv",
        fluidRow(
          box(title = "Convergence Curve – Best Silhouette per Iteration", width = 12,
              status = "primary", solidHeader = TRUE,
            tags$p("Available only for the iterative algorithms (GA, PSO, CSA, GWO); PAM has no convergence curve.", class = "text-muted"),
            plotOutput("plot_conv", height = "420px"))
        )
      ),

      # ── TAB 5: Validation Analysis ───────────────────────────────────────────
      tabItem(tabName = "valid",
        fluidRow(
          box(title = "Cluster Validity Indices (final k)", width = 12, status = "primary", solidHeader = TRUE,
            DT::dataTableOutput("tbl_valid"),
            tags$p("Silhouette: higher is better | Davies–Bouldin: lower is better | Calinski–Harabasz: higher is better",
                   class = "text-muted")
          )
        ),
        fluidRow(
          box(title = "Index Sweep across k (GA)", width = 12, status = "info", solidHeader = TRUE,
            plotOutput("plot_sweep", height = "380px"))
        )
      ),

      # ── TAB 6: Full Tables ──────────────────────────────────────────────
      tabItem(tabName = "tables",
        tabBox(width = 12,
          tabPanel("Mahalanobis", DT::dataTableOutput("tbl_t1"), downloadButton("dl_t1", "Download CSV")),
          tabPanel("VIF", DT::dataTableOutput("tbl_t2"), downloadButton("dl_t2", "Download CSV")),
          tabPanel("Silhouette k=2..10", DT::dataTableOutput("tbl_t3"), downloadButton("dl_t3", "Download CSV")),
          tabPanel("Membership", DT::dataTableOutput("tbl_t4"), downloadButton("dl_t4", "Download CSV")),
          tabPanel("Cluster Profile", DT::dataTableOutput("tbl_t5"), downloadButton("dl_t5", "Download CSV")),
          tabPanel("Validation", DT::dataTableOutput("tbl_t6"), downloadButton("dl_t6", "Download CSV"))
        )
      )

    )  # end tabItems
  )
)

# ─────────────────────────────────────────────────────────────────────────────
# SERVER
# ─────────────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  rv <- reactiveValues(
    df           = NULL,
    data_num     = NULL,
    data_scaled  = NULL,
    DIST         = NULL,
    DMAT         = NULL,
    XMAT         = NULL,
    provinsi     = NULL,
    col_names    = NULL,
    mahal_tab    = NULL,
    vif_tab      = NULL,
    results      = NULL,
    sweep_tab    = NULL,
    conv_list    = NULL,
    valid_tab    = NULL,
    log_msgs     = character(0),
    running      = FALSE,
    prog_val     = 0,
    prog_msg     = ""
  )

  log_add <- function(msg) {
    rv$log_msgs <- c(rv$log_msgs, paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg))
  }

  # ── Upload & render raw table ─────────────────────────────────────────────
  observeEvent(input$file_csv, {
    req(input$file_csv)
    rv$df <- read.csv(input$file_csv$datapath,
                      header    = input$header_cb,
                      sep       = input$sep_cb,
                      check.names = FALSE)
    log_add(paste0("Data loaded: ", nrow(rv$df), " rows, ", ncol(rv$df), " columns"))
  })

  output$tbl_raw <- DT::renderDataTable({
    req(rv$df)
    DT::datatable(rv$df, options = list(scrollX = TRUE, pageLength = 8), rownames = FALSE)
  })

  output$col_name_ui <- renderUI({
    req(rv$df)
    selectInput("col_name", "Observation name/label column:",
                choices  = names(rv$df),
                selected = names(rv$df)[1])
  })

  output$col_drop_ui <- renderUI({
    req(rv$df)
    checkboxGroupInput("cols_drop", "Columns to exclude (optional):",
                       choices = names(rv$df)[-1])
  })

  # ── Preprocessing ─────────────────────────────────────────────────────────
  observeEvent(input$btn_preprocess, {
    req(rv$df, input$col_name)
    df <- rv$df

    # Label column
    prov_col <- input$col_name
    rv$provinsi <- df[[prov_col]]

    # Drop columns
    drop_cols <- c(prov_col, input$cols_drop)
    data_num  <- df[, !names(df) %in% drop_cols]
    data_num  <- data_num[, sapply(data_num, is.numeric), drop = FALSE]
    rv$data_num <- data_num
    rv$col_names <- colnames(data_num)
    log_add(paste0("Numeric columns: ", paste(rv$col_names, collapse = ", ")))

    # NA check
    na_count <- sum(is.na(data_num))
    if (na_count > 0) {
      log_add(paste0("Warning: ", na_count, " NA values found – incomplete rows removed"))
      keep <- complete.cases(data_num)
      rv$provinsi <- rv$provinsi[keep]
      data_num   <- data_num[keep, ]
      rv$data_num <- data_num
    }

    # Mahalanobis
    center <- colMeans(data_num); covm <- cov(data_num)
    mahal  <- mahalanobis(as.matrix(data_num), center, covm)
    cutoff <- qchisq(0.999, df = ncol(data_num))
    rv$mahal_tab <- data.frame(
      Observation      = rv$provinsi,
      Mahalanobis_D2 = round(mahal, 4),
      Cutoff         = round(cutoff, 4),
      Outlier        = mahal > cutoff
    )
    log_add(paste0("Outlier detection done. Outliers: ", sum(rv$mahal_tab$Outlier)))

    # VIF
    vif_vals <- vif_calc(data_num)
    rv$vif_tab <- data.frame(Variable = rv$col_names, VIF = round(vif_vals, 4))

    # Standardization
    data_scaled     <- scale(data_num)
    rownames(data_scaled) <- rv$provinsi
    rv$data_scaled  <- data_scaled
    rv$DIST         <- dist(data_scaled)
    rv$DMAT         <- as.matrix(rv$DIST)
    rv$XMAT         <- as.matrix(data_scaled)
    log_add("z-score standardization done. Ready for clustering.")
  })

  output$tbl_mahal <- DT::renderDataTable({
    req(rv$mahal_tab)
    DT::datatable(rv$mahal_tab, rownames = FALSE, options = list(pageLength = 8)) |>
      DT::formatStyle("Outlier",
        backgroundColor = DT::styleEqual(c(TRUE, FALSE), c("#f8d7da", "#d4edda")))
  })

  output$tbl_vif <- DT::renderDataTable({
    req(rv$vif_tab)
    DT::datatable(rv$vif_tab, rownames = FALSE, options = list(pageLength = 8)) |>
      DT::formatStyle("VIF",
        backgroundColor = DT::styleInterval(c(5, 10), c("#d4edda", "#fff3cd", "#f8d7da")))
  })

  # ── Run Clustering ────────────────────────────────────────────────────────
  observeEvent(input$btn_run, {
    req(rv$data_scaled, rv$DIST, rv$DMAT)

    n        <- nrow(rv$data_scaled)
    k        <- input$k_val
    algos    <- input$algo_sel
    seed_v   <- input$seed_val
    ga_pop   <- input$ga_pop
    ga_iter  <- input$ga_iter
    swarm_sz <- input$swarm
    gen_iter <- input$gen_iter

    rv$log_msgs <- character(0)
    rv$results  <- NULL
    rv$conv_list <- list()
    rv$running  <- TRUE

    withProgress(message = "Running clustering...", value = 0, {

      make_progress <- function(label) {
        function(cur, total, msg = "") {
          setProgress(value = cur / total, message = label, detail = msg)
        }
      }

      results   <- list()
      conv_list <- list()

      if ("PAM" %in% algos) {
        log_add("Running PAM (baseline)...")
        incProgress(0, message = "PAM baseline...")
        t0  <- proc.time()["elapsed"]
        res <- run_PAM(k, rv$data_scaled, rv$DMAT, rv$DIST)
        runtime <- as.numeric(proc.time()["elapsed"] - t0)
        lab <- labels_of_medoids(rv$DMAT, res$med)
        ei  <- extra_indices(rv$XMAT, lab)
        results[["PAM"]] <- list(med = res$med, sil = res$sil, labels = lab,
                                  DB = ei["DB"], CH = ei["CH"], runtime = runtime)
        log_add(sprintf("PAM done | Silhouette = %.4f | %.2fs", res$sil, runtime))
      }

      if ("GA" %in% algos) {
        log_add("Running GA...")
        t0  <- proc.time()["elapsed"]
        res <- run_GA(k, n, rv$DMAT, rv$DIST, seed = seed_v,
                      popSize = ga_pop, maxiter = ga_iter, trace = TRUE,
                      progress_fn = make_progress("GA"))
        runtime <- as.numeric(proc.time()["elapsed"] - t0)
        lab <- labels_of_medoids(rv$DMAT, res$med)
        ei  <- extra_indices(rv$XMAT, lab)
        results[["GA"]] <- list(med = res$med, sil = res$sil, labels = lab,
                                 DB = ei["DB"], CH = ei["CH"], runtime = runtime)
        conv_list[["GA"]] <- res$conv
        log_add(sprintf("GA done | Silhouette = %.4f | %.2fs", res$sil, runtime))
      }

      if ("PSO" %in% algos) {
        log_add("Running PSO...")
        t0  <- proc.time()["elapsed"]
        res <- run_PSO(k, n, rv$DMAT, rv$DIST, seed = seed_v,
                       s = swarm_sz, maxit = gen_iter, trace = TRUE,
                       progress_fn = make_progress("PSO"))
        runtime <- as.numeric(proc.time()["elapsed"] - t0)
        lab <- labels_of_medoids(rv$DMAT, res$med)
        ei  <- extra_indices(rv$XMAT, lab)
        results[["PSO"]] <- list(med = res$med, sil = res$sil, labels = lab,
                                  DB = ei["DB"], CH = ei["CH"], runtime = runtime)
        conv_list[["PSO"]] <- res$conv
        log_add(sprintf("PSO done | Silhouette = %.4f | %.2fs", res$sil, runtime))
      }

      if ("CSA" %in% algos) {
        log_add("Running CSA (Crow Search)...")
        t0  <- proc.time()["elapsed"]
        res <- run_CSA(k, n, rv$DMAT, rv$DIST, seed = seed_v,
                       flock = swarm_sz, num_iter = gen_iter, trace = TRUE,
                       progress_fn = make_progress("CSA"))
        runtime <- as.numeric(proc.time()["elapsed"] - t0)
        lab <- labels_of_medoids(rv$DMAT, res$med)
        ei  <- extra_indices(rv$XMAT, lab)
        results[["CSA"]] <- list(med = res$med, sil = res$sil, labels = lab,
                                  DB = ei["DB"], CH = ei["CH"], runtime = runtime)
        conv_list[["CSA"]] <- res$conv
        log_add(sprintf("CSA done | Silhouette = %.4f | %.2fs", res$sil, runtime))
      }

      if ("GWO" %in% algos) {
        log_add("Running GWO (Grey Wolf Optimizer)...")
        t0  <- proc.time()["elapsed"]
        res <- run_GWO(k, n, rv$DMAT, rv$DIST, seed = seed_v,
                       num_wolves = swarm_sz, num_iter = gen_iter, trace = TRUE,
                       progress_fn = make_progress("GWO"))
        runtime <- as.numeric(proc.time()["elapsed"] - t0)
        lab <- labels_of_medoids(rv$DMAT, res$med)
        ei  <- extra_indices(rv$XMAT, lab)
        results[["GWO"]] <- list(med = res$med, sil = res$sil, labels = lab,
                                  DB = ei["DB"], CH = ei["CH"], runtime = runtime)
        conv_list[["GWO"]] <- res$conv
        log_add(sprintf("GWO done | Silhouette = %.4f | %.2fs", res$sil, runtime))
      }

      rv$results   <- results
      rv$conv_list <- conv_list

      # Validation table
      val_df <- data.frame(
        Algorithm  = names(results),
        Silhouette = sapply(results, function(r) round(r$sil, 4)),
        DBI        = sapply(results, function(r) round(r$DB, 4)),
        CH         = sapply(results, function(r) round(r$CH, 2)),
        Runtime_s  = sapply(results, function(r) round(r$runtime, 3))
      )
      rv$valid_tab <- val_df
      log_add("All algorithms finished!")

      # Sweep k = 2..10 with GA in sweep mode
      if (input$mode_sel == "sweep" && "GA" %in% algos) {
        log_add("Sweep k = 2..10 with GA...")
        sweep_res <- data.frame(k = 2:10, Silhouette = NA, DBI = NA, CH = NA)
        for (i in seq_along(2:10)) {
          ki <- (2:10)[i]
          setProgress(i/9, message = paste0("Sweep GA k=", ki, "..."))
          r <- run_GA(ki, n, rv$DMAT, rv$DIST, seed = seed_v,
                      popSize = ga_pop, maxiter = min(ga_iter, 100))
          lab <- labels_of_medoids(rv$DMAT, r$med)
          ei  <- extra_indices(rv$XMAT, lab)
          s   <- sil_of_medoids(rv$DMAT, rv$DIST, r$med)$sil
          sweep_res$Silhouette[i] <- round(s, 4)
          sweep_res$DBI[i]        <- round(ei["DB"], 4)
          sweep_res$CH[i]         <- round(ei["CH"], 2)
          log_add(sprintf("  k=%d | Sil=%.4f | DB=%.4f | CH=%.2f", ki, s, ei["DB"], ei["CH"]))
        }
        rv$sweep_tab <- sweep_res
      }
    })

    rv$running <- FALSE
  })

  # ── Output: log ───────────────────────────────────────────────────────────
  output$log_out <- renderText({
    req(length(rv$log_msgs) > 0)
    paste(tail(rv$log_msgs, 20), collapse = "\n")
  })

  output$progress_ui <- renderUI({
    if (rv$running) {
      tagList(
        tags$div(class = "progress",
          tags$div(class = "progress-bar progress-bar-striped active",
                   style = paste0("width:", rv$prog_val, "%"),
                   paste0(rv$prog_val, "%"))),
        tags$p(rv$prog_msg, class = "progress-text")
      )
    }
  })

  # ── Info boxes ────────────────────────────────────────────────────────────
  output$info_boxes <- renderUI({
    req(rv$results)
    sil_vals <- sapply(rv$results, function(r) r$sil)
    best_algo <- names(which.max(sil_vals))
    best_sil  <- max(sil_vals)
    n_obs     <- nrow(rv$data_scaled)

    fluidRow(
      infoBox("Observations", n_obs, icon = icon("users"), color = "blue"),
      infoBox("Number of clusters", input$k_val, icon = icon("object-group"), color = "green"),
      infoBox("Best algorithm", best_algo, icon = icon("trophy"), color = "yellow"),
      infoBox("Best Silhouette", round(best_sil, 4), icon = icon("star"), color = "purple")
    )
  })

  # ── Plot silhouette ───────────────────────────────────────────────────────
  output$plot_sil <- renderPlot({
    req(rv$results)
    sil_vals <- sapply(rv$results, function(r) r$sil)
    best_algo <- names(which.max(sil_vals))
    best_res  <- rv$results[[best_algo]]
    lab <- best_res$labels
    sw  <- silhouette(lab, rv$DIST)

    sw_df <- data.frame(
      Observation = rownames(rv$data_scaled),
      Cluster   = factor(lab),
      SilWidth  = sw[, 3]
    )
    sw_df <- sw_df[order(sw_df$Cluster, sw_df$SilWidth), ]
    sw_df$Order <- seq_len(nrow(sw_df))

    ggplot(sw_df, aes(x = Order, y = SilWidth, fill = Cluster)) +
      geom_bar(stat = "identity", width = 0.85) +
      geom_hline(yintercept = mean(sw_df$SilWidth), linetype = "dashed",
                 color = "#c0392b", linewidth = 0.8) +
      scale_fill_brewer(palette = "Set2") +
      labs(title = paste0("Silhouette Plot – ", best_algo, " (k=", input$k_val, ")"),
           subtitle = paste0("Mean silhouette = ", round(mean(sw_df$SilWidth), 4)),
           x = "Observation", y = "Silhouette Width", fill = "Cluster") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
  })

  # ── Comparison plot ─────────────────────────────────────────────────────
  output$plot_compare <- renderPlot({
    req(rv$results)
    sil_vals <- sapply(rv$results, function(r) r$sil)
    df_comp  <- data.frame(
      Algorithm  = names(sil_vals),
      Silhouette = as.numeric(sil_vals)
    )
    df_comp$Algorithm <- factor(df_comp$Algorithm, levels = df_comp$Algorithm[order(df_comp$Silhouette)])

    ggplot(df_comp, aes(x = Algorithm, y = Silhouette, fill = Algorithm)) +
      geom_col(width = 0.6, show.legend = FALSE) +
      geom_text(aes(label = round(Silhouette, 4)), hjust = -0.1, size = 4) +
      scale_fill_brewer(palette = "Set1") +
      coord_flip(ylim = c(min(0, min(df_comp$Silhouette) - 0.05),
                           max(df_comp$Silhouette) + 0.1)) +
      labs(title = "Silhouette Comparison by Algorithm",
           x = NULL, y = "Silhouette Coefficient") +
      theme_minimal(base_size = 12)
  })

  # ── Membership table ─────────────────────────────────────────────────────
  output$tbl_member <- DT::renderDataTable({
    req(rv$results)
    sil_vals  <- sapply(rv$results, function(r) r$sil)
    best_algo <- names(which.max(sil_vals))
    best_lab  <- rv$results[[best_algo]]$labels
    best_med  <- rv$results[[best_algo]]$med
    med_names <- rownames(rv$data_scaled)[best_med]

    df_out <- data.frame(
      No       = seq_len(nrow(rv$data_scaled)),
      Observation = rownames(rv$data_scaled),
      Cluster  = best_lab,
      IsMedoid = rownames(rv$data_scaled) %in% med_names
    )

    DT::datatable(df_out, rownames = FALSE,
      options  = list(pageLength = 15, scrollX = TRUE)) |>
      DT::formatStyle("Cluster",
        backgroundColor = DT::styleEqual(
          sort(unique(best_lab)),
          RColorBrewer::brewer.pal(max(3, length(unique(best_lab))), "Set2")[seq_len(length(unique(best_lab)))]
        )) |>
      DT::formatStyle("IsMedoid",
        fontWeight = DT::styleEqual(TRUE, "bold"),
        color      = DT::styleEqual(TRUE, "#c0392b"))
  })

  # ── Cluster profile plot ────────────────────────────────────────────────────
  output$plot_profile <- renderPlot({
    req(rv$results, rv$data_num)
    sil_vals  <- sapply(rv$results, function(r) r$sil)
    best_algo <- names(which.max(sil_vals))
    best_lab  <- rv$results[[best_algo]]$labels

    profil <- aggregate(rv$data_num, by = list(Cluster = best_lab), FUN = mean)
    profil_long <- pivot_longer(profil, -Cluster, names_to = "Variable", values_to = "Value")
    profil_long$Cluster <- factor(profil_long$Cluster)

    ggplot(profil_long, aes(x = Variable, y = Value, fill = Cluster)) +
      geom_col(position = "dodge", width = 0.65) +
      scale_fill_brewer(palette = "Set2") +
      labs(title = paste0("Mean Cluster Profile – ", best_algo),
           x = "Variable", y = "Mean (original scale)", fill = "Cluster") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  })

  # ── Convergence curve ─────────────────────────────────────────────────────
  output$plot_conv <- renderPlot({
    req(rv$conv_list)
    convs <- rv$conv_list[sapply(rv$conv_list, function(x) !is.null(x) && length(x) > 0)]
    req(length(convs) > 0)

    conv_df <- bind_rows(lapply(names(convs), function(a) {
      v <- convs[[a]]
      data.frame(Algorithm = a, Iteration = seq_along(v), BestSilhouette = cummax(v))
    }))

    ggplot(conv_df, aes(Iteration, BestSilhouette, colour = Algorithm)) +
      geom_line(linewidth = 1) +
      scale_color_brewer(palette = "Set1") +
      labs(title = paste0("Convergence Curve (k = ", input$k_val, ")"),
           x = "Iteration / Generation", y = "Best Silhouette Coefficient",
           colour = "Algorithm") +
      theme_minimal(base_size = 13)
  })

  # ── Validation table ────────────────────────────────────────────────────────
  output$tbl_valid <- DT::renderDataTable({
    req(rv$valid_tab)
    DT::datatable(rv$valid_tab, rownames = FALSE,
      options = list(pageLength = 10)) |>
      DT::formatStyle("Silhouette",
        background  = DT::styleColorBar(range(rv$valid_tab$Silhouette, na.rm = TRUE), "#aed6f1"),
        backgroundSize = "100% 80%", backgroundRepeat = "no-repeat",
        backgroundPosition = "center")
  })

  # ── Plot sweep ────────────────────────────────────────────────────────────
  output$plot_sweep <- renderPlot({
    req(rv$sweep_tab)
    sw <- rv$sweep_tab |> pivot_longer(-k, names_to = "Index", values_to = "Value")

    ggplot(sw, aes(k, Value, colour = Index)) +
      geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
      scale_x_continuous(breaks = 2:10) +
      scale_color_brewer(palette = "Set1") +
      facet_wrap(~Index, scales = "free_y") +
      labs(title = "Validity Indices across k (GA)", x = "Number of clusters (k)", y = "Value") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none")
  })

  # ── Full-table tab ─────────────────────────────────────────────────────
  output$tbl_t1 <- DT::renderDataTable({ req(rv$mahal_tab); DT::datatable(rv$mahal_tab, rownames = FALSE) })
  output$tbl_t2 <- DT::renderDataTable({ req(rv$vif_tab);   DT::datatable(rv$vif_tab,   rownames = FALSE) })

  output$tbl_t3 <- DT::renderDataTable({
    req(rv$sweep_tab)
    DT::datatable(rv$sweep_tab, rownames = FALSE,
      options = list(pageLength = 10))
  })

  output$tbl_t4 <- DT::renderDataTable({
    req(rv$results)
    sil_vals  <- sapply(rv$results, function(r) r$sil)
    best_algo <- names(which.max(sil_vals))
    best_lab  <- rv$results[[best_algo]]$labels
    df_out <- data.frame(Observation = rownames(rv$data_scaled), Cluster = best_lab)
    DT::datatable(df_out, rownames = FALSE, options = list(pageLength = 15))
  })

  output$tbl_t5 <- DT::renderDataTable({
    req(rv$results, rv$data_num)
    sil_vals  <- sapply(rv$results, function(r) r$sil)
    best_algo <- names(which.max(sil_vals))
    best_lab  <- rv$results[[best_algo]]$labels
    profil <- aggregate(rv$data_num, by = list(Cluster = best_lab), FUN = mean)
    profil$Mean_all_variables <- rowMeans(profil[, -1])
    DT::datatable(round(profil, 3), rownames = FALSE)
  })

  output$tbl_t6 <- DT::renderDataTable({
    req(rv$valid_tab)
    DT::datatable(rv$valid_tab, rownames = FALSE)
  })

  # ── Download handlers ─────────────────────────────────────────────────────
  output$dl_t1 <- downloadHandler(filename = "mahalanobis.csv",
    content = function(f) write.csv(rv$mahal_tab, f, row.names = FALSE))
  output$dl_t2 <- downloadHandler(filename = "vif.csv",
    content = function(f) write.csv(rv$vif_tab, f, row.names = FALSE))
  output$dl_t3 <- downloadHandler(filename = "silhouette_sweep.csv",
    content = function(f) write.csv(rv$sweep_tab, f, row.names = FALSE))
  output$dl_t4 <- downloadHandler(filename = "cluster_membership.csv",
    content = function(f) {
      req(rv$results)
      sil_vals  <- sapply(rv$results, function(r) r$sil)
      best_algo <- names(which.max(sil_vals))
      best_lab  <- rv$results[[best_algo]]$labels
      df_out <- data.frame(Observation = rownames(rv$data_scaled), Cluster = best_lab)
      write.csv(df_out, f, row.names = FALSE)
    })
  output$dl_t5 <- downloadHandler(filename = "cluster_profile.csv",
    content = function(f) {
      req(rv$results, rv$data_num)
      sil_vals  <- sapply(rv$results, function(r) r$sil)
      best_algo <- names(which.max(sil_vals))
      profil <- aggregate(rv$data_num, by = list(Cluster = rv$results[[best_algo]]$labels), FUN = mean)
      write.csv(round(profil, 4), f, row.names = FALSE)
    })
  output$dl_t6 <- downloadHandler(filename = "validity_indices.csv",
    content = function(f) write.csv(rv$valid_tab, f, row.names = FALSE))
}

# ─────────────────────────────────────────────────────────────────────────────
shinyApp(ui, server)
