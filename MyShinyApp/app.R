# =============================================================================
# SHINY APP: K-MEDOIDS CLUSTERING DENGAN OPTIMASI METAHEURISTIK
# Metode: GA, PSO, CSA (Crow Search), GWO (Grey Wolf Optimizer)
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
    tambahan <- setdiff(seq_len(n), idx)
    idx <- c(idx, sample(tambahan, k - length(idx)))
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

fitness_sil <- function(sol, k, n, DMAT, DIST, maximize = TRUE) {
  s <- sil_of_medoids(DMAT, DIST, decode_medoid(sol, k, n))$sil
  if (maximize) s else -s
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
  g <- ga(
    type    = "real-valued",
    fitness = function(sol) fitness_sil(sol, k, n, DMAT, DIST, maximize = TRUE),
    lower   = rep(1, k), upper = rep(n, k),
    popSize = popSize, maxiter = maxiter,
    pcrossover = 0.8, pmutation = 0.1, elitism = 2,
    monitor = if (!is.null(progress_fn)) {
      function(obj) {
        if (!is.null(progress_fn)) {
          iter <- obj@iter
          progress_fn(iter, maxiter, paste0("GA Generasi ", iter, "/", maxiter))
        }
      }
    } else FALSE,
    seed = seed
  )
  list(
    med  = decode_medoid(as.numeric(g@solution[1, ]), k, n),
    conv = if (trace) as.numeric(g@summary[, "max"]) else NULL,
    sil  = max(g@summary[, "max"])
  )
}

run_PSO <- function(k, n, DMAT, DIST, seed = 42, s = 30, maxit = 100,
                    trace = FALSE, progress_fn = NULL) {
  set.seed(seed)
  iter_count <- 0
  ctrl <- list(
    maxit = maxit, s = s, w = 0.729, c.p = 1.494, c.g = 1.494,
    trace = if (trace) 1 else 0, trace.stats = trace, REPORT = 1
  )
  p <- psoptim(
    par    = sample(seq_len(n), k),
    fn     = function(sol) {
      iter_count <<- iter_count + 1
      if (!is.null(progress_fn) && iter_count %% k == 0) {
        prog_iter <- ceiling(iter_count / (s * k))
        progress_fn(min(prog_iter, maxit), maxit,
                    paste0("PSO Iterasi ~", min(prog_iter, maxit), "/", maxit))
      }
      fitness_sil(sol, k, n, DMAT, DIST, maximize = FALSE)
    },
    lower  = rep(1, k), upper = rep(n, k),
    control = ctrl
  )
  best_sil <- -p$value
  list(
    med  = decode_medoid(p$par, k, n),
    conv = if (trace && !is.null(p$stats$error)) -as.numeric(p$stats$error) else NULL,
    sil  = best_sil
  )
}

run_CSA <- function(k, n, DMAT, DIST, seed = 42, flock = 30, num_iter = 100,
                    fl = 2.0, AP = 0.1, trace = FALSE, progress_fn = NULL) {
  set.seed(seed)
  pos     <- lapply(seq_len(flock), function(i) decode_medoid(sample(n, k), k, n))
  memory  <- pos
  sil_mem <- sapply(memory, function(m) sil_of_medoids(DMAT, DIST, m)$sil)
  best    <- which.max(sil_mem)
  best_med <- memory[[best]]; best_sil <- sil_mem[best]
  conv <- numeric(0)
  for (it in seq_len(num_iter)) {
    if (!is.null(progress_fn)) progress_fn(it, num_iter, paste0("CSA Iterasi ", it, "/", num_iter))
    for (i in seq_len(flock)) {
      j <- sample(setdiff(seq_len(flock), i), 1)
      if (runif(1) >= AP) {
        new_med <- integer(k)
        for (d in seq_len(k))
          new_med[d] <- if (runif(1) < fl / (fl + 1)) memory[[j]][sample(k, 1)] else pos[[i]][d]
      } else new_med <- sample(n, k)
      new_med <- decode_medoid(new_med, k, n)
      s <- sil_of_medoids(DMAT, DIST, new_med)$sil
      pos[[i]] <- new_med
      if (s > sil_mem[i]) {
        memory[[i]] <- new_med; sil_mem[i] <- s
        if (s > best_sil) { best_sil <- s; best_med <- new_med }
      }
    }
    conv <- c(conv, best_sil)
  }
  list(med = best_med, conv = if (trace) conv else NULL, sil = best_sil)
}

run_GWO <- function(k, n, DMAT, DIST, seed = 42, num_wolves = 30, num_iter = 100,
                    trace = FALSE, progress_fn = NULL) {
  set.seed(seed)
  pos <- matrix(runif(num_wolves * k, 1, n), nrow = num_wolves)
  fit <- apply(pos, 1, function(p) sil_of_medoids(DMAT, DIST, decode_medoid(p, k, n))$sil)
  ord <- order(fit, decreasing = TRUE)
  alpha <- pos[ord[1], ]; beta <- pos[ord[2], ]; delta <- pos[ord[3], ]
  best_med <- decode_medoid(alpha, k, n); best_sil <- fit[ord[1]]
  conv <- numeric(0)
  for (it in seq_len(num_iter)) {
    if (!is.null(progress_fn)) progress_fn(it, num_iter, paste0("GWO Iterasi ", it, "/", num_iter))
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
    fit <- apply(pos, 1, function(p) sil_of_medoids(DMAT, DIST, decode_medoid(p, k, n))$sil)
    ord <- order(fit, decreasing = TRUE)
    alpha <- pos[ord[1], ]; beta <- pos[ord[2], ]; delta <- pos[ord[3], ]
    if (fit[ord[1]] > best_sil) {
      best_sil <- fit[ord[1]]
      best_med <- decode_medoid(alpha, k, n)
    }
    conv <- c(conv, best_sil)
  }
  list(med = best_med, conv = if (trace) conv else NULL, sil = best_sil)
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
  skin = "blue",
  dashboardHeader(title = span(icon("leaf"), " K-Medoids Metaheuristik")),

  dashboardSidebar(
    sidebarMenu(
      menuItem("📂 Upload & Preprocessing", tabName = "preproc", icon = icon("database")),
      menuItem("⚙️ Konfigurasi & Jalankan",  tabName = "run",     icon = icon("cogs")),
      menuItem("📊 Hasil Clustering",         tabName = "results", icon = icon("chart-bar")),
      menuItem("📈 Kurva Konvergensi",         tabName = "conv",    icon = icon("line-chart")),
      menuItem("🔍 Analisis Validasi",          tabName = "valid",   icon = icon("check-circle")),
      menuItem("📋 Tabel Lengkap",              tabName = "tables",  icon = icon("table"))
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
          box(title = "Preview Data Asli", width = 12, status = "primary", solidHeader = TRUE,
            fileInput("file_csv", "Pilih file CSV:", accept = ".csv"),
            checkboxInput("header_cb", "Baris pertama sebagai header", value = TRUE),
            selectInput("sep_cb", "Separator:", choices = c("Koma ,"=",","Titik-koma ;"=";","Tab"="\t"), selected = ","),
            uiOutput("col_name_ui"),
            uiOutput("col_drop_ui"),
            actionButton("btn_preprocess", "Proses Data", icon = icon("play"), class = "btn-primary btn-lg")
          ),
          box(title = "Panduan Format Data", width = 6, status = "info", solidHeader = TRUE,
            tags$ul(
              tags$li("Kolom pertama sebaiknya adalah ", tags$b("nama observasi"), " (provinsi, kota, dll.)"),
              tags$li("Kolom numerik akan digunakan sebagai variabel clustering"),
              tags$li("Kolom non-numerik akan diabaikan otomatis"),
              tags$li("Nilai kosong (NA) akan dideteksi dan dilaporkan"),
              tags$li("Data akan di-standarisasi ", tags$b("z-score"), " sebelum clustering")
            ),
            tags$hr(),
            tags$p("Format CSV contoh (IKLH Provinsi Indonesia 2023):"),
            tags$code("Provinsi,IKU,IKA,IKTL,IKL"),
            tags$br(),
            tags$code("Aceh,72.3,68.1,80.2,74.5"),
            tags$br(), tags$code("...dst")
          )
        ),
        fluidRow(
          box(title = "Preview Data Asli", width = 12, status = "primary",
            DT::dataTableOutput("tbl_raw"))
        ),
        fluidRow(
          box(title = "Mahalanobis Distance (Deteksi Outlier)", width = 6, status = "warning", solidHeader = TRUE,
            DT::dataTableOutput("tbl_mahal")),
          box(title = "VIF – Multikolinearitas", width = 6, status = "warning", solidHeader = TRUE,
            DT::dataTableOutput("tbl_vif"),
            tags$p(tags$b("Interpretasi:"), " VIF < 5 = tidak ada multikolinearitas serius", class = "text-muted")
          )
        )
      ),

      # ── TAB 2: Konfigurasi & Jalankan ─────────────────────────────────────
      tabItem(tabName = "run",
        fluidRow(
          box(title = "Konfigurasi Algoritma", width = 5, status = "primary", solidHeader = TRUE,
            sliderInput("k_val", "Jumlah Klaster (k):", min = 2, max = 10, value = 3, step = 1),
            checkboxGroupInput("algo_sel", "Pilih Algoritma:",
              choices  = c("PAM (baseline)"="PAM", "GA"="GA", "PSO"="PSO", "CSA"="CSA", "GWO"="GWO"),
              selected = c("PAM","GA","PSO","CSA","GWO")),
            numericInput("seed_val", "Random Seed:", value = 42, min = 1),
            tags$hr(),
            tags$h5("Parameter GA"),
            sliderInput("ga_pop",  "Ukuran Populasi:", min = 10, max = 200, value = 50, step = 10),
            sliderInput("ga_iter", "Maks. Generasi:",  min = 50, max = 500, value = 100, step = 50),
            tags$hr(),
            tags$h5("Parameter PSO/CSA/GWO"),
            sliderInput("swarm", "Ukuran Swarm/Flock:", min = 10, max = 100, value = 30, step = 10),
            sliderInput("gen_iter", "Maks. Iterasi:",   min = 50, max = 500, value = 100, step = 50)
          ),
          box(title = "Mode Analisis", width = 4, status = "info", solidHeader = TRUE,
            radioButtons("mode_sel", "Pilih mode:",
              choices = c("Satu nilai k (cepat)" = "single",
                          "Sweep k = 2..10 (komprehensif)" = "sweep"),
              selected = "single"),
            tags$p("Mode sweep akan menjalankan semua algoritma untuk setiap k dari 2 hingga 10.",
                   class = "text-muted"),
            tags$hr(),
            tags$h5("Kolom label medoid"),
            tags$p("Kolom nama (kolom pertama data) akan digunakan sebagai label medoid.", class = "text-muted")
          ),
          box(title = "Jalankan", width = 3, status = "success", solidHeader = TRUE,
            actionButton("btn_run", "▶ Mulai Clustering", icon = icon("rocket"),
                         class = "btn-success btn-lg", style = "width:100%; margin-bottom:10px"),
            tags$hr(),
            uiOutput("progress_ui"),
            verbatimTextOutput("log_out", placeholder = TRUE)
          )
        )
      ),

      # ── TAB 3: Hasil Clustering ────────────────────────────────────────────
      tabItem(tabName = "results",
        fluidRow(
          uiOutput("info_boxes")
        ),
        fluidRow(
          box(title = "Plot Silhouette", width = 6, status = "primary", solidHeader = TRUE,
            plotOutput("plot_sil", height = "350px")),
          box(title = "Perbandingan Skor Silhouette", width = 6, status = "primary", solidHeader = TRUE,
            plotOutput("plot_compare", height = "350px"))
        ),
        fluidRow(
          box(title = "Keanggotaan Klaster (Algoritma Terbaik)", width = 12, status = "success", solidHeader = TRUE,
            DT::dataTableOutput("tbl_member"))
        ),
        fluidRow(
          box(title = "Profil Rata-rata Per Klaster", width = 12, status = "warning", solidHeader = TRUE,
            plotOutput("plot_profile", height = "350px"))
        )
      ),

      # ── TAB 4: Kurva Konvergensi ───────────────────────────────────────────
      tabItem(tabName = "conv",
        fluidRow(
          box(title = "Kurva Konvergensi – Best Silhouette per Iterasi", width = 12,
              status = "primary", solidHeader = TRUE,
            tags$p("Hanya tersedia untuk algoritma iteratif (GA, PSO, CSA, GWO). PAM tidak memiliki kurva konvergensi.", class = "text-muted"),
            plotOutput("plot_conv", height = "420px"))
        )
      ),

      # ── TAB 5: Analisis Validasi ───────────────────────────────────────────
      tabItem(tabName = "valid",
        fluidRow(
          box(title = "Indeks Validasi Klaster (k final)", width = 12, status = "primary", solidHeader = TRUE,
            DT::dataTableOutput("tbl_valid"),
            tags$p("Silhouette: tinggi = baik | Davies-Bouldin: rendah = baik | Calinski-Harabasz: tinggi = baik",
                   class = "text-muted")
          )
        ),
        fluidRow(
          box(title = "Sweep Indeks Lintas k (GA)", width = 12, status = "info", solidHeader = TRUE,
            plotOutput("plot_sweep", height = "380px"))
        )
      ),

      # ── TAB 6: Tabel Lengkap ──────────────────────────────────────────────
      tabItem(tabName = "tables",
        tabBox(width = 12,
          tabPanel("Mahalanobis", DT::dataTableOutput("tbl_t1"), downloadButton("dl_t1", "Unduh CSV")),
          tabPanel("VIF", DT::dataTableOutput("tbl_t2"), downloadButton("dl_t2", "Unduh CSV")),
          tabPanel("Silhouette k=2..10", DT::dataTableOutput("tbl_t3"), downloadButton("dl_t3", "Unduh CSV")),
          tabPanel("Keanggotaan", DT::dataTableOutput("tbl_t4"), downloadButton("dl_t4", "Unduh CSV")),
          tabPanel("Profil Klaster", DT::dataTableOutput("tbl_t5"), downloadButton("dl_t5", "Unduh CSV")),
          tabPanel("Validasi", DT::dataTableOutput("tbl_t6"), downloadButton("dl_t6", "Unduh CSV"))
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
    log_add(paste0("Data dimuat: ", nrow(rv$df), " baris, ", ncol(rv$df), " kolom"))
  })

  output$tbl_raw <- DT::renderDataTable({
    req(rv$df)
    DT::datatable(rv$df, options = list(scrollX = TRUE, pageLength = 8), rownames = FALSE)
  })

  output$col_name_ui <- renderUI({
    req(rv$df)
    selectInput("col_name", "Kolom nama/label observasi:",
                choices  = names(rv$df),
                selected = names(rv$df)[1])
  })

  output$col_drop_ui <- renderUI({
    req(rv$df)
    checkboxGroupInput("cols_drop", "Kolom yang dikecualikan (opsional):",
                       choices = names(rv$df)[-1])
  })

  # ── Preprocessing ─────────────────────────────────────────────────────────
  observeEvent(input$btn_preprocess, {
    req(rv$df, input$col_name)
    df <- rv$df

    # Label kolom
    prov_col <- input$col_name
    rv$provinsi <- df[[prov_col]]

    # Drop kolom
    drop_cols <- c(prov_col, input$cols_drop)
    data_num  <- df[, !names(df) %in% drop_cols]
    data_num  <- data_num[, sapply(data_num, is.numeric), drop = FALSE]
    rv$data_num <- data_num
    rv$col_names <- colnames(data_num)
    log_add(paste0("Kolom numerik: ", paste(rv$col_names, collapse = ", ")))

    # NA check
    na_count <- sum(is.na(data_num))
    if (na_count > 0) {
      log_add(paste0("Peringatan: ", na_count, " nilai NA ditemukan – baris bermasalah dihapus"))
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
      Observasi      = rv$provinsi,
      Mahalanobis_D2 = round(mahal, 4),
      Cutoff         = round(cutoff, 4),
      Outlier        = mahal > cutoff
    )
    log_add(paste0("Deteksi outlier selesai. Outlier: ", sum(rv$mahal_tab$Outlier)))

    # VIF
    vif_vals <- vif_calc(data_num)
    rv$vif_tab <- data.frame(Variabel = rv$col_names, VIF = round(vif_vals, 4))

    # Standarisasi
    data_scaled     <- scale(data_num)
    rownames(data_scaled) <- rv$provinsi
    rv$data_scaled  <- data_scaled
    rv$DIST         <- dist(data_scaled)
    rv$DMAT         <- as.matrix(rv$DIST)
    rv$XMAT         <- as.matrix(data_scaled)
    log_add("Standarisasi z-score selesai. Siap untuk clustering.")
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

    withProgress(message = "Menjalankan clustering...", value = 0, {

      make_progress <- function(label) {
        function(cur, total, msg = "") {
          setProgress(value = cur / total, message = label, detail = msg)
        }
      }

      results   <- list()
      conv_list <- list()

      if ("PAM" %in% algos) {
        log_add("Menjalankan PAM (baseline)...")
        incProgress(0, message = "PAM baseline...")
        t0  <- proc.time()["elapsed"]
        res <- run_PAM(k, rv$data_scaled, rv$DMAT, rv$DIST)
        runtime <- as.numeric(proc.time()["elapsed"] - t0)
        lab <- labels_of_medoids(rv$DMAT, res$med)
        ei  <- extra_indices(rv$XMAT, lab)
        results[["PAM"]] <- list(med = res$med, sil = res$sil, labels = lab,
                                  DB = ei["DB"], CH = ei["CH"], runtime = runtime)
        log_add(sprintf("PAM selesai | Silhouette = %.4f | %.2fs", res$sil, runtime))
      }

      if ("GA" %in% algos) {
        log_add("Menjalankan GA...")
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
        log_add(sprintf("GA selesai | Silhouette = %.4f | %.2fs", res$sil, runtime))
      }

      if ("PSO" %in% algos) {
        log_add("Menjalankan PSO...")
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
        log_add(sprintf("PSO selesai | Silhouette = %.4f | %.2fs", res$sil, runtime))
      }

      if ("CSA" %in% algos) {
        log_add("Menjalankan CSA (Crow Search)...")
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
        log_add(sprintf("CSA selesai | Silhouette = %.4f | %.2fs", res$sil, runtime))
      }

      if ("GWO" %in% algos) {
        log_add("Menjalankan GWO (Grey Wolf Optimizer)...")
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
        log_add(sprintf("GWO selesai | Silhouette = %.4f | %.2fs", res$sil, runtime))
      }

      rv$results   <- results
      rv$conv_list <- conv_list

      # Tabel validasi
      val_df <- data.frame(
        Algoritma  = names(results),
        Silhouette = sapply(results, function(r) round(r$sil, 4)),
        DBI        = sapply(results, function(r) round(r$DB, 4)),
        CH         = sapply(results, function(r) round(r$CH, 2)),
        Runtime_s  = sapply(results, function(r) round(r$runtime, 3))
      )
      rv$valid_tab <- val_df
      log_add("Semua algoritma selesai!")

      # Sweep k=2..10 dengan GA jika mode sweep
      if (input$mode_sel == "sweep" && "GA" %in% algos) {
        log_add("Sweep k=2..10 dengan GA...")
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
      infoBox("Observasi", n_obs, icon = icon("users"), color = "blue"),
      infoBox("Jumlah Klaster", input$k_val, icon = icon("object-group"), color = "green"),
      infoBox("Algoritma Terbaik", best_algo, icon = icon("trophy"), color = "yellow"),
      infoBox("Silhouette Terbaik", round(best_sil, 4), icon = icon("star"), color = "purple")
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
      Observasi = rownames(rv$data_scaled),
      Cluster   = factor(lab),
      SilWidth  = sw[, 3]
    ) |> arrange(Cluster, SilWidth) |>
      mutate(Order = row_number())

    ggplot(sw_df, aes(x = Order, y = SilWidth, fill = Cluster)) +
      geom_bar(stat = "identity", width = 0.85) +
      geom_hline(yintercept = mean(sw_df$SilWidth), linetype = "dashed",
                 color = "#c0392b", linewidth = 0.8) +
      scale_fill_brewer(palette = "Set2") +
      labs(title = paste0("Silhouette Plot – ", best_algo, " (k=", input$k_val, ")"),
           subtitle = paste0("Rata-rata silhouette = ", round(mean(sw_df$SilWidth), 4)),
           x = "Observasi", y = "Silhouette Width", fill = "Klaster") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
  })

  # ── Plot perbandingan ─────────────────────────────────────────────────────
  output$plot_compare <- renderPlot({
    req(rv$results)
    sil_vals <- sapply(rv$results, function(r) r$sil)
    df_comp  <- data.frame(
      Algoritma  = names(sil_vals),
      Silhouette = as.numeric(sil_vals)
    )
    df_comp$Algoritma <- factor(df_comp$Algoritma, levels = df_comp$Algoritma[order(df_comp$Silhouette)])

    ggplot(df_comp, aes(x = Algoritma, y = Silhouette, fill = Algoritma)) +
      geom_col(width = 0.6, show.legend = FALSE) +
      geom_text(aes(label = round(Silhouette, 4)), hjust = -0.1, size = 4) +
      scale_fill_brewer(palette = "Set1") +
      coord_flip(ylim = c(min(0, min(df_comp$Silhouette) - 0.05),
                           max(df_comp$Silhouette) + 0.1)) +
      labs(title = "Perbandingan Silhouette per Algoritma",
           x = NULL, y = "Silhouette Coefficient") +
      theme_minimal(base_size = 12)
  })

  # ── Tabel keanggotaan ─────────────────────────────────────────────────────
  output$tbl_member <- DT::renderDataTable({
    req(rv$results)
    sil_vals  <- sapply(rv$results, function(r) r$sil)
    best_algo <- names(which.max(sil_vals))
    best_lab  <- rv$results[[best_algo]]$labels
    best_med  <- rv$results[[best_algo]]$med
    med_names <- rownames(rv$data_scaled)[best_med]

    df_out <- data.frame(
      No       = seq_len(nrow(rv$data_scaled)),
      Observasi = rownames(rv$data_scaled),
      Klaster  = best_lab,
      IsMedoid = rownames(rv$data_scaled) %in% med_names
    )

    DT::datatable(df_out, rownames = FALSE,
      options  = list(pageLength = 15, scrollX = TRUE)) |>
      DT::formatStyle("Klaster",
        backgroundColor = DT::styleEqual(
          sort(unique(best_lab)),
          RColorBrewer::brewer.pal(max(3, length(unique(best_lab))), "Set2")[seq_len(length(unique(best_lab)))]
        )) |>
      DT::formatStyle("IsMedoid",
        fontWeight = DT::styleEqual(TRUE, "bold"),
        color      = DT::styleEqual(TRUE, "#c0392b"))
  })

  # ── Plot profil klaster ────────────────────────────────────────────────────
  output$plot_profile <- renderPlot({
    req(rv$results, rv$data_num)
    sil_vals  <- sapply(rv$results, function(r) r$sil)
    best_algo <- names(which.max(sil_vals))
    best_lab  <- rv$results[[best_algo]]$labels

    profil <- aggregate(rv$data_num, by = list(Klaster = best_lab), FUN = mean)
    profil_long <- pivot_longer(profil, -Klaster, names_to = "Variabel", values_to = "Nilai")
    profil_long$Klaster <- factor(profil_long$Klaster)

    ggplot(profil_long, aes(x = Variabel, y = Nilai, fill = Klaster)) +
      geom_col(position = "dodge", width = 0.65) +
      scale_fill_brewer(palette = "Set2") +
      labs(title = paste0("Profil Rata-rata Klaster – ", best_algo),
           x = "Variabel", y = "Rata-rata (skala asli)", fill = "Klaster") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  })

  # ── Kurva konvergensi ─────────────────────────────────────────────────────
  output$plot_conv <- renderPlot({
    req(rv$conv_list)
    convs <- rv$conv_list[sapply(rv$conv_list, function(x) !is.null(x) && length(x) > 0)]
    req(length(convs) > 0)

    conv_df <- bind_rows(lapply(names(convs), function(a) {
      v <- convs[[a]]
      data.frame(Algoritma = a, Iterasi = seq_along(v), BestSilhouette = cummax(v))
    }))

    ggplot(conv_df, aes(Iterasi, BestSilhouette, colour = Algoritma)) +
      geom_line(linewidth = 1) +
      scale_color_brewer(palette = "Set1") +
      labs(title = paste0("Kurva Konvergensi (k = ", input$k_val, ")"),
           x = "Iterasi / Generasi", y = "Best Silhouette Coefficient",
           colour = "Algoritma") +
      theme_minimal(base_size = 13)
  })

  # ── Tabel validasi ────────────────────────────────────────────────────────
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
    sw <- rv$sweep_tab |> pivot_longer(-k, names_to = "Indeks", values_to = "Nilai")

    ggplot(sw, aes(k, Nilai, colour = Indeks)) +
      geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
      scale_x_continuous(breaks = 2:10) +
      scale_color_brewer(palette = "Set1") +
      facet_wrap(~Indeks, scales = "free_y") +
      labs(title = "Indeks Validasi Lintas k (GA)", x = "Jumlah Klaster (k)", y = "Nilai") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none")
  })

  # ── Tabel tab lengkap ─────────────────────────────────────────────────────
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
    df_out <- data.frame(Observasi = rownames(rv$data_scaled), Klaster = best_lab)
    DT::datatable(df_out, rownames = FALSE, options = list(pageLength = 15))
  })

  output$tbl_t5 <- DT::renderDataTable({
    req(rv$results, rv$data_num)
    sil_vals  <- sapply(rv$results, function(r) r$sil)
    best_algo <- names(which.max(sil_vals))
    best_lab  <- rv$results[[best_algo]]$labels
    profil <- aggregate(rv$data_num, by = list(Klaster = best_lab), FUN = mean)
    profil$Rata_rata_total <- rowMeans(profil[, -1])
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
  output$dl_t4 <- downloadHandler(filename = "keanggotaan_klaster.csv",
    content = function(f) {
      req(rv$results)
      sil_vals  <- sapply(rv$results, function(r) r$sil)
      best_algo <- names(which.max(sil_vals))
      best_lab  <- rv$results[[best_algo]]$labels
      df_out <- data.frame(Observasi = rownames(rv$data_scaled), Klaster = best_lab)
      write.csv(df_out, f, row.names = FALSE)
    })
  output$dl_t5 <- downloadHandler(filename = "profil_klaster.csv",
    content = function(f) {
      req(rv$results, rv$data_num)
      sil_vals  <- sapply(rv$results, function(r) r$sil)
      best_algo <- names(which.max(sil_vals))
      profil <- aggregate(rv$data_num, by = list(Klaster = rv$results[[best_algo]]$labels), FUN = mean)
      write.csv(round(profil, 4), f, row.names = FALSE)
    })
  output$dl_t6 <- downloadHandler(filename = "validasi_indeks.csv",
    content = function(f) write.csv(rv$valid_tab, f, row.names = FALSE))
}

# ─────────────────────────────────────────────────────────────────────────────
shinyApp(ui, server)
