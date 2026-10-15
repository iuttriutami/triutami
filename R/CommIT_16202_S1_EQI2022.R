# =============================================================================
# CommIT 16202 - SUPPLEMENTARY S1: replication on EQI 2022 (Reviewer 2: generalization)
# Same encoding, repair, fitness, optimizers and budget as the 2023 analysis.
#   - Mahalanobis screening, PAM baseline, exhaustive global optimum (k = 2..5)
#   - equal-budget GA/PSO/CSA/GWO/RS: 3,000 evaluations, 30 runs, k = 3..5
#   - agreement of the 2022 and 2023 partitions (ARI on provinces present in both years)
# Output: /kaggle/working/commIT_supplementary/supp_S1_*.csv + supp_placeholders_S1.csv
# =============================================================================
PATH_CSV_2023 <- "/kaggle/input/datasets/triutamiiut/2023-eqi/Komponen Penyusun Indeks Kualitas Lingkungan Hidup Menurut Provinsi 2023.csv"
PATH_CSV_2022 <- ""   # leave "" to search /kaggle/input for a CSV with "2022" in its name
SEED <- 42L; N_RUN <- 30L; BUDGET <- 3000L
EB_K <- 3:5; EXH_K <- 2:5
DIR_OUT <- if (dir.exists("/kaggle/working")) "/kaggle/working/commIT_supplementary" else "commIT_supplementary"
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

for (p in c("cluster", "GA", "pso")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(library(cluster))
T0 <- proc.time()[["elapsed"]]
TOK <- new.env(); tok <- function(k, v) assign(k, as.character(v), envir = TOK)
f4 <- function(x) sprintf("%.4f", x); f3 <- function(x) sprintf("%.3f", x)
fpct <- function(x) sprintf("%.0f%%", 100 * x)
say <- function(...) cat(sprintf(...), "\n")
wcsv <- function(x, f) write.csv(x, file.path(DIR_OUT, f), row.names = FALSE)
find_csv <- function(path, pattern) {
  if (nzchar(path) && file.exists(path)) return(path)
  cand <- list.files("/kaggle/input", pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  cand <- cand[grepl(pattern, basename(cand)) & grepl("Kualitas|EQI|IKLH", basename(cand), ignore.case = TRUE)]
  if (!length(cand)) stop("CSV not found for pattern '", pattern, "': add the dataset via 'Add Input'.")
  message("Using: ", cand[1]); cand[1]
}
PATH_CSV_2023 <- find_csv(PATH_CSV_2023, "2023")
PATH_CSV_2022 <- find_csv(PATH_CSV_2022, "2022")

# ---- DATA ---------------------------------------------------------------------------
read_eqi <- function(path) {
  hdr <- grep("Indeks Kualitas Udara", readLines(path, warn = FALSE))[1]   # handles raw BPS downloads with title rows
  if (is.na(hdr)) stop("No 'Indeks Kualitas Udara' header found in ", basename(path))
  df <- read.csv(path, header = TRUE, check.names = FALSE, skip = hdr - 1, na.strings = c("", "NA", "-"))
  names(df)[1] <- "Provinsi"
  df <- df[!toupper(trimws(df$Provinsi)) %in% c("INDONESIA", ""), , drop = FALSE]
  nm <- names(df)
  pick <- function(ok) { j <- which(ok); if (length(j) != 1) stop("Column not identified in ", basename(path), ": ", paste(nm[j], collapse = "; ")); j }
  j <- c(AQI  = pick(grepl("Kualitas Udara", nm)),
         WQI  = pick(grepl("Kualitas Air", nm) & !grepl("Laut", nm)),
         LCQI = pick(grepl("Tutupan Lahan", nm)),
         LQI  = pick(grepl("Kualitas Lahan", nm)))
  X <- as.data.frame(lapply(df[, j], function(v) as.numeric(gsub(",", ".", as.character(v)))))
  names(X) <- names(j)
  prov <- toupper(trimws(gsub("\\s+", " ", as.character(df$Provinsi))))
  keep <- stats::complete.cases(X)
  if (any(!keep)) message("Dropped rows with missing components: ", paste(prov[!keep], collapse = ", "))
  list(prov = prov[keep], X = X[keep, ])
}
set_data <- function(d) {        # sets the globals used by the shared helper functions
  assign("provinsi", d$prov, envir = globalenv()); assign("n", nrow(d$X), envir = globalenv())
  Z <- scale(d$X); assign("DIST", dist(Z), envir = globalenv())
  assign("DMAT", as.matrix(dist(Z)), envir = globalenv())
}

# ---- SHARED HELPERS (verbatim from the main script) -----------------------------------
labels_of <- function(med, D = DMAT) max.col(-D[, med, drop = FALSE], ties.method = "first")
fast_sil_lab <- function(lab, D = DMAT) {
  k <- max(lab); nn <- length(lab)
  if (length(unique(lab)) < 2) return(-1)
  M <- matrix(0, nn, k); M[cbind(seq_len(nn), lab)] <- 1
  sz <- colSums(M); S <- D %*% M
  idx <- cbind(seq_len(nn), lab)
  a <- S[idx] / pmax(sz[lab] - 1, 1)
  avg <- sweep(S, 2, pmax(sz, 1), "/"); avg[, sz == 0] <- Inf; avg[idx] <- Inf
  b <- avg[cbind(seq_len(nn), max.col(-avg, ties.method = "first"))]
  s <- (b - a) / pmax(a, b); s[sz[lab] == 1] <- 0
  mean(s)
}
sil_med <- function(med, D = DMAT) fast_sil_lab(labels_of(med, D), D)
decode_medoid <- function(sol, k) {                 # identical to v9 encoding/repair
  idx <- unique(pmax(1L, pmin(n, as.integer(round(sol)))))
  if (length(idx) < k) idx <- c(idx, sample(setdiff(seq_len(n), idx), k - length(idx)))
  idx[seq_len(k)]
}
ari <- function(a, b) {                              # Hubert-Arabie adjusted Rand index
  tab <- table(a, b); c2 <- function(x) sum(x * (x - 1) / 2)
  s <- c2(tab); sa <- c2(rowSums(tab)); sb <- c2(colSums(tab)); tot <- c2(length(a))
  e <- sa * sb / tot; mx <- (sa + sb) / 2
  if (mx == e) return(1); (s - e) / (mx - e)
}


# ---- COUNTED FITNESS ---------------------------------------------------------------
new_counter <- function(k, budget, target = NA_real_) {
  st <- new.env(parent = emptyenv())
  st$calls <- 0L; st$best <- -Inf; st$best_med <- NULL; st$hit <- NA_integer_
  st$trace <- numeric(min(budget, 20000L)); st$t_budget <- NA_real_
  st$t0 <- proc.time()[["elapsed"]]
  st$fitness <- function(sol) {
    if (st$calls >= budget) stop(structure(class = c("FitnessBudgetReached", "error", "condition"),
                                           list(message = "budget reached", call = NULL)))
    med <- decode_medoid(sol, k)
    s <- sil_med(med)
    st$calls <- st$calls + 1L
    if (s > st$best + 1e-12) { st$best <- s; st$best_med <- med }
    if (is.na(st$hit) && !is.na(target) && s >= target - 1e-9) st$hit <- st$calls
    if (st$calls > length(st$trace)) st$trace <- c(st$trace, numeric(length(st$trace)))
    st$trace[st$calls] <- st$best
    if (st$calls == budget) st$t_budget <- proc.time()[["elapsed"]] - st$t0
    s
  }
  st
}

# ---- OPTIMIZERS (same encoding, repair and fitness for all) --------------------------
# GA / PSO: the CRAN packages GA::ga and pso::psoptim, exactly as in v9.
# CSA: discrete variant of Askarzadeh's rule (see manuscript Sec. II-D).
# GWO: continuous positions decoded by the shared rounding/repair decoder.
run_GA <- function(f, k, pop = 50, maxiter = 100000L)
  GA::ga(type = "real-valued", fitness = function(s) f$fitness(s), lower = rep(1, k),
         upper = rep(n, k), popSize = pop, maxiter = maxiter, run = maxiter,
         pcrossover = 0.8, pmutation = 0.1, elitism = 2, monitor = FALSE)
run_PSO <- function(f, k, pop = 30, maxiter = 100000L)
  pso::psoptim(par = sample(seq_len(n), k), fn = function(s) -f$fitness(s),
               lower = rep(1, k), upper = rep(n, k),
               control = list(maxit = maxiter, s = pop, w = 0.729, c.p = 1.494, c.g = 1.494, trace = 0))
run_CSA <- function(f, k, pop = 30, maxiter = 100000L, fl = 2, AP = 0.1) {
  pos <- lapply(seq_len(pop), function(i) decode_medoid(sample(n, k), k))
  mem <- pos; fmem <- vapply(mem, f$fitness, numeric(1))
  for (it in seq_len(maxiter)) for (i in seq_len(pop)) {
    j <- sample(setdiff(seq_len(pop), i), 1)
    if (runif(1) >= AP) {
      new <- integer(k)
      for (d in seq_len(k)) new[d] <- if (runif(1) < fl / (fl + 1)) mem[[j]][sample(k, 1)] else pos[[i]][d]
    } else new <- sample(n, k)
    new <- decode_medoid(new, k); s <- f$fitness(new); pos[[i]] <- new
    if (s > fmem[i]) { mem[[i]] <- new; fmem[i] <- s }
  }
}
run_GWO <- function(f, k, pop = 30, maxiter = 100000L, a_iters = 100L) {
  P <- matrix(runif(pop * k, 1, n), nrow = pop); fit <- apply(P, 1, f$fitness)
  for (it in seq_len(maxiter)) {
    o <- order(fit, decreasing = TRUE); L <- list(P[o[1], ], P[o[2], ], P[o[3], ])
    a <- max(0, 2 - it * (2 / a_iters))
    for (i in seq_len(pop)) {
      for (d in seq_len(k)) {
        Xl <- vapply(1:3, function(l) { A <- 2 * a * runif(1) - a; C <- 2 * runif(1)
          L[[l]][d] - A * abs(C * L[[l]][d] - P[i, d]) }, numeric(1))
        P[i, d] <- mean(Xl)
      }
      P[i, ] <- pmin(pmax(P[i, ], 1), n)
    }
    fit <- apply(P, 1, f$fitness)
  }
}
run_RS <- function(f, k, ...) repeat f$fitness(sample(n, k))
RUN <- list(GA = run_GA, PSO = run_PSO, CSA = run_CSA, GWO = run_GWO, RS = run_RS)
# nominal (v9) settings: GA 50 x 200, PSO 30 x 200, CSA 30 x 100, GWO 30 x 100
NOMINAL <- list(GA = list(pop = 50, maxiter = 200L), PSO = list(pop = 30, maxiter = 200L),
                CSA = list(pop = 30, maxiter = 100L), GWO = list(pop = 30, maxiter = 100L))

run_one <- function(alg, k, seed, budget = BUDGET, pop = NULL, maxiter = 100000L, target = NA_real_) {
  set.seed(seed)
  f <- new_counter(k, budget, target)
  args <- list(f = f, k = k, maxiter = maxiter); if (!is.null(pop)) args$pop <- pop
  t0 <- proc.time()[["elapsed"]]
  tryCatch(do.call(RUN[[alg]], args), FitnessBudgetReached = function(e) invisible(NULL))
  el <- proc.time()[["elapsed"]] - t0
  if (budget < 1e6 && f$calls != budget) stop(sprintf("%s k=%d seed=%d: %d calls", alg, k, seed, f$calls))
  list(best = f$best, med = f$best_med, calls = f$calls, hit = f$hit,
       secs = if (is.na(f$t_budget)) el else f$t_budget, trace = f$trace[seq_len(f$calls)])
}


exhaustive <- function(k) {
  cm <- utils::combn(n, k)
  sc <- vapply(seq_len(ncol(cm)), function(j) sil_med(cm[, j]), numeric(1))
  b <- which.max(sc)
  list(best = sc[b], med = cm[, b], n_opt = sum(sc > sc[b] - 1e-12), lab = labels_of(cm[, b]))
}
d23 <- read_eqi(PATH_CSV_2023)
stopifnot(length(d23$prov) == 38L)

# ---- S1. REPLICATION ON EQI 2022 ----------------------------------------------------
{
  d22 <- read_eqi(PATH_CSV_2022)
  say("EQI 2022: %d provinces", length(d22$prov)); print(d22$prov)
  # 2023 reference partitions (exhaustive)
  set_data(d23); g23 <- lapply(2:3, exhaustive); names(g23) <- c("2", "3")
  lab23 <- lapply(g23, `[[`, "lab"); prov23 <- provinsi
  # 2022
  set_data(d22)
  D2 <- mahalanobis(as.matrix(d22$X), colMeans(d22$X), cov(d22$X))
  mah <- data.frame(Province = provinsi, D2 = round(D2, 4), Outlier = D2 > qchisq(0.999, 4))[order(-D2), ]
  wcsv(mah, "supp_S1_2022_mahalanobis.csv")
  pam <- sapply(2:6, function(k) sil_med(cluster::pam(DIST, k = k, diss = TRUE)$id.med))
  ex <- lapply(EXH_K, exhaustive); names(ex) <- EXH_K
  glob <- sapply(ex, `[[`, "best")
  part_txt <- function(lab) { sz <- table(lab); paste(sapply(names(sz)[order(sz)], function(c)
    if (sz[[c]] <= 3) paste(provinsi[lab == as.integer(c)], collapse = " + ") else sprintf("%d provinces", sz[[c]])), collapse = " | ") }
  tabS1 <- data.frame(k = EXH_K, PAM_SC = pam[EXH_K - 1], global_SC = glob,
                      n_optimal_sets = sapply(ex, `[[`, "n_opt"),
                      sizes = sapply(ex, function(e) paste(sort(tabulate(e$lab)), collapse = "-")),
                      partition = sapply(ex, function(e) part_txt(e$lab)))
  wcsv(tabS1, "supp_S1_2022_exhaustive_vs_PAM.csv"); print(tabS1, row.names = FALSE)
  prof <- aggregate(d22$X, by = list(Cluster = ex[["3"]]$lab), FUN = mean); prof$n <- as.integer(table(ex[["3"]]$lab))
  wcsv(cbind(prof[1], round(prof[-1], 2)), "supp_S1_2022_profile_k3.csv")
  # agreement with 2023 on provinces with unchanged boundaries
  # Papua and Papua Barat lost territory to the four new provinces in 2022, so they are excluded
  common <- setdiff(intersect(provinsi, prov23), c("PAPUA", "PAPUA BARAT"))
  ariS <- sapply(c("2", "3"), function(k) ari(ex[[k]]$lab[match(common, provinsi)], lab23[[k]][match(common, prov23)]))
  wcsv(data.frame(k = 2:3, n_common = length(common), ARI_2022_vs_2023 = ariS), "supp_S1_ARI_2022_vs_2023.csv")
  # equal-budget benchmark 2022
  eb <- do.call(rbind, lapply(EB_K, function(k) do.call(rbind, lapply(seq_len(N_RUN), function(r)
    do.call(rbind, lapply(c("GA", "PSO", "CSA", "GWO", "RS"), function(a) {
      o <- run_one(a, k, SEED + r, target = glob[[as.character(k)]])
      data.frame(k = k, run = r, algorithm = a, best_SC = round(o$best, 10), success = !is.na(o$hit), evals_to_opt = o$hit)
    }))))))
  wcsv(eb, "supp_S1_2022_equal_budget_runs.csv")
  sr <- aggregate(cbind(success, best_SC) ~ algorithm + k, data = eb, FUN = mean)
  names(sr) <- c("algorithm", "k", "success_rate", "mean_SC")
  sr <- sr[order(sr$k, match(sr$algorithm, c("GA", "PSO", "CSA", "GWO", "RS"))), ]
  wcsv(sr, "supp_S1_2022_equal_budget_summary.csv"); print(sr, row.names = FALSE)
  # tests (rounded values; Holm-adjusted Wilcoxon of CSA against each optimizer)
  wx <- do.call(rbind, lapply(EB_K, function(k) do.call(rbind, lapply(c("GA", "PSO", "GWO"), function(a) {
    x <- eb$best_SC[eb$k == k & eb$algorithm == "CSA"]; y <- eb$best_SC[eb$k == k & eb$algorithm == a]
    p <- if (all(x == y)) 1 else suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE, correct = FALSE)$p.value)
    data.frame(k = k, comparison = paste("CSA vs", a), mean_diff = mean(x - y), p_raw = p) }))))
  wx$p_holm <- ave(wx$p_raw, wx$k, FUN = function(p) p.adjust(p, "holm"))
  wcsv(wx, "supp_S1_2022_wilcoxon_CSA.csv"); print(wx, row.names = FALSE)
  for (k in EXH_K) { tok(sprintf("S22_GLOBAL_K%d", k), f4(glob[[as.character(k)]])); tok(sprintf("S22_PAM_K%d", k), f4(pam[k - 1])) }
  tok("S22_N", n); tok("S22_K2_PART", tabS1$partition[tabS1$k == 2]); tok("S22_K3_PART", tabS1$partition[tabS1$k == 3])
  tok("S22_ARI_K2", f3(ariS[["2"]])); tok("S22_ARI_K3", f3(ariS[["3"]])); tok("S22_NCOMMON", length(common))
  for (i in seq_len(nrow(sr))) tok(sprintf("S22_SR_%s_K%d", sr$algorithm[i], sr$k[i]), fpct(sr$success_rate[i]))
  tok("S22_MAHAL_OUT", if (any(mah$Outlier)) paste(mah$Province[mah$Outlier], collapse = ", ") else "none")
}


ph <- data.frame(key = ls(TOK), value = unlist(mget(ls(TOK), envir = TOK)))
write.csv(ph, file.path(DIR_OUT, "supp_placeholders_S1.csv"), row.names = FALSE, fileEncoding = "UTF-8")
writeLines(capture.output(sessionInfo()), file.path(DIR_OUT, "supp_S1_sessionInfo.txt"))
cat("\n=== ARI 2022 vs 2023 ===\n"); print(read.csv(file.path(DIR_OUT, "supp_S1_ARI_2022_vs_2023.csv")))
say("DONE in %.1f min -> %s", (proc.time()[["elapsed"]] - T0) / 60, normalizePath(DIR_OUT))
